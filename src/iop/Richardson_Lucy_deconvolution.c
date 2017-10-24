/*
    This file is part of darktable,
    copyright (c) 2009--2010 johannes hanika.
    copyright (c) 2014-2016 Roman Lebedev.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "bauhaus/bauhaus.h"
#include "common/richardson_lucy_deconvolution.h"

#ifdef HAVE_FFTW3_OMP
#include <fftw3.h>
#endif


DT_MODULE_INTROSPECTION(1, dt_iop_rlucydec_params_t)


typedef struct dt_iop_rlucydec_params_t
{
  int blur_type;
  int quality;
  float artifacts_damping;
  int deblur_strength;
  int blur_width;
  int blur_strength;
  int refine;
  float mask[4];
  int use_mask;
  float backvsmask_ratio;
} dt_iop_rlucydec_params_t;

typedef struct dt_iop_rlucydec_gui_data_t
{
	int drawing;
	float pt_down_x, pt_down_y, pt_last_x, pt_last_y;

  GtkWidget *bt_draw_focus_zone;
  GtkWidget *cmb_blur_type;
  GtkWidget *sl_quality;
  GtkWidget *sl_artifacts_damping;
  GtkWidget *sl_deblur_strength;
  GtkWidget *sl_blur_width;
  GtkWidget *sl_blur_strength;
  GtkWidget *chk_refine;
  GtkWidget *chk_use_mask;
  GtkWidget *sl_backvsmask_ratio;
} dt_iop_rlucydec_gui_data_t;

typedef struct dt_iop_rlucydec_params_t dt_iop_rlucydec_data_t;

const char *name()
{
  return _("Richardson-Lucy deconvolution");
}

int groups()
{
  return IOP_GROUP_TONE;
}

int flags()
{
  return IOP_FLAGS_SUPPORTS_BLENDING;
}

#define PREAMBLE                                        \
  cairo_save (cr);                                      \
  const gint s = MIN (w, h);                            \
  cairo_translate (cr, x + (w / 2.0) - (s / 2.0),       \
                   y + (h / 2.0) - (s / 2.0));          \
  cairo_scale (cr, s, s);                               \
  cairo_push_group (cr);                                \
  cairo_set_source_rgba (cr, 1.0, 1.0, 1.0, 1.0);       \
  cairo_set_line_cap (cr, CAIRO_LINE_CAP_ROUND);        \
  cairo_set_line_width (cr, 0.1);

#define POSTAMBLE                                               \
  cairo_pop_group_to_source (cr);                               \
  cairo_paint_with_alpha (cr, flags & CPF_ACTIVE ? 1.0 : 0.5);  \
  cairo_restore (cr);

static void _cairo_paint_tool_draw_focus_zone(cairo_t *cr,
                                          const gint x, const gint y, const gint w, const gint h, const gint flags)
{
  PREAMBLE;
  
  cairo_set_line_width(cr, 0.1);
  cairo_rectangle(cr, 0., 0., 1., 1.);
  cairo_stroke(cr);

  POSTAMBLE;
}


static void rlucy_blur_type_callback(GtkComboBox *combo, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  
  dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;

  p->blur_type = dt_bauhaus_combobox_get((GtkWidget *)combo);
  
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void rlucy_quality_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;
  
  p->quality = dt_bauhaus_slider_get(slider);
  
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void rlucy_artifacts_damping_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;
  
  p->artifacts_damping = dt_bauhaus_slider_get(slider);
  
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void rlucy_deblur_strength_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;
  
  p->deblur_strength = dt_bauhaus_slider_get(slider);
  
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void rlucy_blur_width_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;
  
  p->blur_width = dt_bauhaus_slider_get(slider);
  
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void rlucy_blur_strength_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;
  
  p->blur_strength = dt_bauhaus_slider_get(slider);
  
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void rlucy_chk_refine_callback(GtkWidget *widget, dt_iop_module_t *self)
{
	if(self->dt->gui->reset) return;
	dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;
	
	p->refine = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget));
	
	dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void rlucy_chk_use_mask_callback(GtkWidget *widget, dt_iop_module_t *self)
{
	if(self->dt->gui->reset) return;
	dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;
	
	p->use_mask = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget));
	
	dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void rlucy_backvsmask_ratio_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;
  
  p->backvsmask_ratio = dt_bauhaus_slider_get(slider);
  
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void rlucy_draw_focus_zone_callback(GtkToggleButton *togglebutton, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;

  dt_iop_request_focus(self);
  
  dt_control_queue_redraw_center();
  
}

static float get_zoom_scale (dt_develop_t *develop)
{
  const dt_dev_zoom_t zoom = dt_control_get_dev_zoom ();
  const int closeup = dt_control_get_dev_closeup ();
  return dt_dev_get_zoom_scale (develop, zoom, closeup ? 2 : 1, 1);
}

static void get_point_scale(struct dt_iop_module_t *module, float x, float y, float *pt_x, float *pt_y, float *scale)
{
  float pzx, pzy;
  dt_dev_get_pointer_zoom_pos(darktable.develop, x, y, &pzx, &pzy);
  pzx += 0.5f;
  pzy += 0.5f;
  float wd = darktable.develop->preview_pipe->backbuf_width;
  float ht = darktable.develop->preview_pipe->backbuf_height;
  float pts[2] = { pzx * wd, pzy * ht };
  dt_dev_distort_backtransform_plus(darktable.develop, darktable.develop->preview_pipe,
                                    module->priority + 1, 9999999, pts, 1);
  dt_dev_distort_backtransform_plus(darktable.develop, darktable.develop->preview_pipe,
                                    0, module->priority - 1, pts, 1);
  float nx = pts[0] / darktable.develop->preview_pipe->iwidth;
  float ny = pts[1] / darktable.develop->preview_pipe->iheight;

  *scale = darktable.develop->preview_pipe->iscale / get_zoom_scale(module->dev);
  *pt_x = (nx * darktable.develop->pipe->iwidth);
  *pt_y = (ny * darktable.develop->pipe->iheight);
}

static void distord_rec(struct dt_iop_module_t *self, dt_develop_t *dev, dt_dev_pixelpipe_t *pipe, 
    float from_scale, float to_scale, int pmin, int pmax, 
    float *pt_x, float *pt_y, float *pt_w, float *pt_h)
{
  float x = *pt_x;
  float y = *pt_y;
  float w = *pt_w;
  float h = *pt_h;
  
  float buffer[4] = {x / from_scale, y / from_scale, (x+w) / from_scale, (y+h) / from_scale};
  int len = 2;
  
  if (pmin < self->priority && pmax > self->priority)
  {
    dt_dev_distort_transform_plus (dev, pipe, pmin, self->priority - 1, buffer, len);
    dt_dev_distort_transform_plus (dev, pipe, self->priority + 1, pmax, buffer, len);
  }
  else
    dt_dev_distort_transform_plus (dev, pipe, pmin, pmax, buffer, len);

  *pt_x = buffer[0] * to_scale;
  *pt_y = buffer[1] * to_scale;
  *pt_w = buffer[2] * to_scale - *pt_x;
  *pt_h = buffer[3] * to_scale - *pt_y;

}

void gui_focus(struct dt_iop_module_t *self, gboolean in)
{
  dt_iop_rlucydec_gui_data_t *g = (dt_iop_rlucydec_gui_data_t *)self->gui_data;

  if (!in)
  {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->bt_draw_focus_zone), FALSE);
  }
}

void commit_params(struct dt_iop_module_t *self, dt_iop_params_t *params, dt_dev_pixelpipe_t *pipe,
    								dt_dev_pixelpipe_iop_t *piece)
{
  memcpy(piece->data, params, sizeof(dt_iop_rlucydec_params_t));
}

void init_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  piece->data = malloc(sizeof(dt_iop_rlucydec_data_t));
  self->commit_params(self, self->default_params, pipe, piece);
}

void cleanup_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  free(piece->data);
  piece->data = NULL;
}

void gui_update(struct dt_iop_module_t *self)
{
  dt_iop_rlucydec_gui_data_t *g = (dt_iop_rlucydec_gui_data_t *)self->gui_data;
  dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;

  dt_bauhaus_combobox_set(g->cmb_blur_type, p->blur_type);
  
  dt_bauhaus_slider_set(g->sl_quality, p->quality);
  dt_bauhaus_slider_set(g->sl_artifacts_damping, p->artifacts_damping);
  dt_bauhaus_slider_set(g->sl_deblur_strength, p->deblur_strength);
  dt_bauhaus_slider_set(g->sl_blur_width, p->blur_width);
  dt_bauhaus_slider_set(g->sl_blur_strength, p->blur_strength);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->chk_refine), p->refine);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->chk_use_mask), p->use_mask);
  dt_bauhaus_slider_set(g->sl_backvsmask_ratio, p->backvsmask_ratio);

}

void init(dt_iop_module_t *module)
{
	module->data = NULL;
  module->params = calloc(1, sizeof(dt_iop_rlucydec_params_t));
  module->default_params = calloc(1, sizeof(dt_iop_rlucydec_params_t));
  module->default_enabled = 0;
//  module->priority = 164; // module order created by iop_dependencies.py, do not edit! // from exposure
  module->priority = 582; // module order created by iop_dependencies.py, do not edit! // from bilat
  module->params_size = sizeof(dt_iop_rlucydec_params_t);
  module->gui_data = NULL;
  
  dt_iop_rlucydec_params_t tmp = {0};
  
  tmp.blur_type = blur_type_kaiser;
  tmp.quality = 50;
  tmp.artifacts_damping = 0.00003f;
  tmp.deblur_strength = 20;
  tmp.blur_width = 11;
  tmp.blur_strength = 8;
  tmp.refine = 0;
  tmp.backvsmask_ratio = 0.f;

  memcpy(module->params, &tmp, sizeof(dt_iop_rlucydec_params_t));
  memcpy(module->default_params, &tmp, sizeof(dt_iop_rlucydec_params_t));

}

void cleanup(dt_iop_module_t *module)
{
  free(module->params);
  module->params = NULL;
}

void gui_init(struct dt_iop_module_t *self)
{
  const int bs = DT_PIXEL_APPLY_DPI(14);
  
  self->gui_data = malloc(sizeof(dt_iop_rlucydec_gui_data_t));
  dt_iop_rlucydec_gui_data_t *g = (dt_iop_rlucydec_gui_data_t *)self->gui_data;
  dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;

  self->widget = GTK_WIDGET(gtk_box_new(GTK_ORIENTATION_VERTICAL, DT_BAUHAUS_SPACE));

  g->drawing = 0;
  g->pt_down_x = g->pt_down_y = g->pt_last_x = g->pt_last_y = 0.f;

  g->cmb_blur_type = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->cmb_blur_type, NULL, _("blur type"));
  dt_bauhaus_combobox_add(g->cmb_blur_type, _("auto"));
  dt_bauhaus_combobox_add(g->cmb_blur_type, _("gaussian"));
  dt_bauhaus_combobox_add(g->cmb_blur_type, _("kaiser"));
  g_object_set(g->cmb_blur_type, "tooltip-text", _("blur type."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->cmb_blur_type), "value-changed", G_CALLBACK(rlucy_blur_type_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), GTK_WIDGET(g->cmb_blur_type), TRUE, TRUE, 0);

  g->sl_quality = dt_bauhaus_slider_new_with_range(self, 1, 100, 1, p->quality, 0);
  dt_bauhaus_widget_set_label(g->sl_quality, _("quality"), _("quality"));
//  dt_bauhaus_slider_set_format(g->sl_quality, "%.01f");
  g_object_set(g->sl_quality, "tooltip-text", _("quality of the refining\n"
      "(ie the total number of iterations to compute)."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_quality), "value-changed", G_CALLBACK(rlucy_quality_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_quality, TRUE, TRUE, 0);

  g->sl_artifacts_damping = dt_bauhaus_slider_new_with_range(self, 0.00003, 1.0, 0.00001, p->artifacts_damping, 6);
  dt_bauhaus_widget_set_label(g->sl_artifacts_damping, _("artifacts damping"), _("artifacts damping"));
//  dt_bauhaus_slider_set_format(g->sl_artifacts_damping, "%.01f");
  g_object_set(g->sl_artifacts_damping, "tooltip-text", _("the noise and artifacts reduction factor lambda.\n"
      "increase it if smudges, noise or ringing appear."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_artifacts_damping), "value-changed", G_CALLBACK(rlucy_artifacts_damping_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_artifacts_damping, TRUE, TRUE, 0);

  g->sl_blur_width = dt_bauhaus_slider_new_with_range(self, 3, 11.0, 1, p->blur_width, 0);
  dt_bauhaus_widget_set_label(g->sl_blur_width, _("blur width"), _("blur width"));
//  dt_bauhaus_slider_set_format(g->sl_blur_width, "%.01f");
  g_object_set(g->sl_blur_width, "tooltip-text", _("width of the blur in px."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_blur_width), "value-changed", G_CALLBACK(rlucy_blur_width_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_blur_width, TRUE, TRUE, 0);

  g->sl_blur_strength = dt_bauhaus_slider_new_with_range(self, 1, 50.0, 1, p->blur_strength, 0);
  dt_bauhaus_widget_set_label(g->sl_blur_strength, _("blur strength"), _("blur strength"));
//  dt_bauhaus_slider_set_format(g->sl_blur_strength, "%.01f");
  g_object_set(g->sl_blur_strength, "tooltip-text", _("standard deviation of the blur kernel."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_blur_strength), "value-changed", G_CALLBACK(rlucy_blur_strength_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_blur_strength, TRUE, TRUE, 0);

	g->chk_refine = gtk_check_button_new_with_label("refine");
	g_object_set(G_OBJECT(g->chk_refine), "tooltip-text", _("refine"), (char *)NULL);
	g_signal_connect(G_OBJECT(g->chk_refine), "toggled", G_CALLBACK(rlucy_chk_refine_callback), self);

	gtk_box_pack_start(GTK_BOX(self->widget), g->chk_refine, TRUE, TRUE, 0);
	
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->chk_refine), p->refine);
	
  g->sl_deblur_strength = dt_bauhaus_slider_new_with_range(self, 0, 100.0, 1, p->deblur_strength, 0);
  dt_bauhaus_widget_set_label(g->sl_deblur_strength, _("deblur strength"), _("deblur strength"));
//  dt_bauhaus_slider_set_format(g->sl_deblur_strength, "%.01f");
  g_object_set(g->sl_deblur_strength, "tooltip-text", _("number of debluring iterations to perform."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_deblur_strength), "value-changed", G_CALLBACK(rlucy_deblur_strength_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_deblur_strength, TRUE, TRUE, 0);

  GtkWidget *hbox_mask = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);

  g->bt_draw_focus_zone = dtgtk_togglebutton_new(_cairo_paint_tool_draw_focus_zone, CPF_STYLE_FLAT | CPF_DO_NOT_USE_BORDER);
  g_object_set(G_OBJECT(g->bt_draw_focus_zone), "tooltip-text", _("activates focus zone drawing"), (char *)NULL);
  g_signal_connect(G_OBJECT(g->bt_draw_focus_zone), "toggled", G_CALLBACK(rlucy_draw_focus_zone_callback), self);
  gtk_widget_set_size_request(GTK_WIDGET(g->bt_draw_focus_zone), bs, bs);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->bt_draw_focus_zone), FALSE);
  gtk_box_pack_end(GTK_BOX(hbox_mask), g->bt_draw_focus_zone, FALSE, FALSE, 0);

  g->chk_use_mask = gtk_check_button_new_with_label("use mask");
  g_object_set(G_OBJECT(g->chk_use_mask), "tooltip-text", _("use mask"), (char *)NULL);
  g_signal_connect(G_OBJECT(g->chk_use_mask), "toggled", G_CALLBACK(rlucy_chk_use_mask_callback), self);

  gtk_box_pack_start(GTK_BOX(hbox_mask), g->chk_use_mask, TRUE, TRUE, 0);
  
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->chk_use_mask), p->use_mask);

  gtk_box_pack_start(GTK_BOX(self->widget), hbox_mask, TRUE, TRUE, 0);

	g->sl_backvsmask_ratio = dt_bauhaus_slider_new_with_range(self, 0.0, 1.0, 0.01, p->backvsmask_ratio, 3);
  dt_bauhaus_widget_set_label(g->sl_backvsmask_ratio, _("back vs. mask ratio"), _("back vs. mask ratio"));
//  dt_bauhaus_slider_set_format(g->sl_backvsmask_ratio, "%.01f");
  g_object_set(g->sl_backvsmask_ratio, "tooltip-text", _("when a mask is used, the ratio  of weights of the whole image / the masked zone.\n"
         "0 means only the masked zone is used,\n"
      "1 means the masked zone is ignored and only the whole image is taken.\n"
      "0 runs faster, 1 runs much slower."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_backvsmask_ratio), "value-changed", G_CALLBACK(rlucy_backvsmask_ratio_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_backvsmask_ratio, TRUE, TRUE, 0);

}

void init_global(dt_iop_module_so_t *module)
{
  module->data = NULL;
  
  // FIXME: this should go into dt startup
#ifdef HAVE_FFTW3_OMP
#ifdef _OPENMP
  
	fftwf_init_threads();
	
#endif
#endif
  
}

void cleanup_global(dt_iop_module_so_t *module)
{
  module->data = NULL;
  
  // FIXME: this should go into dt cleanup
#ifdef HAVE_FFTW3_OMP
#ifdef _OPENMP
  
  fftwf_cleanup_threads();
	
#endif
#endif
  
}

void gui_cleanup(struct dt_iop_module_t *self)
{
  free(self->gui_data);
  self->gui_data = NULL;
}

void gui_post_expose(struct dt_iop_module_t *self, cairo_t *cr, int32_t width, int32_t height,
                     int32_t pointerx, int32_t pointery)
{
  dt_develop_t *dev = self->dev;
  dt_iop_rlucydec_gui_data_t *g = (dt_iop_rlucydec_gui_data_t *)self->gui_data;
  dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;
  if (!g) return;
  if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(g->bt_draw_focus_zone)))  return;

  const float bb_width = dev->preview_pipe->backbuf_width;
  const float bb_height = dev->preview_pipe->backbuf_height;

  const float scale = MAX (bb_width, bb_height);
  if (bb_width < 1.0 || bb_height < 1.0)
    return;

  // You're not supposed to understand this
  const float zoom_x = dt_control_get_dev_zoom_x ();
  const float zoom_y = dt_control_get_dev_zoom_y ();
  const float zoom_scale = get_zoom_scale (dev);

  double dashes = DT_PIXEL_APPLY_DPI(5.0) * (1.0 / (scale * zoom_scale));

  float x, y, w, h;

  if (g->drawing)
  {
    x = MIN(g->pt_down_x, g->pt_last_x);
    y = MIN(g->pt_down_y, g->pt_last_y);
    w = fabs(g->pt_down_x - g->pt_last_x);
    h = fabs(g->pt_down_y - g->pt_last_y);
  }
  else
  {
    x = p->mask[0];
    y = p->mask[1];
    w = p->mask[2];
    h = p->mask[3];

//  printf("gui_post_expose x=%f, y=%f, w=%f, h=%f\n", x, y, w, h);
//  printf("gui_post_expose x=%f, y=%f, w=%f, h=%f\n", p->mask[0], p->mask[1], p->mask[2], p->mask[3]);
  }

  const float iscale = dev->preview_pipe->iscale;

  float from_scale = iscale;
  float to_scale = 1.0 / scale;

  distord_rec(self, dev, dev->preview_pipe, from_scale, to_scale, 0, 9999999, &x, &y, &w, &h);
  if (w <= 0. || h <= 0.) return;

  cairo_save(cr);

  // setup CAIRO coordinate system
  cairo_translate (cr, 0.5 * width, 0.5 * height); // origin @ center of view
  cairo_scale     (cr, zoom_scale, zoom_scale);    // the zoom
  cairo_translate (cr, -bb_width * (0.5 + zoom_x), -bb_height * (0.5 + zoom_y));
  cairo_scale (cr, scale, scale);

  cairo_set_line_width(cr, dashes / 2.0);
  cairo_rectangle(cr, x, y, w, h);
  cairo_set_source_rgb(cr, .7, .7, .7);
  cairo_stroke(cr);

  cairo_restore(cr);

//	printf("gui_post_expose x=%f, y=%f, w=%f, h=%f\n", x, y, w, h);

}

int mouse_moved(struct dt_iop_module_t *self, double x, double y, double pressure, int which)
{
	dt_iop_rlucydec_gui_data_t *g = (dt_iop_rlucydec_gui_data_t *)self->gui_data;
	if (g->drawing)
	{
		float pt_x, pt_y;
	  float scale;

	  get_point_scale(self, x, y, &pt_x, &pt_y, &scale);

	  g->pt_last_x = pt_x;
	  g->pt_last_y = pt_y;

		dt_control_queue_redraw_center();
		
		return 1;
	}

	return 0;
}

int button_pressed(struct dt_iop_module_t *self, double x, double y, double pressure, int which, int type,
                   uint32_t state)
{
	dt_iop_rlucydec_gui_data_t *g = (dt_iop_rlucydec_gui_data_t *)self->gui_data;
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(g->bt_draw_focus_zone)))
	{
		float pt_x, pt_y;
	  float scale;
	  
	  get_point_scale(self, x, y, &pt_x, &pt_y, &scale);

	  g->pt_down_x = g->pt_last_x = pt_x;
	  g->pt_down_y = g->pt_last_y = pt_y;

		g->drawing = 1;
		
		return 1;
	}
	
	return 0;
}

int button_released(struct dt_iop_module_t *self, double x, double y, int which, uint32_t state)
{
	dt_iop_rlucydec_gui_data_t *g = (dt_iop_rlucydec_gui_data_t *)self->gui_data;
	dt_iop_rlucydec_params_t *p = (dt_iop_rlucydec_params_t *)self->params;
	if (g->drawing)
	{
		float pt_x, pt_y;
	  float scale;

	  get_point_scale(self, x, y, &pt_x, &pt_y, &scale);

	  g->pt_last_x = pt_x;
	  g->pt_last_y = pt_y;

	  //-----------------------------

	  pt_x = MIN(g->pt_down_x, g->pt_last_x);
	  pt_y = MIN(g->pt_down_y, g->pt_last_y);
	  float pt_w = fabs(g->pt_down_x - g->pt_last_x);
	  float pt_h = fabs(g->pt_down_y - g->pt_last_y);

    float from_scale = 1.0; //iscale;
    float to_scale = 1.0 ;/// scale;

//	  printf("button_released x=%f, y=%f, w=%f, h=%f\n", pt_x, pt_y, pt_w, pt_h);
	  distord_rec(self, self->dev, self->dev->preview_pipe, from_scale, to_scale, 0, 9999999, &pt_x, &pt_y, &pt_w, &pt_h);
//	  printf("button_released x=%f, y=%f, w=%f, h=%f\n\n", pt_x, pt_y, pt_w, pt_h);
	  
	  p->mask[0] = pt_x;
	  p->mask[1] = pt_y;
	  p->mask[2] = pt_w;
	  p->mask[3] = pt_h;
	  
	  g->pt_down_x = g->pt_last_x = 0.;
	  g->pt_down_y = g->pt_last_y = 0.;

	  dt_dev_add_history_item (darktable.develop, self, TRUE);
	  
	  //--------------------------------
	  
	  g->drawing = 0;
		
		return 1;
	}
	
	return 0;
}

static void check_nan(const float* im, const int size)
{
	int i_nan = 0;
	
	for (int i = 0; i < size; i++) if ( isnan(im[i]) ) i_nan++;
	
	if (i_nan > 0) printf("Richardson-Lucy deconvolution nan: %i\n", i_nan);
}

void process_internal(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid, void *const ovoid,
    const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out, const int use_sse)
{
  const dt_iop_rlucydec_data_t *const p = (const dt_iop_rlucydec_data_t *const)piece->data;
  dt_develop_t *dev = self->dev;

  const dt_iop_colorspace_type_t cst = dt_iop_module_colorspace(self);
  const int ch = piece->colors;

  float *im_tmp = NULL;
  float *kernel = NULL;

  int blur_type = p->blur_type;
  int quality = p->quality;
  float artifacts_damping = p->artifacts_damping;
  int deblur_strength = p->deblur_strength;
  int blur_width = p->blur_width;
  int blur_strength = p->blur_strength;
  int refine = p->refine;
  int *mask = NULL;
  float backvsmask_ratio = p->backvsmask_ratio;

  float pt_x = p->mask[0];
  float pt_y = p->mask[1];
  float pt_w = p->mask[2];
  float pt_h = p->mask[3];
  int rl_mask[4] = {0};

  if (p->use_mask)
  {
    float from_scale = 1.0;
    float to_scale = roi_in->scale / piece->iscale;

    distord_rec(self, dev, piece->pipe, from_scale, to_scale, 0, self->priority - 1, &pt_x, &pt_y, &pt_w, &pt_h);

    //  printf("process_internal pt_x=%f, pt_y=%f, pt_w=%f, pt_h=%f\n", pt_x, pt_y, pt_w, pt_h);
    //  printf("process_internal roi_in->x=%i, roi_in->y=%i, roi_in->width=%i, roi_in->height=%i\n", roi_in->x, roi_in->y, roi_in->width, roi_in->height);

    pt_x -= (float)roi_in->x;
    pt_y -= (float)roi_in->y;

    //  printf("process_internal x=%f, y=%f, w=%f, h=%f\n", pt_x, pt_y, pt_w, pt_h);
    if (pt_w > 0. && pt_h > 0.)
    {
      rl_mask[0] = MAX(pt_x, 0.);
      rl_mask[1] = MAX(pt_y, 0.);
      if (pt_x < 0.)
        rl_mask[2] = pt_w + pt_x;
      else
        rl_mask[2] = pt_w;
      if (pt_y < 0.)
        rl_mask[3] = pt_h + pt_y;
      else
        rl_mask[3] = pt_h;
      rl_mask[2] = MIN(rl_mask[2], roi_in->width-rl_mask[0]);
      rl_mask[3] = MIN(rl_mask[3], roi_in->height-rl_mask[1]);

      if (rl_mask[2] > 0 && rl_mask[3] > 0)
      {
        mask = rl_mask;
        //  		printf("process_internal x=%i, y=%i, w=%i, h=%i\n", rl_mask[0], rl_mask[1], rl_mask[2], rl_mask[3]);
      }
    }
  }

  if (cst == iop_cs_rgb)
  {
    kernel = (float*)dt_alloc_align(64, blur_width * blur_width * ch * sizeof(float));
    if (im_tmp && kernel)
    {
      memcpy(ovoid, ivoid, roi_in->width * roi_in->height * ch * sizeof(float));

      richardson_lucy_build_kernel(im_tmp, roi_out->width, roi_out->height, 1, 
          kernel, blur_width, blur_strength, 
          blur_type, quality, artifacts_damping);
      
      richardson_lucy(im_tmp, roi_out->width, roi_out->height, 1, 
          kernel, blur_width, 
          quality, artifacts_damping, deblur_strength, 
          refine, mask, backvsmask_ratio);
    }
  }
  else if (cst == iop_cs_Lab)
  {
    memcpy(ovoid, ivoid, roi_in->width * roi_in->height * ch * sizeof(float));

    im_tmp = (float*)dt_alloc_align(64, roi_out->width * roi_out->height * sizeof(float));
    kernel = (float*)dt_alloc_align(64, blur_width * blur_width * sizeof(float));
    if (im_tmp && kernel)
    {
      float *im_src = (float*)ivoid;

      for (int i = 0; i < roi_out->width * roi_out->height; i++)
      {
        im_tmp[i] = im_src[i*ch];
      }

      richardson_lucy_build_kernel(im_tmp, roi_out->width, roi_out->height, 1, 
          kernel, blur_width, blur_strength, 
          blur_type, quality, artifacts_damping);
      
      richardson_lucy(im_tmp, roi_out->width, roi_out->height, 1, 
          kernel, blur_width, 
          quality, artifacts_damping, deblur_strength, 
          refine, mask, backvsmask_ratio);
      
      float *im_dest = (float*)ovoid;

      for (int i = 0; i < roi_out->width * roi_out->height; i++)
      {
        im_dest[i*ch] = im_tmp[i];
      }
    }
  }
  else
  {
    memcpy(ovoid, ivoid, roi_in->width * roi_in->height * ch * sizeof(float));
  }

  check_nan((float*)ovoid, roi_out->width*roi_out->height*ch);

  if (im_tmp) dt_free_align(im_tmp);
  if (kernel) dt_free_align(kernel);

}

void process(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid, void *const ovoid,
             const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
	process_internal(self, piece, ivoid, ovoid, roi_in, roi_out, 0);
}

#if defined(__SSE__x)
void process_sse2(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid, void *const ovoid,
    				const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
	process_internal(self, piece, ivoid, ovoid, roi_in, roi_out, 1);
}
#endif

// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
