/*
    This file is part of darktable,
    copyright (c) 2018 edgardo hoszowski.

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


/*
   based on G'MIC filter "Freaky B&W"
   http://gmic.eu/

#@gui Freaky B&amp;W : fx_freaky_bw, fx_freaky_bw_preview
#@gui : Strength (%) = float(90,0,100)
#@gui : Oddness (%) = float(20,0,100)
#@gui : Brightness (%) = float(0,-100,100)
#@gui : Contrast (%) = float(0,-100,100)
#@gui : Gamma (%) = float(0,-100,100)
#@gui : sep = separator(), Preview type = choice("Full","Forward horizontal","Forward vertical","Backward
horizontal","Backward vertical","Duplicate top","Duplicate left","Duplicate bottom","Duplicate right")
#@gui : sep = separator(), note = note("<small>Author: <i>David Tschumperl&#233;</i>.      Latest update:
<i>09/30/2015</i>.</small>")
fx_freaky_bw :
  -repeat $! -l[$>] -split_opacity -l[0]
    -to_rgb

    # Estimate gradient field of B&W result.
    --expand_xy 1,0 -channels. 0,4
    -f. ">if (c!=4,i,
           Rx = i(x+1,y,0,0) - i(x,y,0,0);
           Ry = i(x,y+1,0,0) - i(x,y,0,0);
           Rn = Rx^2 + Ry^2;
           Gx = i(x+1,y,0,1) - i(x,y,0,1);
           Gy = i(x,y+1,0,1) - i(x,y,0,1);
           Gn = Gx^2 + Gy^2;
           Bx = i(x+1,y,0,2) - i(x,y,0,2);
           By = i(x,y+1,0,2) - i(x,y,0,2);
           Bn = Bx^2 + By^2;
           n = 1e-5 + max(Rn,Gn,Bn)^"{$2%}";
           val = 0;
          if (Rn>=Gn && Rn>=Bn,
            i(x,y,0,3) = Rx/n; val=Ry/n,
          if (Gn>=Rn && Gn>=Bn,
            i(x,y,0,3) = Gx/n; val=Gy/n,
            i(x,y,0,3) = Bx/n; val=By/n));
          val
         )"
    -channels. 3,4
    -luminance[0] ia={0,ia}

    # Estimate laplacian of final image.
    -s. c
    -f.. "i - i(x-1,y,0,0)"
    -f. "i - i(x,y-1,0,0)"
    -+[-2,-1]

    # Reconstruct image from laplacian.
    -fft.
    100%,100%,1,1,'cos(2*x*pi/w)+cos(2*y*pi/h)' -*. 2 --. 4
    -=. 1 -/[-3,-2] . -rm.
    -=.. 0 -=. 0
    -ifft[-2,-1] -rm.
    -shrink_xy. 1 -+. $ia -n. 0,255

    # Merge result with original color image.
    -j[0] [1],0,0,0,0,{$1%} -rm.
    -adjust_colors ${3-5}
  -endl -a c -endl -done

fx_freaky_bw_preview :
  -gui_split_preview "-fx_freaky_bw $*",$-1

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "bauhaus/bauhaus.h"
#include "common/colorspaces_inline_conversions.h"
#include "develop/imageop_math.h"

#ifdef HAVE_FFTW3
#include <fftw3.h>
#endif

#ifdef HAVE_FFTW3_OMP
#include <fftw3.h>
#endif

DT_MODULE_INTROSPECTION(1, dt_iop_fbw_params_t)

#define FBW_PREVIEW_LVL_MIN -3.0f
#define FBW_PREVIEW_LVL_MAX 3.0f

#ifdef HAVE_FFTW3

typedef struct
{
  float *in_src;
  fftwf_complex *out_src;

  fftwf_plan plan_src;
  fftwf_plan plan_inv;

  int width_src;
  int height_src;

  int width_dest;
  int height_dest;

  int width_fft;
  int height_fft;
  int width_fft_complex;
  int height_fft_complex;

  dt_pthread_mutex_t *fftw3_lock;
} fbw_fft_t;

#else

typedef float fftwf_complex[2];

typedef struct
{
  float *in_src;
  fftwf_complex *out_src;

  int width_src;
  int height_src;

  int width_dest;
  int height_dest;

  int width_fft;
  int height_fft;
  int width_fft_complex;
  int height_fft_complex;
} fbw_fft_t;

#endif

typedef enum dt_iop_fbw_drag_types_t {
  dt_iop_fbw_lvlbar_drag_left = 1,
  dt_iop_fbw_lvlbar_drag_middle = 2,
  dt_iop_fbw_lvlbar_drag_right = 3
} dt_iop_fbw_drag_types_t;

typedef enum dt_iop_fbw_bw_methods_t { dt_iop_fbw_bw_mix = 0, dt_iop_fbw_bw_max = 1 } dt_iop_fbw_bw_methods_t;

typedef struct dt_iop_fbw_params_t
{
  int bw_method;
  float oddness;
  float red;
  float green;
  float blue;
  float levels[3];
} dt_iop_fbw_params_t;

typedef struct dt_iop_fbw_gui_data_t
{
  dt_pthread_mutex_t lock;
  uint64_t hash;

  float l_min;
  float l_max;
  float l_mid;

  gboolean is_dragging;

  int auto_levels;
  float levels[3];
  float lvlbar_mouse_x, lvlbar_mouse_y;

  GtkWidget *cmb_bw_method;
  GtkWidget *sl_oddness;
  GtkWidget *sl_red;
  GtkWidget *sl_green;
  GtkWidget *sl_blue;
  GtkWidget *vbox_rgb;
  GtkWidget *levels_bar;
  GtkWidget *bt_auto_levels;
} dt_iop_fbw_gui_data_t;

typedef struct dt_iop_fbw_params_t dt_iop_fbw_data_t;

typedef struct dt_iop_fbw_global_data_t
{
  dt_pthread_mutex_t fftw3_lock;
} dt_iop_fbw_global_data_t;

const char *name()
{
  return _("fbw_eh");
}

int groups()
{
  return IOP_GROUP_COLOR;
}

int flags()
{
  return IOP_FLAGS_SUPPORTS_BLENDING;
}

#define PREAMBLE                                                                                                  \
  cairo_save(cr);                                                                                                 \
  const gint s = MIN(w, h);                                                                                       \
  cairo_translate(cr, x + (w / 2.0) - (s / 2.0), y + (h / 2.0) - (s / 2.0));                                      \
  cairo_scale(cr, s, s);                                                                                          \
  cairo_push_group(cr);                                                                                           \
  cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);                                                                  \
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);                                                                   \
  cairo_set_line_width(cr, 0.1);

#define POSTAMBLE                                                                                                 \
  cairo_pop_group_to_source(cr);                                                                                  \
  cairo_paint_with_alpha(cr, flags &CPF_ACTIVE ? 1.0 : 0.5);                                                      \
  cairo_restore(cr);

static void _fbw_cairo_paint_auto_levels(cairo_t *cr, const gint x, const gint y, const gint w, const gint h,
                                         const gint flags, void *data)
{
  PREAMBLE;

  cairo_move_to(cr, .1, 0.3);
  cairo_line_to(cr, .1, 1.);
  cairo_stroke(cr);

  cairo_move_to(cr, .5, 0.1);
  cairo_line_to(cr, .5, 1.);
  cairo_stroke(cr);

  cairo_move_to(cr, .9, 0.3);
  cairo_line_to(cr, .9, 1.);
  cairo_stroke(cr);

  cairo_move_to(cr, 0., 1.0);
  cairo_line_to(cr, 1.0, 1.0);
  cairo_stroke(cr);

  POSTAMBLE;
}

#undef PREAMBLE
#undef POSTAMBLE

static void clamp_minmax(float levels_old[3], float levels_new[3])
{
  // left or right has changed
  if((levels_old[0] != levels_new[0] || levels_old[2] != levels_new[2]) && levels_old[1] == levels_new[1])
  {
    // if old left and right are the same just use the new values
    if(levels_old[2] != levels_old[0])
    {
      // set the new value but keep the middle proportional
      const float left = MAX(levels_new[0], FBW_PREVIEW_LVL_MIN);
      const float right = MIN(levels_new[2], FBW_PREVIEW_LVL_MAX);

      const float percentage = (levels_old[1] - levels_old[0]) / (levels_old[2] - levels_old[0]);
      levels_new[1] = left + (right - left) * percentage;
      levels_new[0] = left;
      levels_new[2] = right;
    }
  }

  // if all zero make it gray
  if(levels_new[0] == 0.f && levels_new[1] == 0.f && levels_new[2] == 0.f)
  {
    levels_new[0] = -1.5f;
    levels_new[1] = 0.f;
    levels_new[2] = 1.5f;
  }

  // check the range
  if(levels_new[2] < levels_new[0] + 0.05f * 2.f) levels_new[2] = levels_new[0] + 0.05f * 2.f;
  if(levels_new[1] < levels_new[0] + 0.05f) levels_new[1] = levels_new[0] + 0.05f;
  if(levels_new[1] > levels_new[2] - 0.05f) levels_new[1] = levels_new[2] - 0.05f;

  {
    // set the new value but keep the middle proportional
    const float left = MAX(levels_new[0], FBW_PREVIEW_LVL_MIN);
    const float right = MIN(levels_new[2], FBW_PREVIEW_LVL_MAX);

    const float percentage = (levels_new[1] - levels_new[0]) / (levels_new[2] - levels_new[0]);
    levels_new[1] = left + (right - left) * percentage;
    levels_new[0] = left;
    levels_new[2] = right;
  }
}

static void show_hide_controls(dt_iop_module_t *self, dt_iop_fbw_gui_data_t *d, dt_iop_fbw_params_t *p)
{
  switch(p->bw_method)
  {
    case dt_iop_fbw_bw_mix:
      gtk_widget_show(GTK_WIDGET(d->vbox_rgb));
      break;
    default:
    case dt_iop_fbw_bw_max:
      gtk_widget_hide(GTK_WIDGET(d->vbox_rgb));
      break;
  }
}

// levels bar

#define FBW_LVLBAR_INSET DT_PIXEL_APPLY_DPI(5)

static float mouse_x_to_levels(const float mouse_x, const float width)
{
  return (mouse_x * ((FBW_PREVIEW_LVL_MAX - FBW_PREVIEW_LVL_MIN) / width)) + FBW_PREVIEW_LVL_MIN;
}

static float levels_to_mouse_x(const float levels, const float width)
{
  return ((levels - FBW_PREVIEW_LVL_MIN) * (width / (FBW_PREVIEW_LVL_MAX - FBW_PREVIEW_LVL_MIN)));
}

static int mouse_x_is_over_levels(const float mouse_x, const float levels, const float width)
{
  const float arrw = DT_PIXEL_APPLY_DPI(7.0f) * .5f;
  const float middle = levels_to_mouse_x(levels, width);

  return (mouse_x > middle - arrw && mouse_x < middle + arrw);
}

static int mouse_x_to_levels_index(const float mouse_x, const float levels[3], const float width)
{
  int levels_index = -1;

  const float mouse_x_left = levels_to_mouse_x(levels[0], width);
  const float mouse_x_middle = levels_to_mouse_x(levels[1], width);
  const float mouse_x_right = levels_to_mouse_x(levels[2], width);

  if(mouse_x <= mouse_x_left + (mouse_x_middle - mouse_x_left) / 2.f)
    levels_index = 0;
  else if(mouse_x <= mouse_x_middle + (mouse_x_right - mouse_x_middle) / 2.f)
    levels_index = 1;
  else
    levels_index = 2;

  return levels_index;
}

static void levels_update(const float levels[3], dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;

  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;

  float levels_old[3] = { p->levels[0], p->levels[1], p->levels[2] };

  p->levels[0] = levels[0];
  p->levels[1] = levels[1];
  p->levels[2] = levels[2];

  clamp_minmax(levels_old, p->levels);

  gtk_widget_queue_draw(g->levels_bar);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static gboolean levelsbar_leave_notify(GtkWidget *widget, GdkEventCrossing *event, dt_iop_module_t *self)
{
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;

  g->lvlbar_mouse_x = -1;
  g->lvlbar_mouse_y = -1;

  gtk_widget_queue_draw(g->levels_bar);

  return TRUE;
}

static gboolean levelsbar_button_press(GtkWidget *widget, GdkEventButton *event, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return TRUE;

  dt_iop_request_focus(self);

  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;

  const int inset = FBW_LVLBAR_INSET;

  GtkAllocation allocation;
  gtk_widget_get_allocation(widget, &allocation);
  const float width = allocation.width - 2 * inset;

  if(event->button == 1 && event->type == GDK_2BUTTON_PRESS)
  {
    // reset values
    const float levels[3] = { FBW_PREVIEW_LVL_MIN, 0.f, FBW_PREVIEW_LVL_MAX };
    levels_update(levels, self);
  }
  else if(event->button == 1)
  {
    // left slider
    if(mouse_x_is_over_levels(g->lvlbar_mouse_x, p->levels[0], width))
    {
      g->is_dragging = dt_iop_fbw_lvlbar_drag_left;
    }
    // middle slider
    else if(mouse_x_is_over_levels(g->lvlbar_mouse_x, p->levels[1], width))
    {
      g->is_dragging = dt_iop_fbw_lvlbar_drag_middle;
    }
    // right slider
    else if(mouse_x_is_over_levels(g->lvlbar_mouse_x, p->levels[2], width))
    {
      g->is_dragging = dt_iop_fbw_lvlbar_drag_right;
    }
    else
    {
      const int lvl_idx = mouse_x_to_levels_index(g->lvlbar_mouse_x, p->levels, width);
      if(lvl_idx >= 0)
      {
        float levels[3] = { p->levels[0], p->levels[1], p->levels[2] };
        levels[lvl_idx] = mouse_x_to_levels(g->lvlbar_mouse_x, width);
        levels_update(levels, self);
      }
    }
  }

  return TRUE;
}

static gboolean levelsbar_button_release(GtkWidget *widget, GdkEventButton *event, dt_iop_module_t *self)
{
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;
  if(event->button == 1)
  {
    g->is_dragging = 0;
  }
  return TRUE;
}

static gboolean levelsbar_scrolled(GtkWidget *widget, GdkEventScroll *event, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return TRUE;

  dt_iop_request_focus(self);

  int delta_y;
  if(dt_gui_get_scroll_unit_deltas(event, NULL, &delta_y))
  {
    dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;
    dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;

    const int inset = FBW_LVLBAR_INSET;

    GtkAllocation allocation;
    gtk_widget_get_allocation(widget, &allocation);
    const float width = allocation.width - 2 * inset;

    const int lvl_idx = mouse_x_to_levels_index(g->lvlbar_mouse_x, p->levels, width);
    if(lvl_idx >= 0)
    {
      float levels[3] = { p->levels[0], p->levels[1], p->levels[2] };
      levels[lvl_idx] = CLAMP(levels[lvl_idx] - (0.05 * delta_y), FBW_PREVIEW_LVL_MIN, FBW_PREVIEW_LVL_MAX);
      levels_update(levels, self);
    }
  }

  return TRUE;
}

static gboolean levelsbar_motion_notify(GtkWidget *widget, GdkEventMotion *event, dt_iop_module_t *self)
{
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  const int inset = FBW_LVLBAR_INSET;

  GtkAllocation allocation;
  gtk_widget_get_allocation(widget, &allocation);
  const float width = allocation.width - 2 * inset;
  const float height = allocation.height - 2 * inset;

  /* record mouse position within control */
  g->lvlbar_mouse_x = CLAMP(event->x - inset, 0, width);
  g->lvlbar_mouse_y = CLAMP(event->y - inset, 0, height);

  float levels[3] = { p->levels[0], p->levels[1], p->levels[2] };

  if(g->is_dragging == dt_iop_fbw_lvlbar_drag_left)
  {
    levels[0] = mouse_x_to_levels(g->lvlbar_mouse_x, width);
    levels_update(levels, self);
  }
  else if(g->is_dragging == dt_iop_fbw_lvlbar_drag_middle)
  {
    levels[1] = mouse_x_to_levels(g->lvlbar_mouse_x, width);
    levels_update(levels, self);
  }
  else if(g->is_dragging == dt_iop_fbw_lvlbar_drag_right)
  {
    levels[2] = mouse_x_to_levels(g->lvlbar_mouse_x, width);
    levels_update(levels, self);
  }

  gtk_widget_queue_draw(g->levels_bar);

  return TRUE;
}

static gboolean levelsbar_draw(GtkWidget *widget, cairo_t *crf, dt_iop_module_t *self)
{
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  const int inset = FBW_LVLBAR_INSET;
  const float arrw = DT_PIXEL_APPLY_DPI(7.0f);

  GdkRGBA color;
  GtkStyleContext *context = gtk_widget_get_style_context(widget);
  gtk_style_context_get_color(context, gtk_widget_get_state_flags(widget), &color);

  GtkAllocation allocation;
  gtk_widget_get_allocation(widget, &allocation);
  float width = allocation.width;
  float height = allocation.height;

  cairo_surface_t *cst = dt_cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
  cairo_t *cr = cairo_create(cst);

  // translate and scale
  width -= 2.f * inset;
  height -= 2.f * inset;
  cairo_save(cr);

  // draw backgrownd
  cairo_pattern_t *gradient = NULL;
  gradient = cairo_pattern_create_linear(0, 0, width, height);
  if(gradient != NULL)
  {
    cairo_pattern_add_color_stop_rgb(gradient, 0, 0., 0., 0.);
    cairo_pattern_add_color_stop_rgb(gradient, 1, .5, .5, .5);

    cairo_set_line_width(cr, 0.1);
    cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
    cairo_set_source(cr, gradient);
    cairo_rectangle(cr, inset, inset - DT_PIXEL_APPLY_DPI(2), width, height + 2. * DT_PIXEL_APPLY_DPI(2));
    cairo_fill(cr);
    cairo_stroke(cr);
    cairo_pattern_destroy(gradient);
  }

  cairo_set_antialias(cr, CAIRO_ANTIALIAS_DEFAULT);
  cairo_restore(cr);

  /* render control points handles */
  cairo_set_source_rgba(cr, color.red, color.green, color.blue, 1.);
  cairo_set_line_width(cr, DT_PIXEL_APPLY_DPI(1.));

  // draw arrows
  for(int i = 0; i < 3; i++)
  {
    const float levels_value = p->levels[i];
    const float middle = levels_to_mouse_x(levels_value, width);
    const gboolean is_under_mouse
        = g->lvlbar_mouse_x >= 0.f && mouse_x_to_levels_index(g->lvlbar_mouse_x, p->levels, width) == i;
    const int is_dragging = (g->is_dragging == dt_iop_fbw_lvlbar_drag_left && i == 0)
                            || (g->is_dragging == dt_iop_fbw_lvlbar_drag_middle && i == 1)
                            || (g->is_dragging == dt_iop_fbw_lvlbar_drag_right && i == 2);

    cairo_move_to(cr, inset + middle, height + (2 * inset) - 1);
    cairo_rel_line_to(cr, -arrw * .5f, 0);
    cairo_rel_line_to(cr, arrw * .5f, -arrw);
    cairo_rel_line_to(cr, arrw * .5f, arrw);
    cairo_close_path(cr);

    if(is_under_mouse || is_dragging)
      cairo_fill(cr);
    else
      cairo_stroke(cr);
  }

  /* push mem surface into widget */
  cairo_destroy(cr);
  cairo_set_source_surface(crf, cst, 0, 0);
  cairo_paint(crf);
  cairo_surface_destroy(cst);

  return TRUE;
}

static void bw_method_callback(GtkComboBox *combo, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;

  const int reset = darktable.gui->reset;
  darktable.gui->reset = 1;

  p->bw_method = dt_bauhaus_combobox_get((GtkWidget *)combo);

  show_hide_controls(self, g, p);

  darktable.gui->reset = reset;

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void oddness_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  p->oddness = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void red_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  p->red = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void green_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  p->green = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void blue_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  p->blue = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void auto_levels_callback(GtkToggleButton *togglebutton, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;

  const int reset = darktable.gui->reset;
  darktable.gui->reset = 1;

  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;

  if(!isnan(g->l_min) && !isnan(g->l_max) && !isnan(g->l_mid))
  {
    p->levels[0] = g->l_min;
    // p->levels[1] = g->l_mid;
    p->levels[2] = g->l_max;
    p->levels[1] = ((p->levels[2] - p->levels[0]) / 2.f) + p->levels[0];
  }

  gtk_toggle_button_set_active(togglebutton, FALSE);

  darktable.gui->reset = reset;

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

void commit_params(struct dt_iop_module_t *self, dt_iop_params_t *params, dt_dev_pixelpipe_t *pipe,
                   dt_dev_pixelpipe_iop_t *piece)
{
  memcpy(piece->data, params, sizeof(dt_iop_fbw_params_t));
}

void init_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  piece->data = malloc(sizeof(dt_iop_fbw_data_t));
  self->commit_params(self, self->default_params, pipe, piece);
}

void cleanup_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  free(piece->data);
  piece->data = NULL;
}

void gui_update(struct dt_iop_module_t *self)
{
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  dt_bauhaus_combobox_set(g->cmb_bw_method, p->bw_method);
  dt_bauhaus_slider_set(g->sl_oddness, p->oddness);
  dt_bauhaus_slider_set(g->sl_red, p->red);
  dt_bauhaus_slider_set(g->sl_green, p->green);
  dt_bauhaus_slider_set(g->sl_blue, p->blue);

  gtk_widget_queue_draw(GTK_WIDGET(g->levels_bar));

  show_hide_controls(self, g, p);

  dt_pthread_mutex_lock(&g->lock);
  g->l_min = NAN;
  g->l_max = NAN;
  g->l_mid = NAN;
  g->hash = 0;
  dt_pthread_mutex_unlock(&g->lock);
}

void init_global(dt_iop_module_so_t *module)
{
  dt_iop_fbw_global_data_t *gd = (dt_iop_fbw_global_data_t *)malloc(sizeof(dt_iop_fbw_global_data_t));
  module->data = gd;

  dt_pthread_mutex_init(&gd->fftw3_lock, NULL);

// FIXME: this should go into dt startup
#ifdef HAVE_FFTW3_OMP
#ifdef _OPENMP

  fftwf_init_threads();

#endif
#endif
}

void cleanup_global(dt_iop_module_so_t *module)
{
  dt_iop_fbw_global_data_t *gd = (dt_iop_fbw_global_data_t *)module->data;

  dt_pthread_mutex_destroy(&gd->fftw3_lock);

// FIXME: this should go into dt cleanup
#ifdef HAVE_FFTW3_OMP
#ifdef _OPENMP

  fftwf_cleanup_threads();
  fftwf_plan_with_nthreads(dt_get_num_threads());

#endif
#endif

  free(module->data);
  module->data = NULL;
}

void init(dt_iop_module_t *module)
{
  module->data = NULL;
  module->params = calloc(1, sizeof(dt_iop_fbw_params_t));
  module->default_params = calloc(1, sizeof(dt_iop_fbw_params_t));
  module->default_enabled = 0;
  module->priority = 161; // module order created by iop_dependencies.py, do not edit! // from exposure
  module->params_size = sizeof(dt_iop_fbw_params_t);
  module->gui_data = NULL;

  dt_iop_fbw_params_t tmp = { 0 };

  tmp.bw_method = dt_iop_fbw_bw_mix;
  tmp.oddness = 15.f;

  tmp.red = 0.22248840f;
  tmp.green = 0.71690369f;
  tmp.blue = 0.06060791f;

  tmp.levels[0] = FBW_PREVIEW_LVL_MIN;
  tmp.levels[1] = 0.f;
  tmp.levels[2] = FBW_PREVIEW_LVL_MAX;

  memcpy(module->params, &tmp, sizeof(dt_iop_fbw_params_t));
  memcpy(module->default_params, &tmp, sizeof(dt_iop_fbw_params_t));
}

void cleanup(dt_iop_module_t *module)
{
  free(module->params);
  module->params = NULL;
}

void gui_init(struct dt_iop_module_t *self)
{
  const int bs = DT_PIXEL_APPLY_DPI(14);

  self->gui_data = malloc(sizeof(dt_iop_fbw_gui_data_t));
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  dt_pthread_mutex_init(&g->lock, NULL);

  g->l_min = NAN;
  g->l_max = NAN;
  g->l_mid = NAN;
  g->hash = 0;

  g->auto_levels = 0;
  g->levels[0] = FBW_PREVIEW_LVL_MIN;
  g->levels[1] = 0.f;
  g->levels[2] = FBW_PREVIEW_LVL_MAX;

  g->is_dragging = 0;
  g->lvlbar_mouse_x = -1;
  g->lvlbar_mouse_y = -1;

  self->widget = GTK_WIDGET(gtk_box_new(GTK_ORIENTATION_VERTICAL, DT_BAUHAUS_SPACE));

  g->cmb_bw_method = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->cmb_bw_method, NULL, _("b&w conversion method"));
  dt_bauhaus_combobox_add(g->cmb_bw_method, _("rgb mix"));
  dt_bauhaus_combobox_add(g->cmb_bw_method, _("rgb max"));
  g_object_set(g->cmb_bw_method, "tooltip-text", _("b&w conversion method."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->cmb_bw_method), "value-changed", G_CALLBACK(bw_method_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->cmb_bw_method, TRUE, TRUE, 0);

  g->sl_oddness = dt_bauhaus_slider_new_with_range(self, 0.0, 100.0, 1.0, p->oddness, 2);
  dt_bauhaus_widget_set_label(g->sl_oddness, _("oddness"), _("oddness"));
  g_object_set(g->sl_oddness, "tooltip-text", _("oddness."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_oddness), "value-changed", G_CALLBACK(oddness_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_oddness, TRUE, TRUE, 0);

  g->vbox_rgb = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);

  g->sl_red = dt_bauhaus_slider_new_with_range(self, -2.0, 2.0, 0.005, p->red, 3);
  dt_bauhaus_widget_set_label(g->sl_red, _("red"), _("red"));
  g_object_set(g->sl_red, "tooltip-text", _("red."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_red), "value-changed", G_CALLBACK(red_callback), self);

  gtk_box_pack_start(GTK_BOX(g->vbox_rgb), g->sl_red, TRUE, TRUE, 0);

  g->sl_green = dt_bauhaus_slider_new_with_range(self, -2.0, 2.0, 0.005, p->green, 3);
  dt_bauhaus_widget_set_label(g->sl_green, _("green"), _("green"));
  g_object_set(g->sl_green, "tooltip-text", _("green."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_green), "value-changed", G_CALLBACK(green_callback), self);

  gtk_box_pack_start(GTK_BOX(g->vbox_rgb), g->sl_green, TRUE, TRUE, 0);

  g->sl_blue = dt_bauhaus_slider_new_with_range(self, -2.0, 2.0, 0.005, p->blue, 3);
  dt_bauhaus_widget_set_label(g->sl_blue, _("blue"), _("blue"));
  g_object_set(g->sl_blue, "tooltip-text", _("blue."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_blue), "value-changed", G_CALLBACK(blue_callback), self);

  gtk_box_pack_start(GTK_BOX(g->vbox_rgb), g->sl_blue, TRUE, TRUE, 0);

  gtk_box_pack_start(GTK_BOX(self->widget), g->vbox_rgb, TRUE, TRUE, 0);

  GtkWidget *prev_lvl = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);

  g->levels_bar = gtk_drawing_area_new();

  gtk_widget_set_tooltip_text(g->levels_bar, _("adjust levels."));
  g_signal_connect(G_OBJECT(g->levels_bar), "draw", G_CALLBACK(levelsbar_draw), self);
  g_signal_connect(G_OBJECT(g->levels_bar), "motion-notify-event", G_CALLBACK(levelsbar_motion_notify), self);
  g_signal_connect(G_OBJECT(g->levels_bar), "leave-notify-event", G_CALLBACK(levelsbar_leave_notify), self);
  g_signal_connect(G_OBJECT(g->levels_bar), "button-press-event", G_CALLBACK(levelsbar_button_press), self);
  g_signal_connect(G_OBJECT(g->levels_bar), "button-release-event", G_CALLBACK(levelsbar_button_release), self);
  g_signal_connect(G_OBJECT(g->levels_bar), "scroll-event", G_CALLBACK(levelsbar_scrolled), self);
  gtk_widget_add_events(GTK_WIDGET(g->levels_bar), GDK_POINTER_MOTION_MASK | GDK_POINTER_MOTION_HINT_MASK
                                                       | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK
                                                       | GDK_LEAVE_NOTIFY_MASK | GDK_SCROLL_MASK
                                                       | GDK_SMOOTH_SCROLL_MASK);
  gtk_widget_set_size_request(g->levels_bar, -1, DT_PIXEL_APPLY_DPI(5));

  g->bt_auto_levels
      = dtgtk_togglebutton_new(_fbw_cairo_paint_auto_levels, CPF_STYLE_FLAT | CPF_DO_NOT_USE_BORDER, NULL);
  g_object_set(G_OBJECT(g->bt_auto_levels), "tooltip-text", _("auto levels."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->bt_auto_levels), "toggled", G_CALLBACK(auto_levels_callback), self);
  gtk_widget_set_size_request(GTK_WIDGET(g->bt_auto_levels), bs, bs);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->bt_auto_levels), FALSE);

  gtk_box_pack_end(GTK_BOX(prev_lvl), g->bt_auto_levels, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(prev_lvl), GTK_WIDGET(g->levels_bar), TRUE, TRUE, 0);

  gtk_box_pack_start(GTK_BOX(self->widget), prev_lvl, TRUE, TRUE, 0);

  gtk_widget_show_all(g->vbox_rgb);
  gtk_widget_set_no_show_all(g->vbox_rgb, TRUE);
}

void gui_cleanup(struct dt_iop_module_t *self)
{
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;
  if(g)
  {
    dt_pthread_mutex_destroy(&g->lock);
  }

  free(self->gui_data);
  self->gui_data = NULL;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

static void check_nan(const float *im, const int size, const char *str)
{
  int i_nan = 0;

  for(int i = 0; i < size; i++)
    if(isnan(im[i])) i_nan++;

  if(i_nan > 0) printf("freaky b&w nan: %i %s\n", i_nan, str);
}

//-----------------------------------------------------------------------

#ifdef HAVE_FFTW3

static void fft(fbw_fft_t *fft_fbw, float *image_src, const int width, const int height,
                dt_pthread_mutex_t *fftw3_lock)
{
  fft_fbw->fftw3_lock = fftw3_lock;

  fft_fbw->width_src = width;
  fft_fbw->height_src = height;

  fft_fbw->width_dest = width;
  fft_fbw->height_dest = height;

  fft_fbw->width_fft = width;
  fft_fbw->height_fft = height;

  fft_fbw->width_fft_complex = fft_fbw->width_fft / 2 + 1;
  fft_fbw->height_fft_complex = fft_fbw->height_fft;

  dt_pthread_mutex_lock(fft_fbw->fftw3_lock);

  fft_fbw->in_src = (float *)fftwf_malloc(sizeof(float) * fft_fbw->width_fft * fft_fbw->height_fft);
  fft_fbw->out_src = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * fft_fbw->width_fft_complex
                                                   * fft_fbw->height_fft_complex);

  fft_fbw->plan_src = fftwf_plan_dft_r2c_2d(fft_fbw->height_fft, fft_fbw->width_fft, fft_fbw->in_src,
                                            fft_fbw->out_src, FFTW_ESTIMATE);
  fft_fbw->plan_inv = fftwf_plan_dft_c2r_2d(fft_fbw->height_fft, fft_fbw->width_fft, fft_fbw->out_src,
                                            fft_fbw->in_src, FFTW_ESTIMATE);

  dt_pthread_mutex_unlock(fft_fbw->fftw3_lock);

  memset(fft_fbw->in_src, 0, sizeof(float) * fft_fbw->width_fft * fft_fbw->height_fft);
  memset(fft_fbw->out_src, 0, sizeof(fftwf_complex) * fft_fbw->width_fft_complex * fft_fbw->height_fft_complex);

  const int w = fft_fbw->width_src;
  const int h = fft_fbw->height_src;
  const int wf = fft_fbw->width_fft;
  float *fft_in_src = fft_fbw->in_src;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(fft_in_src, image_src) schedule(static)
#endif
  for(int y = 0; y < h; y++)
  {
    float *in_src = fft_in_src + y * wf;
    const float *src = image_src + y * w;

    for(int x = 0; x < w; x++)
    {
      in_src[x] = src[x];
    }
  }

  fftwf_execute(fft_fbw->plan_src);
}

static void ifft(fbw_fft_t *fft_fbw, float *image_dest)
{
  const float scale = 1.0 / (fft_fbw->width_fft * fft_fbw->height_fft);

  memset(fft_fbw->in_src, 0, sizeof(float) * fft_fbw->width_fft * fft_fbw->height_fft);

  fftwf_execute(fft_fbw->plan_inv);

  const int w = fft_fbw->width_dest;
  const int h = fft_fbw->height_dest;
  const int wf = fft_fbw->width_fft;
  float *in_src = fft_fbw->in_src;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(in_src, image_dest) schedule(static)
#endif
  for(int y = 0; y < h; y++)
  {
    float *dest = image_dest + y * w;
    float *out_inv = in_src + y * wf;

    for(int x = 0; x < w; x++)
    {
      dest[x] = out_inv[x] * scale;
    }
  }
}

static void fft_free(fbw_fft_t *fft_fbw)
{
  dt_pthread_mutex_lock(fft_fbw->fftw3_lock);

  if(fft_fbw->plan_src) fftwf_destroy_plan(fft_fbw->plan_src);
  if(fft_fbw->plan_inv) fftwf_destroy_plan(fft_fbw->plan_inv);

  if(fft_fbw->in_src) fftwf_free(fft_fbw->in_src);
  if(fft_fbw->out_src) fftwf_free(fft_fbw->out_src);

  dt_pthread_mutex_unlock(fft_fbw->fftw3_lock);

  fft_fbw->plan_src = NULL;
  fft_fbw->plan_inv = NULL;

  fft_fbw->in_src = NULL;
  fft_fbw->out_src = NULL;
}

#else

static void fft(fbw_fft_t *fft_fbw, float *image_src, const int width, const int height)
{
  fft_fbw->width_src = width;
  fft_fbw->height_src = height;

  fft_fbw->width_dest = width;
  fft_fbw->height_dest = height;

  fft_fbw->width_fft = width;
  fft_fbw->height_fft = height;

  fft_fbw->width_fft_complex = fft_fbw->width_fft / 2 + 1;
  fft_fbw->height_fft_complex = fft_fbw->height_fft;

  fft_fbw->in_src = (float *)dt_alloc_align(64, sizeof(float) * fft_fbw->width_fft * fft_fbw->height_fft);
  fft_fbw->out_src = (fftwf_complex *)dt_alloc_align(64, sizeof(fftwf_complex) * fft_fbw->width_fft_complex
                                                             * fft_fbw->height_fft_complex);

  memset(fft_fbw->in_src, 0, sizeof(float) * fft_fbw->width_fft * fft_fbw->height_fft);
  memset(fft_fbw->out_src, 0, sizeof(fftwf_complex) * fft_fbw->width_fft_complex * fft_fbw->height_fft_complex);

  const int w = fft_fbw->width_src;
  const int h = fft_fbw->height_src;
  const int wf = fft_fbw->width_fft;
  float *fft_in_src = fft_fbw->in_src;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(fft_in_src, image_src) schedule(static)
#endif
  for(int y = 0; y < h; y++)
  {
    float *in_src = fft_in_src + y * wf;
    const float *src = image_src + y * w;

    for(int x = 0; x < w; x++)
    {
      in_src[x] = src[x];
    }
  }
}

static void ifft(fbw_fft_t *fft_fbw, float *image_dest)
{
  memset(fft_fbw->in_src, 0, sizeof(float) * fft_fbw->width_fft * fft_fbw->height_fft);

  const int w = fft_fbw->width_dest;
  const int h = fft_fbw->height_dest;
  const int wf = fft_fbw->width_fft;
  float *in_src = fft_fbw->in_src;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(in_src, image_dest) schedule(static)
#endif
  for(int y = 0; y < h; y++)
  {
    float *dest = image_dest + y * w;
    float *out_inv = in_src + y * wf;

    for(int x = 0; x < w; x++)
    {
      dest[x] = out_inv[x];
    }
  }

  if(fft_fbw->in_src) dt_free_align(fft_fbw->in_src);
  if(fft_fbw->out_src) dt_free_align(fft_fbw->out_src);
}

static void fft_free(fbw_fft_t *fft_fbw)
{
  if(fft_fbw->in_src) dt_free_align(fft_fbw->in_src);
  if(fft_fbw->out_src) dt_free_align(fft_fbw->out_src);

  fft_fbw->in_src = NULL;
  fft_fbw->out_src = NULL;
}

#endif

static inline float sRGBtoRGB(float sval)
{
  return (sval <= 0.04045f ? sval / (12.92f) : powf((sval + 0.055f) / (1.055f), 2.4f));
}

static inline float RGBtosRGB(float val)
{
  return (val <= 0.0031308f ? val * 12.92f : 1.055f * powf(val, 0.416667f) - 0.055f);
}

static inline float rgb_luminance(const float r, const float g, const float b, const float rgb[3])
{
  return (r * rgb[0] + g * rgb[1] + b * rgb[2]);
}

static void image_to_output(float *img_src, const int width, const int height, const int ch, float *img_dest)
{
  const int stride = width * height;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_dest) schedule(static)
#endif
  for(int i = 0; i < stride; i++)
  {
    img_dest[i * ch] = img_dest[i * ch + 1] = img_dest[i * ch + 2] = sRGBtoRGB(img_src[i]);
  }
}

static void pad_image_mix(float *img_src, const int width, const int height, const int ch, float *img_dest,
                          const int pad_w, const int pad_h, const float rgb[3])
{
  const int iwidth = width + pad_w * 2;
  const int iheight = height + pad_h * 2;

  memset(img_dest, 0, iwidth * iheight * sizeof(float));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_dest, rgb) schedule(static)
#endif
  for(int y = 0; y < height; y++)
  {
    float *s = img_src + y * width * ch;
    float *d = img_dest + (y + pad_h) * iwidth + pad_w;

    for(int x = 0; x < width; x++)
    {
      d[x] = (rgb_luminance(RGBtosRGB(s[x * ch + 0]), RGBtosRGB(s[x * ch + 1]), RGBtosRGB(s[x * ch + 2]), rgb));
    }
  }
}

static void pad_image_max(float *img_src, const int width, const int height, const int ch, float *img_dest,
                          const int pad_w, const int pad_h)
{
  const int iwidth = width + pad_w * 2;
  const int iheight = height + pad_h * 2;

  memset(img_dest, 0, iwidth * iheight * ch * sizeof(float));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_dest) schedule(static)
#endif
  for(int y = 0; y < height; y++)
  {
    float *s = img_src + y * width * ch;
    float *d = img_dest + (y + pad_h) * iwidth * ch + pad_w * ch;

    for(int x = 0; x < width; x++)
    {
      for(int c = 0; c < ch; c++)
      {
        d[x * ch + c] = RGBtosRGB(s[x * ch + c]);
      }
    }
  }
}

static void unpad_image(float *img_src, const int width, const int height, float *img_dest, const int pad_w,
                        const int pad_h)
{
  const int iwidth = width + pad_w * 2;

  memset(img_dest, 0, width * height * sizeof(float));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_dest) schedule(static)
#endif
  for(int y = 0; y < height; y++)
  {
    float *s = img_src + (y + pad_h) * iwidth + pad_w;
    float *d = img_dest + y * width;

    memcpy(d, s, width * sizeof(float));
  }
}

static void adjust_brcogm(float *img_src, const int width, const int height, const int ch, const float oddness,
                          const float image_scale)
{
  const int stride = width * height * ch;

  const float brightness = 0;
  float contrast = 0.f;

  if(oddness < .5f)
    contrast = -powf(oddness, 2.f) * (1.f - image_scale) * .25f;
  else
    contrast = powf(oddness, 2.f) * (1.f - image_scale) * .25f;

  if(contrast != 0.f || brightness != 0.f)
  {
    const float br = brightness;
    const float cn = expf(contrast);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src) schedule(static)
#endif
    for(int i = 0; i < stride; i++)
    {
      img_src[i] = (img_src[i] - .5f);
      img_src[i] *= cn;
      img_src[i] = (img_src[i] + .5f);

      img_src[i] += br;
    }
  }
}

static void normalize(float *img_src, const int width, const int height, const float L, const float H)
{

  float min = INFINITY;
  float max = -INFINITY;
  const int stride = width * height;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src) schedule(static) reduction(min : min) reduction(max : max)
#endif
  for(int i = 0; i < stride; i++)
  {
    min = MIN(min, img_src[i]);
    max = MAX(max, img_src[i]);
  }

  if(min == max) return;

  const float mult = (H - L) / (max - min);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, min) schedule(static)
#endif
  for(int i = 0; i < stride; i++)
  {
    img_src[i] = ((img_src[i] - min) * mult) + L;
  }
}

static void process_stats(const float *const img_src, const int width, const int height, const int ch,
                          const int pad_w, const int pad_h, float levels[3])
{
  const int ch1 = (ch == 4) ? 3 : ch;

  float l_max = -INFINITY;
  float l_min = INFINITY;
  float l_sum = 0.f;

#ifdef _OPENMP
#pragma omp parallel for default(none) schedule(static) reduction(+ : l_sum) reduction(min : l_min)               \
                                                                      reduction(max : l_max)
#endif
  for(int y = pad_h; y < height - pad_h; y++)
  {
    const float *const s = img_src + y * width * ch;

    for(int x = pad_w; x < width - pad_w; x++)
    {
      for(int c = 0; c < ch1; c++)
      {
        const float l = s[x * ch + c];
        l_max = MAX(l_max, l);
        l_min = MIN(l_min, l);
        l_sum += l;
      }
    }
  }

  levels[0] = l_min;
  levels[2] = l_max;
  levels[1] = (l_sum / (float)((height - 2 * pad_h) * (width - 2 * pad_w)));
}

static void adjust_levels(float *img_src, const int width, const int height, const int ch, const float levels[3])
{
  const int size = width * height * ch;

  const float left = levels[0];
  const float middle = levels[1];
  const float right = levels[2];

  if(left == FBW_PREVIEW_LVL_MIN && middle == 0.f && right == FBW_PREVIEW_LVL_MAX) return;

  const float delta = (right - left) / 2.0f;
  const float mid = left + delta;
  const float tmp = (middle - mid) / delta;
  const float in_inv_gamma = pow(10, tmp);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src) schedule(static)
#endif
  for(int i = 0; i < size; i += ch)
  {
    for(int c = 0; c < 1; c++)
    {
      float L_in = img_src[i + c];

      if(L_in <= left)
      {
        img_src[i + c] = 0.f;
      }
      else
      {
        const float percentage = (L_in - left) / (right - left);
        img_src[i + c] = powf(percentage, in_inv_gamma);
      }
    }
  }
}

static void gradient_rgb_mix(float *img_src, float *img_grx, float *img_gry, const int width, const int height,
                             const int pad_w, const int pad_h, const float _oddness, const float image_scale)
{
  const float oddness = _oddness * sqrtf(image_scale);

  memset(img_grx, 0, width * height * sizeof(float));
  memset(img_gry, 0, width * height * sizeof(float));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_grx, img_gry) schedule(static)
#endif
  for(int y = 1; y < height - 1; y++)
  {
    float *s0 = img_src + y * width;
    float *s1 = img_src + (y + 1) * width;
    float *dx = img_grx + y * width;
    float *dy = img_gry + y * width;

    for(int x = 1; x < width - 1; x++)
    {
      float grx = s0[x + 1] - s0[x];
      float gry = s1[x] - s0[x];
      float grn = sqrtf((grx * grx + gry * gry));

      if(grn > 0.f)
      {
        float n = powf(grn, oddness);


        dx[x] = grx / n;
        dy[x] = gry / n;
      }
    }
  }
}

static void gradient_rgb_max(float *img_src, float *img_grx, float *img_gry, const int width, const int height,
                             const int pad_w, const int pad_h, const float _oddness, const float image_scale)
{
  const int ch = 4;
  const float oddness = _oddness * sqrtf(image_scale);

  memset(img_grx, 0, width * height * sizeof(float));
  memset(img_gry, 0, width * height * sizeof(float));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_grx, img_gry) schedule(static)
#endif
  for(int y = 1; y < height - 1; y++)
  {
    float *s0 = img_src + y * width * ch;
    float *s1 = img_src + (y + 1) * width * ch;
    float *dx = img_grx + y * width;
    float *dy = img_gry + y * width;

    for(int x = 1; x < width - 1; x++)
    {
      float RGBx[3] = { 0 };
      float RGBy[3] = { 0 };
      float RGBn[3] = { 0 };

      for(int c = 0; c < 3; c++)
      {
        RGBx[c] = s0[(x + 1) * ch + c] - s0[x * ch + c];
        RGBy[c] = s1[x * ch + c] - s0[x * ch + c];

        RGBn[c] = RGBx[c] * RGBx[c] + RGBy[c] * RGBy[c];
      }

      int max_bn = 0;
      if(RGBn[0] > RGBn[1])
        max_bn = (RGBn[0] > RGBn[2]) ? 0 : 2;
      else
        max_bn = (RGBn[1] > RGBn[2]) ? 1 : 2;

      float n = 1e-5 + powf(RGBn[max_bn], oddness);

      dx[x] = RGBx[max_bn] / n;
      dy[x] = RGBy[max_bn] / n;
    }
  }
}

static void estimate_laplacian(float *img_grx, float *img_gry, float *img_dest, const int width, const int height,
                               const int pad_w, const int pad_h)
{
  const int stride = width * height;

  float *img_bp = NULL;

  img_bp = dt_alloc_align(64, width * height * sizeof(float));
  if(img_bp == NULL) goto cleanup;

  memset(img_bp, 0, stride * sizeof(float));
  memset(img_dest, 0, stride * sizeof(float));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_grx) schedule(static)
#endif
  for(int y = 0; y < height; y++)
  {
    float *s0 = img_grx + y * width;

    for(int x = width - 1; x > 0; x--)
    {
      s0[x] -= s0[x - 1];
    }
  }

  for(int y = height - 1; y > 0; y--)
  {
    float *s0 = img_gry + y * width;
    float *s1 = img_gry + (y - 1) * width;

    for(int x = 0; x < width; x++)
    {
      s0[x] -= s1[x];
    }
  }

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_grx, img_gry, img_dest) schedule(static)
#endif
  for(int i = 0; i < stride; i++)
  {
    img_dest[i] = (img_grx[i] + img_gry[i]);
  }

cleanup:

  if(img_bp) dt_free_align(img_bp);
}

static void recontruct_laplacian(float *img_src, float *img_dest, const int width, const int height,
                                 dt_pthread_mutex_t *fftw3_lock)
{
  fbw_fft_t fft_fbw = { 0 };

  fft(&fft_fbw, img_src, width, height, fftw3_lock);

  const int width_fft = fft_fbw.width_fft;
  const int height_fft = fft_fbw.height_fft;

  const float piw = (2.f * (float)M_PI / (float)width_fft);
  const float pih = (2.f * (float)M_PI / (float)height_fft);

  const int w = fft_fbw.width_fft_complex;
  const int h = fft_fbw.height_fft_complex;
  fftwf_complex *out_src = fft_fbw.out_src;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(out_src) schedule(static)
#endif
  for(int y = 0; y < h; y++)
  {
    fftwf_complex *d = out_src + y * w;

    const float cos_y = cos(pih * y);

    for(int x = 0; x < w; x++)
    {
      const float cos_x = cos(piw * x);

      if(x > 0 || y > 0)
      {
        float cos_xy = (cos_x + cos_y - 2.f) * 2.f;

        d[x][0] /= cos_xy;
        d[x][1] /= cos_xy;
      }
    }
  }

  fft_fbw.out_src[0][0] = 0.f;
  fft_fbw.out_src[0][1] = 0.f;

  ifft(&fft_fbw, img_dest);

  fft_free(&fft_fbw);
}

static void fbw_process(float *img_src, float *img_dest, const int width, const int height, const int ch,
                        const int bw_method, const float oddness, const float *_levels, const float red,
                        const float green, const float blue, const float image_scale, dt_iop_fbw_gui_data_t *g,
                        dt_dev_pixelpipe_iop_t *piece, struct dt_iop_module_t *self)
{
  dt_iop_fbw_global_data_t *gd = (dt_iop_fbw_global_data_t *)self->data;

  float *img_padded = NULL;
  float *img_grx = NULL;
  float *img_gry = NULL;

  const int pad_w = MIN(width, height) / 2;
  const int pad_h = pad_w;

  const int iwidth = width + pad_w * 2;
  const int iheight = height + pad_h * 2;

  float l_min = 0.f;
  float l_max = 1.f;
  float l_mid = .5f;

  float tmp_l_min = NAN;
  float tmp_l_max = NAN;
  float tmp_l_mid = NAN;

  const float rgb[3] = { red, green, blue };

  img_grx = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if(img_grx == NULL) goto cleanup;

  img_gry = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if(img_gry == NULL) goto cleanup;

  img_padded = dt_alloc_align(64, iwidth * iheight * ((bw_method == dt_iop_fbw_bw_mix) ? 1 : ch) * sizeof(float));
  if(img_padded == NULL) goto cleanup;

  if(bw_method == dt_iop_fbw_bw_mix)
    pad_image_mix(img_src, width, height, ch, img_padded, pad_w, pad_h, rgb);
  else
    pad_image_max(img_src, width, height, ch, img_padded, pad_w, pad_h);

  // get min and max from full image, as DT_DEV_PIXELPIPE_FULL may be cropped
  if(self->dev->gui_attached && g && piece->pipe->type == DT_DEV_PIXELPIPE_FULL)
  {
    dt_pthread_mutex_lock(&g->lock);
    const uint64_t hash = g->hash;
    dt_pthread_mutex_unlock(&g->lock);

    // note that the case 'hash == 0' on first invocation in a session implies that g->l_min/lmax
    // is NAN which initiates special handling below to avoid inconsistent results. in all
    // other cases we make sure that the preview pipe has left us with proper readings for
    // l_min/lmax. if data are not yet there we need to wait (with timeout).
    if(hash != 0 && !dt_dev_sync_pixelpipe_hash(self->dev, piece->pipe, 0, self->priority, &g->lock, &g->hash))
      dt_control_log(_("fbw inconsistent output"));

    dt_pthread_mutex_lock(&g->lock);
    tmp_l_min = g->l_min;
    tmp_l_max = g->l_max;
    tmp_l_mid = g->l_mid;
    dt_pthread_mutex_unlock(&g->lock);
  }

  // in all other cases we calculate l_min/lmax here
  if(isnan(tmp_l_min) || isnan(tmp_l_max) || isnan(tmp_l_mid))
  {
    float levels_tmp[3] = { 0 };
    process_stats(img_padded, iwidth, iheight, ((bw_method == dt_iop_fbw_bw_mix) ? 1 : ch), pad_w, pad_h,
                  levels_tmp);
    l_min = levels_tmp[0];
    l_max = levels_tmp[2];
    l_mid = levels_tmp[1];
  }
  else
  {
    l_min = tmp_l_min;
    l_max = tmp_l_max;
    l_mid = tmp_l_mid;
  }

  // PREVIEW pixelpipe stores l_min/lmax
  if(self->dev->gui_attached && g && piece->pipe->type == DT_DEV_PIXELPIPE_PREVIEW)
  {
    uint64_t hash = dt_dev_hash_plus(self->dev, piece->pipe, 0, self->priority);
    dt_pthread_mutex_lock(&g->lock);
    g->l_min = l_min;
    g->l_max = l_max;
    g->l_mid = l_mid;
    g->hash = hash;
    dt_pthread_mutex_unlock(&g->lock);
  }

  if(bw_method == dt_iop_fbw_bw_mix)
    gradient_rgb_mix(img_padded, img_grx, img_gry, iwidth, iheight, pad_w, pad_h, oddness, image_scale);
  else
    gradient_rgb_max(img_padded, img_grx, img_gry, iwidth, iheight, pad_w, pad_h, oddness, image_scale);

  float *img_l = img_padded;

  estimate_laplacian(img_grx, img_gry, img_l, iwidth, iheight, pad_w, pad_h);

  recontruct_laplacian(img_l, img_l, iwidth, iheight, &(gd->fftw3_lock));

  float *img_unpad_gr = img_grx;

  unpad_image(img_l, width, height, img_unpad_gr, pad_w, pad_h);

  {
    float levels_tmp[3] = { 0 };
    process_stats(img_unpad_gr, width, height, 1, 0, 0, levels_tmp);

    printf("%f\t%f\t%f\t%f\t%f\n", levels_tmp[0], levels_tmp[2], levels_tmp[1], image_scale, oddness);
  }

  normalize(img_unpad_gr, width, height, l_min, l_max);

  adjust_levels(img_unpad_gr, width, height, 1, _levels);

  if(0) adjust_brcogm(img_unpad_gr, width, height, 1, oddness, image_scale);

  image_to_output(img_unpad_gr, width, height, ch, img_dest);

cleanup:

  if(img_grx) dt_free_align(img_grx);
  if(img_gry) dt_free_align(img_gry);
  if(img_padded) dt_free_align(img_padded);
}

#undef FBW_PREVIEW_LVL_MIN
#undef FBW_PREVIEW_LVL_MAX

void process_internal(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid,
                      void *const ovoid, const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  const dt_iop_fbw_data_t *const p = (const dt_iop_fbw_data_t *const)piece->data;
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;
  const float image_scale = roi_in->scale / piece->iscale;

  fbw_process((float *)ivoid, (float *)ovoid, roi_in->width, roi_in->height, piece->colors, p->bw_method,
              p->oddness / 100.f, p->levels, p->red, p->green, p->blue, image_scale, g, piece, self);

  if(piece->pipe->mask_display & DT_DEV_PIXELPIPE_DISPLAY_MASK)
    dt_iop_alpha_copy(ivoid, ovoid, roi_out->width, roi_out->height);

  check_nan((float *)ovoid, roi_in->width * roi_in->height * piece->colors, "process_internal");
}

void process(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid,
             void *const ovoid, const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  process_internal(self, piece, ivoid, ovoid, roi_in, roi_out);
}

// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
