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
#include "common/LocalLaplacian_eh.h"

DT_MODULE_INTROSPECTION(1, dt_iop_loclaplab_params_t)

#define LOCLAP_MAX_SCALES 8

typedef struct dt_iop_loclaplab_params_t
{
  float alpha;
  float beta;
  float scales[LOCLAP_MAX_SCALES];
} dt_iop_loclaplab_params_t;

typedef struct dt_iop_loclaplab_gui_data_t
{
  GtkWidget *sl_alpha;
  GtkWidget *sl_beta;
//  GtkWidget *sl_scales[LOCLAP_MAX_SCALES];
} dt_iop_loclaplab_gui_data_t;

typedef struct dt_iop_loclaplab_params_t dt_iop_loclaplab_data_t;

const char *name()
{
  return _("loclaplab_eh");
}

int groups()
{
  return IOP_GROUP_TONE;
}

int flags()
{
  return IOP_FLAGS_SUPPORTS_BLENDING;
}


void commit_params(struct dt_iop_module_t *self, dt_iop_params_t *params, dt_dev_pixelpipe_t *pipe,
    dt_dev_pixelpipe_iop_t *piece)
{
  memcpy(piece->data, params, sizeof(dt_iop_loclaplab_params_t));
}

void init_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  piece->data = malloc(sizeof(dt_iop_loclaplab_data_t));
  self->commit_params(self, self->default_params, pipe, piece);
}

void cleanup_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  free(piece->data);
  piece->data = NULL;
}

void gui_update(struct dt_iop_module_t *self)
{
  dt_iop_loclaplab_gui_data_t *g = (dt_iop_loclaplab_gui_data_t *)self->gui_data;
  dt_iop_loclaplab_params_t *p = (dt_iop_loclaplab_params_t *)self->params;

  dt_bauhaus_slider_set(g->sl_alpha, p->alpha);
  dt_bauhaus_slider_set(g->sl_beta, p->beta);
  //  for (int i = 0; i < LOCLAP_MAX_SCALES; i++)
  //  	dt_bauhaus_slider_set(g->sl_scales[i], p->scales[i]);
}

void init(dt_iop_module_t *module)
{
  module->data = NULL;
  module->params = calloc(1, sizeof(dt_iop_loclaplab_params_t));
  module->default_params = calloc(1, sizeof(dt_iop_loclaplab_params_t));
  module->default_enabled = 0;
  module->priority = 588; // module order created by iop_dependencies.py, do not edit!
  module->params_size = sizeof(dt_iop_loclaplab_params_t);
  module->gui_data = NULL;

  dt_iop_loclaplab_params_t tmp = {0};

  memcpy(module->params, &tmp, sizeof(dt_iop_loclaplab_params_t));
  memcpy(module->default_params, &tmp, sizeof(dt_iop_loclaplab_params_t));

}

void cleanup(dt_iop_module_t *module)
{
  free(module->params);
  module->params = NULL;
}

static void alpha_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  dt_iop_loclaplab_params_t *p = (dt_iop_loclaplab_params_t *)self->params;

  p->alpha = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void beta_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  dt_iop_loclaplab_params_t *p = (dt_iop_loclaplab_params_t *)self->params;

  p->beta = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}
/*
static void scales_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  dt_iop_loclaplab_params_t *p = (dt_iop_loclaplab_params_t *)self->params;
  dt_iop_loclaplab_gui_data_t *g = (dt_iop_loclaplab_gui_data_t *)self->gui_data;

  for (int i = 0; i < LOCLAP_MAX_SCALES; i++)
  {
  	if (slider == g->sl_scales[i])
  	{
  		p->scales[i] = dt_bauhaus_slider_get(slider);
  		break;
  	}
  }

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}
 */
void gui_init(struct dt_iop_module_t *self)
{
  self->gui_data = malloc(sizeof(dt_iop_loclaplab_gui_data_t));
  dt_iop_loclaplab_gui_data_t *g = (dt_iop_loclaplab_gui_data_t *)self->gui_data;

  self->widget = GTK_WIDGET(gtk_box_new(GTK_ORIENTATION_VERTICAL, DT_BAUHAUS_SPACE));

  g->sl_alpha = dt_bauhaus_slider_new_with_range(self, -5.0, 5.0, 0.05, 0.0, 3);
  dt_bauhaus_widget_set_label(g->sl_alpha, _("alpha"), _("alpha"));
  g_object_set(g->sl_alpha, "tooltip-text", _("specifies how much of an effect to apply.\n"
      " 0 gives no effect,\n1 produces a moderate amount of contrast"
      " enhancement,\n-1 produces a piecewise flattening."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_alpha), "value-changed", G_CALLBACK(alpha_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_alpha, TRUE, TRUE, 0);

  g->sl_beta = dt_bauhaus_slider_new_with_range(self, -5.0, 5.0, 0.05, 0.0, 3);
  dt_bauhaus_widget_set_label(g->sl_beta, _("beta"), _("beta"));
  g_object_set(g->sl_beta, "tooltip-text", _("specifies how the effect should change with respect to"
      " scale.\n0 applies the same effect at all scales,\n1 applies the effect only at fine scales,\n"
      "-1 applies the effectonly at coarse scales."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_beta), "value-changed", G_CALLBACK(beta_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_beta, TRUE, TRUE, 0);

  /*  for (int i = 0; i < LOCLAP_MAX_SCALES; i++)
  {
    g->sl_scales[i] = dt_bauhaus_slider_new_with_range(self, -5.0, 5.0, 0.05, 0.0, 3);
    dt_bauhaus_widget_set_label(g->sl_scales[i], _("scale"), _("scale"));
    g_object_set(g->sl_scales[i], "tooltip-text", _("scale."), (char *)NULL);
    g_signal_connect(G_OBJECT(g->sl_scales[i]), "value-changed", G_CALLBACK(scales_callback), self);

    gtk_box_pack_start(GTK_BOX(self->widget), g->sl_scales[i], TRUE, TRUE, 0);
  }*/

}

void gui_cleanup(struct dt_iop_module_t *self)
{
  free(self->gui_data);
  self->gui_data = NULL;
}

void process_internal(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid, void *const ovoid,
    const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out, const int use_sse)
{
  const dt_iop_loclaplab_data_t *const p = (const dt_iop_loclaplab_data_t *const)piece->data;

  const dt_iop_colorspace_type_t cst = dt_iop_module_colorspace(self);
  const int stride = roi_in->width * roi_in->height;
  const int ch = piece->colors;

  // just to be sure...
  if (cst != iop_cs_Lab) return;

  memcpy(ovoid, ivoid, roi_in->width * roi_in->height * ch * sizeof(float));

  float *image_L = NULL;
  image_L = dt_alloc_align(64, roi_in->width * roi_in->height * sizeof(float));

  float *in = (float*)ivoid;

  for (int i = 0; i < stride; i++, in += ch)
    image_L[i] = *in;

  loclap_LocalLaplacian(image_L, image_L, roi_out->width, roi_out->height, p->alpha, p->beta, NULL/*p->scales*/, (roi_in->scale / piece->iscale), use_sse);

  float *out = (float*)ovoid;

  for (int i = 0; i < stride; i++, out += ch)
    *out = image_L[i];

  if (image_L) dt_free_align(image_L);

}

void process(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid, void *const ovoid,
    const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  process_internal(self, piece, ivoid, ovoid, roi_in, roi_out, 0);
}

#if defined(__SSE__)
void process_sse2(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid, void *const ovoid,
    const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  process_internal(self, piece, ivoid, ovoid, roi_in, roi_out, 1);
}
#endif

// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
