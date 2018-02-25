/*
    This file is part of darktable,
    copyright (c) 2017 edgardo hoszowski.

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
#@gui : sep = separator(), Preview type = choice("Full","Forward horizontal","Forward vertical","Backward horizontal","Backward vertical","Duplicate top","Duplicate left","Duplicate bottom","Duplicate right")
#@gui : sep = separator(), note = note("<small>Author: <i>David Tschumperl&#233;</i>.      Latest update: <i>09/30/2015</i>.</small>")
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

typedef struct dt_iop_fbw_params_t
{
  float oddness;
  float brightness;
  float contrast;
  float gamma;
  float beta;
} dt_iop_fbw_params_t;

typedef struct dt_iop_fbw_gui_data_t
{
  GtkWidget *sl_oddness;
  GtkWidget *sl_brightness;
  GtkWidget *sl_contrast;
  GtkWidget *sl_gamma;
  GtkWidget *sl_beta;
} dt_iop_fbw_gui_data_t;

typedef struct dt_iop_fbw_params_t dt_iop_fbw_data_t;

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

  dt_bauhaus_slider_set(g->sl_oddness, p->oddness);
  dt_bauhaus_slider_set(g->sl_brightness, p->brightness);
  dt_bauhaus_slider_set(g->sl_contrast, p->contrast);
  dt_bauhaus_slider_set(g->sl_gamma, p->gamma);
  dt_bauhaus_slider_set(g->sl_beta, p->beta);

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
  fftwf_plan_with_nthreads(dt_get_num_threads());
  
#endif
#endif
  
}

void init(dt_iop_module_t *module)
{
  module->data = NULL;
  module->params = calloc(1, sizeof(dt_iop_fbw_params_t));
  module->default_params = calloc(1, sizeof(dt_iop_fbw_params_t));
  module->default_enabled = 0;
//  module->priority = 588; // module order created by iop_dependencies.py, do not edit! // from bilat (local contrast)
  module->priority = 161; // module order created by iop_dependencies.py, do not edit! // from exposure
//  module->priority = 955; // module order created by iop_dependencies.py, do not edit! // from borders
//  module->priority = 147; // module order created by iop_dependencies.py, do not edit!
//  module->priority = 632; // module order created by iop_dependencies.py, do not edit! // from monochrome
//  module->priority = 294; // module order created by iop_dependencies.py, do not edit! // basecurve
  module->params_size = sizeof(dt_iop_fbw_params_t);
  module->gui_data = NULL;

  dt_iop_fbw_params_t tmp = {0};

  tmp.oddness = 10.f;
  tmp.brightness = 0.f;
  tmp.contrast = 0.f;
  tmp.gamma = 0.f;
  
  memcpy(module->params, &tmp, sizeof(dt_iop_fbw_params_t));
  memcpy(module->default_params, &tmp, sizeof(dt_iop_fbw_params_t));

}

void cleanup(dt_iop_module_t *module)
{
  free(module->params);
  module->params = NULL;
}

static void oddness_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  p->oddness = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void brightness_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  p->brightness = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void contrast_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  p->contrast = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void gamma_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  p->gamma = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void beta_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  p->beta = dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

void gui_init(struct dt_iop_module_t *self)
{
  self->gui_data = malloc(sizeof(dt_iop_fbw_gui_data_t));
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;
  dt_iop_fbw_params_t *p = (dt_iop_fbw_params_t *)self->params;

  self->widget = GTK_WIDGET(gtk_box_new(GTK_ORIENTATION_VERTICAL, DT_BAUHAUS_SPACE));

  g->sl_oddness = dt_bauhaus_slider_new_with_range(self, 0.0, 100.0, 1.0, p->oddness, 2);
  dt_bauhaus_widget_set_label(g->sl_oddness, _("oddness"), _("oddness"));
  g_object_set(g->sl_oddness, "tooltip-text", _("oddness."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_oddness), "value-changed", G_CALLBACK(oddness_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_oddness, TRUE, TRUE, 0);

  g->sl_brightness = dt_bauhaus_slider_new_with_range(self, -100.0, 100.0, 1.0, p->brightness, 2);
  dt_bauhaus_widget_set_label(g->sl_brightness, _("brightness"), _("brightness"));
  g_object_set(g->sl_brightness, "tooltip-text", _("brightness."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_brightness), "value-changed", G_CALLBACK(brightness_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_brightness, TRUE, TRUE, 0);

  g->sl_contrast = dt_bauhaus_slider_new_with_range(self, -100.0, 100.0, 1.0, p->contrast, 2);
  dt_bauhaus_widget_set_label(g->sl_contrast, _("contrast"), _("contrast"));
  g_object_set(g->sl_contrast, "tooltip-text", _("contrast."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_contrast), "value-changed", G_CALLBACK(contrast_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_contrast, TRUE, TRUE, 0);

  g->sl_gamma = dt_bauhaus_slider_new_with_range(self, -100.0, 100.0, 1.0, p->gamma, 2);
  dt_bauhaus_widget_set_label(g->sl_gamma, _("gamma"), _("gamma"));
  g_object_set(g->sl_gamma, "tooltip-text", _("gamma."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_gamma), "value-changed", G_CALLBACK(gamma_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_gamma, TRUE, TRUE, 0);

  g->sl_beta = dt_bauhaus_slider_new_with_range(self, 0.0, 10000.0, 1.0, 0.0, 8);
  dt_bauhaus_widget_set_label(g->sl_beta, _("beta"), _("beta"));
  g_object_set(g->sl_beta, "tooltip-text", _("beta."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_beta), "value-changed", G_CALLBACK(beta_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_beta, TRUE, TRUE, 0);


}

void gui_cleanup(struct dt_iop_module_t *self)
{
  free(self->gui_data);
  self->gui_data = NULL;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

static void check_nan(const float* im, const int size, const char *str)
{
  int i_nan = 0;

  for (int i = 0; i < size; i++) if ( isnan(im[i]) ) i_nan++;

  if (i_nan > 0) printf("freaky b&w nan: %i %s\n", i_nan, str);
}

//-----------------------------------------------------------------------

#ifdef HAVE_FFTW3

// from:
// https://github.com/jeremyfix/FFTConvolution
//

// Code adapted from gsl/fft/factorize.c
static void fftw3_factorize(const int n, int *n_factors, int factors[], int * implemented_factors)
{
  int nf = 0;
  int ntest = n;
  int factor;
  int i = 0;

  if (n == 0)
  {
    printf("Length n must be positive integer\n");
    return ;
  }

  if (n == 1)
  {
    factors[0] = 1;
    *n_factors = 1;
    return ;
  }

  /* deal with the implemented factors */
  while (implemented_factors[i] && ntest != 1)
  {
    factor = implemented_factors[i];
    while ((ntest % factor) == 0)
    {
      ntest = ntest / factor;
      factors[nf] = factor;
      nf++;
    }
    i++;
  }

  // Ok that's it
  if(ntest != 1)
  {
    factors[nf] = ntest;
    nf++;
  }

  /* check that the factorization is correct */
  {
    int product = 1;

    for (i = 0; i < nf; i++)
    {
      product *= factors[i];
    }

    if (product != n)
    {
      printf("factorization failed");
    }
  }

  *n_factors = nf;
}

static int fftw3_is_optimal(int n, int * implemented_factors)
{
  // We check that n is not a multiple of 4*4*4*2
  if(n % 4*4*4*2 == 0)
    return 0;

  int nf = 0;
  int factors[64];
  int i = 0;
  
  fftw3_factorize(n, &nf, factors,implemented_factors);

  // We just have to check if the last factor belongs to GSL_FACTORS
  while(implemented_factors[i])
  {
    if(factors[nf-1] == implemented_factors[i])
      return 1;
    ++i;
  }
  
  return 0;
}

static int fftw3_find_closest_factor(int n, int * implemented_factor)
{
  int j;
  if (fftw3_is_optimal(n,implemented_factor))
    return n;
  else
  {
    j = n+1;
    while (!fftw3_is_optimal(j,implemented_factor))
      ++j;
    return j;
  }
}

static void fft(fbw_fft_t *fft_fbw, float *image_src, const int width, const int height, const float image_scale)
{
  int FFTW_FACTORS[7] = {13,11,7,5,3,2,0}; // end with zero to detect the end of the array
  
  fft_fbw->width_src = width;
  fft_fbw->height_src = height;

  fft_fbw->width_dest = width;
  fft_fbw->height_dest = height;

  fft_fbw->width_fft = width;
  fft_fbw->height_fft = height;

if(0)  fft_fbw->width_fft = fftw3_find_closest_factor(fft_fbw->width_fft, FFTW_FACTORS);
if(0)  fft_fbw->height_fft = fftw3_find_closest_factor(fft_fbw->height_fft, FFTW_FACTORS);

  fft_fbw->width_fft_complex = fft_fbw->width_fft/2+1;
  fft_fbw->height_fft_complex = fft_fbw->height_fft;

  dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
  
  fft_fbw->in_src = (float*)fftwf_malloc(sizeof(float) * fft_fbw->width_fft * fft_fbw->height_fft);
  fft_fbw->out_src = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * fft_fbw->width_fft_complex * fft_fbw->height_fft_complex);
  
  fft_fbw->plan_src = fftwf_plan_dft_r2c_2d(fft_fbw->height_fft, fft_fbw->width_fft, fft_fbw->in_src, fft_fbw->out_src, FFTW_ESTIMATE);
  fft_fbw->plan_inv = fftwf_plan_dft_c2r_2d(fft_fbw->height_fft, fft_fbw->width_fft, fft_fbw->out_src, fft_fbw->in_src, FFTW_ESTIMATE);

  dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);

  memset(fft_fbw->in_src, 0, sizeof(float) * fft_fbw->width_fft * fft_fbw->height_fft);
  memset(fft_fbw->out_src, 0, sizeof(fftwf_complex) * fft_fbw->width_fft_complex * fft_fbw->height_fft_complex);

  const int w = fft_fbw->width_src;
  const int h = fft_fbw->height_src;
  const int wf = fft_fbw->width_fft;
  float *fft_in_src = fft_fbw->in_src;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(fft_in_src, image_src) schedule(static)
#endif
  for (int y = 0; y < h; y++)
  {
    float *in_src = fft_in_src + y * wf;
    const float *src = image_src + y * w;

    for (int x = 0; x < w; x++)
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
  for (int y = 0; y < h; y++)
  {
    float *dest = image_dest + y * w;
    float *out_inv = in_src + y * wf;

    for (int x = 0; x < w; x++)
    {
      dest[x] = out_inv[x] * scale;
    }
  }

}

static void fft_free(fbw_fft_t *fft_fbw)
{

  dt_pthread_mutex_lock(&darktable.plugin_threadsafe);

  if (fft_fbw->plan_src) fftwf_destroy_plan(fft_fbw->plan_src);
  if (fft_fbw->plan_inv) fftwf_destroy_plan(fft_fbw->plan_inv);
  
  if (fft_fbw->in_src) fftwf_free(fft_fbw->in_src);
  if (fft_fbw->out_src) fftwf_free(fft_fbw->out_src);

  dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);

  fft_fbw->plan_src = NULL;
  fft_fbw->plan_inv = NULL;
  
  fft_fbw->in_src = NULL;
  fft_fbw->out_src = NULL;

}


#else

static void fft(fbw_fft_t *fft_fbw, float *image_src, const int width, const int height, const float image_scale)
{
  fft_fbw->width_src = width;
  fft_fbw->height_src = height;

  fft_fbw->width_dest = width;
  fft_fbw->height_dest = height;

  fft_fbw->width_fft = width;
  fft_fbw->height_fft = height;

  fft_fbw->width_fft_complex = fft_fbw->width_fft/2+1;
  fft_fbw->height_fft_complex = fft_fbw->height_fft;

  fft_fbw->in_src = (float*)dt_alloc_align(64, sizeof(float) * fft_fbw->width_fft * fft_fbw->height_fft);
  fft_fbw->out_src = (fftwf_complex*)dt_alloc_align(64, sizeof(fftwf_complex) * fft_fbw->width_fft_complex * fft_fbw->height_fft_complex);
  
  memset(fft_fbw->in_src, 0, sizeof(float) * fft_fbw->width_fft * fft_fbw->height_fft);
  memset(fft_fbw->out_src, 0, sizeof(fftwf_complex) * fft_fbw->width_fft_complex * fft_fbw->height_fft_complex);

  const int w = fft_fbw->width_src;
  const int h = fft_fbw->height_src;
  const int wf = fft_fbw->width_fft;
  float *fft_in_src = fft_fbw->in_src;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(fft_in_src, image_src) schedule(static)
#endif
  for (int y = 0; y < h; y++)
  {
    float *in_src = fft_in_src + y * wf;
    const float *src = image_src + y * w;

    for (int x = 0; x < w; x++)
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
  for (int y = 0; y < h; y++)
  {
    float *dest = image_dest + y * w;
    float *out_inv = in_src + y * wf;

    for (int x = 0; x < w; x++)
    {
      dest[x] = out_inv[x];
    }
  }

  if (fft_fbw->in_src) dt_free_align(fft_fbw->in_src);
  if (fft_fbw->out_src) dt_free_align(fft_fbw->out_src);

}

#endif

static inline float sRGBtoRGB(float sval, const float beta)
{
  return (sval<=0.04045f?sval/(12.92f):powf((sval + 0.055f)/(1.055f),2.4f));
}

static inline float RGBtosRGB(float val, const float beta)
{
  return (val<=0.0031308f?val*12.92f:1.055f*powf(val,0.416667f) - 0.055f);
}

static inline float rgb_luminance(const float r, const float g, const float b)
{
  return (r * 0.22248840f + g * 0.71690369f + b * 0.06060791f);
}

static void image_to_output(float *img_src, const int width, const int height, const int ch, float *img_dest, 
    const dt_iop_colorspace_type_t cst)
{
  const int stride = width * height;
  
  if (cst == iop_cs_rgb)
  {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_dest) schedule(static)
#endif
    for (int i = 0; i < stride; i++)
    {
      img_dest[i*ch] = img_dest[i*ch+1] = img_dest[i*ch+2] = sRGBtoRGB(img_src[i], 1.f);
    }
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_dest) schedule(static)
#endif
    for (int i = 0; i < stride; i++)
    {
      float sRGB[3] = { img_src[i], img_src[i], img_src[i] };
      float XYZ[3] = {0};
      
      dt_sRGB_to_XYZ(sRGB, XYZ);
      dt_XYZ_to_Lab(XYZ, img_dest + i*ch);
    }
  }
}

static void pad_image(float *img_src, const int width, const int height, const int ch, float *img_dest, 
    const int pad_w, const int pad_h, 
    float *l_min, float *l_max,
    const dt_iop_colorspace_type_t cst, const float beta)
{
  const int iwidth = width + pad_w*2;
  const int iheight = height + pad_h*2;

  float min = INFINITY;
  float max = -INFINITY;
      
  memset(img_dest, 0, iwidth * iheight * sizeof(float));

  if (cst == iop_cs_rgb)
  {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_dest) schedule(static) reduction(min : min) reduction(max : max)
#endif
    for (int y = 0; y < height; y++)
    {
      float *s = img_src + y * width * ch;
      float *d = img_dest + (y + pad_h) * iwidth + pad_w;
      
      for (int x = 0; x < width; x++)
      {
        const float l = ( rgb_luminance(RGBtosRGB(s[x*ch+0], 1.f), RGBtosRGB(s[x*ch+1], 1.f), RGBtosRGB(s[x*ch+2], 1.f)) );
        min = MIN(min, l);
        max = MAX(max, l);
        
        d[x] = l;
      }
    }
  }
  else
  {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_dest) schedule(static) reduction(min : min) reduction(max : max)
#endif
    for (int y = 0; y < height; y++)
    {
      float *s = img_src + y * width * ch;
      float *d = img_dest + (y + pad_h) * iwidth + pad_w;
      
      for (int x = 0; x < width; x++)
      {
        float XYZ[3] = {0};
        
        dt_Lab_to_XYZ(s + x*ch, XYZ);
        dt_XYZ_to_sRGB(XYZ, d + x);
        
        const float l = d[x];
        min = MIN(min, l);
        max = MAX(max, l);
      }
    }
  }
  
  *l_min = min;
  *l_max = max;
  
}

static void unpad_image(float *img_src, const int width, const int height, float *img_dest, 
    const int pad_w, const int pad_h)
{
  const int iwidth = width + pad_w*2;

  memset(img_dest, 0, width * height * sizeof(float));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_dest) schedule(static)
#endif
  for (int y = 0; y < height; y++)
  {
    float *s = img_src + (y + pad_h) * iwidth + pad_w;
    float *d = img_dest + y * width;
    
    memcpy(d, s, width * sizeof(float));
  }

}

static void adjust_brcogm(float *img_src, const int width, const int height, const int ch,
    const float brightness, const float contrast, const float gamma)
{
  const int stride = width * height * ch;
  
  if (gamma != 0.f || contrast != 0.f || brightness != 0.f)
  {
    const float gm = powf(10.f, -gamma);
    const float cn = expf(contrast / .25f);
    const float br = brightness * 2.f;
    
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src) schedule(static)
#endif
    for (int i = 0; i < stride; i++)
    {
      if (img_src[i] < 0.f)
        img_src[i] = -powf(-img_src[i], -gm);
      else
        img_src[i] = powf(img_src[i], gm);
      
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
  for (int i = 0; i < stride; i++)
  {
    min = MIN(min, img_src[i]);
    max = MAX(max, img_src[i]);
  }
  
  if (min == max) return;
  
  const float mult = (H-L) / (max - min);
  
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, min) schedule(static)
#endif
  for (int i = 0; i < stride; i++)
  {
    img_src[i] = ((img_src[i] - min) * mult) + L;
  }
  
}

static void gradient_rgb(float *img_src, float *img_grx, float *img_gry, float *img_grn, const int width, const int height, 
    const int pad_w, const int pad_h, 
    const float _oddness, const float image_scale, const float beta)
{
  const float oddness = _oddness * sqrtf(image_scale);
  
  memset(img_grx, 0, width * height * sizeof(float));
  memset(img_gry, 0, width * height * sizeof(float));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_src, img_grx, img_gry, img_grn) schedule(static)
#endif
  for (int y = 1; y < height-1; y++)
  {
    float *s0 = img_src + y * width;
    float *s1 = img_src + (y + 1) * width;
    float *dx = img_grx + y * width;
    float *dy = img_gry + y * width;
    float *dn = img_grn + y * width;
    
    for (int x = 1; x < width-1; x++)
    {
      float grx = s0[x + 1] - s0[x];
      float gry = s1[x] - s0[x];
      float grn = sqrtf(grx*grx + gry*gry);

      const float n = 1e-15 + powf(grn, oddness);
      
      dx[x] = grx / n;
      dy[x] = gry / n;
      dn[x] = grn;
    }
  }

}

static void estimate_laplacian(float *img_grx, float *img_gry, float *img_grn, float *img_dest, const int width, const int height, 
    const int pad_w, const int pad_h, 
    const float beta, const float image_scale, const float oddness)
{
  const int stride = width * height;

  float *img_bp = NULL;
  
  img_bp = dt_alloc_align(64, width * height * sizeof(float));
  if (img_bp == NULL) goto cleanup;

  memset(img_bp, 0, stride * sizeof(float));
  memset(img_dest, 0, stride * sizeof(float));

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_grx) schedule(static)
#endif
  for (int y = 0; y < height; y++)
  {
    float *s0 = img_grx + y * width;

    for (int x = width-1; x > 0; x--)
    {
      s0[x] -= s0[x-1];
    }
  }

  for (int y = height-1; y > 0; y--)
  {
    float *s0 = img_gry + y * width;
    float *s1 = img_gry + (y - 1) * width;

    for (int x = 0; x < width; x++)
    {
      s0[x] -= s1[x];
    }
  }

  
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_grx, img_gry, img_grn, img_dest) schedule(static)
#endif
  for (int i = 0; i < stride; i++)
  {
    img_dest[i] = (img_grx[i] + img_gry[i]);
  }

/*
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_grx, img_gry, f) schedule(static)
#endif
  for (int y = 0; y < height; y++)
  {
    float *dx = img_grx + y * width;
    float *dy = img_gry + y * width;
    float *d = f + y * width;
    
    for (int x = 0; x < width; x++)
    {
      d[x] = (dx[x] + dy[x]);
    }
  }
*/

/*
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_grn) schedule(static)
#endif
  for (int y = pad_h; y < height-pad_h; y++)
  {
    float *d = img_grn + y * width;
    
    for (int x = pad_w; x < width-pad_w; x++)
    {
      d[x] = 0.f;
    }
  }
*/
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_bp, img_grn) schedule(static)
#endif
  for (int y = 2; y < height-2; y++)
  {
    float *d = img_bp + y * width;
    float *bi0 = img_grn + y * width;
    float *bi1 = img_grn + (y + 1) * width;
    float *bi_1 = img_grn + (y - 1) * width;

    for (int x = 2; x < width-2; x++)
    {
      d[x] = -4.f*bi0[x] + bi0[x+1] + bi0[x-1] + bi_1[x] + bi1[x];
    }
  }

  
/*  for (int y = 0; y < height; y++)
  {
    float *d = img_dest + y * width;
    float *f1 = f + y * width;
    float *bi = img_bp + y * width;

    for (int x = 0; x < width; x++)
    {
      d[x] = (f1[x] - bi[x]);
    }
  }
*/
  

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(img_bp, img_grn, img_dest) schedule(static)
#endif
  for (int i = 0; i < stride; i++)
  {
    img_dest[i] -= (beta * img_bp[i]);
  }

cleanup:
  
  if (img_bp) dt_free_align(img_bp);
  

}

static void recontruct_laplacian(float *img_src, float *img_grn, float *img_dest, const int width, const int height, 
    const float beta, const float oddness, const float image_scale)
{
  fbw_fft_t fft_fbw = {0};
  
  fft(&fft_fbw, img_src, width, height, image_scale);
  
  const int width_fft = fft_fbw.width_fft;
  const int height_fft = fft_fbw.height_fft;
  
  const float piw = (2.f * (float)M_PI / (float)width_fft);
  const float pih = (2.f * (float)M_PI / (float)height_fft);
  
  const int w = fft_fbw.width_fft_complex;
  const int h = fft_fbw.height_fft_complex;
  fftwf_complex *out_src = fft_fbw.out_src;
  
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(out_src, img_grn) schedule(static)
#endif
  for (int y = 0; y < h; y++)
  {
    fftwf_complex *d = out_src + y * w;

    const float cos_y = cos(pih*y);

    for (int x = 0; x < w; x++)
    {
      const float cos_x = cos(piw*x);
      
      float cos_xy = 1.f;
      if (x > 0 || y > 0) cos_xy = (cos_x + cos_y - 2.f) * 2.f;
      
      d[x][0] /= cos_xy;
      d[x][1] /= cos_xy;
    }
  }

  fft_fbw.out_src[0][0] = 0.f;
  fft_fbw.out_src[0][1] = 0.f;

  ifft(&fft_fbw, img_dest);

  fft_free(&fft_fbw);

}

static void fbw_process(float *img_src, float *img_dest, const int width, const int height, const int ch, 
    const float oddness, const float brightness, const float contrast, const float gamma,
    const float beta, const float image_scale, const dt_iop_colorspace_type_t cst)
{
  
  float *img_padded = NULL;
  float *img_grx = NULL;
  float *img_gry = NULL;
  float *img_grn = NULL;
  
  const int pad_w = MIN(width, height)/2;
  const int pad_h = pad_w;
  
  const int iwidth = width + pad_w*2;
  const int iheight = height + pad_h*2;
  
  float l_min = 0.f;
  float l_max = 1.f;
  
  img_grx = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_grx == NULL) goto cleanup;
  
  img_gry = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_gry == NULL) goto cleanup;
  
  img_grn = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_grn == NULL) goto cleanup;
  
  img_padded = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_padded == NULL) goto cleanup;
  
  pad_image(img_src, width, height, ch, img_padded, pad_w, pad_h, &l_min, &l_max, cst, beta);

  gradient_rgb(img_padded, img_grx, img_gry, img_grn, iwidth, iheight, pad_w, pad_h, oddness, image_scale, beta);

  float *img_l = img_padded;
  
  estimate_laplacian(img_grx, img_gry, img_grn, img_l, iwidth, iheight, pad_w, pad_h, beta, image_scale, oddness);

  recontruct_laplacian(img_l, img_grn, img_l, iwidth, iheight, beta, oddness, image_scale);

  float *img_unpad_gr = img_grx;
  
  unpad_image(img_l, width, height, img_unpad_gr, pad_w, pad_h);
  
  normalize(img_unpad_gr, width, height, l_min, l_max);

  adjust_brcogm(img_unpad_gr, width, height, 1, brightness, contrast, gamma);
  
  image_to_output(img_unpad_gr, width, height, ch, img_dest, cst);

cleanup:

  if (img_grx) dt_free_align(img_grx);
  if (img_gry) dt_free_align(img_gry);
  if (img_grn) dt_free_align(img_grn);
  if (img_padded) dt_free_align(img_padded);
  
}

void process_internal(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid, void *const ovoid,
    const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out, const int use_sse)
{
  const dt_iop_fbw_data_t *const p = (const dt_iop_fbw_data_t *const)piece->data;
  const float image_scale = roi_in->scale / piece->iscale;
  const dt_iop_colorspace_type_t cst = dt_iop_module_colorspace(self);

  fbw_process((float *)ivoid, (float *)ovoid, roi_in->width, roi_in->height, piece->colors, 
                  p->oddness/100.f, p->brightness/100.f, p->contrast/100.f, p->gamma/100.f,
                  p->beta, image_scale, cst);

  if(piece->pipe->mask_display & DT_DEV_PIXELPIPE_DISPLAY_MASK) dt_iop_alpha_copy(ivoid, ovoid, roi_out->width, roi_out->height);
  
    check_nan((float *)ovoid, roi_in->width * roi_in->height * piece->colors, "process_internal");
           
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
