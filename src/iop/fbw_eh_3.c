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

#ifdef HAVE_FFTW3
#include <fftw3.h>
#endif

#ifdef HAVE_FFTW3_OMP
#include <fftw3.h>
#endif


DT_MODULE_INTROSPECTION(1, dt_iop_fbw_params_t)


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

} fbw_fftw3_t;

typedef struct dt_iop_fbw_params_t
{
  float oddness;
  float beta;
  float brightness;
  float contrast;
  float gamma;
} dt_iop_fbw_params_t;

typedef struct dt_iop_fbw_gui_data_t
{
  GtkWidget *sl_oddness;
  GtkWidget *sl_brightness;
  GtkWidget *sl_contrast;
  GtkWidget *sl_gamma;
  GtkWidget *sl_beta;
  
  float l_min;
  float l_max;
  uint64_t hash;
  dt_pthread_mutex_t lock;
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

  dt_pthread_mutex_lock(&g->lock);
  g->l_min = NAN;
  g->l_max = NAN;
  g->hash = 0;
  dt_pthread_mutex_unlock(&g->lock);

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

  dt_pthread_mutex_init(&g->lock, NULL);
  g->l_min = NAN;
  g->l_max = NAN;
  g->hash = 0;

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
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;
  dt_pthread_mutex_destroy(&g->lock);

  free(self->gui_data);
  self->gui_data = NULL;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

//#define RNG_UP_RGB 255.f
#define RNG_UP_RGB 1.f

static void check_nan(const float* im, const int size)
{
  int i_nan = 0;

  for (int i = 0; i < size; i++) if ( isnan(im[i]) ) i_nan++;

  if (i_nan > 0) printf("freaky b&w nan: %i\n", i_nan);
}

static void print_min_max(float *img_src, const int width, const int height, const int ch, float *_min, float *_max, const char *str)
{
/*  const int ch1 = (ch==4) ? 3: ch;
  
  float min = INFINITY;
  float max = -INFINITY;
  const int stride = width * height * ch;
  
  for (int i = 0; i < stride; i+=ch)
  {
    for (int c = 0; c < ch1; c++)
    {
      min = MIN(min, img_src[i+c]);
      max = MAX(max, img_src[i+c]);
    }
  }
  
  if (_min) *_min = min;
  if (_max) *_max = max;
  
  printf("%s: min=%f, max=%f\n", str, min, max);
  */
}

//-----------------------------------------------------------------------

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
  
static void fft(fbw_fftw3_t *fftw3_fbw, float *image_src, const int width, const int height, const float image_scale)
{
  int FFTW_FACTORS[7] = {13,11,7,5,3,2,0}; // end with zero to detect the end of the array
  
  fftw3_fbw->width_src = width;
  fftw3_fbw->height_src = height;

  fftw3_fbw->width_dest = width;
  fftw3_fbw->height_dest = height;

  fftw3_fbw->width_fft = width;
  fftw3_fbw->height_fft = height;

if(0)  fftw3_fbw->width_fft = fftw3_find_closest_factor(fftw3_fbw->width_fft, FFTW_FACTORS);
if(0)  fftw3_fbw->height_fft = fftw3_find_closest_factor(fftw3_fbw->height_fft, FFTW_FACTORS);

//  fftw3_fbw->width_fft /= image_scale;
//  fftw3_fbw->height_fft /= image_scale;

  fftw3_fbw->width_fft_complex = fftw3_fbw->width_fft/2+1;
  fftw3_fbw->height_fft_complex = fftw3_fbw->height_fft;

  dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
  
  fftw3_fbw->in_src = (float*)fftwf_malloc(sizeof(float) * fftw3_fbw->width_fft * fftw3_fbw->height_fft);
  fftw3_fbw->out_src = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * fftw3_fbw->width_fft_complex * fftw3_fbw->height_fft_complex);
  
  fftw3_fbw->plan_src = fftwf_plan_dft_r2c_2d(fftw3_fbw->height_fft, fftw3_fbw->width_fft, fftw3_fbw->in_src, fftw3_fbw->out_src, FFTW_ESTIMATE);
  fftw3_fbw->plan_inv = fftwf_plan_dft_c2r_2d(fftw3_fbw->height_fft, fftw3_fbw->width_fft, fftw3_fbw->out_src, fftw3_fbw->in_src, FFTW_ESTIMATE);

  dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);

  memset(fftw3_fbw->in_src, 0, sizeof(float) * fftw3_fbw->width_fft * fftw3_fbw->height_fft);
  memset(fftw3_fbw->out_src, 0, sizeof(fftwf_complex) * fftw3_fbw->width_fft_complex * fftw3_fbw->height_fft_complex);

  for (int y = 0; y < fftw3_fbw->height_src; y++)
  {
    float *in_src = fftw3_fbw->in_src + y * fftw3_fbw->width_fft;
    const float *src = image_src + y * fftw3_fbw->width_src;

    for (int x = 0; x < fftw3_fbw->width_src; x++)
    {
      in_src[x] = src[x];
    }
  }

//  dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
  
  fftwf_execute(fftw3_fbw->plan_src);
  
//  dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);

}

static void ifft(fbw_fftw3_t *fftw3_fbw, float *image_dest)
{
  const float scale = 1.0 / (fftw3_fbw->width_fft * fftw3_fbw->height_fft);
  
  memset(fftw3_fbw->in_src, 0, sizeof(float) * fftw3_fbw->width_fft * fftw3_fbw->height_fft);

//  dt_pthread_mutex_lock(&darktable.plugin_threadsafe);
  
  fftwf_execute(fftw3_fbw->plan_inv);
  
//  dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);

  for (int y = 0; y < fftw3_fbw->height_dest; y++)
  {
    float *dest = image_dest + y * fftw3_fbw->width_dest;
    float *out_inv = fftw3_fbw->in_src + y * fftw3_fbw->width_fft;

    for (int x = 0; x < fftw3_fbw->width_dest; x++)
    {
      dest[x] = out_inv[x] * scale;
    }
  }

  dt_pthread_mutex_lock(&darktable.plugin_threadsafe);

  if (fftw3_fbw->plan_src) fftwf_destroy_plan(fftw3_fbw->plan_src);
  if (fftw3_fbw->plan_inv) fftwf_destroy_plan(fftw3_fbw->plan_inv);
  
  if (fftw3_fbw->in_src) fftwf_free(fftw3_fbw->in_src);
  if (fftw3_fbw->out_src) fftwf_free(fftw3_fbw->out_src);

  dt_pthread_mutex_unlock(&darktable.plugin_threadsafe);
  
}

static inline float sRGBtoRGB(float sval, const float beta)
{
  return (sval<=0.04045f?sval/(12.92f):powf((sval + 0.055f)/(1.055f),2.4f));
}

static inline float RGBtosRGB(float val, const float beta)
{
  return (val<=0.0031308f?val*12.92f:1.055f*powf(val,0.416667f) - 0.055f);
}
/*
static inline float sRGBtoRGB(float sval, const float beta)
{
  return (sval);
}

static inline float RGBtosRGB(float val, const float beta)
{
  return (val);
}
*/
static void image_srgb_to_rgb_bw(float *img_src, const int width, const int height, float *img_dest, const float beta)
{

  const int stride = width * height;
  for (int i = 0; i < stride; i++)
  {
    img_dest[i*4] = img_dest[i*4+1] = img_dest[i*4+2] = sRGBtoRGB(img_src[i] * (1.f / RNG_UP_RGB), 1.f);
  }

}

static void symetrize_image(
  float *img_dest
, float *img_src
, const int width
, const int height
, const int ch
, const unsigned p_borderSize_w
, const unsigned p_borderSize_h
, const float _beta
){
  const float beta = 1;
  
  const unsigned ws = width;
  const unsigned wd = width;
  const unsigned hs = height;
  const unsigned hd = height;

  //! Top and bottom
  for (unsigned y = 0; y < p_borderSize_h; y++)
  {
    float *im2_td = img_dest + y * wd * ch;
    float *im2_ts = img_src + (2 * p_borderSize_h - y - 1) * ws * ch;
    float *im2_bd = img_dest + (hd - y - 1) * wd * ch;
    float *im2_bs = img_src + (hs - 2 * p_borderSize_h + y) * ws * ch;

    for (unsigned x = 0; x < wd*ch; x++)
    {
      im2_td[x] = im2_ts[x] * beta;
      im2_bd[x] = im2_bs[x] * beta;
    }
  }

  //! Right and left
  for (unsigned y = 0; y < hd; y++)
  {
    float *im2_rd = img_dest + y * wd * ch;
    float *im2_rs = img_src + y * ws * ch;
    float *im2_ld = img_dest + y * wd * ch;
    float *im2_ls = img_dest + y * wd * ch;

    for (unsigned x = 0; x < p_borderSize_w; x++)
    {
      for (int c = 0; c < ch; c++)
      {
        im2_rd[(x + wd - p_borderSize_w)*ch + c] = im2_rs[(ws - p_borderSize_w - 1 - x)*ch + c] * beta;
        im2_ld[x*ch + c] = im2_ls[(p_borderSize_w * 2 - 1 - x)*ch + c] * beta;
      }
    }
  }

}

static void image_srgb_to_lab_bw(float *img_src, const int width, const int height, float *img_dest)
{

  const int stride = width * height;
  for (int i = 0; i < stride; i++)
  {
    float sRGB[3] = { img_src[i] * (1.f / RNG_UP_RGB), img_src[i] * (1.f / RNG_UP_RGB), img_src[i] * (1.f / RNG_UP_RGB) };
    float XYZ[3] = {0};
    
    dt_sRGB_to_XYZ(sRGB, XYZ);
    dt_XYZ_to_Lab(XYZ, img_dest + i*4);
  }

}

static void image_rgb_to_srgb_padded(float *img_src, const int width, const int height, float *img_dest, 
    const int pad_w, const int pad_h, 
                                      const float beta)
{
  const int iwidth = width + pad_w*2;
  const int iheight = height + pad_h*2;

  memset(img_dest, 0, iwidth * iheight * 3 * sizeof(float));

  for (int y = 0; y < height; y++)
  {
    float *s = img_src + y * width * 4;
    float *d = img_dest + (y + pad_h) * iwidth * 3 + pad_w * 3;
    
    for (int x = 0; x < width; x++)
    {
      for (int c = 0; c < 3; c++)
      {
        d[x*3+c] = RGBtosRGB(s[x*4+c], 1.f) * RNG_UP_RGB;
      }
    }
  }

if(0)  symetrize_image(img_dest, img_dest, iwidth, iheight, 3, pad_w, pad_h, beta);
  
}

static void image_lab_to_srgb_padded(float *img_src, const int width, const int height, float *img_dest, 
    const int pad_w, const int pad_h,
                                      const float beta)
{
  const int iwidth = width + pad_w*2;
  const int iheight = height + pad_h*2;

  memset(img_dest, 0, iwidth * iheight * 3 * sizeof(float));

  for (int y = 0; y < height; y++)
  {
    float *s = img_src + y * width * 4;
    float *d = img_dest + (y + pad_h) * iwidth * 3 + pad_w * 3;
    
    for (int x = 0; x < width; x++)
    {
      float XYZ[3] = {0};
      
      dt_Lab_to_XYZ(s + x*4, XYZ);
      dt_XYZ_to_sRGB(XYZ, d + x*3);
      
      for (int c = 0; c < 3; c++)
      {
        d[x*3+c] *= RNG_UP_RGB;
      }
    }
  }

}

static void unpad_image(float *img_src, const int width, const int height, float *img_dest, 
    const int pad_w, const int pad_h)
{
  const int iwidth = width + pad_w*2;

  memset(img_dest, 0, width * height * sizeof(float));

  for (int y = 0; y < height; y++)
  {
    float *s = img_src + (y + pad_h) * iwidth + pad_w;
    float *d = img_dest + y * width;
    
    memcpy(d, s, width * sizeof(float));
  }

}

static void normalize(float *img_src, const int width, const int height, const float L, const float H)
{
  
  float min = INFINITY;
  float max = -INFINITY;
  const int stride = width * height;
  
  for (int i = 0; i < stride; i++)
  {
    min = MIN(min, img_src[i]);
    max = MAX(max, img_src[i]);
  }
  
  if (min == max) return;
  
  const float mult = (H-L) / (max - min);
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] = ((img_src[i] - min) * mult) + L;
  }
  
}

static void gradient_rgb(float *img_src, float *img_dest_x, float *img_dest_y, float *img_dest_n, const int width, const int height,
    const int pad_w, const int pad_h,
    const float oddness, const float image_scale, const float beta)
{
  const int ch = 3;
//  float norm = 0.f;
//  float scale = beta;//sqrtf(image_scale);
  
  memset(img_dest_x, 0, width * height * sizeof(float));
  memset(img_dest_y, 0, width * height * sizeof(float));

  for (int y = 0; y < height-1; y++)
  {
    float *s0 = img_src + y * width * ch;
    float *s1 = img_src + (y + 1) * width * ch;
    float *dx = img_dest_x + y * width;
    float *dy = img_dest_y + y * width;
    float *dn = img_dest_n + y * width;
    
    for (int x = 0; x < width-1; x++)
    {
      float RGBx[3] = {0};
      float RGBy[3] = {0};
      float RGBn[3] = {0};

      for (int c = 0; c < 3; c++)
      {
        RGBx[c] = s0[(x+1)*ch + c] - s0[x*ch + c];
        RGBy[c] = s1[x*ch + c] - s0[x*ch + c];
        
//        RGBx[c] *= 1./((width-pad_w)*(height-pad_h));
//        RGBy[c] *= 1./((width-pad_w)*(height-pad_h));
        
        RGBn[c] = (RGBx[c]*RGBx[c] + RGBy[c]*RGBy[c]);
      }
      
      int max_bn = 0;
      if (RGBn[0] > RGBn[1])
        max_bn = (RGBn[0] > RGBn[2]) ? 0: 2;
      else
        max_bn = (RGBn[1] > RGBn[2]) ? 1: 2;
      
//      if (x >= pad_w-1 && x < width-pad_w-1 && y >= pad_h-1 && y < height-pad_h-1)
//      norm += sqrtf(RGBn[max_bn]);

      float n;
      
//      float n = 1e-5 + powf(RGBn[max_bn], oddness/(image_scale));
//      if (x >= pad_w-1 && x < width-pad_w-1 && y >= pad_h-1 && y < height-pad_h-1)
        n = 1e-15 + powf(RGBn[max_bn], oddness /* *sqrtf(image_scale) */ );
//      else
//        n = 1.f;
      
      dx[x] = RGBx[max_bn] / n;
      dy[x] = RGBy[max_bn] / n;
      dn[x] = RGBn[max_bn];
//      dn[x] = n;
    }
  }

/*  for (int y = 0; y < height-1; y++)
  {
    float *dx = img_dest_x + y * width;
    float *dy = img_dest_y + y * width;
    
    for (int x = 0; x < width-1; x++)
    {
//      if (x >= pad_w-1 && x < width-pad_w-1 && y >= pad_h-1 && y < height-pad_h-1)
      {
      dx[x] *= norm/image_scale;
      dy[x] *= norm/image_scale;
      }
    }
  }
  */
}

static void estimate_laplacian(float *img_src_x, float *img_src_y, float *img_src_n, float *img_dest, const int width, const int height, 
    const int _pad_w, const int _pad_h, 
    const float _beta, const float image_scale)
{
  const int stride = width * height;

  memset(img_dest, 0, stride * sizeof(float));

  for (int y = 0; y < height; y++)
  {
    float *s0 = img_src_x + y * width;

    for (int x = width-1; x > 0; x--)
    {
      s0[x] -= s0[x-1];
//      s0[x] -= s0[x-1] * sqrtf( (float)((width-x)*(height-y)) / (float)(stride) );
    }
  }

  for (int y = height-1; y > 0; y--)
  {
    float *s0 = img_src_y + y * width;
    float *s1 = img_src_y + (y - 1) * width;

    for (int x = 0; x < width; x++)
    {
      s0[x] -= s1[x];
//      s0[x] -= s1[x] * sqrtf( (float)((width-x)*(height-y)) / (float)(stride) );
    }
  }
/*
  const float beta = 1.f + (1.f - (image_scale)) / 100.f;
//  const float beta = _beta;
  
  const int pad_w = _pad_w;
  const int pad_h = _pad_h;
  
  float *img_dest_x = img_src_x;
  float *img_dest_y = img_src_y;

  for (int y = 0; y < pad_h; y++)
  {
    float *dx = img_dest_x + y * width;
    float *dy = img_dest_y + y * width;
    
    for (int x = 0; x < width; x++)
    {
      dx[x] *= beta;
      dy[x] *= beta;
    }
  }

  for (int y = height-pad_h; y < height; y++)
  {
    float *dx = img_dest_x + y * width;
    float *dy = img_dest_y + y * width;
    
    for (int x = 0; x < width; x++)
    {
      dx[x] *= beta;
      dy[x] *= beta;
    }
  }

  for (int y = 0; y < height; y++)
  {
    float *dx = img_dest_x + y * width;
    float *dy = img_dest_y + y * width;
    
    for (int x = 0; x < pad_w; x++)
    {
      dx[x] *= beta;
      dy[x] *= beta;
    }
  }

  for (int y = 0; y < height; y++)
  {
    float *dx = img_dest_x + y * width;
    float *dy = img_dest_y + y * width;
    
    for (int x = width-pad_w; x < width; x++)
    {
      dx[x] *= beta;
      dy[x] *= beta;
    }
  }
*/
  
/*  for (int i = 0; i < stride; i++)
  {
    img_dest[i] = (img_src_x[i] + img_src_y[i]);
  }
*/
  for (int y = 0; y < height; y++)
  {
    float *dx = img_src_x + y * width;
    float *dy = img_src_y + y * width;
//    float *dn = img_src_n + y * width;
    float *d = img_dest + y * width;
    
    for (int x = 0; x < width; x++)
    {
      d[x] = (dx[x] + dy[x]);
//      d[x] = (dx[x] + dy[x]) / (dn[x]/(1.f+(float)(x*y)));
//      d[x] = (dx[x] + dy[x]) * ( 1e-15 + (float)(((float)(((float)x*(float)width)+((float)y*(float)height)))   /*/(float)(height*width)*/)  );
//      d[x] = (dx[x] * ( (float)(width-x) ) + dy[x] * (  (float)(height-y) ) );
      
//      d[x] *= ( 1e-15 + (float)(((float)(((float)x*(float)width)+((float)y*(float)height)))   /*/(float)(height*width)*/)  );
    }
  }
  
//  symetrize_image(img_dest, img_dest, width, height, 1, pad_w, pad_h, beta);
  

}

static void recontruct_laplacian(float *img_src, float *img_grn, float *img_dest, const int width, const int height, 
    const float beta, const float oddness, const float image_scale)
{
  fbw_fftw3_t fftw3_fbw = {0};
//  fbw_fftw3_t fftw3_fbw_grn = {0};
  
  fft(&fftw3_fbw, img_src, width, height, image_scale);
//  fft(&fftw3_fbw_grn, img_grn, width, height, image_scale);
  
  const int width_fft = fftw3_fbw.width_fft;
  const int height_fft = fftw3_fbw.height_fft;
  
  const float piw = (2.f * (float)M_PI / (float)width_fft);
  const float pih = (2.f * (float)M_PI / (float)height_fft);
  
  for (int y = 0; y < fftw3_fbw.height_fft_complex; y++)
  {
    fftwf_complex *d = fftw3_fbw.out_src + y * fftw3_fbw.width_fft_complex;
//    fftwf_complex *grn = fftw3_fbw_grn.out_src + y * fftw3_fbw_grn.width_fft_complex;
//    float *grn = img_grn + y * width;

    const float cos_y = cos(pih*y);

    for (int x = 0; x < fftw3_fbw.width_fft_complex; x++)
    {
      const float cos_x = cos(piw*x);
      
      float cos_xy = 1.f;
      if (x > 0 || y > 0) cos_xy = (cos_x + cos_y - 2.f) * 2.f;
      
//      cos_xy *= powf(powf(image_scale*image_scale, oddness), grn[x]);
//      cos_xy *= powf((image_scale*image_scale), grn[x]);
      
      d[x][0] /= cos_xy;
      d[x][1] /= cos_xy;
      
//      d[x][0] /= cos_xy / grn[x][0];
//      d[x][1] /= cos_xy / grn[x][1];
    }
  }

  fftw3_fbw.out_src[0][0] = 0.f;
  fftw3_fbw.out_src[0][1] = 0.f;

//  ifft(&fftw3_fbw_grn, img_dest);
  ifft(&fftw3_fbw, img_dest);

}

static void adjust_brcogm2(float *img_src, float *img_grn, const int width, const int height, const int ch, 
    const int pad_w, const int pad_h,
    const float brightness, const float contrast, const float gamma,
    const float beta, const float oddness, const float image_scale)
{
//  const int stride = width * height * ch;
//  const float image_scale = 1.f - _image_scale;
  /*
//  const float gm = powf(10.f, -gamma);
  
  const float gm = powf(10.f, -sqrtf(image_scale) / 10.f); // 6 .06
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] = powf(img_src[i], gm);
  }
  */
//  const float cn = expf(contrast / .25f);
/*  
  const float cn = expf(image_scale * image_scale * .25f * .25f); // 12  .12
//  const float cn = (image_scale); // 12  .12
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] = (img_src[i] - .5f);
    img_src[i] *= cn;
    img_src[i] = (img_src[i] + .5f);
  }
  
//  const float br = brightness * 2.f;
  
  const float br = image_scale * image_scale * .25f * .25f;
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] += br;
  }
  
  
  const float cn = expf(beta);
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] *= cn;
  }
  */
  
  const int iwidth = width + pad_w*2;
  
  for (int y = 0; y < height; y++)
  {
    float *s = img_src + y * width * ch;
    float *d = img_grn + (y + pad_h) * iwidth + pad_w;
    
    for (int x = 0; x < width; x++)
    {
      for (int c = 0; c < ch; c++)
//        s[x*ch+c] *= powf(d[x], (1.f-image_scale)*oddness);
//        s[x*ch+c] /= powf(d[x], -beta);
      s[x*ch+c] += (d[x]*image_scale*image_scale);
    }
  }
  
}

static void adjust_brcogm(float *img_src, const int width, const int height, const int ch,
    const float brightness, const float contrast, const float gamma)
{
  const int stride = width * height * ch;
  
  const float gm = powf(10.f, -gamma);
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] = powf(img_src[i], gm);
  }
  
  const float cn = expf(contrast / .25f);
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] = (img_src[i] - .5f);
    img_src[i] *= cn;
    img_src[i] = (img_src[i] + .5f);
  }
  
  const float br = brightness * 2.f;
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] += br;
  }
  
}

static inline float rgb_luminance(const float r, const float g, const float b)
{
  return (r * 0.22248840f + g * 0.71690369f + b * 0.06060791f);
}

static void luminance_minmax(float *img_src, const int width, const int height,
    float *l_min, float *l_max,
    const dt_iop_colorspace_type_t cst, const float beta)
{
  const int stride = width * height * 4;
  
  float min = INFINITY;
  float max = -INFINITY;
      
  if (cst == iop_cs_rgb)
  {
    for (int i = 0; i < stride; i+=4)
    {
      const float l = RGBtosRGB( rgb_luminance(img_src[i+0], img_src[i+1], img_src[i+2]), 1.f ) * RNG_UP_RGB;
      min = MIN(min, l);
      max = MAX(max, l);
    }
  }
  else
  {
    for (int i = 0; i < stride; i+=4)
    {
      const float l = img_src[i];
      min = MIN(min, l);
      max = MAX(max, l);
    }
  }
  
  *l_min = min;
  *l_max = max;
}

static void fbw_process(float *img_src, float *img_dest, const int width, const int height, const int ch, 
    const float oddness, const float brightness, const float contrast, const float gamma,
    const float beta, const float image_scale, const dt_iop_colorspace_type_t cst,
    float l_min, float l_max)
{
  
  float *img_padded = NULL;
  float *img_grx = NULL;
  float *img_gry = NULL;
  float *img_grn = NULL;
  
/*  int tmp = MAX(width, height);
  int n = 1;
  while (powf(2, n) < tmp) n++;
  tmp = powf(2, n);
  
  const int pad_w = (tmp-width)/2;
  const int pad_h = (tmp-height)/2;*/
  const int pad_w = MIN(width, height)/2;
  const int pad_h = pad_w;
//  const int pad_w = width;
//  const int pad_h = height;
//  const int pad_w = 1;
//  const int pad_h = 1;
  
  const int iwidth = width + pad_w*2;
  const int iheight = height + pad_h*2;
  
  float min = 0.f;
  float max = 0.f;
  
  img_grx = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_grx == NULL) goto cleanup;
  
  img_gry = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_gry == NULL) goto cleanup;
  
  img_grn = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_grn == NULL) goto cleanup;
  
  img_padded = dt_alloc_align(64, iwidth * iheight * 3 * sizeof(float));
  if (img_padded == NULL) goto cleanup;
  
  printf("\nimage_scale=%f\n", image_scale);
  print_min_max(img_src, width, height, 4, NULL, NULL, "ivoid");
  
  luminance_minmax(img_src, width, height, &l_min, &l_max, cst, beta);
  
  if (cst == iop_cs_rgb)
    image_rgb_to_srgb_padded(img_src, width, height, img_padded, pad_w, pad_h, beta);
  else
    image_lab_to_srgb_padded(img_src, width, height, img_padded, pad_w, pad_h, beta);
  
/*  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      for (int c = 0; c < 3; c++)
      {
        int y1 = y + iheight - height;
        int x1 = x + iwidth - width;
        
        img_dest[y*width*ch + x*ch + c] = img_padded[y1*iwidth*3 + x1*3 + c];
      }
    }
  }
  goto cleanup;
  */
  printf("l_min=%f, l_max=%f\n", l_min, l_max);
  
  print_min_max(img_padded, iwidth, iheight, 3, NULL, NULL, "img_padded");
  
/*
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
*/

  gradient_rgb(img_padded, img_grx, img_gry, img_grn, iwidth, iheight, pad_w, pad_h, oddness, image_scale, beta);

  print_min_max(img_grx, iwidth, iheight, 1, &min, &max, "img_grx");
//  normalize(img_grx, iwidth, iheight, min*(1./iwidth), max*(1./iwidth));
  print_min_max(img_gry, iwidth, iheight, 1, &min, &max, "img_gry");
//  normalize(img_gry, iwidth, iheight, min*(1./iheight), max*(1./iheight));
 
/*
  -luminance[0] ia={0,ia}
*/
  
/*
  # Estimate laplacian of final image.
  -s. c
  -f.. "i - i(x-1,y,0,0)"
  -f. "i - i(x,y-1,0,0)"
  -+[-2,-1]
*/
  float *img_l = img_padded;
  
  
  estimate_laplacian(img_grx, img_gry, img_grn, img_l, iwidth, iheight, pad_w, pad_h, beta, image_scale);

  print_min_max(img_grx, iwidth, iheight, 1, NULL, NULL, "img_grx 2");
  print_min_max(img_gry, iwidth, iheight, 1, NULL, NULL, "img_gry 2");
  print_min_max(img_l, iwidth, iheight, 1, &min, &max, "img_l");
  
//  normalize(img_l, iwidth, iheight, min*image_scale, max*image_scale);
  
/*
  # Reconstruct image from laplacian.
  -fft.
  100%,100%,1,1,'cos(2*x*pi/w)+cos(2*y*pi/h)' -*. 2 --. 4
  -=. 1 -/[-3,-2] . -rm.
  -=.. 0 -=. 0
  -ifft[-2,-1] -rm.
*/
if(1)  recontruct_laplacian(img_l, img_grn, img_l, iwidth, iheight, beta, oddness, image_scale);

/*  const int stride = iwidth * iheight;
  for (int i = 0; i < stride; i++)
    img_l[i] *= image_scale;
*/
  print_min_max(img_l, iwidth, iheight, 1, NULL, NULL, "recontruct_laplacian");

/*
  -shrink_xy. 1 -+. $ia -n. 0,255
*/
  float *img_unpad_gr = img_grx;
  
//  img_l = img_grn;
//  float *img_unpad_gr = img_gry;
  
  unpad_image(img_l, width, height, img_unpad_gr, pad_w, pad_h);
  
  print_min_max(img_unpad_gr, width, height, 1, NULL, NULL, "unpad_image");

/*  printf("avg_l=%f\n", avg_l);
  printf("beta=%f\n", beta);

  const int stride = width * height;
  for (int i = 0; i < stride; i++)
    img_unpad_gr[i] += avg_l;
*/
  
  if(0)  adjust_brcogm2(img_unpad_gr, img_grn, width, height, 1, pad_w, pad_h, 0, 0, 0, beta, oddness, image_scale);

if(0)
  normalize(img_unpad_gr, width, height, 0.f, RNG_UP_RGB);
else
  normalize(img_unpad_gr, width, height, l_min, l_max);

  print_min_max(img_unpad_gr, width, height, 1, NULL, NULL, "normalize");

/*
  # Merge result with original color image.
  -j[0] [1],0,0,0,0,{$1%} -rm.
  -adjust_colors ${3-5}
  -endl -a c -endl -done
*/

  if(0)  adjust_brcogm2(img_unpad_gr, img_grn, width, height, 1, pad_w, pad_h, 0, 0, 0, beta, oddness, image_scale);

  if(1)  adjust_brcogm(img_unpad_gr, width, height, 1, brightness, contrast, gamma);
  
  if (cst == iop_cs_rgb)
    image_srgb_to_rgb_bw(img_unpad_gr, width, height, img_dest, beta);
  else
    image_srgb_to_lab_bw(img_unpad_gr, width, height, img_dest);

  print_min_max(img_dest, width, height, 4, NULL, NULL, "ovoid");
  
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
  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;

  float l_min;
  float l_max;
  float tmp_l_min = NAN;
  float tmp_l_max = NAN;

  if(self->dev->gui_attached && g && piece->pipe->type == DT_DEV_PIXELPIPE_FULL)
  {
    dt_pthread_mutex_lock(&g->lock);
    const uint64_t hash = g->hash;
    dt_pthread_mutex_unlock(&g->lock);

    // note that the case 'hash == 0' on first invocation in a session implies that g->l_min / l_max
    // is NAN which initiates special handling below to avoid inconsistent results. in all
    // other cases we make sure that the preview pipe has left us with proper readings for
    // l_min / l_max. if data are not yet there we need to wait (with timeout).
    if(hash != 0 && !dt_dev_sync_pixelpipe_hash(self->dev, piece->pipe, 0, self->priority, &g->lock, &g->hash))
      dt_control_log(_("freaky b&w: inconsistent output"));

    dt_pthread_mutex_lock(&g->lock);
    tmp_l_min = g->l_min;
    tmp_l_max = g->l_max;
    dt_pthread_mutex_unlock(&g->lock);
  }

  // in all other cases we calculate l_mim / l_max here
  if(isnan(tmp_l_min) || isnan(tmp_l_max))
  {
    luminance_minmax((float *)ivoid, roi_in->width, roi_in->height, &l_min, &l_max, cst, p->beta);
  }
  else
  {
    l_min = tmp_l_min;
    l_max = tmp_l_max;
  }

  // PREVIEW pixelpipe stores l_min / l_max
  if(self->dev->gui_attached && g && piece->pipe->type == DT_DEV_PIXELPIPE_PREVIEW)
  {
    uint64_t hash = dt_dev_hash_plus(self->dev, piece->pipe, 0, self->priority);
    dt_pthread_mutex_lock(&g->lock);
    g->l_min = l_min;
    g->l_max = l_max;
    g->hash = hash;
    dt_pthread_mutex_unlock(&g->lock);
  }

  fbw_process((float *)ivoid, (float *)ovoid, roi_in->width, roi_in->height, piece->colors, 
                  p->oddness/100.f, p->brightness/100.f, p->contrast/100.f, p->gamma/100.f,
                  p->beta, image_scale, cst, l_min, l_max);

  check_nan((float *)ovoid, roi_in->width * roi_in->height * piece->colors);
  
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
