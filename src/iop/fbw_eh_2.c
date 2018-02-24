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
  
//  float avg_l;
//  uint64_t hash;
//  dt_pthread_mutex_t lock;
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
  
/*  dt_pthread_mutex_lock(&g->lock);
  g->avg_l = NAN;
  g->hash = 0;
  dt_pthread_mutex_unlock(&g->lock);
*/  
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
//  module->priority = 161; // module order created by iop_dependencies.py, do not edit! // from exposure
//  module->priority = 955; // module order created by iop_dependencies.py, do not edit! // from borders
//  module->priority = 147; // module order created by iop_dependencies.py, do not edit!
  module->priority = 632; // module order created by iop_dependencies.py, do not edit! // from monochrome
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

/*  dt_pthread_mutex_init(&g->lock, NULL);
  g->avg_l = NAN;
  g->hash = 0;
*/
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

  g->sl_beta = dt_bauhaus_slider_new_with_range(self, 0.0, 10000.0, 1.0, 0.0, 6);
  dt_bauhaus_widget_set_label(g->sl_beta, _("beta"), _("beta"));
  g_object_set(g->sl_beta, "tooltip-text", _("beta."), (char *)NULL);
  g_signal_connect(G_OBJECT(g->sl_beta), "value-changed", G_CALLBACK(beta_callback), self);

  gtk_box_pack_start(GTK_BOX(self->widget), g->sl_beta, TRUE, TRUE, 0);


}

void gui_cleanup(struct dt_iop_module_t *self)
{
//  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;
//  dt_pthread_mutex_destroy(&g->lock);

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

static void print_min_max(float *img_src, const int width, const int height, const int ch, const char *str)
{
  const int ch1 = (ch==4) ? 3: ch;
  
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
  
  printf("%s: min=%f, max=%f\n", str, min, max);
  
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
/*
static inline float sRGBtoRGB(const float sval)
{
  return (sval<=0.04045f?sval/12.92f:powf((sval + 0.055f)/(1.055f),2.4f));
}

static inline float RGBtosRGB(const float val)
{
  return (val<=0.0031308f?val*12.92f:1.055f*powf(val,0.416667f) - 0.055f);
}
*/
static inline float sRGBtoRGB(const float sval)
{
  return (sval);
}

static inline float RGBtosRGB(const float val)
{
  return (val);
}

static void image_srgb_to_rgb_bw(float *img_src, const int width, const int height, float *img_dest)
{

  const int stride = width * height;
  for (int i = 0; i < stride; i++)
  {
    img_dest[i*4] = img_dest[i*4+1] = img_dest[i*4+2] = sRGBtoRGB(img_src[i] * (1.f / RNG_UP_RGB));
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
/*
static void image_rgb_to_srgb_padded(float *img_src, const int width, const int height, float *img_dest, const int pad, 
                                      const float beta)
{
  const int ch_src = 4;
  const int ch_dest = 3;
  
  const int iwidth = width + pad;
  const int iheight = height + pad;

//  memset(img_dest, 0, width_dest * height_dest * ch * sizeof(float));

  // copy image
  for (int y = 0; y < height; y++)
  {
    float *s = img_src + y * width * ch_src;
    float *d = img_dest + (y + pad) * iwidth * ch_dest + pad * ch_dest;
    
    for (int x = 0; x < width; x++)
    {
      for (int c = 0; c < ch_dest; c++)
      {
        d[x*ch_dest+c] = RGBtosRGB(s[x*ch_src+c]) * RNG_UP_RGB;
      }
    }
  }
  
  // pad left and right
  for (int y = 0; y < iheight; y++)
  {
    const float *const src_l = img_dest + y * iwidth * ch_dest;
    float *dest_l = img_dest + y * iwidth * ch_dest;
    
    const float *const src_r = img_dest + y * iwidth * ch_dest + (iwidth-pad-1) * ch_dest;
    float *dest_r = img_dest + y * iwidth * ch_dest + (iwidth-pad) * ch_dest;
    
    for (int x = 0; x < pad; x++)
    {
      for (int c = 0; c < ch_dest; c++)
      {
        dest_l[x*ch_dest + c] = src_l[c];
        dest_r[x*ch_dest + c] = src_r[c];
      }
    }
  }

  // pad top and bottom
  const float *const src_t = img_dest + pad * iwidth * ch_dest;
  const float *const src_b = img_dest + (iheight - pad - 1) * iwidth * ch_dest;
  
  for (int y = 0; y < pad; y++)
  {
    float *dest_t = img_dest + y * iwidth * ch_dest;
    float *dest_b = img_dest + (y + iheight - pad) * iwidth * ch_dest;
    
    memcpy(dest_t, src_t, iwidth * ch_dest * sizeof(float));
    memcpy(dest_b, src_b, iwidth * ch_dest * sizeof(float));
  }

}
*/

static void image_rgb_to_srgb_padded(float *img_src, const int width, const int height, float *img_dest, const int pad, 
                                      const float beta)
{
  const int iwidth = width + pad;
  const int iheight = height + pad;

  memset(img_dest, 0, iwidth * iheight * 3 * sizeof(float));

  for (int y = 0; y < height; y++)
  {
    float *s = img_src + y * width * 4;
    float *d = img_dest + (y + pad) * iwidth * 3 + pad * 3;
    
    for (int x = 0; x < width; x++)
    {
      for (int c = 0; c < 3; c++)
      {
        d[x*3+c] = RGBtosRGB(s[x*4+c]) * RNG_UP_RGB;
      }
    }
  }

}

static void image_lab_to_srgb_padded(float *img_src, const int width, const int height, float *img_dest, const int pad, 
                                      const float beta)
{
  const int iwidth = width + pad;
  const int iheight = height + pad;

  memset(img_dest, 0, iwidth * iheight * 3 * sizeof(float));

  for (int y = 0; y < height; y++)
  {
    float *s = img_src + y * width * 4;
    float *d = img_dest + (y + pad) * iwidth * 3 + pad * 3;
    
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

static void unpad_image(float *img_src, const int width, const int height, float *img_dest, const int pad)
{
  const int iwidth = width + pad;

  memset(img_dest, 0, width * height * sizeof(float));

  for (int y = 0; y < height; y++)
  {
    float *s = img_src + (y + pad) * iwidth + pad;
    float *d = img_dest + y*width;
    
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

static inline float normalize_xy(const int xy, const int xy_L, const int xy_H, const int xy_min, const int xy_max)
{
  const float mult = (xy_H-xy_L) / (xy_max - xy_min);
  return ((xy - xy_min) * mult) + xy_L;
}

/*
static inline float rgb_luminance(const float r, const float g, const float b)
{
  return (r * 0.22248840f + g * 0.71690369f + b * 0.06060791f);
}

static float avg_luminance_rgb(float *img_src, const int width, const int height)
{
  float avg = 0.f;
      
  const int stride = width * height * 4;
  for (int i = 0; i < stride; i+=4)
  {
    avg += RGBtosRGB( rgb_luminance(img_src[i+0], img_src[i+1], img_src[i+2]) ) * RNG_UP_RGB;
  }

  avg /= (float)(width * height);
  
  return avg;
}
*/
static void gradient_rgb(float *img_src, float *img_dest_x, float *img_dest_y, float *img_dest_n, const int width, const int height,
    const float oddness, const float image_scale, const float beta)
{
  const int ch = 3;
  
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
        
//        RGBx[c] *= powf(s0[x*ch + c], -image_scale);
//        RGBy[c] *= powf(s0[x*ch + c], -image_scale);
        
        RGBn[c] = RGBx[c]*RGBx[c] + RGBy[c]*RGBy[c];
      }
      
      int max_bn = 0;
      if (RGBn[0] > RGBn[1])
        max_bn = (RGBn[0] > RGBn[2]) ? 0: 2;
      else
        max_bn = (RGBn[1] > RGBn[2]) ? 1: 2;
      
      float n = 1e-5 + powf(RGBn[max_bn], oddness);
      
//      n /= expf(RGBn[max_bn]*image_scale * .25f);
//      n *= powf(10.f, RGBn[max_bn]*image_scale);
      
      dx[x] = RGBx[max_bn] / n;
      dy[x] = RGBy[max_bn] / n;
      dn[x] = s0[(x+1)*ch + max_bn];
    }
  }

}

static void estimate_laplacian(float *img_src_x, float *img_src_y, float *img_dest, const int width, const int height, 
    const float beta, const float image_scale)
{
  const int stride = width * height;

  memset(img_dest, 0, stride * sizeof(float));

  for (int y = 0; y < height; y++)
  {
    float *s0 = img_src_x + y * width;

    for (int x = width-1; x > 0; x--)
    {
      s0[x] -= s0[x-1];
    }
  }

  for (int y = height-1; y > 0; y--)
  {
    float *s0 = img_src_y + y * width;
    float *s1 = img_src_y + (y - 1) * width;

    for (int x = 0; x < width; x++)
    {
      s0[x] -= s1[x];
    }
  }

  for (int i = 0; i < stride; i++)
  {
    img_dest[i] = (img_src_x[i] + img_src_y[i]);
  }

}

 void _recontruct_laplacian(float *img_src, float *img_grn, float *img_dest, const int width, const int height, 
    const int x_from, const int x_to, const int y_from, const int y_to, 
    const float beta, const float image_scale)
{
  fbw_fftw3_t fftw3_fbw = {0};
  
  fft(&fftw3_fbw, img_src, width, height, image_scale);
  
  const int width_fft = fftw3_fbw.width_fft;
  const int height_fft = fftw3_fbw.height_fft;
  
//  const int width_fft = x_to-x_from;
//  const int height_fft = y_to-y_from;
  
  const float piw = (2.f * (float)M_PI / (float)width_fft);
  const float pih = (2.f * (float)M_PI / (float)height_fft);
  
/*  printf("x_from=%i, x_to=%i, y_from=%i, y_to=%i\n", x_from, x_to, y_from, y_to);

  printf("norm x %i %f\n", 0, normalize_xy(0, x_from, x_to, 0, fftw3_fbw.width_fft_complex));
  printf("norm x %i %f\n", fftw3_fbw.width_fft, normalize_xy(fftw3_fbw.width_fft, x_from, x_to, 0, fftw3_fbw.width_fft));
  
  printf("norm y %i %f\n", 0, normalize_xy(0, y_from, y_to, 0, fftw3_fbw.height_fft_complex));
  printf("norm y %i %f\n", fftw3_fbw.height_fft, normalize_xy(fftw3_fbw.height_fft, y_from, y_to, 0, fftw3_fbw.height_fft));
  */
  for (int y = 0; y < fftw3_fbw.height_fft_complex; y++)
  {
    fftwf_complex *d = fftw3_fbw.out_src + y * fftw3_fbw.width_fft_complex;

    const float yy = y;
//    const float yy = normalize_xy(y, y_from, y_to, 0, fftw3_fbw.height_fft);
//    const float yy = y * (y_to / fftw3_fbw.height_fft);
    
    const float cos_y = cos(pih*yy);

    for (int x = 0; x < fftw3_fbw.width_fft_complex; x++)
    {
      const float xx = x;
//      const float xx = normalize_xy(x, x_from, x_to, 0, fftw3_fbw.width_fft);
//      const float xx = x * (x_to / fftw3_fbw.width_fft);
      
      const float cos_x = cos(piw*xx);
      
      float cos_xy = 1.f;
      if (x > 0 || y > 0) cos_xy = (cos_x + cos_y - 2.f) * 2.f;
      
      d[x][0] /= cos_xy;
      d[x][1] /= cos_xy;
    }
  }

  fftw3_fbw.out_src[0][0] = 0.f;
  fftw3_fbw.out_src[0][1] = 0.f;

  ifft(&fftw3_fbw, img_dest);

}

static void recontruct_laplacian(float *img_src, float *img_dest, const int width, const int height, 
    const float beta, const float image_scale)
{
  fbw_fftw3_t fftw3_fbw = {0};
  
  fft(&fftw3_fbw, img_src, width, height, image_scale);
  
  const int width_fft = fftw3_fbw.width_fft;
  const int height_fft = fftw3_fbw.height_fft;
  
  
//  const float piw = (2.f * (float)M_PI / (float)width_fft);
//  const float pih = (2.f * (float)M_PI / (float)height_fft);
  
  for (int y = 0; y < fftw3_fbw.height_fft_complex; y++)
  {
    fftwf_complex *d = fftw3_fbw.out_src + y * fftw3_fbw.width_fft_complex;

//    const float yy = y;
    
//    const float cos_y = cos(pih*yy);

    for (int x = 0; x < fftw3_fbw.width_fft_complex; x++)
    {
//      const float xx = x;
      
//      const float cos_x = cos(piw*xx);
      
      float cos_xy = 1.f;
//      if (x > 0 || y > 0) cos_xy = (cos_x + cos_y - 2.f) * 2.f;
      if (x > 0 || y > 0) cos_xy = (2*cos(M_PI*x/(width_fft-1))-2) + (2*cos(M_PI*y/(height_fft-1)) - 2);
      
      d[x][0] /= cos_xy;
      d[x][1] /= cos_xy;
    }
  }

//  fftw3_fbw.out_src[0][0] = 0.f;
//  fftw3_fbw.out_src[0][1] = 0.f;

  ifft(&fftw3_fbw, img_dest);

}

static void adjust_brcogm(float *img_src, const int width, const int height, const int ch,
    const float brightness, const float contrast, const float gamma,
    const float beta, const float _image_scale)
{
  const int stride = width * height * ch;
  const float image_scale = 1.f - _image_scale;
  
//  const float gm = powf(10.f, -gamma);
  
  const float gm = powf(10.f, -sqrtf(image_scale) / 10.f); // 6 .06
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] = powf(img_src[i], gm);
  }
  
//  const float cn = expf(contrast / .25f);
  
  const float cn = expf(image_scale / .25f); // 12  .12
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] = (img_src[i] - .5f);
    img_src[i] *= cn;
    img_src[i] = (img_src[i] + .5f);
  }
  
//  const float br = brightness * 2.f;
  
  const float br = brightness;
  
  for (int i = 0; i < stride; i++)
  {
    img_src[i] += br;
  }
  
  
}

static void adjust_brcogm2(float *img_src, const int width, const int height, const int ch,
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

 void fbw_process(float *img_src, float *img_dest, const int width, const int height, const int ch, 
    const float oddness, const float brightness, const float contrast, const float gamma,
    const float beta, const float avg_l, const float image_scale,
    const int x_from, const int x_to, const int y_from, const int y_to)
{
  
  float *img_padded = NULL;
  float *img_grx = NULL;
  float *img_gry = NULL;
  float *img_grn = NULL;
  
  const int pad = 1;
  
  const int iwidth = width + pad*2;
  const int iheight = height + pad*2;
  
  img_grx = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_grx == NULL) goto cleanup;
  
  img_gry = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_gry == NULL) goto cleanup;
  
  img_grn = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_grn == NULL) goto cleanup;
  
  img_padded = dt_alloc_align(64, iwidth * iheight * 3 * sizeof(float));
  if (img_padded == NULL) goto cleanup;
  
  print_min_max(img_src, width, height, 4, "ivoid");
  
  if (0)
  image_rgb_to_srgb_padded(img_src, width, height, img_padded, pad, beta);
  else
  image_lab_to_srgb_padded(img_src, width, height, img_padded, pad, beta);
  
  print_min_max(img_padded, iwidth, iheight, 3, "img_padded");
  
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

if(0)  adjust_brcogm(img_padded, iwidth, iheight, 3, 0, 0, 0, beta, image_scale);
  
  gradient_rgb(img_padded, img_grx, img_gry, img_grn, iwidth, iheight, oddness, image_scale, beta);

  print_min_max(img_grx, iwidth, iheight, 1, "img_grx");
  print_min_max(img_gry, iwidth, iheight, 1, "img_gry");
  
//  normalize(img_grx, iwidth, iheight, -10.f, 10.f);
//  normalize(img_gry, iwidth, iheight, -10.f, 10.f);
  
  printf("image_scale=%f\n", image_scale);
  
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
  
if(1)  estimate_laplacian(img_grx, img_gry, img_l, iwidth, iheight, beta, image_scale);

  print_min_max(img_grx, iwidth, iheight, 1, "img_grx 2");
  print_min_max(img_gry, iwidth, iheight, 1, "img_gry 2");
  print_min_max(img_l, iwidth, iheight, 1, "img_l");
  
  if(0){
  float *boundary_image = img_grn;
  
  for (int y = 2; y < iheight-2; y++)
  {
    float *d = boundary_image + y * iwidth;
    
    for (int x = 2; x < iwidth-2; x++)
    {
      d[x] = 0.f;
    }
  }
  
  float *f_bp = img_grx;
  
  memset(f_bp, 0, iwidth * iheight * sizeof(float));

  for (int y = 1; y < iheight-1; y++)
  {
    float *bp = f_bp + y * iwidth;
    float *b_i = boundary_image + y * iwidth;
    float *b_i_up = boundary_image + (y - 1) * iwidth;
    float *b_i_dn = boundary_image + (y + 1) * iwidth;

    for (int x = 1; x < iwidth-1; x++)
    {
      bp[x] = (-4.f * b_i[x] + b_i[x+1] + b_i[x-1] + b_i_up[x] + b_i_dn[x]);// * image_scale; 
    }
  }

  float *ff_bp = img_grx;

  for (int y = 0; y < 1; y++)
  {
    float *bp_t = ff_bp + y * iwidth;
    float *l_t = img_l + y * iwidth;
    float *bp_b = ff_bp + (iheight - y - 1) * iwidth;
    float *l_b = img_l + (iheight - y - 1) * iwidth;

    for (int x = 0; x < iwidth; x++)
    {
      l_t[x] += bp_t[x];
      l_b[x] += bp_b[x];
    }
  }
  for (int y = 0; y < iheight; y++)
  {
    float *bp = ff_bp + y * iwidth;
    float *l = img_l + y * iwidth;
   
    for (int x = 0; x < 1; x++)
    {
      l[x] += bp[x];
      l[iwidth - x - 1] += bp[iwidth - x - 1];
    }
  }

  }
/*
  # Reconstruct image from laplacian.
  -fft.
  100%,100%,1,1,'cos(2*x*pi/w)+cos(2*y*pi/h)' -*. 2 --. 4
  -=. 1 -/[-3,-2] . -rm.
  -=.. 0 -=. 0
  -ifft[-2,-1] -rm.
*/
if(1)  _recontruct_laplacian(img_l, img_grn, img_l, iwidth, iheight, x_from, x_to, y_from, y_to, beta, image_scale);

/*  const int stride = iwidth * iheight;
  for (int i = 0; i < stride; i++)
    img_l[i] *= image_scale;
*/
  print_min_max(img_l, iwidth, iheight, 1, "recontruct_laplacian");

/*
  -shrink_xy. 1 -+. $ia -n. 0,255
*/
  float *img_unpad_gr = img_grx;
  
//  img_l = img_grx;
//  float *img_unpad_gr = img_gry;
  
  unpad_image(img_l, width, height, img_unpad_gr, pad);
  
  print_min_max(img_unpad_gr, width, height, 1, "unpad_image");

/*  printf("avg_l=%f\n", avg_l);
  printf("beta=%f\n", beta);

  const int stride = width * height;
  for (int i = 0; i < stride; i++)
    img_unpad_gr[i] += avg_l;
*/
  
  normalize(img_unpad_gr, width, height, 0.f, RNG_UP_RGB);

  print_min_max(img_unpad_gr, width, height, 1, "normalize");

/*
  # Merge result with original color image.
  -j[0] [1],0,0,0,0,{$1%} -rm.
  -adjust_colors ${3-5}
  -endl -a c -endl -done
*/
//  adjust_brcogm(img_unpad_gr, width, height, 1, image_scale/10.f, image_scale/10.f, 0.f);
if(0)  adjust_brcogm(img_unpad_gr, width, height, 1, 0, 0, 0, beta, image_scale);

if (0)
image_srgb_to_rgb_bw(img_unpad_gr, width, height, img_dest);
else
image_srgb_to_lab_bw(img_unpad_gr, width, height, img_dest);

if(0)  adjust_brcogm2(img_dest, width, height, ch, brightness, contrast, gamma);
  
  print_min_max(img_dest, width, height, 4, "ovoid");
  
cleanup:

  if (img_grx) dt_free_align(img_grx);
  if (img_gry) dt_free_align(img_gry);
  if (img_grn) dt_free_align(img_grn);
  if (img_padded) dt_free_align(img_padded);
  

}

 void _fbw_process(float *img_src, float *img_dest, const int width, const int height, const int ch, 
    const float oddness, const float brightness, const float contrast, const float gamma,
    const float beta, const float avg_l, const float image_scale,
    const int x_from, const int x_to, const int y_from, const int y_to)
{
  
  float *img_padded = NULL;
  float *img_grx = NULL;
  float *img_gry = NULL;
  float *img_grxx = NULL;
  float *img_gryy = NULL;
  float *img_grn = NULL;
  float *f = NULL;
  float *f_bp = NULL;
  
  const int pad = 1;
  
  const int iwidth = width + pad*2;
  const int iheight = height + pad*2;
  
  img_grx = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_grx == NULL) goto cleanup;
  
  img_gry = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_gry == NULL) goto cleanup;
  
  img_grxx = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_grxx == NULL) goto cleanup;
  
  img_gryy = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_gryy == NULL) goto cleanup;
  
  img_grn = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (img_grn == NULL) goto cleanup;
  
  f = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (f == NULL) goto cleanup;
  
  f_bp = dt_alloc_align(64, iwidth * iheight * sizeof(float));
  if (f_bp == NULL) goto cleanup;
  
  img_padded = dt_alloc_align(64, iwidth * iheight * 3 * sizeof(float));
  if (img_padded == NULL) goto cleanup;
  
  float *img_dest_x = img_grx;
  float *img_dest_y = img_gry;
  float *img_dest_n = img_grn;
  
  memset(img_dest_x, 0, width * height * sizeof(float));
  memset(img_dest_y, 0, width * height * sizeof(float));

  for (int y = 1; y < height-1; y++)
  {
    float *s0 = img_src + y * width * ch;
    float *s1 = img_src + (y + 1) * width * ch;
    float *dx = img_dest_x + y * width;
    float *dy = img_dest_y + y * width;
    float *dn = img_dest_n + y * width;
    
    for (int x = 1; x < width-1; x++)
    {
      float RGBx[3] = {0};
      float RGBy[3] = {0};
      float RGBn[3] = {0};

      for (int c = 0; c < 3; c++)
      {
        RGBx[c] = s0[(x+1)*ch + c] - s0[x*ch + c];
        RGBy[c] = s1[x*ch + c] - s0[x*ch + c];
        
        RGBn[c] = RGBx[c]*RGBx[c] + RGBy[c]*RGBy[c];
      }
      
      int max_bn = 0;
      if (RGBn[0] > RGBn[1])
        max_bn = (RGBn[0] > RGBn[2]) ? 0: 2;
      else
        max_bn = (RGBn[1] > RGBn[2]) ? 1: 2;
      
      float n = 1e-5 + powf(RGBn[max_bn], oddness);
      n = 1.f;
      
      dx[x] = RGBx[max_bn] / n;
      dy[x] = RGBy[max_bn] / n;
      dn[x] = s0[(x+1)*ch + max_bn];
    }
  }

  memset(img_grxx, 0, width * height * sizeof(float));
  memset(img_gryy, 0, width * height * sizeof(float));

  for (int y = 1; y < height-1; y++)
  {
    float *gx = img_grx + y * width;
    float *gy = img_gry + y * width;
    float *gy1 = img_gry + (y + 1) * width;
    float *gxx = img_grxx + y * width;
    float *gyy = img_gryy + (y + 1) * width;
    
    for (int x = 1; x < width-1; x++)
    {
      gyy[x] = gy1[x] - gy[x];
      gxx[x+1] = gx[x+1] - gx[x];
    }
  }
  
  memset(f, 0, width * height * sizeof(float));
  for (int i = 0; i < width*height; i++)
    f[i] = img_grxx[i] + img_gryy[i];
  
  
  float *boundary_image = img_grn;
  
  for (int y = 2; y < height-1; y++)
  {
    float *d = boundary_image + y * width;
    
    for (int x = 2; x < width-1; x++)
    {
      d[x] = 0.f;
    }
  }
  
  memset(f_bp, 0, width * height * sizeof(float));

  for (int y = 2; y < height-1; y++)
  {
    float *bp = f_bp + y * width;
    float *b_i = boundary_image + y * width;
    float *b_i_up = boundary_image + (y - 1) * width;
    float *b_i_dn = boundary_image + (y + 1) * width;

    for (int x = 2; x < width-1; x++)
    {
      bp[x] = -4.f * b_i[x] + b_i[x+1] + b_i[x-1] + b_i_up[x] + b_i_dn[x]; 
    }
  }

  for (int i = 0; i < width*height; i++)
    f[i] -= f_bp[i];

  recontruct_laplacian(f, img_grn, width, height, beta, image_scale);
  
  normalize(img_grn, width, height, 0.f, RNG_UP_RGB);
  
  for (int i = 0; i < width*height; i++)
    img_dest[i*ch + 0] = img_dest[i*ch + 1] = img_dest[i*ch + 2] = img_grn[i];

cleanup:

  if (img_grx) dt_free_align(img_grx);
  if (img_gry) dt_free_align(img_gry);
  if (img_grxx) dt_free_align(img_grxx);
  if (img_gryy) dt_free_align(img_gryy);
  if (img_grn) dt_free_align(img_grn);
  if (img_padded) dt_free_align(img_padded);
  if (f) dt_free_align(f);
  if (f_bp) dt_free_align(f_bp);
  

}

void process_internal(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid, void *const ovoid,
    const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out, const int use_sse)
{
  const dt_iop_fbw_data_t *const p = (const dt_iop_fbw_data_t *const)piece->data;
/*  dt_iop_fbw_gui_data_t *g = (dt_iop_fbw_gui_data_t *)self->gui_data;

  float avg_l;
  float tmp_avg_l = NAN;

  if(self->dev->gui_attached && g && piece->pipe->type == DT_DEV_PIXELPIPE_FULL)
  {
    dt_pthread_mutex_lock(&g->lock);
    const uint64_t hash = g->hash;
    dt_pthread_mutex_unlock(&g->lock);

    // note that the case 'hash == 0' on first invocation in a session implies that g->avg_l
    // is NAN which initiates special handling below to avoid inconsistent results. in all
    // other cases we make sure that the preview pipe has left us with proper readings for
    // avg_l. if data are not yet there we need to wait (with timeout).
    if(hash != 0 && !dt_dev_sync_pixelpipe_hash(self->dev, piece->pipe, 0, self->priority, &g->lock, &g->hash))
      dt_control_log(_("freaky b&w: inconsistent output"));

    dt_pthread_mutex_lock(&g->lock);
    tmp_avg_l = g->avg_l;
    dt_pthread_mutex_unlock(&g->lock);
  }

  // in all other cases we calculate avg_l here
  if(isnan(tmp_avg_l))
  {
    avg_l = avg_luminance_rgb((float *)ivoid, roi_in->width, roi_in->height);
  }
  else
  {
    avg_l = tmp_avg_l;
  }

  // PREVIEW pixelpipe stores avg_l
  if(self->dev->gui_attached && g && piece->pipe->type == DT_DEV_PIXELPIPE_PREVIEW)
  {
    uint64_t hash = dt_dev_hash_plus(self->dev, piece->pipe, 0, self->priority);
    dt_pthread_mutex_lock(&g->lock);
    g->avg_l = avg_l;
    g->hash = hash;
    dt_pthread_mutex_unlock(&g->lock);
  }
*/
  const float image_scale = roi_in->scale / piece->iscale;
  
  const int x_from = roi_in->x / image_scale;
  const int x_to = x_from + roi_in->width / image_scale;
  const int y_from = roi_in->y / image_scale;
  const int y_to = y_from + roi_in->height / image_scale;
  
//  printf("x_from=%i, x_to=%i, y_from=%i, y_to=%i\n", x_from, x_to, y_from, y_to);
//  printf("roi_in->width=%i, x_to-x_from=%i, roi_in->height=%i, y_to-y_from=%i\n", roi_in->width, x_to-x_from, roi_in->height, y_to-y_from);
  
  fbw_process((float *)ivoid, (float *)ovoid, roi_in->width, roi_in->height, piece->colors, 
                  p->oddness/100.f, p->brightness/100.f, p->contrast/100.f, p->gamma/100.f,
                  p->beta, /* avg_l */ 0.f, image_scale,
                  x_from, x_to, y_from, y_to);

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
