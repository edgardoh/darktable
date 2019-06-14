/*
    This file is part of darktable,
    copyright (c) 2016 Roman Lebedev.

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

#include "common/color_picker.h"
#include "common/darktable.h"
#include "common/colorspaces_inline_conversions.h"
#include "common/iop_order.h"
#include "develop/format.h"
#include "develop/imageop.h"
#include "develop/imageop_math.h"
/*
#define QS_SWAP(a, b) \
{ \
  float tmp = a; \
  a = b; \
  b = tmp; \
} \

static int _qs_partition(float *arr, const int left, const int right)
{
  const int mid = left + (right - left) / 2;
  const int pivot = arr[mid];
  // move the mid point value to the front.
  QS_SWAP(arr[mid], arr[left]);
  int i = left + 1;
  int j = right;
  while (i <= j)
  {
    while(i <= j && arr[i] <= pivot) i++;
    while(i <= j && arr[j] > pivot) j--;

    if (i < j)
    {
      QS_SWAP(arr[i], arr[j]);
    }
  }
  QS_SWAP(arr[i - 1], arr[left]);
  return i - 1;
}

static void _quicksort(float *arr, const int left, const int right, const int sz)
{
  if (left >= right) return;

  int part = _qs_partition(arr, left, right);

  _quicksort(arr, left, part - 1, sz);
  _quicksort(arr, part + 1, right, sz);
}

void quicksort(float *arr, const int sz)
{
  _quicksort(arr, 0, sz - 1, sz);
}
*/

void swap(float array[], int left, int right)
{
	float temp;
	temp = array[left];
	array[left] = array[right];
	array[right] = temp;
}

int partition(float array[], int left, int right, int pivot_index)
{
	float pivot_value = array[pivot_index];
	int store_index = left;
	int i;

//	printf("%d %d %d\n", left, right, pivot_index);
	swap(array, pivot_index, right);
	for (i = left; i < right; i++)
		if (array[i] <= pivot_value) {
			swap(array, i, store_index);
			++store_index;
		}
	swap(array, store_index, right);
	return store_index;
}

void quicksort(float array[], int left, int right)
{
	int pivot_index = left;
	int pivot_new_index;

	loop:
	//printf("sorting %d to %d\n", left, right);
	if (right > left) {
		pivot_new_index = partition(array, left, right, pivot_index);
	//	printf("L\n");
		quicksort(array, left, pivot_new_index - 1);
//		printf("R\n");
		pivot_index = left = pivot_new_index + 1;
		goto loop;
	}
}

void color_picker_quicksort(float *arr, const int sz)
{
  quicksort(arr, 0, sz - 1);
}

static void color_picker_helper_1ch_gap(const float *const pixel, const dt_iop_roi_t *roi, const int *const box,
                                        float *const picked_color, float *const picked_color_min, float *const picked_color_max,
                                        const dt_iop_colorspace_type_t image_cst, const dt_iop_colorspace_type_t cst_to,
                                        const dt_iop_order_iccprofile_info_t *const profile_info, const dt_color_picker_channels_t channel)
{
printf("color_picker_helper_1ch_gap\n");
  const int width = roi->width;
  const size_t size_w = (box[2] - box[0]);
  const size_t size_h = (box[3] - box[1]);
  const size_t size = size_w * size_h;
  float *img = dt_alloc_align(64, size * sizeof(float));
  const float w = 1.0f / (float)size;

//  float *img_tmp = img;
  for(size_t j = box[1]; j < box[3]; j++)
  {
    for(size_t i = box[0]; i < box[2]; i++)
    {
      const size_t k = 4 * (width * j + i);
      float Lab[3] = { pixel[k], pixel[k + 1], pixel[k + 2] };

      if(channel == DT_COLOR_PICKER_CHANNEL_GREY)
      {
        Lab[0] = (profile_info) ? dt_ioppr_get_rgb_matrix_luminance(Lab, profile_info) : dt_camera_rgb_luminance(Lab);
      }
      else
      {
        if(cst_to == iop_cs_LCh)
          dt_Lab_2_LCH(pixel + k, Lab);
        else if(cst_to == iop_cs_HSL)
          dt_RGB_2_HSL(pixel + k, Lab);

        Lab[0] = Lab[channel];
      }

      img[(j - box[1]) * size_w + (i - box[0])] = Lab[0];

//      picked_color[0] += Lab[0];

//      *img_tmp = Lab[0];
//      img_tmp++;
    }
  }

//  picked_color[0] *= w;

  printf("color_picker_helper_1ch_gap 1\n");

  color_picker_quicksort(img, size);

  for(size_t i = 0; i < size - 1; i++)
  {
    if(img[i] > img[i + 1]) printf("[color_picker_helper_1ch_gap] ERROR in quicksort %f > %f\n", img[i], img[i + 1]);
  }

  printf("color_picker_helper_1ch_gap 2, %f : %f\n", img[0], img[size - 1]);

  float dist = img[size - 1] - img[0];
  picked_color[0] = 0.f;
  picked_color_min[0] = img[0];
  picked_color_max[0] = img[size - 1];

  for(size_t i = 0; i < size - 1; i++)
  {
    picked_color[0] += img[i];

    if(dist < img[i + 1] - img[i])
    {
      dist = img[i + 1] - img[i];
      picked_color_min[0] = img[i];
      picked_color_max[0] = img[i + 1];
    }
  }

  picked_color[0] = (picked_color[0] + img[size - 1]) * w;
  picked_color[1] = picked_color[2] = picked_color[0];
  picked_color_min[1] = picked_color_min[2] = picked_color_min[0];
  picked_color_max[1] = picked_color_max[2] = picked_color_max[0];

  printf("color_picker_helper_1ch_gap 2, %f : %f\n", picked_color_min[0], picked_color_max[0]);

  dt_free_align(img);
}

static void color_picker_helper_1ch(const float *const pixel, const dt_iop_roi_t *roi, const int *const box,
                                        float *const picked_color, float *const picked_color_min, float *const picked_color_max,
                                        const dt_iop_colorspace_type_t image_cst, const dt_iop_colorspace_type_t cst_to,
                                        const dt_iop_order_iccprofile_info_t *const profile_info, const dt_color_picker_channels_t channel)
{
  const int width = roi->width;

  const size_t size = ((box[3] - box[1]) * (box[2] - box[0]));

  const float w = 1.0f / (float)size;

  picked_color[0] = 0.f;
  picked_color_min[0] = INFINITY;
  picked_color_max[0] = -INFINITY;

  for(size_t j = box[1]; j < box[3]; j++)
  {
    for(size_t i = box[0]; i < box[2]; i++)
    {
      const size_t k = 4 * (width * j + i);
      float Lab[3] = { pixel[k], pixel[k + 1], pixel[k + 2] };

      if(channel == DT_COLOR_PICKER_CHANNEL_GREY)
      {
        Lab[0] = (profile_info) ? dt_ioppr_get_rgb_matrix_luminance(Lab, profile_info) : dt_camera_rgb_luminance(Lab);
      }
      else
      {
        if(cst_to == iop_cs_LCh)
          dt_Lab_2_LCH(pixel + k, Lab);
        else if(cst_to == iop_cs_HSL)
          dt_RGB_2_HSL(pixel + k, Lab);

        Lab[0] = Lab[channel];
      }

      picked_color[0] += Lab[0];
      picked_color_min[0] = fminf(picked_color_min[0], Lab[0]);
      picked_color_max[0] = fmaxf(picked_color_max[0], Lab[0]);
    }
  }

  picked_color[0] *= w;
  picked_color[1] = picked_color[2] = picked_color[0];
  picked_color_min[1] = picked_color_min[2] = picked_color_min[0];
  picked_color_max[1] = picked_color_max[2] = picked_color_max[0];
}

static void color_picker_helper_4ch_seq(const dt_iop_buffer_dsc_t *dsc, const float *const pixel,
                                        const dt_iop_roi_t *roi, const int *const box, float *const picked_color,
                                        float *const picked_color_min, float *const picked_color_max,
                                        const dt_iop_colorspace_type_t cst_to)
{
  const int width = roi->width;

  const size_t size = ((box[3] - box[1]) * (box[2] - box[0]));

  const float w = 1.0f / (float)size;

  // code path for small region, especially for color picker point mode
  for(size_t j = box[1]; j < box[3]; j++)
  {
    for(size_t i = box[0]; i < box[2]; i++)
    {
      const size_t k = 4 * (width * j + i);
      float Lab[3] = { pixel[k], pixel[k + 1], pixel[k + 2] };
      if(cst_to == iop_cs_LCh)
        dt_Lab_2_LCH(pixel + k, Lab);
      if(cst_to == iop_cs_HSL)
        dt_RGB_2_HSL(pixel + k, Lab);
      picked_color[0] += w * Lab[0];
      picked_color[1] += w * Lab[1];
      picked_color[2] += w * Lab[2];
      picked_color_min[0] = fminf(picked_color_min[0], Lab[0]);
      picked_color_min[1] = fminf(picked_color_min[1], Lab[1]);
      picked_color_min[2] = fminf(picked_color_min[2], Lab[2]);
      picked_color_max[0] = fmaxf(picked_color_max[0], Lab[0]);
      picked_color_max[1] = fmaxf(picked_color_max[1], Lab[1]);
      picked_color_max[2] = fmaxf(picked_color_max[2], Lab[2]);
    }
  }
}

static void color_picker_helper_4ch_parallel(const dt_iop_buffer_dsc_t *dsc, const float *const pixel,
                                             const dt_iop_roi_t *roi, const int *const box,
                                             float *const picked_color, float *const picked_color_min,
                                             float *const picked_color_max, const dt_iop_colorspace_type_t cst_to)
{
  const int width = roi->width;

  const size_t size = ((box[3] - box[1]) * (box[2] - box[0]));

  const float w = 1.0f / (float)size;

  const int numthreads = dt_get_num_threads();

  float *const mean = malloc((size_t)3 * numthreads * sizeof(float));
  float *const mmin = malloc((size_t)3 * numthreads * sizeof(float));
  float *const mmax = malloc((size_t)3 * numthreads * sizeof(float));

  for(int n = 0; n < 3 * numthreads; n++)
  {
    mean[n] = 0.0f;
    mmin[n] = INFINITY;
    mmax[n] = -INFINITY;
  }

#ifdef _OPENMP
#pragma omp parallel default(none) \
  dt_omp_firstprivate(w, cst_to, pixel, width, box, mean, mmin, mmax)
#endif
  {
    const int tnum = dt_get_thread_num();

    float *const tmean = mean + 3 * tnum;
    float *const tmmin = mmin + 3 * tnum;
    float *const tmmax = mmax + 3 * tnum;

#ifdef _OPENMP
#pragma omp for schedule(static) collapse(2)
#endif
    for(size_t j = box[1]; j < box[3]; j++)
    {
      for(size_t i = box[0]; i < box[2]; i++)
      {
        const size_t k = 4 * (width * j + i);
        float Lab[3] = { pixel[k], pixel[k + 1], pixel[k + 2] };
        if(cst_to == iop_cs_LCh)
          dt_Lab_2_LCH(pixel + k, Lab);
        if(cst_to == iop_cs_HSL)
          dt_RGB_2_HSL(pixel + k, Lab);
        tmean[0] += w * Lab[0];
        tmean[1] += w * Lab[1];
        tmean[2] += w * Lab[2];
        tmmin[0] = fminf(tmmin[0], Lab[0]);
        tmmin[1] = fminf(tmmin[1], Lab[1]);
        tmmin[2] = fminf(tmmin[2], Lab[2]);
        tmmax[0] = fmaxf(tmmax[0], Lab[0]);
        tmmax[1] = fmaxf(tmmax[1], Lab[1]);
        tmmax[2] = fmaxf(tmmax[2], Lab[2]);
      }
    }
  }

  for(int n = 0; n < numthreads; n++)
  {
    for(int k = 0; k < 3; k++)
    {
      picked_color[k] += mean[3 * n + k];
      picked_color_min[k] = fminf(picked_color_min[k], mmin[3 * n + k]);
      picked_color_max[k] = fmaxf(picked_color_max[k], mmax[3 * n + k]);
    }
  }

  free(mmax);
  free(mmin);
  free(mean);
}

static void color_picker_helper_4ch(const dt_iop_buffer_dsc_t *dsc, const float *const pixel,
                                    const dt_iop_roi_t *roi, const int *const box, float *const picked_color,
                                    float *const picked_color_min, float *const picked_color_max, const dt_iop_colorspace_type_t cst_to)
{
  const size_t size = ((box[3] - box[1]) * (box[2] - box[0]));

  if(size > 100) // avoid inefficient multi-threading in case of small region size (arbitrary limit)
    return color_picker_helper_4ch_parallel(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max, cst_to);
  else
    return color_picker_helper_4ch_seq(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max, cst_to);
}

static void color_picker_helper_bayer_seq(const dt_iop_buffer_dsc_t *const dsc, const float *const pixel,
                                          const dt_iop_roi_t *const roi, const int *const box,
                                          float *const picked_color, float *const picked_color_min,
                                          float *const picked_color_max)
{
  const int width = roi->width;
  const uint32_t filters = dsc->filters;

  uint32_t weights[4] = { 0u, 0u, 0u, 0u };

  // code path for small region, especially for color picker point mode
  for(size_t j = box[1]; j < box[3]; j++)
  {
    for(size_t i = box[0]; i < box[2]; i++)
    {
      const int c = FC(j + roi->y, i + roi->x, filters);
      const size_t k = width * j + i;

      const float v = pixel[k];

      picked_color[c] += v;
      picked_color_min[c] = fminf(picked_color_min[c], v);
      picked_color_max[c] = fmaxf(picked_color_max[c], v);
      weights[c]++;
    }
  }

  // and finally normalize data. For bayer, there is twice as much green.
  for(int c = 0; c < 4; c++)
  {
    picked_color[c] = weights[c] ? (picked_color[c] / (float)weights[c]) : 0.0f;
  }
}

static void color_picker_helper_bayer_parallel(const dt_iop_buffer_dsc_t *const dsc, const float *const pixel,
                                               const dt_iop_roi_t *const roi, const int *const box,
                                               float *const picked_color, float *const picked_color_min,
                                               float *const picked_color_max)
{
  const int width = roi->width;
  const uint32_t filters = dsc->filters;

  uint32_t weights[4] = { 0u, 0u, 0u, 0u };

  const int numthreads = dt_get_num_threads();

  float *const msum = malloc((size_t)4 * numthreads * sizeof(float));
  float *const mmin = malloc((size_t)4 * numthreads * sizeof(float));
  float *const mmax = malloc((size_t)4 * numthreads * sizeof(float));
  uint32_t *const cnt = malloc((size_t)4 * numthreads * sizeof(uint32_t));

  for(int n = 0; n < 4 * numthreads; n++)
  {
    msum[n] = 0.0f;
    mmin[n] = INFINITY;
    mmax[n] = -INFINITY;
    cnt[n] = 0u;
  }

#ifdef _OPENMP
#pragma omp parallel default(none) \
  dt_omp_firstprivate(pixel, width, roi, filters, box, msum, mmin, mmax, cnt)
#endif
  {
    const int tnum = dt_get_thread_num();

    float *const tsum = msum + 4 * tnum;
    float *const tmmin = mmin + 4 * tnum;
    float *const tmmax = mmax + 4 * tnum;
    uint32_t *const tcnt = cnt + 4 * tnum;

#ifdef _OPENMP
#pragma omp for schedule(static) collapse(2)
#endif
    for(size_t j = box[1]; j < box[3]; j++)
    {
      for(size_t i = box[0]; i < box[2]; i++)
      {
        const int c = FC(j + roi->y, i + roi->x, filters);
        const size_t k = width * j + i;

        const float v = pixel[k];

        tsum[c] += v;
        tmmin[c] = fminf(tmmin[c], v);
        tmmax[c] = fmaxf(tmmax[c], v);
        tcnt[c]++;
      }
    }
  }

  for(int n = 0; n < numthreads; n++)
  {
    for(int c = 0; c < 4; c++)
    {
      picked_color[c] += msum[4 * n + c];
      picked_color_min[c] = fminf(picked_color_min[c], mmin[4 * n + c]);
      picked_color_max[c] = fmaxf(picked_color_max[c], mmax[4 * n + c]);
      weights[c] += cnt[4 * n + c];
    }
  }

  free(cnt);
  free(mmax);
  free(mmin);
  free(msum);

  // and finally normalize data. For bayer, there is twice as much green.
  for(int c = 0; c < 4; c++)
  {
    picked_color[c] = weights[c] ? (picked_color[c] / (float)weights[c]) : 0.0f;
  }
}

static void color_picker_helper_bayer(const dt_iop_buffer_dsc_t *dsc, const float *const pixel,
                                      const dt_iop_roi_t *roi, const int *const box, float *const picked_color,
                                      float *const picked_color_min, float *const picked_color_max)
{
  const size_t size = ((box[3] - box[1]) * (box[2] - box[0]));

  if(size > 100) // avoid inefficient multi-threading in case of small region size (arbitrary limit)
    return color_picker_helper_bayer_parallel(dsc, pixel, roi, box, picked_color, picked_color_min,
                                              picked_color_max);
  else
    return color_picker_helper_bayer_seq(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max);
}

static void color_picker_helper_xtrans_seq(const dt_iop_buffer_dsc_t *const dsc, const float *const pixel,
                                           const dt_iop_roi_t *const roi, const int *const box,
                                           float *const picked_color, float *const picked_color_min,
                                           float *const picked_color_max)
{
  const int width = roi->width;
  const uint8_t(*const xtrans)[6] = (const uint8_t(*const)[6])dsc->xtrans;

  uint32_t weights[3] = { 0u, 0u, 0u };

  // code path for small region, especially for color picker point mode
  for(size_t j = box[1]; j < box[3]; j++)
  {
    for(size_t i = box[0]; i < box[2]; i++)
    {
      const int c = FCxtrans(j, i, roi, xtrans);
      const size_t k = width * j + i;

      const float v = pixel[k];

      picked_color[c] += v;
      picked_color_min[c] = fminf(picked_color_min[c], v);
      picked_color_max[c] = fmaxf(picked_color_max[c], v);
      weights[c]++;
    }
  }

  // and finally normalize data.
  // X-Trans RGB weighting averages to 2:5:2 for each 3x3 cell
  for(int c = 0; c < 3; c++)
  {
    picked_color[c] /= (float)weights[c];
  }
}

static void color_picker_helper_xtrans_parallel(const dt_iop_buffer_dsc_t *const dsc, const float *const pixel,
                                                const dt_iop_roi_t *const roi, const int *const box,
                                                float *const picked_color, float *const picked_color_min,
                                                float *const picked_color_max)
{
  const int width = roi->width;
  const uint8_t(*const xtrans)[6] = (const uint8_t(*const)[6])dsc->xtrans;

  uint32_t weights[3] = { 0u, 0u, 0u };

  const int numthreads = dt_get_num_threads();

  float *const msum = malloc((size_t)3 * numthreads * sizeof(float));
  float *const mmin = malloc((size_t)3 * numthreads * sizeof(float));
  float *const mmax = malloc((size_t)3 * numthreads * sizeof(float));
  uint32_t *const cnt = malloc((size_t)3 * numthreads * sizeof(uint32_t));

  for(int n = 0; n < 3 * numthreads; n++)
  {
    msum[n] = 0.0f;
    mmin[n] = INFINITY;
    mmax[n] = -INFINITY;
    cnt[n] = 0u;
  }

#ifdef _OPENMP
#pragma omp parallel default(none) \
  dt_omp_firstprivate(pixel, width, roi, xtrans, box, cnt, msum, mmin, mmax)
#endif
  {
    const int tnum = dt_get_thread_num();

    float *const tsum = msum + 3 * tnum;
    float *const tmmin = mmin + 3 * tnum;
    float *const tmmax = mmax + 3 * tnum;
    uint32_t *const tcnt = cnt + 3 * tnum;

#ifdef _OPENMP
#pragma omp for schedule(static) collapse(2)
#endif
    for(size_t j = box[1]; j < box[3]; j++)
    {
      for(size_t i = box[0]; i < box[2]; i++)
      {
        const int c = FCxtrans(j, i, roi, xtrans);
        const size_t k = width * j + i;

        const float v = pixel[k];

        tsum[c] += v;
        tmmin[c] = fminf(tmmin[c], v);
        tmmax[c] = fmaxf(tmmax[c], v);
        tcnt[c]++;
      }
    }
  }

  for(int n = 0; n < numthreads; n++)
  {
    for(int c = 0; c < 3; c++)
    {
      picked_color[c] += msum[3 * n + c];
      picked_color_min[c] = fminf(picked_color_min[c], mmin[3 * n + c]);
      picked_color_max[c] = fmaxf(picked_color_max[c], mmax[3 * n + c]);
      weights[c] += cnt[3 * n + c];
    }
  }

  free(cnt);
  free(mmax);
  free(mmin);
  free(msum);

  // and finally normalize data.
  // X-Trans RGB weighting averages to 2:5:2 for each 3x3 cell
  for(int c = 0; c < 3; c++)
  {
    picked_color[c] /= (float)weights[c];
  }
}

static void color_picker_helper_xtrans(const dt_iop_buffer_dsc_t *dsc, const float *const pixel,
                                       const dt_iop_roi_t *roi, const int *const box, float *const picked_color,
                                       float *const picked_color_min, float *const picked_color_max)
{
  const size_t size = ((box[3] - box[1]) * (box[2] - box[0]));

  if(size > 100) // avoid inefficient multi-threading in case of small region size (arbitrary limit)
    return color_picker_helper_xtrans_parallel(dsc, pixel, roi, box, picked_color, picked_color_min,
                                               picked_color_max);
  else
    return color_picker_helper_xtrans_seq(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max);
}

void dt_color_picker_helper(const dt_iop_buffer_dsc_t *dsc, const float *const pixel, const dt_iop_roi_t *roi,
                            const int *const box, float *const picked_color, float *const picked_color_min,
                            float *const picked_color_max, const dt_iop_colorspace_type_t image_cst, const dt_iop_colorspace_type_t picker_cst)
{
  printf("[dt_color_picker_helper]\n");
  if((dsc->channels == 4u) && ((image_cst == picker_cst) || (picker_cst == iop_cs_NONE)))
    color_picker_helper_4ch(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max, picker_cst);
  else if(dsc->channels == 4u && image_cst == iop_cs_Lab && picker_cst == iop_cs_LCh)
    color_picker_helper_4ch(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max, picker_cst);
  else if(dsc->channels == 4u && image_cst == iop_cs_rgb && picker_cst == iop_cs_HSL)
    color_picker_helper_4ch(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max, picker_cst);
  else if(dsc->channels == 1u && dsc->filters && dsc->filters != 9u)
    color_picker_helper_bayer(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max);
  else if(dsc->channels == 1u && dsc->filters && dsc->filters == 9u)
    color_picker_helper_xtrans(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max);
  else
    dt_unreachable_codepath();
}

void dt_color_picker_helper_1c_gap(const dt_iop_buffer_dsc_t *dsc, const float *const pixel, const dt_iop_roi_t *roi,
                            const int *const box, float *const picked_color, float *const picked_color_min,
                            float *const picked_color_max, const dt_iop_colorspace_type_t image_cst, const dt_iop_colorspace_type_t picker_cst,
                            const dt_iop_order_iccprofile_info_t *const profile_info,
                            const dt_color_picker_channels_t channel)
{
  printf("[dt_color_picker_helper_1c_gap] channel=%i\n", channel);
  if((dsc->channels == 4u) && ((image_cst == picker_cst) || (picker_cst == iop_cs_NONE)))
    color_picker_helper_1ch_gap(pixel, roi, box, picked_color, picked_color_min, picked_color_max, image_cst, picker_cst, profile_info, channel);
  else if(dsc->channels == 4u && image_cst == iop_cs_Lab && picker_cst == iop_cs_LCh)
    color_picker_helper_1ch_gap(pixel, roi, box, picked_color, picked_color_min, picked_color_max, image_cst, picker_cst, profile_info, channel);
  else if(dsc->channels == 4u && image_cst == iop_cs_rgb && picker_cst == iop_cs_HSL)
    color_picker_helper_1ch_gap(pixel, roi, box, picked_color, picked_color_min, picked_color_max, image_cst, picker_cst, profile_info, channel);
  else if(dsc->channels == 1u && dsc->filters && dsc->filters != 9u)
    color_picker_helper_bayer(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max);
  else if(dsc->channels == 1u && dsc->filters && dsc->filters == 9u)
    color_picker_helper_xtrans(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max);
  else
    dt_unreachable_codepath();
}

void dt_color_picker_helper_1c(const dt_iop_buffer_dsc_t *dsc, const float *const pixel, const dt_iop_roi_t *roi,
                            const int *const box, float *const picked_color, float *const picked_color_min,
                            float *const picked_color_max, const dt_iop_colorspace_type_t image_cst, const dt_iop_colorspace_type_t picker_cst,
                            const dt_iop_order_iccprofile_info_t *const profile_info,
                            const dt_color_picker_channels_t channel)
{
  printf("[dt_color_picker_helper_1c] channel=%i\n", channel);
  if((dsc->channels == 4u) && ((image_cst == picker_cst) || (picker_cst == iop_cs_NONE)))
    color_picker_helper_1ch(pixel, roi, box, picked_color, picked_color_min, picked_color_max, image_cst, picker_cst, profile_info, channel);
  else if(dsc->channels == 4u && image_cst == iop_cs_Lab && picker_cst == iop_cs_LCh)
    color_picker_helper_1ch(pixel, roi, box, picked_color, picked_color_min, picked_color_max, image_cst, picker_cst, profile_info, channel);
  else if(dsc->channels == 4u && image_cst == iop_cs_rgb && picker_cst == iop_cs_HSL)
    color_picker_helper_1ch(pixel, roi, box, picked_color, picked_color_min, picked_color_max, image_cst, picker_cst, profile_info, channel);
  else if(dsc->channels == 1u && dsc->filters && dsc->filters != 9u)
    color_picker_helper_bayer(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max);
  else if(dsc->channels == 1u && dsc->filters && dsc->filters == 9u)
    color_picker_helper_xtrans(dsc, pixel, roi, box, picked_color, picked_color_min, picked_color_max);
  else
    dt_unreachable_codepath();
}

// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
