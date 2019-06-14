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

#pragma once

struct dt_iop_buffer_dsc_t;
struct dt_iop_roi_t;
struct dt_iop_order_iccprofile_info_t;
enum dt_iop_colorspace_type_t;

typedef enum dt_color_picker_channels_t
{
  DT_COLOR_PICKER_CHANNEL_0 = 0,
  DT_COLOR_PICKER_CHANNEL_1 = 1,
  DT_COLOR_PICKER_CHANNEL_2 = 2,
  DT_COLOR_PICKER_CHANNEL_GREY = 3
} dt_color_picker_channels_t;

void dt_color_picker_helper(const struct dt_iop_buffer_dsc_t *dsc, const float *const pixel,
                            const struct dt_iop_roi_t *roi, const int *const box, float *const picked_color,
                            float *const picked_color_min, float *const picked_color_max,
                            const enum dt_iop_colorspace_type_t image_cst, const enum dt_iop_colorspace_type_t picker_cst);

void dt_color_picker_helper_1c_gap(const struct dt_iop_buffer_dsc_t *dsc, const float *const pixel, const struct dt_iop_roi_t *roi,
                            const int *const box, float *const picked_color, float *const picked_color_min,
                            float *const picked_color_max, const enum dt_iop_colorspace_type_t image_cst, const enum dt_iop_colorspace_type_t picker_cst,
                            const struct dt_iop_order_iccprofile_info_t *const profile_info,
                            const dt_color_picker_channels_t channel);

void dt_color_picker_helper_1c(const struct dt_iop_buffer_dsc_t *dsc, const float *const pixel, const struct dt_iop_roi_t *roi,
                            const int *const box, float *const picked_color, float *const picked_color_min,
                            float *const picked_color_max, const enum dt_iop_colorspace_type_t image_cst, const enum dt_iop_colorspace_type_t picker_cst,
                            const struct dt_iop_order_iccprofile_info_t *const profile_info,
                            const dt_color_picker_channels_t channel);

// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
