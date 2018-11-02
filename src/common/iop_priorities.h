/*
    This file is part of darktable,
    copyright (c) 2018 edgardo hoszowski.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DT_IOP_PRIORITIES_H
#define DT_IOP_PRIORITIES_H

#include "develop/imageop.h"

typedef struct dt_iop_priority_entry_v0_t
{
  int priority;
  float iop_order;
  char operation[20];
  int default_enabled;
} dt_iop_priority_entry_v0_t;

typedef struct dt_iop_priority_entry_t
{
  int priority;
  float iop_order;
  char operation[20];
} dt_iop_priority_entry_t;

typedef struct _former_iop_priorities_t
{
  int priority;
  int multi_priority;
  char operation[20];
} _former_iop_priorities_t;

typedef struct dt_iop_priority_rule_t
{
  char op_prev[20];
  char op_next[20];
} dt_iop_priority_rule_t;

GList *dt_get_iop_priorities_v0();
int dt_iop_priorities_get_default_output_cst_v0(const char *op_name);

int dt_iop_priorities_get_default_input_colorspace(const char *op_name);
int dt_iop_priorities_get_default_output_colorspace(const char *op_name);

gint dt_sort_iop_so_by_priority(gconstpointer a, gconstpointer b);
gint dt_sort_iop_by_order(gconstpointer a, gconstpointer b);

// void dt_iop_priorities_set_so_iop_priorities();
GList *dt_iop_priorities_load_iop_priorities();
GList *dt_iop_priorities_load_iop_priority_rules();

int dt_iop_priorities_check_iop_priorities();

int dt_iop_priorities_get_iop_priority(const char *op_name);
int dt_iop_priorities_get_new_iop_multi_priority(dt_develop_t *dev, const char *op_name);

/** returns the iop_order before module_next if module can be moved */
float dt_get_iop_order_before_iop(dt_develop_t *dev, dt_iop_module_t *module, dt_iop_module_t *module_next,
                                  const int validate_order, const int log_error);
/** returns the iop_order after module_prev if module can be moved */
float dt_get_iop_order_after_iop(dt_develop_t *dev, dt_iop_module_t *module, dt_iop_module_t *module_prev,
                                 const int validate_order, const int log_error);

int dt_move_iop_before(dt_develop_t *dev, dt_iop_module_t *module, dt_iop_module_t *module_next,
                       const int validate_order, const int log_error);
int dt_move_iop_after(dt_develop_t *dev, dt_iop_module_t *module, dt_iop_module_t *module_prev,
                      const int validate_order, const int log_error);
float get_iop_default_order(const char *op_name);
dt_iop_module_t *dt_iop_priorities_get_last_instance(dt_develop_t *dev, dt_iop_module_t *module);

// void dt_iop_priorities_rebuild_order(dt_develop_t *dev);

void dt_iop_priorities_read_pipe(dt_develop_t *dev, const int imgid);
void dt_iop_priorities_write_pipe(dt_develop_t *dev, const int imgid);

/** transform image from cst_from to cst_to */
void dt_iop_transform_image_colorspace(struct dt_iop_module_t *self, float *image, const int width,
                                       const int height, const int cst_from, const int cst_to, int *converted_cst);
int dt_iop_transform_image_colorspace_cl(struct dt_iop_module_t *self, const int devid, cl_mem dev_img,
                                         const int width, const int height, const int cst_from, const int cst_to,
                                         int *converted_cst);

int dt_iop_priorities_check_priorities(dt_develop_t *dev, const char *msg);

dt_iop_module_t *dt_iop_priorities_get_iop_by_op(dt_develop_t *dev, const char *op_name);

#endif
