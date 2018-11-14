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

typedef struct dt_iop_priority_entry_t
{
  int priority;
  float iop_order;
  char operation[20];
} dt_iop_priority_entry_t;

typedef struct dt_iop_priority_rule_t
{
  char op_prev[20];
  char op_next[20];
} dt_iop_priority_rule_t;

/** returns the curren priorities version */
int get_priorities_version();
/** returns a list of dt_iop_priority_entry_t with the version 1 of iop priorities */
GList *dt_ioppr_get_iop_priorities_v1();
/** returns a list of dt_iop_priority_entry_t with the current version of iop priorities */
GList *dt_ioppr_get_iop_priorities();
/** returns the dt_iop_priority_entry_t of priorities_list with operation = op_name */
dt_iop_priority_entry_t *dt_ioppr_get_iop_priority_entry(GList *priorities_list, const char *op_name);
/** returns the iop_order from priorities_list list with operation = op_name */
float dt_ioppr_get_iop_order(GList *priorities_list, const char *op_name);

/** migrates *_priorities_list from _old_version to _new_version, it can be a prioty or module list, depending on mode */
int dt_ioppr_legacy_priorities(GList **_iop_list, GList **_priorities_list, const int _old_version, const int _new_version);

/** returns 1 if there's a module_so whithout a priority defined */
int dt_ioppr_check_so_iop_priorities(GList *iop_list, GList *priorities_list);

/* returns a list of dt_iop_priority_rule_t with the current priority rules */
GList *dt_ioppr_get_iop_priority_rules();

/** returns the default input/output colorspace for a module based on priorities_list */
int dt_ioppr_module_default_input_colorspace(GList *priorities_list, const char *op_name);
int dt_ioppr_module_default_output_colorspace(GList *priorities_list, const char *op_name);

/** returns a duplicate of iop_priorities */
GList *dt_ioppr_priorities_copy_deep(GList *iop_priorities);

/** sort two modules by (iop_order, priority) */
gint dt_sort_iop_by_order(gconstpointer a, gconstpointer b);

/** returns the iop_order before module_next if module can be moved */
float dt_ioppr_get_iop_order_before_iop(GList *iop_list, dt_iop_module_t *module, dt_iop_module_t *module_next,
                                  const int validate_order, const int log_error);
/** returns the iop_order after module_prev if module can be moved */
float dt_ioppr_get_iop_order_after_iop(GList *iop_list, dt_iop_module_t *module, dt_iop_module_t *module_prev,
                                 const int validate_order, const int log_error);

int dt_ioppr_move_iop_before(GList **_iop_list, dt_iop_module_t *module, dt_iop_module_t *module_next,
                       const int validate_order, const int log_error);
int dt_ioppr_move_iop_after(GList **_iop_list, dt_iop_module_t *module, dt_iop_module_t *module_prev,
                      const int validate_order, const int log_error);

/** transform image from cst_from to cst_to */
void dt_iop_transform_image_colorspace(struct dt_iop_module_t *self, float *image, const int width,
                                       const int height, const int cst_from, const int cst_to, int *converted_cst);
int dt_iop_transform_image_colorspace_cl(struct dt_iop_module_t *self, const int devid, cl_mem dev_img,
                                         const int width, const int height, const int cst_from, const int cst_to,
                                         int *converted_cst);

// for debug only
int dt_ioppr_check_priorities(dt_develop_t *dev, const int imgid, const char *msg);
void dt_ioppr_print_module_priorities(GList *history_list, const char *msg);
void dt_ioppr_print_history_priorities(GList *iop_list, const char *msg);
void dt_ioppr_print_priorities(GList *priorities_list, const char *msg);

#endif
