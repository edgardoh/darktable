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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
// clang-format off
#include "common/darktable.h"
#include "common/debug.h"
#include "common/file_location.h"
#include "common/iop_priorities.h"
#include "common/colorspaces_inline_conversions.h"
// clang-format on

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DT_IOP_PRIORITIES_VERSION (3)

static inline void _ioppr_insert_iop_after(GList **_priorities_list, const char *op_new, const char *op_previous);
static inline void _ioppr_insert_iop_before(GList **_priorities_list, const char *op_new, const char *op_next);
static inline void _ioppr_move_iop_after(GList **_priorities_list, const char *op_current, const char *op_prev, const int dont_move);
static inline void _ioppr_move_iop_before(GList **_priorities_list, const char *op_current, const char *op_next, const int dont_move);


// migrates *_priorities_list from old_version to the next version (version + 1)
static int _ioppr_legacy_priorities_step(GList **_priorities_list, const int old_version, const int dont_move)
{
  int new_version = -1;
  
  // version 1 --> 2
  if(old_version == 1)
  {
    _ioppr_insert_iop_after(_priorities_list, "dummy1", "demosaic");
    _ioppr_insert_iop_before(_priorities_list, "dummy2", "colorin");
    _ioppr_insert_iop_after(_priorities_list, "dummy3", "colorin");
    _ioppr_insert_iop_before(_priorities_list, "dummy4", "colorout");
    _ioppr_insert_iop_after(_priorities_list, "dummy5", "colorout");
    _ioppr_insert_iop_before(_priorities_list, "dummy6", "gamma");
    _ioppr_move_iop_before(_priorities_list, "lens", "exposure", dont_move);
    
    new_version = 2;
  }
  // version 2 --> 3
  else if(old_version == 2)
  {
    _ioppr_move_iop_after(_priorities_list, "tonecurve", "levels", dont_move);
    
    new_version = 3;
  }

  if(new_version <= 0)
    fprintf(stderr, "[_ioppr_legacy_priorities_step] missing step migrating from version %i\n", old_version);
  
  return new_version;
}

int get_priorities_version()
{
  return DT_IOP_PRIORITIES_VERSION;
}

static GList *_ioppr_buil_priorities_list(const dt_iop_priority_entry_t prior_entry[])
{
  GList *priorities = NULL;
  
  int i = 0;
  while(prior_entry[i].priority >= 0)
  {
    dt_iop_priority_entry_t *prior = calloc(1, sizeof(dt_iop_priority_entry_t));

    prior->priority = i + 1;
    prior->iop_order = (float)prior->priority;
    snprintf(prior->operation, sizeof(prior->operation), "%s", prior_entry[i].operation);

    priorities = g_list_append(priorities, prior);
    i++;
  }
  
  return priorities;
}

GList *dt_ioppr_get_iop_priorities_v1()
{
  GList *priorities = NULL;

  const dt_iop_priority_entry_t prior_entry[] = { { 0, 0.f, "rawprepare" },
                                                  { 0, 0.f, "invert" },
                                                  { 0, 0.f, "temperature" },
                                                  { 0, 0.f, "highlights" },
                                                  { 0, 0.f, "cacorrect" },
                                                  { 0, 0.f, "hotpixels" },
                                                  { 0, 0.f, "rawdenoise" },
                                                  { 0, 0.f, "demosaic" },
                                                  { 0, 0.f, "denoiseprofile" },
                                                  { 0, 0.f, "tonemap" },
                                                  { 0, 0.f, "exposure" },
                                                  { 0, 0.f, "spots" },
                                                  { 0, 0.f, "retouch" },
                                                  { 0, 0.f, "lens" },
                                                  { 0, 0.f, "ashift" },
                                                  { 0, 0.f, "liquify" },
                                                  { 0, 0.f, "rotatepixels" },
                                                  { 0, 0.f, "scalepixels" },
                                                  { 0, 0.f, "flip" },
                                                  { 0, 0.f, "clipping" },
                                                  { 0, 0.f, "graduatednd" },
                                                  { 0, 0.f, "basecurve" },
                                                  { 0, 0.f, "bilateral" },
                                                  { 0, 0.f, "profile_gamma" },
                                                  { 0, 0.f, "hazeremoval" },
                                                  { 0, 0.f, "colorin" },
                                                  { 0, 0.f, "colorreconstruct" },
                                                  { 0, 0.f, "colorchecker" },
                                                  { 0, 0.f, "defringe" },
                                                  { 0, 0.f, "equalizer" },
                                                  { 0, 0.f, "vibrance" },
                                                  { 0, 0.f, "colorbalance" },
                                                  { 0, 0.f, "colorize" },
                                                  { 0, 0.f, "colortransfer" },
                                                  { 0, 0.f, "colormapping" },
                                                  { 0, 0.f, "bloom" },
                                                  { 0, 0.f, "nlmeans" },
                                                  { 0, 0.f, "globaltonemap" },
                                                  { 0, 0.f, "shadhi" },
                                                  { 0, 0.f, "atrous" },
                                                  { 0, 0.f, "bilat" },
                                                  { 0, 0.f, "colorzones" },
                                                  { 0, 0.f, "lowlight" },
                                                  { 0, 0.f, "monochrome" },
                                                  { 0, 0.f, "filmic" },
                                                  { 0, 0.f, "colisa" },
                                                  { 0, 0.f, "zonesystem" },
                                                  { 0, 0.f, "tonecurve" },
                                                  { 0, 0.f, "levels" },
                                                  { 0, 0.f, "relight" },
                                                  { 0, 0.f, "colorcorrection" },
                                                  { 0, 0.f, "sharpen" },
                                                  { 0, 0.f, "lowpass" },
                                                  { 0, 0.f, "highpass" },
                                                  { 0, 0.f, "grain" },
                                                  { 0, 0.f, "colorcontrast" },
                                                  { 0, 0.f, "colorout" },
                                                  { 0, 0.f, "channelmixer" },
                                                  { 0, 0.f, "soften" },
                                                  { 0, 0.f, "vignette" },
                                                  { 0, 0.f, "splittoning" },
                                                  { 0, 0.f, "velvia" },
                                                  { 0, 0.f, "clahe" },
                                                  { 0, 0.f, "finalscale" },
                                                  { 0, 0.f, "overexposed" },
                                                  { 0, 0.f, "rawoverexposed" },
                                                  { 0, 0.f, "borders" },
                                                  { 0, 0.f, "watermark" },
                                                  { 0, 0.f, "dither" },
                                                  { 0, 0.f, "gamma" },
                                                  { -1, 0.f, "\0" }
  };

  priorities = _ioppr_buil_priorities_list(prior_entry);
  
  return priorities;
}

// returns the first priority entry that matches operation == op_name
dt_iop_priority_entry_t *dt_ioppr_get_iop_priority_entry(GList *priorities_list, const char *op_name)
{
  dt_iop_priority_entry_t *priority_entry = NULL;
  
  GList *priorities = g_list_first(priorities_list);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)(priorities->data);

    if(strcmp(prior->operation, op_name) == 0)
    {
      priority_entry = prior;
      break;
    }

    priorities = g_list_next(priorities);
  }

  return priority_entry;
}

// returns the iop_order asociated with the priority entry that matches operation == op_name
float dt_ioppr_get_iop_order(GList *priorities_list, const char *op_name)
{
  float iop_order = FLT_MAX;
  dt_iop_priority_entry_t *priority = dt_ioppr_get_iop_priority_entry(priorities_list, op_name);
  
  if(priority)
    iop_order = priority->iop_order;

  return iop_order;
}

// sets the priority on each entry of priorities_list based on the order of the list
// if rebuild_order == 1 it also sets the iop_order
static void _ioppr_rebuild_priorities(GList *priorities_list, const int rebuild_order)
{
  int i = 1;
  GList *priorities = g_list_first(priorities_list);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)(priorities->data);

    prior->priority = i;
    if(rebuild_order)
      prior->iop_order = (float)prior->priority;

    i++;
    priorities = g_list_next(priorities);
  }
}

// insert op_new before op_next on *_priorities_list
// it sets the iop_order but not the priority on op_new
static void _ioppr_insert_priority_before(GList **_priorities_list, const char *op_new, const char *op_next)
{
  GList *priorities_list = *_priorities_list;
  
  // check that the new priority don't exists on the list
  if(dt_ioppr_get_iop_priority_entry(priorities_list, op_new) == NULL)
  {
    // create a new priority entry
    dt_iop_priority_entry_t *priority_new = (dt_iop_priority_entry_t*)calloc(1, sizeof(dt_iop_priority_entry_t));
    snprintf(priority_new->operation, sizeof(priority_new->operation), "%s", op_new);
    
    //_ioppr_insert_priority_after(&priorities_list, new_entry, op_previous);
    
    // search for the previous one
    int position = 0;
    int found = 0;
    dt_iop_priority_entry_t *priority_prev = NULL;
    dt_iop_priority_entry_t *priority_next = NULL;
    GList *priorities = g_list_first(priorities_list);
    while(priorities)
    {
      dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)(priorities->data);
      if(strcmp(prior->operation, op_next) == 0)
      {
        priority_next = prior;
        found = 1;
        break;
      }
      priority_prev = prior;
      position++;
      
      priorities = g_list_next(priorities);
    }
    
    if(found)
    {
      // set the iop_order
      priority_new->iop_order = priority_prev->iop_order + (priority_next->iop_order - priority_prev->iop_order) / 2.f;
      
      // insert it on the proper order
      priorities_list = g_list_insert(priorities_list, priority_new, position);
    }
    else
      fprintf(stderr, "[_ioppr_insert_priority_before] module %s don't exists on priority list\n", op_next);
  }
  else
    fprintf(stderr, "[_ioppr_insert_priority_before] module %s already exists on priority list\n", op_new);

  *_priorities_list = priorities_list;
}

// insert priority_new after op_prev on *_priorities_list
// it updates the iop_order but not the priority on priority_new
static void _ioppr_insert_priority_after(GList **_priorities_list, const char *op_new, const char *op_prev)
{
  GList *priorities_list = *_priorities_list;
  
  // inserting after op_prev is the same as moving before the very next one after op_prev
  dt_iop_priority_entry_t *prior_next = NULL;
  GList *priorities = g_list_last(priorities_list);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)priorities->data;
    if(strcmp(prior->operation, op_prev) == 0) break;

    prior_next = prior;
    priorities = g_list_previous(priorities);
  }
  if(prior_next == NULL)
  {
    fprintf(
        stderr,
        "[_ioppr_insert_priority_after] can't find module previous to %s while moving %s after it\n",
        op_prev, op_new);
  }
  else
    _ioppr_insert_priority_before(&priorities_list, op_new, prior_next->operation);

  *_priorities_list = priorities_list;
}

static inline void _ioppr_insert_iop_after(GList **_priorities_list, const char *op_new, const char *op_previous)
{
  /*
  GList *priorities_list = *_priorities_list;

  // check that the new priority don't exists on the list
  if(dt_ioppr_get_iop_priority_entry(priorities_list, op_new) == NULL)
  {
    // create a new priority entry
    dt_iop_priority_entry_t *new_entry = (dt_iop_priority_entry_t*)calloc(1, sizeof(dt_iop_priority_entry_t));
    snprintf(new_entry->operation, sizeof(new_entry->operation), "%s", op_new);
    
    _ioppr_insert_priority_after(&priorities_list, new_entry, op_previous);
  }
  else
    fprintf(stderr, "[_ioppr_insert_iop_after] module %s already exists on priority list\n", op_new);

  *_priorities_list = priorities_list;*/
  _ioppr_insert_priority_after(_priorities_list, op_new, op_previous);
}

static inline void _ioppr_insert_iop_before(GList **_priorities_list, const char *op_new, const char *op_next)
{
  /*
  GList *priorities_list = *_priorities_list;

  // check that the new priority don't exists on the list
  if(dt_ioppr_get_iop_priority_entry(priorities_list, op_new) == NULL)
  {
    // create a new priority entry
    dt_iop_priority_entry_t *new_entry = (dt_iop_priority_entry_t*)calloc(1, sizeof(dt_iop_priority_entry_t));
    snprintf(new_entry->operation, sizeof(new_entry->operation), "%s", op_new);
    
    _ioppr_insert_priority_before(&priorities_list, new_entry, op_next);
  }
  else
    fprintf(stderr, "[_ioppr_insert_iop_before] module %s already exists on priority list\n", op_new);

  *_priorities_list = priorities_list;*/
  _ioppr_insert_priority_before(_priorities_list, op_new, op_next);
}

static void _ioppr_move_priority_before(GList **_priorities_list, const char *op_current, const char *op_next, const int dont_move)
{
  if(dont_move) return;
  
  GList *priorities_list = *_priorities_list;
  
  int position = 0;
  int found = 0;
  dt_iop_priority_entry_t *priority_prev = NULL;
  dt_iop_priority_entry_t *priority_next = NULL;
  dt_iop_priority_entry_t *priority_current = NULL;
  GList *priorities_current = NULL;
  
  // search for the current one
  GList *priorities = g_list_first(priorities_list);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)(priorities->data);
    if(strcmp(prior->operation, op_current) == 0)
    {
      priorities_current = priorities;
      priority_current = prior;
      found = 1;
      break;
    }
    
    priorities = g_list_next(priorities);
  }

  if(found)
  {
    // remove it from the list
    priorities_list = g_list_remove_link(priorities_list, priorities_current);
  }
  else
    fprintf(stderr, "[_ioppr_move_iop_before] current module %s don't exists on priority list\n", op_current);

  // search for the previous and next one
  if(found)
  {
    found = 0;
    priorities = g_list_first(priorities_list);
    while(priorities)
    {
      dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)(priorities->data);
      if(strcmp(prior->operation, op_next) == 0)
      {
        priority_next = prior;
        found = 1;
        break;
      }
      priority_prev = prior;
      position++;
      
      priorities = g_list_next(priorities);
    }
  }
  
  if(found)
  {
    // set the iop_order
    priority_current->iop_order = priority_prev->iop_order + (priority_next->iop_order - priority_prev->iop_order) / 2.f;
    
    // insert it on the proper order
    priorities_list = g_list_insert(priorities_list, priority_current, position);
  }
  else
    fprintf(stderr, "[_ioppr_move_iop_before] next module %s don't exists on priority list\n", op_next);

  *_priorities_list = priorities_list;
}

static void _ioppr_move_priority_after(GList **_priorities_list, const char *op_current, const char *op_prev, const int dont_move)
{
  if(dont_move) return;
  
  GList *priorities_list = *_priorities_list;
  
  // moving after op_prev is the same as moving before the very next one after op_prev
  dt_iop_priority_entry_t *prior_next = NULL;
  GList *priorities = g_list_last(priorities_list);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)priorities->data;
    if(strcmp(prior->operation, op_prev) == 0) break;

    prior_next = prior;
    priorities = g_list_previous(priorities);
  }
  if(prior_next == NULL)
  {
    fprintf(
        stderr,
        "[_ioppr_move_iop_after] can't find module previous to %s while moving %s after it\n",
        op_prev, op_current);
  }
  else
    _ioppr_move_iop_before(&priorities_list, op_current, prior_next->operation, dont_move);

  *_priorities_list = priorities_list;
}

static inline void _ioppr_move_iop_before(GList **_priorities_list, const char *op_current, const char *op_next, const int dont_move)
{
  _ioppr_move_priority_before(_priorities_list, op_current, op_next, dont_move);
}

static inline void _ioppr_move_iop_after(GList **_priorities_list, const char *op_current, const char *op_prev, const int dont_move)
{
  _ioppr_move_priority_after(_priorities_list, op_current, op_prev, dont_move);
}

/*
  migrates *_iop_list and *_priorities_list from version _old_version to _new_version
  version == 0 means current version
  for both list the result will be:
  priority matching the current version
  iop_order matching _new_version
  list sorted by (iop_order, priority)
*/
int dt_ioppr_legacy_priorities(GList **_iop_list, GList **_priorities_list, const int _old_version, const int _new_version)
{
  GList *iop_list = *_iop_list;
  GList *priorities_list = *_priorities_list;
  int old_version = (_old_version == 0) ? get_priorities_version(): _old_version;
  const int new_version = (_new_version == 0) ? get_priorities_version(): _new_version;
  int old_version_priority = old_version;
  
  // first we will deal with the priorities
  // we'll just start from scratch
  g_list_free_full(priorities_list, free);
  priorities_list = dt_ioppr_get_iop_priorities_v1();
  old_version_priority = 1;
  
  // we have a list of version 1 of priorities
  // now we want to get a list but for version new_version
  // this will end with a list with priority and iop_order for new_version
  // sorted by (iop_order, priority)
  // this can change both priority and iop_order if new modules are added or existing ones moved
  while(old_version_priority < new_version && old_version_priority > 0)
  {
    old_version_priority = _ioppr_legacy_priorities_step(&priorities_list, old_version_priority/*, DT_IOPPR_MIGRATE_PRIORITY*/, FALSE);
  }
  _ioppr_rebuild_priorities(priorities_list, TRUE);

  // now deal with the modules
  // we will replace the iop_order with the migrated priorities so they match version new_version
  // priority already has to be the current version, so we don't need to worry about it
  // we'll change iop_order only for base modules, multi-instances will be flaged as unused with FLT_MAX
  GList *modules = g_list_first(iop_list);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)(modules->data);

    if(mod->multi_priority == 0)
    {
      const float iop_order = dt_ioppr_get_iop_order(priorities_list, mod->op);
      if(iop_order > 0.f)
      {
        mod->iop_order = iop_order;
      }
      else
      {
        // this can happen if the module was created after version new_version
        // in this case the next step will fix it
        mod->iop_order = FLT_MAX;
      }
    }
    // muti-instances will be set by read history
    else
    {
      mod->iop_order = FLT_MAX;
    }
    
    modules = g_list_next(modules);
  }
  // now we have a list of modules for version new_version

  // now we want to add any module created after version new_version
  // but we don't want to change the iop_order, that has to stay at version new_version
  // so we'll insert the new modules and set the iop_order to match the position on the pipe that they have on current version
  // priority will be set to the new value (for version new_version), but iop_order will not be changed
  // for existing modules
  while(old_version_priority < get_priorities_version() && old_version_priority > 0)
  {
    old_version_priority = _ioppr_legacy_priorities_step(&priorities_list, old_version_priority/*, DT_IOPPR_MIGRATE_PRIORITY*/, TRUE);
  }
  _ioppr_rebuild_priorities(priorities_list, FALSE);
  
  // now that we have a list of priorities for version new_version but with all new modules
  // we take care of the iop_order of new modules on iop list
  modules = g_list_first(iop_list);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)(modules->data);

    if(mod->multi_priority == 0 && mod->iop_order == FLT_MAX)
    {
      const float iop_order = dt_ioppr_get_iop_order(priorities_list, mod->op);
      if(iop_order > 0.f)
      {
        mod->iop_order = iop_order;
      }
      else
        fprintf(stderr, "[dt_ioppr_legacy_priorities] can't find iop_order for module %s\n", mod->op);
    }
    
    modules = g_list_next(modules);
  }
  // we need to set the right order
  iop_list = g_list_sort(iop_list, dt_sort_iop_by_order);

  *_iop_list = iop_list;
  *_priorities_list = priorities_list;

  return new_version;
}

GList *dt_ioppr_get_iop_priorities()
{
  GList *priorities = dt_ioppr_get_iop_priorities_v1();
  
  int old_version = 1;
  while(old_version < get_priorities_version() && old_version > 0)
  {
    old_version = _ioppr_legacy_priorities_step(&priorities, old_version/*, DT_IOPPR_MIGRATE_PRIORITY*/, FALSE);
  }
  
  // and rebuild the priorities
  _ioppr_rebuild_priorities(priorities, TRUE);

  if(old_version != get_priorities_version())
  {
    fprintf(stderr, "[dt_ioppr_get_iop_priorities] error migrating priorities, current verion is %i, migrated version is %i\n", 
        get_priorities_version(), old_version);
  }
  
  return priorities;
}

// check if all modules on iop_list have a priority defined for version
int dt_ioppr_check_so_iop_priorities(GList *iop_list, GList *priorities_list)
{
  int priority_missing = 0;

  // check if all the modules have their priority assigned
  GList *modules = g_list_first(iop_list);
  while(modules)
  {
    dt_iop_module_so_t *mod = (dt_iop_module_so_t *)(modules->data);
    
    dt_iop_priority_entry_t *entry = dt_ioppr_get_iop_priority_entry(priorities_list, mod->op);
    if(entry == NULL)
    {
      priority_missing = 1;
      fprintf(stderr, "[dt_ioppr_check_so_iop_priorities] missing priority for module %s\n", mod->op);
    }
    modules = g_list_next(modules);
  }

  return priority_missing;
}

void dt_ioppr_print_module_priorities(GList *iop_list, const char *msg)
{
  GList *modules = g_list_first(iop_list);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)(modules->data);
    
    fprintf(stderr, "[%s] module %s %s priority=%i, multi_priority=%i, iop_order=%f\n", msg, mod->op, mod->multi_name, mod->priority, mod->multi_priority, mod->iop_order);

    modules = g_list_next(modules);
  }
}

void dt_ioppr_print_history_priorities(GList *history_list, const char *msg)
{
  GList *history = g_list_first(history_list);
  while(history)
  {
    dt_dev_history_item_t *hist = (dt_dev_history_item_t *)(history->data);
    
    fprintf(stderr, "[%s] module %s %s multi_priority=%i, iop_order=%f\n", msg, hist->op_name, hist->multi_name, hist->multi_priority, hist->iop_order);

    history = g_list_next(history);
  }
}

void dt_ioppr_print_priorities(GList *priorities_list, const char *msg)
{
  GList *priorities = g_list_first(priorities_list);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)(priorities->data);
    
    fprintf(stderr, "[%s] priority %s priority=%i, iop_order=%f\n", msg, prior->operation, prior->priority, prior->iop_order);

    priorities = g_list_next(priorities);
  }
}

int dt_ioppr_module_default_input_colorspace(GList *priorities_list, const char *op_name)
{
  int cst = iop_cs_gamma_rgb;
  const float iop_order_module = dt_ioppr_get_iop_order(priorities_list, op_name);
  
  if(iop_order_module < 0.f)
    fprintf(stderr, "[dt_ioppr_module_default_input_colorspace] invalid default iop_order for module %s(%f)\n", op_name, iop_order_module);
  else if(iop_order_module <= dt_ioppr_get_iop_order(priorities_list, "demosaic"))
    cst = iop_cs_RAW;
  else if(iop_order_module <= dt_ioppr_get_iop_order(priorities_list, "colorin"))
    cst = iop_cs_linear_rgb;
  else if(iop_order_module <= dt_ioppr_get_iop_order(priorities_list, "colorout"))
    cst = iop_cs_Lab;
  else if(iop_order_module <= dt_ioppr_get_iop_order(priorities_list, "gamma"))
    cst = iop_cs_gamma_rgb;
  else
    fprintf(stderr, "[dt_ioppr_module_default_input_colorspace] invalid default iop_order for module %s(%f)\n", op_name, iop_order_module);

  return cst;
}

int dt_ioppr_module_default_output_colorspace(GList *priorities_list, const char *op_name)
{
  int cst = iop_cs_gamma_rgb;
  const float iop_order_module = dt_ioppr_get_iop_order(priorities_list, op_name);
  
  if(iop_order_module < 0.f)
    fprintf(stderr, "[dt_ioppr_module_default_output_colorspace] invalid default iop_order for module %s(%f)\n", op_name, iop_order_module);
  else if(iop_order_module < dt_ioppr_get_iop_order(priorities_list, "demosaic"))
    cst = iop_cs_RAW;
  else if(iop_order_module < dt_ioppr_get_iop_order(priorities_list, "colorin"))
    cst = iop_cs_linear_rgb;
  else if(iop_order_module < dt_ioppr_get_iop_order(priorities_list, "colorout"))
    cst = iop_cs_Lab;
  else if(iop_order_module <= dt_ioppr_get_iop_order(priorities_list, "gamma"))
    cst = iop_cs_gamma_rgb;
  else
    fprintf(stderr, "[dt_ioppr_module_default_output_colorspace] invalid default iop_order for module %s(%f)\n", op_name, iop_order_module);

  return cst;
}

static void *_dup_priority_entry(const void *src, gpointer data)
{
  dt_iop_priority_entry_t *scr_entry = (dt_iop_priority_entry_t *)src;
  dt_iop_priority_entry_t *new_entry = malloc(sizeof(dt_iop_priority_entry_t));
  memcpy(new_entry, scr_entry, sizeof(dt_iop_priority_entry_t));
  return (void *)new_entry;
}

GList *dt_ioppr_priorities_copy_deep(GList *iop_priorities)
{
  return (GList *)g_list_copy_deep(iop_priorities, _dup_priority_entry, NULL);
}

// helper to sort a GList by (iop_order, priority)
gint dt_sort_iop_by_order(gconstpointer a, gconstpointer b)
{
  const dt_iop_module_t *am = (const dt_iop_module_t *)a;
  const dt_iop_module_t *bm = (const dt_iop_module_t *)b;
  if(am->iop_order > bm->iop_order) return 1;
  if(am->iop_order < bm->iop_order) return -1;
  return (am->priority - bm->priority);
}

GList *dt_ioppr_get_iop_priority_rules()
{
  GList *rules = NULL;

  const dt_iop_priority_rule_t rule_entry[] = { { "rawprepare", "invert" },
                                                { "invert", "temperature" },
                                                { "temperature", "highlights" },
                                                { "highlights", "cacorrect" },
                                                { "cacorrect", "hotpixels" },
                                                { "hotpixels", "rawdenoise" },
                                                { "rawdenoise", "demosaic" },
                                                { "demosaic", "colorin" },
                                                { "colorin", "colorout" },
                                                { "colorout", "gamma" },
                                                { "\0", "\0" } };

  int i = 0;
  while(rule_entry[i].op_prev[0])
  {
    dt_iop_priority_rule_t *rule = calloc(1, sizeof(dt_iop_priority_rule_t));

    snprintf(rule->op_prev, sizeof(rule->op_prev), "%s", rule_entry[i].op_prev);
    snprintf(rule->op_next, sizeof(rule->op_next), "%s", rule_entry[i].op_next);
    // printf("[_get_iop_priority_rules] %s %s\n", rule->op_prev, rule->op_next);
    rules = g_list_append(rules, rule);
    i++;
  }

  return rules;
}

int dt_iop_priorities_get_new_iop_multi_priority(dt_develop_t *dev, const char *op_name)
{
  int multi_priority_new = -1;
  GList *modules = g_list_first(dev->iop);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)(modules->data);

    if(strcmp(mod->op, op_name) == 0)
    {
      multi_priority_new = MAX(multi_priority_new, mod->multi_priority);
    }
    modules = g_list_next(modules);
  }
  return (multi_priority_new + 1);
}

dt_iop_module_t *dt_iop_priorities_get_last_instance(dt_develop_t *dev, dt_iop_module_t *module)
{
  dt_iop_module_t *last = NULL;
  GList *modules = g_list_last(dev->iop);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
    if(strcmp(mod->op, module->op) == 0)
    {
      last = mod;
      break;
    }
    modules = g_list_previous(modules);
  }
  return last;
}

// if module can be placed before than module_next on the pipe
// it returns the new iop_order
// if it cannot be placed it returns -1.f
// this assums that the order is always positive
float dt_ioppr_get_iop_order_before_iop(GList *iop_list, dt_iop_module_t *module, dt_iop_module_t *module_next,
                                  const int validate_order, const int log_error)
{
  float iop_order = -1.f;
  int module_found = 0;
  int rule_found = 0;

  // this can happen when migrating from an old version
  // let's handle it here
  if(module->iop_order == module_next->iop_order)
  {
    // we need to find the module before this two
    GList *modules = g_list_first(iop_list);
    dt_iop_module_t *prev = NULL;
    while(modules)
    {
      dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;

      if(mod == module || mod == module_next)
      {
        module_found = 1;
        break;
      }
      prev = mod;
      modules = g_list_next(modules);
    }

    if(module_found)
    {
      if(prev)
      {
        // calculate new iop_order
        iop_order = module_next->iop_order - (module_next->iop_order - prev->iop_order) / 2.f;
        /* fprintf(stderr,
               "[dt_ioppr_get_iop_order_before_iop] 1-calculated new iop_order=%f "
            "for %s(%f) between %s(%f) and %s(%f)\n",
               iop_order, module->op, module->iop_order, module_next->op, module_next->iop_order, prev->op,
               prev->iop_order); */
      }
      else
      {
        // those are the first modules, this can't happen
        if(log_error)
          fprintf(stderr, "[dt_ioppr_get_iop_order_before_iop] %s(%f) and %s(%f) are the first modules!!!\n", module->op,
                  module->iop_order, module_next->op, module_next->iop_order);
      }
    }
  }
  // module is before on the pipe
  // move it up
  else if(module->iop_order < module_next->iop_order)
  {
    GList *modules = g_list_first(iop_list);
    dt_iop_module_t *prev = NULL;
    while(modules)
    {
      dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
      if(mod == module)
      {
        module_found = 1; // we found the module
      }
      else if(module_found)
      {
        // we reach the next module, all good
        // calculate the new iop_order
        if(mod == module_next)
        {
          // this is already the previous module!
          if(module == prev)
          {
            if(log_error)
              fprintf(stderr, "[dt_ioppr_get_iop_order_before_iop] %s(%f) is already previous to %s(%f)\n", module->op,
                      module->iop_order, prev->op, prev->iop_order);
            break;
          }

          // calculate new iop_order
          iop_order = module_next->iop_order - (module_next->iop_order - prev->iop_order) / 2.f;
          /* fprintf(stderr,
                 "[dt_ioppr_get_iop_order_before_iop] 2-calculated new iop_order=%f "
              "for %s(%f) between %s(%f) and %s(%f)\n",
                 iop_order, module->op, module->iop_order, module_next->op, module_next->iop_order, mod->op,
                 mod->iop_order); */
          break;
        }
        else if(validate_order)
        {
          // check if module can be moved around this one
          if(mod->flags() & IOP_FLAGS_FENCE)
          {
            if(log_error)
              fprintf(stderr, "[dt_ioppr_get_iop_order_before_iop] can't move %s(%f) pass %s(%f)\n", module->op,
                      module->iop_order, mod->op, mod->iop_order);
            break;
          }

          // is there a rule about swapping this two?
          // int rule_found = 0;
          GList *rules = g_list_first(darktable.iop_priority_rules);
          while(rules)
          {
            dt_iop_priority_rule_t *rule = (dt_iop_priority_rule_t *)rules->data;

            if(strcmp(module->op, rule->op_prev) == 0 && strcmp(mod->op, rule->op_next) == 0)
            {
              if(log_error)
                fprintf(stderr,
                        "[dt_ioppr_get_iop_order_before_iop] found rule %s %s while moving %s(%f) before %s(%f)\n",
                        rule->op_prev, rule->op_next, module->op, module->iop_order, module_next->op,
                        module_next->iop_order);
              rule_found = 1;
              break;
            }

            rules = g_list_next(rules);
          }
          if(rule_found) break;
        }
      }

      prev = mod;
      modules = g_list_next(modules);
    }
  }
  else
  // module is next on the pipe
  // move it down
  {
    GList *modules = g_list_last(iop_list);
    int next_module_found = 0;
    while(modules)
    {
      dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
      if(mod == module)
      {
        module_found = 1; // we found the module
        modules = g_list_previous(modules);
        continue;
      }
      if(module_found)
      {
        if(next_module_found)
        {
          // we reach the module before module_next, all good
          // calculate the new iop_order
          iop_order = mod->iop_order + (module_next->iop_order - mod->iop_order) / 2.f;
          /* fprintf(stderr,
                 "[dt_ioppr_get_iop_order_before_iop] 3-calculated new iop_order=%f "
              "for %s %s(%f) between %s %s(%f) and %s %s(%f)\n",
                 iop_order, module->op, module->multi_name, module->iop_order,
                 mod->op, mod->multi_name, mod->iop_order,
                 module_next->op, module_next->multi_name, module_next->iop_order); */

          break;
        }

        if(mod == module_next)
        {
          next_module_found = 1; // we found the next module
          // modules = g_list_previous(modules);
          // continue;
        }
        if(validate_order)
        {
          // check if module can be moved around this one
          if(mod->flags() & IOP_FLAGS_FENCE)
          {
            if(log_error)
              fprintf(stderr, "[dt_ioppr_get_iop_order_before_iop] can't move %s(%f) pass %s(%f)\n", module->op,
                      module->iop_order, mod->op, mod->iop_order);
            break;
          }

          // is there a rule about swapping this two?
          // int rule_found = 0;
          GList *rules = g_list_first(darktable.iop_priority_rules);
          while(rules)
          {
            dt_iop_priority_rule_t *rule = (dt_iop_priority_rule_t *)rules->data;

            if(strcmp(mod->op, rule->op_prev) == 0 && strcmp(module->op, rule->op_next) == 0)
            {
              if(log_error)
                fprintf(stderr,
                        "[dt_ioppr_get_iop_order_before_iop] found rule %s %s while moving %s(%f) before %s(%f)\n",
                        rule->op_prev, rule->op_next, module->op, module->iop_order, module_next->op,
                        module_next->iop_order);
              rule_found = 1;
              break;
            }

            rules = g_list_next(rules);
          }
          if(rule_found) break;
        }
      }

      modules = g_list_previous(modules);
    }
    if(!next_module_found && !rule_found)
      if(log_error)
        fprintf(stderr, "[dt_ioppr_get_iop_order_before_iop] can't find next module %s(%f) while moving %s(%f)\n",
                module_next->op, module_next->iop_order, module->op, module->iop_order);
  }

  if(!module_found) fprintf(stderr, "[dt_ioppr_get_iop_order_before_iop] can't find module %s\n", module->op);

  return iop_order;
}

// if module can be placed after than module_prev on the pipe
// it returns the new iop_order
// if it cannot be placed it returns -1.f
// this assums that the order is always positive
float dt_ioppr_get_iop_order_after_iop(GList *iop_list, dt_iop_module_t *module, dt_iop_module_t *module_prev,
                                 const int validate_order, const int log_error)
{
  // dt_develop_t *dev = module->dev;
  float iop_order = -1.f;

  // moving after module_prev is the same as moving before the very next one after module_prev
  GList *modules = g_list_last(iop_list);
  dt_iop_module_t *module_next = NULL;
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
    if(mod == module_prev) break;

    module_next = mod;
    modules = g_list_previous(modules);
  }
  if(module_next == NULL)
  {
    fprintf(
        stderr,
        "[dt_ioppr_get_iop_order_after_iop] can't find module previous to %s %s(%f) while moving %s %s(%f) after it\n",
        module_prev->op, module_prev->multi_name, module_prev->iop_order, module->op, module->multi_name,
        module->iop_order);
  }
  else
    iop_order = dt_ioppr_get_iop_order_before_iop(iop_list, module, module_next, validate_order, log_error);

  return iop_order;
}

// changes the module->iop_order so it comes before in the pipe than module_next
// sort dev->iop to reflect the changes
// return 1 if iop_order is changed, 0 otherwise
int dt_ioppr_move_iop_before(GList **_iop_list, dt_iop_module_t *module, dt_iop_module_t *module_next,
                       const int validate_order, const int log_error)
{
  GList *iop_list = *_iop_list;
  int moved = 0;

  // dt_ioppr_check_priorities(dev, "dt_ioppr_move_iop_before begin");
  // printf("[dt_ioppr_move_iop_before] moving %s(%f) before %s(%f)\n", module->op, module->iop_order, module_next->op,
  // module_next->iop_order);

  const float iop_order = dt_ioppr_get_iop_order_before_iop(iop_list, module, module_next, validate_order, log_error);
  // printf("[dt_ioppr_move_iop_before] new iop_order for %s is %f\n", module->op, iop_order);

  if(iop_order >= 0.f)
  {
    module->iop_order = iop_order;
    iop_list = g_list_sort(iop_list, dt_sort_iop_by_order);
    moved = 1;
  }
  else if(log_error)
    fprintf(stderr, "[dt_ioppr_move_iop_before] module %s is already before %s\n", module->op, module_next->op);

  // dt_ioppr_check_priorities(dev, "dt_ioppr_move_iop_before end");

  *_iop_list = iop_list;
  
  return moved;
}

// changes the module->iop_order so it comes after in the pipe than module_prev
// sort dev->iop to reflect the changes
// return 1 if iop_order is changed, 0 otherwise
int dt_ioppr_move_iop_after(GList **_iop_list, dt_iop_module_t *module, dt_iop_module_t *module_prev,
                      const int validate_order, const int log_error)
{
  // printf("[dt_ioppr_move_iop_after] moving %s(%f) after %s(%f)\n", module->op, module->iop_order, module_prev->op,
  // module_prev->iop_order);
  GList *iop_list = *_iop_list;
  int moved = 0;

  // dt_ioppr_check_priorities(dev, "dt_ioppr_move_iop_after begin");

  const float iop_order = dt_ioppr_get_iop_order_after_iop(iop_list, module, module_prev, validate_order, log_error);
  // printf("[dt_ioppr_move_iop_after] new iop_order for %s is %f\n", module->op, iop_order);
  if(iop_order >= 0.f)
  {
    module->iop_order = iop_order;
    iop_list = g_list_sort(iop_list, dt_sort_iop_by_order);
    moved = 1;
  }
  else if(log_error)
    fprintf(stderr, "[dt_ioppr_move_iop_after] module %s is already after %s\n", module->op, module_prev->op);

  // dt_ioppr_check_priorities(dev, "dt_ioppr_move_iop_after end");

  *_iop_list = iop_list;
  
  return moved;
}

static dt_dev_history_item_t *_search_history_by_op_multi_priority_ext(dt_develop_t *dev, const char *op_name,
                                                                       const int multi_priority)
{
  int cnt = 0;
  dt_dev_history_item_t *ret_hist = NULL;
  GList *history = g_list_first(dev->history);
  while(history)
  {
    dt_dev_history_item_t *hist = (dt_dev_history_item_t *)history->data;

    // always look at active history
    if(cnt >= dev->history_end) break;

    if(strcmp(op_name, hist->module->op) == 0 && multi_priority == hist->multi_priority)
    {
      ret_hist = hist;
      break;
    }

    cnt++;
    history = g_list_next(history);
  }
  return ret_hist;
}

static inline dt_dev_history_item_t *_search_history_by_op_multi_priority(dt_develop_t *dev,
                                                                          dt_iop_module_t *module)
{
  return _search_history_by_op_multi_priority_ext(dev, module->op, module->multi_priority);
}

static GList *_get_fence_modules_list(GList *iop_list)
{
  GList *fences = NULL;
  GList *modules = g_list_first(iop_list);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;

    if(mod->flags() & IOP_FLAGS_FENCE)
    {
      fences = g_list_append(fences, mod);
    }

    modules = g_list_next(modules);
  }
  return fences;
}

static void _ioppr_check_rules(GList *iop_list, const int imgid, const char *msg)
{
  GList *modules = NULL;

  // check for IOP_FLAGS_FENCE on each module
  // create a list of fences modules
  GList *fences = _get_fence_modules_list(iop_list);

  // check if each module is between the fences
  modules = g_list_first(iop_list);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
    if(mod->iop_order == FLT_MAX)
    {
      modules = g_list_next(modules);
      continue;
    }

    dt_iop_module_t *fence_prev = NULL;
    dt_iop_module_t *fence_next = NULL;

    GList *mod_fences = g_list_first(fences);
    while(mod_fences)
    {
      dt_iop_module_t *mod_fence = (dt_iop_module_t *)mod_fences->data;

      // mod should be before this fence
      if(mod->priority < mod_fence->priority)
      {
        if(fence_next == NULL)
          fence_next = mod_fence;
        else if(mod_fence->priority < fence_next->priority)
          fence_next = mod_fence;
      }
      // mod should be after this fence
      else if(mod->priority > mod_fence->priority)
      {
        if(fence_prev == NULL)
          fence_prev = mod_fence;
        else if(mod_fence->priority > fence_prev->priority)
          fence_prev = mod_fence;
      }

      mod_fences = g_list_next(mod_fences);
    }

    // now check if mod is between the fences
    if(fence_next && mod->iop_order > fence_next->iop_order)
    {
      fprintf(stderr, "[_ioppr_check_rules] found fence %s %s module %s %s(%f) is after %s %s(%f) image %i (%s)\n",
              fence_next->op, fence_next->multi_name, mod->op, mod->multi_name, mod->iop_order, fence_next->op,
              fence_next->multi_name, fence_next->iop_order, imgid, msg);
    }
    if(fence_prev && mod->iop_order < fence_prev->iop_order)
    {
      fprintf(stderr, "[_ioppr_check_rules] found fence %s %s module %s %s(%f) is before %s %s(%f) image %i (%s)\n",
              fence_prev->op, fence_prev->multi_name, mod->op, mod->multi_name, mod->iop_order, fence_prev->op,
              fence_prev->multi_name, fence_prev->iop_order, imgid, msg);
    }


    modules = g_list_next(modules);
  }

  // for each module check if it doesn't break a rule
  modules = g_list_first(iop_list);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
    if(mod->iop_order == FLT_MAX)
    {
      modules = g_list_next(modules);
      continue;
    }

    // we have a module, now check each rule
    GList *rules = g_list_first(darktable.iop_priority_rules);
    while(rules)
    {
      dt_iop_priority_rule_t *rule = (dt_iop_priority_rule_t *)rules->data;

      // mod must be before rule->op_next
      if(strcmp(mod->op, rule->op_prev) == 0)
      {
        // check if there's a rule->op_next module before mod
        GList *modules_prev = g_list_previous(modules);
        while(modules_prev)
        {
          dt_iop_module_t *mod_prev = (dt_iop_module_t *)modules_prev->data;

          if(strcmp(mod_prev->op, rule->op_next) == 0)
          {
            fprintf(stderr, "[_ioppr_check_rules] found rule %s %s module %s %s(%f) is after %s %s(%f) image %i (%s)\n",
                    rule->op_prev, rule->op_next, mod->op, mod->multi_name, mod->iop_order, mod_prev->op,
                    mod_prev->multi_name, mod_prev->iop_order, imgid, msg);
          }

          modules_prev = g_list_previous(modules_prev);
        }
      }
      // mod must be after rule->op_prev
      else if(strcmp(mod->op, rule->op_next) == 0)
      {
        // check if there's a rule->op_prev module after mod
        GList *modules_next = g_list_next(modules);
        while(modules_next)
        {
          dt_iop_module_t *mod_next = (dt_iop_module_t *)modules_next->data;

          if(strcmp(mod_next->op, rule->op_prev) == 0)
          {
            fprintf(stderr, "[_ioppr_check_rules] found rule %s %s module %s %s(%f) is before %s %s(%f) image %i (%s)\n",
                    rule->op_prev, rule->op_next, mod->op, mod->multi_name, mod->iop_order, mod_next->op,
                    mod_next->multi_name, mod_next->iop_order, imgid, msg);
          }

          modules_next = g_list_next(modules_next);
        }
      }

      rules = g_list_next(rules);
    }

    modules = g_list_next(modules);
  }

  if(fences) g_list_free(fences);
}

int dt_ioppr_check_priorities(dt_develop_t *dev, const int imgid, const char *msg)
{
  int priorities_ok = 1;

  // check if gamma is the last iop
  {
    GList *modules = g_list_last(dev->iop);
    while(modules)
    {
      dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
      if(mod->iop_order == FLT_MAX)
      {
        if(mod->enabled)
        {
          priorities_ok = 0;
          fprintf(stderr, "[dt_ioppr_check_priorities] module not used but enabled!! %s %s(%f) image %i (%s)\n",
                  mod->op, mod->multi_name, mod->iop_order,imgid, msg);
        }
        modules = g_list_previous(dev->iop);
      }
      else
        break;
    }
    if(modules)
    {
      dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;

      if(strcmp(mod->op, "gamma") != 0)
      {
        priorities_ok = 0;
        fprintf(stderr, "[dt_ioppr_check_priorities] gamma is not the last iop, last is %s %s(%f) image %i (%s)\n",
                mod->op, mod->multi_name, mod->iop_order,imgid, msg);
      }
    }
    else
    {
      // fprintf(stderr, "[dt_ioppr_check_priorities] dev->iop is empty image %i (%s)\n",imgid, msg);
    }
  }

  // check if there's duplicate or out-of-order iop_order
  {
    dt_iop_module_t *mod_prev = NULL;
    GList *modules = g_list_first(dev->iop);
    while(modules)
    {
      dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
      if(mod->iop_order != FLT_MAX)
      {
        if(mod_prev)
        {
          if(mod->iop_order < mod_prev->iop_order)
          {
            priorities_ok = 0;
            fprintf(stderr,
                    "[dt_ioppr_check_priorities] module %s %s(%f) should be after %s %s(%f) image %i (%s)\n",
                    mod->op, mod->multi_name, mod->iop_order, mod_prev->op, mod_prev->multi_name,
                    mod_prev->iop_order, imgid, msg);
          }
          else if(mod->iop_order == mod_prev->iop_order)
          {
            priorities_ok = 0;
            fprintf(
                stderr,
                "[dt_ioppr_check_priorities] module %s %s(%i)(%f) and %s %s(%i)(%f) has the same order image %i (%s)\n",
                mod->op, mod->multi_name, mod->multi_priority, mod->iop_order, mod_prev->op, 
                mod_prev->multi_name, mod_prev->multi_priority, mod_prev->iop_order, imgid, msg);
          }
        }
      }
      mod_prev = mod;
      modules = g_list_next(modules);
    }
  }

  _ioppr_check_rules(dev->iop, imgid, msg);

  return priorities_ok;
}

void dt_iop_transform_image_colorspace(struct dt_iop_module_t *self, float *image, const int width,
                                       const int height, const int cst_from, const int cst_to, int *converted_cst)
{
  if(cst_from == cst_to)
  {
    *converted_cst = cst_to;
    return;
  }

  int transformed = 0;
  const int ch = 4;
  const int im_size = width * height * ch;

  if(cst_from == iop_cs_Lab)
  {
    if(cst_to == iop_cs_linear_rgb)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(image) schedule(static)
#endif
      for(int i = 0; i < im_size; i += ch)
      {
        float XYZ[3];
        dt_Lab_to_XYZ(image + i, XYZ);
        dt_XYZ_to_RGB(XYZ, image + i);
      }
      printf("[dt_iop_transform_image_colorspace] transforming %s %s from lab to linear rgb\n", self->op,
             self->multi_name);
      transformed = 1;
    }
    else if(cst_to == iop_cs_gamma_rgb)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(image) schedule(static)
#endif
      for(int i = 0; i < im_size; i += ch)
      {
        float XYZ[3];
        dt_Lab_to_XYZ(image + i, XYZ);
        dt_XYZ_to_sRGB(XYZ, image + i);
      }
      printf("[dt_iop_transform_image_colorspace] transforming %s %s from lab to gamma rgb\n", self->op,
             self->multi_name);
      transformed = 1;
    }
  }
  else if(cst_from == iop_cs_linear_rgb)
  {
    if(cst_to == iop_cs_Lab)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(image) schedule(static)
#endif
      for(int i = 0; i < im_size; i += ch)
      {
        float XYZ[3];
        dt_RGB_to_XYZ(image + i, XYZ);
        dt_XYZ_to_Lab(XYZ, image + i);
      }
      printf("[dt_iop_transform_image_colorspace] transforming %s %s from linear rgb to lab\n", self->op,
             self->multi_name);
      transformed = 1;
    }
    else if(cst_to == iop_cs_gamma_rgb)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(image) schedule(static)
#endif
      for(int i = 0; i < im_size; i += ch)
      {
        dt_RGB_to_sRGB(image + i, image + i);
      }
      printf("[dt_iop_transform_image_colorspace] transforming %s %s from linear rgb to gamma rgb\n", self->op,
             self->multi_name);
      transformed = 1;
    }
  }
  else if(cst_from == iop_cs_gamma_rgb)
  {
    if(cst_to == iop_cs_Lab)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(image) schedule(static)
#endif
      for(int i = 0; i < im_size; i += ch)
      {
        float XYZ[3];
        dt_sRGB_to_XYZ(image + i, XYZ);
        dt_XYZ_to_Lab(XYZ, image + i);
      }
      printf("[dt_iop_transform_image_colorspace] transforming %s %s from gamma rgb to lab\n", self->op,
             self->multi_name);
      transformed = 1;
    }
    else if(cst_to == iop_cs_linear_rgb)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(image) schedule(static)
#endif
      for(int i = 0; i < im_size; i += ch)
      {
        dt_sRGB_to_RGB(image + i, image + i);
      }
      printf("[dt_iop_transform_image_colorspace] transforming %s %s from gamma rgb to linear rgb\n", self->op,
             self->multi_name);
      transformed = 1;
    }
  }

  if(!transformed)
  {
    fprintf(stderr,
            "[dt_iop_transform_image_colorspace] %s %s invalid colorspace transformation cst_from=%i, cst_to=%i\n",
            self->op, self->multi_name, cst_from, cst_to);

    *converted_cst = cst_from;
  }
  else
    *converted_cst = cst_to;
}

int dt_iop_transform_image_colorspace_cl(struct dt_iop_module_t *self, const int devid, cl_mem dev_img,
                                         const int width, const int height, const int cst_from, const int cst_to,
                                         int *converted_cst)
{
  cl_int err = CL_SUCCESS;
  if(cst_from == cst_to)
  {
    *converted_cst = cst_to;
    return (err == CL_SUCCESS) ? TRUE : FALSE;
  }

  const int ch = 4;
  float *src_buffer = NULL;

  const cl_mem in_colort = dt_opencl_alloc_device_buffer(devid, width * height * ch * sizeof(float));
  if(in_colort == NULL)
  {
    fprintf(stderr, "[dt_iop_transform_image_colorspace_cl] error allocating memory for color transformation\n");
    err = CL_MEM_OBJECT_ALLOCATION_FAILURE;
    goto cleanup;
  }

  src_buffer = dt_alloc_align(64, width * height * ch * sizeof(float));
  if(src_buffer == NULL)
  {
    fprintf(stderr, "[dt_iop_transform_image_colorspace_cl] error allocating memory for color transformation 1\n");
    err = CL_MEM_OBJECT_ALLOCATION_FAILURE;
    goto cleanup;
  }

  {
    size_t origin[] = { 0, 0, 0 };
    size_t region[] = { width, height, 1 };
    err = dt_opencl_enqueue_copy_image_to_buffer(devid, dev_img, in_colort, origin, region, 0);
    if(err != CL_SUCCESS)
    {
      fprintf(stderr,
              "[dt_iop_transform_image_colorspace_cl] error allocating memory for color transformation 2\n");
      goto cleanup;
    }
  }

  err = dt_opencl_read_buffer_from_device(devid, (void *)src_buffer, in_colort, 0,
                                          (size_t)width * height * ch * sizeof(float), CL_TRUE);
  if(err != CL_SUCCESS)
  {
    fprintf(stderr, "[dt_iop_transform_image_colorspace_cl] error allocating memory for color transformation 3\n");
    goto cleanup;
  }

  // just call the CPU version for now
  dt_iop_transform_image_colorspace(self, src_buffer, width, height, cst_from, cst_to, converted_cst);

  err = dt_opencl_write_buffer_to_device(devid, src_buffer, in_colort, 0, width * height * ch * sizeof(float),
                                         TRUE);
  if(err != CL_SUCCESS)
  {
    fprintf(stderr, "[dt_iop_transform_image_colorspace_cl] error allocating memory for color transformation 4\n");
    goto cleanup;
  }

  {
    size_t origin[] = { 0, 0, 0 };
    size_t region[] = { width, height, 1 };
    err = dt_opencl_enqueue_copy_buffer_to_image(devid, in_colort, dev_img, 0, origin, region);
    if(err != CL_SUCCESS)
    {
      fprintf(stderr,
              "[dt_iop_transform_image_colorspace_cl] error allocating memory for color transformation 5\n");
      goto cleanup;
    }
  }

cleanup:
  if(src_buffer) dt_free_align(src_buffer);
  if(in_colort) dt_opencl_release_mem_object(in_colort);

  return (err == CL_SUCCESS) ? TRUE : FALSE;
}

