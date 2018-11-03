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

GList *dt_get_iop_priorities_v0()
{
  GList *priorities = NULL;

  const dt_iop_priority_entry_v0_t prior_entry[]
      = { /*{ 0, 0.f, "rawspeed", 1 },*/
          { 0, 0.f, "rawprepare", 1 },     { 0, 0.f, "invert", 0 },
          { 0, 0.f, "temperature", 1 },    { 0, 0.f, "highlights", 1 },
          { 0, 0.f, "cacorrect", 0 },      { 0, 0.f, "hotpixels", 0 },
          { 0, 0.f, "rawdenoise", 0 },     { 0, 0.f, "demosaic", 1 },
          { 0, 0.f, "denoiseprofile", 0 }, { 0, 0.f, "tonemap", 0 },
          { 0, 0.f, "exposure", 0 },       { 0, 0.f, "spots", 0 },
          { 0, 0.f, "retouch", 0 },        { 0, 0.f, "lens", 0 },
          { 0, 0.f, "ashift", 0 },         { 0, 0.f, "liquify", 0 },
          { 0, 0.f, "rotatepixels", 0 },   { 0, 0.f, "scalepixels", 0 },
          { 0, 0.f, "flip", 1 },           { 0, 0.f, "graduatednd", 0 },
          { 0, 0.f, "basecurve", 0 },      { 0, 0.f, "bilateral", 0 },
          { 0, 0.f, "profile_gamma", 0 },  { 0, 0.f, "hazeremoval", 0 },
          { 0, 0.f, "colorin", 1 },        { 0, 0.f, "colorreconstruct", 0 },
          { 0, 0.f, "colorchecker", 0 },   { 0, 0.f, "defringe", 0 },
          { 0, 0.f, "equalizer", 0 },      { 0, 0.f, "vibrance", 0 },
          { 0, 0.f, "colorbalance", 0 },   { 0, 0.f, "clipping", 0 },
          { 0, 0.f, "colorize", 0 },       { 0, 0.f, "colortransfer", 0 },
          { 0, 0.f, "colormapping", 0 },   { 0, 0.f, "bloom", 0 },
          { 0, 0.f, "nlmeans", 0 },        { 0, 0.f, "globaltonemap", 0 },
          { 0, 0.f, "shadhi", 0 },         { 0, 0.f, "atrous", 0 },
          { 0, 0.f, "bilat", 0 },          { 0, 0.f, "colorzones", 0 },
          { 0, 0.f, "lowlight", 0 },       { 0, 0.f, "monochrome", 0 },
          { 0, 0.f, "colisa", 0 },         { 0, 0.f, "zonesystem", 0 },
          { 0, 0.f, "tonecurve", 0 },      { 0, 0.f, "levels", 0 },
          { 0, 0.f, "relight", 0 },        { 0, 0.f, "colorcorrection", 0 },
          { 0, 0.f, "sharpen", 0 },        { 0, 0.f, "lowpass", 0 },
          { 0, 0.f, "highpass", 0 },       { 0, 0.f, "grain", 0 },
          { 0, 0.f, "colorcontrast", 0 },  { 0, 0.f, "colorout", 1 },
          { 0, 0.f, "channelmixer", 0 },   { 0, 0.f, "soften", 0 },
          { 0, 0.f, "vignette", 0 },       { 0, 0.f, "splittoning", 0 },
          { 0, 0.f, "velvia", 0 },         { 0, 0.f, "clahe", 0 },
          { 0, 0.f, "finalscale", 1 },     { 0, 0.f, "overexposed", 1 },
          { 0, 0.f, "rawoverexposed", 1 }, { 0, 0.f, "borders", 0 },
          { 0, 0.f, "watermark", 0 },      { 0, 0.f, "dither", 0 },
          { 0, 0.f, "gamma", 1 },          { -1, 0.f, "\0", 0 }
        };

  int i = 0;
  while(prior_entry[i].priority >= 0)
  {
    dt_iop_priority_entry_v0_t *prior = calloc(1, sizeof(dt_iop_priority_entry_v0_t));

    prior->priority = i + 1;
    prior->iop_order = (float)i;
    snprintf(prior->operation, sizeof(prior->operation), "%s", prior_entry[i].operation);
    prior->default_enabled = prior_entry[i].default_enabled;

    priorities = g_list_append(priorities, prior);
    i++;
  }

  return priorities;
}

GList *_get_iop_priorities()
{
  GList *priorities = NULL;

  const dt_iop_priority_entry_t prior_entry[] = { /*{ 0, 0.f, "rawspeed" },*/
                                                  { 0, 0.f, "rawprepare" },
                                                  { 0, 0.f, "invert" },
                                                  { 0, 0.f, "temperature" },
                                                  { 0, 0.f, "highlights" },
                                                  { 0, 0.f, "cacorrect" },
                                                  { 0, 0.f, "hotpixels" },
                                                  { 0, 0.f, "rawdenoise" },
                                                  { 0, 0.f, "demosaic" },
                                                  { 0, 0.f, "dummy1" },
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
                                                  { 0, 0.f, "graduatednd" },
                                                  { 0, 0.f, "basecurve" },
                                                  { 0, 0.f, "bilateral" },
                                                  { 0, 0.f, "profile_gamma" },
                                                  { 0, 0.f, "hazeremoval" },
                                                  { 0, 0.f, "dummy2" },
                                                  { 0, 0.f, "colorin" },
                                                  { 0, 0.f, "dummy3" },
                                                  { 0, 0.f, "colorreconstruct" },
                                                  { 0, 0.f, "colorchecker" },
                                                  { 0, 0.f, "defringe" },
                                                  { 0, 0.f, "equalizer" },
                                                  { 0, 0.f, "vibrance" },
                                                  { 0, 0.f, "colorbalance" },
                                                  { 0, 0.f, "clipping" },
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
                                                  { 0, 0.f, "dummy4" },
                                                  { 0, 0.f, "colorout" },
                                                  { 0, 0.f, "dummy5" },
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
                                                  { 0, 0.f, "dummy6" },
                                                  { 0, 0.f, "gamma" },
                                                  { -1, 0.f, "\0" }
  };

  int i = 0;
  while(prior_entry[i].priority >= 0)
  {
    dt_iop_priority_entry_t *prior = calloc(1, sizeof(dt_iop_priority_entry_t));

    prior->priority = i + 1;
    prior->iop_order = (float)i;
    snprintf(prior->operation, sizeof(prior->operation), "%s", prior_entry[i].operation);

    priorities = g_list_append(priorities, prior);
    i++;
  }

  return priorities;
}

int dt_iop_priorities_get_default_output_cst_v0(const char *op_name)
{
  int cst = 0;
  int priority_demosaic = -1;
  int priority_colorin = -1;
  int priority_colorout = -1;
  int priority_gamma = -1;
  int priority_module = -1;

  GList *prior_vo = dt_get_iop_priorities_v0();

  GList *priorities = g_list_first(prior_vo);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)priorities->data;

    if(strcmp(prior->operation, op_name) == 0) priority_module = prior->priority;
    if(strcmp(prior->operation, "demosaic") == 0) priority_demosaic = prior->priority;
    if(strcmp(prior->operation, "colorin") == 0) priority_colorin = prior->priority;
    if(strcmp(prior->operation, "colorout") == 0) priority_colorout = prior->priority;
    if(strcmp(prior->operation, "gamma") == 0) priority_gamma = prior->priority;

    priorities = g_list_next(priorities);
  }

  g_list_free_full(prior_vo, free);

  if(priority_module < 0.f)
    fprintf(stderr, "[dt_iop_priorities_get_default_output_cst_v0] can't find default output cst for module %s\n",
            op_name);
  else if(priority_module < priority_demosaic)
    cst = iop_cs_RAW;
  else if(priority_module < priority_colorin)
    cst = iop_cs_linear_rgb;
  else if(priority_module < priority_colorout)
    cst = iop_cs_Lab;
  else if(priority_module <= priority_gamma)
    cst = iop_cs_gamma_rgb;
  else
    fprintf(stderr, "[dt_iop_priorities_get_default_output_cst_v0] can't find default output cst for module %s\n",
            op_name);

  return cst;
}

int dt_iop_priorities_get_default_input_colorspace(const char *op_name)
{
  // printf("[dt_iop_priorities_get_default_input_colorspace] begin\n");
  int cst = 0;
  int priority_demosaic = -1;
  int priority_colorin = -1;
  int priority_colorout = -1;
  int priority_gamma = -1;
  int priority_module = -1;

  GList *priorities = g_list_first(darktable.iop_priorities);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)(priorities->data);

    if(strcmp(prior->operation, op_name) == 0) priority_module = prior->priority;
    if(strcmp(prior->operation, "demosaic") == 0) priority_demosaic = prior->priority;
    if(strcmp(prior->operation, "colorin") == 0) priority_colorin = prior->priority;
    if(strcmp(prior->operation, "colorout") == 0) priority_colorout = prior->priority;
    if(strcmp(prior->operation, "gamma") == 0) priority_gamma = prior->priority;

    priorities = g_list_next(priorities);
  }

  if(priority_module < 0.f)
    fprintf(stderr,
            "[dt_iop_priorities_get_default_input_colorspace] can't find default input cst for module %s\n",
            op_name);
  else if(priority_module <= priority_demosaic)
    cst = iop_cs_RAW;
  else if(priority_module <= priority_colorin)
    cst = iop_cs_linear_rgb;
  else if(priority_module <= priority_colorout)
    cst = iop_cs_Lab;
  else if(priority_module <= priority_gamma)
    cst = iop_cs_gamma_rgb;
  else
    fprintf(stderr,
            "[dt_iop_priorities_get_default_input_colorspace] can't find default input cst for module %s\n",
            op_name);

  return cst;
}

int dt_iop_priorities_get_default_output_colorspace(const char *op_name)
{
  int cst = 0;
  int priority_demosaic = -1;
  int priority_colorin = -1;
  int priority_colorout = -1;
  int priority_gamma = -1;
  int priority_module = -1;

  GList *priorities = g_list_first(darktable.iop_priorities);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)(priorities->data);

    if(strcmp(prior->operation, op_name) == 0) priority_module = prior->priority;
    if(strcmp(prior->operation, "demosaic") == 0) priority_demosaic = prior->priority;
    if(strcmp(prior->operation, "colorin") == 0) priority_colorin = prior->priority;
    if(strcmp(prior->operation, "colorout") == 0) priority_colorout = prior->priority;
    if(strcmp(prior->operation, "gamma") == 0) priority_gamma = prior->priority;

    priorities = g_list_next(priorities);
  }

  if(priority_module < 0.f)
    fprintf(stderr,
            "[dt_iop_priorities_get_default_output_colorspace] can't find default output cst for module %s(%i)\n",
            op_name, priority_module);
  else if(priority_module < priority_demosaic)
    cst = iop_cs_RAW;
  else if(priority_module < priority_colorin)
    cst = iop_cs_linear_rgb;
  else if(priority_module < priority_colorout)
    cst = iop_cs_Lab;
  else if(priority_module <= priority_gamma)
    cst = iop_cs_gamma_rgb;
  else
    fprintf(stderr,
            "[dt_iop_priorities_get_default_output_colorspace] can't find default output cst for module %s(%i)\n",
            op_name, priority_module);

  return cst;
}

gint dt_sort_iop_so_by_priority(gconstpointer a, gconstpointer b)
{
  const dt_iop_module_so_t *am = (const dt_iop_module_so_t *)a;
  const dt_iop_module_so_t *bm = (const dt_iop_module_so_t *)b;
  return (am->priority - bm->priority);
}

gint dt_sort_iop_by_order(gconstpointer a, gconstpointer b)
{
  const dt_iop_module_t *am = (const dt_iop_module_t *)a;
  const dt_iop_module_t *bm = (const dt_iop_module_t *)b;
  if(am->iop_order > bm->iop_order) return 1;
  if(am->iop_order < bm->iop_order) return -1;
  return (am->priority - bm->priority);
}

// loads config file iop_priorities and set the priority in darktable.iop list
// iop_priorities is a list of operation names in the order that must be executed
// by default on the pipe
GList *_load_iop_priorities_from_file()
{
  GList *priorities = NULL;
  FILE *f = 0;
  char datadir[PATH_MAX] = { 0 };
  dt_loc_get_user_config_dir(datadir, sizeof(datadir));
  char iop_priorities[PATH_MAX] = { 0 };
  snprintf(iop_priorities, sizeof(iop_priorities), "%s/iop_priorities", datadir);

#define LINE_SIZE 1023

  char line[LINE_SIZE + 1];

  // read priorities from file and update darktable.iop list
  f = g_fopen(iop_priorities, "rb");
  if(f)
  {
    int read = 0;
    int i = 1;
    while(!feof(f))
    {
      read = fscanf(f, "%s", line);
      if(read > 0 && line[0] != '#' && line[0] != ' ' && line[0] != '\0')
      {
        dt_iop_priority_entry_t *prior = calloc(1, sizeof(dt_iop_priority_entry_t));

        prior->priority = i;
        prior->iop_order = (float)i;
        snprintf(prior->operation, sizeof(prior->operation), "%s", line);
        // prior->default_enabled = prior->default_enabled;

        priorities = g_list_append(priorities, prior);
        i++;
        /*        int module_found = 0;
                GList *modules = g_list_first(darktable.iop);
                while(modules)
                {
                  dt_iop_module_so_t *mod = (dt_iop_module_so_t *)(modules->data);
                  if(strcmp(mod->op, line) == 0)
                  {
                    mod->priority = i++;
                    module_found = 1;
                    break;
                  }
                  modules = g_list_next(modules);
                }
                if(!module_found) fprintf(stderr, "[dt_iop_priorities_set_so_iop_priorities] module %s is not
           installed!!!\n", line);*/
      }
    }
  }

  if(f) fclose(f);

#undef LINE_SIZE

  return priorities;
}

#if 0
// loads config file iop_priorities and set the priority in darktable.iop list
// iop_priorities is a list of operation names in the order that must be executed
// by default on the pipe
void dt_iop_priorities_set_so_iop_priorities()
{
  printf("[dt_iop_priorities_set_so_iop_priorities] begin\n");
  GList *loaded_priorities = _load_iop_priorities_from_file();
  if(loaded_priorities == NULL) loaded_priorities = _get_iop_priorities();

  GList *priorities = g_list_first(loaded_priorities);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)(priorities->data);

    int module_found = 0;
    GList *modules = g_list_first(darktable.iop);
    while(modules)
    {
      dt_iop_module_so_t *mod = (dt_iop_module_so_t *)(modules->data);
      if(strcmp(mod->op, prior->operation) == 0)
      {
        mod->priority = prior->priority;
        module_found = 1;
        break;
      }
      modules = g_list_next(modules);
    }
    if(!module_found)
      fprintf(stderr, "[dt_iop_priorities_set_so_iop_priorities] module %s is not installed!!!\n",
              prior->operation);

    priorities = g_list_next(priorities);
  }

  /*  FILE *f = 0;
    char datadir[PATH_MAX] = { 0 };
    dt_loc_get_user_config_dir(datadir, sizeof(datadir));
    char iop_priorities[PATH_MAX] = { 0 };
    snprintf(iop_priorities, sizeof(iop_priorities), "%s/iop_priorities", datadir);

  #define LINE_SIZE 1023

    char line[LINE_SIZE + 1];

    // read priorities from file and update darktable.iop list
    f = g_fopen(iop_priorities, "rb");
    if(f)
    {
      int read = 0;
      int i = 1;
      while(!feof(f))
      {
        read = fscanf(f, "%s", line);
        if(read > 0 && line[0] != '#' && line[0] != ' ' && line[0] != '\0')
        {
          int module_found = 0;
          GList *modules = g_list_first(darktable.iop);
          while(modules)
          {
            dt_iop_module_so_t *mod = (dt_iop_module_so_t *)(modules->data);
            if(strcmp(mod->op, line) == 0)
            {
              mod->priority = i++;
              module_found = 1;
              break;
            }
            modules = g_list_next(modules);
          }
          if(!module_found) fprintf(stderr, "[dt_iop_priorities_set_so_iop_priorities] module %s is not
  installed!!!\n", line);
        }
      }
    }

    if(f) fclose(f);
  */
  // and sort the iop list
  darktable.iop = g_list_sort(darktable.iop, dt_sort_iop_so_by_priority);

  // now check if all the modules have their priority assigned
  GList *modules = g_list_first(darktable.iop);
  while(modules)
  {
    dt_iop_module_so_t *mod = (dt_iop_module_so_t *)(modules->data);
    if(mod->priority < 0)
    {
      fprintf(stderr, "[dt_iop_priorities_set_so_iop_priorities] missing priority for module %s\n", mod->op);
    }
    modules = g_list_next(modules);
  }

  if(loaded_priorities) g_list_free_full(loaded_priorities, free);
  /*
  #undef LINE_SIZE*/
  printf("[dt_iop_priorities_set_so_iop_priorities] end\n");
}
#endif

int dt_iop_priorities_check_iop_priorities()
{
  int priority_missing = 0;

  // check if all the modules have their priority assigned
  GList *modules = g_list_first(darktable.iop);
  while(modules)
  {
    dt_iop_module_so_t *mod = (dt_iop_module_so_t *)(modules->data);
    if(mod->priority < 0)
    {
      priority_missing = 1;
      fprintf(stderr, "[dt_iop_priorities_check_iop_priorities] missing priority for module %s\n", mod->op);
    }
    modules = g_list_next(modules);
  }

  return priority_missing;
}

// loads config file iop_priority_rules into a list and returns it
// iop_priority_rules is a list of pairs of operation names
// that indicates that the fist one must always come before in the pipe than the second one
// so colorin colorout tell that colorin always come before colorout in the pipe
// any module can be moved arround those two
GList *_load_iop_priority_rules_from_file()
{
  GList *rules = NULL;
  FILE *f = 0;
  char datadir[PATH_MAX] = { 0 };
  dt_loc_get_user_config_dir(datadir, sizeof(datadir));
  char iop_priority_rules[PATH_MAX] = { 0 };
  snprintf(iop_priority_rules, sizeof(iop_priority_rules), "%s/iop_priority_rules", datadir);

#define LINE_SIZE 1023

  char line1[LINE_SIZE + 1];
  char line2[LINE_SIZE + 1];

  // read priority rules from file and add it to a list
  f = g_fopen(iop_priority_rules, "rb");
  if(f)
  {
    int read = 0;
    while(!feof(f))
    {
      read = fscanf(f, "%s %s", line1, line2);
      if(read > 1 && line1[0] != '#' && line1[0] != ' ' && line1[0] != '\0')
      {
        dt_iop_priority_rule_t *rule = calloc(1, sizeof(dt_iop_priority_rule_t));

        snprintf(rule->op_prev, sizeof(rule->op_prev), "%s", line1);
        snprintf(rule->op_next, sizeof(rule->op_next), "%s", line2);
        // printf("[_load_iop_priority_rules_from_file] %s %s\n", rule->op_prev, rule->op_next);
        rules = g_list_append(rules, rule);
      }
    }
  }

  if(f) fclose(f);

#undef LINE_SIZE

  return rules;
}

GList *_get_iop_priority_rules()
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

GList *dt_iop_priorities_load_iop_priorities()
{
  GList *priorities = _load_iop_priorities_from_file();
  if(priorities == NULL) priorities = _get_iop_priorities();

  /*  GList *rules_list = g_list_first(rules);
    while(rules_list)
    {
      dt_iop_priority_rule_t *rule = (dt_iop_priority_rule_t *)rules_list->data;

      printf("[dt_iop_priorities_load_iop_priority_rules] %s %s\n", rule->op_prev, rule->op_next);

      rules_list = g_list_next(rules_list);
    }
  */
  return priorities;
}

GList *dt_iop_priorities_load_iop_priority_rules()
{
  GList *rules = _load_iop_priority_rules_from_file();
  if(rules == NULL) rules = _get_iop_priority_rules();

  /*  GList *rules_list = g_list_first(rules);
    while(rules_list)
    {
      dt_iop_priority_rule_t *rule = (dt_iop_priority_rule_t *)rules_list->data;

      printf("[dt_iop_priorities_load_iop_priority_rules] %s %s\n", rule->op_prev, rule->op_next);

      rules_list = g_list_next(rules_list);
    }
  */
  return rules;
}

int dt_iop_priorities_get_iop_priority(const char *op_name)
{
  int priority = -1;

  GList *priorities = g_list_first(darktable.iop_priorities);
  while(priorities)
  {
    dt_iop_priority_entry_t *prior = (dt_iop_priority_entry_t *)(priorities->data);

    if(strcmp(prior->operation, op_name) == 0)
    {
      priority = prior->priority;
      break;
    }

    priorities = g_list_next(priorities);
  }

  if(priority < 0)
    fprintf(stderr, "[dt_iop_priorities_get_iop_priority] can't find priority for module %s\n", op_name);

  return priority;
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
float dt_get_iop_order_before_iop(dt_develop_t *dev, dt_iop_module_t *module, dt_iop_module_t *module_next,
                                  const int validate_order, const int log_error)
{
  // dt_develop_t *dev = module->dev;
  float iop_order = -1.f;
  int module_found = 0;
  int rule_found = 0;

  // this can happen when migrating from an old version
  // let's handle it here
  if(module->iop_order == module_next->iop_order)
  {
    // we need to find the module before this two
    GList *modules = g_list_first(dev->iop);
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
               "[dt_get_iop_order_before_iop] 1-calculated new iop_order=%f "
            "for %s(%f) between %s(%f) and %s(%f)\n",
               iop_order, module->op, module->iop_order, module_next->op, module_next->iop_order, prev->op,
               prev->iop_order); */
      }
      else
      {
        // those are the first modules, this can't happen
        if(log_error)
          fprintf(stderr, "[dt_get_iop_order_before_iop] %s(%f) and %s(%f) are the first modules!!!\n", module->op,
                  module->iop_order, module_next->op, module_next->iop_order);
      }
    }
  }
  // module is before on the pipe
  // move it up
  else if(module->iop_order < module_next->iop_order)
  {
    GList *modules = g_list_first(dev->iop);
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
              fprintf(stderr, "[dt_get_iop_order_before_iop] %s(%f) is already previous to %s(%f)\n", module->op,
                      module->iop_order, prev->op, prev->iop_order);
            break;
          }

          // calculate new iop_order
          iop_order = module_next->iop_order - (module_next->iop_order - prev->iop_order) / 2.f;
          /* fprintf(stderr,
                 "[dt_get_iop_order_before_iop] 2-calculated new iop_order=%f "
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
              fprintf(stderr, "[dt_get_iop_order_before_iop] can't move %s(%f) pass %s(%f)\n", module->op,
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
                        "[dt_get_iop_order_before_iop] found rule %s %s while moving %s(%f) before %s(%f)\n",
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
    GList *modules = g_list_last(dev->iop);
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
                 "[dt_get_iop_order_before_iop] 3-calculated new iop_order=%f "
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
              fprintf(stderr, "[dt_get_iop_order_before_iop] can't move %s(%f) pass %s(%f)\n", module->op,
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
                        "[dt_get_iop_order_before_iop] found rule %s %s while moving %s(%f) before %s(%f)\n",
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
        fprintf(stderr, "[dt_get_iop_order_before_iop] can't find next module %s(%f) while moving %s(%f)\n",
                module_next->op, module_next->iop_order, module->op, module->iop_order);
  }

  if(!module_found) fprintf(stderr, "[dt_get_iop_order_before_iop] can't find module %s\n", module->op);

  return iop_order;
}

// if module can be placed after than module_prev on the pipe
// it returns the new iop_order
// if it cannot be placed it returns -1.f
// this assums that the order is always positive
float dt_get_iop_order_after_iop(dt_develop_t *dev, dt_iop_module_t *module, dt_iop_module_t *module_prev,
                                 const int validate_order, const int log_error)
{
  // dt_develop_t *dev = module->dev;
  float iop_order = -1.f;

  // moving after module_prev is the same as moving before the very next one after module_prev
  GList *modules = g_list_last(dev->iop);
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
        "[dt_get_iop_order_after_iop] can't find module previous to %s %s(%f) while moving %s %s(%f) after it\n",
        module_prev->op, module_prev->multi_name, module_prev->iop_order, module->op, module->multi_name,
        module->iop_order);
  }
  else
    iop_order = dt_get_iop_order_before_iop(dev, module, module_next, validate_order, log_error);

  return iop_order;
}

// changes the module->iop_order so it comes before in the pipe than module_next
// sort dev->iop to reflect the changes
// return 1 if iop_order is changed, 0 otherwise
int dt_move_iop_before(dt_develop_t *dev, dt_iop_module_t *module, dt_iop_module_t *module_next,
                       const int validate_order, const int log_error)
{
  // dt_develop_t *dev = module->dev;
  int moved = 0;

  // dt_iop_priorities_check_priorities(dev, "dt_move_iop_before begin");
  // printf("[dt_move_iop_before] moving %s(%f) before %s(%f)\n", module->op, module->iop_order, module_next->op,
  // module_next->iop_order);

  const float iop_order = dt_get_iop_order_before_iop(dev, module, module_next, validate_order, log_error);
  // printf("[dt_move_iop_before] new iop_order for %s is %f\n", module->op, iop_order);

  if(iop_order >= 0.f)
  {
    module->iop_order = iop_order;
    dev->iop = g_list_sort(dev->iop, dt_sort_iop_by_order);
    moved = 1;
  }
  else if(log_error)
    fprintf(stderr, "[dt_move_iop_before] module %s is already before %s\n", module->op, module_next->op);

  // dt_iop_priorities_check_priorities(dev, "dt_move_iop_before end");

  return moved;
}

// changes the module->iop_order so it comes after in the pipe than module_prev
// sort dev->iop to reflect the changes
// return 1 if iop_order is changed, 0 otherwise
int dt_move_iop_after(dt_develop_t *dev, dt_iop_module_t *module, dt_iop_module_t *module_prev,
                      const int validate_order, const int log_error)
{
  // printf("[dt_move_iop_after] moving %s(%f) after %s(%f)\n", module->op, module->iop_order, module_prev->op,
  // module_prev->iop_order);
  // dt_develop_t *dev = module->dev;
  int moved = 0;

  // dt_iop_priorities_check_priorities(dev, "dt_move_iop_after begin");

  /*  // moving after module_prev is the same as moving before the very next one after module_prev
    GList *modules = g_list_last(dev->iop);
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
      fprintf(stderr, "[dt_move_iop_after] can't find module previous to %s while moving %s\n", module_prev->op,
              module->op);
      return 0;
    }

    const float iop_order = dt_get_iop_order_before_iop(module, module_next, validate_order);*/
  const float iop_order = dt_get_iop_order_after_iop(dev, module, module_prev, validate_order, log_error);
  // printf("[dt_move_iop_after] new iop_order for %s is %f\n", module->op, iop_order);
  if(iop_order >= 0.f)
  {
    module->iop_order = iop_order;
    dev->iop = g_list_sort(dev->iop, dt_sort_iop_by_order);
    moved = 1;
  }
  else if(log_error)
    fprintf(stderr, "[dt_move_iop_after] module %s is already after %s\n", module->op, module_prev->op);

  // dt_iop_priorities_check_priorities(dev, "dt_move_iop_after end");

  return moved;
}

float get_iop_default_order(const char *op_name)
{
  float iop_order = -1.f;

  GList *modules = g_list_first(darktable.iop);
  while(modules)
  {
    dt_iop_module_so_t *mod = (dt_iop_module_so_t *)(modules->data);
    if(strcmp(mod->op, op_name) == 0)
    {
      iop_order = (float)mod->priority;
      break;
    }
    modules = g_list_next(modules);
  }

  return iop_order;
}
/*
// updates iop_order on history for a given module
static void _iop_priorities_update_history_iop_order(dt_develop_t *dev, dt_iop_module_t *module, const float
iop_order)
{
  GList *history = g_list_first(dev->history);
  while(history)
  {
    dt_dev_history_item_t *hist = (dt_dev_history_item_t *)(history->data);
    if(strcmp(hist->op_name, module->op) == 0 && hist->multi_priority == module->multi_priority)
    {
      hist->iop_order = iop_order;
    }
    history = g_list_next(history);
  }
}

// rebuilds the iop_order on the history and dev->iop
// this assume that all entries in history for a given instance
// have the same iop_order, so it should be called after compress history
void dt_iop_priorities_rebuild_order(dt_develop_t *dev)
{
  // we want, when possible, that iop_order is the defauls order, so iqual to priority
  // that's only possible for the first instance of a multi-instance module
  // the rest of the multi-instances will be evenly distributed

  float last_iop_order = 0.f;
  float increment = 1.f;

  // go through all modules in dev and darktable
  // they are in iop_order order, so will read it in synch
  GList *dev_modules = g_list_first(dev->iop);
  GList *darktable_modules = g_list_first(darktable.iop);
  while(dev_modules && darktable_modules)
  {
    dt_iop_module_t *dev_mod = (dt_iop_module_t *)dev_modules->data;
    dt_iop_module_so_t *dt_mod = (dt_iop_module_so_t *)(darktable_modules->data);

    // this is the first apperance of the module on the pipe
    // we can use the default iop_order
    if(dev_mod->priority == dt_mod->priority)
    {
      printf("[dt_iop_priorities_rebuild_order] 1st apperance dev %s(%i) dt %s(%i)\n", dev_mod->op,
dev_mod->priority, dt_mod->op, dt_mod->priority);
      last_iop_order = (float)dt_mod->priority;
      if(dev_mod->iop_order != last_iop_order)
      {
        _iop_priorities_update_history_iop_order(dev, dev_mod, last_iop_order);
        dev_mod->iop_order = last_iop_order;
      }
      // reset the increment for the next out-of-order module
      increment = 1.f;

      // move to the next module on both lists
      dev_modules = g_list_next(dev_modules);
      darktable_modules = g_list_next(darktable_modules);
    }
    else
    {
      if((dev_mod->flags() & (IOP_FLAGS_HIDDEN | IOP_FLAGS_ONE_INSTANCE | IOP_FLAGS_NO_HISTORY_STACK)) ==
(IOP_FLAGS_HIDDEN | IOP_FLAGS_ONE_INSTANCE | IOP_FLAGS_NO_HISTORY_STACK))
      {
        printf("[dt_iop_priorities_rebuild_order] IOP_FLAGS_HIDDEN dev %s(%i) dt %s(%i)\n", dev_mod->op,
dev_mod->priority, dt_mod->op, dt_mod->priority);
        // this modules still can have the default iop_order
        // but there's no need to update history
        last_iop_order = (float)dt_mod->priority;
        if(dev_mod->iop_order != last_iop_order)
          dev_mod->iop_order = last_iop_order;
        // reset the increment for the next out-of-order module
        increment = 1.f;

        // just need to move dev list
        dev_modules = g_list_next(dev_modules);
      }
      else
      {
        printf("[dt_iop_priorities_rebuild_order] out-of-order dev %s(%i) dt %s(%i)\n", dev_mod->op,
dev_mod->priority, dt_mod->op, dt_mod->priority);
        // this module is out-of-order
        // it means that it was moved or is a multi-instance
        // we will increment by .5 / number of modules until we reach a module that
        // is in order
        increment *= .5f;
        last_iop_order += increment;
        if(dev_mod->iop_order != last_iop_order)
        {
          _iop_priorities_update_history_iop_order(dev, dev_mod, last_iop_order);
          dev_mod->iop_order = last_iop_order;
        }
      }
    }


  }

  if(dev_modules==NULL && darktable_modules!=NULL)
  fprintf(stderr, "[dt_iop_priorities_rebuild_order] there are more modules on dev than in darktable!\n");
  if(dev_modules!=NULL && darktable_modules==NULL)
  fprintf(stderr, "[dt_iop_priorities_rebuild_order] there are more modules on darktable than in dev!\n");
}
*/

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

static dt_iop_pipe_entry_t *_search_pipe_by_op_multi_priority(GList *main_pipe_list, dt_iop_module_t *module)
{
  dt_iop_pipe_entry_t *ret_pi = NULL;
  GList *pipe = g_list_first(main_pipe_list);
  while(pipe)
  {
    dt_iop_pipe_entry_t *pi = (dt_iop_pipe_entry_t *)pipe->data;

    if(strcmp(module->op, pi->operation) == 0 && module->multi_priority == pi->multi_priority)
    {
      ret_pi = pi;
      break;
    }

    pipe = g_list_next(pipe);
  }
  return ret_pi;
}

static GList *_get_fence_modules_list(dt_develop_t *dev)
{
  GList *fences = NULL;
  GList *modules = g_list_first(dev->iop);
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

static void _check_rules(dt_develop_t *dev, const char *msg)
{
  GList *modules = NULL;

  // check for IOP_FLAGS_FENCE on each module
  // create a list of fences modules
  GList *fences = _get_fence_modules_list(dev);

  // check if each module is between the fences
  modules = g_list_first(dev->iop);
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
      fprintf(stderr, "[_check_rules] found fence %s %s module %s %s(%f) is after %s %s(%f) (%s)\n",
              fence_next->op, fence_next->multi_name, mod->op, mod->multi_name, mod->iop_order, fence_next->op,
              fence_next->multi_name, fence_next->iop_order, msg);
    }
    if(fence_prev && mod->iop_order < fence_prev->iop_order)
    {
      fprintf(stderr, "[_check_rules] found fence %s %s module %s %s(%f) is before %s %s(%f) (%s)\n",
              fence_prev->op, fence_prev->multi_name, mod->op, mod->multi_name, mod->iop_order, fence_prev->op,
              fence_prev->multi_name, fence_prev->iop_order, msg);
    }


    modules = g_list_next(modules);
  }

  // for each module check if it doesn't break a rule
  modules = g_list_first(dev->iop);
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
            fprintf(stderr, "[_check_rules] found rule %s %s module %s %s(%f) is after %s %s(%f) (%s)\n",
                    rule->op_prev, rule->op_next, mod->op, mod->multi_name, mod->iop_order, mod_prev->op,
                    mod_prev->multi_name, mod_prev->iop_order, msg);
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
            fprintf(stderr, "[_check_rules] found rule %s %s module %s %s(%f) is before %s %s(%f) (%s)\n",
                    rule->op_prev, rule->op_next, mod->op, mod->multi_name, mod->iop_order, mod_next->op,
                    mod_next->multi_name, mod_next->iop_order, msg);
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

#if 0
void dt_iop_priorities_read_pipe(dt_develop_t *dev, const int imgid)
{
  // printf("dt_iop_priorities_read_pipe begin\n");
  dt_iop_priorities_check_priorities(dev, "dt_iop_priorities_read_pipe begin");

  sqlite3_stmt *stmt;

  GList *modules = NULL;

  const int log_error = 1;

  int demosaic_priority = -1;
  int colorin_priority = -1;
  int colorout_priority = -1;

  float demosaic_iop_order = -1.f;
  float colorin_iop_order = -1.f;
  float colorout_iop_order = -1.f;
  float gamma_iop_order = -1.f;

  dt_iop_module_t *demosaic_module = NULL;
  dt_iop_module_t *colorin_module = NULL;
  dt_iop_module_t *colorout_module = NULL;
  dt_iop_module_t *gamma_module = NULL;

  int demosaic_priority_old = -1;
  int colorin_priority_old = -1;
  int colorout_priority_old = -1;
  int gamma_priority_old = -1;

  GList *main_pipe_list = NULL;

  // read the last state of the pipe and update the iop_order on dev->iop
  // this will ensure that the pipe order is correct even for modules that are not in the history
  // main.pipe has only entries that are not in history but enabled anyway
  DT_DEBUG_SQLITE3_PREPARE_V2(dt_database_get(darktable.db),
                              "SELECT operation, priority, multi_priority, iop_order "
                              "FROM main.pipe WHERE imgid = ?1",
                              -1, &stmt, NULL);
  DT_DEBUG_SQLITE3_BIND_INT(stmt, 1, imgid);
  // dev->history_end = 0;
  int order_changed = 0;
  while(sqlite3_step(stmt) == SQLITE_ROW)
  {
    const char *opname = (const char *)sqlite3_column_text(stmt, 0);
    const int priority = sqlite3_column_int(stmt, 1);
    const int multi_priority = sqlite3_column_int(stmt, 2);
    const float iop_order = sqlite3_column_double(stmt, 3);

    // printf("[dt_iop_priorities_read_pipe] read module %s\n", opname);

    modules = g_list_first(dev->iop);
    while(modules)
    {
      dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
      if(strcmp(mod->op, opname) == 0 && mod->multi_priority == multi_priority)
      {
        if(mod->iop_order != iop_order)
        {
          /* fprintf(stderr,
                  "[dt_iop_priorities_read_pipe] order changed for %s %s (%i) pipe %f current %f image %i\n",
                  mod->op, mod->multi_name, mod->multi_priority, iop_order, mod->iop_order, imgid); */
          order_changed = 1;
          mod->iop_order = iop_order;
        }

        // save it to the list, it will be needed later
        dt_iop_pipe_entry_t *former = calloc(1, sizeof(dt_iop_pipe_entry_t));
        former->priority = priority;
        former->multi_priority = multi_priority;
        snprintf(former->operation, sizeof(former->operation), "%s", opname);

        main_pipe_list = g_list_append(main_pipe_list, former);

        break;
      }
      modules = g_list_next(modules);
    }

    if(strcmp(opname, "demosaic") == 0) demosaic_priority_old = priority;
    if(strcmp(opname, "colorin") == 0) colorin_priority_old = priority;
    if(strcmp(opname, "colorout") == 0) colorout_priority_old = priority;
    if(strcmp(opname, "gamma") == 0) gamma_priority_old = priority;
  }
  sqlite3_finalize(stmt);

  if(order_changed) dev->iop = g_list_sort(dev->iop, dt_sort_iop_by_order);
  if(order_changed) printf("[dt_iop_priorities_read_pipe] order has changed\n");
  // dt_iop_priorities_check_priorities(dev, "dt_iop_priorities_read_pipe 1");

  if(demosaic_priority_old < 0 || colorin_priority_old < 0 || colorout_priority_old < 0 || gamma_priority_old < 0)
  {
    if(demosaic_priority_old < 0)
      fprintf(stderr, "[dt_iop_priorities_read_pipe] can't read old priority for module demosaic on image %i\n",
              imgid);
    if(colorin_priority_old < 0)
      fprintf(stderr, "[dt_iop_priorities_read_pipe] can't read old priority for module colorin on image %i\n",
              imgid);
    if(colorout_priority_old < 0)
      fprintf(stderr, "[dt_iop_priorities_read_pipe] can't read old priority for module colorout on image %i\n",
              imgid);
    if(gamma_priority_old < 0)
      fprintf(stderr, "[dt_iop_priorities_read_pipe] can't read old priority for module gamma on image %i\n",
              imgid);
  }

  // now we have to check if the priorities has changed
  // we will use demosaic, colorin, colorout and gamma for this
  int priorities_changed = 0;
  modules = g_list_first(dev->iop);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;

    if(strcmp(mod->op, "demosaic") == 0)
    {
      if(mod->priority != demosaic_priority_old) priorities_changed = 1;
      demosaic_priority = mod->priority;
      demosaic_iop_order = mod->iop_order;
      demosaic_module = mod;
    }
    if(strcmp(mod->op, "colorin") == 0)
    {
      if(mod->priority != colorin_priority_old) priorities_changed = 1;
      colorin_priority = mod->priority;
      colorin_iop_order = mod->iop_order;
      colorin_module = mod;
    }
    if(strcmp(mod->op, "colorout") == 0)
    {
      if(mod->priority != colorout_priority_old) priorities_changed = 1;
      colorout_priority = mod->priority;
      colorout_iop_order = mod->iop_order;
      colorout_module = mod;
    }
    if(strcmp(mod->op, "gamma") == 0)
    {
      if(mod->priority != gamma_priority_old) priorities_changed = 1;
      gamma_iop_order = mod->iop_order;
      gamma_module = mod;
    }

    modules = g_list_next(modules);
  }

  // priorities has changed since the history was saved
  // this means that modules that are not in the history or in the main.pipe table
  // can be out of order
  // since the user may have changed the order of any module we don't have a way
  // to know the proper place for this modules
  // but at least we can place them in the right colorspace
  if(priorities_changed || order_changed)
  {
    if(priorities_changed) printf("[dt_iop_priorities_read_pipe] priorities has changed\n");

    // go through all modules and check if they are in the right colorspace
    // but only if they are not in the history or main.pipe
    modules = g_list_first(dev->iop);
    while(modules)
    {
      dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;

      int found_history = (_search_history_by_op_multi_priority(dev, mod) != NULL);

      int found_pipe = (_search_pipe_by_op_multi_priority(main_pipe_list, mod) != NULL);

      int reset_list = 0;

      // module is not in history or pipe, let's check if it is in the right colorspace
      if(!found_history && !found_pipe)
      {
        // let's check where it should be
        if(mod->priority < demosaic_priority) // should be raw
        {
          if(mod->iop_order >= demosaic_iop_order)
          {
            if(dt_move_iop_before(dev, mod, demosaic_module, 0, log_error))
              reset_list = 1;
            else
              fprintf(stderr, "[dt_iop_priorities_read_pipe] can't move module %s(%f) before demosaic(%f)\n",
                      mod->op, mod->iop_order, demosaic_iop_order);
          }
        }
        else if(mod->priority < colorin_priority) // should be rgb
        {
          if(mod->iop_order >= colorin_iop_order)
          {
            if(dt_move_iop_before(dev, mod, colorin_module, 0, log_error))
              reset_list = 1;
            else
              fprintf(stderr, "[dt_iop_priorities_read_pipe] can't move module %s(%f) before colorin(%f)\n",
                      mod->op, mod->iop_order, colorin_iop_order);
          }
          else if(mod->iop_order <= demosaic_iop_order)
          {
            if(dt_move_iop_after(dev, mod, demosaic_module, 0, log_error))
              reset_list = 1;
            else
              fprintf(stderr, "[dt_iop_priorities_read_pipe] can't move module %s(%f) after demosaic(%f)\n",
                      mod->op, mod->iop_order, demosaic_iop_order);
          }
        }
        else if(mod->priority < colorout_priority) // should be lab
        {
          if(mod->iop_order >= colorout_iop_order)
          {
            if(dt_move_iop_before(dev, mod, colorout_module, 0, log_error))
              reset_list = 1;
            else
              fprintf(stderr, "[dt_iop_priorities_read_pipe] can't move module %s(%f) before colorout(%f)\n",
                      mod->op, mod->iop_order, colorout_iop_order);
          }
          else if(mod->iop_order <= colorin_iop_order)
          {
            if(dt_move_iop_after(dev, mod, colorin_module, 0, log_error))
              reset_list = 1;
            else
              fprintf(stderr, "[dt_iop_priorities_read_pipe] can't move module %s(%f) after colorin(%f)\n",
                      mod->op, mod->iop_order, colorin_iop_order);
          }
        }
        else // should be rgb
        {
          if(mod->iop_order >= gamma_iop_order)
          {
            if(dt_move_iop_before(dev, mod, gamma_module, 0, log_error))
              reset_list = 1;
            else
              fprintf(stderr, "[dt_iop_priorities_read_pipe] can't move module %s(%f) before gamma(%f)\n", mod->op,
                      mod->iop_order, gamma_iop_order);
          }
          else if(mod->iop_order <= colorout_iop_order)
          {
            if(dt_move_iop_after(dev, mod, colorout_module, 0, log_error))
              reset_list = 1;
            else
              fprintf(stderr, "[dt_iop_priorities_read_pipe] can't move module %s(%f) after colorout(%f)\n",
                      mod->op, mod->iop_order, colorout_iop_order);
          }
        }
      } // if(!found_history && !found_pipe)

      if(reset_list)
        modules = g_list_first(dev->iop);
      else
        modules = g_list_next(modules);
    } // while(modules)
  }   // if(priorities_changed)

  if(main_pipe_list) g_list_free_full(main_pipe_list, free);

  dt_iop_priorities_check_priorities(dev, "dt_iop_priorities_read_pipe end");
  // printf("dt_iop_priorities_read_pipe end\n");
}
#endif

typedef struct _iop_priorities_read_pipe_t
{
  int order_changed;
  int priorities_changed;

  int demosaic_priority;
  int colorin_priority;
  int colorout_priority;

  float demosaic_iop_order;
  float colorin_iop_order;
  float colorout_iop_order;
  float gamma_iop_order;

  dt_iop_module_t *demosaic_module;
  dt_iop_module_t *colorin_module;
  dt_iop_module_t *colorout_module;
  dt_iop_module_t *gamma_module;

  int demosaic_priority_old;
  int colorin_priority_old;
  int colorout_priority_old;
  int gamma_priority_old;

  GList *fences_list;
} _iop_priorities_read_pipe_t;

static void _read_pipe(dt_develop_t *dev, _iop_priorities_read_pipe_t *params)
{
  // read all non-history modules
  GList *main_pipe = g_list_first(dev->main_pipe_list);
  while(main_pipe)
  {
    dt_iop_pipe_entry_t *pipe_entry = (dt_iop_pipe_entry_t *)main_pipe->data;
    /*
        const char *opname = (const char *)sqlite3_column_text(stmt, 0);
        const int priority = sqlite3_column_int(stmt, 1);
        const int multi_priority = sqlite3_column_int(stmt, 2);
        const float iop_order = sqlite3_column_double(stmt, 3);
    */
    // printf("[dt_iop_priorities_read_pipe] read module %s\n", opname);

    // let's make sure is not in history
    const int history_found
        = (_search_history_by_op_multi_priority_ext(dev, pipe_entry->operation, pipe_entry->multi_priority)
           != NULL);
    if(!history_found)
    {
      // for each non-history module change the iop_order on dev->iop
      GList *modules = g_list_first(dev->iop);
      while(modules)
      {
        dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
        if(strcmp(mod->op, pipe_entry->operation) == 0 && mod->multi_priority == pipe_entry->multi_priority)
        {
          if(mod->iop_order != pipe_entry->iop_order)
          {
            // let's make sure no other module has the same iop_order
            GList *order_modules = g_list_first(dev->iop);
            while(order_modules)
            {
              dt_iop_module_t *order_mod = (dt_iop_module_t *)order_modules->data;

              // if we find one change the iop_order
              if(order_mod->iop_order == pipe_entry->iop_order)
              {
                // but only if is not in history
                // we don't want to mess with the current edit
                const int found_history
                    = (_search_pipe_by_op_multi_priority(dev->main_pipe_list, order_mod) != NULL
                       || _search_history_by_op_multi_priority(dev, order_mod) != NULL);

                if(!found_history)
                {
                  // this goes after the new one
                  if(order_mod->priority > mod->priority)
                  {
                    dt_iop_module_t *order_mod_next = NULL;
                    GList *order_modules_next = g_list_next(order_modules);
                    if(order_modules_next) order_mod_next = (dt_iop_module_t *)order_modules_next->data;
                    if(order_mod_next)
                      order_mod->iop_order += (order_mod_next->iop_order - order_mod->iop_order) / 2.f;
                    else
                      order_mod->iop_order += 1.f;
                  }
                  // this goes before the new one
                  else
                  {
                    dt_iop_module_t *order_mod_prev = NULL;
                    GList *order_modules_prev = g_list_previous(order_modules);
                    if(order_modules_prev) order_mod_prev = (dt_iop_module_t *)order_modules_prev->data;
                    if(order_mod_prev)
                      order_mod->iop_order -= (order_mod->iop_order - order_mod_prev->iop_order) / 2.f;
                    else
                      order_mod->iop_order *= .5f;
                  }
                }
              }

              order_modules = g_list_next(order_modules);
            }

            /* fprintf(stderr,
                    "[dt_iop_priorities_read_pipe] order changed for %s %s (%i) pipe %f current %f image %i\n",
                    mod->op, mod->multi_name, mod->multi_priority, iop_order, mod->iop_order, imgid); */
            params->order_changed = 1;
            mod->iop_order = pipe_entry->iop_order;
          }

          break;
        }
        modules = g_list_next(modules);
      }
    }

    if(strcmp(pipe_entry->operation, "demosaic") == 0) params->demosaic_priority_old = pipe_entry->priority;
    if(strcmp(pipe_entry->operation, "colorin") == 0) params->colorin_priority_old = pipe_entry->priority;
    if(strcmp(pipe_entry->operation, "colorout") == 0) params->colorout_priority_old = pipe_entry->priority;
    if(strcmp(pipe_entry->operation, "gamma") == 0) params->gamma_priority_old = pipe_entry->priority;

    main_pipe = g_list_next(main_pipe);
  }

  if(params->order_changed) dev->iop = g_list_sort(dev->iop, dt_sort_iop_by_order);
  // dt_iop_priorities_check_priorities(dev, "dt_iop_priorities_read_pipe 1");
}

static void _check_changed_priorities(dt_develop_t *dev, _iop_priorities_read_pipe_t *params)
{
  params->demosaic_priority = -1;
  params->colorin_priority = -1;
  params->colorout_priority = -1;

  params->demosaic_iop_order = -1.f;
  params->colorin_iop_order = -1.f;
  params->colorout_iop_order = -1.f;
  params->gamma_iop_order = -1.f;

  GList *modules = g_list_first(dev->iop);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;

    if(strcmp(mod->op, "demosaic") == 0)
    {
      if(mod->priority != params->demosaic_priority_old) params->priorities_changed = 1;
      params->demosaic_priority = mod->priority;
      params->demosaic_iop_order = mod->iop_order;
      params->demosaic_module = mod;
    }
    if(strcmp(mod->op, "colorin") == 0)
    {
      if(mod->priority != params->colorin_priority_old) params->priorities_changed = 1;
      params->colorin_priority = mod->priority;
      params->colorin_iop_order = mod->iop_order;
      params->colorin_module = mod;
    }
    if(strcmp(mod->op, "colorout") == 0)
    {
      if(mod->priority != params->colorout_priority_old) params->priorities_changed = 1;
      params->colorout_priority = mod->priority;
      params->colorout_iop_order = mod->iop_order;
      params->colorout_module = mod;
    }
    if(strcmp(mod->op, "gamma") == 0)
    {
      if(mod->priority != params->gamma_priority_old) params->priorities_changed = 1;
      params->gamma_iop_order = mod->iop_order;
      params->gamma_module = mod;
    }

    modules = g_list_next(modules);
  }
}

static int _relocate_non_history_module(dt_develop_t *dev, GList *modules, dt_iop_module_t *mod,
                                        _iop_priorities_read_pipe_t *params, const int log_error)
{
  int relocated = 0;

  // printf("[_relocate_non_history_module] relocating %s %s(%f)\n", mod->op, mod->multi_name, mod->iop_order);

  // let's check where it should be
  dt_iop_module_t *mod_prev = NULL;
  dt_iop_module_t *mod_next = NULL;

  // first check if it is in the right colorspace
  if(mod->priority < params->demosaic_priority) // should be raw
  {
    if(mod->iop_order >= params->demosaic_iop_order)
    {
      if(mod_next == NULL || params->demosaic_module->priority < mod_next->priority)
        mod_next = params->demosaic_module;
    }
  }
  else if(mod->priority < params->colorin_priority) // should be linear rgb
  {
    if(mod->iop_order >= params->colorin_iop_order)
    {
      if(mod_next == NULL || params->colorin_module->priority < mod_next->priority)
        mod_next = params->colorin_module;
    }
    else if(mod->iop_order <= params->demosaic_iop_order)
    {
      if(mod_prev == NULL || params->demosaic_module->priority > mod_prev->priority)
        mod_prev = params->demosaic_module;
    }
  }
  else if(mod->priority < params->colorout_priority) // should be lab
  {
    if(mod->iop_order >= params->colorout_iop_order)
    {
      if(mod_next == NULL || params->colorout_module->priority < mod_next->priority)
        mod_next = params->colorout_module;
    }
    else if(mod->iop_order <= params->colorin_iop_order)
    {
      if(mod_prev == NULL || params->colorin_module->priority > mod_prev->priority)
        mod_prev = params->colorin_module;
    }
  }
  else // should be gamma rgb
  {
    if(mod->iop_order >= params->gamma_iop_order)
    {
      if(mod_next == NULL || params->gamma_module->priority < mod_next->priority) mod_next = params->gamma_module;
    }
    else if(mod->iop_order <= params->colorout_iop_order)
    {
      if(mod_prev == NULL || params->colorout_module->priority > mod_prev->priority)
        mod_prev = params->colorout_module;
    }
  }

  // now check if it doesn't break any rules

  // check if mod is between the modules with IOP_FLAGS_FENCE flag
  dt_iop_module_t *fence_prev = NULL;
  dt_iop_module_t *fence_next = NULL;

  GList *mod_fences = g_list_first(params->fences_list);
  while(mod_fences)
  {
    dt_iop_module_t *mod_fence = (dt_iop_module_t *)mod_fences->data;

    // mod should be before this fence
    if(mod->priority < mod_fence->priority)
    {
      if(fence_next == NULL || mod_fence->priority < fence_next->priority) fence_next = mod_fence;
    }
    // mod should be after this fence
    else if(mod->priority > mod_fence->priority)
    {
      if(fence_prev == NULL || mod_fence->priority > fence_prev->priority) fence_prev = mod_fence;
    }

    mod_fences = g_list_next(mod_fences);
  }

  // now check if mod is between the fences
  if(fence_next && mod->iop_order > fence_next->iop_order)
  {
    if(mod_next == NULL || fence_next->priority < mod_next->priority) mod_next = fence_next;
  }
  if(fence_prev && mod->iop_order < fence_prev->iop_order)
  {
    if(mod_prev == NULL || fence_prev->priority > mod_prev->priority) mod_prev = fence_prev;
  }

  // for each module check if it doesn't break a rule
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
        dt_iop_module_t *rule_mod_prev = (dt_iop_module_t *)modules_prev->data;

        if(rule_mod_prev->iop_order != FLT_MAX)
        {
          if(strcmp(rule_mod_prev->op, rule->op_next) == 0)
          {
            if(mod_next == NULL || rule_mod_prev->priority < mod_next->priority) mod_next = rule_mod_prev;
          }
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
        dt_iop_module_t *rule_mod_next = (dt_iop_module_t *)modules_next->data;

        if(rule_mod_next->iop_order != FLT_MAX)
        {
          if(strcmp(rule_mod_next->op, rule->op_prev) == 0)
          {
            if(mod_prev == NULL || rule_mod_next->priority > mod_prev->priority) mod_prev = rule_mod_next;
          }
        }

        modules_next = g_list_next(modules_next);
      }
    }

    rules = g_list_next(rules);
  }

  // now let's move mod if needed
  if(mod_next)
  {
    if(dt_move_iop_before(dev, mod, mod_next, 0, log_error))
      relocated = 1;
    else
      fprintf(stderr, "[dt_iop_priorities_read_pipe] can't move module %s %s(%f) before %s %s(%f)\n", mod->op,
              mod->multi_name, mod->iop_order, mod_next->op, mod_next->multi_name, mod_next->iop_order);
  }
  else if(mod_prev)
  {
    if(dt_move_iop_after(dev, mod, mod_prev, 0, log_error))
      relocated = 1;
    else
      fprintf(stderr, "[dt_iop_priorities_read_pipe] can't move module %s %s(%f) after %s %s(%f)\n", mod->op,
              mod->multi_name, mod->iop_order, mod_prev->op, mod_prev->multi_name, mod_prev->iop_order);
  }

  return relocated;
}

void dt_iop_priorities_read_pipe(dt_develop_t *dev, const int imgid)
{
  sqlite3_stmt *stmt;

  if(dev->main_pipe_list)
  {
    g_list_free_full(dev->main_pipe_list, free);
    dev->main_pipe_list = NULL;
  }

  DT_DEBUG_SQLITE3_PREPARE_V2(dt_database_get(darktable.db),
                              "SELECT operation, priority, multi_priority, iop_order "
                              "FROM main.pipe WHERE imgid = ?1",
                              -1, &stmt, NULL);
  DT_DEBUG_SQLITE3_BIND_INT(stmt, 1, imgid);
  while(sqlite3_step(stmt) == SQLITE_ROW)
  {
    // save it to the list
    dt_iop_pipe_entry_t *pipe_entry = calloc(1, sizeof(dt_iop_pipe_entry_t));

    snprintf(pipe_entry->operation, sizeof(pipe_entry->operation), "%s",
             (const char *)sqlite3_column_text(stmt, 0));
    pipe_entry->priority = sqlite3_column_int(stmt, 1);
    pipe_entry->multi_priority = sqlite3_column_int(stmt, 2);
    pipe_entry->iop_order = sqlite3_column_double(stmt, 3);

    dev->main_pipe_list = g_list_append(dev->main_pipe_list, pipe_entry);
  }
  sqlite3_finalize(stmt);
}

void dt_iop_priorities_relocate_non_history_modules(dt_develop_t *dev)
{
  // printf("dt_iop_priorities_relocate_non_history_modules begin\n");
  dt_times_t start;
  dt_get_times(&start);

  // dt_iop_priorities_check_priorities(dev, "dt_iop_priorities_relocate_non_history_modules begin");

  // printf("[dt_iop_priorities_relocate_non_history_modules] history_end=%i\n", dev->history_end);

  // GList *modules = NULL;

  const int log_error = 1;
  /*
    int params->demosaic_priority = -1;
    int params->colorin_priority = -1;
    int params->colorout_priority = -1;

    float params->demosaic_iop_order = -1.f;
    float params->colorin_iop_order = -1.f;
    float params->colorout_iop_order = -1.f;
    float params->gamma_iop_order = -1.f;

    dt_iop_module_t *params->demosaic_module = NULL;
    dt_iop_module_t *params->colorin_module = NULL;
    dt_iop_module_t *params->colorout_module = NULL;
    dt_iop_module_t *params->gamma_module = NULL;

    int params->demosaic_priority_old = -1;
    int params->colorin_priority_old = -1;
    int params->colorout_priority_old = -1;
    int params->gamma_priority_old = -1;

    GList *main_pipe_list = NULL;
  */

  _iop_priorities_read_pipe_t _params = { 0 };
  _iop_priorities_read_pipe_t *params = &_params;

  // int order_changed = 0;

  params->demosaic_priority_old = -1;
  params->colorin_priority_old = -1;
  params->colorout_priority_old = -1;
  params->gamma_priority_old = -1;

  // read the last state of the pipe and update the iop_order on dev->iop
  // this will ensure that the pipe order is correct even for modules that are not in the history
  // main.pipe has only entries that are not in history but enabled anyway
  _read_pipe(dev, params);

  // check if we have all the priorities we need
  if(params->demosaic_priority_old < 0 || params->colorin_priority_old < 0 || params->colorout_priority_old < 0
     || params->gamma_priority_old < 0)
  {
    if(params->demosaic_priority_old < 0)
      fprintf(stderr,
              "[dt_iop_priorities_relocate_non_history_modules] can't read old priority for module demosaic\n");
    if(params->colorin_priority_old < 0)
      fprintf(stderr,
              "[dt_iop_priorities_relocate_non_history_modules] can't read old priority for module colorin\n");
    if(params->colorout_priority_old < 0)
      fprintf(stderr,
              "[dt_iop_priorities_relocate_non_history_modules] can't read old priority for module colorout\n");
    if(params->gamma_priority_old < 0)
      fprintf(stderr,
              "[dt_iop_priorities_relocate_non_history_modules] can't read old priority for module gamma\n");
  }

  // now we have to check if the priorities has changed
  // we will use demosaic, colorin, colorout and gamma for this
  _check_changed_priorities(dev, params);

  // priorities has changed since the history was saved
  // this means that modules that are not in the history or in the main.pipe table
  // can be out of order
  // since the user may have changed the order of any module we don't have a way
  // to know the proper place for this modules
  // but at least we can place them in the right colorspace and ensure all the rules
  if(params->priorities_changed || params->order_changed)
  {
    // if(params->priorities_changed) printf("[dt_iop_priorities_relocate_non_history_modules] priorities has
    // changed\n");
    // if(params->order_changed) printf("[dt_iop_priorities_relocate_non_history_modules] order has changed\n");

    // create a list of fences modules
    params->fences_list = _get_fence_modules_list(dev);

    // go through all modules and check if they are in the right place in the pipe
    // but only if they are not in the history or main.pipe
    GList *modules = g_list_first(dev->iop);
    while(modules)
    {
      dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;

      int reset_list = 0;

      if(mod->iop_order != FLT_MAX)
      {
        const int found_history = (_search_pipe_by_op_multi_priority(dev->main_pipe_list, mod) != NULL
                                   || _search_history_by_op_multi_priority(dev, mod) != NULL);

        // module is not in history or pipe, let's check if it is in the right place in the pipe
        // and reposition it if needed
        if(!found_history)
        {
          reset_list = _relocate_non_history_module(dev, modules, mod, params, log_error);
          /*
                    if(reset_list)
                    {
                      printf("[dt_iop_priorities_relocate_non_history_modules] relocated %s %s(%f)\n", mod->op,
             mod->multi_name, mod->iop_order);
                    }*/
        }
        /*        else
                  printf("[dt_iop_priorities_relocate_non_history_modules] found in history %s %s(%f)\n", mod->op,
           mod->multi_name, mod->iop_order);
        */
      }

      // reposition changes the order of the pipe
      // so we may need to start over
      if(reset_list)
        modules = g_list_first(dev->iop);
      else
        modules = g_list_next(modules);
    }
  }

  // if(params->main_pipe_list) g_list_free_full(params->main_pipe_list, free);
  if(params->fences_list) g_list_free(params->fences_list);

  dt_iop_priorities_check_priorities(dev, "dt_iop_priorities_relocate_non_history_modules end");
  /*
    {
      GList *modules = g_list_first(dev->iop);
      while(modules)
      {
        dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;

        printf("[dt_iop_priorities_relocate_non_history_modules] %s %s(%f)\n", mod->op, mod->multi_name,
    mod->iop_order);

        modules = g_list_next(modules);
      }
    }
  */
  dt_show_times(&start, "[dt_iop_priorities_relocate_non_history_modules] read pipe", NULL);

  // printf("dt_iop_priorities_read_pipe end\n");
}

void dt_iop_priorities_write_pipe(dt_develop_t *dev, const int imgid)
{
  dt_iop_priorities_check_priorities(dev, "dt_iop_priorities_write_pipe begin");
  // printf("[dt_iop_priorities_write_pipe] begin\n");

  sqlite3_stmt *stmt;

  // clear the pipe first
  DT_DEBUG_SQLITE3_PREPARE_V2(dt_database_get(darktable.db), "DELETE FROM main.pipe WHERE imgid = ?1", -1, &stmt,
                              NULL);
  DT_DEBUG_SQLITE3_BIND_INT(stmt, 1, imgid);
  // sqlite3_step(stmt);
  if(sqlite3_step(stmt) != SQLITE_DONE)
  {
    fprintf(stderr, "[dt_iop_priorities_write_pipe] error deleting pipe for image %i\n", imgid);
    fprintf(stderr, "%s\n", sqlite3_errmsg(dt_database_get(darktable.db)));
  }
  sqlite3_finalize(stmt);

  // go through all modules and write the default_enabled that not exists on history
  GList *modules = g_list_first(dev->iop);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)modules->data;
    if(mod->default_enabled)
    {
      /*
            int history_found = 0;
            // always save this modules, will be used to check if priorities has changed
            if(strcmp(mod->op, "demosaic") != 0 && strcmp(mod->op, "colorin") != 0 && strcmp(mod->op, "colorout")
         != 0
               && strcmp(mod->op, "gamma") != 0)
            {
              history_found = (_search_history_by_op_multi_priority(dev, mod) != NULL);
            }

            if(!history_found)
      */
      {
        DT_DEBUG_SQLITE3_PREPARE_V2(
            dt_database_get(darktable.db),
            "INSERT INTO main.pipe (imgid, operation, priority, multi_priority, iop_order) "
            "VALUES (?1, ?2, ?3, ?4, ?5)",
            -1, &stmt, NULL);
        DT_DEBUG_SQLITE3_BIND_INT(stmt, 1, imgid);
        DT_DEBUG_SQLITE3_BIND_TEXT(stmt, 2, mod->op, -1, SQLITE_TRANSIENT);
        DT_DEBUG_SQLITE3_BIND_INT(stmt, 3, mod->priority);
        DT_DEBUG_SQLITE3_BIND_INT(stmt, 4, mod->multi_priority);
        DT_DEBUG_SQLITE3_BIND_DOUBLE(stmt, 5, mod->iop_order);
        // sqlite3_step(stmt);
        if(sqlite3_step(stmt) != SQLITE_DONE)
        {
          fprintf(stderr, "[dt_iop_priorities_write_pipe] error inserting operation %s on pipe for image %i\n",
                  mod->op, imgid);
          fprintf(stderr, "%s\n", sqlite3_errmsg(dt_database_get(darktable.db)));
        }
        sqlite3_finalize(stmt);
      }
    }
    modules = g_list_next(modules);
  }
  // printf("[dt_iop_priorities_write_pipe] end\n");
  dt_iop_priorities_check_priorities(dev, "dt_iop_priorities_write_pipe end");
}

int dt_iop_priorities_check_priorities(dt_develop_t *dev, const char *msg)
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
          fprintf(stderr, "[dt_iop_priorities_check_priorities] module not used but enabled!! %s %s(%f) (%s)\n",
                  mod->op, mod->multi_name, mod->iop_order, msg);
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
        fprintf(stderr, "[dt_iop_priorities_check_priorities] gamma is not the last iop, last is %s %s(%f) (%s)\n",
                mod->op, mod->multi_name, mod->iop_order, msg);
      }
    }
    else
    {
      // fprintf(stderr, "[dt_iop_priorities_check_priorities] dev->iop is empty (%s)\n", msg);
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
                    "[dt_iop_priorities_check_priorities] module %s %s(%f) should be after %s %s(%f) (%s)\n",
                    mod->op, mod->multi_name, mod->iop_order, mod_prev->op, mod_prev->multi_name,
                    mod_prev->iop_order, msg);
          }
          else if(mod->iop_order == mod_prev->iop_order)
          {
            priorities_ok = 0;
            fprintf(
                stderr,
                "[dt_iop_priorities_check_priorities] module %s %s(%f) and %s %s(%f) has the same order (%s)\n",
                mod->op, mod->multi_name, mod->iop_order, mod_prev->op, mod_prev->multi_name, mod_prev->iop_order,
                msg);
          }
        }
      }
      mod_prev = mod;
      modules = g_list_next(modules);
    }
  }

  _check_rules(dev, msg);

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
  //  printf("[dt_iop_transform_image_colorspace_cl] transforming %s %s from %i to %i\n", self->op,
  //  self->multi_name,
  //         cst_from, cst_to);
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

dt_iop_module_t *dt_iop_priorities_get_iop_by_op(dt_develop_t *dev, const char *op_name)
{
  dt_iop_module_t *mod_ret = NULL;
  GList *modules = g_list_first(dev->iop);
  while(modules)
  {
    dt_iop_module_t *mod = (dt_iop_module_t *)(modules->data);

    if(strcmp(mod->op, op_name) == 0)
    {
      mod_ret = mod;
      break;
    }
    modules = g_list_next(modules);
  }
  return mod_ret;
}
