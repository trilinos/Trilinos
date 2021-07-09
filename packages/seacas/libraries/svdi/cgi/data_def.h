/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef DATA_DEF_H
#define DATA_DEF_H
/* ------------------------------------------------------------*/
/* >> DATA STRUCTURE DEFINITIONS                               */
/*-------------------------------------------------------------*/

#include "cgi.h"

/* define a structure to store SVDI attributes */
/* --keeps track of SVDI attribs which have been set for each surface */
/* --these values are stored in CGI terms, not SVDI/NDC terms */
typedef struct
{
  int   fg_color;   /* foreground color for vector */
  int   bg_color;   /* background color for vector */
  float intensity;  /* current intensity */
  int   line_type;  /* current line style */
  float line_width; /* current line width */
  float char_x;     /* current char box size in x */
  float char_y;     /* current char box size in y */
  float fg_rgb[3];  /* foreground color for raster */
  float bg_rgb[3];  /* background color for raster */
} vdi_attrib_struct;

/* state list structure for this device */
/* -- stores state of each surface for this device */
typedef struct
{

  /* cstate.i - definitions for cgi part of state lists for devices
   * 8 Sep 1989 - last date modified
   * Pat McGee, jpm@lanl.gov
   */

  int          cgi_inited;                /* 0 for uninitialized, 1 for inited */
  int          err_flag[MAX_ERROR_CLASS]; /* on or off for each of 7 classes */
  error_report err_queue[ERROR_LIST_SIZE];
  int          err_head_ptr; /* head of error queue */
  int          err_count;    /* count of errors in queue */

  /* ----------control, negotiation and errors---------- */
  /* control state list */
  /* -- control state: */
  int pic_dirty;        /* 0 = clean, 1 = dirty */
  int color_index_prec; /* -1 for bad */

  /* -- VDC to Device mapping: */
  point vdc1;    /* corner 1 of vdc extent */
  point vdc2;    /* corner 2 of vdc extent */
  point vp1;     /* corner 1 of viewport (mode:fraction view surf)*/
  point vp2;     /* corner 1 of viewport (mode:fraction view surf)*/
  point eff_vp1; /* corner 1 of effective viewport */
  point eff_vp2; /* corner 2 of effective viewport */

  /* -- clipping: */
  cenum clip_indicator;    /* on, off  */
  point clip_rect1;        /* corner 1 of requested clip rectangle */
  point clip_rect2;        /* corner 2 of requested clip rectangle */
  cenum ds_clip_indicator; /* off, draw surf clip rect, viewport */

  /* ----------output and attributes---------- */
  /* attribute state lists */
  /* line attributes state list */
  int   line_type;
  float line_width;
  int   line_color[3]; /* use elem 0 for indexed, all 3 for direct*/
  cenum line_csm;      /* csm in which line color was last specified */

  /* marker attributes state list */
  int   mark_type;
  float mark_size;
  int   mark_color[3]; /* use elem 0 for indexed, all 3 for direct*/
  cenum mark_csm;      /* csm in which mark color was last specified */

  /* text attributes state list */
  int   text_color[3]; /* use elem 0 for indexed, all 3 for direct*/
  cenum text_csm;      /* csm in which text color was last specified */
  float char_height;

  /* fill attributes state list */
  int   fill_color[3];  /* use elem 0 for indexed, all 3 for direct*/
  cenum fill_csm;       /* csm in which fill color was last specified */
  int   interior_style; /* 0 = hollow, 1 = solid */

  /* output and attribute control state list */
  cenum csm; /* current color mode:
                  0=direct, 1=indexed, -1=bad*/
  rgb color_table[COLOR_TABLE_SIZE];

  /* end cstate.i */

  int next_free_state; /* for free state list */
  int this_index;      /* index of statelist in array */

  /* SVDI dependent stuff */
  float             xscale;     /* x scale to map VDC to NDC viewport */
  float             yscale;     /* y scale to map VDC to NDC viewport */
  float             xoffset;    /* x offset to map VDC to NDC v.p. */
  float             yoffset;    /* y offset to map VDC to NDC v.p. */
  int               color_set;  /* TRUE = color table set by user */
  vdi_attrib_struct vdi_attrib; /* current vdi attributes */

  /* other info */
  point eff_clip_rect1; /* corner 1 of effective clip rectangle */
  point eff_clip_rect2; /* corner 2 of effective clip rectangle */
  point clipmin;        /* min corner of actual clip region - NDC */
  point clipmax;        /* max corner of actual clip region - NDC */
  int   clip_on;        /* logical - TRUE = clipping is on */
  int   fg_index;       /* for mapping foreground and background */
  int   bg_index;       /*  indices */

  /* device dependent stuff - batch device */
  /* -- BUFFER_SIZE is 0 if interactive device */
  char filename[100];       /* name of file */
  int  file_d;              /* file descriptor - -1 if file not open */
  char buffer[BUFFER_SIZE]; /* for buffering data for output to file */
  int  buff_ptr;            /* keep track of buffer space */

  /* device dependent stuff - interactive device  */
  /* a state list exists for each logical input device as identified
   * by its input class and its input device index.  we support only
   * locator, and 1 input device index.
   */
  int   input_dev_class; /* we only support locator */
  int   input_dev_index; /* we only support 1 */
  cenum input_dev_state; /* ready or released */

} surf_statelist;

/* description table */
/* -- stores description table for this device  - one per device */
typedef struct
{

  /* cdescr.i - definitions for cgi description table */
  /* 9 Sep 1989 - last date modified
   */

  /* ----------control, negotiation and errors---------- */
  /* device identity description table  */
  cenum dev_class; /* device class: IN,OUT,INOUT */
  char  dev_id[4]; /* device identification */

  /* output device description table  */
  cenum  copy_class;       /* hard, soft */
  cenum  display_type;     /* vector,raster,other */
  cenum  bcolor_cap;       /* background color capability */
  cenum  dynamic_mod_bg;   /* dynamic mod accept for bg color */
  cenum  dynamic_mod_map;  /* dynamic mod accept for mapping */
  dpoint dc1;              /* device bottom left corner */
  dpoint dc2;              /* device upper right corner */
  float  draw_surf_width;  /* measured in millimeters */
  float  draw_surf_height; /* measured in millimeters */
  cenum  pix_loc;          /* pixel location */

  /* ----------output and attributes---------- */
  /* attribute description tables */
  /* primitive support description table */
  int   max_pts_polyline;      /* -1,128..n  */
  int   max_pts_disj_polyline; /* -1,0..n    */
  int   max_pts_polygon;       /* -1,0..n    */
  int   max_pts_polygon_set;   /* -1,0..n    */
  int   max_pts_polymarker;    /* -1,128..n  */
  int   max_pts_cellarray;     /* -1,0..n    */
  int   max_pts_closed_fig;    /* -1,0..n    */
  int   max_chars_text;        /* 80..n      */
  cenum cellarray_fill_cap;    /* none,global,local */
  cenum cellarray_align_cap;   /* axis,skewed */
  cenum compound_text_cap;     /* none,global,local */
  cenum compound_fig_cap;      /* none,global,local */

  /* line description table  */
  int linewidth_nominal; /* nominal scaled line width ( >0 ) */
  int linewidth_min;     /* minimum scaled line width ( >0 ) */
  int linewidth_max;     /* maximum scaled line width ( >0 ) */
  int line_types[5];     /* list of available line types */
  int line_widths[5];    /* list of available scale line widths */

  /* marker description table */
  int mark_nominal;  /* nominal scaled marker size ( >0 ) */
  int mark_min;      /* minimum scaled marker size ( >0 ) */
  int mark_max;      /* maximum scaled marker size ( >0 ) */
  int mark_types[1]; /* list of available marker types */
  int mark_sizes[1]; /* list of available scale marker sizes */

  /* text description table */
  /* fill description table */
  cenum interior_styles[2]; /* list of available interior styles */

  /* attribute and control description table */
  int   num_simul_colors;  /* 2..n */
  int   num_avail_colors;  /* 2..n */
  int   num_avail_int;     /* if direct, use all 3 */
  cenum csm_avail;         /* indexed only, indexed and direct */
  cenum dynamic_mod_ct;    /* irg, cbs, imm */
  cenum color_overwrite;   /* irg, cbs, imm */
  cenum monochrome_device; /* no,yes */

  /* end cdescr.i */

  /* SVDI part of description table */
  /* --used to store default SVDI setup and to reinitialize */
  int   set_flag;            /* TRUE = SVDI defaults set, FALSE = not */
  int   svdi_type;           /* 0 = SVDI, 1 = SVDI+raster */
  float xndc_max;            /* max x ndc value */
  float yndc_max;            /* max y ndc value */
  float att_array[14];       /* output status of attributes */
                             /* --stuff stored in CGI terms */
  int   num_cols;            /* number of colors */
  int   col_mode;            /* color mode */
  int   index_array[256];    /* color index array */
  float color_array[256][3]; /* color values */

} dev_descrip_table;

#endif
