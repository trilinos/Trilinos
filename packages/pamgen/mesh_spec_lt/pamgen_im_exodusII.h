// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

/*****************************************************************************
 *
 * exodusII.h - Exodus II include file, for general use
 *
 * author - Sandia National Laboratories
 *          
 * environment - UNIX
 *
 * exit conditions - 
 *
 * revision history - 
 *
 *****************************************************************************/


#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0 
#endif

#ifndef IM_EXODUS_II_HDR
#define IM_EXODUS_II_HDR

/*
 * need following extern if this include file is used in a C++ program, to
 * keep the C++ compiler from mangling the function names.
 */
#ifdef __cplusplus
extern "C" {
#endif

  /*
   * The following are miscellaneous constants used in the EXODUS II API.
   */

#define IM_EX_NOCLOBBER            0 /* Don't overwrite existing database, default */
#define IM_EX_CLOBBER              1
#define IM_EX_NORMAL_MODEL         2 /* disable mods that permit storage of larger models */
#define IM_EX_LARGE_MODEL          4 /* enable mods that permit storage of larger models */
#define IM_EX_NETCDF4              8 /* use the hdf5-based netcdf4 output */
#define IM_EX_NOSHARE             16 /* Do not open netcdf file in "share" mode */
#define IM_EX_SHARE               32 /* Do open netcdf file in "share" mode */

#define IM_EX_READ                 0
#define IM_EX_WRITE                1

#define IM_EX_INQ_FILE_TYPE        1               /* inquire EXODUS II file type*/
#define IM_EX_INQ_API_VERS         2               /* inquire API version number */
#define IM_EX_INQ_DB_VERS          3               /* inquire database version   */
                                                /*   number                   */
#define IM_EX_INQ_TITLE            4               /* inquire database title     */
#define IM_EX_INQ_DIM              5               /* inquire number of          */
                                                /*   dimensions               */
#define IM_EX_INQ_NODES            6               /* inquire number of nodes    */
#define IM_EX_INQ_ELEM             7               /* inquire number of elements */
#define IM_EX_INQ_ELEM_BLK         8               /* inquire number of element  */
                                                /*   blocks                   */
#define IM_EX_INQ_NODE_SETS        9               /* inquire number of node sets*/
#define IM_EX_INQ_NS_NODE_LEN      10              /* inquire length of node set */
                                                /*   node list                */
#define IM_EX_INQ_SIDE_SETS        11              /* inquire number of side sets*/
#define IM_EX_INQ_SS_NODE_LEN      12              /* inquire length of side set */
                                                /*   node list                */
#define IM_EX_INQ_SS_ELEM_LEN      13              /* inquire length of side set */
                                                /*   element list             */
#define IM_EX_INQ_QA               14              /* inquire number of QA       */
                                                /*   records                  */
#define IM_EX_INQ_INFO             15              /* inquire number of info     */
                                                /*   records                  */
#define IM_EX_INQ_TIME             16              /* inquire number of time     */
                                                /*   steps in the database    */
#define IM_EX_INQ_EB_PROP          17              /* inquire number of element  */
                                                /*   block properties         */
#define IM_EX_INQ_NS_PROP          18              /* inquire number of node set */
                                                /*   properties               */
#define IM_EX_INQ_SS_PROP          19              /* inquire number of side set */
#define IM_EX_INQ_NS_DF_LEN        20              /* inquire length of node set */
                                                /*   distribution factor  list*/
#define IM_EX_INQ_SS_DF_LEN        21              /* inquire length of node set */
                                                /*   distribution factor  list*/
#define IM_EX_INQ_LIB_VERS         22              /* inquire API Lib vers number*/
#define IM_EX_INQ_EM_PROP          23              /* inquire number of element  */
                                                /*   map properties           */
#define IM_EX_INQ_NM_PROP          24              /* inquire number of node     */
                                                /*   map properties           */
#define IM_EX_INQ_ELEM_MAP         25              /* inquire number of element  */
                                                /*   maps                     */
#define IM_EX_INQ_NODE_MAP         26              /* inquire number of node     */
                                                /*   maps                     */

  /*   properties               */
#define IM_EX_ELEM_BLOCK           1               /* element block property code*/
#define IM_EX_NODE_SET             2               /* node set property code     */
#define IM_EX_SIDE_SET             3               /* side set property code     */
#define IM_EX_ELEM_MAP             4               /* element map property code  */
#define IM_EX_NODE_MAP             5               /* node map property code     */

  /*   max string lengths; constants that are used as netcdf dimensions must be
       of type long       */
#define MAX_STR_LENGTH          32L
#define MAX_VAR_NAME_LENGTH     20
#define MAX_LINE_LENGTH         80L
#define MAX_ERR_LENGTH          256

  /*   for netCDF 3.4, we estimate the size of the header; 
       if estimate is larger than this max, set the estimate to this max;
       I've never measured a header larger than 20K   */
#define MAX_HEADER_SIZE         30000

  /* routines for file initialization i/o */


  extern int im_ex_get_coord_names (int    exoid,
				 char **coord_names);
  extern int im_ex_get_coord (int exoid,
			   void *x_coor,
			   void *y_coor,
			   void *z_coor);
 
  extern int im_ex_get_ids (int  exoid, int obj_type, int *ids);
  extern int im_ex_get_elem_blk_ids (int  exoid, int *ids);
  extern int im_ex_get_elem_blk_parent_mesh(int  exoid, int *ids);
  extern int im_ex_get_ns_parent_mesh(int  exoid,  int *ids);
  extern int im_ex_get_ss_parent_mesh(int  exoid,  int *ids);
  extern int im_ex_get_elem_block (int   exoid,
				int   elem_blk_id,
				char *elem_type,
				int  *num_elem_this_blk, 
				int  *num_nodes_per_elem,
				int  *num_attr);

  extern int im_ex_get_elem_conn (int   exoid,
			       int   elem_blk_id,
			       int  *connect);


  extern int im_ex_get_elem_num_map (int  exoid,int *elem_map);

  extern int im_ex_get_info (int exoid, char **info);

  extern int im_ex_get_init (int   exoid,
			  char *title,
			  int  *num_dim,
			  int  *num_nodes,
			  int  *num_elem, 
			  int  *num_elem_blk,
			  int  *num_node_sets,
			  int  *num_side_sets);

  extern int im_ex_get_map (int  exoid, int *elem_map);


  extern int im_ex_get_node_map (int   exoid,
			      int   map_id,
			      int  *node_map);

  extern int im_ex_get_node_num_map (int  exoid,
				  int *node_map);

  extern int im_ex_get_node_set_param (int  exoid,
				    int  node_set_id,
				    int *num_nodes_in_set,
				    int *num_df_in_set);

  extern int im_ex_get_node_set (int   exoid,
			      int   node_set_id,
			      int  *node_set_node_list);

  extern int im_ex_get_node_set_ids (int  exoid,
				  int *ids);


  extern int im_ex_get_qa (int exoid,
			char *qa_record[][4]);

  extern int im_ex_get_side_set_param (int  exoid,
				    int  side_set_id,
				    int *num_side_in_set, 
				    int *num_dist_fact_in_set);
  extern int im_ex_get_side_set (int   exoid,
			      int   side_set_id,
			      int  *side_set_elem_list, 
			      int  *side_set_side_list);

  extern int im_ex_get_side_set_node_list(int exoid,
				       int side_set_id,
				       int *side_set_node_cnt_list,
				       int *side_set_node_list);

  extern int im_ex_get_side_set_ids (int  exoid,
				  int *ids);
  extern int im_ex_inquire(int, int, int*, float*, char*);


  
  /* ERROR CODE DEFINITIONS AND STORAGE                                       */
  extern int exerrval;            /* shared error return value                */
  extern int exoptval;            /* error reporting flag (default is quiet)  */
  
#ifdef __cplusplus
}                               /* close brackets on extern "C" declaration */
#endif

#endif

/* im_ex_opts function codes - codes are OR'ed into exopts                     */
#define IM_EX_VERBOSE      1       /* verbose mode message flag                */
#define IM_EX_DEBUG        2       /* debug mode def                           */
#define IM_EX_ABORT        4       /* abort mode flag def                      */

/* Exodus error return codes - exerrval return values:                      */
#define IM_EX_MEMFAIL       1000   /* memory allocation failure flag def       */
#define IM_EX_BADFILEMODE   1001   /* bad file mode def                        */
#define IM_EX_BADFILEID     1002   /* bad file id def                          */
#define IM_EX_WRONGFILETYPE 1003   /* wrong file type for function             */
#define IM_EX_LOOKUPFAIL    1004   /* id table lookup failed                   */
#define IM_EX_BADPARAM      1005   /* bad parameter passed                     */
#define IM_EX_NULLENTITY   -1006   /* null entity found                        */
#define IM_EX_MSG          -1000   /* message print code - no error implied    */
#define IM_EX_PRTLASTMSG   -1001   /* print last error message msg code        */

#include "pamgen_im_exodusII_ext.h"
