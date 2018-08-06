/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef NEM_SPREAD_H
#define NEM_SPREAD_H

#include "globals.h"
#include "nem_spread.h"        // for NemSpread, etc
#include "pe_str_util_const.h" // for strip_string, token_compare, etc
#include "rf_allo.h"
#include "rf_io_const.h"

#define UTIL_NAME "nem_spread"
#define VER_STR "6.16 (2016/08/25)"

extern void   check_exodus_error(int, const char *);
extern double second(void);

template <typename T, typename INT> class NemSpread
{
public:
  void   load_mesh();
  int    check_inp();
  void   read_coord(int exoid, int max_name_length);
  void   read_elem_blk_ids(int mesh_exoid, int max_name_length);
  void   read_elem_blk(int exoid);
  void   extract_elem_blk();
  void   extract_global_element_ids(INT global_ids[], size_t Num_Elem, int iproc);
  void   extract_global_node_ids(INT global_ids[], size_t Num_Node, int iproc);
  size_t extract_elem_connect(INT elem_blk[], int icurrent_elem_blk, size_t istart_elem,
                              size_t iend_elem, int *local_ielem_blk, int iproc);
  void extract_elem_attr(T *elem_attr, int icurrent_elem_blk, size_t istart_elem, size_t iend_elem,
                         int natt_p_elem, int iproc);
  void find_elem_block(INT *proc_elem_blk, int iproc, int proc_for);
  void read_node_set_ids(int mesh_exoid, INT /*num_nodes_in_node_set*/[], INT /*num_df_in_nsets*/[],
                         int max_name_length);
  void read_side_set_ids(int mesh_exoid, INT /*num_elem_in_ssets*/[], INT /*num_df_in_ssets*/[],
                         int max_name_length);
  void read_node_sets(int exoid, INT * /*num_nodes_in_node_set*/, INT * /*num_df_in_nsets*/);
  void read_side_sets(int exoid, INT * /*num_elem_in_ssets*/, INT * /*num_df_in_ssets*/);

  void read_nodal_vars(int mesh_exoid);

  void load_lb_info();

  int write_var_param(int mesh_exoid, int max_name_length, int num_glob, char **gv_names,
                      int num_node, char **nv_names, int num_elem, char **ev_names, int *local_ebtt,
                      int num_nset, char **ns_names, int *local_nstt, int num_sset, char **ss_names,
                      int *local_sstt);

  void write_parExo_data(int mesh_exoid, int max_name_length, int iproc, INT *Num_Nodes_In_NS,
                         INT *Num_Elems_In_SS, INT *Num_Elems_In_EB);

  void write_var_timestep(int exoid, int proc, int time_step, INT *eb_ids_global,
                          INT *ss_ids_global, INT *ns_ids_global);

  void process_lb_data(INT *Integer_Vector, int indx);
  void read_proc_init(int lb_exoid, int proc_info[], int **proc_ids_ptr);
  void read_lb_init(int lb_exoid, INT *Int_Space, INT *Int_Node_Num, INT *Bor_Node_Num,
                    INT *Ext_Node_Num, INT *Int_Elem_Num, INT *Bor_Elem_Num, INT *Node_Comm_Num,
                    INT *Elem_Comm_Num, char *Title);

  void read_cmap_params(int lb_exoid, INT *Node_Comm_Num, INT *Elem_Comm_Num, INT *Num_N_Comm_Maps,
                        INT *Num_E_Comm_Maps, ELEM_COMM_MAP<INT> **E_Comm_Map,
                        NODE_COMM_MAP<INT> **N_Comm_Map, INT *cmap_max_size, INT **comm_vec);

  void create_elem_types();
  void read_mesh_param();
  void read_restart_params();
  void read_restart_data();

  int read_var_param(int exoid, int max_name_length);

  int read_vars(int exoid, int index, INT *eb_ids, INT *eb_cnts, INT ***eb_map_ptr,
                INT **eb_cnts_local, INT *ss_ids, INT *ss_cnts, INT *ns_ids, INT *ns_cnts);
  int read_elem_vars(int exoid, int index, INT *eb_ids, INT *eb_cnts, INT ***eb_map_ptr,
                     INT **eb_cnts_local);
  int read_elem_vars_1(int exoid, int index, INT *eb_ids, INT *eb_cnts, INT ***eb_map_ptr,
                       INT **eb_cnts_local, int iblk, int eb_offset, INT *local_offset);
  int read_sset_vars(int exoid, int index, INT *ss_ids, INT *ss_cnts);
  int read_sset_vars_1(int exoid, int index, INT *ss_ids, INT *ss_cnts, int iset);
  int read_nset_vars(int exoid, int index, INT *ns_ids, INT *ns_cnts);
  int read_nset_vars_1(int exoid, int index, INT *ns_ids, INT *ns_cnts, int iset);
  int read_nodal_vars(int exoid, int index);
  int compare_mesh_param(int exoid);

  int  int64db;
  int  int64api;
  bool force64db; /* Store all ints as 64-bit on output databases. */

  int                    io_ws;
  Restart_Description<T> Restart_Info;
  Globals<T, INT>        globals;

  /*-----------------------------------------------------------------------------
   *
   *  variable definitions used for reading and setting up the mesh
   *  information on the processors.  These variables are only used in
   *  the program in the file, el_exoII_io.c.  Use of these variables
   *  outside this file constitutes an error.
   *
   *----------------------------------------------------------------------------*/

  INT *Node_Set_Ids;           /* Vector of node set ids             *
                                * (ptr to vector of length, 	     *
                                *  Num_Node_Set)    		     */
  INT *Side_Set_Ids;           /* Side set ids  		     *
                                * (ptr to vector of length, 	     *
                                *  Num_Side_Set)      		     */
  INT *Num_Elem_In_Blk;        /* Number of elements in each element *
                                * block (ptr to vector of length,    *
                                *        Num_Elem_Blk)     	     */
  INT *Num_Nodes_Per_Elem;     /* Number of nodes per element in each*
                                * elem block (ptr to vector of       *
                                * length, Num_Elem_Blk)              */
  INT *Num_Attr_Per_Elem;      /* Number of attributes per element in*
                                * each block (ptr to vector of       *
                                * length, Num_Elem_Blk)              */
  INT *Elem_Blk_Ids;           /* Element block id's of each element *
                                * block (ptr to vector of length,    *
                                *        Num_Elem_Blk)               */
  char **Elem_Blk_Types;       /* Element block types for each       *
                                * element block (ptr to vector of    *
                                * length, Num_Elem_Blk, of vectors of*
                                * char of length MAX_STR_LENGTH+1)   */
  char **Elem_Blk_Names;       /* Element block names for each       *
                                * element block (ptr to vector of    *
                                * length, Num_Elem_Blk, of vectors of*
                                * char of length MAX_STR_LENGTH+1)   */
  char **Node_Set_Names;       /* Nodeset names for each             *
                                * nodeset (ptr to vector of          *
                                * length, Num_Node_Set, of vectors of*
                                * char of length MAX_STR_LENGTH+1)   */
  char **Side_Set_Names;       /* Sideset names for each             *
                                * sideset (ptr to vector of          *
                                * length, Num_Side_Set, of vectors of*
                                * char of length MAX_STR_LENGTH+1)   */
  char ***Elem_Blk_Attr_Names; /* Element block attribute names for each
                                * attribute in each element block
                                * ragged array (size varies for each
                                * element block          	     */

  int *GM_Elem_Types; /* This is a list of the element      *
                       * in the entire global mesh. It is   *
                       * stored on Processor 0 to           *
                       * facilitate the reading of the side *
                       * set distribution factors.          */

  char *Coord_Name[3]{}; /* The name(s) of the coordinate axes.   */

  int  Proc_Info[6]{};
  int *Proc_Ids;

  NemSpread()
      : int64db(0), int64api(0), force64db(false), io_ws(0), Node_Set_Ids(nullptr),
        Side_Set_Ids(nullptr), Num_Elem_In_Blk(nullptr), Num_Nodes_Per_Elem(nullptr),
        Num_Attr_Per_Elem(nullptr), Elem_Blk_Ids(nullptr), Elem_Blk_Types(nullptr),
        Elem_Blk_Names(nullptr), Node_Set_Names(nullptr), Side_Set_Names(nullptr),
        Elem_Blk_Attr_Names(nullptr), GM_Elem_Types(nullptr), Proc_Ids(nullptr)
  {
    Coord_Name[0] = Coord_Name[1] = Coord_Name[2] = nullptr;
    Proc_Info[0] = Proc_Info[1] = Proc_Info[2] = Proc_Info[3] = Proc_Info[4] = Proc_Info[5] = 0;
  }

  ~NemSpread() { safe_free((void **)&GM_Elem_Types); }
};

#endif
