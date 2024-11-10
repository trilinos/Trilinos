/*
 * Copyright(C) 1999-2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include <array>

#include "globals.h"
#include "pe_str_util_const.h" // for strip_string, token_compare, etc
#include "rf_allo.h"
#include "rf_io_const.h"

#define UTIL_NAME "nem_spread"
#define VER_STR   "7.04 (2024/02/28)"

extern void   check_exodus_error(int, const char *);
extern double second();

template <typename T, typename INT> class NemSpread
{
public:
  void   load_mesh();
  int    check_inp();
  void   read_coord(int exoid, int max_name_length);
  void   read_elem_blk_ids(int mesh_exoid, int max_name_length);
  void   read_elem_blk(int exoid);
  void   extract_elem_blk();
  void   extract_global_element_ids(const std::vector<INT> &global_ids, size_t Num_Elem, int iproc);
  void   extract_global_node_ids(const std::vector<INT> &global_ids, size_t Num_Node, int iproc);
  size_t extract_elem_connect(INT elem_blk[], int icurrent_elem_blk, size_t istart_elem,
                              size_t iend_elem, int *local_ielem_blk, int iproc);
  void extract_elem_attr(T *elem_attr, int icurrent_elem_blk, size_t istart_elem, size_t iend_elem,
                         int natt_p_elem, int iproc);
  void find_elem_block(std::vector<INT> &proc_elem_blk, int iproc, int proc_for);
  void read_node_set_ids(int mesh_exoid, std::vector<INT> &num_nodes_in_node_set,
                         std::vector<INT> &num_df_in_nsets, int max_name_length);
  void read_side_set_ids(int mesh_exoid, std::vector<INT> &num_elem_in_ssets,
                         std::vector<INT> &num_df_in_ssets, int max_name_length);
  void read_node_sets(int exoid, std::vector<INT> &num_nodes_in_node_set,
                      std::vector<INT> &num_df_in_nsets);
  void read_side_sets(int exoid, std::vector<INT> &num_elem_in_ssets,
                      std::vector<INT> &num_df_in_ssets);

  void read_nodal_vars(int mesh_exoid);

  void load_lb_info();

  int write_var_param(int mesh_exoid, int max_name_length, int num_glob, char **gv_names,
                      int num_node, char **nv_names, int num_elem, char **ev_names, int *local_ebtt,
                      int num_nset, char **ns_names, int *local_nstt, int num_sset, char **ss_names,
                      int *local_sstt);

  void write_parExo_data(int mesh_exoid, int max_name_length, int iproc,
                         std::vector<INT> &Num_Nodes_In_NS, std::vector<INT> &Num_Elems_In_SS,
                         std::vector<INT> &Num_Elems_In_EB);

  void write_var_timestep(int exoid, int proc, int time_step, INT *eb_ids_global,
                          INT *ss_ids_global, INT *ns_ids_global);

  void process_lb_data(INT *Integer_Vector, int indx);
  void read_proc_init(int lb_exoid, std::array<int, 6> &proc_info, std::vector<int> &proc_ids);
  void read_lb_init(int lb_exoid, std::vector<INT> &Int_Space, std::vector<INT> &Int_Node_Num,
                    std::vector<INT> &Bor_Node_Num, std::vector<INT> &Ext_Node_Num,
                    std::vector<INT> &Int_Elem_Num, std::vector<INT> &Bor_Elem_Num,
                    std::vector<INT> &Node_Comm_Num, std::vector<INT> &Elem_Comm_Num, char *Title);

  void read_cmap_params(int lb_exoid, std::vector<ELEM_COMM_MAP<INT>> &E_Comm_Map,
                        std::vector<NODE_COMM_MAP<INT>> &N_Comm_Map, INT *cmap_max_size);

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

  int  int64db{0};
  int  int64api{0};
  bool force64db{false}; /* Store all ints as 64-bit on output databases. */

  int                    io_ws{0};
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

  std::vector<INT> Node_Set_Ids;        /* Vector of node set ids             *
                                         * (ptr to vector of length,          *
                                         *  Num_Node_Set)                     */
  std::vector<INT> Side_Set_Ids;        /* Side set ids                       *
                                         * (ptr to vector of length,          *
                                         *  Num_Side_Set)                     */
  std::vector<INT> Num_Elem_In_Blk;     /* Number of elements in each element *
                                         * block (ptr to vector of length,    *
                                         *        Num_Elem_Blk)               */
  std::vector<INT> Num_Nodes_Per_Elem;  /* Number of nodes per element in each*
                                         * elem block (ptr to vector of       *
                                         * length, Num_Elem_Blk)              */
  std::vector<INT> Num_Attr_Per_Elem;   /* Number of attributes per element in*
                                         * each block (ptr to vector of       *
                                         * length, Num_Elem_Blk)              */
  std::vector<INT> Elem_Blk_Ids;        /* Element block id's of each element *
                                         * block (ptr to vector of length,    *
                                         *        Num_Elem_Blk)               */
  char **Elem_Blk_Types{nullptr};       /* Element block types for each       *
                                         * element block (ptr to vector of    *
                                         * length, Num_Elem_Blk, of vectors of*
                                         * char of length MAX_STR_LENGTH+1)   */
  char **Elem_Blk_Names{nullptr};       /* Element block names for each       *
                                         * element block (ptr to vector of    *
                                         * length, Num_Elem_Blk, of vectors of*
                                         * char of length MAX_STR_LENGTH+1)   */
  char **Node_Set_Names{nullptr};       /* Nodeset names for each             *
                                         * nodeset (ptr to vector of          *
                                         * length, Num_Node_Set, of vectors of*
                                         * char of length MAX_STR_LENGTH+1)   */
  char **Side_Set_Names{nullptr};       /* Sideset names for each             *
                                         * sideset (ptr to vector of          *
                                         * length, Num_Side_Set, of vectors of*
                                         * char of length MAX_STR_LENGTH+1)   */
  char ***Elem_Blk_Attr_Names{nullptr}; /* Element block attribute names for each
                                         * attribute in each element block
                                         * ragged array (size varies for each
                                         * element block                      */

  int *GM_Elem_Types{nullptr}; /* This is a list of the element      *
                                * in the entire global mesh. It is   *
                                * stored on Processor 0 to           *
                                * facilitate the reading of the side *
                                * set distribution factors.          */

  char *Coord_Name[3]{}; /* The name(s) of the coordinate axes.   */

  std::array<int, 6> Proc_Info{};
  std::vector<int>   Proc_Ids{};

  NemSpread() { Coord_Name[0] = Coord_Name[1] = Coord_Name[2] = nullptr; }

  ~NemSpread() { safe_free((void **)&GM_Elem_Types); }
};
