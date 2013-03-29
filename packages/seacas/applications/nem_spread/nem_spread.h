#ifndef NEM_SPREAD_H
#define NEM_SPREAD_H

#include <rf_allo.h>
#include "rf_io_const.h"
#include "globals.h"

#define UTIL_NAME "nem_spread"
#define VER_STR   "6.10 (2013/03/29)"

extern void check_exodus_error (int, const char *);
extern double second               (void);

template <typename T, typename INT>
  class NemSpread {
 public:

  void load_mesh();
  int check_inp();
  void read_coord(int mesh_exoid, int max_name_length);
  void read_elem_blk_ids(int mesh_exoid, int max_name_length);
  void read_elem_blk(int mesh_exoid);
  void extract_elem_blk();
  void extract_global_element_ids(INT global_ids[], size_t Num_Elem, int iproc);
  void extract_global_node_ids(INT global_ids[], size_t Num_Node, int iproc);
  size_t extract_elem_connect(INT elem_blk[], int icurrent_elem_blk,
			   size_t istart_elem, size_t iend_elem,
			   int *local_ielem_blk, int indx);
  void extract_elem_attr(T *elem_attr, int icurrent_elem_blk,
			 size_t istart_elem, size_t iend_elem,
			 int natt_p_elem, int indx);
  void find_elem_block(INT *proc_elem_blk, int indx, int proc_for);
  void read_node_set_ids(int mesh_exoid, INT [], INT [], int max_name_length);
  void read_side_set_ids(int mesh_exoid, INT [], INT [], int max_name_length);
  void read_node_sets(int mesh_exoid, INT *, INT *);
  void read_side_sets(int mesh_exoid, INT *, INT *);

  void read_nodal_vars(int mesh_exoid);
  void read_aux_nodal_vars(int mesh_exoid);
  void read_exo_nv(int exoid, int itotal_nodes, int n_glob_nodes, int proc_id,
		   int num_proc, int n_gnv_to_read, INT *gnode, int time_index,
		   INT *gnv_index, T **glob_var_arr);

  void load_lb_info();
  
  int write_var_param(int mesh_exoid, int max_name_length, 
		      int num_glob, char **gv_names,
		      int num_node, char **nv_names, 
		      int num_elem, char **ev_names, int *local_ebtt,
		      int num_nset, char **ns_names, int *local_nstt,
		      int num_sset, char **ss_names, int *local_sstt);

  void write_parExo_data(int mesh_exoid, int max_name_length,
			 int iproc,
			 INT *Num_Nodes_In_NS,
			 INT *Num_Elems_In_SS,
			 INT *Num_Elems_In_EB);

  void write_var_timestep(int exoid, int proc, int time_step, 
			  INT *eb_ids_global, INT *ss_ids_global, INT *ns_ids_global);
  
  void comm_lb_data (int iproc, int proc_for,
		     INT *Integer_Vector_Ptr, int read_length);
  void recv_lb_data (INT *Integer_Vector, int indx);
  void process_lb_data (INT *Integer_Vector, int index);
  void read_proc_init(int lb_exoid, int proc_info[], int **proc_ids_ptr);
  void read_lb_init (int exoid, INT *Int_Space, 
		     INT *Int_Node_Num, INT *Bor_Node_Num,
		     INT *Ext_Node_Num, INT *Int_Elem_Num,
		     INT *Bor_Elem_Num, INT *Node_Comm_Num,
		     INT *Elem_Comm_Num, char *Title);

  void read_cmap_params(int exoid, INT *Node_Comm_Num,
			INT *Elem_Comm_Num, INT *Num_N_Comm_Maps,
			INT *Num_E_Comm_Maps, ELEM_COMM_MAP<INT> **E_Comm_Map,
			NODE_COMM_MAP<INT> **N_Comm_Map, INT *cmap_max_size,
			INT **comm_vec);

  void create_elem_types();
  void read_mesh_param ();
  void read_restart_params();
  void read_restart_data ();

  int read_var_param (int exoid, int max_name_length);

  int broadcast_var_param (Restart_Description<T> *restart, int max_name_length);
  int read_vars(int exoid, int index, INT *eb_ids,
		INT *eb_cnts, INT ***eb_map_ptr, INT **eb_cnts_local,
		INT *ss_ids, INT *ss_cnts, INT *ns_ids, INT *ns_cnts);
  int read_elem_vars(int exoid, int index, INT *eb_ids,
		     INT *eb_cnts, INT ***eb_map_ptr,
		     INT **eb_cnts_local);
  int read_elem_vars_1(int exoid, int index, INT *eb_ids,
		       INT *eb_cnts, INT ***eb_map_ptr,
		       INT **eb_cnts_local, int iblk,
		       int eb_offset, INT *local_offset);
  int read_sset_vars(int exoid, int index, INT *ss_ids, INT *ss_cnts);
  int read_sset_vars_1(int exoid, int index, INT *ss_ids,
		       INT *ss_cnts, int iblk);
  int read_nset_vars(int exoid, int index, INT *ns_ids,
		     INT *ns_cnts);
  int read_nset_vars_1(int exoid, int index, INT *ns_ids,
		       INT *ns_cnts, int iblk);
  int read_nodal_vars_1 (int exoid, int index);
  int read_nodal_vars (int exoid, int index);
  int compare_mesh_param (int exoid);

  int int64db;
  int int64api;
  bool force64db; /* Store all ints as 64-bit on output databases. */
  
  int io_ws;
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

  INT *Node_Set_Ids;       /* Vector of node set ids             *
			    * (ptr to vector of length, 	     *
			    *  Num_Node_Set)    		     */
  INT *Side_Set_Ids;       /* Side set ids  		     *
			    * (ptr to vector of length, 	     *
			    *  Num_Side_Set)      		     */
  INT *Num_Elem_In_Blk;    /* Number of elements in each element *
			    * block (ptr to vector of length,    *
			    *        Num_Elem_Blk)     	     */
  INT *Num_Nodes_Per_Elem; /* Number of nodes per element in each*
			    * elem block (ptr to vector of       *
			    * length, Num_Elem_Blk)              */
  INT *Num_Attr_Per_Elem;  /* Number of attributes per element in*
			    * each block (ptr to vector of       *
			    * length, Num_Elem_Blk)              */
  INT *Elem_Blk_Ids;       /* Element block id's of each element *
			    * block (ptr to vector of length,    *
			    *        Num_Elem_Blk)               */
  char **Elem_Blk_Types;   /* Element block types for each       *
			    * element block (ptr to vector of    *
			    * length, Num_Elem_Blk, of vectors of*
			    * char of length MAX_STR_LENGTH+1)   */
  char **Elem_Blk_Names;   /* Element block names for each       *
			    * element block (ptr to vector of    *
			    * length, Num_Elem_Blk, of vectors of*
			    * char of length MAX_STR_LENGTH+1)   */
  char **Node_Set_Names;   /* Nodeset names for each             *
			    * nodeset (ptr to vector of          *
			    * length, Num_Node_Set, of vectors of*
			    * char of length MAX_STR_LENGTH+1)   */
  char **Side_Set_Names;   /* Sideset names for each             *
			    * sideset (ptr to vector of          *
			    * length, Num_Side_Set, of vectors of*
			    * char of length MAX_STR_LENGTH+1)   */
  char ***Elem_Blk_Attr_Names;  /* Element block attribute names for each 
				 * attribute in each element block
				 * ragged array (size varies for each
				 * element block          	     */

  int *GM_Elem_Types;      /* This is a list of the element      *
			    * in the entire global mesh. It is   *
			    * stored on Processor 0 to           *
			    * facilitate the reading of the side *
			    * set distribution factors.          */

  char *Coord_Name[3];    /* The name(s) of the coordinate axes.   */

  int     Proc_Info[6];
  int    *Proc_Ids;

  NemSpread() :
    int64db(0), int64api(0), force64db(false), io_ws(0)
    {
      Coord_Name[0] = Coord_Name[1] = Coord_Name[2] = NULL;
      Proc_Info[0] = Proc_Info[1] = Proc_Info[2] = Proc_Info[3] = Proc_Info[4] = Proc_Info[5] = 0;
    }

    ~NemSpread()
      {
	safe_free((void**) &GM_Elem_Types);
      }
};

#endif
