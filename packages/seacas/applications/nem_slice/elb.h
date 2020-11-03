/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef _EXOIILB_CONST_H_
#define _EXOIILB_CONST_H_

#include "elb_elem.h"
#include <cstdio>
#include <exodusII.h>
#include <string>
#include <vector>

#define ELB_VERSION "4.19"
#define UTIL_NAME "nem_slice"
#define ELB_FALSE 0
#define ELB_TRUE 1

/* Macro for maximum value */
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

/*
 * Constants for memory allocation of graph structures. The smaller these
 * values, the more memory-efficient the code will be. Larger values
 * will likely speed execution and prevent swap thrashing.
 */
#define SURND_ALLOC 8
#define ADJ_ALLOC 8
#define MEM_CHUNK_SIZE 16 /* Value MUST be >= 2 */
#define MEM_GROWTH 1.5

#define MAX_INP_LINE 10240

#if defined(__GNUC__) && __GNUC__ >= 7 && !__INTEL_COMPILER
#define FALL_THROUGH [[gnu::fallthrough]]
#else
#define FALL_THROUGH ((void)0)
#endif /* __GNUC__ >= 7 */

template <typename INT> void vec_free(std::vector<INT> &V)
{
  V.clear();
  V.shrink_to_fit();
}

/* Prototype for timing function */
extern double get_time();

/* Structure used for the description of the machine for which the
 * load balance is to be constructed. */
struct Machine_Description
{
  int type{-1};
  int num_dims{-1};
  int dim[3]{};
  int num_boxes{-1};     /* added for cluster type machines */
  int procs_per_box{-1}; /* added for cluster type machines, if there is only
                        one box, then this is the same as num_procs */
  int num_procs{-1};

  Machine_Description() { dim[0] = dim[1] = dim[2] = -1; }
};

/* Structure used for the description of what type of load balance is
 * to be performed. */
template <typename INT> struct LB_Description
{
  int         type{-1};
  int         ignore_z{0};
  int         refine{-1};
  int         num_sects{-1};
  int         cnctd_dom{-1};
  int         outfile{-1};
  std::string file{};

  /* Calculated quantities */
  int *vertex2proc{nullptr};

  /* Nodal */
  std::vector<std::vector<INT>> int_nodes;
  std::vector<std::vector<INT>> bor_nodes;
  std::vector<std::vector<INT>> ext_nodes;
  std::vector<std::vector<INT>> ext_procs;

  /* Elemental */
  std::vector<std::vector<std::vector<INT>>> born_procs{};
  std::vector<std::vector<INT>>              int_elems;
  std::vector<std::vector<INT>>              bor_elems;
  std::vector<std::vector<INT>>              e_cmap_elems;
  std::vector<std::vector<INT>>              e_cmap_sides;
  std::vector<std::vector<INT>>              e_cmap_procs;
  std::vector<std::vector<INT>>              e_cmap_neigh;

  LB_Description() {}
};

/* Structure for the problem description. */
struct Problem_Description
{
  int    type{-1};
  int    read_coords{-1};
  int    coarse_flag{-1};
  int    alloc_graph{-1};
  size_t num_vertices{0};
  int    vis_out{-1};
  int    skip_checks{-1};     /* put in to skip some error checks for some meshes  */
  int    face_adj{-1};        /* true if using face definition of adjacencies      */
  int    partial_adj{0};      /* true if allowing partial (3/4) of nodes to */
                              /* determine adjancencies */
  int   global_mech{-1};      /* true if looking for mechanisms in original mesh   */
  int   local_mech{-1};       /* true if looking for mechanisms in subdivided mesh */
  int   find_cnt_domains{-1}; /* true if finding number of connected domains in a graph */
  int   mech_add_procs{-1};   /* adds processors in cases of mechanisms       */
  int   dsd_add_procs{-1};    /* adds processors in cases of disconnected subdomains */
  int   no_sph{-1};
  int   fix_columns{0}; /* detect, fix vertical column partitioning */
  char *groups{nullptr};
  int * group_no{nullptr};
  int   num_groups{-1};
  int   int64db{0};  /* integer types for output mesh database */
  int   int64api{0}; /* integer types for exodus api calls */

  Problem_Description() {}
};

/* Structure for parameters needed for the Eigensolver in Chaco */
struct Solver_Description
{
  double tolerance{-1.0};
  int    rqi_flag{-1};
  int    vmax{-1};

  Solver_Description() {}
};

/* Structure used to store information about the weighting scheme, if
 * any, that is to be used. */
template <typename INT> struct Weight_Description
{
  int type{-1};   /* See weight type below for possible types */
  int ow_read{0}; /* 1 if element block settings overwrite exodus file read */

  std::string exo_filename{};
  std::string exo_varname{};

  int exo_tindx{-1};
  int exo_vindx{-1};

  /* Variable parameters */
  int nvals{0};

  /* vectors to hold element block weights */
  std::vector<INT> elemblk{};     /* Id of element block */
  std::vector<INT> elemblk_wgt{}; /* Weight of that element block */

  /* vector to indicate if weight value has already been overwritten */
  std::vector<INT> ow{};

  std::vector<int>   vertices{};
  std::vector<float> edges{};

  Weight_Description<INT>() {}
};

/* Structure used to store information about the FEM mesh */
template <typename INT> struct Mesh_Description
{
  size_t              num_nodes{0};
  size_t              num_elems{0};
  size_t              num_dims{0};
  size_t              num_el_blks{0};
  std::vector<INT>    eb_cnts{};
  std::vector<INT>    eb_ids{};
  std::vector<INT>    eb_npe{};
  std::vector<E_Type> eb_type{};
  size_t              num_node_sets{0};
  size_t              num_side_sets{0};
  size_t              max_np_elem{0};
  size_t              ns_list_len{0};
  char                title[MAX_LINE_LENGTH + 1]{};
  float *             coords{nullptr};
  E_Type *            elem_type{nullptr};
  INT **              connect;

  Mesh_Description() : connect(nullptr) {}
};

/* Structure for handling meshes with spheres */
struct Sphere_Info
{
  size_t num{0};
  int *  adjust{nullptr};
  int *  begin{nullptr};
  int *  end{nullptr};

  Sphere_Info() {}
};

/* Structure used to store various information about the graph */
template <typename INT> struct Graph_Description
{
  size_t                        nadj{0};
  int                           max_nsur{0};
  std::vector<INT>              adj{};
  std::vector<INT>              start{};
  std::vector<std::vector<INT>> sur_elem;
  Graph_Description<INT>() {}
};

/* Various constants */
#define NODAL 0
#define ELEMENTAL 1

#define UTIL_NAME "nem_slice"

/* Load balance types */
#define MULTIKL 0
#define SPECTRAL 1
#define INERTIAL 2
#define LINEAR 3
#define RANDOM 4
#define SCATTERED 5
#define INFILE 6
#define KL_REFINE 7
#define NO_REFINE 8
#define NUM_SECTS 9
#define CNCT_DOM 10
#define OUTFILE 11
#define ZPINCH 12
#define BRICK 13
#define ZOLTAN_RCB 14
#define ZOLTAN_RIB 15
#define ZOLTAN_HSFC 16
#define IGNORE_Z 17

/* Machine types */
#define MESH 0
#define HCUBE 1
#define HYPERCUBE 2
#define CLUSTER 3

/* Solver options */
#define TOLER 0
#define USE_RQI 1
#define VMAX 2

/* ISSUES options */

#define LOCAL_ISSUES 0
#define GLOBAL_ISSUES 1

/* Weighting options */
/*
 * NOTE: the position of NO_WEIGHT, READ_EXO, EL_BLK, and EWGT_ON
 * should not be changed. These are the possible values for the
 * "type" variable in the Weight struct. They need to b 0, 1, 2, & 4
 * to allow bit masking for the type. The other variables are not
 * currently used in the type, but are needed since they appear
 * on the command line.
 */
#define NO_WEIGHT 0
#define READ_EXO 1
#define EL_BLK 2
#define VAR_INDX 3
#define EDGE_WGT 4
#define TIME_INDX 5
#define VAR_NAME 6

#endif /* _EXOIILB_CONST_H_ */
