/*
 * Copyright(C) 1999-2020, 2022, 2023, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include "elb_elem.h"
#include "vector_data.h"
#include <cstdio>
#include <exodusII.h>
#include <string>
#include <vector>

#define ELB_VERSION "5.04 (2024/08/19)"
#define UTIL_NAME   "nem_slice"
#define ELB_FALSE   0
#define ELB_TRUE    1

/* Macro for maximum value */
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

/*
 * Constants for memory allocation of graph structures. The smaller these
 * values, the more memory-efficient the code will be. Larger values
 * will likely speed execution and prevent swap thrashing.
 */
#define SURND_ALLOC    8
#define ADJ_ALLOC      8
#define MEM_CHUNK_SIZE 16 /* Value MUST be >= 2 */
#define MEM_GROWTH     1.5

#define MAX_INP_LINE 10240

#if (__cplusplus >= 201703L)
#define FALL_THROUGH [[fallthrough]]
#elif defined(__GNUC__) && __GNUC__ >= 7 && !__INTEL_COMPILER
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
  std::vector<std::vector<INT>> int_nodes{};
  std::vector<std::vector<INT>> bor_nodes{};
  std::vector<std::vector<INT>> ext_nodes{};
  std::vector<std::vector<INT>> ext_procs{};

  /* Elemental */
  std::vector<std::vector<std::vector<INT>>> born_procs{};
  std::vector<std::vector<INT>>              int_elems{};
  std::vector<std::vector<INT>>              bor_elems{};
  std::vector<std::vector<INT>>              e_cmap_elems{};
  std::vector<std::vector<INT>>              e_cmap_sides{};
  std::vector<std::vector<INT>>              e_cmap_procs{};
  std::vector<std::vector<INT>>              e_cmap_neigh{};

  LB_Description() = default;
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
                              /* determine adjacencies */
  int   global_mech{-1};      /* true if looking for mechanisms in original mesh   */
  int   local_mech{-1};       /* true if looking for mechanisms in subdivided mesh */
  int   find_cnt_domains{-1}; /* true if finding number of connected domains in a graph */
  int   mech_add_procs{-1};   /* adds processors in cases of mechanisms       */
  int   dsd_add_procs{-1};    /* adds processors in cases of disconnected subdomains */
  int   no_sph{-1};
  int   fix_columns{0}; /* detect, fix vertical column partitioning */
  char *groups{nullptr};
  std::vector<int> group_no{};
  int              num_groups{-1};
  int              int64db{0};  /* integer types for output mesh database */
  int              int64api{0}; /* integer types for exodus api calls */

  Problem_Description() = default;
};

/* Structure for parameters needed for the Eigensolver in Chaco */
struct Solver_Description
{
  double tolerance{-1.0};
  int    rqi_flag{-1};
  int    vmax{-1};

  Solver_Description() = default;
};

/* Structure used to store information about the weighting scheme, if
 * any, that is to be used. */
struct Weight_Description
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
  std::vector<int> elemblk{};     /* Id of element block */
  std::vector<int> elemblk_wgt{}; /* Weight of that element block */

  /* vector to indicate if weight value has already been overwritten */
  std::vector<int> ow{};

  std::vector<int>   vertices{};
  std::vector<float> edges{};

  Weight_Description() = default;
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
  std::vector<float>  coords{};
  std::vector<E_Type> elem_type{};
  INT               **connect;

  Mesh_Description() : connect(nullptr) {}
};

/* Structure for handling meshes with spheres */
struct Sphere_Info
{
  size_t           num{0};
  std::vector<int> adjust{};
  std::vector<int> begin{};
  std::vector<int> end{};

  Sphere_Info() = default;
};

/* Structure used to store various information about the graph */
template <typename INT> struct Graph_Description
{
  size_t                        nadj{0};
  int                           max_nsur{0};
  std::vector<INT>              adj{};
  std::vector<INT>              start{};
  std::vector<std::vector<INT>> sur_elem;
  Graph_Description() = default;
};

/* Various constants */
enum DecompType { NODAL, ELEMENTAL };

#define UTIL_NAME "nem_slice"

/* Load balance types */
enum Balance {
  MULTIKL,
  SPECTRAL,
  INERTIAL,
  LINEAR,
  RANDOM,
  SCATTERED,
  INFILE,
  KL_REFINE,
  NO_REFINE,
  NUM_SECTS,
  CNCT_DOM,
  OUTFILE,
  ZPINCH,
  BRICK,
  ZOLTAN_RCB,
  ZOLTAN_RIB,
  ZOLTAN_HSFC,
  IGNORE_Z
};

/* Machine types */
enum MachineType { MESH, HCUBE, HYPERCUBE, CLUSTER };

/* Solver options */
enum SolverOptions { TOLER, USE_RQI, VMAX };

/* ISSUES options */
enum Issues { LOCAL_ISSUES, GLOBAL_ISSUES };

/* Weighting options */
/*
 * NOTE: the position of NO_WEIGHT, READ_EXO, EL_BLK, and EWGT_ON
 * should not be changed. These are the possible values for the
 * "type" variable in the Weight struct. They need to b 0, 1, 2, & 4
 * to allow bit masking for the type. The other variables are not
 * currently used in the type, but are needed since they appear
 * on the command line.
 */
enum WeightingOptions { NO_WEIGHT, READ_EXO, EL_BLK, VAR_INDX, EDGE_WGT, TIME_INDX, VAR_NAME };
