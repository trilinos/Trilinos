// Copyright(C) 2021, 2022, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifdef USE_ZOLTAN
#include "zoltan.h"       // for Zoltan_Set_Param, etc
#include "zoltan_types.h" // for ZOLTAN_ID_PTR*, ZOLTAN_OK, etc
#include <mpi.h>          // for MPI_Finalize, etc
#endif

#include <algorithm>
#include <fmt/format.h>
#include <random>
#include <vector>

#include "Cell.h"
#include "Grid.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Sort.h"

extern unsigned int debug_level;

namespace {
#ifdef USE_ZOLTAN
#define STRINGIFY(x) #x
#define TOSTRING(x)  STRINGIFY(x)
#define ZCHECK(funcall)                                                                            \
  do {                                                                                             \
    ierr = (funcall);                                                                              \
    if (ierr == ZOLTAN_FATAL) {                                                                    \
      fmt::print(stderr, "Error returned from {} ({}:{})\n", TOSTRING(funcall), __FILE__,          \
                 __LINE__);                                                                        \
      goto End;                                                                                    \
    }                                                                                              \
  } while (0)

  /*****************************************************************************/
  /***** Global data structure used by Zoltan callbacks.                   *****/
  /***** Could implement Zoltan callbacks without global data structure,   *****/
  /***** but using the global data structure makes implementation quick.   *****/
  struct
  {
    size_t ndot; /* Length of x, y, z, and part (== # of elements) */
    int   *vwgt; /* vertex weights */
    float *x;    /* x-coordinates */
    float *y;    /* y-coordinates */
    float *z;    /* z-coordinates */
  } Zoltan_Data;

  /*****************************************************************************/
  /***** ZOLTAN CALLBACK FUNCTIONS *****/
  int zoltan_num_dim(void * /*data*/, int *ierr)
  {
    /* Return dimensionality of coordinate data.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */
    *ierr = ZOLTAN_OK;
    if (Zoltan_Data.z != nullptr) {
      return 3;
    }
    if (Zoltan_Data.y != nullptr) {
      return 2;
    }
    return 1;
  }

  int zoltan_num_obj(void * /*data*/, int *ierr)
  {
    /* Return number of objects.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */
    *ierr = ZOLTAN_OK;
    return Zoltan_Data.ndot;
  }

  void zoltan_obj_list(void * /*data*/, int /*ngid_ent*/, int /*nlid_ent*/, ZOLTAN_ID_PTR gids,
                       ZOLTAN_ID_PTR /*lids*/, int wdim, float *wgts, int *ierr)
  {
    /* Return list of object IDs.
     * Return only global IDs; don't need local IDs since running in serial.
     * gids are array indices for coordinate and vwgts arrays.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */
    for (size_t i = 0; i < Zoltan_Data.ndot; i++) {
      gids[i] = i;
      if (wdim != 0) {
        wgts[i] = static_cast<float>(Zoltan_Data.vwgt[i]);
      }
    }

    *ierr = ZOLTAN_OK;
  }

  void zoltan_geom(void * /*data*/, int /*ngid_ent*/, int /*nlid_ent*/, int nobj,
                   const ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR /*lids*/, int ndim, double *geom,
                   int *ierr)
  {
    /* Return coordinates for objects.
     * gids are array indices for coordinate arrays.
     * Using global data structure Zoltan_Data, initialized in ZOLTAN_RCB_assign.
     */

    for (int i = 0; i < nobj; i++) {
      size_t j       = gids[i];
      geom[i * ndim] = Zoltan_Data.x[j];
      if (ndim > 1) {
        geom[i * ndim + 1] = Zoltan_Data.y[j];
      }
      if (ndim > 2) {
        geom[i * ndim + 2] = Zoltan_Data.z[j];
      }
    }

    *ierr = ZOLTAN_OK;
  }
#endif
} // namespace

void decompose_grid(Grid &grid, int ranks, const std::string &method)
{
  if (ranks == 1) {
    return;
  }

  if (method == "CYCLIC" || method == "LINEAR" || method == "RANDOM") {
    std::vector<int> rank_vec(grid.size());

    // Generate a "cyclic" vector 0, 1, 2, ... ranks-1, 0, 1, ...
    std::generate(rank_vec.begin(), rank_vec.end(),
                  [n = 0, ranks]() mutable { return n++ % ranks; });

    if (method == "LINEAR") {
      Ioss::sort(rank_vec.begin(), rank_vec.end());
    }
    else if (method == "RANDOM") {
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::shuffle(rank_vec.begin(), rank_vec.end(), std::default_random_engine(seed));
    }

    size_t k = 0;
    for (size_t j = 0; j < grid.JJ(); j++) {
      for (size_t i = 0; i < grid.II(); i++) {
        auto &cell = grid.get_cell(i, j);
        cell.set_rank(Loc::C, rank_vec[k++]);
      }
    }
    return;
  }

#ifdef USE_ZOLTAN
  // Below here are Zoltan decompositions...
  std::vector<float> x(grid.size());
  std::vector<float> y(grid.size());
  std::vector<int>   w(grid.size());

  // Get coordinates of each cell centroid... (actually origin of cell)
  size_t idx = 0;
  for (size_t j = 0; j < grid.JJ(); j++) {
    for (size_t i = 0; i < grid.II(); i++) {
      const auto &cell = grid.get_cell(i, j);
      x[idx]           = cell.m_offX;
      y[idx]           = cell.m_offY;

      const auto &element_blocks = cell.region()->get_element_blocks();
      for (const auto *block : element_blocks) {
        w[idx] += block->entity_count();
      }
      idx++;
    }
  }

  /* Copy mesh data and pointers into structure accessible from callback fns. */
  Zoltan_Data.ndot = grid.size();
  Zoltan_Data.vwgt = Data(w);
  Zoltan_Data.x    = Data(x);
  Zoltan_Data.y    = Data(y);
  Zoltan_Data.z    = nullptr;

  /* Initialize Zoltan */
  int           argc   = 0;
  char        **argv   = nullptr;
  ZOLTAN_ID_PTR zgids  = nullptr;
  ZOLTAN_ID_PTR zlids  = nullptr; /* Useful output from Zoltan_LB_Partition */
  int          *zprocs = nullptr; /* Useful output from Zoltan_LB_Partition */
  int          *zparts = nullptr; /* Useful output from Zoltan_LB_Partition */

  int   ierr = 0;
  float ver  = 0.0;
  Zoltan_Initialize(argc, argv, &ver);
  struct Zoltan_Struct *zz = Zoltan_Create(Ioss::ParallelUtils::comm_world());

  /* Register Callback functions */
  /* Using global Zoltan_Data; could register it here instead as data field. */
  ZCHECK(Zoltan_Set_Fn(zz, ZOLTAN_NUM_GEOM_FN_TYPE,
                       reinterpret_cast<ZOLTAN_VOID_FN *>(zoltan_num_dim), nullptr));
  ZCHECK(Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE,
                       reinterpret_cast<ZOLTAN_VOID_FN *>(zoltan_num_obj), nullptr));
  ZCHECK(Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE,
                       reinterpret_cast<ZOLTAN_VOID_FN *>(zoltan_obj_list), nullptr));
  ZCHECK(Zoltan_Set_Fn(zz, ZOLTAN_GEOM_MULTI_FN_TYPE,
                       reinterpret_cast<ZOLTAN_VOID_FN *>(zoltan_geom), nullptr));

  /* Set parameters for Zoltan */
  ZCHECK(Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0"));
  {
    std::string str = fmt::format("{}", ranks);
    ZCHECK(Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTITIONS", str.c_str()));
    ZCHECK(Zoltan_Set_Param(zz, "LB_METHOD", method.c_str()));
  }
  ZCHECK(Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0"));
  ZCHECK(Zoltan_Set_Param(zz, "REMAP", "0"));
  ZCHECK(Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTITION_ASSIGNMENTS"));
  if (Zoltan_Data.vwgt != nullptr) {
    ZCHECK(Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"));
  }
  ZCHECK(Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1"));

  /* Call partitioner */
  {
    fmt::print(" Using Zoltan version {:.2}, method {}\n", static_cast<double>(ver), method);
    int           zngid_ent = 0;
    int           znlid_ent = 0; /* Useful output from Zoltan_LB_Partition */
    int           znobj     = 0;
    ZOLTAN_ID_PTR dummy1    = nullptr;
    ZOLTAN_ID_PTR dummy2    = nullptr; /* Empty output from Zoltan_LB_Partition */
    int           dummy0    = 0;
    int          *dummy3    = nullptr;
    int          *dummy4    = nullptr;
    int           changes   = 0;

    ZCHECK(Zoltan_LB_Partition(zz, &changes, &zngid_ent, &znlid_ent, &dummy0, &dummy1, &dummy2,
                               &dummy3, &dummy4, &znobj, &zgids, &zlids, &zprocs, &zparts));

    /* Sanity check */
    if (grid.size() != static_cast<size_t>(znobj)) {
      fmt::print(stderr, "Sanity check failed; ndot {} != znobj {}.\n", grid.size(),
                 static_cast<size_t>(znobj));
      goto End;
    }

    idx = 0;
    for (size_t j = 0; j < grid.JJ(); j++) {
      for (size_t i = 0; i < grid.II(); i++) {
        auto &cell = grid.get_cell(i, j);
        cell.set_rank(Loc::C, zparts[idx++]);
      }
    }
  }

End:
  /* Clean up */
  Zoltan_LB_Free_Part(&zgids, &zlids, &zprocs, &zparts);
  Zoltan_Destroy(&zz);
  if (ierr != 0) {
    MPI_Finalize();
    exit(-1);
  }
#endif
}
