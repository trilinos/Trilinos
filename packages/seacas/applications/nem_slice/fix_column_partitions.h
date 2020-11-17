/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef _ELB_FIX_COLUMN_PARTITIONS_H
#define _ELB_FIX_COLUMN_PARTITIONS_H

#include "elb.h" // for LB_Description<INT>, etc

/*! @brief If the mesh is columnar, ensure each column is fully in on partition
  @param lb    Load balancing or partitioning information (may be modified)
  @param mesh  Description of the mesh
  @param graph Description of the adjacency graph

  **** ASSUMES COLUMNS ARE STRICTLY VERTICAL, i.e., NO LATERAL FACE ****
  **** HAS A Z-COMPONENT IN ITS NORMAL ****

*/

template <typename INT>
int fix_column_partitions(LB_Description<INT> *lb, Mesh_Description<INT> const *const mesh,
                          Graph_Description<INT> const *const graph);

#endif // _ELB_FIX_COLUMN_PARTITIONS_H
