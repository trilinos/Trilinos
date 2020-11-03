// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#ifndef EJOIN_MAPPING_H
#define EJOIN_MAPPING_H

#include "EJ_CodeTypes.h"

template <typename INT>
void eliminate_omitted_nodes(RegionVector &part_mesh, std::vector<INT> &global_node_map,
                             std::vector<INT> &local_node_map, bool fill_global);

template <typename INT>
void build_reverse_node_map(Ioss::Region &output_region, RegionVector &part_mesh,
                            std::vector<INT> &global_node_map, std::vector<INT> &local_node_map);

template <typename INT>
void build_local_element_map(RegionVector &part_mesh, std::vector<INT> &local_element_map);

template <typename INT>
void generate_element_ids(RegionVector &part_mesh, const std::vector<INT> &local_element_map,
                          std::vector<INT> &global_element_map);
#endif
