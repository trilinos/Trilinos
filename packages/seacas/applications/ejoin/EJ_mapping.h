// Copyright(C) 1999-2020, 2022, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include "EJ_CodeTypes.h"

using IO_map =
    std::map<Ioss::GroupingEntity *, std::vector<std::pair<Ioss::GroupingEntity *, size_t>>>;

template <typename INT>
void eliminate_omitted_nodes(const RegionVector &part_mesh, std::vector<INT> &global_node_map,
                             std::vector<INT> &local_node_map, bool fill_global);

template <typename INT>
void build_reverse_node_map(Ioss::Region &output_region, const RegionVector &part_mesh,
                            std::vector<INT> &global_node_map, std::vector<INT> &local_node_map);

template <typename INT>
std::vector<INT> build_local_element_map(const RegionVector &part_mesh, Ioss::Region &output_region,
                                         const IO_map &output_input_map);

template <typename INT>
void generate_element_ids(const RegionVector &part_mesh, const std::vector<INT> &local_element_map,
                          std::vector<INT> &global_element_map);

template <typename INT>
void select_nodeset_nodes(const RegionVector &part_mesh, std::vector<INT> &local_node_map,
                          const Omissions &nset_match);
