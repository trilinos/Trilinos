// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "EJ_mapping.h"
#include "Ioss_ElementBlock.h"   // for ElementBlock
#include "Ioss_GroupingEntity.h" // for GroupingEntity
#include "Ioss_NodeBlock.h"      // for NodeBlock
#include "Ioss_Property.h"       // for Property
#include "Ioss_Region.h"         // for Region, etc
#include "Ioss_SmartAssert.h"
#include <algorithm> // for sort, unique
#include <cstddef>   // for size_t
#include <fmt/ostream.h>
#include <numeric>
#include <utility> // for make_pair, pair

namespace {
  bool entity_is_omitted(Ioss::GroupingEntity *block)
  {
    bool omitted = false;
    if (block->property_exists("omitted")) {
      omitted = (block->get_property("omitted").get_int() == 1);
    }
    return omitted;
  }
} // namespace

template <typename INT>
void eliminate_omitted_nodes(RegionVector &part_mesh, std::vector<INT> &global_node_map,
                             std::vector<INT> &local_node_map, bool fill_global)
{
  size_t offset     = 0;
  size_t j          = 0;
  size_t part_count = part_mesh.size();
  for (size_t p = 0; p < part_count; p++) {
    bool has_omissions        = part_mesh[p]->get_property("block_omission_count").get_int() > 0;
    Ioss::NodeBlock *nb       = part_mesh[p]->get_node_blocks()[0];
    size_t           loc_size = nb->entity_count();
    if (has_omissions) {
      // If there are any omitted element blocks for this part, don't
      // map the nodes that are only connected to omitted element
      // blocks.
      std::vector<char> node_status;
      nb->get_field_data("node_connectivity_status", node_status);
      for (size_t i = 0; i < node_status.size(); i++) {
        if (node_status[i] != 1) {
          local_node_map[offset + i] = j;
          if (fill_global) {
            global_node_map.push_back(j + 1);
          }
          j++;
        }
        else {
          local_node_map[offset + i] = -1;
        }
      }
    }
    else {
      for (size_t i = 0; i < loc_size; i++) {
        local_node_map[offset + i] = j;
        if (fill_global) {
          global_node_map.push_back(j + 1);
        }
        j++;
      }
    }
    offset += loc_size;
  }
}

template void eliminate_omitted_nodes(RegionVector &part_mesh, std::vector<int> &global_node_map,
                                      std::vector<int> &local_node_map, bool fill_global);
template void eliminate_omitted_nodes(RegionVector &        part_mesh,
                                      std::vector<int64_t> &global_node_map,
                                      std::vector<int64_t> &local_node_map, bool fill_global);

template <typename INT>
void build_reverse_node_map(Ioss::Region & /*global*/, RegionVector &part_mesh,
                            std::vector<INT> &global_node_map, std::vector<INT> &local_node_map)
{
  // Instead of using <set> and <map>, consider using a sorted vector...
  // Append all local node maps to the global node map.
  // Sort the global node map
  // Remove duplicates.
  // Position within map is now the map...
  // When building the local-part node to global id, use binary_search...

  size_t part_count = part_mesh.size();

  // Global node map and count.
  std::vector<std::vector<int>> global_nodes(part_count);

  size_t tot_size = 0;
  for (size_t p = 0; p < part_count; p++) {
    Ioss::NodeBlock *nb       = part_mesh[p]->get_node_blocks()[0];
    size_t           loc_size = nb->entity_count();
    tot_size += loc_size;
    global_nodes[p].resize(loc_size);
  }
  global_node_map.resize(tot_size);

  size_t offset            = 0;
  bool   any_omitted_nodes = false;
  for (size_t p = 0; p < part_count; p++) {
    Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
    nb->get_field_data("ids", global_nodes[p]);

    // If there are any omitted element blocks for this part, set
    // the global id of any nodes that are only connected to omitted
    // element blocks to 0.
    bool has_omissions = part_mesh[p]->get_property("block_omission_count").get_int() > 0;
    if (has_omissions) {
      std::vector<char> node_status;
      nb->get_field_data("node_connectivity_status", node_status);
      for (size_t i = 0; i < node_status.size(); i++) {
        if (node_status[i] == 1) {
          any_omitted_nodes  = true;
          global_nodes[p][i] = 0;
        }
      }
    }

    std::copy(global_nodes[p].begin(), global_nodes[p].end(), &global_node_map[offset]);
    offset += global_nodes[p].size();
  }

  // Now, sort the global_node_map array and remove duplicates...
  Ioss::Utils::uniquify(global_node_map);

  // If any omitted nodes, remove them from the global_node_map.
  // The id will be 0
  if (any_omitted_nodes) {
    auto pos = std::remove(global_node_map.begin(), global_node_map.end(), 0);
    global_node_map.erase(pos, global_node_map.end());
  }

  size_t output_node_count = global_node_map.size();

  // See whether the node numbers are contiguous.  If so, we can map
  // the nodes back to their original location. Since the nodes are
  // sorted and there are no duplicates, we just need to see if the id
  // at global_node_map.size() == global_node_map.size();
  size_t max_id = global_node_map[output_node_count - 1];

  bool is_contiguous = max_id == output_node_count;
  fmt::print("Node map {} contiguous.\n", (is_contiguous ? "is" : "is not"));

  // Create the map that maps from a local part node to the
  // global map. This combines the mapping local part node to
  // 'global id' and then 'global id' to global position. The
  // mapping is now a direct lookup instead of a lookup followed by
  // a reverse map.
  auto cur_pos = global_node_map.begin();
  for (size_t p = 0; p < part_count; p++) {
    size_t noffset    = part_mesh[p]->get_property("node_offset").get_int();
    size_t node_count = global_nodes[p].size();
    for (size_t i = 0; i < node_count; i++) {
      INT global_node = global_nodes[p][i];

      if (global_node > 0) {
        if (cur_pos == global_node_map.end() || *cur_pos != global_node) {
          auto iter = std::lower_bound(global_node_map.begin(), global_node_map.end(), global_node);
          if (iter == global_node_map.end()) {
            INT n = global_node;
            fmt::print("{:n}\n", n);
            SMART_ASSERT(iter != global_node_map.end());
          }
          cur_pos = iter;
        }
        size_t nodal_value          = cur_pos - global_node_map.begin();
        local_node_map[noffset + i] = nodal_value;
        ++cur_pos;
      }
      else {
        local_node_map[noffset + i] = -1;
      }
    }
  }

  // Update the nodal ids to give a unique, non-repeating set.  If contiguous, then
  // there is nothing to do.  If not contiguous, then need to determine if there are any
  // repeats (id reuse) and if so, generate a new id for the repeated uses.
  if (!is_contiguous) {
    bool repeat_found = false;
    INT  id_last      = global_node_map[0];
    for (size_t i = 1; i < output_node_count; i++) {
      if (global_node_map[i] == id_last) {
        global_node_map[i] = ++max_id;
        repeat_found       = true;
      }
      else {
        id_last = global_node_map[i];
      }
    }
    if (repeat_found) {
      fmt::print("Duplicate node ids were found. Their ids have been renumbered to remove "
                 "duplicates.\n");
    }
  }
}

template void build_reverse_node_map(Ioss::Region &global, RegionVector &part_mesh,
                                     std::vector<int> &global_node_map,
                                     std::vector<int> &local_node_map);
template void build_reverse_node_map(Ioss::Region &global, RegionVector &part_mesh,
                                     std::vector<int64_t> &global_node_map,
                                     std::vector<int64_t> &local_node_map);

template <typename INT>
void build_local_element_map(RegionVector &part_mesh, std::vector<INT> &local_element_map)
{
  size_t global = 0;
  size_t offset = 0;
  for (auto &p : part_mesh) {

    const Ioss::ElementBlockContainer &         ebs = p->get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator i   = ebs.begin();

    while (i != ebs.end()) {
      Ioss::ElementBlock *eb       = *i++;
      size_t              num_elem = eb->entity_count();
      if (entity_is_omitted(eb)) {
        // Fill local_element_map with -1 for the omitted elements.
        for (size_t j = 0; j < num_elem; j++) {
          local_element_map[offset + j] = -1;
        }
      }
      else {
        for (size_t j = 0; j < num_elem; j++) {
          local_element_map[offset + j] = global++;
        }
      }
      offset += num_elem;
    }
  }
}

template void build_local_element_map(RegionVector &part_mesh, std::vector<int> &local_element_map);
template void build_local_element_map(RegionVector &        part_mesh,
                                      std::vector<int64_t> &local_element_map);

template <typename INT>
void generate_element_ids(RegionVector &part_mesh, const std::vector<INT> &local_element_map,
                          std::vector<INT> &global_element_map)
{
  // Follow same logic as 'build_local_element_map' to ensure elements
  // are processed in same order.

  // Many models do not use the element number map at all, so they
  // will have a 1..numel map.  If all parts have that, then we don't
  // want to do any fancy duplicate removal and other processing, just
  // output a 1..numel map for the output mesh...  We still generate
  // the global_element_map, but check whether any of the part blocks
  // have a non-1..numel map...
  bool   has_map = false;
  size_t offset  = 0;
  for (auto &p : part_mesh) {
    const Ioss::ElementBlockContainer &         ebs = p->get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator i   = ebs.begin();

    while (i != ebs.end()) {
      Ioss::ElementBlock *eb       = *i++;
      INT                 num_elem = eb->entity_count();
      if (!entity_is_omitted(eb)) {
        std::vector<INT> part_ids;
        eb->get_field_data("ids", part_ids);

        if (!has_map) {
          INT eb_offset = eb->get_offset();
          for (INT j = 0; j < num_elem; j++) {
            if (part_ids[j] != eb_offset + j + 1) {
              has_map = true;
              break;
            }
          }
        }

        for (INT j = 0; j < num_elem; j++) {
          INT gpos = local_element_map[offset + j];
          if (gpos >= 0) {
            global_element_map[gpos] = part_ids[j];
          }
        }
      }
      offset += num_elem;
    }
  }
  // Check for duplicates...
  // NOTE: Used to use an indexed sort here, but if there was a
  // duplicate id, it didn't really care whether part 1 or part N's
  // index came first which causes really screwy element maps.
  // Instead, lets sort a vector containing pairs of <id, index> where
  // the index will always? increase for increasing part numbers...
  if (has_map) {
    std::vector<std::pair<INT, INT>> index(global_element_map.size());
    for (size_t i = 0; i < index.size(); i++) {
      index[i] = std::make_pair(global_element_map[i], (INT)i);
    }

    std::sort(index.begin(), index.end());

    INT max_id = index[index.size() - 1].first + 1;

    size_t beg = 0;
    for (size_t i = 1; i < index.size(); i++) {
      if (index[beg].first == index[i].first) {
        // Duplicate found... Assign it a new id greater than any
        // existing id...  (What happens if we exceed INT_MAX?)
        global_element_map[index[i].second] = max_id++;
        // Keep 'beg' the same in case multiple duplicate of this value.
      }
      else {
        beg = i;
      }
    }
  }
  else {
    INT one = 1;
    std::iota(global_element_map.begin(), global_element_map.end(), one);
  }
}

template void generate_element_ids(RegionVector &          part_mesh,
                                   const std::vector<int> &local_element_map,
                                   std::vector<int> &      global_element_map);
template void generate_element_ids(RegionVector &              part_mesh,
                                   const std::vector<int64_t> &local_element_map,
                                   std::vector<int64_t> &      global_element_map);
