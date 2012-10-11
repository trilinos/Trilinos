// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Ioss_Map.h>
#include <Ioss_Utils.h>
#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {
  typedef std::vector<Ioss::IdPair>::const_iterator RMapI;
}

Ioss::Map::Map() :
  entityCount(-1), sequentialG2L(true), entityReordered(false)
{}


void Ioss::Map::build_reverse_map(ReverseMapContainer *Map, const int64_t *ids,
				 int64_t num_to_get, int64_t offset, int processor)
{
  // Stored as a sorted vector of <global_id, local_id> pairs...
  // To build incrementally:
  // 0. PRE: reverseElementMap is sorted, size >= 0.
  // 1. Build vector of current ids. -- new_ids
  // 2. Sort that vector.
  // 3. Copy reverseElementMap to old_ids, empty reverseElementMap.
  // 4. Merge old_ids and new_ids to reverseElementMap.
  // 5. Check for duplicate global_ids...

  // Build a vector containing the current ids...
  ReverseMapContainer new_ids(num_to_get);
  for (int64_t i=0; i < num_to_get; i++) {
    int64_t local_id = offset + i + 1;
    new_ids[i] = std::make_pair(ids[i], local_id);
  }

  // Sort that vector...
  std::sort(new_ids.begin(), new_ids.end(), IdPairCompare());

  int64_t new_id_min = new_ids.empty() ? 0 : new_ids.front().first;
  int64_t old_id_max = Map->empty() ? 0 : Map->back().first;
  if (new_id_min > old_id_max) {
    Map->insert(Map->end(), new_ids.begin(), new_ids.end());
  } else {
    // Copy reverseElementMap to old_ids, empty reverseElementMap.
    ReverseMapContainer old_ids;
    old_ids.swap(*Map);
    assert(Map->empty());
    
    // Merge old_ids and new_ids to reverseElementMap.
    Map->reserve(old_ids.size() + new_ids.size());
    std::merge(old_ids.begin(), old_ids.end(),
	       new_ids.begin(), new_ids.end(),
	       std::inserter(*Map, Map->begin()), IdPairCompare());
    
    // Check for duplicate ids...
    verify_no_duplicate_ids(*Map, processor);
  }

}

void Ioss::Map::build_reverse_map(int processor)
{
  if (map[0] == 1) {
    build_reverse_map(&reverse, &map[1], map.size()-1, 0, processor);
  }
}

void Ioss::Map::build_reorder_map(int64_t start, int64_t count)
{
  // Note: To further add confusion, the reorderEntityaMap is 0-based
  // and the reverseEntityMap and entityMap are 1-baed. This is
  // just a consequence of how they are intended to be used...
  //
  // start is based on a 0-based array -- start of the reorderMap to build.
      
  if (reorder.empty())
    reorder.resize(map.size()-1);
      
  int64_t my_end = start+count;
  for (int64_t i=start; i < my_end; i++) {
    int64_t global_id = map[i+1];
    int64_t orig_local_id = global_to_local(global_id) - 1;
	
    // If we assume that partial output is not being used (it
    // currently isn't in Sierra), then the reordering should only be
    // a permutation of the original ordering within this entity block...
    assert(orig_local_id >= start && orig_local_id <= my_end);
    reorder[i] = orig_local_id;
  }
}
    
void Ioss::Map::build_reverse_map(ReverseMapContainer *Map, const int *ids,
				 int num_to_get, int offset, int processor)
{
  // Stored as a sorted vector of <global_id, local_id> pairs...
  // To build incrementally:
  // 0. PRE: reverseElementMap is sorted, size >= 0.
  // 1. Build vector of current ids. -- new_ids
  // 2. Sort that vector.
  // 3. Copy reverseElementMap to old_ids, empty reverseElementMap.
  // 4. Merge old_ids and new_ids to reverseElementMap.
  // 5. Check for duplicate global_ids...

  // Build a vector containing the current ids...
  ReverseMapContainer new_ids(num_to_get);
  for (int i=0; i < num_to_get; i++) {
    int local_id = offset + i + 1;
    new_ids[i] = std::make_pair(ids[i], local_id);
  }

  // Sort that vector...
  std::sort(new_ids.begin(), new_ids.end(), IdPairCompare());

  int new_id_min = new_ids.empty() ? 0 : new_ids.front().first;
  int old_id_max = Map->empty() ? 0 : Map->back().first;
  if (new_id_min > old_id_max) {
    Map->insert(Map->end(), new_ids.begin(), new_ids.end());
  } else {
    // Copy reverseElementMap to old_ids, empty reverseElementMap.
    ReverseMapContainer old_ids;
    old_ids.swap(*Map);
    assert(Map->empty());
    
    // Merge old_ids and new_ids to reverseElementMap.
    Map->reserve(old_ids.size() + new_ids.size());
    std::merge(old_ids.begin(), old_ids.end(),
	       new_ids.begin(), new_ids.end(),
	       std::inserter(*Map, Map->begin()), IdPairCompare());
    
    // Check for duplicate ids...
    verify_no_duplicate_ids(*Map, processor);
  }
}

void Ioss::Map::verify_no_duplicate_ids(std::vector<IdPair> &reverse_map, int processor)
{
  // Check for duplicate ids...
  std::vector<IdPair>::iterator dup = std::adjacent_find(reverse_map.begin(),
							 reverse_map.end(),
							 IdPairEqual());

  if (dup != reverse_map.end()) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Duplicate global id detected on processor "
	   << processor << ".\n"
	   << "       Global id " << (*dup).first
	   << " assigned to local entities "
	   << (*dup).second << " and "
	   << (*++dup).second << ".\n";
    IOSS_ERROR(errmsg);
  }
}

bool Ioss::Map::is_sequential(const MapContainer &the_map)
{
  // Assumes the_map runs from [1..size) Slot zero will contain -1 if the
  // vector is sequential; 1 if not sequential, and 0 if it has not
  // yet been determined...
  // Once the the_map has been determined to be sequential/not-sequential,
  // slot zero is set appropriately.
  // 'sequential' is defined here to mean i==the_map[i] for all 0<i<the_map.size()

  // Check slot zero...
  if (the_map[0] == -1)
    return true;
  else if (the_map[0] ==  1)
    return false;
  else {
    size_t size = the_map.size();
    for (size_t i=1; i < size; i++)
      if (the_map[i] != (int64_t)i) {
	MapContainer &new_map = const_cast<MapContainer&>(the_map);
	new_map[0] = 1;
	return false;
      }
    MapContainer &new_map = const_cast<MapContainer&>(the_map);
    new_map[0] = -1;
    return true;
  }
}

// Node and Element mapping functions.  The ExodusII database
// stores ids in a local-id system (1..NUMNP), (1..NUMEL) but
// Sierra wants entities in a global system. These routines
// take care of the mapping from local <-> global

int64_t Ioss::Map::local_to_global(int64_t /* local */)  const
{
#if 0
  assert(local <= nodeCount && local > 0);
  const MapContainer &node_map = get_node_map();
  int64_t global = node_map[local];
  return global;
#endif
  return -1;
}

int64_t Ioss::Map::global_to_local(int64_t global, bool must_exist) const
{
  assert(entityCount >= 0);
  int64_t local = global;
  if (!sequentialG2L) {
    std::pair<RMapI, RMapI> iter = std::equal_range(reverse.begin(), reverse.end(),
						    global, Ioss::IdPairCompare());
    if (iter.first != iter.second)
      local = (iter.first)->second;
    else
      local = 0;
    assert(!must_exist || iter.first != iter.second);
  } else if (!must_exist && global > entityCount) {
    local = 0;
  }
  if (local > entityCount || (local <= 0  && must_exist)) {
    std::ostringstream errmsg;
    errmsg << "Entity (element, face, edge, node) with global id equal to " << global
	   << " returns a local id of " << local
	   << " which is invalid. This should not happen, please report.\n";
    IOSS_ERROR(errmsg);
      }
  return local;
}
