/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Map.h>

#include <assert.h>
#include <Ioss_Utils.h>
#include <string>
#include <algorithm>
#include <sstream>

Ioss::Map::Map() :
  entityCount(-1), sequentialG2L(true), entityReordered(false)
{}

void Ioss::Map::build_reverse_map(ReverseMapContainer *Map, const int *ids,
				 int num_to_get, int offset,
				 const std::string& type, int processor)
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
    verify_no_duplicate_ids(*Map, type, processor);
  }

}

void Ioss::Map::verify_no_duplicate_ids(std::vector<IdPair> &reverse_map,
				  const std::string &type, int processor)
{
  // Check for duplicate ids...
  std::vector<IdPair>::iterator dup = std::adjacent_find(reverse_map.begin(),
							 reverse_map.end(),
							 IdPairEqual());

  if (dup != reverse_map.end()) {
    std::ostringstream errmsg;
    errmsg << "\nERROR: Duplicate global " << type << " id detected on processor "
	   << processor << ".\n"
	   << "       Global id " << (*dup).first
	   << " assigned to local " << type << "s "
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
      if (the_map[i] != (int)i) {
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

int Ioss::Map::local_to_global(int /* local */)  const
{
#if 0
  assert(local <= nodeCount && local > 0);
  const std::vector<int> node_map = get_node_map();
  int global = node_map[local];
  return global;
#endif
  return -1;
}

typedef std::vector<Ioss::IdPair>::const_iterator RMapI;
int Ioss::Map::global_to_local(int global, bool must_exist) const
{
  assert(entityCount >= 0);
  int local = global;
  if (!sequentialG2L) {
    std::pair<RMapI, RMapI> iter = std::equal_range(reverseMap.begin(), reverseMap.end(),
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
    errmsg << "Entity with global id equal to " << global
	   << " returns a local id of " << local
	   << " which is invalid. This should not happen, please report.\n";
    IOSS_ERROR(errmsg);
      }
  return local;
}
