// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//
//     * Neither the name of NTESS nor the names of its
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

#include <Ioss_Field.h> // for Field, etc
#include <Ioss_Map.h>
#include <Ioss_Sort.h>
#include <Ioss_Utils.h> // for IOSS_ERROR
#include <algorithm>    // for adjacent_find, lower_bound, etc
#include <cassert>      // for assert
#include <cstddef>      // for size_t
#include <iterator>     // for insert_iterator, inserter
#include <sstream>      // for operator<<, basic_ostream, etc
#include <string>       // for char_traits, operator<<, etc
#include <sys/types.h>  // for ssize_t
#include <utility>      // for pair, make_pair
#include <vector>       // for vector, vector<>::iterator, etc

namespace {
  // Determines whether the input map is sequential (m_map[i] == i)
  bool is_sequential(const Ioss::MapContainer &the_map)
  {
    // Assumes the_map runs from [1..size) Slot zero will contain -1 if the
    // vector is sequential; 1 if not sequential, and 0 if it has not
    // yet been determined...
    // Once the the_map has been determined to be sequential/not-sequential,
    // slot zero is set appropriately.
    // 'sequential' is defined here to mean i==the_map[i] for all
    // 0<i<the_map.size()

    // Check slot zero...
    if (the_map[0] == -1) {
      return true;
    }
    if (the_map[0] == 1) {
      return false;
    }

    Ioss::MapContainer &new_map = const_cast<Ioss::MapContainer &>(the_map);
    size_t              size    = the_map.size();
    for (size_t i = 1; i < size; i++) {
      if (the_map[i] != static_cast<int64_t>(i)) {
        new_map[0] = 1;
        return false;
      }
    }
    new_map[0] = -1;
    return true;
  }

  // map global to local ids

  // Class to support storing global/local element id map in sorted vector...
  class IdPairCompare
  {
  public:
    IdPairCompare() = default;
    bool operator()(const Ioss::IdPair &lhs, const Ioss::IdPair &rhs) const
    {
      return key_less(lhs.first, rhs.first);
    }
    bool operator()(const Ioss::IdPair &lhs, const Ioss::IdPair::first_type &k) const
    {
      return key_less(lhs.first, k);
    }
    // Assignment operator
    // Copy constructor
  private:
    bool key_less(const Ioss::IdPair::first_type &k1, const Ioss::IdPair::first_type &k2) const
    {
      return k1 < k2;
    }
  };

  class IdPairEqual
  {
  public:
    IdPairEqual() = default;
    bool operator()(const Ioss::IdPair &lhs, const Ioss::IdPair &rhs) const
    {
      return key_equal(lhs.first, rhs.first);
    }
    // Assignment operator
    // Copy constructor
  private:
    bool key_equal(const Ioss::IdPair::first_type &k1, const Ioss::IdPair::first_type &k2) const
    {
      return k1 == k2;
    }
  };

  using RMapI = std::vector<Ioss::IdPair>::const_iterator;

  template <typename INT>
  void map_implicit_data_internal(INT *ids, size_t count, const Ioss::MapContainer &map,
                                  size_t offset)
  {
    // Map the "local" ids (offset+1..offset+count) to the global ids. The local
    // ids are implicit
    if (is_sequential(map)) {
      for (size_t i = 0; i < count; i++) {
        ids[i] = offset + 1 + i;
      }
    }
    else {
      for (size_t i = 0; i < count; i++) {
        ids[i] = map[offset + 1 + i];
      }
    }
  }

} // namespace

void Ioss::Map::release_memory()
{
  IOSS_FUNC_ENTER(m_);
  MapContainer().swap(m_map);
  MapContainer().swap(m_reorder);
  ReverseMapContainer().swap(m_reverse);
}

void Ioss::Map::build_reverse_map()
{
  if (m_map[0] == 1) {
    build_reverse_map(m_map.size() - 1, 0);
  }
}

void Ioss::Map::build_reverse_map(int64_t num_to_get, int64_t offset)
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
  IOSS_FUNC_ENTER(m_);
  ReverseMapContainer new_ids(num_to_get);
  for (int64_t i = 0; i < num_to_get; i++) {
    int64_t local_id = offset + i + 1;
    new_ids[i]       = std::make_pair(m_map[local_id], local_id);

    if (m_map[local_id] <= 0) {
      std::ostringstream errmsg;
      errmsg << "\nERROR: " << m_entityType << " map detected non-positive global id "
             << m_map[local_id] << " for " << m_entityType << " with local id " << local_id
             << " on processor " << m_myProcessor << ".\n";
      IOSS_ERROR(errmsg);
    }
  }

  // Sort that vector...
  Ioss::qsort(new_ids);

  int64_t new_id_min = new_ids.empty() ? 0 : new_ids.front().first;
  int64_t old_id_max = m_reverse.empty() ? 0 : m_reverse.back().first;
  if (new_id_min > old_id_max) {
    m_reverse.insert(m_reverse.end(), new_ids.begin(), new_ids.end());
  }
  else {
    // Copy reverseElementMap to old_ids, empty reverseElementMap.
    ReverseMapContainer old_ids;
    old_ids.swap(m_reverse);
    assert(m_reverse.empty());

    // Merge old_ids and new_ids to reverseElementMap.
    m_reverse.reserve(old_ids.size() + new_ids.size());
    std::merge(old_ids.begin(), old_ids.end(), new_ids.begin(), new_ids.end(),
               std::inserter(m_reverse, m_reverse.begin()), IdPairCompare());
  }
// Check for duplicate ids...
#ifndef NDEBUG
  verify_no_duplicate_ids(m_reverse);
#endif
}

void Ioss::Map::verify_no_duplicate_ids(std::vector<Ioss::IdPair> &reverse_map)
{
  // Check for duplicate ids...
  auto dup = std::adjacent_find(reverse_map.begin(), reverse_map.end(), IdPairEqual());

  if (dup != reverse_map.end()) {
    auto               other = dup + 1;
    std::ostringstream errmsg;
    errmsg << "\nERROR: Duplicate " << m_entityType << " global id detected on processor "
           << m_myProcessor << ", filename '" << m_filename << "'.\n"
           << "       Global id " << (*dup).first << " assigned to local " << m_entityType << "s "
           << (*dup).second << " and " << (*other).second << ".\n";
    IOSS_ERROR(errmsg);
  }
}

template void Ioss::Map::set_map(int *ids, size_t count, size_t offset);
template void Ioss::Map::set_map(int64_t *ids, size_t count, size_t offset);

template <typename INT> void Ioss::Map::set_map(INT *ids, size_t count, size_t offset)
{
  IOSS_FUNC_ENTER(m_);
  for (size_t i = 0; i < count; i++) {
    ssize_t local_id = offset + i + 1;
    m_map[local_id]  = ids[i];
    if (local_id != ids[i]) {
      m_map[0] = 1;
    }
    if (ids[i] <= 0) {
      std::ostringstream errmsg;
      errmsg << "\nERROR: " << m_entityType << " mapping routines detected non-positive global id "
             << ids[i] << " for local id " << local_id << " on processor " << m_myProcessor
             << ", filename '" << m_filename << "'.\n";
      IOSS_ERROR(errmsg);
    }
  }
}

void Ioss::Map::reverse_map_data(void *data, const Ioss::Field &field, size_t count) const
{
  IOSS_FUNC_ENTER(m_);
  assert(!m_map.empty());
  if (!is_sequential(m_map)) {
    if (field.get_type() == Ioss::Field::INTEGER) {
      int *connect = static_cast<int *>(data);
      for (size_t i = 0; i < count; i++) {
        int global_id = connect[i];
        connect[i]    = global_to_local__(global_id, true);
      }
    }
    else {
      int64_t *connect = static_cast<int64_t *>(data);
      for (size_t i = 0; i < count; i++) {
        int64_t global_id = connect[i];
        connect[i]        = global_to_local__(global_id, true);
      }
    }
  }
}

void Ioss::Map::map_data(void *data, const Ioss::Field &field, size_t count) const
{
  IOSS_FUNC_ENTER(m_);
  if (!is_sequential(m_map)) {
    if (field.get_type() == Ioss::Field::INTEGER) {
      int *datum = static_cast<int *>(data);
      for (size_t i = 0; i < count; i++) {
        datum[i] = m_map[datum[i]];
      }
    }
    else {
      int64_t *datum = static_cast<int64_t *>(data);
      for (size_t i = 0; i < count; i++) {
        datum[i] = m_map[datum[i]];
      }
    }
  }
}

void Ioss::Map::map_implicit_data(void *data, const Ioss::Field &field, size_t count,
                                  size_t offset) const
{
  IOSS_FUNC_ENTER(m_);
  if (field.get_type() == Ioss::Field::INTEGER) {
    map_implicit_data_internal(static_cast<int *>(data), count, m_map, offset);
  }
  else {
    map_implicit_data_internal(static_cast<int64_t *>(data), count, m_map, offset);
  }
}

template size_t Ioss::Map::map_field_to_db_scalar_order(double *             variables,
                                                        std::vector<double> &db_var,
                                                        size_t begin_offset, size_t count,
                                                        size_t stride, size_t offset);
template size_t Ioss::Map::map_field_to_db_scalar_order(int *variables, std::vector<double> &db_var,
                                                        size_t begin_offset, size_t count,
                                                        size_t stride, size_t offset);
template size_t Ioss::Map::map_field_to_db_scalar_order(int64_t *            variables,
                                                        std::vector<double> &db_var,
                                                        size_t begin_offset, size_t count,
                                                        size_t stride, size_t offset);

template <typename T>
size_t Ioss::Map::map_field_to_db_scalar_order(T *variables, std::vector<double> &db_var,
                                               size_t begin_offset, size_t count, size_t stride,
                                               size_t offset)
{
  IOSS_FUNC_ENTER(m_);
  size_t num_out = 0;
  if (!m_reorder.empty()) {
    size_t k = offset;
    for (size_t j = begin_offset; j < count * stride; j += stride) {
      // Map to storage location.
      ssize_t where = m_reorder[k++] - offset;
      if (where >= 0) {
        assert(where < (ssize_t)count);
        db_var[where] = variables[j];
        num_out++;
      }
    }
  }
  else {
    size_t k = 0;
    for (size_t j = begin_offset; j < count * stride; j += stride) {
      // Map to storage location.
      db_var[k++] = variables[j];
    }
    num_out = count;
  }
  return num_out;
}

void Ioss::Map::build_reorder_map(int64_t start, int64_t count)
{
  // This routine builds a map that relates the current node id order
  // to the original node ordering in affect at the time the file was
  // created. That is, the node map used to define the topology of the
  // model.  Now, if there are changes in node ordering at the
  // application level, we build the node reorder map to map the
  // current order into the original order.  An added complication is
  // that this is more than just a reordering... It may be that the
  // application has 'ghosted' nodes that it doesnt want put out on
  // the database, so the reorder map must handle a node that is not
  // in the original mesh and map that to an invalid value (currently
  // using -1 as invalid value...)

  // Note: To further add confusion, the reorder map is 0-based
  // and the reverse map and 'map' are 1-baed. This is
  // just a consequence of how they are intended to be used...
  //
  // start is based on a 0-based array -- start of the reorderMap to build.

  IOSS_FUNC_ENTER(m_);
  if (m_reorder.empty()) {
    // See if actually need a reorder map first...
    bool    need_reorder_map = false;
    int64_t my_end           = start + count;
    for (int64_t i = start; i < my_end; i++) {
      int64_t global_id     = m_map[i + 1];
      int64_t orig_local_id = global_to_local__(global_id) - 1;

      // The reordering should only be a permutation of the original
      // ordering within this entity block...
      assert(orig_local_id >= start && orig_local_id <= my_end);
      if (i != orig_local_id) {
        need_reorder_map = true;
        break;
      }
    }
    if (need_reorder_map) {
      int64_t map_size = m_map.size() - 1;
      m_reorder.resize(map_size);
      // If building a partial reorder map, assume all entries
      // are a direct 1-1 and then let the partial fills overwrite
      // if needed.
      if (start > 0 || my_end < map_size) {
        for (size_t i = 0; i < m_reorder.size(); i++) {
          m_reorder[i] = i;
        }
      }
    }
    else {
      return;
    }
  }

  int64_t my_end = start + count;
  for (int64_t i = start; i < my_end; i++) {
    int64_t global_id     = m_map[i + 1];
    int64_t orig_local_id = global_to_local__(global_id) - 1;

    // The reordering should only be a permutation of the original
    // ordering within this entity block...
    assert(orig_local_id >= start && orig_local_id <= my_end);
    m_reorder[i] = orig_local_id;
  }
}

// Node and Element mapping functions.  The ExodusII database
// stores ids in a local-id system (1..NUMNP), (1..NUMEL) but
// Sierra wants entities in a global system. These routines
// take care of the mapping from local <-> global

int64_t Ioss::Map::global_to_local(int64_t global, bool must_exist) const
{
  IOSS_FUNC_ENTER(m_);
  return global_to_local__(global, must_exist);
}

int64_t Ioss::Map::global_to_local__(int64_t global, bool must_exist) const
{
  int64_t local = global;
  if (m_map[0] == 1) {
    auto iter = std::lower_bound(m_reverse.begin(), m_reverse.end(), global, IdPairCompare());
    if (iter != m_reverse.end() && iter->first == global) {
      local = iter->second;
    }
    else {
      local = 0;
    }
  }
  else if (!must_exist && global > static_cast<int64_t>(m_map.size()) - 1) {
    local = 0;
  }
  if (local > static_cast<int64_t>(m_map.size()) - 1 || (local <= 0 && must_exist)) {
    std::ostringstream errmsg;
    errmsg << "ERROR: Ioss Mapping routines detected " << m_entityType
           << " with global id equal to " << global << " returns a local id of " << local
           << " which is invalid\n"
           << "on processor " << m_myProcessor << ", filename '" << m_filename << "'.\n"
           << "This should not happen, please report.\n";
    IOSS_ERROR(errmsg);
  }
  return local;
}
