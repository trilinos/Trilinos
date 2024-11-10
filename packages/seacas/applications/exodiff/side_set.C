// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ED_SystemInterface.h" // for SystemInterface, etc
#include "exodusII.h"           // for ex_set, etc
#include "iqsort.h"             // for index_qsort
#include "side_set.h"
#include "smart_assert.h" // for SMART_ASSERT
#include <cstdlib>        // for exit
#include <vector>         // for vector

template <typename INT> Side_Set<INT>::Side_Set() : Exo_Entity() {}

template <typename INT> Side_Set<INT>::Side_Set(int file_id, size_t id) : Exo_Entity(file_id, id)
{
  SMART_ASSERT((int)id != EX_INVALID_ID);
}

template <typename INT>
Side_Set<INT>::Side_Set(int file_id, size_t id, size_t ns, size_t ndf)
    : Exo_Entity(file_id, id, ns), num_dist_factors(ndf)
{
  SMART_ASSERT(id > 0);
}

template <typename INT> Side_Set<INT>::~Side_Set()
{
  SMART_ASSERT(Check_State());

  delete[] elmts;
  delete[] sides;
  delete[] sideIndex;
  delete[] dfIndex;
  delete[] dist_factors;
}

template <typename INT> EXOTYPE Side_Set<INT>::exodus_type() const { return EX_SIDE_SET; }

template <typename INT> void Side_Set<INT>::entity_load_params()
{
  std::vector<ex_set> sets(1);
  sets[0].id                       = id_;
  sets[0].type                     = EX_SIDE_SET;
  sets[0].entry_list               = nullptr;
  sets[0].extra_list               = nullptr;
  sets[0].distribution_factor_list = nullptr;

  int err = ex_get_sets(fileId, 1, Data(sets));

  if (err < 0) {
    Error(fmt::format("{}: Failed to get sideset parameters for sideset {}. !  Aborting...\n",
                      __func__, id_));
  }

  numEntity        = sets[0].num_entry;
  num_dist_factors = sets[0].num_distribution_factor;
}

template <typename INT> void Side_Set<INT>::apply_map(const std::vector<INT> &elmt_map)
{
  SMART_ASSERT(!elmt_map.empty());
  if (elmts != nullptr) {
    delete[] elmts;
    elmts = nullptr;
    delete[] sides;
    sides = nullptr;
    delete[] sideIndex;
    sideIndex = nullptr;
  }
  load_sides(elmt_map);
}

template <typename INT> void Side_Set<INT>::load_sides(const std::vector<INT> &elmt_map) const
{
  if ((elmts == nullptr || sides == nullptr) && numEntity > 0) {
    elmts     = new INT[numEntity];
    sides     = new INT[numEntity];
    sideIndex = new INT[numEntity];

    int err = ex_get_set(fileId, EX_SIDE_SET, id_, elmts, sides);

    if (err < 0) {
      Error(fmt::format("{}: Failed to read side set {}!  Aborting...\n", __func__, id_));
    }

    if (!elmt_map.empty()) {
      for (size_t i = 0; i < numEntity; i++) {
        elmts[i] = 1 + elmt_map[elmts[i] - 1];
      }
    }

    if (interFace.ssmap_flag) {
      for (size_t i = 0; i < numEntity; i++) {
        sideIndex[i] = i;
        elmts[i]     = elmts[i] * 8 + sides[i];
      }

      index_qsort(elmts, sideIndex, numEntity);

      // Recover elmts...
      for (size_t i = 0; i < numEntity; i++) {
        elmts[i] = elmts[i] / 8;
      }
    }
    else {
      for (size_t i = 0; i < numEntity; i++) {
        sideIndex[i] = i;
      }
    }
    SMART_ASSERT(Check_State());
  }
}

template <typename INT> void Side_Set<INT>::load_df() const
{
  if (elmts == nullptr) {
    std::vector<INT> tmp;
    load_sides(tmp);
  }

  if (dist_factors != nullptr) {
    return; // Already loaded.
  }

  dfIndex = new INT[numEntity + 1];
  SMART_ASSERT(dfIndex != nullptr);
  std::vector<int> count(numEntity);

  // Handle the sierra "universal side set" which only has a single df per face...
  if (num_dist_factors == numEntity) {
    for (size_t i = 0; i < numEntity; i++) {
      count[i] = 1;
    }
  }
  else {
    int err = ex_get_side_set_node_count(fileId, id_, Data(count));
    if (err < 0) {
      Error(fmt::format("{}: Failed to read side set node count for sideset {}!  Aborting...\n",
                        __func__, id_));
    }
  }

  // Convert raw counts to index...
  size_t index = 0;
  for (size_t i = 0; i < numEntity; i++) {
    dfIndex[i] = index;
    index += count[i];
  }
  dfIndex[numEntity] = index;

  // index value should now equal df count for this sideset...
  if (index != num_dist_factors) {
    Error(fmt::format("{}: Mismatch in distribution factor count for sideset {}, "
                      "file says there should be {},\n\t\tbut ex_get_side_set_node_count says "
                      "there should be {}!  Aborting...\n",
                      __func__, id_, num_dist_factors, index));
  }
  SMART_ASSERT(index == num_dist_factors);
  dist_factors = new double[index];
  int err      = ex_get_set_dist_fact(fileId, EX_SIDE_SET, id_, dist_factors);
  if (err < 0) {
    Error(fmt::format(
        "{}: Failed to read side set distribution factors for sideset {}!  Aborting...\n", __func__,
        id_));
  }
}

template <typename INT> const INT *Side_Set<INT>::Elements() const
{
  std::vector<INT> tmp;
  load_sides(tmp);
  return elmts;
}

template <typename INT> const INT *Side_Set<INT>::Sides() const
{
  std::vector<INT> tmp;
  load_sides(tmp);
  return sides;
}

template <typename INT> std::pair<INT, INT> Side_Set<INT>::Side_Id(size_t position) const
{
  std::vector<INT> tmp;
  load_sides(tmp);
  SMART_ASSERT(position < numEntity);
  return std::make_pair(elmts[sideIndex[position]], sides[sideIndex[position]]);
}

template <typename INT> size_t Side_Set<INT>::Side_Index(size_t position) const
{
  std::vector<INT> tmp;
  load_sides(tmp);
  SMART_ASSERT(position < numEntity);
  return sideIndex[position];
}

template <typename INT> const double *Side_Set<INT>::Distribution_Factors() const
{
  if (dist_factors == nullptr) {
    load_df();
  }
  return dist_factors;
}

template <typename INT> void Side_Set<INT>::Free_Distribution_Factors() const
{
  if (dist_factors) {
    delete[] dist_factors;
    dist_factors = nullptr;
  }
}

template <typename INT>
std::pair<INT, INT> Side_Set<INT>::Distribution_Factor_Range(size_t side) const
{
  if (dfIndex == nullptr) {
    load_df();
  }
  if (dfIndex == nullptr) {
    Error(fmt::format("{}: Failed to get distribution factors for sideset {}!  Aborting...\n",
                      __func__, id_));
  }
  size_t side_index = sideIndex[side];
  return std::make_pair(dfIndex[side_index], dfIndex[side_index + 1]);
}

template <typename INT> int Side_Set<INT>::Check_State() const
{
  SMART_ASSERT(id_ >= EX_INVALID_ID);
  SMART_ASSERT(!(id_ == EX_INVALID_ID && numEntity > 0));
  SMART_ASSERT(!(id_ == EX_INVALID_ID && num_dist_factors > 0));
  SMART_ASSERT(!(id_ == EX_INVALID_ID && elmts));
  SMART_ASSERT(!(id_ == EX_INVALID_ID && sides));
  SMART_ASSERT(!(id_ == EX_INVALID_ID && dist_factors));

  SMART_ASSERT(!(elmts && !sides));
  SMART_ASSERT(!(!elmts && sides));

  return 1;
}

template class Side_Set<int>;
template class Side_Set<int64_t>;
