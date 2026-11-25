// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_UNIT_TEST_UTILS_AKRI_UNITCUBITJOURNALUTILS_HPP_
#define KRINO_KRINO_UNIT_TEST_UTILS_AKRI_UNITCUBITJOURNALUTILS_HPP_

#include <stk_math/StkVector.hpp>
#include <array>
#include <string>
#include <vector>

namespace krino {

// What a CubitJournal wrote, read back for tests: the nodes, the elements (as indices into the
// nodes) and the vertices, in the order written.
struct ParsedCubitJournal
{
  std::vector<stk::math::Vector3d> nodes;
  std::vector<std::array<int,2>> edges;
  std::vector<std::array<int,3>> tris;
  std::vector<std::array<int,4>> tets;
  std::vector<stk::math::Vector3d> vertices;
};

ParsedCubitJournal parse_cubit_journal(const std::string & journal);
double total_area_of_tris(const ParsedCubitJournal & journal);
double total_length_of_edges(const ParsedCubitJournal & journal);

}

#endif /* KRINO_KRINO_UNIT_TEST_UTILS_AKRI_UNITCUBITJOURNALUTILS_HPP_ */
