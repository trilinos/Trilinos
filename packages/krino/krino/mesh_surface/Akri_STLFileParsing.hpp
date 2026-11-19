// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_MESH_SURFACE_AKRI_STLFILEPARSING_HPP_
#define KRINO_KRINO_MESH_SURFACE_AKRI_STLFILEPARSING_HPP_

#include <stk_math/StkVector.hpp>

#include <array>
#include <istream>

namespace krino
{

// Reads the facets of an ASCII or binary STL file in order: the header first (open), then one valid
// facet at a time, scaled, with facets whose normal is not finite or that are degenerate skipped
// and the vertex order flipped when sign is -1.  Errors are thrown.
class STLFacetReader
{
public:
  STLFacetReader(std::istream & input, const int sign, const stk::math::Vector3d & scale);

  void read_header();
  bool is_ascii() const { return myIsAscii; }
  // Binary files only: the number of facets the header says follow.
  unsigned num_binary_facets() const { return myNumBinaryFacets; }

  // The next valid facet, or false at the end of the file.
  bool read_next_valid_facet(std::array<stk::math::Vector3d, 3> & facetPts);

private:
  std::istream & myInput;
  int mySign;
  stk::math::Vector3d myScale;
  bool myIsAscii{true};
  unsigned myNumBinaryFacets{0};
  unsigned myNumBinaryFacetsRead{0};
};

} // namespace krino

#endif /* KRINO_KRINO_MESH_SURFACE_AKRI_STLFILEPARSING_HPP_ */
