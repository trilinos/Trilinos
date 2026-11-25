// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SURFACE_AKRI_FACETCOORDINATES_HPP_
#define KRINO_KRINO_SURFACE_AKRI_FACETCOORDINATES_HPP_

#include <Akri_Facet.hpp>
#include <Akri_FacetAppender.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <cstddef>
#include <vector>

namespace krino
{

// Facets (segments in 2D, triangles in 3D) each given by the coordinates of its own vertices, with
// no shared vertices: a facet soup, as contouring writes.
class FacetCoordinates : public FacetAppender
{
public:
  explicit FacetCoordinates(const unsigned dim) : myDim(dim)
  {
    STK_ThrowRequireMsg(
        dim == 2 || dim == 3, "FacetCoordinates supports dimension 2 or 3, not " << dim << ".");
  }

  unsigned spatial_dimension() const { return myDim; }
  // The vertices of facet f are coordinates()[dim*f] to coordinates()[dim*f+dim-1], in order.
  const std::vector<stk::math::Vector3d> & coordinates() const { return myCoords; }
  // The facets as FACETs (Facet2d or Facet3d) with copies of their vertices' coordinates.
  template <class FACET> std::vector<FACET> make_facets() const
  {
    STK_ThrowRequireMsg(FACET::DIM == int(myDim),
        "FacetCoordinates of dimension " << myDim << " made into facets of dimension " << FACET::DIM
                                         << ".");
    std::vector<FACET> facets;
    facets.reserve(num_facets());
    for (size_t f = 0; f < num_facets(); ++f)
    {
      if constexpr (FACET::DIM == 3)
        facets.emplace_back(myCoords[3 * f], myCoords[3 * f + 1], myCoords[3 * f + 2]);
      else
        facets.emplace_back(myCoords[2 * f], myCoords[2 * f + 1]);
    }
    return facets;
  }

  virtual void
  emplace_back_2d(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1) override
  {
    STK_ThrowAssert(2 == myDim);
    myCoords.push_back(x0);
    myCoords.push_back(x1);
  }
  virtual void emplace_back_3d(const stk::math::Vector3d & x0,
      const stk::math::Vector3d & x1,
      const stk::math::Vector3d & x2) override
  {
    STK_ThrowAssert(3 == myDim);
    myCoords.push_back(x0);
    myCoords.push_back(x1);
    myCoords.push_back(x2);
  }
  virtual void clear() override { myCoords.clear(); }
  virtual size_t num_facets() const override { return myCoords.size() / myDim; }

private:
  unsigned myDim;
  std::vector<stk::math::Vector3d> myCoords;
};

} // namespace krino

#endif /* KRINO_KRINO_SURFACE_AKRI_FACETCOORDINATES_HPP_ */
