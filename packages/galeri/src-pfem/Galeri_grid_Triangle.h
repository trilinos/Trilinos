// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! 
 * \file Galeri_grid_Triangle.h
 *
 * \brief Class for grid triangles.
 *
 * \author Marzio Sala, ETHZ
 *
 * \date Last modified on Aug-06
 */

#ifndef GALERI_GRID_TRIANGLE_H
#define GALERI_GRID_TRIANGLE_H

#include "Galeri_grid_Element.h"
#include "Galeri_grid_Segment.h"

namespace Galeri {
namespace grid {

/*!
 * \class Triangle
 *
 * \brief Class for grid triangles.
 *
 * A triangle is composed by 3 vertices, and by 3 components, defined as
 * Galeri::grid::Point's.
 */ 
class Triangle : public Element
{
  public:
    Triangle()
    {
      setLabel("Galeri::grid::Triangle");
      setNumVertices(3);
      setNumComponents(3);
      Segment component;
      for (int i = 0; i < 3; ++i)
        setComponent(i, component);
    }
};

} // namespace grid
} // namespace Galeri
#endif
