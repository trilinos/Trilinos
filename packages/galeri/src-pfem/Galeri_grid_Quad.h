// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! 
 * \file Galeri_grid_Quad.h
 *
 * \brief Class for grid quadrilaterals.
 *
 * \author Marzio Sala, ETHZ
 *
 * \date Last modified on Aug-06
 */

#ifndef GALERI_GRID_QUAD_H
#define GALERI_GRID_QUAD_H

#include "Galeri_grid_Element.h"
#include "Galeri_grid_Segment.h"

namespace Galeri {
namespace grid {

/*!
 * \class Quad
 *
 * \brief Class for grid quadrilaterals.
 *
 * A quad is composed by 4 vertices, and by 4 components, defined as
 * Galeri::grid::Segment's.
 */ 
class Quad : public Element
{
  public:
    Quad()
    {
      setLabel("Galeri::grid::Quad");
      setNumVertices(4);
      setNumComponents(4);
      Segment component;
      for (int i = 0; i < 4; ++i)
        setComponent(i, component);
    }
}; 

} // namespace grid
} // namespace Galeri
#endif
