// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! 
 * \file Galeri_grid_Segment.h
 *
 * \brief Class for grid segment.
 *
 * \author Marzio Sala, ETHZ
 *
 * \date Last modified on Aug-06
 */

#ifndef GALERI_GRID_SEGMENT_H
#define GALERI_GRID_SEGMENT_H

#include "Galeri_grid_Element.h"
#include "Galeri_grid_Point.h"

namespace Galeri {
namespace grid {

/*!
 * \class Segment
 *
 * \brief Class for grid segments.
 *
 * A grid segment is composed by 2 vertices, and 2 components, defined as
 * Galeri::grid::Point's.
 */ 
class Segment : public Element
{
  public:
    Segment()
    {
      setLabel("Galeri::grid::Segment");
      setNumVertices(2);
      setNumComponents(2);
      Point component;
      for (int i = 0; i < 2; ++i)
        setComponent(i, component);
    }
};

} // namespace grid
} // namespace Galeri
#endif
