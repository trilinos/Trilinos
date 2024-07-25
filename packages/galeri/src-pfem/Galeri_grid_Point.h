// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! 
 * \file Galeri_grid_Point.h
 *
 * \brief Class for grid points.
 *
 * \author Marzio Sala, ETHZ
 *
 * \date Last modified on Aug-06
 */

#ifndef GALERI_GRID_POINT_H
#define GALERI_GRID_POINT_H

#include "Galeri_grid_Element.h"

namespace Galeri {
namespace grid {

/*!
 * \class Point
 *
 * \brief Class for grid points.
 *
 * This is the simplest geometrical object.
 */ 
class Point : public Element
{
  public:
    //! default constructor
    Point()
    {
      setLabel("Galeri::grid::Point");
      setNumVertices(1);
      setNumComponents(0);
    }
}; 

} // namespace grid
} // namespace Galeri
#endif
