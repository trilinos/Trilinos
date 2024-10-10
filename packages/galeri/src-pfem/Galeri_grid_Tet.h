// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! 
 * \file Galeri_grid_Tet.h
 *
 * \brief Class for grid tetrahedra.
 *
 * \author Marzio Sala, ETHZ
 *
 * \date Last modified on Aug-06
 */

#ifndef GALERI_GRID_TET_H
#define GALERI_GRID_TET_H

#include "Galeri_grid_Element.h"
#include "Galeri_grid_Segment.h"

namespace Galeri {
namespace grid {

/*!
 * \class Tet
 *
 * \brief Class for grid tetrahedra.
 *
 * A tet is composed by four vertices, and the four components are
 * Galeri::grid::Triangle's.
 */ 
class Tet : public Element
{
  public:
    Tet()
    {
      setLabel("Galeri::grid::Tet");
      setNumVertices(4);
      setNumComponents(4);
      Triangle component;
      for (int i = 0; i < 4; ++i)
        setComponent(i, component);
    }
};

} // namespace grid
} // namespace Galeri
#endif
