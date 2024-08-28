// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! 
 * \file Galeri_grid_Hex.h
 *
 * \brief Class for grid hexahedra.
 *
 * \author Marzio Sala, ETHZ
 *
 * \date Last modified on Aug-06
 */

#ifndef GALERI_GRID_HEX_H
#define GALERI_GRID_HEX_H

#include "Teuchos_Assert.hpp"

#include "Galeri_grid_Element.h"
#include "Galeri_grid_Segment.h"

namespace Galeri {
namespace grid {

/*!
 * \class Hex
 *
 * \brief Class for grid hexahedra.
 *
 * A hexahedron is composed by 8 vertices, and by 8 components, defined as
 * Galeri::grid::Quad's.
 */ 
class Hex : public Element
{
  public:
    Hex()
    {
      setLabel("Galeri::grid::Hex");
      setNumVertices(8);
      setNumComponents(6);
      Quad component;
      for (int i = 0; i < 6; ++i)
        setComponent(i, component);
    }
};

} // namespace grid
} // namespace Galeri
#endif
