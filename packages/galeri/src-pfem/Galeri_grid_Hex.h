// @HEADER
// ************************************************************************
//
//            Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
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

#include "Teuchos_TestForException.hpp"

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


