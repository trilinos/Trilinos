// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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


