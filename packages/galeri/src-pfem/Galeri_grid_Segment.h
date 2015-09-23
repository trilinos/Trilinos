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
