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

#ifndef GALERI_GRID_GENERATOR_H
#define GALERI_GRID_GENERATOR_H

#include "Galeri_grid_Loadable.h"
#include "Galeri_grid_Quad.h"

#include "Teuchos_RefCountPtr.hpp"

namespace Galeri {
namespace grid {

class Generator 
{
public:
  static void
  getSquareWithTriangles(Epetra_Comm& comm, 
                         const int numGlobalElementsX, const int numGlobalElementsY,
                         const int numDomainsX, const int numDomainsY,
                         Galeri::grid::Loadable& domain, Galeri::grid::Loadable& boundary);

  static void
  getSquareWithQuads(Epetra_Comm& comm, 
                    const int numGlobalElementsX, const int numGlobalElementsY,
                    const int numDomainsX, const int numDomainsY,
                    Galeri::grid::Loadable& domain, Galeri::grid::Loadable& boundary);

  static void
  getSquare(Epetra_Comm& comm, 
            const int numGlobalElementsX, const int numGlobalElementsY,
            const int numDomainsX, const int numDomainsY,
            Galeri::grid::Loadable& domain, Galeri::grid::Loadable& boundary,
            const string what);

  static void
  getCubeWithHexs(Epetra_Comm& comm, 
                  const int numGlobalElementsX, const int numGlobalElementsY, const int numGlobalElementsZ,
                  const int numDomainsX, const int numDomainsY, const int numDomainsZ,
                  Galeri::grid::Loadable& domain, Galeri::grid::Loadable& boundary);

}; // class Generator

} // namespace grid
} // namespace Galeri

#endif
