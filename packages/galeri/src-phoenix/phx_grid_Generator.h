// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef PHX_GRID_GENERATOR_H
#define PHX_GRID_GENERATOR_H

#include "phx_grid_Loadable.h"
#include "phx_grid_Quad.h"

#include "Teuchos_RefCountPtr.hpp"

namespace phx {
namespace grid {

class Generator 
{
public:
  static void
  getSquareWithTriangles(Epetra_Comm& comm, 
                         const int numGlobalElementsX, const int numGlobalElementsY,
                         const int numDomainsX, const int numDomainsY,
                         phx::grid::Loadable& domain, phx::grid::Loadable& boundary);

  static void
  getSquareWithQuads(Epetra_Comm& comm, 
                    const int numGlobalElementsX, const int numGlobalElementsY,
                    const int numDomainsX, const int numDomainsY,
                    phx::grid::Loadable& domain, phx::grid::Loadable& boundary);

  static void
  getSquare(Epetra_Comm& comm, 
            const int numGlobalElementsX, const int numGlobalElementsY,
            const int numDomainsX, const int numDomainsY,
            phx::grid::Loadable& domain, phx::grid::Loadable& boundary,
            const string what);

  static void
  getCubeWithHexs(Epetra_Comm& comm, 
                  const int numGlobalElementsX, const int numGlobalElementsY, const int numGlobalElementsZ,
                  const int numDomainsX, const int numDomainsY, const int numDomainsZ,
                  phx::grid::Loadable& domain, phx::grid::Loadable& boundary);

}; // class Generator

} // namespace grid
} // namespace phx

#endif
