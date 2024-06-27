// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
            const std::string what);

  static void
  getCubeWithHexs(Epetra_Comm& comm, 
                  const int numGlobalElementsX, const int numGlobalElementsY, const int numGlobalElementsZ,
                  const int numDomainsX, const int numDomainsY, const int numDomainsZ,
                  Galeri::grid::Loadable& domain, Galeri::grid::Loadable& boundary);

}; // class Generator

} // namespace grid
} // namespace Galeri

#endif
