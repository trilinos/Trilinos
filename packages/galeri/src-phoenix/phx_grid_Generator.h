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
