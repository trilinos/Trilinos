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
  static map<string, RefCountPtr<Loadable> >
  getSquareWithQuad(Epetra_Comm& comm, const int nx, const int ny,
                    const int mx, const int my)
  {
    double lx = 1.0;
    double ly = 1.0;

    int numGlobalElements = nx * ny;
    int numGlobalVertices = (nx + 1) * (ny + 1);

    double deltax = lx / nx;
    double deltay = ly / ny;

    RefCountPtr<Epetra_Map> elementMap = rcp(new Epetra_Map(numGlobalElements, 0, comm));
    RefCountPtr<phx::grid::Element> quad = rcp(new phx::grid::Quad);
    RefCountPtr<phx::grid::Loadable> domain = rcp(new Loadable(elementMap, quad));

    for (int i = 0; i < numGlobalElements; ++i)
    {
      int ix = i % nx;
      int iy = i / nx;

      int base = ix + iy * (nx + 1);
      domain->setGlobalConnectivity(i, 0, base);
      domain->setGlobalConnectivity(i, 1, base + 1);
      domain->setGlobalConnectivity(i, 2, base + 2 + nx);
      domain->setGlobalConnectivity(i, 3, base + 1 + nx);
    }

    domain->freezeConnectivity();

    for (int i = 0; i < numGlobalVertices; ++i)
    {
      int ix = i % (nx + 1);
      int iy = i / (nx + 1);

      domain->setGlobalCoordinates(i, 0, ix * deltax);
      domain->setGlobalCoordinates(i, 1, iy * deltay);
    }

    domain->freezeCoordinates();

    cout << *domain;
    // now build boundary faces
    int numGlobalBoundaries = nx;
    RefCountPtr<Epetra_Map> boundaryMap = rcp(new Epetra_Map(numGlobalBoundaries, 0, comm));
    RefCountPtr<phx::grid::Element> segment = rcp(new phx::grid::Segment);
    RefCountPtr<phx::grid::Loadable> boundary = rcp(new Loadable(boundaryMap, segment));

    for (int i = 0; i < nx; ++i)
    {
      boundary->setGlobalConnectivity(i, 0, i);
      boundary->setGlobalConnectivity(i, 1, i + 1);
    }

    boundary->freezeConnectivity();

    for (int i = 0; i < nx + 1; ++i)
    {
      boundary->setGlobalCoordinates(i, 0, i * deltax);
      boundary->setGlobalCoordinates(i, 1, 0.0);
    }

    boundary->freezeCoordinates();

    map<string, RefCountPtr<Loadable> > patches;
    patches["domain"] = domain;
    patches["boundary"] = boundary;

    return(patches);
  }


}; // class Generator

} // namespace grid
} // namespace phx

#endif
