#ifndef PHX_GRID_TET_H
#define PHX_GRID_TET_H

#include "Teuchos_TestForException.hpp"

#include "phx_grid_Element.h"
#include "phx_grid_Segment.h"

namespace phx {
namespace grid {

class Tet : public Element
{
  public:
    Tet()
    {
      setLabel("phx::grid::Tet");
      setNumVertices(4);
      setNumComponents(4);
      Triangle component;
      for (int i = 0; i < 4; ++i)
        setComponent(i, component);
    }
};

} // namespace grid
} // namespace phx
#endif

