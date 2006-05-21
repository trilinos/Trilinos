#ifndef PHX_GRID_TRIANGLE_H
#define PHX_GRID_TRIANGLE_H

#include "Teuchos_TestForException.hpp"

#include "phx_grid_Element.h"
#include "phx_grid_Segment.h"

namespace phx {
namespace grid {

class Triangle : public Element
{
  public:
    Triangle()
    {
      setLabel("phx::grid::Triangle");
      setNumVertices(3);
      setNumComponents(3);
      Segment component;
      for (int i = 0; i < 3; ++i)
        setComponent(i, component);
    }
};

} // namespace grid
} // namespace phx
#endif
