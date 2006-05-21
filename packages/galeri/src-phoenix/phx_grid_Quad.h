#ifndef PHX_GRID_QUAD_H
#define PHX_GRID_QUAD_H

#include "Teuchos_TestForException.hpp"

#include "phx_grid_Element.h"
#include "phx_grid_Segment.h"

namespace phx {
namespace grid {

class Quad : public Element
{
  public:
    Quad()
    {
      setLabel("phx::grid::Quad");
      setNumVertices(4);
      setNumComponents(4);
      Segment component;
      for (int i = 0; i < 4; ++i)
        setComponent(i, component);
    }
}; 

} // namespace grid
} // namespace phx
#endif
