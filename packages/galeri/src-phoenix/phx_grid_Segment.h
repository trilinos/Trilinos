#ifndef PHX_GRID_SEGMENT_H
#define PHX_GRID_SEGMENT_H

#include "phx_grid_Element.h"
#include "phx_grid_Point.h"

namespace phx {
namespace grid {

class Segment : public Element
{
  public:
    Segment()
    {
      setLabel("phx::grid::Segment");
      setNumVertices(2);
      setNumComponents(2);
      Point component;
      for (int i = 0; i < 2; ++i)
        setComponent(i, component);
    }
};

} // namespace grid
} // namespace phx
#endif
