#ifndef PHX_GRID_QUAD_H
#define PHX_GRID_QUAD_H

#include "phx_grid_Element.h"
#include "phx_grid_Segment.h"

namespace phx {
namespace grid {

class Quad : public Element
{
  public:
    Quad(const int NumDimensions)
    {
      setLabel("GridQuad");
      setNumVertices(4);
      setNumDimensions(NumDimensions);
      setNumComponents(4);
      Segment Component(NumDimensions);
      for (int i = 0; i < 3; ++i)
        setComponent(i, Component);
    }
}; 

} // namespace grid
} // namespace phx
#endif
