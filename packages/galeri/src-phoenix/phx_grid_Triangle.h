#ifndef PHX_GRID_TRIANGLE_H
#define PHX_GRID_TRIANGLE_H

#include "phx_grid_Element.h"
#include "phx_grid_Segment.h"

namespace phx {
namespace grid {

class Triangle : public Element
{
  public:
    Triangle(const int NumDimensions)
    {
      setLabel("GridTriangle");
      setNumVertices(3);
      setNumDimensions(NumDimensions);
      setNumComponents(3);
      Segment Component(NumDimensions);
      for (int i = 0; i < 3; ++i)
        setComponent(i, Component);
    }
};

} // namespace grid
} // namespace phx
#endif
