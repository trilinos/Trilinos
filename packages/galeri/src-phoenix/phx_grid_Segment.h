#ifndef PHX_GRID_SEGMENT_H
#define PHX_GRID_SEGMENT_H

#include "phx_grid_Element.h"
#include "phx_grid_Point.h"

namespace phx {
namespace grid {

class Segment : public Element
{
  public:
    Segment(const int NumDimensions)
    {
      setLabel("GridSegment");
      setNumVertices(2);
      setNumDimensions(NumDimensions);
      setNumComponents(2);
      Point Component(NumDimensions);
      for (int i = 0; i < 2; ++i)
        setComponent(i, Component);
    }
};

} // namespace grid
} // namespace phx
#endif
