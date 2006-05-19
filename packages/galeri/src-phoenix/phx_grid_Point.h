#ifndef PHX_GRID_POINT_H
#define PHX_GRID_POINT_H

#include "phx_grid_Element.h"

namespace phx {
namespace grid {

class Point : public Element
{
  public:
    Point(const int NumDimensions)
    {
      setLabel("GridPoint");
      setNumVertices(1);
      setNumDimensions(NumDimensions);
      setNumComponents(0);
    }
}; 

} // namespace grid
} // namespace phx
#endif
