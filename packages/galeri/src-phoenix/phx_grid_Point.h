#ifndef PHX_GRID_POINT_H
#define PHX_GRID_POINT_H

#include "Teuchos_TestForException.hpp"

#include "phx_grid_Element.h"

namespace phx {
namespace grid {

class Point : public Element
{
  public:
    Point()
    {
      setLabel("phx::grid::Point");
      setNumVertices(1);
      setNumComponents(0);
    }
}; 

} // namespace grid
} // namespace phx
#endif
