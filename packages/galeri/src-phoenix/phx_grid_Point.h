#ifndef PHX_GRID_POINT_H
#define PHX_GRID_POINT_H

#include "Teuchos_TestForException.hpp"

#include "phx_grid_Element.h"

namespace phx {
namespace grid {

class Point : public Element
{
  public:
    Point(const int numDimensions)
    {
      TEST_FOR_EXCEPTION(numDimensions < 1, std::out_of_range,
                         "numDimensions = " << numDimensions << ", should be > 0");

      setLabel("phx::grid::Point");
      setNumVertices(1);
      setNumDimensions(numDimensions);
      setNumComponents(0);
    }
}; 

} // namespace grid
} // namespace phx
#endif
