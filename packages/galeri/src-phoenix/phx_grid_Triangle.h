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
    Triangle(const int numDimensions)
    {
      TEST_FOR_EXCEPTION(numDimensions > 1, std::out_of_range,
                         "numDimensions = " << numDimensions << ", should be > 1");

      setLabel("phx::grid::Triangle");
      setNumVertices(3);
      setNumDimensions(numDimensions);
      setNumComponents(3);
      Segment component(numDimensions);
      for (int i = 0; i < 3; ++i)
        setComponent(i, component);
    }
};

} // namespace grid
} // namespace phx
#endif
