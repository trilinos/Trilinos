#ifndef PHX_GRID_HEX_H
#define PHX_GRID_HEX_H

#include "Teuchos_TestForException.hpp"

#include "phx_grid_Element.h"
#include "phx_grid_Segment.h"

namespace phx {
namespace grid {

class Hex : public Element
{
  public:
    Hex(const int numDimensions)
    {
      TEST_FOR_EXCEPTION(numDimensions < 3, std::out_of_range,
                         "numDimensions = " << numDimensions << ", should be > 2");

      setLabel("phx::grid::Hex");
      setNumVertices(3);
      setNumDimensions(numDimensions);
      setNumComponents(6);
      Quad component(numDimensions);
      for (int i = 0; i < 6; ++i)
        setComponent(i, component);
    }
};

} // namespace grid
} // namespace phx
#endif


