#include <Akri_Triangle.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_WindingNumber.hpp>
#include <vector>

#include <Akri_WindingNumberCalcs.hpp>

namespace krino {

double compute_facet_winding_number(const stk::math::Vector3d & facetX0, const stk::math::Vector3d & facetX1, const stk::math::Vector3d & facetX2, const stk::math::Vector3d & queryLoc)
{
    return compute_exact_facet_winding_number(facetX0, facetX1, facetX2, queryLoc);
}

double compute_faceted_surface_winding_number(const std::vector<std::array<stk::math::Vector3d,3>> & surfFacets, const stk::math::Vector3d & queryLoc)
{
  double windingNumber = 0;
  for (const auto & facetCoords : surfFacets)
      windingNumber += compute_facet_winding_number(facetCoords[0], facetCoords[1], facetCoords[2], queryLoc);
  return windingNumber;
}

}
