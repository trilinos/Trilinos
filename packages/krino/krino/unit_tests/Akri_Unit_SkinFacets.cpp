#include <gtest/gtest.h>
#include <stk_math/StkVector.hpp>
#include <Akri_WindingNumber.hpp>
#include <Akri_Facet.hpp>
#include <Akri_SearchTree.hpp>
#include <Akri_MeshSurface.hpp>
#include <Akri_OutputUtils.hpp>
#include <stk_util/environment/EnvData.hpp>

namespace krino {

double compute_facet_area(const std::vector<Facet3d> & facets)
{
  double area = 0.;
  for (const auto & facet : facets)
    area += facet.facet_area();
  return area;
}

double compute_facet_size(const std::vector<Facet3d> & facets)
{
  const double area = compute_facet_area(facets);
  const double equilateralTriSide = std::sqrt(area/facets.size()) * 2.*std::pow(3.,-0.25);
  return equilateralTriSide;
}

double compute_winding_number(const std::vector<const Facet3d *> & windingNumberFacets, const stk::math::Vector3d &x)
{
  double windingNumber = 0;
  for (const auto & facet : windingNumberFacets)
      windingNumber += compute_facet_winding_number(facet->facet_vertex(0), facet->facet_vertex(1), facet->facet_vertex(2), x);
  return windingNumber;
}

int compute_sign_from_winding_number(const std::vector<const Facet3d *> & windingNumberFacets, const stk::math::Vector3d &x)
{
  const double windingNumber = compute_winding_number(windingNumberFacets, x);
  return (windingNumber > 0.5) ? -1 : 1;
}

double safe_point_signed_distance(const std::vector<const Facet3d *> & windingNumberFacets, SearchTree<const Facet3d*> & facetTree, const stk::math::Vector3d &x, const bool alwaysUseWinding)
{
  std::vector<const Facet3d*> nearestFacets;
  facetTree.find_closest_entities( x, nearestFacets, 0. );

  double distanceSqrMag = std::numeric_limits<double>::max();
  bool haveNeg = false;
  bool havePos = false;
  for (const auto * facet : nearestFacets)
  {
      if (facet->point_distance_sign(x) < 0)
        haveNeg = true;
      else
        havePos = true;
      distanceSqrMag = std::min(distanceSqrMag, facet->point_distance_squared(x));
  }
  const double distanceMag = std::sqrt(distanceSqrMag);
  if (alwaysUseWinding || (haveNeg && havePos))
  {
    return distanceMag * compute_sign_from_winding_number(windingNumberFacets, x);
  }
  if (haveNeg)
    return -distanceMag;
  return distanceMag;
}

bool is_facet_on_exterior(const std::vector<const Facet3d *> & windingNumberFacets, SearchTree<const Facet3d*> & facetTree, const Facet3d & facet, const double defeatureSize, const double tol)
{
  const stk::math::Vector3d facetCentroid = facet.centroid();
  const stk::math::Vector3d facetNormal = facet.facet_normal();

  const double centroidDisplacedDist = safe_point_signed_distance(windingNumberFacets, facetTree, facetCentroid+defeatureSize*facetNormal, false);

  if (centroidDisplacedDist > tol*defeatureSize)
    return true;

  return false;
}

std::vector<const Facet3d *> prune_interior_facets(const std::vector<const Facet3d *> & windingNumberFacets, const std::vector<Facet3d> & facetsToConsider, const double defeatureSize, const double tol)
{
  using FACET = Facet3d;
  SearchTree<const FACET*> facetTree( windingNumberFacets, FACET::get_centroid, FACET::insert_into_bounding_box );

  std::vector<const FACET *> retainedFacets;
  for (const auto & facet : facetsToConsider)
    if (is_facet_on_exterior(windingNumberFacets, facetTree, facet, defeatureSize, tol))
      retainedFacets.push_back(&facet);
  return retainedFacets;
}

std::vector<const Facet3d *> keep_boundary_facets(const std::vector<const Facet3d *> & windingNumberFacets, const std::vector<Facet3d> & facetsToConsider, const double defeatureSize)
{
  std::vector<const Facet3d *> retainedFacets;
  for (const auto & facet : facetsToConsider)
  {
      const stk::math::Vector3d facetCentroid = facet.centroid();
      const stk::math::Vector3d facetNormal = facet.facet_normal();

      const double windingNumberPos = compute_winding_number(windingNumberFacets, facetCentroid+0.25*defeatureSize*facetNormal);
      const double windingNumberNeg = compute_winding_number(windingNumberFacets, facetCentroid-0.25*defeatureSize*facetNormal);

      if (windingNumberPos<0.5 && windingNumberNeg>0.5)
          retainedFacets.push_back(&facet);
  }
  return retainedFacets;
}

template<typename FACET>
std::vector<const FACET *> facets_to_pointers(const std::vector<FACET> & srcFacets)
{
  std::vector<const FACET *> facetPtrs;
  facetPtrs.reserve(srcFacets.size());
  for (const auto & facet : srcFacets)
    facetPtrs.push_back(&facet);
  return facetPtrs;
}

bool does_file_exist(const std::string& filename)
{
  if (std::ifstream(filename))
    return true;
  return false;
}



TEST(skin_facets,demo)
{
  const std::string inputStlFilename = "geom.stl";
  if (1 != stk::EnvData::parallel_size() || !does_file_exist(inputStlFilename))
    return;

  stk::math::Vector3d scaleVec{1., 1., 1.};
  STLSurface surf("ls", sierra::Diag::sierraTimer(), inputStlFilename, 1, scaleVec);

  using FACET = Facet3d;

  BoundingBox surfBBox = surf.get_bounding_box();

  surf.build_local_facets(surfBBox);

  const std::vector<FACET> & localFacets = surf.get_facets();
  const double facetSize = compute_facet_size(localFacets);
  std::cout << "Facet size " << facetSize << std::endl;

  std::vector<const FACET *> allFacets = facets_to_pointers(localFacets);

  const double defeatureSize = facetSize;

  std::vector<const FACET *> windingNumberFacets = prune_interior_facets(allFacets, surf.get_facets(), defeatureSize, 0.5);
  write_stl("skin.stl", windingNumberFacets);

  const bool computeTighterSkinFacets = false;
  if (computeTighterSkinFacets)
  {
    std::vector<const FACET *> retainedFacets = keep_boundary_facets(windingNumberFacets, surf.get_facets(), defeatureSize);

    write_stl("tighterSkin.stl", retainedFacets);
  }
}

}
