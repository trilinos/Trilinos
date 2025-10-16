#include <Akri_ClosestPointRedistance.hpp>
#include <Akri_ContourElement.hpp>
#include <Akri_Faceted_Surface.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_util/diag/Timer.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalBoundingBox.hpp>
#include <Akri_Sign.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

ClosestPointRedistance::ClosestPointRedistance(const stk::mesh::BulkData & mesh,
    const FieldRef& coordinates,
    const FieldRef& distance)
: ClosestPointRedistance(mesh, coordinates, distance, sierra::Diag::sierraTimer())
{
}

ClosestPointRedistance::ClosestPointRedistance(const stk::mesh::BulkData & mesh,
    const FieldRef& coordinates,
    const FieldRef& distance,
    stk::diag::Timer & parentTimer)
: myMesh(mesh), myCoordinates(coordinates), myDistance(distance), myTimer("Closest Point Redistance", parentTimer)
{
}

void ClosestPointRedistance::build_facets_for_elements(const stk::mesh::BulkData & mesh,
    const FieldRef coordinates,
    const FieldRef distance,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const double elementLengthScale,
    FacetedSurfaceBase & facets)
{
  facets.clear();

  for ( auto && elem : elementsToIntersect )
  {
    ContourElement lsElem( mesh, elem, coordinates, distance );
    lsElem.compute_subelement_decomposition(elementLengthScale);

    lsElem.build_subelement_facets( facets );
  }
}

stk::mesh::Selector ClosestPointRedistance::elements_to_contour_and_redistance(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & activeElementSelector,
    const FieldRef distance)
{
  stk::mesh::Selector selector = activeElementSelector & stk::mesh::selectField(distance) & (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part());
  return selector;
}

void ClosestPointRedistance::build_facets_for_selected_elements(const stk::mesh::BulkData & mesh,
    const FieldRef coordinates,
    const FieldRef distance,
    const stk::mesh::Selector & activeElementSelector,
    const double elementLengthScale,
    FacetedSurfaceBase & facets)
{
  const stk::mesh::Selector contourSelector = elements_to_contour_and_redistance(mesh, activeElementSelector, distance);
  std::vector<stk::mesh::Entity> elems;
  stk::mesh::get_selected_entities(contourSelector, mesh.buckets(stk::topology::ELEMENT_RANK), elems);

  build_facets_for_elements(mesh, coordinates, distance, elems, elementLengthScale, facets);
}

void ClosestPointRedistance::build_isosurface_facets(const stk::mesh::Selector & activeElementSelector, const double elementLengthScale, FacetedSurfaceBase & facets) const
{
  stk::diag::TimeBlock timer__(myTimer);

  build_facets_for_selected_elements(myMesh, myCoordinates, myDistance, activeElementSelector, elementLengthScale, facets);
}

void ClosestPointRedistance::redistance(const stk::mesh::Selector & activeElementSelector, const double elementLengthScale, const double narrowBandDistance, const bool doEnforceSignAwayFromInterface) const
{
  std::unique_ptr<FacetedSurfaceBase> facets = FacetedSurfaceBase::build(myMesh.mesh_meta_data().spatial_dimension());
  build_isosurface_facets(activeElementSelector, elementLengthScale, *facets);

  redistance_using_facets(activeElementSelector, *facets, narrowBandDistance, doEnforceSignAwayFromInterface, elementLengthScale);
}

static double get_sign_change_tolerance(const double doEnforceSignAwayFromInterface, const double elementLengthScale)
{
  // If this is too large, then sharp edges can propagate incorrect signs
  // (even through walls, etc).
  // If this is too small, then a phase can't disappear because the sign
  // preservation will prevent it even if the subelement contouring process
  // neglects it.  So this should be slightly larger than the tolerance in
  // compute_subelement_decomposition.
  if (!doEnforceSignAwayFromInterface)
    return 0.;
  const double signChangeTol = 5.e-4*elementLengthScale;
  return signChangeTol;
}

void ClosestPointRedistance::redistance_using_facets(const stk::mesh::Selector & activeElementSelector,
    FacetedSurfaceBase & facets,
    const double narrowBandDistance,
    const bool doEnforceSignAwayFromInterface,
    const double elementLengthScale) const
{
  stk::diag::TimeBlock timer__(myTimer);

  const stk::mesh::Selector redistanceSelector = elements_to_contour_and_redistance(myMesh, activeElementSelector, myDistance);
  const BoundingBox nodeBBox = compute_nodal_bbox(myMesh, redistanceSelector, myCoordinates);

  facets.prepare_to_compute(nodeBBox, narrowBandDistance);

  redistance_selected_element_nodes_using_prepared_facets(activeElementSelector, facets, narrowBandDistance, doEnforceSignAwayFromInterface, elementLengthScale);
}

static double facet_signed_distance(const FacetedSurfaceBase & facets, const stk::math::Vector3d & x, const int previousSign, const double narrowBandSize, const bool doEnforceSign)
{
  if (doEnforceSign)
  {
    return previousSign * facets.point_unsigned_distance(x, narrowBandSize, narrowBandSize);
  }
  return facets.truncated_point_signed_distance(x, narrowBandSize, previousSign*narrowBandSize);
}

static void redistance_node_using_facets(const FacetedSurfaceBase & facets,
    const stk::math::Vector3d & nodeCoords,
    const double narrowBandDistance,
    const bool doEnforceSignAwayFromInterface,
    const double signChangeTol,
    double & nodeDist)
{
  const int previousSign = sign(nodeDist);
  const bool doEnforceSign = doEnforceSignAwayFromInterface && (std::abs(nodeDist) > signChangeTol);

  nodeDist = facet_signed_distance(facets, nodeCoords, previousSign, narrowBandDistance, doEnforceSign);
}

void ClosestPointRedistance::redistance_given_nodes_using_facets(FacetedSurfaceBase & facets,
    const std::vector<stk::mesh::Entity> & nodesToRedistance,
    const double narrowBandDistance,
    const bool doEnforceSignAwayFromInterface,
    const double elementLengthScale) const
{
  stk::diag::TimeBlock timer__(myTimer);

  const unsigned dim = myMesh.mesh_meta_data().spatial_dimension();
  const BoundingBox nodeBbox = compute_nodal_bbox(myMesh, myCoordinates, nodesToRedistance );
  const double signChangeTol = get_sign_change_tolerance(doEnforceSignAwayFromInterface, elementLengthScale);

  facets.prepare_to_compute(nodeBbox, narrowBandDistance);

  for ( auto && node : nodesToRedistance )
  {
    double & nodeDist = get_scalar_field(myMesh, myDistance, node);
    const stk::math::Vector3d nodeCoords = get_vector_field(myMesh, myCoordinates, node, dim);
    redistance_node_using_facets(facets, nodeCoords, narrowBandDistance, doEnforceSignAwayFromInterface, signChangeTol, nodeDist);
  }
}

void ClosestPointRedistance::redistance_selected_element_nodes_using_prepared_facets(const stk::mesh::Selector & activeElementSelector,
    const FacetedSurfaceBase & facets,
    const double narrowBandDistance,
    const bool doEnforceSignAwayFromInterface,
    const double elementLengthScale) const
{
  const unsigned dim = myMesh.mesh_meta_data().spatial_dimension();
  const stk::mesh::Selector redistanceSelector = elements_to_contour_and_redistance(myMesh, activeElementSelector, myDistance);
  const double signChangeTol = get_sign_change_tolerance(doEnforceSignAwayFromInterface, elementLengthScale);

  for (auto * bucketPtr : myMesh.get_buckets(stk::topology::NODE_RANK, redistanceSelector))
  {
    const double * coordsData = field_data<double>(myCoordinates , *bucketPtr);
    double * distData = field_data<double>(myDistance , *bucketPtr);

    for (size_t i = 0; i < bucketPtr->size(); ++i)
    {
      const stk::math::Vector3d nodeCoords(coordsData+i*dim, dim);
      redistance_node_using_facets(facets, nodeCoords, narrowBandDistance, doEnforceSignAwayFromInterface, signChangeTol, distData[i]);
    }
  }
}

}


