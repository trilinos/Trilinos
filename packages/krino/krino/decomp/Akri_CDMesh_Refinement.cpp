// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_RefinementDrivers.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CDMesh_Refinement.hpp>
#include <Akri_CDMesh_Utils.hpp>
#include <Akri_CurvatureErrorEstimator.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodeMarker.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_RefinementManager.hpp>
#include <Akri_RefinementMarkerUtils.hpp>
#include <Akri_RefinementTreeElementParts.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <Akri_RefinementSupport.hpp>
#include <Akri_Sign.hpp>
#include <Akri_Surface.hpp>

namespace krino {

Refinement_Marker determine_interface_refinement_marker_with_error_estimate(const stk::mesh::BulkData & mesh,
  const bool isElementOnInterface,
  const int interfaceMinRefineLevel,
  const int interfaceMaxRefineLevel,
  const int elementRefineLevel,
  const bool isDefaultCoarsen,
  const CurvatureErrorEstimator & errorEstimator,
  const stk::mesh::Entity elem)

{
  Refinement_Marker marker = determine_refinement_marker(isElementOnInterface, interfaceMinRefineLevel, elementRefineLevel, isDefaultCoarsen);

  if (interfaceMaxRefineLevel > interfaceMinRefineLevel && isElementOnInterface)
  {
    // Error estimate can be expensive so only evaluate it if the marker depends on it
    const bool doRefineIfErrorLargerThanTol = marker != Refinement_Marker::REFINE && elementRefineLevel < interfaceMaxRefineLevel;
    const bool doUnrefineIfErrorSmallerThanTol = marker != Refinement_Marker::COARSEN && elementRefineLevel > interfaceMinRefineLevel && isDefaultCoarsen;
    const bool markerDependsOnErrorEst = doRefineIfErrorLargerThanTol || doUnrefineIfErrorSmallerThanTol;
    if (markerDependsOnErrorEst)
    {
      const double elemErrEst = errorEstimator.estimate_element_error(mesh, elem);

      if (doRefineIfErrorLargerThanTol && elemErrEst > errorEstimator.refine_tolerance())
        marker = Refinement_Marker::REFINE;
      else if (doUnrefineIfErrorSmallerThanTol && elemErrEst < errorEstimator.unrefine_tolerance())
        marker = Refinement_Marker::COARSEN;
    }
  }
  return marker;
}

void
mark_possible_cut_elements_for_adaptivity(const stk::mesh::BulkData& mesh,
      const RefinementManager & refinement,
      const InterfaceGeometry & interfaceGeometry,
      const RefinementSupport & refinementSupport,
      const int refinementIterCount)
{
  const int minRefineLevel = refinementSupport.get_interface_minimum_refinement_level();
  const bool isDefaultCoarsen = should_continue_to_coarsen(refinementIterCount, minRefineLevel);
  const std::vector<stk::mesh::Entity> possiblyCutElems = interfaceGeometry.get_possibly_cut_elements(mesh);
  initialize_marker(mesh, refinement, isDefaultCoarsen);
  mark_elements_to_reach_refinement_level(mesh, refinement, possiblyCutElems, minRefineLevel, isDefaultCoarsen);
}

void
mark_elements_that_intersect_interval(const stk::mesh::BulkData& mesh,
      const RefinementManager & refinement,
      const InterfaceGeometry & interfaceGeometry,
      const std::array<double,2> refinementInterval,
      const int numRefineLevels,
      const bool isDefaultCoarsen)
{
  // Used to refine elements with 1 to many interfaces within the interfaceGeometry that each are to be refined with the same interval to the same refinement level
  // If isDefaultCoarsen=false, then the refinement will be cumulative (since default will be to do nothing)
  initialize_marker(mesh, refinement, isDefaultCoarsen);
  interfaceGeometry.prepare_to_intersect_elements(mesh);
  constexpr bool isDefaultCoarsenForEachLS = false;
  std::vector<stk::mesh::Entity> elementsInInterval;
  for (auto & surfaceIdentifier : interfaceGeometry.get_surface_identifiers())
  {
    interfaceGeometry.fill_elements_that_intersect_distance_interval(mesh, surfaceIdentifier, refinementInterval, elementsInInterval);
    mark_elements_to_reach_refinement_level(mesh, refinement, elementsInInterval, numRefineLevels, isDefaultCoarsenForEachLS);
  }
}

void
mark_interface_elements_for_adaptivity(const stk::mesh::BulkData& mesh,
    const RefinementManager & refinement,
    const InterfaceGeometry & interfaceGeometry,
    const RefinementSupport & refinementSupport,
    const FieldRef coordsField,
    const int numRefinements)
{
  // This refinement strategy cuts elements by the user-specified number of adapt levels
  // before the conformal decomposition.

  const FieldRef elementMarkerField = refinement.get_marker_field_and_sync_to_host();

  const stk::mesh::Selector locallyOwnedNotParent = mesh.mesh_meta_data().locally_owned_part() & !refinement.parent_part();
  const int interfaceMinRefineLevel = refinementSupport.get_interface_minimum_refinement_level();
  const int interfaceMaxRefineLevel = refinementSupport.get_interface_maximum_refinement_level();
  std::vector<stk::mesh::Entity> entities;
  std::vector<stk::mesh::Entity> children;

  const NodeToCapturedDomainsMap nodesToCapturedDomains;
  const std::vector<IntersectionPoint> edgeIntersections = interfaceGeometry.get_edge_intersection_points(mesh, nodesToCapturedDomains);
  const bool isDefaultCoarsen = should_continue_to_coarsen(numRefinements, interfaceMinRefineLevel);

  FieldRef nodeMarkerField = refinementSupport.get_nonconforming_refinement_node_marker_field();
  mark_nodes_near_edge_intersections(mesh, edgeIntersections, nodeMarkerField);

  std::unique_ptr<CurvatureErrorEstimator> errorEstimator;
  if (interfaceMaxRefineLevel > interfaceMinRefineLevel)
    errorEstimator = build_curvature_error_estimator(refinementSupport, interfaceGeometry, refinementSupport.get_timer());
  if (errorEstimator)
    errorEstimator->precompute(mesh, coordsField, !refinement.parent_part(), nodeMarkerField);

  stk::mesh::get_selected_entities( locallyOwnedNotParent, mesh.buckets( stk::topology::ELEMENT_RANK ), entities );

  auto elementMarkerData = elementMarkerField.data<int, stk::mesh::ReadWrite>();
  for( auto&& elem : entities )
  {
    bool hasNearbyCrossing = element_has_marked_node(mesh, elem, nodeMarkerField);

    const int elementRefineLevel = refinement.fully_refined_level(elem);

    Refinement_Marker refMarker;

    if (errorEstimator)
      refMarker = determine_interface_refinement_marker_with_error_estimate(mesh, hasNearbyCrossing, interfaceMinRefineLevel, interfaceMaxRefineLevel, elementRefineLevel, isDefaultCoarsen, *errorEstimator, elem);
    else
      refMarker = determine_refinement_marker(hasNearbyCrossing, interfaceMinRefineLevel, elementRefineLevel, isDefaultCoarsen);

    elementMarkerData.entity_values(elem)() = static_cast<int>(refMarker);
  }

  write_refinement_level_sizes(mesh, refinement, coordsField, entities, interfaceMaxRefineLevel);

}

static unsigned determine_parent_element_part_ordinal_based_on_child_element_part_ordinals(const stk::mesh::BulkData& mesh, const Phase_Support & phaseSupport, const std::set<unsigned> & childElemPartOrdinals)
{
  STK_ThrowRequire(!childElemPartOrdinals.empty());
  const stk::mesh::Part & firstChildElementPart = mesh.mesh_meta_data().get_part(*childElemPartOrdinals.begin());

  if (childElemPartOrdinals.size() > 1 || phaseSupport.is_nonconformal(firstChildElementPart))
  {
    stk::mesh::Part & parentElementPart = phaseSupport.get_nonconforming_part(firstChildElementPart);
    return parentElementPart.mesh_meta_data_ordinal();
  }
  return firstChildElementPart.mesh_meta_data_ordinal();
}

std::vector<std::pair<stk::mesh::Entity,unsigned>> get_owned_adaptivity_parents_and_their_element_part(const stk::mesh::BulkData& mesh, const RefinementManager & refinement, const Phase_Support & phaseSupport)
{
  return get_owned_refinement_parents_and_their_element_part(mesh, refinement,
      [&mesh, &phaseSupport](const std::set<unsigned> & childElemPartOrdinals)
      { return determine_parent_element_part_ordinal_based_on_child_element_part_ordinals(mesh, phaseSupport, childElemPartOrdinals); });
}

}
