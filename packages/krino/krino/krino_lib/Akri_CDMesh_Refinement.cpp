// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AdaptivityHelpers.hpp>
#include <Akri_CDFEM_Snapper.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CDMesh_Refinement.hpp>
#include <Akri_CDMesh_Utils.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_RefinementInterface.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <Akri_RefinementSupport.hpp>

namespace krino {

namespace {
  void set_refine_if_not_reached_max_refine(const RefinementInterface & refinement,
    stk::mesh::Entity elem,
    const int interface_max_refine_level,
    FieldRef elem_marker)
  {
    int & marker = *field_data<int>(elem_marker, elem);
    const int refine_level = refinement.fully_refined_level(elem);

    if (refine_level >= interface_max_refine_level )
    {
      marker = static_cast<int>(Refinement_Marker::NOTHING);
    }
    else
    {
      marker = static_cast<int>(Refinement_Marker::REFINE);
    }
  }
}

static bool node_is_snapped_to_interface(stk::mesh::Entity node,
    const InterfaceID & interface,
    const std::unordered_map<stk::mesh::Entity, std::set<InterfaceID>> & nodesSnappedInterfaces)
{
  auto nodeSnappedInterfaces = nodesSnappedInterfaces.find(node);
  const bool nodeIsSnapped = nodeSnappedInterfaces != nodesSnappedInterfaces.end() &&
      nodeSnappedInterfaces->second.find(interface) != nodeSnappedInterfaces->second.end();
  return nodeIsSnapped;
}

void
refine_edges_with_multiple_unsnapped_crossings(const stk::mesh::BulkData& mesh,
    const RefinementInterface & refinement,
    const std::vector<IntersectionPoint> & edgeIntersections,
    const int interface_max_refine_level,
    FieldRef elem_marker_field,
    const std::unordered_map<stk::mesh::Entity, std::set<InterfaceID>> & nodesSnappedInterfaces)
{
  // Refine any edge with multiple crossings that are not snapped to the nodes
  std::set<std::array<stk::mesh::Entity,2>> edgesWithUnsnappedCrossings;
  for (auto && edgeIntersection : edgeIntersections)
  {
    const EdgeIntersection edge(edgeIntersection);

    const InterfaceID & interface = edge.interface;
    const bool node1IsSnapped = node_is_snapped_to_interface(edge.nodes[0], interface, nodesSnappedInterfaces);
    const bool node2IsSnapped = node_is_snapped_to_interface(edge.nodes[1], interface, nodesSnappedInterfaces);
    if (!node1IsSnapped && !node2IsSnapped)
    {
      const auto insertion = edgesWithUnsnappedCrossings.insert(edge.nodes);
      const bool alreadyInSet = !insertion.second;
      if (alreadyInSet)
      {
        std::vector<stk::mesh::Entity> edge_elems;
        stk::mesh::get_entities_through_relations(mesh, stk::mesh::EntityVector{edge.nodes[0], edge.nodes[1]},
            stk::topology::ELEMENT_RANK, edge_elems);
        for (auto && elem : edge_elems)
        {
          set_refine_if_not_reached_max_refine(refinement, elem, interface_max_refine_level, elem_marker_field);
        }
      }
    }
  }
}

void
refine_edges_with_nodes_with_multiple_snapped_interfaces(const stk::mesh::BulkData& mesh,
    const RefinementInterface & refinement,
    const std::vector<IntersectionPoint> & edgeIntersections,
    const int interface_max_refine_level,
    FieldRef elem_marker_field,
    const std::unordered_map<stk::mesh::Entity, std::set<InterfaceID>> & node_snapped_interfaces)
{
  for (auto && edgeIntersection : edgeIntersections)
  {
    const EdgeIntersection edge(edgeIntersection);

    auto node1_snapped_interfaces = node_snapped_interfaces.find(edge.nodes[0]);
    auto node2_snapped_interfaces = node_snapped_interfaces.find(edge.nodes[1]);
    const unsigned num_node1_interfaces = (node1_snapped_interfaces != node_snapped_interfaces.end()) ?
        node1_snapped_interfaces->second.size() : 0;
    const unsigned num_node2_interfaces = (node2_snapped_interfaces != node_snapped_interfaces.end()) ?
        node2_snapped_interfaces->second.size() : 0;
    if (num_node1_interfaces > 1 || num_node2_interfaces > 1)
    {
      std::vector<stk::mesh::Entity> edge_elems;
      stk::mesh::get_entities_through_relations(mesh, stk::mesh::EntityVector{edge.nodes[0], edge.nodes[1]},
          stk::topology::ELEMENT_RANK, edge_elems);
      for (auto && elem : edge_elems)
      {
        set_refine_if_not_reached_max_refine(refinement, elem, interface_max_refine_level, elem_marker_field);
      }
    }
  }
}

void
determine_which_interfaces_snap_to_each_node_and_unsnappable_nodes(const stk::mesh::BulkData& mesh,
    const RefinementInterface & refinement,
    const std::vector<IntersectionPoint> & edgeIntersections,
    const std::vector<InterfaceID> & active_interface_ids,
    const CDFEM_Snapper & snapper,
    FieldRef node_marker_field,
    std::unordered_map<stk::mesh::Entity, std::set<InterfaceID>> & node_snapped_interfaces,
    std::set<stk::mesh::Entity> & unsnappable_nodes)
{
  for (auto && interface : active_interface_ids)
  {
    stk::mesh::field_fill(0, node_marker_field);
    for (auto && edgeIntersection : edgeIntersections)
    {
      const EdgeIntersection edge(edgeIntersection);

      int & node1_marker = *field_data<int>(node_marker_field, edge.nodes[0]);
      int & node2_marker = *field_data<int>(node_marker_field, edge.nodes[1]);
      if (snapper.snap_lo(edge.crossingLocation))
      {
        if (parts_are_compatible_for_snapping(mesh, edge.nodes[0], edge.nodes[1]))
          node1_marker = 2;
        else
          node1_marker = std::max(node1_marker, 1);
      }
      if (snapper.snap_hi(edge.crossingLocation))
      {
        if (parts_are_compatible_for_snapping(mesh, edge.nodes[1], edge.nodes[0]))
          node2_marker = 2;
        else
          node2_marker = std::max(node2_marker, 1);
      }
    }
    stk::mesh::parallel_max(mesh, {&node_marker_field.field()});

    for(const auto & b_ptr : mesh.buckets( stk::topology::NODE_RANK ))
    {
      stk::mesh::Bucket & b = *b_ptr;
      int * node_marker = field_data<int>(node_marker_field, b);
      for(unsigned i=0; i < b.size(); ++i)
      {
        if (node_marker[i] == 1)
        {
          unsnappable_nodes.insert(b[i]);
        }
        else if (node_marker[i] == 2)
        {
          node_snapped_interfaces[b[i]].insert(interface);
        }
      }
    }
  }
}

void
resolve_fine_features(const stk::mesh::BulkData& mesh,
    const RefinementInterface & refinement,
    const std::vector<Surface_Identifier> & surfaceIdentifiers,
    const std::vector<IntersectionPoint> & edgeIntersections,
    const int interface_max_refine_level,
    FieldRef elem_marker_field,
    FieldRef node_marker_field)
{
  if (!CDFEM_Support::is_active(mesh.mesh_meta_data()))
    return;

  STK_ThrowRequireMsg(false, "Unfinished capability resolve_fine_features");
  const std::vector<InterfaceID> activeInterfaceIds; // = cdmesh.active_interface_ids(surfaceIdentifiers); // FIXME: This should look at the actual edge intersections, right?

  const auto & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());

  std::unordered_map<stk::mesh::Entity, std::set<InterfaceID>> node_snapped_interfaces;
  std::set<stk::mesh::Entity> unsnappable_nodes;

  // What should this look like for regular level set usage?
  determine_which_interfaces_snap_to_each_node_and_unsnappable_nodes(mesh, refinement, edgeIntersections, activeInterfaceIds, cdfemSupport.get_snapper(), node_marker_field, node_snapped_interfaces, unsnappable_nodes);

  // Attempt at minimally aggressive. Works pretty well and has has lower element counts.
  refine_edges_with_multiple_unsnapped_crossings(mesh, refinement, edgeIntersections, interface_max_refine_level, elem_marker_field, node_snapped_interfaces);
  //refine_edges_with_unsnappable_nodes(interface_max_refine_level, elem_marker_field, refine_level_field, transition_element_field, unsnappable_nodes);
  refine_edges_with_nodes_with_multiple_snapped_interfaces(mesh, refinement, edgeIntersections, interface_max_refine_level, elem_marker_field, node_snapped_interfaces);
}

void
mark_nearest_node_on_cut_edges(const stk::mesh::BulkData& mesh,
    const std::vector<IntersectionPoint> & edgeIntersections,
    FieldRef node_marker_field)
{
  const double overlap = 0.25;

  stk::mesh::field_fill(0, node_marker_field);

  for (auto && edgeIntersection : edgeIntersections)
  {
    const EdgeIntersection edge(edgeIntersection);

    if (edge.crossingLocation < 0.5+overlap)
      *field_data<int>(node_marker_field, edge.nodes[0]) = 1;
    if (edge.crossingLocation > 0.5-overlap)
      *field_data<int>(node_marker_field, edge.nodes[1]) = 1;
  }
  stk::mesh::parallel_max(mesh, {&node_marker_field.field()});
}

static void initialize_marker(const stk::mesh::BulkData& mesh,
      const RefinementInterface & refinement,
      const bool isDefaultCoarsen)
{
  const FieldRef elementMarkerField = refinement.get_marker_field_and_sync_to_host();
  const int initialVal = isDefaultCoarsen ? static_cast<int>(Refinement_Marker::COARSEN) : static_cast<int>(Refinement_Marker::NOTHING);
  stk::mesh::field_fill(initialVal, elementMarkerField);
}

int determine_refinement_marker(const bool isElementIndicated, const int interfaceMinRefineLevel, const int elementRefineLevel, const bool isDefaultCoarsen)
{
  auto marker = isDefaultCoarsen ? Refinement_Marker::COARSEN : Refinement_Marker::NOTHING;
  const int targetRefineLevel = isElementIndicated ? interfaceMinRefineLevel : 0;
  if (elementRefineLevel < targetRefineLevel)
    marker = Refinement_Marker::REFINE;
  else if (elementRefineLevel == targetRefineLevel)
    marker = Refinement_Marker::NOTHING;
  return static_cast<int>(marker);
}

static void mark_given_elements(const stk::mesh::BulkData& mesh,
      const RefinementInterface & refinement,
      const std::vector<stk::mesh::Entity> elementsToMark,
      const int minRefineLevel,
      const bool isDefaultCoarsen)
{
  const FieldRef elementMarkerField = refinement.get_marker_field_and_sync_to_host();
  constexpr bool doMarkElement = true;

  for( auto&& elem : elementsToMark )
  {
    int & marker = *field_data<int>(elementMarkerField, elem);
    const int elementRefineLevel = refinement.fully_refined_level(elem);

    marker = determine_refinement_marker(doMarkElement, minRefineLevel, elementRefineLevel, isDefaultCoarsen);
  }
}

static bool should_continue_to_coarsen(const int refinementIterCount, const int targetRefineLevel)
{
  // Stop coarsening after num_levels of refinement to avoid infinite looping
  // from the interface position moving between elements because of snapping changes
  // with refinement
  return refinementIterCount < targetRefineLevel;
}

void
mark_possible_cut_elements_for_adaptivity(const stk::mesh::BulkData& mesh,
      const RefinementInterface & refinement,
      const InterfaceGeometry & interfaceGeometry,
      const RefinementSupport & refinementSupport,
      const int refinementIterCount)
{
  const int minRefineLevel = refinementSupport.get_interface_minimum_refinement_level();
  const bool isDefaultCoarsen = should_continue_to_coarsen(refinementIterCount, minRefineLevel);
  const std::vector<stk::mesh::Entity> possiblyCutElems = interfaceGeometry.get_possibly_cut_elements(mesh);
  initialize_marker(mesh, refinement, isDefaultCoarsen);
  mark_given_elements(mesh, refinement, possiblyCutElems, minRefineLevel, isDefaultCoarsen);
}

void
mark_elements_that_intersect_interval(const stk::mesh::BulkData& mesh,
      const RefinementInterface & refinement,
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
    mark_given_elements(mesh, refinement, elementsInInterval, numRefineLevels, isDefaultCoarsenForEachLS);
  }
}

double compute_edge_length(const FieldRef coordsField, const stk::topology elemTopology, const stk::mesh::Entity * const elemNodes, const unsigned iEdge)
{
  std::array<stk::mesh::Entity,3> edgeNodes;

  elemTopology.edge_nodes(elemNodes, iEdge, edgeNodes.data());

  const stk::math::Vector3d edge_node1_coords(field_data<double>(coordsField, edgeNodes[0]), elemTopology.dimension());
  const stk::math::Vector3d edge_node2_coords(field_data<double>(coordsField, edgeNodes[1]), elemTopology.dimension());
  return (edge_node2_coords-edge_node1_coords).length();
}

void write_refinement_level_sizes(const stk::mesh::BulkData& mesh,
    const RefinementInterface & refinement,
    const FieldRef coordsField,
    const std::vector<stk::mesh::Entity> & elements,
    const int interfaceMaxRefineLevel)
{
  const int maxNumRefinementLevels = interfaceMaxRefineLevel+1;
  std::vector<double> min_edge_lengths(maxNumRefinementLevels, std::numeric_limits<double>::max());
  std::vector<double> max_edge_lengths(maxNumRefinementLevels, 0.0);

  for( auto&& elem : elements )
  {
    if (!refinement.is_transition(elem))
    {
      const int elementRefineLevel = refinement.fully_refined_level(elem);
      STK_ThrowRequire(elementRefineLevel < maxNumRefinementLevels);
      const stk::topology elemTopology = mesh.bucket(elem).topology();
      const stk::mesh::Entity * const elemNodes = mesh.begin(elem, stk::topology::NODE_RANK);

      for (unsigned iEdge = 0; iEdge < elemTopology.num_edges(); ++iEdge)
      {
        const double length = compute_edge_length(coordsField, elemTopology, elemNodes, iEdge);

        min_edge_lengths[elementRefineLevel] = std::min(min_edge_lengths[elementRefineLevel], length);
        max_edge_lengths[elementRefineLevel] = std::max(max_edge_lengths[elementRefineLevel], length);
      }
    }
  }

  std::vector<double> localVals = min_edge_lengths;
  stk::all_reduce_min(mesh.parallel(), localVals.data(), min_edge_lengths.data(), min_edge_lengths.size());
  localVals = max_edge_lengths;
  stk::all_reduce_max(mesh.parallel(), localVals.data(), max_edge_lengths.data(), max_edge_lengths.size());

  for (int i=0; i<maxNumRefinementLevels; ++i)
  {
    krinolog << "Min and Max sizes for refinement level " << i <<  " " << min_edge_lengths[i] << " " << max_edge_lengths[i] << stk::diag::dendl;
  }
}

bool element_has_marked_node(const stk::mesh::BulkData& mesh, stk::mesh::Entity elem, const FieldRef nodeMarkerField)
{
  const unsigned num_nodes = mesh.num_nodes(elem);
  const stk::mesh::Entity * const elemNodes = mesh.begin_nodes(elem);
  for(unsigned i=0; i < num_nodes; ++i)
  {
    const int nodeMarker = *field_data<int>(nodeMarkerField, elemNodes[i]);
    if(nodeMarker)
      return true;
  }
  return false;
}

void
mark_interface_elements_for_adaptivity(const stk::mesh::BulkData& mesh,
    const RefinementInterface & refinement,
    const InterfaceGeometry & interfaceGeometry,
    const RefinementSupport & refinementSupport,
    const FieldRef coordsField,
    const int numRefinements)
{
  // This refinement strategy cuts elements by the user-specified number of adapt levels
  // before the conformal decomposition.

  const FieldRef elementMarkerField = refinement.get_marker_field_and_sync_to_host();

  const stk::mesh::Selector locally_owned_selector(mesh.mesh_meta_data().locally_owned_part());
  const int interfaceMinRefineLevel = refinementSupport.get_interface_minimum_refinement_level();
  const int interfaceMaxRefineLevel = refinementSupport.get_interface_maximum_refinement_level();
  std::vector<stk::mesh::Entity> entities;
  std::vector<stk::mesh::Entity> children;

  const NodeToCapturedDomainsMap nodesToCapturedDomains;
  const std::vector<IntersectionPoint> edgeIntersections = interfaceGeometry.get_edge_intersection_points(mesh, nodesToCapturedDomains);
  const bool isDefaultCoarsen = should_continue_to_coarsen(numRefinements, interfaceMinRefineLevel);

  FieldRef nodeMarkerField = refinementSupport.get_nonconforming_refinement_node_marker_field();
  mark_nearest_node_on_cut_edges(mesh, edgeIntersections, nodeMarkerField);

  stk::mesh::get_selected_entities( locally_owned_selector, mesh.buckets( stk::topology::ELEMENT_RANK ), entities );
  for( auto&& elem : entities )
  {
    bool hasCrossing = element_has_marked_node(mesh, elem, nodeMarkerField);

    const int elementRefineLevel = refinement.fully_refined_level(elem);

    int & marker = *field_data<int>(elementMarkerField, elem);
    marker = determine_refinement_marker(hasCrossing, interfaceMinRefineLevel, elementRefineLevel, isDefaultCoarsen);
  }

  write_refinement_level_sizes(mesh, refinement, coordsField, entities, interfaceMaxRefineLevel);

  if (interfaceMinRefineLevel > interfaceMaxRefineLevel)
  {
    resolve_fine_features(mesh, refinement, interfaceGeometry.get_surface_identifiers(), edgeIntersections, interfaceMaxRefineLevel, elementMarkerField, nodeMarkerField);
  }
}

void
refine_edges_with_unsnappable_nodes(const stk::mesh::BulkData& mesh,
    const RefinementInterface & refinement,
    const std::vector<IntersectionPoint> & edgeIntersections,
    const CDFEM_Snapper & snapper,
    const int interface_max_refine_level,
    FieldRef elem_marker_field,
    std::set<stk::mesh::Entity> unsnappable_nodes)
{
  for (auto && edgeIntersection : edgeIntersections)
  {
    const EdgeIntersection edge(edgeIntersection);

    if ((snapper.snap_lo(edge.crossingLocation) && unsnappable_nodes.find(edge.nodes[0]) != unsnappable_nodes.end()) ||
        (snapper.snap_hi(edge.crossingLocation) && unsnappable_nodes.find(edge.nodes[1]) != unsnappable_nodes.end()))
    {
      krinolog << "Refining unsnappable edge  " << mesh.identifier(edge.nodes[0]) << " " << mesh.identifier(edge.nodes[1]) << " " << debug_output(mesh, edgeIntersection) << stk::diag::dendl;
      std::vector<stk::mesh::Entity> edge_elems;
      stk::mesh::get_entities_through_relations(mesh, stk::mesh::EntityVector{edge.nodes[0], edge.nodes[1]},
          stk::topology::ELEMENT_RANK, edge_elems);
      for (auto && elem : edge_elems)
      {
        set_refine_if_not_reached_max_refine(refinement, elem, interface_max_refine_level, elem_marker_field);
      }
    }
  }
}

static unsigned determine_parent_element_part_ordinal_based_on_child_element_part_ordinals(const stk::mesh::BulkData& mesh, const Phase_Support & phaseSupport, const std::set<unsigned> & childElemPartOrdinals)
{
  STK_ThrowRequire(!childElemPartOrdinals.empty());
  const stk::mesh::Part & firstChildElementPart = mesh.mesh_meta_data().get_part(*childElemPartOrdinals.begin());

  if (childElemPartOrdinals.size() > 1 || phaseSupport.is_nonconformal(&firstChildElementPart))
  {
    const stk::mesh::Part * parentElementPart = phaseSupport.find_nonconformal_part(firstChildElementPart);
    STK_ThrowAssert(nullptr != parentElementPart);
    return parentElementPart->mesh_meta_data_ordinal();
  }
  return firstChildElementPart.mesh_meta_data_ordinal();
}

static void fill_sorted_owned_leaf_children_and_their_element_part_and_their_unique_parents(const stk::mesh::BulkData& mesh,
    const RefinementInterface & refinement,
    std::vector<std::pair<stk::mesh::EntityId,unsigned>> & sortedChildIdsAndElementPart,
    std::vector<stk::mesh::Entity> & nextIterationOfParentsToProcess)
{
  const stk::mesh::Selector adaptivityLeafChild = mesh.mesh_meta_data().locally_owned_part() & refinement.child_part() & !refinement.parent_part();
  const auto & buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, adaptivityLeafChild);

  for( auto&& bucket : buckets )
  {
    for (auto && elem : *bucket)
    {
      if (refinement.is_child(elem))  //FIXME: shouldnt be needed, but it is in percept right now
      {
        sortedChildIdsAndElementPart.emplace_back(mesh.identifier(elem), find_element_part(mesh, elem).mesh_meta_data_ordinal());
        nextIterationOfParentsToProcess.push_back(refinement.get_parent(elem));
      }
    }
  }
  std::sort(sortedChildIdsAndElementPart.begin(), sortedChildIdsAndElementPart.end());
  stk::util::sort_and_unique(nextIterationOfParentsToProcess);
}

static const unsigned * find_child_element_part_ordinal(const std::vector<std::pair<stk::mesh::EntityId,unsigned>> & sortedChildIdsAndElementPart, stk::mesh::EntityId childId)
{
  auto compare = [](const std::pair<stk::mesh::EntityId,unsigned> & idAndElementPart, const stk::mesh::EntityId id) { return idAndElementPart.first < id; };
  auto it = std::lower_bound(sortedChildIdsAndElementPart.begin(), sortedChildIdsAndElementPart.end(), childId, compare);
  if (it != sortedChildIdsAndElementPart.end() && it->first == childId)
    return &(it->second);
  return nullptr;
}

static std::set<unsigned> get_child_element_part_ordinals(const std::vector<std::pair<stk::mesh::EntityId,unsigned>> & sortedChildIdsAndElementPart,
    const std::vector<stk::mesh::EntityId> & childElemIds)
{
  std::set<unsigned> childElemPartOrdinals;
  for (auto&& childElemId : childElemIds)
  {
    const unsigned * childElemOrd = find_child_element_part_ordinal(sortedChildIdsAndElementPart, childElemId);
    if (nullptr == childElemOrd)
    {
      // not ready to process, some children not yet available
      childElemPartOrdinals.clear();
      break;
    }
    childElemPartOrdinals.insert(*childElemOrd);
  }
  return childElemPartOrdinals;
}

struct RelocatedChildElemWithElemPart {
  int parentElemOwningProc;
  stk::mesh::EntityId childElemId;
  stk::mesh::EntityId parentElemId;
  unsigned childElemPartOrdinal;
};

static
void pack_child_element_data_for_parent_element_processor(const stk::mesh::BulkData & mesh,
    const std::vector<RelocatedChildElemWithElemPart> & childElemDataToCommunicate,
    stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto && data : childElemDataToCommunicate)
    {
      commSparse.send_buffer(data.parentElemOwningProc).pack(data.childElemId);
      commSparse.send_buffer(data.parentElemOwningProc).pack(data.parentElemId);
      commSparse.send_buffer(data.parentElemOwningProc).pack(data.childElemPartOrdinal);
    }
  });
}

static
void unpack_child_element_data(const stk::mesh::BulkData & mesh,
    std::vector<std::pair<stk::mesh::EntityId,unsigned>> & iterationChildIdsAndElementPart,
    std::vector<stk::mesh::Entity> & parentsToProcess,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId childElemId;
      commSparse.recv_buffer(procId).unpack(childElemId);
      stk::mesh::EntityId parentElemId;
      commSparse.recv_buffer(procId).unpack(parentElemId);
      unsigned childElemPartOrdinal;
      commSparse.recv_buffer(procId).unpack(childElemPartOrdinal);
      iterationChildIdsAndElementPart.emplace_back(childElemId, childElemPartOrdinal);
      const stk::mesh::Entity parentElem = mesh.get_entity(stk::topology::ELEMENT_RANK, parentElemId);
      parentsToProcess.push_back(parentElem);
    }
  });
}

static void communicate_elements_to_parent_element_processor(const stk::mesh::BulkData& mesh,
    const std::vector<RelocatedChildElemWithElemPart> & childElemDataToCommunicate,
    std::vector<std::pair<stk::mesh::EntityId,unsigned>> & iterationChildIdsAndElementPart,
    std::vector<stk::mesh::Entity> & parentsToProcess)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_child_element_data_for_parent_element_processor(mesh, childElemDataToCommunicate, commSparse);
  unpack_child_element_data(mesh, iterationChildIdsAndElementPart, parentsToProcess, commSparse);
}

static void update_children_and_their_element_part_and_their_unique_parents(const stk::mesh::BulkData& mesh,
    const RefinementInterface & refinement,
    const std::vector<std::pair<stk::mesh::Entity,unsigned>> & iterationChildrenAndElementPart,
    std::vector<std::pair<stk::mesh::EntityId,unsigned>> & sortedChildIdsAndElementPart,
    std::vector<stk::mesh::Entity> & parentsToProcess)
{
  parentsToProcess.clear();

  std::vector<std::pair<stk::mesh::EntityId,unsigned>> iterationChildIdsAndElementPart;
  iterationChildIdsAndElementPart.reserve(iterationChildrenAndElementPart.size());

  std::vector<RelocatedChildElemWithElemPart> childElemDataToCommunicate;

  for (auto && childAndElementPart : iterationChildrenAndElementPart)
  {
    const stk::mesh::Entity child = childAndElementPart.first;
    const stk::mesh::EntityId childId = mesh.identifier(child);
    iterationChildIdsAndElementPart.emplace_back(childId, childAndElementPart.second);
    const auto parentIdAndOwnerRank = refinement.get_parent_id_and_parallel_owner_rank(child);
    if (parentIdAndOwnerRank.second == mesh.parallel_rank())
      parentsToProcess.push_back(mesh.get_entity(stk::topology::ELEMENT_RANK, parentIdAndOwnerRank.first));
    else
      childElemDataToCommunicate.push_back(RelocatedChildElemWithElemPart{parentIdAndOwnerRank.second, childId, parentIdAndOwnerRank.first, childAndElementPart.second});
  }

  communicate_elements_to_parent_element_processor(mesh, childElemDataToCommunicate, iterationChildIdsAndElementPart, parentsToProcess);

  std::sort(iterationChildIdsAndElementPart.begin(), iterationChildIdsAndElementPart.end());
  stk::util::insert_keep_sorted(iterationChildIdsAndElementPart, sortedChildIdsAndElementPart, std::less<std::pair<stk::mesh::EntityId,unsigned>>());
  stk::util::sort_and_unique(parentsToProcess);
}

std::vector<std::pair<stk::mesh::Entity,unsigned>> get_owned_adaptivity_parents_and_their_element_part(const stk::mesh::BulkData& mesh, const RefinementInterface & refinement, const Phase_Support & phaseSupport)
{
  // This used to be pretty simple.  We just recursed down to the leaves of the tree and checked the block
  // parts of the conforming elements.  But now the tree may be spread across processors, so we need
  // to communicate this information up the tree.

  std::vector<std::pair<stk::mesh::Entity,unsigned>> parentsAndElementPart;

  std::vector<std::pair<stk::mesh::EntityId,unsigned>> sortedChildIdsAndElementPart;
  std::vector<stk::mesh::Entity> parentsToProcess;

  fill_sorted_owned_leaf_children_and_their_element_part_and_their_unique_parents(mesh, refinement, sortedChildIdsAndElementPart, parentsToProcess);

  std::vector<stk::mesh::EntityId> childElemIds;
  std::vector<std::pair<stk::mesh::Entity,unsigned>> iterationChildrenAndElementPart;

  while (stk::is_true_on_any_proc(mesh.parallel(), !parentsToProcess.empty()))
  {
    iterationChildrenAndElementPart.clear();

    for( auto&& parent : parentsToProcess )
    {
      refinement.fill_child_element_ids(parent, childElemIds);
      const std::set<unsigned> childElemPartOrdinals = get_child_element_part_ordinals(sortedChildIdsAndElementPart, childElemIds);
      if (!childElemPartOrdinals.empty())
      {
        const unsigned elementPartOrdinal = determine_parent_element_part_ordinal_based_on_child_element_part_ordinals(mesh, phaseSupport, childElemPartOrdinals);
        parentsAndElementPart.emplace_back(parent, elementPartOrdinal);

        if (refinement.is_child(parent))
          iterationChildrenAndElementPart.emplace_back(parent, elementPartOrdinal);
      }
    }

    update_children_and_their_element_part_and_their_unique_parents(mesh, refinement, iterationChildrenAndElementPart, sortedChildIdsAndElementPart, parentsToProcess);
  }
  return parentsAndElementPart;
}

}
