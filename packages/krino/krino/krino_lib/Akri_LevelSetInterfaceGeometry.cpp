// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Element_Cutter.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_ElementCutterUtils.hpp>
#include <Akri_LevelSetInterfaceGeometry.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_SnapInfo.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <Akri_MasterElementDeterminer.hpp>
#include <Akri_Surface_Identifier.hpp>

namespace krino {

std::unique_ptr<Element_Cutter> create_element_cutter_with_error_handling(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges,
    const std::vector<bool> & areParentEdgesAreOrientedSameAsElementEdges,
    const Phase_Support & phaseSupport,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker)
{
  if (krinolog.shouldPrint(LOG_DEBUG))
    krinolog << "Building cutting planes for element global_id=" << mesh.identifier(element) << "\n";

  stk::topology elementTopology = mesh.bucket(element).topology();
  const MasterElement & masterElement = MasterElementDeterminer::getMasterElement(elementTopology);
  const bool oneLSPerPhase = phaseSupport.has_one_levelset_per_phase();

  std::unique_ptr<Element_Cutter> cutter;

  try
  {
    cutter = create_element_cutter(oneLSPerPhase, masterElement, elementParentEdges, areParentEdgesAreOrientedSameAsElementEdges, intersectingPlanesDiagonalPicker);
  }
  catch(const std::exception & err)
  {
    krinolog << "Error constructing cutting surfaces for Mesh_Element " << mesh.identifier(element) << stk::diag::push << stk::diag::dendl;
    for(unsigned i=0; i < elementParentEdges.size(); ++i)
    {
      krinolog << "Edge " << i;
      if(elementParentEdges[i]) {
        krinolog << " node ids = ";
        for(auto && node : elementParentEdges[i]->get_nodes()) krinolog << mesh.identifier(node) << " ";
        krinolog << stk::diag::push << stk::diag::dendl;
        krinolog << *elementParentEdges[i] << stk::diag::pop << stk::diag::dendl;
      }
      else
      {
        krinolog << " No parent edge." << stk::diag::dendl;
      }
    }
    krinolog << stk::diag::pop << stk::diag::dendl;
    krinolog << err.what() << stk::diag::dendl;
    throw err;
  }

  return cutter;
}

static stk::math::Vector3d get_centroid(const std::vector<stk::math::Vector3d> & elemNodesCoords)
{
  stk::math::Vector3d centroid = stk::math::Vector3d::ZERO;
  for(auto && nodeCoords : elemNodesCoords)
  {
    centroid += nodeCoords;
  }
  centroid *= 1./elemNodesCoords.size();
  return centroid;
}

LevelSetElementCutter::LevelSetElementCutter(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const ParentEdgeMap & parentEdges,
    const Phase_Support & phaseSupport,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker)
{
  fill_element_parent_edges(mesh, element, parentEdges, myParentEdges, myParentEdgesAreOrientedSameAsElementEdges);
  myElementInterfaceCutter = create_element_cutter_with_error_handling(mesh, element, myParentEdges, myParentEdgesAreOrientedSameAsElementEdges, phaseSupport, intersectingPlanesDiagonalPicker);
}

int LevelSetElementCutter::interface_sign_for_uncrossed_element(const InterfaceID interface, const std::vector<stk::math::Vector3d> & elemNodesCoords) const
{
  return myElementInterfaceCutter->sign_at_position(interface, get_centroid(elemNodesCoords));
}

std::pair<int, double> LevelSetElementCutter::interface_edge_crossing_sign_and_position(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const
{
  const double sign0 = myElementInterfaceCutter->sign_at_position(interface, edgeNodeCoords[0]);
  const double sign1 = myElementInterfaceCutter->sign_at_position(interface, edgeNodeCoords[1]);
  if (sign0 == -sign1)
    return {sign1, myElementInterfaceCutter->interface_crossing_position(interface, edgeNodeCoords)};
  return {0, -1.};
}

std::string LevelSetElementCutter::visualize(const stk::mesh::BulkData & mesh) const
{
  std::ostringstream os;

  std::vector<ElementIntersection> interiorIntersections;
  myElementInterfaceCutter->fill_interior_intersections(interiorIntersections);

  os << "Interior intersections " << interiorIntersections.size() << ": \n";
  for (auto && intersection : interiorIntersections)
  {
    os << "Interior intersection at " << intersection.parametricCoords << " with domains { ";
    for (int domain : intersection.sortedDomains) os<< domain << " ";
        os << "}\n";
  }

  os << "Parent edges: \n";
  for (auto && edge : myParentEdges)
  {
    os << "Edge with nodes " << mesh.identifier(edge->get_parent_nodes().first) << " and " << mesh.identifier(edge->get_parent_nodes().second) << ":\n";
    if (edge) os << *edge;
  }

  os << myElementInterfaceCutter->visualize() << "\n";

  return os.str();
}

std::unique_ptr<ElementCutter> LevelSetInterfaceGeometry::build_element_cutter(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker) const
{
  std::unique_ptr<ElementCutter> cutter;
  cutter.reset( new LevelSetElementCutter(mesh, element, myParentEdges, myPhaseSupport, intersectingPlanesDiagonalPicker) );
  return cutter;
}

PhaseTag LevelSetInterfaceGeometry::get_starting_phase(const ElementCutter * cutter) const
{
  const LevelSetElementCutter * LSCutter = dynamic_cast<const LevelSetElementCutter *>(cutter);
  STK_ThrowRequire(LSCutter);
  return LSCutter->get_starting_phase(myPhaseSupport, myLSFields);
}

PhaseTag LevelSetElementCutter::get_starting_phase(const Phase_Support & phaseSupport, const std::vector<LS_Field> & LSFields) const
{
  PhaseTag phase;

  // For elements with no interfaces, this will be the uncrossed phase of the mesh element.
  // For elements with interface, this is the starting phase that will be inherited by the
  // subelements and then incrementally updated as we process the interfaces.
  if (LSFields.size() > 1 && myElementInterfaceCutter->is_one_ls_per_phase())
  {
    int LSphase = myElementInterfaceCutter->get_starting_phase_for_cutting_surfaces();
    if (LSphase == -1)
    {
      const auto & phasesPresent = get_phases_present_on_edges(myParentEdges);
      STK_ThrowRequire(!phasesPresent.empty() && (myElementInterfaceCutter->get_num_cutting_surfaces() > 0 || 1 == phasesPresent.size()));
      LSphase = *phasesPresent.begin();
    }
    phase.add(LSFields[LSphase].identifier, -1);
  }
  else
  {
    const int num_ls = LSFields.size();
    for(int ls_index=0; ls_index < num_ls; ++ls_index)
    {
      const InterfaceID interface(ls_index,ls_index);
      int sign = 0;
      bool shouldSetPhase = false;
      for(auto && edge : myParentEdges)
      {
        if (edge)
        {
          shouldSetPhase = true;
          if (edge->have_crossing(InterfaceID(ls_index,ls_index)))
          {
            shouldSetPhase = false;
            break;
          }
          else
            sign = edge->get_crossing_sign(interface);
        }
      }
      if (shouldSetPhase)
        phase.add(LSFields[ls_index].identifier, sign);
    }
  }

  return phase;
}

static std::vector<int> find_next_phase_candidates(const std::vector<InterfaceID> & interfaces,
    const std::vector<int> & pathSoFar)
{
  const int currentPhase = pathSoFar.back();
  std::vector<int> nextPhaseCandidates;
  for (auto && interface : interfaces)
  {
    if (interface.first_ls() == currentPhase)
    {
      if (std::find(pathSoFar.begin(), pathSoFar.end(), interface.second_ls()) == pathSoFar.end())
        nextPhaseCandidates.push_back(interface.second_ls());
    }
    else if (interface.second_ls() == currentPhase)
    {
      if (std::find(pathSoFar.begin(), pathSoFar.end(), interface.first_ls()) == pathSoFar.end())
        nextPhaseCandidates.push_back(interface.first_ls());
    }
  }
  return nextPhaseCandidates;
}

static std::vector<int> shortest_path_to_end(const std::vector<int> & pathSoFar,
    const std::vector<InterfaceID> & interfaces,
    const std::set<int> & endPhases)
{
  if (endPhases.count(pathSoFar.back()) > 0)
    return pathSoFar;

  const std::vector<int> nextPhaseCandidates = find_next_phase_candidates(interfaces, pathSoFar);
  if (nextPhaseCandidates.empty())
  {
    return {};
  }

  std::vector<int> shortestPath;
  size_t shortestPathSize = std::numeric_limits<size_t>::max();
  for (int nextPhase : nextPhaseCandidates)
  {
    std::vector<int> path = pathSoFar;
    path.push_back(nextPhase);
    const auto fullPath = shortest_path_to_end(path, interfaces, endPhases);
    if (!fullPath.empty() && fullPath.size() < shortestPathSize)
    {
      shortestPath = fullPath;
      shortestPathSize = fullPath.size();
    }
  }
  return shortestPath;
}

std::vector<int> shortest_path_from_begin_to_end(const std::vector<InterfaceID> & interfaces,
    const int beginPhase,
    const std::set<int> & endPhases)
{
  std::vector<int> startPath;
  startPath.push_back(beginPhase);
  return shortest_path_to_end(startPath, interfaces, endPhases);
}

static bool captures_interface(const std::vector<int> * sortedDomains, const InterfaceID & interface)
{
  if (sortedDomains == nullptr)
    return false;
  if (interface.first_ls() == interface.second_ls())
    return first_sorted_vector_of_domains_contains_all_domains_in_second_vector(*sortedDomains, {interface.first_ls()});
  return first_sorted_vector_of_domains_contains_all_domains_in_second_vector(*sortedDomains, {interface.first_ls(),interface.second_ls()});
}

static bool interface_has_uncaptured_edge_intersection(const LevelSetElementCutter & cutter,
    const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains,
    const InterfaceID & interface)
{
  STK_ThrowRequire(elemNodesCoords.size() == 4 || elemNodesCoords.size() == 3);
  const stk::topology topology = (elemNodesCoords.size() == 4)? stk::topology::TETRAHEDRON_4 : stk::topology::TRIANGLE_3_2D;

  const unsigned numEdges = topology.num_edges();
  for(unsigned i=0; i < numEdges; ++i)
  {
    const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(topology, i);
    const int n0 = edgeNodeOrdinals[0];
    const int n1 = edgeNodeOrdinals[1];
    if (!captures_interface(elemNodesSnappedDomains[n0], interface) && !captures_interface(elemNodesSnappedDomains[n1], interface))
    {
      const std::array<stk::math::Vector3d,2> edgeNodeCoords{elemNodesCoords[n0], elemNodesCoords[n1]};
      if (cutter.have_crossing(interface, edgeNodeCoords))
        return true;
    }
  }
  return false;
}

static std::vector<InterfaceID>
get_sorted_cutting_interfaces_with_uncaptured_intersections(const LevelSetElementCutter & cutter,
    const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains)
{
  std::set<InterfaceID> interfacesWithUncapturedCrossings;
    for (auto && interface : cutter.get_sorted_cutting_interfaces())
      if (interface_has_uncaptured_edge_intersection(cutter, elemNodesCoords, elemNodesSnappedDomains, interface))
        interfacesWithUncapturedCrossings.insert(interface);

  cutter.add_interfaces_with_uncaptured_intersection_within_element(elemNodesCoords, elemNodesSnappedDomains, interfacesWithUncapturedCrossings);

  std::vector<InterfaceID> interfaces(interfacesWithUncapturedCrossings.begin(), interfacesWithUncapturedCrossings.end());
  return interfaces;
}

static int get_interface_index(const std::vector<InterfaceID> & sortedInterfaces, const InterfaceID interface)
{
  const auto iter = std::lower_bound(sortedInterfaces.begin(), sortedInterfaces.end(), interface);
  return std::distance(sortedInterfaces.begin(), iter);
}

std::vector<int>
LevelSetElementCutter::get_interface_signs_based_on_crossings(const std::vector<stk::math::Vector3d> & elemNodesCoords,
    const std::vector<const std::vector<int> *> & elemNodesSnappedDomains) const
{
  const auto allInterfaces = get_sorted_cutting_interfaces();
  std::vector<int> interfaceSigns(allInterfaces.size(), 0);

  const auto intersectingInterfaces = get_sorted_cutting_interfaces_with_uncaptured_intersections(*this, elemNodesCoords, elemNodesSnappedDomains);
  const stk::math::Vector3d centroid = get_centroid(elemNodesCoords);

  const bool oneLSPerPhase = myElementInterfaceCutter->is_one_ls_per_phase();
  if (oneLSPerPhase)
  {
    std::set<int> subPhases;
    const int ownerStartPhase = get_starting_phase_for_cutting_surfaces();

    for (auto && interface : allInterfaces)
    {
      if (std::binary_search(intersectingInterfaces.begin(), intersectingInterfaces.end(), interface))
      {
        subPhases.insert(interface.first_ls());
        subPhases.insert(interface.second_ls());
      }
      else
      {
        interfaceSigns[get_interface_index(allInterfaces, interface)] = -2;
      }
    }
    if (subPhases.empty())
    {
      subPhases.insert(myElementInterfaceCutter->get_ls_per_interface_phase_at_location(centroid));
    }

    if (subPhases.count(ownerStartPhase) == 0)
    {
      std::vector<int> fixPath = shortest_path_from_begin_to_end(allInterfaces, ownerStartPhase, subPhases);
      STK_ThrowRequireMsg(!fixPath.empty(), "Cannot fix starting phase.");
      for (unsigned i=1; i<fixPath.size(); ++i)
      {
        const InterfaceID interface(fixPath[i-1], fixPath[i]);
        const int sign = (fixPath[i] > fixPath[i-1]) ? 1 : -1;
        interfaceSigns[get_interface_index(allInterfaces, interface)] = sign;
      }
    }
  }
  else
  {
    for (auto && interface : allInterfaces)
      if (!std::binary_search(intersectingInterfaces.begin(), intersectingInterfaces.end(), interface))
        interfaceSigns[get_interface_index(allInterfaces, interface)] = sign_at_position(interface, centroid);
  }

  return interfaceSigns;
}

void LevelSetElementCutter::update_edge_crossings(const unsigned iEdge, const std::vector<std::vector<double>> & nodesIsovar)
{
  STK_ThrowRequire(iEdge < myParentEdges.size());
  CDFEM_Parent_Edge * parentEdge = const_cast<CDFEM_Parent_Edge *>(myParentEdges[iEdge]);
  STK_ThrowRequire(parentEdge);
  const bool oneLSPerPhase = myElementInterfaceCutter->is_one_ls_per_phase();
  if (myParentEdgesAreOrientedSameAsElementEdges[iEdge])
  {
    parentEdge->find_crossings(oneLSPerPhase, nodesIsovar);
  }
  else
  {
    const std::vector<std::vector<double>> orientedIsovar = {nodesIsovar[1], nodesIsovar[0]};
    parentEdge->find_crossings(oneLSPerPhase, orientedIsovar);
  }
}

static std::function<bool(const std::array<unsigned,4> &)>
temporary_build_always_true_diagonal_picker()
{
  auto diagonalPicker =
  [](const std::array<unsigned,4> & faceNodes)
  {
    return true;
  };
  return diagonalPicker;
}

typedef std::function<bool(const CDFEM_Parent_Edge & edge)> ParentEdgeFilter;

static ParentEdgeFilter keep_owned_edges_filter(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & parentElementSelector)
{
  const auto filter = [&mesh, &parentElementSelector](const CDFEM_Parent_Edge & edge)
  {
    const std::pair<stk::mesh::Entity,stk::mesh::Entity> edgeNodes = edge.get_parent_nodes();
    std::vector<stk::mesh::Entity> edgeElems;
    stk::mesh::get_entities_through_relations(mesh, stk::mesh::EntityVector{edgeNodes.first, edgeNodes.second}, stk::topology::ELEMENT_RANK, edgeElems);
    {
    bool foundOwnedElement = false;
    for (auto && edgeElem : edgeElems)
      if (parentElementSelector(mesh.bucket(edgeElem)) && mesh.parallel_owner_rank(edgeElem) == mesh.parallel_rank())
        foundOwnedElement = true;
    STK_ThrowRequire(foundOwnedElement);
    }
    const int parallelRank = mesh.parallel_rank(); // Assumes local proc owns at least one selected element of edge
    for (auto && edgeElem : edgeElems)
      if (mesh.parallel_owner_rank(edgeElem) < parallelRank && parentElementSelector(mesh.bucket(edgeElem)))
        return false;
    return true;
  };
  return filter;
}

static ParentEdgeFilter keep_all_edges_filter()
{
  const auto filter = [](const CDFEM_Parent_Edge & edge)
  {
    return true;
  };
  return filter;
}

static void append_intersection_points_from_filtered_parent_edges(std::vector<IntersectionPoint> & intersectionPoints,
    const ParentEdgeMap & parentEdges,
    const IntersectionPointFilter & intersectionPointFilter,
    const ParentEdgeFilter edgeFilter)
{
  const bool intersectionPointIsOwned = true;
  std::vector<int> intersectionPointSortedDomains;
  for (auto && mapEntry : parentEdges)
  {
    const CDFEM_Parent_Edge & edge = mapEntry.second;
    const auto & edgeCrossings = edge.get_crossings();
    if (!edgeCrossings.empty() && edgeFilter(edge))
    {
      const std::pair<stk::mesh::Entity,stk::mesh::Entity> edgeNodes = edge.get_parent_nodes();
      for (auto & crossing : edgeCrossings)
      {
        const InterfaceID & interface = crossing.first;
        const double location = crossing.second;
        interface.fill_sorted_domains(intersectionPointSortedDomains);
        std::vector<stk::mesh::Entity> intersectionPointNodes{edgeNodes.first, edgeNodes.second};
        if (intersectionPointFilter(intersectionPointNodes, intersectionPointSortedDomains))
          intersectionPoints.emplace_back(intersectionPointIsOwned, intersectionPointNodes, std::vector<double>{1.-location, location}, intersectionPointSortedDomains);
      }
    }
  }
}

static void append_intersection_points_from_owned_parent_edges(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & parentElementSelector,
    std::vector<IntersectionPoint> & intersectionPoints,
    const ParentEdgeMap & parentEdges,
    const IntersectionPointFilter & intersectionPointFilter)
{
  append_intersection_points_from_filtered_parent_edges(intersectionPoints, parentEdges, intersectionPointFilter, keep_owned_edges_filter(mesh, parentElementSelector));
}

static void append_intersection_points_from_all_parent_edges(std::vector<IntersectionPoint> & intersectionPoints,
    const ParentEdgeMap & parentEdges,
    const IntersectionPointFilter & intersectionPointFilter)
{
  append_intersection_points_from_filtered_parent_edges(intersectionPoints, parentEdges, intersectionPointFilter, keep_all_edges_filter());
}

static bool determine_if_mesh_has_higher_order_midside_nodes_with_level_set_locally(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields)
{
  // Assume we can check master element on one bucket and get the right answer

  for ( auto && bucket_ptr : mesh.get_buckets(stk::topology::ELEMENT_RANK, phaseSupport.get_all_decomposed_blocks_selector()) )
  {
    if (bucket_ptr->topology() != bucket_ptr->topology().base())
    {
      for (auto && LSField : LSFields)
      {
        stk::topology fieldTopology = MasterElementDeterminer::get_field_topology(*bucket_ptr, LSField.isovar);
        if (fieldTopology != fieldTopology.base())
          return true;
      }
    }
  }
  return false;
}

LevelSetInterfaceGeometry::LevelSetInterfaceGeometry(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields)
: myActivePart(activePart),
  myCdfemSupport(cdfemSupport),
  myPhaseSupport(phaseSupport),
  myLSFields(LSFields)
{
  set_ls_fields(LSFields);
}

void LevelSetInterfaceGeometry::set_ls_fields(const std::vector<LS_Field> & lsFields)
{
  myLSFields = lsFields;
  for (auto && lsField : myLSFields)
    mySurfaceIdentifiers.push_back(lsField.identifier);
  set_parent_element_selector();
}

static stk::mesh::Selector get_active_and_any_levelset_selector(const stk::mesh::Part & activePart, const std::vector<LS_Field> & LSFields)
{
  stk::mesh::Selector levelsetSelector;
  if (!LSFields.empty())
  {
    levelsetSelector = stk::mesh::selectField(LSFields[0].isovar.field());
    for (size_t i=1; i<LSFields.size(); ++i)
      levelsetSelector |= stk::mesh::selectField(LSFields[i].isovar.field());
  }
  return levelsetSelector & activePart;
}

void LevelSetInterfaceGeometry::set_parent_element_selector()
{
  myParentElementSelector = (is_cdfem_use_case(myPhaseSupport)) ?
    get_cdfem_parent_element_selector(myActivePart, myCdfemSupport, myPhaseSupport) :
    get_active_and_any_levelset_selector(myActivePart, myLSFields);
}

void LevelSetInterfaceGeometry::build_parent_edges_for_mesh(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  const bool addHigherOrderMidSideNodes = determine_if_mesh_has_higher_order_midside_nodes_with_level_set_locally(mesh, myPhaseSupport, myLSFields);
  myParentsToChildMapper.build_map(mesh, myActivePart, myCdfemSupport, addHigherOrderMidSideNodes);
  const bool shouldLinearizeEdges = false;
  myParentEdges = build_parent_edges(mesh, myParentsToChildMapper, shouldLinearizeEdges, myParentElementSelector, myPhaseSupport, myLSFields);
}

void LevelSetInterfaceGeometry::build_parent_edges_for_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  const bool addHigherOrderMidSideNodes = determine_if_mesh_has_higher_order_midside_nodes_with_level_set_locally(mesh, myPhaseSupport, myLSFields);
  myParentsToChildMapper.build_map(mesh, myActivePart, myCdfemSupport, addHigherOrderMidSideNodes);
  const bool shouldLinearizeEdges = false;
  myParentEdges = build_parent_edges_using_elements(mesh, myParentsToChildMapper, shouldLinearizeEdges, elementsToIntersect, myParentElementSelector, myPhaseSupport, myLSFields);
}

void LevelSetInterfaceGeometry::prepare_to_decompose_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  build_parent_edges_for_mesh(mesh, nodesToCapturedDomains);
}

void LevelSetInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh) const
{
  const NodeToCapturedDomainsMap emptyNodesToCapturedDomains;
  prepare_to_intersect_elements(mesh, emptyNodesToCapturedDomains);
}

void LevelSetInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  build_parent_edges_for_mesh(mesh, nodesToCapturedDomains);
}

void LevelSetInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  build_parent_edges_for_elements(mesh, elementsToIntersect, nodesToCapturedDomains);
}

static double compute_edge_length(const stk::mesh::BulkData & mesh, const std::array<stk::mesh::Entity,2> & edgeNodes)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  const stk::math::Vector3d x0(field_data<double>(coordsField, edgeNodes[0]), dim);
  const stk::math::Vector3d x1(field_data<double>(coordsField, edgeNodes[1]), dim);
  return (x1-x0).length();
}

static bool edge_is_possibly_cut(const stk::mesh::BulkData & mesh, const std::vector<LS_Field> & LSFields, const std::array<stk::mesh::Entity,2> & edgeNodes)
{
  const double edgeLength = compute_edge_length(mesh, edgeNodes);

  for(auto && LSField : LSFields)
  {
    FieldRef distVar = LSField.isovar;
    const double * dist0Ptr = field_data<double>(distVar, edgeNodes[0]);
    const double * dist1Ptr = field_data<double>(distVar, edgeNodes[1]);
    if (dist0Ptr && dist1Ptr && std::abs(*dist0Ptr) < edgeLength && std::abs(*dist1Ptr) < edgeLength)
      return true;
  }
  return false;
}

static bool element_has_possibly_cut_edge(const stk::mesh::BulkData & mesh, const std::vector<LS_Field> & LSFields, stk::mesh::Entity elem)
{
  stk::topology elemTopology = mesh.bucket(elem).topology();
  const StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};

  const unsigned numEdges = elemTopology.num_edges();
  for(unsigned i=0; i < numEdges; ++i)
  {
    const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(elemTopology, i);
    if (edge_is_possibly_cut(mesh, LSFields, {{elemNodes[edgeNodeOrdinals[0]], elemNodes[edgeNodeOrdinals[1]]}}))
      return true;
  }
  return false;
}

std::vector<stk::mesh::Entity> LevelSetInterfaceGeometry::get_active_elements_that_may_be_cut_by_levelsets(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const std::vector<LS_Field> & LSFields)
{
  std::vector<stk::mesh::Entity> possibleCutElements;

  const stk::mesh::Selector activeLocallyOwned = activePart & mesh.mesh_meta_data().locally_owned_part();

  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeLocallyOwned))
    for(const auto & elem : *bucketPtr)
      if (element_has_possibly_cut_edge(mesh, LSFields, elem))
        possibleCutElements.push_back(elem);

  return possibleCutElements;
}

std::vector<stk::mesh::Entity> LevelSetInterfaceGeometry::get_possibly_cut_elements(const stk::mesh::BulkData & mesh) const
{
  return get_active_elements_that_may_be_cut_by_levelsets(mesh, myActivePart, myLSFields);
}

static void fill_node_levelset(const stk::mesh::BulkData & mesh, const LS_Field & LSField, const stk::mesh::Entity elem, std::vector<double> & nodeLS)
{
  nodeLS.clear();
  for (auto node : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
  {
    const double * lsPtr = field_data<double>(LSField.isovar, node);
    if (lsPtr)
      nodeLS.push_back(*lsPtr);
  }
}

static bool element_intersects_distance_interval(const stk::mesh::BulkData & mesh, const LS_Field & LSField, const stk::mesh::Entity elem, const std::array<double,2> & loAndHi, std::vector<double> & elemNodeDistWorkspace)
{
  // Note that levelset is assumed to signed distance here
  fill_node_levelset(mesh, LSField, elem, elemNodeDistWorkspace);
  return (elemNodeDistWorkspace.size() == mesh.num_nodes(elem) && InterfaceGeometry::element_with_nodal_distance_intersects_distance_interval(elemNodeDistWorkspace, loAndHi));
}

void LevelSetInterfaceGeometry::fill_active_elements_that_intersect_levelset_interval(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const LS_Field & lsField, const std::array<double,2> loAndHi, std::vector<stk::mesh::Entity> & elementsThaIntersectInterval)
{
  std::vector<double> elementNodeDist;

  const stk::mesh::Selector activeLocallyOwned = activePart & mesh.mesh_meta_data().locally_owned_part();

  elementsThaIntersectInterval.clear();
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeLocallyOwned))
    for(const auto & elem : *bucketPtr)
      if (element_intersects_distance_interval(mesh, lsField, elem, loAndHi, elementNodeDist))
        elementsThaIntersectInterval.push_back(elem);
}

const LS_Field & LevelSetInterfaceGeometry::get_ls_field_with_identifier(const Surface_Identifier surfaceIdentifier) const
{
  auto iter = std::find(mySurfaceIdentifiers.begin(), mySurfaceIdentifiers.end(), surfaceIdentifier);
  STK_ThrowRequire(iter != mySurfaceIdentifiers.end());
  return myLSFields[std::distance(mySurfaceIdentifiers.begin(), iter)];
}

void LevelSetInterfaceGeometry::fill_elements_that_intersect_distance_interval(const stk::mesh::BulkData & mesh, const Surface_Identifier surfaceIdentifier, const std::array<double,2> loAndHi, std::vector<stk::mesh::Entity> & elementsThaIntersectInterval) const
{
  fill_active_elements_that_intersect_levelset_interval(mesh, myActivePart, get_ls_field_with_identifier(surfaceIdentifier), loAndHi, elementsThaIntersectInterval);
}

bool LevelSetInterfaceGeometry::have_enough_levelsets_to_have_interior_intersections_or_multiple_crossings() const
{
  const unsigned minNumLSForInteriorIntersectionsOrMultipleElementCrossings = myPhaseSupport.has_one_levelset_per_phase() ? 3 : 2;
  return myLSFields.size() >= minNumLSForInteriorIntersectionsOrMultipleElementCrossings;
}

bool LevelSetInterfaceGeometry::snapped_elements_may_have_new_intersections() const
{
  return have_enough_levelsets_to_have_interior_intersections_or_multiple_crossings();
}

std::vector<IntersectionPoint> LevelSetInterfaceGeometry::get_edge_intersection_points(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  std::vector<IntersectionPoint> intersectionPoints;
  prepare_to_intersect_elements(mesh, nodesToCapturedDomains);
  const IntersectionPointFilter intersectionPointFilter = keep_all_intersection_points_filter();
  append_intersection_points_from_all_parent_edges(intersectionPoints, myParentEdges, intersectionPointFilter);
  return intersectionPoints;
}

void LevelSetInterfaceGeometry::append_element_intersection_points(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints) const
{
  prepare_to_intersect_elements(mesh, elementsToIntersect, nodesToCapturedDomains);
  append_intersection_points_from_owned_parent_edges(mesh, myParentElementSelector, intersectionPoints, myParentEdges, intersectionPointFilter);
  if (have_enough_levelsets_to_have_interior_intersections_or_multiple_crossings())
  {
    append_intersection_points_from_within_elements_and_owned_faces(mesh, myParentElementSelector, elementsToIntersect, *this, intersectionPointFilter, intersectionPoints);
  }
}

static int get_domain_for_uncut_element(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity element,
    const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges,
    const bool oneLSPerPhase)
{
  if (oneLSPerPhase)
  {
    const auto edgePhases = get_phases_present_on_edges_and_interior(elementParentEdges);
    if (edgePhases.size() == 1)
      return *edgePhases.begin();
  }
  return -1;
}

void LevelSetInterfaceGeometry::store_phase_for_uncut_elements(const stk::mesh::BulkData & mesh) const
{
  std::vector<stk::mesh::Entity> elementsToIntersect;
  stk::mesh::get_entities( mesh, stk::topology::ELEMENT_RANK, elementsToIntersect, false);

  const bool addHigherOrderMidSideNodes = determine_if_mesh_has_higher_order_midside_nodes_with_level_set_locally(mesh, myPhaseSupport, myLSFields);
  myParentsToChildMapper.build_map(mesh, myActivePart, myCdfemSupport, addHigherOrderMidSideNodes);
  const bool shouldLinearizeEdges = true;
  ParentEdgeMap parentEdges = build_parent_edges_using_elements(mesh, myParentsToChildMapper, shouldLinearizeEdges, elementsToIntersect, myParentElementSelector, myPhaseSupport, myLSFields);

  const bool oneLSPerPhase = myLSFields.size() > 1 && myPhaseSupport.has_one_levelset_per_phase();
  std::vector<const CDFEM_Parent_Edge *> elementParentEdges;
  std::vector<bool> areParentEdgesOrientedSameAsElementEdges;

  const std::vector<stk::mesh::Entity> parentElements = get_owned_parent_elements(mesh, myParentElementSelector);
  for (auto element : parentElements)
  {
    fill_element_parent_edges(mesh, element, parentEdges, elementParentEdges, areParentEdgesOrientedSameAsElementEdges);
    const int elementDomain = get_domain_for_uncut_element(mesh,
        element,
        elementParentEdges,
        oneLSPerPhase);
    if (elementDomain >= 0)
      myUncutElementPhases[element] = elementDomain;
  }
}

static NodeToCapturedDomainsMap store_and_communicate_new_snap_node_domains(const stk::mesh::BulkData & mesh, const std::vector<IntersectionPoint> & intersectionPoints, const std::vector<SnapInfo> & snapInfos)
{
  std::vector<stk::mesh::Entity> snapNodes;
  NodeToCapturedDomainsMap newSnapnodesToCapturedDomains;
  for (auto && snapInfo : snapInfos)
  {
    if (snapInfo.get_owner() == mesh.parallel_rank())
    {
      const size_t intersectionPointIndex = snapInfo.get_intersection_point_index();
      stk::mesh::Entity snapNode = mesh.get_entity(stk::topology::NODE_RANK, snapInfo.get_node_global_id());
      snapNodes.push_back(snapNode);
      const IntersectionPoint & intersectionPoint = intersectionPoints[intersectionPointIndex];

      newSnapnodesToCapturedDomains[snapNode] = intersectionPoint.get_sorted_domains();
    }
  }

  communicate_node_captured_domains_for_given_nodes(mesh, snapNodes, newSnapnodesToCapturedDomains);

  return newSnapnodesToCapturedDomains;
}

static std::vector<stk::math::Vector3d> get_node_parametric_coords_after_snapping(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    stk::mesh::Entity snapNode,
    const stk::math::Vector3d snapNodeLocation)
{
  stk::topology elementTopology = mesh.bucket(element).topology();
  const MasterElement & masterElement = MasterElementDeterminer::getMasterElement(elementTopology);
  const double * elemNodeParamCoords = masterElement.nodal_parametric_coordinates();
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  const StkMeshEntities elementNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());

  std::vector<stk::math::Vector3d> nodeLocations;
  fill_element_node_coordinates(mesh, element, coordsField, nodeLocations);

  std::vector<stk::math::Vector3d> nodesCoords;
  for (unsigned n=0; n<elementNodes.size(); ++n)
  {
    if (elementNodes[n] == snapNode)
    {
      const stk::math::Vector3d snapNodeParamCoords = get_parametric_coordinates_of_point(nodeLocations, snapNodeLocation);
      nodesCoords.push_back(snapNodeParamCoords);
    }
    else
    {
      const stk::math::Vector3d nodeCoords(&elemNodeParamCoords[n*dim],dim);
      nodesCoords.push_back(nodeCoords);
    }
  }

  return nodesCoords;
}

static double truncate_to_maintain_positive_shape_function(const double shapeFcn, const double baseShapeFcn)
{
  if (shapeFcn >= baseShapeFcn)
    return 1.;
  return baseShapeFcn/(baseShapeFcn-shapeFcn);
}

static stk::math::Vector3d find_point_within_deformed_and_undeformed_tet(const std::vector<stk::math::Vector3d> & deformedElementParamCoords,
    const int lnn)
{
  stk::topology topology = stk::topology::TETRAHEDRON_4;
  std::array<int,4> permutations{0, 1, 2, 4};
  std::array<int,4> permuteNodes;
  topology.permutation_node_ordinals(permutations[lnn], permuteNodes.data());
  STK_ThrowAssert(permuteNodes[0] == lnn);
  const stk::math::Vector3d & pt = deformedElementParamCoords[lnn];
  const stk::math::Vector3d oppositePt = 1./3.*(deformedElementParamCoords[permuteNodes[1]] + deformedElementParamCoords[permuteNodes[2]] + deformedElementParamCoords[permuteNodes[3]]);
  double fraction = 1.0;
  fraction = std::min(fraction, truncate_to_maintain_positive_shape_function(pt[0], oppositePt[0]));
  fraction = std::min(fraction, truncate_to_maintain_positive_shape_function(pt[1], oppositePt[1]));
  fraction = std::min(fraction, truncate_to_maintain_positive_shape_function(pt[2], oppositePt[2]));
  fraction = std::min(fraction, truncate_to_maintain_positive_shape_function(1.-pt[0]-pt[1]-pt[2], 1.-oppositePt[0]-oppositePt[1]-oppositePt[2]));
  const double centroidWt = 0.25*fraction;
  return centroidWt*pt + (1.-centroidWt)*oppositePt;
}

static stk::math::Vector3d find_point_within_deformed_and_undeformed_tri(const std::vector<stk::math::Vector3d> & deformedElementParamCoords,
    const int lnn)
{
  stk::topology topology = stk::topology::TRIANGLE_3_2D;
  std::array<int,3> permutations{0, 2, 1};
  std::array<int,3> permuteNodes;
  topology.permutation_node_ordinals(permutations[lnn], permuteNodes.data());
  STK_ThrowAssert(permuteNodes[0] == lnn);
  const stk::math::Vector3d & pt = deformedElementParamCoords[lnn];
  const stk::math::Vector3d oppositePt = 0.5*(deformedElementParamCoords[permuteNodes[1]] + deformedElementParamCoords[permuteNodes[2]]);
  double fraction = 1.0;
  fraction = std::min(fraction, truncate_to_maintain_positive_shape_function(pt[0], oppositePt[0]));
  fraction = std::min(fraction, truncate_to_maintain_positive_shape_function(pt[1], oppositePt[1]));
  fraction = std::min(fraction, truncate_to_maintain_positive_shape_function(1.-pt[0]-pt[1], 1.-oppositePt[0]-oppositePt[1]));
  const double centroidWt = 1./3.*fraction;
  return centroidWt*pt + (1.-centroidWt)*oppositePt;
}

static int get_node_of_element(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    stk::mesh::Entity node)
{
  const StkMeshEntities elementNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  for (unsigned n=0; n<elementNodes.size(); ++n)
    if (elementNodes[n] == node)
      return n;
  STK_ThrowRequire(false);
  return -1;
}

static stk::math::Vector3d find_point_within_deformed_and_undeformed_element(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    stk::mesh::Entity snapNode,
    const std::vector<stk::math::Vector3d> & deformedElementParamCoords)
{
  const int lnn = get_node_of_element(mesh, element, snapNode);

  stk::topology elementTopology = mesh.bucket(element).topology();
  if (elementTopology.base() == stk::topology::TETRAHEDRON_4)
    return find_point_within_deformed_and_undeformed_tet(deformedElementParamCoords, lnn);
  STK_ThrowRequire(elementTopology.base() == stk::topology::TRIANGLE_3_2D);
  return find_point_within_deformed_and_undeformed_tri(deformedElementParamCoords, lnn);
}

static std::vector<const std::vector<int>*> get_node_snap_domains(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    stk::mesh::Entity snapNode,
    const std::vector<int> & snapNodeDomains,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  const StkMeshEntities elementNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  std::vector<const std::vector<int>*> nodeDomains;
  nodeDomains.reserve(elementNodes.size());
  for (auto node : StkMeshEntities{mesh.begin_nodes(element), mesh.end_nodes(element)})
  {
    if (node == snapNode)
    {
      nodeDomains.push_back(&snapNodeDomains);
    }
    else
    {
      const auto iter = nodesToCapturedDomains.find(node);
      if (iter == nodesToCapturedDomains.end())
        nodeDomains.push_back(nullptr);
      else
        nodeDomains.push_back(&(iter->second));
    }
  }

  return nodeDomains;
}

static bool domain_is_common_to_all_entities(const int domain, const std::vector<const std::vector<int>*> & entitiesDomains)
{
  for (auto && domains : entitiesDomains)
    if (!std::binary_search(domains->begin(), domains->end(), domain))
      return false;
  return true;
}

static std::set<int> get_common_domains(std::vector<const std::vector<int>*> entitiesDomains)
{
  std::set<int> commonDomains;

  if (entitiesDomains.empty())
    return commonDomains;

  for (auto & entityDomains : entitiesDomains)
    if (!entityDomains)
      return commonDomains;

  for (int domain : *entitiesDomains[0])
    if (domain_is_common_to_all_entities(domain, entitiesDomains))
      commonDomains.insert(domain);

  return commonDomains;
}

static bool will_have_uncaptured_edge_intersection_after_snapping(const Element_Cutter & cutter, stk::topology elementTopology, const std::vector<stk::math::Vector3d> nodeCoords,
    const std::vector<const std::vector<int>*> & nodeDomains)
{
  std::vector<InterfaceID> interfacesWithCuttingSurface;
  cutter.fill_interfaces_with_cutting_surface(interfacesWithCuttingSurface);

  const unsigned numEdges = elementTopology.num_edges();
  for(unsigned i=0; i < numEdges; ++i)
  {
    const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(elementTopology(), i);
    const std::vector<int>* node0Domains = nodeDomains[edgeNodeOrdinals[0]];
    const std::vector<int>* node1Domains = nodeDomains[edgeNodeOrdinals[1]];

    for (auto && interface : interfacesWithCuttingSurface)
    {
      if (!captures_interface(node0Domains, interface) && !captures_interface(node1Domains, interface))
      {
        const std::array<stk::math::Vector3d,2> edgeNodeCoords{nodeCoords[edgeNodeOrdinals[0]], nodeCoords[edgeNodeOrdinals[1]]};
        if (cutter.have_crossing(interface, edgeNodeCoords))
          return true;
      }
    }
  }
  return false;
}

static int get_domain_of_element_if_it_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
    const ParentsToChildMapper parentsToChildMapper,
    const stk::mesh::Selector & parentElementSelector,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields,
    stk::mesh::Entity element,
    stk::mesh::Entity snapNode,
    const stk::math::Vector3d & snapNodeLocation,
    const std::vector<const std::vector<int>*> & elementNodeDomains)
{
  const auto diagonalPicker = temporary_build_always_true_diagonal_picker();
  const bool oneLSPerPhase = true;

  const bool shouldLinearizeEdges = true;
  ParentEdgeMap parentEdges = build_parent_edges_using_elements(mesh, parentsToChildMapper, shouldLinearizeEdges, {element}, parentElementSelector, phaseSupport, LSFields);
  if (parentEdges.empty())
    return -1; // not parent element

  std::vector<const CDFEM_Parent_Edge *> elementParentEdges;
  std::vector<bool> areParentEdgesOrientedSameAsElementEdges;
  fill_element_parent_edges(mesh, element, parentEdges, elementParentEdges, areParentEdgesOrientedSameAsElementEdges);

  const auto edgePhases = get_phases_present_on_edges_and_interior(elementParentEdges);
  if (edgePhases.size() == 1)
    return *edgePhases.begin();

  stk::topology elementTopology = mesh.bucket(element).topology();
  const MasterElement & masterElement = MasterElementDeterminer::getMasterElement(elementTopology);
  std::unique_ptr<Element_Cutter> elementCutter = create_element_cutter(oneLSPerPhase, masterElement, elementParentEdges, areParentEdgesOrientedSameAsElementEdges, diagonalPicker);

  const std::vector<stk::math::Vector3d> nodeParamCoordsAfterSnapping = get_node_parametric_coords_after_snapping(mesh, element, snapNode, snapNodeLocation);
  if (!will_have_uncaptured_edge_intersection_after_snapping(*elementCutter, elementTopology, nodeParamCoordsAfterSnapping, elementNodeDomains))
  {
    const stk::math::Vector3d evaluationPt = find_point_within_deformed_and_undeformed_element(mesh, element, snapNode, nodeParamCoordsAfterSnapping);
    return elementCutter->get_ls_per_interface_phase_at_location(evaluationPt);
  }

  return -1;
}

static int get_sign_of_uncrossed_edges(const InterfaceID interface, const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges)
{
  int uncrossedSign = 0;
  for (auto && edge : elementParentEdges)
  {
    if (edge->have_crossing(interface))
    {
      return 0;
    }
    else
    {
      const int edgeSign = edge->get_crossing_sign(interface);
      if (edgeSign == -uncrossedSign)
        return 0;
      uncrossedSign = edgeSign;
    }
  }
  return uncrossedSign;
}

static int get_sign_of_element_if_it_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
    const ParentsToChildMapper parentsToChildMapper,
    const stk::mesh::Selector & parentElementSelector,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields,
    stk::mesh::Entity element,
    stk::mesh::Entity snapNode,
    const stk::math::Vector3d & snapNodeLocation,
    const InterfaceID & interface)
{
  const auto diagonalPicker = temporary_build_always_true_diagonal_picker();
  const bool oneLSPerPhase = false;

  const bool shouldLinearizeEdges = true;
  ParentEdgeMap parentEdges = build_parent_edges_using_elements(mesh, parentsToChildMapper, shouldLinearizeEdges, {element}, parentElementSelector, phaseSupport, LSFields);
  if (parentEdges.empty())
    return -1; // not parent element

  std::vector<const CDFEM_Parent_Edge *> elementParentEdges;
  std::vector<bool> areParentEdgesOrientedSameAsElementEdges;
  fill_element_parent_edges(mesh, element, parentEdges, elementParentEdges, areParentEdgesOrientedSameAsElementEdges);

  const int uncrossedSign = get_sign_of_uncrossed_edges(interface, elementParentEdges);
  if (uncrossedSign != 0)
    return uncrossedSign;

  stk::topology elementTopology = mesh.bucket(element).topology();
  const MasterElement & masterElement = MasterElementDeterminer::getMasterElement(elementTopology);
  std::unique_ptr<Element_Cutter> elementCutter = create_element_cutter(oneLSPerPhase, masterElement, elementParentEdges, areParentEdgesOrientedSameAsElementEdges, diagonalPicker);

  const std::vector<stk::math::Vector3d> nodeParamCoordsAfterSnapping = get_node_parametric_coords_after_snapping(mesh, element, snapNode, snapNodeLocation);
  const stk::math::Vector3d evaluationPt = find_point_within_deformed_and_undeformed_element(mesh, element, snapNode, nodeParamCoordsAfterSnapping);
  return elementCutter->sign_at_position(interface, evaluationPt);
}

static void set_domains_for_element_if_it_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
    const ParentsToChildMapper parentsToChildMapper,
    const stk::mesh::Selector & parentElementSelector,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields,
    stk::mesh::Entity element,
    stk::mesh::Entity snapNode,
    const std::vector<int> & snapNodeDomains,
    const stk::math::Vector3d & snapNodeLocation,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    ElementToDomainMap & elementsToDomain )
{
  auto iter = elementsToDomain.lower_bound(element);
  if (iter == elementsToDomain.end() || iter->first != element)
  {
    const std::vector<const std::vector<int>*> nodeDomains = get_node_snap_domains(mesh, element, snapNode, snapNodeDomains, nodesToCapturedDomains);

    const std::set<int> commonDomains = get_common_domains(nodeDomains);

    int elementDomain = -1;
    if (commonDomains.size() == 1)
    {
      elementDomain = *commonDomains.begin();
    }
    else if (commonDomains.size() > 1)
    {
      elementDomain = get_domain_of_element_if_it_will_be_uncut_after_snapping(mesh, parentsToChildMapper, parentElementSelector, phaseSupport, LSFields, element, snapNode, snapNodeLocation, nodeDomains);
    }

    if (elementDomain >= 0)
        elementsToDomain.emplace_hint(iter, element, elementDomain);
  }
}

static void set_sign_for_element_if_it_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
    const ParentsToChildMapper parentsToChildMapper,
    const stk::mesh::Selector & parentElementSelector,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields,
    stk::mesh::Entity element,
    stk::mesh::Entity snapNode,
    const std::vector<int> & snapNodeDomains,
    const stk::math::Vector3d & snapNodeLocation,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    ElementToDomainMap & elementsToDomain )
{
  auto iter = elementsToDomain.lower_bound(element);
  if (iter == elementsToDomain.end() || iter->first != element)
  {
    const std::vector<const std::vector<int>*> nodeDomains = get_node_snap_domains(mesh, element, snapNode, snapNodeDomains, nodesToCapturedDomains);

    const std::set<int> commonDomains = get_common_domains(nodeDomains);

    int elementSign = 0;
    if (commonDomains.size() == 1)
    {
      const int lsIndex = *commonDomains.begin();
      const InterfaceID interface(lsIndex,lsIndex);
      elementSign = get_sign_of_element_if_it_will_be_uncut_after_snapping(mesh, parentsToChildMapper, parentElementSelector, phaseSupport, LSFields, element, snapNode, snapNodeLocation, interface);
    }

    if (elementSign != 0)
        elementsToDomain.emplace_hint(iter, element, elementSign);
  }
}

void LevelSetInterfaceGeometry::store_phase_for_elements_that_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
      const std::vector<IntersectionPoint> & intersectionPoints,
      const std::vector<SnapInfo> & snapInfos,
      const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  const bool oneLSPerPhase = myLSFields.size() > 1 && myPhaseSupport.has_one_levelset_per_phase();
  if (!oneLSPerPhase && myLSFields.size() > 1)
    return; //FIXME: Fix for more than one ls per interface

  const NodeToCapturedDomainsMap newSnapnodesToCapturedDomains = store_and_communicate_new_snap_node_domains(mesh, intersectionPoints, snapInfos);

  for (auto && snapInfo : snapInfos)
  {
    stk::mesh::Entity snapNode = mesh.get_entity(stk::topology::NODE_RANK, snapInfo.get_node_global_id());
    const stk::math::Vector3d snapLocation = snapInfo.get_snap_location();
    const auto & newSnapNodeDomains = newSnapnodesToCapturedDomains.at(snapNode);
    for (auto elem : StkMeshEntities{mesh.begin_elements(snapNode), mesh.end_elements(snapNode)})
    {
      if (mesh.bucket(elem).owned() && mesh.bucket(elem).member(myActivePart))
      {
        if (oneLSPerPhase)
          set_domains_for_element_if_it_will_be_uncut_after_snapping(mesh, myParentsToChildMapper, myParentElementSelector, myPhaseSupport, myLSFields, elem, snapNode, newSnapNodeDomains, snapLocation, nodesToCapturedDomains, myUncutElementPhases);
        else
          set_sign_for_element_if_it_will_be_uncut_after_snapping(mesh, myParentsToChildMapper, myParentElementSelector, myPhaseSupport, myLSFields, elem, snapNode, newSnapNodeDomains, snapLocation, nodesToCapturedDomains, myUncutElementPhases);
      }
    }
  }
}

}
