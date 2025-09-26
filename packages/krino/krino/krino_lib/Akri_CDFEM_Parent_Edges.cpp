// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh_Utils.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MasterElement.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <Akri_ParentsToChildMapper.hpp>
#include <Akri_Phase_Support.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/HashEntityAndEntityKey.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <unordered_map>

namespace krino
{

std::set<int> get_phases_present_on_edges(const std::vector<const CDFEM_Parent_Edge *> & parentEdges)
{
  std::set<int> phasesPresent;
  for(auto && parentEdge : parentEdges)
    if(parentEdge)
      for (int phase : parentEdge->get_edge_phases())
        phasesPresent.insert(phase);
  return phasesPresent;
}

static void add_interface_phases(const InterfaceID & interface, std::set<int> & phases)
{
  phases.insert(interface.first_ls());
  phases.insert(interface.second_ls());
}

static std::set<InterfaceID> get_all_interfaces_including_fake(const std::vector<const CDFEM_Parent_Edge *> & parentEdges)
{
  std::set<InterfaceID> allInterfacesIncludingFake;
  for(auto && parentEdge : parentEdges)
    if(parentEdge)
      for (auto && crossing : parentEdge->get_crossings_including_fake())
        allInterfacesIncludingFake.insert(crossing.first);
  return allInterfacesIncludingFake;
}

bool phase_has_interfaces_to_all_other_phases(const int phase, const std::set<int> & phasesPresent, const std::set<InterfaceID> & interfaces)
{
  for (int otherPhase : phasesPresent)
    if (interfaces.find(InterfaceID(phase,otherPhase)) == interfaces.end())
      return false;
  return true;
}

void add_phases_possibly_present_on_interior(const std::vector<const CDFEM_Parent_Edge *> & parentEdges, std::set<int> & phasesPresent)
{
  const std::set<InterfaceID> allInterfacesIncludingFake = get_all_interfaces_including_fake(parentEdges);

  std::set<int> phasesPossiblyPresent;
  for (auto && interface : allInterfacesIncludingFake)
    add_interface_phases(interface, phasesPossiblyPresent);

  // This should be unit tested.  What about order dependency?
  for (int phase : phasesPossiblyPresent)
    if (phasesPresent.find(phase) == phasesPresent.end())
      if (phase_has_interfaces_to_all_other_phases(phase, phasesPresent, allInterfacesIncludingFake))
        phasesPresent.insert(phase);
}

std::set<int> get_phases_present_on_edges_and_interior(const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges)
{
  std::set<int> phasesPresent = get_phases_present_on_edges(elementParentEdges);
  add_phases_possibly_present_on_interior(elementParentEdges, phasesPresent);
  return phasesPresent;
}

static
std::set<InterfaceID> get_interfaces_with_any_real_edge_crossings(const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges)
{
  std::set<InterfaceID> interfacesPresent;
  for(auto && parentEdge : elementParentEdges)
    if(parentEdge)
      for(const auto & crossing : parentEdge->get_crossings())
        interfacesPresent.insert(crossing.first);
  return interfacesPresent;
}

std::set<InterfaceID> get_interfaces_present(const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges)
{
  return get_interfaces_with_any_real_edge_crossings(elementParentEdges);
}

void fill_element_parent_edges(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity elem,
    const ParentEdgeMap & parentEdges,
    std::vector<const CDFEM_Parent_Edge *> & elementParentEdges,
    std::vector<bool> & areParentEdgesAreOrientedSameAsElementEdges)
{
  const stk::topology stk_topology = mesh.bucket(elem).topology();
  const unsigned numEdges = stk_topology.num_edges();
  const stk::mesh::Entity * const nodes = mesh.begin_nodes(elem);
  elementParentEdges.assign(numEdges, nullptr);
  areParentEdgesAreOrientedSameAsElementEdges.assign(numEdges, true);
  for(unsigned i=0; i < numEdges; ++i)
  {
    const unsigned * edge_node_ordinals = get_edge_node_ordinals(stk_topology, i);
    const CDFEM_Parent_Edge * parent_edge =
        find_parent_edge(mesh, parentEdges, nodes[edge_node_ordinals[0]], nodes[edge_node_ordinals[1]]);
    elementParentEdges[i] = parent_edge;
    if(parent_edge && parent_edge->get_parent_nodes().first != nodes[edge_node_ordinals[0]])
      areParentEdgesAreOrientedSameAsElementEdges[i] = false;
  }
}

void fill_face_nodes_and_parent_edges(const stk::topology & elementTopology,
    const int iFace,
    const std::vector<stk::mesh::Entity> & elementNodes,
    const std::vector<const CDFEM_Parent_Edge *> & elementParentEdges,
    const std::vector<bool> & areParentEdgesOrientedSameAsElementEdges,
    std::vector<stk::mesh::Entity> & faceNodes,
    std::vector<const CDFEM_Parent_Edge *> & faceParentEdges,
    std::vector<bool> & areParentEdgesOrientedSameAsFaceEdges)
{
  STK_ThrowAssert(elementTopology == stk::topology::TETRAHEDRON_4);
  constexpr std::array<std::array<int,3>,4> faceEdges = {{ {{0,4,3}}, {{1,5,4}}, {{3,5,2}}, {{2,1,0}} }};
  constexpr std::array<std::array<bool,3>,4> isFaceEdgeOrientedSameAsElementEdge = {{ {{true,true,false}}, {{true,true,false}}, {{true,false,true}}, {{false,false,false}} }};

  faceNodes.resize(3);
  elementTopology.face_nodes(elementNodes, iFace, faceNodes.data());

  constexpr int numFaceEdges = 3;
  faceParentEdges.resize(numFaceEdges);
  areParentEdgesOrientedSameAsFaceEdges.resize(3);
  for (int i=0; i<numFaceEdges; ++i)
  {
    faceParentEdges[i] = elementParentEdges[faceEdges[iFace][i]];
    areParentEdgesOrientedSameAsFaceEdges[i] = isFaceEdgeOrientedSameAsElementEdge[iFace][i] ?
        areParentEdgesOrientedSameAsElementEdges[faceEdges[iFace][i]] :
        !areParentEdgesOrientedSameAsElementEdges[faceEdges[iFace][i]];
  }
}

static bool in_block_decomposed_by_ls(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const stk::mesh::Entity node_or_elem,
    const LS_Field & lsField )
{
  STK_ThrowAssert(phaseSupport.is_cdfem_use_case());

  for(auto * partPtr : mesh.bucket(node_or_elem).supersets())
  {
    if (partPtr->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
      stk::io::is_part_io_part(*partPtr) &&
      !phaseSupport.is_nonconformal(*partPtr) &&
      phaseSupport.level_set_is_used_by_nonconformal_part(lsField.identifier, &phaseSupport.find_nonconformal_part(*partPtr)))
      return true;
  }
  return false;
}

static bool has_levelset_registered(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity node_or_elem,
    const LS_Field & lsField )
{
  const stk::mesh::Selector fieldSelector(lsField.isovar.field());
  return fieldSelector(mesh.bucket(node_or_elem));
}

static bool
has_io_part_containing_phase(const Phase_Support & phase_support, const stk::mesh::PartVector & parts, const PhaseTag & phase)
{
  for(stk::mesh::PartVector::const_iterator part_iter = parts.begin(); part_iter != parts.end(); ++part_iter)
  {
    const stk::mesh::Part * const part = *part_iter;
    if (part->primary_entity_rank() != stk::topology::ELEMENT_RANK || // limit ourselves to volume parts
        !(stk::io::is_part_io_part(*part) ||
        phase_support.is_nonconformal(*part)))
      continue;

    const PhaseTag & iopart_phase = phase_support.get_iopart_phase(*part);
    if (iopart_phase.contain(phase))
    {
      return true;
    }
  }

  return false;
}

static bool node_touches_alive_block(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    stk::mesh::Entity node,
    const LS_Field & lsField )
{
  const CDFEM_Inequality_Spec * death_spec = lsField.deathPtr;
  if (nullptr == death_spec) return false;

  const PhaseTag & alive_phase = death_spec->get_active_phase();
  const stk::mesh::PartVector & node_parts = mesh.bucket(node).supersets();
  return has_io_part_containing_phase(phaseSupport, node_parts, alive_phase);
}

static bool node_touches_dead_block(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    stk::mesh::Entity node,
    const LS_Field & lsField )
{
  const CDFEM_Inequality_Spec * death_spec = lsField.deathPtr;
  if (nullptr == death_spec) return false;

  const PhaseTag & dead_phase = death_spec->get_deactivated_phase();
  const stk::mesh::PartVector & node_parts = mesh.bucket(node).supersets();
  return has_io_part_containing_phase(phaseSupport, node_parts, dead_phase);
}

static bool node_has_real_ls_value(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    stk::mesh::Entity node,
    const LS_Field & lsField )
{
  if( node_touches_dead_block(mesh, phaseSupport, node, lsField) &&
      !node_touches_alive_block(mesh, phaseSupport, node, lsField) )
    return false;

  return (phaseSupport.is_cdfem_use_case()) ?
      in_block_decomposed_by_ls(mesh, phaseSupport, node, lsField) :
      has_levelset_registered(mesh, node, lsField);
}

static void debug_print_edge_info(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & lsFields,
    const std::vector<stk::mesh::Entity> & edge_nodes,
    const std::vector<std::vector<double> > & nodes_isovar)
{
  const int num_nodes = nodes_isovar.size();
  if(krinolog.shouldPrint(LOG_DEBUG))
  {
    const auto old_precision = krinolog.getStream().precision();
    krinolog.getStream().precision(16);
    const int num_ls = lsFields.size();
    krinolog << stk::diag::dendl << "CDFEM_Parent_Edge::find_crossings():" << "\n";
    for ( int n = 0; n < num_nodes; ++n )
    {
      krinolog << "  Node: " << mesh.identifier(edge_nodes[n]) << ", in_block_decomposed_by_ls = { ";
      for ( int ls_index = 0; ls_index < num_ls; ++ls_index ) krinolog << in_block_decomposed_by_ls(mesh, phaseSupport, edge_nodes[n], lsFields[ls_index]) << " ";
      krinolog << "}, has_real_ls_value = { ";
      for ( int ls_index = 0; ls_index < num_ls; ++ls_index ) krinolog << node_has_real_ls_value(mesh, phaseSupport, edge_nodes[n], lsFields[ls_index]) << " ";
      krinolog << "}, ls = { ";
      for ( int ls_index = 0; ls_index < num_ls; ++ls_index ) krinolog << nodes_isovar[n][ls_index] << " ";
      krinolog << "}";
      krinolog << stk::diag::dendl;
    }
    krinolog.getStream().precision(old_precision);
  }
}

double ls_node_value(const stk::mesh::BulkData & mesh,
  const Phase_Support & phaseSupport,
  const LS_Field & lsField,
  const stk::mesh::Entity node)
{
  const CDFEM_Inequality_Spec * death_spec = lsField.deathPtr;
  FieldRef isovar = lsField.isovar;

  const bool node_has_ls = node_has_real_ls_value(mesh, phaseSupport, node, lsField);
  double nodeLS = 0.0;
  if (node_has_ls)
  {
    if (nullptr != death_spec && isovar.entity_rank() == stk::topology::ELEMENT_RANK)
    {
      // currently requires aura to work correctly in parallel
      STK_ThrowAssertMsg(mesh.is_automatic_aura_on(), "Capability requires aura.");
      bool have_pos_elem = false;
      bool have_neg_elem = false;
      const unsigned num_node_elems = mesh.num_elements(node);
      const stk::mesh::Entity* node_elems = mesh.begin_elements(node);
      for (unsigned node_elem_index=0; node_elem_index<num_node_elems; ++node_elem_index)
      {
	stk::mesh::Entity node_elem = node_elems[node_elem_index];
	const double * isoptr = field_data<double>(isovar, node_elem);
	if (nullptr != isoptr)
	{
	  if (*isoptr - lsField.isoval < 0.)
	    have_neg_elem = true;
	  else
	    have_pos_elem = true;
	}
      }
      nodeLS = (have_pos_elem) ? (have_neg_elem ? 0.0 : 1.0) : -1.0;
    }
    else
    {
      const double * isoptr = field_data<double>(isovar, node);
      STK_ThrowRequireMsg(nullptr != isoptr, "Isovar " << isovar.name() << " missing on node " << debug_entity(mesh, node));
      nodeLS = *isoptr - lsField.isoval;
    }
  }
  else if (nullptr != death_spec && node_touches_dead_block(mesh, phaseSupport, node, lsField))
  {
    const bool dead_is_positive = death_spec->get_deactivated_phase().contain(lsField.identifier,+1);
    if (dead_is_positive) nodeLS = 1.0;
  }

  if (nullptr != death_spec && node_touches_dead_block(mesh, phaseSupport, node, lsField))
  {
    const bool dead_is_positive = death_spec->get_deactivated_phase().contain(lsField.identifier, +1);
    if (dead_is_positive && nodeLS < 0.0)
    {
      if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "Setting node " << mesh.identifier(node) << " to zero to enforce irreversibility, ls = " << nodeLS << "\n";
      nodeLS = 0.0;
    }
    else if (!dead_is_positive && nodeLS >= 0.0)
    {
      if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << "Setting node " << mesh.identifier(node) << " to -REAL_MIN to enforce irreversibility, ls = " << nodeLS << "\n";
      nodeLS = -std::numeric_limits<double>::min();
    }
  }

  return nodeLS;
}

static void edge_ls_node_values(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & lsFields,
    const std::vector<stk::mesh::Entity> & edgeNodes,
    std::vector<std::vector<double> >& nodesIsovar)
{
  const unsigned numNodes = edgeNodes.size();
  const unsigned numLS = lsFields.size();

  nodesIsovar.resize(numNodes);
  for (unsigned n = 0; n < numNodes; ++n)
  {
    nodesIsovar[n].resize(numLS);
    for (unsigned lsIndex = 0; lsIndex < numLS; ++lsIndex)
      nodesIsovar[n][lsIndex] = ls_node_value(mesh, phaseSupport, lsFields[lsIndex], edgeNodes[n]);
  }
}

static CDFEM_Parent_Edge &
build_parent_edge(const stk::mesh::BulkData & mesh, ParentEdgeMap & parentEdges, const ParentsToChildMapper & parentsToChildMapper, const bool linearizeEdge, stk::mesh::Entity node0, stk::mesh::Entity node1)
{
  const stk::mesh::EntityId id0 = mesh.identifier(node0);
  const stk::mesh::EntityId id1 = mesh.identifier(node1);
  const ParentEdgeKey edge_key(id0, id1);
  CDFEM_Parent_Edge & edge = parentEdges[edge_key];

  std::vector<stk::mesh::Entity> edgeNodes;
  std::vector<double> edgeNodePositions;

  if(!edge.valid())
  {
    const std::array<stk::mesh::Entity,2> nodes = (id0 < id1) ? std::array<stk::mesh::Entity,2>{node0, node1} : std::array<stk::mesh::Entity,2>{node1, node0};
    if (linearizeEdge)
      fill_linear_edge_nodes_and_positions(nodes[0], nodes[1], edgeNodes, edgeNodePositions);
    else
      fill_edge_nodes_and_positions(mesh, nodes[0], nodes[1], parentsToChildMapper, edgeNodes, edgeNodePositions);
    edge = CDFEM_Parent_Edge(edgeNodes, edgeNodePositions);
  }
  return edge;
}

static void find_parent_edge_crossings(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & lsFields,
    ParentEdgeMap & parentEdges)
{
  std::vector<std::vector<double>> nodes_isovar;

  for (auto && map_entry : parentEdges)
  {
    CDFEM_Parent_Edge & edge = map_entry.second;
    const std::vector<stk::mesh::Entity> edge_nodes = edge.get_nodes();
    edge_ls_node_values(mesh, phaseSupport, lsFields, edge_nodes, nodes_isovar);
    debug_print_edge_info(mesh, phaseSupport, lsFields, edge_nodes, nodes_isovar);
    edge.find_crossings(phaseSupport.has_one_levelset_per_phase(), nodes_isovar);
  }
  if(krinolog.shouldPrint(LOG_DEBUG)) krinolog << stk::diag::dendl;
}

stk::mesh::Selector get_decomposed_cdfem_parent_element_selector(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
{
  stk::mesh::Selector parentElementSelector =
      (cdfemSupport.get_parent_part() | (activePart & !cdfemSupport.get_child_part())) &
      phaseSupport.get_all_decomposed_blocks_selector();

  return parentElementSelector;
}

stk::mesh::Selector get_potential_cdfem_parent_element_selector(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport)
{
  stk::mesh::Selector parentElementSelector =
      (cdfemSupport.get_parent_part() | (activePart & !cdfemSupport.get_child_part()));;

  return parentElementSelector;
}

std::vector<stk::mesh::Entity> get_owned_parent_elements(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & parentElementSelector)
{
  const stk::mesh::Selector ownedParentElementSelector = parentElementSelector & mesh.mesh_meta_data().locally_owned_part();
  std::vector<stk::mesh::Entity> parentElements;
  stk::mesh::get_selected_entities( ownedParentElementSelector, mesh.buckets(stk::topology::ELEMENT_RANK), parentElements, false );
  return parentElements;
}

ParentEdgeMap
build_parent_edges(const stk::mesh::BulkData & mesh,
    const ParentsToChildMapper & parentsToChildMapper,
    const bool shouldLinearizeEdges,
    const stk::mesh::Selector & parentElementSelector,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields)
{
  std::vector<stk::mesh::Entity> elements;
  stk::mesh::get_entities( mesh, stk::topology::ELEMENT_RANK, mesh.mesh_meta_data().locally_owned_part(), elements, false);

  return build_parent_edges_using_elements(mesh, parentsToChildMapper, shouldLinearizeEdges, elements, parentElementSelector, phaseSupport, LSFields);
}

ParentEdgeMap
build_parent_edges_using_elements(const stk::mesh::BulkData & mesh,
    const ParentsToChildMapper & parentsToChildMapper,
    const bool shouldLinearizeEdges,
    const std::vector<stk::mesh::Entity> & elements,
    const stk::mesh::Selector & parentElementSelector,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields)
{
  ParentEdgeMap parentEdges;

  for (auto && elem : elements)
  {
    if (parentElementSelector(mesh.bucket(elem)))
    {
      const stk::topology topology = mesh.bucket(elem).topology();
      const stk::mesh::Entity* elem_nodes = mesh.begin_nodes(elem);
      const unsigned numEdges = topology.num_edges();

      for (unsigned iedge = 0; iedge < numEdges; ++iedge)
      {
        const unsigned * edge_node_ordinals = get_edge_node_ordinals(topology, iedge);
        stk::mesh::Entity node0 = elem_nodes[edge_node_ordinals[0]];
        stk::mesh::Entity node1 = elem_nodes[edge_node_ordinals[1]];

        build_parent_edge(mesh, parentEdges, parentsToChildMapper, shouldLinearizeEdges, node0, node1);
      }
    }
  }

  find_parent_edge_crossings(mesh, phaseSupport, LSFields, parentEdges);

  return parentEdges;
}

const CDFEM_Parent_Edge *
find_parent_edge(const stk::mesh::BulkData & mesh, const ParentEdgeMap & parentEdges, stk::mesh::Entity node0, stk::mesh::Entity node1)
{
  const stk::mesh::EntityId id0 = mesh.identifier(node0);
  const stk::mesh::EntityId id1 = mesh.identifier(node1);
  const ParentEdgeKey edge_key(id0, id1);

  ParentEdgeMap::const_iterator it = parentEdges.find(edge_key);

  if (it != parentEdges.end())
    return &(it->second);

  return nullptr;
}

} // namespace
