// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <functional>
#include <sstream>
#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <optional>
#include <math.h>
#include <set>

#include <Akri_AuxMetaData.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CubitJournal.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Element.hpp>
#include <Akri_ElementCutterUtils.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MultipleEdgeCrossings.hpp>
#include <Akri_ParametricElementLibrary.hpp>
#include <Akri_PhaseTagRegistry.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_ProlongationData.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_SharpEdgeRecoveryBySwaps.hpp>
#include <Akri_SharpFeatureIntersectionPoints.hpp>
#include <Akri_SharpFeatureLeafUtils.hpp>
#include <Akri_SubElement.hpp>
#include <Akri_SubElementNodeAncestry.hpp>
#include <Akri_SubElementNodeOrder.hpp>
#include <stk_math/StkVector.hpp>

namespace krino{

ElementObj::ElementObj(ParametricElementLibrary & parametricElemLibrary, const stk::mesh::BulkData & stkMesh, stk::mesh::Entity elemEntity)
    : myParametricElemLib(parametricElemLibrary),
      my_parametric_elem(parametricElemLibrary.get_parametric_element(stkMesh.bucket(elemEntity).topology())),
      my_phase(&PhaseTagRegistry::empty_phase()),
      my_entity(elemEntity),
      my_entityId(stkMesh.identifier(elemEntity)),
      my_prolongation_data(NULL)
{
}

ElementObj::ElementObj(ParametricElementLibrary & parametricElemLibrary, const stk::topology elem_topology, const NodeVec & nodes)
   :  myParametricElemLib(parametricElemLibrary),
      my_parametric_elem(parametricElemLibrary.get_parametric_element(elem_topology)),
      my_nodes(nodes),
      my_phase(&PhaseTagRegistry::empty_phase()),
      my_entity(),
      my_entityId(0),
      my_prolongation_data(NULL)
{
}

ElementObj::~ElementObj() {}

void
ElementObj::integration_locations(
    std::vector<stk::math::Vector3d> & intg_pt_locations,
    const ParametricElement & me)
{
  const unsigned num_intg_pts = me.num_intg_pts();
  const double * intg_pt_loc_ptr = me.intg_pt_locations();
  const unsigned dim = me.topology_dimension();

  intg_pt_locations.resize(num_intg_pts);
  for(unsigned i=0; i<num_intg_pts; ++i) intg_pt_locations[i] = stk::math::Vector3d(intg_pt_loc_ptr+i*dim, dim);
}

void
ElementObj::integration_weights(
    std::vector<double> & intg_weights, // includes both gauss point weight and detJ
    const int numCoordDims,
    const std::vector<double> & mesh_coords,
    const ParametricElement & me,
    const ParametricElement & mesh_me)
{
  const unsigned num_intg_pts = me.num_intg_pts();
  const unsigned dim = me.topology_dimension();

  std::vector<double> det_J(num_intg_pts);
  double det_J_error;

  if ( me.get_topology() == mesh_me.get_topology() )
  {
    mesh_me.determinant( numCoordDims, 1, mesh_coords.data(), det_J.data(), &det_J_error );
  }
  else
  {
    const int num_coord_dofs = mesh_me.get_topology().num_nodes();
    const double * intg_pt_loc_ptr = me.intg_pt_locations();
    std::vector<double> d_shapef_coords(num_coord_dofs*dim*num_intg_pts);
    mesh_me.shape_fcn_deriv(num_intg_pts, intg_pt_loc_ptr, d_shapef_coords.data());
    mesh_me.determinant(
      numCoordDims,           // Number of coordinate dimensions
      num_intg_pts,           // Number of target points
      num_coord_dofs,         // Number of coord shape functions
      d_shapef_coords.data(), // Mesh shape function derivatives
      1,                      // Number of elements
      mesh_coords.data(),        // Mesh coordinate values
      det_J.data(),              // Determinant of the transformation Jacobian for each element (output)
      &det_J_error );         // Determinant error (output)
  }

  const double * intg_wt_ptr = me.intg_weights();

  intg_weights.resize(num_intg_pts);
  for(unsigned i=0; i<num_intg_pts; ++i) intg_weights[i] = intg_wt_ptr[i]*det_J[i];
}

void
ElementObj::gather_nodal_field(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef field, std::vector<double> & result, unsigned dim, unsigned npe)
{
  if (npe == 0)
  {
    npe = mesh.bucket(element).topology().num_nodes();
  }
  result.resize(npe*dim);

  const stk::mesh::Entity* elem_nodes = mesh.begin_nodes(element);

  auto fieldData = field.data<double>();
  for (unsigned n=0; n<npe; ++n)
  {
    auto node_data = fieldData.entity_values(elem_nodes[n]);
    for (unsigned d=0; d<dim; ++d)
    {
      result[n*dim+d] = node_data(stk::mesh::ComponentIdx(d));
    }
  }
}

double
ElementObj::volume(ParametricElementLibrary & parametricElemLibrary, const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coords_field)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<double> coords;
  ElementObj::gather_nodal_field(mesh, element, coords_field, coords, dim);

  const ParametricElement& parametric_elem = parametricElemLibrary.get_parametric_element(mesh.bucket(element).topology());

  std::vector<double> intg_weights;
  ElementObj::integration_weights( intg_weights, dim, coords, parametric_elem );

  double vol = 0.0;
  for ( double intg_wt : intg_weights )
  {
    vol += intg_wt;
  }
  return vol;
}

PhaseTag
ElementObj::update_phase(const std::vector<Surface_Identifier> & surfaceIDs, const PhaseTag & startPhase, const InterfaceID interface_key, const int sign)
{
  PhaseTag phase = startPhase;
  if (sign == -1 || sign == 1)
  {
    if (interface_key.is_single_ls())
    {
      phase.add(surfaceIDs[interface_key.first_ls()],sign);
    }
    else
    {
      const int other_ls_id = (sign < 0) ? interface_key.second_ls() : interface_key.first_ls();
      if (phase.contain(surfaceIDs[other_ls_id], -1))
      {
        const int ls_id = (sign < 0) ? interface_key.first_ls() : interface_key.second_ls();

        phase.clear();
        phase.add(surfaceIDs[ls_id], -1);
      }
    }
  }

  return phase;
}

PhaseTag
ElementObj::update_phase(const std::vector<Surface_Identifier> & surfaceIDs, const PhaseTag & startPhase, const std::vector<InterfaceID> & interfaces, const std::vector<int> & interfaceSigns)
{
  PhaseTag phase = startPhase;
  const int numInterfaces = interfaces.size();
  STK_ThrowRequire((int)interfaceSigns.size() == numInterfaces);

  if (numInterfaces > 0)
  {
    bool badConfig = false;
    int iterCount = 0;
    while (true)
    {
      PhaseTag iterationPhase = phase;
      for (int i=0; i<numInterfaces; ++i)
        phase = update_phase(surfaceIDs, phase, interfaces[i], interfaceSigns[i]);
      if (phase == iterationPhase)
        break;

      if (iterCount++ > numInterfaces)
      {
        krinolog << "BAD INTERFACE CROSSING CONFIG WITH STARTING PHASE: " << startPhase << stk::diag::dendl;
        for (int i=0; i<numInterfaces; ++i)
          krinolog << "  " << interfaces[i] << " " << interfaceSigns[i] << stk::diag::dendl;
        badConfig = true;
        break;
      }
    }

    if (badConfig)
    {
      phase = startPhase;
      for (int i=0; i<numInterfaces; ++i)
        phase = update_phase(surfaceIDs, phase, interfaces[i], interfaceSigns[i]);
    }
  }

  return phase;
}

void ElementObj::fill_node_owner_coords(const Mesh_Element * owner, std::vector<stk::math::Vector3d> & nodeOwnerCoords) const
{
  nodeOwnerCoords.clear();
  nodeOwnerCoords.reserve(my_nodes.size());

  for ( auto && node : my_nodes )
    nodeOwnerCoords.push_back(node->owner_coords(owner));
}

stk::math::Vector3d
ElementObj::compute_local_coords_from_owner_coordinates(const Mesh_Element * owner, const stk::math::Vector3d & ptOwnerCoords) const
{
  std::vector<stk::math::Vector3d> nodeOwnerCoords;
  fill_node_owner_coords(owner, nodeOwnerCoords);

  STK_ThrowAssert((spatial_dim() == 2 && (nodeOwnerCoords.size() == 3 || nodeOwnerCoords.size() == 6)) ||
      (spatial_dim() == 3 && (nodeOwnerCoords.size() == 4 || nodeOwnerCoords.size() == 10)));
  if (spatial_dim() == 2 && nodeOwnerCoords.size() == 6)
    nodeOwnerCoords.resize(3);
  else if (spatial_dim() == 3 && nodeOwnerCoords.size() == 10)
    nodeOwnerCoords.resize(4);

  return get_parametric_coordinates_of_point(nodeOwnerCoords, ptOwnerCoords);
}

void ElementObj::add_subelement(std::unique_ptr<SubElement> subelem)
{
  my_subelements.emplace_back(std::move(subelem));
}

bool ElementObj::fill_path_to_leaf(const SubElement * leaf, std::vector<ElementObj *> & path)
{
  path.push_back(this);
  for (auto && subelem : my_subelements)
    if (subelem.get() == leaf || subelem->fill_path_to_leaf(leaf, path)) return true;
  path.pop_back();
  return false;
}

void ElementObj::remove_leaf_subelements(const std::set<const SubElement *> & leaves)
{
  std::set<const SubElement *> emptied;
  for (auto && subelem : my_subelements)
  {
    if (leaves.count(subelem.get()) != 0 || !subelem->have_subelements()) continue;
    subelem->remove_leaf_subelements(leaves);
    if (!subelem->have_subelements()) emptied.insert(subelem.get());
  }
  std::erase_if(my_subelements,
      [&](const std::unique_ptr<SubElement> & subelem)
      { return leaves.count(subelem.get()) != 0 || emptied.count(subelem.get()) != 0; });
}

void ElementObj::replace_leaf_subelements(const std::vector<const SubElement *> & oldLeaves,
    std::vector<std::unique_ptr<SubElement>> newLeaves)
{
  STK_ThrowRequire(!oldLeaves.empty());
  std::vector<ElementObj *> commonPath;
  for (auto * leaf : oldLeaves)
  {
    std::vector<ElementObj *> path;
    STK_ThrowRequireMsg(fill_path_to_leaf(leaf, path), "Leaf to replace is not below the element.");
    if (commonPath.empty())
      commonPath = path;
    else
      commonPath.erase(
          std::mismatch(commonPath.begin(), commonPath.end(), path.begin(), path.end()).first,
          commonPath.end());
  }
  ElementObj * ancestor = commonPath.back();
  ancestor->remove_leaf_subelements(
      std::set<const SubElement *>(oldLeaves.begin(), oldLeaves.end()));
  for (auto && leaf : newLeaves)
    ancestor->add_subelement(std::move(leaf));
}

void ElementObj::get_subelements( std::vector<const SubElement *> & subelems ) const
{
  if ( !my_subelements.empty() )
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->get_subelements(subelems);
    }
    return;
  }

  const SubElement * subelem = dynamic_cast<const SubElement *>( this );
  if (NULL != subelem) subelems.push_back(subelem);
}

void ElementObj::get_subelements( std::vector<SubElement *> & subelems )
{
  if ( !my_subelements.empty() )
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->get_subelements(subelems);
    }
    return;
  }

  SubElement * subelem = dynamic_cast<SubElement *>( this );
  if (NULL != subelem) subelems.push_back(subelem);
}

bool
ElementObj::have_refined_edges() const
{ /* %TRACE% */  /* %TRACE% */

  const stk::topology Top = topology();
  const int num_edges = Top.num_edges();

  for ( int edge = 0; edge < num_edges; edge++ )
  {
    const unsigned * const lnn = get_edge_node_ordinals(Top, edge);

    const int num_edge_nodes = Top.edge_topology(edge).num_nodes();
    STK_ThrowRequire(2 == num_edge_nodes || 3 == num_edge_nodes);

    if ((2 == num_edge_nodes &&
         NULL != SubElementNode::common_child({my_nodes[lnn[0]], my_nodes[lnn[1]]})) ||
        (3 == num_edge_nodes &&
         (NULL != SubElementNode::common_child({my_nodes[lnn[0]], my_nodes[lnn[2]]}) ||
          NULL != SubElementNode::common_child({my_nodes[lnn[1]], my_nodes[lnn[2]]}))))
    {
      return true;
    }
  }
  return false;
}

void ElementObj::cut_interior_intersection_point(CDMesh & /*mesh*/, const stk::math::Vector3d & /*pCoords*/, const std::vector<int> & /*sortedDomains*/)
{
  throw std::runtime_error("Incorrect usage of ElementObj.  The type of element cannot cut_interior_intersection_point.");
}

void ElementObj::cut_face_intersection_point(const int /*iFace*/, const stk::math::Vector3d & /*pCoords*/, const std::vector<int> & /*sortedDomains*/)
{
  throw std::runtime_error("Incorrect usage of ElementObj.  The type of element cannot cut_face_intersection_point.");
}

bool ElementObj::captures_intersection_point_domains(const std::vector<int> & intersectionPointDomains) const
{
  for (auto && node : my_nodes)
    if (node->captures_intersection_point_domains(intersectionPointDomains))
      return true;
  return false;
}

void
ElementObj::prolongate_fields(const CDMesh & mesh) const
{
  const FieldSet & elementFields = mesh.get_element_fields();
  if (elementFields.empty()) return;

  const ProlongationElementData * prolong_element = nullptr;
  const SubElement * subelem = dynamic_cast<const SubElement * >(this);
  if (nullptr == subelem)
    prolong_element = mesh.fetch_prolong_element(entityId());
  else
    prolong_element = mesh.fetch_prolong_element(subelem->get_owner().entityId());

  //I think that adaptivity currently can lead to prolong_element == NULL: STK_ThrowAssert(NULL != prolong_element);

  for(auto && field : elementFields)
  {
    if (has_field_data(field, entity()))
    {
      auto val = field.data<double, stk::mesh::ReadWrite>().entity_values(entity());
      const unsigned fieldLength = field.length();
      const double * owner_val = (nullptr == prolong_element) ? nullptr : prolong_element->get_field_data(field);
      if (nullptr == owner_val)
        for (unsigned i=0; i<fieldLength; ++i) val(stk::mesh::ScalarIdx(i)) = 0.;
      else
        for (unsigned i=0; i<fieldLength; ++i) val(stk::mesh::ScalarIdx(i)) = owner_val[i];
    }
  }
}

const ParametricElement*
ElementObj::get_evaluation_parametric_element(const FieldRef field) const
{ /* %TRACE% */  /* %TRACE% */
  // Supports Q1Q1, Q2Q2, and Q2Q1 cases
  const ParametricElement* calc_parametric_elem = &parametric_elem();
  const unsigned full_npe = parametric_elem().get_topology().num_nodes();
  STK_ThrowAssert(field.type_is<double>());

  for ( unsigned n = 0; n < full_npe; n++ )
  {
    if (!has_field_data(field, my_nodes[n]->entity()))
    {
      calc_parametric_elem = &myParametricElemLib.get_parametric_element(parametric_elem().get_topology().base());
      STK_ThrowRequire(n >= calc_parametric_elem->num_nodes());
      break;
    }
  }
  return calc_parametric_elem;
}

void
Mesh_Element::determine_decomposed_elem_phase(const CDMesh & mesh, const std::vector<Surface_Identifier> & surfaceIDs)
{
  if(have_subelements())
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->determine_decomposed_elem_phase(mesh, surfaceIDs);
    }
    // Phase for Mesh_Element with subelements is left empty
    return;
  }
}

static std::function<bool(const std::array<unsigned,4> &)>
cut_element_intersecting_planes_diagonal_picker(const NodeVec & cutElemNodes)
{
  auto intersectingPlanesDiagonalPicker =
  [&cutElemNodes](const std::array<unsigned,4> & faceNodes)
  {
    return ElementObj::determine_diagonal_for_internal_quad_of_cut_tet_from_owner_nodes(cutElemNodes[faceNodes[0]], cutElemNodes[faceNodes[1]], cutElemNodes[faceNodes[2]], cutElemNodes[faceNodes[3]]);
  };
  return intersectingPlanesDiagonalPicker;
}

std::function<bool(const std::array<unsigned,4> &)>
Mesh_Element::get_diagonal_picker() const
{
  return cut_element_intersecting_planes_diagonal_picker(my_nodes);
}


bool
Mesh_Element::is_single_coincident() const
{
  std::vector<const SubElement *> conformal_subelems;
  get_subelements(conformal_subelems);
  if(conformal_subelems.size() != 1) return false;

  const SubElement * subelem = conformal_subelems[0];
  if(subelem->topology().num_nodes() != coord_topology().num_nodes()) return false;
  for (auto && node : get_nodes())
    if (std::find(subelem->get_nodes().begin(), subelem->get_nodes().end(), node) == subelem->get_nodes().end())
      return false;
  return true;
}

Mesh_Element::~Mesh_Element() {}

Mesh_Element::Mesh_Element(CDMesh & mesh,
                           stk::mesh::Entity elemEntity)
    : ElementObj(mesh.aux_meta().parametric_element_library(), mesh.stk_bulk(), elemEntity),
      my_subelement_order(0),
      my_have_interface(false)
{
  // set subelement specs
  const stk::topology me_topology = my_parametric_elem.get_topology();

  std::tie(my_subelement_topology, my_subelement_order) = determine_subelement_topology(me_topology);
  STK_ThrowRequire(stk::topology::INVALID_TOPOLOGY != my_subelement_topology);

  const unsigned npe_coords = me_topology.num_nodes();
  const stk::mesh::Entity* elem_nodes = mesh.stk_bulk().begin_nodes(my_entity);

  my_nodes.reserve(npe_coords);
  const unsigned npe_base = me_topology.base().num_nodes();
  for (unsigned i=0; i<npe_base; ++i)
  {
    const SubElementNode * node = mesh.create_mesh_node(this, i, elem_nodes[i]);
    my_nodes.push_back(node);
  }
  if (npe_coords > npe_base)
  {
    my_nodes.resize(npe_coords);
    STK_ThrowRequireMsg(npe_coords == npe_base + me_topology.num_edges(), "Unexpected topology");
    for (unsigned edge_i=0; edge_i<me_topology.num_edges(); ++edge_i) // higher order nodes
    {
      const unsigned * edge_lnn = get_edge_node_ordinals(me_topology, edge_i);
      const unsigned inode = edge_lnn[2];
      STK_ThrowAssert(inode >= npe_base && inode<npe_coords);
      my_nodes[inode] = mesh.create_midside_node(this, my_nodes[edge_lnn[0]], my_nodes[edge_lnn[1]], elem_nodes[inode]);
    }
  }
}

std::pair<stk::topology, unsigned>
Mesh_Element::determine_subelement_topology(stk::topology elem_topology)
{
  switch(elem_topology())
  {
    case stk::topology::TRIANGLE_3_2D:
        return std::pair<stk::topology, unsigned>(stk::topology::TRIANGLE_3_2D, 1);
    case stk::topology::TRIANGLE_6_2D:
        return std::pair<stk::topology, unsigned>(stk::topology::TRIANGLE_6_2D, 2);
    case stk::topology::TRIANGLE_3:
    case stk::topology::SHELL_TRIANGLE_3:
    case stk::topology::SHELL_TRIANGLE_3_ALL_FACE_SIDES:
        return std::pair<stk::topology, unsigned>(stk::topology::TRIANGLE_3, 1); // Using side topology instead of shell topology for shells so that it looks just like lower dimensional element wrt to sides
    case stk::topology::TRIANGLE_6:
    case stk::topology::SHELL_TRIANGLE_6:
    case stk::topology::SHELL_TRIANGLE_6_ALL_FACE_SIDES:
        return std::pair<stk::topology, unsigned>(stk::topology::TRIANGLE_6, 2); // Using side topology instead of shell topology for shells so that it looks just like lower dimensional element wrt to sides
    case stk::topology::TETRAHEDRON_4:
        return std::pair<stk::topology, unsigned>(stk::topology::TETRAHEDRON_4, 1);
    case stk::topology::TETRAHEDRON_10:
        return std::pair<stk::topology, unsigned>(stk::topology::TETRAHEDRON_10, 2);
    default:
        return std::pair<stk::topology, unsigned>(stk::topology::INVALID_TOPOLOGY, 0);
  }
}

stk::math::Vector3d
Mesh_Element::get_node_parametric_coords( const int lnn ) const
{
  const double * nodal_parametric_coordinates = my_parametric_elem.nodal_parametric_coordinates();
  const int dim = spatial_dim();
  return stk::math::Vector3d(&nodal_parametric_coordinates[lnn*dim],dim);
}

static int get_local_node_number( const NodeVec & nodes, const SubElementNode * node )
{
  for (size_t iNode=0; iNode<nodes.size(); ++iNode)
    if (nodes[iNode] == node)
      return iNode;
  STK_ThrowRequireMsg(false, "Failed to find local node.");
  return -1;
}

stk::math::Vector3d
Mesh_Element::get_node_parametric_coords( const SubElementNode * node ) const
{
  STK_ThrowAssertMsg(!get_nodes().empty(), "Attempt to use get_node_parametric_coords before NodeVec filled.");
  return get_node_parametric_coords(get_local_node_number(get_nodes(), node));
}

void Mesh_Element::find_child_coordinates_at_owner_coordinates(const stk::math::Vector3d & ownerCoordinates, const ElementObj *& child, stk::math::Vector3d & childPCoords) const
{
  if (!have_subelements())
  {
    child = this;
    childPCoords = ownerCoordinates;
    return;
  }

  std::vector<const SubElement *> conformal_subelems;
  get_subelements(conformal_subelems);

  double minSqrDist = std::numeric_limits<double>::max();
  for ( auto&& subelement : conformal_subelems )
  {
    const stk::math::Vector3d currentChildPCoords = subelement->compute_local_coords_from_owner_coordinates(this, ownerCoordinates);
    const double currentChildSqrDist = compute_parametric_square_distance(currentChildPCoords);
    if (currentChildSqrDist < minSqrDist)
    {
      minSqrDist = currentChildSqrDist;
      child = subelement;
      childPCoords = currentChildPCoords;
    }
  }
}

stk::math::Vector3d
Mesh_Element::coordinates( const stk::math::Vector3d & p_coords ) const
{ /* %TRACE% */  /* %TRACE% */

  const int npeCoords = my_nodes.size();
  std::vector<double> shapeFcn(npeCoords);
  my_parametric_elem.shape_fcn(1, p_coords.data(), shapeFcn.data());

  stk::math::Vector3d coords(stk::math::Vector3d::ZERO);
  for ( int n = 0; n < npeCoords; n++ )
    coords += shapeFcn[n] * my_nodes[n]->coordinates();

  return coords;
}

static double mean_ratio_of_simplex_for_journal(const std::vector<stk::math::Vector3d> & coords)
{
  return (coords.size() == 4) ? MeanRatioQualityMetric::tet_mean_ratio(coords) : MeanRatioQualityMetric::tri2d_mean_ratio(coords);
}

static void add_simplex_to_journal(CubitJournal & journal, const std::vector<stk::math::Vector3d> & coords)
{
  if (coords.size() == 4)
    journal.add_tet({{coords[0], coords[1], coords[2], coords[3]}});
  else
    journal.add_tri({{coords[0], coords[1], coords[2]}});
}

void Mesh_Element::add_to_journal(CubitJournal & journal,
    const std::function<std::string(const PhaseTag &)> & phaseName,
    const double sliverMeanRatio) const
{
  const std::string prefix = "element_" + std::to_string(entityId());
  const NodeVec parentNodes = vertex_nodes();
  const size_t numVertices = parentNodes.size();
  auto vertex_coordinates = [numVertices](const NodeVec & nodes)
  {
    std::vector<stk::math::Vector3d> coords;
    for (size_t i=0; i<numVertices; ++i)
      coords.push_back(nodes[i]->coordinates());
    return coords;
  };

  journal.set_group(prefix + "_parent", "grey");
  std::ostringstream parentComment;
  parentComment << "Element " << entityId() << " with nodes";
  for (auto * node : parentNodes)
    parentComment << " " << node->entityId();
  journal.add_comment(parentComment.str());
  add_simplex_to_journal(journal, vertex_coordinates(parentNodes));

  std::vector<const ElementObj *> leaves;
  std::vector<const SubElement *> subelements;
  get_subelements(subelements);
  leaves.assign(subelements.begin(), subelements.end());
  if (leaves.empty())
    leaves.push_back(this);

  static const std::array<std::string,6> phaseColors{{"blue", "green", "cyan", "yellow", "orange", "purple"}};
  std::vector<std::string> phaseGroups;
  for (auto * leaf : leaves)
  {
    const std::vector<stk::math::Vector3d> coords = vertex_coordinates(leaf->get_nodes());
    const double meanRatio = mean_ratio_of_simplex_for_journal(coords);
    if (meanRatio < sliverMeanRatio)
    {
      journal.set_group(prefix + "_slivers", "red");
    }
    else
    {
      const std::string group = prefix + "_" + phaseName(leaf->get_phase());
      auto iter = std::find(phaseGroups.begin(), phaseGroups.end(), group);
      if (iter == phaseGroups.end())
        iter = phaseGroups.insert(phaseGroups.end(), group);
      journal.set_group(group, phaseColors[std::distance(phaseGroups.begin(), iter) % phaseColors.size()]);
    }
    std::ostringstream leafComment;
    leafComment << "leaf with mean ratio " << meanRatio << ", phase " << leaf->get_phase();
    journal.add_comment(leafComment.str());
    add_simplex_to_journal(journal, coords);
  }
}

std::string Mesh_Element::visualize(const CDMesh & mesh) const
{
  std::ostringstream os;
  if (myCutter)
    os << "Debugging mesh element " << entityId() << " with cutter:\n" << myCutter->visualize(mesh.stk_bulk()) << "\n";
  else
    os << "Debugging mesh element " << entityId() << " without cutter.\n";
  return os.str();
}

double
Mesh_Element::interface_crossing_position(const InterfaceID interface, const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const
{
  STK_ThrowRequire(get_cutter());
  const auto [crossingSign, position] = get_cutter()->interface_edge_crossing_sign_and_position(interface, edgeNodeCoords);
  STK_ThrowRequireMsg(crossingSign!=0, "Request for interface_crossing_position on edge without crossing.");
  return position;
}

static ElementIntersectionPointFilter build_element_intersection_filter(const NodeVec & nodes)
{
  auto filter =
  [&nodes](const std::vector<int> & intersectionPointDomains)
  {
    for (auto && node : nodes)
      if (node->captures_intersection_point_domains(intersectionPointDomains))
        return false;
    return true;
  };
  return filter;
}

void
Mesh_Element::fill_face_interior_intersections(const NodeVec & faceNodes, std::vector<ElementIntersection> & faceIntersectionPoints) const
{
  STK_ThrowRequire(get_cutter() && faceNodes.size() == 3);
  const std::array<stk::math::Vector3d,3> faceNodeOwnerCoords = {{faceNodes[0]->owner_coords(this), faceNodes[1]->owner_coords(this), faceNodes[2]->owner_coords(this)}};
  const ElementIntersectionPointFilter intersectionPointFilter = build_element_intersection_filter(faceNodes);
  get_cutter()->fill_tetrahedron_face_interior_intersections(faceNodeOwnerCoords, intersectionPointFilter, faceIntersectionPoints);
}

std::pair<int, double>
Mesh_Element::interface_edge_crossing_sign_and_position(const InterfaceID interface, const SubElementNode * node1, const SubElementNode * node2) const
{
  STK_ThrowRequire(get_cutter());
  std::array<stk::math::Vector3d,2> edgeNodeCoords{node1->owner_coords(this), node2->owner_coords(this)};
  return get_cutter()->interface_edge_crossing_sign_and_position(interface, edgeNodeCoords);
}

void Mesh_Element::fill_multiple_edge_crossing_positions(const InterfaceID interface,
    const SubElementNode * node1,
    const SubElementNode * node2,
    std::vector<double> & positions) const
{
  STK_ThrowRequire(get_cutter());
  get_cutter()->fill_multiple_edge_crossing_positions(
      interface, {{node1->coordinates(), node2->coordinates()}}, positions);
}

bool
Mesh_Element::is_prolonged() const
{ /* %TRACE% */  /* %TRACE% */
  for (auto && node : my_nodes)
  {
    if (!(node->is_prolonged())) return false;
  }
  return true;
}

int Mesh_Element::get_interface_index(const InterfaceID interface) const
{
  STK_ThrowAssert(have_interface(interface));
  const auto iter = std::lower_bound(myCuttingInterfaces.begin(), myCuttingInterfaces.end(), interface);
  return std::distance(myCuttingInterfaces.begin(), iter);
}

int Mesh_Element::get_interface_sign_for_uncrossed_subelement(const InterfaceID interface, const std::vector<stk::math::Vector3d> & elemNodeCoords) const
{
   return get_cutter()->interface_sign_for_uncrossed_element(interface, elemNodeCoords);
}


bool
Mesh_Element::triangulate(const CDMesh & mesh, const InterfaceGeometry & interfaceGeometry)
{ /* %TRACE% */  /* %TRACE% */

  const auto & decomposedBlocksSelector = mesh.get_phase_support().get_all_decomposed_blocks_selector();
  const bool inDecomposedBlock = decomposedBlocksSelector(mesh.stk_bulk().bucket(entity()));
  if (inDecomposedBlock)
  {
    if (get_phase().empty())
    {
      create_cutter(mesh, interfaceGeometry);

      if (!myCuttingInterfaces.empty())
        my_have_interface = true;

      set_registered_phase(mesh.get_registered_phase(interfaceGeometry.get_starting_phase(myCutter.get())));

      // The cutter is only used for elements that have interfaces, so free it for the (typically many) uncut elements.
      if (myCuttingInterfaces.empty())
        myCutter.reset();
    }
    else
    {
      my_have_interface = false; // uncut element with phase already set
    }
  }

  if (!have_subelements() && (have_interface() || have_refined_edges()))
  {
    create_base_subelement();
  }

  return false;
}

const SubElementNode * get_node_matching_entity(const std::vector<const SubElementNode *> & nodes, stk::mesh::Entity stkNode)
{
  for (auto && node : nodes)
    if (stkNode == node->entity())
      return node;
  return nullptr;
}

static IntersectionPointFilter
keep_all_intersecion_points_filter()
{
  auto filter =
  [](const std::vector<stk::mesh::Entity> & /*intersectionPointNodes*/, const std::vector<int> & /*intersectionPointSortedDomains*/)
  {
    return true;
  };
  return filter;
}

stk::math::Vector3d Mesh_Element::get_intersection_point_parametric_coordinates(const IntersectionPoint & intersectionPoint) const
{
  stk::math::Vector3d pCoords(stk::math::Vector3d::ZERO);
  const auto & intersectionPointNodes = intersectionPoint.get_nodes();
  const auto & intersectionPointWeights = intersectionPoint.get_weights();
  const size_t numIntersectionPointNodes = intersectionPointNodes.size();
  for (size_t iNode=0; iNode<numIntersectionPointNodes; ++iNode)
  {
    const SubElementNode * node = get_node_matching_entity(my_nodes, intersectionPointNodes[iNode]);
    STK_ThrowAssert(node);
    pCoords += intersectionPointWeights[iNode] * node->owner_coords(this);
  }
  return pCoords;
}

std::vector<int>
Mesh_Element::get_interface_signs_based_on_crossings(const NodeVec & nodes) const
{
  std::vector<stk::math::Vector3d> nodeCoords;
  std::vector<const std::vector<int>*> nodeDomains;

  nodeCoords.clear();
  nodeDomains.clear();
  nodeCoords.reserve(nodes.size());
  nodeDomains.reserve(nodes.size());
  for (auto && node : nodes)
  {
    nodeCoords.push_back(node->owner_coords(this));
    nodeDomains.push_back(&node->get_sorted_node_domains());
  }
  return get_cutter()->get_interface_signs_based_on_crossings(nodeCoords, nodeDomains);
}

void
Mesh_Element::cut_interior_intersection_points(ParametricElementLibrary & parametricElemLibrary, CDMesh & mesh)
{ /* %TRACE% */  /* %TRACE% */
  if (!have_interface())
    return;

  std::vector<IntersectionPoint> intersectionPoints;
  stk::topology elementTopology = mesh.stk_bulk().bucket(entity()).topology();
  const ParametricElement & parametricElement = parametricElemLibrary.get_parametric_element(elementTopology);
  const std::vector<stk::mesh::Entity> elementNodes(mesh.stk_bulk().begin_nodes(entity()), mesh.stk_bulk().end_nodes(entity()));
  append_intersection_points_from_element_interior(parametricElement, elementNodes, *myCutter, keep_all_intersecion_points_filter(), intersectionPoints);

  for (auto && intersectionPoint : intersectionPoints)
  {
    const stk::math::Vector3d pCoords = get_intersection_point_parametric_coordinates(intersectionPoint);

    const ElementObj * containingElem = nullptr;
    stk::math::Vector3d containingElemPCoords;
    find_child_coordinates_at_owner_coordinates(pCoords, containingElem, containingElemPCoords);
    STK_ThrowAssert(containingElem);

    ElementObj * elem = const_cast<ElementObj *>(containingElem);
    if (!elem->captures_intersection_point_domains(intersectionPoint.get_sorted_domains()))
      elem->cut_interior_intersection_point(mesh, containingElemPCoords, intersectionPoint.get_sorted_domains());
  }

  for (auto && subelem : my_subelements)
    subelem->cut_face_interior_intersection_points(mesh);

  std::vector<const SubElement *> leafSubElements;
  get_subelements(leafSubElements);
  for (auto && leafSubElement : leafSubElements)
  {
    SubElement * subElem = const_cast<SubElement *>(leafSubElement);
    const std::vector<int> interfaceSigns = get_interface_signs_based_on_crossings(subElem->get_nodes());
    subElem->set_interface_signs(interfaceSigns);
  }
}

void Mesh_Element::cut_multiple_edge_crossings(CDMesh & mesh)
{
  if (!have_interface() || !get_cutter()) return;

  // Split the leaves one crossing at a time, until no leaf edge has a crossing that its end signs
  // do not show: a second crossing of an edge is then on the edge between the first and its end.
  for (bool cutThisPass = true; cutThisPass;)
  {
    cutThisPass = false;
    std::vector<SubElement *> leafSubElements;
    get_subelements(leafSubElements);
    for (auto * leafSubElement : leafSubElements)
      cutThisPass |= leafSubElement->cut_first_multiple_edge_crossing(mesh);
  }
}

void
Mesh_Element::create_cutter(const CDMesh & mesh, const InterfaceGeometry & interfaceGeometry)
{ /* %TRACE% */  /* %TRACE% */
  const auto intersectingPlanesDiagonalPicker = get_diagonal_picker();
  myCutter = interfaceGeometry.build_element_cutter(mesh.stk_bulk(), entity(), intersectingPlanesDiagonalPicker);
  myCuttingInterfaces = myCutter->get_sorted_cutting_interfaces();
}

void Mesh_Element::split_leaves_at_node(const NodeVec & entityNodes, const SubElementNode * node)
{
  std::vector<SubElement *> leaves;
  get_subelements(leaves);
  for (auto * leaf : leaves)
  {
    const NodeVec & leafNodes = leaf->get_nodes();
    std::vector<int> entityNodeOrdinals;
    for (auto * entityNode : entityNodes)
    {
      const auto iter = std::find(leafNodes.begin(), leafNodes.end(), entityNode);
      if (iter == leafNodes.end()) break;
      entityNodeOrdinals.push_back(std::distance(leafNodes.begin(), iter));
    }
    if (entityNodeOrdinals.size() != entityNodes.size()) continue;

    if (entityNodes.size() == leafNodes.size())
    {
      leaf->split_interior_at_node(node);
    }
    else if (entityNodes.size() == 3)
    {
      const int oppositeNode =
          6 - entityNodeOrdinals[0] - entityNodeOrdinals[1] - entityNodeOrdinals[2];
      leaf->split_face_at_node(oppositeNode, node);
    }
    else
    {
      STK_ThrowRequire(entityNodes.size() == 2);
      leaf->bisect_edge(entityNodeOrdinals[0], entityNodeOrdinals[1], node);
    }
  }
}

void Mesh_Element::insert_intersection_point(CDMesh & mesh, const IntersectionPoint & reducedPoint)
{
  STK_ThrowRequire(is_cut_by_rank());

  // Without subelements, the element itself is the leaf the point is in.
  const stk::math::Vector3d ownerCoords =
      get_intersection_point_parametric_coordinates(reducedPoint);
  const ElementObj * containingElem = nullptr;
  stk::math::Vector3d leafPCoords;
  find_child_coordinates_at_owner_coordinates(ownerCoords, containingElem, leafPCoords);
  STK_ThrowRequire(containingElem);

  // The tolerance was applied on this element, when the point was reduced: on the leaf, a weight is
  // zero only to roundoff.
  constexpr double onLeafEntityTol = 1.e-10;
  // The leaves are linear, as is the element itself without them, on its vertices (the first nodes
  // of a quadratic element).
  const NodeVec leafNodes = (containingElem == this) ? vertex_nodes() : containingElem->get_nodes();
  const int numLeafNodes = leafNodes.size();
  STK_ThrowRequire(numLeafNodes <= 4);
  std::vector<int> entityOrdinals;
  std::vector<double> entityWeights;
  fill_entity_of_weights_above_tolerance(simplex_weights_at_parametric_coordinates(leafPCoords, numLeafNodes),
      onLeafEntityTol,
      entityOrdinals,
      entityWeights);
  STK_ThrowRequire(!entityOrdinals.empty());
  NodeVec entityNodes;
  for (const int ordinal : entityOrdinals)
    entityNodes.push_back(leafNodes[ordinal]);

  const SubElementNode * node = nullptr;
  if (entityNodes.size() == 1)
    node = entityNodes[0];
  else if (entityNodes.size() == 2)
    node = mesh.create_edge_node(this, entityNodes[0], entityNodes[1], entityWeights[1]);
  else
    node = mesh.create_child_internal_or_face_node(this, entityNodes, entityWeights);
  node->add_node_domains_counting_repeats(reducedPoint.get_sorted_domains());

  if (entityNodes.size() > 1)
  {
    if (!have_subelements()) create_base_subelement();
    split_leaves_at_node(entityNodes, node);
  }
}

void Mesh_Element::reset_leaf_interface_signs()
{
  if (!get_cutter()) return;
  std::vector<SubElement *> leaves;
  get_subelements(leaves);
  for (auto * leaf : leaves)
    leaf->set_interface_signs(get_interface_signs_based_on_crossings(leaf->get_nodes()));
}

std::vector<stk::mesh::Entity> Mesh_Element::vertex_entities(const stk::mesh::BulkData & mesh) const
{
  const stk::mesh::Entity * elemNodes = mesh.begin_nodes(entity());
  return std::vector<stk::mesh::Entity>(elemNodes, elemNodes + vertex_topology().num_nodes());
}

std::vector<NodeVec> Mesh_Element::leaf_node_lists() const
{
  if (!have_subelements()) return {vertex_nodes()};
  std::vector<const SubElement *> leaves;
  get_subelements(leaves);
  std::vector<NodeVec> leavesNodes;
  for (auto * leaf : leaves)
    leavesNodes.push_back(leaf->get_nodes());
  return leavesNodes;
}

std::vector<IntersectionPoint>
Mesh_Element::get_corner_intersection_points(const CDMesh & mesh) const
{
  if (!have_interface() || !get_cutter() ||
      !get_cutter()->might_have_interior_or_face_intersections())
    return {};

  const stk::mesh::BulkData & stkMesh = mesh.stk_bulk();
  const ParametricElement & parametricElement =
      get_parametric_elem_library().get_parametric_element(vertex_topology());
  const std::vector<stk::mesh::Entity> elementNodes = vertex_entities(stkMesh);
  std::vector<IntersectionPoint> interiorPoints;
  append_intersection_points_from_element_interior(parametricElement,
      elementNodes,
      *myCutter,
      keep_all_intersection_points_filter(),
      interiorPoints);

  // The corners are the points of codimension equal to the spatial dimension: the corners of the
  // surfaces in 3D; in 2D, also the points where three phases meet.  A corner is not dropped for
  // the domains of the nodes near it, as an edge crossing is: another corner of the same surface
  // has the same domains.  A corner at a node is captured by it when inserted.
  const int cornerCodimension = stkMesh.mesh_meta_data().spatial_dimension();
  return reduced_intersection_points_of_codimension(
      stkMesh, interiorPoints, cornerCodimension, mesh.get_snapper().get_edge_tolerance());
}

std::vector<ElementIntersection>
Mesh_Element::get_uncaptured_sharp_edge_crossings_of_tet_face(const NodeVec & faceNodes, const double edgeTol) const
{
  constexpr int sharpEdgeCodimension = 2;
  const ElementIntersectionPointFilter keepAll = [](const std::vector<int> &) { return true; };
  const std::array<stk::math::Vector3d, 3> faceNodeOwnerCoords{{faceNodes[0]->owner_coords(this),
      faceNodes[1]->owner_coords(this),
      faceNodes[2]->owner_coords(this)}};
  const std::array<const std::vector<int> *, 3> faceNodesDomains{{&faceNodes[0]->get_sorted_node_domains(),
      &faceNodes[1]->get_sorted_node_domains(),
      &faceNodes[2]->get_sorted_node_domains()}};
  std::vector<ElementIntersection> faceIntersections;
  get_cutter()->fill_tetrahedron_face_interior_intersections(
      faceNodeOwnerCoords, keepAll, faceIntersections);
  std::erase_if(faceIntersections,
      [&](const ElementIntersection & faceIntersection)
      {
        return sharpEdgeCodimension != sharp_feature_codimension(faceIntersection.sortedDomains) ||
            tet_face_point_is_captured(faceNodesDomains,
                faceIntersection.parametricCoords,
                faceIntersection.sortedDomains,
                edgeTol);
      });
  return faceIntersections;
}

std::vector<IntersectionPoint> Mesh_Element::sharp_edge_crossing_points_of_leaf_faces(
    const std::vector<NodeVec> & leavesNodes,
    const std::vector<stk::mesh::Entity> & elementNodes,
    const double edgeTol) const
{
  std::vector<IntersectionPoint> crossings;
  for (const auto & face : unique_faces_of_tets(leavesNodes))
  {
    const NodeVec faceNodes(face.nodes.begin(), face.nodes.end());
    const std::array<stk::math::Vector3d, 3> faceNodeOwnerCoords{
        {faceNodes[0]->owner_coords(this), faceNodes[1]->owner_coords(this), faceNodes[2]->owner_coords(this)}};
    for (auto & faceIntersection : get_uncaptured_sharp_edge_crossings_of_tet_face(faceNodes, edgeTol))
    {
      const stk::math::Vector3d & faceCoords = faceIntersection.parametricCoords;
      const stk::math::Vector3d ownerCoords = (1. - faceCoords[0] - faceCoords[1]) * faceNodeOwnerCoords[0] +
          faceCoords[0] * faceNodeOwnerCoords[1] + faceCoords[1] * faceNodeOwnerCoords[2];
      crossings.emplace_back(true,
          elementNodes,
          simplex_weights_at_parametric_coordinates(ownerCoords, elementNodes.size()),
          faceIntersection.sortedDomains);
    }
  }
  return crossings;
}

std::vector<IntersectionPoint> Mesh_Element::points_of_leaf_nodes_with_domains(const stk::mesh::BulkData & stkMesh,
    const std::vector<NodeVec> & leavesNodes,
    const std::vector<stk::mesh::Entity> & elementNodes) const
{
  constexpr double roundoffTol = 1.e-12;
  std::vector<IntersectionPoint> nodesWithDomains;
  std::set<const SubElementNode *> nodesDone;
  for (auto & leafNodes : leavesNodes)
  {
    for (auto * node : leafNodes)
    {
      if (node->get_sorted_node_domains().empty() || !nodesDone.insert(node).second) continue;
      const stk::math::Vector3d & ownerCoords = node->owner_coords(this);
      nodesWithDomains.push_back(reduced_intersection_point_in_canonical_form(stkMesh,
          IntersectionPoint(true,
              elementNodes,
              simplex_weights_at_parametric_coordinates(ownerCoords, elementNodes.size()),
              node->get_sorted_node_domains()),
          roundoffTol));
    }
  }
  return nodesWithDomains;
}

std::vector<IntersectionPoint> Mesh_Element::get_edge_intersection_points(const CDMesh & mesh) const
{
  if (!have_interface() || !get_cutter() ||
      !get_cutter()->might_have_interior_or_face_intersections() ||
      vertex_topology() != stk::topology::TETRAHEDRON_4)
    return {};

  const stk::mesh::BulkData & stkMesh = mesh.stk_bulk();
  const std::vector<stk::mesh::Entity> elementNodes = vertex_entities(stkMesh);
  const double edgeTol = mesh.get_snapper().get_edge_tolerance();
  constexpr int sharpEdgeCodimension = 2;
  const std::vector<NodeVec> leavesNodes = leaf_node_lists();

  std::vector<IntersectionPoint> points = reduced_intersection_points_of_codimension(stkMesh,
      sharp_edge_crossing_points_of_leaf_faces(leavesNodes, elementNodes, edgeTol),
      sharpEdgeCodimension,
      edgeTol);

  // A crossing within the edge tolerance of a node of a leaf whose domains capture it is that node,
  // found again: the sharp edge through the node crosses a leaf face beside it, as in a sliver.
  // Cutting it would make a node beside the one there, and its reduction onto the entity of the
  // element it is near would move it off the sharp edge.
  const std::vector<IntersectionPoint> nodesWithDomains =
      points_of_leaf_nodes_with_domains(stkMesh, leavesNodes, elementNodes);
  std::erase_if(points,
      [&](const IntersectionPoint & point)
      {
        return std::any_of(nodesWithDomains.begin(),
            nodesWithDomains.end(),
            [&](const IntersectionPoint & node) { return is_captured_by_intersection_point(node, point, edgeTol); });
      });
  return points;
}

// The largest distance between the mesh nodes the edge's nodes come from: the size of the entity of
// the parent elements that the edge is in, the same in every element with the edge.
static double
extent_of_mesh_ancestors_of_edge(const SubElementNode * node0, const SubElementNode * node1)
{
  NodeSet ancestors;
  node0->get_ancestors(ancestors);
  node1->get_ancestors(ancestors);
  double extent = 0.;
  for (auto * ancestor0 : ancestors)
    for (auto * ancestor1 : ancestors)
      extent = std::max(extent, (ancestor1->coordinates() - ancestor0->coordinates()).length());
  return extent;
}

void Mesh_Element::append_multiple_crossing_points_of_leaf_edge(const SubElementNode * node0,
    const SubElementNode * node1,
    const std::vector<stk::mesh::Entity> & elementNodes,
    const double edgeTol,
    const bool sharpFeaturesAreCut,
    std::vector<IntersectionPoint> & crossings) const
{
  // The tolerances are fractions of the extent of the edge's mesh ancestors, the same in every
  // element with the edge.
  const double extent = extent_of_mesh_ancestors_of_edge(node0, node1);
  std::vector<int> sortedDomains;
  std::vector<double> positions;
  for (auto && interface : get_sorted_cutting_interfaces())
  {
    interface.fill_sorted_domains(sortedDomains);
    if (edge_is_chord_of_interface(node0->get_sorted_node_domains(), node1->get_sorted_node_domains(), sortedDomains))
      continue;
    fill_multiple_edge_crossing_positions(interface, node0, node1, positions);
    if (positions.empty()) continue;
    const double atEndTol =
        multiple_crossing_end_tolerance(sharpFeaturesAreCut && get_cutter()->interface_has_sharp_features(interface), edgeTol) * extent;
    for (const double position : multiple_crossing_positions_to_insert(positions,
             node0->coordinates(),
             node1->coordinates(),
             node0->get_sorted_node_domains(),
             node1->get_sorted_node_domains(),
             sortedDomains,
             edgeTol * extent,
             atEndTol))
    {
      const stk::math::Vector3d ownerCoords =
          (1. - position) * node0->owner_coords(this) + position * node1->owner_coords(this);
      crossings.emplace_back(true,
          elementNodes,
          simplex_weights_at_parametric_coordinates(ownerCoords, elementNodes.size()),
          sortedDomains);
    }
  }
}

std::vector<IntersectionPoint>
Mesh_Element::get_multiple_edge_crossing_points(const CDMesh & mesh) const
{
  if (!have_interface() || !get_cutter() || !is_cut_by_rank()) return {};

  const stk::mesh::BulkData & stkMesh = mesh.stk_bulk();
  const std::vector<stk::mesh::Entity> elementNodes = vertex_entities(stkMesh);
  const double edgeTol = mesh.get_snapper().get_edge_tolerance();
  const stk::topology topology = vertex_topology();
  std::array<unsigned, 2> edgeNodeOrdinals;

  std::set<std::pair<const SubElementNode *, const SubElementNode *>> edgesDone;
  std::vector<IntersectionPoint> crossings;
  for (auto & leafNodes : leaf_node_lists())
  {
    for (unsigned iEdge = 0; iEdge < topology.num_edges(); ++iEdge)
    {
      topology.edge_node_ordinals(iEdge, edgeNodeOrdinals.data());
      const SubElementNode * node0 = leafNodes[edgeNodeOrdinals[0]];
      const SubElementNode * node1 = leafNodes[edgeNodeOrdinals[1]];
      if (edgesDone.insert(std::minmax(node0, node1)).second)
        append_multiple_crossing_points_of_leaf_edge(
            node0, node1, elementNodes, edgeTol, mesh.get_cdfem_support().cuts_sharp_features(), crossings);
    }
  }

  constexpr int surfaceCrossingCodimension = 1;
  constexpr double roundoffTol = 1.e-12;
  return reduced_intersection_points_of_codimension(
      stkMesh, crossings, surfaceCrossingCodimension, roundoffTol);
}

void Mesh_Element::insert_intersection_points(
    CDMesh & mesh, std::vector<IntersectionPoint> points, const double repeatTol)
{
  if (points.empty()) return;

  sort_intersection_points_in_canonical_order(mesh.stk_bulk(), points);
  remove_repeated_intersection_points(points, repeatTol);

  for (auto & point : points)
    insert_intersection_point(mesh, point);
}

using TetNodes = std::array<const SubElementNode *, 4>;
using SortedTriNodes = std::array<const SubElementNode *, 3>;

// The nodes of a face in a canonical order, by ancestry: the same in every leaf with the face and on
// every processor.  The swaps find a face's crossings with its nodes in this order, so that whether
// it is crossed does not depend on the order of its nodes in a leaf, even to roundoff: otherwise a
// swap accepted because its new faces were found uncrossed could find them crossed, and be undone,
// without end.  That roundoff can change the decomposition, so the order does not depend on the ids
// where they are not parallel consistent.
static SortedTriNodes canonical_tri_nodes(const bool globalIdsAreParallelConsistent,
    const SubElementNode * node0,
    const SubElementNode * node1,
    const SubElementNode * node2)
{
  SortedTriNodes nodes{{node0, node1, node2}};
  std::sort(nodes.begin(),
      nodes.end(),
      [globalIdsAreParallelConsistent](const SubElementNode * a, const SubElementNode * b)
      {
        return less_conditionally_by_coordinates_because_it_impacts_structure(
            globalIdsAreParallelConsistent, a->get_ancestry(), b->get_ancestry());
      });
  return nodes;
}

// The nodes of a face sorted by id, as a key: only to find the face again.
static SortedTriNodes tri_nodes_key(
    const SubElementNode * node0, const SubElementNode * node1, const SubElementNode * node2)
{
  SortedTriNodes nodes{{node0, node1, node2}};
  std::sort(nodes.begin(),
      nodes.end(),
      [](const SubElementNode * a, const SubElementNode * b)
      { return less_by_global_id_because_it_doesnt_impact_structure(a->get_ancestry(), b->get_ancestry()); });
  return nodes;
}

// The new leaves take the parent side ids of the old leaves' sides they are on; their interface
// signs are set after the swaps.
void Mesh_Element::replace_tet_leaves(
    const std::vector<const SubElement *> & oldLeaves, const std::vector<TetNodes> & tets)
{
  std::map<SortedTriNodes, int> parentSideIds;
  for (auto * leaf : oldLeaves)
  {
    const NodeVec & leafNodes = leaf->get_nodes();
    for (int iSide = 0; iSide < 4; ++iSide)
    {
      const auto & side = tetSideNodeOrdinals[iSide];
      parentSideIds[tri_nodes_key(leafNodes[side[0]], leafNodes[side[1]], leafNodes[side[2]])] =
          leaf->parent_side_id(iSide);
    }
  }

  std::vector<std::unique_ptr<SubElement>> newLeaves;
  for (auto & tet : tets)
  {
    std::vector<int> sideIds;
    for (auto & side : tetSideNodeOrdinals)
    {
      const auto iter =
          parentSideIds.find(tri_nodes_key(tet[side[0]], tet[side[1]], tet[side[2]]));
      sideIds.push_back(iter == parentSideIds.end() ? -1 : iter->second);
    }
    newLeaves.push_back(
        std::make_unique<SubElement_Tet_4>(NodeVec(tet.begin(), tet.end()), sideIds, this));
  }
  replace_leaf_subelements(oldLeaves, std::move(newLeaves));
}

void Mesh_Element::recover_sharp_edges_in_tet_by_swaps(const CDMesh & mesh)
{
  if (!have_interface() || !get_cutter() ||
      !get_cutter()->might_have_interior_or_face_intersections() ||
      vertex_topology() != stk::topology::TETRAHEDRON_4 || !have_subelements())
    return;

  const double edgeTol = mesh.get_snapper().get_edge_tolerance();
  const bool globalIdsAreParallelConsistent = mesh.get_cdfem_support().get_global_ids_are_parallel_consistent();
  using Recovery = SharpEdgeRecoveryBySwaps<const SubElementNode *>;
  const Recovery recovery(edgeTol,
      [](const SubElementNode * node) { return node->coordinates(); },
      [this](const SubElementNode * node) { return node->owner_coords(this); },
      [](const SubElementNode * node, const std::vector<int> & sortedDomains)
      { return node_domains_capture_point_domains(node->get_sorted_node_domains(), sortedDomains); },
      [globalIdsAreParallelConsistent](
          const SubElementNode * node0, const SubElementNode * node1, const SubElementNode * node2)
      { return canonical_tri_nodes(globalIdsAreParallelConsistent, node0, node1, node2); },
      [this, edgeTol](const Recovery::Face & faceNodes)
      {
        std::vector<Recovery::FaceCrossing> crossings;
        for (auto & crossing :
            get_uncaptured_sharp_edge_crossings_of_tet_face(NodeVec(faceNodes.begin(), faceNodes.end()), edgeTol))
          crossings.push_back({crossing.parametricCoords, crossing.sortedDomains});
        return crossings;
      });

  // Each swap replaces its leaves in place, as they are found, so that the order of the leaves, and
  // so of the swaps after it, is as without the class.
  while (true)
  {
    std::vector<SubElement *> leaves;
    get_subelements(leaves);
    std::vector<TetNodes> tets;
    for (auto * leaf : leaves)
    {
      const NodeVec & leafNodes = leaf->get_nodes();
      tets.push_back({{leafNodes[0], leafNodes[1], leafNodes[2], leafNodes[3]}});
    }
    const std::optional<Recovery::Swap> swap = recovery.find_swap(tets);
    if (!swap) break;
    std::vector<const SubElement *> oldLeaves;
    for (const size_t iTet : swap->oldTets)
      oldLeaves.push_back(leaves[iTet]);
    replace_tet_leaves(oldLeaves, swap->newTets);
  }
}

void
Mesh_Element::create_base_subelement()
{ /* %TRACE% */  /* %TRACE% */

  stk::topology baseTopology = coord_topology();

  const unsigned num_sides = baseTopology.num_sides();

  std::vector<int> parent_side_ids(num_sides);
  for (unsigned i=0; i<num_sides; ++i)
  {
    parent_side_ids[i] = i;
  }

  std::unique_ptr<SubElement> base_subelement;
  if (stk::topology::TETRAHEDRON_4 == baseTopology)
  {
    base_subelement = std::make_unique<SubElement_Tet_4>( my_nodes, parent_side_ids, this);
  }
  else if (stk::topology::TETRAHEDRON_10 == baseTopology)
  {
    // purposely use lower order base subelement
    std::vector<const SubElementNode *> sub_nodes(my_nodes.begin(), my_nodes.begin()+4);
    base_subelement = std::make_unique<SubElement_Tet_4>( sub_nodes, parent_side_ids, this);
  }
  else if (stk::topology::TRIANGLE_3_2D == baseTopology || stk::topology::TRIANGLE_3 == baseTopology)
  {
    base_subelement = std::make_unique<SubElement_Tri_3>( my_nodes, parent_side_ids, this);
  }
  else if (stk::topology::TRIANGLE_6_2D == baseTopology || stk::topology::TRIANGLE_6 == baseTopology)
  {
    // purposely use lower order base subelement
    std::vector<const SubElementNode *> sub_nodes(my_nodes.begin(), my_nodes.begin()+3);
    base_subelement = std::make_unique<SubElement_Tri_3>( sub_nodes, parent_side_ids, this);
  }
  STK_ThrowErrorMsgIf(!base_subelement, "Elements with topology " << baseTopology.name() << " not supported for CDFEM.");

  base_subelement->initialize_interface_signs();
  add_subelement(std::move(base_subelement));
}

void
Mesh_Element::determine_node_signs(const CDMesh & mesh, const InterfaceID interface_key)
{
  if (have_interface(interface_key))
  {
    for(auto && subelement : my_subelements)
    {
      subelement->SubElement::determine_node_signs(mesh, interface_key);
    }
  }
}

void
Mesh_Element::determine_node_scores(const CDMesh & mesh, const InterfaceID interface_key)
{
  if (have_interface(interface_key))
  {
    for(auto && subelement : my_subelements)
    {
      subelement->SubElement::determine_node_scores(mesh, interface_key);
    }
  }
}

void
Mesh_Element::decompose(CDMesh & mesh, const InterfaceID interface_key)
{
  if (have_interface(interface_key))
  {
    for(auto && subelement : my_subelements)
    {
      subelement->decompose(mesh, interface_key);
    }
  }
}

void
Mesh_Element::handle_hanging_children(CDMesh & mesh, const InterfaceID & interface)
{
  if (!have_subelements() && have_edges_with_children())
  {
    create_base_subelement();
  }

  for(auto && subelement : my_subelements)
  {
    subelement->handle_hanging_children(mesh, interface);
  }
}

void
Mesh_Element::build_quadratic_subelements(CDMesh & mesh)
{
  if (2 == subelement_order())
  {
    for (auto && subelement : my_subelements)
    {
      subelement->build_quadratic_subelements(mesh);
    }
  }
}

bool
ElementObj::will_cutting_quad_from_0to2_cut_largest_angle(const SubElementNode * n0, const SubElementNode * n1, const SubElementNode * n2, const SubElementNode * n3)
{
  return krino::will_cutting_quad_from_0to2_cut_largest_angle(n0->coordinates(), n1->coordinates(), n2->coordinates(), n3->coordinates());
}

bool
ElementObj::have_edges_with_children() const
{ /* %TRACE% */  /* %TRACE% */
  const int num_edges = topology().num_edges();

  // Iterate edges looking for any common children of the edge nodes
  for ( int edge = 0; edge < num_edges; ++edge )
  {
    const unsigned * edge_node_ordinals = get_edge_node_ordinals(topology(), edge);

    const SubElementNode * node0 = my_nodes[edge_node_ordinals[0]];
    const SubElementNode * node1 = my_nodes[edge_node_ordinals[1]];
    const SubElementNode * child = SubElementNode::common_child({node0, node1});
    if( child )
    {
      return true;
    }
  }
  return false;
}

bool
ElementObj::determine_diagonal_for_internal_quad_of_cut_tet_from_owner_nodes(const SubElementNode * pn0, const SubElementNode * pn1, const SubElementNode * pn2, const SubElementNode * pn3)
{
  // true:  connect child of pn0 and pn2 (n0) to child of pn1 and pn3 (n2)
  // false: connect child of pn1 and pn2 (n1) to child of pn0 and pn3 (n3)

  // This method is used for the cutter.  For CUT_QUADS_BY_GLOBAL_IDENTIFIER, use logic that is consistent with the way elements are cut.
  // For other strategies, still use id based logic because this is only being used for the cutter, and it will be removed when
  // we use the actual facets instead of the cutter.

  stk::mesh::EntityId pn0_id = pn0->entityId();
  stk::mesh::EntityId pn1_id = pn1->entityId();
  stk::mesh::EntityId pn2_id = pn2->entityId();
  stk::mesh::EntityId pn3_id = pn3->entityId();

  STK_ThrowAssert(pn0_id != 0 && pn1_id != 0 && pn2_id != 0 && pn3_id != 0);

  if ((pn0_id < pn1_id) == (pn2_id < pn3_id))
  {
    return true;
  }
  return false;
}

bool
ElementObj::determine_diagonal_for_internal_quad_of_cut_tet_from_edge_nodes(const Simplex_Generation_Method simplexMethod, const SubElementNode * n0, const SubElementNode * n1, const SubElementNode * n2, const SubElementNode * n3,
    const bool face0, const bool face1, const bool face2, const bool face3)
{
  // true: connect n0 and n2, false: connect n1 and n3

  if (simplexMethod == CUT_QUADS_BY_GLOBAL_IDENTIFIER)
  {
    SubElementNodeAncestry ancestry0 = n0->get_ancestry();
    SubElementNodeAncestry ancestry1 = n1->get_ancestry();
    SubElementNodeAncestry ancestry2 = n2->get_ancestry();
    SubElementNodeAncestry ancestry3 = n3->get_ancestry();

    const auto lessById = [](const SubElementNodeAncestry & a, const SubElementNodeAncestry & b)
    { return less_by_global_id_because_it_is_requested(a, b); };
    SubElementNodeAncestry ancestry02 = std::min(ancestry0, ancestry2, lessById);
    SubElementNodeAncestry ancestry13 = std::min(ancestry1, ancestry3, lessById);

    return lessById(ancestry02, ancestry13);
  }

  // Make sure that diagonal can be picked either way and still avoid Steiner point
  const int caseId = (face0 ? 1 : 0) + (face1 ? 2 : 0) + (face2 ? 4 : 0) + (face3 ? 8 : 0);
  STK_ThrowRequire(caseId == 3 || caseId == 5 || caseId == 10 || caseId == 12);

  return will_cutting_quad_from_0to2_cut_largest_angle(n0,n1,n2,n3);
}

} // namespace krino
