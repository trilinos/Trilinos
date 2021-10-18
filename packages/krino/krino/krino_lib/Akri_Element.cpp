// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Element.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_SubElement.hpp>
#include <Akri_SubElementNodeAncestry.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_ElementCutterUtils.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_ProlongationData.hpp>
#include <Akri_MeshHelpers.hpp>
#include <math.h>

#include "../interface_geometry_interface/Akri_InterfaceGeometry.hpp"
#include <Akri_MasterElementDeterminer.hpp>
namespace krino{

ElementObj::ElementObj(const stk::mesh::BulkData & stkMesh, stk::mesh::Entity elemEntity)
    : my_master_elem(MasterElementDeterminer::getMasterElement(stkMesh.bucket(elemEntity).topology())),
      my_entity(elemEntity),
      my_entityId(stkMesh.identifier(elemEntity)),
      my_prolongation_data(NULL)
{
}

ElementObj::ElementObj(const stk::topology elem_topology, const NodeVec & nodes)
   :  my_master_elem(MasterElementDeterminer::getMasterElement(elem_topology)),
      my_nodes(nodes),
      my_entity(),
      my_entityId(0),
      my_prolongation_data(NULL)
{
}

ElementObj::~ElementObj() {}

void
ElementObj::integration_locations(
    std::vector<Vector3d> & intg_pt_locations,
    const MasterElement & me)
{
  const unsigned num_intg_pts = me.num_intg_pts();
  const double * intg_pt_loc_ptr = me.intg_pt_locations();
  const unsigned dim = me.topology_dimension();

  intg_pt_locations.resize(num_intg_pts);
  for(unsigned i=0; i<num_intg_pts; ++i) intg_pt_locations[i] = Vector3d(intg_pt_loc_ptr+i*dim, dim);
}

void
ElementObj::integration_weights(
    std::vector<double> & intg_weights, // includes both gauss point weight and detJ
    const int numCoordDims,
    const std::vector<double> & mesh_coords,
    const MasterElement & me,
    const MasterElement & mesh_me)
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

  for (unsigned n=0; n<npe; ++n)
  {
    const double * node_data = field_data<double>(field, elem_nodes[n]);
    for (unsigned d=0; d<dim; ++d)
    {
      result[n*dim+d] = node_data[d];
    }
  }
}

double
ElementObj::volume(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coords_field)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<double> coords;
  ElementObj::gather_nodal_field(mesh, element, coords_field, coords, dim);

  const MasterElement& master_elem = MasterElementDeterminer::getMasterElement(mesh.bucket(element).topology());

  std::vector<double> intg_weights;
  ElementObj::integration_weights( intg_weights, dim, coords, master_elem );

  double vol = 0.0;
  for ( double intg_wt : intg_weights )
  {
    vol += intg_wt;
  }
  return vol;
}

PhaseTag
ElementObj::update_phase(const CDMesh & mesh, const PhaseTag & startPhase, const InterfaceID interface_key, const int sign)
{
  PhaseTag phase = startPhase;
  if (sign == -1 || sign == 1)
  {
    if (interface_key.is_single_ls())
    {
      phase.add(mesh.ls_field(interface_key.first_ls()).identifier,sign);
    }
    else
    {
      const int other_ls_id = (sign < 0) ? interface_key.second_ls() : interface_key.first_ls();
      if (phase.contain(mesh.ls_field(other_ls_id).identifier, -1))
      {
        const int ls_id = (sign < 0) ? interface_key.first_ls() : interface_key.second_ls();

        phase.clear();
        phase.add(mesh.ls_field(ls_id).identifier, -1);
      }
    }
  }

  return phase;
}

PhaseTag
ElementObj::update_phase(const CDMesh & mesh, const PhaseTag & startPhase, const std::vector<InterfaceID> & interfaces, const std::vector<int> & interfaceSigns)
{
  PhaseTag phase = startPhase;
  const int numInterfaces = interfaces.size();
  ThrowRequire((int)interfaceSigns.size() == numInterfaces);

  if (numInterfaces > 0)
  {
    bool badConfig = false;
    int iterCount = 0;
    while (true)
    {
      PhaseTag iterationPhase = phase;
      for (int i=0; i<numInterfaces; ++i)
        phase = update_phase(mesh, phase, interfaces[i], interfaceSigns[i]);
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
        phase = update_phase(mesh, phase, interfaces[i], interfaceSigns[i]);
    }
  }

  return phase;
}

void ElementObj::fill_node_owner_coords(const Mesh_Element * owner, std::vector<Vector3d> & nodeOwnerCoords) const
{
  nodeOwnerCoords.clear();
  nodeOwnerCoords.reserve(my_nodes.size());

  for ( auto && node : my_nodes )
    nodeOwnerCoords.push_back(node->owner_coords(owner));
}

Vector3d
ElementObj::compute_local_coords_from_owner_coordinates(const Mesh_Element * owner, const Vector3d & ptOwnerCoords) const
{
  std::vector<Vector3d> nodeOwnerCoords;
  fill_node_owner_coords(owner, nodeOwnerCoords);

  ThrowAssert((spatial_dim() == 2 && (nodeOwnerCoords.size() == 3 || nodeOwnerCoords.size() == 6)) ||
      (spatial_dim() == 3 && (nodeOwnerCoords.size() == 4 || nodeOwnerCoords.size() == 10)));
  if (spatial_dim() == 2 && nodeOwnerCoords.size() == 6)
    nodeOwnerCoords.resize(3);
  else if (spatial_dim() == 3 && nodeOwnerCoords.size() == 10)
    nodeOwnerCoords.resize(4);

  return get_parametric_coordinates_of_point(nodeOwnerCoords, ptOwnerCoords);
}

void
ElementObj::integration_weights(
    std::vector<double> & intg_weights, // includes both gauss point weight and detJ
    const MasterElement & me ) const
{
  const unsigned dim = spatial_dim();
  const unsigned nnodes = my_nodes.size();
  std::vector<double> flat_coords(nnodes*dim);
  for ( unsigned n = 0; n < nnodes; n++ )
  {
    const Vector3d & node_coords = my_nodes[n]->coordinates();
    for (unsigned d=0; d<dim; ++d)
    {
      flat_coords[n*dim+d] = node_coords[d];
    }
  }

  integration_weights(intg_weights, dim, flat_coords, me, my_master_elem);
}

double
Mesh_Element::volume() const
{ /* %TRACE% */  /* %TRACE% */

  std::vector<double> intg_weights;
  integration_weights( intg_weights, my_master_elem );

  double vol = 0.0;

  for ( double intg_weight : intg_weights )
  {
    vol += intg_weight;
  }
  return vol;
}

void ElementObj::add_subelement(std::unique_ptr<SubElement> subelem)
{
  my_subelements.emplace_back(std::move(subelem));
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
    ThrowRequire(2 == num_edge_nodes || 3 == num_edge_nodes);

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

void ElementObj::cut_interior_intersection_point(CDMesh & mesh, const Vector3d & pCoords, const std::vector<int> & sortedDomains)
{
  throw std::runtime_error("Incorrect usage of ElementObj.  The type of element cannot cut_interior_intersection_point.");
}

void ElementObj::cut_face_intersection_point(const int iFace, const Vector3d & pCoords, const std::vector<int> & sortedDomains)
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
  const CDMesh* old_mesh = mesh.get_old_mesh();
  if (nullptr != old_mesh)
  {
    const FieldSet & element_fields = mesh.get_element_fields();
    if (element_fields.empty()) return;

    const ProlongationElementData * prolong_element = NULL;
    const SubElement * subelem = dynamic_cast<const SubElement * >(this);
    if (NULL == subelem)
    {
      prolong_element = old_mesh->fetch_prolong_element(entityId());
    }
    else
    {
      prolong_element = old_mesh->fetch_prolong_element(subelem->get_owner().entityId());
    }

    //I think that adaptivity currently can lead to prolong_element == NULL: ThrowAssert(NULL != prolong_element);

    for(FieldSet::const_iterator it = element_fields.begin(); it != element_fields.end(); ++it)
    {
      const FieldRef field = *it;

      double * val = field_data<double>(field, entity());
      if (NULL == val) continue;

      const unsigned field_length = field.length();

      const double * owner_val = (NULL == prolong_element) ? NULL : prolong_element->get_field_data(field);
      if (NULL == owner_val)
      {
        std::fill(val, val+field_length, 0.);
      }
      else
      {
        std::copy(owner_val, owner_val+field_length, val);
      }
    }
  }
}

const MasterElement*
ElementObj::get_evaluation_master_element(const FieldRef field) const
{ /* %TRACE% */  /* %TRACE% */
  // Supports Q1Q1, Q2Q2, and Q2Q1 cases
  const MasterElement* calc_master_elem = &master_elem();
  const unsigned full_npe = master_elem().get_topology().num_nodes();
  ThrowAssert(field.type_is<double>());

  for ( unsigned n = 0; n < full_npe; n++ )
  {
    double * data = field_data<double>(field, my_nodes[n]->entity());
    if (nullptr == data)
    {
      calc_master_elem = &MasterElementDeterminer::getMasterElement(master_elem().get_topology().base());
      ThrowRequire(n >= calc_master_elem->num_nodes());
      break;
    }
  }
  return calc_master_elem;
}

void
ElementObj::evaluate_prolongation_field(const CDMesh & mesh, const FieldRef field, const unsigned field_length, const Vector3d & p_coords, double * result) const
{ /* %TRACE% */  /* %TRACE% */

  // Figuring out the field master element here is actually quite hard since the entity may not exist any more.
  // We'll assume that the field master element is the master_elem or the one with topology master_elem->get_topology().base().
  // This will handle the Q2Q1 case.

  for (unsigned i=0; i<field_length; ++i) result[i] = 0.0;
  FieldRef initial_field;

  const MasterElement* calc_master_elem = &master_elem();
  const int full_npe = master_elem().get_topology().num_nodes();
  std::vector<const double *> node_data(full_npe, nullptr);
  const std::vector<double> zeros(field_length,0.0);

  for ( int n = 0; n < full_npe; n++ )
  {
    const ProlongationNodeData * prolong_data = mesh.fetch_prolong_node(my_nodes[n]->entityId());
    if (nullptr != prolong_data) node_data[n] = prolong_data->get_field_data(field);
    if (node_data[n] == nullptr)
    {
      if (!initial_field.valid())
      {
        initial_field = mesh.get_cdfem_support().get_initial_prolongation_field( field );
      }
      if (initial_field.valid())
      {
        node_data[n] = prolong_data->get_field_data(initial_field);
      }
    }
    if (node_data[n] == nullptr)
    {
      calc_master_elem = &MasterElementDeterminer::getMasterElement(master_elem().get_topology().base());
      node_data[n] = zeros.data();
    }
  }

  const int npe = calc_master_elem->get_topology().num_nodes();
  std::vector<double> shapefcn (npe,0.);
  calc_master_elem->shape_fcn(1, p_coords.data(), shapefcn.data());

  for ( int n = 0; n < npe; n++ )
  {
    ThrowRequire(nullptr != node_data[n]);
    for (unsigned i=0; i<field_length; ++i) result[i] += shapefcn[n]*node_data[n][i];
  }
}

void
Mesh_Element::determine_decomposed_elem_phase(const CDMesh & mesh)
{
  if(have_subelements())
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->determine_decomposed_elem_phase(mesh);
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
  if(subelem->topology() != coord_topology()) return false;
  if (get_nodes() != subelem->get_nodes()) return false;
  return true;
}

Mesh_Element::~Mesh_Element() {}

Mesh_Element::Mesh_Element(CDMesh & mesh,
                           stk::mesh::Entity elemEntity)
    : ElementObj(mesh.stk_bulk(), elemEntity),
      my_subelement_order(0),
      my_have_interface(false)
{
  // set subelement specs
  const stk::topology me_topology = my_master_elem.get_topology();

  std::tie(my_subelement_topology, my_subelement_order) = determine_subelement_topology(me_topology);
  ThrowRequire(stk::topology::INVALID_TOPOLOGY != my_subelement_topology);

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
    ThrowRequireMsg(npe_coords == npe_base + me_topology.num_edges(), "Unexpected topology");
    for (unsigned edge_i=0; edge_i<me_topology.num_edges(); ++edge_i) // higher order nodes
    {
      const unsigned * edge_lnn = get_edge_node_ordinals(me_topology, edge_i);
      const unsigned inode = edge_lnn[2];
      ThrowAssert(inode >= npe_base && inode<npe_coords);
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
    case stk::topology::TETRAHEDRON_4:
        return std::pair<stk::topology, unsigned>(stk::topology::TETRAHEDRON_4, 1);
    case stk::topology::TETRAHEDRON_10:
        return std::pair<stk::topology, unsigned>(stk::topology::TETRAHEDRON_10, 2);
    default:
        return std::pair<stk::topology, unsigned>(stk::topology::INVALID_TOPOLOGY, 0);
  }
}

Vector3d
Mesh_Element::get_node_parametric_coords( const int lnn ) const
{
  const double * nodal_parametric_coordinates = my_master_elem.nodal_parametric_coordinates();
  const int dim = spatial_dim();
  return Vector3d(&nodal_parametric_coordinates[lnn*dim],dim);
}

static int get_local_node_number( const NodeVec & nodes, const SubElementNode * node )
{
  for (size_t iNode=0; iNode<nodes.size(); ++iNode)
    if (nodes[iNode] == node)
      return iNode;
  ThrowRequireMsg(false, "Failed to find local node.");
  return -1;
}

Vector3d
Mesh_Element::get_node_parametric_coords( const SubElementNode * node ) const
{
  ThrowAssertMsg(!get_nodes().empty(), "Attempt to use get_node_parametric_coords before NodeVec filled.");
  return get_node_parametric_coords(get_local_node_number(get_nodes(), node));
}

static double compute_parametric_square_distance(const Vector3d childPCoords)
{
  double dist2 = 0.0;
  if (childPCoords[0] < 0.0) dist2 += childPCoords[0]*childPCoords[0];
  if (childPCoords[1] < 0.0) dist2 += childPCoords[1]*childPCoords[1];
  if (childPCoords[2] < 0.0) dist2 += childPCoords[2]*childPCoords[2];
  const double zeta = 1.0 - childPCoords[0] - childPCoords[1] - childPCoords[2];
  if (zeta < 0.0) dist2 += zeta*zeta;
  return dist2;
}

void Mesh_Element::find_child_coordinates_at_owner_coordinates(const Vector3d & ownerCoordinates, const ElementObj *& child, Vector3d & childPCoords) const
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
    const Vector3d currentChildPCoords = subelement->compute_local_coords_from_owner_coordinates(this, ownerCoordinates);
    const double currentChildSqrDist = compute_parametric_square_distance(currentChildPCoords);
    if (currentChildSqrDist < minSqrDist)
    {
      minSqrDist = currentChildSqrDist;
      child = subelement;
      childPCoords = currentChildPCoords;
    }
  }
}

Vector3d
Mesh_Element::coordinates( const Vector3d & p_coords ) const
{ /* %TRACE% */  /* %TRACE% */

  const int npeCoords = my_nodes.size();
  std::vector<double> shapeFcn(npeCoords);
  my_master_elem.shape_fcn(1, p_coords.data(), shapeFcn.data());

  Vector3d coords(Vector3d::ZERO);
  for ( int n = 0; n < npeCoords; n++ )
    coords += shapeFcn[n] * my_nodes[n]->coordinates();

  return coords;
}

std::string Mesh_Element::visualize(const CDMesh & mesh) const
{
  std::ostringstream os;
  os << "Debugging mesh element " << entityId() << " with cutter:\n" << myCutter->visualize(mesh.stk_bulk()) << "\n";
  return os.str();
}

double
Mesh_Element::interface_crossing_position(const InterfaceID interface, const Segment3d & edge) const
{
  ThrowRequire(get_cutter());
  return get_cutter()->interface_crossing_position(interface, edge);
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
Mesh_Element::fill_face_interior_intersections(const NodeVec & faceNodes, const InterfaceID & interface1, const InterfaceID & interface2, std::vector<ElementIntersection> & faceIntersectionPoints) const
{
  ThrowRequire(get_cutter() && faceNodes.size() == 3);
  const std::array<Vector3d,3> faceNodeOwnerCoords = {{faceNodes[0]->owner_coords(this), faceNodes[1]->owner_coords(this), faceNodes[2]->owner_coords(this)}};
  const ElementIntersectionPointFilter intersectionPointFilter = build_element_intersection_filter(faceNodes);
  get_cutter()->fill_tetrahedron_face_interior_intersections(faceNodeOwnerCoords, interface1, interface2, intersectionPointFilter, faceIntersectionPoints);
}

int
Mesh_Element::interface_node_sign(const InterfaceID interface, const SubElementNode * node) const
{
  ThrowRequire(get_cutter());
  return get_cutter()->sign_at_position(interface, node->owner_coords(this));
}

double
Mesh_Element::interface_crossing_position(const InterfaceID interface, const SubElementNode * node1, const SubElementNode * node2) const
{
  ThrowRequire(get_cutter());
  Segment3d edge(node1->owner_coords(this), node2->owner_coords(this));
  return get_cutter()->interface_crossing_position(interface, edge);
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
  ThrowAssert(have_interface(interface));
  const auto iter = std::lower_bound(myCuttingInterfaces.begin(), myCuttingInterfaces.end(), interface);
  return std::distance(myCuttingInterfaces.begin(), iter);
}

bool
Mesh_Element::triangulate(const CDMesh & mesh, const InterfaceGeometry & interfaceGeometry)
{ /* %TRACE% */  /* %TRACE% */

  const auto & decomposedBlocksSelector = mesh.get_phase_support().get_all_decomposed_blocks_selector();
  const bool inDecomposedBlock = decomposedBlocksSelector(mesh.stk_bulk().bucket(entity()));
  if (inDecomposedBlock)
  {
    create_cutter(mesh, interfaceGeometry);

    if (!myCuttingInterfaces.empty())
      my_have_interface = true;

    if (my_phase.empty())
      my_phase = interfaceGeometry.get_starting_phase(myCutter.get());
    else
      my_have_interface = false; // uncut element with phase already set
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
  [](const std::vector<stk::mesh::Entity> & intersectionPointNodes, const std::vector<int> & intersectionPointSortedDomains)
  {
    return true;
  };
  return filter;
}

Vector3d Mesh_Element::get_intersection_point_parametric_coordinates(const IntersectionPoint & intersectionPoint) const
{
  Vector3d pCoords(Vector3d::ZERO);
  const auto & intersectionPointNodes = intersectionPoint.get_nodes();
  const auto & intersectionPointWeights = intersectionPoint.get_weights();
  const size_t numIntersectionPointNodes = intersectionPointNodes.size();
  for (size_t iNode=0; iNode<numIntersectionPointNodes; ++iNode)
  {
    const SubElementNode * node = get_node_matching_entity(my_nodes, intersectionPointNodes[iNode]);
    ThrowAssert(node);
    pCoords += intersectionPointWeights[iNode] * node->owner_coords(this);
  }
  return pCoords;
}

std::vector<int>
Mesh_Element::get_interface_signs_based_on_crossings(const NodeVec & nodes) const
{
  std::vector<Vector3d> nodeCoords;
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
Mesh_Element::cut_interior_intersection_points(CDMesh & mesh)
{ /* %TRACE% */  /* %TRACE% */
  if (!have_interface())
    return;

  std::vector<IntersectionPoint> intersectionPoints;
  stk::topology elementTopology = mesh.stk_bulk().bucket(entity()).topology();
  const MasterElement & masterElement = MasterElementDeterminer::getMasterElement(elementTopology);
  const std::vector<stk::mesh::Entity> elementNodes(mesh.stk_bulk().begin_nodes(entity()), mesh.stk_bulk().end_nodes(entity()));
  append_intersection_points_from_element_interior(masterElement, elementNodes, *myCutter, keep_all_intersecion_points_filter(), intersectionPoints);

  for (auto && intersectionPoint : intersectionPoints)
  {
    const Vector3d pCoords = get_intersection_point_parametric_coordinates(intersectionPoint);

    const ElementObj * containingElem = nullptr;
    Vector3d containingElemPCoords;
    find_child_coordinates_at_owner_coordinates(pCoords, containingElem, containingElemPCoords);
    ThrowAssert(containingElem);

    ElementObj * elem = const_cast<ElementObj *>(containingElem);
    if (!elem->captures_intersection_point_domains(intersectionPoint.get_sorted_domains()))
      elem->cut_interior_intersection_point(mesh, containingElemPCoords, intersectionPoint.get_sorted_domains());
  }

  const std::vector<InterfaceID> interfaces = get_sorted_cutting_interfaces();
  for (size_t i1=0; i1<interfaces.size(); ++i1)
  {
    for (size_t i2=i1+1; i2<interfaces.size(); ++i2)
    {
      for (auto && subelem : my_subelements)
        subelem->cut_face_interior_intersection_points(mesh, interfaces[i1], interfaces[i2]);
    }
  }

  std::vector<const SubElement *> leafSubElements;
  get_subelements(leafSubElements);
  for (auto && leafSubElement : leafSubElements)
  {
    SubElement * subElem = const_cast<SubElement *>(leafSubElement);
    const std::vector<int> interfaceSigns = get_interface_signs_based_on_crossings(subElem->get_nodes());
    subElem->set_interface_signs(interfaceSigns);
  }
}

void
Mesh_Element::create_cutter(const CDMesh & mesh, const InterfaceGeometry & interfaceGeometry)
{ /* %TRACE% */  /* %TRACE% */
  const auto intersectingPlanesDiagonalPicker = get_diagonal_picker();
  myCutter = interfaceGeometry.build_element_cutter(mesh.stk_bulk(), entity(), intersectingPlanesDiagonalPicker);
  myCuttingInterfaces = myCutter->get_sorted_cutting_interfaces();
}

void
Mesh_Element::create_base_subelement()
{ /* %TRACE% */  /* %TRACE% */

  stk::topology base_topology = coord_topology();

  const unsigned num_sides = base_topology.num_sides();

  std::vector<int> parent_side_ids(num_sides);
  for (unsigned i=0; i<num_sides; ++i)
  {
    parent_side_ids[i] = i;
  }

  std::unique_ptr<SubElement> base_subelement;
  if (stk::topology::TETRAHEDRON_4 == base_topology)
  {
    base_subelement = std::make_unique<SubElement_Tet_4>( my_nodes, parent_side_ids, this);
  }
  else if (stk::topology::TETRAHEDRON_10 == base_topology)
  {
    // purposely use lower order base subelement
    std::vector<const SubElementNode *> sub_nodes(my_nodes.begin(), my_nodes.begin()+4);
    base_subelement = std::make_unique<SubElement_Tet_4>( sub_nodes, parent_side_ids, this);
  }
  else if (stk::topology::TRIANGLE_3_2D == base_topology)
  {
    base_subelement = std::make_unique<SubElement_Tri_3>( my_nodes, parent_side_ids, this);
  }
  else if (stk::topology::TRIANGLE_6_2D == base_topology)
  {
    // purposely use lower order base subelement
    std::vector<const SubElementNode *> sub_nodes(my_nodes.begin(), my_nodes.begin()+3);
    base_subelement = std::make_unique<SubElement_Tri_3>( sub_nodes, parent_side_ids, this);
  }
  ThrowErrorMsgIf(!base_subelement, "Elements with topology " << base_topology.name() << " not supported for CDFEM.");

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

int
ElementObj::evaluate_quad(const SubElementNode * n0, const SubElementNode * n1, const SubElementNode * n2, const SubElementNode * n3)
{ /* %TRACE% */  /* %TRACE% */
  return evaluate_quad(n0->coordinates(), n1->coordinates(), n2->coordinates(), n3->coordinates());
}

int
ElementObj::evaluate_quad(const Vector3d & x0, const Vector3d & x1, const Vector3d & x2, const Vector3d & x3)
{ /* %TRACE% */  /* %TRACE% */
  // given 4 angles
  // angle0 - angle subtended by (x3-x0) and (x1-x0)
  // angle1 - angle subtended by (x0-x1) and (x2-x1)
  // angle2 - angle subtended by (x1-x2) and (x3-x2)
  // angle3 - angle subtended by (x2-x3) and (x0-x3)
  // returns -1 if the max of angle 0 and angle 2 is significant larger than the max of angle 1 and angle 3
  // returns +1 if the opposite is true
  // returns 0 if neither is true (the max angles are basically the same) OR one of the sides is degenerate

  const Vector3d side0 = x1 - x0;
  const Vector3d side1 = x2 - x1;
  const Vector3d side2 = x3 - x2;
  const Vector3d side3 = x0 - x3;

  const double side_len0 = side0.length();
  const double side_len1 = side1.length();
  const double side_len2 = side2.length();
  const double side_len3 = side3.length();

  // here, meas0 = -cos(angle0) * side_len0*side_len1*side_len2*side_len3
  const double meas0 = Dot(side3,side0) * side_len1*side_len2;
  const double meas1 = Dot(side0,side1) * side_len2*side_len3;
  const double meas2 = Dot(side1,side2) * side_len3*side_len0;
  const double meas3 = Dot(side2,side3) * side_len0*side_len1;

  if ( krinolog.shouldPrint(LOG_DEBUG) )
  {
    krinolog << "Lengths: " << side_len0 << ", " << side_len1 << ", " << side_len2 << ", " << side_len3 << "\n";
    krinolog << "Angle measures: " << meas0 << ", " << meas1 << ", " << meas2 << ", " << meas3 << "\n";
    krinolog << "Angles: " << std::acos(-meas0/(side_len0*side_len1*side_len2*side_len3))*180.0/3.141592654 << ", "
             << std::acos(-meas1/(side_len0*side_len1*side_len2*side_len3))*180.0/3.141592654 << ", "
             << std::acos(-meas2/(side_len0*side_len1*side_len2*side_len3))*180.0/3.141592654 << ", "
             << std::acos(-meas3/(side_len0*side_len1*side_len2*side_len3))*180.0/3.141592654 << "\n";
  }

  const double meas02 = std::max(meas0,meas2);
  const double meas13 = std::max(meas1,meas3);

  const double tol = std::numeric_limits<float>::epsilon()*(side_len0*side_len1*side_len2*side_len3);

  if (meas02 > (meas13 + tol))
  {
    return -1;
  }
  else if (meas13 > (meas02 + tol))
  {
    return +1;
  }
  else
  {
    return 0;
  }
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

  ThrowAssert(pn0_id != 0 && pn1_id != 0 && pn2_id != 0 && pn3_id != 0);

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

    SubElementNodeAncestry ancestry02 = std::min(ancestry0, ancestry2);
    SubElementNodeAncestry ancestry13 = std::min(ancestry1, ancestry3);

    return (ancestry02 < ancestry13);
  }

  // Make sure that diagonal can be picked either way and still avoid Steiner point
  const int caseId = (face0 ? 1 : 0) + (face1 ? 2 : 0) + (face2 ? 4 : 0) + (face3 ? 8 : 0);
  ThrowRequire(caseId == 3 || caseId == 5 || caseId == 10 || caseId == 12);

  return (evaluate_quad(n0,n1,n2,n3) == -1);
}

} // namespace krino
