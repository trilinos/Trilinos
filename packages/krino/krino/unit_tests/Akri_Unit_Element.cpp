// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <stk_mesh/base/GetEntities.hpp>

#include <Akri_Unit_Single_Element_Fixtures.hpp>

#include <Akri_CDFEM_Parent_Edge.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_SubElement.hpp>
#include <Akri_Element.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_LevelSetInterfaceGeometry.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

template <int TOPO>
class Mesh_Element_Fixture : public ::testing::Test
{
public:
  Mesh_Element_Fixture() :
    elem_fixture(static_cast<stk::topology::topology_t>(TOPO)),
    krino_mesh(elem_fixture.stk_fixture.bulk_data())
  {
    elem_fixture.generate_mesh();
    check_entity_counts();
    Phase_Support::get(stk_meta()).force_cdfem_use_case_for_minimal_unit_tests();
    Phase_Support::get(stk_meta()).set_one_levelset_per_phase(false);
    const NodeToCapturedDomainsMap nodesToCapturedDomains;
    interfaceGeometry = std::make_unique<LevelSetInterfaceGeometry>(AuxMetaData::get(stk_meta()).active_part(), CDFEM_Support::get(stk_meta()), Phase_Support::get(stk_meta()));
    interfaceGeometry->prepare_to_decompose_elements(krino_mesh.stk_bulk(), nodesToCapturedDomains);
  }
  virtual ~Mesh_Element_Fixture() {};
  void check_entity_counts()
  {
    std::vector<stk::mesh::Entity> entities;
    stk::mesh::get_entities(stk_bulk(), stk::topology::ELEMENT_RANK, entities);
    ASSERT_EQ(1u, entities.size());
    stk::mesh::get_entities(stk_bulk(), stk::topology::NODE_RANK, entities);
    ASSERT_EQ(elem_fixture.my_topology.num_nodes(), entities.size());
  }
  stk::mesh::BulkData & stk_bulk() { return elem_fixture.stk_fixture.bulk_data(); }
  stk::mesh::MetaData & stk_meta() { return elem_fixture.stk_fixture.meta_data(); }
  stk::mesh::Entity elem() { return elem_fixture.my_elem; }

  unsigned node_ordinal(const stk::mesh::Entity node)
  {
    const stk::mesh::Entity * const elem_nodes = stk_bulk().begin_nodes(elem());
    const unsigned num_nodes = stk_bulk().num_nodes(elem());
    const stk::mesh::Entity * found = std::find(elem_nodes, elem_nodes+num_nodes, node);
    return found - elem_nodes;
  }

  void generate_nonconformal_elements() {
    krino_mesh.clear();
    krino_mesh.generate_nonconformal_elements();
  }
  void triangulate()
  {
    Mesh_Element & meshElem = get_mesh_element();
    meshElem.create_cutter(krino_mesh, *interfaceGeometry); // update cutter with for edges that now have crossings found
    get_mesh_element().create_cutter(krino_mesh, *interfaceGeometry);
    krino_mesh.triangulate(*interfaceGeometry);
  }

  LevelSetElementCutter & get_cutter()
  {
    Mesh_Element & meshElem = get_mesh_element();
    meshElem.create_cutter(krino_mesh, *interfaceGeometry);
    LevelSetElementCutter * cutter = dynamic_cast<LevelSetElementCutter *>(meshElem.get_cutter());
    STK_ThrowRequire(cutter);
    return *cutter;
  }

  void find_edge_crossings(const std::vector<double> & node_LS_values)
  {
    std::vector<std::vector<double> > nodesIsovar(2);
    LevelSetElementCutter & cutter = get_cutter();
    for (unsigned iEdge=0; iEdge<get_mesh_element().topology().num_edges(); ++iEdge)
    {
      const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(elem_fixture.my_topology, iEdge);
      nodesIsovar[0].assign(1, node_LS_values[edgeNodeOrdinals[0]]);
      nodesIsovar[1].assign(1, node_LS_values[edgeNodeOrdinals[1]]);
      cutter.update_edge_crossings(iEdge, nodesIsovar);
    }
  }

  void
  find_edge_crossings(const std::vector<std::vector<double> > & node_LS_values)
  {
    std::vector<std::vector<double> > nodesIsovar(2);
    LevelSetElementCutter & cutter = get_cutter();
    for (unsigned iEdge=0; iEdge<get_mesh_element().topology().num_edges(); ++iEdge)
    {
      const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(elem_fixture.my_topology, iEdge);
      nodesIsovar[0] = node_LS_values[edgeNodeOrdinals[0]];
      nodesIsovar[1] = node_LS_values[edgeNodeOrdinals[1]];
      cutter.update_edge_crossings(iEdge, nodesIsovar);
    }
  }

  void generate_mesh_element_and_cutter(const std::vector<double> & node_LS_values)
  {
    generate_nonconformal_elements();
    find_edge_crossings(node_LS_values);
    triangulate();
  }

  void generate_mesh_element_and_cutter(const std::vector<std::vector<double> > & node_LS_values)
  {
    generate_nonconformal_elements();
    find_edge_crossings(node_LS_values);
    triangulate();
  }

  Mesh_Element & get_mesh_element() {
    return *krino_mesh.elements.front();
  }

  std::array<const SubElementNode *,2> get_edge_nodes(const unsigned edgeOrdinal)
  {
    const auto & nodes = get_mesh_element().get_nodes();
    const unsigned * edge_node_ordinals = get_edge_node_ordinals(elem_fixture.my_topology, edgeOrdinal);
    return {{nodes[edge_node_ordinals[0]], nodes[edge_node_ordinals[1]]}};
  }

  std::array<stk::math::Vector3d,2> get_edge_segment(const unsigned edge_ord)
  {
    const auto & edgeNodes = get_edge_nodes(edge_ord);
    return std::array<stk::math::Vector3d,2>{get_mesh_element().get_node_parametric_coords(edgeNodes[0]),
                                  get_mesh_element().get_node_parametric_coords(edgeNodes[1])};
  }

  SingleElementFixture elem_fixture;
  CDMesh krino_mesh;
  std::unique_ptr<LevelSetInterfaceGeometry> interfaceGeometry;
};

typedef Mesh_Element_Fixture<stk::topology::TRI_3_2D> Mesh_Element_Tri3;
TEST_F(Mesh_Element_Tri3, generate)
{
  generate_nonconformal_elements();
  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_EQ(elem(), mesh_elem.entity());
  const NodeVec & mesh_nodes = mesh_elem.get_nodes();
  ASSERT_EQ(3u, mesh_nodes.size());
  const stk::mesh::Entity * const elem_nodes = stk_bulk().begin_nodes(elem());
  for(unsigned i=0; i < mesh_nodes.size(); ++i)
  {
    EXPECT_EQ(mesh_nodes[i]->entity(), elem_nodes[i]);
  }
}

typedef Mesh_Element_Tri3 Mesh_Element_Tri3_One_LS;
TEST_F(Mesh_Element_Tri3_One_LS, All_Nodes_Positive)
{
  const InterfaceID iface(0,0);
  krino_mesh.add_interface_id(iface);

  const std::vector<double> node_LS_values(3, 1.);

  generate_mesh_element_and_cutter(node_LS_values);
  krino_mesh.determine_node_signs(iface);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_FALSE(mesh_elem.have_interface());

  for(unsigned i=0; i < mesh_elem.topology().num_edges(); ++i)
  {
    EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(InterfaceID(0,0), get_edge_segment(i)));
  }
}

TEST_F(Mesh_Element_Tri3_One_LS, All_Nodes_Negative)
{
  const InterfaceID iface(0,0);
  krino_mesh.add_interface_id(iface);

  const std::vector<double> node_LS_values(3, -1.);

  generate_mesh_element_and_cutter(node_LS_values);
  krino_mesh.determine_node_signs(iface);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_FALSE(mesh_elem.have_interface());

  for(unsigned i=0; i < mesh_elem.topology().num_edges(); ++i)
  {
    EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(InterfaceID(0,0), get_edge_segment(i)));
  }
}

TEST_F(Mesh_Element_Tri3_One_LS, Node0_Pos_Node1_2_Neg)
{
  const InterfaceID iface(0,0);
  krino_mesh.add_interface_id(iface);

  std::vector<double> node_LS_values(3);
  node_LS_values[0] = 1.;
  node_LS_values[1] = -0.5;
  node_LS_values[2] = -0.5;

  generate_mesh_element_and_cutter(node_LS_values);
  krino_mesh.determine_node_signs(iface);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_TRUE(mesh_elem.have_interface());
  const NodeVec & mesh_nodes = mesh_elem.get_nodes();
  EXPECT_EQ(3u, mesh_nodes.size());
  EXPECT_EQ(+1, mesh_nodes[0]->get_node_sign());
  EXPECT_EQ(-1, mesh_nodes[1]->get_node_sign());
  EXPECT_EQ(-1, mesh_nodes[2]->get_node_sign());

  // NOTE: Edge 0 is from nodes 0->1, Edge 1 is 1->2, Edge 2 is 2->0
  EXPECT_DOUBLE_EQ(2./3., mesh_elem.interface_crossing_position(iface, get_edge_segment(0)));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface, get_edge_segment(1)));
  EXPECT_DOUBLE_EQ(1./3., mesh_elem.interface_crossing_position(iface, get_edge_segment(2)));
}

TEST_F(Mesh_Element_Tri3_One_LS, Node0_Snapped_Node1_2_Pos)
{
  const InterfaceID iface(0,0);
  krino_mesh.add_interface_id(iface);

  std::vector<double> node_LS_values(3);
  node_LS_values[0] = 0;
  node_LS_values[1] = 1.;
  node_LS_values[2] = 1.;

  generate_mesh_element_and_cutter(node_LS_values);
  krino_mesh.determine_node_signs(iface);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_FALSE(mesh_elem.have_interface());

  // NOTE: Edge 0 is from nodes 0->1, Edge 1 is 1->2, Edge 2 is 2->0
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface, get_edge_segment(0)));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface, get_edge_segment(1)));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface, get_edge_segment(2)));
}

TEST_F(Mesh_Element_Tri3_One_LS, Node0_Snapped_Node1_2_Neg)
{
  const InterfaceID iface(0,0);
  krino_mesh.add_interface_id(iface);

  std::vector<double> node_LS_values(3);
  node_LS_values[0] = 0;
  node_LS_values[1] = -1.;
  node_LS_values[2] = -1.;

  generate_mesh_element_and_cutter(node_LS_values);
  krino_mesh.determine_node_signs(iface);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_TRUE(mesh_elem.have_interface());
  const NodeVec & mesh_nodes = mesh_elem.get_nodes();
  EXPECT_EQ(3u, mesh_nodes.size());
  EXPECT_EQ( 0, mesh_nodes[0]->get_node_sign());
  EXPECT_EQ(-1, mesh_nodes[1]->get_node_sign());
  EXPECT_EQ(-1, mesh_nodes[2]->get_node_sign());

  // NOTE: Edge 0 is from nodes 0->1, Edge 1 is 1->2, Edge 2 is 2->0
  // NOTE: Because we snap-away from nodes, crossing will be epsilon away from nodes.
  EXPECT_NEAR(0., mesh_elem.interface_crossing_position(iface, get_edge_segment(0)), std::numeric_limits<float>::epsilon());
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface, get_edge_segment(1)));
  EXPECT_NEAR(1., mesh_elem.interface_crossing_position(iface, get_edge_segment(2)), std::numeric_limits<float>::epsilon());
}

TEST_F(Mesh_Element_Tri3_One_LS, Node0_Snapped_Node1_Pos_Node_2_Neg)
{
  const InterfaceID iface(0,0);
  krino_mesh.add_interface_id(iface);

  std::vector<double> node_LS_values(3);
  node_LS_values[0] = 0;
  // Crossing position 0.608351703529745 for edge 1
  // gives a test case that fails to intersect exactly with node 0
  node_LS_values[1] = 1.;
  node_LS_values[2] = -(1/0.60835170352974499152765019971412 - 1.);

  generate_mesh_element_and_cutter(node_LS_values);
  krino_mesh.determine_node_signs(iface);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_TRUE(mesh_elem.have_interface());
  const NodeVec & mesh_nodes = mesh_elem.get_nodes();
  EXPECT_EQ(3u, mesh_nodes.size());
  EXPECT_EQ( 0, mesh_nodes[0]->get_node_sign());
  EXPECT_EQ(+1, mesh_nodes[1]->get_node_sign());
  EXPECT_EQ(-1, mesh_nodes[2]->get_node_sign());

  // NOTE: Edge 0 is from nodes 0->1, Edge 1 is 1->2, Edge 2 is 2->0
  // NOTE: Because we snap-away from nodes, crossing will be epsilon away from nodes.
  EXPECT_DOUBLE_EQ(node_LS_values[1]/(node_LS_values[1] - node_LS_values[2]), mesh_elem.interface_crossing_position(iface, get_edge_segment(1)));
  EXPECT_NEAR(1., mesh_elem.interface_crossing_position(iface, get_edge_segment(2)), std::numeric_limits<float>::epsilon());
}

TEST_F(Mesh_Element_Tri3_One_LS, Node0_1_Snapped_Node2_Pos)
{
  const InterfaceID iface(0,0);
  krino_mesh.add_interface_id(iface);

  std::vector<double> node_LS_values(3);
  node_LS_values[0] = 0;
  node_LS_values[1] = 0;
  node_LS_values[2] = 1.;

  generate_mesh_element_and_cutter(node_LS_values);
  krino_mesh.determine_node_signs(iface);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_FALSE(mesh_elem.have_interface());

  // NOTE: Edge 0 is from nodes 0->1, Edge 1 is 1->2, Edge 2 is 2->0
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface, get_edge_segment(0)));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface, get_edge_segment(1)));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface, get_edge_segment(2)));
}

TEST_F(Mesh_Element_Tri3_One_LS, Node0_1_Snapped_Node2_Neg)
{
  const InterfaceID iface(0,0);
  krino_mesh.add_interface_id(iface);

  std::vector<double> node_LS_values(3);
  node_LS_values[0] = 0;
  node_LS_values[1] = 0;
  node_LS_values[2] = -1.;

  generate_mesh_element_and_cutter(node_LS_values);
  krino_mesh.determine_node_signs(iface);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_TRUE(mesh_elem.have_interface());
  const NodeVec & mesh_nodes = mesh_elem.get_nodes();
  EXPECT_EQ(3u, mesh_nodes.size());
  EXPECT_EQ( 0, mesh_nodes[0]->get_node_sign());
  EXPECT_EQ( 0, mesh_nodes[1]->get_node_sign());
  EXPECT_EQ(-1, mesh_nodes[2]->get_node_sign());

  // NOTE: Edge 0 is from nodes 0->1, Edge 1 is 1->2, Edge 2 is 2->0
  // TODO: What should intersections be for snapped edge?
  // I don't think we should ever be looking for the interface crossing position of it
  // NOTE: Because we snap-away from nodes, crossing will be epsilon away from nodes.
  EXPECT_NEAR(0., mesh_elem.interface_crossing_position(iface, get_edge_segment(1)), std::numeric_limits<float>::epsilon());
  EXPECT_NEAR(1., mesh_elem.interface_crossing_position(iface, get_edge_segment(2)), std::numeric_limits<float>::epsilon());
}

typedef Mesh_Element_Tri3 Mesh_Element_Tri3_Three_LS;
TEST_F(Mesh_Element_Tri3_Three_LS, All_Phase2_Unsnapped)
{
  const LS_Field ls1("LS1", Surface_Identifier(1));
  const LS_Field ls2("LS2", Surface_Identifier(2));
  const LS_Field ls3("LS3", Surface_Identifier(3));
  interfaceGeometry->set_ls_fields({ls1, ls2, ls3});

  const InterfaceID iface01(0,1);
  const InterfaceID iface02(0,2);
  const InterfaceID iface12(1,2);
  krino_mesh.add_interface_id(iface01);
  krino_mesh.add_interface_id(iface02);
  krino_mesh.add_interface_id(iface12);
  Phase_Support::get(stk_meta()).set_one_levelset_per_phase(true);

  std::vector<std::vector<double> > node_LS_values(3);
  node_LS_values[0].resize(3);
  node_LS_values[1].resize(3);
  node_LS_values[2].resize(3);
  // LS 0
  node_LS_values[0][0] = 1;
  node_LS_values[1][0] = 1;
  node_LS_values[2][0] = 1;
  // LS 1
  node_LS_values[0][1] = 1;
  node_LS_values[1][1] = 1;
  node_LS_values[2][1] = 1;
  // LS 2
  node_LS_values[0][2] = -1;
  node_LS_values[1][2] = -1;
  node_LS_values[2][2] = -1;

  generate_mesh_element_and_cutter(node_LS_values);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_FALSE(mesh_elem.have_interface());
}

TEST_F(Mesh_Element_Tri3_Three_LS, One_Interface)
{
  const LS_Field ls1("LS1", Surface_Identifier(1));
  const LS_Field ls2("LS2", Surface_Identifier(2));
  const LS_Field ls3("LS3", Surface_Identifier(3));
  interfaceGeometry->set_ls_fields({ls1, ls2, ls3});
  const InterfaceID iface01(0,1);
  const InterfaceID iface02(0,2);
  const InterfaceID iface12(1,2);
  krino_mesh.add_interface_id(iface01);
  krino_mesh.add_interface_id(iface02);
  krino_mesh.add_interface_id(iface12);
  Phase_Support::get(stk_meta()).set_one_levelset_per_phase(true);

  std::vector<std::vector<double> > node_LS_values(3);
  node_LS_values[0].resize(3);
  node_LS_values[1].resize(3);
  node_LS_values[2].resize(3);
  // LS 0
  node_LS_values[0][0] = 0.038335;
  node_LS_values[1][0] = 0.014437;
  node_LS_values[2][0] = -0.01288;
  // LS 1
  node_LS_values[0][1] = 0.037250;
  node_LS_values[1][1] = 0.070753;
  node_LS_values[2][1] = 0.021996;
  // LS 2
  node_LS_values[0][2] = 0.01;
  node_LS_values[1][2] = 0.01;
  node_LS_values[2][2] = 0.01;

  generate_mesh_element_and_cutter(node_LS_values);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_TRUE(mesh_elem.have_interface());
  const NodeVec & mesh_nodes = mesh_elem.get_nodes();
  EXPECT_EQ(3u, mesh_nodes.size());

  const stk::math::Vector3d node0_coords = mesh_nodes[0]->owner_coords(&mesh_elem);
  const stk::math::Vector3d node1_coords = mesh_nodes[1]->owner_coords(&mesh_elem);
  const stk::math::Vector3d node2_coords = mesh_nodes[2]->owner_coords(&mesh_elem);
  const std::array<stk::math::Vector3d,2> edge0{node0_coords, node1_coords};
  const std::array<stk::math::Vector3d,2> edge1{node1_coords, node2_coords};
  const std::array<stk::math::Vector3d,2> edge2{node2_coords, node0_coords};

  krino_mesh.determine_node_signs(iface01);

  if (mesh_elem.have_interface(iface01))
  {
    EXPECT_EQ(+1, mesh_nodes[0]->get_node_sign());
    EXPECT_EQ(+1, mesh_nodes[1]->get_node_sign());
    EXPECT_EQ(+1, mesh_nodes[2]->get_node_sign());
  }

  krino_mesh.determine_node_signs(iface02);

  EXPECT_TRUE(mesh_elem.have_interface(iface02));
  EXPECT_EQ(+1, mesh_nodes[0]->get_node_sign());
  EXPECT_EQ(+1, mesh_nodes[1]->get_node_sign());
  EXPECT_EQ(-1, mesh_nodes[2]->get_node_sign());

  krino_mesh.determine_node_signs(iface12);

  if (mesh_elem.have_interface(iface12))
  {
    EXPECT_EQ(+1, mesh_nodes[0]->get_node_sign());
    EXPECT_EQ(+1, mesh_nodes[1]->get_node_sign());
    EXPECT_EQ(+1, mesh_nodes[2]->get_node_sign());
  }

  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface01, edge0));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface01, edge1));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface01, edge2));

  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface12, edge0));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface12, edge1));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface12, edge2));

  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface02, edge0));
}

TEST_F(Mesh_Element_Tri3_Three_LS, Handle_Hanging_Child)
{
  const LS_Field ls1("LS1", Surface_Identifier(1));
  const LS_Field ls2("LS2", Surface_Identifier(2));
  const LS_Field ls3("LS3", Surface_Identifier(3));
  interfaceGeometry->set_ls_fields({ls1, ls2, ls3});
  const InterfaceID iface01(0,1);
  const InterfaceID iface02(0,2);
  const InterfaceID iface12(1,2);
  krino_mesh.add_interface_id(iface01);
  krino_mesh.add_interface_id(iface02);
  krino_mesh.add_interface_id(iface12);
  Phase_Support::get(stk_meta()).set_one_levelset_per_phase(true);

  std::vector<std::vector<double> > node_LS_values(3);
  node_LS_values[0].resize(3);
  node_LS_values[1].resize(3);
  node_LS_values[2].resize(3);
  // LS 0
  node_LS_values[0][0] = 0.04;
  node_LS_values[1][0] = 0.015;
  node_LS_values[2][0] = 0.025;
  // LS 1
  node_LS_values[0][1] = 0.;
  node_LS_values[1][1] = 0.02;
  node_LS_values[2][1] = 0.02;
  // LS 2
  node_LS_values[0][2] = 0.01;
  node_LS_values[1][2] = 0.01;
  node_LS_values[2][2] = 0.01;

  generate_mesh_element_and_cutter(node_LS_values);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_TRUE(mesh_elem.have_interface());
  const NodeVec & mesh_nodes = mesh_elem.get_nodes();
  EXPECT_EQ(3u, mesh_nodes.size());

  const stk::math::Vector3d node0_coords = mesh_nodes[0]->owner_coords(&mesh_elem);
  const stk::math::Vector3d node1_coords = mesh_nodes[1]->owner_coords(&mesh_elem);
  const stk::math::Vector3d node2_coords = mesh_nodes[2]->owner_coords(&mesh_elem);
  const std::array<stk::math::Vector3d,2> edge0{node0_coords, node1_coords};
  const std::array<stk::math::Vector3d,2> edge1{node1_coords, node2_coords};
  const std::array<stk::math::Vector3d,2> edge2{node2_coords, node0_coords};

  krino_mesh.determine_node_signs(iface01);

  if (mesh_elem.have_interface(iface01))
  {
    EXPECT_EQ(+1, mesh_nodes[0]->get_node_sign());
    EXPECT_EQ(+1, mesh_nodes[1]->get_node_sign());
    EXPECT_EQ(+1, mesh_nodes[2]->get_node_sign());
  }

  krino_mesh.determine_node_signs(iface02);

  if (mesh_elem.have_interface(iface02))
  {
    EXPECT_EQ(+1, mesh_nodes[0]->get_node_sign());
    EXPECT_EQ(+1, mesh_nodes[1]->get_node_sign());
    EXPECT_EQ(+1, mesh_nodes[2]->get_node_sign());
  };

  krino_mesh.determine_node_signs(iface12);

  EXPECT_TRUE(mesh_elem.have_interface(iface12));
  EXPECT_EQ(-1, mesh_nodes[0]->get_node_sign());
  EXPECT_EQ(+1, mesh_nodes[1]->get_node_sign());
  EXPECT_EQ(+1, mesh_nodes[2]->get_node_sign());

  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface01, edge0));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface01, edge1));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface01, edge2));

  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface02, edge0));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface02, edge1));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface02, edge2));

  EXPECT_DOUBLE_EQ(0.5, mesh_elem.interface_crossing_position(iface12, edge0));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface12, edge1));
  EXPECT_DOUBLE_EQ(0.5, mesh_elem.interface_crossing_position(iface12, edge2));
}

TEST_F(Mesh_Element_Tri3_Three_LS, Zero_Crossings_For_Phases_Present_Bug)
{
  const LS_Field ls1("LS1", Surface_Identifier(1));
  const LS_Field ls2("LS2", Surface_Identifier(2));
  const LS_Field ls3("LS3", Surface_Identifier(3));
  interfaceGeometry->set_ls_fields({ls1, ls2, ls3});
  const InterfaceID iface01(0,1);
  const InterfaceID iface02(0,2);
  const InterfaceID iface12(1,2);
  krino_mesh.add_interface_id(iface01);
  krino_mesh.add_interface_id(iface02);
  krino_mesh.add_interface_id(iface12);
  Phase_Support::get(stk_meta()).set_one_levelset_per_phase(true);

  std::vector<std::vector<double> > node_LS_values(4);
  node_LS_values[0].resize(3);
  node_LS_values[1].resize(3);
  node_LS_values[2].resize(3);
  node_LS_values[3].resize(3);
  // Node 0
  node_LS_values[0][0] = -0.00765969;
  node_LS_values[0][1] = 0.532721;
  node_LS_values[0][2] = 0.;
  // Node 1
  node_LS_values[1][0] = -0.000100754;
  node_LS_values[1][1] = 0.;
  node_LS_values[1][2] = 0.;
  // Node 2
  node_LS_values[2][0] = -0.00666939;
  node_LS_values[2][1] = 0.337202;
  node_LS_values[2][2] = 0.;
  // Node 3 (Mid-node for edge 1)
  node_LS_values[3][0] = -2.76e-6;
  node_LS_values[3][1] = -0.00614775;
  node_LS_values[3][2] = 0.;

  generate_mesh_element_and_cutter(node_LS_values);

  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_FALSE(mesh_elem.have_interface());
  const NodeVec & mesh_nodes = mesh_elem.get_nodes();
  EXPECT_EQ(3u, mesh_nodes.size());

  const stk::math::Vector3d node0_coords = mesh_nodes[0]->owner_coords(&mesh_elem);
  const stk::math::Vector3d node1_coords = mesh_nodes[1]->owner_coords(&mesh_elem);
  const stk::math::Vector3d node2_coords = mesh_nodes[2]->owner_coords(&mesh_elem);
  const std::array<stk::math::Vector3d,2> edge0{node0_coords, node1_coords};
  const std::array<stk::math::Vector3d,2> edge1{node1_coords, node2_coords};
  const std::array<stk::math::Vector3d,2> edge2{node2_coords, node0_coords};

  krino_mesh.determine_node_signs(iface01);

  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface01, edge0));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface01, edge1));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface01, edge2));

  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface02, edge0));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface02, edge1));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface02, edge2));

  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface12, edge0));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface12, edge1));
  EXPECT_ANY_THROW(mesh_elem.interface_crossing_position(iface12, edge2));
}

class Mesh_Element_Tet4 : public Mesh_Element_Fixture<stk::topology::TET_4>
{
protected:
  void test_node_signs(const std::vector<int> & gold)
  {
    const NodeVec & mesh_nodes = get_mesh_element().get_nodes();
    ASSERT_EQ(4u, gold.size());
    ASSERT_EQ(4u, mesh_nodes.size());
    for (int n=0; n<4; ++n)
    {
      EXPECT_EQ(gold[n], mesh_nodes[n]->get_node_sign());
    }
  }

  void test_node_score_rank(const std::vector<unsigned> & goldNodeRankHiToLow)
  {
    const NodeVec & mesh_nodes = get_mesh_element().get_nodes();
    ASSERT_EQ(4u, goldNodeRankHiToLow.size());
    ASSERT_EQ(4u, mesh_nodes.size());
    for (int n1=0; n1<4; ++n1)
    {
      for (int n2=n1+1; n2<4; ++n2)
      {
        if (goldNodeRankHiToLow[n1] > goldNodeRankHiToLow[n2])
          EXPECT_LT(mesh_nodes[n1]->get_node_score(), mesh_nodes[n2]->get_node_score());
        else if (goldNodeRankHiToLow[n1] < goldNodeRankHiToLow[n2])
          EXPECT_GT(mesh_nodes[n1]->get_node_score(), mesh_nodes[n2]->get_node_score());
        else
          EXPECT_EQ(mesh_nodes[n1]->get_node_score(), mesh_nodes[n2]->get_node_score());
      }
    }
  }

  template<typename VEC>
  VEC permute(const VEC & input, const int permutation)
  {
    VEC output;
    output.resize(input.size());
    stk::topology topo = stk::topology::TETRAHEDRON_4;
    topo.permutation_nodes(input.data(), permutation, output.data());
    return output;
  }

  void set_ls_and_test_node_signs_and_scores(const std::vector<double> & nodeLsValues, const std::vector<int> & goldNodeSigns, const std::vector<unsigned> & goldNodeRankHiToLo, const int permutation)
  {
    const InterfaceID iface(0,0);
    krino_mesh.add_interface_id(iface);
    Phase_Support::get(stk_meta()).set_one_levelset_per_phase(false);

    const auto permutedNodeLsValues = permute(nodeLsValues, permutation);
    const auto permutedGoldNodeSigns = permute(goldNodeSigns, permutation);
    const auto permutedGoldNodeRankHiToLow = permute(goldNodeRankHiToLo, permutation);

    generate_mesh_element_and_cutter(permutedNodeLsValues);
    krino_mesh.determine_node_signs(iface);

    EXPECT_TRUE(get_mesh_element().have_interface());

    test_node_signs(permutedGoldNodeSigns);

    krino_mesh.decompose_edges(iface);
    krino_mesh.determine_node_scores(iface);

    test_node_score_rank(permutedGoldNodeRankHiToLow);
  }

  void set_ls_and_test_node_signs_and_scores_for_all_permutations(const std::vector<double> & nodeLsValues, const std::vector<int> & goldNodeSigns, const std::vector<unsigned> & goldNodeRankHiToLo)
  {
    for (int iperm=0; iperm<12; ++iperm)
    {
      set_ls_and_test_node_signs_and_scores(nodeLsValues, goldNodeSigns, goldNodeRankHiToLo, iperm);
    }
  }
};

TEST_F(Mesh_Element_Tet4, generate)
{
  generate_nonconformal_elements();
  Mesh_Element & mesh_elem = get_mesh_element();
  EXPECT_EQ(elem(), mesh_elem.entity());
  const NodeVec & mesh_nodes = mesh_elem.get_nodes();
  ASSERT_EQ(4u, mesh_nodes.size());
  const stk::mesh::Entity * const elem_nodes = stk_bulk().begin_nodes(elem());
  for(unsigned i=0; i < mesh_nodes.size(); ++i)
  {
    EXPECT_EQ(mesh_nodes[i]->entity(), elem_nodes[i]);
  }
}

TEST_F(Mesh_Element_Tet4, OneInterfaceCheckNodeScore_allPermutations)
{
  set_ls_and_test_node_signs_and_scores_for_all_permutations({1., -0.2, -0.5, -0.8}, {+1, -1, -1, -1}, {0,3,2,1});
}

TEST_F(Mesh_Element_Tet4, OneInterfaceCheckNodeScore_onlyQuadNodesGetScore)
{
  set_ls_and_test_node_signs_and_scores_for_all_permutations({0.2, -0.2, -0.5, -0.8}, {+1, -1, -1, -1}, {0,3,2,1});
}

TEST_F(Mesh_Element_Tet4, OneInterfaceCheckNodeScore_ScoreBasedOnAngleNotPosition)
{
  krino_mesh.get_cdfem_support().set_simplex_generation_method(CUT_QUADS_BY_LARGEST_ANGLE);
  set_ls_and_test_node_signs_and_scores({-1.0, -1.01, -1.02, 1.0}, {-1, -1, -1, +1}, {1,3,2,0}, 0);
  set_ls_and_test_node_signs_and_scores({-1.0, -1.01, -1.02, 1.0}, {-1, -1, -1, +1}, {1,3,2,0}, 3);
  set_ls_and_test_node_signs_and_scores({-1.0, -1.01, -1.02, 1.0}, {-1, -1, -1, +1}, {1,3,2,0}, 6);
}

}
