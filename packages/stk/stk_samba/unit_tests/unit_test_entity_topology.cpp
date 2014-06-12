#include <gtest/gtest.h>

#include <samba/entity_topology.hpp>

#include <iostream>
#include <sstream>
#include <string>

TEST(samba, entity_topology_print)
{
  using namespace samba;

  std::ostringstream out;

  for (entity_topology t = {0}; is_valid_topology(t); ++t)
  {
    print_topology_detail(out, t);
  }

  std::string expected =
"{entity_topology:node}:\n\
           dimension: 0\n\
        num_vertices: 0\n\
           num_nodes: 0\n\
           num_edges: 0\n\
           num_faces: 0\n\
           num_sides: 0\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:particle}:\n\
           dimension: 1\n\
        num_vertices: 1\n\
           num_nodes: 1\n\
           num_edges: 0\n\
           num_faces: 0\n\
           num_sides: 0\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:line_2}:\n\
           dimension: 1\n\
        num_vertices: 2\n\
           num_nodes: 2\n\
           num_edges: 0\n\
           num_faces: 0\n\
           num_sides: 0\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:line_3}:\n\
           dimension: 1\n\
        num_vertices: 2\n\
           num_nodes: 3\n\
           num_edges: 0\n\
           num_faces: 0\n\
           num_sides: 0\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:beam_2}:\n\
           dimension: 2\n\
        num_vertices: 2\n\
           num_nodes: 2\n\
           num_edges: 1\n\
           num_faces: 0\n\
           num_sides: 1\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:beam_3}:\n\
           dimension: 2\n\
        num_vertices: 2\n\
           num_nodes: 3\n\
           num_edges: 1\n\
           num_faces: 0\n\
           num_sides: 1\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:shell_line_2}:\n\
           dimension: 2\n\
        num_vertices: 2\n\
           num_nodes: 2\n\
           num_edges: 2\n\
           num_faces: 0\n\
           num_sides: 2\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:shell_line_3}:\n\
           dimension: 2\n\
        num_vertices: 2\n\
           num_nodes: 3\n\
           num_edges: 2\n\
           num_faces: 0\n\
           num_sides: 2\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:triangle_3}:\n\
           dimension: 2\n\
        num_vertices: 3\n\
           num_nodes: 3\n\
           num_edges: 3\n\
           num_faces: 0\n\
           num_sides: 3\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:triangle_4}:\n\
           dimension: 2\n\
        num_vertices: 3\n\
           num_nodes: 4\n\
           num_edges: 3\n\
           num_faces: 0\n\
           num_sides: 3\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:triangle_6}:\n\
           dimension: 2\n\
        num_vertices: 3\n\
           num_nodes: 6\n\
           num_edges: 3\n\
           num_faces: 0\n\
           num_sides: 3\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:quadrilateral_4}:\n\
           dimension: 2\n\
        num_vertices: 4\n\
           num_nodes: 4\n\
           num_edges: 4\n\
           num_faces: 0\n\
           num_sides: 4\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:quadrilateral_8}:\n\
           dimension: 2\n\
        num_vertices: 4\n\
           num_nodes: 8\n\
           num_edges: 4\n\
           num_faces: 0\n\
           num_sides: 4\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:quadrilateral_9}:\n\
           dimension: 2\n\
        num_vertices: 4\n\
           num_nodes: 9\n\
           num_edges: 4\n\
           num_faces: 0\n\
           num_sides: 4\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:pentagon_5}:\n\
           dimension: 2\n\
        num_vertices: 5\n\
           num_nodes: 5\n\
           num_edges: 5\n\
           num_faces: 0\n\
           num_sides: 5\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:hexagon_6}:\n\
           dimension: 2\n\
        num_vertices: 6\n\
           num_nodes: 6\n\
           num_edges: 6\n\
           num_faces: 0\n\
           num_sides: 6\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:shell_triangle_3}:\n\
           dimension: 3\n\
        num_vertices: 3\n\
           num_nodes: 3\n\
           num_edges: 3\n\
           num_faces: 0\n\
           num_sides: 3\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:shell_triangle_6}:\n\
           dimension: 3\n\
        num_vertices: 3\n\
           num_nodes: 6\n\
           num_edges: 3\n\
           num_faces: 0\n\
           num_sides: 3\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:shell_quadrilateral_4}:\n\
           dimension: 3\n\
        num_vertices: 4\n\
           num_nodes: 4\n\
           num_edges: 4\n\
           num_faces: 0\n\
           num_sides: 4\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:shell_quadrilateral_8}:\n\
           dimension: 3\n\
        num_vertices: 4\n\
           num_nodes: 8\n\
           num_edges: 4\n\
           num_faces: 0\n\
           num_sides: 4\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:shell_quadrilateral_9}:\n\
           dimension: 3\n\
        num_vertices: 4\n\
           num_nodes: 9\n\
           num_edges: 4\n\
           num_faces: 0\n\
           num_sides: 4\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:tetrahedron_4}:\n\
           dimension: 3\n\
        num_vertices: 4\n\
           num_nodes: 4\n\
           num_edges: 6\n\
           num_faces: 4\n\
           num_sides: 4\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:tetrahedron_8}:\n\
           dimension: 3\n\
        num_vertices: 4\n\
           num_nodes: 8\n\
           num_edges: 6\n\
           num_faces: 4\n\
           num_sides: 4\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:tetrahedron_10}:\n\
           dimension: 3\n\
        num_vertices: 4\n\
           num_nodes: 10\n\
           num_edges: 6\n\
           num_faces: 4\n\
           num_sides: 4\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:tetrahedron_11}:\n\
           dimension: 3\n\
        num_vertices: 4\n\
           num_nodes: 11\n\
           num_edges: 6\n\
           num_faces: 4\n\
           num_sides: 4\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:pyramid_5}:\n\
           dimension: 3\n\
        num_vertices: 5\n\
           num_nodes: 5\n\
           num_edges: 8\n\
           num_faces: 5\n\
           num_sides: 5\n\
   homogeneous_sides: 0\n\
\n\
{entity_topology:pyramid_13}:\n\
           dimension: 3\n\
        num_vertices: 5\n\
           num_nodes: 13\n\
           num_edges: 8\n\
           num_faces: 5\n\
           num_sides: 5\n\
   homogeneous_sides: 0\n\
\n\
{entity_topology:pyramid_14}:\n\
           dimension: 3\n\
        num_vertices: 5\n\
           num_nodes: 14\n\
           num_edges: 8\n\
           num_faces: 5\n\
           num_sides: 5\n\
   homogeneous_sides: 0\n\
\n\
{entity_topology:wedge_6}:\n\
           dimension: 3\n\
        num_vertices: 6\n\
           num_nodes: 6\n\
           num_edges: 9\n\
           num_faces: 5\n\
           num_sides: 5\n\
   homogeneous_sides: 0\n\
\n\
{entity_topology:wedge_15}:\n\
           dimension: 3\n\
        num_vertices: 6\n\
           num_nodes: 15\n\
           num_edges: 9\n\
           num_faces: 5\n\
           num_sides: 5\n\
   homogeneous_sides: 0\n\
\n\
{entity_topology:wedge_18}:\n\
           dimension: 3\n\
        num_vertices: 6\n\
           num_nodes: 18\n\
           num_edges: 9\n\
           num_faces: 5\n\
           num_sides: 5\n\
   homogeneous_sides: 0\n\
\n\
{entity_topology:hexahedron_8}:\n\
           dimension: 3\n\
        num_vertices: 8\n\
           num_nodes: 8\n\
           num_edges: 12\n\
           num_faces: 6\n\
           num_sides: 6\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:hexahedron_20}:\n\
           dimension: 3\n\
        num_vertices: 8\n\
           num_nodes: 20\n\
           num_edges: 12\n\
           num_faces: 6\n\
           num_sides: 6\n\
   homogeneous_sides: 1\n\
\n\
{entity_topology:hexahedron_27}:\n\
           dimension: 3\n\
        num_vertices: 8\n\
           num_nodes: 27\n\
           num_edges: 12\n\
           num_faces: 6\n\
           num_sides: 6\n\
   homogeneous_sides: 1\n\
\n";

  EXPECT_EQ(expected, out.str());
}

TEST(samba, entity_topology_side_nodes)
{
  using namespace samba;

  entity_topology hex8 = entity_topology::hex_8();
  connectivity_ordinal const* const null_rels = NULL;

  EXPECT_NE(null_rels, side_nodes(hex8, 1 /*valid side-id*/));
  EXPECT_EQ(null_rels, side_nodes(hex8, 666 /*bad side-id*/));

  EXPECT_NE(null_rels, edge_nodes(hex8, 1 /*valid edge-id*/));
  EXPECT_EQ(null_rels, edge_nodes(hex8, 666 /*bad edge-id*/));

  entity_topology hex8_side_topo = side_topology(hex8, 1 /*side id*/);
  EXPECT_EQ(entity_topology::quad_4(), hex8_side_topo);
  connectivity_ordinal const* const side_nodes = samba::side_nodes(hex8, 1 /*side id*/);
  connectivity_ordinal const* const face_nodes = samba::face_nodes(hex8, 1 /*side id*/);
  unsigned num_nodes_in_side = num_nodes(hex8_side_topo);
  EXPECT_TRUE(std::equal(side_nodes, side_nodes + num_nodes_in_side, face_nodes));

  // Test topology with non-homogenoues sides
  entity_topology mids14 = entity_topology::pyramid_14();
  entity_topology t_side_topo = side_topology(mids14, 1 /*side id*/);
  entity_topology q_side_topo = side_topology(mids14, 4 /*side id*/);
  EXPECT_EQ(entity_topology::triangle_6(), t_side_topo);
  EXPECT_EQ(entity_topology::quad_9(),     q_side_topo);
}
