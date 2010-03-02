/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <unit_tests/UnitTestGridMeshFixture.hpp>

#include <unit_tests/stk_utest_macros.hpp>
#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

using namespace stk;
using namespace stk::mesh;

/*
The following fixture creates the mesh below on proc 0
1-16 Quadrilateral<4>
17-41 Nodes
42-45 ShellLine<2>

17---18---19---20---21
|  1 |  2 |  3 || 4  |
22---23---24---25---26
|  5 |  6 |  7 ||  8 |
27---28---29---30---31
|  9 | 10 | 11 || 12 |
32---33---34---35---36
| 13 | 14 | 15 || 16 |
37---38---39---40---41
*/

GridMeshFixture::GridMeshFixture(ParallelMachine pm) : m_node_ids(),
                                                       m_quad_face_ids(),
                                                       m_shell_face_ids(),
                                                       m_meta_data(NULL),
                                                       m_bulk_data(NULL)
{
  m_meta_data = new MetaData( fem_entity_type_names() );
  m_quad_part = &(m_meta_data->declare_part("quad_part", Face));
  m_shell_part = &(m_meta_data->declare_part("shell_part", Face));
  set_cell_topology<shards::Quadrilateral<4> >(*m_quad_part);
  set_cell_topology<shards::ShellLine<2> >(*m_shell_part);
  // Add shells from nodes 20 to 40
  // shells are considered to be a Face

  m_meta_data->commit();

  const size_t num_entities_per_bucket = 100;
  m_bulk_data = new BulkData( *m_meta_data, pm, num_entities_per_bucket );

  PartVector no_parts;

  generate_grid();
  m_bulk_data->modification_end();
}

GridMeshFixture::~GridMeshFixture()
{
  // do not delete face part because the meta data owns it
  delete m_bulk_data;
  delete m_meta_data;
}

void GridMeshFixture::generate_grid()
{
  const unsigned num_nodes = 25;
  const unsigned num_quad_faces = 16;
  const unsigned num_shell_faces = 4;
  const unsigned p_rank = m_bulk_data->parallel_rank();
  std::vector<Entity*> all_entities;

  // we don't want anything on any of the processes expect rank 0
  if (p_rank != 0) {
    return;
  }

  m_node_ids.resize( num_nodes );
  m_quad_face_ids.resize( num_quad_faces );
  m_shell_face_ids.resize( num_shell_faces );

  // assign ids, quads, nodes, then shells
  // (we need this order to be this way in order for our connectivity setup to work)
  {
    unsigned curr_id = 1;
    for (unsigned  i = 0 ; i < num_quad_faces; ++i, ++curr_id) {
      m_quad_face_ids[i] = curr_id;
    }
    for (unsigned  i = 0 ; i < num_nodes; ++i, ++curr_id) {
      m_node_ids[i] = curr_id;
    }
    for (unsigned  i = 0 ; i < num_shell_faces; ++i, ++curr_id) {
      m_shell_face_ids[i] = curr_id;
    }
  }

  // declare entities such that entity_id - 1 is the index of the
  // entity in the all_entities vector
  std::vector<Entity*> quad_faces;
  std::vector<Entity*> nodes;
  std::vector<Entity*> shell_faces;
  {
    const PartVector no_parts;

    // declare faces
    PartVector face_parts;
    face_parts.push_back(m_quad_part);
    for (unsigned i = 0; i < num_quad_faces; ++i) {
      Entity& new_face = m_bulk_data->declare_entity(Face, m_quad_face_ids[i], face_parts);
      quad_faces.push_back(&new_face);
      all_entities.push_back(&new_face);
      switch (m_quad_face_ids[i]) {
      case 6:
      case 7:
      case 10:
      case 11:
      case 14:
      case 15:
        m_closure.push_back(&new_face);
      default:
        break;
      }
    }

    // declare nodes
    for (unsigned i = 0; i < num_nodes; ++i) {
      Entity& new_node = m_bulk_data->declare_entity(Node, m_node_ids[i], no_parts);
      nodes.push_back(&new_node);
      all_entities.push_back(&new_node);
      switch (m_node_ids[i]) {
      case 23:
      case 24:
      case 25:
      case 28:
      case 29:
      case 30:
      case 33:
      case 34:
      case 35:
      case 38:
      case 39:
      case 40:
        m_closure.push_back(&new_node);
      default:
        break;
      }
    }

    // declare shell faces
    PartVector shell_parts;
    shell_parts.push_back(m_shell_part);
    for (unsigned i = 0; i < num_shell_faces; ++i) {
      Entity& new_shell = m_bulk_data->declare_entity(Face, m_shell_face_ids[i], shell_parts);
      all_entities.push_back(&new_shell);
      shell_faces.push_back(&new_shell);
    }
  }

  // sort the closure
  std::sort(m_closure.begin(), m_closure.end(), stk::mesh::EntityLess());

  // declare relationships
  {
    // declare quad relationships
    for (unsigned i = 0; i < num_quad_faces; ++i) {
      Entity& face = *(quad_faces[i]);
      unsigned face_id = m_quad_face_ids[i];
      unsigned row = (face_id - 1) / 4;

      unsigned node_id = num_quad_faces + face_id + row;
      int chg_list[4] = {0, 5, 1, -5}; // (right-hand rule) counterclockwise
      for (unsigned chg_itr = 0; chg_itr < 4; ++chg_itr) {
        node_id += chg_list[chg_itr];
        Entity& node = *(all_entities[node_id - 1]);
        m_bulk_data->declare_relation( face , node , chg_itr);
      }
    }

    // declare shell relationships
    unsigned node_list[5] = {20, 25, 30, 35, 40};
    for (unsigned i = 0; i < num_shell_faces; ++i) {
      Entity& shell = *(shell_faces[i]);
      Entity& node1 = *all_entities[node_list[i] - 1];
      Entity& node2 = *all_entities[node_list[i+1] - 1];
      m_bulk_data->declare_relation(shell, node1, 0);
      m_bulk_data->declare_relation(shell, node2, 1);
    }
  }
}
