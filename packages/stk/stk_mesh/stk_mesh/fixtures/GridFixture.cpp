/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/fixtures/GridFixture.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>

/*
The following fixture creates the mesh below
1-16 Quadrilateral<4>
17-41 Nodes

17---18---19---20---21
|  1 |  2 |  3 |  4 |
22---23---24---25---26
|  5 |  6 |  7 |  8 |
27---28---29---30---31
|  9 | 10 | 11 | 12 |
32---33---34---35---36
| 13 | 14 | 15 | 16 |
37---38---39---40---41
*/

namespace stk {
namespace mesh {
namespace fixtures {


GridFixture::GridFixture(stk::ParallelMachine pm)
  : m_spatial_dimension(2)
  , m_meta_data( fem::entity_rank_names(m_spatial_dimension) )
  , m_bulk_data( m_meta_data , pm )
  , m_fem(  m_meta_data, m_spatial_dimension )
  , m_quad_part( declare_part<shards::Quadrilateral<4> >(m_meta_data, "quad_part") )
  , m_dead_part( declare_part(m_meta_data, "dead_part"))
{}

GridFixture::~GridFixture()
{ }

void GridFixture::generate_grid()
{
  const unsigned num_nodes = 25;
  const unsigned num_quad_faces = 16;
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned p_size = m_bulk_data.parallel_size();
  const EntityRank element_rank = fem::element_rank(m_fem);
  std::vector<Entity*> all_entities;

  // assign ids, quads, nodes, then shells
  // (we need this order to be this way in order for our connectivity setup to  work)
  std::vector<unsigned> quad_face_ids(num_quad_faces);
  std::vector<unsigned> node_ids(num_nodes);
  {
    unsigned curr_id = 1;
    for (unsigned  i = 0 ; i < num_quad_faces; ++i, ++curr_id) {
      quad_face_ids[i] = curr_id;
    }
    for (unsigned  i = 0 ; i < num_nodes; ++i, ++curr_id) {
      node_ids[i] = curr_id;
    }
  }

  // Note:  This block of code would normally be replaced with a call to stk_io
  // to generate the mesh.

  // declare entities such that entity_id - 1 is the index of the
  // entity in the all_entities vector
  {
    const PartVector no_parts;
    const unsigned first_quad = (p_rank * num_quad_faces) / p_size;
    const unsigned end_quad = ((p_rank + 1) * num_quad_faces) / p_size;

    // declare faces
    PartVector face_parts;
    face_parts.push_back(&m_quad_part);
    const unsigned num_nodes_per_quad = 4;
    // (right-hand rule) counterclockwise:
    const int stencil_for_4x4_quad_mesh[num_nodes_per_quad] = {0, 5, 1, -5};
    for (unsigned i = first_quad; i < end_quad; ++i) {

      unsigned face_id = quad_face_ids[i];
      unsigned row = (face_id - 1) / num_nodes_per_quad;

      Entity& face = m_bulk_data.declare_entity(element_rank, face_id, face_parts);

      unsigned node_id = num_quad_faces + face_id + row;

      for (unsigned chg_itr = 0; chg_itr < num_nodes_per_quad; ++chg_itr) {
        node_id += stencil_for_4x4_quad_mesh[chg_itr];
        Entity& node = m_bulk_data.declare_entity(fem::NODE_RANK, node_id, no_parts);
        m_bulk_data.declare_relation( face , node , chg_itr);
      }
    }
  }
}

} // fixtures
} // mesh
} // stk

