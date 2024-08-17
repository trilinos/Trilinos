// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/GridFixture.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <vector>                       // for vector
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc

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
  : m_spatial_dimension(2),
    m_bulk_data_ptr(stk::unit_test_util::build_mesh(2, pm)),
    m_bulk_data(*m_bulk_data_ptr),
    m_fem_meta(m_bulk_data.mesh_meta_data()),
    m_quad_part( m_fem_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D) ),
    m_dead_part( m_fem_meta.declare_part("dead_part"))
{
}

GridFixture::~GridFixture()
{ }

void GridFixture::generate_grid()
{
  const unsigned num_nodes = 25;
  const unsigned num_quad_faces = 16;
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned p_size = m_bulk_data.parallel_size();

  for (unsigned i_rank = 0; i_rank < p_size; ++i_rank) {
    fill_node_map(num_nodes, num_quad_faces, i_rank);
  }

  std::vector<Entity> all_entities;

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
      Entity face = m_bulk_data.declare_element(face_id, face_parts);

      unsigned node_id = num_quad_faces + face_id + row;

      for (unsigned chg_itr = 0; chg_itr < num_nodes_per_quad; ++chg_itr) {
        node_id += stencil_for_4x4_quad_mesh[chg_itr];
        Entity node = m_bulk_data.declare_node(node_id, no_parts);
        m_bulk_data.declare_relation( face , node , chg_itr);
        DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, node_id, node);
      }
    }
  }
}

void GridFixture::fill_node_map(unsigned num_nodes, unsigned num_quad_faces, int p_rank)
{
  const unsigned p_size = m_bulk_data.parallel_size();
  std::vector<Entity> all_entities;

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

  // Iterate entities such that entity_id - 1 is the index of the
  // entity in the all_entities vector
  {
    const PartVector no_parts;
    const unsigned first_quad = (p_rank * num_quad_faces) / p_size;
    const unsigned end_quad = ((p_rank + 1) * num_quad_faces) / p_size;

    const unsigned num_nodes_per_quad = 4;
    // (right-hand rule) counterclockwise:
    const int stencil_for_4x4_quad_mesh[num_nodes_per_quad] = {0, 5, 1, -5};
    for (unsigned i = first_quad; i < end_quad; ++i) {

      unsigned face_id = quad_face_ids[i];
      unsigned row = (face_id - 1) / num_nodes_per_quad;
      unsigned node_id = num_quad_faces + face_id + row;

      for (unsigned chg_itr = 0; chg_itr < num_nodes_per_quad; ++chg_itr) {
        node_id += stencil_for_4x4_quad_mesh[chg_itr];
        AddToNodeProcsMMap(m_nodes_to_procs, node_id, p_rank);
      }
    }
  }
}

namespace simple_fields {

GridFixture::GridFixture(stk::ParallelMachine pm)
  : m_spatial_dimension(2),
    m_bulk_data_ptr(stk::unit_test_util::build_mesh(2, pm)),
    m_bulk_data(*m_bulk_data_ptr),
    m_fem_meta(m_bulk_data.mesh_meta_data()),
    m_quad_part( m_fem_meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D) ),
    m_dead_part( m_fem_meta.declare_part("dead_part"))
{
}

GridFixture::~GridFixture()
{ }

void GridFixture::generate_grid()
{
  const unsigned num_nodes = 25;
  const unsigned num_quad_faces = 16;
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned p_size = m_bulk_data.parallel_size();

  for (unsigned i_rank = 0; i_rank < p_size; ++i_rank) {
    fill_node_map(num_nodes, num_quad_faces, i_rank);
  }

  std::vector<Entity> all_entities;

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
      Entity face = m_bulk_data.declare_element(face_id, face_parts);

      unsigned node_id = num_quad_faces + face_id + row;

      for (unsigned chg_itr = 0; chg_itr < num_nodes_per_quad; ++chg_itr) {
        node_id += stencil_for_4x4_quad_mesh[chg_itr];
        Entity node = m_bulk_data.declare_node(node_id, no_parts);
        m_bulk_data.declare_relation( face , node , chg_itr);
        stk::mesh::fixtures::DoAddNodeSharings(m_bulk_data, m_nodes_to_procs, node_id, node);
      }
    }
  }
}

void GridFixture::fill_node_map(unsigned num_nodes, unsigned num_quad_faces, int p_rank)
{
  const unsigned p_size = m_bulk_data.parallel_size();
  std::vector<Entity> all_entities;

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

  // Iterate entities such that entity_id - 1 is the index of the
  // entity in the all_entities vector
  {
    const PartVector no_parts;
    const unsigned first_quad = (p_rank * num_quad_faces) / p_size;
    const unsigned end_quad = ((p_rank + 1) * num_quad_faces) / p_size;

    const unsigned num_nodes_per_quad = 4;
    // (right-hand rule) counterclockwise:
    const int stencil_for_4x4_quad_mesh[num_nodes_per_quad] = {0, 5, 1, -5};
    for (unsigned i = first_quad; i < end_quad; ++i) {

      unsigned face_id = quad_face_ids[i];
      unsigned row = (face_id - 1) / num_nodes_per_quad;
      unsigned node_id = num_quad_faces + face_id + row;

      for (unsigned chg_itr = 0; chg_itr < num_nodes_per_quad; ++chg_itr) {
        node_id += stencil_for_4x4_quad_mesh[chg_itr];
        stk::mesh::fixtures::AddToNodeProcsMMap(m_nodes_to_procs, node_id, p_rank);
      }
    }
  }
}

} // namespace simple_fields

} // fixtures
} // mesh
} // stk

