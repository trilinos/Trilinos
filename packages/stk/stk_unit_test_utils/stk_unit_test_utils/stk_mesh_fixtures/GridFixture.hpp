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

#ifndef STK_MESH_FIXTURES_GRID_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_GRID_MESH_FIXTURE_HPP

#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include "stk_unit_test_utils/stk_mesh_fixtures/FixtureNodeSharing.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace fixtures {

class GridFixture
{
public:
  GridFixture(stk::ParallelMachine pm);

  ~GridFixture();

  MetaData& fem_meta() { return m_fem_meta; }
  BulkData& bulk_data() { return m_bulk_data; }

  Part* quad_part() const { return & m_quad_part; }
  Part* dead_part() const { return & m_dead_part; }

  void generate_grid();

  const unsigned m_spatial_dimension;

  std::shared_ptr<BulkData> m_bulk_data_ptr;
  BulkData&  m_bulk_data;
  MetaData&  m_fem_meta;
  Part &    m_quad_part;
  Part &    m_dead_part;

private:
  NodeToProcsMMap m_nodes_to_procs;

  void fill_node_map(unsigned num_nodes, unsigned num_quad_faces, int p_rank);
};

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
GridFixture
{
public:
  GridFixture(stk::ParallelMachine pm);

  ~GridFixture();

  MetaData& fem_meta() { return m_fem_meta; }
  BulkData& bulk_data() { return m_bulk_data; }

  Part* quad_part() const { return & m_quad_part; }
  Part* dead_part() const { return & m_dead_part; }

  void generate_grid();

  const unsigned m_spatial_dimension;

  std::shared_ptr<BulkData> m_bulk_data_ptr;
  BulkData&  m_bulk_data;
  MetaData&  m_fem_meta;
  Part &    m_quad_part;
  Part &    m_dead_part;

private:
  NodeToProcsMMap m_nodes_to_procs;

  void fill_node_map(unsigned num_nodes, unsigned num_quad_faces, int p_rank);
};

} // namespace simple_fields

} // fixtures
} // mesh
} // stk

#endif

