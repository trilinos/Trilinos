/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_UNITTEST_RING_MESH_FIXTURE_HPP
#define STK_MESH_UNITTEST_RING_MESH_FIXTURE_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

class RingMeshFixture {
public:
  stk::mesh::MetaData   m_meta_data;
  stk::mesh::BulkData   m_bulk_data;
  stk::mesh::PartVector m_edge_parts ;
  stk::mesh::Part     & m_edge_part_extra ;
  const size_t          m_num_edge_per_proc ;
  std::vector<stk::mesh::EntityId> m_node_ids , m_edge_ids ;

  RingMeshFixture( stk::ParallelMachine pm ,
                   unsigned num_edge_per_proc = 10 ,
                   bool use_edge_parts = false );

  ~RingMeshFixture();

  // Testing for a simple loop of mesh entities:
  // node[i] : edge[i] : node[ ( i + 1 ) % node.size() ]
  void generate_mesh( bool generate_aura = true );

  void test_shift_loop( bool generate_aura );

private:

   RingMeshFixture();
   RingMeshFixture( const RingMeshFixture & );
   RingMeshFixture & operator = ( const RingMeshFixture & );
};

#endif
