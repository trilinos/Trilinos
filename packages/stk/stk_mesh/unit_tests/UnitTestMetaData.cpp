/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stddef.h>                     // for NULL, size_t
#include <exception>                    // for exception
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names
#include <stk_mesh/baseImpl/EntityRepository.hpp>  // for EntityRepository
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <string>                       // for string, operator==, etc
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for AssertHelper
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/NamedPair.hpp"
namespace stk { namespace mesh { class Part; } }






using stk::mesh::MetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityRank;
using std::cout;
using std::endl;

//----------------------------------------------------------------------

namespace {

STKUNIT_UNIT_TEST( UnitTestMetaData, testMetaData )
{
  //Test functions in MetaData.cpp
  const int spatial_dimension = 3;
  MetaData metadata_committed(spatial_dimension);
  MetaData metadata_not_committed(spatial_dimension);
  MetaData metadata(spatial_dimension);

  Part &pa = metadata.declare_part( std::string("a") , 0 );
  Part &pb = metadata.declare_part( std::string("b") , 0 );
  Part &pc = metadata.declare_part( std::string("c") , 0 );
  Part &pd = metadata.declare_part( std::string("d") , 0 );
  Part &pe = metadata.declare_part( std::string("e") , 0 );
  PartVector part_vector;
  metadata_committed.commit();

  //test get_part with part that does not exist
  std::string test_string = "this_part_does_not_exist";
  STKUNIT_ASSERT_THROW( metadata_committed.get_part(test_string,"test_throw"),std::runtime_error);

  //test get_part with valid part
  STKUNIT_ASSERT( metadata.get_part(std::string("a"),"do_not_throw"));



  part_vector.push_back(& pa);
  part_vector.push_back(& pb);
  part_vector.push_back(& pc);
  part_vector.push_back(& pd);

  //Test declare_part_subset
  STKUNIT_ASSERT_THROW(  metadata.declare_part_subset( pe, pe), std::runtime_error);

  metadata.commit();
}

STKUNIT_UNIT_TEST( UnitTestMetaData, rankHigherThanDefined )
{
  //Test function entity_rank_name in MetaData.cpp
  const int spatial_dimension = 3;
  const std::vector<std::string> & rank_names = stk::mesh::entity_rank_names();
  MetaData metadata(spatial_dimension, rank_names);

  const std::string& i_name2 =  metadata.entity_rank_name( stk::topology::EDGE_RANK );

  STKUNIT_ASSERT( i_name2 == rank_names[stk::topology::EDGE_RANK] );

  EntityRank one_rank_higher_than_defined = static_cast<EntityRank>(rank_names.size());

  STKUNIT_ASSERT_THROW(
    metadata.entity_rank_name( one_rank_higher_than_defined ),
    std::runtime_error
                        );
}

STKUNIT_UNIT_TEST( UnitTestMetaData, testEntityRepository )
{
  static const size_t spatial_dimension = 3;

  //Test Entity repository - covering EntityRepository.cpp/hpp
  stk::mesh::MetaData meta ( spatial_dimension );
  stk::mesh::Part & part = meta.declare_part("another part");

  meta.commit();

  stk::mesh::BulkData bulk ( meta , MPI_COMM_WORLD , 100 );
  std::vector<stk::mesh::Part *>  add_part;
  add_part.push_back ( &part );

  int rank = stk::parallel_machine_rank( MPI_COMM_WORLD );
  int size = stk::parallel_machine_size( MPI_COMM_WORLD );
  PartVector tmp(1);

  bulk.modification_begin();

  int id_base = 0;
  for ( id_base = 0 ; id_base < 97 ; ++id_base )
  {
    int new_id = size * id_base + rank;
    bulk.declare_entity( stk::topology::NODE_RANK , new_id+1 , add_part );
  }

  int new_id = size * (++id_base) + rank;
  stk::mesh::Entity elem  = bulk.declare_entity( stk::topology::ELEMENT_RANK , new_id+1 , add_part );

  //new_id = size * (++id_base) + rank;
  // stk::mesh::Entity elem2  = bulk.declare_entity( stk::topology::ELEMENT_RANK , new_id+1 , add_part );

  stk::mesh::impl::EntityRepository &e = bulk.get_entity_repository();

  bulk.entity_comm_clear(bulk.entity_key(elem));

  bulk.entity_comm_clear_ghosting(bulk.entity_key(elem));

  const stk::mesh::Ghosting & ghost = bulk.shared_aura();

  bulk.modification_end();

  STKUNIT_ASSERT_FALSE(bulk.entity_comm_erase(bulk.entity_key(elem), ghost));

  const stk::mesh::EntityCommInfo comm_info( ghost.ordinal() , 0 );

  STKUNIT_ASSERT_FALSE(bulk.entity_comm_erase(bulk.entity_key(elem), comm_info));

  STKUNIT_ASSERT(bulk.entity_comm_insert(elem, comm_info));

  //Checking internal_create_entity.
  //   Hey, this doesn't seem to test much! -- PGX
  e.internal_create_entity( stk::mesh::EntityKey( stk::topology::ELEM_RANK, 2 ));
  e.internal_create_entity( stk::mesh::EntityKey( stk::topology::ELEM_RANK, 5 ));
  e.internal_create_entity( stk::mesh::EntityKey( stk::topology::ELEM_RANK, 7 ));

  //Checking get_entity with invalid key - no rank or id
  {
    int ok = 0 ;
    try {

      stk::mesh::Entity elem3 = e.get_entity(stk::mesh::EntityKey());
      if(bulk.is_valid(elem3)){
        // CAROL FIXME
      }

    }
    catch( const std::exception & x ) {
      ok = 1 ;
      std::cout << "UnitTestMetaData CORRECTLY caught error for : "
                << x.what()
                << std::endl ;
    }
    if ( ! ok ) {
      throw std::runtime_error("UnitTestMetaData FAILED to catch error for get_entity - invalid key");
    }
  }
}

STKUNIT_UNIT_TEST( UnitTestMetaData, noEntityTypes )
{
  //MetaData constructor fails because there are no entity types:
  std::vector<std::string> wrong_names(1, "foo");
  STKUNIT_ASSERT_THROW(
    MetaData metadata(3 /*dim*/, wrong_names),
    std::runtime_error
    );
}
STKUNIT_UNIT_TEST( UnitTestMetaData, declare_part_with_rank )
{
  //MetaData constructor fails because there are no entity types:
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);
  metadata.declare_part("foo");
  STKUNIT_ASSERT_NO_THROW(metadata.declare_part("foo",stk::topology::EDGE_RANK));
  STKUNIT_ASSERT_NO_THROW(metadata.declare_part("foo",stk::topology::EDGE_RANK));

  // Should throw because we're trying to change rank
  STKUNIT_ASSERT_THROW(metadata.declare_part("foo",stk::topology::FACE_RANK),std::runtime_error);

  // Should not throw since we did not provide rank
  metadata.declare_part("foo");
}

STKUNIT_UNIT_TEST( UnitTestMetaData, declare_attribute_no_delete )
{
  //Coverage of declare_attribute_no_delete in MetaData.hpp
  const CellTopologyData * singleton = NULL;
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);
  Part &pa = metadata.declare_part( std::string("a") , stk::topology::NODE_RANK );
  metadata.declare_attribute_no_delete( pa, singleton);
  metadata.commit();
}

}
//----------------------------------------------------------------------





