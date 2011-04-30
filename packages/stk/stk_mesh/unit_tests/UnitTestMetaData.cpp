/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>

#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>

#include <stk_mesh/base/FieldRelation.hpp>
#include <stk_mesh/base/PartRelation.hpp>

#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

using stk::mesh::MetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::PartRelation;
using stk::mesh::EntityRank;
using std::cout;
using std::endl;

//----------------------------------------------------------------------

namespace {
int stencil_test_function( unsigned  from_type ,
                           unsigned  to_type ,
                           unsigned  identifier )
{
  return 0;
}



STKUNIT_UNIT_TEST( UnitTestMetaData, testMetaData )
{
  //Test functions in MetaData.cpp
  const int spatial_dimension = 3;
  const std::vector<std::string> & rank_names = stk::mesh::fem::entity_rank_names(spatial_dimension);
  MetaData metadata_committed(rank_names);
  MetaData metadata_not_committed(rank_names);
  MetaData metadata(rank_names);

  Part &pa = metadata.declare_part( std::string("a") , 0 );
  Part &pb = metadata.declare_part( std::string("b") , 0 );
  Part &pc = metadata.declare_part( std::string("c") , 0 );
  Part &pd = metadata.declare_part( std::string("d") , 0 );
  Part &pe = metadata.declare_part( std::string("e") , 0 );
  Part &pf = metadata.declare_part( std::string("f") , 0 );
  Part &pg = metadata.declare_part( std::string("g") , 0 );
  Part &ph = metadata.declare_part( std::string("h") , 0 );
  PartVector part_vector;
  metadata_committed.commit();

  //test get_part with part that does not exist
  std::string test_string = "this_part_does_not_exist";
  STKUNIT_ASSERT_THROW( metadata_committed.get_part(test_string,"test_throw"),std::runtime_error);

  //test get_part with valid part
  STKUNIT_ASSERT( metadata.get_part(std::string("a"),"do_not_throw"));

  //test declare part
  metadata.declare_part_relation( pe,stencil_test_function, pg);

  part_vector.push_back(& pa);
  part_vector.push_back(& pb);
  part_vector.push_back(& pc);
  part_vector.push_back(& pd);

  //Part * const intersection_part = &
  metadata.declare_part(part_vector);

  //Test declare_part_subset
  STKUNIT_ASSERT_THROW(  metadata.declare_part_subset( pe, pe), std::runtime_error);

  //Test declare_part_relation with parts that are not subsets of each other
  STKUNIT_ASSERT_THROW(  metadata.declare_part_relation( pg,stencil_test_function, ph), std::logic_error);

  //Test declare_part_relation with a NULL stencil function
  STKUNIT_ASSERT_THROW(  metadata.declare_part_relation( pe,NULL, pe), std::runtime_error);

  //Test declare_part_relation with parts that are subsets of each other
  metadata.declare_part_subset( pd, pf);
  STKUNIT_ASSERT_THROW(  metadata.declare_part_relation( pd,stencil_test_function, pf), std::runtime_error);

  metadata.commit();
}

STKUNIT_UNIT_TEST( UnitTestMetaData, rankHigherThanDefined )
{
  //Test function entity_rank_name in MetaData.cpp
  const int spatial_dimension = 3;
  const std::vector<std::string> & rank_names = stk::mesh::fem::entity_rank_names(spatial_dimension);
  MetaData metadata(rank_names);
  int i = 2;

  const std::string& i_name2 =  metadata.entity_rank_name( i );

  STKUNIT_ASSERT( i_name2 == rank_names[i] );

  EntityRank one_rank_higher_than_defined = rank_names.size();

  STKUNIT_ASSERT_THROW(
    metadata.entity_rank_name( one_rank_higher_than_defined ),
    std::runtime_error
                        );
}

STKUNIT_UNIT_TEST( UnitTestMetaData, testEntityRepository )
{
  static const size_t spatial_dimension = 3;

  //Test Entity repository - covering EntityRepository.cpp/hpp
  stk::mesh::MetaData meta ( stk::mesh::fem::entity_rank_names(spatial_dimension) );
  stk::mesh::Part & part = meta.declare_part("another part");

  meta.commit();

  stk::mesh::BulkData bulk ( meta , MPI_COMM_WORLD , 100 );
  std::vector<stk::mesh::Part *>  add_part;
  add_part.push_back ( &part );

  int  size , rank;
  rank = stk::parallel_machine_rank( MPI_COMM_WORLD );
  size = stk::parallel_machine_size( MPI_COMM_WORLD );
  PartVector tmp(1);

  bulk.modification_begin();

  int id_base = 0;
  for ( id_base = 0 ; id_base < 97 ; ++id_base )
  {
    int new_id = size * id_base + rank;
    bulk.declare_entity( 0 , new_id+1 , add_part );
  }

  int new_id = size * (++id_base) + rank;
  stk::mesh::Entity & elem  = bulk.declare_entity( 3 , new_id+1 , add_part );

  //new_id = size * (++id_base) + rank;
  // stk::mesh::Entity & elem2  = bulk.declare_entity( 3 , new_id+1 , add_part );

  stk::mesh::impl::EntityRepository e;

  e.comm_clear( elem );

  e.comm_clear_ghosting( elem );

  const stk::mesh::Ghosting & ghost = bulk.shared_aura();

  bulk.modification_end();

  STKUNIT_ASSERT_FALSE(e.erase_ghosting(elem, ghost));

  const stk::mesh::EntityCommInfo comm_info( ghost.ordinal() , 0 );

  STKUNIT_ASSERT_FALSE(e.erase_comm_info(elem, comm_info));

  STKUNIT_ASSERT(e.insert_comm_info(elem, comm_info));

  //Checking internal_create_entity

  e.internal_create_entity( stk::mesh::EntityKey( 3, 2 ));
  e.internal_create_entity( stk::mesh::EntityKey( 3, 5 ));
  e.internal_create_entity( stk::mesh::EntityKey( 3, 7 ));

  //Checking get_entity with invalid key - no rank or id
  {
    int ok = 0 ;
    try {

      stk::mesh::Entity * elem3 = e.get_entity(stk::mesh::EntityKey());
      if(elem3){
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
  std::vector<std::string> empty_names;
  STKUNIT_ASSERT_THROW(
    MetaData metadata(empty_names),
    std::runtime_error
    );
}
STKUNIT_UNIT_TEST( UnitTestMetaData, declare_part_with_rank )
{
  //MetaData constructor fails because there are no entity types:
  const int spatial_dimension = 3;
  MetaData metadata(stk::mesh::fem::entity_rank_names(spatial_dimension));
  metadata.declare_part("foo");
  STKUNIT_ASSERT_NO_THROW(metadata.declare_part("foo",1));
  STKUNIT_ASSERT_NO_THROW(metadata.declare_part("foo",1));

  // Should throw because we're trying to change rank
  STKUNIT_ASSERT_THROW(metadata.declare_part("foo",2),std::runtime_error);

  // Should not throw since we did not provide rank
  metadata.declare_part("foo");
}

STKUNIT_UNIT_TEST( UnitTestMetaData, declare_attribute_no_delete )
{
  //Coverage of declare_attribute_no_delete in MetaData.hpp
  const CellTopologyData * singleton = NULL;
  const int spatial_dimension = 3;
  MetaData metadata(stk::mesh::fem::entity_rank_names(spatial_dimension));
  Part &pa = metadata.declare_part( std::string("a") , 0 );
  metadata.declare_attribute_no_delete( pa, singleton);
  metadata.commit();
}

}
//----------------------------------------------------------------------





