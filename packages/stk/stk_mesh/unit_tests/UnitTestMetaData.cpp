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

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>

#include <stk_mesh/base/FieldRelation.hpp>
#include <stk_mesh/base/PartRelation.hpp>
#include <stk_mesh/fem/Stencils.hpp>

using stk::mesh::MetaData;
using stk::mesh::TopologicalMetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::PartRelation;
using stk::mesh::EntityRank;
using std::cout;
using std::endl;

namespace stk {
namespace mesh {

class UnitTestMetaData {
public:
  UnitTestMetaData();
  ~UnitTestMetaData() {}

  const int spatial_dimension;
  MetaData metadata_committed;
  MetaData metadata_not_committed;
  MetaData metadata;
  TopologicalMetaData top;
  Part &pa;
  Part &pb;
  Part &pc;
  Part &pd;
  Part &pe;
  Part &pf;
  Part &pg;
  Part &ph;
  Part &partLeft; 
  PartVector part_vector;
  const CellTopologyData * singleton;

};


UnitTestMetaData::UnitTestMetaData()
  : spatial_dimension(3)
  , metadata_committed( TopologicalMetaData::entity_rank_names(spatial_dimension) )
  , metadata_not_committed( TopologicalMetaData::entity_rank_names(spatial_dimension) )
  , metadata( TopologicalMetaData::entity_rank_names(spatial_dimension) )
  , top( metadata_committed, spatial_dimension )
  , pa( metadata.declare_part( std::string("a") , 0 ) )
  , pb( metadata.declare_part( std::string("b") , 0 ) )
  , pc( metadata.declare_part( std::string("c") , 0 ) )
  , pd( metadata.declare_part( std::string("d") , 0 ) )
  , pe( metadata.declare_part( std::string("e") , 0 ) )
  , pf( metadata.declare_part( std::string("f") , 0 ) )
  , pg( metadata.declare_part( std::string("g") , 0 ) )
  , ph( metadata.declare_part( std::string("h") , 0 ) )
  , partLeft (top.declare_part<shards::Tetrahedron<4> >( "block_left" ))
  , part_vector()
  , singleton( NULL )

{
 metadata_committed.commit();

}


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
  UnitTestMetaData umeta;

  STKUNIT_ASSERT_THROW(umeta.metadata_not_committed.assert_committed("test throw"),std::logic_error);

  STKUNIT_ASSERT_THROW(umeta.metadata_committed.assert_not_committed("test throw"),std::logic_error);

  STKUNIT_ASSERT_THROW(umeta.metadata_not_committed.assert_same_mesh_meta_data("test throw",umeta.metadata_committed),std::logic_error);

  //test get_part with part that does not exist
  std::string test_string = "this_part_does_not_exist";
  STKUNIT_ASSERT_THROW(umeta.metadata_committed.get_part(test_string,"test_throw"),std::runtime_error);

  //test get_part with valid part
  STKUNIT_ASSERT(umeta.metadata_committed.get_part("block_left","do_not_throw"));

  //test declare part
  umeta.metadata.declare_part_relation(umeta.pe,stencil_test_function,umeta.pg);

  umeta.part_vector.push_back(&umeta.pa);
  umeta.part_vector.push_back(&umeta.pb);
  umeta.part_vector.push_back(&umeta.pc);
  umeta.part_vector.push_back(&umeta.pd);

  //Part * const intersection_part = &
  umeta.metadata.declare_part(umeta.part_vector);
 
  //Test declare_part_subset
  STKUNIT_ASSERT_THROW( umeta.metadata.declare_part_subset(umeta.pe,umeta.pe), std::runtime_error);
 
  //Test declare_part_relation with parts that are not subsets of each other
  STKUNIT_ASSERT_THROW( umeta.metadata.declare_part_relation(umeta.pg,stencil_test_function,umeta.ph), std::runtime_error);

  //Test declare_part_relation with a NULL stencil function
  STKUNIT_ASSERT_THROW( umeta.metadata.declare_part_relation(umeta.pe,NULL,umeta.pe), std::runtime_error);

  //Test declare_part_relation with parts that are subsets of each other
  umeta.metadata.declare_part_subset(umeta.pd,umeta.pf);
  STKUNIT_ASSERT_THROW( umeta.metadata.declare_part_relation(umeta.pd,stencil_test_function,umeta.pf), std::runtime_error);
    
  umeta.metadata.commit();
}

STKUNIT_UNIT_TEST( UnitTestMetaData, rankHigherThanDefined )
{
  //Test function entity_rank_name in MetaData.cpp
  UnitTestMetaData umeta;
  int i = 2;

  const std::string& i_name2 = umeta.metadata_committed.entity_rank_name( i );

  const std::vector<std::string> & type_names = TopologicalMetaData::entity_rank_names(umeta.spatial_dimension);
   
  STKUNIT_ASSERT( i_name2 == type_names[i] );
   
  EntityRank one_rank_higher_than_defined = type_names.size();
    
  STKUNIT_ASSERT_THROW( 
    umeta.metadata_committed.entity_rank_name( one_rank_higher_than_defined ),
    std::runtime_error
    );
}

STKUNIT_UNIT_TEST( UnitTestMetaData, testEntityRepository )
{
  //Test Entity repository - covering EntityRepository.cpp/hpp

  stk::mesh::MetaData meta ( stk::mesh::fem_entity_rank_names() );
  stk::mesh::Part & part = meta.declare_part( "another part");

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
  UnitTestMetaData umeta;
  std::vector<std::string> empty_names;
  STKUNIT_ASSERT_THROW(
      MetaData metadata5(empty_names),
      std::runtime_error
      );
}

STKUNIT_UNIT_TEST( UnitTestMetaData, declare_attribute_no_delete )
{
  //Coverage of declare_attribute_no_delete in MetaData.hpp 
  UnitTestMetaData umeta;
  umeta.metadata.declare_attribute_no_delete(umeta.pa, umeta.singleton);
}

}
//----------------------------------------------------------------------

}
}


