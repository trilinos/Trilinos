/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <unit_tests/UnitTestFixture.hpp>

namespace sunit {


  Stk_Mesh_Fixture::Stk_Mesh_Fixture( const std::vector<std::string> & entity_type_names )
    : m_MetaData( entity_type_names )
    , m_BulkData( m_MetaData , MPI_COMM_WORLD )
  { }

  Stk_Mesh_Fixture::~Stk_Mesh_Fixture() {}

  std::vector<stk::mesh::Part*> getPartVector(
      const sunit::Stk_Mesh_Fixture & fix,
      const std::vector<std::string> & names
      )
  {
    std::vector<stk::mesh::Part*> partVector;
    const stk::mesh::MetaData & metaData = fix.get_MetaData();
    const std::vector<std::string>::const_iterator it_begin = names.begin();
    std::vector<std::string>::const_iterator it = it_begin;
    const std::vector<std::string>::const_iterator it_end = names.end();
    for ( ; it != it_end ; ++it ) {
      stk::mesh::Part *part = metaData.get_part(*it);
      partVector.push_back(part);
    }
    return partVector;
  }


  stk::mesh::Part* getPart(
      const sunit::Stk_Mesh_Fixture & fix,
      const std::string name
      )
  {
    std::vector<std::string> names;
    names.push_back(name);
    std::vector<stk::mesh::Part*> partVector = getPartVector(fix,names);
    return partVector[0];
  }


  const stk::mesh::Bucket & getBucketContainingEntity(
      const sunit::Stk_Mesh_Fixture & fix,
      stk::mesh::EntityType ent_type,
      stk::mesh::EntityId ent_id
      )
  {
    const stk::mesh::BulkData & bulkData = fix.get_BulkData();
    return bulkData.get_entity(ent_type,ent_id)->bucket();
  }


  ExampleFixture::~ExampleFixture() {}

  ExampleFixture::ExampleFixture()
    : Stk_Mesh_Fixture( std::vector<std::string>( 1 , std::string("MyEntityType") ) )
  {
    // Create Parts and commit:
    stk::mesh::MetaData & metaData = get_NonconstMetaData();
    std::string myPartName;
    stk::mesh::EntityType myRank = 0;

    myPartName = "PartA";
    metaData.declare_part(myPartName,myRank);

    myPartName = "PartB";
    metaData.declare_part(myPartName,myRank);

    myPartName = "PartC";
    metaData.declare_part(myPartName,myRank);

    myPartName = "PartD";
    metaData.declare_part(myPartName,myRank);

    metaData.commit();

    // Create Entities and assign to parts:
    {
      stk::mesh::BulkData & bulkData = get_NonconstBulkData();
      bulkData.modification_begin();
      stk::mesh::EntityType ent_type = 0; // rank
      stk::mesh::EntityId ent_id = 1; // Unique ID
      std::vector<stk::mesh::Part*> partMembership;

      // Entity1 is contained in PartA
      partMembership.clear();
      partMembership.push_back(metaData.get_part("PartA"));
      {
        bulkData.declare_entity(ent_type, ent_id, partMembership);
      }
      ++ent_id;

      // Entity2 is contained in PartA and PartB
      partMembership.clear();
      partMembership.push_back(metaData.get_part("PartA"));
      partMembership.push_back(metaData.get_part("PartB"));
      {
        bulkData.declare_entity(ent_type, ent_id, partMembership);
      }
      ++ent_id;

      // Entity3 is contained in PartB and PartC
      partMembership.clear();
      partMembership.push_back(metaData.get_part("PartB"));
      partMembership.push_back(metaData.get_part("PartC"));
      {
        bulkData.declare_entity(ent_type, ent_id, partMembership);
      }
      ++ent_id;

      // Entity4 is contained in PartC
      partMembership.clear();
      partMembership.push_back(metaData.get_part("PartC"));
      {
        bulkData.declare_entity(ent_type, ent_id, partMembership);
      }
      ++ent_id;

      // Entity5 is not contained in any Part
      partMembership.clear();
      {
        bulkData.declare_entity(ent_type, ent_id, partMembership);
      }
      ++ent_id;
    }
  }

  const stk::mesh::Bucket & getExampleBucket(
      const sunit::Stk_Mesh_Fixture & fix,
      stk::mesh::EntityId ent_id
      )
  {
    stk::mesh::EntityType ent_type = 0;
    const stk::mesh::BulkData & bulkData = fix.get_BulkData();
    return bulkData.get_entity(ent_type,ent_id)->bucket();
  }



  stk::mesh::Part* getExamplePart(
      const sunit::Stk_Mesh_Fixture & fix,
      const std::string name
      )
  {
    stk::mesh::Part* part;
    if (name == "PartU") {
      part = &(fix.get_MetaData().universal_part());
    } else {
      std::vector<std::string> names;
      names.push_back(name);
      std::vector<stk::mesh::Part*> partVector = getPartVector(fix,names);
      part = partVector[0];
    }
    return part;
  }



} // namespace sunit

