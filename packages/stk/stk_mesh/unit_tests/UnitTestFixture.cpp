
#include <unit_tests/UnitTestFixture.hpp>

namespace sunit {


  Stk_Mesh_Fixture::Stk_Mesh_Fixture()
  {
  }

  Stk_Mesh_Fixture::~Stk_Mesh_Fixture()
  {
  }

  void Stk_Mesh_Fixture::createMetaData(const std::vector<std::string> & entity_type_names) 
  {
    m_MetaData = Teuchos::rcp(new stk::mesh::MetaData(entity_type_names));
  }


  void Stk_Mesh_Fixture::createBulkData()
  {
    // Create BulkData:
    // Note:  MetaData.commit() bust have been called prior to creating BulkData.
    m_BulkData = Teuchos::rcp(new stk::mesh::BulkData(*m_MetaData,MPI_COMM_WORLD));
  }


  Teuchos::RCP<const stk::mesh::BulkData> Stk_Mesh_Fixture::get_BulkData() const
  {
    return m_BulkData;
  }


  Teuchos::RCP<stk::mesh::BulkData> Stk_Mesh_Fixture::get_NonconstBulkData() 
  {
    return m_BulkData;
  }


  Teuchos::RCP<const stk::mesh::MetaData> Stk_Mesh_Fixture::get_MetaData() const
  {
    return m_MetaData;
  }


  Teuchos::RCP<stk::mesh::MetaData> Stk_Mesh_Fixture::get_NonconstMetaData() 
  {
    return m_MetaData;
  }
  


  Teuchos::RCP<Stk_Mesh_Fixture> stk_Mesh_Fixture()
  {
    return Teuchos::rcp(new Stk_Mesh_Fixture);
  }


  std::vector<stk::mesh::Part*> getPartVector(
      const sunit::Stk_Mesh_Fixture & fix,
      const std::vector<std::string> & names
      )
  {
    std::vector<stk::mesh::Part*> partVector;
    Teuchos::RCP<const stk::mesh::MetaData> metaData = fix.get_MetaData();
    const std::vector<std::string>::const_iterator it_begin = names.begin();
    std::vector<std::string>::const_iterator it = it_begin;
    const std::vector<std::string>::const_iterator it_end = names.end();
    for ( ; it != it_end ; ++it ) {
      stk::mesh::Part *part = metaData->get_part(*it);
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
    Teuchos::RCP<const stk::mesh::BulkData> bulkData = fix.get_BulkData();
    return bulkData->get_entity(ent_type,ent_id)->bucket();
  }


  Teuchos::RCP<Stk_Mesh_Fixture> getExampleFixture()
  {
    // Create Mesh Meta Data:
    Teuchos::RCP<Stk_Mesh_Fixture> fix = sunit::stk_Mesh_Fixture();
    std::vector<std::string> entity_type_names;
    entity_type_names.push_back("MyEntityType");
    fix->createMetaData(entity_type_names);

    // Create Parts and commit:
    Teuchos::RCP<stk::mesh::MetaData> metaData = fix->get_NonconstMetaData();
    std::string myPartName;
    stk::mesh::EntityType myRank = 0;

    myPartName = "PartA";
    metaData->declare_part(myPartName,myRank);

    myPartName = "PartB";
    metaData->declare_part(myPartName,myRank);

    myPartName = "PartC";
    metaData->declare_part(myPartName,myRank);

    myPartName = "PartD";
    metaData->declare_part(myPartName,myRank);

    metaData->commit();

    // Create BulkData
    fix->createBulkData();

    // Create Entities and assign to parts:
    {
      Teuchos::RCP<stk::mesh::BulkData> bulkData = fix->get_NonconstBulkData();
      stk::mesh::EntityType ent_type = 0; // rank
      stk::mesh::EntityId ent_id = 1; // Unique ID
      std::vector<stk::mesh::Part*> partMembership;

      // Entity1 is contained in PartA
      partMembership.clear();
      partMembership.push_back(metaData->get_part("PartA"));
      {
        bulkData->declare_entity(ent_type, ent_id, partMembership);
      }
      ++ent_id;

      // Entity2 is contained in PartA and PartB
      partMembership.clear();
      partMembership.push_back(metaData->get_part("PartA"));
      partMembership.push_back(metaData->get_part("PartB"));
      {
        bulkData->declare_entity(ent_type, ent_id, partMembership);
      }
      ++ent_id;

      // Entity3 is contained in PartB and PartC
      partMembership.clear();
      partMembership.push_back(metaData->get_part("PartB"));
      partMembership.push_back(metaData->get_part("PartC"));
      {
        bulkData->declare_entity(ent_type, ent_id, partMembership);
      }
      ++ent_id;

      // Entity4 is contained in PartC
      partMembership.clear();
      partMembership.push_back(metaData->get_part("PartC"));
      {
        bulkData->declare_entity(ent_type, ent_id, partMembership);
      }
      ++ent_id;

      // Entity5 is not contained in any Part
      partMembership.clear();
      {
        bulkData->declare_entity(ent_type, ent_id, partMembership);
      }
      ++ent_id;
    }
    return fix;
  }

  const stk::mesh::Bucket & getExampleBucket(
      const sunit::Stk_Mesh_Fixture & fix,
      stk::mesh::EntityId ent_id
      )
  {
    stk::mesh::EntityType ent_type = 0;
    Teuchos::RCP<const stk::mesh::BulkData> bulkData = fix.get_BulkData();
    return bulkData->get_entity(ent_type,ent_id)->bucket();
  }



  stk::mesh::Part* getExamplePart(
      const sunit::Stk_Mesh_Fixture & fix,
      const std::string name
      )
  {
    stk::mesh::Part* part;
    if (name == "PartU") {
      part = &(fix.get_MetaData()->universal_part());
    } else {
      std::vector<std::string> names;
      names.push_back(name);
      std::vector<stk::mesh::Part*> partVector = getPartVector(fix,names);
      part = partVector[0];
    }
    return part;
  }



} // namespace sunit

