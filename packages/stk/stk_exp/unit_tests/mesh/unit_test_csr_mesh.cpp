#ifndef __IBMCPP__
#include <sierra/mesh/csr/csr_mesh.hpp>
#include <sierra/mesh/fixture/csr_mesh_factory.hpp>

#include <sierra/mesh/fixture/hex_fixture.hpp>

#include <gtest/gtest.h>

#include <iostream>

TEST( csr_mesh, basic)
{
  using namespace sierra;
  using namespace sierra::mesh;
  using namespace sierra::mesh::fixture;
  using namespace sierra::mesh::details;

  //create a mesh of 2x3x4 elements, which is 3x4x5 nodes. (total of 24 elements, 60 nodes)
  unsigned nx=2, ny=3, nz=4;
  hex_fixture fixture(nx,ny,nz);

  fixture.generate_mesh();

  size_t num_entities = 0;
  BOOST_FOREACH( bucket_key bucket, fixture.m_mesh.get_buckets() ) {
//    std::cout << bucket << std::endl;
//    std::cout << "\tParts: ";
//    BOOST_FOREACH( part_key part, fixture.m_mesh.get_parts(bucket)) {
//      std::cout << part.base() << " ";
//    }
//    std::cout << std::endl;
//    std::cout << "\tNum Entities: " << fixture.m_mesh.num_entities(bucket) << std::endl;
    num_entities += fixture.m_mesh.num_entities(bucket);
//    BOOST_FOREACH( entity_key entity, fixture.m_mesh.get_entities(bucket)) {
//      std::cout << "\tNum Out relations: " << fixture.m_mesh.num_out_relations(entity) << std::endl;
//    }
  }

  EXPECT_EQ(num_entities, (size_t)(60+24));
  boost::shared_ptr<csr_mesh> meshptr = csr_mesh_factory::create_from_modifiable(fixture.m_mesh);
  csr_mesh& mesh = *meshptr;

  //The csr-mesh should have the same number of buckets as the modifiable-mesh:
  EXPECT_EQ(boost::distance(fixture.m_mesh.get_buckets()),
            boost::distance(mesh.get_buckets()));

  size_t csr_num_entities = 0;
  csr_mesh::bucket_range buckets = mesh.get_buckets();
  BOOST_FOREACH( bucket_key bucket, buckets ) {
//    std::cout << bucket << std::endl;
//    std::cout << "\tParts: ";
//    const csr_mesh::part_set& parts = mesh.get_parts(bucket);
//    BOOST_FOREACH( part_key part, parts) {
//      std::cout << part.base() << " ";
//    }
//    std::cout << std::endl;
//    std::cout << "\tNum Entities: " << mesh.num_entities(bucket) << std::endl;
    csr_mesh::entity_descriptor_range entities = mesh.get_entities(bucket);
    csr_num_entities += boost::distance(entities);
//    BOOST_FOREACH( entity_descriptor entity, entities) {
//      std::cout << "\tNum Node relations: " << mesh.num_relations(entity,entity_rank(0)) << std::endl;
//      //std::cout << "\tNum Element relations: " << mesh.num_relations(entity,entity_rank(4)) << std::endl;
//    }
  }

  EXPECT_EQ(csr_num_entities, num_entities);
}

#endif
