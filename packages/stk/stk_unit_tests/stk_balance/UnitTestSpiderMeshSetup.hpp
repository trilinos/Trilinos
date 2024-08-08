#ifndef UNITTESTSPIDERMESHSETUP_HPP
#define UNITTESTSPIDERMESHSETUP_HPP

#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_io/WriteMesh.hpp"

inline
void make_mesh_non_spider_no_volume_elements(stk::mesh::BulkData & bulk)
{
  std::string meshDesc = "2,1,BEAM_2,1,7\n"
                         "2,2,BEAM_2,2,7\n"
                         "2,3,BEAM_2,3,7\n"
                         "2,4,BEAM_2,4,7\n"
                         "3,5,BEAM_2,5,7\n"
                         "3,6,BEAM_2,6,7\n";

  std::vector<double> coordinates {
    0,0,1, 0,0,0, 1,0,1, 1,0,0, 2,0,1, 2,0,0,
    1,1,0.5
  };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

inline
void make_mesh_non_spider_not_enough_legs(stk::mesh::BulkData & bulk)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                         "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                         "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                         "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                         "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                         "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                         "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                         "2,9,BEAM_2,8,37\n"
                         "2,10,BEAM_2,9,37\n"
                         "2,11,BEAM_2,12,37\n"
                         "3,12,BEAM_2,13,37\n"
                         "3,13,BEAM_2,16,37\n";

  std::vector<double> coordinates {
    0,1,1, 0,0,1, 0,0,0, 0,1,0,
    1,1,1, 1,0,1, 1,0,0, 1,1,0,
    2,1,1, 2,0,1, 2,0,0, 2,1,0,
    3,1,1, 3,0,1, 3,0,0, 3,1,0,
    4,1,1, 4,0,1, 4,0,0, 4,1,0,
    5,1,1, 5,0,1, 5,0,0, 5,1,0,
    6,1,1, 6,0,1, 6,0,0, 6,1,0,
    7,1,1, 7,0,1, 7,0,0, 7,1,0,
    8,1,1, 8,0,1, 8,0,0, 8,1,0,
    2,2,0.5
  };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

inline
void make_mesh_one_spider_no_body_element(stk::mesh::BulkData & bulk)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                         "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                         "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                         "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                         "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                         "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                         "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                         "2,9,BEAM_2,5,37\n"
                         "2,10,BEAM_2,8,37\n"
                         "2,11,BEAM_2,9,37\n"
                         "2,12,BEAM_2,12,37\n"
                         "3,13,BEAM_2,13,37\n"
                         "3,14,BEAM_2,16,37\n";

  std::vector<double> coordinates {
    0,1,1, 0,0,1, 0,0,0, 0,1,0,
    1,1,1, 1,0,1, 1,0,0, 1,1,0,
    2,1,1, 2,0,1, 2,0,0, 2,1,0,
    3,1,1, 3,0,1, 3,0,0, 3,1,0,
    4,1,1, 4,0,1, 4,0,0, 4,1,0,
    5,1,1, 5,0,1, 5,0,0, 5,1,0,
    6,1,1, 6,0,1, 6,0,0, 6,1,0,
    7,1,1, 7,0,1, 7,0,0, 7,1,0,
    8,1,1, 8,0,1, 8,0,0, 8,1,0,
    2,2,0.5
  };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

inline
void make_mesh_one_spider_particle_body(stk::mesh::BulkData & bulk)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                         "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                         "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                         "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                         "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                         "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                         "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                         "2,9,BEAM_2,5,37\n"
                         "2,10,BEAM_2,8,37\n"
                         "2,11,BEAM_2,9,37\n"
                         "2,12,BEAM_2,12,37\n"
                         "3,13,BEAM_2,13,37\n"
                         "3,14,BEAM_2,16,37\n"
                         "3,15,PARTICLE,37\n";

  std::vector<double> coordinates {
    0,1,1, 0,0,1, 0,0,0, 0,1,0,
    1,1,1, 1,0,1, 1,0,0, 1,1,0,
    2,1,1, 2,0,1, 2,0,0, 2,1,0,
    3,1,1, 3,0,1, 3,0,0, 3,1,0,
    4,1,1, 4,0,1, 4,0,0, 4,1,0,
    5,1,1, 5,0,1, 5,0,0, 5,1,0,
    6,1,1, 6,0,1, 6,0,0, 6,1,0,
    7,1,1, 7,0,1, 7,0,0, 7,1,0,
    8,1,1, 8,0,1, 8,0,0, 8,1,0,
    2,2,0.5
  };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

inline
void make_mesh_one_spider_beam_body(stk::mesh::BulkData & bulk)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                         "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                         "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                         "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                         "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                         "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                         "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                         "2,9,BEAM_2,5,37\n"
                         "2,10,BEAM_2,8,37\n"
                         "2,11,BEAM_2,9,37\n"
                         "2,12,BEAM_2,12,37\n"
                         "3,13,BEAM_2,13,37\n"
                         "3,14,BEAM_2,16,37\n"
                         "3,15,BEAM_2,37,38\n";

  std::vector<double> coordinates {
    0,1,1, 0,0,1, 0,0,0, 0,1,0,
    1,1,1, 1,0,1, 1,0,0, 1,1,0,
    2,1,1, 2,0,1, 2,0,0, 2,1,0,
    3,1,1, 3,0,1, 3,0,0, 3,1,0,
    4,1,1, 4,0,1, 4,0,0, 4,1,0,
    5,1,1, 5,0,1, 5,0,0, 5,1,0,
    6,1,1, 6,0,1, 6,0,0, 6,1,0,
    7,1,1, 7,0,1, 7,0,0, 7,1,0,
    8,1,1, 8,0,1, 8,0,0, 8,1,0,
    2,2,0.5,
    2,3,0.5
  };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

inline
void make_mesh_compound_spider_beam_body(stk::mesh::BulkData & bulk)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                         "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                         "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                         "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                         "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                         "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                         "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                         "2,9,BEAM_2,5,37\n"
                         "2,10,BEAM_2,8,37\n"
                         "2,11,BEAM_2,9,37\n"
                         "2,12,BEAM_2,12,37\n"
                         "3,13,BEAM_2,13,37\n"
                         "3,14,BEAM_2,16,37\n"
                         "0,15,BEAM_2,21,38\n"
                         "0,16,BEAM_2,24,38\n"
                         "0,17,BEAM_2,25,38\n"
                         "0,18,BEAM_2,28,38\n"
                         "1,19,BEAM_2,29,38\n"
                         "1,20,BEAM_2,32,38\n"
                         "3,21,BEAM_2,37,38\n";

  std::vector<double> coordinates {
    0,1,1, 0,0,1, 0,0,0, 0,1,0,
    1,1,1, 1,0,1, 1,0,0, 1,1,0,
    2,1,1, 2,0,1, 2,0,0, 2,1,0,
    3,1,1, 3,0,1, 3,0,0, 3,1,0,
    4,1,1, 4,0,1, 4,0,0, 4,1,0,
    5,1,1, 5,0,1, 5,0,0, 5,1,0,
    6,1,1, 6,0,1, 6,0,0, 6,1,0,
    7,1,1, 7,0,1, 7,0,0, 7,1,0,
    8,1,1, 8,0,1, 8,0,0, 8,1,0,
    2,2,0.5,
    6,2,0.5
  };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

inline
void make_mesh_two_spiders_particle_body(stk::mesh::BulkData & bulk)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                         "1,3,HEX_8,9,10,11,12,13,14,15,16\n"
                         "1,4,HEX_8,13,14,15,16,17,18,19,20\n"
                         "2,5,HEX_8,17,18,19,20,21,22,23,24\n"
                         "2,6,HEX_8,21,22,23,24,25,26,27,28\n"
                         "3,7,HEX_8,25,26,27,28,29,30,31,32\n"
                         "3,8,HEX_8,29,30,31,32,33,34,35,36\n"
                         "2,9,BEAM_2,5,37\n"
                         "2,10,BEAM_2,8,37\n"
                         "2,11,BEAM_2,9,37\n"
                         "2,12,BEAM_2,12,37\n"
                         "3,13,BEAM_2,13,37\n"
                         "3,14,BEAM_2,16,37\n"
                         "3,15,PARTICLE,37\n"
                         "0,16,BEAM_2,21,38\n"
                         "0,17,BEAM_2,24,38\n"
                         "0,18,BEAM_2,25,38\n"
                         "0,19,BEAM_2,28,38\n"
                         "1,20,BEAM_2,29,38\n"
                         "1,21,BEAM_2,32,38\n"
                         "1,22,PARTICLE,38\n";

  std::vector<double> coordinates {
    0,1,1, 0,0,1, 0,0,0, 0,1,0,
    1,1,1, 1,0,1, 1,0,0, 1,1,0,
    2,1,1, 2,0,1, 2,0,0, 2,1,0,
    3,1,1, 3,0,1, 3,0,0, 3,1,0,
    4,1,1, 4,0,1, 4,0,0, 4,1,0,
    5,1,1, 5,0,1, 5,0,0, 5,1,0,
    6,1,1, 6,0,1, 6,0,0, 6,1,0,
    7,1,1, 7,0,1, 7,0,0, 7,1,0,
    8,1,1, 8,0,1, 8,0,0, 8,1,0,
    2,2,0.5,
    6,2,0.5
  };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

inline
void write_serial_cube_mesh_with_spider(unsigned meshSize, bool addParticleBody, const std::string & fileName)
{
  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    stk::mesh::MeshBuilder builder(MPI_COMM_SELF);
    builder.set_spatial_dimension(3);
    builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
    auto bulk = builder.create();
    stk::mesh::MetaData & meta = bulk->mesh_meta_data();

    stk::mesh::Part & block2Part = meta.declare_part_with_topology("block_2", stk::topology::BEAM_2);
    stk::mesh::Part & block3Part = meta.declare_part_with_topology("block_3", stk::topology::PARTICLE);
    stk::io::put_io_part_attribute(block2Part);
    stk::io::put_io_part_attribute(block3Part);

    unsigned newNodeId = (meshSize+1) * (meshSize+1) * (meshSize+1) + 1;
    unsigned newElemId = meshSize * meshSize * meshSize + 1;
    stk::io::fill_mesh("generated:" + std::to_string(meshSize) + "x" + std::to_string(meshSize) +
                       "x" + std::to_string(meshSize), *bulk);

    double maxCoord = meshSize;
    const auto & coordinates = *static_cast<const stk::mesh::Field<double>*>(meta.coordinate_field());

    bulk->modification_begin();

    stk::mesh::Entity spiderNode = bulk->declare_node(newNodeId);
    double * spiderCoords = stk::mesh::field_data(coordinates, spiderNode);
    spiderCoords[0] = maxCoord + 1.0;
    spiderCoords[1] = maxCoord / 2.0;
    spiderCoords[2] = maxCoord / 2.0;

    const stk::mesh::BucketVector ownedBuckets = bulk->get_buckets(stk::topology::NODE_RANK, meta.locally_owned_part());
    const stk::mesh::ConstPartVector elemParts {&block2Part};

    const auto allNodes = stk::mesh::get_entities(*bulk, stk::topology::NODE_RANK, meta.locally_owned_part());
    stk::mesh::EntityVector edgeNodes;
    for (const stk::mesh::Entity & node : allNodes) {
      const double * coords = stk::mesh::field_data(coordinates, node);
      if (std::abs(coords[0] - maxCoord) < 0.1) {
        if (bulk->num_elements(node) <= 2) {
          edgeNodes.push_back(node);
        }
      }
    }

    for (stk::mesh::Entity node : edgeNodes) {
      stk::mesh::Entity elem = bulk->declare_element(newElemId++, elemParts);
      bulk->declare_relation(elem, node, 0);
      bulk->declare_relation(elem, spiderNode, 1);
    }

    if (addParticleBody) {
      const stk::mesh::ConstPartVector bodyParts {&block3Part};
      stk::mesh::Entity spiderBody = bulk->declare_element(newElemId++, bodyParts);
      bulk->declare_relation(spiderBody, spiderNode, 0);
    }

    bulk->modification_end();

    stk::io::write_mesh(fileName, *bulk);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

#endif // UNITTESTSPIDERMESHSETUP_HPP
