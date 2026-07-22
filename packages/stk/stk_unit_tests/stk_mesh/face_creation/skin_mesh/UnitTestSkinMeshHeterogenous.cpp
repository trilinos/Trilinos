#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

#include "SetupKeyholeMesh.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/QuadFixture.hpp"  // for QuadFixture
#include "stk_unit_test_utils/stk_mesh_fixtures/degenerate_mesh.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/heterogeneous_mesh.hpp"

#include "UnitTestSkinMeshUseCaseUtils.hpp"

namespace {

using namespace stk::mesh::impl;
using namespace stk::mesh;
using stk::unit_test_util::build_mesh;

TEST(ElementGraph, heterogeneous_mesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) <= 17)
  {
    std::string fileName("hetero.g");
    {
      std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA);
      stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
      stk::mesh::BulkData& bulk_data = *bulkPtr;
      stk::mesh::fixtures::VectorFieldType & node_coord =
          meta_data.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
      stk::mesh::put_field_on_mesh( node_coord , meta_data.universal_part() , 3, nullptr);
      stk::io::set_field_output_type(node_coord, stk::io::FieldOutputType::VECTOR_3D);

      stk::mesh::fixtures::heterogeneous_mesh_meta_data( meta_data , node_coord );
      meta_data.commit();

      stk::mesh::fixtures::heterogeneous_mesh_bulk_data( bulk_data , node_coord );
      if (stk::parallel_machine_rank(comm) == 0)
      {
        stk::io::write_mesh(fileName, bulk_data);
      }
    }

    stk::parallel_machine_barrier(comm);

    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, comm, stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& bulk_data = *bulkPtr;

    for(int i=1;i<=17;i++)
    {
      std::ostringstream os;
      os << "skin_" << i;
      std::string part_name = os.str();
      stk::mesh::Part &tmp = meta_data.declare_part(part_name, meta_data.side_rank());
      stk::io::put_io_part_attribute(tmp);
    }

    stk::mesh::Part& skin = meta_data.declare_part("skin", meta_data.side_rank());

    stk::unit_test_util::read_from_serial_file_and_decompose(fileName, bulk_data, "RIB");
    unlink(fileName.c_str());
    EXPECT_NO_FATAL_FAILURE(ElemGraphTestUtils::skin_boundary(bulk_data, meta_data.locally_owned_part(), {&skin}));
    std::vector<size_t> mesh_counts;
    stk::mesh::comm_mesh_counts(bulk_data, mesh_counts);
    EXPECT_EQ(23u, mesh_counts[meta_data.side_rank()]);

    std::vector<std::pair<stk::mesh::EntityId, int>> id_and_num_faces = {
      {7, 1},
      {8, 1},
      {9, 1},
      {10, 2},
      {11, 1},
      {4, 2},
      {5, 2},
      {6, 2},
      {1, 3},
      {2, 2},
      {3, 4},
      {15, 0},
      {16, 0},
      {17, 1},
      {12, 0},
      {13, 0},
      {14, 1}
    };

    bulk_data.modification_begin();

    for(size_t i=0;i<id_and_num_faces.size();++i)
    {
      stk::mesh::EntityId id = id_and_num_faces[i].first;
      int gold_num_faces = id_and_num_faces[i].second;
      stk::mesh::Entity elem = bulk_data.get_entity(stk::topology::ELEM_RANK, id);
      if(bulk_data.is_valid(elem))
      {
        int num_faces = bulk_data.num_sides(elem);
        const stk::mesh::Entity *faces = bulk_data.begin_faces(elem);
        for(int j=0;j<num_faces;++j)
        {
          std::ostringstream os;
          os << "skin_" << bulk_data.identifier(elem);
          stk::mesh::PartVector add_parts;
          add_parts.push_back(meta_data.get_part(os.str()));
          bulk_data.change_entity_parts(faces[j], add_parts, {});

        }
        EXPECT_EQ(gold_num_faces, num_faces) << "element " << id << " has topology " << bulk_data.bucket(elem).topology() << " with num faces " << num_faces << " not same as gold value " << gold_num_faces << std::endl;
      }
    }

    bulk_data.modification_begin();

    //stk::io::write_mesh("heter.g", bulk_data, comm);
  }
}

}
