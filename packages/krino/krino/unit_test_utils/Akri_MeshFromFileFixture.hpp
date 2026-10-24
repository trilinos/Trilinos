#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_MESHFROMFILEFIXTURE_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_MESHFROMFILEFIXTURE_HPP_
#include <gtest/gtest.h>
#include <memory>

#include <Akri_AuxMetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshFromFile.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_OutputUtils.hpp>

namespace krino {

class MeshFromFileFixture : public ::testing::Test
{
public:
  bool read_mesh_if_present_and_supported(const std::string & inputMeshName)
  {
    if ((is_parallel_io_enabled() || stk::parallel_machine_size(myComm) == 1) && does_file_exist(inputMeshName))
    {
      myMeshFromFile = std::make_unique<MeshFromFile>(inputMeshName, myComm, "rib");
      myMeshFromFile->populate_mesh();
      activate_all_entities(myMeshFromFile->bulk_data(), get_aux_meta().active_part());
      return true;
    }
    return false;
  }
  bool is_test_active() const { return myMeshFromFile.get() != nullptr; }
  stk::mesh::BulkData & get_mesh() { return myMeshFromFile->bulk_data(); }
  const stk::mesh::BulkData & get_mesh() const { return myMeshFromFile->bulk_data(); }
  stk::ParallelMachine get_comm() const { return myComm; }
  int parallel_size() const { return stk::parallel_machine_size(myComm); }
  int parallel_rank() const { return stk::parallel_machine_rank(myComm); }
  AuxMetaData & get_aux_meta() { return AuxMetaData::get(myMeshFromFile->meta_data()); }
  const AuxMetaData & get_aux_meta() const { return AuxMetaData::get(myMeshFromFile->meta_data()); }
  CoordinatesFieldRef get_coordinates_field() { const CoordinatesFieldRef coordsField(myMeshFromFile->meta_data().coordinate_field(), myMeshFromFile->meta_data().spatial_dimension()); return coordsField; }

protected:

  stk::ParallelMachine myComm{MPI_COMM_WORLD};
  std::unique_ptr<MeshFromFile> myMeshFromFile;
};

}



#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_MESHFROMFILEFIXTURE_HPP_ */
