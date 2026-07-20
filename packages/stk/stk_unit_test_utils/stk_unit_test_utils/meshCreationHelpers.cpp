#include <stk_unit_test_utils/meshCreationHelpers.hpp>
#include <stddef.h>                    // for size_t
#include <gtest/gtest.h>
#include <stk_io/DatabasePurpose.hpp>  // for DatabasePurpose::READ_MESH, etc
#include <stk_io/StkMeshIoBroker.hpp>  // for StkMeshIoBroker
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/Field.hpp>     // for Field
#include <stk_mesh/base/MetaData.hpp>  // for MetaData, put_field
#include <stk_topology/topology.hpp>   // for topology, etc
#include "stk_mesh/base/BulkData.hpp"  // for BulkData
#include "stk_mesh/base/ForEachEntity.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"

#include <limits>

namespace stk
{
namespace unit_test_util
{

size_t write_mesh_data__field_1__field_2__field_3(const std::string & filename, MPI_Comm /*communicator*/, stk::mesh::BulkData & bulk, stk::io::StkMeshIoBroker & stkIo)
{
    stkIo.set_bulk_data(bulk);
    size_t index = stkIo.add_mesh_database("generated:1x1x4", stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &meta = bulk.mesh_meta_data();
    stk::mesh::Field<double> & field1 = meta.declare_field<double>(stk::topology::ELEM_RANK, "field_1", 1);
    stk::mesh::Field<double> & field2 = meta.declare_field<double>(stk::topology::ELEM_RANK, "field_2", 1);
    stk::mesh::Field<double> & field3 = meta.declare_field<double>(stk::topology::ELEM_RANK, "field_3", 1);

    double fieldValues[] = {1.0, 2.0, 3.0};
    stk::mesh::put_field_on_mesh(field1, meta.universal_part(), fieldValues);
    stk::mesh::put_field_on_mesh(field2, meta.universal_part(), fieldValues + 1);
    stk::mesh::put_field_on_mesh(field3, meta.universal_part(), fieldValues + 2);
    stkIo.populate_bulk_data();

    size_t results_output_index = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    stkIo.add_field(results_output_index, field1);
    stkIo.add_field(results_output_index, field2);
    stkIo.add_field(results_output_index, field3);

    stkIo.write_output_mesh(results_output_index);
    return results_output_index;
}

void create_mesh_without_time_steps(const std::string & filename, MPI_Comm communicator)
{
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(communicator, stk::mesh::BulkData::AUTO_AURA);
    stk::io::StkMeshIoBroker stkIoWriter(communicator);
    write_mesh_data__field_1__field_2__field_3(filename, communicator, *bulk, stkIoWriter);
}

void create_mesh_with__field_1__field_2__field_3(const std::string & filename, MPI_Comm communicator)
{
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(communicator, stk::mesh::BulkData::AUTO_AURA);
    stk::io::StkMeshIoBroker stkIo(communicator);
    size_t results_output_index = write_mesh_data__field_1__field_2__field_3(filename, communicator, *bulk, stkIo);
    double time = 0.0;
    stkIo.begin_output_step(results_output_index, time);
    stkIo.write_defined_output_fields(results_output_index);
    stkIo.end_output_step(results_output_index);
}

void scale_to_unit_bbox(stk::mesh::BulkData& mesh)
{
  double maxX = 0.0, maxY = 0.0, maxZ = 0.0;

  auto* coordField = mesh.mesh_meta_data().coordinate_field();
  STK_ThrowRequireMsg(coordField != nullptr, "Error, coordinate_field is nullptr");

  const unsigned spatialDim = mesh.mesh_meta_data().spatial_dimension();
  stk::mesh::field_data_execute<double,stk::mesh::ReadOnly>(*coordField,
  [&](auto& coordData) {
    stk::mesh::for_each_entity_run(mesh, stk::topology::NODE_RANK, mesh.mesh_meta_data().universal_part(),
      [&](const stk::mesh::BulkData&, stk::mesh::Entity node) {
        auto coord = coordData.entity_values(node);
        const bool nonNegative = coord(0_comp) >= 0.0 && coord(1_comp) >= 0.0 &&
          (spatialDim==2 || coord(2_comp) >= 0.0);
        STK_ThrowRequireMsg(nonNegative, "Currently only supporting non-negative coords.");

        maxX = std::max(maxX, coord(0_comp));
        maxY = std::max(maxY, coord(1_comp));
        if (spatialDim == 3) {
          maxZ = std::max(maxZ, coord(2_comp));
        }
      });
  });

  const double xfactor = 1.0/maxX;
  const double yfactor = 1.0/maxY;
  const double zfactor = spatialDim==3 ? 1.0/maxZ : 1.0;

  stk::mesh::field_data_execute<double,stk::mesh::ReadWrite>(*coordField,
  [&](auto& coordData) {
    stk::mesh::for_each_entity_run(mesh, stk::topology::NODE_RANK, mesh.mesh_meta_data().universal_part(),
      [&](const stk::mesh::BulkData&, stk::mesh::Entity node) {
        auto coord = coordData.entity_values(node);
        coord(0_comp) *= xfactor;
        coord(1_comp) *= yfactor;
        if (spatialDim == 3) {
          coord(2_comp) *= zfactor;
        }
      });
  });
}

void add_field_to_mesh(stk::mesh::BulkData& mesh,
                       const std::string& fieldName,
                       const FieldEvaluator& fieldEval,
                       const stk::mesh::EntityRank fieldRank,
                       const Ioss::Field::RoleType fieldRole)
{
  STK_ThrowRequireMsg(fieldRank==stk::topology::NODE_RANK,"Currently add_field_to_mesh, taking FieldEvaluator, only supports nodal fields.");
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  meta.enable_late_fields();
  stk::mesh::Field<double> &field = meta.declare_field<double>(fieldRank, fieldName, 1);
  stk::io::set_field_role(field, fieldRole);
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), nullptr);
  set_node_field(mesh, field, fieldEval);
}

std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_field(const stk::ParallelMachine comm,
                       const std::string& fieldName,
                       const double fieldValue,
                       const stk::mesh::EntityRank fieldRank,
                       const Ioss::Field::RoleType fieldRole,
                       const std::string generatedMeshString,
                       const int numCopies,
                       const int spatialDim)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(spatialDim, comm);

  stk::mesh::Field<double> &field=
          bulk->mesh_meta_data().declare_field<double>(fieldRank, fieldName, 1);
  stk::io::set_field_role(field, fieldRole);
  std::vector<double> fieldValues(numCopies, fieldValue);
  stk::mesh::put_field_on_mesh(field, bulk->mesh_meta_data().universal_part(), 1, numCopies, fieldValues.data());
  stk::io::fill_mesh(generatedMeshString, *bulk);
  stk::mesh::Part& edgePart = bulk->mesh_meta_data().declare_part("edges_1", stk::topology::EDGE_RANK);
  stk::mesh::create_edges(*bulk, bulk->mesh_meta_data().universal_part(), &edgePart);

  if (bulk->mesh_meta_data().get_part("nodeset_1") != nullptr) {
    EXPECT_EQ(bulk->mesh_meta_data().get_part("nodeset_1")->primary_entity_rank(), stk::topology::NODE_RANK);
  }
  if (bulk->mesh_meta_data().get_part("edges_1") != nullptr) {
    EXPECT_EQ(bulk->mesh_meta_data().get_part("edges_1")->primary_entity_rank(), stk::topology::EDGE_RANK);
  }
  if (bulk->mesh_meta_data().get_part("surface_1") != nullptr) {
    EXPECT_EQ(bulk->mesh_meta_data().get_part("surface_1")->primary_entity_rank(), stk::topology::FACE_RANK);
  }
  if (bulk->mesh_meta_data().get_part("block_1") != nullptr) {
    EXPECT_EQ(bulk->mesh_meta_data().get_part("block_1")->primary_entity_rank(), stk::topology::ELEM_RANK);
  }

  return bulk;
}

std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_field(const stk::ParallelMachine comm,
                       const std::string& fieldName,
                       const FieldEvaluator& fieldEval,
                       const stk::mesh::EntityRank fieldRank,
                       const Ioss::Field::RoleType fieldRole,
                       const std::string generatedMeshString,
                       const int numCopies,
                       const int spatialDim)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = create_mesh_with_field(comm, fieldName, 0.0, fieldRank, fieldRole, generatedMeshString, numCopies, spatialDim);
  STK_ThrowRequireMsg(fieldRank==stk::topology::NODE_RANK,"Currently create_mesh_with_field, taking FieldEvaluator, only supports nodal fields.");
  const stk::mesh::FieldBase* field = bulk->mesh_meta_data().get_field(stk::topology::NODE_RANK, fieldName);
  STK_ThrowRequireMsg(field != nullptr, "Failed to retrieve field '"<<fieldName<<"'.");
  set_node_field(*bulk, *field, fieldEval);
  return bulk;
}

}
}
