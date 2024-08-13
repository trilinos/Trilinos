// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "GeneratedMeshToFile.hpp"
#include <iosfwd>                               // for ostream
#include <string>                               // for operator+, string
#include <vector>                               // for vector
#include "Ioss_Field.h"                         // for Field, etc
#include "mpi.h"                                // for ompi_communicator_t
#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"           // for StkMeshIoBroker
#include "stk_mesh/base/BulkData.hpp"           // for BulkData, etc
#include "stk_mesh/base/CoordinateSystems.hpp"  // for Cartesian
#include "stk_mesh/base/Field.hpp"              // for Field
#include "stk_mesh/base/MetaData.hpp"           // for MetaData, put_field
#include "stk_unit_test_utils/ioUtils.hpp"      // for FieldValueSetter
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{

GeneratedMeshToFile::GeneratedMeshToFile(stk::ParallelMachine comm,
                                         stk::mesh::BulkData::AutomaticAuraOption auraOption)
  : bulkPtr(build_mesh(3, comm, auraOption)),
    bulk(*bulkPtr),
    meta(bulk.mesh_meta_data())
{
}

void GeneratedMeshToFile::write_mesh()
{
  broker.write_output_mesh(outputFileIndex);
}

void GeneratedMeshToFile::setup_mesh(const std::string &meshSizeSpec,
                                     const std::string &outputFileName,
                                     bool useBigIds)
{
  broker.set_bulk_data(bulk);
  broker.add_mesh_database("generated:" + meshSizeSpec, stk::io::READ_MESH);
  broker.create_input_mesh();
  broker.populate_bulk_data();
  if (useBigIds) {
    broker.property_add(Ioss::Property("INTEGER_SIZE_API", 8));
    broker.property_add(Ioss::Property("INTEGER_SIZE_DB", 8));
  }
  outputFileIndex = broker.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);
}

GeneratedMeshToFileWithTransientFields::GeneratedMeshToFileWithTransientFields(stk::ParallelMachine comm,
                                                                               stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                                                               const std::string& fieldBaseName,
                                                                               stk::topology::rank_t rank)
  : GeneratedMeshToFile(comm, auraOption),
    fieldRank(rank),
    scalarField(meta.declare_field<double>(fieldRank, fieldBaseName+"_scalar", 1)),
    vectorField(meta.declare_field<double>(fieldRank, fieldBaseName+"_vector", 1))
{
  stk::mesh::put_field_on_mesh(scalarField, meta.universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(vectorField, meta.universal_part(), 3, nullptr);

  stk::io::set_field_output_type(vectorField, stk::io::FieldOutputType::VECTOR_3D);
}

void GeneratedMeshToFileWithTransientFields::write_mesh_with_field(const std::vector<double>& timeSteps,
                                                                   const FieldValueSetter &fieldValueSetter,
                                                                   const std::string& globalVariableName)
{
  broker.add_field(outputFileIndex, scalarField);
  broker.add_field(outputFileIndex, vectorField);
  broker.add_global(outputFileIndex, globalVariableName+"_double", Ioss::Field::DOUBLE);
  broker.add_global(outputFileIndex, globalVariableName+"_int", Ioss::Field::INTEGER);
  int vectorLength = 3;
  broker.add_global(outputFileIndex, globalVariableName+"_real_vec", vectorLength, Ioss::Field::REAL);

  for(unsigned step = 0; step < timeSteps.size(); step++)
  {
    double time = timeSteps[step];
    int timestep = step+1;

    fieldValueSetter.populate_field(bulk, &scalarField, step, time);
    fieldValueSetter.populate_field(bulk, &vectorField, step, time);
    broker.begin_output_step(outputFileIndex, time);
    broker.write_defined_output_fields(outputFileIndex);
    broker.write_global(outputFileIndex, globalVariableName+"_double", timeSteps[step]);
    broker.write_global(outputFileIndex, globalVariableName+"_int", timestep);
    std::vector<double> values(vectorLength, timeSteps[step]);
    broker.write_global(outputFileIndex, globalVariableName+"_real_vec", values);
    broker.end_output_step(outputFileIndex);
  }
}

namespace simple_fields {

GeneratedMeshToFile::GeneratedMeshToFile(stk::ParallelMachine comm,
                                         stk::mesh::BulkData::AutomaticAuraOption auraOption)
  : bulkPtr(build_mesh(3, comm, auraOption)),
    bulk(*bulkPtr),
    meta(bulk.mesh_meta_data())
{
}

void GeneratedMeshToFile::write_mesh()
{
  broker.write_output_mesh(outputFileIndex);
}

void GeneratedMeshToFile::setup_mesh(const std::string &meshSizeSpec,
                                     const std::string &outputFileName,
                                     bool useBigIds)
{
  broker.set_bulk_data(bulk);
  broker.add_mesh_database("generated:" + meshSizeSpec, stk::io::READ_MESH);
  broker.create_input_mesh();
  broker.populate_bulk_data();
  if (useBigIds) {
    broker.property_add(Ioss::Property("INTEGER_SIZE_API", 8));
    broker.property_add(Ioss::Property("INTEGER_SIZE_DB", 8));
  }
  outputFileIndex = broker.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);
}

GeneratedMeshToFileWithTransientFields::GeneratedMeshToFileWithTransientFields(stk::ParallelMachine comm,
                                                                               stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                                                               const std::string& fieldBaseName,
                                                                               stk::topology::rank_t rank)
  : GeneratedMeshToFile(comm, auraOption),
    fieldRank(rank),
    scalarField(meta.declare_field<double>(fieldRank, fieldBaseName+"_scalar", 1)),
    vectorField(meta.declare_field<double>(fieldRank, fieldBaseName+"_vector", 1))
{
  stk::mesh::put_field_on_mesh(scalarField, meta.universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(vectorField, meta.universal_part(), 3, nullptr);

  stk::io::set_field_output_type(vectorField, stk::io::FieldOutputType::VECTOR_3D);
}

void GeneratedMeshToFileWithTransientFields::write_mesh_with_field(const std::vector<double>& timeSteps,
                                                                   const FieldValueSetter &fieldValueSetter,
                                                                   const std::string& globalVariableName)
{
  broker.add_field(outputFileIndex, scalarField);
  broker.add_field(outputFileIndex, vectorField);
  broker.add_global(outputFileIndex, globalVariableName+"_double", Ioss::Field::DOUBLE);
  broker.add_global(outputFileIndex, globalVariableName+"_int", Ioss::Field::INTEGER);
  int vectorLength = 3;
  broker.add_global(outputFileIndex, globalVariableName+"_real_vec", vectorLength, Ioss::Field::REAL);

  for(unsigned step = 0; step < timeSteps.size(); step++)
  {
    double time = timeSteps[step];
    int timestep = step+1;

    fieldValueSetter.populate_field(bulk, &scalarField, step, time);
    fieldValueSetter.populate_field(bulk, &vectorField, step, time);
    broker.begin_output_step(outputFileIndex, time);
    broker.write_defined_output_fields(outputFileIndex);
    broker.write_global(outputFileIndex, globalVariableName+"_double", timeSteps[step]);
    broker.write_global(outputFileIndex, globalVariableName+"_int", timestep);
    std::vector<double> values(vectorLength, timeSteps[step]);
    broker.write_global(outputFileIndex, globalVariableName+"_real_vec", values);
    broker.end_output_step(outputFileIndex);
  }
}

} // namespace simple_fields

}
}
