// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <gtest/gtest.h>
#include "ioUtils.hpp"
#include <stddef.h>                                  // for size_t
#include <unistd.h>                                  // for unlink
#include <stk_io/StkIoUtils.hpp>
#include <string>                                    // for allocator, etc
#include <vector>                                    // for vector
#include "GeneratedMeshToFile.hpp"
#include "TextMeshToFile.hpp"
#include "Ioss_Property.h"                           // for Property
#include "mpi.h"                                     // for MPI_COMM_SELF, etc
#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"                // for StkMeshIoBroker
#include "stk_mesh/base/BulkData.hpp"                // for BulkData, etc
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/FieldBase.hpp"               // for FieldBase, etc
#include "stk_mesh/base/GetEntities.hpp"             // for get_entities
#include "stk_mesh/base/Types.hpp"                   // for EntityRank, etc
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_util/parallel/Parallel.hpp"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{

void generated_mesh_to_file_in_serial(const std::string &meshSizeSpec, const std::string &fileName)
{
  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    GeneratedMeshToFile gMesh(MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA);

    gMesh.setup_mesh(meshSizeSpec, fileName);
    gMesh.write_mesh();
  }
}

void text_mesh_to_file_in_serial(const std::string& meshDesc, const std::string& fileName)
{
  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    TextMeshToFile tMesh(MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA);

    tMesh.setup_mesh(meshDesc, fileName);
    tMesh.write_mesh();
  }
}

void generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(const std::string &meshSizeSpec, stk::mesh::BulkData &mesh, const std::string &decompositionMethod)
{
    // meshSizeSpec should NOT include generated:, just "2x2x1" for example.
    // decomposition methods: "linear", "rcb", "rib", "hsfc", "block", "cyclic", "random", "kway", "geom_kway", "metis_sfc"
    const std::string tempFilename = "exodus_" + meshSizeSpec + ".e";
    generated_mesh_to_file_in_serial(meshSizeSpec,tempFilename);

    read_from_serial_file_and_decompose(tempFilename, mesh, decompositionMethod);
    unlink(tempFilename.c_str());
}

MeshFromFile::MeshFromFile(const MPI_Comm& c)
  : m_comm(c),
    m_empty(true),
    bulk(stk::mesh::MeshBuilder(m_comm).create()),
    meta(bulk->mesh_meta_data()),
    broker()
{
}

void
MeshFromFile::fill_from_serial(const std::string& fileName)
{
  broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));
  stk::io::fill_mesh_preexisting(broker, fileName, *bulk);
  m_empty = false;
}

void
MeshFromFile::fill_from_parallel(const std::string& baseName)
{
  stk::io::fill_mesh_preexisting(broker, baseName, *bulk);
  m_empty = false;
}


TransientVerifier::TransientVerifier(const MPI_Comm& c)
  : m_comm(c),
    m_epsilon(std::numeric_limits<double>::epsilon())
{ }

void
TransientVerifier::verify_num_transient_fields(const MeshFromFile& mesh, unsigned expectedNumFields) const
{
  stk::io::FieldNameToPartVector fieldNamePartVector = mesh.broker.get_nodal_var_names();
  EXPECT_EQ(fieldNamePartVector.size(), expectedNumFields);
}

void
TransientVerifier::verify_time_steps(const MeshFromFile& mesh, const std::vector<double>& expectedTimeSteps) const
{
  EXPECT_EQ(expectedTimeSteps, mesh.broker.get_time_steps());
}

void
TransientVerifier::verify_global_variables_at_each_time_step(MeshFromFile& mesh,
                                                             const std::string& globalVariableName,
                                                             const std::vector<double>& expectedTimeSteps) const
{
  verify_time_steps(mesh, expectedTimeSteps);
  verify_global_variable_names(mesh, globalVariableName);

  int index = 0;
  for(double timeStep : expectedTimeSteps) {
    index++;
    mesh.broker.read_defined_input_fields(index);

    verify_global_int(mesh, globalVariableName, index);

    verify_global_double(mesh, globalVariableName, timeStep);
    verify_global_real_vec(mesh, globalVariableName, timeStep);
  }
}

void
TransientVerifier::verify_sideset_orientation(const MeshFromFile& mesh,
                                              int expectedProc,
                                              const stk::mesh::EntityId expectedId,
                                              const stk::mesh::ConnectivityOrdinal expectedOrdinal) const
{
  const stk::mesh::BulkData& bulk = *mesh.bulk;
  std::vector<const stk::mesh::SideSet *> sidesets = bulk.get_sidesets();

  if (stk::parallel_machine_rank(m_comm) == expectedProc) {
    ASSERT_EQ(sidesets.size(), 1u);
    ASSERT_EQ(sidesets[0]->size(), 1u);

    const stk::mesh::SideSetEntry sideSetEntry = (*sidesets[0])[0];
    EXPECT_EQ(bulk.identifier(sideSetEntry.element), expectedId);
    EXPECT_EQ(sideSetEntry.side, expectedOrdinal);
  }
  else {
    ASSERT_EQ(sidesets.size(), 1u);
    EXPECT_EQ(sidesets[0]->size(), 0u);
  }
}

void
TransientVerifier::compare_entity_rank_names(const MeshFromFile& meshA, const MeshFromFile& meshB) const
{
  EXPECT_EQ(meshA.meta.entity_rank_names(), meshB.meta.entity_rank_names());
}

void
TransientVerifier::verify_transient_field_names(const MeshFromFile& mesh, const std::string& fieldBaseName) const
{
  const stk::mesh::FieldVector transientFields = stk::io::get_transient_fields(mesh.meta, stk::topology::NODE_RANK);
  verify_transient_field_name(transientFields[0], fieldBaseName+"_scalar");
  verify_transient_field_name(transientFields[1], fieldBaseName+"_vector");
}

void
TransientVerifier::verify_transient_fields(MeshFromFile& mesh) const
{
  stk::io::StkMeshIoBroker& broker = mesh.broker;

  const stk::mesh::FieldVector transientFields = stk::io::get_transient_fields(mesh.meta);
  const std::vector<double> timeSteps = broker.get_time_steps();

  for(int iStep=0; iStep<broker.get_num_time_steps(); iStep++)
  {
    double readTime = broker.read_defined_input_fields_at_step(iStep+1, nullptr);
    EXPECT_EQ(timeSteps[iStep], readTime);

    for(stk::mesh::FieldBase* field : transientFields)
      verify_transient_field_values(*mesh.bulk, field, readTime);
  }
}

void
TransientVerifier::verify_decomp(MeshFromFile& mesh, const stk::mesh::EntityIdProcVec& expectedDecomp) const
{
  for (const stk::mesh::EntityIdProc& idProc : expectedDecomp) {
    const stk::mesh::Entity entity = mesh.bulk->get_entity(stk::topology::ELEMENT_RANK, idProc.first);
    EXPECT_EQ(mesh.bulk->parallel_owner_rank(entity), idProc.second);
  }
}

void
TransientVerifier::verify_global_variable_names(const MeshFromFile& mesh, const std::string& baseName) const
{
  std::vector<std::string> goldNames = {baseName+"_double",
                                        baseName+"_int",
                                        baseName+"_real_vec"};
  std::vector<std::string> names;
  mesh.broker.get_global_variable_names(names);
  EXPECT_EQ(names, goldNames);
}

void
TransientVerifier::verify_global_double(const MeshFromFile& mesh, const std::string& variable, double goldValue) const
{
  double value;
  EXPECT_TRUE(mesh.broker.get_global(variable+"_double", value));
  EXPECT_NEAR(value, goldValue, m_epsilon);
}

void
TransientVerifier::verify_global_int(const MeshFromFile& mesh, const std::string& variable, int goldValue) const
{
  int value;
  EXPECT_TRUE(mesh.broker.get_global(variable+"_int", value));
  EXPECT_EQ(value, goldValue);
}

void
TransientVerifier::verify_global_real_vec(const MeshFromFile& mesh, const std::string& variable, double goldValue) const
{
  std::vector<double> values;
  EXPECT_TRUE(mesh.broker.get_global(variable+"_real_vec", values));

  EXPECT_EQ(values.size(), 3u);
  for (double value : values) {
    EXPECT_NEAR(value, goldValue, m_epsilon);
  }
}

void
TransientVerifier::verify_transient_field_name(stk::mesh::FieldBase* field, const std::string& fieldName) const
{
  EXPECT_EQ(0, strcasecmp(fieldName.c_str(), field->name().c_str()))
      << fieldName << "   " << field->name();
}

void
TransientVerifier::verify_transient_field_values(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase* field, double timeStep) const
{
  const stk::mesh::BucketVector & entityBuckets = bulk.get_buckets(field->entity_rank(),bulk.mesh_meta_data().locally_owned_part());

  for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];

    for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
      stk::mesh::Entity entity = entityBucket[entityIndex];

      double * data = static_cast<double*> (stk::mesh::field_data(*field, entity));
      unsigned numEntriesPerEntity = stk::mesh::field_scalars_per_entity(*field, entity);

      for(unsigned i=0; i<numEntriesPerEntity; i++) {
        EXPECT_EQ(i + 100*timeStep + static_cast<double>(bulk.identifier(entity)), data[i]);
      }
    }
  }
}

void generated_mesh_with_transient_data_to_file_in_serial(const std::string &meshSizeSpec,
                                                          const std::string &fileName,
                                                          const std::string& fieldName,
                                                          stk::topology::rank_t fieldRank,
                                                          const std::string& globalVariableName,
                                                          const std::vector<double>& timeSteps,
                                                          const FieldValueSetter &fieldValueSetter)
{
    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
    {
        GeneratedMeshToFileWithTransientFields gMesh(MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA, fieldName,
                                                     fieldRank);

        gMesh.setup_mesh(meshSizeSpec, fileName);
        gMesh.write_mesh_with_field(timeSteps, fieldValueSetter, globalVariableName);
    }
}

void read_from_serial_file_and_decompose(const std::string& fileName, stk::mesh::BulkData &mesh, const std::string &decompositionMethod)
{
    stk::io::StkMeshIoBroker broker;
    broker.set_bulk_data(mesh);
    broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", decompositionMethod));
    broker.add_mesh_database(fileName, stk::io::READ_MESH);
    broker.create_input_mesh();
    broker.populate_bulk_data();
}

void IdAndTimeFieldValueSetter::populate_field(stk::mesh::BulkData &bulk, stk::mesh::FieldBase* field, const unsigned step, const double time) const
{
    stk::mesh::EntityRank fieldRank = field->entity_rank();

    std::vector<stk::mesh::Entity> entities;
    stk::mesh::get_entities(bulk, fieldRank, entities);

    stk::mesh::FieldVector allTransientFields = stk::io::get_transient_fields(bulk.mesh_meta_data());

    for(stk::mesh::FieldBase * transientField : allTransientFields)
    {
        for(size_t i = 0; i < entities.size(); i++)
        {
            unsigned numEntriesPerEntity = stk::mesh::field_scalars_per_entity(*transientField, entities[i]);
            double value = 100.0 * time + static_cast<double>(bulk.identifier(entities[i]));
            double *data = static_cast<double*> (stk::mesh::field_data(*transientField, entities[i]));
            for(unsigned j=0; j<numEntriesPerEntity; j++)
                data[j] = value + j;
        }
    }
}

namespace simple_fields {

void generated_mesh_to_file_in_serial(const std::string &meshSizeSpec, const std::string &fileName)
{
  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    stk::unit_test_util::GeneratedMeshToFile gMesh(MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA);

    gMesh.setup_mesh(meshSizeSpec, fileName);
    gMesh.write_mesh();
  }
}

void text_mesh_to_file_in_serial(const std::string& meshDesc, const std::string& fileName)
{
  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    stk::unit_test_util::TextMeshToFile tMesh(MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA);

    tMesh.setup_mesh(meshDesc, fileName);
    tMesh.write_mesh();
  }
}

void generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(const std::string &meshSizeSpec, stk::mesh::BulkData &mesh, const std::string &decompositionMethod)
{
    // meshSizeSpec should NOT include generated:, just "2x2x1" for example.
    // decomposition methods: "linear", "rcb", "rib", "hsfc", "block", "cyclic", "random", "kway", "geom_kway", "metis_sfc"
    const std::string tempFilename = "exodus_" + meshSizeSpec + ".e";
    stk::unit_test_util::generated_mesh_to_file_in_serial(meshSizeSpec,tempFilename);

    stk::unit_test_util::read_from_serial_file_and_decompose(tempFilename, mesh, decompositionMethod);
    unlink(tempFilename.c_str());
}

MeshFromFile::MeshFromFile(const MPI_Comm& c)
  : m_comm(c),
    m_empty(true),
    bulk(stk::mesh::MeshBuilder(m_comm).create()),
    meta(bulk->mesh_meta_data()),
    broker()
{
}

void
MeshFromFile::fill_from_serial(const std::string& fileName)
{
  broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));
  stk::io::fill_mesh_preexisting(broker, fileName, *bulk);
  m_empty = false;
}

void
MeshFromFile::fill_from_parallel(const std::string& baseName)
{
  stk::io::fill_mesh_preexisting(broker, baseName, *bulk);
  m_empty = false;
}


TransientVerifier::TransientVerifier(const MPI_Comm& c)
  : m_comm(c),
    m_epsilon(std::numeric_limits<double>::epsilon())
{ }

void
TransientVerifier::verify_num_transient_fields(const MeshFromFile& mesh, unsigned expectedNumFields) const
{
  stk::io::FieldNameToPartVector fieldNamePartVector = mesh.broker.get_nodal_var_names();
  EXPECT_EQ(fieldNamePartVector.size(), expectedNumFields);
}

void
TransientVerifier::verify_time_steps(const MeshFromFile& mesh, const std::vector<double>& expectedTimeSteps) const
{
  EXPECT_EQ(expectedTimeSteps, mesh.broker.get_time_steps());
}

void
TransientVerifier::verify_global_variables_at_each_time_step(MeshFromFile& mesh,
                                                             const std::string& globalVariableName,
                                                             const std::vector<double>& expectedTimeSteps) const
{
  verify_time_steps(mesh, expectedTimeSteps);
  verify_global_variable_names(mesh, globalVariableName);

  int index = 0;
  for(double timeStep : expectedTimeSteps) {
    index++;
    mesh.broker.read_defined_input_fields(index);

    verify_global_int(mesh, globalVariableName, index);

    verify_global_double(mesh, globalVariableName, timeStep);
    verify_global_real_vec(mesh, globalVariableName, timeStep);
  }
}

void
TransientVerifier::verify_sideset_orientation(const MeshFromFile& mesh,
                                              int expectedProc,
                                              const stk::mesh::EntityId expectedId,
                                              const stk::mesh::ConnectivityOrdinal expectedOrdinal) const
{
  const stk::mesh::BulkData& bulk = *mesh.bulk;
  std::vector<const stk::mesh::SideSet *> sidesets = bulk.get_sidesets();

  if (stk::parallel_machine_rank(m_comm) == expectedProc) {
    ASSERT_EQ(sidesets.size(), 1u);
    ASSERT_EQ(sidesets[0]->size(), 1u);

    const stk::mesh::SideSetEntry sideSetEntry = (*sidesets[0])[0];
    EXPECT_EQ(bulk.identifier(sideSetEntry.element), expectedId);
    EXPECT_EQ(sideSetEntry.side, expectedOrdinal);
  }
  else {
    ASSERT_EQ(sidesets.size(), 1u);
    EXPECT_EQ(sidesets[0]->size(), 0u);
  }
}

void
TransientVerifier::compare_entity_rank_names(const MeshFromFile& meshA, const MeshFromFile& meshB) const
{
  EXPECT_EQ(meshA.meta.entity_rank_names(), meshB.meta.entity_rank_names());
}

void
TransientVerifier::verify_transient_field_names(const MeshFromFile& mesh, const std::string& fieldBaseName) const
{
  const stk::mesh::FieldVector transientFields = stk::io::get_transient_fields(mesh.meta, stk::topology::NODE_RANK);
  verify_transient_field_name(transientFields[0], fieldBaseName+"_scalar");
  verify_transient_field_name(transientFields[1], fieldBaseName+"_vector");
}

void
TransientVerifier::verify_transient_fields(MeshFromFile& mesh) const
{
  stk::io::StkMeshIoBroker& broker = mesh.broker;

  const stk::mesh::FieldVector transientFields = stk::io::get_transient_fields(mesh.meta);
  const std::vector<double> timeSteps = broker.get_time_steps();

  for(int iStep=0; iStep<broker.get_num_time_steps(); iStep++)
  {
    double readTime = broker.read_defined_input_fields_at_step(iStep+1, nullptr);
    EXPECT_EQ(timeSteps[iStep], readTime);

    for(stk::mesh::FieldBase* field : transientFields)
      verify_transient_field_values(*mesh.bulk, field, readTime);
  }
}

void
TransientVerifier::verify_decomp(MeshFromFile& mesh, const stk::mesh::EntityIdProcVec& expectedDecomp) const
{
  for (const stk::mesh::EntityIdProc& idProc : expectedDecomp) {
    const stk::mesh::Entity entity = mesh.bulk->get_entity(stk::topology::ELEMENT_RANK, idProc.first);
    EXPECT_EQ(mesh.bulk->parallel_owner_rank(entity), idProc.second);
  }
}

void
TransientVerifier::verify_global_variable_names(const MeshFromFile& mesh, const std::string& baseName) const
{
  std::vector<std::string> goldNames = {baseName+"_double",
                                        baseName+"_int",
                                        baseName+"_real_vec"};
  std::vector<std::string> names;
  mesh.broker.get_global_variable_names(names);
  EXPECT_EQ(names, goldNames);
}

void
TransientVerifier::verify_global_double(const MeshFromFile& mesh, const std::string& variable, double goldValue) const
{
  double value;
  EXPECT_TRUE(mesh.broker.get_global(variable+"_double", value));
  EXPECT_NEAR(value, goldValue, m_epsilon);
}

void
TransientVerifier::verify_global_int(const MeshFromFile& mesh, const std::string& variable, int goldValue) const
{
  int value;
  EXPECT_TRUE(mesh.broker.get_global(variable+"_int", value));
  EXPECT_EQ(value, goldValue);
}

void
TransientVerifier::verify_global_real_vec(const MeshFromFile& mesh, const std::string& variable, double goldValue) const
{
  std::vector<double> values;
  EXPECT_TRUE(mesh.broker.get_global(variable+"_real_vec", values));

  EXPECT_EQ(values.size(), 3u);
  for (double value : values) {
    EXPECT_NEAR(value, goldValue, m_epsilon);
  }
}

void
TransientVerifier::verify_transient_field_name(stk::mesh::FieldBase* field, const std::string& fieldName) const
{
  EXPECT_EQ(0, strcasecmp(fieldName.c_str(), field->name().c_str()))
      << fieldName << "   " << field->name();
}

void
TransientVerifier::verify_transient_field_values(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase* field, double timeStep) const
{
  const stk::mesh::BucketVector & entityBuckets = bulk.get_buckets(field->entity_rank(),bulk.mesh_meta_data().locally_owned_part());

  for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
    stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];

    for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
      stk::mesh::Entity entity = entityBucket[entityIndex];

      double * data = static_cast<double*> (stk::mesh::field_data(*field, entity));
      unsigned numEntriesPerEntity = stk::mesh::field_scalars_per_entity(*field, entity);

      for(unsigned i=0; i<numEntriesPerEntity; i++) {
        EXPECT_EQ(i + 100*timeStep + static_cast<double>(bulk.identifier(entity)), data[i]);
      }
    }
  }
}

void generated_mesh_with_transient_data_to_file_in_serial(const std::string &meshSizeSpec,
                                                          const std::string &fileName,
                                                          const std::string& fieldName,
                                                          stk::topology::rank_t fieldRank,
                                                          const std::string& globalVariableName,
                                                          const std::vector<double>& timeSteps,
                                                          const FieldValueSetter &fieldValueSetter)
{
    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
    {
        stk::unit_test_util::GeneratedMeshToFileWithTransientFields gMesh(MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA, fieldName,
                                                     fieldRank);

        gMesh.setup_mesh(meshSizeSpec, fileName);
        gMesh.write_mesh_with_field(timeSteps, fieldValueSetter, globalVariableName);
    }
}

void read_from_serial_file_and_decompose(const std::string& fileName, stk::mesh::BulkData &mesh, const std::string &decompositionMethod)
{
    stk::io::StkMeshIoBroker broker;
    broker.set_bulk_data(mesh);
    broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", decompositionMethod));
    broker.add_mesh_database(fileName, stk::io::READ_MESH);
    broker.create_input_mesh();
    broker.populate_bulk_data();
}

} // namespace simple_fields

} // namespace unit_test_util
} // namespace stk

