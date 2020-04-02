#include <stk_unit_test_utils/unittestMeshUtils.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "stk_unit_test_utils/GenerateALefRAMesh.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_io/WriteMesh.hpp"
#include <stk_io/StkMeshIoBroker.hpp>

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_mesh/base/FieldParallel.hpp"
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include <iostream>
#include <limits>

class StkBalanceRunner
{
public:
  StkBalanceRunner(MPI_Comm c)
    : m_comm(c)
  { }

  void set_filename(const std::string& name)
  {
    m_options.m_inFile = name;
  }

  void set_output_dir(const std::string& name)
  {
    m_options.outputDirectory = name;
  }

  void set_app_type_defaults(stk::balance::AppTypeDefaults defaults)
  {
    m_options.set_app_type_default(defaults);
  }

  void set_initial_decomp_method(const std::string& method)
  {
    m_options.set_initial_decomp_method(method);
  }

  void set_decomp_method(const std::string& method)
  {
    m_options.set_decomp_method(method);
  }

  void run_end_to_end() const
  {
    stk::balance::run_stk_rebalance(m_options, m_comm);
  }

private:
  MPI_Comm m_comm;
  stk::balance::ParsedOptions m_options;
};

class MeshFromFile
{
public:
  MeshFromFile(const MPI_Comm& c)
    : m_comm(c),
      m_empty(true),
      bulk(meta, m_comm)
  { }

  void fill_from_serial(const std::string& fileName)
  {
    stk::io::fill_mesh_with_auto_decomp(fileName, bulk, broker);
    m_empty = false;
  }

  void fill_from_parallel(const std::string& baseName)
  {
    stk::io::fill_mesh_preexisting(broker, baseName, bulk);
    m_empty = false;
  }

  bool is_empty() const { return m_empty; }

private:
  MPI_Comm m_comm;
  bool m_empty;

public:
  stk::mesh::MetaData meta;
  stk::mesh::BulkData bulk;
  stk::io::StkMeshIoBroker broker;
};

class TransientWriter
{
public:
  TransientWriter(MPI_Comm c, const std::string& name)
    : m_comm(c),
      m_fileBaseName(name),
      m_fieldName("myFieldName"),
      m_varName("TestTime")
  { }

  void set_field_name(const std::string& name)
  {
    m_fieldName = name;
  }

  void set_global_variable_name(const std::string& name)
  {
    m_varName = name;
  }

  void set_time_steps(const std::vector<double>& steps)
  {
    m_timeSteps = steps;
  }

  void write_static_mesh(const std::string& meshDesc) const
  {
    stk::unit_test_util::generated_mesh_to_file_in_serial(meshDesc, m_fileBaseName);
  }

  void write_transient_mesh(const std::string& meshDesc) const
  {
    stk::unit_test_util::generated_mesh_with_transient_data_to_file_in_serial(
          meshDesc, m_fileBaseName, m_fieldName, m_varName, m_timeSteps, m_fieldSetter);
  }

  void write_two_element_mesh_with_sideset(stk::unit_test_util::ElementOrdering elemOrdering) const
  {
    if (stk::parallel_machine_rank(m_comm) == 0) {
      stk::mesh::MetaData meta(3);
      stk::mesh::BulkData bulk(meta, MPI_COMM_SELF);

      stk::unit_test_util::create_AB_mesh_with_sideset_and_field(
            bulk, stk::unit_test_util::LEFT, elemOrdering, "dummyField");
      stk::io::write_mesh_with_fields(m_fileBaseName, bulk, 1, 1.0);
    }
  }

private:
  MPI_Comm m_comm;
  std::string m_fileBaseName;
  std::string m_fieldName;
  std::string m_varName;
  std::vector<double> m_timeSteps;
  stk::unit_test_util::IdAndTimeFieldValueSetter m_fieldSetter;
};

class TransientVerifier
{
public:
  TransientVerifier(const MPI_Comm& c)
    : m_comm(c),
      m_epsilon(std::numeric_limits<double>::epsilon())
  { }

  void verify_num_transient_fields(const MeshFromFile& mesh, unsigned expectedNumFields) const
  {
    stk::io::FieldNameToPartVector fieldNamePartVector = mesh.broker.get_nodal_var_names();
    EXPECT_EQ(fieldNamePartVector.size(), expectedNumFields);
  }

  void verify_time_steps(const MeshFromFile& mesh, const std::vector<double>& expectedTimeSteps) const
  {
    EXPECT_EQ(expectedTimeSteps, mesh.broker.get_time_steps());
  }

  void verify_global_variables_at_each_time_step(MeshFromFile& mesh,
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

  void verify_sideset_orientation(const MeshFromFile& mesh,
                                  const stk::mesh::EntityId expectedId,
                                  const stk::mesh::ConnectivityOrdinal expectedOrdinal) const
  {
    std::vector<const stk::mesh::SideSet *> sidesets = mesh.bulk.get_sidesets();

    if (stk::parallel_machine_rank(m_comm) == 0) {
      ASSERT_EQ(sidesets.size(), 1u);
      ASSERT_EQ(sidesets[0]->size(), 1u);

      const stk::mesh::SideSetEntry sideSetEntry = (*sidesets[0])[0];
      EXPECT_EQ(mesh.bulk.identifier(sideSetEntry.element), expectedId);
      EXPECT_EQ(sideSetEntry.side, expectedOrdinal);
    }
    else {
      ASSERT_EQ(sidesets.size(), 1u);
      EXPECT_EQ(sidesets[0]->size(), 0u);
    }
  }

  void compare_entity_rank_names(const MeshFromFile& meshA, const MeshFromFile& meshB) const
  {
    EXPECT_EQ(meshA.meta.entity_rank_names(), meshB.meta.entity_rank_names());
  }

  void verify_transient_field_names(const MeshFromFile& mesh, const std::string& fieldBaseName) const
  {
    const stk::mesh::FieldVector transientFields = stk::io::get_transient_fields(mesh.meta, stk::topology::NODE_RANK);
    verify_transient_field_name(transientFields[0], fieldBaseName+"_scalar");
    verify_transient_field_name(transientFields[1], fieldBaseName+"_vector");
  }

  void verify_transient_fields(MeshFromFile& mesh) const
  {
    stk::io::StkMeshIoBroker& broker = mesh.broker;

    const stk::mesh::FieldVector transientFields = stk::io::get_transient_fields(mesh.meta);
    const std::vector<double> timeSteps = broker.get_time_steps();

    for(int iStep=0; iStep<broker.get_num_time_steps(); iStep++)
    {
      double readTime = broker.read_defined_input_fields_at_step(iStep+1, nullptr);
      EXPECT_EQ(timeSteps[iStep], readTime);

      for(stk::mesh::FieldBase* field : transientFields)
        verify_transient_field_values(mesh.bulk, field, readTime);
    }
  }

private:
  void verify_global_variable_names(const MeshFromFile& mesh, const std::string& baseName) const
  {
    std::vector<std::string> goldNames = {baseName+"_double",
                                          baseName+"_int",
                                          baseName+"_real_vec"};
    std::vector<std::string> names;
    mesh.broker.get_global_variable_names(names);
    EXPECT_EQ(names, goldNames);
  }

  void verify_global_double(const MeshFromFile& mesh, const std::string& variable, double goldValue) const
  {
    double value;
    EXPECT_TRUE(mesh.broker.get_global(variable+"_double", value));
    EXPECT_NEAR(value, goldValue, m_epsilon);
  }

  void verify_global_int(const MeshFromFile& mesh, const std::string& variable, int goldValue) const
  {
    int value;
    EXPECT_TRUE(mesh.broker.get_global(variable+"_int", value));
    EXPECT_EQ(value, goldValue);
  }

  void verify_global_real_vec(const MeshFromFile& mesh, const std::string& variable, double goldValue) const
  {
    std::vector<double> values;
    EXPECT_TRUE(mesh.broker.get_global(variable+"_real_vec", values));

    EXPECT_EQ(values.size(), 3u);
    for (double value : values) {
      EXPECT_NEAR(value, goldValue, m_epsilon);
    }
  }

  void verify_transient_field_name(stk::mesh::FieldBase* field, const std::string& fieldName) const
  {
    EXPECT_EQ(0, strcasecmp(fieldName.c_str(), field->name().c_str()))
        << fieldName << "   " << field->name();
  }

  void verify_transient_field_values(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase* field, double timeStep) const
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

  MPI_Comm m_comm;
  const double m_epsilon;
};

class TransientFieldBalance : public stk::unit_test_util::MeshFixture
{
public:
  TransientFieldBalance()
    : fileBaseName("transient.e"),
      balanceRunner(get_comm()),
      writer(get_comm(), fileBaseName),
      verifier(get_comm()),
      m_initialMeshPtr(new MeshFromFile(get_comm())),
      m_balancedMeshPtr(new MeshFromFile(get_comm()))
  {
    balanceRunner.set_filename(fileBaseName);
    balanceRunner.set_output_dir(".");
    balanceRunner.set_app_type_defaults(stk::balance::SD_DEFAULTS);
    balanceRunner.set_decomp_method("rcb");
  }

  MeshFromFile& get_initial_mesh()
  {
    if (m_initialMeshPtr->is_empty()) read_initial_mesh();
    return *m_initialMeshPtr;
  }

  MeshFromFile& get_balanced_mesh()
  {
    if (m_balancedMeshPtr->is_empty()) read_balanced_mesh();
    return *m_balancedMeshPtr;
  }

  void cleanup_files()
  {
    unlink_serial_file(fileBaseName);
    unlink_parallel_file(fileBaseName);
  }

private:
  void read_initial_mesh()
  {
    m_initialMeshPtr->fill_from_serial(fileBaseName);
  }

  void read_balanced_mesh()
  {
    m_balancedMeshPtr->fill_from_parallel(fileBaseName);
  }

  void unlink_parallel_file(const std::string& baseName)
  {
    std::string parallelName = stk::io::construct_parallel_filename(baseName, get_parallel_size(), get_parallel_rank());
    unlink_serial_file(parallelName);
  }

  void unlink_serial_file(const std::string& fileName)
  {
    unlink(fileName.c_str());
  }

protected:
  const std::string fileBaseName;
  StkBalanceRunner balanceRunner;
  TransientWriter writer;
  const TransientVerifier verifier;

private:
  const std::unique_ptr<MeshFromFile> m_initialMeshPtr;
  const std::unique_ptr<MeshFromFile> m_balancedMeshPtr;
};

TEST_F(TransientFieldBalance, verifyStaticDataTransfer)
{
  if (get_parallel_size() == 2 ||
      get_parallel_size() == 4 ||
      get_parallel_size() == 8 ||
      get_parallel_size() == 16) {

    writer.write_static_mesh("1x4x4");

    MeshFromFile& initialMesh = get_initial_mesh();
    verifier.verify_time_steps(initialMesh, {});
    verifier.verify_num_transient_fields(initialMesh, 0u);

    balanceRunner.run_end_to_end();

    MeshFromFile& balancedMesh = get_balanced_mesh();
    verifier.verify_time_steps(balancedMesh, {});
    verifier.verify_num_transient_fields(balancedMesh, 0u);

    cleanup_files();
  }
}

TEST_F(TransientFieldBalance, verifyNumberOfSteps)
{
  if (get_parallel_size() != 2) return;

  const std::vector<double> timeSteps = {0.0, 1.0, 2.0, 3.0, 4.0};

  writer.set_time_steps(timeSteps);
  writer.write_transient_mesh("1x1x20");

  MeshFromFile& initialMesh = get_initial_mesh();
  verifier.verify_time_steps(initialMesh, timeSteps);
  verifier.verify_num_transient_fields(initialMesh, 2u);

  balanceRunner.run_end_to_end();

  MeshFromFile& balancedMesh = get_balanced_mesh();
  verifier.verify_time_steps(balancedMesh, timeSteps);

  cleanup_files();
}

TEST_F(TransientFieldBalance, verifyGlobalVariable)
{
  if (get_parallel_size() != 2) return;

  const std::vector<double> timeSteps = {0.0, 1.0, 2.0, 3.0, 4.0};
  const std::string globalVariableName = "test_time";

  writer.set_time_steps(timeSteps);
  writer.set_global_variable_name(globalVariableName);
  writer.write_transient_mesh("1x1x20");

  MeshFromFile& initialMesh = get_initial_mesh();
  verifier.verify_time_steps(initialMesh, timeSteps);
  verifier.verify_num_transient_fields(initialMesh, 2u);

  balanceRunner.run_end_to_end();

  MeshFromFile& balancedMesh = get_balanced_mesh();
  verifier.verify_global_variables_at_each_time_step(balancedMesh, globalVariableName, timeSteps);

  cleanup_files();
}

TEST_F(TransientFieldBalance, verifyTransientDataTransferOnFourProcessors)
{
  if (get_parallel_size() != 4) return;

  const std::vector<double> timeSteps = {0.0, 1.0, 2.0, 3.0, 4.0};
  const std::string fieldName = "myField";

  writer.set_time_steps(timeSteps);
  writer.set_field_name(fieldName);
  writer.write_transient_mesh("1x4x4");

  MeshFromFile& initialMesh = get_initial_mesh();
  verifier.verify_time_steps(initialMesh, timeSteps);
  verifier.verify_num_transient_fields(initialMesh, 2u);
  verifier.verify_transient_field_names(initialMesh, fieldName);

  balanceRunner.run_end_to_end();

  MeshFromFile& balancedMesh = get_balanced_mesh();
  verifier.verify_time_steps(balancedMesh, timeSteps);
  verifier.verify_num_transient_fields(balancedMesh, 2u);
  verifier.verify_transient_field_names(balancedMesh, fieldName);

  verifier.compare_entity_rank_names(initialMesh, balancedMesh);

  verifier.verify_transient_fields(balancedMesh);

  cleanup_files();
}

TEST_F(TransientFieldBalance, verifyTransientDataTransferWithSidesets)
{
  if (get_parallel_size() != 2) return;

  const stk::mesh::EntityId expectedId = 1;
  const stk::mesh::ConnectivityOrdinal expectedOrdinal = 5;

  writer.write_two_element_mesh_with_sideset(stk::unit_test_util::INCREASING);

  MeshFromFile& initialMesh = get_initial_mesh();
  verifier.verify_sideset_orientation(initialMesh, expectedId, expectedOrdinal);

  balanceRunner.run_end_to_end();

  MeshFromFile& balancedMesh = get_balanced_mesh();
  verifier.verify_sideset_orientation(balancedMesh, expectedId, expectedOrdinal);

  cleanup_files();
}

TEST_F(TransientFieldBalance, verifyTransientDataTransferWithSidesetsOnMovedElements)
{
  if (get_parallel_size() != 2) return;

  const stk::mesh::EntityId expectedId = 2;
  const stk::mesh::ConnectivityOrdinal expectedOrdinal = 5;

  writer.write_two_element_mesh_with_sideset(stk::unit_test_util::DECREASING);

  MeshFromFile& initialMesh = get_initial_mesh();
  verifier.verify_sideset_orientation(initialMesh, expectedId, expectedOrdinal);

  balanceRunner.set_initial_decomp_method("BLOCK");
  balanceRunner.run_end_to_end();

  MeshFromFile& balancedMesh = get_balanced_mesh();
  verifier.verify_sideset_orientation(balancedMesh, expectedId, expectedOrdinal);

  cleanup_files();
}

TEST(TestTransientFieldBalance, verifyThrowIfInputFileEqualsOutputFile)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    std::string serialMeshName = "sixteen_hex_transient.e";
    std::string parallelOutputMeshName = "sixteen_hex_transient.e";

    stk::unit_test_util::generated_mesh_to_file_in_serial("1x4x4", serialMeshName);

    stk::balance::BasicZoltan2Settings rcbOptions;
    EXPECT_THROW(stk::balance::run_stk_balance_with_settings(parallelOutputMeshName, serialMeshName, MPI_COMM_WORLD, rcbOptions), std::logic_error);
    unlink(serialMeshName.c_str());
  }
}

TEST(TestTransientFieldBalance, verifyThrowIfInputFileEqualsDotSlashOutputFile)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    std::string serialMeshName = "sixteen_hex_transient.e";
    std::string parallelOutputMeshName = "./sixteen_hex_transient.e";

    stk::unit_test_util::generated_mesh_to_file_in_serial("1x4x4", serialMeshName);

    stk::balance::BasicZoltan2Settings rcbOptions;
    EXPECT_THROW(stk::balance::run_stk_balance_with_settings(parallelOutputMeshName, serialMeshName, MPI_COMM_WORLD, rcbOptions), std::logic_error);
    unlink(serialMeshName.c_str());
  }
}
