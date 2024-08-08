#include "gtest/gtest.h"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_balance/setup/LifeCycle.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include "stk_unit_test_utils/GeneratedMeshToFile.hpp"
#include <vector>
#include <unistd.h>

class Args
{
public:
  Args(const std::vector<std::string> & arguments)
    : m_stringArgs(arguments),
      m_argc(m_stringArgs.size()),
      m_argv(arguments.empty() ? nullptr : new const char*[m_argc])
  {
    for (int i = 0; i < m_argc; ++i) {
      m_argv[i] = m_stringArgs[i].c_str();
    }
  }

  ~Args()
  {
    delete [] m_argv;
  }

  int argc() { return m_argc; }
  const char** argv() { return m_argv; }

private:
  const std::vector<std::string> m_stringArgs;
  int m_argc;
  const char** m_argv;
};

class TestLifeCycle : public stk::unit_test_util::MeshFixture
{
protected:
  TestLifeCycle()
    : m_execName("stk_balance"),
      m_inFile("mesh.g")
  {
  }

  void build_serial_mesh() {
    stk::unit_test_util::generated_mesh_to_file_in_serial("1x1x4", m_inFile);
    MPI_Barrier(get_comm());
  }

  void build_parallel_mesh() {
    stk::unit_test_util::GeneratedMeshToFile gMesh(get_comm(), stk::mesh::BulkData::AUTO_AURA);
    const bool useBigIds = false;
    gMesh.setup_mesh("1x1x4", m_inFile, useBigIds);
    gMesh.write_mesh();
  }

  void build_parallel_mesh_with_big_ids() {
    stk::unit_test_util::GeneratedMeshToFile gMesh(get_comm(), stk::mesh::BulkData::AUTO_AURA);
    const bool useBigIds = true;
    gMesh.setup_mesh("1x1x4", m_inFile, useBigIds);
    gMesh.write_mesh();
  }

  void build_serial_mesh_with_transient_field_data() {
    std::vector<double> transientTimeSteps = {0.0, 1.0, 2.0};
    std::string transientFieldName = "transient_field";
    std::string globalVariableName = "global_variable";
    stk::unit_test_util::generated_mesh_with_transient_data_to_file_in_serial("1x1x4",
                                                                                             m_inFile,
                                                                                             transientFieldName,
                                                                                             stk::topology::NODE_RANK,
                                                                                             globalVariableName,
                                                                                             transientTimeSteps,
                                                                                             stk::unit_test_util::IdAndTimeFieldValueSetter());
    MPI_Barrier(get_comm());
  }

  void build_parallel_mesh_with_transient_field_data() {
    std::vector<double> transientTimeSteps = {0.0, 1.0, 2.0};
    std::string transientFieldName = "transient_field";
    std::string globalVariableName = "global_variable";
    stk::unit_test_util::GeneratedMeshToFileWithTransientFields gMesh(get_comm(),
                                                                                     stk::mesh::BulkData::AUTO_AURA,
                                                                                     transientFieldName,
                                                                                     stk::topology::NODE_RANK);
    const bool useBigIds = false;
    gMesh.setup_mesh("1x1x4", m_inFile, useBigIds);
    gMesh.write_mesh_with_field(transientTimeSteps,
                                stk::unit_test_util::IdAndTimeFieldValueSetter(),
                                globalVariableName);
  }

  std::vector<const char*> build_command_line(const std::vector<std::string>& options) const {
    std::vector<const char*> args = {m_execName.c_str()};

    for (const std::string & option : options) {
      args.push_back(option.c_str());
    }

    return args;
  }

  void clean_up_temporary_files(const stk::balance::LifeCycle & lifeCycle) {
    MPI_Barrier(get_comm());

    if (get_parallel_rank() == 0) {
      const stk::balance::BalanceSettings & balanceSettings = lifeCycle.get_balance_settings();
      unlink(balanceSettings.get_log_filename().c_str());

      for (unsigned i = 0; i < balanceSettings.get_num_input_processors(); ++i) {
        unlink(stk::io::construct_filename_for_serial_or_parallel(balanceSettings.get_input_filename(),
                                                                  balanceSettings.get_num_input_processors(), i).c_str());
      }

      if (lifeCycle.exit_code() == stk::balance::LifeCycleStatus::SUCCESS) {
        for (unsigned i = 0; i < balanceSettings.get_num_output_processors(); ++i) {
          unlink(stk::io::construct_filename_for_serial_or_parallel(balanceSettings.get_output_filename(),
                                                                    balanceSettings.get_num_output_processors(), i).c_str());
        }
      }
    }
  }

  std::string m_execName;
  std::string m_inFile;

private:
};

TEST_F(TestLifeCycle, balance)
{
  build_serial_mesh();
  Args args({m_execName, m_inFile});
  stk::balance::LifeCycle lifeCycle(get_comm(), args.argc(), args.argv());
  lifeCycle.run();

  EXPECT_EQ(lifeCycle.exit_code(), stk::balance::LifeCycleStatus::SUCCESS);

  clean_up_temporary_files(lifeCycle);
}

TEST_F(TestLifeCycle, balanceWithTransientFieldData)
{
  build_serial_mesh_with_transient_field_data();
  Args args({m_execName, m_inFile});
  stk::balance::LifeCycle lifeCycle(get_comm(), args.argc(), args.argv());
  lifeCycle.run();

  EXPECT_EQ(lifeCycle.exit_code(), stk::balance::LifeCycleStatus::SUCCESS);

  clean_up_temporary_files(lifeCycle);
}

TEST_F(TestLifeCycle, rebalance)
{
  build_parallel_mesh();
  Args args({m_execName, m_inFile, "--rebalance-to=4"});
  stk::balance::LifeCycle lifeCycle(get_comm(), args.argc(), args.argv());
  lifeCycle.run();

  EXPECT_EQ(lifeCycle.exit_code(), stk::balance::LifeCycleStatus::SUCCESS);

  clean_up_temporary_files(lifeCycle);
}

TEST_F(TestLifeCycle, rebalanceWithBigIds)
{
  build_parallel_mesh_with_big_ids();
  Args args({m_execName, m_inFile, "--rebalance-to=4"});
  stk::balance::LifeCycle lifeCycle(get_comm(), args.argc(), args.argv());
  lifeCycle.run();

  if (get_parallel_size() == 4) {
    EXPECT_EQ(lifeCycle.exit_code(), stk::balance::LifeCycleStatus::REBALANCE_CORRUPTION_ERROR);
  }
  else {
    EXPECT_EQ(lifeCycle.exit_code(), stk::balance::LifeCycleStatus::SUCCESS);
  }

  clean_up_temporary_files(lifeCycle);
}

TEST_F(TestLifeCycle, rebalanceWithTransientFieldData)
{
  build_parallel_mesh_with_transient_field_data();
  Args args({m_execName, m_inFile, "--rebalance-to=4"});
  stk::balance::LifeCycle lifeCycle(get_comm(), args.argc(), args.argv());
  lifeCycle.run();

  if (get_parallel_size() == 4) {
    EXPECT_EQ(lifeCycle.exit_code(), stk::balance::LifeCycleStatus::REBALANCE_CORRUPTION_ERROR);
  }
  else {
    EXPECT_EQ(lifeCycle.exit_code(), stk::balance::LifeCycleStatus::SUCCESS);
  }

  clean_up_temporary_files(lifeCycle);
}



