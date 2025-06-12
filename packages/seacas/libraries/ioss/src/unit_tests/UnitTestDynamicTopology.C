// Copyright(C) 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <string>
#include <vector>

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include <algorithm>
#include <cctype> // std::tolower
#include <fstream>
#include <functional>
#include <iostream>
#include <random>

#include <unistd.h> // for unlink

#include "Ionit_Initializer.h"
#include "Ioss_ChangeSet.h"
#include "Ioss_ChangeSetFactory.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include "Ioss_DynamicTopology.h"
#include "Ioss_DynamicTopologyBroker.h"
#include "Ioss_DynamicTopologyFileControl.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_Field.h" // for Field, etc
#include "Ioss_FileInfo.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_Utils.h"

#include "exodus/Ioex_DatabaseIO.h"

namespace {
  std::string get_many_block_mesh_desc(unsigned numBlocks)
  {
    std::ostringstream    oss;
    std::vector<unsigned> elementIds(numBlocks);
    std::iota(elementIds.begin(), elementIds.end(), 1);

    unsigned proc = 0;
    for (unsigned i = 0; i < numBlocks; ++i) {
      unsigned elemId      = elementIds[i];
      unsigned firstNodeId = i * 4 + 1;
      oss << proc << "," << elemId << ",HEX_8,";
      for (unsigned node = firstNodeId; node < firstNodeId + 8; ++node) {
        oss << node << ",";
      }
      unsigned blockId = i + 1;
      oss << "block_" << blockId;

      if (i < numBlocks - 1) {
        oss << "\n";
      }

      proc++;
    }

    oss << "|coordinates:";

    std::vector<double> planeCoords = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};

    for (double coord : planeCoords) {
      oss << coord << ",";
    }

    for (unsigned i = 1; i <= numBlocks; ++i) {
      for (unsigned point = 0; point < 4; ++point) {
        planeCoords[3 * point + 2] += 1;
      }

      for (double coord : planeCoords) {
        oss << coord << ",";
      }
    }

    return oss.str();
  }

  void define_model(const Ioss::Region &i_region, Ioss::Region &o_region)
  {
    Ioss::DatabaseIO *o_database = o_region.get_database();

    o_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

    Ioss::NodeBlock *i_nb        = i_region.get_node_blocks()[0];
    int64_t          spatial_dim = 3;
    int64_t          num_nodes   = i_nb->entity_count();
    Ioss::NodeBlock *o_nb = new Ioss::NodeBlock(o_database, "nodeblock_1", num_nodes, spatial_dim);
    o_region.add(o_nb);

    for (Ioss::ElementBlock *i_eb : i_region.get_element_blocks()) {
      Ioss::ElementBlock *o_eb = new Ioss::ElementBlock(
          o_database, i_eb->name(), i_eb->topology()->name(), i_eb->entity_count());
      o_eb->property_add(i_eb->get_property("id"));
      o_region.add(o_eb);
    }

    o_region.end_mode(Ioss::STATE_DEFINE_MODEL);
  }

  void write_model(const Ioss::Region &i_region, Ioss::Region &o_region)
  {
    Ioss::NodeBlock *i_nb = i_region.get_node_blocks()[0];
    Ioss::NodeBlock *o_nb = o_region.get_node_blocks()[0];

    o_region.begin_mode(Ioss::STATE_MODEL);
    std::vector<double> coordinates;
    std::vector<int>    node_ids;
    i_nb->get_field_data("ids", node_ids);
    i_nb->get_field_data("mesh_model_coordinates", coordinates);

    o_nb->put_field_data("ids", node_ids);
    o_nb->put_field_data("mesh_model_coordinates", coordinates);

    for (Ioss::ElementBlock *i_eb : i_region.get_element_blocks()) {
      Ioss::ElementBlock *o_eb = o_region.get_element_block(i_eb->name());
      std::vector<int>    elem_ids;
      std::vector<int>    connectivity;

      i_eb->get_field_data("ids", elem_ids);
      i_eb->get_field_data("connectivity", connectivity);

      o_eb->put_field_data("ids", elem_ids);
      o_eb->put_field_data("connectivity", connectivity);
    }

    o_region.end_mode(Ioss::STATE_MODEL);
  }

  void define_transient(Ioss::Region &o_region, const std::string &elemFieldName)
  {
    o_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

    for (Ioss::ElementBlock *o_eb : o_region.get_element_blocks()) {
      size_t      num_elem = o_eb->entity_count();
      std::string storage  = "scalar";

      Ioss::Field field(elemFieldName, Ioss::Field::REAL, storage, 1, Ioss::Field::Field::TRANSIENT,
                        num_elem);
      o_eb->field_add(field);
    }
    o_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  }

  int write_transient(Ioss::Region &o_region, const std::string &elemFieldName, const double time)
  {
    o_region.begin_mode(Ioss::STATE_TRANSIENT);
    int step = o_region.add_state(time);
    o_region.begin_state(step);

    for (Ioss::ElementBlock *o_eb : o_region.get_element_blocks()) {
      size_t num_elem = o_eb->entity_count();

      std::vector<double> field_data(num_elem);
      std::vector<int>    elem_ids;

      o_eb->get_field_data("ids", elem_ids);
      for (size_t i = 0; i < elem_ids.size(); i++) {
        field_data[i] = (double)elem_ids[i] + 100 * time;
      }

      o_eb->put_field_data(elemFieldName, field_data);
    }

    o_region.end_state(step);
    o_region.end_mode(Ioss::STATE_TRANSIENT);

    return step;
  }

  class Observer : public Ioss::DynamicTopologyObserver
  {
  public:
    Observer(Ioss::Region &inputRegion_, const std::string &elemFieldName_,
             const Ioss::FileControlOption fileControlOption_)
        : Ioss::DynamicTopologyObserver(nullptr), inputRegion(inputRegion_),
          elemFieldName(elemFieldName_), fileControlOption(fileControlOption_)
    {
    }

    virtual ~Observer() {}

    void define_model() override { ::define_model(inputRegion, *get_region()); }

    void write_model() override { ::write_model(inputRegion, *get_region()); }

    void define_transient() override { ::define_transient(*get_region(), elemFieldName); }

    Ioss::FileControlOption get_control_option() const { return fileControlOption; }

  private:
    Observer();

    Ioss::Region           &inputRegion;
    const std::string       elemFieldName;
    Ioss::FileControlOption fileControlOption;
  };

  struct OutputParams
  {
    OutputParams(const std::string &outFile_) : outFile(outFile_) {}

    OutputParams(const std::string &outFile_, const std::string &elemFieldName_)
        : outFile(outFile_), elemFieldName(elemFieldName_)
    {
    }

    void set_data(const std::vector<double> &output_times_, const std::vector<bool> &output_steps_,
                  const std::vector<bool> &modification_steps_)
    {
      ASSERT_EQ(output_times_.size(), output_steps_.size());
      ASSERT_EQ(output_times_.size(), modification_steps_.size());

      size_t numSteps = output_times_.size();
      for (size_t i = 1; i < numSteps; i++) {
        // Monotone increasing
        ASSERT_TRUE(output_times_[i] > output_times_[i - 1]);
      }

      output_times       = output_times_;
      output_steps       = output_steps_;
      modification_steps = modification_steps_;
    }

    void set_data(const std::vector<bool> &output_steps_,
                  const std::vector<bool> &modification_steps_)
    {
      ASSERT_EQ(output_steps_.size(), modification_steps_.size());

      size_t numSteps = output_steps_.size();
      for (size_t i = 0; i < numSteps; i++) {
        output_times.push_back((double)i);
      }
      output_steps       = output_steps_;
      modification_steps = modification_steps_;
    }

    struct OutputParams &add(const double time, const bool do_output, const bool do_modification)
    {
      size_t numSteps = output_times.size();
      if (numSteps > 0) {
        // Monotone increasing
        EXPECT_TRUE(time > output_times[numSteps - 1]);
      }

      output_times.push_back(time);
      output_steps.push_back(do_output);
      modification_steps.push_back(do_modification);

      return *this;
    }

    void clear()
    {
      output_times.clear();
      output_steps.clear();
      modification_steps.clear();
    }

    void set_cyclic_count(unsigned cyclicCount_) { cyclicCount = cyclicCount_; }

    unsigned            cyclicCount{0};
    std::string         outFile{"file.g"};
    std::string         elemFieldName{"elem_field"};
    std::vector<double> output_times;
    std::vector<bool>   output_steps;
    std::vector<bool>   modification_steps;
  };

  void do_output(Ioss::Region &o_region, const OutputParams &params, size_t step, double &minTime,
                 int &maxStep, bool &doneOutputAfterModification)
  {
    if (params.output_steps[step]) {
      if (!doneOutputAfterModification) {
        minTime = params.output_times[step];
      }

      write_transient(o_region, params.elemFieldName, params.output_times[step]);

      auto min_result = o_region.get_min_time();
      EXPECT_EQ(1, min_result.first);
      EXPECT_NEAR(minTime, min_result.second, 1.0e-6);

      auto max_result = o_region.get_max_time();
      EXPECT_EQ(maxStep, max_result.first);
      EXPECT_NEAR(params.output_times[step], max_result.second, 1.0e-6);

      maxStep++;
      doneOutputAfterModification = true;
    }
  }

  void run_topology_change(const Ioss::Region &i_region, Ioss::Region &o_region,
                           const OutputParams &params)
  {
    auto observer = o_region.get_mesh_modification_observer();

    define_model(i_region, o_region);
    write_model(i_region, o_region);

    define_transient(o_region, params.elemFieldName);

    auto numSteps = params.output_steps.size();

    int maxStep = 1;

    double minTime = numSteps > 0 ? params.output_times[0] : 0.0;

    bool doneOutputAfterModification = true;

    for (size_t i = 0; i < numSteps; i++) {
      if (params.modification_steps[i]) {
        observer->set_topology_modification(Ioss::TOPOLOGY_UNKNOWN);
        maxStep                     = 1;
        doneOutputAfterModification = false;
      }

      do_output(o_region, params, i, minTime, maxStep, doneOutputAfterModification);
    }
  }

  void cleanup_linear_multi_files(const std::string &outFile, int numOutputs = 1)
  {
    Ioss::ParallelUtils util(Ioss::ParallelUtils::comm_world());

    for (int i = 1; i <= numOutputs; i++) {
      std::string baseFile =
          Ioss::DynamicTopologyFileControl::get_linear_database_filename(outFile, i);
      std::string parallelFile =
          Ioss::Utils::decode_filename(baseFile, util.parallel_rank(), util.parallel_size());
      unlink(parallelFile.c_str());
    }
  }

  std::shared_ptr<Ioss::Region> construct_region(const OutputParams    &params,
                                                 Ioss::PropertyManager &propertyManager,
                                                 Ioss::DatabaseUsage    db_usage,
                                                 const std::string     &name)
  {
    std::string outFile = params.outFile;
    if (params.cyclicCount > 0) {
      outFile = Ioss::DynamicTopologyFileControl::get_cyclic_database_filename(
          params.outFile, params.cyclicCount, 0);
    }
    propertyManager.add(Ioss::Property("base_filename", params.outFile));
    Ioss::DatabaseIO *database = Ioss::IOFactory::create(
        "exodus", outFile, db_usage, Ioss::ParallelUtils::comm_world(), propertyManager);

    auto region = std::make_shared<Ioss::Region>(database, name);
    region->set_file_cyclic_count(params.cyclicCount);
    EXPECT_TRUE(database != nullptr);
    EXPECT_TRUE(database->ok(true));

    return region;
  }

  void run_multi_file_topology_change(const OutputParams &params)
  {
    Ioss::Init::Initializer io;
    Ioss::ParallelUtils     util(Ioss::ParallelUtils::comm_world());

    int numBlocks = util.parallel_size();

    std::string meshDesc = get_many_block_mesh_desc(numBlocks);

    Ioss::PropertyManager propertyManager;

    Ioss::DatabaseIO *i_database = Ioss::IOFactory::create(
        "textmesh", meshDesc, Ioss::READ_MODEL, Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region i_region(i_database, "input_model");
    EXPECT_TRUE(i_database != nullptr);
    EXPECT_TRUE(i_database->ok(true));

    auto o_region = construct_region(params, propertyManager, Ioss::WRITE_RESULTS, "output_model");

    auto fileControlOption = Ioss::FileControlOption::CONTROL_AUTO_MULTI_FILE;
    auto observer = std::make_shared<Observer>(i_region, params.elemFieldName, fileControlOption);
    o_region->register_mesh_modification_observer(observer);

    run_topology_change(i_region, *o_region, params);
  }

  TEST(TestDynamicWrite, multi_file_simple_topology_modification)
  {
    std::string outFile("multiFileManyBlocks.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    std::vector<bool> output_steps{true, true, true, true, true, true};
    std::vector<bool> modification_steps{false, true, false, true, true, false};

    params.set_data(output_steps, modification_steps);

    cleanup_linear_multi_files(outFile, params.output_times.size());
    run_multi_file_topology_change(params);
    cleanup_linear_multi_files(outFile, params.output_times.size());
  }

  void cleanup_cyclic_multi_files(const std::string &outFile, unsigned cyclicCount = 3)
  {
    Ioss::ParallelUtils util(Ioss::ParallelUtils::comm_world());

    for (unsigned i = 1; i <= cyclicCount; i++) {
      std::string baseFile =
          Ioss::DynamicTopologyFileControl::get_cyclic_database_filename(outFile, cyclicCount, i);
      std::string parallelFile =
          Ioss::Utils::decode_filename(baseFile, util.parallel_rank(), util.parallel_size());
      unlink(parallelFile.c_str());
    }
  }

  TEST(TestDynamicWrite, multi_file_cyclic_topology_modification)
  {
    std::string outFile("cyclicMultiFileManyBlocks.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    std::vector<double> output_times{0.0, 0.5, 1.5, 1.75, 2.0, 3.0};
    std::vector<bool>   output_steps{true, true, true, true, true, true};
    std::vector<bool>   modification_steps{false, true, false, true, true, false};

    params.set_data(output_times, output_steps, modification_steps);
    params.set_cyclic_count(3);

    cleanup_cyclic_multi_files(outFile, params.cyclicCount);
    run_multi_file_topology_change(params);
    cleanup_cyclic_multi_files(outFile, params.cyclicCount);
  }

  void fill_internal_file_change_set_gold_names(const int                 numChangeSets,
                                                std::vector<std::string> &gold_names,
                                                std::vector<std::string> &gold_full_names)
  {
    gold_names.clear();
    gold_full_names.clear();

    for (int i = 1; i <= numChangeSets; i++) {
      std::string setName = Ioss::DynamicTopologyFileControl::get_internal_file_change_set_name(i);

      gold_names.push_back(setName);
      gold_full_names.push_back("/" + setName);
    }
  }

  void test_internal_file_change_set_names(Ioss::DatabaseIO *database)
  {
    Ioss::NameList names      = database->internal_change_set_describe(false);
    Ioss::NameList full_names = database->internal_change_set_describe(true);

    std::vector<std::string> gold_names;
    std::vector<std::string> gold_full_names;

    fill_internal_file_change_set_gold_names(database->num_internal_change_set(), gold_names,
                                             gold_full_names);

    EXPECT_EQ(gold_names, names);
    EXPECT_EQ(gold_full_names, full_names);
  }

  void cleanup_single_file(const std::string &outFile)
  {
    Ioss::ParallelUtils util(Ioss::ParallelUtils::comm_world());

    std::string file =
        Ioss::Utils::decode_filename(outFile, util.parallel_rank(), util.parallel_size());
    unlink(file.c_str());
  }

  void run_single_file_simple_topology_change(const OutputParams &params)
  {
    Ioss::Init::Initializer io;
    Ioss::ParallelUtils     util(Ioss::ParallelUtils::comm_world());

    int numBlocks = util.parallel_size();

    std::string meshDesc = get_many_block_mesh_desc(numBlocks);

    Ioss::PropertyManager propertyManager;

    Ioss::DatabaseIO *i_database = Ioss::IOFactory::create(
        "textmesh", meshDesc, Ioss::READ_MODEL, Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region i_region(i_database, "input_model");
    EXPECT_TRUE(i_database != nullptr);
    EXPECT_TRUE(i_database->ok(true));

    propertyManager.add(Ioss::Property("ENABLE_FILE_GROUPS", 1));
    propertyManager.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND_GROUP));
    Ioss::DatabaseIO *o_database =
        Ioss::IOFactory::create("exodus", params.outFile, Ioss::WRITE_RESULTS,
                                Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region o_region(o_database, "output_model");
    EXPECT_TRUE(o_database != nullptr);
    EXPECT_TRUE(o_database->ok(true));
    EXPECT_TRUE(o_database->supports_internal_change_set());

    auto fileControlOption = Ioss::FileControlOption::CONTROL_AUTO_GROUP_FILE;
    auto observer = std::make_shared<Observer>(i_region, params.elemFieldName, fileControlOption);
    o_region.register_mesh_modification_observer(observer);

    run_topology_change(i_region, o_region, params);
    test_internal_file_change_set_names(o_database);
  }

  TEST(TestDynamicWrite, single_file_simple_topology_modification)
  {
    std::string outFile("singleFileManyBlocks.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    params.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, true, false)
        .add(3.0, true, true)
        .add(4.0, true, true)
        .add(5.0, true, false);

    cleanup_single_file(outFile);
    run_single_file_simple_topology_change(params);
    cleanup_single_file(outFile);
  }

  TEST(TestDynamicWrite, single_file_groups_not_enabled)
  {
    Ioss::Init::Initializer io;
    Ioss::ParallelUtils     util(Ioss::ParallelUtils::comm_world());

    int numBlocks = util.parallel_size();
    if (numBlocks > 1)
      GTEST_SKIP();

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1"
                           "|coordinates:0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1";

    Ioss::PropertyManager propertyManager;

    Ioss::DatabaseIO *i_database = Ioss::IOFactory::create(
        "textmesh", meshDesc, Ioss::READ_MODEL, Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region i_region(i_database, "input_model");
    EXPECT_TRUE(i_database != nullptr);
    EXPECT_TRUE(i_database->ok(true));

    std::string outFile("singleFileGroupsNotEnabled.g");
    std::string elemFieldName = "elem_field";
    cleanup_single_file(outFile);

    // Need the line below to allow this to pass
    // propertyManager.add(Ioss::Property("ENABLE_FILE_GROUPS", 1));
    propertyManager.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND_GROUP));
    Ioss::DatabaseIO *o_database = Ioss::IOFactory::create(
        "exodus", outFile, Ioss::WRITE_RESULTS, Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region o_region(o_database, "output_model");
    EXPECT_TRUE(o_database != nullptr);
    EXPECT_TRUE(o_database->ok(true));

    auto fileControlOption = Ioss::FileControlOption::CONTROL_AUTO_GROUP_FILE;
    auto observer          = std::make_shared<Observer>(i_region, elemFieldName, fileControlOption);
    EXPECT_THROW(o_region.register_mesh_modification_observer(observer), std::runtime_error);
    cleanup_single_file(outFile);
  }

  TEST(TestDynamicWrite, create_subgroup_with_file_reopen)
  {
    std::string outFile("subgroupManyBlocks.g");
    std::string elemFieldName = "elem_field";

    Ioss::Init::Initializer io;
    Ioss::ParallelUtils     util(Ioss::ParallelUtils::comm_world());

    std::string file1 =
        Ioss::Utils::decode_filename(outFile, util.parallel_rank(), util.parallel_size());
    unlink(file1.c_str());

    int numBlocks = util.parallel_size();
    if (numBlocks > 1)
      GTEST_SKIP();

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1"
                           "|coordinates:0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1";

    Ioss::PropertyManager propertyManager;

    Ioss::DatabaseIO *i_database = Ioss::IOFactory::create(
        "textmesh", meshDesc, Ioss::READ_MODEL, Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region i_region(i_database, "input_model");
    EXPECT_TRUE(i_database != nullptr);
    EXPECT_TRUE(i_database->ok(true));

    {
      propertyManager.add(Ioss::Property("ENABLE_FILE_GROUPS", 1));
      Ioss::DatabaseIO *o_database =
          Ioss::IOFactory::create("exodus", outFile, Ioss::WRITE_RESULTS,
                                  Ioss::ParallelUtils::comm_world(), propertyManager);

      Ioex::BaseDatabaseIO *ex_database = dynamic_cast<Ioex::BaseDatabaseIO *>(o_database);
      EXPECT_TRUE(nullptr != ex_database);

      Ioss::Region o_region(o_database, "output_model");
      EXPECT_TRUE(o_database != nullptr);
      EXPECT_TRUE(o_database->ok(true));
      EXPECT_TRUE(o_database->supports_internal_change_set());

      ex_database->create_subgroup("GROUP_1");
    }

    {
      propertyManager.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND_GROUP));
      Ioss::DatabaseIO *o_database =
          Ioss::IOFactory::create("exodus", outFile, Ioss::WRITE_RESULTS,
                                  Ioss::ParallelUtils::comm_world(), propertyManager);

      Ioex::BaseDatabaseIO *ex_database = dynamic_cast<Ioex::BaseDatabaseIO *>(o_database);
      EXPECT_TRUE(nullptr != ex_database);

      Ioss::Region o_region(o_database, "output_model");
      EXPECT_TRUE(o_database != nullptr);
      EXPECT_TRUE(o_database->ok(true));

      EXPECT_TRUE(ex_database->supports_group());

      // Group pointer is automatically at first child
      ex_database->create_subgroup("GROUP_2");

      Ioss::NameList names      = ex_database->groups_describe(false);
      Ioss::NameList full_names = ex_database->groups_describe(true);

      std::vector<std::string> gold_names{"/", "GROUP_1", "GROUP_2"};
      std::vector<std::string> gold_full_names{"/", "/GROUP_1", "/GROUP_1/GROUP_2"};

      EXPECT_EQ(gold_names, names);
      EXPECT_EQ(gold_full_names, full_names);
    }

    unlink(file1.c_str());
  }

  TEST(TestDynamicWrite, create_subgroup_with_file_persistence_and_child_group)
  {
    std::string outFile("subgroupManyBlocks.g");
    std::string elemFieldName = "elem_field";

    Ioss::Init::Initializer io;
    Ioss::ParallelUtils     util(Ioss::ParallelUtils::comm_world());

    std::string file1 =
        Ioss::Utils::decode_filename(outFile, util.parallel_rank(), util.parallel_size());
    unlink(file1.c_str());

    int numBlocks = util.parallel_size();
    if (numBlocks > 1)
      GTEST_SKIP();

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1"
                           "|coordinates:0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1";

    Ioss::PropertyManager propertyManager;

    Ioss::DatabaseIO *i_database = Ioss::IOFactory::create(
        "textmesh", meshDesc, Ioss::READ_MODEL, Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region i_region(i_database, "input_model");
    EXPECT_TRUE(i_database != nullptr);
    EXPECT_TRUE(i_database->ok(true));

    {
      propertyManager.add(Ioss::Property("ENABLE_FILE_GROUPS", 1));
      propertyManager.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND_GROUP));
      Ioss::DatabaseIO *o_database =
          Ioss::IOFactory::create("exodus", outFile, Ioss::WRITE_RESULTS,
                                  Ioss::ParallelUtils::comm_world(), propertyManager);

      Ioex::BaseDatabaseIO *ex_database = dynamic_cast<Ioex::BaseDatabaseIO *>(o_database);
      EXPECT_TRUE(nullptr != ex_database);

      Ioss::Region o_region(o_database, "output_model");
      EXPECT_TRUE(o_database != nullptr);
      EXPECT_TRUE(o_database->ok(true));

      EXPECT_TRUE(ex_database->supports_group());

      ex_database->create_subgroup("GROUP_1");

      // Group pointer is at "GROUP_1" ... "GROUP_2" is a child
      ex_database->create_subgroup("GROUP_2");

      Ioss::NameList names      = ex_database->groups_describe(false);
      Ioss::NameList full_names = ex_database->groups_describe(true);

      std::vector<std::string> gold_names{"/", "GROUP_1", "GROUP_2"};
      std::vector<std::string> gold_full_names{"/", "/GROUP_1", "/GROUP_1/GROUP_2"};

      EXPECT_EQ(gold_names, names);
      EXPECT_EQ(gold_full_names, full_names);
    }

    unlink(file1.c_str());
  }

  TEST(TestDynamicWrite, create_subgroup_with_file_persistence_and_no_child_group)
  {
    std::string outFile("subgroupManyBlocks.g");
    std::string elemFieldName = "elem_field";

    Ioss::Init::Initializer io;
    Ioss::ParallelUtils     util(Ioss::ParallelUtils::comm_world());

    std::string file1 =
        Ioss::Utils::decode_filename(outFile, util.parallel_rank(), util.parallel_size());
    unlink(file1.c_str());

    int numBlocks = util.parallel_size();
    if (numBlocks > 1)
      GTEST_SKIP();

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1"
                           "|coordinates:0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1";

    Ioss::PropertyManager propertyManager;

    Ioss::DatabaseIO *i_database = Ioss::IOFactory::create(
        "textmesh", meshDesc, Ioss::READ_MODEL, Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region i_region(i_database, "input_model");
    EXPECT_TRUE(i_database != nullptr);
    EXPECT_TRUE(i_database->ok(true));

    {
      propertyManager.add(Ioss::Property("ENABLE_FILE_GROUPS", 1));
      propertyManager.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND_GROUP));
      Ioss::DatabaseIO *o_database =
          Ioss::IOFactory::create("exodus", outFile, Ioss::WRITE_RESULTS,
                                  Ioss::ParallelUtils::comm_world(), propertyManager);

      Ioex::BaseDatabaseIO *ex_database = dynamic_cast<Ioex::BaseDatabaseIO *>(o_database);
      EXPECT_TRUE(nullptr != ex_database);

      Ioss::Region o_region(o_database, "output_model");
      EXPECT_TRUE(o_database != nullptr);
      EXPECT_TRUE(o_database->ok(true));

      EXPECT_TRUE(ex_database->supports_group());

      ex_database->create_subgroup("GROUP_1");

      // Group pointer is reset to root group
      EXPECT_TRUE(ex_database->open_root_group());
      ex_database->create_subgroup("GROUP_2");

      Ioss::NameList names      = ex_database->groups_describe(false);
      Ioss::NameList full_names = ex_database->groups_describe(true);

      std::vector<std::string> gold_names{"/", "GROUP_1", "GROUP_2"};
      std::vector<std::string> gold_full_names{"/", "/GROUP_1", "/GROUP_2"};

      EXPECT_EQ(gold_names, names);
      EXPECT_EQ(gold_full_names, full_names);
    }

    unlink(file1.c_str());
  }

  void run_topology_change_with_multiple_output(const Ioss::Region &i_region,
                                                Ioss::Region &o_region1, Ioss::Region &o_region2,
                                                const OutputParams &params1,
                                                const OutputParams &params2)
  {
    ASSERT_EQ(params1.modification_steps, params2.modification_steps);

    auto observer1 = o_region1.get_mesh_modification_observer();
    auto observer2 = o_region2.get_mesh_modification_observer();

    define_model(i_region, o_region1);
    write_model(i_region, o_region1);
    define_transient(o_region1, params1.elemFieldName);

    define_model(i_region, o_region2);
    write_model(i_region, o_region2);
    define_transient(o_region2, params2.elemFieldName);

    auto numSteps = params1.output_steps.size();

    int maxStep1 = 1;
    int maxStep2 = 1;

    double minTime1 = numSteps > 0 ? params1.output_times[0] : 0.0;
    double minTime2 = numSteps > 0 ? params2.output_times[0] : 0.0;

    bool doneOutputAfterModification1 = true;
    bool doneOutputAfterModification2 = true;

    for (size_t i = 0; i < numSteps; i++) {
      if (params1.modification_steps[i]) {
        EXPECT_TRUE(params2.modification_steps[i]);

        observer1->set_topology_modification(Ioss::TOPOLOGY_UNKNOWN);
        maxStep1 = 1;
        maxStep2 = 1;

        EXPECT_EQ(Ioss::TOPOLOGY_UNKNOWN, observer1->get_topology_modification());
        EXPECT_EQ(Ioss::TOPOLOGY_UNKNOWN, observer2->get_topology_modification());

        doneOutputAfterModification1 = false;
        doneOutputAfterModification2 = false;
      }

      do_output(o_region1, params1, i, minTime1, maxStep1, doneOutputAfterModification1);
      do_output(o_region2, params2, i, minTime2, maxStep2, doneOutputAfterModification2);
    }
  }

  void run_single_file_simple_topology_change_with_multiple_output(const std::string  &model,
                                                                   const OutputParams &params1,
                                                                   const OutputParams &params2)
  {
    Ioss::Init::Initializer io;
    Ioss::ParallelUtils     util(Ioss::ParallelUtils::comm_world());

    auto broker = Ioss::DynamicTopologyBroker::broker();
    broker->register_model(model);

    int numBlocks = util.parallel_size();

    std::string meshDesc = get_many_block_mesh_desc(numBlocks);

    Ioss::PropertyManager propertyManager;

    Ioss::DatabaseIO *i_database = Ioss::IOFactory::create(
        "textmesh", meshDesc, Ioss::READ_MODEL, Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region i_region(i_database, "input_model");
    EXPECT_TRUE(i_database != nullptr);
    EXPECT_TRUE(i_database->ok(true));

    propertyManager.add(Ioss::Property("ENABLE_FILE_GROUPS", 1));
    propertyManager.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND_GROUP));

    Ioss::DatabaseIO *o_database1 =
        Ioss::IOFactory::create("exodus", params1.outFile, Ioss::WRITE_RESULTS,
                                Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region o_region1(o_database1, "region1");
    EXPECT_TRUE(o_database1 != nullptr);
    EXPECT_TRUE(o_database1->ok(true));
    EXPECT_TRUE(o_database1->supports_internal_change_set());

    auto fileControlOption = Ioss::FileControlOption::CONTROL_AUTO_GROUP_FILE;
    auto observer1 = std::make_shared<Observer>(i_region, params1.elemFieldName, fileControlOption);
    broker->register_observer(model, observer1, o_region1);

    Ioss::DatabaseIO *o_database2 =
        Ioss::IOFactory::create("exodus", params2.outFile, Ioss::WRITE_RESULTS,
                                Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region o_region2(o_database2, "region2");
    EXPECT_TRUE(o_database2 != nullptr);
    EXPECT_TRUE(o_database2->ok(true));
    EXPECT_TRUE(o_database2->supports_internal_change_set());

    auto observer2 = std::make_shared<Observer>(i_region, params2.elemFieldName, fileControlOption);
    broker->register_observer(model, observer2, o_region2);

    run_topology_change_with_multiple_output(i_region, o_region1, o_region2, params1, params2);

    test_internal_file_change_set_names(o_database1);
    test_internal_file_change_set_names(o_database2);
  }

  TEST(TestDynamicWrite, single_file_simple_topology_modification_with_multiple_output)
  {
    std::string outFile1("singleFileManyBlocks1.g");
    std::string outFile2("singleFileManyBlocks2.g");
    std::string elemFieldName = "elem_field";
    std::string model         = "multiple-output";

    OutputParams params1(outFile1, elemFieldName);

    params1.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, false, false)
        .add(3.0, true, true)
        .add(4.0, true, false)
        .add(5.0, true, true);

    OutputParams params2(outFile2, elemFieldName);

    params2.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, true, false)
        .add(3.0, false, true)
        .add(4.0, true, false)
        .add(5.0, true, true);

    cleanup_single_file(outFile1);
    cleanup_single_file(outFile2);
    run_single_file_simple_topology_change_with_multiple_output(model, params1, params2);
    cleanup_single_file(outFile1);
    cleanup_single_file(outFile2);
  }

  TEST(TestDynamicWrite, same_model_triggers_same_modification_for_all_observers)
  {
    Ioss::Init::Initializer io;
    Ioss::ParallelUtils     util(Ioss::ParallelUtils::comm_world());

    std::string outFile1("sameModelManyBlocks1.g");
    std::string outFile2("sameModelManyBlocks2.g");
    std::string elemFieldName("elem_field");
    std::string model("same-model");

    auto broker = Ioss::DynamicTopologyBroker::broker();
    broker->register_model(model);

    std::string file1 =
        Ioss::Utils::decode_filename(outFile1, util.parallel_rank(), util.parallel_size());
    std::string file2 =
        Ioss::Utils::decode_filename(outFile2, util.parallel_rank(), util.parallel_size());

    unlink(file1.c_str());
    unlink(file2.c_str());

    int numBlocks = util.parallel_size();
    if (numBlocks > 1)
      GTEST_SKIP();

    std::string meshDesc = get_many_block_mesh_desc(numBlocks);

    Ioss::PropertyManager propertyManager;

    Ioss::DatabaseIO *i_database = Ioss::IOFactory::create(
        "textmesh", meshDesc, Ioss::READ_MODEL, Ioss::ParallelUtils::comm_world(), propertyManager);
    Ioss::Region i_region(i_database, "input_model");
    EXPECT_TRUE(i_database != nullptr);
    EXPECT_TRUE(i_database->ok(true));

    {
      propertyManager.add(Ioss::Property("ENABLE_FILE_GROUPS", 1));
      propertyManager.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND_GROUP));

      Ioss::DatabaseIO *o_database1 =
          Ioss::IOFactory::create("exodus", outFile1, Ioss::WRITE_RESULTS,
                                  Ioss::ParallelUtils::comm_world(), propertyManager);
      Ioss::Region o_region1(o_database1, "region1");
      EXPECT_TRUE(o_database1 != nullptr);
      EXPECT_TRUE(o_database1->ok(true));
      EXPECT_TRUE(o_database1->supports_internal_change_set());

      auto fileControlOption = Ioss::FileControlOption::CONTROL_AUTO_GROUP_FILE;
      auto observer1 = std::make_shared<Observer>(i_region, elemFieldName, fileControlOption);
      broker->register_observer(model, observer1, o_region1);

      Ioss::DatabaseIO *o_database2 =
          Ioss::IOFactory::create("exodus", outFile2, Ioss::WRITE_RESULTS,
                                  Ioss::ParallelUtils::comm_world(), propertyManager);
      Ioss::Region o_region2(o_database2, "region2");
      EXPECT_TRUE(o_database2 != nullptr);
      EXPECT_TRUE(o_database2->ok(true));
      EXPECT_TRUE(o_database2->supports_internal_change_set());

      auto observer2 = std::make_shared<Observer>(i_region, elemFieldName, fileControlOption);
      broker->register_observer(model, observer2, o_region2);

      EXPECT_EQ(Ioss::TOPOLOGY_SAME, observer1->get_topology_modification());
      EXPECT_EQ(Ioss::TOPOLOGY_SAME, observer2->get_topology_modification());

      observer1->set_topology_modification(Ioss::TOPOLOGY_UNKNOWN);

      EXPECT_EQ(Ioss::TOPOLOGY_UNKNOWN, observer1->get_topology_modification());
      EXPECT_EQ(Ioss::TOPOLOGY_UNKNOWN, observer2->get_topology_modification());
    }

    unlink(file1.c_str());
    unlink(file2.c_str());
  }

  void test_single_file_simple_topology_change_data(Ioss::Region      &i_region,
                                                    const std::string &elemFieldName, int gold_step,
                                                    double gold_time)
  {
    i_region.begin_state(gold_step);
    for (Ioss::ElementBlock *i_eb : i_region.get_element_blocks()) {
      size_t num_elem = i_eb->entity_count();

      std::vector<double> field_data(num_elem);
      std::vector<int>    elem_ids;

      i_eb->get_field_data(elemFieldName, field_data);
      i_eb->get_field_data("ids", elem_ids);

      for (size_t i = 0; i < elem_ids.size(); i++) {
        double gold_value = (double)elem_ids[i] + 100 * gold_time;
        EXPECT_NEAR(gold_value, field_data[i], 1.0e-6);
      }
    }
  }

  void read_and_test_single_file_simple_topology_change(const OutputParams &params)
  {
    Ioss::PropertyManager propertyManager;
    auto i_region = construct_region(params, propertyManager, Ioss::READ_RESTART, "input_model");
    Ioss::DatabaseIO *i_database = i_region->get_database();

    test_internal_file_change_set_names(i_database);

    EXPECT_TRUE(i_database->supports_internal_change_set());

    auto numSteps = params.output_steps.size();

    int numMods = 0;

    int maxStep = 1;

    double minTime = numSteps > 0 ? params.output_times[0] : 0.0;
    double maxTime = numSteps > 0 ? params.output_times[0] : 0.0;

    bool doneOutputAfterModification = true;

    Ioss::NameList names = i_database->internal_change_set_describe(false);

    for (size_t i = 0; i < numSteps; i++) {
      maxTime = params.output_times[i];

      for (size_t j = i + 1; j < numSteps; j++) {
        if (params.modification_steps[j]) {
          maxTime = params.output_times[j - 1];
          break;
        }

        maxTime = params.output_times[j];
      }

      if (params.modification_steps[i]) {
        numMods++;
        maxStep = 1;

        EXPECT_TRUE(i_region->load_internal_change_set_mesh(names[numMods]));

        doneOutputAfterModification = false;
      }

      if (params.output_steps[i]) {
        if (!doneOutputAfterModification) {
          minTime = params.output_times[i];
        }
        auto min_result = i_region->get_min_time();
        EXPECT_EQ(1, min_result.first);
        EXPECT_NEAR(minTime, min_result.second, 1.0e-6);
        test_single_file_simple_topology_change_data(*i_region, params.elemFieldName, 1, minTime);

        if ((((i + 1) < numSteps) && params.modification_steps[i + 1]) || (i == (numSteps - 1))) {
          auto max_result = i_region->get_max_time();
          EXPECT_EQ(maxStep, max_result.first);
          EXPECT_NEAR(maxTime, max_result.second, 1.0e-6);

          test_single_file_simple_topology_change_data(*i_region, params.elemFieldName, maxStep,
                                                       maxTime);
        }

        maxStep++;
        doneOutputAfterModification = true;
      }
    }
  }

  TEST(TestDynamicRead, single_file_simple_topology_modification)
  {
    std::string outFile("singleFileManyBlocks.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    params.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, true, false)
        .add(3.0, true, true)
        .add(4.0, true, true)
        .add(5.0, true, false);

    cleanup_single_file(outFile);
    run_single_file_simple_topology_change(params);
    read_and_test_single_file_simple_topology_change(params);
    cleanup_single_file(outFile);
  }

  std::tuple<std::string, int, double> read_and_locate_db_state(const OutputParams &params,
                                                                double              targetTime)
  {
    Ioss::PropertyManager propertyManager;
    auto region = construct_region(params, propertyManager, Ioss::READ_RESTART, "input_model");
    return region->locate_db_state(targetTime);
  }

  void run_single_file_locate_db_time_state(const std::string &outFile, const OutputParams &params,
                                            double targetTime, const std::string &goldSet,
                                            int goldState, double goldTime)
  {
    cleanup_single_file(outFile);
    run_single_file_simple_topology_change(params);

    std::string changeSet;
    int         nearestState;
    double      nearestTime;

    std::tie(changeSet, nearestState, nearestTime) = read_and_locate_db_state(params, targetTime);

    EXPECT_EQ(goldSet, changeSet);
    EXPECT_EQ(goldState, nearestState);
    EXPECT_EQ(goldTime, nearestTime);

    cleanup_single_file(outFile);
  }

  TEST(TestDynamicRead, single_file_locate_db_time_state)
  {
    std::string outFile("singleFileManyBlocksPositiveTime.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    params.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, true, false)
        .add(3.0, true, true)
        .add(4.0, true, true)
        .add(5.0, true, false);

    double targetTime = 3.5;

    std::string goldSet   = Ioss::DynamicTopologyFileControl::get_internal_file_change_set_name(3);
    int         goldState = 1;
    double      goldTime  = 3.0;

    run_single_file_locate_db_time_state(outFile, params, targetTime, goldSet, goldState, goldTime);
  }

  TEST(TestDynamicRead, single_file_locate_db_time_state_all_negative_time)
  {
    std::string outFile("singleFileManyBlocksNegativeTime.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    params.add(-5.0, true, false)
        .add(-4.0, true, true)
        .add(-3.0, true, false)
        .add(-2.0, true, true)
        .add(-1.0, true, true)
        .add(-0.0, true, false);

    double targetTime = -1.5;

    std::string goldSet   = Ioss::DynamicTopologyFileControl::get_internal_file_change_set_name(4);
    int         goldState = 1;
    double      goldTime  = -1.0;

    run_single_file_locate_db_time_state(outFile, params, targetTime, goldSet, goldState, goldTime);
  }

  void run_multi_file_locate_db_time_state(const std::string &outFile, const OutputParams &params,
                                           double targetTime, const std::string &goldFile,
                                           int goldState, double goldTime)
  {
    if (params.cyclicCount > 0) {
      cleanup_cyclic_multi_files(outFile, params.cyclicCount);
    }
    else {
      cleanup_linear_multi_files(outFile, params.output_times.size());
    }

    run_multi_file_topology_change(params);

    std::string file;
    int         nearestState;
    double      nearestTime;

    std::tie(file, nearestState, nearestTime) = read_and_locate_db_state(params, targetTime);

    EXPECT_EQ(goldFile, file);
    EXPECT_EQ(goldState, nearestState);
    EXPECT_EQ(goldTime, nearestTime);

    if (params.cyclicCount > 0) {
      cleanup_cyclic_multi_files(outFile, params.cyclicCount);
    }
    else {
      cleanup_linear_multi_files(outFile, params.output_times.size());
    }
  }

  TEST(TestDynamicRead, linear_multi_file_locate_db_time_state)
  {
    std::string outFile("linearMultiFileManyBlocksPositiveTime.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    params.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, true, false)
        .add(3.0, true, true)
        .add(4.0, true, true)
        .add(5.0, true, false);

    double targetTime = 3.5;

    std::string goldFile =
        Ioss::DynamicTopologyFileControl::get_linear_database_filename(outFile, 3);
    int    goldState = 1;
    double goldTime  = 3.0;

    run_multi_file_locate_db_time_state(outFile, params, targetTime, goldFile, goldState, goldTime);
  }

  TEST(TestDynamicRead, cyclic_multi_file_locate_db_time_state)
  {
    std::string outFile("cyclicMultiFileManyBlocksPositiveTime.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    std::vector<double> output_times{0.0, 0.5, 1.5, 1.75, 2.0, 3.0};
    std::vector<bool>   output_steps{true, true, true, true, true, true};
    std::vector<bool>   modification_steps{false, true, false, true, true, false};
    //                                      -A     -B    -B     -C    -B    -B

    params.set_data(output_times, output_steps, modification_steps);
    params.set_cyclic_count(3);

    double targetTime = 1.95;

    std::string goldFile = Ioss::DynamicTopologyFileControl::get_cyclic_database_filename(
        outFile, params.cyclicCount, 2);
    int    goldState = 1;
    double goldTime  = 2.0;

    run_multi_file_locate_db_time_state(outFile, params, targetTime, goldFile, goldState, goldTime);
  }

  std::tuple<std::string, int, double> read_and_locate_db_max_time(const OutputParams &params)
  {
    Ioss::PropertyManager propertyManager;

    std::string outFile = params.outFile;

    if (params.cyclicCount > 0) {
      outFile = Ioss::DynamicTopologyFileControl::get_cyclic_database_filename(
          params.outFile, params.cyclicCount, 0);
      propertyManager.add(Ioss::Property("base_filename", params.outFile));
    }

    Ioss::DatabaseIO *db = Ioss::IOFactory::create(
        "exodus", outFile, Ioss::READ_RESTART, Ioss::ParallelUtils::comm_world(), propertyManager);

    Ioss::Region region(db, "input_model");
    region.set_file_cyclic_count(params.cyclicCount);

    EXPECT_TRUE(db != nullptr);
    EXPECT_TRUE(db->ok(true));

    return region.get_db_max_time();
  }

  TEST(TestDynamicRead, single_file_locate_db_max_time)
  {
    std::string outFile("singleFileManyBlocksMaxTime.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    params.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, true, false)
        .add(3.0, true, true)
        .add(4.0, true, true)
        .add(5.0, true, false);

    cleanup_single_file(outFile);
    run_single_file_simple_topology_change(params);

    std::string maxSet;
    int         maxState;
    double      maxTime;

    std::tie(maxSet, maxState, maxTime) = read_and_locate_db_max_time(params);

    std::string goldSet   = Ioss::DynamicTopologyFileControl::get_internal_file_change_set_name(4);
    int         goldState = 2;
    double      goldTime  = 5.0;

    EXPECT_EQ(goldSet, maxSet);
    EXPECT_EQ(goldState, maxState);
    EXPECT_EQ(goldTime, maxTime);

    cleanup_single_file(outFile);
  }

  TEST(TestDynamicRead, linear_multi_file_locate_db_max_time)
  {
    std::string outFile("linearMultiFileManyBlocksMaxTime.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    params.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, true, false)
        .add(3.0, true, true)
        .add(4.0, true, true)
        .add(5.0, true, false);

    cleanup_linear_multi_files(outFile, params.output_times.size());
    run_multi_file_topology_change(params);

    std::string maxFile;
    int         maxState;
    double      maxTime;

    std::tie(maxFile, maxState, maxTime) = read_and_locate_db_max_time(params);

    std::string goldFile =
        Ioss::DynamicTopologyFileControl::get_linear_database_filename(outFile, 4);
    int    goldState = 2;
    double goldTime  = 5.0;

    EXPECT_EQ(goldFile, maxFile);
    EXPECT_EQ(goldState, maxState);
    EXPECT_EQ(goldTime, maxTime);

    cleanup_linear_multi_files(outFile, params.output_times.size());
  }

  TEST(TestDynamicRead, cyclic_multi_file_locate_db_max_time)
  {
    std::string outFile("cyclicMultiFileManyBlocksMaxTime.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    std::vector<double> output_times{0.0, 0.5, 1.5, 1.75, 2.0, 3.0};
    std::vector<bool>   output_steps{true, true, true, true, true, true};
    std::vector<bool>   modification_steps{false, true, false, true, true, false};
    //                                      -A     -B    -B     -C    -B    -B

    params.set_data(output_times, output_steps, modification_steps);
    params.set_cyclic_count(3);

    cleanup_cyclic_multi_files(outFile, params.cyclicCount);
    run_multi_file_topology_change(params);

    std::string maxFile;
    int         maxState;
    double      maxTime;

    std::tie(maxFile, maxState, maxTime) = read_and_locate_db_max_time(params);

    std::string goldFile = Ioss::DynamicTopologyFileControl::get_cyclic_database_filename(
        outFile, params.cyclicCount, 2); // -B
    int    goldState = 2;
    double goldTime  = 3.0;

    EXPECT_EQ(goldFile, maxFile);
    EXPECT_EQ(goldState, maxState);
    EXPECT_EQ(goldTime, maxTime);

    cleanup_cyclic_multi_files(outFile, params.cyclicCount);
  }

  std::tuple<std::string, int, double> read_and_locate_db_min_time(const OutputParams &params)
  {
    Ioss::PropertyManager propertyManager;
    auto region = construct_region(params, propertyManager, Ioss::READ_RESTART, "input_model");

    return region->get_db_min_time();
  }

  TEST(TestDynamicRead, single_file_locate_db_min_time)
  {
    std::string outFile("singleFileManyBlocksMinTime.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    params.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, true, false)
        .add(3.0, true, true)
        .add(4.0, true, true)
        .add(5.0, true, false);

    cleanup_single_file(outFile);
    run_single_file_simple_topology_change(params);

    std::string minSet;
    int         minState;
    double      minTime;

    std::tie(minSet, minState, minTime) = read_and_locate_db_min_time(params);

    std::string goldSet   = Ioss::DynamicTopologyFileControl::get_internal_file_change_set_name(1);
    int         goldState = 1;
    double      goldTime  = 0.0;

    EXPECT_EQ(goldSet, minSet);
    EXPECT_EQ(goldState, minState);
    EXPECT_EQ(goldTime, minTime);

    cleanup_single_file(outFile);
  }

  TEST(TestDynamicRead, linear_multi_file_locate_db_min_time)
  {
    std::string outFile("linearMultiFileManyBlocksMinTime.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    params.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, true, false)
        .add(3.0, true, true)
        .add(4.0, true, true)
        .add(5.0, true, false);

    cleanup_linear_multi_files(outFile, params.output_times.size());
    run_multi_file_topology_change(params);

    std::string minFile;
    int         minState;
    double      minTime;

    std::tie(minFile, minState, minTime) = read_and_locate_db_min_time(params);

    std::string goldFile =
        Ioss::DynamicTopologyFileControl::get_linear_database_filename(outFile, 1);
    int    goldState = 1;
    double goldTime  = 0.0;

    EXPECT_EQ(goldFile, minFile);
    EXPECT_EQ(goldState, minState);
    EXPECT_EQ(goldTime, minTime);

    cleanup_linear_multi_files(outFile, params.output_times.size());
  }

  TEST(TestDynamicRead, cyclic_multi_file_locate_db_min_time)
  {
    std::string outFile("cyclicMultiFileManyBlocksMinTime.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    std::vector<double> output_times{0.0, 0.5, 1.5, 1.75, 2.0, 3.0};
    std::vector<bool>   output_steps{true, true, true, true, true, true};
    std::vector<bool>   modification_steps{false, true, false, true, true, false};
    //                                      -A     -B    -B     -C    -B    -B

    params.set_data(output_times, output_steps, modification_steps);
    params.set_cyclic_count(3);

    cleanup_cyclic_multi_files(outFile, params.cyclicCount);
    run_multi_file_topology_change(params);

    std::string minFile;
    int         minState;
    double      minTime;

    std::tie(minFile, minState, minTime) = read_and_locate_db_min_time(params);

    std::string goldFile = Ioss::DynamicTopologyFileControl::get_cyclic_database_filename(
        outFile, params.cyclicCount, 1);
    int    goldState = 1;
    double goldTime  = 0.0;

    EXPECT_EQ(goldFile, minFile);
    EXPECT_EQ(goldState, minState);
    EXPECT_EQ(goldTime, minTime);

    cleanup_cyclic_multi_files(outFile, params.cyclicCount);
  }

  unsigned get_num_change_sets(const OutputParams &params)
  {
    auto numSteps = params.output_steps.size();

    unsigned numMods   = 0;
    int      numModInc = 0;

    int currentModStep = -1;

    for (size_t i = 0; i < numSteps; i++) {
      if (params.modification_steps[i]) {
        currentModStep = i;
        numMods++;
      }

      if (params.output_steps[i] && (currentModStep < 0)) {
        numModInc = 1;
      }
    }

    numMods += numModInc;

    return numMods;
  }

  void read_and_test_single_file_topology_change_set(const OutputParams &params)
  {
    Ioss::PropertyManager propertyManager;
    auto i_region = construct_region(params, propertyManager, Ioss::READ_RESTART, "input_model");

    auto changeSet = Ioss::ChangeSetFactory::create(i_region.get());
    changeSet->populate_change_sets();

    EXPECT_EQ(Ioss::CHANGE_SET_INTERNAL_FILES, changeSet->database_format());

    std::vector<std::string> gold_names;
    std::vector<std::string> gold_full_names;

    unsigned numMods = get_num_change_sets(params);

    fill_internal_file_change_set_gold_names(numMods, gold_names, gold_full_names);

    EXPECT_EQ(gold_names, changeSet->names());
  }

  TEST(TestChangeSet, single_file_simple_topology_modification)
  {
    std::string outFile("singleFileManyBlocksChangeSet.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    params.add(0.0, true, false)
        .add(1.0, true, true)
        .add(2.0, true, false)
        .add(3.0, true, true)
        .add(4.0, true, true)
        .add(5.0, true, false);

    cleanup_single_file(outFile);
    run_single_file_simple_topology_change(params);
    read_and_test_single_file_topology_change_set(params);
    cleanup_single_file(outFile);
  }

  void read_and_test_cyclic_multi_file_topology_change_set(const OutputParams &params)
  {
    ASSERT_TRUE(params.cyclicCount > 0);

    Ioss::PropertyManager propertyManager;
    auto i_region = construct_region(params, propertyManager, Ioss::READ_RESTART, "input_model");

    auto changeSet = Ioss::ChangeSetFactory::create(i_region.get());
    changeSet->populate_change_sets();

    EXPECT_EQ(Ioss::CHANGE_SET_CYCLIC_MULTI_FILES, changeSet->database_format());

    std::vector<std::string> gold_names;

    unsigned numMods = get_num_change_sets(params);
    if (numMods > params.cyclicCount) {
      numMods = params.cyclicCount;
    }

    for (unsigned i = 1; i <= numMods; i++) {
      gold_names.push_back(Ioss::DynamicTopologyFileControl::get_cyclic_database_filename(
          params.outFile, params.cyclicCount, i));
    }

    EXPECT_EQ(gold_names, changeSet->names());
  }

  TEST(TestChangeSet, multi_file_cyclic_topology_modification)
  {
    std::string outFile("cyclicMultiFileManyBlocksChangeSet.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    std::vector<double> output_times{0.0, 0.5, 1.5, 1.75, 2.0, 3.0};
    std::vector<bool>   output_steps{true, true, true, true, true, true};
    std::vector<bool>   modification_steps{false, true, false, true, true, false};

    params.set_data(output_times, output_steps, modification_steps);
    params.set_cyclic_count(3);

    cleanup_cyclic_multi_files(outFile, params.cyclicCount);
    run_multi_file_topology_change(params);
    read_and_test_cyclic_multi_file_topology_change_set(params);
    cleanup_cyclic_multi_files(outFile, params.cyclicCount);
  }

  void read_and_test_linear_multi_file_topology_change_set(const OutputParams &params)
  {
    ASSERT_TRUE(params.cyclicCount == 0);

    Ioss::PropertyManager propertyManager;
    auto i_region = construct_region(params, propertyManager, Ioss::READ_RESTART, "input_model");

    auto changeSet = Ioss::ChangeSetFactory::create(i_region.get());
    changeSet->populate_change_sets();

    EXPECT_EQ(Ioss::CHANGE_SET_LINEAR_MULTI_FILES, changeSet->database_format());

    std::vector<std::string> gold_names;

    unsigned numMods = get_num_change_sets(params);

    for (unsigned i = 1; i <= numMods; i++) {
      gold_names.push_back(
          Ioss::DynamicTopologyFileControl::get_linear_database_filename(params.outFile, i));
    }

    EXPECT_EQ(gold_names, changeSet->names());
  }

  TEST(TestChangeSet, multi_file_linear_topology_modification)
  {
    std::string outFile("linearMultiFileManyBlocksChangeSet.g");
    std::string elemFieldName = "elem_field";

    OutputParams params(outFile, elemFieldName);

    std::vector<double> output_times{0.0, 0.5, 1.5, 1.75, 2.0, 3.0};
    std::vector<bool>   output_steps{true, true, true, true, true, true};
    std::vector<bool>   modification_steps{false, true, false, true, true, false};

    params.set_data(output_times, output_steps, modification_steps);

    cleanup_linear_multi_files(outFile, params.output_times.size());
    run_multi_file_topology_change(params);
    read_and_test_linear_multi_file_topology_change_set(params);
    cleanup_linear_multi_files(outFile, params.output_times.size());
  }
} // namespace
