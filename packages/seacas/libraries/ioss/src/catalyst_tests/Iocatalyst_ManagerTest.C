// Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_ParallelUtils.h>
#include <catalyst/Iocatalyst_CatalystManager.h>
#include <catalyst_tests/Iocatalyst_DatabaseIOTest.h>
#include <catalyst_tests/Iocatalyst_LoggingTest.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>

using namespace Iocatalyst;

class ManagerTest : public ::testing::Test
{
protected:
  Ioss::PropertyManager               props;
  Ioss::ParallelUtils                 putils;
  CatalystManager::CatalystProps      catalystProps;
  CatalystManager::CatalystPipelineID id;
  conduit_cpp::Node                   n;
  void                                reset() { CatalystManager::getInstance().reset(); }
  void                                initialize()
  {
    id            = CatalystManager::getInstance().initialize(props, putils);
    catalystProps = CatalystManager::getInstance().getCatalystProps(id);
  }

  void compareConduit(const conduit_cpp::Node &n, const conduit_cpp::Node &m)
  {
    EXPECT_EQ(n.to_string(), m.to_string());
  }

  void checkExecuteProps(const conduit_cpp::Node &n, CatalystManager::CatalystProps &p, int state,
                         double time)
  {
    auto csp = CatalystManager::getInstance().getCatStatePath();
    EXPECT_EQ(n[csp + CatalystManager::TIMESTEP].as_int(), state - 1);
    EXPECT_EQ(n[csp + CatalystManager::CYCLE].as_int(), state - 1);
    EXPECT_EQ(n[csp + CatalystManager::TIME].as_double(), time);
    auto sn = std::to_string(p.catalystPipelineID);
    EXPECT_EQ(n[csp + CatalystManager::PIPELINES + CatalystManager::FS + sn].as_string(), sn);
  }

  void checkExecuteData(const conduit_cpp::Node &n, std::string &channel_name, int state,
                        double time, const conduit_cpp::Node &m)
  {
    auto ccp = CatalystManager::getInstance().getCatChannelsPath();
    auto ip  = ccp + channel_name + CatalystManager::FS;
    EXPECT_EQ(n[ip + CatalystManager::TYPE].as_string(), CatalystManager::IOSS);
    auto              dp = ip + CatalystManager::FS + CatalystManager::DATA;
    conduit_cpp::Node d;
    d.set_node(n[dp]);
    d.remove(CatalystManager::TIMESTEP);
    d.remove(CatalystManager::CYCLE);
    d.remove(CatalystManager::TIME);
    compareConduit(d, m);
    EXPECT_DOUBLE_EQ(n[dp + CatalystManager::FS + CatalystManager::TIME].as_double(), time);
    EXPECT_EQ(n[dp + CatalystManager::FS + CatalystManager::CYCLE].as_int(), state - 1);
    EXPECT_EQ(n[dp + CatalystManager::FS + CatalystManager::TIMESTEP].as_int(), state - 1);
  }
};

TEST_F(LoggingTest, LoggingDefault)
{
  Ioss::ParallelUtils putils;
  CatalystManager::getInstance().writeToCatalystLogFile(putils, props);
  EXPECT_FALSE(isFileExists(CatalystLogging::getDefaultLogFileName().c_str()));
}

TEST_F(LoggingTest, LoggingNotEnabled)
{
  Ioss::ParallelUtils putils;
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", false));
  CatalystManager::getInstance().writeToCatalystLogFile(putils, props);
  EXPECT_FALSE(isFileExists(CatalystLogging::getDefaultLogFileName().c_str()));
}

TEST_F(LoggingTest, LoggingEnabled)
{
  Ioss::ParallelUtils putils;
  props.add(Ioss::Property("CATALYST_LOGGING_ENABLED", true));
  props.add(Ioss::Property("CATALYST_LOGGING_STRING_PROP", "foo"));
  props.add(Ioss::Property("CATALYST_LOGGING_INTEGER_PROP", 6));
  props.add(Ioss::Property("CATALYST_LOGGING_REAL_PROP", 3.7556));
  CatalystManager::getInstance().writeToCatalystLogFile(putils, props);
  EXPECT_TRUE(isFileExists(CatalystLogging::getDefaultLogFileName().c_str()));
}

TEST_F(ManagerTest, CatalystPipelineID)
{
  initialize();
  EXPECT_EQ(catalystProps.catalystPipelineID, 0);

  initialize();
  EXPECT_EQ(catalystProps.catalystPipelineID, 1);

  initialize();
  EXPECT_EQ(catalystProps.catalystPipelineID, 2);
}

TEST_F(ManagerTest, CATALYST_BLOCK_PARSE_JSON_STRING)
{
  std::string jsonScript = "{foo: 12}";
  props.add(Ioss::Property("CATALYST_BLOCK_PARSE_JSON_STRING", jsonScript));
  initialize();
  EXPECT_EQ(catalystProps.catalystBlockJSON, jsonScript);
}

TEST_F(ManagerTest, PHACTORI_JSON_SCRIPT)
{
  std::string   jsonFileName = "jsonFile.json";
  std::string   jsonScript   = "{foo: 12}";
  std::ofstream jsonFile;
  jsonFile.open(jsonFileName);
  jsonFile << jsonScript;
  jsonFile.close();
  props.add(Ioss::Property(CatalystManager::PHACTORI_JSON_SCRIPT, jsonFileName));
  initialize();
  EXPECT_EQ(catalystProps.catalystBlockJSON, jsonScript);
  remove(jsonFileName.c_str());
}

TEST_F(ManagerTest, CATALYST_SCRIPT)
{
  initialize();
  EXPECT_EQ(catalystProps.catalystPythonFilename,
            CatalystManager::getInstance().getCatalystPythonDriverPath());

  std::string catalystFileName = "/path/to/file/catalystFile.txt";
  props.add(Ioss::Property(CatalystManager::CATALYST_SCRIPT, catalystFileName));
  initialize();
  EXPECT_EQ(catalystProps.catalystPythonFilename, catalystFileName);
}

TEST_F(ManagerTest, CATALYST_SCRIPT_EXTRA_FILE)
{
  std::string extraFileName = "extraFileName.txt";
  props.add(Ioss::Property(CatalystManager::CATALYST_SCRIPT_EXTRA_FILE, extraFileName));
  initialize();
  EXPECT_EQ(catalystProps.catalystScriptExtraFile, extraFileName);
}

TEST_F(ManagerTest, CATALYST_BLOCK_PARSE_INPUT_DECK_NAME)
{
  std::string inputDeckName = "contact.i";
  props.add(Ioss::Property(CatalystManager::CATALYST_BLOCK_PARSE_INPUT_DECK_NAME, inputDeckName));
  initialize();
  EXPECT_EQ(catalystProps.catalystInputDeckName, inputDeckName);
}

TEST_F(ManagerTest, CATALYST_ENABLE_LOGGING)
{
  initialize();
  EXPECT_FALSE(catalystProps.enableLogging);

  props.add(Ioss::Property(CatalystManager::CATALYST_ENABLE_LOGGING, true));
  initialize();
  EXPECT_TRUE(catalystProps.enableLogging);
}

TEST_F(ManagerTest, CATALYST_DEBUG_LEVEL)
{
  initialize();
  EXPECT_EQ(catalystProps.debugLevel, 0);

  props.add(Ioss::Property(CatalystManager::CATALYST_DEBUG_LEVEL, 3));
  initialize();
  EXPECT_EQ(catalystProps.debugLevel, 3);
}

TEST_F(ManagerTest, CATALYST_OUTPUT_DIRECTORY)
{
  initialize();
  EXPECT_EQ(catalystProps.catalystOutputDirectory,
            CatalystManager::getInstance().CATALYST_OUTPUT_DEFAULT);

  std::string catalystOutputDirectory = "catalyst";
  props.add(Ioss::Property(CatalystManager::CATALYST_OUTPUT_DIRECTORY, catalystOutputDirectory));
  initialize();
  EXPECT_EQ(catalystProps.catalystOutputDirectory, catalystOutputDirectory);
}

TEST_F(ManagerTest, CATALYST_INPUT_NAME)
{
  initialize();
  EXPECT_EQ(catalystProps.catalystInputName, CatalystManager::getInstance().CATALYST_INPUT_DEFAULT);

  std::string catalystInputName = "mesh";
  props.add(Ioss::Property(CatalystManager::CATALYST_INPUT_NAME, catalystInputName));
  initialize();
  EXPECT_EQ(catalystProps.catalystInputName, catalystInputName);
}

TEST_F(ManagerTest, CATALYST_MULTI_INPUT_PIPELINE_NAME)
{
  initialize();
  EXPECT_FALSE(catalystProps.enableCatalystMultiInputPipeline);

  std::string catalystMultiInputPipelineName = "multi";
  props.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME,
                           catalystMultiInputPipelineName));
  initialize();
  EXPECT_EQ(catalystProps.catalystMultiInputPipelineName, catalystMultiInputPipelineName);
  EXPECT_TRUE(catalystProps.enableCatalystMultiInputPipeline);
}

TEST_F(ManagerTest, InitializeConduitDefault)
{
  reset();
  compareConduit(CatalystManager::getInstance().getInitializeConduit(), n);
}

TEST_F(ManagerTest, InitializeConduitCatalystFile)
{
  reset();
  std::string catalystFileName = "/path/to/file/catalystFile.txt";
  props.add(Ioss::Property(CatalystManager::CATALYST_SCRIPT, catalystFileName));
  initialize();
  CatalystManager::getInstance().addScriptProps(n, catalystProps);
  compareConduit(CatalystManager::getInstance().getInitializeConduit(), n);
}

TEST_F(ManagerTest, InitializeConduitPhactoriJSON)
{
  reset();
  std::string js = "some json";
  props.add(Ioss::Property(CatalystManager::CATALYST_BLOCK_PARSE_JSON_STRING, js));
  props.add(Ioss::Property(CatalystManager::CATALYST_INPUT_NAME, "data"));
  props.add(Ioss::Property(CatalystManager::CATALYST_SCRIPT_EXTRA_FILE, "extra"));
  props.add(Ioss::Property(CatalystManager::CATALYST_BLOCK_PARSE_INPUT_DECK_NAME, "adagio"));
  props.add(Ioss::Property(CatalystManager::CATALYST_OUTPUT_DIRECTORY, "temp"));
  props.add(Ioss::Property(CatalystManager::CATALYST_ENABLE_LOGGING, true));
  props.add(Ioss::Property(CatalystManager::CATALYST_DEBUG_LEVEL, 11));
  initialize();
  CatalystManager::getInstance().addScriptProps(n, catalystProps);
  compareConduit(CatalystManager::getInstance().getInitializeConduit(), n);
}

TEST_F(ManagerTest, InitializeConduitMultipleScripts)
{
  reset();
  std::string catalystFileName = "/path/to/file/catalystFile.txt";
  props.add(Ioss::Property(CatalystManager::CATALYST_SCRIPT, catalystFileName));
  initialize();

  Ioss::PropertyManager propsOne;
  std::string           otherFile = "/path/to/other/file";
  props.add(Ioss::Property(CatalystManager::CATALYST_SCRIPT, otherFile));
  auto idOne       = CatalystManager::getInstance().initialize(propsOne, putils);
  auto catPropsOne = CatalystManager::getInstance().getCatalystProps(idOne);

  Ioss::PropertyManager propsTwo;
  std::string           js = "json";
  props.add(Ioss::Property(CatalystManager::CATALYST_BLOCK_PARSE_JSON_STRING, js));
  idOne            = CatalystManager::getInstance().initialize(propsOne, putils);
  auto catPropsTwo = CatalystManager::getInstance().getCatalystProps(idOne);

  CatalystManager::getInstance().addScriptProps(n, catalystProps);
  CatalystManager::getInstance().addScriptProps(n, catPropsOne);
  CatalystManager::getInstance().addScriptProps(n, catPropsTwo);
  compareConduit(CatalystManager::getInstance().getInitializeConduit(), n);
}

TEST_F(ManagerTest, ExecuteConduitOneScript)
{
  reset();
  std::string catalystFileName = "/path/to/file/catalystFile.txt";
  props.add(Ioss::Property(CatalystManager::CATALYST_SCRIPT, catalystFileName));
  initialize();

  conduit_cpp::Node m;
  int               state = 2;
  double            time  = 10.4;
  m["some/data"]          = 32;
  CatalystManager::getInstance().addExecuteProps(n, catalystProps, state, time);
  CatalystManager::getInstance().addExecuteData(n, catalystProps.catalystInputName, state, time, m);
  checkExecuteProps(n, catalystProps, state, time);
  checkExecuteData(n, catalystProps.catalystInputName, state, time, m);
}

TEST_F(ManagerTest, ExecuteConduitInputName)
{
  reset();
  props.add(Ioss::Property(CatalystManager::CATALYST_INPUT_NAME, "dataset"));
  initialize();

  int               state = 10;
  double            time  = 4.5;
  conduit_cpp::Node m;
  m["other/data"] = 90;
  CatalystManager::getInstance().addExecuteProps(n, catalystProps, state, time);
  CatalystManager::getInstance().addExecuteData(n, catalystProps.catalystInputName, state, time, m);
  checkExecuteProps(n, catalystProps, state, time);
  checkExecuteData(n, catalystProps.catalystInputName, state, time, m);
}

TEST_F(ManagerTest, ManagerStateDefault)
{
  reset();
  EXPECT_EQ(CatalystManager::getInstance().getManagerState(), CatalystManager::mInit);
}

TEST_F(ManagerTest, InvalidIDGetCatalystProps)
{
  reset();
  EXPECT_THROW(CatalystManager::getInstance().getCatalystProps(1), std::runtime_error);
}

TEST_F(ManagerTest, ManagerExecuteStateChange)
{
  reset();
  initialize();
  EXPECT_EQ(CatalystManager::getInstance().getManagerState(), CatalystManager::mInit);

  conduit_cpp::Node m;
  CatalystManager::getInstance().execute(catalystProps.catalystPipelineID, 2, 10.2, m);
  EXPECT_EQ(CatalystManager::getInstance().getManagerState(), CatalystManager::mExecute);

  EXPECT_THROW(CatalystManager::getInstance().initialize(props, putils), std::runtime_error);
}

TEST_F(ManagerTest, ManagerPipelineState)
{
  reset();
  initialize();
  EXPECT_EQ(CatalystManager::getInstance().getPipelineState(catalystProps.catalystPipelineID),
            CatalystManager::pExecute);
}

TEST_F(ManagerTest, ManagerFinalizeStateChange)
{
  reset();
  initialize();
  CatalystManager::getInstance().finalize(catalystProps.catalystPipelineID);
  EXPECT_EQ(CatalystManager::getInstance().getManagerState(), CatalystManager::mFinalize);
  EXPECT_EQ(CatalystManager::getInstance().getPipelineState(catalystProps.catalystPipelineID),
            CatalystManager::pFinalize);
  EXPECT_THROW(CatalystManager::getInstance().initialize(props, putils), std::runtime_error);
  conduit_cpp::Node m;
  EXPECT_THROW(CatalystManager::getInstance().execute(catalystProps.catalystPipelineID, 2, 10.2, m),
               std::runtime_error);
}

TEST_F(ManagerTest, ManagerFinalizeStateChangeMultipleScripts)
{
  reset();
  initialize();

  auto idOne = CatalystManager::getInstance().initialize(props, putils);
  auto idTwo = CatalystManager::getInstance().initialize(props, putils);

  CatalystManager::getInstance().finalize(catalystProps.catalystPipelineID);
  EXPECT_EQ(CatalystManager::getInstance().getPipelineState(catalystProps.catalystPipelineID),
            CatalystManager::pFinalize);
  EXPECT_EQ(CatalystManager::getInstance().getManagerState(), CatalystManager::mExecute);

  conduit_cpp::Node m;
  CatalystManager::getInstance().execute(idOne, 2, 10.2, m);
  CatalystManager::getInstance().finalize(idOne);
  EXPECT_EQ(CatalystManager::getInstance().getPipelineState(idOne), CatalystManager::pFinalize);
  EXPECT_EQ(CatalystManager::getInstance().getManagerState(), CatalystManager::mExecute);
  EXPECT_THROW(CatalystManager::getInstance().execute(idOne, 2, 10.2, m), std::runtime_error);

  CatalystManager::getInstance().finalize(idTwo);
  EXPECT_EQ(CatalystManager::getInstance().getPipelineState(idTwo), CatalystManager::pFinalize);
  EXPECT_EQ(CatalystManager::getInstance().getManagerState(), CatalystManager::mFinalize);
}

TEST_F(ManagerTest, ManagerGetCatDataPath)
{
  EXPECT_EQ(CatalystManager::getInstance().getCatDataPath(props), "catalyst/channels/input/data");
  std::string name = "foo";
  props.add(Ioss::Property(CatalystManager::CATALYST_INPUT_NAME, name));
  EXPECT_EQ(CatalystManager::getInstance().getCatDataPath(props), "catalyst/channels/foo/data");
}

TEST_F(ManagerTest, ManagerMultiInputStateChange)
{
  reset();
  std::string catalystMultiInputPipelineName = "multi";
  props.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME,
                           catalystMultiInputPipelineName));
  initialize();
  EXPECT_EQ(catalystProps.enableCatalystMultiInputPipeline, true);
  EXPECT_EQ(catalystProps.catalystMultiInputPipelineName, catalystMultiInputPipelineName);
  n["my/data"] = 12;
  int    state = 12;
  double time  = 55.3;
  CatalystManager::getInstance().setMultiInputWaitState(id, state, time, n);
  auto p = CatalystManager::getInstance().getCatalystProps(id);
  EXPECT_EQ(p.pipelineState, CatalystManager::pWaitExecute);
  compareConduit(p.data, n);
  EXPECT_EQ(p.state, state);
  EXPECT_DOUBLE_EQ(p.time, time);
}

TEST_F(ManagerTest, ManagerSetMultiInputStateError)
{
  reset();
  initialize();
  EXPECT_THROW(CatalystManager::getInstance().setMultiInputWaitState(id, 2, 7.6, n),
               std::runtime_error);
}

TEST_F(ManagerTest, ManagerCanExecuteMultiInputScriptError)
{
  reset();
  initialize();
  EXPECT_THROW(CatalystManager::getInstance().canExecuteMultiInputScript(id), std::runtime_error);
}

TEST_F(ManagerTest, ManagerCanExecuteMultiInputScriptOneNoState)
{
  reset();
  props.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME, "multi"));
  initialize();
  EXPECT_FALSE(CatalystManager::getInstance().canExecuteMultiInputScript(id));
}

TEST_F(ManagerTest, ManagerCanExecuteMultiInputScriptOne)
{
  reset();
  props.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME, "multi"));
  initialize();
  CatalystManager::getInstance().setMultiInputWaitState(id, 3, 5.6, n);
  EXPECT_TRUE(CatalystManager::getInstance().canExecuteMultiInputScript(id));
}

TEST_F(ManagerTest, ManagerCanExecuteMultiInputScriptFour)
{
  reset();
  std::string name  = "multi";
  int         state = 2;
  double      time  = 8.9;
  props.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME, name));
  initialize();
  CatalystManager::getInstance().setMultiInputWaitState(id, state, time, n);

  Ioss::PropertyManager propsOne;
  CatalystManager::getInstance().initialize(propsOne, putils);
  EXPECT_TRUE(CatalystManager::getInstance().canExecuteMultiInputScript(id));

  Ioss::PropertyManager propsTwo;
  propsTwo.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME, name));
  auto idTwo = CatalystManager::getInstance().initialize(propsTwo, putils);
  EXPECT_FALSE(CatalystManager::getInstance().canExecuteMultiInputScript(id));
  CatalystManager::getInstance().setMultiInputWaitState(idTwo, state, time, n);
  EXPECT_TRUE(CatalystManager::getInstance().canExecuteMultiInputScript(idTwo));

  Ioss::PropertyManager propsThree;
  propsThree.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME, "foo"));
  auto idThree = CatalystManager::getInstance().initialize(propsThree, putils);
  EXPECT_FALSE(CatalystManager::getInstance().canExecuteMultiInputScript(idThree));
  CatalystManager::getInstance().setMultiInputWaitState(idThree, state, time, n);
  EXPECT_TRUE(CatalystManager::getInstance().canExecuteMultiInputScript(idThree));
  EXPECT_TRUE(CatalystManager::getInstance().canExecuteMultiInputScript(id));
}

TEST_F(ManagerTest, ManagerClearAllMultiInputWaitStatesError)
{
  reset();
  initialize();
  EXPECT_THROW(CatalystManager::getInstance().clearAllMultiInputWaitStates(id), std::runtime_error);
}

TEST_F(ManagerTest, ManagerClearAllMultiInputWaitStatesOne)
{
  reset();
  std::string name = "multi";
  props.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME, name));
  initialize();
  n["my/value"] = 13.1;
  CatalystManager::getInstance().setMultiInputWaitState(id, 3, 9.2, n);
  EXPECT_TRUE(CatalystManager::getInstance().canExecuteMultiInputScript(id));
  auto p = CatalystManager::getInstance().getCatalystProps(id);
  compareConduit(p.data, n);
  CatalystManager::getInstance().clearAllMultiInputWaitStates(id);
  EXPECT_FALSE(CatalystManager::getInstance().canExecuteMultiInputScript(id));
  p = CatalystManager::getInstance().getCatalystProps(id);
  compareConduit(p.data, conduit_cpp::Node());
}

TEST_F(ManagerTest, ManagerClearAllMultiInputWaitStatesThree)
{
  reset();
  std::string name  = "multi";
  int         state = 2;
  double      time  = 8.9;
  props.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME, name));
  initialize();
  CatalystManager::getInstance().setMultiInputWaitState(id, state, time, n);

  Ioss::PropertyManager propsOne;
  CatalystManager::getInstance().initialize(propsOne, putils);

  Ioss::PropertyManager propsTwo;
  propsTwo.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME, name));
  auto idTwo = CatalystManager::getInstance().initialize(propsTwo, putils);
  CatalystManager::getInstance().setMultiInputWaitState(idTwo, state, time, n);

  Ioss::PropertyManager propsThree;
  propsThree.add(Ioss::Property(CatalystManager::CATALYST_MULTI_INPUT_PIPELINE_NAME, "foo"));
  auto idThree = CatalystManager::getInstance().initialize(propsThree, putils);
  CatalystManager::getInstance().setMultiInputWaitState(idThree, state, time, n);

  EXPECT_TRUE(CatalystManager::getInstance().canExecuteMultiInputScript(idTwo));
  EXPECT_TRUE(CatalystManager::getInstance().canExecuteMultiInputScript(idThree));

  CatalystManager::getInstance().clearAllMultiInputWaitStates(id);
  EXPECT_FALSE(CatalystManager::getInstance().canExecuteMultiInputScript(idTwo));
  EXPECT_TRUE(CatalystManager::getInstance().canExecuteMultiInputScript(idThree));

  CatalystManager::getInstance().clearAllMultiInputWaitStates(idThree);
  EXPECT_FALSE(CatalystManager::getInstance().canExecuteMultiInputScript(idThree));
}

TEST_F(ManagerTest, ManagerAddExecuteDataThreeInputs)
{
  int               state = 6;
  double            time  = 20.2;
  conduit_cpp::Node m;
  m["other/data"] = 90;

  conduit_cpp::Node m1;
  std::string       m1Channel = "m1";
  m1["m1/data"]               = 100;

  conduit_cpp::Node m2;
  std::string       m2Channel = "m2";
  m2["m2/data"]               = 500;

  CatalystManager::getInstance().addExecuteProps(n, catalystProps, state, time);
  CatalystManager::getInstance().addExecuteData(n, catalystProps.catalystInputName, state, time, m);
  CatalystManager::getInstance().addExecuteData(n, m1Channel, state + 1, time + 1, m1);
  CatalystManager::getInstance().addExecuteData(n, m2Channel, state + 2, time + 2, m2);
  checkExecuteProps(n, catalystProps, state, time);
  checkExecuteData(n, catalystProps.catalystInputName, state, time, m);
  checkExecuteData(n, m1Channel, state + 1, time + 1, m1);
  checkExecuteData(n, m2Channel, state + 2, time + 2, m2);
}
