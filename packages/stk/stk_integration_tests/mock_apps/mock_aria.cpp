
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Constants.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/CommSplitting.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include "StkMesh.hpp"
#include <iostream>
#include <sstream>

class MockAria
{
public:
  MockAria()
  : m_appName("Mock-Aria"),
    m_mesh(),
    m_doneFlagName("time step status"),
    m_doingSendTransfer(false),
    m_doingRecvTransfer(false),
    m_sendFieldName(),
    m_recvFieldName()
  {}

  ~MockAria()
  {
    stk::parallel_machine_finalize();
  }

  double calculate_time_step()
  {
    return 0.1;
  }

  void read_input_and_setup_split_comms(int argc, char** argv)
  {
    m_commWorld = stk::parallel_machine_init(&argc, &argv);
    int myWorldRank = stk::parallel_machine_rank(m_commWorld);
    int numWorldRanks = stk::parallel_machine_size(m_commWorld);

    int defaultColor = stk::coupling::string_to_color(m_appName);
    int color = stk::get_command_line_option(argc, argv, "app-color", defaultColor);
    std::string defaultSyncMode = "Send";
    std::string syncModeString = stk::get_command_line_option<std::string>(argc, argv, "sync-mode", defaultSyncMode);
    m_syncMode = stk::coupling::string_to_sync_mode(syncModeString);

    m_commApp = stk::coupling::split_comm(m_commWorld, color);
    std::pair<int,int> rootRanks = stk::coupling::calc_my_root_and_other_root_ranks(m_commWorld, m_commApp);
    int myAppRank = stk::parallel_machine_rank(m_commApp);
    int numAppRanks = stk::parallel_machine_size(m_commApp);

    {
      std::ostringstream os;
      os << m_appName << ": color="<<color<<", world rank: " << myWorldRank<<" out of " << numWorldRanks
                      <<", app rank: " << myAppRank << " out of " << numAppRanks << std::endl;
      os << m_appName << ": my root-rank: " << rootRanks.first << ", other app's root-rank: " << rootRanks.second;
      std::cout << os.str() << std::endl;
    }
  }

  void communicate_initial_setup()
  {
    m_myInfo = create_sync_info();

    m_myInfo.set_value(stk::coupling::AppName, m_appName);
    m_myInfo.set_value(stk::coupling::InitialTime, 0.);
    double timeStep = calculate_time_step();
    m_myInfo.set_value(stk::coupling::TimeStep, timeStep);
    constexpr int numberOfSteps = 5;
    m_myInfo.set_value(stk::coupling::FinalTime, numberOfSteps * m_myInfo.get_value<double>(stk::coupling::TimeStep));
    m_myInfo.set_value(stk::coupling::TimeSyncMode, m_syncMode);

    m_otherInfo = m_myInfo.exchange(m_commWorld, m_commApp);

    stk::coupling::check_sync_mode_consistency(m_myInfo, m_otherInfo);

    {
      std::ostringstream os;
      os << m_appName << ": other app 'm_appName': " << m_otherInfo.get_value<std::string>(stk::coupling::AppName);
      if(stk::parallel_machine_rank(m_commApp) == 0) std::cout << os.str() << std::endl;
    }

    {
      ThrowRequireMsg(stk::coupling::check_consistency<double>(m_myInfo, m_otherInfo, stk::coupling::InitialTime, m_syncMode),
                       m_appName << ": initial time is inconsistent with " << m_otherInfo.get_value<std::string>(stk::coupling::AppName));
      ThrowRequireMsg(m_syncMode == stk::coupling::Send || m_otherInfo.has_value<double>(stk::coupling::TimeStep), 
                       m_appName << ": other app ("<< m_otherInfo.get_value<std::string>(stk::coupling::AppName)<<") doesn't have time step");
      ThrowRequireMsg(m_syncMode == stk::coupling::Send || m_otherInfo.has_value<double>(stk::coupling::FinalTime), 
                       m_appName << ": other app ("<< m_otherInfo.get_value<std::string>(stk::coupling::AppName)<<") doesn't have final time");

      m_currentTime = 0.0;
      m_finalTime = numberOfSteps*timeStep;
      m_iWasToldToStop = false;
    }

    m_myInfo.set_value(stk::coupling::CurrentTime, m_myInfo.get_value<double>(stk::coupling::InitialTime));
  }

  void setup_fields_and_transfers()
  {
  }

  void communicate_time_step_info()
  {
    double timeStep = calculate_time_step();

    m_myInfo.set_value(stk::coupling::TimeStep, timeStep);
    m_myInfo.set_value(stk::coupling::CurrentTime, m_currentTime);
    m_myInfo.set_value(m_doneFlagName, false);

    m_otherInfo = m_myInfo.exchange(m_commWorld, m_commApp);
    m_iWasToldToStop = m_otherInfo.get_value<bool>(m_doneFlagName, false);
  }

  void perform_transfers()
  {
  }

  void do_physics_solve()
  {
  }

  void update_current_time()
  {
    const double timeStep = calculate_time_step();
    m_currentTime += timeStep;
    m_myInfo.set_value(stk::coupling::CurrentTime, m_currentTime);
    m_currentTime = stk::coupling::choose_value(m_myInfo, m_otherInfo, stk::coupling::CurrentTime, m_syncMode);

    std::ostringstream os;
    os << m_appName << ": "<<stk::coupling::CurrentTime<<": " << m_currentTime << ", time step: " << timeStep;
    if (stk::parallel_machine_rank(m_commApp) == 0) std::cout << os.str() << std::endl;
  }

  bool time_to_stop()
  {
    bool timeToStop = m_iWasToldToStop || m_currentTime >= m_finalTime;
    return timeToStop;
  }

  void communicate_finish()
  {
    if (!m_iWasToldToStop) {
      m_myInfo.set_value(m_doneFlagName, true);
      m_myInfo.exchange(m_commWorld, m_commApp);
    }
  }

private:
  stk::coupling::SyncInfo create_sync_info()
  {
    return stk::coupling::SyncInfo(m_appName);
  }

  std::string m_appName;
  std::shared_ptr<mock::StkMesh> m_mesh;
  std::string m_doneFlagName;

  stk::ParallelMachine m_commWorld;
  stk::ParallelMachine m_commApp;

  stk::coupling::SyncInfo m_myInfo;
  stk::coupling::SyncInfo m_otherInfo;

  stk::coupling::SyncMode m_syncMode;
  double m_currentTime;
  double m_finalTime;
  bool m_iWasToldToStop;
  int m_step;
  bool m_doingSendTransfer;
  bool m_doingRecvTransfer;
  std::string m_sendFieldName;
  std::string m_recvFieldName;
};

int main(int argc, char** argv)
{
  MockAria app;
  app.read_input_and_setup_split_comms(argc, argv);
  app.communicate_initial_setup();
  app.setup_fields_and_transfers();

  do {
    app.communicate_time_step_info();
    if (app.time_to_stop()) break;
    app.perform_transfers();
    app.do_physics_solve();
    app.update_current_time();
  } while (!app.time_to_stop());

  app.communicate_finish();

  return 0;
}
