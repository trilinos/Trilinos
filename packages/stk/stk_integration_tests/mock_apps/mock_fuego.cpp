
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

class MockFuego
{
public:
  MockFuego()
  : m_appName("Mock-Fuego"),
    m_mesh(),
    m_doneFlagName("time step status"),
    m_iAmRootRank(false),
    m_doingSendTransfer(false),
    m_sendFieldName()
  {}

  ~MockFuego()
  {
    stk::parallel_machine_finalize();
  }

  void read_input_and_setup_split_comms(int argc, char** argv)
  {
    m_commWorld = stk::parallel_machine_init(&argc, &argv);
    int myWorldRank = stk::parallel_machine_rank(m_commWorld);
    int numWorldRanks = stk::parallel_machine_size(m_commWorld);

    {
      std::ostringstream os;
      os << m_appName << ": my world rank is: " << myWorldRank << " out of " << numWorldRanks;
      std::cout << os.str() << std::endl;
    }
    int defaultColor = stk::coupling::string_to_color(m_appName);
    int color = stk::get_command_line_option(argc, argv, "app-color", defaultColor);

    m_commApp = stk::coupling::split_comm(m_commWorld, color);
    std::pair<int,int> rootRanks = stk::coupling::calc_my_root_and_other_root_ranks(m_commWorld, m_commApp);
    int myAppRank = stk::parallel_machine_rank(m_commApp);
    m_iAmRootRank = myAppRank == 0;
    int numAppRanks = stk::parallel_machine_size(m_commApp);

    {
      std::ostringstream os;
      os << m_appName << ": color="<<color<<", my app rank is: " << myAppRank << " out of " << numAppRanks << std::endl;
      os << m_appName << ": my root-rank: " << rootRanks.first << ", other app's root-rank: " << rootRanks.second;
      std::cout << os.str() << std::endl;
    }
  }

  void communicate_initial_setup()
  {
    m_myInfo = create_sync_info();

    m_myInfo.set_value(stk::coupling::AppName, m_appName);

    m_otherInfo = m_myInfo.exchange(m_commWorld, m_commApp);

    {
      std::ostringstream os;
      os << m_appName << ": other app 'app_name': " << m_otherInfo.get_value<std::string>(stk::coupling::AppName);
      if(m_iAmRootRank) std::cout << os.str() << std::endl;
    }
  }

  double compute_time_step()
  {
    return 0.005;
  }

  void do_physics_solve()
  {
  }

  void setup_fields_and_transfers()
  {
    constexpr int numberOfSteps = 10;

    double initialTime = 0.0;
    double timeStep = compute_time_step();
    m_step = 0;
    m_finalTime = initialTime + numberOfSteps * timeStep;
  }

  void perform_transfers()
  {
  }

  bool time_to_stop()
  {
    bool timeToStop = (m_currentTime >= m_finalTime) || m_otherInfo.get_value<bool>(m_doneFlagName, false);
  
    if (timeToStop) {
      std::ostringstream os;
      os << m_appName << " finished, final time: " << m_finalTime << std::endl;
      if(m_iAmRootRank) std::cout << os.str();
    }
    return timeToStop;
  }

  void communicate_time_step_info()
  {
    m_myInfo.set_value("step", m_step);
    m_myInfo.set_value(stk::coupling::CurrentTime, m_currentTime);
    m_myInfo.set_value(stk::coupling::FinalTime, m_finalTime);

    m_otherInfo = m_myInfo.exchange(m_commWorld, m_commApp);
  }

  void update_current_time()
  {
    m_currentTime = m_step*compute_time_step();

    std::ostringstream os;
    os << m_appName << ": "<<stk::coupling::CurrentTime<<": " << m_currentTime;
    if (m_iAmRootRank) std::cout << os.str() << std::endl;

    ++m_step;
  }

  void communicate_finish()
  {
    m_myInfo = create_sync_info();
    m_myInfo.set_value(m_doneFlagName, true);
    m_myInfo.exchange(m_commWorld, m_commApp);
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
  bool m_iAmRootRank;

  stk::coupling::SyncInfo m_myInfo;
  stk::coupling::SyncInfo m_otherInfo;

  double m_currentTime;
  double m_finalTime;
  int m_step;
  bool m_doingSendTransfer;
  std::string m_sendFieldName;
};

int main(int argc, char** argv)
{
  MockFuego app;
  app.read_input_and_setup_split_comms(argc, argv);
  app.communicate_initial_setup();
  app.setup_fields_and_transfers();

  do {
    app.communicate_time_step_info();
    app.perform_transfers();
    app.do_physics_solve();
    app.update_current_time();
  } while (!app.time_to_stop());

  app.communicate_finish();

  return 0;
}
