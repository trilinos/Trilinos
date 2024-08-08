
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Constants.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/CouplingVersions.hpp>
#include <stk_util/parallel/CouplingVersions_impl.hpp>
#include <stk_util/Version.hpp>
#include "MockUtils.hpp"
#include "StkMesh.hpp"
#include "MockMeshUtils.hpp"
#include <iostream>
#include <sstream>

class MockFuego
{
public:
  MockFuego()
  : m_appName("Mock-Fuego"),
    m_mesh(),
    m_doneFlagName("time step status"),
    m_splitComms(),
    m_otherColor(),
    m_iAmRootRank(false),
    m_myInfo(),
    m_otherInfo(),
    m_currentTime(),
    m_finalTime(),
    m_step(),
    m_sendFieldName()
  {}

  ~MockFuego()
  {
    stk::parallel_machine_finalize();
  }

  void read_input_and_setup_split_comms(int argc, char** argv)
  {
    MPI_Comm commWorld = stk::parallel_machine_init(&argc, &argv);
    int myWorldRank = stk::parallel_machine_rank(commWorld);
    int numWorldRanks = stk::parallel_machine_size(commWorld);

    {
      std::ostringstream os;
      os << m_appName << ": STK version: " << stk::version_string() 
         << " (Coupling Version: " << stk::util::get_common_coupling_version() << ")"
         << ": my world rank is: " << myWorldRank << " out of " << numWorldRanks;
      std::cout << os.str() << std::endl;
    }
    int defaultColor = stk::coupling::string_to_color(m_appName);
    int color = stk::get_command_line_option(argc, argv, "app-color", defaultColor);
    int coupling_version_override = stk::get_command_line_option(argc, argv, "stk_coupling_version", STK_MAX_COUPLING_VERSION);
    const std::string defaultFileName = "generated:1x1x4|sideset:x";
    std::string meshFileName = stk::get_command_line_option(argc, argv, "mesh", defaultFileName);

    stk::util::impl::set_error_on_reset(false);

    m_splitComms = stk::coupling::SplitComms(commWorld, color);
    m_splitComms.set_free_comms_in_destructor(true);
    stk::util::impl::set_coupling_version(commWorld, coupling_version_override);
    MPI_Comm splitComm = m_splitComms.get_split_comm();
    const std::vector<int>& otherColors = m_splitComms.get_other_colors();
    if (otherColors.size() != 1) {
      if (otherColors.empty()) {
        std::cout << m_appName << " No other colors, not running in MPMD mode." <<std::endl;
      }
      else {
        mock_utils::exchange_and_print_info(m_splitComms, m_appName, color);
      }
      return;
    }
    m_otherColor = otherColors[0];
    stk::coupling::PairwiseRanks rootRanks = m_splitComms.get_pairwise_root_ranks(m_otherColor);
    int myAppRank = stk::parallel_machine_rank(splitComm);
    m_iAmRootRank = myAppRank == 0;
    int numAppRanks = stk::parallel_machine_size(splitComm);

    {
      std::ostringstream os;
      os << m_appName << ": STK Version: " << stk::version_string()
         << " (Coupling Version: " << stk::util::get_common_coupling_version() << ")" << std::endl;
      os << m_appName << ": color="<<color<<", my app rank is: " << myAppRank << " out of " << numAppRanks << std::endl;
      os << m_appName << ": my root-rank: " << rootRanks.localColorRoot << ", other app's root-rank: " << rootRanks.otherColorRoot;
      std::cout << os.str() << std::endl;
    }

    std::vector<std::string> fieldNames = {"not-set-yet1", "not-set-yet2"};
    mock_utils::read_mesh(splitComm, meshFileName, "surface_1", fieldNames, m_mesh);
  }

  void communicate_initial_setup()
  {
    m_myInfo = create_sync_info();

    m_myInfo.set_value(stk::coupling::AppName, m_appName);

    m_otherInfo = m_myInfo.exchange(m_splitComms, m_otherColor);

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

    m_otherInfo = m_myInfo.exchange(m_splitComms, m_otherColor);
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
    m_myInfo.exchange(m_splitComms, m_otherColor);
  }

  unsigned get_number_of_other_coupled_apps() const
  {
    return m_splitComms.get_other_colors().size();
  }

private:
  stk::coupling::SyncInfo create_sync_info()
  {
    return stk::coupling::SyncInfo(m_appName);
  }

  std::string m_appName;
  std::shared_ptr<mock::StkMesh> m_mesh;
  std::string m_doneFlagName;

  stk::coupling::SplitComms m_splitComms;
  int m_otherColor;
  bool m_iAmRootRank;

  stk::coupling::SyncInfo m_myInfo;
  stk::coupling::SyncInfo m_otherInfo;

  double m_currentTime;
  double m_finalTime;
  int m_step;
  std::string m_sendFieldName;
};

int main(int argc, char** argv)
{
  MockFuego app;
  app.read_input_and_setup_split_comms(argc, argv);
  if (app.get_number_of_other_coupled_apps() != 1) {
    return 0;
  }

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
