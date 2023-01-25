
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Constants.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_coupling/Version.hpp>
#include <stk_transfer/ReducedDependencyGeometricTransfer.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/Version.hpp>
#include <stk_util/parallel/CouplingVersions_impl.hpp>
#include <stk_util/parallel/CouplingVersions.hpp>
#include "MockUtils.hpp"
#include "SparcMesh.hpp"
#include "MockMeshUtils.hpp"
#include "SparcSendAdapter.hpp"
#include "EmptyRecvAdapter.hpp"
#include "SendInterpolate.hpp"
#include <iostream>
#include <sstream>

class MockSparc
{
  using SendTransfer = stk::transfer::ReducedDependencyGeometricTransfer<
                                 mock::SendInterpolate<mock::SparcSendAdapter, mock::EmptyRecvAdapter>>;

public:
  MockSparc()
    : m_appName("Mock-Sparc"),
      m_mesh(),
      m_doneFlagName("time step status"),
      m_splitComms(),
      m_otherColor(),
      m_iAmRootRank(false),
      m_myInfo(),
      m_otherInfo(),
      m_timeSyncMode(),
      m_initialTime(),
      m_timeStep(),
      m_finalTime(),
      m_currentTime(),
      m_isTimeToStop(),
      m_doingSendTransfer(false),
      m_sendFieldName()
  { }

  ~MockSparc()
  {
    stk::parallel_machine_finalize();
  }

  void read_input_and_setup_split_comms(int argc, char** argv)
  {
    MPI_Comm commWorld = stk::parallel_machine_init(&argc, &argv);

    int defaultColor = stk::coupling::string_to_color(m_appName);
    int color = stk::get_command_line_option(argc, argv, "app-color", defaultColor);
    int coupling_version_override = stk::get_command_line_option(argc, argv, "stk_coupling_version", STK_MAX_COUPLING_VERSION);
    const std::string defaultFileName = "generated:1x1x4|sideset:x";
    std::string meshFileName = stk::get_command_line_option(argc, argv, "mesh", defaultFileName);

    stk::util::impl::set_coupling_version(coupling_version_override);
    stk::util::impl::set_error_on_reset(false);
    m_splitComms = stk::coupling::SplitComms(commWorld, color);
    MPI_Comm splitComm = m_splitComms.get_split_comm();
    int myAppRank = stk::parallel_machine_rank(splitComm);
    m_iAmRootRank = myAppRank == 0;

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

    {
      stk::coupling::PairwiseRanks rootRanks = m_splitComms.get_pairwise_root_ranks(otherColors[0]);
      int numAppRanks = stk::parallel_machine_size(splitComm);
      int myWorldRank = stk::parallel_machine_rank(commWorld);
      int numWorldRanks = stk::parallel_machine_size(commWorld);
 
      if (m_iAmRootRank) {
        std::ostringstream os;
        os << m_appName << ": STK version: " << stk::version_string() 
           << " (Coupling Version: " << stk::util::get_common_coupling_version() << ")"<<std::endl;
        os << m_appName << ": my world rank is: " << myWorldRank << " out of " << numWorldRanks
           <<", app rank: " << myAppRank << " out of " << numAppRanks << std::endl;
        os << m_appName << ": my root-rank: " << rootRanks.localColorRoot << ", other app's root-rank: " << rootRanks.otherColorRoot;
        std::cout << os.str() << std::endl;
      }
    }

    std::vector<std::string> fieldNames = {"sparc-traction", "heat-transfer-coefficient"};
    mock_utils::read_mesh(splitComm, meshFileName, fieldNames, m_mesh);

    // TODO: put timeSyncMode in a command-line arg like mock-aria
    m_timeSyncMode = stk::coupling::Send; // This is usually Send, but could be Receive for SPARC-SPARC MPMD coupling
    constexpr int numberOfSteps = 5;
    m_initialTime = 0.0;
    m_timeStep = 0.1;
    m_finalTime = numberOfSteps * m_timeStep;
    m_currentTime = m_initialTime;
    m_isTimeToStop = false;
  }

  stk::coupling::SyncInfo perform_exchange(const stk::coupling::SyncInfo& info)
  {
    return info.exchange(m_splitComms, m_otherColor);
  }

  void communicate_and_check_initial_setup_compatibility()
  {
    m_myInfo = create_sync_info();

    m_myInfo.set_value(stk::coupling::AppName, m_appName); // SPARC currently does not exchange this, but this will be helpful
    m_myInfo.set_value(stk::coupling::TimeSyncMode, m_timeSyncMode);
    // Assuming that these times and time steps are associated with coupling (i.e. TimeStep is the coupled time step)
    m_myInfo.set_value(stk::coupling::InitialTime, m_initialTime);
    m_myInfo.set_value(stk::coupling::TimeStep, m_timeStep);
    m_myInfo.set_value(stk::coupling::FinalTime, m_finalTime);

    m_otherInfo = perform_exchange(m_myInfo);

    {
      std::ostringstream os;
      os << m_appName << ": other app 'app_name': " << m_otherInfo.get_value<std::string>(stk::coupling::AppName);
      if (m_iAmRootRank) std::cout << os.str() << std::endl;
    }

    std::string otherAppName = m_otherInfo.get_value<std::string>(stk::coupling::AppName, "none");

    if (otherAppName=="Mock-Salinas") {
      m_doingSendTransfer = true;
      m_sendFieldName = "sparc-traction";
    }
    if (otherAppName=="Mock-Aria") {
      m_doingSendTransfer = true;
      m_sendFieldName = "heat-transfer-coefficient";
    }

    {
      if (m_doingSendTransfer) {
        std::ostringstream os;
        os << m_appName << ": will do send-transfer (field='"<<m_sendFieldName<<"') "
           <<" to other app: "<<otherAppName<<std::endl;
        if (m_iAmRootRank) std::cout << os.str() << std::endl;
      }
    }

    stk::coupling::check_sync_mode_consistency(m_myInfo, m_otherInfo);

    {
      ThrowRequireMsg(stk::coupling::check_consistency<double>(m_myInfo, m_otherInfo, stk::coupling::InitialTime, m_timeSyncMode),
          m_appName << ": initial time is inconsistent with " << m_otherInfo.get_value<std::string>(stk::coupling::AppName));
      ThrowRequireMsg(m_timeSyncMode == stk::coupling::Send || m_otherInfo.has_value<double>(stk::coupling::TimeStep), 
          m_appName << ": other app ("<< m_otherInfo.get_value<std::string>(stk::coupling::AppName)<<") doesn't have time step");
      ThrowRequireMsg(m_timeSyncMode == stk::coupling::Send || m_otherInfo.has_value<double>(stk::coupling::FinalTime), 
          m_appName << ": other app ("<< m_otherInfo.get_value<std::string>(stk::coupling::AppName)<<") doesn't have final time");
    }
  }

  void check_field_sizes(std::vector<std::pair<std::string,int>> sendFields,
                         std::vector<std::pair<std::string,int>> recvFields)
  {
    ThrowRequireMsg(sendFields.size() == recvFields.size(), "Number of send-fields ("
       <<sendFields.size()<<") doesn't match number of recv-fields ("<<recvFields.size()
       <<")");
    for (unsigned i=0; i<sendFields.size(); ++i) {
      ThrowRequireMsg(sendFields[i].second == recvFields[i].second,
        "Send-field size ("<<sendFields[i].first<<","<<sendFields[i].second<<") "
        <<"doesn't match Recv-field size ("<<recvFields[i].first<<","<<recvFields[i].second<<")");
    }
  }

  void setup_fields_and_transfers()
  {
    if (!m_doingSendTransfer) { return; }

    std::vector<std::pair<std::string,int>> mySendFields;
    std::vector<std::pair<std::string,int>> myRecvFields;
    mySendFields.push_back(std::make_pair(m_sendFieldName, m_mesh->get_field_size()));

    stk::coupling::SyncInfo info = create_sync_info();
    info.set_value("SendFields", mySendFields);
    info.set_value("RecvFields", myRecvFields);

    stk::coupling::SyncInfo otherInfo = perform_exchange(info);
    std::vector<std::pair<std::string,int>> otherSendFields = otherInfo.get_value<std::vector<std::pair<std::string,int>>>("SendFields");
    std::vector<std::pair<std::string,int>> otherRecvFields = otherInfo.get_value<std::vector<std::pair<std::string,int>>>("RecvFields");

    check_field_sizes(mySendFields, otherRecvFields);
    check_field_sizes(otherSendFields, myRecvFields);

    m_mesh->set_sparc_field_values(m_sendFieldName, 4.4);
    std::shared_ptr<mock::SparcSendAdapter> sendAdapter =
       std::make_shared<mock::SparcSendAdapter>(m_splitComms.get_parent_comm(), *m_mesh, m_sendFieldName);
    std::shared_ptr<mock::EmptyRecvAdapter> recvAdapter;
    m_sendTransfer.reset(new SendTransfer(sendAdapter, recvAdapter, "MockSparcSendTransfer", m_splitComms.get_parent_comm()));

    m_sendTransfer->coarse_search();
    m_sendTransfer->communication();
    m_sendTransfer->local_search();

    // todo Create Fields
    // todo Construct transfers
    // todo Loop transfers and initialize (In SPARC, we loop time solvers, and then loop transfers that are driven by each time solver)
    // Each transfer (t) initialize looks like the following:
    //   t->coarse_search();
    //   t->communication();
    //   t->local_search();
    //   CheckFieldSize() - we would like STK to take care of this - this should check for field size consistency
    //   and ideally some way of checking that field transfers are in the same order between apps.
    //   Also, possibly printing this info including field names
    //   This will allow the user to at least see in the output what fields are being sent and received
    //
    // Loop transfers and apply if input file requests pre-solve apply (used for initial condition transfer)
  }

  void communicate_time_step_info()
  {
    m_myInfo = create_sync_info();

    m_myInfo.set_value(m_doneFlagName, m_isTimeToStop);
    m_myInfo.set_value(stk::coupling::TimeStep, m_timeStep);
    m_myInfo.set_value(stk::coupling::CurrentTime, m_currentTime);
    m_myInfo.set_value(stk::coupling::FinalTime, m_finalTime);

    m_otherInfo = perform_exchange(m_myInfo);
  }

  bool time_to_stop()
  {
    return m_isTimeToStop || m_otherInfo.get_value<bool>(m_doneFlagName, false);
  }

  void agree_with_other_app_on_timestep()
  {
    m_timeStep = stk::coupling::choose_value(m_myInfo, m_otherInfo, stk::coupling::TimeStep, m_timeSyncMode);
    m_currentTime = stk::coupling::choose_value(m_myInfo, m_otherInfo, stk::coupling::CurrentTime, m_timeSyncMode);
  }

  void physics_inner_subcycling_loop()
  {
    // Solve physics inner subcycling loop (i.e. from currentTime to currentTime + timeStep
    double sparcTime = m_currentTime;
    const double nextSyncTime = sparcTime + m_timeStep;
    do {
      double sparcTimeStep = 1.0e-3; // arbitrarily set in this example
      if (sparcTime < nextSyncTime && sparcTime + sparcTimeStep > nextSyncTime) {
        sparcTimeStep = nextSyncTime - sparcTime;
      }

      // Solve physics over one time step (do nothing here in mock app)
      // this is where we may add logic to do transfers at the nonlinear iteration level

      sparcTime += sparcTimeStep;

      // isTimeToStop could be set to true if solution convergence achieved, or NaN detected, or other possible reasons
      m_isTimeToStop = false;

    } while (sparcTime < nextSyncTime);
  }

  void perform_transfers()
  {
    if (m_doingSendTransfer) {
      m_sendTransfer->apply();
    }
    //   for(auto t : transfers)
    //     t->apply();
  }

  void compute_my_timestep_and_decide_if_i_want_to_stop()
  {
    m_currentTime += m_timeStep;

    if (m_currentTime >= m_finalTime) {
      m_isTimeToStop = true;
    }
    else if (m_currentTime + m_timeStep > m_finalTime)
    {
      m_timeStep = m_finalTime - m_currentTime;
    }

    {
      if (m_iAmRootRank) std::cout << m_appName << ": "<<stk::coupling::CurrentTime<<": " << m_currentTime << ", final time: " << m_finalTime << ", isTimeToStop: " << m_isTimeToStop << std::endl;
    }
  }

  void communicate_finish()
  {
    m_myInfo = create_sync_info();
    m_myInfo.set_value(m_doneFlagName, true);
    m_otherInfo = perform_exchange(m_myInfo);
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

  const std::string m_appName;
  std::shared_ptr<mock::SparcMesh> m_mesh;
  const std::string m_doneFlagName;

  stk::coupling::SplitComms m_splitComms;
  int m_otherColor;
  bool m_iAmRootRank;

  stk::coupling::SyncInfo m_myInfo;
  stk::coupling::SyncInfo m_otherInfo;

  stk::coupling::SyncMode m_timeSyncMode;
  double m_initialTime;
  double m_timeStep;
  double m_finalTime;
  double m_currentTime;
  bool m_isTimeToStop;

  std::shared_ptr<SendTransfer> m_sendTransfer;
  bool m_doingSendTransfer;
  std::string m_sendFieldName;
};

int main(int argc, char** argv)
{
  MockSparc app;
  app.read_input_and_setup_split_comms(argc, argv);
  if (app.get_number_of_other_coupled_apps() != 1) {
    return 0;
  }

  app.communicate_and_check_initial_setup_compatibility();
  app.setup_fields_and_transfers();

  do {
    app.communicate_time_step_info();
    if (app.time_to_stop()) break;

    app.agree_with_other_app_on_timestep();
    app.physics_inner_subcycling_loop();
    app.perform_transfers();

    app.compute_my_timestep_and_decide_if_i_want_to_stop();
  } while (!app.time_to_stop());

  app.communicate_finish();

  return 0;
}
