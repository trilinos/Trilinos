
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Constants.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_transfer/ReducedDependencyGeometricTransfer.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/Version.hpp>
#include <stk_util/parallel/CouplingVersions_impl.hpp>
#include <stk_util/parallel/CouplingVersions.hpp>
#include "MockUtils.hpp"
#include "StkMesh.hpp"
#include "MockMeshUtils.hpp"
#include "StkSendAdapter.hpp"
#include "StkRecvAdapter.hpp"
#include "SendInterpolate.hpp"
#include "RecvInterpolate.hpp"
#include <iostream>
#include <sstream>

class MockAria
{
  using SendTransfer = stk::transfer::ReducedDependencyGeometricTransfer<
                                 mock::SendInterpolate<mock::StkSendAdapter, mock::StkRecvAdapter>>;
  using RecvTransfer = stk::transfer::ReducedDependencyGeometricTransfer<
                                 mock::RecvInterpolate<mock::StkSendAdapter, mock::StkRecvAdapter>>;
public:
  MockAria()
  : m_appName("Mock-Aria"),
    m_mesh(),
    m_doneFlagName("time step status"),
    m_splitComms(),
    m_otherColor(),
    m_iAmRootRank(false),
    m_myInfo(),
    m_otherInfo(),
    m_syncMode(),
    m_currentTime(),
    m_finalTime(),
    m_iWasToldToStop(),
    m_iWantToStop(false),
    m_sendTransfer(),
    m_recvTransfer1(),
    m_recvTransfer2(),
    m_doingSendTransfer(false),
    m_doingRecvTransfer(false),
    m_wrongTransferOrder(false),
    m_sendFieldName(),
    m_recvFieldName1(),
    m_recvFieldName2()
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
    MPI_Comm commWorld = stk::parallel_machine_init(&argc, &argv);
    int myWorldRank = stk::parallel_machine_rank(commWorld);
    int numWorldRanks = stk::parallel_machine_size(commWorld);

    int defaultColor = stk::coupling::string_to_color(m_appName);
    int color = stk::get_command_line_option(argc, argv, "app-color", defaultColor);

    int coupling_version_override = stk::get_command_line_option(argc, argv, "stk_coupling_version", STK_MAX_COUPLING_VERSION);

    const std::string defaultFileName = "generated:1x1x4|sideset:x";
    std::string meshFileName = stk::get_command_line_option(argc, argv, "mesh", defaultFileName);

    const std::string defaultPartName = "surface_1";
    std::string partName = stk::get_command_line_option(argc, argv, "part-name", defaultPartName);

    m_wrongTransferOrder = stk::get_command_line_option(argc, argv, "wrong-transfer-order", false);

    stk::util::impl::set_error_on_reset(false);
    std::string defaultSyncMode = "Send";
    std::string syncModeString = stk::get_command_line_option<std::string>(argc, argv, "sync-mode", defaultSyncMode);
    m_syncMode = stk::coupling::string_to_sync_mode(syncModeString);

    m_splitComms = stk::coupling::SplitComms(commWorld, color);
    m_splitComms.set_free_comms_in_destructor(true);
    stk::util::impl::set_coupling_version(commWorld, coupling_version_override);
    MPI_Comm splitComm = m_splitComms.get_split_comm();
    int myAppRank = stk::parallel_machine_rank(splitComm);
    int numAppRanks = stk::parallel_machine_size(splitComm);
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
    stk::coupling::PairwiseRanks rootRanks = m_splitComms.get_pairwise_root_ranks(m_otherColor);

    if (m_iAmRootRank) {
      std::ostringstream os;
      os << m_appName << ": STK version: " << stk::version_string()
         << " (Coupling Version: " << stk::util::get_common_coupling_version() << ") ("<<coupling_version_override<<")" << std::endl;
      os << m_appName << ", color="<<color<<", world rank: " << myWorldRank<<" out of " << numWorldRanks
                      <<", app rank: " << myAppRank << " out of " << numAppRanks << std::endl;
      os << m_appName << ": my root-rank: " << rootRanks.localColorRoot << ", other app's root-rank: " << rootRanks.otherColorRoot;
      std::cout << os.str() << std::endl;
    }

    std::vector<std::string> fieldNames = {"reference-temperature","heat-transfer-coefficient1", "heat-transfer-coefficient2"};
    mock_utils::read_mesh(splitComm, meshFileName, partName, fieldNames, m_mesh);
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

    m_otherInfo = m_myInfo.exchange(m_splitComms, m_otherColor);

    stk::coupling::check_sync_mode_consistency(m_myInfo, m_otherInfo);

    {
      std::ostringstream os;
      os << m_appName << ": other app 'm_appName': " << m_otherInfo.get_value<std::string>(stk::coupling::AppName);
      if(m_iAmRootRank) std::cout << os.str() << std::endl;
    }

    {
      STK_ThrowRequireMsg(stk::coupling::check_consistency<double>(m_myInfo, m_otherInfo, stk::coupling::InitialTime, m_syncMode),
                       m_appName << ": initial time is inconsistent with " << m_otherInfo.get_value<std::string>(stk::coupling::AppName));
      STK_ThrowRequireMsg(m_syncMode == stk::coupling::Send || m_otherInfo.has_value<double>(stk::coupling::TimeStep), 
                       m_appName << ": other app ("<< m_otherInfo.get_value<std::string>(stk::coupling::AppName)<<") doesn't have time step");
      STK_ThrowRequireMsg(m_syncMode == stk::coupling::Send || m_otherInfo.has_value<double>(stk::coupling::FinalTime), 
                       m_appName << ": other app ("<< m_otherInfo.get_value<std::string>(stk::coupling::AppName)<<") doesn't have final time");

      m_currentTime = 0.0;
      m_finalTime = numberOfSteps*timeStep;
      m_iWasToldToStop = false;
    }

    std::string otherAppName = m_otherInfo.get_value<std::string>(stk::coupling::AppName, "none");
    if (otherAppName == "Mock-Salinas") {
      m_doingSendTransfer = true;
      m_sendFieldName = "reference-temperature";
    }
    if (otherAppName == "Mock-Sparc") {
      m_doingRecvTransfer = true;
      m_recvFieldName1 = "heat-transfer-coefficient1";
      m_recvFieldName2 = "heat-transfer-coefficient2";
    }

    {
      std::ostringstream os;
      if (m_doingSendTransfer) {
        os << m_appName << ": will send-transfer (field='"<<m_sendFieldName<<"') "
           << " to other app: "<<otherAppName<<std::endl;
      }

      if (m_doingRecvTransfer) {
        os << m_appName << ": will recv-transfer (field='"<<m_recvFieldName1<<","<<m_recvFieldName2<<"') "
           <<" from other app: "<<otherAppName<<std::endl;
      }

      if (m_iAmRootRank) std::cout << os.str() << std::endl;
    }

    m_myInfo.set_value(stk::coupling::CurrentTime, m_myInfo.get_value<double>(stk::coupling::InitialTime));
  }

  void check_field_sizes(std::vector<std::pair<std::string,int>> sendFields,
                         std::vector<std::pair<std::string,int>> recvFields)
  {
    STK_ThrowRequireMsg(sendFields.size() == recvFields.size(), "Number of send-fields ("
       <<sendFields.size()<<") doesn't match number of recv-fields ("<<recvFields.size()
       <<")");
    for (unsigned i=0; i<sendFields.size(); ++i) {
      STK_ThrowRequireMsg(sendFields[i].second == recvFields[i].second,
        "Send-field size ("<<sendFields[i].first<<","<<sendFields[i].second<<") "
        <<"doesn't match Recv-field size ("<<recvFields[i].first<<","<<recvFields[i].second<<")");
    }   
  }

  void setup_fields_and_transfers()
  {
    if (!m_doingSendTransfer && !m_doingRecvTransfer) { return; }

    std::vector<std::pair<std::string,int>> mySendFields;
    std::vector<std::pair<std::string,int>> myRecvFields;
    if (m_doingSendTransfer) {
      mySendFields.push_back(std::make_pair(m_sendFieldName, m_mesh->get_field_size()));
    }
    if (m_doingRecvTransfer) {
      myRecvFields.push_back(std::make_pair(m_recvFieldName1, m_mesh->get_field_size()));
      myRecvFields.push_back(std::make_pair(m_recvFieldName2, m_mesh->get_field_size()));
    }

    stk::coupling::SyncInfo info = create_sync_info();
    info.set_value("SendFields", mySendFields);
    info.set_value("RecvFields", myRecvFields);

    MPI_Comm parentComm = m_splitComms.get_parent_comm();
    stk::coupling::SyncInfo otherInfo = info.exchange(m_splitComms, m_otherColor);
    std::vector<std::pair<std::string,int>> otherSendFields = otherInfo.get_value<std::vector<std::pair<std::string,int>>>("SendFields");
    std::vector<std::pair<std::string,int>> otherRecvFields = otherInfo.get_value<std::vector<std::pair<std::string,int>>>("RecvFields");

    check_field_sizes(mySendFields, otherRecvFields);
    check_field_sizes(otherSendFields, myRecvFields);

    if (m_doingSendTransfer) {
      m_mesh->set_stk_field_values(m_sendFieldName, 9.9);
      std::shared_ptr<mock::StkSendAdapter> sendAdapter =
         std::make_shared<mock::StkSendAdapter>(parentComm, *m_mesh, m_sendFieldName);
      std::shared_ptr<mock::StkRecvAdapter> nullRecvAdapter;
      m_sendTransfer.reset(new SendTransfer(sendAdapter, nullRecvAdapter, "MockAriaSendTransfer", parentComm));

      m_sendTransfer->coarse_search();
      m_sendTransfer->communication();
      m_sendTransfer->local_search();
    }

    if (m_doingRecvTransfer) {
      std::shared_ptr<mock::StkSendAdapter> nullSendAdapter;
      std::shared_ptr<mock::StkRecvAdapter> recvAdapter1 =
         std::make_shared<mock::StkRecvAdapter>(parentComm, *m_mesh, m_recvFieldName1);
      std::shared_ptr<mock::StkRecvAdapter> recvAdapter2 =
         std::make_shared<mock::StkRecvAdapter>(parentComm, *m_mesh, m_recvFieldName2);
      m_recvTransfer1.reset(new RecvTransfer(nullSendAdapter, recvAdapter1, "MockAriaRecvTransfer1", parentComm));
      m_recvTransfer2.reset(new RecvTransfer(nullSendAdapter, recvAdapter2, "MockAriaRecvTransfer2", parentComm));

      m_recvTransfer1->coarse_search();
      m_recvTransfer1->communication();
      m_recvTransfer1->local_search();

      m_recvTransfer2->coarse_search();
      m_recvTransfer2->communication();
      m_recvTransfer2->local_search();
    }
  }

  void communicate_time_step_info()
  {
    double timeStep = calculate_time_step();

    m_myInfo.set_value(stk::coupling::TimeStep, timeStep);
    m_myInfo.set_value(stk::coupling::CurrentTime, m_currentTime);
    m_myInfo.set_value(m_doneFlagName, false);

    m_otherInfo = m_myInfo.exchange(m_splitComms, m_otherColor);
    m_iWasToldToStop = m_otherInfo.get_value<bool>(m_doneFlagName, false);
  }

  void perform_transfers()
  {
    if (m_doingSendTransfer) {
      m_sendTransfer->apply();
    }
    if (m_doingRecvTransfer) {
      m_mesh->set_stk_field_values(m_recvFieldName1, 0.0);
      m_mesh->set_stk_field_values(m_recvFieldName2, 0.0);
      m_recvTransfer1->meshb()->called_update_values = false;
      m_recvTransfer2->meshb()->called_update_values = false;

      m_recvTransfer1->apply();
      m_recvTransfer2->apply();

      STK_ThrowRequire(m_recvTransfer1->meshb()->called_update_values);
      STK_ThrowRequire(m_recvTransfer2->meshb()->called_update_values);
      m_recvTransfer1->meshb()->called_update_values = false;
      m_recvTransfer2->meshb()->called_update_values = false;

      const double expectedField1Value = 4.4;
      const bool values1Match = m_mesh->verify_stk_field_values(m_recvFieldName1, expectedField1Value);
      STK_ThrowRequireMsg(values1Match, "Mock-Aria error, field1-values are not correct after transfer");

      const double expectedField2Value = 8.8;
      const bool values2Match = m_mesh->verify_stk_field_values(m_recvFieldName2, expectedField2Value);
      STK_ThrowRequireMsg(values2Match, "Mock-Aria error, field2-values are not correct after transfer");
    }
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

    {
      std::ostringstream os;
      os << m_appName << ": "<<stk::coupling::CurrentTime<<": " << m_currentTime << ", time step: " << timeStep;
      if (m_iAmRootRank) std::cout << os.str() << std::endl;
    }
  }

  bool time_to_stop()
  {
    bool timeToStop = m_iWasToldToStop || m_currentTime >= m_finalTime || m_iWantToStop;
    return timeToStop;
  }

  void communicate_finish()
  {
    if (!m_iWasToldToStop) {
      m_myInfo.set_value(m_doneFlagName, true);
      m_myInfo.exchange(m_splitComms, m_otherColor);
    }
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

  stk::coupling::SyncMode m_syncMode;
  double m_currentTime;
  double m_finalTime;
  bool m_iWasToldToStop;
  bool m_iWantToStop;
  std::shared_ptr<SendTransfer> m_sendTransfer;
  std::shared_ptr<RecvTransfer> m_recvTransfer1;
  std::shared_ptr<RecvTransfer> m_recvTransfer2;
  bool m_doingSendTransfer;
  bool m_doingRecvTransfer;
  bool m_wrongTransferOrder;
  std::string m_sendFieldName;
  std::string m_recvFieldName1;
  std::string m_recvFieldName2;
};

int main(int argc, char** argv)
{
  MockAria app;
  app.read_input_and_setup_split_comms(argc, argv);
  if (app.get_number_of_other_coupled_apps() != 1) {
    return 0;
  }

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
