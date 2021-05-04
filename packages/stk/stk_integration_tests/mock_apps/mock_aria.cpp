
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
#include "MockUtils.hpp"
#include "StkMesh.hpp"
#include "StkSendAdapter.hpp"
#include "StkRecvAdapter.hpp"
#include "EmptySendAdapter.hpp"
#include "EmptyRecvAdapter.hpp"
#include "SendInterpolate.hpp"
#include "RecvInterpolate.hpp"
#include <iostream>
#include <sstream>

class MockAria
{
  using SendTransfer = stk::transfer::ReducedDependencyGeometricTransfer<
                                 mock::SendInterpolate<mock::StkSendAdapter, mock::EmptyRecvAdapter>>;
  using RecvTransfer = stk::transfer::ReducedDependencyGeometricTransfer<
                                 mock::RecvInterpolate<mock::EmptySendAdapter, mock::StkRecvAdapter>>;
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
    m_step(),
    m_sendTransfer(),
    m_recvTransfer(),
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
    MPI_Comm commWorld = stk::parallel_machine_init(&argc, &argv);
    int myWorldRank = stk::parallel_machine_rank(commWorld);
    int numWorldRanks = stk::parallel_machine_size(commWorld);

    int defaultColor = stk::coupling::string_to_color(m_appName);
    int color = stk::get_command_line_option(argc, argv, "app-color", defaultColor);
    std::string defaultSyncMode = "Send";
    std::string syncModeString = stk::get_command_line_option<std::string>(argc, argv, "sync-mode", defaultSyncMode);
    m_syncMode = stk::coupling::string_to_sync_mode(syncModeString);

    m_splitComms = stk::coupling::SplitComms(commWorld, color);
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

    {
      std::ostringstream os;
      os << m_appName << ": STK version: " << stk::version_string()
         << " (Coupling Version: " << stk::coupling::version() << ")" << std::endl;
      os << m_appName << ", color="<<color<<", world rank: " << myWorldRank<<" out of " << numWorldRanks
                      <<", app rank: " << myAppRank << " out of " << numAppRanks << std::endl;
      os << m_appName << ": my root-rank: " << rootRanks.localColorRoot << ", other app's root-rank: " << rootRanks.otherColorRoot;
      std::cout << os.str() << std::endl;
    }

    m_mesh.reset(new mock::StkMesh(splitComm));
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

    std::string otherAppName = m_otherInfo.get_value<std::string>(stk::coupling::AppName, "none");
    if (otherAppName == "Mock-Salinas") {
      m_doingSendTransfer = true;
      m_sendFieldName = "reference-temperature";
    }
    if (otherAppName == "Mock-Sparc") {
      m_doingRecvTransfer = true;
      m_recvFieldName = "heat-transfer-coefficient";
    }

    {
      std::ostringstream os;
      if (m_doingSendTransfer) {
        os << m_appName << ": will send-transfer (field='"<<m_sendFieldName<<"') "
           << " to other app: "<<otherAppName<<std::endl;
      }

      if (m_doingRecvTransfer) {
        os << m_appName << ": will recv-transfer (field='"<<m_recvFieldName<<"') "
           <<" from other app: "<<otherAppName<<std::endl;
      }

      if (m_iAmRootRank) std::cout << os.str() << std::endl;
    }

    m_myInfo.set_value(stk::coupling::CurrentTime, m_myInfo.get_value<double>(stk::coupling::InitialTime));
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
    if (!m_doingSendTransfer && !m_doingRecvTransfer) { return; }

    std::vector<std::pair<std::string,int>> mySendFields;
    std::vector<std::pair<std::string,int>> myRecvFields;
    if (m_doingSendTransfer) {
      mySendFields.push_back(std::make_pair(m_sendFieldName, m_mesh->get_field_size()));
    }
    if (m_doingRecvTransfer) {
      myRecvFields.push_back(std::make_pair(m_recvFieldName, m_mesh->get_field_size()));
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
      m_mesh->set_stk_field_value(m_mesh->get_stk_source_entity_key(), m_sendFieldName, 9.9);
      std::shared_ptr<mock::StkSendAdapter> sendAdapter =
         std::make_shared<mock::StkSendAdapter>(parentComm, *m_mesh, m_sendFieldName);
      std::shared_ptr<mock::EmptyRecvAdapter> recvAdapter;
      m_sendTransfer.reset(new SendTransfer(sendAdapter, recvAdapter, "MockAriaSendTransfer", parentComm));

      m_sendTransfer->coarse_search();
      m_sendTransfer->communication();
      m_sendTransfer->local_search();
    }

    if (m_doingRecvTransfer) {
      std::shared_ptr<mock::EmptySendAdapter> sendAdapter;
      std::shared_ptr<mock::StkRecvAdapter> recvAdapter =
         std::make_shared<mock::StkRecvAdapter>(parentComm, *m_mesh, m_recvFieldName);
      m_recvTransfer.reset(new RecvTransfer(sendAdapter, recvAdapter, "MockAriaRecvTransfer", parentComm));

      m_recvTransfer->coarse_search();
      m_recvTransfer->communication();
      m_recvTransfer->local_search();
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
      if (stk::parallel_machine_rank(m_splitComms.get_split_comm()) == m_mesh->owning_rank()) {
        m_mesh->set_stk_field_value(m_mesh->get_stk_dest_entity_key(), m_recvFieldName, 0.0);
      }
      m_recvTransfer->apply();
      ThrowRequire(m_recvTransfer->meshb()->called_update_values);

      {
        if (stk::parallel_machine_rank(m_splitComms.get_split_comm()) == m_mesh->owning_rank()) {
          std::ostringstream os;
          os << m_appName << " transfer: recvd '"<<m_recvFieldName<<"' value: "
            << m_mesh->get_stk_field_value(m_mesh->get_stk_dest_entity_key(), m_recvFieldName) << std::endl;
          std::cout << os.str();
        }
      }
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
    bool timeToStop = m_iWasToldToStop || m_currentTime >= m_finalTime;
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
  int m_step;
  std::shared_ptr<SendTransfer> m_sendTransfer;
  std::shared_ptr<RecvTransfer> m_recvTransfer;
  bool m_doingSendTransfer;
  bool m_doingRecvTransfer;
  std::string m_sendFieldName;
  std::string m_recvFieldName;
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
