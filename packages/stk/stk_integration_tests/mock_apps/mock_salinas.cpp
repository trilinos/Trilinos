
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Constants.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_transfer/ReducedDependencyGeometricTransfer.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/Version.hpp>
#include <stk_util/parallel/CouplingVersions.hpp>
#include <stk_util/parallel/CouplingVersions_impl.hpp>
#include "MockUtils.hpp"
#include "StkMesh.hpp"
#include "MockMeshUtils.hpp"
#include "StkSendAdapter.hpp"
#include "StkRecvAdapter.hpp"
#include "RecvInterpolate.hpp"
#include <iostream>
#include <sstream>

class MockSalinas
{
  using RecvTransfer = stk::transfer::ReducedDependencyGeometricTransfer<
                                    mock::RecvInterpolate<mock::StkSendAdapter, mock::StkRecvAdapter>>;
public:
  MockSalinas()
    : m_appName("Mock-Salinas"),
      m_mesh(),
      m_doneFlagName("time step status"),
      m_splitComms(),
      m_otherColor(),
      m_otherAppName(),
      m_iAmRootRank(false),
      m_myInfo(),
      m_otherInfo(),
      m_currentTime(),
      m_finalTime(),
      m_doingRecvTransfer(false),
      m_recvFieldName()
  { }

  ~MockSalinas() = default;

  void read_input_and_setup_split_comms(int argc, char** argv)
  {
    MPI_Comm commWorld = stk::initialize(&argc, &argv);
    int myWorldRank = stk::parallel_machine_rank(commWorld);
    int numWorldRanks = stk::parallel_machine_size(commWorld);

    {
      std::ostringstream os;
      os << m_appName
         << ": STK version: " << stk::version_string()
         << ", Coupling Version: " << stk::util::get_common_coupling_version();
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

    MPI_Comm splitComm = m_splitComms.get_split_comm();
    stk::coupling::PairwiseRanks rootRanks = m_splitComms.get_pairwise_root_ranks(otherColors[0]);
    int myAppRank = stk::parallel_machine_rank(splitComm);
    int numAppRanks = stk::parallel_machine_size(splitComm);
    m_iAmRootRank = myAppRank == 0;

    {
      std::ostringstream os;
      os << m_appName << ": color="<<color<<", my world rank is: " << myWorldRank << " out of " << numWorldRanks
         << ", my app rank is: " << myAppRank << " out of " << numAppRanks << std::endl;
      os << m_appName << ": my root-rank: " << rootRanks.localColorRoot << ", other app's root-rank: " << rootRanks.otherColorRoot;
      std::cout << os.str() << std::endl;
    }

    std::vector<std::string> fieldNames = {"traction", "temperature"};
    mock_utils::read_mesh(splitComm, meshFileName, "surface_1", fieldNames, m_mesh);
    m_currentTime = 0.0;
    m_finalTime = 1.0;
  }

  void communicate_initial_setup()
  {
    m_myInfo = create_sync_info();
    m_myInfo.set_value(stk::coupling::AppName, m_appName);

    m_otherInfo = m_myInfo.exchange(m_splitComms, m_otherColor);

    {
      std::ostringstream os;
      os << m_appName << ": other app 'app_name': " << m_otherInfo.get_value<std::string>(stk::coupling::AppName);
      if (m_iAmRootRank) std::cout << os.str() << std::endl;
    }

    m_otherAppName = m_otherInfo.get_value<std::string>(stk::coupling::AppName, "none");

    if (m_otherAppName=="Mock-Sparc") {
      m_doingRecvTransfer = true;
      m_recvFieldName = "traction";
    }
    if (m_otherAppName=="Mock-Aria") {
      m_doingRecvTransfer = true;
      m_recvFieldName = "temperature";
    }

    {
      if (m_doingRecvTransfer) {
        std::ostringstream os;
        os << m_appName << ": will recv-transfer (field='"<<m_recvFieldName<<"') "
           <<" from other app: "<<m_otherAppName<<std::endl;
        if (m_iAmRootRank) std::cout << os.str() << std::endl;
      }
    }
  }

  bool time_to_stop()
  {
    return m_currentTime >= m_finalTime;
  }

  void communicate_time_step_info()
  {
    m_myInfo = create_sync_info();
    m_otherInfo = m_myInfo.exchange(m_splitComms, m_otherColor);
  }

  bool other_app_is_done()
  {
    bool done = m_otherInfo.get_value<bool>(m_doneFlagName, false);

    {
      if (done && m_iAmRootRank) std::cout << m_appName << " ending because other app ending" << std::endl;
    }

    return done;
  }

  void get_time_from_other_app()
  {
    m_currentTime = m_otherInfo.get_value<double>(stk::coupling::CurrentTime);
    m_finalTime = m_otherInfo.get_value<double>(stk::coupling::FinalTime);

    {
      std::ostringstream os;
      os << m_appName << ": "<<stk::coupling::CurrentTime<<": " << m_currentTime << ", final time: " << m_finalTime;
      if(m_iAmRootRank) std::cout << os.str() << std::endl;
    }
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
    if (!m_doingRecvTransfer) { return; }

    std::vector<std::pair<std::string,int>> mySendFields;
    std::vector<std::pair<std::string,int>> myRecvFields;
    myRecvFields.push_back(std::make_pair(m_recvFieldName, m_mesh->get_field_size()));

    stk::coupling::SyncInfo info = create_sync_info();
    info.set_value("SendFields", mySendFields);
    info.set_value("RecvFields", myRecvFields);

    stk::coupling::SyncInfo otherInfo = info.exchange(m_splitComms, m_otherColor);
    std::vector<std::pair<std::string,int>> otherSendFields = otherInfo.get_value<std::vector<std::pair<std::string,int>>>("SendFields");
    std::vector<std::pair<std::string,int>> otherRecvFields = otherInfo.get_value<std::vector<std::pair<std::string,int>>>("RecvFields");

    check_field_sizes(mySendFields, otherRecvFields);
    check_field_sizes(otherSendFields, myRecvFields);

    MPI_Comm pairwiseComm = m_splitComms.get_pairwise_comm(m_otherColor);
    std::shared_ptr<mock::StkSendAdapter> nullSendAdapter;
    std::shared_ptr<mock::StkRecvAdapter> recvAdapter =
       std::make_shared<mock::StkRecvAdapter>(pairwiseComm, *m_mesh, m_recvFieldName);
    m_recvTransfer.reset(new RecvTransfer(nullSendAdapter, recvAdapter, "MockSalinasRecvTransfer", pairwiseComm));

    m_recvTransfer->coarse_search();
    m_recvTransfer->communication();
    m_recvTransfer->local_search();
  }

  void perform_transfers()
  {
    if (m_doingRecvTransfer) {
      m_mesh->set_stk_field_values(m_recvFieldName, 0.0);
      m_recvTransfer->apply();
      STK_ThrowRequire(m_recvTransfer->meshb()->called_update_values);

      const double expectedFieldValue = (m_otherAppName=="Mock-Sparc" ? 4.4 : 9.9);
      const bool valuesMatch = m_mesh->verify_stk_field_values(m_recvFieldName, expectedFieldValue);
      STK_ThrowRequireMsg(valuesMatch, "Mock-Salinas error, field-values are not correct after transfer");
    }
  }

  void do_physics_solve()
  {
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
  std::shared_ptr<mock::StkMesh> m_mesh;
  const std::string m_doneFlagName;

  stk::coupling::SplitComms m_splitComms;
  int m_otherColor;
  std::string m_otherAppName;
  bool m_iAmRootRank;

  stk::coupling::SyncInfo m_myInfo;
  stk::coupling::SyncInfo m_otherInfo;

  double m_currentTime;
  double m_finalTime;
  std::shared_ptr<RecvTransfer> m_recvTransfer;
  bool m_doingRecvTransfer;
  std::string m_recvFieldName;
};

int main(int argc, char** argv)
{
  {
    MockSalinas app;
    app.read_input_and_setup_split_comms(argc, argv);

    if (app.get_number_of_other_coupled_apps() == 1) {
      app.communicate_initial_setup();
      app.setup_fields_and_transfers();

      while (!app.time_to_stop()) {
        app.communicate_time_step_info();

        if (app.other_app_is_done()) break;

        app.get_time_from_other_app();
        app.do_physics_solve();
        app.perform_transfers();
      }
    }
  }

  stk::finalize();

  return 0;
}
