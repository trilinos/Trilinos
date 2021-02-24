
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Constants.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/CommSplitting.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_transfer/ReducedDependencyGeometricTransfer.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include "StkMesh.hpp"
#include "StkRecvAdapter.hpp"
#include "EmptySendAdapter.hpp"
#include "RecvInterpolate.hpp"
#include <iostream>
#include <sstream>

class MockSalinas
{
  using RecvTransfer = stk::transfer::ReducedDependencyGeometricTransfer<
                                    mock::RecvInterpolate<mock::EmptySendAdapter, mock::StkRecvAdapter>>;
public:
  MockSalinas()
    : m_appName("Mock-Salinas"),
      m_mesh(),
      m_doneFlagName("time step status"),
      m_iAmRootRank(false),
      m_doingRecvTransfer(false),
      m_recvFieldName()
  { }

  ~MockSalinas()
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
    int numAppRanks = stk::parallel_machine_size(m_commApp);
    m_iAmRootRank = myWorldRank == rootRanks.first;

    {
      std::ostringstream os;
      os << m_appName << ": color="<<color<<", my app rank is: " << myAppRank << " out of " << numAppRanks << std::endl;
      os << m_appName << ": my root-rank: " << rootRanks.first << ", other app's root-rank: " << rootRanks.second;
      std::cout << os.str() << std::endl;
    }

    m_mesh.reset(new mock::StkMesh(m_commApp));
    m_currentTime = 0.0;
    m_finalTime = 1.0;
  }

  void communicate_initial_setup()
  {
    m_myInfo = create_sync_info();
    m_myInfo.set_value(stk::coupling::AppName, m_appName);

    m_otherInfo = m_myInfo.exchange(m_commWorld, m_commApp);

    {
      std::ostringstream os;
      os << m_appName << ": other app 'app_name': " << m_otherInfo.get_value<std::string>(stk::coupling::AppName);
      if (m_iAmRootRank) std::cout << os.str() << std::endl;
    }

    std::string otherAppName = m_otherInfo.get_value<std::string>(stk::coupling::AppName, "none");

    if (otherAppName=="Mock-Sparc") {
      m_doingRecvTransfer = true;
      m_recvFieldName = "traction";
    }
    if (otherAppName=="Mock-Aria") {
      m_doingRecvTransfer = true;
      m_recvFieldName = "temperature";
    }

    {
      if (m_doingRecvTransfer) {
        std::ostringstream os;
        os << m_appName << ": will recv-transfer (field='"<<m_recvFieldName<<"') "
           <<" from other app: "<<otherAppName<<std::endl;
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
    m_otherInfo = m_myInfo.exchange(m_commWorld, m_commApp);
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

  void setup_fields_and_transfers()
  {
    if (!m_doingRecvTransfer) { return; }
    std::shared_ptr<mock::EmptySendAdapter> sendAdapter;
    std::shared_ptr<mock::StkRecvAdapter> recvAdapter =
       std::make_shared<mock::StkRecvAdapter>(m_commWorld, *m_mesh, m_recvFieldName);
    m_recvTransfer.reset(new RecvTransfer(sendAdapter, recvAdapter, "MockSalinasRecvTransfer", m_commWorld));

    m_recvTransfer->coarse_search();
    m_recvTransfer->communication();
    m_recvTransfer->local_search();
  }

  void perform_transfers()
  {
    if (m_doingRecvTransfer) {
      if (stk::parallel_machine_rank(m_commApp) == m_mesh->owning_rank()) {
        m_mesh->set_stk_field_value(m_mesh->get_stk_dest_entity_key(), m_recvFieldName, 0.0);
      }
      m_recvTransfer->apply();
      ThrowRequire(m_recvTransfer->meshb()->called_update_values);
      if (stk::parallel_machine_rank(m_commApp) == m_mesh->owning_rank()) {
        std::ostringstream os;
        os << m_appName << " transfer: recvd '"<<m_recvFieldName<<"' value: "
          << m_mesh->get_stk_field_value(m_mesh->get_stk_dest_entity_key(), m_recvFieldName) << std::endl;
        std::cout << os.str();
      }
    }
  }

  void do_physics_solve()
  {
  }

private:
  stk::coupling::SyncInfo create_sync_info()
  {
    return stk::coupling::SyncInfo(m_appName);
  }

  const std::string m_appName;
  std::shared_ptr<mock::StkMesh> m_mesh;
  const std::string m_doneFlagName;

  stk::ParallelMachine m_commWorld;
  stk::ParallelMachine m_commApp;
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
  MockSalinas app;

  app.read_input_and_setup_split_comms(argc, argv);
  app.communicate_initial_setup();
  app.setup_fields_and_transfers();

  while (!app.time_to_stop()) {
    app.communicate_time_step_info();

    if (app.other_app_is_done()) break;

    app.get_time_from_other_app();
    app.do_physics_solve();
    app.perform_transfers();
  }

  return 0;
}
