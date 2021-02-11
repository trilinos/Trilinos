
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Constants.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/CommSplitting.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <iostream>
#include <sstream>

std::string app_name()
{
  return "Mock-Aria";
}

double calculate_time_step()
{
  return 0.1;
}

int main(int argc, char** argv)
{
  stk::ParallelMachine commWorld = stk::parallel_machine_init(&argc, &argv);
  int myWorldRank = stk::parallel_machine_rank(commWorld);
  int numWorldRanks = stk::parallel_machine_size(commWorld);

  {
    std::ostringstream os;
    os << app_name() << ": my world rank is: " << myWorldRank << " out of " << numWorldRanks;
    std::cout << os.str() << std::endl;
  }

  int defaultColor = stk::coupling::string_to_color(app_name());
  int color = stk::get_command_line_option(argc, argv, "app-color", defaultColor);
  std::string defaultSyncMode = "Send";
  std::string syncModeString = stk::get_command_line_option<std::string>(argc, argv, "sync-mode", defaultSyncMode);
  stk::coupling::SyncMode syncMode = stk::coupling::string_to_sync_mode(syncModeString);

  stk::ParallelMachine commApp = stk::coupling::split_comm(commWorld, color);
  std::pair<int,int> rootRanks = stk::coupling::calc_my_root_and_other_root_ranks(commWorld, commApp);
  int myAppRank = stk::parallel_machine_rank(commApp);
  int numAppRanks = stk::parallel_machine_size(commApp);

  {
    std::ostringstream os;
    os << app_name() << ": color="<<color<<", my app rank is: " << myAppRank << " out of " << numAppRanks << std::endl;
    os << app_name() << ": my root-rank: " << rootRanks.first << ", other app's root-rank: " << rootRanks.second;
    std::cout << os.str() << std::endl;
  }

  stk::coupling::SyncInfo myInfo;

  myInfo.set_value(stk::coupling::AppName, app_name());
  myInfo.set_value(stk::coupling::InitialTime, 0.);
  double timeStep = calculate_time_step();
  myInfo.set_value(stk::coupling::TimeStep, timeStep);
  constexpr int numberOfSteps = 5;
  myInfo.set_value(stk::coupling::FinalTime, numberOfSteps * myInfo.get_value<double>(stk::coupling::TimeStep));
  myInfo.set_value(stk::coupling::TimeSyncMode, syncMode);

  if(myWorldRank == rootRanks.first) std::cout << "Exchanging Info" << std::endl;
  stk::coupling::SyncInfo otherInfo = myInfo.exchange(commWorld, commApp);

  check_sync_mode_consistency(myInfo, otherInfo);

  {
    std::ostringstream os;
    os << app_name() << ": other app 'app_name': " << otherInfo.get_value<std::string>(stk::coupling::AppName);
    if(myWorldRank == rootRanks.first) std::cout << os.str() << std::endl;
  }

  {
    ThrowRequireMsg(stk::coupling::check_consistency<double>(myInfo, otherInfo, stk::coupling::InitialTime), 
                     app_name() << ": initial time is inconsistent with " << otherInfo.get_value<std::string>(stk::coupling::AppName));
    ThrowRequireMsg(otherInfo.has_value<double>(stk::coupling::TimeStep), 
                     app_name() << ": other app ("<< otherInfo.get_value<std::string>(stk::coupling::AppName)<<") doesn't have time step");
    ThrowRequireMsg(otherInfo.has_value<double>(stk::coupling::FinalTime), 
                     app_name() << ": other app ("<< otherInfo.get_value<std::string>(stk::coupling::AppName)<<") doesn't have final time");
  }

  if(myWorldRank == rootRanks.first) std::cout << "Setting Current Time" << std::endl;
  stk::coupling::copy_value<double>(myInfo, stk::coupling::InitialTime, myInfo, "current time");
  
  double currentTime = 0.0;
  std::string timeStepStatus = "proceed";
  for(int step = 0; step < numberOfSteps; ++step) {
    // Update State and Do Transfer Using stk_transfer
    timeStep = calculate_time_step();
    if (step == numberOfSteps-2) timeStepStatus = "end";

    myInfo.set_value(stk::coupling::TimeStep, timeStep); 
    myInfo.set_value("current time", currentTime); 
    myInfo.set_value("time step status", timeStepStatus);
    otherInfo = myInfo.exchange(commWorld, commApp);
    currentTime = stk::coupling::choose_value(myInfo, otherInfo, "current time", syncMode);

    std::ostringstream os;
    os << app_name() << ": current time: " << currentTime << ", time step: " << timeStep << ". I want to " << myInfo.get_value<std::string>("time step status");
    if(myWorldRank == rootRanks.first) std::cout << os.str() << std::endl;

    if(myInfo.get_value<std::string>("time step status") == "end") break;
    if(otherInfo.has_value<std::string>("time step status") && otherInfo.get_value<std::string>("time step status") == "end") break;

    currentTime += timeStep;
  }
  

  stk::parallel_machine_finalize();
  return 0;
}
