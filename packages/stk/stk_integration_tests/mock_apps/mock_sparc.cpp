
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Constants.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/CommSplitting.hpp>
#include <stk_coupling/ConfigurationInfo.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <iostream>
#include <sstream>

std::string app_name()
{
  return "Mock-Sparc";
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

  stk::coupling::ConfigurationInfo myInfo;

  myInfo.set_value(stk::coupling::AppName, app_name());
  myInfo.set_value(stk::coupling::TimeSyncMode, stk::coupling::Send);

  constexpr int numberOfSteps = 5;
  myInfo.set_value(stk::coupling::InitialTime, 0.);
  myInfo.set_value(stk::coupling::TimeStep, 0.1);
  myInfo.set_value(stk::coupling::FinalTime, numberOfSteps * myInfo.get_value<double>(stk::coupling::TimeStep));

  const stk::coupling::ConfigurationInfo otherInfo = myInfo.exchange(commWorld, commApp);

  {
    std::ostringstream os;
    os << app_name() << ": other app 'app_name': " << otherInfo.get_value<std::string>(stk::coupling::AppName);
    if(myWorldRank == rootRanks.first) std::cout << os.str() << std::endl;
  }

  stk::coupling::check_sync_mode_consistency(myInfo, otherInfo);

  if(myWorldRank == rootRanks.first) std::cout << "Exchanging Info" << std::endl;
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
  
  double currentTime = myInfo.get_value<double>("current time");
  double timeStep = myInfo.get_value<double>(stk::coupling::TimeStep);
  for(int step = 0; step < numberOfSteps; ++step) {
    // Update State and Do Transfer Using stk_transfer
    currentTime += timeStep;
    myInfo.set_value("time step status", (step == numberOfSteps-1) ? "end" : "proceed");
    myInfo.set_value("current time", currentTime);

    std::cout<<app_name()<<" setting current time to "<<currentTime<<std::endl;

    const stk::coupling::ConfigurationInfo otherInfo = myInfo.exchange(commWorld, commApp);
    currentTime = stk::coupling::choose_value(myInfo, otherInfo, "current time", stk::coupling::Send);

    std::ostringstream os;
    os << app_name() << ": current time: " << currentTime << ", I want to " << myInfo.get_value<std::string>("time step status");
    if(myWorldRank == rootRanks.first) std::cout << os.str() << std::endl;

    if(myInfo.get_value<std::string>("time step status") == "end") break;
    if(otherInfo.has_value<std::string>("time step status") && otherInfo.get_value<std::string>("time step status") == "end") break;
  }

  stk::parallel_machine_finalize();
  return 0;
}
