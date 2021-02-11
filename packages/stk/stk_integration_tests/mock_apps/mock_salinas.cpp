
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
  return "Mock-Salinas";
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

  stk::coupling::SyncInfo myInfo;

  myInfo.set_value(stk::coupling::AppName, app_name());

  myInfo.set_value(stk::coupling::InitialTime, 0.0);
  myInfo.set_value("current time", 0.0);
  myInfo.set_value(stk::coupling::TimeStep, 0.0);
  myInfo.set_value(stk::coupling::FinalTime, 0.0);

  stk::coupling::SyncInfo otherInfo = myInfo.exchange(commWorld, commApp);

  {
    std::ostringstream os;
    os << app_name() << ": other app 'app_name': " << otherInfo.get_value<std::string>(stk::coupling::AppName);
    if(myWorldRank == rootRanks.first) std::cout << os.str() << std::endl;
  }

  double currentTime = 0.0;
  double finalTime = 1.0;
  while (currentTime < finalTime) {

    otherInfo = myInfo.exchange(commWorld, commApp);

    currentTime = otherInfo.get_value<double>("current time");
    finalTime = otherInfo.get_value<double>(stk::coupling::FinalTime);

    std::ostringstream os;
    os << app_name() << ": current time: " << currentTime << ", final time: " << finalTime;
    if(myWorldRank == rootRanks.first) std::cout << os.str() << std::endl;
    if(otherInfo.has_value<std::string>("time step status") && otherInfo.get_value<std::string>("time step status") == "end") {
      if(myWorldRank == rootRanks.first) std::cout << app_name() << " ending because other app ending" << std::endl;
      break;
    }
  }

  stk::parallel_machine_finalize();
  return 0;
}
