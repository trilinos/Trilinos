
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
  return "Mock-Fuego";
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

  stk::coupling::ConfigurationInfo otherInfo = myInfo.exchange(commWorld, commApp);

  {
    std::ostringstream os;
    os << app_name() << ": other app 'app_name': " << otherInfo.get_value<std::string>(stk::coupling::AppName);
    if(myWorldRank == rootRanks.first) std::cout << os.str() << std::endl;
  }


  constexpr int numberOfSteps = 10;

  double initialTime = 0.0;
  double timeStep = 0.005;
  double finalTime = numberOfSteps * timeStep;
  for(int step = 0; step <= numberOfSteps; ++step) {
    double currentTime = initialTime + step*timeStep;

    myInfo.set_value("step", step);
    myInfo.set_value("current time", currentTime);
    myInfo.set_value(stk::coupling::FinalTime, finalTime);

    otherInfo = myInfo.exchange(commWorld, commApp);

    std::ostringstream os;
    os << app_name() << ": current time: " << currentTime;
    if(myWorldRank == rootRanks.first) std::cout << os.str() << std::endl;
  }
  
  {
    std::ostringstream os;
    os << app_name() << " finished, final time: " << finalTime << std::endl;
    if(myWorldRank == rootRanks.first) std::cout << os.str();
  }

  stk::parallel_machine_finalize();
  return 0;
}
