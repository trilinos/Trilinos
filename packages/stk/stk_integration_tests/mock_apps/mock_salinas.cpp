
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/CommSplitting.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include "ConfigurationInfo.hpp"
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
    os<<app_name()<<": my world rank is: "<<myWorldRank<<" out of "<<numWorldRanks;
    std::cout<<os.str()<<std::endl;
  }
  bool foundColor = false;
  int color = -1;
  std::tie(foundColor, color) = stk::coupling::get_app_color(argc, argv);
  ThrowRequireMsg(foundColor, "Need to specify '--app-color <int>' for mpmd run");

  stk::ParallelMachine commApp = stk::coupling::split_comm(commWorld, color);
  std::pair<int,int> rootRanks = stk::coupling::calc_my_root_and_other_root_ranks(commWorld, commApp);
  int myAppRank = stk::parallel_machine_rank(commApp);
  int numAppRanks = stk::parallel_machine_size(commApp);

  {
    std::ostringstream os;
    os<<app_name()<<": my app rank is: "<<myAppRank<<" out of "<<numAppRanks<<std::endl;
    os<<app_name()<<": my root-rank: "<<rootRanks.first<<", other app's root-rank: "<<rootRanks.second;
    std::cout<<os.str()<<std::endl;
  }

  mock_coupling::ConfigurationInfo myConfigInfo;

  myConfigInfo.svals["app_name"] = app_name();

  const mock_coupling::ConfigurationInfo otherConfigInfo = myConfigInfo.exchange(commWorld, commApp);

  {
    std::ostringstream os;
    os<<app_name()<<": other app 'app_name': "<<otherConfigInfo.svals.at("app_name");
    if(myWorldRank == rootRanks.first) std::cout<<os.str()<<std::endl;
  }

  stk::parallel_machine_finalize();
  return 0;
}
