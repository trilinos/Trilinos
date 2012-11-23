#ifndef stk_percept_MeshDifference_hpp
#define stk_percept_MeshDifference_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>

#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/RunEnvironment.hpp>


namespace stk
{
  namespace percept
  {

    class MeshDifference
    {
    public:
      MeshDifference() :mesh_opt_1(""),mesh_opt_2(""),print_all_field_diffs(0){}
      void run(int argc,  char** argv);

      void process_options(RunEnvironment& re);
      // command line
      std::string mesh_opt_1, mesh_opt_2;
      int print_all_field_diffs;
    };

      
  }//namespace percept
}//namespace stk

#endif
