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
      MeshDifference() {}
      void run(int argc,  char** argv);

      void process_options(RunEnvironment& re);
      // command line
      std::string mesh_opt_1, mesh_opt_2;
    };

      
  }//namespace percept
}//namespace stk

#endif
