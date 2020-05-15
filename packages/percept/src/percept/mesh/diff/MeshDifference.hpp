// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_MeshDifference_hpp
#define percept_MeshDifference_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>

#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/RunEnvironment.hpp>


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

#endif
