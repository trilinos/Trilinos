// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <iostream>

#include <boost/program_options.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

namespace stk_examples {

  void example_io_1( stk::ParallelMachine comm,
                     const std::string& in_filename,
                     const std::string & out_filename );

  void example_io_2( stk::ParallelMachine comm,
                     const std::string& in_filename,
                     const std::string & out_filename );

  void use_case_5_write_mesh( stk::ParallelMachine , 
                              const std::string & out_filename );
} // namespace stk_examples

namespace bopt = boost::program_options;

int main(int argc, char** argv)
{
  //----------------------------------
  // Broadcast argc and argv to all processors.

  stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

  stk::BroadcastArg b_arg(comm, argc, argv);

  //----------------------------------
  // Process the broadcast command line arguments
  
  bopt::options_description desc("options");

  desc.add_options()
    ("help,h",        "produce help message")
    ("usecase",      bopt::value<std::string>(), "use case, valid = all, 1, 2, requires -mesh option" )
    ("mesh",         bopt::value<std::string>(), "mesh file" )
    ("directory,d",  bopt::value<std::string>(), "working directory" )
    ("output-log,o", bopt::value<std::string>(), "output log path" )
    ("runtest,r",    bopt::value<std::string>(), "runtest pid file" );

  stk::get_options_description().add(desc);
    
  bopt::variables_map &vm = stk::get_variables_map();  
  try {
    bopt::store(bopt::parse_command_line(b_arg.m_argc, b_arg.m_argv, desc), vm);
    bopt::notify(vm);
  }
  catch (std::exception &x) {
    stk::RuntimeDoomedSymmetric() << x.what();
    std::exit(1);
  }
  
  if (vm.count("help")) {
    std::cout << desc << "\n";
    std::exit(EXIT_SUCCESS);
  }

  //----------------------------------

  std::string usecase = "all";
  if (vm.count("usecase")) {
    usecase = boost::any_cast<std::string>(vm["usecase"].value());
  }

  if ( vm.count("mesh") ) {
    std::string in_filename = boost::any_cast<std::string>(vm["mesh"].value());

    if (usecase == "all" || usecase == "1") {
      std::string out_filename = in_filename + ".out-1";
      stk_examples::example_io_1(comm, in_filename, out_filename );
    }
    
    if (usecase == "all" || usecase == "2") {
      std::string out_filename = in_filename + ".out-2";
      stk_examples::example_io_2(comm, in_filename, out_filename );
    }
  }
  else {
    std::string in_filename( "example_use_case_5.g" );

    stk_examples::use_case_5_write_mesh( comm , in_filename );

    std::string out_filename ;
    out_filename.append( "out_" );
    out_filename.append( in_filename );

    stk_examples::example_io_1(comm, in_filename, out_filename );
  }

  stk::parallel_machine_finalize();

  return 0;
}

