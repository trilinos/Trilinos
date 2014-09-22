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

#include <stk_mesh/base/Comm.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/PrintTimer.hpp>

bool use_case_5_driver(stk::ParallelMachine  comm);

bool use_case_6_driver(stk::ParallelMachine  comm,
                      const std::string &working_directory,
                      const std::string &range_mesh,
                      const std::string &range_mesh_type,
                      const std::string &domain_mesh,
                      const std::string &domain_mesh_type);

bool use_case_7_driver(stk::ParallelMachine  comm,
                      const std::string &working_directory,
                      const std::string &domain_mesh,
                      const std::string &domain_mesh_type);

bool use_case_8_driver(stk::ParallelMachine  comm,
                      const std::string &working_directory,
                      const std::string &domain_mesh,
                      const std::string &domain_mesh_type);

namespace bopt = boost::program_options;

int main(int argc, char **argv)
{
  stk::diag::Timer timer("Transfer Use Cases",
                          use_case::TIMER_TRANSFER,
                          use_case::timer());
  use_case::timerSet().setEnabledTimerMask(use_case::TIMER_ALL);

  bool status = true;

  std::string range_mesh        = "9x9x9";
  std::string range_mesh_type   = "generated";
  std::string domain_mesh       = "8x8x8";
  std::string domain_filetype   = "generated";

  //----------------------------------
  // Process the broadcast command line arguments

  bopt::options_description desc("Transfer use case options");

  // NOTE: Options --directory --output-log --runtest are handled/defined in UseCaseEnvironment
  desc.add_options()
    ("range_mesh",    bopt::value<std::string>(&range_mesh),
     "range mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. "
     "See GeneratedMesh documentation for more options. Use 'gears' to generate the gears mesh." )
    ("domain_mesh",   bopt::value<std::string>(&domain_mesh),
     "domain mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. "
     "See GeneratedMesh documentation for more options. Use 'gears' to generate the gears mesh." )
    ("use_case_5",   "transfer use case 5 -- node (range) to node    (domain) copy     search." )
    ("use_case_6",   "transfer use case 6 -- node (range) to node    (domain) copy     search." )
    ("use_case_7",   "transfer use case 7 -- node (range) to node    (domain) copy     search." )
    ("use_case_8",   "transfer use case 8 -- node (range) to node    (domain) copy     search." )
    ("offset",       bopt::value<double>()->default_value(0.1), "transfer use case 3 offset" )
    ("scale",        bopt::value<double>()->default_value(0.0), "transfer use case 3 scale." )
    ;

  stk::get_options_description().add(desc);

  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);
  const std::string working_directory = use_case_environment.m_workingDirectory;

  bopt::variables_map &vm = stk::get_variables_map();

  stk::ParallelMachine comm = use_case_environment.m_comm;

  if (vm.count("use_case_5")) {
     status = use_case_5_driver(comm);
  }
  if (vm.count("use_case_6")) {
     status = status && use_case_6_driver(comm, working_directory, range_mesh, range_mesh_type, domain_mesh, domain_filetype);
  }
  if (vm.count("use_case_7")) {
     status = status && use_case_7_driver(comm, working_directory, domain_mesh, domain_filetype);
  }
  if (vm.count("use_case_8")) {
     status = status && use_case_8_driver(comm, working_directory, domain_mesh, domain_filetype);
  }

  timer.stop();

  const bool collective_result = use_case::print_status(comm, status);
  const int return_code = collective_result ? 0 : -1;
  return return_code;
}
