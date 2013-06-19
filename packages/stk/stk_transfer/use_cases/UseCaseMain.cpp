/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */      
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */      
/*  United States Government.                                             */      
/*------------------------------------------------------------------------*/

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

  timer.stop();

  const bool collective_result = use_case::print_status(comm, status);
  const int return_code = collective_result ? 0 : -1;
  return return_code;
}
