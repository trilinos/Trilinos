/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

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

