/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <strings.h>

#include <boost/program_options.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>

#include <stk_util/environment/ProgramOptions.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

namespace stk_linsys_usecases {

bool use_case_1_driver( MPI_Comm comm ,
                        const std::string &mesh_filename );

bool use_case_2_driver( MPI_Comm comm ,
                        const std::string &mesh_filename );

bool use_case_3_driver( MPI_Comm comm ,
                        const std::string &mesh_filename );

bool use_case_4_driver( MPI_Comm comm ,
                        const std::string &mesh_filename );

bool use_case_5_driver( MPI_Comm comm ,
                        const std::string &mesh_filename,
                        const std::string & solver_params );

bool use_case_7_driver( MPI_Comm comm ,
                        const std::string &mesh_filename,
                        const std::string & solver_params );

} // namespace stk_linsys_usecases


namespace {
  const char use_case_1[]  = "use_case_1" ;
  const char use_case_2[]  = "use_case_2" ;
  const char use_case_3[]  = "use_case_3" ;
  const char use_case_4[]  = "use_case_4" ;
  const char use_case_5[]  = "use_case_5" ;
  const char use_case_7[]  = "use_case_7" ;
}

namespace bopt = boost::program_options;

int main(int argc, char** argv)
{
  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);
  stk::ParallelMachine comm = use_case_environment.m_comm;

  //----------------------------------
  // Broadcast argc and argv to all processors.

  stk::BroadcastArg b_arg(comm, argc, argv);

  //----------------------------------
  // Process the broadcast command line arguments

  bopt::options_description desc("options");

  stk::get_options_description().add(desc);

  int nthreads = 1;
  int bucket_size = 1000;

  desc.add_options()
    ("help,h",        "produce help message")
    ("mesh,m",        bopt::value<std::string>(), "mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. See GeneratedMesh documentation for more options. Use 'gears' to generate the gears mesh." )
    ("directory,d",   bopt::value<std::string>(), "working directory with trailing '/'" )
    ("output-log,o",  bopt::value<std::string>(), "output log path" )
    ("runtest,r",     bopt::value<std::string>(), "runtest pid file" )
    ("threads",       bopt::value<int>(&nthreads)->default_value(1),   "number of threads")
    ("bucket_size",    bopt::value<int>(&bucket_size)->default_value(1000),   "size of buckets used internally in stk::mesh")
    ("tbb",           "Use Threaded Building Blocks algorithm thread runner")
    ("tpi",           "Use Thread Pool Interface algorithm thread runner")
    ("nonthreaded",   "Run algorithms non-threaded [default]")
    ("solver_params",   "Teuchos parameter XML filename")
    ( use_case_1 ,   "use case 1" )
    ( use_case_2 ,   "use case 2" )
    ( use_case_3 ,   "use case 3" )
    ( use_case_4 ,   "use case 4" )
    ( use_case_5 ,   "use case 5" )
    ( use_case_7 ,   "use case 7" );

  bopt::variables_map vm;
  try {
    bopt::store(bopt::parse_command_line(b_arg.m_argc, b_arg.m_argv, desc), vm);
    bopt::notify(vm);
  }
  catch (std::exception &x) {
    std::cout << x.what() << std::endl;
    std::cout << "Test Failed" << std::endl;
    stk::parallel_machine_finalize();
    std::exit(1);
  }


  //----------------------------------

  if (vm.count("help")) {
    std::cout << desc << "\n";
    std::exit(EXIT_SUCCESS);
  }

  std::string thread_runner = "NonThreaded";
  if (vm.count("tbb"))
    thread_runner = "TBB";
  if (vm.count("tpi"))
    thread_runner = "TPI";
  if (vm.count("nonthreaded"))
    thread_runner = "NonThreaded";

  // Used by multiple use cases.
  std::string working_directory = "";
  std::string in_filename = "";
  std::string in_filetype = "exodusii";
  if (vm.count("mesh")) {
    in_filename = boost::any_cast<std::string>(vm["mesh"].value());;

    if (strncasecmp("gen:", in_filename.c_str(), 4) == 0) {
      // Strip off the 'gen:' prefix and set the type to "generated"
      in_filename = in_filename.substr(4, in_filename.size());
      in_filetype = "generated";
    }

    if (strncasecmp("gears:", in_filename.c_str(), 6) == 0) {
      // Strip off the 'gears:' prefix and set the type to "gears"
      in_filename = in_filename.substr(6, in_filename.size());
      in_filetype = "gears";
    }

    if (vm.count("directory")) {
      working_directory = boost::any_cast<std::string>(vm["directory"].value());
    }
  } else {
    std::cout << "OPTION ERROR: The '--mesh <filename>' option is required for all use cases!\n";
    std::exit(EXIT_FAILURE);
  }

  bool success = false;

  if ( vm.count( use_case_1 ) ) {

    if (in_filetype != "generated") {
      std::cout << "OPTION ERROR: For use case 1, only the 'generated' mesh option is supported!\n";
      std::exit(EXIT_FAILURE);
    }

    success = stk_linsys_usecases::use_case_1_driver( comm , in_filename );
  }
  else if ( vm.count( use_case_2 ) ) {

    if (in_filetype != "generated") {
      std::cout << "OPTION ERROR: For use case 2, only the 'generated' mesh option is supported!\n";
      std::exit(EXIT_FAILURE);
    }

    success = stk_linsys_usecases::use_case_2_driver( comm , in_filename );
  }
  else if ( vm.count( use_case_3 ) ) {

    if (in_filetype != "generated") {
      std::cout << "OPTION ERROR: For use case 3, only the 'generated' mesh option is supported!\n";
      std::exit(EXIT_FAILURE);
    }

    success = stk_linsys_usecases::use_case_3_driver( comm , in_filename );
  }
  else if ( vm.count( use_case_4 ) ) {

    if (in_filetype != "generated") {
      std::cout << "OPTION ERROR: For use case 4, only the 'generated' mesh option is supported!\n";
      std::exit(EXIT_FAILURE);
    }

    success = stk_linsys_usecases::use_case_4_driver( comm , in_filename );
  }
  else if ( vm.count( use_case_5 ) ) {

    if (in_filetype != "generated") {
      std::cout << "OPTION ERROR: For use case 5, only the 'generated' mesh option is supported!\n";
      std::exit(EXIT_FAILURE);
    }

    std::string solver_params = "";

    if( vm.count("solver_params")) {
      solver_params = boost::any_cast<std::string>(vm["solver_params"].value());
    }

    success = stk_linsys_usecases::use_case_5_driver( comm , in_filename, solver_params );
  }
  else if ( vm.count( use_case_7 ) ) {

    if (in_filetype != "generated") {
      std::cout << "OPTION ERROR: For use case 7, only the 'generated' mesh option is supported!\n";
      std::exit(EXIT_FAILURE);
    }

    std::string solver_params = "";

    if( vm.count("solver_params")) {
      solver_params = boost::any_cast<std::string>(vm["solver_params"].value());
    }

    success = stk_linsys_usecases::use_case_7_driver( comm , in_filename, solver_params );
  }
  else {
    std::cout << "OPTION ERROR: Missing a use case selection option.  Use --help option to see valid options.\n";
    std::exit(EXIT_FAILURE);
  }

  use_case::print_status(comm, success);

  return 0;
}

