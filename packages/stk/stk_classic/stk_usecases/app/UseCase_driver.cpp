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

namespace stk_classic {
namespace app {

bool use_case_7_driver( MPI_Comm comm ,
                        bool performance_test ,
                        const std::string &mesh_filename ,
                        int num_threads ,
                        const std::string& thread_runner );

bool use_case_blas_driver(MPI_Comm comm,
                        int num_threads,
                        int num_trials,
                        const std::string &working_directory,
                        const std::string &mesh_filename,
                        const std::string &mesh_type,
                        const std::string &thread_runner,
                        int chunk_size,
                        bool performance_test);

bool use_case_14_driver(MPI_Comm comm,
                        int num_threads,
                        int num_trials,
                        const std::string &working_directory,
                        const std::string &mesh_filename,
                        const std::string &mesh_type,
                        const std::string &thread_runner,
                        int chunk_size,
                        bool performance_test);

bool use_case_14a_driver(MPI_Comm comm,
                         int num_threads,
                         int num_trials,
                         const std::string &working_directory,
                         const std::string &mesh_filename,
                         const std::string &mesh_type,
                         const std::string &thread_runner,
                         int chunk_size,
                         bool performance_test);

bool use_case_24_driver( MPI_Comm comm,
                         const std::string &working_directory,
                         const std::string &mesh_filename,
                         const std::string &outputName );

} // namespace app
} // namespace stk_classic


namespace {
  const char use_case_7[]  = "use_case_7" ;
  const char use_case_blas[]  = "use_case_blas" ;
  const char use_case_14[]  = "use_case_14" ;
  const char use_case_14a[] = "use_case_14a" ;
  const char use_case_24[]  = "use_case_24" ;
}

namespace bopt = boost::program_options;

int main(int argc, char** argv)
{
  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);
  stk_classic::ParallelMachine comm = use_case_environment.m_comm;

  //----------------------------------
  // Broadcast argc and argv to all processors.

  stk_classic::BroadcastArg b_arg(comm, argc, argv);

  //----------------------------------
  // Process the broadcast command line arguments

  bopt::options_description desc("options");

  stk_classic::get_options_description().add(desc);

  int num_trials = 0;
  int nthreads = 1;
  int bucket_size = 1000;

  desc.add_options()
    ("help,h",        "produce help message")
    ("performance",   "run performance test [14/14a/blas only]")
    ("mesh,m",        bopt::value<std::string>(), "mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. See GeneratedMesh documentation for more options. Use 'gears' to generate the gears mesh." )
    ("directory,d",   bopt::value<std::string>(), "working directory with trailing '/'" )
    ("output-log,o",  bopt::value<std::string>(), "output log path" )
    ("runtest,r",     bopt::value<std::string>(), "runtest pid file" )
    ("threads",       bopt::value<int>(&nthreads)->default_value(1),   "number of threads [14/14a/blas only]")
    ("trials",        bopt::value<int>(&num_trials)->default_value(1),   "number of trials (execute loops) [14/14a/blas only]")
    ("bucket_size",    bopt::value<int>(&bucket_size)->default_value(1000),   "size of buckets used internally in stk_classic::mesh [14/14a/blas only]")
    ("tbb",           "Use Threaded Building Blocks algorithm thread runner [14/14a/blas only]")
    ("tpi",           "Use Thread Pool Interface algorithm thread runner [14/14a/blas only]")
    ("nonthreaded",   "Run algorithms non-threaded [default] [14/14a/blas only]")
    ( use_case_7 ,   "use case 7" )
    ( use_case_blas , "use case blas (fill, axpby, dot, norm)" )
    ( use_case_14 ,   "use case 14 (hex internal force algorithm)" )
    ( use_case_14a ,  "use case 14a (hex internal force algorithm)" )
    ( use_case_24 ,   "use case 24 " );

  bopt::variables_map vm;
  try {
    bopt::store(bopt::parse_command_line(b_arg.m_argc, b_arg.m_argv, desc), vm);
    bopt::notify(vm);
  }
  catch (std::exception &x) {
    std::cout << x.what() << std::endl;
    std::exit(1);
  }


  //----------------------------------

  if (vm.count("help")) {
    std::cout << desc << "\n";
    std::exit(EXIT_SUCCESS);
  }

  bool run_performance_test = vm.count( "performance" ) > 0;

  std::string thread_runner = "NonThreaded";
  if (vm.count("tbb"))
    thread_runner = "TBB";
  if (vm.count("tpi"))
    thread_runner = "TPI";
//   if (vm.count("dgb"))
//     thread_runner = "DGB";
  if (vm.count("nonthreaded"))
    thread_runner = "NonThreaded";

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

  if ( vm.count( use_case_7 ) ) {

    if (in_filetype != "generated") {
      std::cout << "OPTION ERROR: For use case 7, only the 'generated' mesh option is supported!\n";
      std::exit(EXIT_FAILURE);
    }

    success = stk_classic::app::use_case_7_driver( comm , run_performance_test , in_filename , nthreads, thread_runner );
  }

  else if ( vm.count( use_case_blas ) ) {
    success = stk_classic::app::use_case_blas_driver(comm, nthreads, num_trials, working_directory,
                                             in_filename, in_filetype, thread_runner, bucket_size,
                                             run_performance_test);
  }

  else if ( vm.count( use_case_14 ) ) {
    success = stk_classic::app::use_case_14_driver(comm, nthreads, num_trials, working_directory,
                                           in_filename, in_filetype, thread_runner, bucket_size,
                                           run_performance_test);
  }

  else if ( vm.count( use_case_14a ) ) {
    success = stk_classic::app::use_case_14a_driver(comm, nthreads, num_trials, working_directory,
                                            in_filename, in_filetype, thread_runner, bucket_size,
                                            run_performance_test);
  }

  else if ( vm.count( use_case_24 ) ) {
    if (in_filetype == "generated") {
      std::cout << "OPTION ERROR: The 'generated mesh option' option is not supported for use case 24!\n";
      std::exit(EXIT_FAILURE);
    }

    if (in_filetype == "gears") {
      std::cout << "OPTION ERROR: The 'gears mesh option' option is not supported for use case 24!\n";
      std::exit(EXIT_FAILURE);
    }

    success = stk_classic::app::use_case_24_driver( comm, working_directory, in_filename, "1mCube.e" );
  }
  else {
    std::cout << "OPTION ERROR: Missing a use case selection option.  Use --help option to see valid options.\n";
    std::exit(EXIT_FAILURE);
  }

  use_case::print_status(comm, success);

  return 0;
}

