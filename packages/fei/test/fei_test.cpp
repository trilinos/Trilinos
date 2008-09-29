/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <base/fei_iostream.hpp>

#include <base/fei_utils.hpp>
#include <test_utils/fei_test_utils.hpp>

//Each of the following 'test_*.hpp' headers declares a class which
//specializes the 'tester' interface. Each tester specialization
//provides one or more unit-tests. These are used in the function
//execute_unit_tests below in this file.

#include <test_utils/test_Set.hpp>
#include <test_utils/test_Database.hpp>
#include <test_utils/test_EqnBuffer.hpp>
#include <test_utils/test_EqnCommMgr.hpp>
#include <test_utils/test_SSMat.hpp>
#include <test_utils/test_Algebraic.hpp>
#include <test_utils/test_AztecWrappers.hpp>
#include <test_utils/test_misc.hpp>
#include <test_utils/test_Factory.hpp>
#include <test_utils/test_SNL_FEI_Structure.hpp>
#include <test_utils/test_FEI_Implementation.hpp>
#include <test_utils/test_FEI_Impl.hpp>
#include <test_utils/test_Tables.hpp>
#include <test_utils/test_PointBlockMap.hpp>
#include <test_utils/test_VectorSpace.hpp>
#include <test_utils/test_MatrixGraph.hpp>
#include <test_utils/test_Vector.hpp>
#include <test_utils/test_Utils.hpp>
#include <test_utils/test_Matrix.hpp>
#include <test_utils/test_LinearSystem.hpp>
#include <test_utils/test_benchmarks.hpp>
#include <test_utils/test_FEI.hpp>

///////////// end of 'tester' specializations ///////////////


#include <fei_CommUtils.hpp>
#include <snl_fei_Utils.hpp>
#include <fei_ParameterSet.hpp>

#include "test_utils/poisson_beam_mains.hpp"

#undef fei_file
#define fei_file "fei_test.cpp"

#include <fei_ErrMacros.hpp>

void execute_benchmarks(MPI_Comm comm);

int test_library_plugins(MPI_Comm comm);

int execute_named_test(const std::string& testname,
		       int argc, char** argv, MPI_Comm comm);

void read_input_and_execute_fullsystem_tests(const std::string& filename,
                                             int argc, char** argv, MPI_Comm comm);

void execute_unit_tests(const std::string& path, MPI_Comm comm);

void execute_fullsystem_tests(MPI_Comm comm, const std::string& path,
                              fei::ParameterSet& name_numproc_pairs);

#ifndef FEI_SER
int split_four_procs_into_two_groups(MPI_Comm comm,
				     MPI_Comm& newcomm1, MPI_Group& newgroup1,
				     MPI_Comm& newcomm2, MPI_Group& newgroup2);
#endif

//--------------------------------------------------
// main
//--------------------------------------------------
int main(int argc, char** argv) {

  double start_time = fei::utils::cpu_time();

  int numProcs, localProc;
  CHK_ERR( fei_test_utils::initialize_mpi(argc, argv, localProc, numProcs) );

  if (localProc == 0) {
    FEI_COUT << "FEI version: " << fei::utils::version() << FEI_ENDL;
  }

  //Check whether the -test flag was used.

  std::string testname = fei_test_utils::get_arg_value("-test", argc, argv);
  fei_test_utils::broadcast_string(MPI_COMM_WORLD, 0, testname);

  int errcode = 0;

  if (!testname.empty()) {
    errcode = execute_named_test(testname, argc, argv, MPI_COMM_WORLD);

    if (localProc == 0 && errcode == 0) {
      FEI_COUT << localProc << ": FEI test successful" << FEI_ENDL;
    }
  }
  else {
    //...else since testname is empty,
    //Check whether the -i flag was used to specify an input file.
    //(construct_filename constructs a file name using a combination of
    // '-i <file>' and optional '-d <path>' flags.)
    std::string filename;
    if (localProc == 0) {
      filename = fei_test_utils::construct_filename(argc, argv);
    }
    fei_test_utils::broadcast_string(MPI_COMM_WORLD, 0, filename);

    try {
      read_input_and_execute_fullsystem_tests(filename, argc, argv, MPI_COMM_WORLD);
    }
    catch(std::runtime_error& exc) {
      FEI_CERR << "caught fei test error: "<<exc.what() << FEI_ENDL;
      errcode = -1;
    }
  }

  int global_err_code = 0;
  fei::GlobalSum(MPI_COMM_WORLD, errcode, global_err_code);

#ifndef FEI_SER
  if (MPI_Finalize() != MPI_SUCCESS) ERReturn(-1);
#endif

  if (localProc == 0) {
    double elapsedTime = fei::utils::cpu_time() - start_time;
    FEI_COUT << "Proc0 CPU  time: " << elapsedTime << " seconds." << FEI_ENDL;
  }

  return(global_err_code);
}

void read_input_and_execute_fullsystem_tests(const std::string& filename,
                                             int argc, char** argv,
                                             MPI_Comm comm)
{
  //We'll run some 'full-system' FEI tests, if any are contained in the
  //specified filename.
  //First, fill an array with names that we read from filename. These names
  //are input-file-names, followed by an int specifying the number of processors
  //the test should be run on.

  std::vector<std::string> inputFileNames;
  if (!filename.empty()) {
    const char* filename_c_str = filename.c_str();
    fei_test_utils::read_input_file(filename_c_str, comm, inputFileNames);
    fei::ParameterSet name_numproc_pairs;
    fei::utils::parse_strings(inputFileNames, " ", name_numproc_pairs);

    std::string path = fei_test_utils::get_arg_value("-d", argc, argv);

    //make sure every processor has the path string.
    fei_test_utils::broadcast_string(comm, 0, path);

    try {
      execute_fullsystem_tests(comm, path, name_numproc_pairs);
    }
    catch(std::runtime_error& exc) {
      FEI_CERR << "caught fei test error: "<<exc.what() << FEI_ENDL;
      throw;
    }

#ifndef FEI_SER
    //Next, if we're running on 4 procs, create two new communicators each
    //representing 2 procs, and run some 2-proc tests using one of those
    //communicators.
    int numProcs = fei::numProcs(comm);
    int localProc = fei::localProc(comm);
    if (numProcs == 4) {
      fei::Barrier(comm);
      if (localProc == 0) {
        FEI_COUT
         << "*****************************************************************"
         << FEI_ENDL << "   Running tests with partitioned communicators/groups"
         << FEI_ENDL
         << "*****************************************************************"
         << FEI_ENDL << FEI_ENDL;
      }

      std::string path = fei_test_utils::get_arg_value("-d", argc, argv);

      //make sure every processor has the path string.
      fei_test_utils::broadcast_string(comm, 0, path);

      MPI_Comm newcomm1, newcomm2;
      MPI_Group newgroup1, newgroup2;

      //newcomm1 and newgroup1 will represent procs 0 and 1, while
      //newcomm2 and newgroup2 will represent procs 2 and 3.
      split_four_procs_into_two_groups(comm,
                                       newcomm1, newgroup1,
                                       newcomm2, newgroup2);

      if (localProc < 2) {
        execute_fullsystem_tests(newcomm1, path, name_numproc_pairs);
      }

    }
#endif
  }
}

int test_library_plugins(MPI_Comm comm)
{
  int errcode = 0;

  //--------- factory test --------------
  {
    fei::SharedPtr<tester> tst(new test_Factory(comm));

    try {
      errcode = tst->runtests();
    }
    catch(std::exception& exc) {
      FEI_COUT << "test_library_plugins: caught exception: "
          << exc.what() << FEI_ENDL;
      return(-1);
    }
  }

  //--------- vector test --------------
  {
    fei::SharedPtr<tester> tst(new test_Vector(comm));

    try {
      errcode = tst->runtests();
    }
    catch(std::exception& exc) {
      FEI_COUT << "test_library_plugins: caught exception: "
          << exc.what() << FEI_ENDL;
      return(-1);
    }
  }

  //--------- matrix test --------------
  {
    fei::SharedPtr<tester> tst(new test_Matrix(comm));

    try {
      errcode = tst->runtests();
    }
    catch(std::exception& exc) {
      FEI_COUT << "test_library_plugins: caught exception: "
          << exc.what() << FEI_ENDL;
      return(-1);
    }
  }

  return(errcode);
}

int execute_named_test(const std::string& testname,
		       int argc, char** argv, MPI_Comm comm)
{
  int numProcs = fei::numProcs(comm);
  int localProc = fei::localProc(comm);

  std::string path = fei_test_utils::get_arg_value("-d", argc, argv);

  //make sure every processor has the path string.
  fei_test_utils::broadcast_string(MPI_COMM_WORLD, 0, path);

  int errcode = 0;

  fei::Barrier(comm);

  bool testname_recognized = false;

  if (testname == "unit_tests") {
    testname_recognized = true;
    try {
      execute_unit_tests(path, comm);
    }
    catch(std::runtime_error& exc) {
      FEI_CERR << "caught unit-test error: "<<exc.what() << FEI_ENDL;
      return(-1);
    }
  }

  if (testname == "poisson_main") {
    errcode = poisson_main(argc, argv, comm, numProcs, localProc);
    testname_recognized = true;
  }
  if (testname == "poisson3_main") {
    errcode = poisson3_main(argc, argv, comm, numProcs, localProc);
    testname_recognized = true;
  }

  if (testname == "beam_oldfei_main") {
    errcode = beam_oldfei_main(argc, argv, comm, numProcs, localProc);
    testname_recognized = true;
  }
  if (testname == "beam_main") {
    errcode = beam_main(argc, argv, comm, numProcs, localProc);
    testname_recognized = true;
  }

  if (testname == "feiDriver_main") {
    errcode = feiDriver_main(argc, argv, comm, numProcs, localProc);
    testname_recognized = true;
  }

  if (testname == "cFeiTester_main") {
    errcode = cFeiTester_main(argc, argv, comm, numProcs, localProc);
    testname_recognized = true;
  }

  if (testname == "library_plugins") {
    errcode = test_library_plugins(comm);
    testname_recognized = true;
  }

  if (testname == "benchmarks") {
    testname_recognized = true;
    if (localProc == 0) {
      try {
        execute_benchmarks(comm);
      }
      catch(std::runtime_error& exc) {
        FEI_CERR <<"caught exception from benchmarks: "<<exc.what()<<FEI_ENDL;
        errcode = -1;
      }
    }
  }

  fei::Barrier(comm);

  if (!testname_recognized && localProc == 0) {
    FEI_COUT << "fei_test: '-test' argument used, but value ("<<testname
	     << ") not recognized. Valid values are:"<<FEI_ENDL
	     << "    unit_tests"<<FEI_ENDL
	     << "    benchmarks"<<FEI_ENDL
	     << "    library_plugins"<<FEI_ENDL
	     << "    poisson_main"<<FEI_ENDL
	     << "    poisson3_main"<<FEI_ENDL
	     << "    cube_main"<<FEI_ENDL
	     << "    cube3_main"<<FEI_ENDL
	     << "    feiDriver_main"<<FEI_ENDL
	     << "    cFeiTester_main"<<FEI_ENDL << FEI_ENDL;
    return(-1);
  }
  return(errcode);
}

void execute_benchmarks(MPI_Comm comm)
{
  test_benchmarks tst(comm);

  bool test_failed = false;
  if (tst.runtests() != 0) test_failed = false;

  if (test_failed) throw std::runtime_error("unit-test failed");
}

void execute_unit_tests(const std::string& path,
                        MPI_Comm comm)
{
  std::vector<fei::SharedPtr<tester> > testers;

  //Since each of the following classes implement the tester interface,
  //we can simply stick instances of them into an array, then run through
  //the array running the tests on each class instance.

  testers.push_back(fei::SharedPtr<tester>(new test_misc(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_Utils(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_Set(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_Database(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_EqnBuffer(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_EqnCommMgr(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_SSMat(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_PointBlockMap(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_Algebraic(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_AztecWrappers(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_Tables(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_VectorSpace(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_MatrixGraph(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_Vector(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_Matrix(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_Factory(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_SNL_FEI_Structure(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_FEI_Implementation(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_FEI_Impl(comm)));
  testers.push_back(fei::SharedPtr<tester>(new test_LinearSystem(comm)));

  std::vector<fei::SharedPtr<tester> >::const_iterator
    testers_iter = testers.begin(),
    testers_end  = testers.end();

  std::string failed_test_name;
  bool test_failed = false;
  for(; testers_iter != testers_end; ++testers_iter) {
    fei::Barrier(comm);

    fei::SharedPtr<tester> tst = *testers_iter;
    tst->setPath(path);
    if (tst->runtests() != 0) {
      failed_test_name = tst->getName();
      test_failed = true;
    }
  }

  if (test_failed) {
    std::string str1("unit-test failed: ");
    throw std::runtime_error(str1+failed_test_name);
  }
}

void execute_fullsystem_tests(MPI_Comm comm,
                              const std::string& path,
                              fei::ParameterSet& name_numproc_pairs)
{
  test_FEI test_fei(comm);

  if (!path.empty()) {
    test_fei.setPath(path.c_str());
  }
  else {
    test_fei.setPath(".");
  }

  //We'll iterate name_numproc_pairs, passing filenames to
  //test_fei and have it run tests.

  fei::ParameterSet::const_iterator
    n_iter = name_numproc_pairs.begin(),
    n_end = name_numproc_pairs.end();

  for(; n_iter != n_end; ++n_iter) {
    const char* fileName = (*n_iter).getName().c_str();
    int numProcsToUse = (*n_iter).getIntValue();

    if (numProcsToUse != fei::numProcs(comm)) {
      continue;
    }

    fei::Barrier(comm);
    if (fei::localProc(comm) == 0) {
      FEI_COUT << FEI_ENDL << "*****" << FEI_ENDL
           << fileName << FEI_ENDL << "*****"<<FEI_ENDL;
    }
    fei::Barrier(comm);

    test_fei.setFileName(fileName);

    int resultCode = test_fei.runtests();
    if (resultCode < 0) {
      throw std::runtime_error("nonzero resultCode from test_fei.runtests()");
    }
  }
}

#ifndef FEI_SER
int split_four_procs_into_two_groups(MPI_Comm comm,
				     MPI_Comm& newcomm1, MPI_Group& newgroup1,
				     MPI_Comm& newcomm2, MPI_Group& newgroup2)
{
  //This function is hardwired to operate on 4 processors. It will create two new
  //communicators and two new groups, each representing 2 processors. newcomm1
  //and newgroup1 will contain procs 0 and 1, while newcomm2 and newgroup2 will
  //contain procs 2 and 3.

  int numProcs, localProc;
  MPI_Comm_rank(comm, &localProc);
  MPI_Comm_size(comm, &numProcs);

  if (numProcs != 4) {
    return(-1);
  }

  std::vector<int> procs(numProcs);
  for(int i=0; i<numProcs; ++i) {
    procs[i] = i;
  }

  int midpoint = 2;
  int newgroup1_size = 2, newgroup2_size = 2;

  MPI_Group group;
  MPI_Comm_group(comm, &group);

  MPI_Group_incl(group, newgroup1_size, &procs[0], &newgroup1);
  MPI_Group_incl(group, newgroup2_size, &procs[0]+midpoint, &newgroup2);

  MPI_Comm_create(comm, newgroup1, &newcomm1);
  MPI_Comm_create(comm, newgroup2, &newcomm2);

  return(0);
}
#endif


