/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
// This is a simple program to exercise the FEI, for the purposes of
// testing, code tuning and scaling studies.
//
// This program assembles a linear system from a 2D square Poisson
// problem, using 4-node square elements. There is only 1 degree-of-
// freedom per node.
//
// This problem was coded up with the help of Ray Tuminaro.
// Alan Williams 12-20-2000
//
// The input file for this program should provide the following:
//   WHICH_FEI <OLDFEI|fei::FEI_Impl>
//   SOLVER_LIBRARY <library-name> -- e.g., Aztec
//   L <int> -- the global length (num-elements) on a side of the 2D square 
//
//

#include <fei_iostream.hpp>
#include <cmath>
#include <fei_base.hpp>
#include <FEI_Implementation.hpp>

//Now make provision for using any one of several solver libraries. This is
//handled by the code in LibraryFactory.{hpp,cpp}.

#include <test_utils/LibraryFactory.hpp>
#include <fei_LibraryWrapper.hpp>

//And, we need to include some headers for utility classes which are simply
//used for setting up the data for the example problem.

#include <test_utils/Poisson_Elem.hpp>
#include <test_utils/PoissonData.hpp>

#include <test_utils/ElemBlock.hpp>
#include <test_utils/CRSet.hpp>
#include <test_utils/CommNodeSet.hpp>
#include <test_utils/DataReader.hpp>

#include <test_utils/SolnCheck.hpp>

#include <test_utils/fei_test_utils.hpp>
#include <snl_fei_Utils.hpp>
//
//Include definitions of macros to call functions and check the return code.
//
#undef fei_file
#define fei_file "poisson_main.cpp"
#include <fei_ErrMacros.hpp>

//============================================================================
//Here's the main...
//============================================================================
int poisson_main(int argc, char** argv,
                 MPI_Comm comm, int numProcs, int localProc){
  int outputLevel;
  int masterProc = 0, err = 0;

  double start_time = fei::utils::cpu_time();

  std::vector<std::string> stdstrings;
  CHK_ERR( fei_test_utils::get_filename_and_read_input(argc, argv,
						comm, localProc,
						stdstrings) );
  const char** params = NULL;
  int numParams = 0;
  fei::utils::strings_to_char_ptrs(stdstrings, numParams, params);
  
  fei::ParameterSet paramset;
  fei::utils::parse_strings(stdstrings, " ", paramset);

  std::string which_fei;
  std::string solverName;
  int L = 0;

  int errcode = 0;
  errcode += paramset.getStringParamValue("SOLVER_LIBRARY", solverName);
  errcode += paramset.getStringParamValue("WHICH_FEI", which_fei);
  errcode += paramset.getIntParamValue("L", L);

  if (errcode != 0) {
    fei::console_out() << "Failed to find one or more required parameters in input-file."
	     << FEI_ENDL << "Required parameters:"<<FEI_ENDL
	     << "SOLVER_LIBRARY" << FEI_ENDL
	     << "WHICH_FEI" << FEI_ENDL
	     << "L" << FEI_ENDL;
    return(-1);
  }

  if (localProc == 0) {
    int nodes = (L+1)*(L+1);
    int eqns = nodes;
    FEI_COUT << FEI_ENDL;
    FEI_COUT << "========================================================" 
	 << FEI_ENDL;
    FEI_COUT << "Square size     L: " << L << " elements." << FEI_ENDL;
    FEI_COUT << "Global number of elements: " << L*L << FEI_ENDL;
    FEI_COUT << "Global number of nodes: " << nodes << FEI_ENDL;
    FEI_COUT << "Global number of equations: " << eqns <<FEI_ENDL;
    FEI_COUT << "========================================================" 
	 << FEI_ENDL;
  }

  outputLevel = fei_test_utils::whichArg(numParams, params, "outputLevel 1");
  if (outputLevel >= 0) outputLevel = 1;
  if (outputLevel < 0) outputLevel = 0;

  if ((masterProc == localProc)&&(outputLevel>0)) {
    fei_test_utils::print_args(argc, argv);
  }

  if (outputLevel == 1) {
    if (localProc != 0) outputLevel = 0;
  }

  //PoissonData is the object that will be in charge of generating the
  //data to pump into the FEI object.

  PoissonData poissonData(L, numProcs, localProc, outputLevel);

  double start_init_time = fei::utils::cpu_time();

  fei::SharedPtr<fei::Factory> factory;
  fei::SharedPtr<LibraryWrapper> wrapper;
  fei::SharedPtr<FEI> fei;

  if (which_fei == "OLDFEI") {
    try {
      wrapper = fei::create_LibraryWrapper(comm, solverName.c_str());
    }
    catch (std::runtime_error& exc) {
      fei::console_out() << exc.what()<<FEI_ENDL;
      ERReturn(-1);
    }
    fei.reset(new FEI_Implementation(wrapper, comm));
  }
  else if (which_fei == "fei::FEI_Impl") {
    try {
      factory = fei::create_fei_Factory(comm, solverName.c_str());
    }
    catch (std::runtime_error& exc) {
      fei::console_out() << exc.what()<<FEI_ENDL;
      ERReturn(-1);
    }
    fei = factory->createFEI(comm);
  }
  else {
    fei::console_out() << "poisson_main ERROR, value of 'WHICH_FEI' must be 'OLDFEI' or 'fei::FEI_Impl'"<< FEI_ENDL;
    ERReturn(-1);
  }

  const char* feiVersionString;
  CHK_ERR( fei->version(feiVersionString) );

  if (localProc==0) FEI_COUT << feiVersionString << FEI_ENDL;

  //load some parameters.
  CHK_ERR( fei->parameters( numParams, params ) );

  delete [] params;

  if (outputLevel>0 && localProc==0) FEI_COUT << "setSolveType" << FEI_ENDL;
  CHK_ERR( fei->setSolveType(FEI_SINGLE_SYSTEM) );


  int numFields = poissonData.getNumFields();
  int* fieldSizes = poissonData.getFieldSizes();
  int* fieldIDs = poissonData.getFieldIDs();

  if (outputLevel>0 && localProc==0) FEI_COUT << "initFields" << FEI_ENDL;
  CHK_ERR( fei->initFields( numFields, fieldSizes, fieldIDs ) );


  CHK_ERR( init_elem_connectivities(fei.get(), poissonData) );

  CHK_ERR( set_shared_nodes(fei.get(), poissonData) );

  //The FEI_COUT and IOS_... macros are defined in base/fei_iostream.hpp
  FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);

  if (outputLevel>0 && localProc==0) FEI_COUT << "initComplete" << FEI_ENDL;
  CHK_ERR( fei->initComplete() );

  double fei_init_time = fei::utils::cpu_time() - start_init_time;

  //Now the initialization phase is complete. Next we'll do the load phase,
  //which for this problem just consists of loading the element data
  //(element-wise stiffness arrays and load vectors) and the boundary
  //condition data.
  //This simple problem doesn't have an constraint relations, etc.

  double start_load_time = fei::utils::cpu_time();

  CHK_ERR( load_elem_data(fei.get(), poissonData) );

  CHK_ERR( load_BC_data(fei.get(), poissonData) );

  double fei_load_time = fei::utils::cpu_time() - start_load_time;

  //
  //now the load phase is complete, so we're ready to launch the underlying
  //solver and solve Ax=b
  //
  int status;
  if (outputLevel>0 && localProc==0) FEI_COUT << "solve..." << FEI_ENDL;
  double start_solve_time = fei::utils::cpu_time();
  err = fei->solve(status);
  if (err) {
    if (localProc==0) FEI_COUT << "solve returned err: " << err << FEI_ENDL;
  }

  double iTime, lTime, sTime, rTime;
  CHK_ERR(fei->cumulative_cpu_times(iTime, lTime, sTime, rTime) );

  double solve_time = fei::utils::cpu_time() - start_solve_time;

  if (localProc == 0) {
    FEI_COUT << "FEI cpu-times:" << FEI_ENDL
	 << "    init. phase: " << iTime << FEI_ENDL
	 << "    load  phase: " << lTime << FEI_ENDL
	 << "    solve  time: " << sTime << FEI_ENDL;
  }

  double norm = 0.0;
  FEI_COUT.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);
  CHK_ERR( fei->residualNorm(1, 1, &fieldIDs[0], &norm) );
  if (localProc == 0) {
    FEI_COUT << "returned residual norm: " << norm << FEI_ENDL;
  }

  int itersTaken = 0;
  CHK_ERR( fei->iterations(itersTaken) );

  //
  //We oughtta make sure the solution we just computed is correct...
  //

  int numNodes = 0;
  CHK_ERR( fei->getNumLocalNodes(numNodes) );

  double maxErr = 0.0;
  if (numNodes > 0) {
    int lenNodeIDs = numNodes;
    GlobalID* nodeIDs = new GlobalID[lenNodeIDs];
    double* soln = new double[lenNodeIDs];
    if (nodeIDs != NULL && soln != NULL) {
      CHK_ERR( fei->getLocalNodeIDList(numNodes, nodeIDs, lenNodeIDs) );

      int fieldID = 1;
      CHK_ERR( fei->getNodalFieldSolution(fieldID, numNodes, nodeIDs, soln));

      for(int i=0; i<numNodes; i++) {
	int nID = (int)nodeIDs[i];
	double x = (1.0* ((nID-1)%(L+1)))/L;
	double y = (1.0* ((nID-1)/(L+1)))/L;

	double exactSoln = x*x + y*y;
	double error = std::abs(exactSoln - soln[i]);
	if (maxErr < error) maxErr = error;
      }

      delete [] nodeIDs;
      delete [] soln;
    }
    else {
      fei::console_out() << "allocation of nodeIDs or soln failed." << FEI_ENDL; 
    }

  }

#ifndef FEI_SER
  double globalMaxErr = 0.0;
  MPI_Allreduce(&maxErr, &globalMaxErr, 1, MPI_DOUBLE, MPI_MAX, comm);
  maxErr = globalMaxErr;
#endif
  bool testPassed = true;
  if (maxErr > 1.e-6) testPassed = false;

  if (testPassed && localProc == 0) {
    FEI_COUT << "poisson: TEST PASSED, maxErr = " << maxErr << ", iterations: "
	 << itersTaken << FEI_ENDL;
    //This is something the SIERRA runtest tool looks for in test output...
    FEI_COUT << "SIERRA execution successful" << FEI_ENDL;
  }
  if (testPassed == false && localProc == 0) {
    FEI_COUT << "maxErr = " << maxErr << ", TEST FAILED" << FEI_ENDL;
    FEI_COUT << "(Test is deemed to have passed if the maximum difference"
	 << " between the exact and computed solutions is 1.e-6 or less.)"
	 << FEI_ENDL;
  }

  double elapsed_cpu_time = fei::utils::cpu_time() - start_time;

  //The following IOS_... macros are defined in base/fei_macros.h
  FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
  if (localProc==0) {
    FEI_COUT << "Proc0 cpu times (seconds):" << FEI_ENDL
	 << "   FEI initialize:  " << fei_init_time << FEI_ENDL
         << "   FEI load:        " << fei_load_time << FEI_ENDL
         << "      solve:        " << solve_time << FEI_ENDL
         << "Total program time: " << elapsed_cpu_time << FEI_ENDL;
  }

  wrapper.reset();
  fei.reset();
  factory.reset();
  //If Prometheus is being used, we need to make sure that the
  //LibraryWrapper inside factory is deleted before MPI_Finalize() is called.
  //(That's why we call the 'reset' methods on these shared-pointers rather
  //than just letting them destroy themselves when they go out of scope.)

  return(0);
}

