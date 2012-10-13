/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
// This is a simple program to exercise FEI classes,
// for the purposes of testing, code tuning and scaling studies.
//
// This program assembles a linear system from a 2D square Poisson
// problem, using 4-node square elements. There is only 1 degree-of-
// freedom per node, and only one element-block per processor.
//
// This problem was coded up with Ray Tuminaro.
//
// The input file for this program should provide the following:
//   SOLVER_LIBRARY <library-name> -- e.g., Aztec
//   L <int> -- the global length (num-elements) on a side of the 2D square 
//
// Alan Williams 03-13-2002
//

#include <fei_iostream.hpp>
#include <cmath>

//Including the header fei_base.hpp gets us the declaration for
//the various FEI classes.

#include <fei_base.hpp>

//Now make provision for using any one of several solver libraries. This is
//handled by the code in LibraryFactory.[hC].

#include <test_utils/LibraryFactory.hpp>

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
#include <fei_ErrMacros.hpp>


//==============================================================================
//Here's the main...
//==============================================================================
int poisson3_main(int argc, char** argv,
                  MPI_Comm comm, int numProcs, int localProc){
  int masterProc = 0;

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

  std::string solverName;
  int L = 0;
  int outputLevel = 0;

  int errcode = 0;
  errcode += paramset.getStringParamValue("SOLVER_LIBRARY", solverName);
  errcode += paramset.getIntParamValue("L", L);
  paramset.getIntParamValue("outputLevel", outputLevel);

  if (errcode != 0) {
    fei::console_out() << "Failed to find one or more required parameters in input-file."
	     << FEI_ENDL << "Required parameters:"<<FEI_ENDL
	     << "SOLVER_LIBRARY" << FEI_ENDL
	     << "L" << FEI_ENDL;
    return(-1);
  }

  if (localProc == 0) {
    int nodes = (L+1)*(L+1);
    int eqns = nodes;
    FEI_COUT << FEI_ENDL;
    FEI_COUT << "========================================================" 
	 << FEI_ENDL;
    FEI_COUT << "FEI version: " << fei::utils::version() << FEI_ENDL;
    FEI_COUT << "Square size     L: " << L << " elements." << FEI_ENDL;
    FEI_COUT << "Global number of elements: " << L*L << FEI_ENDL;
    FEI_COUT << "Global number of nodes: " << nodes << FEI_ENDL;
    FEI_COUT << "Global number of equations: " << eqns <<FEI_ENDL;
    FEI_COUT << "========================================================" 
	 << FEI_ENDL;
  }

  if ((masterProc == localProc)&&(outputLevel>0)) {
    fei_test_utils::print_args(argc, argv);
  }

  if (outputLevel == 1) {
    if (localProc != 0) outputLevel = 0;
  }

  //PoissonData is the object that will be in charge of generating the
  //data to pump into the FEI objects.

  PoissonData poissonData(L, numProcs, localProc, outputLevel);

  double start_init_time = fei::utils::cpu_time();

  fei::SharedPtr<fei::Factory> factory;
  try {
    factory = fei::create_fei_Factory(comm, solverName.c_str());
  }
  catch (...) {
    FEI_COUT << "library " << solverName << " not available."<<FEI_ENDL;
    return(-1);
  }

  if (factory.get() == NULL) {
    FEI_COUT << "snl_fei::Factory creation failed." << FEI_ENDL;
    return(-1);
  }

  factory->parameters(paramset);

  fei::SharedPtr<fei::VectorSpace> nodeSpace =
    factory->createVectorSpace(comm, "poisson3");

  fei::SharedPtr<fei::VectorSpace> dummy;
  fei::SharedPtr<fei::MatrixGraph> matrixGraph =
    factory->createMatrixGraph(nodeSpace, dummy, "poisson3");

  //load some control parameters.
  matrixGraph->setParameters(paramset);


  int numFields = poissonData.getNumFields();
  int* fieldSizes = poissonData.getFieldSizes();
  int* fieldIDs = poissonData.getFieldIDs();
  int nodeIDType = 0;

  if (outputLevel>0 && localProc==0) FEI_COUT << "defineFields" << FEI_ENDL;
  nodeSpace->defineFields( numFields, fieldIDs, fieldSizes );

  if (outputLevel>0 && localProc==0) FEI_COUT << "defineIDTypes" << FEI_ENDL;
  nodeSpace->defineIDTypes( 1, &nodeIDType );

  CHK_ERR( init_elem_connectivities(matrixGraph.get(), poissonData) );

  CHK_ERR( set_shared_nodes(nodeSpace.get(), poissonData) );


  //The following IOS_... macros are defined in base/fei_macros.h
  FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
  if (outputLevel>0 && localProc==0) FEI_COUT << "initComplete" << FEI_ENDL;
  CHK_ERR( matrixGraph->initComplete() );

  double fei_init_time = fei::utils::cpu_time() - start_init_time;

  //Now the initialization phase is complete. Next we'll do the load phase,
  //which for this problem just consists of loading the element data
  //(element-wise stiffness arrays and load vectors) and the boundary
  //condition data.
  //This simple problem doesn't have any constraint relations, etc.

  double start_load_time = fei::utils::cpu_time();


  fei::SharedPtr<fei::Matrix> mat = factory->createMatrix(matrixGraph);
  fei::SharedPtr<fei::Vector> solnVec = factory->createVector(nodeSpace, true);
  fei::SharedPtr<fei::Vector> rhsVec  = factory->createVector(nodeSpace);
  fei::SharedPtr<fei::LinearSystem> linSys= factory->createLinearSystem(matrixGraph);

  linSys->setMatrix(mat);
  linSys->setSolutionVector(solnVec);
  linSys->setRHS(rhsVec);

  CHK_ERR( linSys->parameters(paramset));
  CHK_ERR( load_elem_data(matrixGraph.get(), mat.get(), rhsVec.get(), poissonData) );

  CHK_ERR( load_BC_data(linSys.get(), poissonData) );

  CHK_ERR( linSys->loadComplete() );

  double fei_load_time = fei::utils::cpu_time() - start_load_time;

  //
  //now the load phase is complete, so we're ready to launch the underlying
  //solver and solve Ax=b
  //

  fei::SharedPtr<fei::Solver> solver = factory->createSolver();

  int status;
  int itersTaken = 0;

  if (outputLevel>0 && localProc==0) FEI_COUT << "solve..." << FEI_ENDL;
  double start_solve_time = fei::utils::cpu_time();

  int err = solver->solve(linSys.get(),
			  NULL, //preconditioningMatrix
			  paramset, itersTaken, status);

  double solve_time = fei::utils::cpu_time() - start_solve_time;

  if (err!=0) {
    if (localProc==0) FEI_COUT << "solve returned err: " << err <<", status: "
			   << status << FEI_ENDL;
  }

  CHK_ERR( solnVec->scatterToOverlap() );

  //
  //We oughtta make sure the solution we just computed is correct...
  //

  int numNodes = nodeSpace->getNumOwnedAndSharedIDs(nodeIDType);

  double maxErr = 0.0;
  if (numNodes > 0) {
    int lenNodeIDs = numNodes;
    GlobalID* nodeIDs = new GlobalID[lenNodeIDs];
    double* soln = new double[lenNodeIDs];
    if (nodeIDs != NULL && soln != NULL) {
      CHK_ERR( nodeSpace->getOwnedAndSharedIDs(nodeIDType, numNodes,
				      nodeIDs, lenNodeIDs) );

      int fieldID = 1;
      CHK_ERR( solnVec->copyOutFieldData(fieldID, nodeIDType,
				      numNodes, nodeIDs, soln));

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

  double elapsed_cpu_time = fei::utils::cpu_time() - start_time;
  int returnValue = 0;

  //The following IOS_... macros are defined in base/fei_macros.h
  FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
  if (localProc==0) {
    FEI_COUT << "Proc0 cpu times (seconds):" << FEI_ENDL
	 << "   FEI initialize:  " << fei_init_time << FEI_ENDL
         << "   FEI load:        " << fei_load_time << FEI_ENDL
         << "      solve:        " << solve_time << FEI_ENDL
         << "Total program time: " << elapsed_cpu_time << FEI_ENDL;

#if defined(FEI_PLATFORM) && defined(FEI_OPT_LEVEL)
    double benchmark = fei_init_time+fei_load_time;
    FEI_OSTRINGSTREAM testname;
    testname << "poisson3_"<<L<<"_"<<solverName<<"_np"<<numProcs<<"_"
	     <<FEI_PLATFORM<<"_"<<FEI_OPT_LEVEL;
 
    double file_value = 0.0;
    bool file_value_available = true;
    try {
      file_value = fei_test_utils::get_file_benchmark("./poisson3_timings.txt",
					       testname.str().c_str());
    }
    catch(std::runtime_error& exc) {
      file_value_available = false;
    }

    if (file_value_available) {
      bool test_passed = fei_test_utils::check_and_cout_test_result(testname.str(),
							     benchmark,
							     file_value, 10);
      returnValue = test_passed ? 0 : 1;
    }
#endif

  }

  if (testPassed && returnValue==0 && localProc == 0) {
    FEI_COUT.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);
    FEI_COUT << "poisson: TEST PASSED, maxErr = " << maxErr << ", iterations: "
	 << itersTaken << FEI_ENDL;
    //This is something the SIERRA runtest tool looks for in test output...
    FEI_COUT << "SIERRA execution successful" << FEI_ENDL;
  }
  if ((testPassed == false || returnValue != 0) && localProc == 0) {
    if (testPassed==true) {
      FEI_COUT << "maxErr = " << maxErr << ", but time-taken outside margin. TEST FAILED" << FEI_ENDL;
    }
    else {
      FEI_COUT << "maxErr = " << maxErr << ", TEST FAILED"<<FEI_ENDL;
    }
    FEI_COUT << "(Test is deemed to have passed if the maximum difference"
	 << " between the exact and computed solutions is 1.e-6 or less, *AND*"
	 << " time-taken matches file-benchmark if available.)"
	 << FEI_ENDL;
  }

  return(returnValue);
}

