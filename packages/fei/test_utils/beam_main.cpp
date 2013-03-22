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

#include <fei_macros.hpp> //includes std stuff like iostream, etc.

//Including the header fei_base.hpp gets us the declaration for
//the various fei:: and snl_fei:: classes.

#include <fei_base.hpp>

//Now make provision for using any one of several solver libraries. This is
//handled by the code in LibraryFactory.{hpp,cpp}.

#include <test_utils/LibraryFactory.hpp>

//And finally, we need to include some 'utils' headers which are used for
//setting up the data for the example problem.

#include <test_utils/HexBeam.hpp>
#include <test_utils/HexBeamCR.hpp>
#include <snl_fei_Utils.hpp>
#include <test_utils/fei_test_utils.hpp>
//
//Include definitions of macros to call functions and check the return code.
//
#undef fei_file
#define fei_file "cube3.cpp"
#include <fei_ErrMacros.hpp>

//==============================================================================
//Here's the main...
//==============================================================================
int beam_main(int argc, char** argv,
	       MPI_Comm comm, int numProcs, int localProc){

  double start_time = fei::utils::cpu_time();

  std::vector<std::string> stdstrings;
  CHK_ERR( fei_test_utils::get_filename_and_read_input(argc, argv,
						comm, localProc,
						stdstrings) );
  fei::ParameterSet params;
  fei::utils::parse_strings(stdstrings, " ", params);

  std::string solverName;
  std::string datasource;
  std::string constraintform;
  int W = 0;
  int D = 0;
  int DofPerNode = 0;

  int errcode = 0;
  errcode += params.getStringParamValue("SOLVER_LIBRARY", solverName);
  errcode += params.getStringParamValue("DATA_SOURCE", datasource);
  errcode += params.getIntParamValue("W", W);
  errcode += params.getIntParamValue("D", D);
  errcode += params.getIntParamValue("DofPerNode", DofPerNode);

  if (errcode != 0) {
    fei::console_out() << "Failed to find one or more required parameters in input-file."
	     << FEI_ENDL << "Required parameters:"<<FEI_ENDL
	     << "SOLVER_LIBRARY" << FEI_ENDL
	     << "DATA_SOURCE" << FEI_ENDL
	     << "CONSTRAINT_FORM" << FEI_ENDL
	     << "W" << FEI_ENDL << "D" << FEI_ENDL << "DofPerNode" << FEI_ENDL;
    return(-1);
  }

  params.getStringParamValue("CONSTRAINT_FORM", constraintform);

  bool slave_constraints = false;
  if ("slaves" == constraintform) {
    slave_constraints = true;
  }

  HexBeam* hexcubeptr = NULL;
  if (datasource == "HexBeam") {
    hexcubeptr = new HexBeam(W, D, DofPerNode,
			     HexBeam::OneD, numProcs, localProc);
  }
  else {
    hexcubeptr = new HexBeamCR(W, D, DofPerNode,
			       HexBeam::OneD, numProcs, localProc);
  }

  HexBeam& hexcube = *hexcubeptr;

  if (localProc == 0) {
    int numCRs = (W+1)*(W+1)* ((numProcs*2)-1);
    if (hexcube.getNumCRs() < 1) numCRs = 0;
    FEI_COUT << FEI_ENDL;
    FEI_COUT << "========================================================" 
	 << FEI_ENDL;
    FEI_COUT << "FEI version: " << fei::utils::version() << FEI_ENDL;
    FEI_COUT << "--------------------------------------------------------"<<FEI_ENDL;
    FEI_COUT << "Size W: " << W << " (num-elements-along-side-of-cube)"<<FEI_ENDL;
    FEI_COUT << "Size D: " << D << " (num-elements-along-depth-of-cube)"<<FEI_ENDL;
    FEI_COUT << "DOF per node: " << DofPerNode <<FEI_ENDL;
    FEI_COUT << "Num local  elements: " << hexcube.localNumElems_ << FEI_ENDL;
    FEI_COUT << "Num global elements: " << hexcube.totalNumElems_ << FEI_ENDL;
    FEI_COUT << "Num local  DOF: " << hexcube.numLocalDOF_ << FEI_ENDL;
    FEI_COUT << "Num global DOF: " << hexcube.numGlobalDOF_ << FEI_ENDL;
    FEI_COUT << "Num global CRs: " << numCRs << FEI_ENDL;
    FEI_COUT << "========================================================" 
	 << FEI_ENDL;
  }

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
    FEI_COUT << "snl_fei::Factory allocation failed." << FEI_ENDL;
    return(-1);
  }

  factory->parameters(params);

  fei::SharedPtr<fei::VectorSpace> nodeSpace =
    factory->createVectorSpace(comm, "cube3");

  fei::SharedPtr<fei::VectorSpace> dummy;
  fei::SharedPtr<fei::MatrixGraph> matrixGraph =
    factory->createMatrixGraph(nodeSpace, dummy, "cube3");

  //load some solver-control parameters.
  nodeSpace->setParameters(params);
  matrixGraph->setParameters(params);

  int fieldID = 0;
  int fieldSize = hexcube.numDofPerNode();
  int nodeIDType = 0;
  int crIDType = 1;

  nodeSpace->defineFields( 1, &fieldID, &fieldSize );
  nodeSpace->defineIDTypes(1, &nodeIDType );
  nodeSpace->defineIDTypes(1, &crIDType );

  CHK_ERR( HexBeam_Functions::
	   init_elem_connectivities( matrixGraph.get(), hexcube ) );

  CHK_ERR( HexBeam_Functions::
	   init_shared_nodes( matrixGraph.get(), hexcube ) );

  int firstLocalCRID = 0;
  if (slave_constraints) {
    CHK_ERR( HexBeam_Functions::
	     init_slave_constraints( matrixGraph.get(), hexcube) );
  }
  else {
    CHK_ERR( HexBeam_Functions::
	     init_constraints( matrixGraph.get(), hexcube, localProc, firstLocalCRID) );
  }

  CHK_ERR( matrixGraph->initComplete() );

  double fei_init_time = fei::utils::cpu_time() - start_init_time;

  if (localProc == 0) {
    FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
    FEI_COUT << "Initialization cpu time:   " << fei_init_time << FEI_ENDL;
  }

  //Now the initialization phase is complete. Next we'll do the load phase,
  //which for this problem just consists of loading the element data
  //(element-wise stiffness arrays and load vectors) and the boundary
  //condition data.

  double fei_creatematrix_start_time = fei::utils::cpu_time();

  fei::SharedPtr<fei::Matrix> mat = factory->createMatrix(matrixGraph);

  double fei_creatematrix_time = fei::utils::cpu_time() - fei_creatematrix_start_time;
  if (localProc == 0) {
    FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
    FEI_COUT << "Create-Matrix cpu time:   " << fei_creatematrix_time << FEI_ENDL;
  }

  double start_load_time = fei::utils::cpu_time();

  fei::SharedPtr<fei::Vector> solnVec = factory->createVector(matrixGraph, true);
  fei::SharedPtr<fei::Vector> rhsVec  = factory->createVector(matrixGraph);
  fei::SharedPtr<fei::LinearSystem> linSys =
    factory->createLinearSystem(matrixGraph);

  CHK_ERR( linSys->parameters(params));

  linSys->setMatrix(mat);
  linSys->setSolutionVector(solnVec);
  linSys->setRHS(rhsVec);

  CHK_ERR( HexBeam_Functions::
	   load_elem_data(matrixGraph.get(), mat.get(), rhsVec.get(), hexcube) );

  //temporarily disable boundary-conditions
  //  CHK_ERR( load_BC_data(linSys, hexcube) );

  if (!slave_constraints) {
    CHK_ERR( HexBeam_Functions::
	     load_constraints(linSys.get(), hexcube, firstLocalCRID) );
  }

  CHK_ERR( linSys->loadComplete() );

  double fei_load_time = fei::utils::cpu_time() - start_load_time;

  if (localProc == 0) {
    //IOS macros are defined in fei_macros.h
    FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
    FEI_COUT << "Coef. loading cpu time:    " << fei_load_time << FEI_ENDL;
    FEI_COUT << "Total assembly wall time:   "
	 << fei_init_time + fei_creatematrix_time + fei_load_time << FEI_ENDL;
  }

  //
  //now the load phase is complete, so we're ready to launch the underlying
  //solver and solve Ax=b
  //

  //First, check whether the params (set from the input .i file) specify
  //a "Trilinos_Solver" string (which may be Amesos...). If not, it won't
  //matter but if so, the factory may switch based on this value, when
  //creating the solver object instance.

  std::string solver_name_value;
  params.getStringParamValue("Trilinos_Solver", solver_name_value);

  const char* charptr_solvername =
    solver_name_value.empty() ? 0 : solver_name_value.c_str();

  fei::SharedPtr<fei::Solver> solver = factory->createSolver(charptr_solvername);

  int status;
  int itersTaken = 0;

  if (localProc==0) FEI_COUT << "solve..." << FEI_ENDL;
  double start_solve_time = fei::utils::cpu_time();

  int err = solver->solve(linSys.get(),
			  NULL, //preconditioningMatrix
			  params,
			  itersTaken,
			  status);

  double solve_time = fei::utils::cpu_time()-start_solve_time;

  if (err!=0) {
    if (localProc==0) FEI_COUT << "solve returned err: " << err <<", status: "
			   << status << FEI_ENDL;
  }

  if (localProc==0) {
    FEI_COUT << " cpu-time in solve: " << solve_time << FEI_ENDL;
  }

  CHK_ERR( solnVec->scatterToOverlap() );

  delete hexcubeptr;

  double elapsed_cpu_time = fei::utils::cpu_time() - start_time;
  int returnValue = 0;

  //The following IOS_... macros are defined in base/fei_macros.h
  FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
  if (localProc==0) {
    FEI_COUT << "Proc0 cpu times (seconds):" << FEI_ENDL
	 << "   FEI initialize:    " << fei_init_time << FEI_ENDL
         << "   FEI create-matrix: " << fei_creatematrix_time << FEI_ENDL
         << "   FEI load:          " << fei_load_time << FEI_ENDL
         << "      solve:          " << solve_time << FEI_ENDL
         << "Total program time:   " << elapsed_cpu_time << FEI_ENDL;

#if defined(FEI_PLATFORM) && defined(FEI_OPT_LEVEL)
    double benchmark = fei_init_time;

    std::string slavestr;
    if (!constraintform.empty()) {
      slavestr = constraintform;
    }
    if (slavestr.size() > 0) slavestr += "_";

    FEI_OSTRINGSTREAM testname_init;
    testname_init << "cube3_init_"<<slavestr<<W<<"_"<<D<<"_"<<DofPerNode<<"_"
		  <<solverName<<"_np"<<numProcs<<"_"
		  <<FEI_PLATFORM<<"_"<<FEI_OPT_LEVEL;
    FEI_OSTRINGSTREAM testname_create;
    testname_create << "cube3_creatematrix_"<<slavestr<<W<<"_"<<D<<"_"<<DofPerNode
		    <<"_"<<solverName<<"_np"<<numProcs<<"_"
		    <<FEI_PLATFORM<<"_"<<FEI_OPT_LEVEL;
    FEI_OSTRINGSTREAM testname_load;
    testname_load << "cube3_load_"<<slavestr<<W<<"_"<<D<<"_"<<DofPerNode<<"_"
		  <<solverName<<"_np"<<numProcs<<"_"
		  <<FEI_PLATFORM<<"_"<<FEI_OPT_LEVEL;

    double file_init, file_create, file_load;
    bool file_benchmarks_available = true;
    try {
      file_init = fei_test_utils::get_file_benchmark("./cube3_timings.txt",
					      testname_init.str().c_str());
      file_create = fei_test_utils::get_file_benchmark("./cube3_timings.txt",
						testname_create.str().c_str());
      file_load = fei_test_utils::get_file_benchmark("./cube3_timings.txt",
					      testname_load.str().c_str());
    }
    catch (std::runtime_error& exc) {
      file_benchmarks_available = false;
    }

    if (file_benchmarks_available) {

      bool init_test_passed =
	fei_test_utils::check_and_cout_test_result(testname_init.str(), benchmark,
					    file_init, 10);

      benchmark = fei_creatematrix_time;
      bool create_test_passed =
	fei_test_utils::check_and_cout_test_result(testname_create.str(), benchmark,
					    file_create, 10);

      benchmark = fei_load_time;
      bool load_test_passed =
	fei_test_utils::check_and_cout_test_result(testname_load.str(), benchmark,
					    file_load, 10);

      returnValue = init_test_passed&&create_test_passed&&load_test_passed ? 0 : 1;
    }
#endif

  }

  bool testPassed = returnValue==0;
  if (testPassed && localProc == 0) {
    //A string for the SIERRA runtest tool to search for in test output...
    FEI_COUT << "FEI test successful" << FEI_ENDL;
  }

  return(returnValue);
}
