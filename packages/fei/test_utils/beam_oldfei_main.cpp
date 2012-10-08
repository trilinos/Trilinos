/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
// This is a simple program to exercise the new FEI (version 3?) classes,
// for the purposes of testing, code tuning and scaling studies.
//

#include <fei_macros.hpp> //simply includes std stuff like iostream, etc.
                             //users could include that stuff themselves rather
                             //than using this header.

//Including the header fei_base.hpp gets us the declaration for
//various fei classes and interfaces.

#include <fei_base.hpp>
#include <FEI_Implementation.hpp>

//Now make provision for using any one of several solver libraries. This is
//handled by the code in LibraryFactory.[hC].

#include <test_utils/LibraryFactory.hpp>

//And finally, we need to include some 'utils' headers which are used for
//setting up the data for the example problem.

#include <test_utils/HexBeam.hpp>
#include <test_utils/HexBeamCR.hpp>

#include <test_utils/fei_test_utils.hpp>
#include <snl_fei_Utils.hpp>

#undef fei_file
#define fei_file "cube_main.cpp"
#include <fei_ErrMacros.hpp>

//==============================================================================
//Here's the main...
//==============================================================================
int beam_oldfei_main(int argc, char** argv,
	      MPI_Comm comm, int numProcs, int localProc){

  std::vector<std::string> stdstrings;
  CHK_ERR( fei_test_utils::get_filename_and_read_input(argc, argv,
						comm, localProc,
						stdstrings) );

  const char** params = NULL;
  int numParams = 0;
  fei::utils::strings_to_char_ptrs(stdstrings, numParams, params);
  
  fei::ParameterSet paramset;
  fei::utils::parse_strings(stdstrings, " ", paramset);

  std::string whichFEI;
  std::string solverName;
  std::string datasource;
  int W = 0;
  int D = 0;
  int DofPerNode = 0;

  int errcode = 0;
  errcode += paramset.getStringParamValue("SOLVER_LIBRARY", solverName);
  errcode += paramset.getStringParamValue("DATA_SOURCE", datasource);
  errcode += paramset.getStringParamValue("WHICH_FEI", whichFEI);
  errcode += paramset.getIntParamValue("W", W);
  errcode += paramset.getIntParamValue("D", D);
  errcode += paramset.getIntParamValue("DofPerNode", DofPerNode);

  if (errcode != 0) {
    fei::console_out() << "Failed to find one or more required parameters in input-file."
	     << FEI_ENDL << "Required parameters:"<<FEI_ENDL
	     << "SOLVER_LIBRARY" << FEI_ENDL
	     << "DATA_SOURCE" << FEI_ENDL
	     << "CONSTRAINT_FORM" << FEI_ENDL
	     << "W" << FEI_ENDL << "D" << FEI_ENDL << "DofPerNode" << FEI_ENDL;
    return(-1);
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
    int numCRs = (W+1)*(W+1)*(numProcs*2)-1;
    if (hexcube.getNumCRs() < 1) numCRs = 0;
    FEI_COUT << FEI_ENDL;
    FEI_COUT << "========================================================" 
	 << FEI_ENDL;
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

  //CHK_ERR( print_cube_data(hexcube, numProcs, localProc) );

  fei::SharedPtr<fei::Factory> factory;
  fei::SharedPtr<LibraryWrapper> wrapper;
  fei::SharedPtr<FEI> fei;

  if (whichFEI == "OLDFEI") {
    try {
      wrapper = fei::create_LibraryWrapper(comm, solverName.c_str());
    }
    catch (std::runtime_error& exc) {
      fei::console_out() << exc.what() << FEI_ENDL;
      ERReturn(-1);
    }
    fei.reset(new FEI_Implementation(wrapper, comm));
  }
  else if (whichFEI == "fei::FEI_Impl") {
    try {
      factory = fei::create_fei_Factory(comm, solverName.c_str());
    }
    catch (std::runtime_error& exc) {
      fei::console_out() << exc.what() << FEI_ENDL;
      ERReturn(-1);
    }
    fei = factory->createFEI(comm);
  }
  else {
    fei::console_out() << "cube ERROR, value of 'WHICH_FEI' must be 'OLDFEI' or 'fei::FEI_Impl'"<< FEI_ENDL;
    ERReturn(-1);
  }

  //load some control parameters.
  CHK_ERR( fei->parameters(numParams, params) );

  delete [] params;

  CHK_ERR( fei->setSolveType(FEI_SINGLE_SYSTEM) );

  int fieldID = 0;
  int fieldSize = hexcube.numDofPerNode();

  CHK_ERR( fei->initFields( 1, &fieldSize, &fieldID ) );

  CHK_ERR( HexBeam_Functions::init_elem_connectivities( fei.get(), hexcube ) );

  CHK_ERR( HexBeam_Functions::init_shared_nodes( fei.get(), hexcube ) );

  int firstLocalCRID;
  CHK_ERR( HexBeam_Functions::
	   init_constraints( fei.get(), hexcube, firstLocalCRID) );

  CHK_ERR( fei->initComplete() );

  double fei_init_time = fei::utils::cpu_time() - start_init_time;

  if (localProc == 0) {
    FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
    FEI_COUT << "Initialization time:   " << fei_init_time << FEI_ENDL;
  }

  //Now the initialization phase is complete. Next we'll do the load phase,
  //which for this problem just consists of loading the element data
  //(element-wise stiffness arrays and load vectors) and the boundary
  //condition data.

  double start_load_time = fei::utils::cpu_time();

  CHK_ERR( HexBeam_Functions::load_elem_data(fei.get(), hexcube) );

  CHK_ERR( HexBeam_Functions::load_constraints(fei.get(), hexcube, firstLocalCRID) );

  //FEI_COUT << "calling load_BC_data"<<FEI_ENDL;
  //CHK_ERR( load_BC_data(fei, hexcube) );

  fei->loadComplete();

  double fei_load_time = fei::utils::cpu_time() - start_load_time;

  delete hexcubeptr;

  if (localProc == 0) {
    //IOS macros are defined in fei_macros.h
    FEI_COUT.setf(IOS_FIXED, IOS_FLOATFIELD);
    FEI_COUT << "Coef. loading  time:    " << fei_load_time << FEI_ENDL;
    FEI_COUT << "Total assembly time:    " << fei_init_time + fei_load_time << FEI_ENDL;
  }

  //
  //now the load phase is complete, so we're ready to launch the underlying
  //solver and solve Ax=b
  //

  int err;
  int status;
  double start_solve_time = fei::utils::cpu_time();

  err = fei->solve(status);

  if (err || status) {
    if (localProc==0) FEI_COUT << "solve returned status: " << status << FEI_ENDL;
  }

  double solve_time = fei::utils::cpu_time() - start_solve_time;

  if (localProc == 0) {
    FEI_COUT << "Solver time:          " << solve_time << FEI_ENDL;
  }

  int returnValue = 0;
  if (localProc == 0) {
#if defined(FEI_PLATFORM) && defined(FEI_OPT_LEVEL)
    double benchmark = fei_init_time;

    FEI_OSTRINGSTREAM testname_init;
    testname_init << "cube_"<<whichFEI<<"_init_"<<W<<"_"<<D<<"_"<<DofPerNode<<"_"
		  <<solverName<<"_np"<<numProcs<<"_"
		  <<FEI_PLATFORM<<"_"<<FEI_OPT_LEVEL;

    FEI_OSTRINGSTREAM testname_load;
    testname_load << "cube_"<<whichFEI<<"_load_"<<W<<"_"<<D<<"_"<<DofPerNode<<"_"
		  <<solverName<<"_np"<<numProcs<<"_"
		  <<FEI_PLATFORM<<"_"<<FEI_OPT_LEVEL;

    double file_init, file_load;
    bool file_benchmarks_available = true;
    try {
      file_init = fei_test_utils::get_file_benchmark("./cube_timings.txt",
					      testname_init.str().c_str());
      file_load = fei_test_utils::get_file_benchmark("./cube_timings.txt",
					      testname_load.str().c_str());
    }
    catch (std::runtime_error& exc) {
      file_benchmarks_available = false;
    }

    if (file_benchmarks_available) {

      bool init_test_passed =
	fei_test_utils::check_and_cout_test_result(testname_init.str(), benchmark,
					    file_init, 10);

      benchmark = fei_load_time;
      bool load_test_passed =
	fei_test_utils::check_and_cout_test_result(testname_load.str(), benchmark,
					    file_load, 10);

      returnValue = init_test_passed&&load_test_passed ? 0 : 1;
    }
#endif
  }

  bool testPassed = returnValue==0;
  if (testPassed && localProc == 0) {
    //A string for the SIERRA runtest tool to search for in test output...
    FEI_COUT << "cube execution successful" << FEI_ENDL;
  }

  fei.reset();

  return(0);
}
