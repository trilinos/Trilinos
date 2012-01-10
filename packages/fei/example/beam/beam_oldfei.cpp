/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


//
// This is a simple program to exercise the FEI interface
// for the purposes of testing, code tuning and scaling studies.
//

#include <fei_iostream.hpp>

//Including the header fei_base.hpp pulls in the declaration for
//various fei classes and interfaces.

#include <fei_base.hpp>
#include <FEI_Implementation.hpp>

//Now make provision for using any one of several solver libraries. This is
//handled by the code in LibraryFactory.{hpp,cpp}.

#include <test_utils/LibraryFactory.hpp>


//we need to include some 'utils' headers which are used for
//setting up the data for the example problem.

#include <test_utils/HexBeam.hpp>
#include <test_utils/HexBeamCR.hpp>

#include <test_utils/fei_test_utils.hpp>
#include <snl_fei_Utils.hpp>

#undef fei_file
#define fei_file "beam_oldfei.cpp"
#include <fei_ErrMacros.hpp>

//==============================================================================
//Here's the main...
//==============================================================================
int main(int argc, char** argv)
{
  MPI_Comm comm;
  int numProcs, localProc;
  CHK_ERR( fei_test_utils::initialize_mpi(argc, argv, localProc, numProcs) );
  comm = MPI_COMM_WORLD;

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
	     << std::endl << "Required parameters:"<<std::endl
	     << "SOLVER_LIBRARY" << std::endl
	     << "DATA_SOURCE" << std::endl
	     << "WHICH_FEI" << std::endl
	     << "W" << std::endl << "D" << std::endl << "DofPerNode" << std::endl;
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
    std::cout << std::endl;
    std::cout << "========================================================" 
	 << std::endl;
    std::cout << "Size W: " << W << " (num-elements-along-side-of-cube)"<<std::endl;
    std::cout << "Size D: " << D << " (num-elements-along-depth-of-cube)"<<std::endl;
    std::cout << "DOF per node: " << DofPerNode <<std::endl;
    std::cout << "Num local  elements: " << hexcube.localNumElems_ << std::endl;
    std::cout << "Num global elements: " << hexcube.totalNumElems_ << std::endl;
    std::cout << "Num local  DOF: " << hexcube.numLocalDOF_ << std::endl;
    std::cout << "Num global DOF: " << hexcube.numGlobalDOF_ << std::endl;
    std::cout << "Num global CRs: " << numCRs << std::endl;
    std::cout << "========================================================" 
	 << std::endl;
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
      fei::console_out() << exc.what() << std::endl;
      ERReturn(-1);
    }
    fei.reset(new FEI_Implementation(wrapper, comm));
  }
  else if (whichFEI == "fei::FEI_Impl") {
    try {
      factory = fei::create_fei_Factory(comm, solverName.c_str());
    }
    catch (std::runtime_error& exc) {
      fei::console_out() << exc.what() << std::endl;
      ERReturn(-1);
    }
    fei = factory->createFEI(comm);
  }
  else {
    fei::console_out() << "beam ERROR, value of 'WHICH_FEI' must be 'OLDFEI' or 'fei::FEI_Impl'"<< std::endl;
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
    std::cout.setf(IOS_FIXED, IOS_FLOATFIELD);
    std::cout << "Initialization time:   " << fei_init_time << std::endl;
  }

  //Now the initialization phase is complete. Next we'll do the load phase,
  //which for this problem just consists of loading the element data
  //(element-wise stiffness arrays and load vectors) and the boundary
  //condition data.

  double start_load_time = fei::utils::cpu_time();

  CHK_ERR( HexBeam_Functions::load_elem_data(fei.get(), hexcube) );

  CHK_ERR( HexBeam_Functions::load_constraints(fei.get(), hexcube, firstLocalCRID) );

  //std::cout << "calling load_BC_data"<<std::endl;
  //CHK_ERR( load_BC_data(fei, hexcube) );

  fei->loadComplete();

  double fei_load_time = fei::utils::cpu_time() - start_load_time;

  delete hexcubeptr;

  if (localProc == 0) {
    //IOS macros are defined in fei_macros.h
    std::cout.setf(IOS_FIXED, IOS_FLOATFIELD);
    std::cout << "Coef. loading  time:    " << fei_load_time << std::endl;
    std::cout << "Total assembly time:    " << fei_init_time + fei_load_time << std::endl;
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
    if (localProc==0) std::cout << "solve returned status: " << status << std::endl;
  }

  double solve_time = fei::utils::cpu_time() - start_solve_time;

  if (localProc == 0) {
    std::cout << "Solver time:          " << solve_time << std::endl;
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
    std::cout << "beam execution successful" << std::endl;
  }

  fei.reset();

#ifndef FEI_SER
  MPI_Finalize();
#endif

  return(0);
}
