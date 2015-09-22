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
// This is a simple program to exercise FEI classes for the
// purposes of testing, code tuning and scaling studies.
//

#include "fei_iostream.hpp" 

//Including the header fei_base.hpp pulls in the declaration for
//various classes in the 'fei' namespace.

#include "fei_base.hpp"


//Make provision for using any one of several solver libraries. This is
//handled by the code in LibraryFactory.{hpp,cpp}.

#include "test_utils/LibraryFactory.hpp"


//we need to include some 'test_utils' headers which are used for
//setting up the data for the hex-beam example problem.

#include "test_utils/HexBeam.hpp"
#include "test_utils/HexBeamCR.hpp"


//fei_test_utils.hpp declares things like the intialize_mpi function.

#include "test_utils/fei_test_utils.hpp"

//
//Include definitions of macros like 'CHK_ERR' to call functions and check
//the return code.
//
#undef fei_file
#define fei_file "beam.cpp"
#include "fei_ErrMacros.hpp"

//==============================================================================
//Here's the main...
//==============================================================================
int main(int argc, char** argv)
{
  MPI_Comm comm;
  int numProcs, localProc;
  CHK_ERR( fei_test_utils::initialize_mpi(argc, argv, localProc, numProcs) );
  comm = MPI_COMM_WORLD;

  double start_time = fei::utils::cpu_time();

  //read input parameters from a file specified on the command-line with
  // '-i file'
  std::vector<std::string> stdstrings;
  CHK_ERR( fei_test_utils::get_filename_and_read_input(argc, argv,
						comm, localProc,
						stdstrings) );
  fei::ParameterSet params;
  //parse the strings from the input file into a fei::ParameterSet object.
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
	     << std::endl << "Required parameters:"<<std::endl
	     << "SOLVER_LIBRARY" << std::endl
	     << "DATA_SOURCE" << std::endl
	     << "CONSTRAINT_FORM" << std::endl
	     << "W" << std::endl << "D" << std::endl << "DofPerNode" << std::endl;

#ifndef FEI_SER
    MPI_Finalize();
#endif

    return(-1);
  }

  params.getStringParamValue("CONSTRAINT_FORM", constraintform);

  bool slave_constraints = (constraintform == "slaves");

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
    //macros std::cout and std::endl are aliases for std::cout and std::endl,
    //defined in fei_iostream.hpp.
    std::cout << std::endl;
    std::cout << "========================================================\n";
    std::cout << "FEI version: " << fei::utils::version() << std::endl;
    std::cout << "--------------------------------------------------------\n";
    std::cout << "Size W: " << W << " (num-elements-along-side-of-cube)\n";
    std::cout << "Size D: " << D << " (num-elements-along-depth-of-cube)\n";
    std::cout << "DOF per node: " << DofPerNode <<"\n";
    std::cout << "Num local  elements: " << hexcube.localNumElems_ << "\n";
    std::cout << "Num global elements: " << hexcube.totalNumElems_ << "\n";
    std::cout << "Num local  DOF: " << hexcube.numLocalDOF_ << "\n";
    std::cout << "Num global DOF: " << hexcube.numGlobalDOF_ << "\n";
    std::cout << "Num global CRs: " << numCRs << "\n";
    std::cout << "========================================================" 
	 << std::endl;
  }

  double start_init_time = fei::utils::cpu_time();

  fei::SharedPtr<fei::Factory> factory;
  try {
    factory = fei::create_fei_Factory(comm, solverName.c_str());
  }
  catch (...) {
    std::cout << "library " << solverName << " not available."<<std::endl;

#ifndef FEI_SER
    MPI_Finalize();
#endif

    return(-1);
  }
  if (factory.get() == NULL) {
    std::cout << "fei::Factory allocation failed." << std::endl;

#ifndef FEI_SER
    MPI_Finalize();
#endif

    return(-1);
  }

  factory->parameters(params);

  fei::SharedPtr<fei::VectorSpace> nodeSpace =
    factory->createVectorSpace(comm, "beam");

  fei::SharedPtr<fei::VectorSpace> dummy;
  fei::SharedPtr<fei::MatrixGraph> matrixGraph =
    factory->createMatrixGraph(nodeSpace, dummy, "beam");

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

  //HexBeam_Functions::init_elem_connectivities demonstrates the process
  //of initializing element-to-node connectivities on the fei::MatrixGraph.
  CHK_ERR( HexBeam_Functions::
	   init_elem_connectivities( matrixGraph.get(), hexcube ) );

  //HexBeam_Functions::init_shared_nodes demonstrates the process of
  //initializing shared nodes on the fei::MatrixGraph.
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
    std::cout.setf(IOS_FIXED, IOS_FLOATFIELD);
    std::cout << "Initialization cpu time:   " << fei_init_time << std::endl;
  }

  //Now the initialization phase is complete. Next we'll do the load phase,
  //which for this problem just consists of loading the element data
  //(element-wise stiffness arrays and load vectors) and the boundary
  //condition data.

  double fei_creatematrix_start_time = fei::utils::cpu_time();

  fei::SharedPtr<fei::Matrix> mat = factory->createMatrix(matrixGraph);

  double fei_creatematrix_time = fei::utils::cpu_time() - fei_creatematrix_start_time;
  if (localProc == 0) {
    std::cout.setf(IOS_FIXED, IOS_FLOATFIELD);
    std::cout << "Create-Matrix cpu time:   " << fei_creatematrix_time << std::endl;
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
    //IOS macros are defined in fei_stdinc.h
    std::cout.setf(IOS_FIXED, IOS_FLOATFIELD);
    std::cout << "Coef. loading cpu time:    " << fei_load_time << std::endl;
    std::cout << "Total assembly wall time:   "
	 << fei_init_time + fei_creatematrix_time + fei_load_time << std::endl;
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
    solver_name_value.empty() ? NULL : solver_name_value.c_str();

  fei::SharedPtr<fei::Solver> solver = factory->createSolver(charptr_solvername);

  int status;
  int itersTaken = 0;

  if (localProc==0) std::cout << "solve..." << std::endl;
  double start_solve_time = fei::utils::cpu_time();

  int err = solver->solve(linSys.get(),
			  NULL, //preconditioningMatrix
			  params,
			  itersTaken,
			  status);

  double solve_time = fei::utils::cpu_time()-start_solve_time;

  if (err!=0 && localProc==0) {
    std::cout << "solve returned err: " << err <<", status: "
			   << status << std::endl;
  }

  if (localProc==0) {
    std::cout << " cpu-time in solve: " << solve_time << std::endl;
  }

  CHK_ERR( solnVec->scatterToOverlap() );

  delete hexcubeptr;

  double elapsed_cpu_time = fei::utils::cpu_time() - start_time;

  //The following IOS_... macros are defined in base/fei_iostream.hpp
  std::cout.setf(IOS_FIXED, IOS_FLOATFIELD);
  if (localProc==0) {
    std::cout << "Proc0 cpu times (seconds):" << std::endl
	 << "   FEI initialize:    " << fei_init_time << std::endl
         << "   FEI create-matrix: " << fei_creatematrix_time << std::endl
         << "   FEI load:          " << fei_load_time << std::endl
         << "      solve:          " << solve_time << std::endl
         << "Total program time:   " << elapsed_cpu_time << std::endl;

  }


#ifndef FEI_SER
    MPI_Finalize();
#endif

  return(0);
}

