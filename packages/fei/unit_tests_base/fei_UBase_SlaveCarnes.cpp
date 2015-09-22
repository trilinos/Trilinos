
/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <FEI_config.h>

#include "fei_CommUtils.hpp"
#include "fei_SharedPtr.hpp"
#include "fei_ParameterSet.hpp"

#include "fei_base.hpp"
#include "fei_mpi.h"
#include "fei_Factory_Trilinos.hpp"

namespace {

TEUCHOS_UNIT_TEST(SlaveConstraint, Carnes)
{

#ifndef FEI_SER
  MPI_Comm comm = MPI_COMM_WORLD;
#else
  int comm = 0;
#endif

  int numProcs = fei::numProcs(comm);
  if (numProcs > 1) {
    out << "SlaveConstraint.Carnes test only runs on 1 or 2 procs."<<std::endl;
    return;
  }

  fei::ParameterSet params;
  fei::SharedPtr<fei::Factory> factory;
  fei::SharedPtr<fei::VectorSpace> vectorNodeSpace;
  fei::SharedPtr<fei::MatrixGraph> vectorMatrixGraph;
  fei::SharedPtr<fei::Matrix> vectorMatrix;
  fei::SharedPtr<fei::Vector> vectorSolution;
  fei::SharedPtr<fei::Vector> vectorRHS;
  fei::SharedPtr<fei::LinearSystem> vectorLinSys;
  fei::SharedPtr<fei::Solver> vectorSolver; 
  
  const int fei_numFields = 1;
  const int Ndim = 1;
  const int Nnpe = 2;

  std::vector<int> fei_fieldIDs(1,0);
  std::vector<int> fei_fieldSizes(1,Ndim);
  const int fei_nodeIDType = 0;

  // UNSPhysics::initialize_FEI_factory()

  std::string solverName = "Trilinos_Amesos";
  
  params.add(fei::Param("Trilinos_Solver","Amesos_Klu"));
  params.add(fei::Param("PrintTiming", true));
  params.add(fei::Param("PrintStatus", true));
  params.add(fei::Param("ComputeTrueResidual", true));
  params.add(fei::Param("ComputeVectorNorms", true));
  params.add(fei::Param("OutputLevel",3));
  params.add(fei::Param("outputLevel",3));
  
  // NOTE for debugging can put a non-empty string below (see list below)
  std::string outputLevel("FULL_LOGS"); // ALL, FULL_LOGS, MATRIX_FILES, BRIEF_LOGS, STATS
  
  factory = fei::SharedPtr<fei::Factory>(new Factory_Trilinos(comm));
  
  params.add(fei::Param("SOLVER_LIBRARY",solverName.c_str()));
  params.add(fei::Param("FEI_OUTPUT_LEVEL",outputLevel.c_str()));
  params.add(fei::Param("INCLUDE_ALL_SLAVE_CONSTRAINTS",true));
  
  factory->parameters(params);

  // UNSPhysics::build_FEI_VectorSpace()

  vectorNodeSpace = factory->createVectorSpace(comm, "vector_node_space");
  
  fei::SharedPtr<fei::VectorSpace> dummy;
  vectorMatrixGraph = factory->createMatrixGraph(vectorNodeSpace, dummy, "vector_matgraph");
  
  vectorNodeSpace->setParameters(params);
  vectorMatrixGraph->setParameters(params);
  
  vectorNodeSpace->defineFields(fei_numFields, &fei_fieldIDs[0], &fei_fieldSizes[0]);
  
  std::vector<int> idTypes(fei_numFields, fei_nodeIDType); // std constructor
  vectorNodeSpace->defineIDTypes(fei_numFields, &idTypes[0]);
  
  // UNSPhysics::build_FEI_MatrixGraph()
    
  std::vector<int> numFieldsPerID (Nnpe,fei_numFields);
  std::vector<int> fieldIDs (Nnpe*fei_numFields);
  
  const int num_elems = 4;
  
  // local IDs are listed with node slowest (field id fastest)
  // NOTE bcarnes 12/2/29: I tried listing them with node fastest, but
  // the FEI did not build the correct vector space in that case
  int count = 0;
  for (int i=0; i<Nnpe; i++) {
    for (int f=0; f<fei_numFields; f++) {
      fieldIDs[count++] = fei_fieldIDs[f];
    }
  }
  
  // this is the patten with which the local nodal info is defined for
  // the single element data structure book-keeping
  int patternID = vectorMatrixGraph->definePattern(
						   Nnpe, fei_nodeIDType,
						   &numFieldsPerID[0],
						   &fieldIDs[0]);
  // TODO generalize to multiple elem blocks?
  const int blockID = 0; // everything so far is written for a single element type
  vectorMatrixGraph->initConnectivityBlock(blockID, num_elems, patternID);
  
  std::vector<int> nodeIDs(Nnpe);

  if (Ndim==1) {

    // connectivity: 1--2--3--4--5
    
    for (int e=0; e<num_elems; e++) {

      nodeIDs[0] = e+1;
      nodeIDs[1] = e+2;

      vectorMatrixGraph->initConnectivity(blockID, e, &nodeIDs[0]);
    }
  }
  
  // build_FEI_MatrixGraph_Init_Constraints(false);
  
  const int numIDs = 2; // don't change this, a node is constrained to itself, FEI is very general
  const int offsetOfSlave = 0; // no field is constrained to other fields

  // first Ndim comps are for slave dofs, last Ndim components are for master dofs
  std::vector<double> constrained_coeffs (2*Ndim,0.0);
  double constraint_rhs = 0;

  std::vector<int> constraintIDtypes(numIDs, fei_nodeIDType);
  std::vector<int> IDs (numIDs, -1);

  int constrained_dim=0;

  const int offsetIntoSlaveField = constrained_dim;

  // Vx = constant
  constrained_coeffs[constrained_dim] = -1.0;

  if (Ndim==1) {

    constraint_rhs = 1.0;

    IDs[0]=5;
    IDs[1]=5;

    vectorMatrixGraph->initSlaveConstraint(
					   numIDs,
					   &constraintIDtypes[0],
					   &IDs[0],
					   &fieldIDs[0],
					   offsetOfSlave,
					   offsetIntoSlaveField,
					   &constrained_coeffs[0],
					   constraint_rhs);
    constraint_rhs = 0.0;

    IDs[0]=1;
    IDs[1]=1;

    vectorMatrixGraph->initSlaveConstraint(
					   numIDs,
					   &constraintIDtypes[0],
					   &IDs[0],
					   &fieldIDs[0],
					   offsetOfSlave,
					   offsetIntoSlaveField,
					   &constrained_coeffs[0],
					   constraint_rhs);
  }

  vectorMatrixGraph->initComplete();
  
  //UNSPhysics::build_FEI_LinSysAndSolver()
  
  vectorMatrix = factory->createMatrix(vectorMatrixGraph);
  
  vectorSolution = factory->createVector(vectorMatrixGraph, true);
  vectorRHS      = factory->createVector(vectorMatrixGraph);
  vectorLinSys   = factory->createLinearSystem(vectorMatrixGraph);
  
  vectorLinSys->parameters(params);
  
  vectorLinSys->setMatrix(vectorMatrix);
  vectorLinSys->setSolutionVector(vectorSolution);
  vectorLinSys->setRHS(vectorRHS);
  
  vectorSolver = factory->createSolver(NULL);

  //TransientDynamics_STK::assemble_primal_linear_system(int tstep, Real alpha, StaticOrDynamicPhysics WHICHPHYSICS)

  vectorRHS->putScalar(0.0);  
  vectorMatrix->putScalar(0.0);
  
  int len = Nnpe*Ndim;
  std::vector<int> indices (len);

  const double h =1.0/(double)num_elems;

  double* elemMat = new double[len*len];
  double** elemMat2D = new double*[len];
  for(int j=0; j<len; ++j) {
    elemMat2D[j] = &(elemMat[j*len]);
  }

  if (Ndim==1) {

    elemMat2D[0][0] = 1./h;
    elemMat2D[0][1] = -1./h;
    elemMat2D[1][0] = -1./h;
    elemMat2D[1][1] = 1./h;

    for (int e=0; e<num_elems; e++) {

      vectorMatrix->getMatrixGraph()->getConnectivityIndices(
							     blockID,
							     e,
							     len, &indices[0],
							     len);
    
      vectorMatrix->sumIn(
			  len, &indices[0], len, &indices[0],
			  elemMat2D, FEI_DENSE_COL);
    }
  }

  vectorLinSys->loadComplete();

  int solverStatus = -1;
  int itersTaken = -1;
  vectorSolver->solve(
		      vectorLinSys.get(),
		      NULL, //preconditioningMatrix
		      params,
		      itersTaken,
		      solverStatus);
  
  vectorSolution->scatterToOverlap();

  if (Ndim == 1) {
    const int numIDs = 5;
    int solnIDs[numIDs] = {1, 2, 3, 4, 5};
    double solnValues[numIDs] = {0, 0, 0, 0, 0};
    vectorSolution->copyOutFieldData(fei_fieldIDs[0], fei_nodeIDType, numIDs, solnIDs, solnValues);
    double tol = 1.e-6;
    TEUCHOS_TEST_FLOATING_EQUALITY(solnValues[0], 0.0, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(solnValues[1], 0.25, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(solnValues[2], 0.5, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(solnValues[3], 0.75, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(solnValues[4], 1.0, tol, out, success);
    double rhsValues[numIDs] = {0, 0, 0, 0, 0};
    vectorRHS->copyOutFieldData(fei_fieldIDs[0], fei_nodeIDType, numIDs, solnIDs, solnValues);
    TEUCHOS_TEST_FLOATING_EQUALITY(rhsValues[0], 0.0, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(rhsValues[1], 0.0, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(rhsValues[2], 0.0, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(rhsValues[3], 0.0, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(rhsValues[4], 0.0, tol, out, success);
  }

  std::cout << "fei solved rezone linsys in " << itersTaken << " iters "
	    << " with status = " << solverStatus << std::endl;  
}

}

