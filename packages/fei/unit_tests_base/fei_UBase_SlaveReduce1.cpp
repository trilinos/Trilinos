#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <FEI_config.h>

#include "fei_iostream.hpp"
#include "fei_base.hpp"
#include <fei_Factory_Trilinos.hpp>
#include "fei_ErrMacros.hpp"


namespace {

TEUCHOS_UNIT_TEST(SlaveReduce1, Rijk_de_Rooij )
{
	using std::cout;
	using std::endl;

	// ======================================================== //
	// 			I N I T I A L I Z E  M P I  D A T A 			//
	// -------------------------------------------------------- //

	int localProc = 0;
#ifndef FEI_SER
  int numProcs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  if (numProcs != 1) {
    return;
  }
#endif

  MPI_Comm comm = MPI_COMM_WORLD;
  bool verbose = (localProc==0);

	// ======================================================== //
	// 			  	C R E A T E  F A C T O R Y 					//
	// -------------------------------------------------------- //


	fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));
	fei::ParameterSet paramset;
	paramset.add(fei::Param("Trilinos_Solver", "Amesos_Klu"));
	paramset.add(fei::Param("FEI_OUTPUT_PATH", "."));
	paramset.add(fei::Param("FEI_OUTPUT_LEVEL","ALL"));
	factory->parameters(paramset);

	// Check whether the factory allocation was succesful
	TEUCHOS_TEST_INEQUALITY (factory.get() , NULL, out, success);

	// Create rowspace and columnspace
	fei::SharedPtr<fei::VectorSpace> nodeSpace = factory->createVectorSpace(comm,NULL); 	//row space

	fei::SharedPtr<fei::VectorSpace> dummy;			// column space

	// Create matrixgraph from this space. Note: if the matrix is
	// symmetric, then the columns space is not required,i.e. NULL
	fei::SharedPtr<fei::MatrixGraph> matrixGraph =
		factory->createMatrixGraph(nodeSpace, dummy, "StiffnessMatrix");

	//set parameters to nodeSpace and matrixGraph
	nodeSpace->setParameters(paramset);
	matrixGraph->setParameters(paramset);

	// ======================================================== //
	// 			   		S E T  U P  P R O B L E M				//
	// -------------------------------------------------------- //

	// Define the fields for this problem and add to nodeSpace. In this
	// case the field is a displacement field with 1dof per node
	const int fieldID = 0;
	const int fieldSize = 1;
	const int nodeIDType = 0;							// node type
	const int blockID  = 0;
  int err = 0;

	nodeSpace->defineFields( 1, &fieldID, &fieldSize );
	nodeSpace->defineIDTypes(1, &nodeIDType );

	// Initialize element connectivities on proc 0

	// Define pattern
	int patternID = matrixGraph->definePattern(2, nodeIDType, fieldID);
	// Initialize connectivity block
   err = matrixGraph->initConnectivityBlock(blockID, 2, patternID);
  TEUCHOS_TEST_EQUALITY(err, 0, out, success);

    // Initialize connectivities
  int nodeIDs[3] = {0, 1, 2};
	err = matrixGraph->initConnectivity(blockID, 0, &nodeIDs[0]);
  TEUCHOS_TEST_EQUALITY(err, 0, out, success);
	err = matrixGraph->initConnectivity(blockID, 1, &nodeIDs[1]);
  TEUCHOS_TEST_EQUALITY(err, 0, out, success);


	// Initialize slave constraint on proc 0
	int numIDs = 2;
	int nodeTypes[2] = {0, 0};
	int slaveNodeIDs[2] = {2, 1}; // Node 2 and 1
	int fieldIDs[2] = {0, 0};
	int offsetOfSlave = 0;
	int offsetIntoSlaveField = 0;

	// Introduce the weights, such that u_1=2*u_2 + rhsValue
	double weights[2] = {-1.0, 2.0};

	double rhsValue = 2.0;

	int ierr = matrixGraph->initSlaveConstraint(numIDs,
							  nodeTypes,
							  slaveNodeIDs,
							  fieldIDs,
							  offsetOfSlave,
							  offsetIntoSlaveField,
							  weights,
							  rhsValue);

    FEI_COUT << " Initialized Slave constraint with ierr = " << ierr << FEI_ENDL;

	// Initialization is complete
	if(verbose)FEI_COUT << "Started init_complete"  << FEI_ENDL;
	TEUCHOS_TEST_EQUALITY(matrixGraph->initComplete(), 0, out, success);
	if(verbose)FEI_COUT << "Finished init_complete"  << FEI_ENDL;


	// Set up linear system
	// Create StiffnessMatrix, DisplacementVector, and ForceVector based on matrixGraph
	fei::SharedPtr<fei::Matrix> StiffnessMatrix = factory->createMatrix(matrixGraph);
	fei::SharedPtr<fei::Vector> DisplacementVector = factory->createVector(matrixGraph, true);
	fei::SharedPtr<fei::Vector> ForceVector  = factory->createVector(matrixGraph);
	fei::SharedPtr<fei::LinearSystem> LinSys = factory->createLinearSystem(matrixGraph);

	// Create some solverparameters for the factory and for linsys
	fei::ParameterSet paramLin;
	paramLin.add(fei::Param("Trilinos_Solver", "Mumps"));
	paramLin.add(fei::Param("PrintTiming", false));
	paramLin.add(fei::Param("PrintStatus", false));
	factory->parameters(paramLin);

	LinSys->parameters(paramset);
	LinSys->parameters(paramLin);

	// Set StiffnessMatrix, DisplacementVector, and ForceVector to LinSys
	LinSys->setMatrix(StiffnessMatrix);
	LinSys->setSolutionVector(DisplacementVector);
	LinSys->setRHS(ForceVector);


	// Assemble Stiffness matrix
	double  elemMat[4] = {1., -1., -1., 1.};		// Element stiffness matrix
	double* elemMat2d[2] = {elemMat, elemMat+2};
	int matSize = 2;								// Size of element matrix
	int indices[2]; 						// DOF indices for the element

	// Loop through the 2 elements
	for(size_t elemID=0; elemID<2; elemID++){
		// Get connectivity indices of this element
		err = StiffnessMatrix->getMatrixGraph()->getConnectivityIndices(blockID, elemID, matSize,indices,matSize);
		TEUCHOS_TEST_EQUALITY(err, 0, out, success);

		// Sum element Matrix:
		err = StiffnessMatrix->sumIn(matSize, indices, matSize, indices, elemMat2d, FEI_DENSE_COL);
		TEUCHOS_TEST_EQUALITY(err, 0, out, success );

	}

	// Clamp node 0 using essential BC
	int nodeID = 0;
	int offsetIntoField=0;
	double prescribedValue =0.;

	err = LinSys->loadEssentialBCs(1,
							  &nodeID,
							  nodeIDType,
							  fieldID,
							  &offsetIntoField,
							  &prescribedValue);
  TEUCHOS_TEST_EQUALITY(err, 0, out, success);

	// Apply loads to nodes 1 and 2
	int *DOFIndices = new int[fieldSize];			// Indices of nodal dof within the rowspace
	double * LoadVector = new double[fieldSize];	// Load vector at the nodes

	// Loop through nodeIDs 1 and 2
	for(int nodeID=1; nodeID<3; nodeID++){
		LoadVector[0] = 10*nodeID;		// Apply load of 10 to node 1, and 20 to node 2

		// Compute the DOFIndices from the nodeID
		err = matrixGraph->getRowSpace()->getGlobalIndices(1,
														&nodeID,
														nodeIDType,
														fieldID,
														DOFIndices);
    TEUCHOS_TEST_EQUALITY(err, 0, out, success);
		// Sum load in forceVector
		err = ForceVector->sumIn(fieldSize, DOFIndices, LoadVector, 0);
		TEUCHOS_TEST_EQUALITY( err, 0, out, success);
	}

	// Complete loading the LinSys
	LinSys->loadComplete();

  fei::SharedPtr<fei::Solver> solver = factory->createSolver(NULL);

  int solverStatus = -1;
  int itersTaken = -1;
  solver->solve(LinSys.get(), NULL, paramset, itersTaken, solverStatus);

  DisplacementVector->scatterToOverlap();

  const int slaveNodeID = 2;
  double slaveSolnValue = 0;
  DisplacementVector->copyOutFieldData(fieldIDs[0], nodeTypes[0], 1, &slaveNodeID, &slaveSolnValue);
  double expectedSolnValue = 50.0;
  double tol = 1.e-6;
  TEUCHOS_TEST_FLOATING_EQUALITY(slaveSolnValue, expectedSolnValue, tol, out, success);
//  StiffnessMatrix->writeToFile("StiffnessMatrix.mtx");
//  ForceVector->writeToFile("ForceVector.vec");
}

}

