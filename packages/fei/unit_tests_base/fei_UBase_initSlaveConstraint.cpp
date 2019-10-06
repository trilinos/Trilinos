
/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_base.hpp>
#include <fei_Factory_Trilinos.hpp>

//The following test is a program submitted with a bug report by 
//a user named Rijk de Rooij at tudelft.nl.
//
TEUCHOS_UNIT_TEST(fei_UBase, slaveMatrix) 
{
  MPI_Comm comm = MPI_COMM_WORLD;

  int localProc = fei::localProc(comm);
  int numProcs  = fei::numProcs(comm);
  if (numProcs != 2) {
    std::cout<<"test Eqns_unit.feiInitSlave only runs on 2 procs. returning."<<std::endl;
    return;
  }

	bool verbose = (localProc==0);
	
	// ======================================================== //
	// 			  	C R E A T E  F A C T O R Y 					//
	// -------------------------------------------------------- //
	
	
	fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));
	fei::ParameterSet paramset;
	paramset.add(fei::Param("FEI_OUTPUT_PATH", "log/"));
	paramset.add(fei::Param("FEI_OUTPUT_LEVEL","ALL"));
	factory->parameters(paramset);
	
	// Check whether the factory allocation was succesful
	if (factory.get() == NULL) {
		std::cout << "fei::Factory allocation failed." << std::endl;

		return;
	}
	
	// Create rowspace and columnspace
	fei::SharedPtr<fei::VectorSpace> nodeSpace =
		factory->createVectorSpace(comm,NULL); 	//row space

	fei::SharedPtr<fei::VectorSpace> dummy;			// column space
	
	// Create matrixgraph from this space. Note: if the matrix is 
	// symmetric, then the columns space is not required,i.e. NULL
	fei::SharedPtr<fei::MatrixGraph> matrixGraph =
		factory->createMatrixGraph(nodeSpace, dummy, "StiffnessMatrix"); 
	
	//set parameters to nodeSpace and matrixGraph
	nodeSpace->setParameters(paramset);
	matrixGraph->setParameters(paramset);
	
	// ======================================================== //
	// 			  	S E T  U P  P R O B L E M 					//
	// -------------------------------------------------------- //
	
	// Define the fields for this problem and add to nodeSpace. In this
	// case the field is a displacement field with 3dof per node
	int fieldID = 0;
	int fieldSize = 6;
	int nodeIDType = 0;							// node type
	int blockID  = 0;
	int crIDType = 2;							// constraint type

	nodeSpace->defineFields( 1, &fieldID, &fieldSize );
	nodeSpace->defineIDTypes(1, &nodeIDType );
	nodeSpace->defineIDTypes(1, &crIDType );
	
	// Initialize element connectivities on proc 0
	if(verbose){
		// Define pattern
		int patternID = matrixGraph->definePattern(2, nodeIDType, fieldID);
		 // Initialize connectivity block
	    int error_code = matrixGraph->initConnectivityBlock(blockID, 2, patternID);
      TEUCHOS_TEST_EQUALITY(0, error_code, out, success);

	    // Initialize connectivities
	    int* nodeIDs = new int[3];
	    nodeIDs[0] = 0; nodeIDs[1] = 1; nodeIDs[2] = 2;
		TEUCHOS_TEST_EQUALITY(0, matrixGraph->initConnectivity(blockID, 0, &nodeIDs[0]) , out, success);
		TEUCHOS_TEST_EQUALITY(0, matrixGraph->initConnectivity(blockID, 1, &nodeIDs[1]) , out, success);
		
		delete [] nodeIDs;
	}
	
	// Initialize slave constraint on proc 0
	if(verbose){
		int numIDs = 2;
		int* nodeTypes= new int[2]; nodeTypes[0]= 0; nodeTypes[1] = 0;
		int* nodeIDs  = new int[2]; nodeIDs[0] = 0; nodeIDs[1] = 1;
		int* fieldIDs = new int[2]; fieldIDs[0] = 0; fieldIDs[1] = 0;
		int offsetOfSlave = 0;
		int offsetIntoSlaveField = 0;
		
		double* weights = new double[12];
		for(int i=0; i<12; i++) weights[i] = 0;
		weights[0] = -1;
		weights[6] = 1;
		
		double rhsValue = 0.0;
		
		int ierr = matrixGraph->initSlaveConstraint(numIDs,
								  nodeTypes,
								  nodeIDs,
								  fieldIDs,
								  offsetOfSlave,
								  offsetIntoSlaveField,
								  weights,
								  rhsValue);
								  
	    out << " Initialized Slave constraint with ierr = " << ierr << FEI_ENDL;
		
		delete [] fieldIDs;
		delete [] nodeTypes;
		delete [] nodeIDs;
		delete [] weights;
	}
	
	// Initialization is complete
	if(verbose)out << "Started init_complete"  << FEI_ENDL;
	TEUCHOS_TEST_EQUALITY(0, matrixGraph->initComplete(), out, success);
	if(verbose)out << "Finished init_complete"  << FEI_ENDL;
	
   // Set up linear system
   // Create StiffnesMatrix, DisplacementVector, and ForceVector based on matrixGraph

  if (verbose)  out << "Factory create Matrices start " << FEI_ENDL;
  fei::SharedPtr<fei::Vector> DisplacementVector = factory->createVector(matrixGraph, true);
  if (verbose)  out << "Factory create DisplacementVector ended " << FEI_ENDL;
  fei::SharedPtr<fei::Vector> ForceVector  = factory->createVector(matrixGraph);
  if (verbose)  out << "Factory create ForceVector ended " << FEI_ENDL;
  fei::SharedPtr<fei::LinearSystem> LinSys = factory->createLinearSystem(matrixGraph);
  if (verbose)  out << "Factory create LinSys ended " << FEI_ENDL;
  fei::SharedPtr<fei::Matrix> StiffnesMatrix = factory->createMatrix(matrixGraph);
  if (verbose)  out << "Factory create StiffnesMatrix ended " << FEI_ENDL;


	factory.reset();
}

