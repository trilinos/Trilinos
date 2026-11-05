// Xpetra include
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

// FROSch include
#include <ShyLU_DDFROSch_config.h>
#include <FROSch_Tools_def.hpp>
#include <FROSch_GeometricOneLevelPreconditioner_decl.hpp>
#include <FROSch_GeometricOneLevelPreconditioner_def.hpp>

// Teuchos include
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp> 
#include <Teuchos_DefaultComm.hpp> 
#include <Teuchos_OrdinalTraits.hpp> 
#include <Teuchos_Array.hpp> 

// Belos
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosXpetraAdapter.hpp>

// Test include 
#include "FROSchTestLaplace.hpp"

#include <iostream>
#include <fstream>

using UN    = unsigned;
using SC    = double;
using LO    = int;
using GO    = FROSch::DefaultGlobalOrdinal;
using NO    = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
using namespace Belos;

typedef MultiVector<SC, LO, GO, NO> multivector_type;
typedef Belos::OperatorT<multivector_type> operatort_type;
typedef Belos::LinearProblem<SC, multivector_type, operatort_type> linear_problem_type;
typedef Belos::SolverFactory<SC, multivector_type, operatort_type> solverfactory_type;
typedef Belos::SolverManager<SC, multivector_type, operatort_type> solver_type;
typedef XpetraOp<SC, LO , GO, NO> xpetraop_type;




int 
main(int argc, char* argv[]) 
{
  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Comm<int> > CommWorld = DefaultPlatform::getDefaultPlatform().getComm();
  RCP<const Comm<int> > CommSelf  = Teuchos::rcp(new MpiComm<int>(MPI_COMM_SELF));

  RCP<ParameterList> parameterList = getParametersFromXmlFile("ParameterList.xml");
  RCP<ParameterList> belosList     = sublist(parameterList,"Belos List");
  RCP<ParameterList> precList      = sublist(parameterList,"Preconditioner List");

  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();   

  int rank = CommWorld->getRank();

  // The test is designed to be executed with 2 ranks!
  assert(CommWorld->getSize() == 2); 


  /* We consider the following grid:
   *
   * +----+----+----+----+
   * | 12 | 13 | 14 | 15 |
   * +----+----+----+----+
   * |  8 |  9 | 10 | 11 |
   * +----+----+----+----+
   * 
   * +----+----+----+----+
   * |  4 |  5 |  6 |  7 | 
   * +----+----+----+----+
   * |  0 |  1 |  2 |  3 |
   * +----+----+----+----+
   *
   * Where the indices from 0 - 7 are stored on rank 0,
   * and the indices from 8 - 15 are stored on rank 1.
   *
   * With the dof enumeration
   *
   * 20---21---22---23---24
   * |    |    |    |    |
   * 15---16---17---18---19
   * |    |    |    |    |
   * 10---11---12---13---14
   * 
   * 10---11---12---13---14
   * |    |    |    |    | 
   * 2----3----5----7----9
   * |    |    |    |    |
   * 0----1----4----6----8
   * */

  // Geometric information:
  unsigned int dim              = 2;
  LO           nVerticesPerCell = 4;
  LO           nFacesPerCell    = 4;

  // Problem size:
  LO           nLocalDofs       = 15;
  GO           nGlobalDofs      = 25;
  GO           nCellsPerRow     = 4;
  double       domainSize       = 1.0;

  GO           nCells           = nCellsPerRow * nCellsPerRow;
  LO           nLocalCells      = 8;
  double       h                = domainSize / ( (double)nCellsPerRow );


  // ---------------------------------------------------------------------------------
  // Step 1: Locally owned index set

  // Remark: this would normally be provided by the finite element library
  Teuchos::RCP<Xpetra::Map<LO, GO, NO>> rowMap;
  {
    Array<GO> localIndices;

    if (rank == 0) 
      localIndices = Array<GO>{0, 1, 2, 3, 4, 5, 6, 7};

    else if (rank == 1)
      localIndices = Array<GO>{8, 9, 10, 11, 12, 13, 14, 15};

    rowMap = Xpetra::MapFactory<LO, GO, NO>::Build(
      Xpetra::UseTpetra, nCells, localIndices, 0,  CommWorld);
  }


  // ---------------------------------------------------------------------------------
  // Step 2: Dual graph

  // create the connectivity map, i,e, the dual_graph
  // Remark: this would normally be provided by the finite element library
  RCP<CrsGraph<LO, GO, NO>> dualGraph = CrsGraphFactory<LO, GO, NO>::Build(rowMap, nFacesPerCell); 
  {
    Array<Array<GO>> neighbours;
    if (rank == 0) {
      neighbours = Array<Array<GO>>{
        {1, 4},
        {0, 2, 5},
        {1, 3, 6},
        {2, 7},
        {0, 5, 8},
        {1, 4, 6, 9},
        {2, 5, 7, 10},
        {3, 6, 11}
      };
    } else if (rank == 1) {
      neighbours = Array<Array<GO>>{
        {9, 12, 4},
        {8, 10, 13, 5},
        {9, 11, 14, 6},
        {10, 15, 7},
        {8, 13},
        {9, 12, 14},
        {10, 13, 15},
        {11, 14}
      };
    }

    for (GO i = 0; i < neighbours.size(); ++i) 
      dualGraph->insertGlobalIndices(rowMap->getGlobalElement(i), neighbours[i]);

    dualGraph->fillComplete();
  }


  // ---------------------------------------------------------------------------------
  // Step 3: Locally relevant and locally owned DoF list

  // Locally relevant DoF list:
  // Remark: this would normally be provided by the finite element library
  RCP<Map<LO, GO, NO>> locallyRelevantDofs;
  {
    Array<GO> indices(nLocalDofs);
    if (rank == 0) 
      for(GO dofIndex = 0; dofIndex < nLocalDofs; ++dofIndex)
        indices[dofIndex] = dofIndex;
    else if (rank == 1)
      for(GO dofIndex = 0; dofIndex < nLocalDofs; ++dofIndex)
        indices[dofIndex] = 10 + dofIndex;

    locallyRelevantDofs= Xpetra::MapFactory<LO, GO, NO>::Build(UseTpetra, nGlobalDofs, indices, 0, CommWorld);
  }

  // Locally owned DoF list:
  // Remark: this would normally be provided by the finite element library
  RCP<Map<LO, GO, NO>> locallyOwnedDofs;
  {
    Array<GO> indices;
    if (rank == 0) 
      for(GO dofIndex = 0; dofIndex < nLocalDofs; ++dofIndex)
        indices.push_back(dofIndex);
    else if (rank == 1)
      for(GO dofIndex = nLocalDofs; dofIndex < nGlobalDofs; ++dofIndex)
        indices.push_back(dofIndex);

    locallyOwnedDofs = Xpetra::MapFactory<LO, GO, NO>::Build(UseTpetra, nGlobalDofs, indices, 0, CommWorld);
  }



  // ---------------------------------------------------------------------------------
  // Step 4: Triangulation description

  // Node Data
  //   Stores a list of vertices present in the triangulation. 
  //   (We use that the dofs correspond to the vertices in this setting)
  RCP<MultiVector<double, LO, GO, NO>> nodesVector = createNodesVector(dim, nCellsPerRow, nLocalDofs, domainSize, locallyRelevantDofs, rank);
       
  // Cell Data
  //   Store a description of each cell present in the triangulation.
  //   Each cell is described by it's vertices.
  RCP<MultiVector<GO, LO, GO, NO>> cellVector = createCellVector(nLocalCells, rowMap, rank);
    
  // Auxillary Data:
  //   Here we store the global cell index. What information is stored in the 
  //   auxillary vector depends on used FEM-Software.
  RCP<MultiVector<GO, LO, GO, NO>> auxillaryVector = createAuxillaryVector(nLocalCells, rowMap, rank);
  


  // ---------------------------------------------------------------------------------
  // Step 5: Assemble the linear system

  // Assemble system Matrix:
  // With the Node Data and the Cell Data we can assemble the system matrix
  // (this is normally done by the FEM-Software)
  RCP<CrsMatrix<SC, LO, GO, NO>> systemMatrix = 
    assembleMatrix(locallyOwnedDofs, 
                   locallyRelevantDofs, 
                   cellVector,
                   Array<GO>{0, 1, 4, 6, 8, 2, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24} /*dofs on boundary*/);

  // Assemble system RHS
  // (this is normally done by the FEM-Software)
  RCP<MultiVector<SC, LO, GO, NO>> systemRHS = 
    assembleRHS(locallyOwnedDofs, 
                   Array<GO>{0, 1, 4, 6, 8, 2, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24}, /*dofs on boundary*/
                   h);

  // Debugging:
  //auto print_out = Teuchos::getFancyOStream (Teuchos::rcpFromRef(std::cout));
  //systemMatrix->describe(*print_out, Teuchos::VERB_EXTREME);


  // ---------------------------------------------------------------------------------
  // Step 6: Communicate between ranks
  
  /* The overlapping domains should look like this:
   * rank 1:
   * 15---16---17---18---19
   * |    |    |    |    |
   * 10---11---12---13---14
   * |    |    |    |    | 
   * 5----6----7----8----9
   * |    |    |    |    |
   * 0----1----2----3----4
   * 
   * rank 0:
   * 15---16---17---18---19
   * |    |    |    |    |
   * 10---11---12---13---14
   * |    |    |    |    | 
   * 2----3----5----7----9
   * |    |    |    |    |
   * 0----1----4----6----8
   */
  
  // convert to Xpetra::Matrix
  RCP<Matrix<SC, LO, GO, NO>> k = Teuchos::rcp(new CrsMatrixWrap<SC, LO, GO, NO>(systemMatrix));

  // Test the intialize function:
  RCP<GeometricOneLevelPreconditioner<SC, LO, GO, NO>> geometricPreconditioner =
    rcp(new GeometricOneLevelPreconditioner<SC, LO, GO, NO>(k.getConst(), dualGraph, precList));


  // Test the communication function:
  geometricPreconditioner->communicateOverlappingTriangulation(nodesVector,
                                                               cellVector,
                                                               auxillaryVector,
                                                               nodesVector,
                                                               cellVector,
                                                               auxillaryVector);

  // Debugging
  //auto print_out = Teuchos::getFancyOStream (Teuchos::rcpFromRef(std::cout));
  //cellVector->describe(*print_out, Teuchos::VERB_EXTREME);


  // ---------------------------------------------------------------------------------
  // Test 1: Communcation between ranks
  // Check that the cellVector is correctly distributed between the ranks 
  // (this is the most complicated one; therefore, it is sufficient to check the cellVector).
  
  // Create the reference:
  {
    RCP<MultiVector<GO, LO, GO, NO>> referenceCellVector = 
      MultiVectorFactory<GO, LO, GO, NO>::Build(cellVector->getMap(), nVerticesPerCell);
    {
      Array<ArrayRCP<GO>> referenceCellVectorData(nVerticesPerCell);
      for (LO i = 0; i < nVerticesPerCell; ++i)
        referenceCellVectorData[i] = referenceCellVector->getDataNonConst(i);

      // The hard coded (distributed) cellVector:
      Array<Array<GO>> referenceCellDataArray;
      if (rank == 0) 
        referenceCellDataArray = 
          Array<Array<GO>>{
            Array<GO>{0, 1, 2, 3},     Array<GO>{1, 4, 3, 5},     Array<GO>{4, 6, 5, 7},     Array<GO>{6, 8, 7, 9}, 
            Array<GO>{2, 3, 10, 11},   Array<GO>{3, 5, 11, 12},   Array<GO>{5, 7, 12, 13},   Array<GO>{7, 9, 13, 14}, 
            Array<GO>{10, 11, 15, 16}, Array<GO>{11, 12, 16, 17}, Array<GO>{12, 13, 17, 18}, Array<GO>{13, 14, 18, 19}
          };
      else if (rank == 1)
        referenceCellDataArray = 
          Array<Array<GO>>{
            Array<GO>{5, 6, 10, 11},   Array<GO>{6, 7, 11, 12},   Array<GO>{7, 8, 12, 13},   Array<GO>{8, 9, 13, 14}, 
            Array<GO>{10, 11, 15, 16}, Array<GO>{11, 12, 16, 17}, Array<GO>{12, 13, 17, 18}, Array<GO>{13, 14, 18, 19}, 
            Array<GO>{0, 1, 5, 6},     Array<GO>{1, 2, 6, 7},     Array<GO>{2, 3, 7, 8},     Array<GO>{3, 4, 8, 9}
          };

      for(LO cellIndex = 0; cellIndex < 12 /*cells per subdomain*/; ++cellIndex ) 
        for(LO cellDofIndex = 0; cellDofIndex < nVerticesPerCell; ++cellDofIndex) 
          referenceCellVectorData[cellDofIndex][cellIndex] = referenceCellDataArray[cellIndex][cellDofIndex];
    }

    //auto print_out = Teuchos::getFancyOStream (Teuchos::rcpFromRef(std::cout));
    //referenceCellVector->describe(*print_out, Teuchos::VERB_EXTREME);

    // create a deep copy, so we do not destroy the original
    RCP<MultiVector<GO, LO, GO, NO>> cellVectorCopy = 
      MultiVectorFactory<GO, LO, GO, NO>::Build(cellVector, Teuchos::Copy);


    // Compare the vector:
    cellVectorCopy->update(-1.0, *referenceCellVector, 1.0);
    double cellVectorError = 0.0;
    for (LO i = 0; i < nVerticesPerCell; ++i)
      cellVectorError += cellVectorCopy->getVector(i)->norm2();

    assert(cellVectorError < 1e-8);
  }


  // ---------------------------------------------------------------------------------
  // Step 7: Build the local system

  LO nDofsOnSubdomain = 20;
 
  // Reorder the output to match the original cell enumeration
  { 
    Array<GO> indexList = extractIndexList(auxillaryVector);

    // create a deep copy:
    RCP<MultiVector<GO, LO, GO, NO>> cellVectorOriginal = 
      MultiVectorFactory<GO, LO, GO, NO>::Build(cellVector, Teuchos::Copy);

    Array<ArrayRCP<GO>>       cellVectorData(nVerticesPerCell);
    Array<ArrayRCP<const GO>> cellVectorDataSource(nVerticesPerCell);
    for (LO i = 0; i < nVerticesPerCell; ++i) {
      cellVectorData[i]       = cellVector->getDataNonConst(i);
      cellVectorDataSource[i] = cellVectorOriginal->getData(i);
    }

    for(GO cellIndex = 0; cellIndex < cellVectorData[0].size(); ++cellIndex ) 
      for(LO cellDofIndex = 0; cellDofIndex < nVerticesPerCell; ++cellDofIndex) 
        cellVectorData[cellDofIndex][cellIndex] = cellVectorDataSource[cellDofIndex][indexList[cellIndex]];
  }

  RCP<Map<LO,GO,NO>> localDofs;
  {
    Array<GO> indices(nDofsOnSubdomain);
    for(GO dofIndex = 0; dofIndex < nDofsOnSubdomain; ++dofIndex)
      indices[dofIndex] = dofIndex;

    localDofs= Xpetra::MapFactory<LO, GO, NO>::Build(UseTpetra, nDofsOnSubdomain, indices, 0, CommSelf);
  }

  Array<GO> localDofsOnBoundary;
  Array<GO> localDofsOnInterface;
  if (rank == 0) {
    localDofsOnBoundary  = Array<GO>{0, 1, 2, 4, 6, 8, 9, 10, 14, 15, 19};
    localDofsOnInterface = Array<GO>{15, 16, 17, 18, 19};
  } else if (rank == 1) {
    localDofsOnBoundary  = Array<GO>{0, 4, 5, 9, 10, 14, 15, 16, 17, 18, 19};
    localDofsOnInterface = Array<GO>{0, 1, 2, 3, 4};
  }


  // ---------------------------------------------------------------------------------
  // Step 8: Initialize the GeometricOneLevelPreconditioner
  
  // create the overlapping map:
  Array<GO> overlappingArray(nDofsOnSubdomain);
  if (rank == 0)
    for(GO i = 0; i < nDofsOnSubdomain; ++i)
      overlappingArray[i] = i; 
  else if (rank == 1)
  {
    overlappingArray[0] = 2;
    overlappingArray[1] = 3;
    overlappingArray[2] = 5;
    overlappingArray[3] = 7;
    for(GO i = 4; i < nDofsOnSubdomain; ++i)
      overlappingArray[i] = i + 5; 
  }
  Teuchos::RCP<Xpetra::Map<LO, GO, NO>> overlappingMap = 
    Xpetra::MapFactory<LO, GO, NO>::Build(Xpetra::UseTpetra, nDofsOnSubdomain, overlappingArray, 0, CommWorld);

  geometricPreconditioner->initialize(overlappingMap);


  // ---------------------------------------------------------------------------------
  // Step 9: Assemble the local system

  RCP<CrsMatrix<SC, LO, GO, NO>> localNeumannMatrix = 
    assembleMatrix(localDofs, 
                   localDofs, 
                   cellVector,
                   localDofsOnBoundary);

  RCP<CrsMatrix<SC, LO, GO, NO>> localInterfaceMatrix = 
    assembleInterfaceMatrix(localDofs,
                            cellVector,
                            localDofsOnInterface);

  // Debugging
  //auto print_out = Teuchos::getFancyOStream (Teuchos::rcpFromRef(std::cout));
  //if (rank == 0) {
  //  localInterfaceMatrix->describe(*print_out, Teuchos::VERB_EXTREME);
  //  localNeumannMatrix->describe(*print_out, Teuchos::VERB_EXTREME);
  //}


  // ---------------------------------------------------------------------------------
  // Step 10: Compute the preconditioner
  
  // convert to Xpetra::Matrix
  RCP<Matrix<SC, LO, GO, NO>> localNeumannMatrixConvert = Teuchos::rcp(new CrsMatrixWrap<SC, LO, GO, NO>(localNeumannMatrix));
  RCP<Matrix<SC, LO, GO, NO>> localInterfaceMatrixConvert = Teuchos::rcp(new CrsMatrixWrap<SC, LO, GO, NO>(localInterfaceMatrix));

  geometricPreconditioner->compute(localNeumannMatrixConvert, localInterfaceMatrixConvert);


  // ---------------------------------------------------------------------------------
  // Step 11: Solve
  
  RCP<MultiVector<SC, LO, GO, NO>> solution = MultiVectorFactory<SC, LO, GO, NO>::Build(locallyOwnedDofs, 1);
  
  // Convert the geometricPreconditioner to a belos preconditioner:
  RCP< Xpetra::Operator<SC, LO, GO, NO> > xpetraOperator = 
    rcp_dynamic_cast<Xpetra::Operator<SC, LO, GO, NO>>(geometricPreconditioner);
  RCP<operatort_type> belosPrec = rcp(new xpetraop_type(xpetraOperator));

  // Set up the linear equation system for Belos
  RCP<operatort_type> belosA = rcp(new xpetraop_type(systemMatrix));
  RCP<linear_problem_type> linear_problem (new linear_problem_type(belosA, solution, systemRHS));
  linear_problem->setProblem(solution, systemRHS);
  linear_problem->setRightPrec(belosPrec);

  // Build the Belos iterative solver
  solverfactory_type solverfactory;
  RCP<solver_type> solver = solverfactory.create(parameterList->get("Belos Solver Type","GMRES"), belosList);
  solver->setProblem(linear_problem);

  // Solve the linear syste
  solver->solve();

  // Debugging
  //auto print_out = Teuchos::getFancyOStream (Teuchos::rcpFromRef(std::cout));
  //solution->describe(*print_out, Teuchos::VERB_EXTREME);


  // ---------------------------------------------------------------------------------
  // Test 2: Iteration number
  // When something goes wrong, the number of Iterations probably goes up.
 
  int numIters = solver->getNumIters();
  assert(numIters == 3);

  
  // ---------------------------------------------------------------------------------
  // Test 3: Verify the solution
  
  // Create the reference solution
  RCP<MultiVector<SC, LO, GO, NO>> referenceSolution = 
    MultiVectorFactory<SC, LO, GO, NO>::Build(locallyOwnedDofs, 1);
  {
    ArrayRCP<SC> referenceSolutionData = referenceSolution->getDataNonConst(0);

    // The hard coded solution:
    Array<SC> referenceSolutionArray;
    if (rank == 0)
      referenceSolutionArray = Array<SC>{0.0, 0.0, 0.0, 0.0482143, 0.0, 0.0602679, 0.0, 0.0482143, 0.0, 0.0, 0.0, 0.0602679, 0.0776786, 0.0602679, 0.0};
    else if (rank == 1)
      referenceSolutionArray = Array<SC>{0.0, 0.0482143, 0.0602679, 0.0482143, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (size_t dofIndex = 0; dofIndex < locallyOwnedDofs->getLocalNumElements(); ++dofIndex )
      referenceSolutionData[dofIndex] = referenceSolutionArray[dofIndex];
  }

  solution->update(-1.0, *referenceSolution, 1.0);
  double solutionError = solution->getVector(0)->norm2();

  // Check that the solution is correct:
  assert(solutionError < 1e-6);

  return(EXIT_SUCCESS);
}
