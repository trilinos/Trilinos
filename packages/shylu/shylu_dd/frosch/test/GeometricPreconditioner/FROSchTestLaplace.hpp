// Teuchos include
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp> 
#include <Teuchos_DefaultComm.hpp> 
#include <Teuchos_OrdinalTraits.hpp> 
#include <Teuchos_Array.hpp> 

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

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using UN    = unsigned;
using SC    = double;
using LO    = int;
using GO    = FROSch::DefaultGlobalOrdinal;
using NO    = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;



bool
contains(Teuchos::Array<GO> array, GO value) {
  for(GO i : array)
    if (i == value)
      return true;
  return false;
}



Array<GO>
extractIndexList(RCP<MultiVector<GO, LO, GO, NO>>  mv) {
  const size_t m = mv->getLocalLength();

  Array<GO> indexList(m);
  for (unsigned int j = 0; j < m; ++j)
    indexList[j] = j;

  {
    auto data = mv->getData(0);
    std::sort(indexList.begin(),
              indexList.end(),
              [data](unsigned int a, unsigned int b) {
                return data[a] < data[b];
              });
  }

  return indexList;;
}


template<typename Number>
Array<Array<Number>>
extractRcpMultiVector(RCP<MultiVector<Number, LO, GO, NO>>  mv) {
  const size_t n = mv->getNumVectors();
  const size_t m = mv->getLocalLength();

  Array<Array<Number>> list(m, Array<Number>(n));
  for (unsigned int i = 0; i < n; ++i)
    {
      auto data = mv->getData(i);

      for (unsigned int j = 0; j < m; ++j)
        list[j][i] = data[j];
    }

  return list;
}



RCP<MultiVector<double, LO, GO, NO>> 
createNodesVector(
  unsigned int dim,
  LO nCellsPerRow,
  LO nLocalDofs,
  GO domainSize,
  RCP<Map<LO,GO,NO> > locallyRelevantDofs,
  int rank) {
  RCP<MultiVector<double, LO, GO, NO>> 
  nodesVector = 
    MultiVectorFactory<double, LO, GO, NO>::Build(locallyRelevantDofs, dim);

  // Create an Array, such we can access its data
  Array<ArrayRCP<double>> nodesVectorData(dim);
  for (unsigned int i = 0; i < dim; ++i)
    nodesVectorData[i] = nodesVector->getDataNonConst(i);

  // Fill node_vector: 
  // (this part is normally provided by the FEM-Software)
  if (rank == 0)
    for(LO dofIndex = 0; dofIndex < nLocalDofs; ++dofIndex) {
      nodesVectorData[0][dofIndex] = (dofIndex % (nCellsPerRow + 1)) / ( (double)nCellsPerRow ); // x-component of the node location
      nodesVectorData[1][dofIndex] = (dofIndex / (nCellsPerRow + 1)) / ( (double)nCellsPerRow ); // y-component of the node location
    }
  else if (rank == 1)
    for(LO dofIndex = 0; dofIndex < nLocalDofs; ++dofIndex) {
      nodesVectorData[0][dofIndex] =  (dofIndex % (nCellsPerRow + 1)) / ( (double)nCellsPerRow );                       // x-component of the node location
      nodesVectorData[1][dofIndex] = ((dofIndex / (nCellsPerRow + 1)) / ( (double)nCellsPerRow )) + (domainSize / 2.0); // y-component of the node location
    }

  return nodesVector;
}
   


RCP<MultiVector<GO, LO, GO, NO>>
createCellVector(LO nLocalCells, RCP<Map<LO, GO, NO> > rowMap, int rank) {
  // TODO: Generalize this function
  LO nVerticesPerCell = 4;

  RCP<MultiVector<GO, LO, GO, NO>> cellVector = 
    MultiVectorFactory<GO, LO, GO, NO>::Build(rowMap, nVerticesPerCell);

  // Create an Array, such we can access its data
  Array<ArrayRCP<GO>> cellVectorData(nVerticesPerCell);
  for (LO i = 0; i < nVerticesPerCell; ++i)
    cellVectorData[i] = cellVector->getDataNonConst(i);

  // Prepare the information:
  // (this part is normally provided by the FEM-Software)
  Array<Array<GO>> cellDataArray;
  if (rank == 0) 
    cellDataArray = 
      Array<Array<GO>>{
        Array<GO>{0,1,2,3},   Array<GO>{1,4,3,5},   Array<GO>{4,6,5,7},   Array<GO>{6,8,7,9},
        Array<GO>{2,3,10,11}, Array<GO>{3,5,11,12}, Array<GO>{5,7,12,13}, Array<GO>{7,9,13,14}
      };
  else if (rank == 1)
    cellDataArray = 
      Array<Array<GO>>{
        Array<GO>{10,11,15,16}, Array<GO>{11,12,16,17}, Array<GO>{12,13,17,18}, Array<GO>{13,14,18,19},
        Array<GO>{15,16,20,21}, Array<GO>{16,17,21,22}, Array<GO>{17,18,22,23}, Array<GO>{18,19,23,24}
      };

  // copy the data into the cell_data_vector:
  for(LO cellIndex = 0; cellIndex < nLocalCells; ++cellIndex ) 
    for(LO cellDofIndex = 0; cellDofIndex < nVerticesPerCell; ++cellDofIndex) 
      cellVectorData[cellDofIndex][cellIndex] = cellDataArray[cellIndex][cellDofIndex];

  return cellVector;
}



RCP<MultiVector<GO, LO, GO, NO>> 
createAuxillaryVector(LO nLocalCells, RCP<Map<LO, GO, NO> > rowMap, int rank) {
  RCP<MultiVector<GO, LO, GO, NO>> auxillaryVector = 
    MultiVectorFactory<GO, LO, GO, NO>::Build(rowMap, 1);

  ArrayRCP<GO> auxillaryVectorData = auxillaryVector->getDataNonConst(0);
  if (rank == 0) 
    for(LO cellIndex = 0; cellIndex < nLocalCells; ++cellIndex ) 
      auxillaryVectorData[cellIndex] = cellIndex;
  else if (rank == 1)
    for(LO cellIndex = 0; cellIndex < nLocalCells; ++cellIndex ) 
      auxillaryVectorData[cellIndex] = cellIndex + nLocalCells;

  return auxillaryVector;
}



RCP<CrsMatrix<SC, LO, GO, NO>> 
assembleMatrix(
    RCP<Map<LO, GO, NO>>              locallyOwnedDofs, 
    RCP<Map<LO, GO, NO>>              locallyRelevantDofs,
    RCP<MultiVector<GO, LO, GO, NO>>  cellVector,
    Array<GO>                         dofsOnBoundary
){
  GO nGlobalDofs      = locallyOwnedDofs->getGlobalNumElements();
  LO nVerticesPerCell = 4;

  RCP<CrsMatrix<SC, LO, GO, NO>> returnMatrix = 
    CrsMatrixFactory<SC, LO, GO, NO>::Build(locallyOwnedDofs, nGlobalDofs);


  Array<Array<SC>> localStiffnessMatrix = 
    {Array<SC>{ 4.0/6.0, -1.0/6.0, -1.0/6.0, -2.0/6.0}, 
     Array<SC>{-1.0/6.0,  4.0/6.0, -1.0/3.0, -1.0/6.0},
     Array<SC>{-1.0/6.0, -1.0/3.0,  4.0/6.0, -1.0/6.0},
     Array<SC>{-2.0/6.0, -1.0/6.0, -1.0/6.0,  4.0/6.0}};

  Array<Array<GO>> cellDataArray = extractRcpMultiVector(cellVector);

  // for simplicity we first create a dense matrix and create a sparse matrix out of that later.
  Array<Array<SC>> fullMatrix(nGlobalDofs, Array<SC>(nGlobalDofs));
  for(size_t cell = 0; cell < cellVector->getLocalLength(); ++cell)
    for(LO i = 0; i < nVerticesPerCell; ++i)
      for(LO j = 0; j < nVerticesPerCell; ++j)
          fullMatrix[cellDataArray[cell][i]][cellDataArray[cell][j]] += localStiffnessMatrix[i][j];

  // Apply boundary conditions
  for(GO i : dofsOnBoundary) 
    for(GO j = 0; j < nGlobalDofs; ++j) 
      if (i != j) {
        fullMatrix[i][j] = 0;
        fullMatrix[j][i] = 0;
      }

  // Create the distributed sparse matrix 
  for (auto globalRow : locallyRelevantDofs->getLocalElementList()) {
    Array<GO> cols;
    Array<SC> vals;
    for (GO i = 0; i < nGlobalDofs; ++i) {
      if (std::abs(fullMatrix[globalRow][i]) > 1e-6) {
        cols.push_back(i);
        vals.push_back(fullMatrix[globalRow][i]);
      }
    }
    returnMatrix->insertGlobalValues(globalRow, cols, vals);
  }
  returnMatrix->fillComplete();

  return returnMatrix;
}



RCP<MultiVector<SC, LO, GO, NO>>
assembleRHS(
    RCP<Map<LO, GO, NO>>              locallyOwnedDofs,
    Array<GO>                         dofsOnBoundary,
    double                            h /* step width */
) {
  RCP<MultiVector<SC, LO, GO, NO>> returnVector = MultiVectorFactory<SC, LO, GO, NO>::Build(locallyOwnedDofs, 1);
  returnVector->putScalar( h * h );

  for (auto gid : dofsOnBoundary) 
    if (locallyOwnedDofs->isNodeGlobalElement(gid))
      returnVector->replaceLocalValue(locallyOwnedDofs->getLocalElement(gid), 0, 0.0); 

  return returnVector;
}



RCP<CrsMatrix<SC, LO, GO, NO>> 
assembleInterfaceMatrix(
    RCP<Map<LO, GO, NO>>              locallyOwnedDofs,
    RCP<MultiVector<GO, LO, GO, NO>>  cellVector,
    Array<GO>                         dofsOnInterface
  ){
  GO nGlobalDofs      = locallyOwnedDofs->getGlobalNumElements();
  LO nVerticesPerCell = 4;

  RCP<CrsMatrix<SC, LO, GO, NO>> returnMatrix = 
    CrsMatrixFactory<SC, LO, GO, NO>::Build(locallyOwnedDofs, nGlobalDofs);

  Array<Array<GO>> cellDataArray = extractRcpMultiVector(cellVector);

  // for simplicity we first create a dense matrix and create a sparse matrix out of that later.
  Array<Array<SC>> fullMatrix(nGlobalDofs, Array<SC>(nGlobalDofs));
  for(size_t cell = 0; cell < cellVector->getLocalLength(); ++cell)
    for(LO i = 0; i < nVerticesPerCell; ++i) {
      if(!contains(dofsOnInterface, cellDataArray[cell][i]))
        continue;
      for(LO j = 0; j < nVerticesPerCell; ++j)
        if(contains(dofsOnInterface, cellDataArray[cell][j]))
          fullMatrix[cellDataArray[cell][i]][cellDataArray[cell][j]] += 0.0625;
    }

  // Create the distributed sparse matrix 
  for (auto globalRow : locallyOwnedDofs->getLocalElementList()) {
    Array<GO> cols;
    Array<SC> vals;
    for (GO i = 0; i < nGlobalDofs; ++i) {
      if (std::abs(fullMatrix[globalRow][i]) > 1e-6) {
        cols.push_back(i);
        vals.push_back(fullMatrix[globalRow][i]);
      }
    }
    returnMatrix->insertGlobalValues(globalRow, cols, vals);
  }
  returnMatrix->fillComplete();

  return returnMatrix;
}

