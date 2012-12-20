// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_DISTANCELAPLACIANFACTORY_DEF_HPP
#define MUELU_DISTANCELAPLACIANFACTORY_DEF_HPP

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_Import.hpp"
#include "Xpetra_ImportFactory.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_Matrix_fwd.hpp"

#include "MueLu_DistanceLaplacianFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  DistanceLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DistanceLaplacianFactory()
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  DistanceLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~DistanceLaplacianFactory() {}

  //---------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DistanceLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "Coordinates");
  }

  //---------------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DistanceLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
    RCP<MultiVector> coordinates = Get< RCP<MultiVector> >(currentLevel, "Coordinates");

    // Import coordinates of columns indices in local rows.
    RCP<const Map> colMap = A->getColMap();
    RCP<MultiVector> coordinatesWithGhostInfo = MultiVectorFactory::Build(colMap,coordinates->getNumVectors());
    RCP<const Map> rowMap = A->getRowMap();
    RCP<const Import> importer = ImportFactory::Build(rowMap, colMap);
    coordinatesWithGhostInfo->doImport(*coordinates, *importer, Xpetra::INSERT);

    // Allocate proper amount of space for new matrix
    ArrayRCP< size_t> numEntriesPerRow(A->getNodeNumRows());
    ArrayView<const LO> indices;
    ArrayView<const SC> values;
    for (size_t i=0; i<A->getNodeNumRows(); ++i) {
      A->getLocalRowView(i, indices, values);
      numEntriesPerRow[i] = values.size();
    }
    RCP<Matrix> DL = rcp(new CrsMatrixWrap(rowMap, numEntriesPerRow, Xpetra::StaticProfile));

    const SC zero=Teuchos::ScalarTraits<SC>::zero();
    ArrayRCP<const SC> xcoord = coordinatesWithGhostInfo->getData(0);
    ArrayRCP<const SC> ycoord, zcoord;
    if (coordinatesWithGhostInfo->getNumVectors() > 1)
      ycoord = coordinatesWithGhostInfo->getData(1);
    else
      ycoord = ArrayRCP<SC>(xcoord.size(),zero);
    if (coordinatesWithGhostInfo->getNumVectors() > 2)
      zcoord = coordinatesWithGhostInfo->getData(2);
    else
      zcoord = ArrayRCP<SC>(xcoord.size(),zero);

    ArrayRCP<SC> valPtr = ArrayRCP<SC>(A->getNodeMaxNumRowEntries(),zero);

    Array<GO> globalIndices(A->getNodeMaxNumRowEntries());
    for (size_t i=0; i<A->getNodeNumRows(); ++i) {
      A->getLocalRowView(i, indices, values);
      SC sum=zero;
      LO diagIndex=0;
      for (LO j=0; j<indices.size(); ++j) {
        globalIndices[j] = colMap->getGlobalElement((indices[j]));
        LO colInd = indices[j];
        if ((size_t)colInd != i) {
          SC xdelta = xcoord[i]-xcoord[colInd];
          SC ydelta = ycoord[i]-ycoord[colInd];
          SC zdelta = zcoord[i]-zcoord[colInd];
          SC distance = Teuchos::ScalarTraits<SC>::squareroot(xdelta*xdelta + ydelta*ydelta + zdelta*zdelta);
          TEUCHOS_TEST_FOR_EXCEPTION(distance==zero, Exceptions::RuntimeError,
                                     "MueLu::DistanceLaplacian: distance between two distinct nodes is zero");
          valPtr[j] = -1*Teuchos::ScalarTraits<SC>::one() / distance;
          sum -= valPtr[j];
        } else {
          diagIndex=j;
        }
      }
      valPtr[diagIndex]=sum;

      DL->insertGlobalValues(rowMap->getGlobalElement(i),
                             globalIndices.view(0,indices.size()),
                             valPtr.view(0,indices.size()));
    } //for (size_t i=0; ...

    DL->fillComplete(A->getDomainMap(),A->getRangeMap());
    Set(currentLevel,"Distance Laplacian", DL);
  } //Build

} // namespace MueLu

#define MUELU_DISTANCELAPLACIANFACTORY_SHORT
#endif // MUELU_DISTANCELAPLACIANFACTORY_DEF_HPP
