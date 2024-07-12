// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_Q2Q1Q2COARSEGRIDFACTORY_DEF_HPP
#define MUELU_Q2Q1Q2COARSEGRIDFACTORY_DEF_HPP

#include <iostream>
#include <cmath>

#include <Teuchos_SerialDenseMatrix.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include "MueLu_Q2Q1Q2CoarseGridFactory_decl.hpp"
#include <MueLu_Level.hpp>

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Q2Q1Q2CoarseGridFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Q2Q1Q2CoarseGridFactory() {
  GetOStream(Runtime1) << "I constructed a Q2Q1Q2CoarseGridFactory object... Nothing else to do here." << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Q2Q1Q2CoarseGridFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Q2Q1Q2CoarseGridFactory() {
  // Should be empty. All destruction should be handled by Level-based get stuff and RCP
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1Q2CoarseGridFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  Input(fineLevel, "VElementList");
  Input(fineLevel, "PElementList");
  Input(fineLevel, "MElementList");

  Input(coarseLevel, "VElementList");
  Input(coarseLevel, "PElementList");
  Input(coarseLevel, "MElementList");

  // currentLevel.DeclareInput(varName_,factory_,this);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1Q2CoarseGridFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  GetOStream(Runtime1) << "Starting 'build' routine...\n";

  // This will create a list of elements on the coarse grid with a
  // predictable structure, as well as modify the fine grid list of
  // elements, if necessary (i.e. if fineLevel.GetLevelID()==0);
  // BuildCoarseGrid(fineLevel,coarseLevel);

  // This will actually build our prolongator P
  return BuildCoarseGrid(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Q2Q1Q2CoarseGridFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCoarseGrid(Level &fineLevel, Level &coarseLevel) const {
  GetOStream(Runtime1) << "starting 'BuildCoarseGrid' routine...\n";

  RCP<Teuchos::SerialDenseMatrix<GO, GO> > fineElementPDOFs = Get<RCP<Teuchos::SerialDenseMatrix<GO, GO> > >(fineLevel, "PElementList");

  GO totalFineElements = fineElementPDOFs->numRows();

  // Compute number of coarse grid elements in total:
  GO totalCoarseElements = totalFineElements / 4;
  LO nCoarseElements     = (int)sqrt(totalCoarseElements);

  // Initialize some counters:
  size_t EdgeCount   = (nCoarseElements + 1) * (nCoarseElements + 1);
  size_t CenterCount = EdgeCount + 2 * nCoarseElements * (nCoarseElements + 1);

  // Initialize arrays of the proper size:
  RCP<Teuchos::SerialDenseMatrix<GO, GO> > coarseElementVDOFs = rcp(new Teuchos::SerialDenseMatrix<GO, GO>(totalCoarseElements, 18));
  RCP<Teuchos::SerialDenseMatrix<GO, GO> > coarseElementPDOFs = rcp(new Teuchos::SerialDenseMatrix<GO, GO>(totalCoarseElements, 4));
  RCP<Teuchos::SerialDenseMatrix<GO, GO> > coarseElementMDOFs = rcp(new Teuchos::SerialDenseMatrix<GO, GO>(totalCoarseElements, 9));

  for (GO coarseElement = 0; coarseElement < totalCoarseElements; coarseElement++) {
    // ***************************************************************
    // This is less of a pain in the ass for magnetics, so I'm
    // going to build the magnetics list. The velocity follows
    // by doubling everything (and adding 1 for the y-components)
    // and the pressure follows by copying the magnetics nodes.
    // ***************************************************************

    // if (coarseElement is on the Bottom Edge)
    if (coarseElement < nCoarseElements) {
      // Bottom nodes
      (*coarseElementMDOFs)(coarseElement, 0) = coarseElement;
      (*coarseElementMDOFs)(coarseElement, 1) = coarseElement + 1;

      // Bottom edge
      (*coarseElementMDOFs)(coarseElement, 4) = EdgeCount++;

    } else {
      // Bottom Nodes
      (*coarseElementMDOFs)(coarseElement, 0) = (*coarseElementMDOFs)(coarseElement - nCoarseElements, 3);
      (*coarseElementMDOFs)(coarseElement, 1) = (*coarseElementMDOFs)(coarseElement - nCoarseElements, 2);

      // Bottom Edge
      (*coarseElementMDOFs)(coarseElement, 4) = (*coarseElementMDOFs)(coarseElement - nCoarseElements, 6);
    }

    // Right and Top Edges -- must be determined before left edge
    (*coarseElementMDOFs)(coarseElement, 5) = EdgeCount++;
    (*coarseElementMDOFs)(coarseElement, 6) = EdgeCount++;

    // if (coarseElement is on the Left Edge)
    if (coarseElement % nCoarseElements == 0) {
      // Top left node
      (*coarseElementMDOFs)(coarseElement, 3) = (*coarseElementMDOFs)(coarseElement, 0) + nCoarseElements + 1;

      // Left Edge
      (*coarseElementMDOFs)(coarseElement, 7) = EdgeCount++;

    } else {
      // Top left node
      (*coarseElementMDOFs)(coarseElement, 3) = (*coarseElementMDOFs)(coarseElement - 1, 2);

      // Left Edge
      (*coarseElementMDOFs)(coarseElement, 7) = (*coarseElementMDOFs)(coarseElement - 1, 5);
    }

    // Top right node -- Must be the last node to be determined!
    (*coarseElementMDOFs)(coarseElement, 2) = (*coarseElementMDOFs)(coarseElement, 3) + 1;

    // Center Node
    (*coarseElementMDOFs)(coarseElement, 8) = CenterCount++;

    // With Magnetics built, Pressure and Velocity follow without effort.
    // First, Velocity:
    (*coarseElementVDOFs)(coarseElement, 0)  = 2 * (*coarseElementMDOFs)(coarseElement, 0);
    (*coarseElementVDOFs)(coarseElement, 1)  = 2 * (*coarseElementMDOFs)(coarseElement, 0) + 1;
    (*coarseElementVDOFs)(coarseElement, 2)  = 2 * (*coarseElementMDOFs)(coarseElement, 1);
    (*coarseElementVDOFs)(coarseElement, 3)  = 2 * (*coarseElementMDOFs)(coarseElement, 1) + 1;
    (*coarseElementVDOFs)(coarseElement, 4)  = 2 * (*coarseElementMDOFs)(coarseElement, 2);
    (*coarseElementVDOFs)(coarseElement, 5)  = 2 * (*coarseElementMDOFs)(coarseElement, 2) + 1;
    (*coarseElementVDOFs)(coarseElement, 6)  = 2 * (*coarseElementMDOFs)(coarseElement, 3);
    (*coarseElementVDOFs)(coarseElement, 7)  = 2 * (*coarseElementMDOFs)(coarseElement, 3) + 1;
    (*coarseElementVDOFs)(coarseElement, 8)  = 2 * (*coarseElementMDOFs)(coarseElement, 4);
    (*coarseElementVDOFs)(coarseElement, 9)  = 2 * (*coarseElementMDOFs)(coarseElement, 4) + 1;
    (*coarseElementVDOFs)(coarseElement, 10) = 2 * (*coarseElementMDOFs)(coarseElement, 5);
    (*coarseElementVDOFs)(coarseElement, 11) = 2 * (*coarseElementMDOFs)(coarseElement, 5) + 1;
    (*coarseElementVDOFs)(coarseElement, 12) = 2 * (*coarseElementMDOFs)(coarseElement, 6);
    (*coarseElementVDOFs)(coarseElement, 13) = 2 * (*coarseElementMDOFs)(coarseElement, 6) + 1;
    (*coarseElementVDOFs)(coarseElement, 14) = 2 * (*coarseElementMDOFs)(coarseElement, 7);
    (*coarseElementVDOFs)(coarseElement, 15) = 2 * (*coarseElementMDOFs)(coarseElement, 7) + 1;
    (*coarseElementVDOFs)(coarseElement, 16) = 2 * (*coarseElementMDOFs)(coarseElement, 8);
    (*coarseElementVDOFs)(coarseElement, 17) = 2 * (*coarseElementMDOFs)(coarseElement, 8) + 1;

    // Lastly, Pressure:
    (*coarseElementPDOFs)(coarseElement, 0) = (*coarseElementMDOFs)(coarseElement, 0);
    (*coarseElementPDOFs)(coarseElement, 1) = (*coarseElementMDOFs)(coarseElement, 1);
    (*coarseElementPDOFs)(coarseElement, 2) = (*coarseElementMDOFs)(coarseElement, 2);
    (*coarseElementPDOFs)(coarseElement, 3) = (*coarseElementMDOFs)(coarseElement, 3);

  }  // Loop over elements

  Set(coarseLevel, "VElementList", coarseElementVDOFs);
  Set(coarseLevel, "PElementList", coarseElementPDOFs);
  Set(coarseLevel, "MElementList", coarseElementMDOFs);

  // coarseLevel.Keep("VElementList",coarseLevel.GetFactoryManager()->GetFactory("VElementList").get());
  // coarseLevel.Keep("PElementList",coarseLevel.GetFactoryManager()->GetFactory("PElementList").get());
  // coarseLevel.Keep("MElementList",coarseLevel.GetFactoryManager()->GetFactory("MElementList").get());

}  // BuildCoarseGrid

}  // namespace MueLu

#define MUELU_Q2Q1Q2COARSEGRIDFACTORY_SHORT
#endif  // MUELU_Q2Q1Q2COARSEGRIDFACTORY_DEF_HPP
