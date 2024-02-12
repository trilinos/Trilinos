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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_GEOINTERPFACTORY_DEF_HPP
#define MUELU_GEOINTERPFACTORY_DEF_HPP

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

#include "MueLu_GeoInterpFactory_decl.hpp"
#include <MueLu_Level.hpp>

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GeoInterpFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GeoInterpFactory() {
  GetOStream(Runtime1) << "I constructed a GeoInterpFactory object... Nothing else to do here." << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GeoInterpFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~GeoInterpFactory() {
  // Should be empty. All destruction should be handled by Level-based get stuff and RCP
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void GeoInterpFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  Input(fineLevel, "A");
  Input(fineLevel, "A00");
  Input(fineLevel, "A10");
  Input(fineLevel, "A20");

  Input(fineLevel, "VElementList");
  Input(fineLevel, "PElementList");
  Input(fineLevel, "MElementList");

  Input(coarseLevel, "VElementList");
  Input(coarseLevel, "PElementList");
  Input(coarseLevel, "MElementList");

  /*
      coarseLevel.DeclareInput("VElementList",coarseLevel.GetFactoryManager()->GetFactory("VElementList").get(),this);
      coarseLevel.DeclareInput("PElementList",coarseLevel.GetFactoryManager()->GetFactory("PElementList").get(),this);
      coarseLevel.DeclareInput("MElementList",coarseLevel.GetFactoryManager()->GetFactory("PElementList").get(),this);

      fineLevel.DeclareInput("VElementList",fineLevel.GetFactoryManager()->GetFactory("VElementList").get(),this);
      fineLevel.DeclareInput("PElementList",fineLevel.GetFactoryManager()->GetFactory("PElementList").get(),this);
      fineLevel.DeclareInput("MElementList",fineLevel.GetFactoryManager()->GetFactory("PElementList").get(),this);
  */

  // currentLevel.DeclareInput(varName_,factory_,this);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void GeoInterpFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  GetOStream(Runtime1) << "Starting 'build' routine...\n";

  // This will create a list of elements on the coarse grid with a
  // predictable structure, as well as modify the fine grid list of
  // elements, if necessary (i.e. if fineLevel.GetLevelID()==0);
  // BuildCoarseGrid(fineLevel,coarseLevel);

  // This will actually build our prolongator P
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void GeoInterpFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level &fineLevel, Level &coarseLevel) const {
  typedef Teuchos::SerialDenseMatrix<GO, GO> SerialDenseMatrixType;

  GetOStream(Runtime1) << "Starting 'BuildP' routine...\n";

  // DEBUG
  // Teuchos::FancyOStream fout(*GetOStream(Runtime1));
  // fineLevel.print(fout,Teuchos::VERB_HIGH);

  // Get finegrid element lists
  RCP<SerialDenseMatrixType> fineElementPDOFs = Get<RCP<SerialDenseMatrixType> >(fineLevel, "PElementList");
  RCP<SerialDenseMatrixType> fineElementVDOFs = Get<RCP<SerialDenseMatrixType> >(fineLevel, "VElementList");
  RCP<SerialDenseMatrixType> fineElementMDOFs = Get<RCP<SerialDenseMatrixType> >(fineLevel, "MElementList");

  // DEBUG
  GetOStream(Runtime1) << "done getting fine level elements...\n";
  GetOStream(Runtime1) << "getting coarse level elements...\n";
  // coarseLevel.print(fout,Teuchos::VERB_HIGH);

  // Get coarse grid element lists
  RCP<Teuchos::SerialDenseMatrix<GO, GO> > coarseElementVDOFs,
      coarseElementPDOFs,
      coarseElementMDOFs;

  coarseLevel.Get("VElementList", coarseElementVDOFs, coarseLevel.GetFactoryManager()->GetFactory("VElementList").get());
  coarseLevel.Get("PElementList", coarseElementPDOFs, coarseLevel.GetFactoryManager()->GetFactory("PElementList").get());
  coarseLevel.Get("MElementList", coarseElementMDOFs, coarseLevel.GetFactoryManager()->GetFactory("MElementList").get());

  GetOStream(Runtime1) << "computing various numbers...\n";
  // Number of elements?
  GO totalFineElements   = fineElementMDOFs->numRows();
  LO nFineElements       = (int)sqrt(totalFineElements);
  GO totalCoarseElements = coarseElementMDOFs->numRows();
  LO nCoarseElements     = (int)sqrt(totalCoarseElements);

  // Set sizes for *COARSE GRID*
  GO nM = (2 * nCoarseElements + 1) * (2 * nCoarseElements + 1);
  GO nV = 2 * nM;
  GO nP = (nCoarseElements + 1) * (nCoarseElements + 1);

  // Get the row maps for the Ps
  RCP<Matrix> fineA00 = Get<RCP<Matrix> >(fineLevel, "A00");
  RCP<Matrix> fineA10 = Get<RCP<Matrix> >(fineLevel, "A10");
  RCP<Matrix> fineA20 = Get<RCP<Matrix> >(fineLevel, "A20");

  GetOStream(Runtime1) << "creating coarse grid maps...\n";

  RCP<const Map> rowMapforPV = fineA00->getRowMap();
  RCP<const Map> rowMapforPP = fineA10->getRowMap();
  RCP<const Map> rowMapforPM = fineA20->getRowMap();

  GO fNV = rowMapforPV->getGlobalNumElements();
  GO fNP = rowMapforPP->getGlobalNumElements();
  GO fNM = rowMapforPM->getGlobalNumElements();

  // Get the comm for the maps
  RCP<const Teuchos::Comm<int> > comm = rowMapforPV->getComm();

  // Create rowMap for P
  RCP<Matrix> FineA         = Factory::Get<RCP<Matrix> >(fineLevel, "A");
  RCP<const Map> rowMapforP = FineA->getRowMap();

  // Create colMaps for the coarse grid
  RCP<const Map> colMapforPV = Xpetra::MapFactory<LO, GO>::createUniformContigMap(Xpetra::UseTpetra, nV, comm);
  RCP<const Map> colMapforPP = Xpetra::MapFactory<LO, GO>::createUniformContigMap(Xpetra::UseTpetra, nP, comm);
  RCP<const Map> colMapforPM = Xpetra::MapFactory<LO, GO>::createUniformContigMap(Xpetra::UseTpetra, nM, comm);

  GetOStream(Runtime1) << "creating coarse grid matrices...\n";
  // Create our final output Ps for the coarseGrid
  size_t maxEntriesPerRowV = 9,  // No overlap of VX and VY
      maxEntriesPerRowP    = 4,
         maxEntriesPerRowM = 9;

  RCP<Matrix> P  = rcp(new CrsMatrixWrap(rowMapforP, maxEntriesPerRowV));
  RCP<Matrix> PV = rcp(new CrsMatrixWrap(rowMapforPV, maxEntriesPerRowV));
  RCP<Matrix> PP = rcp(new CrsMatrixWrap(rowMapforPP, maxEntriesPerRowP));
  RCP<Matrix> PM = rcp(new CrsMatrixWrap(rowMapforPM, maxEntriesPerRowM));

  //*****************************************************************/
  //
  // All 25 fine grid dofs are completely determined by the coarse
  // grid element in which they reside! So I should loop over coarse
  // grid elements and build 25 rows at a time! If that's not
  // ridiculous... I just have to be careful about duplicates on
  // future elements! But duplicates are easy - just the bottom and
  // left edges.
  //
  //
  // Looking at a fine grid patch, define the following Local-Global
  // relationship (magnetics as an example):
  //
  // Bottom-Left Corner:
  // 0 -> (*fineElementMDOFs)(fineElement[0],0)
  //
  // Bottom Edge:
  // 1 -> (*fineElementMDOFs)(fineElement[0],4)
  // 2 -> (*fineElementMDOFs)(fineElement[0],2)
  // 3 -> (*fineElementMDOFs)(fineElement[1],4)
  // 4 -> (*fineElementMDOFs)(fineElement[1],2)
  //
  // Left Edge:
  // 5 -> (*fineElementMDOFs)(fineElement[0],7)
  // 6 -> (*fineElementMDOFs)(fineElement[0],3)
  // 7 -> (*fineElementMDOFs)(fineElement[2],7)
  // 8 -> (*fineElementMDOFs)(fineElement[2],3)
  //
  // All the rest:
  // 9 -> (*fineElementMDOFs)(fineElement[3],0)
  // 10 -> (*fineElementMDOFs)(fineElement[3],1)
  // 11 -> (*fineElementMDOFs)(fineElement[3],2)
  // 12 -> (*fineElementMDOFs)(fineElement[3],3)
  // 13 -> (*fineElementMDOFs)(fineElement[0],5)
  // 14 -> (*fineElementMDOFs)(fineElement[0],6)
  // 15 -> (*fineElementMDOFs)(fineElement[1],5)
  // 16 -> (*fineElementMDOFs)(fineElement[1],6)
  // 17 -> (*fineElementMDOFs)(fineElement[2],5)
  // 18 -> (*fineElementMDOFs)(fineElement[2],6)
  // 19 -> (*fineElementMDOFs)(fineElement[3],5)
  // 20 -> (*fineElementMDOFs)(fineElement[3],6)
  // 21 -> (*fineElementMDOFs)(fineElement[0],8)
  // 22 -> (*fineElementMDOFs)(fineElement[1],8)
  // 23 -> (*fineElementMDOFs)(fineElement[2],8)
  // 24 -> (*fineElementMDOFs)(fineElement[3],8)
  //
  //*****************************************************************/

  size_t nnz = 0;  // Just to make my copy-paste life easier...
  Teuchos::ArrayRCP<GO> colPtrV(maxEntriesPerRowV, 0);
  Teuchos::ArrayRCP<GO> colPtrM(maxEntriesPerRowM, 0);
  Teuchos::ArrayRCP<SC> valPtrM(maxEntriesPerRowM, 0.);

  Teuchos::ArrayRCP<GO> colPtrP(maxEntriesPerRowP, 0);
  Teuchos::ArrayRCP<SC> valPtrP(maxEntriesPerRowP, 0.);

  // About which fine-grid elements do we care?
  GO fineElement[4] = {0, 1, nFineElements, nFineElements + 1};

  GetOStream(Runtime1) << "start building matrices...\n";
  GetOStream(Runtime1) << "nCoarseElements = " << nCoarseElements << std::endl;

  for (GO coarseElement = 0; coarseElement < totalCoarseElements; coarseElement++) {
    // We don't really care what is shared with future elements -
    // we know a priori what needs to be filled for the element
    // we're dealing with, depending of if it it's on the bottom
    // edge, the left edge, or in the middle. Thoses should be the
    // only cases we care about, and they should require a
    // superset of the work required for an interior node.

    // if (CoarseElement is on bottom edge)
    if (coarseElement < nCoarseElements) {
      // fill in the bottom edge of the element patch
      //  FP = 1
      colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 0);
      valPtrM[0] = 0.375;
      colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 1);
      valPtrM[1] = -0.125;
      colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 4);
      valPtrM[2] = 0.75;

      nnz = 3;
      PM->insertGlobalValues((*fineElementMDOFs)(fineElement[0], 4), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0];
      colPtrV[1] = 2 * colPtrM[1];
      colPtrV[2] = 2 * colPtrM[2];

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 8), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0] + 1;
      colPtrV[1] = 2 * colPtrM[1] + 1;
      colPtrV[2] = 2 * colPtrM[2] + 1;

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 9), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      // FPr = 1
      colPtrP[0] = (*coarseElementPDOFs)(coarseElement, 0);
      valPtrP[0] = 0.5;
      colPtrP[1] = (*coarseElementPDOFs)(coarseElement, 1);
      valPtrP[1] = 0.5;

      nnz = 2;
      PP->insertGlobalValues((*fineElementPDOFs)(fineElement[0], 1), colPtrP.view(0, nnz), valPtrP.view(0, nnz));

      // FP = 2
      colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 4);
      valPtrM[0] = 1.0;

      nnz = 1;
      PM->insertGlobalValues((*fineElementMDOFs)(fineElement[0], 1), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0];

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 2), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0] + 1;

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 3), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      // FPr = 2
      colPtrP[0] = (*coarseElementPDOFs)(coarseElement, 1);
      valPtrP[0] = 1.0;

      nnz = 1;
      PP->insertGlobalValues((*fineElementPDOFs)(fineElement[1], 1), colPtrP.view(0, nnz), valPtrP.view(0, nnz));

      // FP = 3
      colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 0);
      valPtrM[0] = -0.125;
      colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 1);
      valPtrM[1] = 0.375;
      colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 4);
      valPtrM[2] = 0.75;

      nnz = 3;
      PM->insertGlobalValues((*fineElementMDOFs)(fineElement[1], 4), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0];
      colPtrV[1] = 2 * colPtrM[1];
      colPtrV[2] = 2 * colPtrM[2];

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[1], 8), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0] + 1;
      colPtrV[1] = 2 * colPtrM[1] + 1;
      colPtrV[2] = 2 * colPtrM[2] + 1;

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[1], 9), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      // FP = 4
      colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 1);
      valPtrM[0] = 1.0;

      nnz = 1;
      PM->insertGlobalValues((*fineElementMDOFs)(fineElement[1], 1), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0];

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[1], 2), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0] + 1;

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[1], 3), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      // if (CoarseElement is on the bottom left corner)
      if (coarseElement == 0) {
        // fill in the bottom left corner
        //  FP = 0
        colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 0);
        valPtrM[0] = 1.0;

        nnz = 1;
        PM->insertGlobalValues((*fineElementMDOFs)(fineElement[0], 0), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

        colPtrV[0] = 2 * colPtrM[0];

        PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 0), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

        colPtrV[0] = 2 * colPtrM[0] + 1;

        PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 1), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

        // FPr = 0
        colPtrP[0] = (*coarseElementPDOFs)(coarseElement, 0);
        valPtrP[0] = 1.0;

        nnz = 1;
        PP->insertGlobalValues((*fineElementPDOFs)(fineElement[0], 0), colPtrP.view(0, nnz), valPtrP.view(0, nnz));

      }  // if (coarseElement is on the bottom left corner)
    }    // if (coarseElement is on the bottom edge)

    // if (CoarseElement is on left edge)
    if (coarseElement % (nCoarseElements) == 0) {
      // fill in the left edge of the element patch
      //  FP = 5
      colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 0);
      valPtrM[0] = 0.375;
      colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 3);
      valPtrM[1] = -0.125;
      colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 7);
      valPtrM[2] = 0.75;

      nnz = 3;
      PM->insertGlobalValues((*fineElementMDOFs)(fineElement[0], 7), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0];
      colPtrV[1] = 2 * colPtrM[1];
      colPtrV[2] = 2 * colPtrM[2];

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 14), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0] + 1;
      colPtrV[1] = 2 * colPtrM[1] + 1;
      colPtrV[2] = 2 * colPtrM[2] + 1;

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 15), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      // FP = 6
      colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 7);
      valPtrM[0] = 1.0;

      nnz = 1;
      PM->insertGlobalValues((*fineElementMDOFs)(fineElement[0], 3), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0];

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 6), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0] + 1;

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 7), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      // FP = 7
      colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 0);
      valPtrM[0] = -0.125;
      colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 3);
      valPtrM[1] = 0.375;
      colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 7);
      valPtrM[2] = 0.75;

      nnz = 3;
      PM->insertGlobalValues((*fineElementMDOFs)(fineElement[2], 7), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0];
      colPtrV[1] = 2 * colPtrM[1];
      colPtrV[2] = 2 * colPtrM[2];

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[2], 14), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0] + 1;
      colPtrV[1] = 2 * colPtrM[1] + 1;
      colPtrV[2] = 2 * colPtrM[2] + 1;

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[2], 15), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      // FP = 8
      colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 3);
      valPtrM[0] = 1.0;

      nnz = 1;
      PM->insertGlobalValues((*fineElementMDOFs)(fineElement[2], 3), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0];

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[2], 6), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      colPtrV[0] = 2 * colPtrM[0] + 1;

      PV->insertGlobalValues((*fineElementVDOFs)(fineElement[2], 7), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

      // FPr = 3
      colPtrP[0] = (*coarseElementPDOFs)(coarseElement, 0);
      valPtrP[0] = 0.5;
      colPtrP[1] = (*coarseElementPDOFs)(coarseElement, 3);
      valPtrP[1] = 0.5;

      nnz = 2;
      PP->insertGlobalValues((*fineElementPDOFs)(fineElement[0], 3), colPtrP.view(0, nnz), valPtrP.view(0, nnz));

      // FPr = 4
      colPtrP[0] = (*coarseElementPDOFs)(coarseElement, 3);
      valPtrP[0] = 1.0;

      nnz = 1;
      PP->insertGlobalValues((*fineElementPDOFs)(fineElement[2], 3), colPtrP.view(0, nnz), valPtrP.view(0, nnz));

    }  // endif (coarseElement is on left edge)

    // fill in the rest of the patch
    //  FP = 9
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 8);
    valPtrM[0] = 1.0;

    nnz = 1;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[3], 0), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 0), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 1), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 10
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 5);
    valPtrM[0] = 1.0;

    nnz = 1;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[3], 1), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 2), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 3), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 11
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 2);
    valPtrM[0] = 1.0;

    nnz = 1;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[3], 2), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 4), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 5), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 12
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 6);
    valPtrM[0] = 1.0;

    nnz = 1;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[3], 3), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 6), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 7), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 13
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 4);
    valPtrM[0] = 0.375;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 6);
    valPtrM[1] = -0.125;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 8);
    valPtrM[2] = 0.75;

    nnz = 3;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[0], 5), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 10), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 11), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 14
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 5);
    valPtrM[0] = -0.125;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 7);
    valPtrM[1] = 0.375;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 8);
    valPtrM[2] = 0.75;

    nnz = 3;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[0], 6), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 12), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 13), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 15
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 1);
    valPtrM[0] = 0.375;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 2);
    valPtrM[1] = -0.125;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 5);
    valPtrM[2] = 0.75;

    nnz = 3;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[1], 5), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[1], 10), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[1], 11), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 16
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 5);
    valPtrM[0] = 0.375;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 7);
    valPtrM[1] = -0.125;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 8);
    valPtrM[2] = 0.75;

    nnz = 3;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[1], 6), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[1], 12), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[1], 13), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 17
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 4);
    valPtrM[0] = -0.125;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 6);
    valPtrM[1] = 0.375;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 8);
    valPtrM[2] = 0.75;

    nnz = 3;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[2], 5), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[2], 10), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[2], 11), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 18
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 2);
    valPtrM[0] = -0.125;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 3);
    valPtrM[1] = 0.375;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 6);
    valPtrM[2] = 0.75;

    nnz = 3;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[2], 6), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[2], 12), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[2], 13), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 19
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 1);
    valPtrM[0] = -0.125;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 2);
    valPtrM[1] = 0.375;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 5);
    valPtrM[2] = 0.75;

    nnz = 3;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[3], 5), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 10), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 11), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 20
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 2);
    valPtrM[0] = 0.375;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 3);
    valPtrM[1] = -0.125;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 6);
    valPtrM[2] = 0.75;

    nnz = 3;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[3], 6), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 12), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 13), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 21
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 0);
    valPtrM[0] = 0.140625;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 1);
    valPtrM[1] = -0.046875;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 2);
    valPtrM[2] = 0.015625;
    colPtrM[3] = (*coarseElementMDOFs)(coarseElement, 3);
    valPtrM[3] = -0.046875;
    colPtrM[4] = (*coarseElementMDOFs)(coarseElement, 4);
    valPtrM[4] = 0.28125;
    colPtrM[5] = (*coarseElementMDOFs)(coarseElement, 5);
    valPtrM[5] = -0.09375;
    colPtrM[6] = (*coarseElementMDOFs)(coarseElement, 6);
    valPtrM[6] = -0.09375;
    colPtrM[7] = (*coarseElementMDOFs)(coarseElement, 7);
    valPtrM[7] = 0.28125;
    colPtrM[8] = (*coarseElementMDOFs)(coarseElement, 8);
    valPtrM[8] = 0.5625;

    nnz = 9;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[0], 8), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];
    colPtrV[3] = 2 * colPtrM[3];
    colPtrV[4] = 2 * colPtrM[4];
    colPtrV[5] = 2 * colPtrM[5];
    colPtrV[6] = 2 * colPtrM[6];
    colPtrV[7] = 2 * colPtrM[7];
    colPtrV[8] = 2 * colPtrM[8];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 16), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;
    colPtrV[3] = 2 * colPtrM[3] + 1;
    colPtrV[4] = 2 * colPtrM[4] + 1;
    colPtrV[5] = 2 * colPtrM[5] + 1;
    colPtrV[6] = 2 * colPtrM[6] + 1;
    colPtrV[7] = 2 * colPtrM[7] + 1;
    colPtrV[8] = 2 * colPtrM[8] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[0], 17), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 22
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 0);
    valPtrM[0] = -0.046875;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 1);
    valPtrM[1] = 0.140625;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 2);
    valPtrM[2] = -0.046875;
    colPtrM[3] = (*coarseElementMDOFs)(coarseElement, 3);
    valPtrM[3] = 0.015625;
    colPtrM[4] = (*coarseElementMDOFs)(coarseElement, 4);
    valPtrM[4] = 0.28125;
    colPtrM[5] = (*coarseElementMDOFs)(coarseElement, 5);
    valPtrM[5] = 0.28125;
    colPtrM[6] = (*coarseElementMDOFs)(coarseElement, 6);
    valPtrM[6] = -0.09375;
    colPtrM[7] = (*coarseElementMDOFs)(coarseElement, 7);
    valPtrM[7] = -0.09375;
    colPtrM[8] = (*coarseElementMDOFs)(coarseElement, 8);
    valPtrM[8] = 0.5625;

    nnz = 9;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[1], 8), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];
    colPtrV[3] = 2 * colPtrM[3];
    colPtrV[4] = 2 * colPtrM[4];
    colPtrV[5] = 2 * colPtrM[5];
    colPtrV[6] = 2 * colPtrM[6];
    colPtrV[7] = 2 * colPtrM[7];
    colPtrV[8] = 2 * colPtrM[8];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[1], 16), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;
    colPtrV[3] = 2 * colPtrM[3] + 1;
    colPtrV[4] = 2 * colPtrM[4] + 1;
    colPtrV[5] = 2 * colPtrM[5] + 1;
    colPtrV[6] = 2 * colPtrM[6] + 1;
    colPtrV[7] = 2 * colPtrM[7] + 1;
    colPtrV[8] = 2 * colPtrM[8] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[1], 17), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 23
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 0);
    valPtrM[0] = -0.046875;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 1);
    valPtrM[1] = 0.015625;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 2);
    valPtrM[2] = -0.046875;
    colPtrM[3] = (*coarseElementMDOFs)(coarseElement, 3);
    valPtrM[3] = 0.140625;
    colPtrM[4] = (*coarseElementMDOFs)(coarseElement, 4);
    valPtrM[4] = -0.09375;
    colPtrM[5] = (*coarseElementMDOFs)(coarseElement, 5);
    valPtrM[5] = -0.09375;
    colPtrM[6] = (*coarseElementMDOFs)(coarseElement, 6);
    valPtrM[6] = 0.28125;
    colPtrM[7] = (*coarseElementMDOFs)(coarseElement, 7);
    valPtrM[7] = 0.28125;
    colPtrM[8] = (*coarseElementMDOFs)(coarseElement, 8);
    valPtrM[8] = 0.5625;

    nnz = 9;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[2], 8), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];
    colPtrV[3] = 2 * colPtrM[3];
    colPtrV[4] = 2 * colPtrM[4];
    colPtrV[5] = 2 * colPtrM[5];
    colPtrV[6] = 2 * colPtrM[6];
    colPtrV[7] = 2 * colPtrM[7];
    colPtrV[8] = 2 * colPtrM[8];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[2], 16), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;
    colPtrV[3] = 2 * colPtrM[3] + 1;
    colPtrV[4] = 2 * colPtrM[4] + 1;
    colPtrV[5] = 2 * colPtrM[5] + 1;
    colPtrV[6] = 2 * colPtrM[6] + 1;
    colPtrV[7] = 2 * colPtrM[7] + 1;
    colPtrV[8] = 2 * colPtrM[8] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[2], 17), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FP = 24
    colPtrM[0] = (*coarseElementMDOFs)(coarseElement, 0);
    valPtrM[0] = 0.015625;
    colPtrM[1] = (*coarseElementMDOFs)(coarseElement, 1);
    valPtrM[1] = -0.046875;
    colPtrM[2] = (*coarseElementMDOFs)(coarseElement, 2);
    valPtrM[2] = 0.140625;
    colPtrM[3] = (*coarseElementMDOFs)(coarseElement, 3);
    valPtrM[3] = -0.046875;
    colPtrM[4] = (*coarseElementMDOFs)(coarseElement, 4);
    valPtrM[4] = -0.09375;
    colPtrM[5] = (*coarseElementMDOFs)(coarseElement, 5);
    valPtrM[5] = 0.28125;
    colPtrM[6] = (*coarseElementMDOFs)(coarseElement, 6);
    valPtrM[6] = 0.28125;
    colPtrM[7] = (*coarseElementMDOFs)(coarseElement, 7);
    valPtrM[7] = -0.09375;
    colPtrM[8] = (*coarseElementMDOFs)(coarseElement, 8);
    valPtrM[8] = 0.5625;

    nnz = 9;
    PM->insertGlobalValues((*fineElementMDOFs)(fineElement[3], 8), colPtrM.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0];
    colPtrV[1] = 2 * colPtrM[1];
    colPtrV[2] = 2 * colPtrM[2];
    colPtrV[3] = 2 * colPtrM[3];
    colPtrV[4] = 2 * colPtrM[4];
    colPtrV[5] = 2 * colPtrM[5];
    colPtrV[6] = 2 * colPtrM[6];
    colPtrV[7] = 2 * colPtrM[7];
    colPtrV[8] = 2 * colPtrM[8];

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 16), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    colPtrV[0] = 2 * colPtrM[0] + 1;
    colPtrV[1] = 2 * colPtrM[1] + 1;
    colPtrV[2] = 2 * colPtrM[2] + 1;
    colPtrV[3] = 2 * colPtrM[3] + 1;
    colPtrV[4] = 2 * colPtrM[4] + 1;
    colPtrV[5] = 2 * colPtrM[5] + 1;
    colPtrV[6] = 2 * colPtrM[6] + 1;
    colPtrV[7] = 2 * colPtrM[7] + 1;
    colPtrV[8] = 2 * colPtrM[8] + 1;

    PV->insertGlobalValues((*fineElementVDOFs)(fineElement[3], 17), colPtrV.view(0, nnz), valPtrM.view(0, nnz));

    // FPr = 5
    colPtrP[0] = (*coarseElementPDOFs)(coarseElement, 0);
    valPtrP[0] = 0.25;
    colPtrP[1] = (*coarseElementPDOFs)(coarseElement, 1);
    valPtrP[1] = 0.25;
    colPtrP[2] = (*coarseElementPDOFs)(coarseElement, 2);
    valPtrP[2] = 0.25;
    colPtrP[3] = (*coarseElementPDOFs)(coarseElement, 3);
    valPtrP[3] = 0.25;

    nnz = 4;
    PP->insertGlobalValues((*fineElementPDOFs)(fineElement[0], 2), colPtrP.view(0, nnz), valPtrP.view(0, nnz));

    // FPr = 6
    colPtrP[0] = (*coarseElementPDOFs)(coarseElement, 1);
    valPtrP[0] = 0.5;
    colPtrP[1] = (*coarseElementPDOFs)(coarseElement, 2);
    valPtrP[1] = 0.5;

    nnz = 2;
    PP->insertGlobalValues((*fineElementPDOFs)(fineElement[1], 2), colPtrP.view(0, nnz), valPtrP.view(0, nnz));

    // FPr = 7
    colPtrP[0] = (*coarseElementPDOFs)(coarseElement, 2);
    valPtrP[0] = 1.0;

    nnz = 1;
    PP->insertGlobalValues((*fineElementPDOFs)(fineElement[3], 2), colPtrP.view(0, nnz), valPtrP.view(0, nnz));

    // FPr = 8
    colPtrP[0] = (*coarseElementPDOFs)(coarseElement, 2);
    valPtrP[0] = 0.5;
    colPtrP[1] = (*coarseElementPDOFs)(coarseElement, 3);
    valPtrP[1] = 0.5;

    nnz = 2;
    PP->insertGlobalValues((*fineElementPDOFs)(fineElement[3], 3), colPtrP.view(0, nnz), valPtrP.view(0, nnz));

    // Update counters:
    if ((coarseElement + 1) % (nCoarseElements) == 0)  // if the end of a row of c.g. elements
    {
      fineElement[0] = fineElement[3] + 1;
      fineElement[1] = fineElement[0] + 1;
      fineElement[2] = fineElement[0] + nFineElements;
      fineElement[3] = fineElement[2] + 1;
    } else {
      fineElement[0] = fineElement[1] + 1;
      fineElement[1] = fineElement[0] + 1;
      fineElement[2] = fineElement[3] + 1;
      fineElement[3] = fineElement[2] + 1;
    }
  }  // END OF BUILD LOOP

  // Loop over V rows
  for (GO VRow = 0; VRow < fNV; VRow++) {
    Teuchos::ArrayView<const LO> colPtr;
    Teuchos::ArrayView<const SC> valPtr;

    PV->getGlobalRowView(VRow, colPtr, valPtr);

    // Can be directly inserted!
    P->insertGlobalValues(VRow, colPtr, valPtr);
  }

  // Loop over P rows
  for (GO PRow = 0; PRow < fNP; PRow++) {
    Teuchos::ArrayView<const LO> colPtr;
    Teuchos::ArrayView<const SC> valPtr;

    // Now do pressure column:
    PP->getGlobalRowView(PRow, colPtr, valPtr);

    Teuchos::ArrayRCP<LO> newColPtr(colPtr.size(), nV);
    for (LO jj = 0; jj < colPtr.size(); jj++) {
      newColPtr[jj] += colPtr[jj];
    }

    // Insert into A
    P->insertGlobalValues(PRow + fNV, newColPtr.view(0, colPtr.size()), valPtr);
  }

  // Loop over M rows
  for (GO MRow = 0; MRow < fNM; MRow++) {
    Teuchos::ArrayView<const LO> colPtr;
    Teuchos::ArrayView<const SC> valPtr;

    // Now do magnetics column:
    PM->getGlobalRowView(MRow, colPtr, valPtr);

    Teuchos::ArrayRCP<LO> newColPtr(colPtr.size(), nV + nP);
    for (LO jj = 0; jj < colPtr.size(); jj++) {
      newColPtr[jj] += colPtr[jj];
    }

    // Insert into A
    P->insertGlobalValues(MRow + fNV + fNP, newColPtr.view(0, colPtr.size()), valPtr);
  }

  // Fill-complete all matrices
  PV->fillComplete(colMapforPV, rowMapforPV);
  PP->fillComplete(colMapforPP, rowMapforPP);
  PM->fillComplete(colMapforPM, rowMapforPM);
  P->fillComplete();

  // Set prolongators on the coarse grid
  Set(coarseLevel, "PV", PV);
  Set(coarseLevel, "PP", PP);
  Set(coarseLevel, "PM", PM);
  Set(coarseLevel, "P", P);

  Set(coarseLevel, "NV", nV);
  Set(coarseLevel, "NP", nP);
  Set(coarseLevel, "NM", nM);

}  // end buildp

}  // namespace MueLu

#define MUELU_GEOINTERPFACTORY_SHORT
#endif  // MUELU_GEOINTERPFACTORY_DEF_HPP
