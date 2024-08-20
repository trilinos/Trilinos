// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// A very simple driver to test basic MueLu functionality,
// most specifically, a geometric linear interpolation routine.

#include <iostream>
#include <fstream>
#include <cmath>

// MueLu main header: include most common header files in one line
#include <MueLu.hpp>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Epetra
#include <Epetra_CrsMatrix.h>

// EpetraExt
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#define XPETRA_ENABLED
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

#include "MueLu_Level.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>

#include <MueLu_GeoInterpFactory_def.hpp>
#include <MueLu_Q2Q1Q2CoarseGridFactory_def.hpp>
#include <MueLu_SaPFactory_def.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_MHDRAPFactory_def.hpp>

#include <MueLu_UseDefaultTypes.hpp>

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>

//#include "EditCopies/MueLu_MHDVankaSmoother_def.hpp"

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;  // reference count pointers

  typedef MueLu::GeoInterpFactory<SC, LO, GO, NO> GeoInterpFactory;
  typedef MueLu::Q2Q1Q2CoarseGridFactory<SC, LO, GO, NO> Q2Q1Q2CoarseGridFactory;
  typedef MueLu::MHDRAPFactory<SC, LO, GO, NO> MHDRAPFactory;

  //
  // MPI initialization using Teuchos
  // *Included because JHU indicated that it shouldn't cause problems
  // if running in serial
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // Initialize a "FancyOStream" to output to standard out (cout)
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);

    // First, we start with an Xpetra::Map
    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

    string fileNameA00 = "./Matrices/J_00_3.mm";
    string fileNameA01 = "./Matrices/J_01_3.mm";
    string fileNameA02 = "./Matrices/J_02_3.mm";
    string fileNameA10 = "./Matrices/J_10_3.mm";
    string fileNameA11 = "./Matrices/J_11_3.mm";
    string fileNameA12 = "./Matrices/J_12_3.mm";
    string fileNameA20 = "./Matrices/J_20_3.mm";
    string fileNameA21 = "./Matrices/J_21_3.mm";
    string fileNameA22 = "./Matrices/J_22_3.mm";

    RCP<Matrix> A00 = Utils::Read(fileNameA00, lib, comm);
    RCP<Matrix> A01 = Utils::Read(fileNameA01, lib, comm);
    RCP<Matrix> A02 = Utils::Read(fileNameA02, lib, comm);
    RCP<Matrix> A10 = Utils::Read(fileNameA10, lib, comm);
    RCP<Matrix> A11 = Utils::Read(fileNameA11, lib, comm);
    RCP<Matrix> A12 = Utils::Read(fileNameA12, lib, comm);
    RCP<Matrix> A20 = Utils::Read(fileNameA20, lib, comm);
    RCP<Matrix> A21 = Utils::Read(fileNameA21, lib, comm);
    RCP<Matrix> A22 = Utils::Read(fileNameA22, lib, comm);

    std::cout << A00->isFillComplete() << std::endl;

    // Now we want to build A;
    RCP<const Map> VMap = A00->getRowMap();
    RCP<const Map> PMap = A10->getRowMap();
    RCP<const Map> MMap = A20->getRowMap();

    size_t nv = VMap->getGlobalNumElements();
    size_t np = PMap->getGlobalNumElements();
    size_t nm = MMap->getGlobalNumElements();

    RCP<const Map> AMap = Xpetra::MapFactory<LO, GO>::createUniformContigMap(Xpetra::UseTpetra, nv + np + nm, comm);

    size_t maxEntriesPerRow = 84;
    RCP<Matrix> A           = rcp(new CrsMatrixWrap(AMap, maxEntriesPerRow));

    Teuchos::ArrayView<const LO> colPtr;
    Teuchos::ArrayView<const SC> valPtr;

    A00->getLocalRowView(7, colPtr, valPtr);

    std::cout << "valPtr.size() = " << valPtr.size() << std::endl;

    A01->getLocalRowView(7, colPtr, valPtr);

    std::cout << "valPtr.size() = " << valPtr.size() << std::endl;

    // Loop over V rows
    for (size_t VRow = 0; VRow < nv; VRow++) {
      Teuchos::ArrayView<const LO> colPtr;
      Teuchos::ArrayView<const SC> valPtr;

      A00->getLocalRowView(VRow, colPtr, valPtr);

      // Can be directly inserted!
      A->insertGlobalValues(VRow, colPtr, valPtr);

      // Now do pressure column:
      A01->getLocalRowView(VRow, colPtr, valPtr);

      Teuchos::ArrayRCP<LO> newColPtr(colPtr.size(), nv);
      for (LO jj = 0; jj < colPtr.size(); jj++) {
        newColPtr[jj] += colPtr[jj];
      }

      // Insert into A
      A->insertGlobalValues(VRow, newColPtr.view(0, colPtr.size()), valPtr);

      // Now do magnetics column:
      A02->getLocalRowView(VRow, colPtr, valPtr);

      newColPtr.clear();
      newColPtr.resize(colPtr.size(), nv + np);
      for (LO jj = 0; jj < colPtr.size(); jj++) {
        newColPtr[jj] += colPtr[jj];
      }

      // Insert into A
      A->insertGlobalValues(VRow, newColPtr.view(0, colPtr.size()), valPtr);
    }

    // Loop over P rows
    for (size_t PRow = 0; PRow < np; PRow++) {
      Teuchos::ArrayView<const LO> colPtr;
      Teuchos::ArrayView<const SC> valPtr;

      A10->getLocalRowView(PRow, colPtr, valPtr);

      // Can be directly inserted!
      A->insertGlobalValues(PRow + nv, colPtr, valPtr);

      // Now do pressure column:
      A11->getLocalRowView(PRow, colPtr, valPtr);

      Teuchos::ArrayRCP<LO> newColPtr(colPtr.size(), nv);
      for (LO jj = 0; jj < colPtr.size(); jj++) {
        newColPtr[jj] += colPtr[jj];
      }

      // Insert into A
      A->insertGlobalValues(PRow + nv, newColPtr.view(0, colPtr.size()), valPtr);

      // Now do magnetics column:
      A12->getLocalRowView(PRow, colPtr, valPtr);

      newColPtr.clear();
      newColPtr.resize(colPtr.size(), nv + np);
      for (LO jj = 0; jj < colPtr.size(); jj++) {
        newColPtr[jj] += colPtr[jj];
      }

      // Insert into A
      A->insertGlobalValues(PRow + nv, newColPtr.view(0, colPtr.size()), valPtr);
    }

    // Loop over M rows
    for (size_t MRow = 0; MRow < nm; MRow++) {
      Teuchos::ArrayView<const LO> colPtr;
      Teuchos::ArrayView<const SC> valPtr;

      A20->getLocalRowView(MRow, colPtr, valPtr);

      // Can be directly inserted!
      A->insertGlobalValues(MRow + nv + np, colPtr, valPtr);

      // Now do pressure column:
      A21->getLocalRowView(MRow, colPtr, valPtr);

      Teuchos::ArrayRCP<LO> newColPtr(colPtr.size(), nv);
      for (LO jj = 0; jj < colPtr.size(); jj++) {
        newColPtr[jj] += colPtr[jj];
      }

      // Insert into A
      A->insertGlobalValues(MRow + nv + np, newColPtr.view(0, colPtr.size()), valPtr);

      // Now do magnetics column:
      A22->getLocalRowView(MRow, colPtr, valPtr);

      newColPtr.clear();
      newColPtr.resize(colPtr.size(), nv + np);
      for (LO jj = 0; jj < colPtr.size(); jj++) {
        newColPtr[jj] += colPtr[jj];
      }

      // Insert into A
      A->insertGlobalValues(MRow + nv + np, newColPtr.view(0, colPtr.size()), valPtr);
    }

    A->fillComplete();

    // Let's read in the element connectivity info:
    GO totalFineElements = 32 * 32;

    RCP<Teuchos::SerialDenseMatrix<GO, GO> > fineGridVElements = rcp(new Teuchos::SerialDenseMatrix<GO, GO>(totalFineElements, 18));
    RCP<Teuchos::SerialDenseMatrix<GO, GO> > fineGridPElements = rcp(new Teuchos::SerialDenseMatrix<GO, GO>(totalFineElements, 4));
    RCP<Teuchos::SerialDenseMatrix<GO, GO> > fineGridMElements = rcp(new Teuchos::SerialDenseMatrix<GO, GO>(totalFineElements, 9));

    std::ifstream VElementFile("./Matrices/elements_0_0");
    std::ifstream PElementFile("./Matrices/elements_1_1");
    std::ifstream MElementFile("./Matrices/elements_2_2");

    for (GO ii = 0; ii < totalFineElements; ii++) {
      for (LO jj = 0; jj < 9; jj++) {
        VElementFile >> (*fineGridVElements)(ii, 2 * jj);
        VElementFile >> (*fineGridVElements)(ii, 2 * jj + 1);

        MElementFile >> (*fineGridMElements)(ii, jj);
      }

      for (LO kk = 0; kk < 4; kk++) {
        PElementFile >> (*fineGridPElements)(ii, kk);
      }
    }
    VElementFile.close();
    PElementFile.close();
    MElementFile.close();

    Hierarchy H;
    H.setDefaultVerbLevel(Teuchos::VERB_NONE);

    RCP<Level> finest = H.GetLevel();
    finest->setDefaultVerbLevel(Teuchos::VERB_NONE);

    finest->Set("A", A);
    finest->Set("A00", A00);
    finest->Set("A01", A01);
    finest->Set("A02", A02);
    finest->Set("A10", A10);
    finest->Set("A11", A11);
    finest->Set("A12", A12);
    finest->Set("A20", A20);
    finest->Set("A21", A21);
    finest->Set("A22", A22);

    // Set finegrid elements
    finest->Set("VElementList", fineGridVElements);
    finest->Set("PElementList", fineGridPElements);
    finest->Set("MElementList", fineGridMElements);

    // Set count stuff for Vanka Smoother
    finest->Set("NV", A00->getGlobalNumRows());
    finest->Set("NP", A10->getGlobalNumRows());
    finest->Set("NM", A20->getGlobalNumRows());

    // Create a GeoInterpFactory
    RCP<GeoInterpFactory> geoInterp                    = rcp(new GeoInterpFactory());
    RCP<Q2Q1Q2CoarseGridFactory> coarseElementFact     = rcp(new Q2Q1Q2CoarseGridFactory());
    RCP<MueLu::MHDRAPFactory<SC, LO, GO, NO> > rapFact = rcp(new MueLu::MHDRAPFactory<SC, LO, GO, NO>());

    RCP<FactoryManager> M = rcp(new FactoryManager());
    M->SetFactory("A", rapFact);
    M->SetFactory("A00", rapFact);
    M->SetFactory("A01", rapFact);
    M->SetFactory("A02", rapFact);
    M->SetFactory("A10", rapFact);
    M->SetFactory("A11", rapFact);
    M->SetFactory("A12", rapFact);
    M->SetFactory("A20", rapFact);
    M->SetFactory("A21", rapFact);
    M->SetFactory("A22", rapFact);
    M->SetFactory("VElementList", coarseElementFact);
    M->SetFactory("PElementList", coarseElementFact);
    M->SetFactory("MElementList", coarseElementFact);
    M->SetFactory("PV", geoInterp);
    M->SetFactory("PP", geoInterp);
    M->SetFactory("PM", geoInterp);
    M->SetFactory("P", geoInterp);

    H.Setup(*M, 0, 3);

    // Utils::Write("./output/BigAMat.mm",*A);
    std::cout << "Hello world!\n";

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
