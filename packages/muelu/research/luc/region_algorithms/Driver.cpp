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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <MueLu.hpp>
#include <MueLu_ParameterListInterpreter.hpp>  // TODO: move into MueLu.hpp
#include <MueLu_Utilities.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_CoarseMapFactory.hpp>
#include <MueLu_TentativePFactory.hpp>

#include <MueLu_CreateXpetraPreconditioner.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>   // => This header defines Belos::MueLuOp
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void createTwoLevelHierarchy(MueLu::Level& fineLevel, MueLu::Level& coarseLevel,
                             Teuchos::RCP<typename Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A);

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using Teuchos::tuple;

  typedef Tpetra::Map<LO, GO, NO> map_type;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> crs_matrix_type;

  typedef typename crs_matrix_type::scalar_type scalar_type;
  typedef typename crs_matrix_type::local_ordinal_type local_ordinal_type;
  typedef typename crs_matrix_type::global_ordinal_type global_ordinal_type;
  typedef typename crs_matrix_type::node_type node_type;
  typedef typename crs_matrix_type::mag_type magnitude_type;

  bool success = false;
  bool verbose = true;

  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    const local_ordinal_type myRank     = comm->getRank();

    // Manage the way output stream works
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fancy->setShowAllFrontMatter(false).setShowProcRank(true);
    Teuchos::FancyOStream& out = *fancy;
    // out.setOutputToRootOnly(0);
    // // Do not print anything
    // Teuchos::oblackholestream blackhole;
    // std::ostream& out = blackhole;

    // =========================================================================
    // Parameters initialization
    // =========================================================================

    // GO nx = 100, ny = 100, nz = 100;
    // Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);  // manage parameters of Xpetra

    std::string xmlFileName = "driver.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'scalingTest.xml'");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    global_ordinal_type gNumCompFineGIDs = 10, lCompFineGIDOffset = 0;
    global_ordinal_type lRegFineGIDOffset = 0;  // gNumRegFineGIDs  = 11,
    local_ordinal_type lNumCompFineGIDs = 0, lNumRegFineGIDs = 0;
    local_ordinal_type lNumRegCoarseGIDs = 0;  // lNumCompCoarseGIDs =  0,
    if (myRank == 0) {
      lNumCompFineGIDs = 7;
      lNumRegFineGIDs  = 7;

      lCompFineGIDOffset = 0;
      lRegFineGIDOffset  = 0;

      // lNumCompCoarseGIDs = 3;
      lNumRegCoarseGIDs = 3;
    } else if (myRank == 1) {
      lNumCompFineGIDs = 3;
      lNumRegFineGIDs  = 4;

      lCompFineGIDOffset = 7;
      lRegFineGIDOffset  = 6;

      // lNumCompCoarseGIDs = 1;
      lNumRegCoarseGIDs = 2;
    }

    // The initial focus is on getting rowMaps as they are also used as rangeMap and potentially
    // domainMaps...
    Array<global_ordinal_type> fineCompRowGIDs(lNumCompFineGIDs);
    for (local_ordinal_type dof = 0; dof < lNumCompFineGIDs; ++dof) {
      fineCompRowGIDs[dof] = lCompFineGIDOffset + dof;
    }
    out << "fineCompRowGIDs: " << fineCompRowGIDs << std::endl;
    RCP<map_type> fineCompRowMap = rcp(new map_type(10, fineCompRowGIDs, 0, comm));

    // Now the columnMap for the initial composite operator is needed
    Array<global_ordinal_type> fineCompColGIDs(lNumCompFineGIDs + 1);
    if (myRank == 0) {
      for (local_ordinal_type dof = 0; dof < lNumCompFineGIDs + 1; ++dof) {
        fineCompColGIDs[dof] = lCompFineGIDOffset + dof;
      }
    } else if (myRank == 1) {
      for (local_ordinal_type dof = 0; dof < lNumCompFineGIDs + 1; ++dof) {
        fineCompColGIDs[dof] = lCompFineGIDOffset - 1 + dof;
      }
    }
    out << "fineCompColGIDs: " << fineCompColGIDs << std::endl;
    RCP<map_type> fineCompColMap = rcp(new map_type(12, fineCompColGIDs, 0, comm));

    // Create the matrix using our maps and assuming at most 3 entries per row
    RCP<crs_matrix_type> compA(new crs_matrix_type(fineCompRowMap, fineCompColMap, 3));

    // Now each row needs to be filled
    const scalar_type one    = static_cast<scalar_type>(1.0);
    const scalar_type two    = static_cast<scalar_type>(2.0);
    const scalar_type negOne = static_cast<scalar_type>(-1.0);
    for (local_ordinal_type lclRow = 0; lclRow < static_cast<local_ordinal_type>(lNumCompFineGIDs); ++lclRow) {
      const global_ordinal_type gblRow = fineCompRowMap->getGlobalElement(lclRow);
      // A(0, 0) = [1]
      if (gblRow == 0) {
        compA->insertGlobalValues(gblRow,
                                  tuple<global_ordinal_type>(gblRow),
                                  tuple<scalar_type>(one));
      }
      // A(N-1, N-2:N-1) = [-1, 1]
      else if (gblRow == gNumCompFineGIDs - 1) {
        compA->insertGlobalValues(gblRow,
                                  tuple<global_ordinal_type>(gblRow - 1, gblRow),
                                  tuple<scalar_type>(negOne, one));
      }
      // A(i, i-1:i+1) = [-1, 2, -1]
      else {
        compA->insertGlobalValues(gblRow,
                                  tuple<global_ordinal_type>(gblRow - 1, gblRow, gblRow + 1),
                                  tuple<scalar_type>(negOne, two, negOne));
      }
    }
    compA->fillComplete(fineCompRowMap, fineCompRowMap);
    // compA->print(out);

    out << std::endl
        << "Now switching to the region matrix" << std::endl;

    Array<global_ordinal_type> fineRegRowGIDs(lNumRegFineGIDs);
    for (local_ordinal_type dof = 0; dof < lNumRegFineGIDs; ++dof) {
      fineRegRowGIDs[dof] = (lRegFineGIDOffset + dof == 6 && myRank == 1) ? 10 : lRegFineGIDOffset + dof;
    }
    out << "fineRegRowGIDs: " << fineRegRowGIDs << std::endl;
    RCP<map_type> fineRegRowMap = rcp(new map_type(11, fineRegRowGIDs, 0, comm));

    // Create the matrix using the simplest constructor as we assume that the region matrix
    // has now off processor entries, hence no column map is required.
    RCP<crs_matrix_type> regA(new crs_matrix_type(fineRegRowMap, 3));

    // Now each row needs to be filled
    for (local_ordinal_type lclRow = 0; lclRow < static_cast<local_ordinal_type>(lNumRegFineGIDs); ++lclRow) {
      const global_ordinal_type prevGblRow = (lclRow > 0) ? fineRegRowMap->getGlobalElement(lclRow - 1) : -1;
      const global_ordinal_type curGblRow  = fineRegRowMap->getGlobalElement(lclRow);
      const global_ordinal_type nextGblRow = (lclRow < lNumRegFineGIDs - 1) ? fineRegRowMap->getGlobalElement(lclRow + 1) : -1;
      // A(0, 0) = [1]
      if (lclRow == 0) {
        if (curGblRow == 0) {  // Apply Dirichlet BC to the left of the mesh
          regA->insertGlobalValues(curGblRow,
                                   tuple<global_ordinal_type>(curGblRow),
                                   tuple<scalar_type>(one));
        } else {  // Here I hard code the column associated with "-1" because I'm lazy...
          regA->insertGlobalValues(curGblRow,
                                   tuple<global_ordinal_type>(curGblRow, nextGblRow),
                                   tuple<scalar_type>(one, negOne));
        }
      }
      // A(N-1, N-2:N-1) = [-1, 1]
      else if (lclRow == static_cast<local_ordinal_type>(lNumRegFineGIDs - 1)) {
        regA->insertGlobalValues(curGblRow,
                                 tuple<global_ordinal_type>(prevGblRow, curGblRow),
                                 tuple<scalar_type>(negOne, one));
      }
      // A(i, i-1:i+1) = [-1, 2, -1]
      else {
        regA->insertGlobalValues(curGblRow,
                                 tuple<global_ordinal_type>(prevGblRow, curGblRow, nextGblRow),
                                 tuple<scalar_type>(negOne, two, negOne));
      }
    }
    regA->fillComplete();

    Tpetra::MatrixMarket::Writer<crs_matrix_type>::
        writeSparseFile("regA.m", regA, "regA",
                        "region representation of the operator.");

    out << std::endl
        << "Forming the prolongator" << std::endl;

    MueLu::Level fineLevel, coarseLevel;
    RCP<Xpetra::Matrix<SC, LO, GO, NO> > matA = MueLu::TpetraCrs_To_XpetraMatrix<SC, LO, GO, NO>(regA);
    createTwoLevelHierarchy<SC, LO, GO, NO>(fineLevel, coarseLevel, matA);

    // Now the prolongator needs to be created
    Array<global_ordinal_type> coarseRegGIDs(lNumRegCoarseGIDs);
    if (myRank == 0) {
      coarseRegGIDs[0] = 0;
      coarseRegGIDs[1] = 3;
      coarseRegGIDs[2] = 6;
    } else if (myRank == 1) {
      coarseRegGIDs[0] = 10;
      coarseRegGIDs[1] = 9;
    }
    out << "coarseRegGIDs: " << coarseRegGIDs << std::endl;
    RCP<map_type> colMapP = rcp(new map_type(5, coarseRegGIDs, 0, comm));
    // RCP<crs_matrix_type> regP (new crs_matrix_type (fineRegRowMap, 1));
    RCP<crs_matrix_type> regP(new crs_matrix_type(fineRegRowMap, colMapP, 1));

    // Now each row needs to be filled
    for (local_ordinal_type lclRow = 0; lclRow < static_cast<local_ordinal_type>(lNumRegFineGIDs); ++lclRow) {
      const global_ordinal_type gblRow = fineRegRowMap->getGlobalElement(lclRow);
      const global_ordinal_type gblCol = colMapP->getGlobalElement((lclRow + 1) / 3);
      regP->insertGlobalValues(gblRow,
                               tuple<global_ordinal_type>(gblCol),
                               tuple<scalar_type>(one));
    }
    regP->fillComplete();

    Tpetra::MatrixMarket::Writer<crs_matrix_type>::
        writeSparseFile("regP.m", regP, "regP",
                        "region representation of the prolongator.");

    out << std::endl
        << "Computing AP = AxP" << std::endl;

    // The number of non-zeros per row is set to zero because it's hard to guess what it will be...
    RCP<crs_matrix_type> regAP = rcp(new crs_matrix_type(fineRegRowMap, 0));
    Tpetra::MatrixMatrix::
        Multiply<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(*regA, false, *regP,
                                                                                  false, *regAP);

    Tpetra::MatrixMarket::Writer<crs_matrix_type>::
        writeSparseFile("regAP.m", regAP, "regAP",
                        "region representation of AP.");

    out << std::endl
        << "Computing Ac = P'xAP" << std::endl;

    // The number of non-zeros per row is set to zero because it's hard to guess what it will be...
    RCP<crs_matrix_type> regAc = rcp(new crs_matrix_type(regP->getDomainMap(), 0));
    Tpetra::MatrixMatrix::
        Multiply<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(*regP, true, *regAP,
                                                                                  false, *regAc);

    Tpetra::MatrixMarket::Writer<crs_matrix_type>::
        writeSparseFile("regAc.m", regAc, "regAc",
                        "region representation of Ac.");

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}  // main

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void createTwoLevelHierarchy(MueLu::Level& fineLevel, MueLu::Level& coarseLevel,
                             Teuchos::RCP<typename Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
  fineLevel.SetFactoryManager(factoryHandler);
  coarseLevel.SetFactoryManager(factoryHandler);

  coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));

  fineLevel.SetLevelID(0);
  coarseLevel.SetLevelID(1);

  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  A->SetFixedBlockSize(1);
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("DofsPerNode", 2);

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), 1);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<UncoupledAggregationFactory> UnCoupledAggFact = rcp(new UncoupledAggregationFactory());
  UnCoupledAggFact->SetFactory("Graph", dropFact);

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", UnCoupledAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
  TentativePFact->SetFactory("Aggregates", UnCoupledAggFact);
  TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  TentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", TentativePFact.get());  // request Ptent
  coarseLevel.Request("Nullspace", TentativePFact.get());
  coarseLevel.Request(*TentativePFact);
  Teuchos::ParameterList paramList;
  paramList.set("tentative: calculate qr", false);
  TentativePFact->SetParameterList(paramList);
  TentativePFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> Ptent;
  coarseLevel.Get("P", Ptent, TentativePFact.get());
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char* argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
