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
#ifndef MUELU_BLOCKEDRAPFACTORY_DEF_HPP
#define MUELU_BLOCKEDRAPFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_RAPFactory_decl.hpp"
#include "MueLu_BlockedRAPFactory_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BlockedRAPFactory()
    : checkAc_(false), repairZeroDiagonals_(false)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(coarseLevel, "R");
    Input(fineLevel,   "A");
    Input(coarseLevel, "P");

    // call DeclareInput of all user-given transfer factories
    for(std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it!=transferFacts_.end(); ++it) {
      (*it)->CallDeclareInput(coarseLevel);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {  //FIXME make fineLevel const!!
    typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsMatrixClass; // TODO move me

    FactoryMonitor m(*this, "Computing Ac (block)", coarseLevel);

    //
    // Inputs: R, A, P
    //

    RCP<Matrix> R = Get< RCP<Matrix> >(coarseLevel, "R");
    RCP<Matrix> A = Get< RCP<Matrix> >(fineLevel,   "A");
    RCP<Matrix> P = Get< RCP<Matrix> >(coarseLevel, "P");

    //
    // Dynamic casts
    //

    RCP<BlockedCrsMatrixClass> bR, bA, bP;

    try {
      /* using rcp_dynamic_cast with throw_on_fail = true */
      bR = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(R, true);
      bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(A, true);
      bP = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(P, true);
    } catch(std::bad_cast e) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::BadCast, "MueLu::BlockedRAPFactory::Build(): matrices R, A and P must be of type BlockedCrsMatrix. " << e.what());
    }

    /*Utils::Write( "A00.m", CrsMatrixWrap(bA->getMatrix(0,0)) );
    Utils::Write( "A11.m", CrsMatrixWrap(bA->getMatrix(1,1)) );
    Utils::Write( "A01.m", CrsMatrixWrap(bA->getMatrix(0,1)) );
    Utils::Write( "A10.m", CrsMatrixWrap(bA->getMatrix(1,0)) );

    Utils::Write( "P00.m", CrsMatrixWrap(bP->getMatrix(0,0)) );
    Utils::Write( "P11.m", CrsMatrixWrap(bP->getMatrix(1,1)) );*/

    //
    // Build Ac = RAP
    //

    // Triple matrix product for BlockedCrsMatrixClass
    TEUCHOS_TEST_FOR_EXCEPTION((bA->Cols() != bP->Rows()) || (bA->Rows() != bR->Cols()), Exceptions::BadCast, "MueLu::BlockedRAPFactory::Build(): block matrix dimensions do not match.");
    RCP<BlockedCrsMatrixClass> bAP = Utils::TwoMatrixMultiplyBlock(bA, false, bP,  false, true, true);
    RCP<BlockedCrsMatrixClass> bAc = Utils::TwoMatrixMultiplyBlock(bR, false, bAP, false, true, true);

    if (checkAc_)
      CheckMainDiagonal(bAc);

    GetOStream(Statistics0, 0) << Utils::PrintMatrixInfo(*bAc, "Ac (blocked)");

    Set<RCP <Matrix> >(coarseLevel, "A", bAc);

    if (transferFacts_.begin() != transferFacts_.end()) {
      SubFactoryMonitor m1(*this, "Projections", coarseLevel);

      // call Build of all user-given transfer factories
      for(std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
        GetOStream(Runtime0, 0) << "Ac: call transfer factory " << (*it).get() << ": " << (*it)->description() << std::endl;
        (*it)->CallBuild(coarseLevel);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CheckMainDiagonal(RCP<BlockedCrsMatrix> & bAc, bool repairZeroDiagonals) {
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > c00 = bAc->getMatrix(0, 0);
    Teuchos::RCP<CrsMatrix> Aout = Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(c00->getRowMap(), c00->getGlobalMaxNumRowEntries(), Xpetra::StaticProfile);

    RCP<Vector> diagVec = VectorFactory::Build(c00->getRowMap());
    c00->getLocalDiagCopy(*diagVec);
    Teuchos::ArrayRCP< Scalar > diagVal = diagVec->getDataNonConst(0);

    // loop over local rows
    for(size_t row=0; row<c00->getNodeNumRows(); row++) {
      // get global row id
      GlobalOrdinal grid = c00->getRowMap()->getGlobalElement(row); // global row id

      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      c00->getLocalRowView(row, indices, vals);

      // just copy all values in output
      Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(), Teuchos::ScalarTraits<GlobalOrdinal>::zero());
      Teuchos::ArrayRCP<Scalar>        valout(indices.size(), Teuchos::ScalarTraits<Scalar>::zero());

      // just copy values
      for(size_t i=0; i<(size_t)indices.size(); i++) {
        GlobalOrdinal gcid = c00->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)
        indout [i] = gcid;
        valout [i] = vals[i];
      }

      Aout->insertGlobalValues(grid, indout.view(0, indout.size()), valout.view(0, valout.size()));
      if(diagVal[row]==0.0 && repairZeroDiagonals) {
        // always overwrite diagonal entry
        Aout->insertGlobalValues(grid, Teuchos::tuple<GlobalOrdinal>(grid), Teuchos::tuple<Scalar>(1.0));
      }
    }

    Aout->fillComplete(c00->getDomainMap(), c00->getRangeMap());

    bAc->setMatrix(0, 0, Aout);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddTransferFactory(const RCP<const FactoryBase>& factory) {
    // check if it's a TwoLevelFactoryBase based transfer factory
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast, "MueLu::RAPFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. This is very strange. (Note: you can remove this exception if there's a good reason for)");
    transferFacts_.push_back(factory);
  }

} //namespace MueLu

#define MUELU_BLOCKEDRAPFACTORY_SHORT
#endif // MUELU_BLOCKEDRAPFACTORY_DEF_HPP

// TODO add plausibility check
// TODO add CheckMainDiagonal for Blocked operator
// Avoid copying block matrix!
// create new empty Operator
