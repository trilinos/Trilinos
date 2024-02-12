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

#ifndef MUELU_MATRIXANALYSISFACTORY_DEF_HPP_
#define MUELU_MATRIXANALYSISFACTORY_DEF_HPP_

#include "MueLu_MatrixAnalysisFactory_decl.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_StridedMap.hpp>  // for nDofsPerNode...
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MatrixAnalysisFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MatrixAnalysisFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MatrixAnalysisFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~MatrixAnalysisFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> MatrixAnalysisFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A to be permuted.");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixAnalysisFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level & /* coarseLevel */) const {
  Input(fineLevel, "A");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixAnalysisFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level & /* coarseLevel */) const {
  FactoryMonitor m(*this, "MatrixAnalysis Factory ", fineLevel);

  Teuchos::RCP<Matrix> A                       = Get<Teuchos::RCP<Matrix> >(fineLevel, "A");
  Teuchos::RCP<const Teuchos::Comm<int> > comm = A->getRangeMap()->getComm();

  // const ParameterList & pL = GetParameterList();

  // General information
  {
    GetOStream(Runtime0) << "~~~~~~~~ GENERAL INFORMATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (Level " << fineLevel.GetLevelID() << ")" << std::endl;
    GetOStream(Runtime0) << "A is a " << A->getRangeMap()->getGlobalNumElements() << " x " << A->getDomainMap()->getGlobalNumElements() << " matrix." << std::endl;

    if (A->IsView("stridedMaps") && Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
      Teuchos::RCP<const StridedMap> strRowMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps"));
      LocalOrdinal blockid                     = strRowMap->getStridedBlockId();
      if (blockid > -1) {
        std::vector<size_t> stridingInfo = strRowMap->getStridingData();
        LO dofsPerNode                   = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);
        GetOStream(Runtime0) << "Strided row: " << dofsPerNode << " dofs per node. Strided block id = " << blockid << std::endl;
      } else {
        GetOStream(Runtime0) << "Row: " << strRowMap->getFixedBlockSize() << " dofs per node" << std::endl;
      }
    }

    if (A->IsView("stridedMaps") && Teuchos::rcp_dynamic_cast<const StridedMap>(A->getColMap("stridedMaps")) != Teuchos::null) {
      Teuchos::RCP<const StridedMap> strDomMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getColMap("stridedMaps"));
      LocalOrdinal blockid                     = strDomMap->getStridedBlockId();
      if (blockid > -1) {
        std::vector<size_t> stridingInfo = strDomMap->getStridingData();
        LO dofsPerNode                   = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);
        GetOStream(Runtime0) << "Strided column: " << dofsPerNode << " dofs per node. Strided block id = " << blockid << std::endl;
      } else {
        GetOStream(Runtime0) << "Column: " << strDomMap->getFixedBlockSize() << " dofs per node" << std::endl;
      }
    }

    GetOStream(Runtime0) << "A is distributed over " << comm->getSize() << " processors" << std::endl;

    // empty processors
    std::vector<LO> lelePerProc(comm->getSize(), 0);
    std::vector<LO> gelePerProc(comm->getSize(), 0);
    lelePerProc[comm->getRank()] = A->getLocalNumEntries();
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, comm->getSize(), &lelePerProc[0], &gelePerProc[0]);
    if (comm->getRank() == 0) {
      for (int i = 0; i < comm->getSize(); i++) {
        if (gelePerProc[i] == 0) {
          GetOStream(Runtime0) << "Proc " << i << " is empty." << std::endl;
        }
      }
    }
  }

  // Matrix diagonal
  {
    GetOStream(Runtime0) << "~~~~~~~~ MATRIX DIAGONAL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (Level " << fineLevel.GetLevelID() << ")" << std::endl;
    Teuchos::RCP<Vector> diagAVec = VectorFactory::Build(A->getRowMap(), true);
    A->getLocalDiagCopy(*diagAVec);
    Teuchos::ArrayRCP<const Scalar> diagAVecData = diagAVec->getData(0);
    GlobalOrdinal lzerosOnDiagonal               = 0;
    GlobalOrdinal gzerosOnDiagonal               = 0;
    GlobalOrdinal lnearlyzerosOnDiagonal         = 0;
    GlobalOrdinal gnearlyzerosOnDiagonal         = 0;
    GlobalOrdinal lnansOnDiagonal                = 0;
    GlobalOrdinal gnansOnDiagonal                = 0;

    for (size_t i = 0; i < diagAVec->getMap()->getLocalNumElements(); ++i) {
      if (diagAVecData[i] == Teuchos::ScalarTraits<Scalar>::zero()) {
        lzerosOnDiagonal++;
      } else if (Teuchos::ScalarTraits<Scalar>::magnitude(diagAVecData[i]) < 1e-6) {
        lnearlyzerosOnDiagonal++;
      } else if (Teuchos::ScalarTraits<Scalar>::isnaninf(diagAVecData[i])) {
        lnansOnDiagonal++;
      }
    }

    // sum up all entries in multipleColRequests over all processors
    MueLu_sumAll(comm, lzerosOnDiagonal, gzerosOnDiagonal);
    MueLu_sumAll(comm, lnearlyzerosOnDiagonal, gnearlyzerosOnDiagonal);
    MueLu_sumAll(comm, lnansOnDiagonal, gnansOnDiagonal);

    if (gzerosOnDiagonal > 0) GetOStream(Runtime0) << "Found " << gzerosOnDiagonal << " zeros on diagonal of A" << std::endl;
    if (gnearlyzerosOnDiagonal > 0) GetOStream(Runtime0) << "Found " << gnearlyzerosOnDiagonal << " entries with absolute value < 1.0e-6 on diagonal of A" << std::endl;
    if (gnansOnDiagonal > 0) GetOStream(Runtime0) << "Found " << gnansOnDiagonal << " entries with NAN or INF on diagonal of A" << std::endl;
  }

  // Diagonal dominance?
  {
    GetOStream(Runtime0) << "~~~~~~~~ DIAGONAL DOMINANCE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (Level " << fineLevel.GetLevelID() << ")" << std::endl;
    // loop over all local rows in matrix A and keep diagonal entries if corresponding
    // matrix rows are not contained in permRowMap
    GlobalOrdinal lnumWeakDiagDomRows   = 0;
    GlobalOrdinal gnumWeakDiagDomRows   = 0;
    GlobalOrdinal lnumStrictDiagDomRows = 0;
    GlobalOrdinal gnumStrictDiagDomRows = 0;
    double worstRatio                   = 99999999.9;
    for (size_t row = 0; row < A->getRowMap()->getLocalNumElements(); row++) {
      GlobalOrdinal grow = A->getRowMap()->getGlobalElement(row);

      size_t nnz = A->getNumEntriesInLocalRow(row);

      // extract local row information from matrix
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      A->getLocalRowView(row, indices, vals);

      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::MatrixAnalysisFactory::Build: number of nonzeros not equal to number of indices? Error.");

      // find column entry with max absolute value
      double norm1    = 0.0;
      double normdiag = 0.0;
      for (size_t j = 0; j < Teuchos::as<size_t>(indices.size()); j++) {
        GO gcol = A->getColMap()->getGlobalElement(indices[j]);
        if (gcol == grow)
          normdiag = Teuchos::as<double>(Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]));
        else
          norm1 += Teuchos::as<double>(Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]));
      }

      if (normdiag >= norm1)
        lnumWeakDiagDomRows++;
      else if (normdiag > norm1)
        lnumStrictDiagDomRows++;

      if (norm1 != 0.0) {
        if (normdiag / norm1 < worstRatio) worstRatio = normdiag / norm1;
      }
    }

    MueLu_sumAll(comm, lnumWeakDiagDomRows, gnumWeakDiagDomRows);
    MueLu_sumAll(comm, lnumStrictDiagDomRows, gnumStrictDiagDomRows);

    GetOStream(Runtime0) << "A has " << gnumWeakDiagDomRows << "/" << A->getRangeMap()->getGlobalNumElements() << " weakly diagonal dominant rows. (" << Teuchos::as<Scalar>(gnumWeakDiagDomRows / A->getRangeMap()->getGlobalNumElements() * 100) << "%)" << std::endl;
    GetOStream(Runtime0) << "A has " << gnumStrictDiagDomRows << "/" << A->getRangeMap()->getGlobalNumElements() << " strictly diagonal dominant rows. (" << Teuchos::as<Scalar>(gnumStrictDiagDomRows / A->getRangeMap()->getGlobalNumElements() * 100) << "%)" << std::endl;

    double gworstRatio;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, comm->getSize(), &worstRatio, &gworstRatio);
    GetOStream(Runtime0) << "The minimum of the ratio of diagonal element/ sum of off-diagonal elements is " << gworstRatio << ". Values about 1.0 are ok." << std::endl;
  }

  SC one = Teuchos::ScalarTraits<Scalar>::one(), zero = Teuchos::ScalarTraits<Scalar>::zero();

  // multiply with one vector
  {
    GetOStream(Runtime0) << "~~~~~~~~ Av with one vector ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (Level " << fineLevel.GetLevelID() << ")" << std::endl;
    Teuchos::RCP<Vector> ones = VectorFactory::Build(A->getDomainMap(), 1);
    ones->putScalar(one);

    Teuchos::RCP<Vector> res1 = VectorFactory::Build(A->getRangeMap(), false);
    A->apply(*ones, *res1, Teuchos::NO_TRANS, one, zero);
    checkVectorEntries(res1, std::string("after applying the one vector to A"));
  }

  {
    GetOStream(Runtime0) << "~~~~~~~~ Av with random vector ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (Level " << fineLevel.GetLevelID() << ")" << std::endl;
    Teuchos::RCP<Vector> randvec = VectorFactory::Build(A->getDomainMap(), 1);
    randvec->randomize();

    Teuchos::RCP<Vector> resrand = VectorFactory::Build(A->getRangeMap(), false);
    A->apply(*randvec, *resrand, Teuchos::NO_TRANS, one, zero);
    checkVectorEntries(resrand, std::string("after applying random vector to A"));
  }

  // apply Jacobi sweep
  {
    GetOStream(Runtime0) << "~~~~~~~~ Damped Jacobi sweep (one vector) ~~~~~~~~~~~~~~~~~~~~~~~ (Level " << fineLevel.GetLevelID() << ")" << std::endl;
    Teuchos::RCP<Vector> ones = VectorFactory::Build(A->getDomainMap(), 1);
    ones->putScalar(one);

    Teuchos::RCP<Vector> res1 = VectorFactory::Build(A->getRangeMap(), false);
    A->apply(*ones, *res1, Teuchos::NO_TRANS, one, zero);
    checkVectorEntries(res1, std::string("after applying one vector to A"));

    Teuchos::RCP<Vector> invDiag = Utilities::GetMatrixDiagonalInverse(*A);
    checkVectorEntries(invDiag, std::string("in invDiag"));

    Teuchos::RCP<Vector> res2 = VectorFactory::Build(A->getRangeMap(), false);
    res2->elementWiseMultiply(0.8, *invDiag, *res1, 0.0);
    checkVectorEntries(res2, std::string("after scaling Av with invDiag (with v the one vector)"));
    res2->update(1.0, *ones, -1.0);

    checkVectorEntries(res2, std::string("after applying one damped Jacobi sweep (with one vector)"));
  }

  // apply Jacobi sweep
  {
    GetOStream(Runtime0) << "~~~~~~~~ Damped Jacobi sweep (random vector) ~~~~~~~~~~~~~~~~~~~~ (Level " << fineLevel.GetLevelID() << ")" << std::endl;
    Teuchos::RCP<Vector> ones = VectorFactory::Build(A->getDomainMap(), 1);
    ones->randomize();

    Teuchos::RCP<Vector> res1 = VectorFactory::Build(A->getRangeMap(), false);
    A->apply(*ones, *res1, Teuchos::NO_TRANS, one, zero);
    checkVectorEntries(res1, std::string("after applying a random vector to A"));

    Teuchos::RCP<Vector> invDiag = Utilities::GetMatrixDiagonalInverse(*A);
    checkVectorEntries(invDiag, std::string("in invDiag"));

    Teuchos::RCP<Vector> res2 = VectorFactory::Build(A->getRangeMap(), false);
    res2->elementWiseMultiply(0.8, *invDiag, *res1, 0.0);
    checkVectorEntries(res2, std::string("after scaling Av with invDiag (with v a random vector)"));
    res2->update(1.0, *ones, -1.0);

    checkVectorEntries(res2, std::string("after applying one damped Jacobi sweep (with v a random vector)"));
  }

  GetOStream(Runtime0) << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ (Level " << fineLevel.GetLevelID() << ")" << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixAnalysisFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::checkVectorEntries(const Teuchos::RCP<Vector> &vec, std::string info) const {
  SC zero                                      = Teuchos::ScalarTraits<Scalar>::zero();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = vec->getMap()->getComm();
  Teuchos::ArrayRCP<const Scalar> vecData      = vec->getData(0);
  GO lzeros                                    = 0;
  GO gzeros                                    = 0;
  GO lnans                                     = 0;
  GO gnans                                     = 0;

  for (size_t i = 0; i < vec->getMap()->getLocalNumElements(); ++i) {
    if (vecData[i] == zero) {
      lzeros++;
    } else if (Teuchos::ScalarTraits<Scalar>::isnaninf(vecData[i])) {
      lnans++;
    }
  }

  // sum up all entries in multipleColRequests over all processors
  MueLu_sumAll(comm, lzeros, gzeros);
  MueLu_sumAll(comm, lnans, gnans);

  if (gzeros > 0) GetOStream(Runtime0) << "Found " << gzeros << " zeros " << info << std::endl;
  if (gnans > 0) GetOStream(Runtime0) << "Found " << gnans << " entries " << info << std::endl;
}

}  // namespace MueLu

#endif /* MUELU_MATRIXANALYSISFACTORY_DEF_HPP_ */
