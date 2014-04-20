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
#ifndef MUELU_TENTATIVEPNEWFACTORY_DEF_HPP
#define MUELU_TENTATIVEPNEWFACTORY_DEF_HPP

// disable clang warnings
#ifdef __clang__
#pragma clang system_header
#endif

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_StridedMapFactory.hpp>

#include "MueLu_TentativePNewFactory_decl.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> TentativePNewFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Aggregates",         Teuchos::null, "Generating factory of the aggregates");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",          Teuchos::null, "Generating factory of the nullspace");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
    validParamList->set< RCP<const FactoryBase> >("CoarseMap",          Teuchos::null, "Generating factory of the coarse map");
    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePNewFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & fineLevel, Level & coarseLevel) const {
    Input(fineLevel, "A");
    Input(fineLevel, "Aggregates");
    Input(fineLevel, "Nullspace");
    Input(fineLevel, "UnAmalgamationInfo");
    Input(fineLevel, "CoarseMap");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePNewFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TentativePNewFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level& fineLevel, Level& coarseLevel) const {

    FactoryMonitor m(*this, "Build", coarseLevel);

    RCP<Matrix>           A             = Get< RCP<Matrix> >          (fineLevel, "A");
    RCP<Aggregates>       aggregates    = Get< RCP<Aggregates> >      (fineLevel, "Aggregates");
    RCP<AmalgamationInfo> amalgInfo     = Get< RCP<AmalgamationInfo> >(fineLevel, "UnAmalgamationInfo");
    RCP<MultiVector>      fineNullspace = Get< RCP<MultiVector> >     (fineLevel, "Nullspace");
    RCP<const Map>        coarseMap     = Get< RCP<const Map> >       (fineLevel, "CoarseMap");

    RCP<const Map> rowMap = A->getRowMap();
    RCP<const Map> colMap = A->getColMap();

    const GO     indexBase = rowMap->getIndexBase();
    const size_t numRows   = rowMap->getNodeNumElements();

    RCP<const Teuchos::Comm<int> > comm = rowMap->getComm();

    typedef Teuchos::ScalarTraits<SC> STS;
    typedef typename STS::magnitudeType Magnitude;
    const SC     zero    = STS::zero();
    const SC     one     = STS::one();
    const LO     INVALID = Teuchos::OrdinalTraits<LO>::invalid();

    const GO     numAggs   = aggregates->GetNumAggregates();
    const bool   localAggs = !aggregates->AggregatesCrossProcessors();
    const size_t NSDim     = fineNullspace->getNumVectors();

    // Create a lookup table to determine the rows (fine DOFs) that belong to a given aggregate.
    // aggStart is a pointer into aggToRowMap
    // aggStart[i]..aggStart[i+1] are indices into aggToRowMap
    // aggToRowMap[aggStart[i]]..aggToRowMap[aggStart[i+1]-1] are the DOFs in aggregate i
    ArrayRCP<LO> aggStart;
    ArrayRCP<GO> aggToRowMap;
    amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);

    // Find largest aggregate size
    LO maxAggSize = 0;
    for (GO i = 0; i < numAggs; i++)
      maxAggSize = std::max(maxAggSize, aggStart[i+1] - aggStart[i]);

    RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(coarseMap, NSDim);

    // Pull out the nullspace vectors so that we can have random access.
    ArrayRCP<ArrayRCP<const SC> > fineNS  (NSDim);
    ArrayRCP<ArrayRCP<SC> >       coarseNS(NSDim);
    for (size_t i = 0; i < NSDim; i++) {
      fineNS[i] = fineNullspace->getData(i);
      if (coarseMap->getNodeNumElements() > 0)
        coarseNS[i] = coarseNullspace->getDataNonConst(i);
    }

    size_t nnzEstimate = numRows * NSDim;

    // Time to construct the matrix and fill in the values
    RCP<Matrix>    Ptentative = rcp(new CrsMatrixWrap(rowMap, coarseMap, 0, Xpetra::StaticProfile));
    RCP<CrsMatrix> PtentCrs   = rcp_dynamic_cast<CrsMatrixWrap>(Ptentative)->getCrsMatrix();

    ArrayRCP<size_t>  iaPtent;
    ArrayRCP<LO>      jaPtent;
    ArrayRCP<SC>     valPtent;

    PtentCrs->allocateAllValues(nnzEstimate, iaPtent, jaPtent, valPtent);

    ArrayView<size_t> ia  = iaPtent();
    ArrayView<LO>     ja  = jaPtent();
    ArrayView<SC>     val = valPtent();

    ia[0] = 0;
    for (size_t i = 1; i <= numRows; i++)
      ia[i] = ia[i-1] + NSDim;

    for (size_t j = 0; j < nnzEstimate; j++) {
      ja [j] = INVALID;
      val[j] = zero;
    }

    for (GO agg = 0; agg < numAggs; agg++) {
      LO aggSize = aggStart[agg+1] - aggStart[agg];

      Xpetra::global_size_t offset = agg*NSDim;

      // Extract the piece of the nullspace corresponding to the aggregate, and
      // put it in the flat array, "localQR" (in column major format) for the
      // QR routine.
      Teuchos::SerialDenseMatrix<LO,SC> localQR(aggSize, NSDim);
      for (size_t j = 0; j < NSDim; j++) {
        bool bIsZeroNSColumn = true;

        for (LO k = 0; k < aggSize; k++) {
          localQR(k,j) = fineNS[j][rowMap->getLocalElement(aggToRowMap[aggStart[agg]+k])];
          if (localQR(k,j) != zero)
            bIsZeroNSColumn = false;
        }
        TEUCHOS_TEST_FOR_EXCEPTION(bIsZeroNSColumn == true, Exceptions::RuntimeError,
                                   "MueLu::TentativePNewFactory::MakeTentative: fine level NS part has a zero column");
      }

      // Calculate QR decomposition (standard)
      // NOTE: Q is stored in localQR and R is stored in coarseNS
      if (aggSize >= Teuchos::as<LO>(NSDim)) {

        if (NSDim == 1) {
          // Only one nullspace vector, calculate Q and R by hand
          Magnitude norm = zero;
          for (size_t k = 0; k < Teuchos::as<size_t>(aggSize); k++)
            norm += STS::magnitude(localQR(k,0)*localQR(k,0));
          norm = Teuchos::ScalarTraits<Magnitude>::squareroot(norm);

          // R = norm
          coarseNS[0][offset] = norm;

          // Q = localQR(:,0)/norm
          for (LO i = 0; i < aggSize; i++)
            localQR(i,0) /= norm;

        } else {
          Teuchos::SerialQRDenseSolver<LO,SC> qrSolver;
          qrSolver.setMatrix(Teuchos::rcp(&localQR, false));
          qrSolver.factor();

          // R = upper triangular part of localQR
          for (size_t j = 0; j < NSDim; j++)
            for (size_t k = 0; k <= j; k++)
              coarseNS[j][offset+k] = localQR(k,j); //TODO is offset+k the correct local ID?!

          // Calculate Q, the tentative prolongator.
          // The Lapack GEQRF call only works for myAggsize >= NSDim
          qrSolver.formQ();
          Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SC> > qFactor = qrSolver.getQ();
          for (size_t j = 0; j < NSDim; j++)
            for (size_t i = 0; i < Teuchos::as<size_t>(aggSize); i++)
              localQR(i,j) = (*qFactor)(i,j);
        }

      } else {
        // Special handling for aggSize < NSDim (i.e. 1pt nodes)
        throw Exceptions::RuntimeError("aggSize < NSDim is not implemented");
      }

      // Process each row in the local Q factor
      // FIXME: What happens if maps are blocked?
      for (LO j = 0; j < aggSize; j++) {
        GO globalRow = aggToRowMap[aggStart[agg]+j];
        LO localRow  = rowMap->getLocalElement(globalRow); // CMS: There has to be an efficient way to do this...

        size_t rowStart = ia[localRow];
        for (size_t k = 0, lnnz = 0; k < NSDim; k++) {
          // Skip zeros (there may be plenty of them, i.e., NSDim > 1 or boundary conditions)
          if (localQR(j,k) != zero) {
            ja [rowStart+lnnz] = offset + k;
            val[rowStart+lnnz] = localQR(j,k);
            lnnz++;
          }
        }
      }
    }

    // Compress storage (remove all INVALID, which happen when we skip zeros)
    // We do that in-place
    size_t ia_tmp = 0, nnz = 0;
    for (size_t i = 0; i < numRows; i++) {
      for (size_t j = ia_tmp; j < ia[i+1]; j++)
        if (ja[j] != INVALID) {
          ja [nnz] = ja [j];
          val[nnz] = val[j];
          nnz++;
        }
      ia_tmp  = ia[i+1];
      ia[i+1] = nnz;
    }
    if (rowMap->lib() == Xpetra::UseTpetra) {
      // - Cannot resize for Epetra, as it checks for same pointers
      // - Need to resize for Tpetra, as it check ().size() == ia[numRows]
      // NOTE: these invalidate ja and val views
      jaPtent .resize(nnz);
      valPtent.resize(nnz);
    }

    GetOStream(Runtime1) << "TentativePNewFactory : aggregates do not cross process boundaries" << std::endl;

    {
      SubFactoryMonitor m2(*this, "fillComplete", coarseLevel);
      PtentCrs->setAllValues(iaPtent, jaPtent, valPtent);
      PtentCrs->expertStaticFillComplete(coarseMap, A->getDomainMap());
    }

    // If available, use striding information of fine level matrix A for range
    // map and coarseMap as domain map; otherwise use plain range map of
    // Ptent = plain range map of A for range map and coarseMap as domain map.
    // NOTE:
    // The latter is not really safe, since there is no striding information
    // for the range map. This is not really a problem, since striding
    // information is always available on the intermedium levels and the
    // coarsest levels.
    if (A->IsView("stridedMaps") == true)
      Ptentative->CreateView("stridedMaps", A->getRowMap("stridedMaps"), coarseMap);
    else
      Ptentative->CreateView("stridedMaps", Ptentative->getRangeMap(),   coarseMap);

    Set(coarseLevel, "Nullspace", coarseNullspace);
    Set(coarseLevel, "P",         Ptentative);

    if (IsPrint(Statistics1)) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      if (!localAggs)
        params->set("printCommInfo",        true);
      GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*Ptentative, "Ptent", params);
    }
  }

} //namespace MueLu

// TODO ReUse: If only P or Nullspace is missing, TentativePNewFactory can be smart and skip part of the computation.

#define MUELU_TENTATIVEPNEWFACTORY_SHORT
#endif // MUELU_TENTATIVEPNEWFACTORY_DEF_HPP
