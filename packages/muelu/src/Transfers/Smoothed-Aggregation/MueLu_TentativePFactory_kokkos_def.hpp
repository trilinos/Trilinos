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
#ifndef MUELU_TENTATIVEPFACTORY_KOKKOS_DEF_HPP
#define MUELU_TENTATIVEPFACTORY_KOKKOS_DEF_HPP

#ifdef HAVE_MUELU_KOKKOS_REFACTOR

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialQRDenseSolver.hpp>

#include "MueLu_TentativePFactory_kokkos_decl.hpp"

#include "MueLu_Aggregates_kokkos.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_CoarseMapFactory_kokkos.hpp"
#include "MueLu_NullspaceFactory_kokkos.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities_kokkos.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const ParameterList> TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Aggregates",         Teuchos::null, "Generating factory of the aggregates");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",          Teuchos::null, "Generating factory of the nullspace");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
    validParamList->set< RCP<const FactoryBase> >("CoarseMap",          Teuchos::null, "Generating factory of the coarse map");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(fineLevel, "A");
    Input(fineLevel, "Aggregates");
    Input(fineLevel, "Nullspace");
    Input(fineLevel, "UnAmalgamationInfo");
    Input(fineLevel, "CoarseMap");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    auto A             = Get< RCP<Matrix> >           (fineLevel, "A");
    auto aggregates    = Get< RCP<Aggregates_kokkos> >(fineLevel, "Aggregates");
    auto amalgInfo     = Get< RCP<AmalgamationInfo> > (fineLevel, "UnAmalgamationInfo");
    auto fineNullspace = Get< RCP<MultiVector> >      (fineLevel, "Nullspace");
    auto coarseMap     = Get< RCP<const Map> >        (fineLevel, "CoarseMap");

    RCP<Matrix>      Ptentative;
    RCP<MultiVector> coarseNullspace;
    if (!aggregates->AggregatesCrossProcessors())
      BuildPuncoupled(A, aggregates, amalgInfo, fineNullspace, coarseMap, Ptentative, coarseNullspace);
    else
      BuildPcoupled  (A, aggregates, amalgInfo, fineNullspace, coarseMap, Ptentative, coarseNullspace);

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
      GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*Ptentative, "Ptent", params);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::
  BuildPuncoupled(RCP<Matrix> A, RCP<Aggregates_kokkos> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                  RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace) const {
    RCP<const Map> rowMap = A->getRowMap();
    RCP<const Map> colMap = A->getColMap();

    const size_t numRows  = rowMap->getNodeNumElements();

    typedef Teuchos::ScalarTraits<SC> STS;
    const SC zero    = STS::zero();
    const SC one     = STS::one();
    const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();

    auto     aggGraph = aggregates->GetGraph();
    auto     aggRows  = aggGraph.row_map;
    auto     aggCols  = aggGraph.entries;
    const GO numAggs  = aggregates->GetNumAggregates();

    // Aggregates map is based on the amalgamated column map
    // We can skip global-to-local conversion if LIDs in row map are
    // same as LIDs in column map
    bool goodMap = isGoodMap(*rowMap, *colMap);
#if 1
    TEUCHOS_TEST_FOR_EXCEPTION(!goodMap,                    Exceptions::RuntimeError, "For now, need matching row and col maps");
    TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() != 1, Exceptions::RuntimeError, "For now, only block size 1");

    // For now, do a simple translation
    // FIXME: only type is correct here, everything else is not
    Kokkos::View<LO*, DeviceType> aggToRowMapLO("agg2row_map", numAggs);
#else
    ArrayRCP<LO> aggStart;
    ArrayRCP<LO> aggToRowMapLO;
    ArrayRCP<GO> aggToRowMapGO;
    if (goodMap) {
      amalgInfo->UnamalgamateAggregatesLO(*aggregates, aggStart, aggToRowMapLO);
      GetOStream(Runtime1) << "Column map is consistent with the row map, good." << std::endl;

    } else {
      amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMapGO);
      GetOStream(Warnings0) << "Column map is not consistent with the row map\n"
                            << "using GO->LO conversion with performance penalty" << std::endl;
    }
#endif

    const size_t NSDim = fineNullspace->getNumVectors();
    coarseNullspace = MultiVectorFactory::Build(coarseMap, NSDim);

    // Pull out the nullspace vectors so that we can have random access.
    auto fineNS   = fineNullspace  ->template getLocalView<DeviceType>();
    auto coarseNS = coarseNullspace->template getLocalView<DeviceType>();

    size_t nnzEstimate = numRows * NSDim;

    typedef typename Matrix::local_matrix_type          local_matrix_type;
    typedef typename local_matrix_type::row_map_type    rows_type;
    typedef typename local_matrix_type::index_type      cols_type;
    typedef typename local_matrix_type::values_type     vals_type;

    typename rows_type::non_const_type rows("Ptent_aux_rows", numRows+1);
    typename cols_type::non_const_type cols("Ptent_aux_cols", nnzEstimate);
    typename vals_type::non_const_type vals("Ptent_aux_vals", nnzEstimate);

    Kokkos::parallel_for("TentativePF:BuildPuncoupled:for1", numRows+1, KOKKOS_LAMBDA(const LO row) {
      rows(row) = row*NSDim;
    });
    Kokkos::parallel_for("TentativePF:BuildUncoupled:for2", nnzEstimate, KOKKOS_LAMBDA(const LO j) {
      cols(j) = INVALID;
      vals(j) = zero;
    });

    // One thread per aggregate
    typename AppendTrait<decltype(fineNS), Kokkos::RandomAccess>::type fineNSRandom = fineNS;
    Kokkos::parallel_for("TentativePF:BuildUncoupled:main_loop", numAggs, KOKKOS_LAMBDA(const GO agg) {
      LO aggSize = aggRows(agg+1) - aggRows(agg);

      Xpetra::global_size_t offset = agg*NSDim;

      // Extract the piece of the nullspace corresponding to the aggregate, and
      // put it in the flat array, "localQR" (in column major format) for the
      // QR routine.
      // FIXME: can I create views in parallel_regions? If not, I will need to create a view with max aggregate outside?
      Kokkos::View<SC**, DeviceType> localQR("localQR", aggSize, NSDim);
      if (goodMap) {
        for (size_t j = 0; j < NSDim; j++)
          for (LO k = 0; k < aggSize; k++)
            localQR(k,j) = fineNSRandom(aggToRowMapLO(aggRows(agg)+k), j);
      } else {
        // FIXME
        throw Exceptions::RuntimeError("Not implemented");
#if 0
        for (size_t j = 0; j < NSDim; j++)
          for (LO k = 0; k < aggSize; k++)
            // FIXME
            localQR(k,j) = fineNS(rowMap->getLocalElement(aggToRowMapGO(aggStart(agg)+k)), j);
#endif
      }

      // Test for zero columns
      for (size_t j = 0; j < NSDim; j++) {
        bool bIsZeroNSColumn = true;

        for (LO k = 0; k < aggSize; k++)
          if (localQR(k,j) != zero)
            bIsZeroNSColumn = false;

        TEUCHOS_TEST_FOR_EXCEPTION(bIsZeroNSColumn == true, Exceptions::RuntimeError,
            "MueLu::TentativePFactory::MakeTentative: fine level NS part has a zero column");
      }

      // Calculate QR decomposition (standard)
      // NOTE: Q is stored in localQR and R is stored in coarseNS
      if (aggSize >= NSDim) {

        if (NSDim == 1) {
          // Only one nullspace vector, calculate Q and R by hand
          typedef Kokkos::ArithTraits<SC>  ATS;
          typedef typename ATS::magnitudeType Magnitude;

          Magnitude norm = ATS::magnitude(zero);
          for (size_t k = 0; k < Teuchos::as<size_t>(aggSize); k++)
            norm += ATS::magnitude(localQR(k,0)*localQR(k,0));
          norm = Kokkos::ArithTraits<Magnitude>::squareroot(norm);

          // R = norm
          coarseNS(offset, 0) = norm;

          // Q = localQR(:,0)/norm
          for (LO i = 0; i < aggSize; i++)
            localQR(i,0) /= norm;

        } else {
#if 1
          throw Exceptions::RuntimeError("Not implemented");
#else
          // FIXME: Need Kokkos QR solver
          Teuchos::SerialQRDenseSolver<LO,SC> qrSolver;
          qrSolver.setMatrix(Teuchos::rcp(&localQR, false));
          qrSolver.factor();

          // R = upper triangular part of localQR
          for (size_t j = 0; j < NSDim; j++)
            for (size_t k = 0; k <= j; k++)
              coarseNS(offset+k,j) = localQR(k,j); //TODO is offset+k the correct local ID?!

          // Calculate Q, the tentative prolongator.
          // The Lapack GEQRF call only works for myAggsize >= NSDim
          qrSolver.formQ();
          Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SC> > qFactor = qrSolver.getQ();
          for (size_t j = 0; j < NSDim; j++)
            for (size_t i = 0; i < Teuchos::as<size_t>(aggSize); i++)
              localQR(i,j) = (*qFactor)(i,j);
#endif
        }

      } else {
        throw; // for now
#if 0
        // Special handling for aggSize < NSDim (i.e. single node aggregates in structural mechanics)

        // The local QR decomposition is not possible in the "overconstrained"
        // case (i.e. number of columns in localQR > number of rows), which
        // corresponds to #DOFs in Aggregate < NSDim. For usual problems this
        // is only possible for single node aggregates in structural mechanics.
        // (Similar problems may arise in discontinuous Galerkin problems...)
        // We bypass the QR decomposition and use an identity block in the
        // tentative prolongator for the single node aggregate and transfer the
        // corresponding fine level null space information 1-to-1 to the coarse
        // level null space part.

        // NOTE: The resulting tentative prolongation operator has
        // (aggSize*DofsPerNode-NSDim) zero columns leading to a singular
        // coarse level operator A.  To deal with that one has the following
        // options:
        // - Use the "RepairMainDiagonal" flag in the RAPFactory (default:
        //   false) to add some identity block to the diagonal of the zero rows
        //   in the coarse level operator A, such that standard level smoothers
        //   can be used again.
        // - Use special (projection-based) level smoothers, which can deal
        //   with singular matrices (very application specific)
        // - Adapt the code below to avoid zero columns. However, we do not
        //   support a variable number of DOFs per node in MueLu/Xpetra which
        //   makes the implementation really hard.

        // R = extended (by adding identity rows) localQR
        for (size_t j = 0; j < NSDim; j++)
          for (size_t k = 0; k < NSDim; k++)
            if (k < as<size_t>(aggSize))
              coarseNS[j][offset+k] = localQR(k,j);
            else
              coarseNS[j][offset+k] = (k == j ? one : zero);

        // Q = I (rectangular)
        for (size_t i = 0; i < as<size_t>(aggSize); i++)
          for (size_t j = 0; j < NSDim; j++)
            localQR(i,j) = (j == i ? one : zero);
#endif
      }

      // Process each row in the local Q factor
      // FIXME: What happens if maps are block maps?
      for (LO j = 0; j < aggSize; j++) {
#if 1
        LO localRow = (goodMap ? aggToRowMapLO(aggRows(agg)+j) : -1);
#else
        LO localRow = (goodMap ? aggToRowMapLO[aggRows(agg)+j] : rowMap->getLocalElement(aggToRowMapGO[aggStart[agg]+j]));
#endif

        size_t rowStart = rows(localRow);
        for (size_t k = 0, lnnz = 0; k < NSDim; k++) {
          // Skip zeros (there may be plenty of them, i.e., NSDim > 1 or boundary conditions)
          if (localQR(j,k) != zero) {
            cols(rowStart+lnnz) = offset + k;
            vals(rowStart+lnnz) = localQR(j,k);
            lnnz++;
          }
        }
      }
    });

#if 0
    // TODO
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
#endif

    GetOStream(Runtime1) << "TentativePFactory : aggregates do not cross process boundaries" << std::endl;

    // Time to construct the matrix and fill in the values
    Ptentative = rcp(new CrsMatrixWrap(rowMap, coarseMap, 0, Xpetra::StaticProfile));
    RCP<CrsMatrix> PtentCrs = rcp_dynamic_cast<CrsMatrixWrap>(Ptentative)->getCrsMatrix();

    // FIXME: Here we actually need to transform a Kokkos CrsMatrix into Xpetra::Matrix
    // For now, simply do a copy-paste
    ArrayRCP<size_t>  iaPtent;
    ArrayRCP<LO>      jaPtent;
    ArrayRCP<SC>     valPtent;

    PtentCrs->allocateAllValues(nnzEstimate, iaPtent, jaPtent, valPtent);

    ArrayView<size_t> ia  = iaPtent();
    ArrayView<LO>     ja  = jaPtent();
    ArrayView<SC>     val = valPtent();

    // Copy values
    // FIXME: need host views (or DualView)
    for (LO i = 0; i < rows.size(); i++)
      ia[i] = rows(i);
    for (LO j = 0; j < cols.size(); j++) {
      ja [j] = cols(j);
      val[j] = vals(j);
    }

    PtentCrs->setAllValues(iaPtent, jaPtent, valPtent);
    PtentCrs->expertStaticFillComplete(coarseMap, A->getDomainMap());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::
  BuildPcoupled(RCP<Matrix> A, RCP<Aggregates_kokkos> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace) const {
    throw Exceptions::RuntimeError("Construction of coupled tentative P is not implemented");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::isGoodMap(const Map& rowMap, const Map& colMap) const {
    ArrayView<const GO> rowElements = rowMap.getNodeElementList();
    ArrayView<const GO> colElements = colMap.getNodeElementList();

    const size_t numElements = rowElements.size();

    bool goodMap = true;
    for (size_t i = 0; i < numElements; i++)
      if (rowElements[i] != colElements[i]) {
        goodMap = false;
        break;
      }

    return goodMap;
  }

} //namespace MueLu

#define MUELU_TENTATIVEPFACTORY_KOKKOS_SHORT
#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_TENTATIVEPFACTORY_KOKKOS_DEF_HPP
