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

  namespace { // anonymous

    template<class LocalOrdinal, class RowType>
    class ScanFunctor {
    public:
      ScanFunctor(RowType rows) : rows_(rows) { }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LocalOrdinal i, LocalOrdinal& upd, const bool& final) const {
        upd += rows_(i);
        if (final)
          rows_(i) = upd;
      }

    private:
      RowType rows_;
    };

  }

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

    auto     aggGraph = aggregates->GetGraph(); // FIXME
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

    size_t nnzEstimate = numRows * NSDim, nnz = 0;

    typedef typename Matrix::local_matrix_type          local_matrix_type;
    typedef typename local_matrix_type::row_map_type    rows_type;
    typedef typename local_matrix_type::index_type      cols_type;
    typedef typename local_matrix_type::values_type     vals_type;

    // Stage 0: initialize auxilary arrays
    // The main thing to notice is initialization of vals with INVALID. These
    // values will later be used to compress the arrays
    typename rows_type::non_const_type rowsAux("Ptent_aux_rows", numRows+1), rows("Ptent_rows");
    typename cols_type::non_const_type colsAux("Ptent_aux_cols", nnzEstimate);
    typename vals_type::non_const_type valsAux("Ptent_aux_vals", nnzEstimate);

    Kokkos::parallel_for("TentativePF:BuildPuncoupled:for1", numRows+1, KOKKOS_LAMBDA(const LO row) {
      rowsAux(row) = row*NSDim;
    });
    Kokkos::parallel_for("TentativePF:BuildUncoupled:for2", nnzEstimate, KOKKOS_LAMBDA(const LO j) {
      colsAux(j) = INVALID;
      valsAux(j) = zero;
    });

    typedef Kokkos::View<int[10], DeviceType> status_type;
    status_type status("status");

    // Stage 1: construct auxilary arrays.
    // The constructed arrays may have gaps in them (vals(j) == INAVLID)
    // Run one thread per aggregate.
    typename AppendTrait<decltype(fineNS), Kokkos::RandomAccess>::type fineNSRandom = fineNS;
    typename AppendTrait<status_type,      Kokkos::Atomic>      ::type statusAtomic = status;
    if (NSDim == 1) {
      // 1D is special, as it is the easiest. We don't even need to the QR,
      // just normalize an array. Plus, no worries abot small aggregates.
      Kokkos::parallel_reduce("TentativePF:BuildUncoupled:main_loop", numAggs, KOKKOS_LAMBDA(const GO agg, size_t& rowNnz) {
        LO aggSize = aggRows(agg+1) - aggRows(agg);

        // Extract the piece of the nullspace corresponding to the aggregate, and
        // put it in the flat array, "localQR" (in column major format) for the
        // QR routine. Trivial in 1D.
        if (goodMap) {
          // Calculate QR by hand
          typedef Kokkos::ArithTraits<SC>  ATS;
          typedef typename ATS::magnitudeType Magnitude;

          Magnitude norm = ATS::magnitude(zero);
          for (size_t k = 0; k < aggSize; k++) {
            Magnitude dnorm = ATS::magnitude(fineNSRandom(aggToRowMapLO(aggRows(agg)+k),0));
            if (dnorm == zero) {
              // This will result in a zero in a row of tentative prolongator, which is bad
              // Terminate the execution
              statusAtomic(1) = true;
              return;
            }
            norm += dnorm*dnorm;
          }

          // R = norm
          coarseNS(agg, 0) = norm;

          // Q = localQR(:,0)/norm
          for (LO k = 0; k < aggSize; k++) {
            LO localRow = aggToRowMapLO(aggRows(agg)+k);
            SC localVal = fineNSRandom(aggToRowMapLO(aggRows(agg)+k),0) / norm;

            size_t rowStart = rowsAux(localRow);
            colsAux(rowStart) = agg;
            valsAux(rowStart) = localVal;

            // Store true number of nonzeros per row
            rows(localRow+1) = 1;
          }

        } else {
          // FIXME: implement non-standard map QR
          // Look at the original TentativeP for how to do that
          statusAtomic(0) = true;
          return;
        }
      }, nnz);

      typename status_type::HostMirror statusHost = Kokkos::create_mirror_view(status);
      for (int i = 0; i < statusHost.size(); i++)
        if (statusHost(i)) {
          std::ostringstream oss;
          oss << "MueLu::TentativePFactory::MakeTentative: ";
          switch(i) {
            case 0: oss << "!goodMap is not implemented";
            case 1: oss << "fine level NS part has a zero column";
          }
          throw Exceptions::RuntimeError(oss.str());
        }

    } else {
      throw Exceptions::RuntimeError("Ignore NSDim > 1 for now");
#if 0
      Kokkos::parallel_reduce("TentativePF:BuildUncoupled:main_loop", numAggs, KOKKOS_LAMBDA(const GO agg, size_t& nnz) {
        LO aggSize = aggRows(agg+1) - aggRows(agg);

        Xpetra::global_size_t offset = agg*NSDim;

        // Extract the piece of the nullspace corresponding to the aggregate, and
        // put it in the flat array, "localQR" (in column major format) for the
        // QR routine.
        // FIXME: can I create views in parallel_regions? If not, I will need to create a view with max aggregate outside?
        // Can I create local variables? Or do I need View of Views
        Kokkos::View<SC**, DeviceType> localQR("localQR", aggSize, NSDim);
        if (goodMap) {
          for (size_t j = 0; j < NSDim; j++)
            for (LO k = 0; k < aggSize; k++)
              localQR(k,j) = fineNSRandom(aggToRowMapLO(aggRows(agg)+k), j);
        } else {
          statusAtomic(0) = true;
          return;
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

          if (bIsZeroNSColumn) {
            statusAtomic(1) = true;
            return;
          }
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
            statusAtomic(2) = true;
            return;
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
          statusAtomic(3) = true;
          return;
#if 0
          // Special handling for aggSize < NSDim (i.e. single node aggregates in structural mechanics)

          // The local QR decomposition is not possible in the "overconstrained"
          // case (i.e. number of columns in localQR > number of rowsAux), which
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
          //   false) to add some identity block to the diagonal of the zero rowsAux
          //   in the coarse level operator A, such that standard level smoothers
          //   can be used again.
          // - Use special (projection-based) level smoothers, which can deal
          //   with singular matrices (very application specific)
          // - Adapt the code below to avoid zero columns. However, we do not
          //   support a variable number of DOFs per node in MueLu/Xpetra which
          //   makes the implementation really hard.

          // R = extended (by adding identity rowsAux) localQR
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

          size_t rowStart = rowsAux(localRow), lnnz = 0;
          for (size_t k = 0; k < NSDim; k++) {
            // Skip zeros (there may be plenty of them, i.e., NSDim > 1 or boundary conditions)
            if (localQR(j,k) != zero) {
              colsAux(rowStart+lnnz) = offset + k;
              valsAux(rowStart+lnnz) = localQR(j,k);
              lnnz++;
            }
          }
          // Store true number of nonzeros per row
          rows(localRow+1) = lnnz;
        }
      }, nnz);

      typename status_type::HostMirror statusHost = Kokkos::create_mirror_view(status);
      for (int i = 0; i < statusHost.size(); i++)
        if (statusHost(i)) {
          std::ostringstream oss;
          oss << "MueLu::TentativePFactory::MakeTentative: ";
          switch(i) {
            case 0: oss << "!goodMap is not implemented";
            case 1: oss << "fine level NS part has a zero column";
            case 2: oss << "NSDim > 1 is not implemented";
            case 3: oss << "aggSize < NSDim is not imlemented";
          }
          throw Exceptions::RuntimeError(oss.str());
        }
#endif
    }

    // Stage 3: compress the arrays
    ScanFunctor<LO,decltype(rows)> scanFunctor(rows);
    Kokkos::parallel_scan("TentativePF:Build:compress_rows", numRows+1, scanFunctor);

    // The real cols and vals are constructed using calculated (not estimated) nnz
    typename cols_type::non_const_type cols("Ptent_cols", nnz);
    typename vals_type::non_const_type vals("Ptent_vals", nnz);
    Kokkos::parallel_for("TentativePF:Build:compress_cols_vals", numRows, KOKKOS_LAMBDA(const LO i) {
      LO rowStart = rows(i);

      size_t lnnz = 0;
      for (LO j = rowsAux(i); j < rowsAux(i+1); j++)
        if (valsAux(j) != INVALID) {
          cols(rowStart+lnnz) = colsAux(j);
          vals(rowStart+lnnz) = valsAux(j);
          lnnz++;
        }
    });

    GetOStream(Runtime1) << "TentativePFactory : aggregates do not cross process boundaries" << std::endl;

    // Stage 4: construct Xpetra::Matrix
    // FIXME: For now, we simply copy-paste arrays. The proper way to do that
    // would be to construct a Kokkos CrsMatrix, and then construct
    // Xpetra::Matrix out of that.
    Ptentative = rcp(new CrsMatrixWrap(rowMap, coarseMap, 0, Xpetra::StaticProfile));
    RCP<CrsMatrix> PtentCrs = rcp_dynamic_cast<CrsMatrixWrap>(Ptentative)->getCrsMatrix();

    ArrayRCP<size_t>  iaPtent;
    ArrayRCP<LO>      jaPtent;
    ArrayRCP<SC>     valPtent;

    PtentCrs->allocateAllValues(nnzEstimate, iaPtent, jaPtent, valPtent);

    ArrayView<size_t> ia  = iaPtent();
    ArrayView<LO>     ja  = jaPtent();
    ArrayView<SC>     val = valPtent();

    // Copy values
    typename rows_type::HostMirror rowsHost = Kokkos::create_mirror_view(rows);
    typename cols_type::HostMirror colsHost = Kokkos::create_mirror_view(cols);
    typename vals_type::HostMirror valsHost = Kokkos::create_mirror_view(vals);
    for (LO i = 0; i < rowsHost.size(); i++)
      ia[i] = rowsHost(i);
    for (LO j = 0; j < colsHost.size(); j++) {
      ja [j] = colsHost(j);
      val[j] = valsHost(j);
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
