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

#include "Kokkos_UnorderedMap.hpp"

#include "MueLu_TentativePFactory_kokkos_decl.hpp"

#include "MueLu_Aggregates_kokkos.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_CoarseMapFactory_kokkos.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_NullspaceFactory_kokkos.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities_kokkos.hpp"

namespace MueLu {

  namespace { // anonymous

    template<class LocalOrdinal, class View>
    class ReduceMaxFunctor{
    public:
      ReduceMaxFunctor(View view) : view_(view) { }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LocalOrdinal &i, LocalOrdinal& vmax) const {
        if (vmax < view_(i))
          vmax = view_(i);
      }

      KOKKOS_INLINE_FUNCTION
      void join (volatile LocalOrdinal& dst, const volatile LocalOrdinal& src) const {
        if (dst < src) {
          dst = src;
        }
      }

      KOKKOS_INLINE_FUNCTION
      void init (LocalOrdinal& dst) const {
        dst = 0;
      }
    private:
      View view_;
    };

    // local QR decomposition
    template<class LOType, class GOType, class SCType,class DeviceType, class NspType, class aggRowsType, class maxAggDofSizeType, class agg2RowMapLOType, class statusType, class rowsType, class rowsAuxType, class colsAuxType, class valsAuxType>
    class LocalQRDecompFunctor {
    private:
      typedef LOType LO;
      typedef GOType GO;
      typedef SCType SC;

      typedef typename DeviceType::execution_space execution_space;
      typedef Kokkos::ArithTraits<SC> ATS;
      typedef typename ATS::magnitudeType Magnitude;

      typedef Kokkos::View<SC**,typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged> shared_matrix;
      typedef Kokkos::View<SC* ,typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged> shared_vector;

    private:

      NspType fineNS;
      NspType coarseNS;
      aggRowsType aggRows;
      maxAggDofSizeType maxAggDofSize; //< maximum number of dofs in aggregate (max size of aggregate * numDofsPerNode)
      agg2RowMapLOType agg2RowMapLO;
      statusType statusAtomic;
      rowsType rows;
      rowsAuxType rowsAux;
      colsAuxType colsAux;
      valsAuxType valsAux;
      bool doQRStep;
    public:
      LocalQRDecompFunctor(NspType fineNS_, NspType coarseNS_, aggRowsType aggRows_, maxAggDofSizeType maxAggDofSize_, agg2RowMapLOType agg2RowMapLO_, statusType statusAtomic_, rowsType rows_, rowsAuxType rowsAux_, colsAuxType colsAux_, valsAuxType valsAux_, bool doQRStep_) :
        fineNS(fineNS_),
        coarseNS(coarseNS_),
        aggRows(aggRows_),
        maxAggDofSize(maxAggDofSize_),
        agg2RowMapLO(agg2RowMapLO_),
        statusAtomic(statusAtomic_),
        rows(rows_),
        rowsAux(rowsAux_),
        colsAux(colsAux_),
        valsAux(valsAux_),
        doQRStep(doQRStep_)
        { }

      KOKKOS_INLINE_FUNCTION
      void operator() ( const typename Kokkos::TeamPolicy<execution_space>::member_type & thread, size_t& nnz) const {
        auto agg = thread.league_rank();

        // size of aggregate: number of DOFs in aggregate
        auto aggSize = aggRows(agg+1) - aggRows(agg);

        const SC one     = ATS::one();
        const SC two     = one + one;
        const SC zero    = ATS::zero();
        const auto zeroM = ATS::magnitude(zero);

        int m = aggSize;
        int n = fineNS.extent(1);

        // calculate row offset for coarse nullspace
        Xpetra::global_size_t offset = agg * n;

        if (doQRStep) {

          // Extract the piece of the nullspace corresponding to the aggregate
          shared_matrix r(thread.team_shmem(), m, n);     // A (initially), R (at the end)
          for (int j = 0; j < n; j++)
            for (int k = 0; k < m; k++)
              r(k,j) = fineNS(agg2RowMapLO(aggRows(agg)+k),j);
#if 0
          printf("A\n");
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
              printf(" %5.3lf ", r(i,j));
            printf("\n");
          }
#endif

          // Calculate QR decomposition (standard)
          shared_matrix q(thread.team_shmem(), m, m);     // Q
          if (m >= n) {
            bool isSingular = false;

            // Initialize Q^T
            auto qt = q;
            for (int i = 0; i < m; i++) {
              for (int j = 0; j < m; j++)
                qt(i,j) = zero;
              qt(i,i) = one;
            }

            for (int k = 0; k < n; k++) {  // we ignore "n" instead of "n-1" to normalize
              // FIXME_KOKKOS: use team
              Magnitude s = zeroM, norm, norm_x;
              for (int i = k+1; i < m; i++)
                s += pow(ATS::magnitude(r(i,k)), 2);
              norm = sqrt(pow(ATS::magnitude(r(k,k)), 2) + s);

              if (norm == zero) {
                isSingular = true;
                break;
              }

              r(k,k) -= norm*one;

              norm_x = sqrt(pow(ATS::magnitude(r(k,k)), 2) + s);
              if (norm_x == zeroM) {
                // We have a single diagonal element in the column.
                // No reflections required. Just need to restor r(k,k).
                r(k,k) = norm*one;
                continue;
              }

              // FIXME_KOKKOS: use team
              for (int i = k; i < m; i++)
                r(i,k) /= norm_x;

              // Update R(k:m,k+1:n)
              for (int j = k+1; j < n; j++) {
                // FIXME_KOKKOS: use team in the loops
                SC si = zero;
                for (int i = k; i < m; i++)
                  si += r(i,k) * r(i,j);
                for (int i = k; i < m; i++)
                  r(i,j) -= two*si * r(i,k);
              }

              // Update Q^T (k:m,k:m)
              for (int j = k; j < m; j++) {
                // FIXME_KOKKOS: use team in the loops
                SC si = zero;
                for (int i = k; i < m; i++)
                  si += r(i,k) * qt(i,j);
                for (int i = k; i < m; i++)
                  qt(i,j) -= two*si * r(i,k);
              }

              // Fix R(k:m,k)
              r(k,k) = norm*one;
              for (int i = k+1; i < m; i++)
                r(i,k) = zero;
            }

#if 0
            // Q = (Q^T)^T
            for (int i = 0; i < m; i++)
              for (int j = 0; j < i; j++) {
                SC tmp  = qt(i,j);
                qt(i,j) = qt(j,i);
                qt(j,i) = tmp;
              }
#endif

            // Build coarse nullspace using the upper triangular part of R
            for (int j = 0; j < n; j++)
              for (int k = 0; k <= j; k++)
                coarseNS(offset+k,j) = r(k,j);

            if (isSingular) {
              statusAtomic(1) = true;
              return;
            }

          } else {
            // Special handling for m < n (i.e. single node aggregates in structural mechanics)

            // The local QR decomposition is not possible in the "overconstrained"
            // case (i.e. number of columns in qr > number of rowsAux), which
            // corresponds to #DOFs in Aggregate < n. For usual problems this
            // is only possible for single node aggregates in structural mechanics.
            // (Similar problems may arise in discontinuous Galerkin problems...)
            // We bypass the QR decomposition and use an identity block in the
            // tentative prolongator for the single node aggregate and transfer the
            // corresponding fine level null space information 1-to-1 to the coarse
            // level null space part.

            // NOTE: The resulting tentative prolongation operator has
            // (m*DofsPerNode-n) zero columns leading to a singular
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
            //
            // FIXME: do we need to check for singularity here somehow? Zero
            // columns would be easy but linear dependency would require proper QR.

            // R = extended (by adding identity rowsAux) qr
            for (int j = 0; j < n; j++)
              for (int k = 0; k < n; k++)
                if (k < m)
                  coarseNS(offset+k,j) = r(k,j);
                else
                  coarseNS(offset+k,j) = (k == j ? one : zero);

            // Q = I (rectangular)
            for (int i = 0; i < m; i++)
              for (int j = 0; j < n; j++)
                q(i,j) = (j == i ? one : zero);
          }

          // Process each row in the local Q factor and fill helper arrays to assemble P
          for (int j = 0; j < m; j++) {
            LO localRow = agg2RowMapLO(aggRows(agg)+j);
            size_t rowStart = rowsAux(localRow);
            size_t lnnz = 0;
            for (int k = 0; k < n; k++) {
              // skip zeros
              if (q(j,k) != zero) {
                colsAux(rowStart+lnnz) = offset + k;
                valsAux(rowStart+lnnz) = q(j,k);
                lnnz++;
              }
            }
            rows(localRow+1) = lnnz;
            nnz += lnnz;
          }

#if 0
          printf("R\n");
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
              printf(" %5.3lf ", coarseNS(i,j));
            printf("\n");
          }

          printf("Q\n");
          for (int i = 0; i < aggSize; i++) {
            for (int j = 0; j < aggSize; j++)
              printf(" %5.3lf ", q(i,j));
            printf("\n");
          }
#endif
        } else {
          /////////////////////////////
          //      "no-QR" option     //
          /////////////////////////////
          // Local Q factor is just the fine nullspace support over the current aggregate.
          // Local R factor is the identity.
          // TODO I have not implemented any special handling for aggregates that are too
          // TODO small to locally support the nullspace, as is done in the standard QR
          // TODO case above.

          for (int j = 0; j < m; j++) {
            LO localRow = agg2RowMapLO(aggRows(agg)+j);
            size_t rowStart = rowsAux(localRow);
            size_t lnnz = 0;
            for (int k = 0; k < n; k++) {
              const SC qr_jk = fineNS(localRow,k);
              // skip zeros
              if (qr_jk != zero) {
                colsAux(rowStart+lnnz) = offset + k;
                valsAux(rowStart+lnnz) = qr_jk;
                lnnz++;
              }
            }
            rows(localRow+1) = lnnz;
            nnz += lnnz;
          }

          for (int j = 0; j < n; j++)
            coarseNS(offset+j,j) = one;

        }

      }

      // amount of shared memory
      size_t team_shmem_size( int team_size ) const {
        if (doQRStep) {
          int m = maxAggDofSize;
          int n = fineNS.extent(1);
          return shared_matrix::shmem_size(m, n) +    // r
            shared_matrix::shmem_size(m, m);     // q
        } else
          return 0;
      }
    };

  } // namespace anonymous

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const ParameterList> TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("tentative: calculate qr");
  SET_VALID_ENTRY("tentative: build coarse coordinates");
#undef SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Aggregates",         Teuchos::null, "Generating factory of the aggregates");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",          Teuchos::null, "Generating factory of the nullspace");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
    validParamList->set< RCP<const FactoryBase> >("CoarseMap",          Teuchos::null, "Generating factory of the coarse map");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",        Teuchos::null, "Generating factory of the coordinates");

    // Make sure we don't recursively validate options for the matrixmatrix kernels
    ParameterList norecurse;
    norecurse.disableRecursiveValidation();
    validParamList->set<ParameterList> ("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {

    const ParameterList& pL = GetParameterList();

    Input(fineLevel, "A");
    Input(fineLevel, "Aggregates");
    Input(fineLevel, "Nullspace");
    Input(fineLevel, "UnAmalgamationInfo");
    Input(fineLevel, "CoarseMap");
    if( fineLevel.GetLevelID() == 0 &&
        fineLevel.IsAvailable("Coordinates", NoFactory::get()) &&     // we have coordinates (provided by user app)
        pL.get<bool>("tentative: build coarse coordinates") ) {       // and we want coordinates on other levels
      bTransferCoordinates_ = true;                                   // then set the transfer coordinates flag to true
      Input(fineLevel, "Coordinates");
    } else if (bTransferCoordinates_) {
      Input(fineLevel, "Coordinates");
    }
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
    RCP<RealValuedMultiVector> fineCoords;
    if(bTransferCoordinates_) {
      fineCoords = Get< RCP<RealValuedMultiVector> >(fineLevel, "Coordinates");
    }

    RCP<Matrix>      Ptentative;
    RCP<MultiVector> coarseNullspace;
    RCP<RealValuedMultiVector> coarseCoords;

    if(bTransferCoordinates_) {
      ArrayView<const GO> elementAList = coarseMap->getNodeElementList();
      GO                  indexBase    = coarseMap->getIndexBase();

      LO blkSize = 1;
      if (rcp_dynamic_cast<const StridedMap>(coarseMap) != Teuchos::null)
        blkSize = rcp_dynamic_cast<const StridedMap>(coarseMap)->getFixedBlockSize();

      Array<GO>           elementList;
      ArrayView<const GO> elementListView;
      if (blkSize == 1) {
        // Scalar system
        // No amalgamation required
        elementListView = elementAList;

      } else {
        auto numElements = elementAList.size() / blkSize;

        elementList.resize(numElements);

        // Amalgamate the map
        for (LO i = 0; i < Teuchos::as<LO>(numElements); i++)
          elementList[i] = (elementAList[i*blkSize]-indexBase)/blkSize + indexBase;

        elementListView = elementList;
      }

      auto uniqueMap      = fineCoords->getMap();
      auto coarseCoordMap = MapFactory::Build(coarseMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                              elementListView, indexBase, coarseMap->getComm());
      coarseCoords = RealValuedMultiVectorFactory::Build(coarseCoordMap, fineCoords->getNumVectors());

      // Create overlapped fine coordinates to reduce global communication
      RCP<RealValuedMultiVector> ghostedCoords = fineCoords;
      if (aggregates->AggregatesCrossProcessors()) {
        auto nonUniqueMap = aggregates->GetMap();
        auto importer     = ImportFactory::Build(uniqueMap, nonUniqueMap);

        ghostedCoords = RealValuedMultiVectorFactory::Build(nonUniqueMap, fineCoords->getNumVectors());
        ghostedCoords->doImport(*fineCoords, *importer, Xpetra::INSERT);
      }

      // The good new is that his graph has already been constructed for the
      // TentativePFactory and was cached in Aggregates. So this is a no-op.
      auto aggGraph = aggregates->GetGraph();
      auto numAggs  = aggGraph.numRows();

      auto fineCoordsView   = fineCoords  ->template getLocalView<DeviceType>();
      auto coarseCoordsView = coarseCoords->template getLocalView<DeviceType>();

      // Fill in coarse coordinates
      {
        SubFactoryMonitor m2(*this, "AverageCoords", coarseLevel);

        const auto dim = fineCoords->getNumVectors();

        typename AppendTrait<decltype(fineCoordsView), Kokkos::RandomAccess>::type fineCoordsRandomView = fineCoordsView;
        for (size_t j = 0; j < dim; j++) {
          Kokkos::parallel_for("MueLu::TentativeP::BuildCoords", Kokkos::RangePolicy<local_ordinal_type, execution_space>(0, numAggs),
                               KOKKOS_LAMBDA(const LO i) {
                                 // A row in this graph represents all node ids in the aggregate
                                 // Therefore, averaging is very easy

                                 auto aggregate = aggGraph.rowConst(i);

                                 double sum = 0.0; // do not use Scalar here (Stokhos)
                                 for (size_t colID = 0; colID < static_cast<size_t>(aggregate.length); colID++)
                                   sum += fineCoordsRandomView(aggregate(colID),j);

                                 coarseCoordsView(i,j) = sum / aggregate.length;
                               });
        }
      }
    }

    if (!aggregates->AggregatesCrossProcessors())
      BuildPuncoupled(coarseLevel, A, aggregates, amalgInfo, fineNullspace, coarseMap, Ptentative, coarseNullspace, coarseLevel.GetLevelID());
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

    if(bTransferCoordinates_) {
      Set(coarseLevel, "Coordinates", coarseCoords);
    }
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
  BuildPuncoupled(Level& coarseLevel, RCP<Matrix> A, RCP<Aggregates_kokkos> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                  RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace, const int levelID) const {
    auto rowMap = A->getRowMap();
    auto colMap = A->getColMap();

    const size_t numRows  = rowMap->getNodeNumElements();
    const size_t NSDim    = fineNullspace->getNumVectors();

    typedef Kokkos::ArithTraits<SC>     ATS;
    typedef typename ATS::magnitudeType Magnitude;
    const SC zero = ATS::zero(), one = ATS::one();

    const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();

    typename Aggregates_kokkos::local_graph_type aggGraph;
    {
      SubFactoryMonitor m2(*this, "Get Aggregates graph", coarseLevel);
      aggGraph = aggregates->GetGraph();
    }
    auto aggRows  = aggGraph.row_map;
    auto aggCols  = aggGraph.entries;

    // Aggregates map is based on the amalgamated column map
    // We can skip global-to-local conversion if LIDs in row map are
    // same as LIDs in column map
    bool goodMap;
    {
      SubFactoryMonitor m2(*this, "Check good map", coarseLevel);
      goodMap = isGoodMap(*rowMap, *colMap);
    }
    // FIXME_KOKKOS: need to proofread later code for bad maps
    TEUCHOS_TEST_FOR_EXCEPTION(!goodMap, Exceptions::RuntimeError,
        "MueLu: TentativePFactory_kokkos: for now works only with good maps "
        "(i.e. \"matching\" row and column maps)");

    // STEP 1: do unamalgamation
    // The non-kokkos version uses member functions from the AmalgamationInfo
    // container class to unamalgamate the data. In contrast, the kokkos
    // version of TentativePFactory does the unamalgamation here and only uses
    // the data of the AmalgamationInfo container class

    // Extract information for unamalgamation
    LO fullBlockSize, blockID, stridingOffset, stridedBlockSize;
    GO indexBase;
    amalgInfo->GetStridingInformation(fullBlockSize, blockID, stridingOffset, stridedBlockSize, indexBase);
    GO globalOffset = amalgInfo->GlobalOffset();

    // Extract aggregation info (already in Kokkos host views)
    auto         procWinner    = aggregates->GetProcWinner()  ->template getLocalView<DeviceType>();
    auto         vertex2AggId  = aggregates->GetVertex2AggId()->template getLocalView<DeviceType>();
    const size_t numAggregates = aggregates->GetNumAggregates();

    int myPID = aggregates->GetMap()->getComm()->getRank();

    // Create Kokkos::View (on the device) to store the aggreate dof sizes
    // Later used to get aggregate dof offsets
    // NOTE: This zeros itself on construction
    typedef typename Aggregates_kokkos::aggregates_sizes_type::non_const_type AggSizeType;
    AggSizeType aggDofSizes;

    if (stridedBlockSize == 1) {
      SubFactoryMonitor m2(*this, "Calc AggSizes", coarseLevel);

      // FIXME_KOKKOS: use ViewAllocateWithoutInitializing + set a single value
      aggDofSizes = AggSizeType("agg_dof_sizes", numAggregates+1);

      auto sizesConst = aggregates->ComputeAggregateSizes();
      Kokkos::deep_copy(Kokkos::subview(aggDofSizes, Kokkos::make_pair(static_cast<size_t>(1), numAggregates+1)), sizesConst);

    } else {
      SubFactoryMonitor m2(*this, "Calc AggSizes", coarseLevel);

      // FIXME_KOKKOS: use ViewAllocateWithoutInitializing + set a single value
      aggDofSizes = AggSizeType("agg_dof_sizes", numAggregates + 1);

      auto nodeMap = aggregates->GetMap()->getLocalMap();
      auto dofMap  = colMap->getLocalMap();

      Kokkos::parallel_for("MueLu:TentativePF:Build:compute_agg_sizes", range_type(0,numAggregates),
        KOKKOS_LAMBDA(const LO agg) {
          auto aggRowView = aggGraph.rowConst(agg);

          size_t size = 0;
          for (LO colID = 0; colID < aggRowView.length; colID++) {
            GO nodeGID = nodeMap.getGlobalElement(aggRowView(colID));

            for (LO k = 0; k < stridedBlockSize; k++) {
              GO dofGID = (nodeGID - indexBase) * fullBlockSize + k + indexBase + globalOffset + stridingOffset;

              if (dofMap.getLocalElement(dofGID) != INVALID)
                size++;
            }
          }
          aggDofSizes(agg+1) = size;
        });
    }

    // Find maximum dof size for aggregates
    // Later used to reserve enough scratch space for local QR decompositions
    LO maxAggSize = 0;
    ReduceMaxFunctor<LO,decltype(aggDofSizes)> reduceMax(aggDofSizes);
    Kokkos::parallel_reduce("MueLu:TentativePF:Build:max_agg_size", range_type(0, aggDofSizes.extent(0)), reduceMax, maxAggSize);

    // parallel_scan (exclusive)
    // The aggDofSizes View then contains the aggregate dof offsets
    Kokkos::parallel_scan("MueLu:TentativePF:Build:aggregate_sizes:stage1_scan", range_type(0,numAggregates+1),
      KOKKOS_LAMBDA(const LO i, LO& update, const bool& final_pass) {
        update += aggDofSizes(i);
        if (final_pass)
          aggDofSizes(i) = update;
      });

    // Create Kokkos::View on the device to store mapping
    // between (local) aggregate id and row map ids (LIDs)
    Kokkos::View<LO*, DeviceType> agg2RowMapLO(Kokkos::ViewAllocateWithoutInitializing("agg2row_map_LO"), numRows);
    {
      SubFactoryMonitor m2(*this, "Create Agg2RowMap", coarseLevel);

      AggSizeType aggOffsets(Kokkos::ViewAllocateWithoutInitializing("aggOffsets"), numAggregates);
      Kokkos::deep_copy(aggOffsets, Kokkos::subview(aggDofSizes, Kokkos::make_pair(static_cast<size_t>(0), numAggregates)));

      Kokkos::parallel_for("MueLu:TentativePF:Build:createAgg2RowMap", range_type(0, vertex2AggId.extent(0)),
        KOKKOS_LAMBDA(const LO lnode) {
          if (procWinner(lnode, 0) == myPID) {
            // No need for atomics, it's one-to-one
            auto aggID = vertex2AggId(lnode,0);

            auto offset = Kokkos::atomic_fetch_add( &aggOffsets(aggID), stridedBlockSize );
            // FIXME: I think this may be wrong
            // We unconditionally add the whole block here. When we calculated
            // aggDofSizes, we did the isLocalElement check. Something's fishy.
            for (LO k = 0; k < stridedBlockSize; k++)
              agg2RowMapLO(offset + k) = lnode*stridedBlockSize + k;
          }
        });
    }

    // STEP 2: prepare local QR decomposition
    // Reserve memory for tentative prolongation operator
    coarseNullspace = MultiVectorFactory::Build(coarseMap, NSDim);

    // Pull out the nullspace vectors so that we can have random access (on the device)
    auto fineNS   = fineNullspace  ->template getLocalView<DeviceType>();
    auto coarseNS = coarseNullspace->template getLocalView<DeviceType>();

    size_t nnz = 0;                       // actual number of nnz

    typedef typename Xpetra::Matrix<SC,LO,GO,NO>::local_matrix_type    local_matrix_type;
    typedef typename local_matrix_type::row_map_type::non_const_type   rows_type;
    typedef typename local_matrix_type::index_type::non_const_type     cols_type;
    typedef typename local_matrix_type::values_type::non_const_type    vals_type;


    // Device View for status (error messages...)
    typedef Kokkos::View<int[10], DeviceType> status_type;
    status_type status("status");

    typename AppendTrait<decltype(fineNS), Kokkos::RandomAccess>::type fineNSRandom = fineNS;
    typename AppendTrait<status_type,      Kokkos::Atomic>      ::type statusAtomic = status;

    rows_type rows;
    cols_type cols;
    vals_type vals;

    const ParameterList& pL = GetParameterList();
    const bool& doQRStep = pL.get<bool>("tentative: calculate qr");
    if (!doQRStep) {
      GetOStream(Runtime1) << "TentativePFactory : bypassing local QR phase" << std::endl;
      if (NSDim>1)
        GetOStream(Warnings0) << "TentativePFactor : for nontrivial nullspace, this may degrade performance" << std::endl;
    }

    if (NSDim == 1) {
      // 1D is special, as it is the easiest. We don't even need to the QR,
      // just normalize an array. Plus, no worries abot small aggregates.  In
      // addition, we do not worry about compression. It is unlikely that
      // nullspace will have zeros. If it does, a prolongator row would be
      // zero and we'll get singularity anyway.
      SubFactoryMonitor m2(*this, "Stage 1 (LocalQR)", coarseLevel);

      nnz = numRows;

      // FIXME_KOKKOS: use ViewAllocateWithoutInitializing + set a single value
      rows = rows_type("Ptent_rows", numRows+1);
      cols = cols_type(Kokkos::ViewAllocateWithoutInitializing("Ptent_cols"), numRows);
      vals = vals_type(Kokkos::ViewAllocateWithoutInitializing("Ptent_vals"), numRows);

      // Set up team policy with numAggregates teams and one thread per team.
      // Each team handles a slice of the data associated with one aggregate
      // and performs a local QR decomposition (in this case real QR is
      // unnecessary).
      const Kokkos::TeamPolicy<execution_space> policy(numAggregates, 1);

      if (doQRStep) {
        Kokkos::parallel_for("MueLu:TentativePF:BuildUncoupled:main_loop", policy,
          KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<execution_space>::member_type &thread) {
            auto agg = thread.league_rank();

            // size of the aggregate (number of DOFs in aggregate)
            LO aggSize = aggRows(agg+1) - aggRows(agg);

            // Extract the piece of the nullspace corresponding to the aggregate, and
            // put it in the flat array, "localQR" (in column major format) for the
            // QR routine. Trivial in 1D.
            auto norm = ATS::magnitude(zero);

            // Calculate QR by hand
            // FIXME: shouldn't there be stridedblock here?
            // FIXME_KOKKOS: shouldn't there be stridedblock here?
            for (decltype(aggSize) k = 0; k < aggSize; k++) {
              auto dnorm = ATS::magnitude(fineNSRandom(agg2RowMapLO(aggRows(agg)+k),0));
              norm += dnorm*dnorm;
            }
            norm = sqrt(norm);

            if (norm == zero) {
              // zero column; terminate the execution
              statusAtomic(1) = true;
              return;
            }

            // R = norm
            coarseNS(agg, 0) = norm;

            // Q = localQR(:,0)/norm
            for (decltype(aggSize) k = 0; k < aggSize; k++) {
              LO localRow = agg2RowMapLO(aggRows(agg)+k);
              SC localVal = fineNSRandom(agg2RowMapLO(aggRows(agg)+k),0) / norm;

              rows(localRow+1) = localRow+1;
              cols(localRow) = agg;
              vals(localRow) = localVal;

            }
          });

        typename status_type::HostMirror statusHost = Kokkos::create_mirror_view(status);
        Kokkos::deep_copy(statusHost, status);
        for (decltype(statusHost.size()) i = 0; i < statusHost.size(); i++)
          if (statusHost(i)) {
            std::ostringstream oss;
            oss << "MueLu::TentativePFactory::MakeTentative: ";
            switch (i) {
              case 0: oss << "!goodMap is not implemented";               break;
              case 1: oss << "fine level NS part has a zero column";      break;
            }
            throw Exceptions::RuntimeError(oss.str());
          }

      } else {
        Kokkos::parallel_for("MueLu:TentativePF:BuildUncoupled:main_loop_noqr", policy,
          KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<execution_space>::member_type &thread) {
            auto agg = thread.league_rank();

            // size of the aggregate (number of DOFs in aggregate)
            LO aggSize = aggRows(agg+1) - aggRows(agg);

            // R = norm
            coarseNS(agg, 0) = one;

            // Q = localQR(:,0)/norm
            for (decltype(aggSize) k = 0; k < aggSize; k++) {
              LO localRow = agg2RowMapLO(aggRows(agg)+k);
              SC localVal = fineNSRandom(agg2RowMapLO(aggRows(agg)+k),0);

              rows(localRow+1) = localRow+1;
              cols(localRow) = agg;
              vals(localRow) = localVal;

            }
          });
      }

    } else { // NSdim > 1
      // FIXME_KOKKOS: This code branch is completely unoptimized.
      // Work to do:
      //   - Optimize QR decomposition
      //   - Remove INVALID usage similarly to CoalesceDropFactory_kokkos by
      //     packing new values in the beginning of each row
      // We do use auxilary view in this case, so keep a second rows view for
      // counting nonzeros in rows

      // NOTE: the allocation (initialization) of these view takes noticeable time
      size_t nnzEstimate = numRows * NSDim;
      rows_type rowsAux("Ptent_aux_rows", numRows+1);
      cols_type colsAux("Ptent_aux_cols", nnzEstimate);
      vals_type valsAux("Ptent_aux_vals", nnzEstimate);
      rows = rows_type("Ptent_rows", numRows+1);
      {
        // Stage 0: fill in views.
        SubFactoryMonitor m2(*this, "Stage 0 (InitViews)", coarseLevel);

        // The main thing to notice is initialization of vals with INVALID. These
        // values will later be used to compress the arrays
        Kokkos::parallel_for("MueLu:TentativePF:BuildPuncoupled:for1", range_type(0, numRows+1),
          KOKKOS_LAMBDA(const LO row) {
            rowsAux(row) = row*NSDim;
          });
        Kokkos::parallel_for("MueLu:TentativePF:BuildUncoupled:for2", range_type(0, nnzEstimate),
          KOKKOS_LAMBDA(const LO j) {
            colsAux(j) = INVALID;
            valsAux(j) = zero;
          });
      }

      {
        SubFactoryMonitor m2 = SubFactoryMonitor(*this, doQRStep ? "Stage 1 (LocalQR)" : "Stage 1 (Fill coarse nullspace and tentative P)", coarseLevel);
        // Set up team policy with numAggregates teams and one thread per team.
        // Each team handles a slice of the data associated with one aggregate
        // and performs a local QR decomposition
        const Kokkos::TeamPolicy<execution_space> policy(numAggregates,1); // numAggregates teams a 1 thread
        LocalQRDecompFunctor<LocalOrdinal, GlobalOrdinal, Scalar, DeviceType, decltype(fineNSRandom),
            decltype(aggDofSizes /*aggregate sizes in dofs*/), decltype(maxAggSize), decltype(agg2RowMapLO),
            decltype(statusAtomic), decltype(rows), decltype(rowsAux), decltype(colsAux),
            decltype(valsAux)>
                localQRFunctor(fineNSRandom, coarseNS, aggDofSizes, maxAggSize, agg2RowMapLO, statusAtomic,
                               rows, rowsAux, colsAux, valsAux, doQRStep);
        Kokkos::parallel_reduce("MueLu:TentativePF:BuildUncoupled:main_qr_loop", policy, localQRFunctor, nnz);
      }

      typename status_type::HostMirror statusHost = Kokkos::create_mirror_view(status);
      Kokkos::deep_copy(statusHost, status);
      for (decltype(statusHost.size()) i = 0; i < statusHost.size(); i++)
        if (statusHost(i)) {
          std::ostringstream oss;
          oss << "MueLu::TentativePFactory::MakeTentative: ";
          switch(i) {
            case 0: oss << "!goodMap is not implemented";               break;
            case 1: oss << "fine level NS part has a zero column";      break;
          }
          throw Exceptions::RuntimeError(oss.str());
        }

      // Compress the cols and vals by ignoring INVALID column entries that correspond
      // to 0 in QR.

      // The real cols and vals are constructed using calculated (not estimated) nnz
      cols = decltype(cols)("Ptent_cols", nnz);
      vals = decltype(vals)("Ptent_vals", nnz);
      {
        // Stage 2: compress the arrays
        SubFactoryMonitor m2(*this, "Stage 2 (CompressRows)", coarseLevel);

        Kokkos::parallel_scan("MueLu:TentativePF:Build:compress_rows", range_type(0,numRows+1),
          KOKKOS_LAMBDA(const LO i, LO& upd, const bool& final) {
            upd += rows(i);
            if (final)
              rows(i) = upd;
          });
      }

      {
        SubFactoryMonitor m2(*this, "Stage 2 (CompressCols)", coarseLevel);

        // FIXME_KOKKOS: this can be spedup by moving correct cols and vals values
        // to the beginning of rows. See CoalesceDropFactory_kokkos for
        // example.
        Kokkos::parallel_for("MueLu:TentativePF:Build:compress_cols_vals", range_type(0,numRows),
          KOKKOS_LAMBDA(const LO i) {
            LO rowStart = rows(i);

            size_t lnnz = 0;
            for (auto j = rowsAux(i); j < rowsAux(i+1); j++)
              if (colsAux(j) != INVALID) {
                cols(rowStart+lnnz) = colsAux(j);
                vals(rowStart+lnnz) = valsAux(j);
                lnnz++;
              }
          });
      }
    }

    GetOStream(Runtime1) << "TentativePFactory : aggregates do not cross process boundaries" << std::endl;

    {
      // Stage 3: construct Xpetra::Matrix
      SubFactoryMonitor m2(*this, "Stage 3 (LocalMatrix+FillComplete)", coarseLevel);

      local_matrix_type lclMatrix = local_matrix_type("A", numRows, coarseMap->getNodeNumElements(), nnz, vals, rows, cols);

      // Managing labels & constants for ESFC
      RCP<ParameterList> FCparams;
      if (pL.isSublist("matrixmatrix: kernel params"))
        FCparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
      else
        FCparams = rcp(new ParameterList);

      // By default, we don't need global constants for TentativeP
      FCparams->set("compute global constants", FCparams->get("compute global constants", false));
      FCparams->set("Timer Label",              std::string("MueLu::TentativeP-") + toString(levelID));

      auto PtentCrs = CrsMatrixFactory::Build(lclMatrix, rowMap, coarseMap, coarseMap, A->getDomainMap());
      Ptentative = rcp(new CrsMatrixWrap(PtentCrs));
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::
  BuildPcoupled(RCP<Matrix> A, RCP<Aggregates_kokkos> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace) const {
    throw Exceptions::RuntimeError("MueLu: Construction of coupled tentative P is not implemented");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::
  isGoodMap(const Map& rowMap, const Map& colMap) const {
    auto rowLocalMap = rowMap.getLocalMap();
    auto colLocalMap = colMap.getLocalMap();

    const size_t numRows = rowLocalMap.getNodeNumElements();
    const size_t numCols = colLocalMap.getNodeNumElements();

    if (numCols < numRows)
      return false;

    size_t numDiff = 0;
    Kokkos::parallel_reduce("MueLu:TentativePF:isGoodMap", range_type(0, numRows),
      KOKKOS_LAMBDA(const LO i, size_t &diff) {
        diff += (rowLocalMap.getGlobalElement(i) != colLocalMap.getGlobalElement(i));
      }, numDiff);

    return (numDiff == 0);
  }

} //namespace MueLu

#define MUELU_TENTATIVEPFACTORY_KOKKOS_SHORT
#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_TENTATIVEPFACTORY_KOKKOS_DEF_HPP
