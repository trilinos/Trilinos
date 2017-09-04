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

      typedef Kokkos::View<SC**,execution_space,Kokkos::MemoryUnmanaged> shared_matrix;
      typedef Kokkos::View<SC* ,execution_space,Kokkos::MemoryUnmanaged> shared_vector;

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
    public:
      LocalQRDecompFunctor(NspType fineNS_, NspType coarseNS_, aggRowsType aggRows_, maxAggDofSizeType maxAggDofSize_, agg2RowMapLOType agg2RowMapLO_, statusType statusAtomic_, rowsType rows_, rowsAuxType rowsAux_, colsAuxType colsAux_, valsAuxType valsAux_) :
        fineNS(fineNS_),
        coarseNS(coarseNS_),
        aggRows(aggRows_),
        maxAggDofSize(maxAggDofSize_),
        agg2RowMapLO(agg2RowMapLO_),
        statusAtomic(statusAtomic_),
        rows(rows_),
        rowsAux(rowsAux_),
        colsAux(colsAux_),
        valsAux(valsAux_)
        { }

      KOKKOS_INLINE_FUNCTION
      void operator() ( const typename Kokkos::TeamPolicy<execution_space>::member_type & thread, size_t& nnz) const {
        auto agg = thread.league_rank();

        // size of aggregate: number of DOFs in aggregate
        auto aggSize = aggRows(agg+1) - aggRows(agg);

        typedef Kokkos::ArithTraits<SC>     ATS;
        SC one  = ATS::one();
        SC zero = ATS::zero();

        // Extract the piece of the nullspace corresponding to the aggregate, and
        // put it in the flat array, "localQR" (in column major format) for the
        // QR routine.
        shared_matrix localQR(thread.team_shmem(),aggSize,fineNS.dimension_1());
        for (size_t j = 0; j < fineNS.dimension_1(); j++)
          for (LO k = 0; k < aggSize; k++)
            localQR(k,j) = fineNS(agg2RowMapLO(aggRows(agg)+k),j);

        // Test for zero columns
        for (size_t j = 0; j < fineNS.dimension_1(); j++) {
          bool bIsZeroNSColumn = true;

          for (LO k = 0; k < aggSize; k++)
            if (localQR(k,j) != zero)
              bIsZeroNSColumn = false;

          if (bIsZeroNSColumn) {
            statusAtomic(1) = true;
            return;
          }
        }
        // calculate row offset for coarse nullspace
        Xpetra::global_size_t offset = agg * fineNS.dimension_1();

        // Reserve shared memory for local QR decomposition and results (r and q)
        shared_matrix q (thread.team_shmem(),aggSize,aggSize);  // memory containing q part in the end

        // Calculate QR decomposition (standard)
        if (aggSize >= decltype(aggSize)( fineNS.dimension_1() )) {

          // Reserve temporary shared memory for local QR decomposition
          shared_matrix r(thread.team_shmem(),aggSize,fineNS.dimension_1()); // memory containing the r part in the end
          shared_matrix z(thread.team_shmem(),aggSize,fineNS.dimension_1()); // helper matrix (containing parts of localQR)
          shared_vector e(thread.team_shmem(),aggSize); //unit vector
          shared_matrix qk(thread.team_shmem(),aggSize,aggSize);  // memory cotaining one householder reflection part
          shared_matrix qt (thread.team_shmem(),aggSize,aggSize); // temporary
          shared_vector x(thread.team_shmem(),aggSize);

          // standard case
          // do local QR decomposition
          matrix_copy(localQR,z);

          typedef typename ATS::magnitudeType Magnitude;
          Magnitude zeroMagnitude = Kokkos::ArithTraits<Magnitude>::zero();

          for(decltype(localQR.dimension_0()) k = 0; k < localQR.dimension_0() && k < localQR.dimension_1()/*-1*/; k++) {
            // extract minor parts from mat
            matrix_clear(r);  // zero out temporary matrix (there is some potential to speed this up by avoiding this)
            matrix_minor(z,r,k);

            // extract k-th column from current minor
            //auto x  = subview(r, Kokkos::ALL (), k);
            for(decltype(x.dimension_0()) i = 0; i < x.dimension_0(); i++)
              x(i) = r(i,k);
            SC   xn = vnorm(x); // calculate 2-norm of current column vector
            if(ATS::magnitude(localQR(k,k)) > zeroMagnitude) xn = -xn;

            // build k-th unit vector
            for(decltype(e.dimension_0()) i = 0; i < e.dimension_0(); i++)
              e(i) = (i==k) ?  one : zero;
            vmadd(e, x, xn); // e = x + xn * e;
            SC en = vnorm(e);
            vdiv(e,en); // scale vector e

            // build Q(k) and Q matrix
            if (k == 0) {
              vmul ( e, q);
              matrix_mul(q, r, z);
            }
            else {
              matrix_clear(qk); // zero out old qk
              vmul ( e, qk);

              matrix_mul ( qk, q, qt);  // TODO can we avoid qt?
              matrix_copy(qt,q);
              matrix_mul(qk, r, z);
            }
          }

          // build R part
          matrix_mul ( q, localQR, r);

          // upper triangular part of R build coarse NS
          for(decltype(fineNS.dimension_1()) j = 0; j < fineNS.dimension_1(); j++)
            for(decltype(j) k = 0; k <= j; k++)
              coarseNS(offset+k,j) = r(k,j);
          // Don't forget to transpose q
          matrix_transpose(q);
        } else {
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
          for (decltype(fineNS.dimension_1()) j = 0; j < fineNS.dimension_1(); j++)
            for (decltype(fineNS.dimension_1()) k = 0; k < fineNS.dimension_1(); k++)
              if (k < decltype(k)(aggSize))
                coarseNS(offset+k,j) = localQR(k,j);
              else
                coarseNS(offset+k,j) = (k == j ? one : zero);

          // Q = I (rectangular)
          for (decltype(fineNS.dimension_1()) i = 0; i < decltype(fineNS.dimension_1()) (aggSize); i++)
            for (decltype(fineNS.dimension_1()) j = 0; j < fineNS.dimension_1(); j++)
              q(i,j) = (j == i ? one : zero);
        }

        // Process each row in the local Q factor and fill helper arrays to assemble P
        for (LO j = 0; j < aggSize; j++) {
          LO localRow = agg2RowMapLO(aggRows(agg)+j);
          size_t rowStart = rowsAux(localRow);
          size_t lnnz = 0;
          for (decltype(fineNS.dimension_1()) k = 0; k < fineNS.dimension_1(); k++) {
            // skip zeros
            if(q(j,k) != zero) {
              colsAux(rowStart+lnnz) = offset + k;
              valsAux(rowStart+lnnz) = q(j,k);
              lnnz++;
            }
          }
          rows(localRow+1) = lnnz;
          nnz += lnnz;
        }

        // debug
        /*printf("R\n");
        for(int i=0; i<aggSize; i++) {
          for(int j=0; j<fineNS.dimension_1(); j++) {
            printf(" %.3g ",coarseNS(i,j));
          }
          printf("\n");
        }
        printf("Q\n");

        for(int i=0; i<aggSize; i++) {
          for(int j=0; j<aggSize; j++) {
            printf(" %.3g ",q(i,j));
          }
          printf("\n");
        }*/
        // end debug

      }

      KOKKOS_INLINE_FUNCTION
      void matrix_clear ( shared_matrix & m) const {
        typedef Kokkos::ArithTraits<SC>     ATS;
        SC zero = ATS::zero();
        for(decltype(m.dimension_0()) i = 0; i < m.dimension_0(); i++) {
          for(decltype(m.dimension_1()) j = 0; j < m.dimension_1(); j++) {
            m(i,j) = zero;
          }
        }
      }

      KOKKOS_INLINE_FUNCTION
      void matrix_copy (const shared_matrix & mi, shared_matrix & mo) const {
        for(decltype(mi.dimension_0()) i = 0; i < mi.dimension_0(); i++) {
          for(decltype(mi.dimension_1()) j = 0; j < mi.dimension_1(); j++) {
            mo(i,j) = mi(i,j);
          }
        }
      }

      KOKKOS_FUNCTION
      void matrix_transpose ( shared_matrix & m) const {
        for(decltype(m.dimension_0()) i = 0; i < m.dimension_0(); i++) {
          for(decltype(m.dimension_1()) j = 0; j < i; j++) {
            SC t = m(i,j);
            m(i,j) = m(j,i);
            m(j,i) = t;
          }
        }
      }

      KOKKOS_FUNCTION
      void matrix_mul ( const shared_matrix & m1, const shared_matrix & m2, shared_matrix & m1m2) const {
        typedef Kokkos::ArithTraits<SC>     ATS;
        SC zero = ATS::zero();
        for(decltype(m1.dimension_0()) i = 0; i < m1.dimension_0(); i++)
          for(decltype(m2.dimension_1()) j = 0; j < m2.dimension_1(); j++) {
            m1m2(i,j) = zero;
            for(decltype(m1.dimension_1()) k = 0; k < m1.dimension_1(); k++)
              m1m2(i,j) += m1(i,k) * m2(k,j);
          }
      }

      KOKKOS_FUNCTION
      void matrix_minor ( const shared_matrix & mat, shared_matrix & matminor, LO d) const {
        typedef Kokkos::ArithTraits<SC>     ATS;
        SC one = ATS::one();

        for (decltype(d) i = 0; i < d; i++) {
          matminor(i,i) = one;
        }
        for (decltype(mat.dimension_0()) i = d; i < mat.dimension_0(); i++) {
          for (decltype(mat.dimension_1()) j=d; j < mat.dimension_1(); j++) {
            matminor(i,j) = mat(i,j);
          }
        }
      }

      /// \brief Build vmul = I - v*v^T
      /// \param v[in] input vector
      ///
      KOKKOS_FUNCTION
      void vmul ( const shared_vector & v, shared_matrix & vmuldata) const {
        typedef Kokkos::ArithTraits<SC>     ATS;
        SC one = ATS::one();
        for(decltype(v.dimension_0()) i = 0; i < v.dimension_0(); i++) {
          for(decltype(v.dimension_0()) j = 0; j < v.dimension_0(); j++) {
            vmuldata(i,j) = -2. * v(i) * v(j);
          }
        }
        for(decltype(v.dimension_0()) i = 0; i < v.dimension_0(); i++) {
          vmuldata(i,i) += one;
        }
      }

      /// \brief Calculate 2-norm of vector
      /// \param v[in] input vector
      /// \ret norm[out] 2-norm of input vector
      KOKKOS_INLINE_FUNCTION
      SC vnorm(const shared_vector & v) const {
        SC sum = 0;
        for(decltype(v.dimension_0()) i = 0; i < v.dimension_0(); i++)
          sum += v(i) * v(i);
        return sqrt(sum);
      }

      /// \brief Scales input vector v by d: v = v/d
      /// \param v[in] input vector
      /// \param d[in] scalar scaling factor
      KOKKOS_INLINE_FUNCTION
      void vdiv(shared_vector & v, SC d) const {
        for(decltype(v.dimension_0()) i = 0; i < v.dimension_0(); i++)
          v(i) = v(i) / d;
      }

      /// \brief Add vector add to v: v = d * v + va
      /// \param v[in] input vector
      /// \param va[in] input vector
      /// \param d[in] scalar scaling factor
      KOKKOS_INLINE_FUNCTION
      void vmadd(shared_vector & v, const shared_vector & va, SC d) const {
        for(decltype(v.dimension_0()) i = 0; i < v.dimension_0(); i++)
          v(i) = va(i) + d * v(i);
      }

      // amout of shared memory
      size_t team_shmem_size( int team_size ) const {
        return 3 * Kokkos::View<double**,Kokkos::MemoryUnmanaged>::shmem_size(maxAggDofSize,fineNS.dimension_1()) + // mat + matminor + z
               3 * Kokkos::View<double**,Kokkos::MemoryUnmanaged>::shmem_size(maxAggDofSize,maxAggDofSize) +  // qk and q and qt
               2 * Kokkos::View<double*,Kokkos::MemoryUnmanaged>::shmem_size(maxAggDofSize); // e, x
      }
    };

  } // namespace anonymous

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
      BuildPuncoupled(coarseLevel, A, aggregates, amalgInfo, fineNullspace, coarseMap, Ptentative, coarseNullspace);
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
  BuildPuncoupled(Level& coarseLevel, RCP<Matrix> A, RCP<Aggregates_kokkos> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                  RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace) const {
    auto rowMap = A->getRowMap();
    auto colMap = A->getColMap();

    const size_t numRows  = rowMap->getNodeNumElements();
    const size_t NSDim    = fineNullspace->getNumVectors();

    typedef Kokkos::ArithTraits<SC>     ATS;
    typedef typename ATS::magnitudeType Magnitude;
    const SC zero = ATS::zero();

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
    auto     procWinner    = aggregates->GetProcWinner()  ->template getLocalView<DeviceType>();
    auto     vertex2AggId  = aggregates->GetVertex2AggId()->template getLocalView<DeviceType>();
    const GO numAggregates = aggregates->GetNumAggregates();

    int myPid = aggregates->GetMap()->getComm()->getRank();

    // Create Kokkos::View (on the device) to store the aggreate dof size
    // Later used to get aggregate dof offsets
    // NOTE: This zeros itself on construction
    typedef typename Aggregates_kokkos::aggregates_sizes_type::non_const_type AggSizeType;
    AggSizeType sizes;

    if (stridedBlockSize == 1) {
      SubFactoryMonitor m2(*this, "Calc AggSizes", coarseLevel);
      // FIXME_KOKKOS: use ViewAllocateWithoutInitializing + set a single value
      sizes = AggSizeType("agg_dof_sizes", numAggregates+1);

      auto sizesConst = aggregates->ComputeAggregateSizes();
      Kokkos::deep_copy(Kokkos::subview(sizes, Kokkos::make_pair(1, numAggregates+1)), sizesConst);

    } else {
      // FIXME_KOKKOS: This branch of code is completely unoptimized.

      sizes = AggSizeType("agg_dof_sizes", numAggregates + 1);
      typename AppendTrait<AggSizeType, Kokkos::Atomic>::type sizesAtomic = sizes; // atomic access

      Teuchos::ArrayView<const GO> nodeGlobalEltsView = aggregates->GetMap()->getNodeElementList();
      Kokkos::View<GO*, DeviceType> nodeGlobalElts("nodeGlobalElts", nodeGlobalEltsView.size());
      {
        SubFactoryMonitor m2(*this, "Create isNodeGlobalElement", coarseLevel);
        // Store overlapping local node ids belonging to aggregates on the current processor in a view
        // This has to be done in serial on the host
        for(size_t i = 0; i < nodeGlobalElts.size(); i++)
          nodeGlobalElts(i) = nodeGlobalEltsView[i];
      }

      // Create an Kokkos::UnorderedMap to store the mapping globalDofIndex -> bool (isNodeGlobalElement)
      // We want to use that information within parallel kernels and cannot call Xpetra::Map routines in the
      // parallel kernels
      typedef Kokkos::UnorderedMap<LO, bool, DeviceType> map_type;
      map_type isNodeGlobalElement(colMap->getNodeNumElements());

      {
        SubFactoryMonitor m2(*this, "Create isNodeGlobalElement", coarseLevel);

        // create a unordered map GID -> isGlobalElement in colMap of A (available from above)
        // This has to be done on the host
        for (decltype(vertex2AggId.dimension_0()) lnode = 0; lnode < vertex2AggId.dimension_0(); lnode++) {
          if(procWinner(lnode,0) == myPid) {
            GO gnodeid = nodeGlobalElts[lnode];
            for (LO k = 0; k< stridedBlockSize; k++) {
              GO gDofIndex = globalOffset + (gnodeid - indexBase) * fullBlockSize + stridingOffset + k + indexBase;
              bool bIsInColumnMap = colMap->isNodeGlobalElement(gDofIndex);
              isNodeGlobalElement.insert(gDofIndex, bIsInColumnMap);
            }
          }
        }
      }

      {
        SubFactoryMonitor m2(*this, "Calc AggSizes", coarseLevel);

        Kokkos::parallel_for("MueLu:TentativePF:Build:aggSizes", range_type(0, vertex2AggId.dimension_0()),
          KOKKOS_LAMBDA(const LO lnode) {
            if (procWinner(lnode,0) == myPid) {
              auto myAgg = vertex2AggId(lnode,0);
              GO gnodeid = nodeGlobalElts(lnode);
              for (LO k = 0; k< stridedBlockSize; k++) {
                GO gDofIndex = globalOffset + (gnodeid - indexBase) * fullBlockSize + stridingOffset + k + indexBase;
                if (isNodeGlobalElement.value_at(isNodeGlobalElement.find(gDofIndex)) == true) {
                  sizesAtomic(myAgg+1)++;
                }
              }
            }
          });
      }
    }

    // Find maximum dof size for aggregates
    // Later used to reserve enough scratch space for local QR decompositions
    LO maxAggSize = 0;
    ReduceMaxFunctor<LO,decltype(sizes)> reduceMax(sizes);
    Kokkos::parallel_reduce("MueLu:TentativePF:Build:max_agg_size", range_type(0, sizes.dimension_0()), reduceMax, maxAggSize);

    // parallel_scan (exclusive)
    // The sizes View then contains the aggregate dof offsets
    Kokkos::parallel_scan("MueLu:TentativePF:Build:aggregate_sizes:stage1_scan", range_type(0,numAggregates+1),
      KOKKOS_LAMBDA(const LO i, LO& update, const bool& final_pass) {
        update += sizes(i);
        if (final_pass)
          sizes(i) = update;
      });

    // Create Kokkos::View on the device to store mapping between (local) aggregate id and row map ids (LIDs)
    // FIXME_KOKKOS: Does it need to be this complicated for scalar?
    Kokkos::View<LO*, DeviceType> agg2RowMapLO(Kokkos::ViewAllocateWithoutInitializing("agg2row_map_LO"), numRows);
    {
      // FIXME_KOKKOS: This chunk of code is completely unoptimized.

      SubFactoryMonitor m2(*this, "Create Agg2RowMap", coarseLevel);

      // We need this initialized to 0, so no ViewAlloWithoutInitializing
      AggSizeType aggDofCount("aggDofCount", numAggregates);

      Kokkos::parallel_for("MueLu:TentativePF:Build:createAgg2RowMap", range_type(0, vertex2AggId.dimension_0()),
        KOKKOS_LAMBDA(const LO lnode) {
          if (stridedBlockSize == 1) {
            if (procWinner(lnode,0) == myPid) {
              auto myAgg = vertex2AggId(lnode,0);
              LO myAggDofStart = sizes(myAgg);
              auto idx = Kokkos::atomic_fetch_add( &aggDofCount(myAgg), 1 );
              agg2RowMapLO(myAggDofStart + idx) = lnode;
            }
          } else {
            if (procWinner(lnode,0) == myPid) {
              auto myAgg = vertex2AggId(lnode,0);
              LO myAggDofStart = sizes(myAgg);
              auto idx = Kokkos::atomic_fetch_add( &aggDofCount(myAgg), stridedBlockSize );
              for (LO k = 0; k < stridedBlockSize; k++) {
                agg2RowMapLO(myAggDofStart + idx + k) = lnode * stridedBlockSize +k;
              }
            }
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

    TEUCHOS_TEST_FOR_EXCEPTION(goodMap == false, Exceptions::RuntimeError,
        "Only works for non-overlapping aggregates (goodMap == true)");

    rows_type rows;
    cols_type cols;
    vals_type vals;

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

      Kokkos::parallel_for("MueLu:TentativePF:BuildUncoupled:main_loop", policy,
        KOKKOS_LAMBDA(const typename Kokkos::TeamPolicy<execution_space>::member_type &thread) {
          auto agg = thread.league_rank();

          // size of the aggregate (number of DOFs in aggregate)
          LO aggSize = aggRows(agg+1) - aggRows(agg);

          // Extract the piece of the nullspace corresponding to the aggregate, and
          // put it in the flat array, "localQR" (in column major format) for the
          // QR routine. Trivial in 1D.
          if (goodMap) {
            // Calculate QR by hand
            auto norm = ATS::magnitude(zero);
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

          } else {
            // FIXME_KOKKOS: implement non-standard map QR
            // Look at the original TentativeP for how to do that
            statusAtomic(0) = true;
            return;
          }
        });

        typename status_type::HostMirror statusHost = Kokkos::create_mirror_view(status);
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
          SubFactoryMonitor m2(*this, "Stage 1 (LocalQR)", coarseLevel);
          // Set up team policy with numAggregates teams and one thread per team.
          // Each team handles a slice of the data associated with one aggregate
          // and performs a local QR decomposition
          const Kokkos::TeamPolicy<execution_space> policy(numAggregates,1); // numAggregates teams a 1 thread
          LocalQRDecompFunctor<LocalOrdinal, GlobalOrdinal, Scalar, DeviceType, decltype(fineNSRandom),
              decltype(sizes /*aggregate sizes in dofs*/), decltype(maxAggSize), decltype(agg2RowMapLO),
              decltype(statusAtomic), decltype(rows), decltype(rowsAux), decltype(colsAux),
              decltype(valsAux)>
                  localQRFunctor(fineNSRandom, coarseNS, sizes, maxAggSize, agg2RowMapLO, statusAtomic,
                                 rows, rowsAux, colsAux, valsAux);
          Kokkos::parallel_reduce("MueLu:TentativePF:BuildUncoupled:main_qr_loop", policy, localQRFunctor, nnz);
        }

        typename status_type::HostMirror statusHost = Kokkos::create_mirror_view(status);
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
          Kokkos::parallel_for("MueLu:TentativePF:Build:compress_cols_vals", numRows,
            KOKKOS_LAMBDA(const LO i) {
              LO rowStart = rows(i);

              size_t lnnz = 0;
              for (typename decltype(cols)::value_type j = rowsAux(i); j < rowsAux(i+1); j++)
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
#if 1
      // FIXME_KOKKOS: this should be gone once Xpetra propagate "local matrix + 4 maps" constructor
      auto PtentCrs = CrsMatrixFactory::Build(rowMap, coarseMap, lclMatrix);
      PtentCrs->resumeFill();  // we need that for rectangular matrices
      PtentCrs->expertStaticFillComplete(coarseMap, A->getDomainMap());
#else
      auto PtentCrs = CrsMatrixFactory::Build(lclMatrix, rowMap, coarseMap, coarseMap, A->getDomainMap());
#endif
      Ptentative = rcp(new CrsMatrixWrap(PtentCrs));
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::
  BuildPcoupled(RCP<Matrix> A, RCP<Aggregates_kokkos> aggregates, RCP<AmalgamationInfo> amalgInfo, RCP<MultiVector> fineNullspace,
                RCP<const Map> coarseMap, RCP<Matrix>& Ptentative, RCP<MultiVector>& coarseNullspace) const {
    throw Exceptions::RuntimeError("Construction of coupled tentative P is not implemented");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool TentativePFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::
  isGoodMap(const Map& rowMap, const Map& colMap) const {
    bool goodMap = true;
    if (rowMap.lib() == Xpetra::UseEpetra) {
      // Epetra version
      // As Xpetra::Map is missing getLocalMap(), implement a serial variant
      ArrayView<const GO> rowElements = rowMap.getNodeElementList();
      ArrayView<const GO> colElements = colMap.getNodeElementList();

      const size_t numElements = rowElements.size();

      for (size_t i = 0; i < numElements; i++)
        if (rowElements[i] != colElements[i]) {
          goodMap = false;
          break;
        }
    } else {
      // Tpetra version

      // FIXME_KOKKOS: move this map check to Tpetra::CrsGraph
      auto rowTMap = Utilities::Map2TpetraMap(rowMap);
      auto colTMap = Utilities::Map2TpetraMap(colMap);

      auto rowLocalMap = rowTMap->getLocalMap();
      auto colLocalMap = colTMap->getLocalMap();

      const size_t numRows = rowLocalMap.getNodeNumElements();
      const size_t numCols = colLocalMap.getNodeNumElements();

      if (numCols < numRows)
        return false;

      size_t numDiff = 0;
      Kokkos::parallel_reduce("MueLu:TentativePF:isGoodMap", range_type(0, numRows),
        KOKKOS_LAMBDA(const size_t i, size_t &diff) {
          diff += (rowLocalMap.getGlobalElement(i) != colLocalMap.getGlobalElement(i));
        }, numDiff);

      goodMap = (numDiff == 0);
    }

    return goodMap;
  }

} //namespace MueLu

#define MUELU_TENTATIVEPFACTORY_KOKKOS_SHORT
#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_TENTATIVEPFACTORY_KOKKOS_DEF_HPP
