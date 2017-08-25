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

    // only limited support for lambdas with CUDA 7.5
    // see https://github.com/kokkos/kokkos/issues/763
    // This functor probably can go away when using CUDA 8.0
    template<class LocalOrdinal, class rowType, class NSDimType>
    class FillRowAuxArrayFunctor {
    public:
      FillRowAuxArrayFunctor(rowType rows, NSDimType nsDim) : rows_(rows), nsDim_(nsDim) { }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LocalOrdinal row) const {
        rows_(row) = row*nsDim_;
      }

    private:
      rowType   rows_;
      NSDimType nsDim_;
    };

    // only limited support for lambdas with CUDA 7.5
    // see https://github.com/kokkos/kokkos/issues/763
    // This functor probably can go away when using CUDA 8.0
    template<class SC, class LO, class colType, class valType>
    class FillColAuxArrayFunctor {
    public:
      FillColAuxArrayFunctor(SC zero, LO invalid, colType cols, valType vals) : zero_(zero), invalid_(invalid), cols_(cols), vals_(vals) { }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO col) const {
        cols_(col) = invalid_;
        vals_(col) = zero_;
      }

    private:
      SC        zero_;
      LO        invalid_;
      colType   cols_;
      valType   vals_;
    };

    // collect aggregate sizes (number of dofs associated with all nodes in aggregate)
    // aggSizesType needs to be atomic
    template<class aggSizesType, class vertex2AggIdType, class procWinnerType, class nodeGlobalEltsType, class isNodeGlobalEltsType, class LOType, class GOType>
    class CollectAggregateSizeFunctor {
    private:
      typedef LOType LO;
      typedef GOType GO;

      aggSizesType         aggSizes;        //< view containint size of aggregate
      vertex2AggIdType     vertex2AggId;    //< view containing vertex2AggId information
      procWinnerType       procWinner;      //< view containing processor ids
      nodeGlobalEltsType   nodeGlobalElts;  //< view containing global node ids of current proc
      isNodeGlobalEltsType isNodeGlobalElement; //< unordered map storing whether (global) node id is owned by current proc

      LO fullBlockSize, blockID, stridingOffset, stridedBlockSize;
      GO indexBase, globalOffset;
      int myPid;
    public:
      CollectAggregateSizeFunctor(aggSizesType aggSizes_, vertex2AggIdType vertex2AggId_, procWinnerType procWinner_, nodeGlobalEltsType nodeGlobalElts_, isNodeGlobalEltsType isNodeGlobalElement_, LO fullBlockSize_, LOType blockID_, LOType stridingOffset_, LOType stridedBlockSize_, GOType indexBase_, GOType globalOffset_, int myPID_ ) :
        aggSizes(aggSizes_),
        vertex2AggId(vertex2AggId_),
        procWinner(procWinner_),
        nodeGlobalElts(nodeGlobalElts_),
        isNodeGlobalElement(isNodeGlobalElement_),
        fullBlockSize(fullBlockSize_),
        blockID(blockID_),
        stridingOffset(stridingOffset_),
        stridedBlockSize(stridedBlockSize_),
        indexBase(indexBase_),
        globalOffset(globalOffset_),
        myPid(myPID_) {
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO lnode) const {
        if(stridedBlockSize == 1) {
          if(procWinner(lnode,0) == myPid) {
            auto myAgg = vertex2AggId(lnode,0);
            aggSizes(myAgg+1)++; // atomic access
          }
        } else {
          if(procWinner(lnode,0) == myPid) {
            auto myAgg = vertex2AggId(lnode,0);
            GO gnodeid = nodeGlobalElts(lnode);
            for (LO k = 0; k< stridedBlockSize; k++) {
              GO gDofIndex = globalOffset + (gnodeid - indexBase) * fullBlockSize + stridingOffset + k + indexBase;
              if(isNodeGlobalElement.value_at(isNodeGlobalElement.find(gDofIndex)) == true) {
                aggSizes(myAgg+1)++; // atomic access
              }
            }
          }
        }
      }
    };

    //      CreateAgg2RowMapLOFunctor<decltype(agg2RowMapLO), decltype(sizes), decltype(vertex2AggId), decltype(procWinner), decltype(nodeGlobalElts), decltype(isNodeGlobalElement), LO, GO> createAgg2RowMap (agg2RowMapLO, sizes, vertex2AggId, procWinner, nodeGlobalElts, isNodeGlobalElement, fullBlockSize, blockID, stridingOffset, stridedBlockSize, indexBase, globalOffset, myPid );

    template<class agg2RowMapType, class aggSizesType, class vertex2AggIdType, class procWinnerType, class nodeGlobalEltsType, class isNodeGlobalEltsType, class LOType, class GOType>
    class CreateAgg2RowMapLOFunctor {
    private:
      typedef LOType LO;
      typedef GOType GO;

      agg2RowMapType       agg2RowMap;      //< view containing row map entries associated with aggregate
      aggSizesType         aggSizes;        //< view containing size of aggregate
      aggSizesType         aggDofCount;     //< view containing current count of "found" vertices in aggregate
      vertex2AggIdType     vertex2AggId;    //< view containing vertex2AggId information
      procWinnerType       procWinner;      //< view containing processor ids
      nodeGlobalEltsType   nodeGlobalElts;  //< view containing global node ids of current proc
      isNodeGlobalEltsType isNodeGlobalElement; //< unordered map storing whether (global) node id is owned by current proc

      LO fullBlockSize, blockID, stridingOffset, stridedBlockSize;
      GO indexBase, globalOffset;
      int myPid;
    public:
      CreateAgg2RowMapLOFunctor(agg2RowMapType agg2RowMap_, aggSizesType aggSizes_, aggSizesType aggDofCount_, vertex2AggIdType vertex2AggId_, procWinnerType procWinner_, nodeGlobalEltsType nodeGlobalElts_, isNodeGlobalEltsType isNodeGlobalElement_, LO fullBlockSize_, LOType blockID_, LOType stridingOffset_, LOType stridedBlockSize_, GOType indexBase_, GOType globalOffset_, int myPID_ ) :
        agg2RowMap(agg2RowMap_),
        aggSizes(aggSizes_),
        aggDofCount(aggDofCount_),
        vertex2AggId(vertex2AggId_),
        procWinner(procWinner_),
        nodeGlobalElts(nodeGlobalElts_),
        isNodeGlobalElement(isNodeGlobalElement_),
        fullBlockSize(fullBlockSize_),
        blockID(blockID_),
        stridingOffset(stridingOffset_),
        stridedBlockSize(stridedBlockSize_),
        indexBase(indexBase_),
        globalOffset(globalOffset_),
        myPid(myPID_) {
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO lnode) const {
        if(stridedBlockSize == 1) {
          if(procWinner(lnode,0) == myPid) {
            auto myAgg = vertex2AggId(lnode,0);
            LO myAggDofStart = aggSizes(myAgg);
            auto idx = Kokkos::atomic_fetch_add( &aggDofCount(myAgg), 1 );
            agg2RowMap(myAggDofStart + idx) = lnode;
          }
        } else {
          if(procWinner(lnode,0) == myPid) {
            auto myAgg = vertex2AggId(lnode,0);
            LO myAggDofStart = aggSizes(myAgg);
            auto idx = Kokkos::atomic_fetch_add( &aggDofCount(myAgg), stridedBlockSize );
            for (LO k = 0; k< stridedBlockSize; k++) {
              agg2RowMap(myAggDofStart + idx + k) = lnode * stridedBlockSize +k;
            }
          }
        }
      }
    };

    // local QR decomposition (scalar case, no real QR decomposition necessary)
    // Use exactly the same template types. The only one that is not really necessary is the maxAggDofSizeType
    // for reserving temporary scratch space
    template<class LOType, class GOType, class SCType,class DeviceType, class NspType, class aggRowsType, class maxAggDofSizeType, class agg2RowMapLOType, class statusType, class rowsType, class rowsAuxType, class colsAuxType, class valsAuxType>
    class LocalScalarQRDecompFunctor {
    private:
      typedef LOType LO;
      typedef GOType GO;
      typedef SCType SC;

      typedef typename DeviceType::execution_space execution_space;

    private:

      NspType fineNSRandom;
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
      LocalScalarQRDecompFunctor(NspType fineNSRandom_, NspType coarseNS_, aggRowsType aggRows_, maxAggDofSizeType maxAggDofSize_, agg2RowMapLOType agg2RowMapLO_, statusType statusAtomic_, rowsType rows_, rowsAuxType rowsAux_, colsAuxType colsAux_, valsAuxType valsAux_) :
        fineNSRandom(fineNSRandom_),
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
      void operator() ( const typename Kokkos::TeamPolicy<execution_space>::member_type & thread, size_t& rowNnz) const {
        auto agg = thread.league_rank();

        // size of aggregate: number of DOFs in aggregate
        LO aggSize = aggRows(agg+1) - aggRows(agg);

        // Extract the piece of the nullspace corresponding to the aggregate, and
        // put it in the flat array, "localQR" (in column major format) for the
        // QR routine. Trivial in 1D.

        // Calculate QR by hand
        typedef Kokkos::ArithTraits<SC>     ATS;
        typedef typename ATS::magnitudeType Magnitude;

        Magnitude norm = Kokkos::ArithTraits<Magnitude>::zero();
        for (decltype(aggSize) k = 0; k < aggSize; k++) {
          Magnitude dnorm = ATS::magnitude(fineNSRandom(agg2RowMapLO(aggRows(agg)+k),0));
          norm += dnorm*dnorm;
        }
        norm = sqrt(norm);

        if (norm == ATS::zero()) {
          // zero column; terminate the execution
          statusAtomic(1) = true;
          return;
        }

        // R = norm
        coarseNS(agg, 0) = norm;

        // Q = localQR(:,0)/norm
        for (LO k = 0; k < aggSize; k++) {
          LO localRow = agg2RowMapLO(aggRows(agg)+k);
          SC localVal = fineNSRandom(agg2RowMapLO(aggRows(agg)+k),0) / norm;

          size_t rowStart = rowsAux(localRow);
          colsAux(rowStart) = agg;
          valsAux(rowStart) = localVal;

          // Store true number of nonzeros per row
          rows(localRow+1) = 1;
          rowNnz          += 1;
        }
      }
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
        for(decltype(m1.dimension_0()) i = 0; i < m1.dimension_0(); i++) {
          for(decltype(m1.dimension_1()) j = 0; j < m1.dimension_1(); j++) {
            m1m2(i,j) = zero;
            for(decltype(m1.dimension_1()) k = 0; k < m1.dimension_1(); k++) {
              m1m2(i,j) += m1(i,k) * m2(k,j);
            }
          }
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


    // only limited support for lambdas with CUDA 7.5
    // see https://github.com/kokkos/kokkos/issues/763
    // This functor probably can go away when using CUDA 8.0
    template<class LO, class rowType, class colType, class valType>
    class CompressArraysFunctor {
    public:
      CompressArraysFunctor(LO invalid, rowType rowsAux, colType colsAux, valType valsAux, rowType rows, colType cols, valType vals) : invalid_(invalid), rowsAux_(rowsAux), colsAux_(colsAux), valsAux_(valsAux), rows_(rows), cols_(cols), vals_(vals) { }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO i) const {
        LO rowStart = rows_(i);

        size_t lnnz = 0;
        for (decltype(rowsAux_(i)) j = rowsAux_(i); j < rowsAux_(i+1); j++)
          if (colsAux_(j) != invalid_) {
            cols_(rowStart+lnnz) = colsAux_(j);
            vals_(rowStart+lnnz) = valsAux_(j);
            lnnz++;
          }
      }

    private:
      LO        invalid_;
      rowType   rowsAux_;
      colType   colsAux_;
      valType   valsAux_;
      rowType   rows_;
      colType   cols_;
      valType   vals_;
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
    RCP<const Map> rowMap = A->getRowMap();
    RCP<const Map> colMap = A->getColMap();

    const size_t numRows  = rowMap->getNodeNumElements();
    const size_t NSDim    = fineNullspace->getNumVectors();

    typedef Teuchos::ScalarTraits<SC> STS;
    const SC zero    = STS::zero();
    const SC one     = STS::one();
    const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();

    auto     aggGraph = aggregates->GetGraph();
    auto     aggRows  = aggGraph.row_map;
    auto     aggCols  = aggGraph.entries;
    //const GO numAggs  = aggregates->GetNumAggregates();

    // Aggregates map is based on the amalgamated column map
    // We can skip global-to-local conversion if LIDs in row map are
    // same as LIDs in column map
    bool goodMap = isGoodMap(*rowMap, *colMap);

    // STEP 1: do unamalgamation
    // In contrast to the non-kokkos version which uses member functions from the AmalgamationInfo container
    // class to unamalgamate the data. The kokkos version of TentativePFacotry does the unamalgamation here
    // and only uses the data of the AmalgamationInfo container class

    // Extract information for unamalgamation
    LO fullBlockSize, blockID, stridingOffset, stridedBlockSize;
    GO indexBase;
    amalgInfo->GetStridingInformation(fullBlockSize, blockID, stridingOffset, stridedBlockSize, indexBase);
    GO globalOffset = amalgInfo->GlobalOffset();

    // Store overlapping local node ids belonging to aggregates on the current processor in a view
    // This has to be done in serial on the host
    Teuchos::ArrayView<const GO> nodeGlobalEltsView = aggregates->GetMap()->getNodeElementList();
    Kokkos::View<GO*, DeviceType> nodeGlobalElts("nodeGlobalElts", nodeGlobalEltsView.size());
    for(size_t i = 0; i < nodeGlobalElts.size(); i++)
      nodeGlobalElts(i) = nodeGlobalEltsView[i];

    // Extract aggregation info (already in Kokkos host views)
    auto procWinner   = aggregates->GetProcWinner()->getHostLocalView();
    auto vertex2AggId = aggregates->GetVertex2AggId()->getHostLocalView();
    const GO numAggregates = aggregates->GetNumAggregates();

    // Create an Kokkos::UnorderedMap to store the mapping globalDofIndex -> bool (isNodeGlobalElement)
    // We want to use that information within parallel kernels and cannot call Xpetra::Map routines in the
    // parallel kernels
    typedef Kokkos::UnorderedMap<LO, bool, DeviceType> map_type;
    map_type isNodeGlobalElement(colMap->getNodeNumElements());

    int myPid = aggregates->GetMap()->getComm()->getRank();

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

    // Create Kokkos::View (on the device) to store the aggreate dof size
    // Later used to get aggregate dof offsets
    // NOTE: This zeros itself on construction
    typedef Kokkos::View<LO*, DeviceType> aggSizeType;
    aggSizeType sizes("agg_dof_sizes", numAggregates + 1);
    //    sizes(0) = 0;

#if 1

    {
      SubFactoryMonitor m2(*this, "Calc AggSizes", coarseLevel);

      // atomic data access for sizes view
      typename AppendTrait<aggSizeType, Kokkos::Atomic>::type sizesAtomic = sizes;

      // loop over all nodes
      CollectAggregateSizeFunctor <decltype(sizesAtomic), decltype(vertex2AggId), decltype(procWinner), decltype(nodeGlobalElts), decltype(isNodeGlobalElement), LO, GO> collectAggregateSizes(sizesAtomic, vertex2AggId, procWinner, nodeGlobalElts, isNodeGlobalElement, fullBlockSize, blockID, stridingOffset, stridedBlockSize, indexBase, globalOffset, myPid );
      Kokkos::parallel_for("MueLu:TentativePF:Build:getaggsizes", range_type(0,vertex2AggId.dimension_0()), collectAggregateSizes);
    }
#else
    // we have to avoid race conditions when parallelizing the loops below
    if(stridedBlockSize == 1) {
      // loop over all nodes
      for (decltype(vertex2AggId.dimension_0()) lnode = 0; lnode < vertex2AggId.dimension_0(); lnode++) {
        if(procWinner(lnode,0) == myPid) {
          auto myAgg = vertex2AggId(lnode,0);
          sizes(myAgg+1)++;
        }
      }
    } else {
      // loop over all nodes
      for (decltype(vertex2AggId.dimension_0()) lnode = 0; lnode < vertex2AggId.dimension_0(); lnode++) {
        if(procWinner(lnode,0) == myPid) {
          auto myAgg = vertex2AggId(lnode,0);
          GO gnodeid = nodeGlobalElts[lnode];

          for (LO k = 0; k< stridedBlockSize; k++) {
            GO gDofIndex = globalOffset + (gnodeid - indexBase) * fullBlockSize + stridingOffset + k + indexBase;
            if(isNodeGlobalElement.value_at(isNodeGlobalElement.find(gDofIndex)) == true) {
              sizes(myAgg+1)++;
            }
          }
        }
      }
    }
#endif

    // Find maximum dof size for aggregates
    // Later used to reserve enough scratch space for local QR decompositions
    // TODO can be done with a parallel reduce
    LO maxAggSize = 0;
    for(decltype(sizes.dimension_0()) i = 0; i < sizes.dimension_0(); i++) {
      if(sizes(i) > maxAggSize) maxAggSize = sizes(i);
    }

    // parallel_scan (exclusive)
    // The Kokkos::View sizes then contains the aggreate Dof offsets
    ScanFunctor<LO,decltype(sizes)> scanFunctorAggSizes(sizes);
    Kokkos::parallel_scan("MueLu:TentativePF:Build:aggregate_sizes:stage1_scan", range_type(0,numAggregates+1), scanFunctorAggSizes);
    // Create Kokkos::View on the device to store mapping between (local) aggregate id and row map ids (LIDs)
    Kokkos::View<LO*, DeviceType> agg2RowMapLO("agg2row_map_LO", numRows); // initialized to 0

#if 1

    {
      SubFactoryMonitor m2(*this, "Create Agg2RowMap", coarseLevel);
      // NOTE: This zeros itself on construction
      aggSizeType aggDofCount("aggDofCount", numAggregates);


      // atomic data access for agg2RowMapLO view
      //typename AppendTrait<decltype(agg2RowMapLO), Kokkos::Atomic>::type agg2RowMapLOAtomic = agg2RowMapLO;

      CreateAgg2RowMapLOFunctor<decltype(agg2RowMapLO), decltype(sizes), decltype(vertex2AggId), decltype(procWinner), decltype(nodeGlobalElts), decltype(isNodeGlobalElement), LO, GO> createAgg2RowMap (agg2RowMapLO, sizes, aggDofCount,vertex2AggId, procWinner, nodeGlobalElts, isNodeGlobalElement, fullBlockSize, blockID, stridingOffset, stridedBlockSize, indexBase, globalOffset, myPid );
      Kokkos::parallel_for("MueLu:TentativePF:Build:createAgg2RowMap", range_type(0,vertex2AggId.dimension_0()), createAgg2RowMap);
    }
#else
    Kokkos::View<LO*, DeviceType> numDofs("countDofsPerAggregate", numAggregates); // helper view. We probably don't need that
    if(stridedBlockSize == 1) {
      // loop over all nodes
      for (decltype(vertex2AggId.dimension_0()) lnode = 0; lnode < vertex2AggId.dimension_0(); lnode++) {
        if(procWinner(lnode,0) == myPid) {
          auto myAgg = vertex2AggId(lnode,0);
          agg2RowMapLO(sizes(myAgg) + numDofs(myAgg)) = lnode;
          numDofs(myAgg)++;
        }
      }
    } else {
      // loop over all nodes
      for (decltype(vertex2AggId.dimension_0()) lnode = 0; lnode < vertex2AggId.dimension_0(); lnode++) {
        if(procWinner(lnode,0) == myPid) {
          auto myAgg = vertex2AggId(lnode,0);
          GO gnodeid = nodeGlobalElts[lnode];

          for (LO k = 0; k< stridedBlockSize; k++) {
            GO gDofIndex = globalOffset + (gnodeid - indexBase) * fullBlockSize + stridingOffset + k + indexBase;
            if(isNodeGlobalElement.value_at(isNodeGlobalElement.find(gDofIndex)) == true) {
              agg2RowMapLO(sizes(myAgg) + numDofs(myAgg)) = lnode * stridedBlockSize + k;
              numDofs(myAgg)++;
            }
          }
        }
      }
    }
#endif

    // STEP 2: prepare local QR decomposition
    // Reserve memory for tentative prolongation operator

    coarseNullspace = MultiVectorFactory::Build(coarseMap, NSDim);

    // Pull out the nullspace vectors so that we can have random access (on the device)
    auto fineNS   = fineNullspace  ->template getLocalView<DeviceType>();
    auto coarseNS = coarseNullspace->template getLocalView<DeviceType>();

    size_t nnzEstimate = numRows * NSDim; // maximum number of possible nnz
    size_t nnz = 0;                       // actual number of nnz

    typedef typename Xpetra::Matrix<SC,LO,GO,NO>::local_matrix_type          local_matrix_type;
    typedef typename local_matrix_type::row_map_type    rows_type;
    typedef typename local_matrix_type::index_type      cols_type;
    typedef typename local_matrix_type::values_type     vals_type;

    // Stage 0: initialize auxiliary arrays
    // The main thing to notice is initialization of vals with INVALID. These
    // values will later be used to compress the arrays
    typename rows_type::non_const_type rowsAux("Ptent_aux_rows", numRows+1),    rows("Ptent_rows", numRows+1);
    typename cols_type::non_const_type colsAux("Ptent_aux_cols", nnzEstimate);
    typename vals_type::non_const_type valsAux("Ptent_aux_vals", nnzEstimate);

#if 1
    // only limited support for lambdas with CUDA 7.5
    // see https://github.com/kokkos/kokkos/issues/763
    FillRowAuxArrayFunctor<LO, decltype(rowsAux), decltype(NSDim)> rauxf(rowsAux, NSDim);
    Kokkos::parallel_for("MueLu:TentativePF:BuildPuncoupled:for1", range_type(0,numRows+1), rauxf);
    FillColAuxArrayFunctor<SC, LO, decltype(colsAux), decltype(valsAux)> cauxf(zero, INVALID, colsAux, valsAux);
    Kokkos::parallel_for("MueLu:TentativePF:BuildUncoupled:for2", range_type(0,nnzEstimate), cauxf);
#else
    Kokkos::parallel_for("MueLu:TentativePF:BuildPuncoupled:for1", range_type(0,numRows+1), KOKKOS_LAMBDA(const LO row) {
      rowsAux(row) = row*NSDim;
    });
    Kokkos::parallel_for("MueLu:TentativePF:BuildUncoupled:for2", nnzEstimate, KOKKOS_LAMBDA(const LO j) {
     colsAux(j) = INVALID;
     valsAux(j) = zero;
    });
#endif

    // Device View for status (error messages...)
    typedef Kokkos::View<int[10], DeviceType> status_type;
    status_type status("status");

    // Stage 1: construct auxiliary arrays.
    // The constructed arrays may have gaps in them (vals(j) == INVALID)
    // Run one thread per aggregate.
    typename AppendTrait<decltype(fineNS), Kokkos::RandomAccess>::type fineNSRandom = fineNS;
    typename AppendTrait<status_type,      Kokkos::Atomic>      ::type statusAtomic = status;

    TEUCHOS_TEST_FOR_EXCEPTION(goodMap == false, Exceptions::RuntimeError, "Only works for non-overlapping aggregates (goodMap == true)");

    {
      SubFactoryMonitor m2(*this, "Stage 1 (LocalQR)", coarseLevel);

      if (NSDim == 1) {

  #if 1
        // only limited support for lambdas with CUDA 7.5
        // see https://github.com/kokkos/kokkos/issues/763
        //
        // Set up team policy with numAggregates teams and one thread per team.
        // Each team handles a slice of the data associated with one aggregate
        // and performs a local QR decomposition (in this case no real QR decomp
        // is necessary)
        const Kokkos::TeamPolicy<execution_space> policy(numAggregates,1);
        Kokkos::parallel_reduce(policy, LocalScalarQRDecompFunctor<LocalOrdinal, GlobalOrdinal, Scalar, DeviceType, decltype(fineNSRandom), decltype(sizes /*aggregate sizes in dofs*/), decltype(maxAggSize), decltype(agg2RowMapLO), decltype(statusAtomic), decltype(rows), decltype(rowsAux), decltype(colsAux), decltype(valsAux)>(fineNSRandom,coarseNS,sizes,maxAggSize,agg2RowMapLO,statusAtomic,rows,rowsAux,colsAux,valsAux),nnz);
  #else
        // Andrey's original implementation as a lambda
        // 1D is special, as it is the easiest. We don't even need to the QR,
        // just normalize an array. Plus, no worries abot small aggregates.
        Kokkos::parallel_reduce("MueLu:TentativePF:BuildUncoupled:main_loop", numAggregates, KOKKOS_LAMBDA(const GO agg, size_t& rowNnz) {
          LO aggSize = aggRows(agg+1) - aggRows(agg);

          // Extract the piece of the nullspace corresponding to the aggregate, and
          // put it in the flat array, "localQR" (in column major format) for the
          // QR routine. Trivial in 1D.
          if (goodMap) {
            // Calculate QR by hand
            typedef Kokkos::ArithTraits<SC>     ATS;
            typedef typename ATS::magnitudeType Magnitude;

            Magnitude norm = ATS::magnitude(zero);
            for (size_t k = 0; k < aggSize; k++) {
              Magnitude dnorm = ATS::magnitude(fineNSRandom(agg2RowMapLO(aggRows(agg)+k),0));
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
            for (LO k = 0; k < aggSize; k++) {
              LO localRow = agg2RowMapLO(aggRows(agg)+k);
              SC localVal = fineNSRandom(agg2RowMapLO(aggRows(agg)+k),0) / norm;

              size_t rowStart = rowsAux(localRow);
              colsAux(rowStart) = agg;
              valsAux(rowStart) = localVal;

              // Store true number of nonzeros per row
              rows(localRow+1) = 1;
              rowNnz          += 1;
            }

          } else {
            // FIXME: implement non-standard map QR
            // Look at the original TentativeP for how to do that
            statusAtomic(0) = true;
            return;
          }
        }, nnz);
  #endif

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

      } else { // NSdim > 1

        // Set up team policy with numAggregates teams and one thread per team.
        // Each team handles a slice of the data associated with one aggregate
        // and performs a local QR decomposition
        //const Kokkos::TeamPolicy<> policy( numAggregates, 1); // numAggregates teams a 1 thread
        const Kokkos::TeamPolicy<execution_space> policy(numAggregates,1);
        Kokkos::parallel_reduce(policy, LocalQRDecompFunctor<LocalOrdinal, GlobalOrdinal, Scalar, DeviceType, decltype(fineNSRandom), decltype(sizes /*aggregate sizes in dofs*/), decltype(maxAggSize), decltype(agg2RowMapLO), decltype(statusAtomic), decltype(rows), decltype(rowsAux), decltype(colsAux), decltype(valsAux)>(fineNSRandom,coarseNS,sizes,maxAggSize,agg2RowMapLO,statusAtomic,rows,rowsAux,colsAux,valsAux),nnz);

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
      }

    } // subtime monitor stage 1

    {
      SubFactoryMonitor m2(*this, "Stage 2 (CompressRows)", coarseLevel);
      // Stage 2: compress the arrays
      ScanFunctor<LO,decltype(rows)> scanFunctor(rows);
      Kokkos::parallel_scan("MueLu:TentativePF:Build:compress_rows", range_type(0,numRows+1), scanFunctor);
    }

    // The real cols and vals are constructed using calculated (not estimated) nnz
    typename cols_type::non_const_type cols("Ptent_cols", nnz);
    typename vals_type::non_const_type vals("Ptent_vals", nnz);

#if 1
    {
      SubFactoryMonitor m2(*this, "Stage 2 (CompressCols)", coarseLevel);

      CompressArraysFunctor<LO, decltype(rows),decltype(cols),decltype(vals)> comprArr(INVALID, rowsAux, colsAux, valsAux, rows, cols, vals);
      Kokkos::parallel_for("MueLu:TentativePF:Build:compress_cols_vals", range_type(0,numRows), comprArr);
    }
#else
    Kokkos::parallel_for("MueLu:TentativePF:Build:compress_cols_vals", numRows, KOKKOS_LAMBDA(const LO i) {
      LO rowStart = rows(i);

      size_t lnnz = 0;
      for (LO j = rowsAux(i); j < rowsAux(i+1); j++)
        if (colsAux(j) != INVALID) {
          cols(rowStart+lnnz) = colsAux(j);
          vals(rowStart+lnnz) = valsAux(j);
          lnnz++;
        }
    });
#endif

    GetOStream(Runtime1) << "TentativePFactory : aggregates do not cross process boundaries" << std::endl;

    // Stage 3: construct Xpetra::Matrix
#if 1
    // create local CrsMatrix
    {
      SubFactoryMonitor m2(*this, "Stage 3 (LocalMatrix+FillComplete)", coarseLevel);
      local_matrix_type lclMatrix = local_matrix_type("A", numRows, coarseMap->getNodeNumElements(), nnz, vals, rows, cols);
      Teuchos::RCP<CrsMatrix> PtentCrs = CrsMatrixFactory::Build(rowMap,coarseMap,lclMatrix);
      PtentCrs->resumeFill();  // we need that for rectangular matrices
      PtentCrs->expertStaticFillComplete(coarseMap, A->getDomainMap());
      Ptentative = rcp(new CrsMatrixWrap(PtentCrs));
    }
#else
    Ptentative = rcp(new CrsMatrixWrap(rowMap, coarseMap, 0, Xpetra::StaticProfile));
    RCP<CrsMatrix> PtentCrs = rcp_dynamic_cast<CrsMatrixWrap>(Ptentative)->getCrsMatrix();

    ArrayRCP<size_t>  iaPtent;
    ArrayRCP<LO>      jaPtent;
    ArrayRCP<SC>     valPtent;

    PtentCrs->allocateAllValues(nnz, iaPtent, jaPtent, valPtent);
    ArrayView<size_t> ia  = iaPtent();
    ArrayView<LO>     ja  = jaPtent();
    ArrayView<SC>     val = valPtent();

    // Copy values
    typename rows_type::HostMirror rowsHost = Kokkos::create_mirror_view(rows);
    typename cols_type::HostMirror colsHost = Kokkos::create_mirror_view(cols);
    typename vals_type::HostMirror valsHost = Kokkos::create_mirror_view(vals);

    // copy data from device to host
    // shouldn't be do anything if we are already on the host
    Kokkos::deep_copy (rowsHost, rows);
    Kokkos::deep_copy (colsHost, cols);
    Kokkos::deep_copy (valsHost, vals);

    for (LO i = 0; i < rowsHost.size(); i++)
      ia[i] = rowsHost(i);
    for (LO j = 0; j < colsHost.size(); j++) {
      ja [j] = colsHost(j);
      val[j] = valsHost(j);
    }
    PtentCrs->setAllValues(iaPtent, jaPtent, valPtent);
    PtentCrs->expertStaticFillComplete(coarseMap, A->getDomainMap());
#endif
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
