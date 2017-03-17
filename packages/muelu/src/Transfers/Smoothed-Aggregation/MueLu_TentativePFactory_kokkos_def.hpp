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
//#include <Teuchos_SerialDenseMatrix.hpp>   // TODO remove me
//#include <Teuchos_SerialQRDenseSolver.hpp> // TODO remove me

#include <Teuchos_LAPACK.hpp> // TODO remove me

//#include <Tsqr_Matrix.hpp>
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

    // collect aggregate sizes (number of dofs associated with all nodes in aggregate)
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

      int myPid;
      LO fullBlockSize, blockID, stridingOffset, stridedBlockSize;
      GO indexBase, globalOffset;
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
      void operator()(const LO aggId) const {
        if(stridedBlockSize == 1) {
          // loop over all nodes
          for (decltype(vertex2AggId.dimension_0()) lnode = 0; lnode < vertex2AggId.dimension_0(); lnode++) {
            if(procWinner(lnode,0) == myPid) {
              auto myAgg = vertex2AggId(lnode,0);
              if(myAgg == aggId)
                aggSizes(myAgg+1)++;
            }
          }
        } else {
          // loop over all nodes
          for (decltype(vertex2AggId.dimension_0()) lnode = 0; lnode < vertex2AggId.dimension_0(); lnode++) {
            if(procWinner(lnode,0) == myPid) {
              auto myAgg = vertex2AggId(lnode,0);
              if(myAgg == aggId) {
                GO gnodeid = nodeGlobalElts(lnode);

                for (LO k = 0; k< stridedBlockSize; k++) {
                  GO gDofIndex = globalOffset + (gnodeid - indexBase) * fullBlockSize + stridingOffset + k + indexBase;
                  if(isNodeGlobalElement.value_at(isNodeGlobalElement.find(gDofIndex)) == true) {
                    aggSizes(myAgg+1)++;
                  }
                }
              }
            }
          }
        }

      }
    };


    template<class agg2RowMapType, class aggSizesType, class vertex2AggIdType, class procWinnerType, class nodeGlobalEltsType, class isNodeGlobalEltsType, class LOType, class GOType>
    class CreateAgg2RowMapLOFunctor {
    private:
      typedef LOType LO;
      typedef GOType GO;

      agg2RowMapType       agg2RowMap;      //< view containing row map entries associated with aggregate
      aggSizesType         aggSizes;        //< view containing size of aggregate
      vertex2AggIdType     vertex2AggId;    //< view containing vertex2AggId information
      procWinnerType       procWinner;      //< view containing processor ids
      nodeGlobalEltsType   nodeGlobalElts;  //< view containing global node ids of current proc
      isNodeGlobalEltsType isNodeGlobalElement; //< unordered map storing whether (global) node id is owned by current proc

      int myPid;
      LO fullBlockSize, blockID, stridingOffset, stridedBlockSize;
      GO indexBase, globalOffset;
    public:
      CreateAgg2RowMapLOFunctor(agg2RowMapType agg2RowMap_, aggSizesType aggSizes_, vertex2AggIdType vertex2AggId_, procWinnerType procWinner_, nodeGlobalEltsType nodeGlobalElts_, isNodeGlobalEltsType isNodeGlobalElement_, LO fullBlockSize_, LOType blockID_, LOType stridingOffset_, LOType stridedBlockSize_, GOType indexBase_, GOType globalOffset_, int myPID_ ) :
        agg2RowMap(agg2RowMap_),
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
      void operator()(const LO aggId) const {
        LO myAggDofSize  = 0;
        LO myAggDofStart = aggSizes(aggId);
        if(stridedBlockSize == 1) {
          // loop over all nodes
          for (decltype(vertex2AggId.dimension_0()) lnode = 0; lnode < vertex2AggId.dimension_0(); lnode++) {
            if(procWinner(lnode,0) == myPid) {
              auto myAgg = vertex2AggId(lnode,0);
              if(myAgg == aggId) {
                agg2RowMap(myAggDofStart + myAggDofSize) = lnode;
                myAggDofSize++;
              }
            }
          }
        } else {
          // loop over all nodes
          for (decltype(vertex2AggId.dimension_0()) lnode = 0; lnode < vertex2AggId.dimension_0(); lnode++) {
            if(procWinner(lnode,0) == myPid) {
              auto myAgg = vertex2AggId(lnode,0);
              if(myAgg == aggId) {
                GO gnodeid = nodeGlobalElts[lnode];

                for (LO k = 0; k< stridedBlockSize; k++) {
                  GO gDofIndex = globalOffset + (gnodeid - indexBase) * fullBlockSize + stridingOffset + k + indexBase;
                  if(isNodeGlobalElement.value_at(isNodeGlobalElement.find(gDofIndex)) == true) {
                    agg2RowMap(myAggDofStart + myAggDofSize) = lnode * stridedBlockSize + k;
                    myAggDofSize++;
                  }
                }
              }
            }
          }
        }
      }
    };


    typedef typename Kokkos::TeamPolicy<>::member_type team_member ;

    // local QR decomposition
    template<class LOType, class GOType, class SCType,class DeviceType, class NspType, class aggRowsType, class maxAggDofSizeType, class agg2RowMapLOType, class statusType>
    class TestFunctor {
    private:
      typedef LOType LO;
      typedef GOType GO;
      typedef SCType SC;

      //typedef Kokkos::DefaultExecutionSpace::scratch_memory_space shared_space;
      //typedef typename MatrixType::execution_space::scratch_memory_space shared_space;
      //typedef typename DeviceType::scratch_memory_space shared_space;
      typedef Kokkos::View<SC**,Kokkos::MemoryUnmanaged> shared_matrix;
      typedef Kokkos::View<SC*,Kokkos::MemoryUnmanaged> shared_vector;

    private:
      //MatrixType kokkosMatrix; //< local matrix part
      //NnzType nnz;             //< View containing number of nonzeros for current row
      //blkSizeType blkSize;     //< block size (or partial block size in strided maps)

      NspType fineNS;
      NspType coarseNS;
      aggRowsType aggRows;
      maxAggDofSizeType maxAggDofSize; //< maximum number of dofs in aggregate (max size of aggregate * numDofsPerNode)
      agg2RowMapLOType agg2RowMapLO;
      statusType statusAtomic;
    public:
      TestFunctor(NspType fineNS_, NspType coarseNS_, aggRowsType aggRows_, maxAggDofSizeType maxAggDofSize_, agg2RowMapLOType agg2RowMapLO_, statusType statusAtomic_) :
        fineNS(fineNS_),
        coarseNS(coarseNS_),
        aggRows(aggRows_),
        maxAggDofSize(maxAggDofSize_),
        agg2RowMapLO(agg2RowMapLO_),
        statusAtomic(statusAtomic_)
        { }

      KOKKOS_INLINE_FUNCTION
      void operator() ( const team_member & thread) const {
        auto agg = thread.league_rank();

        printf("Team no: %i, thread no: %i\n", thread.league_rank(), thread.team_rank());

        LO aggSize = aggRows(agg+1) - aggRows(agg);

        //SC one  = Teuchos::ScalarTraits<SC>::one();
        //SC zero = Teuchos::ScalarTraits<SC>::zero();

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
            if (localQR(k,j) != 0.0)  // TODO scalar traits
              bIsZeroNSColumn = false;

          if (bIsZeroNSColumn) {
            statusAtomic(1) = true;
            return;
          }
        }


#if 0
        shared_matrix mat(thread.team_shmem(),5,3);
        mat(0,0) = 12.0;
        mat(0,1) = -51.0;
        mat(0,2) = 4.0;
        mat(1,0) = 6.0;
        mat(1,1) = 167.0;
        mat(1,2) = -68.0;
        mat(2,0) = -4.0;
        mat(2,1) = 24.0;
        mat(2,2) = -41.0;
        mat(3,0) = -1.0;
        mat(3,1) = 1.0;
        mat(3,2) = 0.0;
        mat(4,0) = 2.0;
        mat(4,1) = 0.0;
        mat(4,2) = 3.0;
#endif

        ///////////////////////////////////
        // mat is the input matrix

        shared_matrix r(thread.team_shmem(),aggSize,fineNS.dimension_1()); // memory containing the r part in the end
        shared_matrix z(thread.team_shmem(),aggSize,fineNS.dimension_1()); // helper matrix (containing parts of localQR)
        shared_vector e(thread.team_shmem(),aggSize); //unit vector

        shared_matrix qk(thread.team_shmem(),aggSize,aggSize);  // memory cotaining one householder reflection part
        shared_matrix q (thread.team_shmem(),aggSize,aggSize);  // memory containing q part in the end
        shared_matrix qt (thread.team_shmem(),aggSize,aggSize); // temporary

        matrix_copy(localQR,z);

        for(decltype(localQR.dimension_0()) k = 0; k < localQR.dimension_0() && k < localQR.dimension_1()/*-1*/; k++) {
          // extract minor parts from mat
          matrix_clear(r);  // zero out temporary matrix (there is some potential to speed this up by avoiding this)
          matrix_minor(z,r,k);

          // extract k-th column from current minor
          auto x  = subview(r, Kokkos::ALL (), k);
          SC   xn = vnorm(x); // calculate 2-norm of current column vector
          if(localQR(k,k) > 0) xn = -xn;

          // build k-th unit vector
          for(decltype(e.dimension_0()) i = 0; i < e.dimension_0(); i++)
            e(i) = (i==k) ?  1.0 :  0.0;

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

        Xpetra::global_size_t offset = agg * fineNS.dimension_1(); // calculate row offset for coarse nullspace

        // upper triangular part of R build coarse NS
        for(size_t j = 0; j < fineNS.dimension_1(); j++)
          for(size_t k = 0; k <= j; k++)
            coarseNS(offset+k,j) = r(k,j);

        /*printf("R\n");
        for(int i=0; i<aggSize; i++) {
          for(int j=0; j<fineNS.dimension_1(); j++) {
            printf(" %.3g ",r(i,j));
          }
          printf("\n");
        }
        printf("Q\n");
        matrix_transpose(q);
        for(int i=0; i<aggSize; i++) {
          for(int j=0; j<aggSize; j++) {
            printf(" %.3g ",q(i,j));
          }
          printf("\n");
        }*/


        /*auto colview0 = subview(matminor, Kokkos::ALL (), 0);
        auto colview1 = subview(matminor, Kokkos::ALL (), 1);
        auto colview2 = subview(matminor, Kokkos::ALL (), 2);

        for(int i=0; i<3; i++) {
          printf(" %.3g . ", colview0(i));
        }
        printf("\n");

        printf("SQRT(4) = %.4g\n", vnorm(colview0));*/
      }

      KOKKOS_INLINE_FUNCTION
      void matrix_clear ( shared_matrix & m) const {
        //SC zero = Teuchos::ScalarTraits<SC>::zero();
        for(decltype(m.dimension_0()) i = 0; i < m.dimension_0(); i++) {
          for(decltype(m.dimension_1()) j = 0; j < m.dimension_1(); j++) {
            m(i,j) = 0.0;
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
        matrix_clear(m1m2);
        for(decltype(m1.dimension_0()) i = 0; i < m1.dimension_0(); i++) {
          for(decltype(m1.dimension_1()) j = 0; j < m1.dimension_1(); j++) {
            for(decltype(m1.dimension_1()) k = 0; k < m1.dimension_1(); k++) {
              m1m2(i,j) += m1(i,k) * m2(k,j);
            }
          }
        }
      }

      KOKKOS_FUNCTION
      void matrix_minor ( const shared_matrix & mat, shared_matrix & matminor, LO d) const {
        //SC one = Teuchos::ScalarTraits<SC>::one();
        for (LO i = 0; i < d; i++) {
          matminor(i,i) = 1.0; //one;
        }
        for (LO i = d; i < mat.dimension_0(); i++) {
          for (LO j=d; j < mat.dimension_1(); j++) {
            matminor(i,j) = mat(i,j);
          }
        }
      }

      /// \brief Build vmul = I - v*v^T
      /// \param v[in] input vector
      ///
      KOKKOS_FUNCTION
      void vmul ( const shared_vector & v, shared_matrix & vmuldata) const {
        //SC one = Teuchos::ScalarTraits<SC>::one();
        for(decltype(v.dimension_0()) i = 0; i < v.dimension_0(); i++) {
          for(decltype(v.dimension_0()) j = 0; j < v.dimension_0(); j++) {
            vmuldata(i,j) = -2 * v(i) * v(j);
          }
        }
        for(decltype(v.dimension_0()) i = 0; i < v.dimension_0(); i++) {
          vmuldata(i,i) += 1;
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
        printf("team size = %i\n", team_size);
        return 3 * Kokkos::View<double**,Kokkos::MemoryUnmanaged>::shmem_size(maxAggDofSize,fineNS.dimension_1()) + // mat + matminor + z
               3 * Kokkos::View<double**,Kokkos::MemoryUnmanaged>::shmem_size(maxAggDofSize,maxAggDofSize) +  // qk and q and qt
               Kokkos::View<double*,Kokkos::MemoryUnmanaged>::shmem_size(maxAggDofSize); // e
      }



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
    const size_t NSDim    = fineNullspace->getNumVectors();

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
    //TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() != 1, Exceptions::RuntimeError, "For now, only block size 1");

    LO fullBlockSize, blockID, stridingOffset, stridedBlockSize;
    GO indexBase;
    amalgInfo->GetStridingInformation(fullBlockSize, blockID, stridingOffset, stridedBlockSize, indexBase);
    std::cout << "fullBlockSize " << fullBlockSize << " blockID " << blockID << " stridingOffset " << stridingOffset << " stridedBlockSize " << stridedBlockSize << std::endl;

    GO globalOffset = amalgInfo->GlobalOffset();


    // do unamalgamation

    // (overlapping) local node Ids belonging to aggregates on current processor
    Teuchos::ArrayView<const GO> nodeGlobalEltsView = aggregates->GetMap()->getNodeElementList();
    Kokkos::View<GO*, DeviceType> nodeGlobalElts("nodeGlobalElts", nodeGlobalEltsView.size());
    for(size_t i = 0; i < nodeGlobalElts.size(); i++)
      nodeGlobalElts(i) = nodeGlobalEltsView[i];

    // TODO we have to create a view out of it
    auto procWinner   = aggregates->GetProcWinner()->getHostLocalView();
    auto vertex2AggId = aggregates->GetVertex2AggId()->getHostLocalView();
    const GO numAggregates = aggregates->GetNumAggregates();

    std::cout << "numAggregates = " << numAggregates << " vertex2AggId.dim0=" << vertex2AggId.dimension_0()<< " vertex2AggId.dim1=" << vertex2AggId.dimension_1() << std::endl;

    for(decltype(vertex2AggId.dimension_0()) i = 0; i < vertex2AggId.dimension_0(); i++) {
      std::cout << nodeGlobalElts[i] << "(" << vertex2AggId(i,0) << ") ";
    }
    std::cout << std::endl;

    typedef Kokkos::UnorderedMap<LO, bool, DeviceType> map_type;
    //map_type isNodeGlobalElement(numAggregates);
    map_type isNodeGlobalElement(colMap->getNodeNumElements());

    int myPid = aggregates->GetMap()->getComm()->getRank();

    for(LO i = 0; i < colMap->getNodeNumElements(); i++)
      std::cout << colMap->getGlobalElement(i) << " ";
    std::cout << std::endl;

    // create a unordered map GID -> isGlobalElement in colMap of A (available from above)
    // This has to be done on the host
    // I hope that the UnorderedMap can be used within kernels
    for (decltype(vertex2AggId.dimension_0()) lnode = 0; lnode < vertex2AggId.dimension_0(); lnode++) {
      auto myAgg = vertex2AggId(lnode,0);
      if(procWinner(lnode,0) == myPid) {
        GO gnodeid = nodeGlobalElts[lnode];
        for (LO k = 0; k< stridedBlockSize; k++) {
          GO gDofIndex = globalOffset + (gnodeid - indexBase) * fullBlockSize + stridingOffset + k + indexBase;
          bool bIsInColumnMap = colMap->isNodeGlobalElement(gDofIndex);
          isNodeGlobalElement.insert(gDofIndex, bIsInColumnMap);
          //std::cout << gnodeid << "->" << gDofIndex << ": " << isNodeGlobalElement.value_at(isNodeGlobalElement.find(gDofIndex)) << " should be " << bIsInColumnMap << std::endl;
        }
      }
    }

    // write parallel kernel to detect size of aggregates
    // TODO: have an outer for loop over all aggregates
    //       each thread loops through all nodes and checks whether node
    //       is owned by current proc and belongs to the aggregate id associated
    //       with the process. This way we avoid the race condition?
    // Attention: aggregate ids are "local"?
    Kokkos::View<LO*, DeviceType> sizes("agg_dof_sizes", numAggregates + 1);
    sizes(0) = 0;

#if 1
    // It's not clear whether the parallel routine is really faster than the serial one
    // Each process has to loop over all nodes. The parallel for loop is over an additional
    // outer for loop over the aggregates which actually is not needed in the serial code.
    CollectAggregateSizeFunctor <decltype(sizes), decltype(vertex2AggId), decltype(procWinner), decltype(nodeGlobalElts), decltype(isNodeGlobalElement), LO, GO> collectAggregateSizes(sizes, vertex2AggId, procWinner, nodeGlobalElts, isNodeGlobalElement, fullBlockSize, blockID, stridingOffset, stridedBlockSize, indexBase, globalOffset, myPid );
    Kokkos::parallel_for("MueLu:TentativePF:Build:getaggsizes", numAggregates, collectAggregateSizes);

    // the alternative would be to avoid race conditions and allow parallel write acces to sizes
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

    // get maximum size of aggregate
    LO maxAggSize = 0;
    for(LO i = 0; i < sizes.dimension_0(); i++) {
      if(sizes(i) > maxAggSize) maxAggSize = sizes(i);
      std::cout << "aggregate " << i << ": " << sizes(i) << std::endl;
    }
    std::cout << "maxAggSize = " << maxAggSize << std::endl;

    // parallel_scan (exclusive)
    ScanFunctor<LO,decltype(sizes)> scanFunctorAggSizes(sizes);
    Kokkos::parallel_scan("MueLu:TentativePF:Build:aggregate_sizes:stage1_scan", numAggregates+1, scanFunctorAggSizes);

    // create "map" aggregate id 2 row map
    // same problem as above
    // TODO add outer parallel loop over all aggregates
    Kokkos::View<LO*, DeviceType> agg2RowMapLO("agg2row_map_LO", numRows); // initialized to 0


#if 1
    CreateAgg2RowMapLOFunctor<decltype(agg2RowMapLO), decltype(sizes), decltype(vertex2AggId), decltype(procWinner), decltype(nodeGlobalElts), decltype(isNodeGlobalElement), LO, GO> createAgg2RowMap (agg2RowMapLO, sizes, vertex2AggId, procWinner, nodeGlobalElts, isNodeGlobalElement, fullBlockSize, blockID, stridingOffset, stridedBlockSize, indexBase, globalOffset, myPid );
    Kokkos::parallel_for("MueLu:TentativePF:Build:createAgg2RowMap", numAggregates, createAgg2RowMap);
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

    /*for(LO i = 0; i < numRows; i++) {
      std::cout << i << " dof " << agg2RowMapLO(i) << std::endl;
    }*/

#else

    // temporarely keep the unamalgamation code until we have a kokkos version of it
    ArrayRCP<LO> aggStart;
    ArrayRCP<LO> array_aggToRowMapLO;
    ArrayRCP<GO> array_aggToRowMapGO;
    if (goodMap) {
      amalgInfo->UnamalgamateAggregatesLO(*aggregates, aggStart, array_aggToRowMapLO);
      GetOStream(Runtime1) << "Column map is consistent with the row map, good." << std::endl;

    } else {
      amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, array_aggToRowMapGO);
      GetOStream(Warnings0) << "Column map is not consistent with the row map\n"
                            << "using GO->LO conversion with performance penalty" << std::endl;
    }

    // TAW: at the momemt only support for matching row and col maps
    TEUCHOS_TEST_FOR_EXCEPTION(!goodMap,                    Exceptions::RuntimeError, "For now, need matching row and col maps");

    Kokkos::View<LO*, DeviceType> agg2RowMapLO("agg2row_map_LO", array_aggToRowMapLO.size());
    for(size_t i = 0; i < array_aggToRowMapLO.size(); i++) {
      agg2RowMapLO(i) = array_aggToRowMapLO[i];
    }
#endif
    // TODO
    // Use TeamPolicy with scratch_memory_space for local QR
    // Something along the lines:
    //
    //   typedef TeamPolicy<ExecutionSpace>::member_type team_t;
    //   struct functor() {
    //     inline unsigned team_shmem_size( int team_size ) const {
    //       return view_type::shmem_size( team_size, 10, 3 );
    //     }
    //     KOKKOS_INLINE_FUNCTION
    //     void operator() (const team_t& team) const {
    //       view_type matrices(team.team_shmem(), team.team_size(), 10, 3);
    //       auto matrix = subview(matrices, team.team_rank(), Kokkos::ALL(), Kokkos::ALL());
    //     }
    //   }

    coarseNullspace = MultiVectorFactory::Build(coarseMap, NSDim);

    // Pull out the nullspace vectors so that we can have random access
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
    typename rows_type::non_const_type rowsAux("Ptent_aux_rows", numRows+1),    rows("Ptent_rows", numRows+1);
    typename cols_type::non_const_type colsAux("Ptent_aux_cols", nnzEstimate);
    typename vals_type::non_const_type valsAux("Ptent_aux_vals", nnzEstimate);

    Kokkos::parallel_for("MueLu:TentativePF:BuildPuncoupled:for1", numRows+1, KOKKOS_LAMBDA(const LO row) {
      rowsAux(row) = row*NSDim;
    });
    Kokkos::parallel_for("MueLu:TentativePF:BuildUncoupled:for2", nnzEstimate, KOKKOS_LAMBDA(const LO j) {
      colsAux(j) = INVALID;
      valsAux(j) = zero;
    });

    typedef Kokkos::View<int[10], DeviceType> status_type;
    status_type status("status");

    // Stage 1: construct auxilary arrays.
    // The constructed arrays may have gaps in them (vals(j) == INVALID)
    // Run one thread per aggregate.
    typename AppendTrait<decltype(fineNS), Kokkos::RandomAccess>::type fineNSRandom = fineNS;
    typename AppendTrait<status_type,      Kokkos::Atomic>      ::type statusAtomic = status;
    if (NSDim == 1) {
      // 1D is special, as it is the easiest. We don't even need to the QR,
      // just normalize an array. Plus, no worries abot small aggregates.
      Kokkos::parallel_reduce("MueLu:TentativePF:BuildUncoupled:main_loop", numAggs, KOKKOS_LAMBDA(const GO agg, size_t& rowNnz) {
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

      typename status_type::HostMirror statusHost = Kokkos::create_mirror_view(status);
      for (int i = 0; i < statusHost.size(); i++)
        if (statusHost(i)) {
          std::ostringstream oss;
          oss << "MueLu::TentativePFactory::MakeTentative: ";
          switch(i) {
            case 0: oss << "!goodMap is not implemented";               break;
            case 1: oss << "fine level NS part has a zero column";      break;
          }
          throw Exceptions::RuntimeError(oss.str());
        }

    } else {

      TEUCHOS_TEST_FOR_EXCEPTION(goodMap == false, Exceptions::RuntimeError, "Only works for non-overlapping aggregates (goodMap == true)");

      // Set up team policy with numAggregates teams and one thread per team.
      // Each team handles a slice of the data associated with one aggregate
      // and performs a local QR decomposition
      const Kokkos::TeamPolicy<> policy( numAggregates /*1*/ /*10*/ , 1); // 10 teams a 1 thread

      Kokkos::parallel_for( policy, TestFunctor<LocalOrdinal, GlobalOrdinal, Scalar, DeviceType, decltype(fineNSRandom), decltype(sizes /*aggregate sizes in dofs*/), decltype(maxAggSize), decltype(agg2RowMapLO), decltype(statusAtomic)>(fineNSRandom,coarseNS,sizes,maxAggSize,agg2RowMapLO,statusAtomic));


      std::cout << *coarseNullspace << std::endl;

      Teuchos::ArrayRCP< const Scalar > data0 = coarseNullspace->getData(0);
      Teuchos::ArrayRCP< const Scalar > data1 = coarseNullspace->getData(1);
      for (size_t i = 0; i < coarseNullspace->getLocalLength(); i++) {
        std::cout << i << "\t" << coarseNS(i,0) << "\t" << coarseNS(i,1) << std::endl;
        std::cout << i << "\t" << data0[i] << "\t" << data1[i] << std::endl;
      }


      throw Exceptions::RuntimeError("Ignore NSDim > 1 for now");


#if 0
      Kokkos::parallel_reduce("MueLu:TentativePF:BuildUncoupled:main_loop", numAggs, KOKKOS_LAMBDA(const GO agg, size_t& nnz) {
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
              localQR(k,j) = fineNSRandom(agg2RowMapLO(aggRows(agg)+k), j);
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
          LO localRow = (goodMap ? agg2RowMapLO(aggRows(agg)+j) : -1);
#else
          LO localRow = (goodMap ? agg2RowMapLO[aggRows(agg)+j] : rowMap->getLocalElement(aggToRowMapGO[aggStart[agg]+j]));
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

    // Stage 2: compress the arrays
    ScanFunctor<LO,decltype(rows)> scanFunctor(rows);
    Kokkos::parallel_scan("MueLu:TentativePF:Build:compress_rows", numRows+1, scanFunctor);

    // The real cols and vals are constructed using calculated (not estimated) nnz
    typename cols_type::non_const_type cols("Ptent_cols", nnz);
    typename vals_type::non_const_type vals("Ptent_vals", nnz);
    Kokkos::parallel_for("MueLu:TentativePF:Build:compress_cols_vals", numRows, KOKKOS_LAMBDA(const LO i) {
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

    // Stage 3: construct Xpetra::Matrix
    // FIXME: For now, we simply copy-paste arrays. The proper way to do that
    // would be to construct a Kokkos CrsMatrix, and then construct
    // Xpetra::Matrix out of that.
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
