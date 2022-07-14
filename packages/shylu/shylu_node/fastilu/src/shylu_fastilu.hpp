//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov) 
// 
// ************************************************************************
//@HEADER

// The struct that iterates over all non-zeros for the FastILU
//  Contact for bugs and complaints - Siva Rajamanickam (srajama@sandia.gov)
//
#ifndef __FAST_ILU_HPP__
#define __FAST_ILU_HPP__

#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <random>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string>

#include <assert.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#include <KokkosKernels_Sorting.hpp>
#include <KokkosKernels_SparseUtils.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosSparse_sptrsv.hpp>
#include <KokkosSparse_trsv.hpp>
#include <shylu_fastutil.hpp>

//whether to print extra debug output at runtime to stdout
//comment out next line to disable
//#define FASTILU_DEBUG_OUTPUT

//forward declaration
template<class Ordinal, class Scalar, class ExecSpace>
class FastILUFunctor;

template<class Ordinal, class Scalar, class ExecSpace>
class FastICFunctor;

template<class Ordinal, class Scalar, class ExecSpace>
class JacobiIterFunctor;

template<class Ordinal, class Scalar, class ExecSpace>
class JacobiIterFunctorT;

template<class Ordinal, class Scalar, class ExecSpace>
class BlockJacobiIterFunctorU;

template<class Ordinal, class Scalar, class ExecSpace>
class BlockJacobiIterFunctorL;

template<class Ordinal, class Scalar, class ExecSpace>
class ParCopyFunctor;

template<class Ordinal, class Scalar, class Real, class ExecSpace>
class ParScalFunctor;

struct NonTranPermScalTag {};
struct    TranPermScalTag {};
template<class Ordinal, class Scalar, class Real, class ExecSpace>
class PermScalFunctor;

template<class Ordinal, class Scalar, class ExecSpace>
class ParInitZeroFunctor;

template<class Ordinal, class Scalar, class ExecSpace>
class MemoryPrimeFunctorN;

template<class Ordinal, class Scalar, class ExecSpace>
class MemoryPrimeFunctorNnzCoo;

template<class Ordinal, class Scalar, class ExecSpace>
class MemoryPrimeFunctorNnzCsr;

template<class Ordinal, class Scalar, class ExecSpace>
class FastILUPrec
{
    public:
        typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Real;
        typedef Kokkos::View<Ordinal *, ExecSpace> OrdinalArray;
        typedef Kokkos::View<Scalar *, ExecSpace> ScalarArray;
        typedef Kokkos::View<Real *, ExecSpace> RealArray;
        typedef Kokkos::View<Ordinal *, typename ExecSpace::array_layout,
                             Kokkos::Serial, Kokkos::MemoryUnmanaged> UMOrdinalArray;
        typedef Kokkos::View<Scalar *, typename ExecSpace::array_layout,
                             Kokkos::Serial, Kokkos::MemoryUnmanaged> UMScalarArray;
        typedef FastILUPrec<Ordinal, Scalar, ExecSpace> FastPrec;

        using HostSpace = Kokkos::HostSpace;
        using MirrorSpace = typename OrdinalArray::host_mirror_space;

        typedef Kokkos::View<Ordinal *, HostSpace> OrdinalArrayHost;
        typedef Kokkos::View<Scalar  *, HostSpace>  ScalarArrayHost;
        typedef typename OrdinalArray::host_mirror_type OrdinalArrayMirror;
        typedef typename ScalarArray::host_mirror_type  ScalarArrayMirror;

        using STS = Kokkos::ArithTraits<Scalar>;
        using RTS = Kokkos::ArithTraits<Real>;

    private:
        double computeTime;
        double applyTime;
        double initTime;

        Ordinal nRows;
        Ordinal guessFlag;
        Ordinal nFact;
        Ordinal nTrisol;
        Ordinal level;
        Ordinal blkSzILU;
        Ordinal blkSz;
        Scalar omega; //Underrelaxation parameter
        Scalar shift; //Manteuffel Shift

        // Metis
        bool useMetis;
        OrdinalArray  permMetis;
        OrdinalArray ipermMetis;
        OrdinalArrayMirror  permMetisHost;
        OrdinalArrayMirror ipermMetisHost;

        //Lower triangular factor (CSR)
        bool sptrsv_KKSpMV; // use Kokkos-Kernels SpMV for Fast SpTRSV
        ScalarArray lVal;
        OrdinalArray lColIdx;
        OrdinalArray lRowMap;
        // mirrors
        ScalarArrayMirror lVal_;
        OrdinalArrayMirror lColIdx_;
        OrdinalArrayMirror lRowMap_;

        //Lower triangular factor, without (unit) diagonals,
        // for TRSV (not SpTRSV)
        ScalarArrayHost lVal_trsv_;
        OrdinalArrayHost lColIdx_trsv_;
        OrdinalArrayHost lRowMap_trsv_;

        //Upper triangular factor (CSC)
        ScalarArray uVal;
        OrdinalArray uColIdx;
        OrdinalArray uRowMap;
        // mirrors
        ScalarArrayMirror uVal_;
        OrdinalArrayMirror uColIdx_;
        OrdinalArrayMirror uRowMap_;

        //Upper triangular factor (CSR)
        ScalarArray utVal;
        OrdinalArray utColIdx;
        OrdinalArray utRowMap;
        // mirrors
        ScalarArrayMirror utVal_;
        OrdinalArrayMirror utColIdx_;
        OrdinalArrayMirror utRowMap_;

        //Upper triangular factor (CSR), with diagonal extracted out
        // for TRSV (not SpTRSV)
        bool doUnitDiag_TRSV; // perform TRSV with unit diagonals
        ScalarArrayHost   dVal_trsv_;
        ScalarArrayHost  utVal_trsv_;
        OrdinalArrayHost utColIdx_trsv_;
        OrdinalArrayHost utRowMap_trsv_;

        //Pointer to the copy of input A.
        // device
        ScalarArray        aValIn;
        OrdinalArray       aRowMapIn;
        OrdinalArray       aColIdxIn;
        // host
        ScalarArrayMirror  aValHost;
        OrdinalArrayMirror aRowMapHost;
        OrdinalArrayMirror aColIdxHost;

        //A matrix in COO format
        ScalarArray aVal;
        OrdinalArray aRowMap;
        OrdinalArray aRowIdx;
        OrdinalArray aColIdx;
        // mirrors
        ScalarArrayMirror aVal_;
        OrdinalArrayMirror aRowMap_;
        OrdinalArrayMirror aRowIdx_;
        OrdinalArrayMirror aColIdx_;
        OrdinalArrayHost   aLvlIdx_;

        //Diagonal scaling factors
        RealArray diagFact;
        ScalarArray diagElems;

        //Temporary vectors for triangular solves
        ScalarArray xOld;
        ScalarArray xTemp;
        ScalarArray onesVector;

        //This will have the continuation initial
        //guess if guessFlag=1
        Teuchos::RCP<FastPrec> initGuessPrec;

        // forward/backwar substitution for standard SpTrsv
        using MemSpace = typename ExecSpace::memory_space;
        using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle <Ordinal, Ordinal, Scalar, ExecSpace, MemSpace, MemSpace >;
        FastILU::SpTRSV sptrsv_algo;

        KernelHandle khL;
        KernelHandle khU;

        //Internal functions
        //Serial Transpose for now.
        //TODO:convert to parallel.
        void transposeU()
        {
            //Count the elements in each row of Ut
            auto temp = OrdinalArrayHost("temp", nRows + 1);
            auto rowPtrs = OrdinalArrayHost("rowPtrs", nRows);
            for (Ordinal i = 0; i <= nRows; i++) 
            {
                temp[i] = 0;
            }
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = uRowMap_[i]; k < uRowMap_[i+1]; k++) 
                {
                    temp[uColIdx_[k]+1]++;

                }
            }
            //Perform an add scan to get the row map for 
            //the transpose
            for (Ordinal i = 0; i <= nRows; i++) 
            {
                utRowMap_[i] = temp[i];
            }
            for (Ordinal i = 1; i <= nRows; i++) 
            {
                utRowMap_[i] += utRowMap_[i-1];
            }
            //Set the row pointers to their initial places;
            for (Ordinal i = 0; i < nRows; i++) 
            {
                rowPtrs[i] = utRowMap_[i];
            }
            //Copy the data
            Kokkos::deep_copy(uVal_, uVal);
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = uRowMap_[i]; k < uRowMap_[i+1]; k++)
                {
                    Ordinal row = uColIdx_[k];
                    Scalar value = uVal_[k];
                    utVal_[rowPtrs[row]] = value;
                    utColIdx_[rowPtrs[row]] = i;
                    rowPtrs[row]++;
                    assert(rowPtrs[row] <= utRowMap_[row + 1]);
                }
            }
            Kokkos::deep_copy(utRowMap, utRowMap_);
            Kokkos::deep_copy(utColIdx, utColIdx_);
            Kokkos::deep_copy(utVal, utVal_);
        }
        //
        void findFills(int levfill, OrdinalArrayMirror aRowMap, OrdinalArrayMirror aColIdx,
                       int *nzl, std::vector<int> &lRowMap, std::vector<int> &lColIdx, std::vector<int> &lLevel,
                       int *nzu, std::vector<int> &uRowMap, std::vector<int> &uColIdx, std::vector<int> &uLevel) {
            using std::vector;
            using std::cout;
            using std::sort;

            Ordinal n = nRows;

            Ordinal row = 0, i = 0;
            vector<int> lnklst(n);
            vector<int> curlev(n);
            vector<int> iwork(n);

            int knzl = 0;
            int knzu = 0;

            lRowMap[0] = 0;
            uRowMap[0] = 0;

            for (i=0; i<n; i++)
            {
                int first, next, j;
                row = (useMetis ? permMetisHost(i) : i);

                /* copy column indices of row into workspace and sort them */
                int len = aRowMap[row+1] - aRowMap[row];
                next = 0;
                for (j=aRowMap[row]; j<aRowMap[row+1]; j++) {
                    iwork[next++] = (useMetis ? ipermMetisHost(aColIdx[j]) : aColIdx[j]);
                }
                // sort column indices in non-descending (ascending) order
                sort(iwork.begin(), iwork.begin() + len);

                /* construct implied linked list for row */
                first = iwork[0];
                curlev[first] = 0;

                for (j=0; j<=len-2; j++)
                {
                    lnklst[iwork[j]] = iwork[j+1];
                    curlev[iwork[j]] = 0;
                }

                lnklst[iwork[len-1]] = n;
                curlev[iwork[len-1]] = 0;

                /* merge with rows in U */
                next = first;
                while (next < i)
                {
                    int oldlst = next;
                    int nxtlst = lnklst[next];
                    int row = next;
                    int ii;

                    /* scan row */
                    for (ii=uRowMap[row]+1; ii<uRowMap[row+1]; /*nop*/)
                    {
                        if (uColIdx[ii] < nxtlst)
                        {
                            /* new fill-in */
                            int newlev = curlev[row] + uLevel[ii] + 1;
                            if (newlev <= levfill)
                            {
                                lnklst[oldlst]  = uColIdx[ii];
                                lnklst[uColIdx[ii]] = nxtlst;
                                oldlst = uColIdx[ii];
                                curlev[uColIdx[ii]] = newlev;
                            }
                            ii++;
                        }
                        else if (uColIdx[ii] == nxtlst)
                        {
                            int newlev;
                            oldlst = nxtlst;
                            nxtlst = lnklst[oldlst];
                            newlev = curlev[row] + uLevel[ii] + 1;
                            //curlev[uColIdx[ii]] = MIN(curlev[uColIdx[ii]], newlev);
                            if (curlev[uColIdx[ii]] > newlev)
                            {
                                curlev[uColIdx[ii]] = newlev;
                            }
                            ii++;
                        }
                        else /* (jau[ii] > nxtlst) */
                        {
                            oldlst = nxtlst;
                            nxtlst = lnklst[oldlst];
                        }
                    }
                    next = lnklst[next];
                }

                /* gather the pattern into L and U */
                /* L (no diagonal) */
                next = first;
                while (next < i)
                {
                    assert(knzl < *nzl);
                    lLevel[knzl] = curlev[next];
                    lColIdx[knzl++] = next;
                    if (knzl >= *nzl)
                    {
                        *nzl = *nzl + n;
                        lColIdx.resize(*nzl);
                        lLevel.resize(*nzl);
                    }
                    next = lnklst[next];
                }
                lRowMap[i+1] = knzl;
                assert(next == i);
                /* U (with diagonal) */
                while (next < n)
                {
                    assert(knzu < *nzu);
                    uLevel[knzu] = curlev[next];
                    uColIdx[knzu++] = next;
                    if (knzu >= *nzu)
                    {
                        *nzu = *nzu + n;
                        uColIdx.resize(*nzu);
                        uLevel.resize(*nzu);
                    }
                    next = lnklst[next];
                }
                uRowMap[i+1] = knzu;
            }

            *nzl = knzl;
            *nzu = knzu;

        }
        //Symbolic ILU code
        //initializes the matrices L and U and readies them
        //according to the level of fill
        void symbolicILU()
        {
            #ifdef FASTILU_INIT_TIMER
            Kokkos::Timer timer;
            #endif
            using std::vector;
            using std::cout;
            using std::stable_sort;
            using std::sort;
            OrdinalArrayMirror ia = aRowMapHost;
            OrdinalArrayMirror ja = aColIdxHost;
            int *nzu;
            int *nzl;
            nzu = new int[1];
            nzl = new int[1];
            Ordinal i;

            //Compute sparsity structure of ILU
            *nzl = aRowMapHost[nRows];
            *nzu = aRowMapHost[nRows];

            *nzl *= (level + 2);
            *nzu *= (level + 2);
            vector<int> ial(nRows+1);
            vector<int> jal(*nzl);
            vector<int> levell(*nzl);
            vector<int> iau(nRows+1);
            vector<int> jau(*nzu);
            vector<int> levelu(*nzu);

            // TODO: if (initGuess & level > 0), call this with (aRowMap_, aColIdx_) and level = 1
            findFills(level, aRowMapHost, aColIdxHost, // input
                      nzl, ial, jal, levell, // output L in CSR
                      nzu, iau, jau, levelu  // output U in CSR
                     );
            int knzl = *nzl;
            int knzu = *nzu;
            #ifdef FASTILU_INIT_TIMER
            std::cout << " findFills time : " << timer.seconds() << std::endl;
            timer.reset();
            #endif

            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "knzl =" << knzl;
            std::cout << "knzu =" << knzu;

            std::cout << "ILU: nnz = "<< knzl + knzu << std::endl;
            std::cout << "Actual nnz for ILU: " << *nzl + *nzu << std::endl;
            #endif


            //Initialize the A matrix that is to be used in the computation
            aRowMap = OrdinalArray("aRowMap", nRows + 1);
            aColIdx = OrdinalArray("aColIdx", knzl + knzu);
            aRowIdx = OrdinalArray("aRowIds", knzl + knzu);
            aRowMap_ = Kokkos::create_mirror(aRowMap);
            aColIdx_ = Kokkos::create_mirror(aColIdx);
            aRowIdx_ = Kokkos::create_mirror(aRowIdx);

            aLvlIdx_ = OrdinalArrayHost("aLvlIdx", knzl + knzu);

            aVal = ScalarArray("aVal", aColIdx.extent(0));
            aVal_ = Kokkos::create_mirror(aVal);

            Ordinal aRowPtr = 0;
            aRowMap_[0] = aRowPtr;
            for (i = 0; i < nRows; i++) 
            {
                #ifdef FASTILU_DEBUG_OUTPUT
                std::cout << "***row:" << i << std::endl;
                #endif
                for(Ordinal k = ial[i]; k < ial[i+1]; k++)
                {
                    #ifdef FASTILU_DEBUG_OUTPUT
                    std::cout << "jal[k]=" << jal[k] << std::endl;
                    #endif
                    aColIdx_[aRowPtr] = jal[k];
                    aRowIdx_[aRowPtr] = i;
                    aLvlIdx_[aRowPtr] = levell[k];
                    aRowPtr++;
                }
                for(Ordinal k = iau[i]; k < iau[i+1]; k++)
                {
                    aColIdx_[aRowPtr] = jau[k];
                    aRowIdx_[aRowPtr] = i;
                    aLvlIdx_[aRowPtr] = levelu[k];
                    aRowPtr++;
                }
                aRowMap_[i+1] = aRowPtr;
            }
            #ifdef FASTILU_INIT_TIMER
            std::cout << " Copy time : " << timer.seconds() << std::endl;
            timer.reset();
            #endif
            // sort based on ColIdx, RowIdx stays the same (do we need this?)
            using host_space = typename HostSpace::execution_space;
            KokkosKernels::sort_crs_matrix<host_space, OrdinalArrayMirror, OrdinalArrayMirror, ScalarArrayMirror>
              (aRowMap_, aColIdx_, aVal_);
            host_space().fence();
            #ifdef FASTILU_INIT_TIMER
            std::cout << " Sort time : " << timer.seconds() << std::endl;
            timer.reset();
            #endif

            Kokkos::deep_copy(aRowMap, aRowMap_);
            Kokkos::deep_copy(aColIdx, aColIdx_);
            Kokkos::deep_copy(aRowIdx, aRowIdx_);
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished initializing A" << std::endl;
            #endif

            //Compute RowMap for L and U. 
            // > form RowMap for L
            lRowMap = OrdinalArray("lRowMap", nRows + 1);
            lRowMap_ = Kokkos::create_mirror(lRowMap);
            Ordinal nnzL = countL();
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished counting L" << std::endl;
            #endif
            
            // > form RowMap for U and Ut
            uRowMap  = OrdinalArray("uRowMap", nRows + 1);
            utRowMap = OrdinalArray("utRowMap", nRows + 1);
            utRowMap_ = Kokkos::create_mirror(utRowMap);
            uRowMap_  = Kokkos::create_mirror(uRowMap);
            Ordinal nnzU = countU();

            // > form RowMap for Ut
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished counting U" << std::endl;
            #endif

            //Allocate memory and initialize pattern for L, U (transpose).
            lColIdx = OrdinalArray("lColIdx", nnzL);
            uColIdx = OrdinalArray("uColIdx", nnzU);
            utColIdx = OrdinalArray("utColIdx", nnzU);

            lVal = ScalarArray("lVal", nnzL);
            uVal = ScalarArray("uVal", nnzU);
            utVal = ScalarArray("utVal", nnzU);

            //Create mirror
            lColIdx_ = Kokkos::create_mirror(lColIdx);
            uColIdx_ = Kokkos::create_mirror(uColIdx);
            utColIdx_ = Kokkos::create_mirror(utColIdx);

            lVal_    = Kokkos::create_mirror(lVal);
            uVal_    = Kokkos::create_mirror(uVal);
            utVal_    = Kokkos::create_mirror(utVal);
            #ifdef FASTILU_INIT_TIMER
            std::cout << " Mirror : " << timer.seconds() << std::endl;
            timer.reset();
            #endif
        }

        void symbolicILU(OrdinalArrayMirror pRowMap_, OrdinalArrayMirror pColIdx_, ScalarArrayMirror pVal_, OrdinalArrayHost pLvlIdx_)
        {
            Ordinal nnzA = 0;
            for (Ordinal k = 0; k < pRowMap_(nRows); k++)  {
                if(pLvlIdx_(k) <= level) {
                   nnzA++;
                }
            }
            //Initialize the A matrix that is to be used in the computation
            aRowMap = OrdinalArray("aRowMap", nRows + 1);
            aColIdx = OrdinalArray("aColIdx", nnzA);
            aRowIdx = OrdinalArray("aRowIds", nnzA);
            aRowMap_ = Kokkos::create_mirror(aRowMap);
            aColIdx_ = Kokkos::create_mirror(aColIdx);
            aRowIdx_ = Kokkos::create_mirror(aRowIdx);

            aVal = ScalarArray("aVal", nnzA);
            aVal_ = Kokkos::create_mirror(aVal);

            Ordinal aRowPtr = 0;
            aRowMap_[0] = aRowPtr;
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for(Ordinal k = pRowMap_(i); k < pRowMap_(i+1); k++)
                {
                    if (pLvlIdx_(k) <= level) {
                        aVal_[aRowPtr] = pVal_[k];
                        aColIdx_[aRowPtr] = pColIdx_[k];
                        aRowIdx_[aRowPtr] = i;
                        aRowPtr++;
                    }
                }
                aRowMap_[i+1] = aRowPtr;
            }
            // sort based on ColIdx, RowIdx stays the same (do we need this?)
            //using host_space = Kokkos::HostSpace::execution_space;
            //KokkosKernels::sort_crs_matrix<host_space, OrdinalArrayMirror, OrdinalArrayMirror, ScalarArrayMirror>
            //  (aRowMap_, aColIdx_, aVal_);
            //host_space().fence();

            Kokkos::deep_copy(aRowMap, aRowMap_);
            Kokkos::deep_copy(aColIdx, aColIdx_);
            Kokkos::deep_copy(aRowIdx, aRowIdx_);
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished initializing A" << std::endl;
            #endif

            //Now allocate memory for L and U. 
            //
            lRowMap = OrdinalArray("lRowMap", nRows + 1);
            lRowMap_ = Kokkos::create_mirror(lRowMap);
            countL();
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished counting L" << std::endl;
            #endif
            
            uRowMap = OrdinalArray("uRowMap", nRows + 1);
            utRowMap = OrdinalArray("utRowMap", nRows + 1);
            utRowMap_ = Kokkos::create_mirror(utRowMap);
            uRowMap_  = Kokkos::create_mirror(uRowMap);
            countU();
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished counting U" << std::endl;
            #endif

            //Allocate memory and initialize pattern for L, U (transpose).
            lColIdx = OrdinalArray("lColIdx", lRowMap_[nRows]);
            uColIdx = OrdinalArray("uColIdx", uRowMap_[nRows]);
            utColIdx = OrdinalArray("utColIdx", uRowMap_[nRows]);

            lVal = ScalarArray("lVal", lRowMap_[nRows]);
            uVal = ScalarArray("uVal", uRowMap_[nRows]);
            utVal = ScalarArray("utVal", uRowMap_[nRows]);

            //Create mirror
            lColIdx_ = Kokkos::create_mirror(lColIdx);
            uColIdx_ = Kokkos::create_mirror(uColIdx);
            utColIdx_ = Kokkos::create_mirror(utColIdx);

            lVal_    = Kokkos::create_mirror(lVal);
            uVal_    = Kokkos::create_mirror(uVal);
            utVal_   = Kokkos::create_mirror(utVal);
        }

        void numericILU()
        {
            const Scalar zero = STS::zero();
            #ifdef FASTILU_TIMER
            Kokkos::Timer Timer;
            #endif
            if (useMetis && level == 0) { // applied only at the first call (level 0)
              // apply column permutation before sorting it
              FastILUPrec_Functor perm_functor(aColIdxIn, ipermMetis);
              Kokkos::RangePolicy<ColPermTag, ExecSpace> perm_policy (0, aColIdxIn.size());
              Kokkos::parallel_for(
                "numericILU::colPerm", perm_policy, perm_functor);
              ExecSpace().fence();
            }
            //Sort each row of A by ColIdx
            KokkosKernels::sort_crs_matrix<ExecSpace, OrdinalArray, OrdinalArray, ScalarArray>(aRowMapIn, aColIdxIn, aValIn);
            ExecSpace().fence();

            //Copy the host matrix into the initialized a;
            //a contains the structure of ILU(k), values of original Ain is copied at level-0
            FastILUPrec_Functor functor(aValIn, aRowMapIn, aColIdxIn, aVal, diagFact, aRowMap, aColIdx, aRowIdx);
            if (useMetis) {
              FastILUPrec_Functor functor_perm(aValIn, aRowMapIn, aColIdxIn, permMetis, aVal, diagFact, aRowMap, aColIdx, aRowIdx);
              Kokkos::RangePolicy<CopySortedValsPermTag, ExecSpace> copy_perm_policy (0, nRows);
              Kokkos::parallel_for(
                "numericILU::copyVals", copy_perm_policy, functor_perm);
            } else {
              Kokkos::RangePolicy<CopySortedValsTag, ExecSpace> copy_policy (0, nRows);
              Kokkos::parallel_for(
                "numericILU::copyVals", copy_policy, functor);
            }
            ExecSpace().fence();
            #ifdef FASTILU_TIMER
            std::cout << "   + copy values  " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif

            // obtain diagonal scaling factor
            Kokkos::RangePolicy<GetDiagsTag, ExecSpace> get_policy (0, nRows);
            Kokkos::parallel_for(
              "numericILU::getDiags", get_policy, functor);
            ExecSpace().fence();
            // apply diagonal scaling
            Kokkos::RangePolicy<DiagScalTag, ExecSpace> scale_policy (0, nRows);
            Kokkos::parallel_for(
              "numericILU::diagScal", scale_policy, functor);
            Kokkos::deep_copy(aVal_, aVal);

            // applyShift
            if (shift != zero) {
                applyManteuffelShift();
            }
            #ifdef FASTILU_TIMER
            std::cout << "   + apply shift/scale  " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif

            Kokkos::deep_copy(aVal, aVal_);
            #ifdef FASTILU_TIMER
            std::cout << "   + deep copy  " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif

            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished diagonal scaling" << std::endl;
            #endif
            fillL();
            #ifdef FASTILU_TIMER
            std::cout << "   + fill L  " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished copying L" << std::endl;
            #endif
            fillU();
            #ifdef FASTILU_TIMER
            std::cout << "   + fill U  " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished copying U" << std::endl;
            std::cout << "nnz L = " << lRowMap_[nRows] << std::endl;
            std::cout << "nnz U = " << uRowMap_[nRows] << std::endl;
            #endif
        }

        //Initialize the rowMap (rowPtr) for L
        int countL()
        {
            lRowMap_[0] = 0;
            for (Ordinal i = 0; i < nRows; i++) 
            {
                Ordinal row_count = 0;
                for (Ordinal k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx_[k];

                    if (row >= col)
                    {
                       row_count++; 
                    }
                }
                lRowMap_[i+1] = lRowMap_[i] + row_count;
            }
            Kokkos::deep_copy(lRowMap, lRowMap_);
            return lRowMap_[nRows];
        }

        //Put the initial guess into L.
        void fillL()
        {
            // extract L out of A, where A contains the structure of ILU(k), and original nonzero values at level-0
            FastILUPrec_Functor functor(aVal, aRowMap, aColIdx, lVal, lRowMap, lColIdx, diagElems);
            Kokkos::RangePolicy<GetLowerTag, ExecSpace> getL_policy (0, nRows);
            Kokkos::parallel_for(
              "numericILU::getLower", getL_policy, functor);
            ExecSpace().fence();

            if ((level > 0) && (guessFlag !=0))
            {
                // overwrite initial values from warmup runs
                OrdinalArray lGRowMap;
                OrdinalArray lGColIdx;
                ScalarArray lGVal;
                ScalarArray gD;
                initGuessPrec->getL(lGRowMap, lGColIdx, lGVal);
                initGuessPrec->getD(gD);
                Kokkos::deep_copy(diagElems, gD);

                // copy LG into L
                FastILUPrec_Functor functorG(lGVal, lGRowMap, lGColIdx, lVal, lRowMap, lColIdx);
                Kokkos::RangePolicy<CopySortedValsTag, ExecSpace> copy_policy (0, nRows);
                Kokkos::parallel_for(
                  "numericILU::copyVals(G)", copy_policy, functorG);
                ExecSpace().fence();
            }
        }

        //Initialize rowMap of U
        int countU()
        {
            // extract U out of A, where A contains the structure of ILU(k), and original nonzero values at level-0
            for (Ordinal i = 0; i <= nRows; i++)
            {
                uRowMap_[i] = 0;
            }
            for(Ordinal i = 0; i < nRows; i++)
            {
                for(Ordinal k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx_[k];
                    if (row <= col)
                    {
                        uRowMap_[col+1]++;
                    }
                }
            }
            for (Ordinal i = 0; i < nRows; i++)
            {
                uRowMap_[i+1] += uRowMap_[i];
            }
            Kokkos::deep_copy(uRowMap, uRowMap_);
            return uRowMap_[nRows];
        }

        //Put initial guess into U
        void fillU()
        {
#if 1
            #ifdef FASTILU_TIMER
            Kokkos::Timer Timer;
            #endif
            // extract U
            Kokkos::deep_copy(utRowMap, uRowMap); // using utRowMap (will get incremented)
            FastILUPrec_Functor functor(aVal, aRowMap, aColIdx, uVal, utRowMap, uColIdx);
            Kokkos::RangePolicy<GetUpperTag, ExecSpace> getU_policy (0, nRows);
            Kokkos::parallel_for(
              "numericILU::getUpper", getU_policy, functor);
            #ifdef FASTILU_TIMER
            std::cout << "   + transpose_matrix  " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif

            // sort
            KokkosKernels::sort_crs_matrix<ExecSpace, OrdinalArray, OrdinalArray, ScalarArray>
              (uRowMap, uColIdx, uVal);
            ExecSpace().fence();
            #ifdef FASTILU_TIMER
            std::cout << "   + sort_matrix  " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif
#else
            using std::vector;
            vector<Ordinal> colPtrs(nRows, 0);
            for(Ordinal i = 0; i < nRows; i++) 
            {
                colPtrs[i] = uRowMap_[i];
            }
            for(Ordinal i = 0; i < nRows; i++)
            {
                for(Ordinal k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx_[k];
                    if (row <= col)
                    {
                        uVal_[colPtrs[col]] = aVal_[k];
                        uColIdx_[colPtrs[col]] = row;
                        colPtrs[col]++;
                        assert(colPtrs[col] <= uRowMap_[col+1]);
                    }
                }
            }
            Kokkos::deep_copy(uColIdx, uColIdx_);
            Kokkos::deep_copy(uVal, uVal_);
#endif

            if ((level > 0) && (guessFlag !=0))
            {
                // overwrite initial values from warmup runs
                ScalarArray gD;
                initGuessPrec->getD(gD);
                Kokkos::deep_copy(diagElems, gD);

                // copy UG into U
                OrdinalArray uGRowMap;
                OrdinalArray uGColIdx;
                ScalarArray uGVal;
                initGuessPrec->getU(uGRowMap, uGColIdx, uGVal);

                Kokkos::RangePolicy<CopySortedValsTag, ExecSpace> copy_policy (0, nRows);
                FastILUPrec_Functor functorG(uGVal, uGRowMap, uGColIdx, uVal, uRowMap, uColIdx);
                Kokkos::parallel_for(
                  "numericILU::copyVals(G)", copy_policy, functorG);
                ExecSpace().fence();
                #ifdef FASTILU_TIMER
                std::cout << "   + merge_sorted  " << Timer.seconds() << std::endl;
                Timer.reset();
                #endif
            }
        }

        void getL(OrdinalArray &lRowMapOut, OrdinalArray &lColIdxOut, ScalarArray &lValOut)
        {
            lRowMapOut = lRowMap;
            lColIdxOut = lColIdx;
            lValOut = lVal;
        }

        void getU(OrdinalArray &uRowMapOut, OrdinalArray &uColIdxOut, ScalarArray &uValOut)
        {
            uRowMapOut = uRowMap;
            uColIdxOut = uColIdx;
            uValOut = uVal;
        }

        void getUt(OrdinalArray &utRowMapOut, OrdinalArray &utColIdxOut, ScalarArray &utValOut)
        {
            utRowMapOut = utRowMap;
            utColIdxOut = utColIdx;
            utValOut = utVal;
        }

        void getD(ScalarArray &diagElemsOut)
        {
            diagElemsOut = diagElems;
        }

        void applyDiagonalScaling()
        {
            const Real one = RTS::one();
            int anext = 0;
            //First fill Aj and extract the diagonal scaling factors
            //Use diag array to store scaling factors since
            //it gets set to the correct value by findFactorPattern anyway.
            auto diagFact_ = Kokkos::create_mirror(diagFact);
            for (int i = 0; i < nRows; i++)
            {
                for(int k = aRowMap_[i]; k < aRowMap_[i+1]; k++) 
                {
                    aRowIdx_[anext++] = i;
                    if (aColIdx_[k] == i)
                    {
                        diagFact_[i] = one/(RTS::sqrt(STS::abs(aVal_[k])));
                    }
                }
            }

            //Now go through each element of A and apply the scaling
            int row;
            int col;
            Real sc1, sc2;

            for (int i = 0; i < nRows; i++) 
            {
                for (int k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    row = aRowIdx_[k];
                    col = aColIdx_[k];

                    sc1 = diagFact_[row];
                    sc2 = diagFact_[col];
                    aVal_[k] = aVal_[k]*sc1*sc2;
                }
            }
            Kokkos::deep_copy(diagFact, diagFact_);
        }

        void applyManteuffelShift()
        {
            const Scalar one = STS::one();
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx_[k];
                    if (row != col)
                    {
                        aVal_[k] = (one/(one + shift))*aVal_[k];
                    }
                }
            }
        }

        void applyD_Perm(ScalarArray &x, ScalarArray &y)
        {
            Kokkos::RangePolicy<NonTranPermScalTag, ExecSpace> scale_policy (0, nRows);
            PermScalFunctor<Ordinal, Scalar, Real, ExecSpace> functor(x, y, diagFact, permMetis);
            ExecSpace().fence();
            Kokkos::parallel_for(
              "numericILU::applyD_iPerm", scale_policy, functor);
            ExecSpace().fence();
        }

        void applyD_iPerm(ScalarArray &x, ScalarArray &y)
        {
            Kokkos::RangePolicy<TranPermScalTag, ExecSpace> scale_policy (0, nRows);
            PermScalFunctor<Ordinal, Scalar, Real, ExecSpace> functor(x, y, diagFact, ipermMetis);
            ExecSpace().fence();
            Kokkos::parallel_for(
              "numericILU::applyD_iPerm", scale_policy, functor);
            ExecSpace().fence();
        }

        void applyD(ScalarArray &x, ScalarArray &y)
        {
            ParScalFunctor<Ordinal, Scalar, Real, ExecSpace> parScal(x, y, diagFact);
            ExecSpace().fence();
            Kokkos::parallel_for(nRows, parScal);
            ExecSpace().fence();
        }

        void applyL(ScalarArray &x, ScalarArray &y)
        {
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(xOld);
            Kokkos::parallel_for(nRows, parInitZero);
            ExecSpace().fence();
#if 0
            JacobiIterFunctor<Ordinal, Scalar, ExecSpace> jacIter(nRows, lRowMap, lColIdx, lVal, x, y, xOld, onesVector); 
#endif 
            BlockJacobiIterFunctorL<Ordinal, Scalar, ExecSpace> jacIter(nRows, blkSz, lRowMap, lColIdx, lVal, x, y, xOld, onesVector);
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopy(xOld, y);
            Ordinal extent = nRows/blkSz;
            if (nRows%blkSz != 0)
            {
                extent++;
            }
            ExecSpace().fence();
            for (Ordinal i = 0; i < nTrisol; i++) 
            {
                //Kokkos::parallel_for(nRows, jacIter);
                Kokkos::parallel_for(extent, jacIter);
                ExecSpace().fence();
                Kokkos::parallel_for(nRows, parCopy);
                ExecSpace().fence();
            }
            return;
        }

        void applyU(ScalarArray &x, ScalarArray &y)
        {
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(xOld);
            Kokkos::parallel_for(nRows, parInitZero);
            ExecSpace().fence();
#if 0
            JacobiIterFunctor<Ordinal, Scalar, ExecSpace> jacIter(nRows, utRowMap, utColIdx, utVal, x, y, xOld, diagElems); 
#endif
            BlockJacobiIterFunctorU<Ordinal, Scalar, ExecSpace> jacIter(nRows, blkSz, utRowMap, utColIdx, utVal, x, y, xOld, diagElems);
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopy(xOld, y);
            Ordinal extent = nRows/blkSz; 
            if (nRows%blkSz != 0)
            {
                extent++;
            }
            ExecSpace().fence();
            for (Ordinal i = 0; i < nTrisol; i++) 
            {
                //Kokkos::parallel_for(nRows, jacIter);
                Kokkos::parallel_for(extent, jacIter);
                ExecSpace().fence();
                Kokkos::parallel_for(nRows, parCopy);
                ExecSpace().fence();
            }
            return;
        }


    public:
        //Constructor
        //TODO: Use a Teuchos::ParameterList object
        FastILUPrec(OrdinalArray &aRowMapIn_, OrdinalArray &aColIdxIn_, ScalarArray &aValIn_, Ordinal nRow_, FastILU::SpTRSV sptrsv_algo_,
                    Ordinal nFact_, Ordinal nTrisol_, Ordinal level_, Scalar omega_, Scalar shift_, Ordinal guessFlag_,
                    Ordinal blkSzILU_, Ordinal blkSz_)
        {
            nRows = nRow_;
            sptrsv_algo = sptrsv_algo_;
            nFact = nFact_;
            nTrisol = nTrisol_;

            useMetis = false;

            computeTime = 0.0;
            applyTime = 0.0;
            initTime = 0.0;
            //icFlag = icFlag_;
            level = level_;

            // mirror & deep-copy the input matrix
            aRowMapIn = aRowMapIn_;
            aColIdxIn = aColIdxIn_;
            aValIn    = aValIn_;
            aRowMapHost = Kokkos::create_mirror(aRowMapIn_);
            aColIdxHost = Kokkos::create_mirror(aColIdxIn_);
            aValHost = Kokkos::create_mirror(aValIn_);
            Kokkos::deep_copy(aRowMapHost, aRowMapIn_);
            Kokkos::deep_copy(aColIdxHost, aColIdxIn_);
            Kokkos::deep_copy(aValHost,    aValIn_);

            omega = omega_;
            guessFlag = guessFlag_;
            shift = shift_;
            blkSzILU = blkSzILU_;
            blkSz = blkSz_;
            doUnitDiag_TRSV = true; // perform TRSV with unit diagonals
            sptrsv_KKSpMV = true;   // use Kokkos-Kernels SpMV for Fast SpTRSV

            const Scalar one = STS::one();
            onesVector = ScalarArray("onesVector", nRow_);
            Kokkos::deep_copy(onesVector, one);

            diagFact = RealArray("diagFact", nRow_);
            diagElems = ScalarArray("diagElems", nRow_);
            xOld = ScalarArray("xOld", nRow_);
            xTemp = ScalarArray("xTemp", nRow_);

            if ((level > 0) && (guessFlag != 0))
            {
                initGuessPrec = Teuchos::rcp(new FastPrec(aRowMapIn_, aColIdxIn_, aValIn_, nRow_, sptrsv_algo_, 3, 5,
                                                          level_-1, omega_, shift_, guessFlag_, blkSzILU_, blkSz_));
            }
        }

        // internal functors
        struct ColPermTag {};
        struct CopySortedValsTag {};
        struct CopySortedValsPermTag {};
        struct CopyValsTag {};
        struct GetDiagsTag {};
        struct DiagScalTag {};
        struct SwapDiagTag {};

        struct GetLowerTag{};
        struct GetUpperTag{};
        struct FastILUPrec_Functor
        {
            FastILUPrec_Functor(ScalarArray aValIn_, OrdinalArray aRowMapIn_, OrdinalArray aColIdxIn_,
                                ScalarArray aVal_, RealArray diagFact_, OrdinalArray aRowMap_, OrdinalArray aColIdx_, OrdinalArray aRowIdx_) :
            aValIn (aValIn_),
            aRowMapIn (aRowMapIn_),
            aColIdxIn (aColIdxIn_),
            aVal (aVal_),
            diagFact (diagFact_),
            aRowMap (aRowMap_),
            aColIdx (aColIdx_),
            aRowIdx (aRowIdx_)
            {}

            // just calling CopySortedValsPermTag
            FastILUPrec_Functor(ScalarArray aValIn_, OrdinalArray aRowMapIn_, OrdinalArray aColIdxIn_, OrdinalArray perm_,
                                ScalarArray aVal_, RealArray diagFact_, OrdinalArray aRowMap_, OrdinalArray aColIdx_, OrdinalArray aRowIdx_) :
            aValIn (aValIn_),
            aRowMapIn (aRowMapIn_),
            aColIdxIn (aColIdxIn_),
            aVal (aVal_),
            diagFact (diagFact_),
            aRowMap (aRowMap_),
            aColIdx (aColIdx_),
            aRowIdx (aRowIdx_),
            iperm (perm_)
            {}

            // just calling CopySortedValsTag, or GetUpperTag
            FastILUPrec_Functor(ScalarArray aValIn_, OrdinalArray aRowMapIn_, OrdinalArray aColIdxIn_,
                                ScalarArray aVal_, OrdinalArray aRowMap_, OrdinalArray aColIdx_) :
            aValIn (aValIn_),
            aRowMapIn (aRowMapIn_),
            aColIdxIn (aColIdxIn_),
            aVal (aVal_),
            aRowMap (aRowMap_),
            aColIdx (aColIdx_)
            {}

            // just calling GetLowerTag
            FastILUPrec_Functor(ScalarArray aVal_, OrdinalArray aRowMap_, OrdinalArray aColIdx_,
                                ScalarArray lVal_, OrdinalArray lRowMap_, OrdinalArray lColIdx_,
                                ScalarArray diagElems_) :
            aVal (aVal_),
            diagElems (diagElems_),
            aRowMap (aRowMap_),
            aColIdx (aColIdx_),
            lVal (lVal_),
            lRowMap (lRowMap_),
            lColIdx (lColIdx_)
            {}

            // just calling SwapDiagTag
            FastILUPrec_Functor(const SwapDiagTag&, ScalarArray  lVal_, OrdinalArray  lRowMap_, OrdinalArray lColIdx_,
                                ScalarArray utVal_, OrdinalArray utRowMap_, OrdinalArray utColIdx_,
                                ScalarArray diagElems_) :
            diagElems (diagElems_),
            lVal (lVal_),
            lRowMap (lRowMap_),
            lColIdx (lColIdx_),
            utVal (utVal_),
            utRowMap (utRowMap_),
            utColIdx (utColIdx_)
            {}

            // just calling ColPerm
            FastILUPrec_Functor(OrdinalArray aColIdx_, OrdinalArray iperm_) :
            aColIdx (aColIdx_),
            iperm (iperm_)
            {}
            

            // ------------------------------------------------
            // functor to load values
            // both matrices are sorted and, "a" (with fills) contains "aIn" (original)
            KOKKOS_INLINE_FUNCTION
            void operator()(const CopySortedValsTag &, const int i) const {
                Ordinal aPtr = aRowMapIn[i];
                for(Ordinal k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    Ordinal col = aColIdx[k];
                    if (col == aColIdxIn[aPtr])
                    {
                        aVal[k] = aValIn[aPtr];
                        aPtr++;
                    } else
                    {
                        aVal[k] = STS::zero();
                    }
                }
            }

            // ------------------------------------------------
            // functor to load values with perm
            // both matrices are sorted and, "a" (with fills) contains "aIn" (original)
            KOKKOS_INLINE_FUNCTION
            void operator()(const CopySortedValsPermTag &, const int i) const {
                Ordinal aPtr = aRowMapIn[iperm[i]];
                for(Ordinal k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    Ordinal col = aColIdx[k];
                    if (col == aColIdxIn[aPtr])
                    {
                        aVal[k] = aValIn[aPtr];
                        aPtr++;
                    } else
                    {
                        aVal[k] = STS::zero();
                    }
                }
            }

            // functor to load values
            KOKKOS_INLINE_FUNCTION
            void operator()(const CopyValsTag &, const int i) const {
                for(Ordinal k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    Ordinal col = aColIdx[k];
                    for(Ordinal aPtr = aRowMapIn[i]; aPtr < aRowMapIn[i+1]; aPtr++)
                    {
                        if (col == aColIdxIn[aPtr])
                        {
                            aVal[k] = aValIn[aPtr];
                            break;
                        }
                    }
                }

            }

            // functor to extract diagonals (inverted)
            KOKKOS_INLINE_FUNCTION
            void operator()(const GetDiagsTag &, const int i) const {
                const Real one = RTS::one();
                for(int k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    aRowIdx[k] = i;
                    if (aColIdx[k] == i)
                    {
                        diagFact[i] = one/(RTS::sqrt(STS::abs(aVal[k])));
                    }
                }
            }

            // functor to swap diagonals
            KOKKOS_INLINE_FUNCTION
            void operator()(const SwapDiagTag &, const int i) const {
                const Scalar one  = STS::one();
                const Scalar zero = STS::zero();
                // zero the diagonal of L. If sorted, this finds it on first iter.
                Ordinal lRowBegin = lRowMap(i);
                Ordinal lRowEnd = lRowMap(i + 1);
                for(Ordinal j = 0; j < lRowEnd - lRowBegin; j++) {
                  Ordinal reversed = lRowEnd - j - 1;
                  if(lColIdx(reversed) == i) {
                    lVal(reversed) = zero;
                    break;
                  }
                }
                // zero the diagonal of Ut. If sorted, this finds it on first iter.
                Ordinal utRowBegin = utRowMap(i);
                Ordinal utRowEnd = utRowMap(i + 1);
                for(Ordinal j = utRowBegin; j < utRowEnd; j++) {
                  if(utColIdx(j) == i) {
                    utVal(j) = zero;
                    break;
                  }
                }
                // invert D
                diagElems[i] = one / diagElems[i];
            }

            // functor to apply diagonal scaling
            KOKKOS_INLINE_FUNCTION
            void operator()(const DiagScalTag &, const int i) const {
                for (int k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    int row = aRowIdx[k];
                    int col = aColIdx[k];

                    Real sc1 = diagFact[row];
                    Real sc2 = diagFact[col];
                    aVal[k] = aVal[k]*sc1*sc2;
                }
            }

            // ----------------------------------------------------------
            // functor to extract L & diagongals
            KOKKOS_INLINE_FUNCTION
            void operator()(const GetLowerTag &, const int i) const {
                Ordinal lPtr = lRowMap[i];
                for (Ordinal k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx[k];
                    if (row >= col)
                    {
                        if (row == col)
                        {
                            diagElems[row] = aVal[k];
                            lVal[lPtr] = STS::one();
                        } else {
                            lVal[lPtr] = aVal[k];
                        }
                        lColIdx[lPtr] = col;
                        lPtr++;
                    }
                }
            }

            // ----------------------------------------------------------
            // functor to extract U
            KOKKOS_INLINE_FUNCTION
            void operator()(const GetUpperTag &, const int i) const {
                for (Ordinal k = aRowMapIn[i]; k < aRowMapIn[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdxIn[k];
                    if (row <= col)
                    {
                        Ordinal pos = Kokkos::atomic_fetch_add(&(aRowMap[col]), 1);
                        aVal(pos) = aValIn[k];
                        aColIdx(pos) = row;
                    }
                }
            }

            // ----------------------------------------------------------
            // functor to apply column permutation
            KOKKOS_INLINE_FUNCTION
            void operator()(const ColPermTag &, const int i) const {
              aColIdx(i) = iperm(aColIdx(i));
            }

            // member variables
            // + input matrix
            ScalarArray    aValIn;
            OrdinalArray   aRowMapIn;
            OrdinalArray   aColIdxIn;
            // + output matrix
            ScalarArray    aVal;
            ScalarArray    diagElems;
            RealArray      diagFact;
            OrdinalArray   aRowMap;
            OrdinalArray   aColIdx;
            OrdinalArray   aRowIdx;
            // + output L matrix
            ScalarArray    lVal;
            OrdinalArray   lRowMap;
            OrdinalArray   lColIdx;
            // + output U matrix
            ScalarArray    utVal;
            OrdinalArray   utRowMap;
            OrdinalArray   utColIdx;
            // permutation
            OrdinalArray   iperm;
        };

        // set Metis pre-ordering
        template<class MetisArrayHost>
        void setMetisPerm(MetisArrayHost permMetis_, MetisArrayHost ipermMetis_)
        {
          Ordinal nRows_ = permMetis_.size();
          permMetis = OrdinalArray("permMetis", nRows_);
          ipermMetis = OrdinalArray("ipermMetis", nRows_);

          permMetisHost = Kokkos::create_mirror(permMetis);
          ipermMetisHost = Kokkos::create_mirror(ipermMetis);
          for (Ordinal i = 0; i < nRows_; i++) {
            permMetisHost(i) = permMetis_(i);
            ipermMetisHost(i) = ipermMetis_(i);
          }
          Kokkos::deep_copy(permMetis, permMetisHost);
          Kokkos::deep_copy(ipermMetis, ipermMetisHost);
          if ((level > 0) && (guessFlag != 0))
          {
            initGuessPrec->setMetisPerm(permMetis_, ipermMetis_);
          }

          useMetis = true;
        }

        //Symbolic Factorization Phase
        void initialize()
        {
            Kokkos::Timer timer;
            #ifdef FASTILU_INIT_TIMER
            Kokkos::Timer timer2;
            #endif
            // call symbolic that generates A with level associated to each nonzero entry
            // then pass that to initialize the initGuessPrec
            #if 1
            symbolicILU();
            #ifdef FASTILU_INIT_TIMER
            double tic = timer2.seconds();
            timer2.reset();
            std::cout << " + initial SymbolicILU (" << level << ") time : " << tic << std::endl;
            #endif
            if ((level > 0) && (guessFlag != 0))
            {
                initGuessPrec->initialize(aRowMap_, aColIdx_, aVal_, aLvlIdx_);
                #ifdef FASTILU_INIT_TIMER
                tic = timer2.seconds();
                timer2.reset();
                std::cout << "  > SymbolicILU (" << level << ") time : " << tic << std::endl;
                #endif
            }
            #else
            if ((level > 0) && (guessFlag != 0))
            {
                initGuessPrec->initialize();
            }
            symbolicILU();
            #endif
            //Allocate memory for the local A.
            //initialize L, U, A patterns
            #ifdef SHYLU_DEBUG
            Ordinal nnzU = uRowMap[nRows];
            MemoryPrimeFunctorN<Ordinal, Scalar, ExecSpace> copyFunc1(aRowMap, lRowMap, uRowMap, diagElems);
            MemoryPrimeFunctorNnzCsr<Ordinal, Scalar, ExecSpace> copyFunc4(uColIdx, uVal); 

            //Make sure all memory resides on the device
            ExecSpace().fence();
            Kokkos::parallel_for(nRows, copyFunc1);
            Kokkos::parallel_for(nnzU, copyFunc4);

            //Note that the following is a temporary measure
            //to ensure that memory resides on the device. 
            Ordinal nnzL = lRowMap[nRows];
            Ordinal nnzA = aRowMap[nRows];
            MemoryPrimeFunctorNnzCoo<Ordinal, Scalar, ExecSpace> copyFunc2(aColIdx, aRowIdx, aVal);
            MemoryPrimeFunctorNnzCsr<Ordinal, Scalar, ExecSpace> copyFunc3(lColIdx, lVal); 

            ExecSpace().fence();
            Kokkos::parallel_for(nRows, copyFunc1);
            Kokkos::parallel_for(nnzA, copyFunc2);
            Kokkos::parallel_for(nnzL, copyFunc3);
            #endif
            double t = timer.seconds();

            #ifdef FASTILU_INIT_TIMER
            std::cout << "Symbolic phase complete." << std::endl;
            std::cout << "Init time: "<< t << "s" << std::endl;
            #endif
            initTime = t;
            return;
        }

        //Symbolic Factorization Phase
        void initialize(OrdinalArrayMirror pRowMap_, OrdinalArrayMirror pColIdx_, ScalarArrayMirror pVal_, OrdinalArrayHost pLvlIdx_)
        {
            Kokkos::Timer timer;
            #ifdef FASTILU_INIT_TIMER
            Kokkos::Timer timer2;
            #endif
            // call symbolic that generates A with level associated to each nonzero entry
            // then pass that to initialize the initGuessPrec
            #if 1
            symbolicILU(pRowMap_, pColIdx_, pVal_, pLvlIdx_);
            #ifdef FASTILU_INIT_TIMER
            double tic = timer2.seconds();
            timer2.reset();
            std::cout << " - initial SymbolicILU (" << level << ") time : " << tic << std::endl;
            #endif
            if ((level > 0) && (guessFlag != 0))
            {
                initGuessPrec->initialize(pRowMap_, pColIdx_, pVal_, pLvlIdx_);
                #ifdef FASTILU_INIT_TIMER
                double tic = timer2.seconds();
                timer2.reset();
                std::cout << "  = SymbolicILU (" << level << ") time : " << tic << std::endl;
                #endif
            }
            #else
            if ((level > 0) && (guessFlag != 0))
            {
                initGuessPrec->initialize();
            }
            symbolicILU();
            #endif
            //Allocate memory for the local A.
            //initialize L, U, A patterns
            #ifdef SHYLU_DEBUG
            Ordinal nnzU = uRowMap[nRows];
            MemoryPrimeFunctorN<Ordinal, Scalar, ExecSpace> copyFunc1(aRowMap, lRowMap, uRowMap, diagElems);
            MemoryPrimeFunctorNnzCsr<Ordinal, Scalar, ExecSpace> copyFunc4(uColIdx, uVal); 

            //Make sure all memory resides on the device
            ExecSpace().fence();
            Kokkos::parallel_for(nRows, copyFunc1);
            Kokkos::parallel_for(nnzU, copyFunc4);

            //Note that the following is a temporary measure
            //to ensure that memory resides on the device. 
            Ordinal nnzL = lRowMap[nRows];
            Ordinal nnzA = aRowMap[nRows];
            MemoryPrimeFunctorNnzCoo<Ordinal, Scalar, ExecSpace> copyFunc2(aColIdx, aRowIdx, aVal);
            MemoryPrimeFunctorNnzCsr<Ordinal, Scalar, ExecSpace> copyFunc3(lColIdx, lVal); 

            ExecSpace().fence();
            Kokkos::parallel_for(nRows, copyFunc1);
            Kokkos::parallel_for(nnzA, copyFunc2);
            Kokkos::parallel_for(nnzL, copyFunc3);
            #endif
            double t = timer.seconds();

            #ifdef FASTILU_INIT_TIMER
            std::cout << " + Symbolic phase complete." << std::endl;
            std::cout << " + Init time: "<< t << "s" << std::endl;
            #endif
            initTime = t;
            return;
        }

        void setValues(ScalarArray& aValIn_)
        {
          this->aValIn = aValIn_;
          this->aValHost = Kokkos::create_mirror(aValIn_);
          Kokkos::deep_copy(this->aValHost, aValIn_);
          if(!initGuessPrec.is_null())
          {
            initGuessPrec->setValues(aValIn_);
          }
        }

        //Actual computation phase.
        //blkSzILU is the chunk size (hard coded).
        //1 gives the best performance on GPUs.
        //
        void compute()
        {
            Kokkos::Timer timer;
            #ifdef FASTILU_TIMER
            std::cout << "  >> compute <<" << std::endl;
            Kokkos::Timer Timer;
            Kokkos::Timer Timer2;
            #endif
            if ((level > 0) && (guessFlag !=0))
            {
                initGuessPrec->compute();
            }
            #ifdef FASTILU_TIMER
            std::cout << "  > initGuess " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif
            numericILU();
            #ifdef FASTILU_TIMER
            std::cout << "  > numericILU " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif
            FastILUFunctor<Ordinal, Scalar, ExecSpace> iluFunctor(aRowMap_[nRows], blkSzILU,
                    aRowMap, aRowIdx, aColIdx, aVal,
                    lRowMap, lColIdx, lVal, uRowMap, uColIdx, uVal, diagElems, omega);
            Ordinal extent = aRowMap_[nRows]/blkSzILU;
            if (aRowMap_[nRows]%blkSzILU != 0)
            {
                extent++;
            }
            //Ordinal extent = aRowMap[nRows];
            //ExecSpace().fence();

            for (int i = 0; i < nFact; i++) 
            {
                Kokkos::parallel_for(extent, iluFunctor);
            }
            ExecSpace().fence();
            #ifdef FASTILU_TIMER
            std::cout << "  > iluFunctor (" << nFact << ") " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif

            // transposee u
            double t = timer.seconds();
            #if 1
            // transpose
            Kokkos::deep_copy(utRowMap, 0);
            KokkosKernels::Impl::transpose_matrix<OrdinalArray, OrdinalArray, ScalarArray, OrdinalArray, OrdinalArray, ScalarArray, OrdinalArray, ExecSpace>
              (nRows, nRows, uRowMap, uColIdx, uVal, utRowMap, utColIdx, utVal);
            // sort, if the triangular solve algorithm requires a sorted matrix.
            // Currently, only Fast does not require this.
            bool sortRequired = sptrsv_algo != FastILU::SpTRSV::Fast;
            if(sortRequired) {
              KokkosKernels::sort_crs_matrix<ExecSpace, OrdinalArray, OrdinalArray, ScalarArray>
                (utRowMap, utColIdx, utVal);
            }
            if (sptrsv_algo == FastILU::SpTRSV::StandardHost) {
                // deep-copy to host
                Kokkos::deep_copy(lColIdx_, lColIdx);
                Kokkos::deep_copy(lVal_, lVal);

                Kokkos::deep_copy(utRowMap_, utRowMap);
                Kokkos::deep_copy(utColIdx_, utColIdx);
                Kokkos::deep_copy(utVal_, utVal);
            }
            #else
            // transposee u on host (need to copy to & from host)
            transposeU();
            #endif
            computeTime = t;
            #ifdef FASTILU_TIMER
            std::cout << "  > transposeU " << Timer.seconds() << std::endl;
            Timer.reset();
            #endif

            if (sptrsv_algo == FastILU::SpTRSV::Standard) {
                #if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
                KokkosSparse::Experimental::SPTRSVAlgorithm algo = KokkosSparse::Experimental::SPTRSVAlgorithm::SPTRSV_CUSPARSE;
                #else
                KokkosSparse::Experimental::SPTRSVAlgorithm algo = KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1;
                #endif
                // setup L solve
                khL.create_sptrsv_handle(algo, nRows, true);
                #if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
                KokkosSparse::Experimental::sptrsv_symbolic(&khL, lRowMap, lColIdx, lVal);
                #else
                KokkosSparse::Experimental::sptrsv_symbolic(&khL, lRowMap, lColIdx);
                #endif
                // setup U solve
                khU.create_sptrsv_handle(algo, nRows, false);
                #if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
                KokkosSparse::Experimental::sptrsv_symbolic(&khU, utRowMap, utColIdx, utVal);
                #else
                KokkosSparse::Experimental::sptrsv_symbolic(&khU, utRowMap, utColIdx);
                #endif
                #ifdef FASTILU_TIMER
                std::cout << "  > sptrsv_symbolic : nnz(L)=" << lColIdx.extent(0) << " nnz(U)=" << utColIdx.extent(0)
                          << ", " << Timer.seconds() << " seconds" << std::endl;
                #endif
            } else if (sptrsv_algo == FastILU::SpTRSV::StandardHost && doUnitDiag_TRSV) {
                // Prepare L for TRSV by removing unit-diagonals
                lVal_trsv_   = ScalarArrayHost ("lVal_trsv",    lRowMap_[nRows]-nRows);
                lColIdx_trsv_ = OrdinalArrayHost("lColIdx_trsv", lRowMap_[nRows]-nRows);
                lRowMap_trsv_ = OrdinalArrayHost("lRowMap_trsv", nRows+1);

                size_t nnzL = 0;
                lRowMap_trsv_(0) = 0;
                for (Ordinal i = 0; i < nRows; i++) {
                    for (Ordinal k = lRowMap_(i); k < lRowMap_[i+1]; k++) {
                        if (lColIdx_(k) != i) {
                            lVal_trsv_(nnzL) = lVal_(k);
                            lColIdx_trsv_(nnzL) = lColIdx_(k);
                            nnzL++;
                        }
                    }
                    lRowMap_trsv_(i+1)=nnzL;
                }

                // Prepare U by extracting and scaling D
                dVal_trsv_     = ScalarArrayHost ("dVal_trsv",     nRows);
                utVal_trsv_    = ScalarArrayHost ("utVal_trsv",    utRowMap_[nRows]-nRows);
                utColIdx_trsv_ = OrdinalArrayHost("utColIdx_trsv", utRowMap_[nRows]-nRows);
                utRowMap_trsv_ = OrdinalArrayHost("utRowMap_trsv", nRows+1);

                size_t nnzU = 0;
                utRowMap_trsv_(0) = 0;
                for (Ordinal i = 0; i < nRows; i++) {
                    for (Ordinal k = utRowMap_(i); k < utRowMap_[i+1]; k++) {
                        if (utColIdx_(k) == i) {
                            dVal_trsv_(i) = utVal_(k);
                        } else {
                            utVal_trsv_(nnzU) = utVal_(k);
                            utColIdx_trsv_(nnzU) = utColIdx_(k);
                            nnzU++;
                        }
                    }
                    utRowMap_trsv_(i+1)=nnzU;
                }
                for (Ordinal i = 0; i < nRows; i++) {
                    for (Ordinal k = utRowMap_trsv_(i); k < utRowMap_trsv_[i+1]; k++) {
                        utVal_trsv_(k) = utVal_trsv_(k) / dVal_trsv_(i);
                    }
                    dVal_trsv_(i) = STS::one() / dVal_trsv_(i);
                }
            } else if (sptrsv_KKSpMV) {
                FastILUPrec_Functor functor(SwapDiagTag(), lVal, lRowMap, lColIdx, utVal, utRowMap, utColIdx, diagElems);
                Kokkos::RangePolicy<SwapDiagTag, ExecSpace> swap_policy (0, nRows);
                Kokkos::parallel_for(
                  "numericILU::swapDiag", swap_policy, functor);
                ExecSpace().fence();
            }
            #ifdef FASTILU_TIMER
            std::cout << "  >> compute done " << Timer2.seconds() << " <<" << std::endl << std::endl;
            #endif
            return;
        }

        //Preconditioner application. Note that this does
        //*not* support multiple right hand sizes.
        void apply(ScalarArray &x, ScalarArray &y)
        {
            Kokkos::Timer timer;
            const Scalar one  = STS::one();

            //required to prevent contamination of the input.
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopyFunctor(xTemp, x);
            ExecSpace().fence();
            Kokkos::parallel_for(nRows, parCopyFunctor);
            ExecSpace().fence();
            //apply D
            if (useMetis) {
                applyD_Perm(x, xTemp);
            } else {
                applyD(x, xTemp);
            }
            if (sptrsv_algo == FastILU::SpTRSV::Standard) {
                // solve with L
                KokkosSparse::Experimental::sptrsv_solve(&khL, lRowMap, lColIdx, lVal, xTemp, y);
                // solve with U
                KokkosSparse::Experimental::sptrsv_solve(&khU, utRowMap, utColIdx, utVal, y, xTemp);
            } else {
                // wrap x and y into 2D views
                typedef Kokkos::View<Scalar **, ExecSpace> Scalar2dArray;
                Scalar2dArray x2d (const_cast<Scalar*>(xTemp.data()), nRows, 1);
                Scalar2dArray y2d (const_cast<Scalar*>(y.data()), nRows, 1);

                if (sptrsv_algo == FastILU::SpTRSV::StandardHost) {

                    // copy x to host
                    auto x_ = Kokkos::create_mirror(x2d);
                    auto y_ = Kokkos::create_mirror(y2d);
                    Kokkos::deep_copy(x_, x2d);

                    if (doUnitDiag_TRSV) {
                        using crsmat_host_t = KokkosSparse::CrsMatrix<Scalar, Ordinal, HostSpace, void, Ordinal>;
                        using graph_host_t  = typename crsmat_host_t::StaticCrsGraphType;

                        // wrap L into crsmat on host
                        graph_host_t static_graphL(lColIdx_trsv_, lRowMap_trsv_);
                        crsmat_host_t crsmatL("CrsMatrix", nRows, lVal_trsv_, static_graphL);

                        // wrap U into crsmat on host
                        graph_host_t static_graphU(utColIdx_trsv_, utRowMap_trsv_);
                        crsmat_host_t crsmatU("CrsMatrix", nRows, utVal_trsv_, static_graphU);

                        // solve with L, unit-diag
                        KokkosSparse::trsv ("L", "N", "U", crsmatL, x_, y_);
                        // solve with D
                        for (Ordinal i = 0; i < nRows; i++) {
                            y_(i, 0) = dVal_trsv_(i) * y_(i, 0);
                        }
                        // solve with U, unit-diag
                        KokkosSparse::trsv ("U", "N", "U", crsmatU, y_, x_);
                    } else {
                        using crsmat_mirror_t = KokkosSparse::CrsMatrix<Scalar, Ordinal, MirrorSpace, void, Ordinal>;
                        using graph_mirror_t  = typename crsmat_mirror_t::StaticCrsGraphType;

                        // wrap L into crsmat on host
                        graph_mirror_t static_graphL(lColIdx_, lRowMap_);
                        crsmat_mirror_t crsmatL("CrsMatrix", nRows, lVal_, static_graphL);

                        // wrap U into crsmat on host
                        graph_mirror_t static_graphU(utColIdx_, utRowMap_);
                        crsmat_mirror_t crsmatU("CrsMatrix", nRows, utVal_, static_graphU);

                        // solve with L
                        KokkosSparse::trsv ("L", "N", "N", crsmatL, x_, y_);
                        // solve with U
                        KokkosSparse::trsv ("U", "N", "N", crsmatU, y_, x_);
                    }
                    // copy x to device
                    Kokkos::deep_copy(x2d, x_);
                } else {
                    if (sptrsv_KKSpMV) {
                        using crsmat_t = KokkosSparse::CrsMatrix<Scalar, Ordinal, ExecSpace, void, Ordinal>;
                        using graph_t  = typename crsmat_t::StaticCrsGraphType;

                        graph_t static_graphL(lColIdx, lRowMap);
                        crsmat_t crsmatL("CrsMatrix", nRows, lVal, static_graphL);

                        graph_t static_graphU(utColIdx, utRowMap);
                        crsmat_t crsmatU("CrsMatrix", nRows, utVal, static_graphU);

                        Scalar2dArray x2d_old (const_cast<Scalar*>(xOld.data()), nRows, 1);

                        // 1) approximately solve, y = L^{-1}*x
                        // functor to copy RHS x into y (for even iteration)
                        ParCopyFunctor<Ordinal, Scalar, ExecSpace> copy_x2y(y, xTemp);
                        // functor to copy RHS x into xold (for odd iteration)
                        ParCopyFunctor<Ordinal, Scalar, ExecSpace> copy_x2xold(xOld, xTemp);

                        // functor to copy x_old to y (final iteration)
                        ParCopyFunctor<Ordinal, Scalar, ExecSpace> copy_xold2y(y, xOld);

                        // xold = zeros
                        ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> initZeroX(xOld);
                        Kokkos::parallel_for(nRows, initZeroX);
                        //Kokkos::deep_copy(x2d_old, zero);
                        for (Ordinal i = 0; i < nTrisol; i++) 
                        {
                            if (i%2 == 0) {
                                // y = y - L*x_old
                                // > y = x
                                Kokkos::parallel_for(nRows, copy_x2y);
                                ExecSpace().fence();
                                // > y = y - L*x_old
                                KokkosSparse::spmv("N", -one, crsmatL, x2d_old, one, y2d);
                            } else {
                                // x_old = x_old - L*y
                                // > x_old = x
                                Kokkos::parallel_for(nRows, copy_x2xold);
                                ExecSpace().fence();
                                // > x_old = x_old - L*y
                                KokkosSparse::spmv("N", -one, crsmatL, y2d, one, x2d_old);

                                if (i == nTrisol-1) {
                                    // y = x_old
                                    Kokkos::parallel_for(nRows, copy_xold2y);
                                    ExecSpace().fence();
                                }
                            }
                        }

                        // 2) approximately solve, x = U^{-1}*y
                        // functor to copy y into x
                        ParCopyFunctor<Ordinal, Scalar, ExecSpace> copy_y2x(xTemp, y);
                        // functor to copy y into xold
                        ParCopyFunctor<Ordinal, Scalar, ExecSpace> copy_y2xold(xOld, y);

                        // functor to scale x
                        ParScalFunctor<Ordinal, Scalar, Scalar, ExecSpace> scal_x(xTemp, xTemp, diagElems);
                        ParScalFunctor<Ordinal, Scalar, Scalar, ExecSpace> scal_xold(xOld, xOld, diagElems);

                        // functor to copy x_old to x (final iteration)
                        ParCopyFunctor<Ordinal, Scalar, ExecSpace> copy_xold2x(xTemp, xOld);

                        // xold = zeros
                        Kokkos::parallel_for(nRows, initZeroX);
                        //Kokkos::deep_copy(x2d_old, zero);
                        for (Ordinal i = 0; i < nTrisol; i++) 
                        {
                            if (i%2 == 0) {
                                // x = y - U*x_old
                                // > x = y
                                Kokkos::parallel_for(nRows, copy_y2x);
                                ExecSpace().fence();
                                // > x = x - U*x_old
                                KokkosSparse::spmv("N", -one, crsmatU, x2d_old, one, x2d);
                                // > scale x = inv(diag(U))*x
                                Kokkos::parallel_for(nRows, scal_x);
                            } else {
                                // xold = y - U*x
                                // > xold = y
                                Kokkos::parallel_for(nRows, copy_y2xold);
                                ExecSpace().fence();
                                // > x = x - U*x_old
                                KokkosSparse::spmv("N", -one, crsmatU, x2d, one, x2d_old);
                                // > scale x = inv(diag(U))*x
                                Kokkos::parallel_for(nRows, scal_xold);

                                if (i == nTrisol-1) {
                                    // x = x_old
                                    Kokkos::parallel_for(nRows, copy_xold2x);
                                    ExecSpace().fence();
                                }
                            }
                        }
                    } else {
                        //apply L^{-1} to xTemp
                        applyL(xTemp, y);
                        //apply U^{-1} to y
                        applyU(y, xTemp);
                    }
                }
            }
            //apply D again (we assume that the scaling is 
            //symmetric for now).
            if (useMetis) {
                applyD_iPerm(xTemp, y);
            } else {
                applyD(xTemp, y);
            }
            double t = timer.seconds();
            applyTime = t;
            return;
        }

        Ordinal getNFact() const
        {
            return nFact;
        }

        std::string getSpTrsvType() const
        {
            if (sptrsv_algo == FastILU::SpTRSV::StandardHost) {
                return "Standard Host";
            } else if (sptrsv_algo == FastILU::SpTRSV::Standard) {
                return "Standard";
            } else if (sptrsv_algo == FastILU::SpTRSV::Fast) {
                return "Fast";
            }
            return "Invalid";
        }

        Ordinal getNTrisol() const
        {
            return nTrisol;
        }

        Ordinal getNRows() const
        {
            return nRows;
        }

        double getComputeTime() const
        {
            return computeTime;
        }

        double getInitializeTime() const
        {
            return initTime;
        }

        double getApplyTime() const
        {
            return applyTime;
        }

        //Compute the L2 norm of the nonlinear residual (A - LU) on sparsity pattern
        void checkILU() const
        {
            Scalar sum = 0.0;
            Scalar sum_diag = 0.0;
            for (int i = 0 ; i < nRows; i++)
            {
                // Each row in A matrix (result)
                for (int k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    Scalar acc_val = aVal[k];
                    Ordinal lptr, uptr;
                    for ( lptr = lRowMap[i], uptr = uRowMap[aColIdx[k]] ;
                            lptr < lRowMap[i+1] && uptr < uRowMap[aColIdx[k]+1] ; )
                    {
                        if (lColIdx[lptr] == uColIdx[uptr])
                        {
                            acc_val -= lVal[lptr] * uVal[uptr];
                            lptr++;
                            uptr++;
                        }
                        else if (lColIdx[lptr] < uColIdx[uptr])
                            lptr++;
                        else
                            uptr++;
                    }
                    sum += acc_val * acc_val;
                }
            }
            
            for (int i = 0; i < nRows; i++) 
            {
                sum_diag += diagElems[i]*diagElems[i];
            }
            
            std::cout << "l2 norm of nonlinear residual = " << RTS::sqrt(STS::abs(sum)) << std::endl;
            std::cout << "l2 norm of diag. of U = " << RTS::sqrt(STS::abs(sum_diag)) << std::endl;
        }

        void checkIC() const
        {
            //Compute the L2 norm of the nonlinear residual (A - LLt) on sparsity pattern
            //
            Scalar sum = 0.0;
            for (int i = 0 ; i < nRows; i++)
            {
                int row = i;
                // Each row in A matrix (result)
                for (int k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {

                    int col = aColIdx[k];

                    if (row >= col) {
                        Scalar acc_val = aVal[k];
                        Ordinal lptr, uptr;
                        for ( lptr = lRowMap[i], uptr = lRowMap[aColIdx[k]] ;
                                lptr < lRowMap[i+1] && uptr < lRowMap[aColIdx[k]+1] ; )
                        {
                            if (lColIdx[lptr] == lColIdx[uptr])
                            {
                                acc_val -= lVal[lptr] * lVal[uptr];
                                lptr++;
                                uptr++;
                            }
                            else if (lColIdx[lptr] < lColIdx[uptr])
                                lptr++;
                            else
                                uptr++;
                        }

                        sum += acc_val * acc_val;
                    }
                }
            }
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "l2 norm of nonlinear residual = " << std::sqrt(sum) << std::endl;
            #endif
        }
        friend class FastILUFunctor<Ordinal, Scalar, ExecSpace>;
        friend class FastICFunctor<Ordinal, Scalar, ExecSpace>;
        friend class JacobiIterFunctor<Ordinal, Scalar, ExecSpace>;
        friend class ParCopyFunctor<Ordinal, Scalar, ExecSpace>;
        friend class JacobiIterFunctorT<Ordinal, Scalar, ExecSpace>;
        friend class ParScalFunctor<Ordinal, Scalar, Real, ExecSpace>;
        friend class PermScalFunctor<Ordinal, Scalar, Real, ExecSpace>;
        friend class MemoryPrimeFunctorN<Ordinal, Scalar, ExecSpace>;
        friend class MemoryPrimeFunctorNnzCoo<Ordinal, Scalar, ExecSpace>;
        friend class MemoryPrimeFunctorNnzCsr<Ordinal, Scalar, ExecSpace>;
        
};

//TODO: find a way to avoid the if condition (only store lower triangular part of A?)
template<class Ordinal, class Scalar, class ExecSpace>
class FastICFunctor
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        using STS = Kokkos::ArithTraits<Scalar>;

        FastICFunctor (Ordinal nNZ, Ordinal bs, ordinal_array_type Ap, ordinal_array_type Ai,
                ordinal_array_type Aj, scalar_array_type Ax, ordinal_array_type Lp,
                ordinal_array_type Li, scalar_array_type Lx, scalar_array_type diag, Scalar omega)
            :
                nnz(nNZ), blk_size(bs), _Ap(Ap), _Ai(Ai), _Aj(Aj),  _Lp(Lp), _Li(Li), _Ax(Ax), _Lx(Lx), _diag(diag), _omega(omega)
        {}
        
        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal blk_index) const
            {
                Ordinal start = blk_index * blk_size;
                Ordinal end = start + blk_size;

                Ordinal nz_index;

                if (end > nnz)
                {
                    end = nnz;
                }

                for (nz_index = start; nz_index < end && nz_index < nnz; nz_index++)
                {
                    Ordinal i = _Ai[nz_index];
                    Ordinal j = _Aj[nz_index];
                    //Ordinal temp = i;

                    Scalar val = _Ax[nz_index];
                    Scalar acc_val = 0.0;
                    Ordinal lptr = _Lp[i];
                    Ordinal ltptr = _Lp[j];
                    Ordinal endpt = j;
                    if (i >= j) { 

                        for ( ; _Li[lptr] < endpt && _Li[ltptr] < endpt; )
                        {
                            if (_Li[lptr] == _Li[ltptr])
                            {
                                acc_val += _Lx[lptr] * _Lx[ltptr];
                                lptr++;
                                ltptr++;
                            }
                            else if (_Li[lptr] < _Li[ltptr])
                            {
                                lptr++;
                            }
                            else 
                            {
                                ltptr++;
                            }
                        }
                        if (i > j) 
                        {
                            val = (val-acc_val) / _diag[j];
                            for ( ; _Li[lptr] < j ; lptr++) ; // dummy loop
                            assert (_Li[lptr] == j);
                            _Lx[lptr] = ((1.0 - _omega)*_Lx[lptr]) + (_omega*val);
                        }
                        else if (i == j)
                        {
                            //_diag[j] =  std::sqrt(val - acc_val);
                            val = STS::sqrt(val - acc_val);
                            _diag[j] = ((1.0 - _omega) * _diag[j]) + (_omega*val); 
                            for ( ; _Li[lptr] < j ; lptr++) ; // dummy loop
                            assert(_Li[lptr]==i);
                            _Lx[lptr] = _diag[j];
                        }
                    }
                }
            }

        Ordinal nnz, blk_size;
        ordinal_array_type _Ap, _Ai, _Aj, _Lp, _Li;
        scalar_array_type _Ax, _Lx, _diag;
        Scalar _omega;
};

template<class Ordinal, class Scalar, class ExecSpace>
class MemoryPrimeFunctorNnzCsr
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        MemoryPrimeFunctorNnzCsr (ordinal_array_type Ai, 
                scalar_array_type Ax)
            :
                _Ai(Ai),  _Ax(Ax)
        {}

        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal index) const
            {
                _Ai[index];
                _Ax[index];
            }

        ordinal_array_type _Ai; 
        scalar_array_type _Ax;
};

template<class Ordinal, class Scalar, class ExecSpace>
class MemoryPrimeFunctorNnzCoo
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        MemoryPrimeFunctorNnzCoo (ordinal_array_type Ai, 
                ordinal_array_type Aj,
                scalar_array_type Ax)
            :
                _Ai(Ai), _Aj(Aj), _Ax(Ax)
        {}

        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal index) const
            {
                /*
                Ordinal v1, v2;
                Scalar v3;
                */
                
                _Ai[index];
                _Aj[index];
                _Ax[index];
            }

        ordinal_array_type _Ai, _Aj; 
        scalar_array_type _Ax;
};

template<class Ordinal, class Scalar, class ExecSpace>
class MemoryPrimeFunctorN
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        MemoryPrimeFunctorN (ordinal_array_type Ap, 
                ordinal_array_type Lp,
                ordinal_array_type Up,
                scalar_array_type diag)
            :
                _Ap(Ap), _Lp(Lp), _Up(Up),
                 _diag(diag)
        {}

        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal index) const
            {
                //bmk: fix unused warnings?
                //does this functor actually have any side effects, or is it just incomplete?
                /*
                Ordinal v1, v2, v3;
                Scalar v4;
                */
                 
                _Ap[index];
                _Lp[index];
                _Up[index];
                _diag[index];
            }

        ordinal_array_type _Ap, _Lp, _Up;
        scalar_array_type _diag;
};

template<class Ordinal, class Scalar, class ExecSpace>
class FastILUFunctor
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        using STS = Kokkos::ArithTraits<Scalar>;

        FastILUFunctor (Ordinal nNZ, Ordinal bs, ordinal_array_type Ap, ordinal_array_type Ai,
                ordinal_array_type Aj, scalar_array_type Ax, ordinal_array_type Lp,
                ordinal_array_type Li, scalar_array_type Lx, ordinal_array_type Up,
                ordinal_array_type Ui, scalar_array_type Ux, scalar_array_type diag, Scalar omega)
            :
                nnz(nNZ), blk_size(bs), _Ap(Ap), _Ai(Ai), _Aj(Aj),  _Lp(Lp), _Li(Li),_Up(Up),
                _Ui(Ui), _Ax(Ax), _Lx(Lx), _Ux(Ux), _diag(diag), _omega(omega)
        {}

        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal blk_index) const
            {
                const Scalar zero = STS::zero();
                const Scalar one = STS::one();

                Ordinal start = blk_index * blk_size;
                Ordinal end = start + blk_size;

                Ordinal nz_index;

                if (end > nnz) 
                {
                    end = nnz;
                }

                for (nz_index = start; nz_index < end && nz_index < nnz; nz_index++)  
                {
                    Ordinal i = _Ai[nz_index];
                    Ordinal j = _Aj[nz_index];
                    Ordinal lCol;
                    Ordinal uCol;
                    Scalar val = _Ax[nz_index];
                    Scalar acc_val = zero;
                    Scalar lAdd = zero;
                    Ordinal lptr = _Lp[i];
                    Ordinal uptr = _Up[j];

                    while ( lptr < _Lp[i+1] && uptr < _Up[j+1] ) 
                    {
                        lCol = _Li[lptr];
                        uCol = _Ui[uptr];
                        lAdd = zero;
                        if (lCol == uCol)
                        {
                            lAdd = _Lx[lptr] * _Ux[uptr];
                            acc_val += lAdd; 
                        }
                        if (lCol <= uCol)
                        {
                            lptr++;
                        }
                        if (lCol >= uCol)
                        {
                            uptr++;
                        }
                    }

                    acc_val -= lAdd;

                    // Place the value into L or U
                    if (i > j) 
                    {
                        val = (val-acc_val) / _Ux[_Up[j+1]-1];
                        _Lx[lptr-1] = ((one - _omega) * _Lx[lptr-1]) + (_omega * val);
                    }
                    else
                    {
                        val = (val-acc_val);
                        if (i == j) _diag[j] = val;
                        _Ux[uptr-1] = ((one - _omega) * _Ux[uptr - 1]) + (_omega * val);
                    }
                }
            }

        Ordinal nnz, blk_size;
        ordinal_array_type _Ap, _Ai, _Aj, _Lp, _Li, _Up, _Ui;
        scalar_array_type _Ax, _Lx, _Ux, _diag;
        Scalar _omega;
};


template<class Ordinal, class Scalar, class ExecSpace>
class BlockJacobiIterFunctorL
{
    public:
        typedef ExecSpace ESpace;
        typedef Kokkos::View<Ordinal*, ExecSpace> OrdinalArray;
        typedef Kokkos::View<Scalar *, ExecSpace> ScalarArray;

        BlockJacobiIterFunctorL (Ordinal n, Ordinal bs, OrdinalArray aI, OrdinalArray aJ,
                ScalarArray aVal, ScalarArray b, ScalarArray xNew, 
                ScalarArray xOld, ScalarArray diag)
            :
                nRow(n), blkSize(bs), aRPtr(aI), aColIdx(aJ), aVal2(aVal), rhs(b), x2(xNew), x1(xOld),
                diagElems(diag)
        {}

        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal blkID) const
            {
                Ordinal idx1 = blkID * blkSize;
                Ordinal idx2 = idx1 + blkSize;
                Scalar val;
                Ordinal row;
                Ordinal col;
                Ordinal k;

                if (idx2 > nRow) 
                {
                    idx2 = nRow;
                }

                for (row = idx1; row < idx2; row++) 
                {
                    val = 0.0;
                    val = rhs[row];
                    for (k = aRPtr[row]; k < aRPtr[row+1]; k++) 
                    {
                        col = aColIdx[k];
                        if (col >= idx1 && col < row)
                        {
                            val -= aVal2[k]*x2[col];
                        }
                        else if (col < idx1 || col > row)
                        {
                            val -= aVal2[k]*x1[col];
                        }
                    }
                    x2[row] = val/diagElems[row];
                }
            }
        Ordinal nRow;
        Ordinal blkSize;
        OrdinalArray aRPtr, aColIdx;
        ScalarArray aVal2, rhs, x2, x1, diagElems;
};


template<class Ordinal, class Scalar, class ExecSpace>
class BlockJacobiIterFunctorU
{
    public:
        typedef ExecSpace ESpace;
        typedef Kokkos::View<Ordinal*, ExecSpace> OrdinalArray;
        typedef Kokkos::View<Scalar *, ExecSpace> ScalarArray;

        BlockJacobiIterFunctorU (Ordinal n, Ordinal bs, OrdinalArray aI, OrdinalArray aJ,
                ScalarArray aVal, ScalarArray b, ScalarArray xNew, 
                ScalarArray xOld, ScalarArray diag)
            :
                nRow(n), blkSize(bs), aRPtr(aI), aColIdx(aJ), aVal2(aVal), rhs(b), x2(xNew), x1(xOld),
                diagElems(diag)
        {}

        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal blkID) const
            {
                Ordinal idx1 = blkID * blkSize;
                Ordinal idx2 = idx1 + blkSize;
                Scalar val;
                Ordinal row;
                Ordinal col;
                Ordinal k;

                if (idx2 > nRow) 
                {
                    idx2 = nRow;
                }

                for (row = idx2 - 1; row >= idx1; row--) 
                {
                    val = 0.0;
                    val = rhs[row];
                    for (k = aRPtr[row]; k < aRPtr[row+1]; k++) 
                    {
                        col = aColIdx[k];
                        if (col < idx2 && col > row)
                        {
                            val -= aVal2[k]*x2[col];
                        }
                        else if (col >= idx2 || col < row)
                        {
                            val -= aVal2[k]*x1[col];
                        }
                    }
                    x2[row] = val/diagElems[row];
                }
            }
        Ordinal nRow;
        Ordinal blkSize;
        OrdinalArray aRPtr, aColIdx;
        ScalarArray aVal2, rhs, x2, x1, diagElems;
};

template<class Ordinal, class Scalar, class ExecSpace>
class JacobiIterFunctor
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        JacobiIterFunctor (Ordinal n, ordinal_array_type aI, ordinal_array_type aJ,
                scalar_array_type aVal, scalar_array_type b, scalar_array_type xNew,
                scalar_array_type xOld, scalar_array_type diag)
            :
                aI_(aI), aJ_(aJ), aVal_(aVal), b_(b), xNew_(xNew), xOld_(xOld), diag_(diag)
        {}

        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal xId) const
            {
                Scalar rowDot = 0.0;
                Ordinal k;

                //The equation is x_{k+1} = D^{-1}b + (I - D^{-1}A)x_{k}
                //The individual updates are x^{k+1}_{i} = b_{i}/d_{i} + x^{k}_{i} - 
                // \sum_{j = 1}^{n} r_{ij} x_{j}^{k}
                xNew_[xId] = b_[xId]/diag_[xId];
                xNew_[xId] += xOld_[xId];

                for (k = aI_[xId]; k < aI_[xId+1]; k++)
                {
                    rowDot += aVal_[k]*xOld_[aJ_[k]];
                }
                xNew_[xId] -= rowDot/diag_[xId];
            }

        ordinal_array_type aI_, aJ_; 
        scalar_array_type aVal_, b_, xNew_, xOld_, diag_;
};


//Parallel copy operation
template<class Ordinal, class Scalar, class ExecSpace>
class ParCopyFunctor
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        ParCopyFunctor (scalar_array_type xDestination, scalar_array_type xSource)
            :
                xDestination_(xDestination), xSource_(xSource)
        {}

        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal xId) const
            {
                xDestination_[xId] = xSource_[xId];
            }

        scalar_array_type xDestination_, xSource_;
};


template<class Ordinal, class Scalar, class ExecSpace>
class JacobiIterFunctorT
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        JacobiIterFunctorT (Ordinal n, ordinal_array_type aI, ordinal_array_type aJ,
                scalar_array_type aVal, scalar_array_type b, scalar_array_type xNew,
                scalar_array_type xOld, scalar_array_type diag)
            :
                aI_(aI), aJ_(aJ), aVal_(aVal), b_(b), xNew_(xNew), xOld_(xOld), diag_(diag), n_(n)
        {}

        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal xId) const
            {
                Ordinal k;

               // xNew_[xId] += b_[xId]/diag_[xId];
                //xNew_[xId] += xOld_[xId];

                Kokkos::atomic_add(&xNew_[xId], b_[xId]/diag_[xId]);
                Kokkos::atomic_add(&xNew_[xId], xOld_[xId]);

                for (k = aI_[xId]; k < aI_[xId+1]; k++)
                {

                    //y[aJ_[k]] += (aVal_[k]*xOld_[aJ_[k]])/diag_[aJ_[k]];
                    Kokkos::atomic_add(&xNew_[aJ_[k]], -(aVal_[k]*xOld_[xId])/diag_[aJ_[k]]);
                }
            }

        ordinal_array_type aI_, aJ_; 
        scalar_array_type aVal_, b_, xNew_, xOld_, diag_;
        Ordinal n_;
};

template<class Ordinal, class Scalar, class Real, class ExecSpace>
class ParScalFunctor
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;
        typedef Kokkos::View<Real *, ExecSpace> real_array_type;

        ParScalFunctor (scalar_array_type x, scalar_array_type y, real_array_type scaleFactors)
            :
                x_(x), y_(y), scaleFactors_(scaleFactors)
        {
        }


        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal xId) const
            {
               y_[xId] = x_[xId]*scaleFactors_[xId];
            }

        scalar_array_type x_, y_;
        real_array_type scaleFactors_;
};


template<class Ordinal, class Scalar, class Real, class ExecSpace>
class PermScalFunctor
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;
        typedef Kokkos::View<Real *, ExecSpace> real_array_type;

        PermScalFunctor (scalar_array_type x, scalar_array_type y, real_array_type scaleFactors, ordinal_array_type perm)
            :
                x_(x), y_(y), scaleFactors_(scaleFactors), perm_(perm)
        {
        }

        KOKKOS_INLINE_FUNCTION
            void operator()(const NonTranPermScalTag &, const Ordinal xId) const
            {
               // y = D*P*x 
               Ordinal row = perm_(xId);
               y_[xId] = scaleFactors_[xId]*x_[row];
            }

        KOKKOS_INLINE_FUNCTION
            void operator()(const TranPermScalTag &, const Ordinal xId) const
            {
               // y = P'*D*x 
               Ordinal row = perm_(xId);
               y_[xId] = x_[row]*scaleFactors_[row];
            }

        scalar_array_type x_, y_;
        real_array_type scaleFactors_;
        ordinal_array_type perm_;
};


template<class Ordinal, class Scalar, class ExecSpace>
class ParInitZeroFunctor
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        ParInitZeroFunctor(scalar_array_type x)
            :
                x_(x) 
    {
    }


        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal xId) const
            {
                x_[xId] = 0.0;
            }

        scalar_array_type x_ ;
};

#endif

