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

#include <assert.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#include <KokkosSparse_sptrsv.hpp>

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

template<class Ordinal, class Scalar, class ExecSpace>
class ParScalFunctor;

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
        typedef Kokkos::View<Ordinal *, ExecSpace> OrdinalArray;
        typedef Kokkos::View<Scalar *, ExecSpace> ScalarArray;
        typedef Kokkos::View<Ordinal *, typename ExecSpace::array_layout,
        Kokkos::Serial, Kokkos::MemoryUnmanaged> UMOrdinalArray;
        typedef Kokkos::View<Scalar *, typename ExecSpace::array_layout,
        Kokkos::Serial, Kokkos::MemoryUnmanaged> UMScalarArray;
        typedef FastILUPrec<Ordinal, Scalar, ExecSpace> FastPrec;

        typedef Kokkos::View<Ordinal *, Kokkos::HostSpace> OrdinalArrayHost;
        typedef Kokkos::View<Scalar  *, Kokkos::HostSpace>  ScalarArrayHost;
        typedef typename OrdinalArray::host_mirror_type OrdinalArrayMirror;
        typedef typename ScalarArray::host_mirror_type  ScalarArrayMirror;

    private:
        double computeTime;
        double applyTime;
        double initTime;

        Ordinal nRows;
        Ordinal guessFlag;
        Ordinal nFact;
        Ordinal nTrisol;
        Ordinal level;
        Ordinal blkSz;
        Scalar omega; //Underrelaxation parameter
        Scalar shift; //Manteuffel Shift

        //Lower triangular factor (CSR)
        ScalarArray lVal;
        OrdinalArray lColIdx;
        OrdinalArray lRowMap;
        // mirrors
        ScalarArrayMirror lVal_;
        OrdinalArrayMirror lColIdx_;
        OrdinalArrayMirror lRowMap_;

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

        //Pointer to the original host copy of A.
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

        //Diagonal scaling factors
        ScalarArray diagFact;
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
        bool standard_sptrsv;
        KernelHandle khL;
        KernelHandle khU;

        //Internal functions
        //Serial Transpose for now.
        //TODO:convert to parallel.
        void transposeU()
        {
            auto utRowMap_ = Kokkos::create_mirror(utRowMap);
            auto utColIdx_ = Kokkos::create_mirror(utColIdx);
            auto utVal_ = Kokkos::create_mirror(utVal);

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
        //Symbolic ILU code
        //initializes the matrices L and U and readies them
        //according to the level of fill
        void symbolicILU()
        {
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
            Ordinal n = nRows;
            *nzl = aRowMapHost[nRows] * (level + 2);
            *nzu = aRowMapHost[nRows] * (level + 2);
            Ordinal i;
            int levfill = level;
            vector<int> lnklst(n);
            vector<int> curlev(n);
            vector<int> levels(*nzu);
            vector<int> iwork(n);
            vector<int> ial(*nzl);
            vector<int> jal(*nzl);
            vector<int> iau(*nzu);
            vector<int> jau(*nzu);

            int knzl = 0;
            int knzu = 0;

            ial[0] = 0;
            iau[0] = 0;

            for (i=0; i<n; i++)
            {
                int first, next, j;

                /* copy column indices of row into workspace and sort them */

                int len = ia[i+1] - ia[i];
                next = 0;
                for (j=ia[i]; j<ia[i+1]; j++)
                    iwork[next++] = ja[j];
                // shell_sort(len, iwork);
                //stable_sort(iwork.begin(), iwork.begin() + len);
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

                    for (ii=iau[row]+1; ii<iau[row+1]; /*nop*/)
                    {
                        if (jau[ii] < nxtlst)
                        {
                            /* new fill-in */
                            int newlev = curlev[row] + levels[ii] + 1;
                            if (newlev <= levfill)
                            {
                                lnklst[oldlst]  = jau[ii];
                                lnklst[jau[ii]] = nxtlst;
                                oldlst = jau[ii];
                                curlev[jau[ii]] = newlev;
                            }
                            ii++;
                        }
                        else if (jau[ii] == nxtlst)
                        {
                            int newlev;
                            oldlst = nxtlst;
                            nxtlst = lnklst[oldlst];
                            newlev = curlev[row] + levels[ii] + 1;
                            //curlev[jau[ii]] = MIN(curlev[jau[ii]], newlev);
                            if (curlev[jau[ii]] > newlev)
                            {
                                curlev[jau[ii]] = newlev;
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

                next = first;
                while (next < i)
                {
                    assert(knzl < *nzl);
                    jal[knzl++] = next;
                    if (knzl >= *nzl)
                    {
                        //(*nzl)++;
                        *nzl = *nzl + nRows;
                        jal.resize(*nzl);
                    }
                    //jal.push_back(next);
                   // knzl++;
                    next = lnklst[next];
                }
                ial[i+1] = knzl;
                assert(next == i);

#if 0
                if (next != i)
                {
                    /*
                       assert(knzu < *nzu);
                       levels[knzu] = 2*n;
                       jau[knzu++] = i;
                       */
                }
#endif

                while (next < n)
                {
                    assert(knzu < *nzu);
                    levels[knzu] = curlev[next];
                    jau[knzu++] = next;
                    if (knzu >= *nzu)
                    {
                        //(*nzu)++;
                        *nzu = *nzu + nRows;
                        jau.resize(*nzu);
                        levels.resize(*nzu);
                    }
                    //jau.push_back(next);
                    //knzu++;
                    next = lnklst[next];
                }
                iau[i+1] = knzu;
            }

            *nzl = knzl;
            *nzu = knzu;

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
                    aRowPtr++;
                }
                for(Ordinal k = iau[i]; k < iau[i+1]; k++)
                {
                    aColIdx_[aRowPtr] = jau[k];
                    aRowIdx_[aRowPtr] = i;
                    aRowPtr++;
                }
                aRowMap_[i+1] = aRowPtr;
            }
            Kokkos::deep_copy(aRowMap, aRowMap_);
            Kokkos::deep_copy(aColIdx, aColIdx_);
            Kokkos::deep_copy(aRowIdx, aRowIdx_);
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished initializing A" << std::endl;
            #endif

            //Now allocate memory for L and U. 
            //
            lRowMap = OrdinalArray("lRowMap", nRows + 1);
            uRowMap = OrdinalArray("uRowMap", nRows + 1);
            utRowMap = OrdinalArray("utRowMap", nRows + 1);
            countL();
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished counting L" << std::endl;
            #endif
            
            countU();
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished counting U" << std::endl;
            #endif

            //Allocate memory and initialize pattern for L, U (transpose).
            lColIdx = OrdinalArray("lColIdx", lRowMap_[nRows]);
            uColIdx = OrdinalArray("uColIdx", uRowMap_[nRows]);
            utColIdx = OrdinalArray("utColIdx", uRowMap_[nRows]);
        }

        void numericILU()
        {
            aVal = ScalarArray("aVal", aColIdx.extent(0));
            aVal_ = Kokkos::create_mirror(aVal);
            //Copy the host matrix into the initialized a;
            Ordinal aHostPtr = 0;
            for (Ordinal i = 0; i < nRows; i++)
            {
                #ifdef SHYLU_DEBUG
                Ordinal check = aHostPtr;
                #endif
                for(Ordinal k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    Ordinal col = aColIdx_[k];
                    #ifdef FASTILU_DEBUG_OUTPUT
                    std::cout << "col =" << col << std::endl;
                    #endif
                    
                    if (col == aColIdxHost[aHostPtr])
                    {
                       aVal_[k] = aValHost[aHostPtr];
                       aHostPtr++;
                    }
                }
                #ifdef SHYLU_DEBUG
                assert((aHostPtr - check) == (aRowMapHost[i+1] - aRowMapHost[i]));
                #endif
            }

            lVal = ScalarArray("lVal", lRowMap_[nRows]);
            uVal = ScalarArray("uVal", uRowMap_[nRows]);
            utVal = ScalarArray("utVal", uRowMap_[nRows]);
            applyDiagonalScaling();
            applyManteuffelShift();

            Kokkos::deep_copy(aVal, aVal_);
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished diagonal scaling" << std::endl;
            #endif
            fillL();
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished copying L" << std::endl;
            #endif
            fillU();
            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "**Finished copying U" << std::endl;
            std::cout << "nnz L = " << lRowMap_[nRows] << std::endl;
            std::cout << "nnz U = " << uRowMap_[nRows] << std::endl;
            #endif
        }

        //Initialize the rowMap (rowPtr) and colIdx arrays for L
        void countL()
        {
            lRowMap_ = Kokkos::create_mirror(lRowMap);

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
        }

        //Put the initial guess into L.
        void fillL()
        {
            lVal_    = Kokkos::create_mirror(lVal);
            lColIdx_ = Kokkos::create_mirror(lColIdx);

            auto diagElems_ = Kokkos::create_mirror(diagElems);
            Ordinal lPtr = 0; 
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx_[k];

                    if (row >= col)
                    {
                        if (row == col) 
                        {
                            diagElems_[row] = aVal_[k];
                            lVal_[lPtr] = 1.0;
                            lColIdx_[lPtr] = aColIdx_[k];
                            lPtr++;
                        }
                        else 
                        {
                            lVal_[lPtr] = aVal_[k];
                            lColIdx_[lPtr] = aColIdx_[k];
                            lPtr++;
                        }
                    }
                }
            }
            assert(lPtr == lRowMap[nRows]);
            Kokkos::deep_copy(diagElems, diagElems_);

            if ((level > 0) && (guessFlag !=0))
            {
                OrdinalArray lGRowMap;
                OrdinalArray lGColIdx;
                ScalarArray lGVal;
                ScalarArray gD;
                Ordinal lGPtr = 0;
                initGuessPrec->getL(lGRowMap, lGColIdx, lGVal);
                initGuessPrec->getD(gD);

                Kokkos::deep_copy(diagElems, gD);

                auto lGColIdx_ = Kokkos::create_mirror(lGColIdx);
                auto lGVal_ = Kokkos::create_mirror(lGVal);
                Kokkos::deep_copy(lGColIdx_, lGColIdx);
                Kokkos::deep_copy(lGVal_, lGVal);
                for (Ordinal i = 0; i < nRows; i++) 
                {
                    #ifdef SHYLU_DEBUG
                    Ordinal check = lGPtr;
                    #endif
                    for (Ordinal k = lRowMap_[i]; k < lRowMap_[i+1]; k++)
                    {
                        //Ordinal row = i;
                        Ordinal col = lColIdx_[k];
                        if (col == lGColIdx_[lGPtr])
                        {
                            lVal_[k] = lGVal_[lGPtr];
                            lGPtr++;
                        }
                    }
                    #ifdef SHYLU_DEBUG
                    Ordinal rowLen = lGPtr - check;
                    assert(rowLen == lGRowMap[i+1] - lGRowMap[i]);
                    #endif
                }
            }
            Kokkos::deep_copy(lColIdx, lColIdx_);
            Kokkos::deep_copy(lVal, lVal_);
        }
        //Initialize rowMap and colIdx arrays of U
        void countU()
        {
            using std::vector;
            vector<Ordinal> colCounts(nRows + 1, 0);
            for(Ordinal i = 0; i < nRows; i++)
            {
                for(Ordinal k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx_[k];
                    if (row <= col)
                    {
                        colCounts[col+1]++;
                    }
                }
            }
            for (Ordinal i = 0; i < nRows; i++)
            {
                colCounts[i+1] += colCounts[i];
            }

            uRowMap_ = Kokkos::create_mirror(uRowMap);
            for (Ordinal i = 0; i <= nRows; i++)
            {
                uRowMap_[i] = colCounts[i];
            }
            Kokkos::deep_copy(uRowMap, uRowMap_);
        }

        //Put initial guess into U
        void fillU()
        {
            uVal_    = Kokkos::create_mirror(uVal);
            uColIdx_ = Kokkos::create_mirror(uColIdx);

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
            if ((level > 0) && (guessFlag !=0))
            {
                OrdinalArray uGRowMap;
                OrdinalArray uGColIdx;
                ScalarArray uGVal;
                ScalarArray gD;
                Ordinal uGPtr = 0;

                initGuessPrec->getU(uGRowMap, uGColIdx, uGVal);
                initGuessPrec->getD(gD);

                Kokkos::deep_copy(diagElems, gD);

                auto uGColIdx_ = Kokkos::create_mirror(uGColIdx);
                auto uGVal_ = Kokkos::create_mirror(uGVal);
                Kokkos::deep_copy(uGColIdx_, uGColIdx);
                Kokkos::deep_copy(uGVal_, uGVal);
                for (Ordinal i = 0; i < nRows; i++) 
                {
                    #ifdef SHYLU_DEBUG
                    Ordinal check = uGPtr;
                    #endif
                    for (Ordinal k = uRowMap_[i]; k < uRowMap_[i+1]; k++)
                    {
                        //unused: Ordinal row = i;
                        Ordinal col = uColIdx_[k];
                        if (col == uGColIdx_[uGPtr])
                        {
                            uVal_[k] = uGVal_[uGPtr];
                            uGPtr++;
                        }
                    }
                    #ifdef SHYLU_DEBUG
                    assert((uGPtr - check) == (uGRowMap[i+1] - uGRowMap[i]));
                    #endif
                }
            }
            Kokkos::deep_copy(uColIdx, uColIdx_);
            Kokkos::deep_copy(uVal, uVal_);
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

        void getD(ScalarArray &diagElemsOut)
        {
            diagElemsOut = diagElems;
        }

        void applyDiagonalScaling()
        {
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
                        diagFact_[i] = 1.0/std::sqrt(std::fabs(aVal_[k]));
                        //diagFact[i] = std::sqrt(std::fabs(aVal[k]));
                        #ifdef FASTILU_DEBUG_OUTPUT
                        std::cout << "diagFact["<<i<<"]="<<aVal_[k]<<std::endl;
                        #endif
                    }
                }
            }

            //Now go through each element of A and apply the scaling
            int row;
            int col;
            double sc1, sc2;

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
            //Scalar shift = 0.05;
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx_[k];
                    if (row != col)
                    {
                        aVal_[k] = (1.0/(1.0 + shift))*aVal_[k];
                    }
                }
            }
        }

        void applyD(ScalarArray &x, ScalarArray &y)
        {
            ParScalFunctor<Ordinal, Scalar, ExecSpace> parScal(nRows, x, y, diagFact);
            ExecSpace().fence();
            Kokkos::parallel_for(nRows, parScal);
            ExecSpace().fence();

        }

        void applyL(ScalarArray &x, ScalarArray &y)
        {
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(nRows, xOld);
            Kokkos::parallel_for(nRows, parInitZero);
            ExecSpace().fence();
#if 0
            JacobiIterFunctor<Ordinal, Scalar, ExecSpace> jacIter(nRows, lRowMap, lColIdx, lVal, x, y, xOld, onesVector); 
#endif 
            BlockJacobiIterFunctorL<Ordinal, Scalar, ExecSpace> jacIter(nRows, blkSz, lRowMap, lColIdx, lVal, x, y, xOld, onesVector);
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopy(nRows, xOld, y);
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
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(nRows, xOld);
            Kokkos::parallel_for(nRows, parInitZero);
            ExecSpace().fence();
#if 0
            JacobiIterFunctor<Ordinal, Scalar, ExecSpace> jacIter(nRows, utRowMap, utColIdx, utVal, x, y, xOld, diagElems); 
#endif
            BlockJacobiIterFunctorU<Ordinal, Scalar, ExecSpace> jacIter(nRows, blkSz, utRowMap, utColIdx, utVal, x, y, xOld, diagElems);
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopy(nRows, xOld, y);
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
        FastILUPrec(OrdinalArray &aRowMapIn, OrdinalArray &aColIdxIn, ScalarArray &aValIn, Ordinal nRow_, bool standard_sptrsv_,
                    Ordinal nFact_, Ordinal nTrisol_, Ordinal level_, Scalar omega_, Scalar shift_, Ordinal guessFlag_, Ordinal blkSz_)
        {
            nRows = nRow_;
            standard_sptrsv = standard_sptrsv_;
            nFact = nFact_;
            nTrisol = nTrisol_;
            computeTime = 0.0;
            applyTime = 0.0;
            initTime = 0.0;
            //icFlag = icFlag_;
            level = level_;

            // mirror & deep-copy the input matrix
            aRowMapHost = Kokkos::create_mirror(aRowMapIn);
            aColIdxHost = Kokkos::create_mirror(aColIdxIn);
            aValHost    = Kokkos::create_mirror(aValIn);
            Kokkos::deep_copy(aRowMapHost, aRowMapIn);
            Kokkos::deep_copy(aColIdxHost, aColIdxIn);
            Kokkos::deep_copy(aValHost,    aValIn);

            omega = omega_;
            guessFlag = guessFlag_;
            shift = shift_;
            blkSz = blkSz_;

            const Scalar one = Kokkos::ArithTraits<Scalar>::one();
            onesVector = ScalarArray("onesVector", nRow_);
            Kokkos::deep_copy(onesVector, one);

            diagFact = ScalarArray("diagFact", nRow_);
            diagElems = ScalarArray("diagElems", nRow_);
            xOld = ScalarArray("xOld", nRow_);
            xTemp = ScalarArray("xTemp", nRow_);

            if (level > 0)
            {
                initGuessPrec = Teuchos::rcp(new FastPrec(aRowMapIn, aColIdxIn, aValIn, nRow_, standard_sptrsv, 3, 5,
                                                          level_-1, omega_, shift_, guessFlag_, blkSz_));
            }
        }

        //Symbolic Factorization Phase
        void initialize()
        {
            Kokkos::Timer timer;
            if ((level > 0) && (guessFlag !=0))
            {
                initGuessPrec->initialize();
            }
            symbolicILU();
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

            #ifdef FASTILU_DEBUG_OUTPUT
            std::cout << "Symbolic phase complete." << std::endl;
            std::cout << "Init time: "<< t << "s" << std::endl;
            #endif
            initTime = t;
            return;
            
        }

        void setValues(ScalarArray& aValsIn)
        {
          this->aValHost = Kokkos::create_mirror(aValsIn);
          Kokkos::deep_copy(this->aValHost, aValsIn);
          if(!initGuessPrec.is_null())
          {
            initGuessPrec->setValues(aValsIn);
          }
        }

        //Actual computation phase.
        //blkSzILU is the chunk size (hard coded).
        //1 gives the best performance on GPUs.
        //
        void compute()
        {
            Kokkos::Timer timer;
            if ((level > 0) && (guessFlag !=0))
            {
                initGuessPrec->compute();
            }
            numericILU();
            Ordinal blkSzILU = 4096;
            FastILUFunctor<Ordinal, Scalar, ExecSpace> iluFunctor(aRowMap_[nRows], blkSzILU,
                    aRowMap, aColIdx, aRowIdx, aVal, 
                    lRowMap, lColIdx, lVal, uRowMap, uColIdx, uVal, diagElems, omega);
            Ordinal extent = aRowMap_[nRows]/blkSzILU;
            if (aRowMap_[nRows]%blkSzILU != 0)
            {
                extent++;
            }
            
            //Ordinal extent = aRowMap[nRows];
            ExecSpace().fence();

            for (int i = 0; i < nFact; i++) 
            {
                Kokkos::parallel_for(extent, iluFunctor);
            }
            //ExecSpace().fence();

            // transposee u on host (need to copy to & from host)
            double t = timer.seconds();
            transposeU();
            computeTime = t;

            if (standard_sptrsv) {
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
            }
            return;
        }

        //Preconditioner application. Note that this does
        //*not* support multiple right hand sizes.
        void apply(ScalarArray &x, ScalarArray &y)
        {
            Kokkos::Timer timer;

            //required to prevent contamination of the input.
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopyFunctor(nRows, xTemp, x);
            ExecSpace().fence();
            Kokkos::parallel_for(nRows, parCopyFunctor);
            ExecSpace().fence();

            //apply D
            applyD(x, xTemp);
            if (standard_sptrsv) {
                // solve with L
                KokkosSparse::Experimental::sptrsv_solve(&khL, lRowMap, lColIdx, lVal, xTemp, y);
                // solve with U
                KokkosSparse::Experimental::sptrsv_solve(&khU, utRowMap, utColIdx, utVal, y, xTemp);
            } else {
                //apply L
                applyL(xTemp, y);
                //apply U or Lt depending on icFlag
                applyU(y, xTemp);
            }
            //apply D again (we assume that the scaling is 
            //symmetric for now).
            applyD(xTemp, y);
            double t = timer.seconds();
            applyTime = t;
            return;
        }

        Ordinal getNFact() const
        {
            return nFact;
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
            
            std::cout << "l2 norm of nonlinear residual = " << std::sqrt(sum) << std::endl;
            std::cout << "l2 norm of diag. of U = " << std::sqrt(sum_diag) << std::endl;
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
        friend class ParScalFunctor<Ordinal, Scalar, ExecSpace>;
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
                            val = std::sqrt(val - acc_val);
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
                    Scalar acc_val = 0.0;
                    Scalar lAdd = 0.0;
                    Ordinal lptr = _Lp[i];
                    Ordinal uptr = _Up[j];

                    while ( lptr < _Lp[i+1] && uptr < _Up[j+1] ) 
                    {
                        lCol = _Li[lptr];
                        uCol = _Ui[uptr];
                        lAdd = 0.0;
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
                        _Lx[lptr-1] = ((1 - _omega) * _Lx[lptr-1]) + (_omega * val);
                    }
                    else
                    {
                        val = (val-acc_val);
                        if (i == j) _diag[j] = val;
                        _Ux[uptr-1] = ((1 - _omega) * _Ux[uptr - 1]) + (_omega * val);
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

        ParCopyFunctor (Ordinal n, scalar_array_type xDestination, scalar_array_type xSource)
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

template<class Ordinal, class Scalar, class ExecSpace>
class ParScalFunctor
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        ParScalFunctor (Ordinal n, scalar_array_type x, scalar_array_type y, scalar_array_type scaleFactors)
            :
                x_(x), y_(y), scaleFactors_(scaleFactors)
        {
        }


        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal xId) const
            {
               // x_[xId] *= scaleFactors_[xId];
               y_[xId] = x_[xId]*scaleFactors_[xId];
            }

        scalar_array_type x_, y_, scaleFactors_;
};


template<class Ordinal, class Scalar, class ExecSpace>
class ParInitZeroFunctor
{
    public:
        typedef ExecSpace execution_space;
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        ParInitZeroFunctor(Ordinal n, scalar_array_type x)
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

