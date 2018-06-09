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
#include <impl/Kokkos_Timer.hpp>

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

        //Upper triangular factor (CSC)
        ScalarArray uVal;
        OrdinalArray uColIdx;
        OrdinalArray uRowMap;

        //Upper triangular factor (CSR)
        ScalarArray utVal;
        OrdinalArray utColIdx;
        OrdinalArray utRowMap;

        //Pointer to the original host copy of A.
        ScalarArray aValHost;
        OrdinalArray aRowMapHost;
        OrdinalArray aColIdxHost;

        //A matrix in COO format
        ScalarArray aVal;
        OrdinalArray aRowMap;
        OrdinalArray aRowIdx;
        OrdinalArray aColIdx;

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

        //Internal functions
        //Serial Transpose for now.
        //TODO:convert to parallel.
        void transposeU()
        {
            //Count the elements in each row of Ut
            auto temp = OrdinalArray("temp", nRows + 1);
            auto rowPtrs = OrdinalArray("rowPtrs", nRows);
            for (Ordinal i = 0; i <= nRows; i++) 
            {
                temp[i] = 0;
            }
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = uRowMap[i]; k < uRowMap[i+1]; k++) 
                {
                    temp[uColIdx[k]+1]++;

                }
            }
            //Perform an add scan to get the row map for 
            //the transpose
            for (Ordinal i = 0; i <= nRows; i++) 
            {
                utRowMap[i] = temp[i];
            }
            for (Ordinal i = 1; i <= nRows; i++) 
            {
                utRowMap[i] += utRowMap[i-1];
            }
            //Set the row pointers to their initial places;
            for (Ordinal i = 0; i < nRows; i++) 
            {
                rowPtrs[i] = utRowMap[i];
            }
            //Copy the data
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = uRowMap[i]; k < uRowMap[i+1]; k++)
                {
                    Ordinal row = uColIdx[k];
                    Scalar value = uVal[k];
                    utVal[rowPtrs[row]] = value;
                    utColIdx[rowPtrs[row]] = i;
                    rowPtrs[row]++;
                    assert(rowPtrs[row] <= utRowMap[row + 1]);
                }
            }
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
            OrdinalArray ia = aRowMapHost;
            OrdinalArray ja = aColIdxHost;
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
            Ordinal aRowPtr = 0;

            aRowMap[0] = 0;
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
                    aColIdx[aRowPtr] = jal[k];
                    aRowIdx[aRowPtr] = i;
                    aRowPtr++;
                }
                for(Ordinal k = iau[i]; k < iau[i+1]; k++)
                {
                    aColIdx[aRowPtr] = jau[k];
                    aRowIdx[aRowPtr] = i;
                    aRowPtr++;
                }
                aRowMap[i+1] = aRowPtr;
            }

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
            lColIdx = OrdinalArray("lColIdx", lRowMap[nRows]);
            uColIdx = OrdinalArray("uColIdx", uRowMap[nRows]);
            utColIdx = OrdinalArray("utColIdx", uRowMap[nRows]);
        }

        void numericILU()
        {
            aVal = ScalarArray("aVal", aColIdx.extent(0));
            //Copy the host matrix into the initialized a;
            Ordinal aHostPtr = 0;
            Ordinal check = 0;
            for (Ordinal i = 0; i < nRows; i++)
            {
                check = aHostPtr;
                for(Ordinal k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    Ordinal col = aColIdx[k];
                    #ifdef FASTILU_DEBUG_OUTPUT
                    std::cout << "col =" << col << std::endl;
                    #endif
                    
                    if (col == aColIdxHost[aHostPtr])
                    {
                       aVal[k] = aValHost[aHostPtr]; 
                       aHostPtr++;
                    }
                }
                assert((aHostPtr - check) == (aRowMapHost[i+1] - aRowMapHost[i]));
            }
            lVal = ScalarArray("lVal", lRowMap[nRows]);
            uVal = ScalarArray("uVal", uRowMap[nRows]);
            utVal = ScalarArray("utVal", uRowMap[nRows]);
            applyDiagonalScaling();
            applyManteuffelShift();
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
            std::cout << "nnz L = " << lRowMap[nRows] << std::endl;
            std::cout << "nnz U = " << uRowMap[nRows] << std::endl;
            #endif
        }

        //Initialize the rowMap (rowPtr) and colIdx arrays for L
        void countL()
        {
            lRowMap[0] = 0;
            for (Ordinal i = 0; i < nRows; i++) 
            {
                Ordinal row_count = 0;
                for (Ordinal k = aRowMap[i]; k < aRowMap[i+1]; k++) 
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx[k];

                    if (row >= col)
                    {
                       row_count++; 
                    }
                }
                lRowMap[i+1] = lRowMap[i] + row_count;
            }
        }

        //Put the initial guess into L.
        void fillL()
        {
            Ordinal lPtr = 0; 
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = aRowMap[i]; k < aRowMap[i+1]; k++) 
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx[k];

                    if (row >= col)
                    {
                        if (row == col) 
                        {
                            diagElems[row] = aVal[k];
                            lVal[lPtr] = 1.0;
                            lColIdx[lPtr] = aColIdx[k];
                            lPtr++;
                        }
                        else 
                        {
                            lVal[lPtr] = aVal[k];
                            lColIdx[lPtr] = aColIdx[k];
                            lPtr++;
                        }
                    }
                }
            }
            assert(lPtr == lRowMap[nRows]);
            if ((level > 0) && (guessFlag !=0))
            {
                OrdinalArray lGRowMap;
                OrdinalArray lGColIdx;
                ScalarArray lGVal;
                ScalarArray gD;
                Ordinal lGPtr = 0;
                Ordinal check = 0;
                Ordinal rowLen;

                initGuessPrec->getL(lGRowMap, lGColIdx, lGVal);
                initGuessPrec->getD(gD);

                for (Ordinal i = 0; i < nRows; i++)
                {
                    diagElems[i] = gD[i];
                }
                for (Ordinal i = 0; i < nRows; i++) 
                {
                    check = lGPtr;
                    for (Ordinal k = lRowMap[i]; k < lRowMap[i+1]; k++)
                    {
                        //Ordinal row = i;
                        Ordinal col = lColIdx[k];

                        if (col == lGColIdx[lGPtr])
                        {
                            lVal[k] = lGVal[lGPtr];
                            lGPtr++;
                        }
                    }
                    rowLen = lGPtr - check;
                    assert(rowLen == lGRowMap[i+1] - lGRowMap[i]);
                }
            }
        }
        //Initialize rowMap and colIdx arrays of U
        void countU()
        {
            using std::vector;
            vector<Ordinal> colCounts(nRows + 1, 0);
            for(Ordinal i = 0; i < nRows; i++)
            {
                for(Ordinal k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx[k];
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

            for (Ordinal i = 0; i <= nRows; i++)
            {
                uRowMap[i] = colCounts[i];
            }
        }

        //Put initial guess into U
        void fillU()
        {
            using std::vector;
            vector<Ordinal> colPtrs(nRows, 0);
            for(Ordinal i = 0; i < nRows; i++) 
            {
                colPtrs[i] = uRowMap[i];
            }
            for(Ordinal i = 0; i < nRows; i++)
            {
                for(Ordinal k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx[k];
                    if (row <= col)
                    {
                        uVal[colPtrs[col]] = aVal[k];
                        uColIdx[colPtrs[col]] = row;
                        colPtrs[col]++;
                        assert(colPtrs[col] <= uRowMap[col+1]);
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
                Ordinal check = 0;
                Ordinal rowLen;

                initGuessPrec->getU(uGRowMap, uGColIdx, uGVal);
                initGuessPrec->getD(gD);

                for (Ordinal i = 0; i < nRows; i++)
                {
                    diagElems[i] = gD[i];
                }
                for (Ordinal i = 0; i < nRows; i++) 
                {
                    check = uGPtr;
                    for (Ordinal k = uRowMap[i]; k < uRowMap[i+1]; k++)
                    {
                        //unused: Ordinal row = i;
                        Ordinal col = uColIdx[k];

                        if (col == uGColIdx[uGPtr])
                        {
                            uVal[k] = uGVal[uGPtr];
                            uGPtr++;
                        }
                    }
                    rowLen = uGPtr - check;
                    assert(rowLen == uGRowMap[i+1] - uGRowMap[i]);
                }
            }

        }

        void getL(OrdinalArray &lRowMap_, OrdinalArray &lColIdx_, ScalarArray &lVal_)
        {
            lRowMap_ = lRowMap;
            lColIdx_ = lColIdx;
            lVal_ = lVal;
        }

        void getU(OrdinalArray &uRowMap_, OrdinalArray &uColIdx_, ScalarArray &uVal_)
        {
            uRowMap_ = uRowMap;
            uColIdx_ = uColIdx;
            uVal_ = uVal;
        }

        void getD(ScalarArray &diagElems_)
        {
            diagElems_ = diagElems;
        }

        void applyDiagonalScaling()
        {
            int anext = 0;
            //First fill Aj and extract the diagonal scaling factors
            //Use diag array to store scaling factors since
            //it gets set to the correct value by findFactorPattern anyway.
            for (int i = 0; i < nRows; i++) 
            {
                for(int k = aRowMap[i]; k < aRowMap[i+1]; k++) 
                {
                    aRowIdx[anext++] = i;
                    if (aColIdx[k] == i) 
                    {
                        diagFact[i] = 1.0/std::sqrt(std::fabs(aVal[k]));
                        //diagFact[i] = std::sqrt(std::fabs(aVal[k]));
                        #ifdef FASTILU_DEBUG_OUTPUT
                        std::cout << "diagFact["<<i<<"]="<<aVal[k]<<std::endl;
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
                for (int k = aRowMap[i]; k < aRowMap[i+1]; k++) 
                {
                    row = aRowIdx[k];
                    col = aColIdx[k];

                    sc1 = diagFact[row];
                    sc2 = diagFact[col];
                    aVal[k] = aVal[k]*sc1*sc2;
                }
            }
        }

        void applyManteuffelShift()
        {
            //Scalar shift = 0.05;
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = aRowMap[i]; k < aRowMap[i+1]; k++)
                {
                    Ordinal row = i;
                    Ordinal col = aColIdx[k];
                    if (row != col)
                    {
                        aVal[k] = (1.0/(1.0 + shift))*aVal[k];
                    }
                }
            }
        }

        void applyD(ScalarArray &x, ScalarArray &y)
        {
            ParScalFunctor<Ordinal, Scalar, ExecSpace> parScal(nRows, x, y, diagFact);
            ExecSpace::fence();
            Kokkos::parallel_for(nRows, parScal);
            ExecSpace::fence();

        }

        void applyL(ScalarArray &x, ScalarArray &y)
        {


            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(nRows, xOld);
            Kokkos::parallel_for(nRows, parInitZero);
            ExecSpace::fence();
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
            ExecSpace::fence();
            for (Ordinal i = 0; i < nTrisol; i++) 
            {
                //Kokkos::parallel_for(nRows, jacIter);
                Kokkos::parallel_for(extent, jacIter);
                ExecSpace::fence();
                Kokkos::parallel_for(nRows, parCopy);
                ExecSpace::fence();
            }
            return;
        }

        void applyU(ScalarArray &x, ScalarArray &y)
        {
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(nRows, xOld);
            Kokkos::parallel_for(nRows, parInitZero);
            ExecSpace::fence();
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
            ExecSpace::fence();
            for (Ordinal i = 0; i < nTrisol; i++) 
            {
                //Kokkos::parallel_for(nRows, jacIter);
                Kokkos::parallel_for(extent, jacIter);
                ExecSpace::fence();
                Kokkos::parallel_for(nRows, parCopy);
                ExecSpace::fence();
            }
            return;
        }


    public:
        //Constructor
        //TODO: Use a Teuchos::ParameterList object
        FastILUPrec(OrdinalArray &aRowMap_, OrdinalArray &aColIdx_, ScalarArray &aVal_, Ordinal nRow_,
                Ordinal nFact_, Ordinal nTrisol_, Ordinal level_, Scalar omega_, Scalar shift_, Ordinal guessFlag_, Ordinal blkSz_)
        {
            nRows = nRow_;
            nFact = nFact_;
            nTrisol = nTrisol_;
            computeTime = 0.0;
            applyTime = 0.0;
            initTime = 0.0;
            //icFlag = icFlag_;
            level = level_;
            aRowMapHost = aRowMap_;
            aColIdxHost = aColIdx_;
            aValHost = aVal_;
            omega = omega_;
            guessFlag = guessFlag_;
            shift = shift_;
            blkSz = blkSz_;

            onesVector = ScalarArray("onesVector", nRow_);
            for (int i = 0; i < nRow_; i++)
            {
                onesVector[i] = 1.0;
            }

            diagFact = ScalarArray("diagFact", nRow_);
            diagElems = ScalarArray("diagElems", nRow_);
            xOld = ScalarArray("xOld", nRow_);
            xTemp = ScalarArray("xTemp", nRow_);

            if (level > 0)
            {
                initGuessPrec = Teuchos::rcp(new FastPrec(aRowMap_, aColIdx_, aVal_, nRow_, 3, 5, 
                            level_ - 1, omega_, shift_, guessFlag_, blkSz_));
            }
        }

        //Symbolic Factorization Phase
        void initialize()
        {
            Kokkos::Impl::Timer timer;
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
            ExecSpace::fence();
            Kokkos::parallel_for(nRows, copyFunc1);
            Kokkos::parallel_for(nnzU, copyFunc4);

            //Note that the following is a temporary measure
            //to ensure that memory resides on the device. 
            Ordinal nnzL = lRowMap[nRows];
            Ordinal nnzA = aRowMap[nRows];
            MemoryPrimeFunctorNnzCoo<Ordinal, Scalar, ExecSpace> copyFunc2(aColIdx, aRowIdx, aVal);
            MemoryPrimeFunctorNnzCsr<Ordinal, Scalar, ExecSpace> copyFunc3(lColIdx, lVal); 

            ExecSpace::fence();
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
          this->aValHost = aValsIn;
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
            Kokkos::Impl::Timer timer;
            if ((level > 0) && (guessFlag !=0))
            {
                initGuessPrec->compute();
            }
            numericILU();
            Ordinal blkSzILU = 4096;
            FastILUFunctor<Ordinal, Scalar, ExecSpace> iluFunctor(aRowMap[nRows], blkSzILU, aRowMap, 
                    aColIdx, aRowIdx, aVal, 
                    lRowMap, lColIdx, lVal, uRowMap, uColIdx, uVal, diagElems, omega);
            Ordinal extent = aRowMap[nRows]/blkSzILU;
            if (aRowMap[nRows]%blkSzILU != 0)
            {
                extent++;
            }
            
            //Ordinal extent = aRowMap[nRows];
            ExecSpace::fence();

            for (int i = 0; i < nFact; i++) 
            {
                Kokkos::parallel_for(extent, iluFunctor);
            }
            //ExecSpace::fence();

            double t = timer.seconds();
            transposeU();
            computeTime = t;
            return;
        }

        //Preconditioner application. Note that this does
        //*not* support multiple right hand sizes.
        void apply(ScalarArray &x, ScalarArray &y)
        {
            Kokkos::Impl::Timer timer;

            //required to prevent contamination of the input.
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopyFunctor(nRows, xTemp, x);
            ExecSpace::fence();
            Kokkos::parallel_for(nRows, parCopyFunctor);
            ExecSpace::fence();

            //apply D
            applyD(x, xTemp);
            //apply L
            applyL(xTemp, y);
            //apply U or Lt depending on icFlag
            applyU(y, xTemp);
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

        FastICFunctor (Ordinal nvtx, ordinal_array_type Ap, ordinal_array_type Ai,
                ordinal_array_type Aj, scalar_array_type Ax, ordinal_array_type Lp,
                ordinal_array_type Li, scalar_array_type Lx, scalar_array_type diag, Scalar omega)
            :
                _Ap(Ap), _Ai(Ai), _Aj(Aj),  _Lp(Lp), _Li(Li), _Ax(Ax), _Lx(Lx), _diag(diag), _omega(omega)
    {
    }
        
        KOKKOS_INLINE_FUNCTION
            void operator()(const Ordinal nz_index) const
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
    {
    }
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
    {
    }
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
    {
    }
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
    {
    }
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
    {
    }

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
    {
    }

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
    {
    } 

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
        typedef Kokkos::View<Ordinal *, ExecSpace> ordinal_array_type;
        typedef Kokkos::View<Scalar *, ExecSpace> scalar_array_type;

        ParCopyFunctor (Ordinal n, scalar_array_type xDestination, scalar_array_type xSource)
            :
                xDestination_(xDestination), xSource_(xSource)
    {
    }


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
    {
    } 

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

