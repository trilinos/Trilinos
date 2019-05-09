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


// The struct that iterates over all non-zeros for the FastIC
//  Contact for bugs and complaints - Siva Rajamanickam (srajama@sandia.gov)
//
#ifndef __FAST_IC_HPP__
#define __FAST_IC_HPP__
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

#include "shylu_fastilu.hpp"

//whether to print extra debug output at runtime to stdout
//comment out next line to disable
//#define FASTIC_DEBUG_OUTPUT

template<class Ordinal, class Scalar, class ExecSpace>
class FastICPrec
{

    public:
        typedef Kokkos::View<Ordinal *, ExecSpace> OrdinalArray;
        typedef Kokkos::View<Scalar *, ExecSpace> ScalarArray;
        typedef Kokkos::View<Ordinal *, typename ExecSpace::array_layout, Kokkos::Serial, 
                Kokkos::MemoryUnmanaged> UMOrdinalArray;
        typedef Kokkos::View<Scalar *, typename ExecSpace::array_layout, Kokkos::Serial, 
                Kokkos::MemoryUnmanaged> UMScalarArray;
        typedef FastICPrec<Ordinal, Scalar, ExecSpace> FastPrec;

    private:
        double computeTime;
        double applyTime;
        double initTime;

        Ordinal guessFlag;
        Ordinal nRows;
        Ordinal nFact;
        Ordinal nTrisol;
        Ordinal level;
        Ordinal blkSz;
        Scalar omega;
        Scalar shift;

        //Lower triangular factor (CSR)
        ScalarArray lVal;
        OrdinalArray lColIdx;
        OrdinalArray lRowMap;

        //Upper triangular factor (CSR)
        ScalarArray ltVal;
        OrdinalArray ltColIdx;
        OrdinalArray ltRowMap;

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
        ScalarArray diagElemsInv;

        //Temporary vectors for triangular solves
        ScalarArray xOld;
        ScalarArray xTemp;
        ScalarArray onesVector;

        Teuchos::RCP<FastPrec> initGuessPrec;


    public:
        FastICPrec(OrdinalArray &aRowMap_, OrdinalArray &aColIdx_, ScalarArray &aVal_, Ordinal nRow_,
                Ordinal nFact_, Ordinal nTrisol_, Ordinal level_, Scalar omega_, 
                Scalar shift_, Ordinal guessFlag_, Ordinal blkSz_)
        {
            nRows = nRow_;
            nFact = nFact_;
            nTrisol = nTrisol_;
            computeTime = 0.0;
            applyTime = 0.0;
            initTime = 0.0;
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
            diagElemsInv = ScalarArray("diagElems", nRow_);
            xOld = ScalarArray("xOld", nRow_);
            xTemp = ScalarArray("xTemp", nRow_);

            if (level > 0)
            {
                initGuessPrec = Teuchos::rcp(new FastPrec(aRowMap_, aColIdx_, aVal_, nRow_, 3, 5, 
                            level_ - 1, omega_, shift_, guessFlag_, blkSz_));
            }

        }

        void initialize()
        {
            Kokkos::Impl::Timer timer;
            //Here we have to allocate memory for L, U.
            if((level > 0) && (guessFlag != 0))
            {
                initGuessPrec->initialize();
            }
            //symbolicILUAndInit(); //Note the diagonal scaling is now embedded in this routine.
            symbolicILU();
            //initialize L, U, A patterns
            #ifdef SHYLU_DEBUG
            MemoryPrimeFunctorN<Ordinal, Scalar, ExecSpace> copyFunc1(aRowMap, lRowMap, lRowMap, diagElems);

            //Make sure all memory resides on the device
            ExecSpace::fence();
            Kokkos::parallel_for(nRows, copyFunc1);

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

            #ifdef FASTIC_DEBUG_OUTPUT
            std::cout << "Symbolic phase complete." << std::endl;
            std::cout << "Init time: "<< t << "s" << std::endl;
            #endif
            initTime = t;
        }

        void getL(OrdinalArray &lRowMap_, OrdinalArray &lColIdx_, ScalarArray &lVal_)
        {
            lRowMap_ = lRowMap;
            lColIdx_ = lColIdx;
            lVal_ = lVal;
        }

        void getD(ScalarArray &diagElems_)
        {
            diagElems_ = diagElems;
        }

        void transposeL()
        {
            //Count the elements in each row of Lt
            auto temp = OrdinalArray("temp", nRows + 1);
            auto rowPtrs = OrdinalArray("rowPtrs", nRows);
            for (Ordinal i = 0; i <= nRows; i++) 
            {
                temp[i] = 0;
            }
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = lRowMap[i]; k < lRowMap[i+1]; k++) 
                {
                    temp[lColIdx[k]+1]++;

                }
            }
            //Perform an add scan to get the row map for 
            //the transpose
            for (Ordinal i = 0; i <= nRows; i++) 
            {
                ltRowMap[i] = temp[i];
            }
            for (Ordinal i = 1; i <= nRows; i++) 
            {
                ltRowMap[i] += ltRowMap[i-1];
            }
            //Set the row pointers to their initial places;
            for (Ordinal i = 0; i < nRows; i++) 
            {
                rowPtrs[i] = ltRowMap[i];
            }
            //Copy the data
            for (Ordinal i = 0; i < nRows; i++) 
            {
                for (Ordinal k = lRowMap[i]; k < lRowMap[i+1]; k++)
                {
                    Ordinal row = lColIdx[k];
                    Scalar value = lVal[k];
                    ltVal[rowPtrs[row]] = value;
                    ltColIdx[rowPtrs[row]] = i;
                    rowPtrs[row]++;
                    assert(rowPtrs[row] <= ltRowMap[row + 1]);
                }
            }
        }

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
            int n = nRows;
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

            #ifdef FASTIC_DEBUG_OUTPUT
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
                #ifdef FASTIC_DEBUG_OUTPUT
                std::cout << "***row:" << i << std::endl;
                #endif
                for(Ordinal k = ial[i]; k < ial[i+1]; k++)
                {
                    #ifdef FASTIC_DEBUG_OUTPUT
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
            //Now allocate memory for L and U. 
            //
            lRowMap = OrdinalArray("lRowMap", nRows + 1);
            ltRowMap = OrdinalArray("ltRowMap", nRows + 1);
            countL();
            //Allocate memory and initialize pattern for L, U (transpose).
            lColIdx = OrdinalArray("lColIdx", lRowMap[nRows]);
            ltColIdx = OrdinalArray("ltColIdx", lRowMap[nRows]);
            #ifdef FASTIC_DEBUG_OUTPUT
            std::cout << "nnz L = " << lRowMap[nRows] << std::endl;
            #endif
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
                    #ifdef FASTIC_DEBUG_OUTPUT
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
            ltVal = ScalarArray("ltVal", lRowMap[nRows]);
            applyDiagonalScaling();
            applyManteuffelShift();
            fillL();
        }

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
                        }
                        lVal[lPtr] = aVal[k];
                        lColIdx[lPtr] = aColIdx[k];
                        lPtr++;
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
                        //unused: Ordinal row = i;
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
                        #ifdef FASTIC_DEBUG_OUTPUT
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

        void applyDD(ScalarArray &x, ScalarArray &y)
        {
            ParScalFunctor<Ordinal, Scalar, ExecSpace> parScal(nRows, x, y, diagElemsInv);
            ExecSpace::fence();
            Kokkos::parallel_for(nRows, parScal);
            ExecSpace::fence();

        }
        void applyD(ScalarArray &x, ScalarArray &y)
        {
            ParScalFunctor<Ordinal, Scalar, ExecSpace> parScal(nRows, x, y, diagFact);
            ExecSpace::fence();
            Kokkos::parallel_for(nRows, parScal);
            ExecSpace::fence();

        }
        void applyLIC(ScalarArray &x, ScalarArray &y)
        {
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(nRows, xOld);
            Kokkos::parallel_for(nRows, parInitZero);
            ExecSpace::fence();
#if 0
            JacobiIterFunctor<Ordinal, Scalar, ExecSpace> jacIter(nRows, lRowMap, lColIdx, lVal, 
                                                                  x, y, xOld, diagElems); 
#endif 
            BlockJacobiIterFunctorL<Ordinal, Scalar, ExecSpace> jacIter(nRows, blkSz, lRowMap, 
                                                                        lColIdx, lVal,
                                                                        x, y, xOld, diagElems);
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopy(nRows, xOld, y);
            Ordinal extent = nRows/blkSz;
            if(nRows%blkSz != 0)
            {
                extent++;
            }
            ExecSpace::fence();
            for (Ordinal i = 0; i < nTrisol; i++) 
            {
                Kokkos::parallel_for(extent, jacIter);
                ExecSpace::fence();
                Kokkos::parallel_for(nRows, parCopy);
                ExecSpace::fence();
            }
            return;
        }

        void applyLT(ScalarArray &x, ScalarArray &y)
        {
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(nRows, xOld);
            Kokkos::parallel_for(nRows, parInitZero);
            ExecSpace::fence();
#if 0
            JacobiIterFunctor<Ordinal, Scalar, ExecSpace> jacIter(nRows, ltRowMap, ltColIdx, ltVal, 
                                                                  x, y, xOld, diagElems); 
#endif
            BlockJacobiIterFunctorU<Ordinal, Scalar, ExecSpace> jacIter(nRows, blkSz, ltRowMap, ltColIdx,
                                                                        ltVal, x, y, xOld, diagElems);
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopy(nRows, xOld, y);
            Ordinal extent = nRows/blkSz;
            if (nRows%blkSz != 0)
            {
                extent++;
            }
            ExecSpace::fence();
            for (Ordinal i = 0; i < nTrisol; i++) 
            {
                Kokkos::parallel_for(extent, jacIter);
                ExecSpace::fence();
                Kokkos::parallel_for(nRows, parCopy);
                ExecSpace::fence();
            }
            return;
        }


#if 0
        void applyLT(ScalarArray &x, ScalarArray &y)
        {
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(nRows, xOld);
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero1(nRows, y);
            Kokkos::parallel_for(nRows, parInitZero);
            ExecSpace::fence();
            //Transpose Jacobi implementation
            JacobiIterFunctorT<Ordinal, Scalar, ExecSpace> jacIterT(nRows, lRowMap, lColIdx, lVal, x, y, xOld, onesVector);
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopy(nRows, xOld, y);

            ExecSpace::fence();

            for (Ordinal i = 0; i < nTrisol; i++) 
            {
                Kokkos::parallel_for(nRows, parInitZero1);
                ExecSpace::fence();
                Kokkos::parallel_for(nRows, jacIterT);
                ExecSpace::fence();
                Kokkos::parallel_for(nRows, parCopy);
                ExecSpace::fence();
            }

            return;
        }
#endif

        void setValues(ScalarArray& aValsIn)
        {
          this->aValHost = aValsIn;
          if(!initGuessPrec.is_null())
          {
            initGuessPrec->setValues(aValsIn);
          }
        }

        void compute()
        {
            Kokkos::Impl::Timer timer;
            if((level > 0) && (guessFlag != 0))
            {
                initGuessPrec->compute();
            }
            numericILU();
            FastICFunctor<Ordinal, Scalar, ExecSpace> icFunctor(nRows, aRowMap, aColIdx, aRowIdx, aVal,
                    lRowMap, lColIdx, lVal, diagElems, omega);
            ExecSpace::fence();

            for (int i = 0; i < nFact; i++) 
            {
                Kokkos::parallel_for(aRowMap[nRows], icFunctor);
            }
            ExecSpace::fence();

            double t = timer.seconds();
            computeTime = t;

            for (int i = 0; i < nRows; i++) 
            {
                diagElemsInv[i] = 1.0/diagElems[i];
            }
            transposeL();
        }

        void apply(ScalarArray &x, ScalarArray &y)
        {

            Kokkos::Impl::Timer timer;
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopyFunctor(nRows, xTemp, x);
            ExecSpace::fence();
            Kokkos::parallel_for(nRows, parCopyFunctor);
            ExecSpace::fence();

            applyD(x, xTemp);
            applyLIC(xTemp, y);
            //applyDD(y, xTemp);
            applyLT(y, xTemp);
            applyD(xTemp, y);

//            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopyFunctor2(nRows, y, xTemp);
//            ExecSpace::fence();
//            Kokkos::parallel_for(nRows, parCopyFunctor2);
//            ExecSpace::fence();

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
            #ifdef FASTIC_DEBUG_OUTPUT
            std::cout << "l2 norm of nonlinear residual = " << std::sqrt(sum) << std::endl;
            #endif
        }
        friend class FastICFunctor<Ordinal, Scalar, ExecSpace>;
        friend class JacobiIterFunctor<Ordinal, Scalar, ExecSpace>;
        friend class ParCopyFunctor<Ordinal, Scalar, ExecSpace>;
        friend class JacobiIterFunctorT<Ordinal, Scalar, ExecSpace>;
        friend class ParScalFunctor<Ordinal, Scalar, ExecSpace>;
        friend class MemoryPrimeFunctorN<Ordinal, Scalar, ExecSpace>;
        friend class MemoryPrimeFunctorNnzCoo<Ordinal, Scalar, ExecSpace>;
        friend class MemoryPrimeFunctorNnzCsr<Ordinal, Scalar, ExecSpace>;
};


#endif
