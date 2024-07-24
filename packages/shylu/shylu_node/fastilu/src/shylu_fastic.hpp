// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include <string>

#include <assert.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

#include "shylu_fastilu.hpp"

//whether to print extra debug output at runtime to stdout
//comment out next line to disable
//#define FASTIC_DEBUG_OUTPUT

template<class Ordinal, class Scalar, class ExecSpace>
class FastICPrec
{

    public:
        typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Real;
        typedef Kokkos::View<Ordinal *, ExecSpace> OrdinalArray;
        typedef Kokkos::View<Scalar *, ExecSpace> ScalarArray;
        typedef Kokkos::View<Real *, ExecSpace> RealArray;
        typedef Kokkos::View<Ordinal *, typename ExecSpace::array_layout, Kokkos::Serial,
                             Kokkos::MemoryUnmanaged> UMOrdinalArray;
        typedef Kokkos::View<Scalar *, typename ExecSpace::array_layout, Kokkos::Serial,
                             Kokkos::MemoryUnmanaged> UMScalarArray;
        typedef FastICPrec<Ordinal, Scalar, ExecSpace> FastPrec;

        typedef Kokkos::View<Ordinal *, Kokkos::HostSpace> OrdinalArrayHost;
        typedef Kokkos::View<Scalar  *, Kokkos::HostSpace>  ScalarArrayHost;
        typedef typename OrdinalArray::host_mirror_type OrdinalArrayMirror;
        typedef typename ScalarArray::host_mirror_type  ScalarArrayMirror;

        typedef Kokkos::RangePolicy<ExecSpace> RangePolicy;

        using STS = Kokkos::ArithTraits<Scalar>;
        using RTS = Kokkos::ArithTraits<Real>;

    private:
        double computeTime;
        double applyTime;
        double initTime;

        Ordinal guessFlag;
        Ordinal nRows;
        Ordinal nFact;
        Ordinal nTrisol;
        Ordinal level;
        Ordinal blkSzIC;
        Ordinal blkSz;
        Scalar omega;
        Scalar shift;

        //Lower triangular factor (CSR)
        ScalarArray lVal;
        OrdinalArray lColIdx;
        OrdinalArray lRowMap;
        // mirrors
        ScalarArrayMirror lVal_;
        OrdinalArrayMirror lColIdx_;
        OrdinalArrayMirror lRowMap_;

        //Upper triangular factor (CSR)
        ScalarArray ltVal;
        OrdinalArray ltColIdx;
        OrdinalArray ltRowMap;

        //Pointer to the original host copy of A.
        ScalarArrayMirror aValHost;
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
        RealArray diagFact;
        ScalarArray diagElems;
        ScalarArray diagElemsInv;
        // mirrors
        ScalarArrayMirror diagElems_;
        ScalarArrayMirror diagElemsInv_;

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
        KernelHandle khLt;

    public:
        FastICPrec(OrdinalArray &aRowMapIn, OrdinalArray &aColIdxIn, ScalarArray &aValIn, Ordinal nRow_, bool standard_sptrsv_,
                   Ordinal nFact_, Ordinal nTrisol_, Ordinal level_, Scalar omega_, Scalar shift_, Ordinal guessFlag_, Ordinal blkSzIC_, Ordinal blkSz_)
        {
            nRows = nRow_;
            standard_sptrsv = standard_sptrsv_;
            nFact = nFact_;
            nTrisol = nTrisol_;
            computeTime = 0.0;
            applyTime = 0.0;
            initTime = 0.0;
            level = level_;

            // mirror & deep-copy the input matrix
            aRowMapHost = Kokkos::create_mirror_view(aRowMapIn);
            aColIdxHost = Kokkos::create_mirror_view(aColIdxIn);
            aValHost    = Kokkos::create_mirror_view(aValIn);
            Kokkos::deep_copy(aRowMapHost, aRowMapIn);
            Kokkos::deep_copy(aColIdxHost, aColIdxIn);
            Kokkos::deep_copy(aValHost,    aValIn);

            omega = omega_;
            guessFlag = guessFlag_;
            shift = shift_;
            blkSzIC = blkSzIC_;
            blkSz = blkSz_;

            const Scalar one = STS::one();
            onesVector = ScalarArray("onesVector", nRow_);
            Kokkos::deep_copy(onesVector, one);

            diagFact = RealArray("diagFact", nRow_);
            diagElems = ScalarArray("diagElems", nRow_);
            diagElemsInv = ScalarArray("diagElems", nRow_);
            xOld = ScalarArray("xOld", nRow_);
            xTemp = ScalarArray("xTemp", nRow_);

            if (level > 0)
            {
                initGuessPrec = Teuchos::rcp(new FastPrec(aRowMapIn, aColIdxIn, aValIn, nRow_, standard_sptrsv, 3, 5,
                                                          level_-1, omega_, shift_, guessFlag_, blkSzIC_, blkSz_));
            }

        }

        void initialize()
        {
            Kokkos::Timer timer;
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
            ExecSpace().fence();
            Kokkos::parallel_for(RangePolicy(0, nRows), copyFunc1);

            //Note that the following is a temporary measure
            //to ensure that memory resides on the device.
            Ordinal nnzL = lRowMap[nRows];
            Ordinal nnzA = aRowMap[nRows];

            MemoryPrimeFunctorNnzCoo<Ordinal, Scalar, ExecSpace> copyFunc2(aColIdx, aRowIdx, aVal);
            MemoryPrimeFunctorNnzCsr<Ordinal, Scalar, ExecSpace> copyFunc3(lColIdx, lVal);

            ExecSpace().fence();
            Kokkos::parallel_for(RangePolicy(0, nRows), copyFunc1);
            Kokkos::parallel_for(RangePolicy(0, nnzA), copyFunc2);
            Kokkos::parallel_for(RangePolicy(0, nnzL), copyFunc3);
            #endif
            double t = timer.seconds();

            #ifdef FASTIC_DEBUG_OUTPUT
            std::cout << "Symbolic phase complete." << std::endl;
            std::cout << "Init time: "<< t << "s" << std::endl;
            #endif
            initTime = t;
        }

        void getL(OrdinalArray &lRowMapArg, OrdinalArray &lColIdxArg, ScalarArray &lValArg)
        {
            lRowMapArg = lRowMap;
            lColIdxArg = lColIdx;
            lValArg = lVal;
        }

        void getD(ScalarArray &diagElemsArg)
        {
            diagElemsArg = diagElems;
        }

        void transposeL()
        {
            auto ltRowMap_ = Kokkos::create_mirror_view(ltRowMap);
            auto ltColIdx_ = Kokkos::create_mirror_view(ltColIdx);
            auto ltVal_ = Kokkos::create_mirror_view(ltVal);

            //Count the elements in each row of Lt
            auto temp = OrdinalArrayHost("temp", nRows + 1);
            auto rowPtrs = OrdinalArrayHost("rowPtrs", nRows);
            for (Ordinal i = 0; i <= nRows; i++)
            {
                temp[i] = 0;
            }
            for (Ordinal i = 0; i < nRows; i++)
            {
                for (Ordinal k = lRowMap_[i]; k < lRowMap_[i+1]; k++)
                {
                    temp[lColIdx_[k]+1]++;

                }
            }
            //Perform an add scan to get the row map for
            //the transpose
            for (Ordinal i = 0; i <= nRows; i++)
            {
                ltRowMap_[i] = temp[i];
            }
            for (Ordinal i = 1; i <= nRows; i++)
            {
                ltRowMap_[i] += ltRowMap_[i-1];
            }
            //Set the row pointers to their initial places;
            for (Ordinal i = 0; i < nRows; i++)
            {
                rowPtrs[i] = ltRowMap_[i];
            }
            //Copy the data
            Kokkos::deep_copy(lVal_, lVal);
            for (Ordinal i = 0; i < nRows; i++)
            {
                for (Ordinal k = lRowMap_[i]; k < lRowMap_[i+1]; k++)
                {
                    Ordinal row = lColIdx_[k];
                    Scalar value = lVal_[k];
                    ltVal_[rowPtrs[row]] = value;
                    ltColIdx_[rowPtrs[row]] = i;
                    rowPtrs[row]++;
                    assert(rowPtrs[row] <= ltRowMap_[row + 1]);
                }
            }
            Kokkos::deep_copy(ltRowMap, ltRowMap_);
            Kokkos::deep_copy(ltColIdx, ltColIdx_);
            Kokkos::deep_copy(ltVal,ltVal_);
        }

        void symbolicILU()
        {
            using std::vector;
            using std::cout;
            using std::stable_sort;
            using std::sort;
            OrdinalArrayMirror ia = aRowMapHost;
            OrdinalArrayMirror ja = aColIdxHost;
            int n = nRows;
            int nzl = aRowMapHost[nRows] * (level + 2);
            int nzu = aRowMapHost[nRows] * (level + 2);
            Ordinal i;
            int levfill = level;
            vector<int> lnklst(n);
            vector<int> curlev(n);
            vector<int> levels(nzu);
            vector<int> iwork(n);
            vector<int> ial(nzl);
            vector<int> jal(nzl);
            vector<int> iau(nzu);
            vector<int> jau(nzu);

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
                    assert(knzl < nzl);
                    jal[knzl++] = next;
                    if (knzl >= nzl)
                    {
                        //(*nzl)++;
                        nzl += nRows;
                        jal.resize(nzl);
                    }
                    //jal.push_back(next);
                   // knzl++;
                    next = lnklst[next];
                }
                ial[i+1] = knzl;
                assert(next == i);

                while (next < n)
                {
                    assert(knzu < nzu);
                    levels[knzu] = curlev[next];
                    jau[knzu++] = next;
                    if (knzu >= nzu)
                    {
                        //(*nzu)++;
                        nzu += nRows;
                        jau.resize(nzu);
                        levels.resize(nzu);
                    }
                    //jau.push_back(next);
                    //knzu++;
                    next = lnklst[next];
                }
                iau[i+1] = knzu;
            }

            nzl = knzl;
            nzu = knzu;

            #ifdef FASTIC_DEBUG_OUTPUT
            std::cout << "knzl =" << knzl;
            std::cout << "knzu =" << knzu;
            std::cout << "ILU: nnz = "<< knzl + knzu << std::endl;
            std::cout << "Actual nnz for ILU: " << nzl + nzu << std::endl;
            #endif
            //Initialize the A matrix that is to be used in the computation
            aRowMap = OrdinalArray("aRowMap", nRows + 1);
            aColIdx = OrdinalArray("aColIdx", knzl + knzu);
            aRowIdx = OrdinalArray("aRowIds", knzl + knzu);
            aRowMap_ = Kokkos::create_mirror_view(aRowMap);
            aColIdx_ = Kokkos::create_mirror_view(aColIdx);
            aRowIdx_ = Kokkos::create_mirror_view(aRowIdx);

            Ordinal aRowPtr = 0;
            aRowMap_[0] = 0;
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

            //Now allocate memory for L and U.
            //
            lRowMap = OrdinalArray("lRowMap", nRows + 1);
            ltRowMap = OrdinalArray("ltRowMap", nRows + 1);
            countL();
            //Allocate memory and initialize pattern for L, U (transpose).
            lColIdx = OrdinalArray("lColIdx", lRowMap_[nRows]);
            ltColIdx = OrdinalArray("ltColIdx", lRowMap_[nRows]);
            #ifdef FASTIC_DEBUG_OUTPUT
            std::cout << "nnz L = " << lRowMap[nRows] << std::endl;
            #endif
        }

        void numericILU()
        {
            aVal = ScalarArray("aVal", aColIdx.extent(0));
            aVal_ = Kokkos::create_mirror_view(aVal);

            //Copy the host matrix into the initialized a;
            for (Ordinal i = 0; i < nRows; i++)
            {
                for(Ordinal k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    Ordinal col = aColIdx_[k];
                    #ifdef FASTIC_DEBUG_OUTPUT
                    std::cout << "col =" << col << std::endl;
                    #endif
                    for(Ordinal aHostPtr = aRowMapHost[i]; aHostPtr < aRowMapHost[i+1]; aHostPtr++)
                    {
                        if (col == aColIdxHost[aHostPtr])
                        {
                            aVal_[k] = aValHost[aHostPtr];
                            break;
                        }
                    }
                }
            }
            lVal = ScalarArray("lVal", lRowMap_[nRows]);
            ltVal = ScalarArray("ltVal", lRowMap_[nRows]);
            applyDiagonalScaling();
            applyManteuffelShift();
            Kokkos::deep_copy(aVal, aVal_);

            fillL();
        }

        void countL()
        {
            lRowMap_ = Kokkos::create_mirror_view(lRowMap);

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

        void fillL()
        {
            lVal_    = Kokkos::create_mirror_view(lVal);
            lColIdx_ = Kokkos::create_mirror_view(lColIdx);
            diagElems_ = Kokkos::create_mirror_view(diagElems);

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
                        }
                        lVal_[lPtr] = aVal_[k];
                        lColIdx_[lPtr] = aColIdx_[k];
                        lPtr++;
                    }
                }
            }
            Kokkos::deep_copy(diagElems, diagElems_);
            assert(lPtr == lRowMap_[nRows]);

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

                auto lGColIdx_ = Kokkos::create_mirror_view(lGColIdx);
                auto lGVal_ = Kokkos::create_mirror_view(lGVal);
                Kokkos::deep_copy(lGColIdx_, lGColIdx);
                Kokkos::deep_copy(lGVal_, lGVal);
                for (Ordinal i = 0; i < nRows; i++)
                {
                    #ifdef SHYLU_DEBUG
                    Ordinal check = lGPtr;
                    #endif
                    for (Ordinal k = lRowMap_[i]; k < lRowMap_[i+1]; k++)
                    {
                        //unused: Ordinal row = i;
                        Ordinal col = lColIdx_[k];

                        if (col == lGColIdx_[lGPtr])
                        {
                            lVal_[k] = lGVal_[lGPtr];
                            lGPtr++;
                        }
                    }
                    #ifdef SHYLU_DEBUG
                    rowLen = lGPtr - check;
                    assert(rowLen == lGRowMap[i+1] - lGRowMap[i]);
                    #endif
                }
            }
            Kokkos::deep_copy(lColIdx, lColIdx_);
            Kokkos::deep_copy(lVal, lVal_);
        }

        void applyDiagonalScaling()
        {
            const Real one = RTS::one();
            int anext = 0;
            //First fill Aj and extract the diagonal scaling factors
            //Use diag array to store scaling factors since
            //it gets set to the correct value by findFactorPattern anyway.
            auto diagFact_ = Kokkos::create_mirror_view(diagFact);
            for (int i = 0; i < nRows; i++)
            {
                for(int k = aRowMap_[i]; k < aRowMap_[i+1]; k++)
                {
                    aRowIdx_[anext++] = i;
                    if (aColIdx_[k] == i)
                    {
                        diagFact_[i] = one/RTS::sqrt(STS::abs(aVal_[k]));
                        #ifdef FASTIC_DEBUG_OUTPUT
                        std::cout << "diagFact["<<i<<"]="<<aVal[k]<<std::endl;
                        #endif
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
            //Scalar shift = 0.05;
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

        void applyDD(ScalarArray &x, ScalarArray &y)
        {
            ParScalFunctor<Ordinal, Scalar, Scalar, ExecSpace> parScal(x, y, diagElemsInv);
            ExecSpace().fence();
            Kokkos::parallel_for(RangePolicy(0, nRows), parScal);
            ExecSpace().fence();

        }
        void applyD(ScalarArray &x, ScalarArray &y)
        {
            ParScalFunctor<Ordinal, Scalar, Real, ExecSpace> parScal(x, y, diagFact);
            ExecSpace().fence();
            Kokkos::parallel_for(RangePolicy(0, nRows), parScal);
            ExecSpace().fence();

        }
        void applyLIC(ScalarArray &x, ScalarArray &y)
        {
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(xOld);
            Kokkos::parallel_for(RangePolicy(0, nRows), parInitZero);
            ExecSpace().fence();
#if 0
            JacobiIterFunctor<Ordinal, Scalar, ExecSpace> jacIter(nRows, lRowMap, lColIdx, lVal,
                                                                  x, y, xOld, diagElems);
#endif
            BlockJacobiIterFunctorL<Ordinal, Scalar, ExecSpace> jacIter(nRows, blkSz, lRowMap,
                                                                        lColIdx, lVal,
                                                                        x, y, xOld, diagElems);
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopy(xOld, y);
            Ordinal extent = nRows/blkSz;
            if(nRows%blkSz != 0)
            {
                extent++;
            }
            ExecSpace().fence();
            for (Ordinal i = 0; i < nTrisol; i++)
            {
                Kokkos::parallel_for(RangePolicy(0, extent), jacIter);
                ExecSpace().fence();
                Kokkos::parallel_for(RangePolicy(0, nRows), parCopy);
                ExecSpace().fence();
            }
            return;
        }

        void applyLT(ScalarArray &x, ScalarArray &y)
        {
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(xOld);
            Kokkos::parallel_for(RangePolicy(0, nRows), parInitZero);
            ExecSpace().fence();
#if 0
            JacobiIterFunctor<Ordinal, Scalar, ExecSpace> jacIter(nRows, ltRowMap, ltColIdx, ltVal,
                                                                  x, y, xOld, diagElems);
#endif
            BlockJacobiIterFunctorU<Ordinal, Scalar, ExecSpace> jacIter(nRows, blkSz, ltRowMap, ltColIdx,
                                                                        ltVal, x, y, xOld, diagElems);
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopy(xOld, y);
            Ordinal extent = nRows/blkSz;
            if (nRows%blkSz != 0)
            {
                extent++;
            }
            ExecSpace().fence();
            for (Ordinal i = 0; i < nTrisol; i++)
            {
                Kokkos::parallel_for(RangePolicy(0, extent), jacIter);
                ExecSpace().fence();
                Kokkos::parallel_for(RangePolicy(0, nRows), parCopy);
                ExecSpace().fence();
            }
            return;
        }


#if 0
        void applyLT(ScalarArray &x, ScalarArray &y)
        {
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero(nRows, xOld);
            ParInitZeroFunctor<Ordinal, Scalar, ExecSpace> parInitZero1(nRows, y);
            Kokkos::parallel_for(RangePolicy(0, nRows), parInitZero);
            ExecSpace().fence();
            //Transpose Jacobi implementation
            JacobiIterFunctorT<Ordinal, Scalar, ExecSpace> jacIterT(nRows, lRowMap, lColIdx, lVal, x, y, xOld, onesVector);
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopy(nRows, xOld, y);

            ExecSpace().fence();

            for (Ordinal i = 0; i < nTrisol; i++)
            {
                Kokkos::parallel_for(RangePolicy(0, nRows), parInitZero1);
                ExecSpace().fence();
                Kokkos::parallel_for(RangePolicy(0, nRows), jacIterT);
                ExecSpace().fence();
                Kokkos::parallel_for(RangePolicy(0, nRows), parCopy);
                ExecSpace().fence();
            }

            return;
        }
#endif

        void setValues(ScalarArray& aValsIn)
        {
          this->aValHost = Kokkos::create_mirror_view(aValsIn);
          Kokkos::deep_copy(this->aValHost, aValsIn);
          if(!initGuessPrec.is_null())
          {
            initGuessPrec->setValues(aValsIn);
          }
        }

        void compute()
        {
            const Scalar one = STS::one();
            Kokkos::Timer timer;
            if((level > 0) && (guessFlag != 0))
            {
                initGuessPrec->compute();
            }
            numericILU();
            FastICFunctor<Ordinal, Scalar, ExecSpace> icFunctor(aRowMap_[nRows], blkSzIC,
                    aRowMap, aColIdx, aRowIdx, aVal, lRowMap, lColIdx, lVal, diagElems, omega);
            ExecSpace().fence();
            Ordinal extent = aRowMap_[nRows]/blkSzIC;
            if (aRowMap_[nRows]%blkSzIC != 0)
            {
                extent++;
            }
            //Ordinal extent = aRowMap[nRows];
            for (int i = 0; i < nFact; i++)
            {
                Kokkos::parallel_for(RangePolicy(0, extent), icFunctor);
            }
            ExecSpace().fence();

            double t = timer.seconds();
            computeTime = t;

            Kokkos::deep_copy(diagElems_, diagElems);
            diagElemsInv_ = Kokkos::create_mirror_view(diagElemsInv);
            for (int i = 0; i < nRows; i++)
            {
                diagElemsInv_[i] = one/diagElems_[i];
            }
            Kokkos::deep_copy(diagElemsInv, diagElemsInv_);
            transposeL();

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
                // setup Lt solve
                khLt.create_sptrsv_handle(algo, nRows, false);
                #if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
                KokkosSparse::Experimental::sptrsv_symbolic(&khLt, ltRowMap, ltColIdx, ltVal);
                #else
                KokkosSparse::Experimental::sptrsv_symbolic(&khLt, ltRowMap, ltColIdx);
                #endif
            }
        }

        void apply(ScalarArray &x, ScalarArray &y)
        {
            Kokkos::Timer timer;
            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopyFunctor(xTemp, x);
            ExecSpace().fence();
            Kokkos::parallel_for(RangePolicy(0, nRows), parCopyFunctor);
            ExecSpace().fence();

            applyD(x, xTemp);
            if (standard_sptrsv) {
                // solve with L
                KokkosSparse::Experimental::sptrsv_solve(&khL, lRowMap, lColIdx, lVal, xTemp, y);
                // solve with Lt
                KokkosSparse::Experimental::sptrsv_solve(&khLt, ltRowMap, ltColIdx, ltVal, y, xTemp);
            } else {
                applyLIC(xTemp, y);
                //applyDD(y, xTemp);
                applyLT(y, xTemp);
            }
            applyD(xTemp, y);

//            ParCopyFunctor<Ordinal, Scalar, ExecSpace> parCopyFunctor2(nRows, y, xTemp);
//            ExecSpace().fence();
//            Kokkos::parallel_for(RangePolicy(0, nRows), parCopyFunctor2);
//            ExecSpace().fence();

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
            return (standard_sptrsv ? "Standard" : "Fast");
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
            Scalar sum = STS::zero();
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
        friend class ParScalFunctor<Ordinal, Scalar, Real, ExecSpace>;
        friend class MemoryPrimeFunctorN<Ordinal, Scalar, ExecSpace>;
        friend class MemoryPrimeFunctorNnzCoo<Ordinal, Scalar, ExecSpace>;
        friend class MemoryPrimeFunctorNnzCsr<Ordinal, Scalar, ExecSpace>;
};


#endif
