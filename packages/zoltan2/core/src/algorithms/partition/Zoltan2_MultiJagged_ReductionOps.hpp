// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_MultiJagged_ReductionOps.hpp
  \brief Contains Teuchos redcution operators for the Multi-jagged algorthm.
 */

#ifndef _ZOLTAN2_MultiJagged_ReductionOps_HPP_
#define _ZOLTAN2_MultiJagged_ReductionOps_HPP_

#include <Teuchos_ReductionOp.hpp>

namespace Teuchos{
template <typename Ordinal, typename T>
class MultiJaggedCombinedReductionOp  : public ValueTypeReductionOp<Ordinal,T>
{
private:
    Ordinal numSum_0, numMin_1, numMin_2;
    std::vector <Ordinal> *partVector;
    Ordinal vectorBegin;
    Ordinal k;
    int reductionType;

public:
    /*! \brief Default Constructor
     */
    MultiJaggedCombinedReductionOp ():numSum_0(0), numMin_1(0),
        numMin_2(0), k(0), partVector(NULL), vectorBegin(0), reductionType(0){}

    /*! \brief Constructor
     *   \param nsum  the count of how many sums will be computed at the
     *             start of the list.
     *   \param nmin  following the sums, this many minimums will be computed.
     *   \param nmax  following the minimums, this many maximums will be computed.
     */
    MultiJaggedCombinedReductionOp (Ordinal nsum, Ordinal nmin1, Ordinal nmin2,
     Ordinal k_):
        numSum_0(nsum), numMin_1(nmin1), numMin_2(nmin2), partVector(NULL),
        vectorBegin(0), k(k_), reductionType(0){}


    MultiJaggedCombinedReductionOp (std::vector <Ordinal> *pVector, Ordinal vBegin,
     Ordinal k_):
        numSum_0(0), numMin_1(0), numMin_2(0), partVector(pVector),
        vectorBegin(vBegin), k(k_), reductionType(1){}


    /*! \brief Implement Teuchos::ValueTypeReductionOp interface
     */
    void reduce( const Ordinal /* count */, const T inBuffer[], T inoutBuffer[]) const
    {
        if (reductionType == 0){
            Ordinal next=0;
            for(Ordinal ii = 0; ii < k ; ++ii){
                for (Ordinal i=0; i < numSum_0; i++, next++)
                    inoutBuffer[next] += inBuffer[next];

                for (Ordinal i=0; i < numMin_1; i++, next++)
                    if (inoutBuffer[next] < inBuffer[next])
                        inoutBuffer[next] = inBuffer[next];

                for (Ordinal i=0; i < numMin_2; i++, next++)
                    if (inoutBuffer[next] > inBuffer[next])
                        inoutBuffer[next] = inBuffer[next];
            }
        }
        else {
            Ordinal next=0;
            for(Ordinal ii = 0; ii < k ; ++ii){
                Ordinal partPartition = (*partVector)[ii + vectorBegin];
                Ordinal tnumSum_ = 2 * partPartition - 1;
                Ordinal tnumMin_1 = partPartition - 1;
                Ordinal tnumMin_2 = tnumMin_1 ;
                for (Ordinal i=0; i < tnumSum_; i++, next++)
                    inoutBuffer[next] += inBuffer[next];

                for (Ordinal i=0; i < tnumMin_1; i++, next++)
                    if (inoutBuffer[next] < inBuffer[next])
                        inoutBuffer[next] = inBuffer[next];

                for (Ordinal i=0; i < tnumMin_2; i++, next++)
                    if (inoutBuffer[next] > inBuffer[next])
                        inoutBuffer[next] = inBuffer[next];
            }
        }
    }
};


template <typename Ordinal, typename T>
class MultiJaggedCombinedMinMaxTotalReductionOp :
     public ValueTypeReductionOp<Ordinal,T>
{
private:
    Ordinal numMin, numMax, numTotal;

public:
    /*! \brief Default Constructor
     */
    MultiJaggedCombinedMinMaxTotalReductionOp ():numMin(0), numMax(0), numTotal(0)
    {}

    /*! \brief Constructor
     *   \param nsum  the count of how many sums will be computed at the
     *             start of the list.
     *   \param nmin  following the sums,  #minimums will be computed.
     *   \param nmax  following the minimums, #maximums will be computed.
     */
    MultiJaggedCombinedMinMaxTotalReductionOp (Ordinal nmin, Ordinal nmax, Ordinal
     nTotal): numMin(nmin), numMax(nmax), numTotal(nTotal){}

    /*! \brief Implement Teuchos::ValueTypeReductionOp interface
     */
    void reduce( const Ordinal /* count */, const T inBuffer[], T inoutBuffer[]) const
    {
        Ordinal next=0;

        for (Ordinal i=0; i < numMin; i++, next++)
            if (inoutBuffer[next] > inBuffer[next])
                inoutBuffer[next] = inBuffer[next];

        for (Ordinal i=0; i < numMax; i++, next++)
            if (inoutBuffer[next] < inBuffer[next])
                inoutBuffer[next] = inBuffer[next];


        for (Ordinal i=0; i < numTotal; i++, next++)
            inoutBuffer[next] += inBuffer[next];
    }
};
} // namespace Teuchos

#endif //_ZOLTAN2_MultiJagged_ReductionOps_HPP_
