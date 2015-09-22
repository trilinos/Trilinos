// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_MAX_INDEX_LESS_THAN_BOUND_HPP
#define RTOPPACK_ROP_MAX_INDEX_LESS_THAN_BOUND_HPP

#include "RTOpPack_ROpMaxIndex.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"


namespace RTOpPack {


/** \brief . */
template<class Scalar>
class ROpMaxIndexLessThanBoundEleWiseReductionOp {
public:
  /** \brief . */
  ROpMaxIndexLessThanBoundEleWiseReductionOp(
    const Scalar &bound = ScalarTraits<Scalar>::zero()
    )
    :bound_(bound)
    {}
  /** \brief . */
  void operator()(const index_type i, const Scalar &v0,
    ScalarIndex<Scalar> &reduct) const
    {
      if(
        v0 <  bound_
        &&
        (
          v0 > reduct.scalar
          ||
          ( v0 == reduct.scalar && i < reduct.index )
          )
        )
      {
        reduct = ScalarIndex<Scalar>(v0, i);
      }
    }
private:
  Scalar bound_;
};


/** \brief Returns the maximum element less than some bound along with
 * its index: <tt>result.scalar = x(k)</tt> and <tt>result.index =
 * k</tt> such that <tt>x(k) >= x(i)</tt> for all <tt>i</tt> where
 * <tt>x(i) < bound</tt> and <tt>k</tt> is the minimum index to break
 * ties.
 *
 * If no element is less than <tt>bound</tt> then <tt>results.index <
 * 0</tt>.
 */
template<class Scalar>
class ROpMaxIndexLessThanBound
  : public ROp_1_CoordVariantScalarReduction<
      Scalar,
      ScalarIndex<Scalar>,
      ROpMaxIndexLessThanBoundEleWiseReductionOp<Scalar>,
      ROpMaxIndexReductObjReductionOp<Scalar> >
{
public:
  /** \brief . */
  ROpMaxIndexLessThanBound(
    const Scalar &bound_in = Teuchos::ScalarTraits<Scalar>::zero()
    )
    {
      this->setOpNameBase("ROpMaxIndexLessThanBound");
      bound(bound_in);
      this->initReductObjValue(
        ScalarIndex<Scalar>(-ScalarTraits<Scalar>::rmax(), -1));
    }
  /** \brief . */
  void bound(const Scalar& bound_in)
    { 
      this->setEleWiseReduction(
        ROpMaxIndexLessThanBoundEleWiseReductionOp<Scalar>(bound_in)
        );
    }
  /** \brief . */
  ScalarIndex<Scalar> operator()(const ReductTarget& reduct_obj ) const
    { return this->getRawVal(reduct_obj); }
};


} // namespace RTOpPack


#endif // RTOPPACK_ROP_MAX_INDEX_LESS_THAN_BOUND_HPP
