/*
// @HEADER
// ***********************************************************************
// 
//    OptiPack: Collection of simple Thyra-based Optimization ANAs
//                 Copyright (2009) Sandia Corporation
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
*/

#ifndef OPTIPACK_DEFAULT_POLY_LINE_SEARCH_POINT_EVALUATOR_HPP
#define OPTIPACK_DEFAULT_POLY_LINE_SEARCH_POINT_EVALUATOR_HPP


#include "OptiPack_LineSearchPointEvaluatorBase.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace OptiPack {


/** \brief Default line search point evaluator using a polynomial linear
 * combination of vectors.
 *
 * This object computes:
 
 \verbatim

   p = sum( alpha^i * vec[i], i = 0...n-1 )

 \endverbatim

 * This allows, for instance, a curvy-linear line search algorithm.
 */
template<typename Scalar>
class DefaultPolyLineSearchPointEvaluator : public LineSearchPointEvaluatorBase<Scalar>
{
public:

  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors/intializers/accessors. */
  //@{

  /** \brief . */
  DefaultPolyLineSearchPointEvaluator();
  
  /** \brief . */
  void initialize(const ArrayView<const RCP<const Thyra::VectorBase<Scalar> > > &vecs);

  //@}

  /** \name Overridden from LineSearchPointEvaluatorBase. */
  //@{

  /** \brief . */
  virtual void computePoint( const ScalarMag &alpha,
    const Ptr<Thyra::VectorBase<Scalar> > &p
    ) const;

  //@}

private:

  Array<RCP<const Thyra::VectorBase<Scalar> > > vecs_;

};


/** \brief Nonmember constructor.
 *
 * \relates DefaultPolyLineSearchPointEvaluator
 */
template<typename Scalar>
const RCP<DefaultPolyLineSearchPointEvaluator<Scalar> >
defaultPolyLineSearchPointEvaluator()
{
  return Teuchos::rcp(new DefaultPolyLineSearchPointEvaluator<Scalar>);
}


//
// Implementations
//


// Constructors/intializers/accessors


template<typename Scalar>
DefaultPolyLineSearchPointEvaluator<Scalar>::DefaultPolyLineSearchPointEvaluator()
{}
  

template<typename Scalar>
void DefaultPolyLineSearchPointEvaluator<Scalar>::initialize(
  const ArrayView<const RCP<const Thyra::VectorBase<Scalar> > > &vecs
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT(vecs.size());
#endif
  vecs_ = vecs;
}


// Overridden from LineSearchPointEvaluatorBase


template<typename Scalar>
void DefaultPolyLineSearchPointEvaluator<Scalar>::computePoint( const ScalarMag &alpha,
  const Ptr<Thyra::VectorBase<Scalar> > &p
  ) const
{
  typedef ScalarTraits<Scalar> ST;
  using Teuchos::as;
  using Thyra::V_V;
  using Thyra::Vp_StV;
  V_V( p, *vecs_[0] );
  if (alpha != ST::zero()) {
    ScalarMag alpha_i = alpha;
    const int n = vecs_.size();
    for (int i = 1; i < n; ++i, alpha_i *= alpha) {
      Vp_StV(p, alpha_i, *vecs_[i]); 
    }
  }
}


} // namespace OptiPack


#endif // OPTIPACK_DEFAULT_POLY_LINE_SEARCH_POINT_EVALUATOR_HPP
