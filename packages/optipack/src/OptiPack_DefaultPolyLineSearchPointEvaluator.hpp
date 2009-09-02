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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
