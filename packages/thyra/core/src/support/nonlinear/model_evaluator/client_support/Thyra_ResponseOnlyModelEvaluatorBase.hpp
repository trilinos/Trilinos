// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_RESPONSE_ONLY_MODEL_EVALUATOR_BASE_HPP
#define THYRA_RESPONSE_ONLY_MODEL_EVALUATOR_BASE_HPP


#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {


/** \brief This base class defines default function implementations
 * appropritate for a response-only model evaluator <tt>(p) -> g(j)</tt>, for
 * <tt>j=0...Ng-1</tt>.
 *
 * The minimum that a subclass must to is to define implementations for
 * <tt>get_p_space()</tt>, <tt>get_g_space()</tt>, <tt>createInArgs()</tt>,
 * <tt>createOutArgsImpl</tt>, and <tt>evalModelImpl()</tt>.
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class ResponseOnlyModelEvaluatorBase : virtual public ModelEvaluatorDefaultBase<Scalar> {
public:

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief Throws exception. */
  RCP<const VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief Returns null. */
  RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief Throws exception. */
  RCP<const VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief Returns this->createInArgs(). */
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief Returns this->createInArgs(). */
  ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief Returns this->createInArgs(). */
  ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;
  /** \brief Thorws exception. */
  RCP<LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief Thorws exception. */
  RCP<LinearOpBase<Scalar> > create_W_op() const;
  /** \brief Thorws exception. */
  RCP<const LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  /** \brief Does nothing and ignores input. */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
    const bool wasSolved
    );

  //@}
  
};


// /////////////////////////////////
// Implementations


// Public functions overridden from ModelEvaulator


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::get_x_space() const
{
  return Teuchos::null;
}


template<class Scalar>
RCP<const Teuchos::Array<std::string> >
ResponseOnlyModelEvaluatorBase<Scalar>::get_p_names(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, this->Np() );
#endif
  return Teuchos::null;
}

template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::get_f_space() const
{
  return Teuchos::null;
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ResponseOnlyModelEvaluatorBase<Scalar>::getNominalValues() const
{ return this->createInArgs(); }


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ResponseOnlyModelEvaluatorBase<Scalar>::getLowerBounds() const
{ return this->createInArgs(); }


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ResponseOnlyModelEvaluatorBase<Scalar>::getUpperBounds() const
{ return this->createInArgs(); }


template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::create_W() const
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"Error, if \'W\' is supported by the ModelEvaluator subclass then"
    " this function create_W() must be overridden by the subclass to return"
    " a non-null object!"
    );
  return Teuchos::null; // Should never be called!
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::create_W_op() const
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"Error, if \'W\' is supported by the ModelEvaluator subclass then"
    " this function create_W_op() must be overridden by the subclass "
    <<this->description()<<" to return a non-null object!"
    );
  return Teuchos::null; // Should never be called!
}


template<class Scalar>
RCP<const LinearOpWithSolveFactoryBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::get_W_factory() const
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"Error, if \'W\' is supported by the ModelEvaluator subclass then"
    " this function get_W_factory() must be overridden by the subclass "
    <<this->description()<<" to return a non-null object!"
    );
  return Teuchos::null; // Should never be called!
}


template<class Scalar>
void ResponseOnlyModelEvaluatorBase<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
  const bool wasSolved
  )
{
  // This final point is just ignored by default!
}


} // namespace Thyra


#endif // THYRA_RESPONSE_ONLY_MODEL_EVALUATOR_BASE_HPP
