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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_STATE_FUNC_MODEL_EVALUATOR_BASE_HPP
#define THYRA_STATE_FUNC_MODEL_EVALUATOR_BASE_HPP

#include "Thyra_ModelEvaluatorDefaultBase.hpp"


namespace Thyra {


/** \brief This base class defines default function implementations
 * appropritate for a set of nonlinear state functions of the form
 * <tt>x -> f(x)</tt>.
 *
 * The minimum that a subclass must do is to define implementations for
 * <tt>get_x_space()</tt>, <tt>get_f_space()</tt>, <tt>createInArgs()</tt>,
 * <tt>createOutArgsImpl</tt>, and <tt>evalModelImpl()</tt>.
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class StateFuncModelEvaluatorBase : virtual public ModelEvaluatorDefaultBase<Scalar> {
public:

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief Throws exception. */
  RCP<const VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief Throws exception. */
  RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief Throws exception. */
  RCP<const VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief Throws exception. */
  Teuchos::ArrayView<const std::string> get_g_names(int j) const;
  /** \brief Returns this->createInArgs(). */
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief Returns this->createInArgs(). */
  ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief Returns this->createInArgs(). */
  ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;
  /** \brief Throws exception. */
  RCP<LinearOpBase<Scalar> > create_W_op() const;
  /** \brief Returns null. */
  RCP<PreconditionerBase<Scalar> > create_W_prec() const;
  /** \brief Returns null . */
  RCP<const LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  /** \brief Ignores input and does nothing. */
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
StateFuncModelEvaluatorBase<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"ModelEvaluator<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::get_p_space(l): "
    "Error, this function was not overridden in *this = \'"<<this->description()<<"\'!"
    );
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}


template<class Scalar>
RCP<const Teuchos::Array<std::string> >
StateFuncModelEvaluatorBase<Scalar>::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"ModelEvaluator<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::get_p_names(l): "
    "Error, this function was not overridden in *this = \'"<<this->description()<<"\'!"
    );
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
StateFuncModelEvaluatorBase<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"ModelEvaluator<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::get_g_space(j): "
    " Error, this function was not overridden in \'"
    <<this->description()<<"\'!"
    );
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}


template<class Scalar>
Teuchos::ArrayView<const std::string>
StateFuncModelEvaluatorBase<Scalar>::get_g_names(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"ModelEvaluator<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::get_g_names(j): "
    "Error, this function was not overridden in *this = \'"<<this->description()<<"\'!"
    );
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::ArrayView<const std::string>(Teuchos::null));
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
StateFuncModelEvaluatorBase<Scalar>::getNominalValues() const
{ return this->createInArgs(); }


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
StateFuncModelEvaluatorBase<Scalar>::getLowerBounds() const
{ return this->createInArgs(); }


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
StateFuncModelEvaluatorBase<Scalar>::getUpperBounds() const
{ return this->createInArgs(); }


template<class Scalar>
RCP<LinearOpBase<Scalar> >
StateFuncModelEvaluatorBase<Scalar>::create_W_op() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"Error, if \'W_op\' is supported by the ModelEvaluator subclass then"
    " this function create_W_op() must be overridden by the subclass "
    <<this->description()<<" to return a non-null object!"
    );
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}


template<class Scalar>
RCP<PreconditionerBase<Scalar> >
StateFuncModelEvaluatorBase<Scalar>::create_W_prec() const
{
  return Teuchos::null;
}


template<class Scalar>
RCP<const LinearOpWithSolveFactoryBase<Scalar> >
StateFuncModelEvaluatorBase<Scalar>::get_W_factory() const
{
  return Teuchos::null;
}


template<class Scalar>
void StateFuncModelEvaluatorBase<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
  const bool wasSolved
  )
{
  // This final point is just ignored by default!
}


} // namespace Thyra


#endif // THYRA_STATE_FUNC_MODEL_EVALUATOR_BASE_HPP
