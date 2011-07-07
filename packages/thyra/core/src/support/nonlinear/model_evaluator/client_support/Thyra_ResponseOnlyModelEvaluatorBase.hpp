// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
