//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_IMPLICITBDF_STEPPER_ERR_WT_VEC_CALC_H
#define Rythmos_IMPLICITBDF_STEPPER_ERR_WT_VEC_CALC_H

#include "Rythmos_ErrWtVecCalcBase.hpp"

namespace Rythmos {

template<class Scalar>
class ImplicitBDFStepperErrWtVecCalc
  : virtual public ErrWtVecCalcBase<Scalar>
{
  public:

    /** \brief . */
    void errWtVecSet(
         Thyra::VectorBase<Scalar>* weight
         ,const Thyra::VectorBase<Scalar>& vector
         ,Scalar relTol
         ,Scalar absTol
         ) const;

    /** \name Overridden from ParameterListAcceptor */
    //@{
    /** \brief . */
    void setParameterList(RCP<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    RCP<Teuchos::ParameterList> getNonconstParameterList();

    /** \brief . */
    RCP<Teuchos::ParameterList> unsetParameterList();

    /** \brief . */
    RCP<const Teuchos::ParameterList> getValidParameters() const;

    //@}

  private:
    RCP<Teuchos::ParameterList> paramList_;
};


template<class Scalar>
void ImplicitBDFStepperErrWtVecCalc<Scalar>::errWtVecSet(
     Thyra::VectorBase<Scalar>* weight
     ,const Thyra::VectorBase<Scalar>& vector
     ,Scalar relTol
     ,Scalar absTol
     ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT(weight==NULL);
  TEST_FOR_EXCEPTION(
      ( ( relTol == ST::zero() ) && ( absTol == ST::zero() ) ), std::logic_error,
      "Error, relTol and absTol cannot both be zero!\n"
      );
  Thyra::VectorBase<Scalar> &w = *weight;
  Thyra::abs(ptr(&w),vector);
  Vt_S(ptr(&w),relTol);
  Vp_S(ptr(&w),absTol);
  reciprocal(ptr(&w),w);
  Vt_StV(ptr(&w),ST::one(),w); // We square w because of how weighted norm_2 is computed.
  // divide by N to get RMS norm
  int N = vector.space()->dim();
  Vt_S(ptr(&w),Scalar(1.0/N));
  // Now you can compute WRMS norm as:
  // Scalar WRMSnorm = norm_2(w,y); // WRMS norm of y with respect to weights w.

  using Teuchos::as;
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"errWtVecSet");

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
    *out << "weight = " << std::endl;
    weight->describe(*out,verbLevel);
  }
}

template<class Scalar>
void ImplicitBDFStepperErrWtVecCalc<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(paramList == Teuchos::null);
  paramList->validateParameters(*this->getValidParameters(),0);
  paramList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
ImplicitBDFStepperErrWtVecCalc<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = paramList_;
  paramList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
ImplicitBDFStepperErrWtVecCalc<Scalar>::getNonconstParameterList()
{
  return(paramList_);
}

template<class Scalar>
RCP<const Teuchos::ParameterList>
ImplicitBDFStepperErrWtVecCalc<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList>
      pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return (validPL);
}

} // namespace Rythmos

#endif // Rythmos_IMPLICITBDF_STEPPER_ERR_WT_VEC_CALC_H

