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

#ifndef Rythmos_FORWARDEULER_STEPPER_H
#define Rythmos_FORWARDEULER_STEPPER_H

#include "Rythmos_Stepper.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"

namespace Rythmos {

/** \brief . */
template<class Scalar>
class ForwardEulerStepper : virtual public Stepper<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    
    /** \brief . */
    ForwardEulerStepper();
    ForwardEulerStepper(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model);
    
    /** \brief . */
    ~ForwardEulerStepper();

    /** \brief . */
    Scalar TakeStep(Scalar dt);
   
    /** \brief . */
    Scalar TakeStep();

    /** \brief . */
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_solution() const;

    /** \brief . */
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_residual() const;
    
    /** \brief . */
    std::string description() const;

    /** \brief . */
    std::ostream& describe(
      std::ostream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ,const std::string          leadingIndent
      ,const std::string          indentSpacer
      ) const;
    
    /// Redefined from InterpolationBufferBase 
    /// Add points to buffer
    bool SetPoints(
      const std::vector<Scalar>& time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
      ,const std::vector<ScalarMag>& accuracy_vec
      );
    
    /// Get values from buffer
    bool GetPoints(
      const std::vector<Scalar>& time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
      ,std::vector<ScalarMag>* accuracy_vec) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBufferBase<Scalar> & IB);

    /// Get interpolation nodes
    bool GetNodes(std::vector<Scalar>* time_vec) const;

    /// Remove interpolation nodes
    bool RemoveNodes(std::vector<Scalar>& time_vec);

    /// Get order of interpolation
    int GetOrder() const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();

  private:

    Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > model_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > solution_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > residual_vector_;
    Scalar t_;

};

template<class Scalar>
ForwardEulerStepper<Scalar>::ForwardEulerStepper(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model_ = model;
  t_ = ST::zero();
  solution_vector_ = model_->getNominalValues().get_x()->clone_v();
  residual_vector_ = Thyra::createMember(model_->get_f_space());
}

template<class Scalar>
ForwardEulerStepper<Scalar>::ForwardEulerStepper()
{
}

template<class Scalar>
ForwardEulerStepper<Scalar>::~ForwardEulerStepper()
{
}

template<class Scalar>
Scalar ForwardEulerStepper<Scalar>::TakeStep()
{
  // print something out about this method not supporting automatic variable step-size
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return(-ST::one());
}

template<class Scalar>
Scalar ForwardEulerStepper<Scalar>::TakeStep(Scalar dt)
{
/*
  Thyra::ModelEvaluatorBase::InArgs<Scalar>   inArgs  = model_->createInArgs();
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>  outArgs = model_->createOutArgs();

  inArgs.set_x(solution_vector_);
  inArgs.set_t(t_+dt);

  outArgs.set_f(residual_vector_);

  model_->evalModel(inArgs,outArgs);
*/
  Thyra::eval_f<Scalar>(*model_,*solution_vector_,t_+dt,&*residual_vector_);

  // solution_vector = solution_vector + dt*residual_vector
  Thyra::Vp_StV(&*solution_vector_,dt,*residual_vector_); 
  t_ += dt;

  return(dt);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ForwardEulerStepper<Scalar>::get_solution() const
{
  return(solution_vector_);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ForwardEulerStepper<Scalar>::get_residual() const
{
  return(residual_vector_);
}

template<class Scalar>
std::string ForwardEulerStepper<Scalar>::description() const
{
  std::string name = "Rythmos::ForwardEulerStepper";
  return(name);
}

template<class Scalar>
std::ostream& ForwardEulerStepper<Scalar>::describe(
      std::ostream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ,const std::string          leadingIndent
      ,const std::string          indentSpacer
      ) const
{
  if (verbLevel == Teuchos::VERB_EXTREME)
  {
    out << description() << "::describe" << std::endl;
    out << "model_ = " << std::endl;
    out << model_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "solution_vector_ = " << std::endl;
    out << solution_vector_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "residual_vector_ = " << std::endl;
    out << residual_vector_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  }
  return(out);
}

template<class Scalar>
bool ForwardEulerStepper<Scalar>::SetPoints(
    const std::vector<Scalar>& time_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
    ,const std::vector<ScalarMag> & accuracy_vec 
    )
{
  return(false);
}

template<class Scalar>
bool ForwardEulerStepper<Scalar>::GetPoints(
    const std::vector<Scalar>& time_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
    ,std::vector<ScalarMag>* accuracy_vec) const
{
  return(false);
}

template<class Scalar>
bool ForwardEulerStepper<Scalar>::SetRange(
    const Scalar& time_lower
    ,const Scalar& time_upper
    ,const InterpolationBufferBase<Scalar>& IB)
{
  return(false);
}

template<class Scalar>
bool ForwardEulerStepper<Scalar>::GetNodes(std::vector<Scalar>* time_vec) const
{
  return(false);
}

template<class Scalar>
bool ForwardEulerStepper<Scalar>::RemoveNodes(std::vector<Scalar>& time_vec) 
{
  return(false);
}

template<class Scalar>
int ForwardEulerStepper<Scalar>::GetOrder() const
{
  return(1);
}

template <class Scalar>
void ForwardEulerStepper<Scalar>::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> ForwardEulerStepper<Scalar>::getParameterList()
{
  return(Teuchos::null);
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> ForwardEulerStepper<Scalar>::unsetParameterList()
{
  return(Teuchos::null);
}


} // namespace Rythmos

#endif //Rythmos_FORWARDEULER_STEPPER_H
