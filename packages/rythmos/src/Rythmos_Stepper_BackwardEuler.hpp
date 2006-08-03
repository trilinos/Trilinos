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

#ifndef Rythmos_STEPPER_BACKWARD_EULER_H
#define Rythmos_STEPPER_BACKWARD_EULER_H

#include "Rythmos_Stepper.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_SingleResidSSDAEModelEvaluator.hpp"

namespace Rythmos {

/** \brief . */
template<class Scalar>
class BackwardEulerStepper : virtual public Stepper<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** \brief . */
    BackwardEulerStepper(
      const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> >  &model
      ,const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> >  &solver
      );

    /** \brief . */
    void setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model);

    /** \brief . */
    void setSolver(const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &solver);

    /** \brief . */
    Scalar TakeStep(Scalar dt);
   
    /** \brief . */
    Scalar TakeStep();

    /** \brief . */
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_solution() const;

    /** \brief . */
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_residual() const;

    /// Redefined from describable
    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(
      Teuchos::FancyOStream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const;

    /// Redefined from InterpolationBuffer 
    /// Add points to buffer
    bool SetPoints(
      const std::vector<Scalar>& time_list
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_list
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_list);
    
    /// Get values from buffer
    bool GetPoints(
      const std::vector<Scalar>& time_list
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_list
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_list
      ,std::vector<ScalarMag>* accuracy_list) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBuffer<Scalar> & IB);

    /// Get interpolation nodes
    bool GetNodes(std::vector<Scalar>* time_list) const;

    /// Remove interpolation nodes
    bool RemoveNodes(std::vector<Scalar>& time_list) const;

    /// Get order of interpolation
    int GetOrder() const;

  private:

    Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > model_;
    Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > solver_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > scaled_x_old_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > f_;
    Scalar t_;
    Scalar t_old_;

    Teuchos::RefCountPtr<Thyra::SingleResidSSDAEModelEvaluator<Scalar> >  neModel_;

};

// ////////////////////////////
// Defintions

template<class Scalar>
BackwardEulerStepper<Scalar>::BackwardEulerStepper(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
  ,const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  setModel(model);
  setSolver(solver);
}

template<class Scalar>
void BackwardEulerStepper<Scalar>::setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model_ = model;
  t_ = ST::zero();
  x_ = model_->getNominalValues().get_x()->clone_v();
  f_ = Thyra::createMember(model_->get_f_space());

  scaled_x_old_ = x_->clone_v();
  
}

template<class Scalar>
void BackwardEulerStepper<Scalar>::setSolver(const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &solver)
{
  solver_ = solver;
}

template<class Scalar>
Scalar BackwardEulerStepper<Scalar>::TakeStep()
{
  // print something out about this method not supporting automatic variable step-size
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return(-ST::one());
}

template<class Scalar>
Scalar BackwardEulerStepper<Scalar>::TakeStep(Scalar dt)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  //
  // Setup the nonlinear equations:
  //
  //   f( (1/dt)* x + (-1/dt)*x_old), x, t ) = 0
  //
  V_StV( &*scaled_x_old_, Scalar(-ST::one()/dt), *x_ );
  t_old_ = t_;
  if(!neModel_.get())
    neModel_ = Teuchos::rcp(new Thyra::SingleResidSSDAEModelEvaluator<Scalar>());
  neModel_->initialize(model_,Scalar(ST::one()/dt),scaled_x_old_,ST::one(),Teuchos::null,t_old_+dt,Teuchos::null);
  if(solver_->getModel().get()!=neModel_.get())
    solver_->setModel(neModel_);
  //
  // Solve the implicit nonlinear system to a tolerance of ???
  //
  Thyra::assign(&*x_,ST::zero());
  solver_->solve(&*x_); // Note that x in input is x_old and on output is the solved x!
  //
  // Update the step
  //
  t_ += dt;

  return(dt);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > BackwardEulerStepper<Scalar>::get_solution() const
{
  return(x_);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > BackwardEulerStepper<Scalar>::get_residual() const
{
  return(f_);
}

template<class Scalar>
std::string BackwardEulerStepper<Scalar>::description() const
{
  std::string name = "Rythmos::BackwardEulerStepper";
  return(name);
}

template<class Scalar>
void BackwardEulerStepper<Scalar>::describe(
      Teuchos::FancyOStream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const
{
  if (verbLevel == Teuchos::VERB_EXTREME)
  {
    out << description() << "::describe:";
    out << "\nmodel_ = " << std::endl;
    model_->describe(out,verbLevel);
    out << "\nsolver_ = " << std::endl;
    solver_->describe(out,verbLevel);
    out << "\nx_ = " << std::endl;
    x_->describe(out,verbLevel);
    out << "\nscaled_x_old_ = " << std::endl;
    scaled_x_old_->describe(out,verbLevel);
    out << "\nf_ = " << std::endl;
    f_->describe(out,verbLevel);
    out << "\nt_ = " << t_;
    out << "\nt_old_ = " << t_old_;
    out << std::endl;
//    out << "neModel_ = " << std::endl;
//    out << neModel_->describe(out,verbLevel) << std::endl;
  }
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::SetPoints(
    const std::vector<Scalar>& time_list
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_list
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_list)
{
  return(false);
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::GetPoints(
    const std::vector<Scalar>& time_list
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_list
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_list
    ,std::vector<ScalarMag>* accuracy_list) const
{
  return(false);
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::SetRange(
    const Scalar& time_lower
    ,const Scalar& time_upper
    ,const InterpolationBuffer<Scalar>& IB)
{
  return(false);
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::GetNodes(std::vector<Scalar>* time_list) const
{
  return(false);
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::RemoveNodes(std::vector<Scalar>& time_list) const
{
  return(false);
}

template<class Scalar>
int BackwardEulerStepper<Scalar>::GetOrder() const
{
  return(1);
}


} // namespace Rythmos

#endif //Rythmos_STEPPER_BACKWARD_EULER_H
