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

#ifndef Rythmos_BACKWARD_EULER_STEPPER_H
#define Rythmos_BACKWARD_EULER_STEPPER_H

#include "Rythmos_Stepper.hpp"
#include "Rythmos_DataStore.hpp"
#include "Rythmos_LinearInterpolator.hpp"
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

    /// Redefined from describable
    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(
      Teuchos::FancyOStream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const;

    /// Redefined from InterpolationBufferBase 
    /// Add points to buffer
    /// This will take the last one or two points in the list and set up to integrate from here.
    bool SetPoints(
      const std::vector<Scalar>& time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
      ,const std::vector<ScalarMag> & accuracy_vec 
      );
    
    /// Get values from buffer
    /// This will interpolate points if t_old_ != t_
    bool GetPoints(
      const std::vector<Scalar>& time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
      ,std::vector<ScalarMag>* accuracy_vec) const;

    /// Fill data in from another interpolation buffer
    /// This will do the same as SetPoints
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBufferBase<Scalar> & IB);

    /// Get interpolation nodes
    /// This will return t_old_ and t_ provided t_old_ != t_
    bool GetNodes(std::vector<Scalar>* time_vec) const;

    /// Remove interpolation nodes
    /// This would allow removal of t_old_ and/or t_ which might be used to reject steps.
    bool RemoveNodes(std::vector<Scalar>& time_vec);

    /// Get order of interpolation
    /// This will return 1.
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
    Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > solver_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > scaled_x_old_;
    Scalar t_;
    Scalar t_old_;
    int numSteps;

    Teuchos::RefCountPtr<Thyra::SingleResidSSDAEModelEvaluator<Scalar> >  neModel_;

    Teuchos::RefCountPtr<Teuchos::ParameterList> parameterList_;

};

// ////////////////////////////
// Defintions

template<class Scalar>
BackwardEulerStepper<Scalar>::BackwardEulerStepper(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
  ,const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  out->precision(15);
  out->setMaxLenLinePrefix(32);
  out->pushLinePrefix("Rythmos::BackwardEulerStepper");
  out->setShowLinePrefix(true);
  out->setTabIndentStr("    ");

  setModel(model);
  setSolver(solver);
  numSteps = 0;
}

template<class Scalar>
void BackwardEulerStepper<Scalar>::setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model)
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"BES::setModel");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "model = " << model->description() << std::endl;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model_ = model;
  t_ = ST::zero();
  t_old_  = t_;
  x_ = model_->getNominalValues().get_x()->clone_v();

  scaled_x_old_ = x_->clone_v();
  
}

template<class Scalar>
void BackwardEulerStepper<Scalar>::setSolver(const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &solver)
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"BES::setSolver");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "solver = " << solver->description() << std::endl;
  solver_ = solver;
}

template<class Scalar>
Scalar BackwardEulerStepper<Scalar>::TakeStep()
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"BES::TakeStep()");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_LOW) )
    *out << "This is not valid for BackwardEulerStepper at this time." << std::endl;
  // print something out about this method not supporting automatic variable step-size
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return(-ST::one());
}

template<class Scalar>
Scalar BackwardEulerStepper<Scalar>::TakeStep(Scalar dt)
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"BES::TakeStep(dt)");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "dt = " << dt << std::endl;
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
  // 10/9/06 tscoffe:  I could use FE as a predictor here.
  Thyra::assign(&*x_,ST::zero());
  solver_->solve(&*x_); // Note that x in input is x_old and on output is the solved x!
  //
  // Update the step
  //
  t_ += dt;
  numSteps++;
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "t_old_ = " << t_old_ << std::endl;
    *out << "t_ = " << t_ << std::endl;
  }

  return(dt);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > BackwardEulerStepper<Scalar>::get_solution() const
{
  return(x_);
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
  if ( (static_cast<int>(verbLevel) == static_cast<int>(Teuchos::VERB_DEFAULT) ) ||
       (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)     )
     )
  {
    out << description() << "::describe:" << std::endl;
    out << "model = " << model_->description() << std::endl;
    out << "solver = " << solver_->description() << std::endl;
    out << "neModel = " << neModel_->description() << std::endl;
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
  {
    out << "t_ = " << t_ << std::endl;
    out << "t_old_ = " << t_old_ << std::endl;
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM))
  {
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH))
  {
    out << "model_ = " << std::endl;
    model_->describe(out,verbLevel);
    out << "solver_ = " << std::endl;
    solver_->describe(out,verbLevel);
    out << "neModel = " << std::endl;
    neModel_->describe(out,verbLevel);
    out << "x_ = " << std::endl;
    x_->describe(out,verbLevel);
    out << "scaled_x_old_ = " << std::endl;
    scaled_x_old_->describe(out,verbLevel);
  }
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::SetPoints(
    const std::vector<Scalar>& time_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
    ,const std::vector<ScalarMag> & accuracy_vec 
    )
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"BES::SetPoints");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "time_vec = " << std::endl;
    for (int i=0 ; i<time_vec.size() ; ++i)
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
  }
  typedef Teuchos::ScalarTraits<Scalar> ST;
  if (time_vec.size() == 0)
  {
    return(false);
  }
  else if (time_vec.size() == 1)
  {
    int n = 0;
    t_ = time_vec[n];
    t_old_ = t_;
    Thyra::V_V(&*x_,*x_vec[n]);
    Thyra::V_V(&*scaled_x_old_,*x_);
  }
  else 
  {
    int n = time_vec.size()-1;
    int nm1 = time_vec.size()-2;
    t_ = time_vec[n];
    t_old_ = time_vec[nm1];
    Thyra::V_V(&*x_,*x_vec[n]);
    Scalar dt = t_ - t_old_;
    Thyra::V_StV(&*scaled_x_old_,Scalar(-ST::one()/dt),*x_vec[nm1]);
  }
  return(true);
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::GetPoints(
    const std::vector<Scalar>& time_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
    ,std::vector<ScalarMag>* accuracy_vec) const
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"BES::GetPoints");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    for (int i=0 ; i<time_vec.size() ; ++i)
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
    *out << "I can interpolate in the interval [" << t_old_ << "," << t_ << "]." << std::endl;
  }
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typename DataStore<Scalar>::DataStoreVector_t ds_nodes;
  typename DataStore<Scalar>::DataStoreVector_t ds_out;
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot_temp; // Teuchos::null
  if (t_old_ != t_)
  {
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      *out << "Passing two points to interpolator:  " << t_old_ << " and " << t_ << std::endl;
    DataStore<Scalar> ds_temp;
    Scalar dt = t_ - t_old_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_temp = scaled_x_old_->clone_v();
    Thyra::Vt_S(&*x_temp,Scalar(-ST::one()*dt));  // undo the scaling
    ds_temp.time = t_old_;
    ds_temp.x = x_temp;
    ds_temp.xdot = xdot_temp;
    ds_temp.accuracy = ScalarMag(dt);
    ds_nodes.push_back(ds_temp);
    ds_temp.time = t_;
    ds_temp.x = x_;
    ds_temp.xdot = xdot_temp;
    ds_temp.accuracy = ScalarMag(dt);
    ds_nodes.push_back(ds_temp);
  }
  else
  {
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      *out << "Passing one point to interpolator:  " << t_ << std::endl;
    DataStore<Scalar> ds_temp;
    ds_temp.time = t_;
    ds_temp.x = x_;
    ds_temp.xdot = xdot_temp;
    ds_temp.accuracy = ScalarMag(ST::zero());
    ds_nodes.push_back(ds_temp);
  }
  LinearInterpolator<Scalar> interpolator;
  bool status = interpolator.interpolate(ds_nodes,time_vec,&ds_out);
  if (!status) return(status);
  std::vector<Scalar> time_out;
  DataStoreVectorToVector(ds_out,&time_out,x_vec,xdot_vec,accuracy_vec);
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "Passing out the interpolated values:" << std::endl;
    for (int i=0; i<time_out.size() ; ++i)
    {
      *out << "time[" << i << "] = " << time_out[i] << std::endl;
      *out << "x_vec[" << i << "] = " << std::endl;
      (*x_vec)[i]->describe(*out,Teuchos::VERB_EXTREME);
      if ( (*xdot_vec)[i] == Teuchos::null)
        *out << "xdot_vec[" << i << "] = Teuchos::null" << std::endl;
      else
      {
        *out << "xdot_vec[" << i << "] = " << std::endl;
        (*xdot_vec)[i]->describe(*out,Teuchos::VERB_EXTREME);
      }
      *out << "accuracy[" << i << "] = " << (*accuracy_vec)[i] << std::endl;
    }
  }
  return(status);
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::SetRange(
    const Scalar& time_lower
    ,const Scalar& time_upper
    ,const InterpolationBufferBase<Scalar>& IB)
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"BES::SetRange");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "time_lower = " << time_lower << std::endl;
    *out << "time_upper = " << time_upper << std::endl;
    *out << "IB = " << IB.description() << std::endl;
  }
  // TODO:
  // get node_list from IB, crop it to [time_lower,time_upper], crop x_vec to same,
  // pass to SetPoints.
  return(false);
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::GetNodes(std::vector<Scalar>* time_vec) const
{
  time_vec->clear();
  time_vec->push_back(t_old_);
  if (numSteps > 0)
    time_vec->push_back(t_);
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"BES::GetNodes");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << this->description() << std::endl;
    for (int i=0 ; i<time_vec->size() ; ++i)
    {
      *out << "time_vec[" << i << "] = " << (*time_vec)[i] << std::endl;
    }
  }
  return(true);
}

template<class Scalar>
bool BackwardEulerStepper<Scalar>::RemoveNodes(std::vector<Scalar>& time_vec) 
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"BES::RemoveNodes");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "time_vec = " << std::endl;
    for (int i=0 ; i<time_vec.size() ; ++i)
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
  }
  // TODO:
  // if any time in time_vec matches t_ or t_old_, then do the following:
  // remove t_old_:  set t_old_ = t_ and set scaled_x_old_ = x_
  // remove t_:  set t_ = t_old_ and set x_ = -dt*scaled_x_old_
  return(false);
}

template<class Scalar>
int BackwardEulerStepper<Scalar>::GetOrder() const
{
  return(1);
}

template <class Scalar>
void BackwardEulerStepper<Scalar>::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  parameterList_ = paramList;
  int outputLevel = parameterList_->get( "outputLevel", int(-1) );
  outputLevel = min(max(outputLevel,-1),4);
  this->setVerbLevel(static_cast<Teuchos::EVerbosityLevel>(outputLevel));
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> BackwardEulerStepper<Scalar>::getParameterList()
{
  return(parameterList_);
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> BackwardEulerStepper<Scalar>::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

} // namespace Rythmos

#endif //Rythmos_BACKWARD_EULER_STEPPER_H
