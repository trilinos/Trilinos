//@HEADER

// ***********************************************************************
//
//                     Rythmos Package
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

#ifndef SIN_COS_MODEL_HPP
#define SIN_COS_MODEL_HPP

#include "Rythmos_ConfigDefs.h"
#include "Rythmos_Types.hpp"

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

using Thyra::ModelEvaluatorBase;

namespace Rythmos {

  /*
   * This is the canonical Sine Cosine differential equation 
   * 
   * \ddot{x} = -x 
   * 
   * with a few enhancements.
   * We start with the exact solution to the differential equation as:
   *
   * x0(t) = a + b*sin((f/L)*t+phi)
   * x1(t) =   b*(f/L)*cos((f/L)*t+phi)
   *
   * Then the form of the model is:
   *
   * d/dt x0(t) = x1(t)
   * d/dt x1(t) = (f/L)*(f/L)*(a-x0(t)) [a=0,f=1,L=1]
   *
   * With Initial conditions:
   *
   * x0(t0=0) = gamma0 [0.0]
   * x1(t0=0) = gamma1 [1.0]
   *
   * We can use gamma0 and gamma1 to solve for phi and b:
   *
   * phi = atan(((f/L)/gamma1)*(gamma0-a))-(f/L)*t0 [0.0]
   * b = gamma1/((f/L)*cos((f/L)*t0+phi)) [1.0]
   *
   * Therefore this model has three model parameters and two initial conditions
   * which effect the exact solution as above.
   *
   * p = (a, f, L)
   *
   * \dot{x}=F(x,t,p)
   * F_0 = x1
   * F_1 = (f/L)^2*(a-x0)
   *
   */

class SinCosModel 
  : public Thyra::StateFuncModelEvaluatorBase<double>,
    public Teuchos::ParameterListAcceptorDefaultBase
{
  public:

  // Constructor
  SinCosModel();

  // Exact solution
  ModelEvaluatorBase::InArgs<double> getExactSolution(double t) const;

  // Exact sensitivity solution
  ModelEvaluatorBase::InArgs<double> getExactSensSolution(int j, double t) const;

  // Set explicit/implicit flag
  void setImplicitFlag(bool implicit);

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<double> > get_x_space() const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<double> > get_f_space() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> getNominalValues() const;
  /** \brief . */
  RCP<Thyra::LinearOpWithSolveBase<double> > create_W() const;
  /** \brief . */
  RCP<Thyra::LinearOpBase<double> > create_W_op() const;
  /** \brief . */
  RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > get_W_factory() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<double> createInArgs() const;

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<double> > get_p_space(int l) const;
  /** \brief . */
  RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<double> > get_g_space(int j) const;

  //@}
  
  /** \name Public functions overridden from ParameterListAcceptor. */
  //@{
  
  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

private:

  /** \brief. */
  void setupInOutArgs_() const;

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<double> &inArgs_bar,
    const ModelEvaluatorBase::OutArgs<double> &outArgs_bar
    ) const;

  //@}

  void calculateCoeffFromIC_();

private:
  int dim_;         // Number of state unknowns (2)
  int Np_;          // Number of parameter vectors (1)
  int np_;          // Number of parameters in this vector (2)
  int Ng_;          // Number of observation functions (1)
  int ng_;          // Number of elements in this observation function (1)
  bool isImplicit_; // false => \dot{x} = f(x,t)    W = beta*df/dx
                    // true =>  F(\dot{x},x,t) = 0  W = alpha*dF/dxdot + beta*dF/dx
  bool haveIC_;     // false => no nominal values are provided (default=true)
  bool acceptModelParams_; // Changes inArgs to require parameters
  mutable bool isInitialized_;
  mutable ModelEvaluatorBase::InArgs<double> inArgs_;
  mutable ModelEvaluatorBase::OutArgs<double> outArgs_;
  mutable ModelEvaluatorBase::InArgs<double> nominalValues_;
  RCP<const Thyra::VectorSpaceBase<double> > x_space_;
  RCP<const Thyra::VectorSpaceBase<double> > f_space_;
  RCP<const Thyra::VectorSpaceBase<double> > p_space_;
  RCP<const Thyra::VectorSpaceBase<double> > g_space_;

  // Parameters for the model:  x_0(t) = a + b*sin(f*t+phi)
  //                            x_1(t) = b*f*cos(f*t+phi)
  double a_; // This is a model parameter 
  double f_; // This is a model parameter
  double L_; // This is a model parameter
  double phi_; // This is a parameter determined from the IC
  double b_; // This is a parameter determined from the IC
  double t0_ic_; // This is the time value where the initial condition is specified
  double x0_ic_; // initial condition for x0
  double x1_ic_; // initial condition for x1
};

// Non-member constructor
RCP<SinCosModel> sinCosModel(bool implicit);
RCP<SinCosModel> sinCosModel();


} // namespace Rythmos 

#endif // SIN_COS_MODEL_HPP
