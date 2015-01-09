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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef LOG_TIME_MODEL_HPP
#define LOG_TIME_MODEL_HPP

#include "Rythmos_ConfigDefs.h"
#include "Rythmos_Types.hpp"

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

using Thyra::ModelEvaluatorBase;

namespace Rythmos {

  /*
   * This is the ordinary differential equation which has a logarithmic
   * in time solution.  This stresses the variable time integration to
   * see if it can increase the time step "quickly enough" to reduce
   * overall computation.
   *
   * \ddot{x} =
   *      (a*t^3*(8*b^2*d + 8*c*d*t + b*sqrt(t)*((7 + 9*c)*d + (-1 + c)*t^4)))/
   *      (2.*(b + Sqrt(t))^2*(d + t^4)^2)
   *
   * where a=1.4, b=0.0001, c=0.1, and d=1.0e-36.
   * The exact solution to the differential equation is
   *
   * x(t) = (a*(b*t^4 + c*t^4.5))/((b + sqrt(t))*(d + t^4))
   *
   * with initial conditions:
   *
   * x(t0=0) = 0
   *
   */

class LogTimeModel
  : public Thyra::StateFuncModelEvaluatorBase<double>,
    public Teuchos::ParameterListAcceptorDefaultBase
{
  public:

  // Constructor
  LogTimeModel();

  // Exact solution
  ModelEvaluatorBase::InArgs<double> getExactSolution(double t) const;

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
  int dim_;         // Number of state unknowns (1)
  int Np_;          // Number of parameter vectors (1)
  int np_;          // Number of parameters in this vector (1)
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

  // Parameters for the model
  double a_; // This is a model parameter
  double b_; // This is a model parameter
  double c_; // This is a model parameter
  double d_; // This is a model parameter
  double t0_ic_; // This is the time value where the initial condition is specified
  double x_ic_; // initial condition for x0
};

// Non-member constructor
RCP<LogTimeModel> logTimeModel(bool implicit);
RCP<LogTimeModel> logTimeModel();


} // namespace Rythmos

#endif // LOG_TIME_MODEL_HPP
