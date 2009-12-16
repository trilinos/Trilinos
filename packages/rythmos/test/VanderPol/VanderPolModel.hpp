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

#ifndef VAN_DER_POL_MODEL_HPP
#define VAN_DER_POL_MODEL_HPP

#include "Rythmos_ConfigDefs.h"

#ifdef Rythmos_ENABLE_Sacado

#include "Rythmos_Types.hpp"

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

using Thyra::ModelEvaluatorBase;

namespace Rythmos {

  /*
   * This is the Van der Pol Equation from 
   * "Solving Differential Equations I, Nonstiff Problems" by
   * E. Hairer, S.P. Norsett, and G. Wanner.  2000, Second edition, Springer.
   * 
   * \ddot{x_1} = x_2
   * \ddot{x_2} = \epsilon*(1-x_1^2)*x_2-x_1, \epsilon > 0
   *
   * Hairer, etal. show the solution x(t) for large \epsilon:
   *
   * log(x_1) - \frac{x_1^2}{2} = \frac{t-t_0}{\epsilon} + C
   *
   * For small \epsilon, the exact solution is more complicated.
   * 
   * There is one parameter in this model, \epsilon.
   *
   */

class VanderPolModel
  : public Thyra::StateFuncModelEvaluatorBase<double>,
    public Teuchos::ParameterListAcceptorDefaultBase
{
  public:

  // Constructor
  VanderPolModel();

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
  void setupInOutArgs_();

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

private:

  int dim_;         // Number of state unknowns (2)
  int Np_;          // Number of parameter vectors (1)
  int np_;          // Number of parameters in this vector (1)
  int Ng_;          // Number of observation functions (0)
  int ng_;          // Number of elements in this observation function (0)
  bool isImplicit_; // false => \dot{x} = f(x,t)    W = beta*df/dx
                    // true =>  F(\dot{x},x,t) = 0  W = alpha*dF/dxdot + beta*dF/dx
  bool haveIC_;     // false => no nominal values are provided (default=true)
  bool acceptModelParams_; // Changes inArgs to require parameters
  bool isInitialized_;
  ModelEvaluatorBase::InArgs<double> inArgs_;
  ModelEvaluatorBase::OutArgs<double> outArgs_;
  ModelEvaluatorBase::InArgs<double> nominalValues_;
  RCP<const Thyra::VectorSpaceBase<double> > x_space_;
  RCP<const Thyra::VectorSpaceBase<double> > f_space_;
  RCP<const Thyra::VectorSpaceBase<double> > p_space_;
  RCP<const Thyra::VectorSpaceBase<double> > g_space_;

  // Parameters for the model:  
  double epsilon_; // This is a model parameter 
  double t0_ic_; // initial time
  double x0_ic_; // initial condition for x0
  double x1_ic_; // initial condition for x1

  template<class ScalarT>
  void eval_f(
    const ArrayView<const ScalarT> &x_dot,
    const ArrayView<const ScalarT> &x,
    const ScalarT &eps,
    const ScalarT &t,
    const ArrayView<ScalarT> &f
    ) const;

};

// Non-member constructors
RCP<VanderPolModel> vanderPolModel();
RCP<VanderPolModel> vanderPolModel(bool implicit);
RCP<VanderPolModel> vanderPolModel(const RCP<ParameterList> &pl);

} // namespace Rythmos 

#endif //Rythmos_ENABLE_Sacado
#endif // VAN_DER_POL_MODEL_HPP
