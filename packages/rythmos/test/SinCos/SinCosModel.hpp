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

#include "Teuchos_ParameterList.hpp"

using Thyra::ModelEvaluatorBase;

namespace Rythmos {

class SinCosModel : public Thyra::StateFuncModelEvaluatorBase<double> 
{
  public:

  // Constructor
  SinCosModel();

  // Exact solution
  ModelEvaluatorBase::InArgs<double> getExactSolution(double t) const;

  // Set explicit/implicit flag
  void setImplicitFlag(bool implicit);

  // Set have/not-have initial condition
  void setHaveIC(bool haveIC);

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

  //@}

private:

  /** \brief. */
  void initialize_();

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
  int dim_;         // 2
  bool isImplicit_; // false => \dot{x} = f(x,t)    W = beta*df/dx
                    // true =>  F(\dot{x},x,t) = 0  W = alpha*dF/dxdot + beta*dF/dx
  bool haveIC_;     // false => no nominal values are provided (default=true)
  bool isInitialized_;
  ModelEvaluatorBase::InArgs<double> inArgs_;
  ModelEvaluatorBase::OutArgs<double> outArgs_;
  ModelEvaluatorBase::InArgs<double> nominalValues_;
  RCP<const Thyra::VectorSpaceBase<double> > x_space_;
  RCP<const Thyra::VectorSpaceBase<double> > f_space_;

};

// Non-member constructor
RCP<SinCosModel> sinCosModel(bool implicit, bool haveIC=true);


} // namespace Rythmos 

#endif // SIN_COS_MODEL_HPP
