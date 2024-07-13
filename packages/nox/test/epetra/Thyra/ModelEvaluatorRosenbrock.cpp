// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "ModelEvaluatorRosenbrock.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_CrsMatrix.h"

// Thyra support
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

// Thyra Epetra support
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"

RosenbrockModelEvaluator::
RosenbrockModelEvaluator(const Teuchos::RCP<const Epetra_Comm>& comm):
  d_(0.0),
  epetra_comm_(comm)
{
  using Teuchos::rcp;

  const int nx = 2;

  map_x_ = rcp(new Epetra_Map(nx,0,*epetra_comm_));
  x_space_ = Thyra::create_VectorSpace(map_x_);

  f_epetra_map_ = map_x_;
  f_space_ = x_space_;

  x0_ = rcp(new Epetra_Vector(*map_x_));
  (*x0_)[0] = -1.2;
  (*x0_)[1] = 1.0;

  x_analytic_ = rcp(new Epetra_Vector(*map_x_));
  (*x_analytic_)[0] = 1.0;
  (*x_analytic_)[1] = 1.0;

  // Initialize the graph for W CrsMatrix object
  W_graph_ = rcp(new Epetra_CrsGraph(::Copy,*map_x_,nx));
  {
    int indices[nx] = { 0, 1 };
    for( int i = 0; i < nx; ++i )
      W_graph_->InsertGlobalIndices(i,nx,indices);
  }
  W_graph_->FillComplete();

}

// Overridden from Thyra::ModelEvaulator

Teuchos::RCP<const Thyra::VectorSpaceBase<double> >
RosenbrockModelEvaluator::get_x_space() const
{
  return x_space_;
}

Teuchos::RCP<const Thyra::VectorSpaceBase<double> >
RosenbrockModelEvaluator::get_f_space() const
{
  return f_space_;
}

Thyra::ModelEvaluatorBase::InArgs<double>
RosenbrockModelEvaluator::getNominalValues() const
{
  nominal_values_ = this->createInArgs();
  nominal_values_.set_x(Thyra::create_Vector(x0_,x_space_));
  nominal_values_.set_x_dot(Teuchos::null);
  nominal_values_.set_alpha(0.0);
  nominal_values_.set_beta(1.0);
  return nominal_values_;
}

Teuchos::RCP<Thyra::LinearOpBase<double> >
RosenbrockModelEvaluator::create_W_op() const
{
  Teuchos::RCP<Epetra_CrsMatrix> W_epetra =
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));

  return Thyra::nonconstEpetraLinearOp(W_epetra);
}

Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> >
RosenbrockModelEvaluator::get_W_factory() const
{
  return W_factory_;
}

// Misc
Teuchos::RCP<const Epetra_Vector>
RosenbrockModelEvaluator::get_analytic_solution() const
{
  return x_analytic_;
}
void
RosenbrockModelEvaluator::
set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> >& W_factory)
{
  W_factory_ = W_factory;
}

Thyra::ModelEvaluatorBase::InArgs<double>
RosenbrockModelEvaluator::createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<double> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_x_dot,true);
  inArgs.setSupports(IN_ARG_t,true);
  inArgs.setSupports(IN_ARG_alpha,true);
  inArgs.setSupports(IN_ARG_beta,true);
  return inArgs;
}

Thyra::ModelEvaluatorBase::OutArgs<double>
RosenbrockModelEvaluator::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<double> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W_op,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_FULL
      ,false // supportsAdjoint
      )
    );
  return outArgs;
}

void
RosenbrockModelEvaluator::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<double>& inArgs,
          const Thyra::ModelEvaluatorBase::OutArgs<double>& outArgs) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  //
  // Get the input arguments
  //

  TEUCHOS_ASSERT(nonnull(inArgs.get_x()));
  const Thyra::ConstDetachedVectorView<double> x(inArgs.get_x());

  Teuchos::RCP<const Thyra::VectorBase<double> > x_dot_rcp = inArgs.get_x_dot();
  double alpha = inArgs.get_alpha();
  double beta = inArgs.get_beta();
  //double time = inArgs.get_t();

  //
  // Get the output arguments
  //
  const Teuchos::RCP< Thyra::VectorBase<double> > f_out = outArgs.get_f();
  const Teuchos::RCP< Thyra::LinearOpBase< double > > W_out = outArgs.get_W_op();

  //
  // Compute the functions
  //
  if(nonnull(f_out)) {
    if (nonnull(x_dot_rcp)) {
      NOX_FUNC_TIME_MONITOR("RosenbrockModelEvaluator::eval f_out");
      Thyra::DetachedVectorView<double> f(f_out);
      const Thyra::ConstDetachedVectorView<double> x_dot(x_dot_rcp);
      f[0] = x_dot[0] + 10.0 * ( x[1] - x[0] * x[0] );
      f[1] = x_dot[1] + 1.0 - x[0];
      std::cout << "Residual Evaluation: Transient" << std::endl;
    }
    else {
      NOX_FUNC_TIME_MONITOR("RosenbrockModelEvaluator::eval f_out");
      const Thyra::DetachedVectorView<double> f(f_out);
      f[0] = 10.0 * ( x[1] - x[0] * x[0] );
      f[1] = 1.0 - x[0];
      std::cout << "Residual Evaluation: Steady-State" << std::endl;
    }
  }
  if(nonnull(W_out)) {
    Teuchos::RCP<Epetra_Operator> W_epetra= Thyra::get_Epetra_Operator(*W_out);
    Teuchos::RCP<Epetra_CrsMatrix> W_epetracrs = rcp_dynamic_cast<Epetra_CrsMatrix>(W_epetra);
    TEUCHOS_ASSERT(nonnull(W_epetracrs));
    Epetra_CrsMatrix& DfDx = *W_epetracrs;

    DfDx.PutScalar(0.0);

    int indexes[2];
    indexes[0] = 0;
    indexes[1] = 1;

    double values[2];

    if (nonnull(x_dot_rcp)) {
      // Augment Jacobian for Pseudo-transient, use identity for V
      // Row [0]
      values[0] = alpha * 1.0 + beta * ( -20.0 * x[0] );
      values[1] =               beta * 10.0;
      DfDx.SumIntoGlobalValues( 0, 2, values, indexes );
      // Row [1]
      values[0] =               beta * (-1.0);
      values[1] = alpha * 1.0 + beta * 0.0;
      DfDx.SumIntoGlobalValues( 1, 2, values, indexes );
      std::cout << "Jacobian Evaluation: Transient, alpha=" << alpha << ", beta=" << beta << std::endl;
    }
    else {
      // Row [0]
      values[0] = -20.0 * x[0];
      values[1] = 10.0;
      DfDx.SumIntoGlobalValues( 0, 2, values, indexes );
      // Row [1]
      values[0] = -1.0;
      values[1] = 0.0;
      DfDx.SumIntoGlobalValues( 1, 2, values, indexes );
      std::cout << "Jacobian Evaluation: Steady-State, alpha=" << alpha << ", beta=" << beta << std::endl;
    }

  }
}
