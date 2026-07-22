// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_THYRA_ROSENBROCK_MODEL_EVALUATOR_HPP
#define NOX_THYRA_ROSENBROCK_MODEL_EVALUATOR_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

/** \brief ModelEvaluator for the Roesenbrock problem

   This example is the "%Rosenbrock function" from Jorge J. More',
   Burton S. Garbow, and Kenneth E. Hillstrom, Testing Unconstrained
   Optimization Software, ACM TOMS, Vol. 7, No. 1, March 1981,
   pp. 14-41.

   It comes originally from H. H. %Rosenbrock, An Automatic Method for
   Finding the Greatest or Least Value of a Function, J. Comput.
   3(1960):175-184.

   The function is defined as
   \f[
   F(x) = \left[
   \begin{array}{c}
   10 (x[2] - x[1]^2) \\
   1 - x[1]
   \end{array}
   \right]
   \f]

   The initial guess is given by
   \f[
   x_0 = \left[
   \begin{array}{c}
   -1.2\\
   1
   \end{array}
   \right]
   \f]

   The solution is
   \f[
   x_* = \left[
   \begin{array}{c}
   1\\
   1
   \end{array}
   \right]
   \f]
 */
class RosenbrockModelEvaluator : public Thyra::StateFuncModelEvaluatorBase<double> {
public:

  RosenbrockModelEvaluator(const Teuchos::RCP<const Epetra_Comm>& comm);

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > get_x_space() const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > get_f_space() const;
  Thyra::ModelEvaluatorBase::InArgs<double> getNominalValues() const;
  Teuchos::RCP<Thyra::LinearOpBase<double> > create_W_op() const;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > get_W_factory() const;
  Thyra::ModelEvaluatorBase::InArgs<double> createInArgs() const;

  //@}

  void set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> >& W_factory);

  Teuchos::RCP<const Epetra_Vector> get_analytic_solution() const;

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  Thyra::ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const;

  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs
    ) const;

  //@}

private:

  double    d_;

  Teuchos::RCP<const Epetra_Comm>  epetra_comm_;

  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > x_space_;
  Teuchos::RCP<const Epetra_Map>   map_x_;

  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > f_space_;
  Teuchos::RCP<const Epetra_Map>   f_epetra_map_;

  Teuchos::RCP<Epetra_Vector> x0_;
  Teuchos::RCP<Epetra_Vector> x_analytic_;

  Teuchos::RCP<Epetra_CrsGraph>  W_graph_;

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > W_factory_;

  mutable Thyra::ModelEvaluatorBase::InArgs<double> nominal_values_;
};

#endif
