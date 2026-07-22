// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MATRIXBASED_LOWS_H
#define MATRIXBASED_LOWS_H

#include <Teuchos_RCP.hpp>
#include <Thyra_MultiVectorBase_decl.hpp>
#include <Thyra_LinearOpWithSolveBase_decl.hpp>




//! MatrixBased_LOWS provides a concrete implementation of LinearOpWithSolve based on an existing matrix
/*!
  * This class imports a given matrix (linear operator) and allows to initialize the solver
  * using a provided Stratimikos parameter list.
  */
class MatrixBased_LOWS : public Thyra::LinearOpWithSolveBase<double>
{
public:
  // Constructor
  MatrixBased_LOWS(
      const Teuchos::RCP<Thyra::LinearOpBase<double>> &matrix);

  //! Destructor
  virtual ~MatrixBased_LOWS();

  //! Overrides Thyra::LinearOWithSolvepBase purely virtual method
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> domain() const;

  //! Overrides Thyra::LinearOpWithSolveBase purely virtual method
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> range() const;

  //! Returns the matrix passed to the constructor
  Teuchos::RCP<Thyra::LinearOpBase<double>> getMatrix();

  //! Initialize the solver from a Stratimikos parameter list
  void initializeSolver(Teuchos::RCP<Teuchos::ParameterList> solverParamList);

  //@}

protected:
  //! Overrides Thyra::LinearOpWithSolveBase purely virtual method
  bool opSupportedImpl(Thyra::EOpTransp M_trans) const;

  //! Overrides Thyra::LinearOpWithSolveBase purely virtual method
  void applyImpl(const Thyra::EOpTransp M_trans,
                  const Thyra::MultiVectorBase<double> &X,
                  const Teuchos::Ptr<Thyra::MultiVectorBase<double>> &Y,
                  const double alpha,
                  const double beta) const;

  //! Overrides Thyra::LinearOpWithSolveBase purely virtual method
  Thyra::SolveStatus<double> solveImpl(
      const Thyra::EOpTransp transp,
      const Thyra::MultiVectorBase<double> &B,
      const Teuchos::Ptr<Thyra::MultiVectorBase<double>> &X,
      const Teuchos::Ptr<const Thyra::SolveCriteria<double>> solveCriteria) const;

  const Teuchos::RCP<Thyra::LinearOpBase<double>> mat_;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double>> solver_;
  //@}

}; 

#endif // MATRIXBASED_LOWS_H