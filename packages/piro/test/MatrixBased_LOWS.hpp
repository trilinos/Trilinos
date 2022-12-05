// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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