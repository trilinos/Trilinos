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

#ifndef PIRO_MATRIXFREEDECORATOR_HPP
#define PIRO_MATRIXFREEDECORATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

#include "Teuchos_RCP.hpp"

namespace Piro {

/** \brief Decorator class that creates a Jacobian (W) operator
 * using matrix-free directional derivatives.
 *
 * This class wraps a model evaluator supporting residual calculation
 * and adds support for the W matrix that implements a Thyra-based version of
 * NOX_Epetra_MatrixFree conforming to the Thyra::LinerarOp interface.
 *
 * This class supports time-dependent problems (characterized by x_dot != null)
 * and uses the input values of the alpha and beta coefficients.
 */

template <typename Scalar>
class MatrixFreeDecorator : public Thyra::ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \name Constructors/initializers */
  //@{
  /** \brief . */
  explicit MatrixFreeDecorator(Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model);
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  //@}

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  /** \brief . */
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

  //@}
};

}

#include "Piro_MatrixFreeDecorator_Def.hpp"

#endif /* PIRO_MATRIXFREEDECORATOR_HPP */
