// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Epetra_ModelEvaluatorInterface.H"

#ifdef HAVE_NOX_EPETRAEXT

#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

// *****************************************************************
// *****************************************************************
NOX::Epetra::ModelEvaluatorInterface::
ModelEvaluatorInterface(const
            Teuchos::RCP<EpetraExt::ModelEvaluator>& m) :
  model_(m)
{
  inargs_ = model_->createInArgs();
}

// *****************************************************************
// *****************************************************************
NOX::Epetra::ModelEvaluatorInterface::~ModelEvaluatorInterface()
{
}

// *****************************************************************
// *****************************************************************
bool NOX::Epetra::ModelEvaluatorInterface::
computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  EpetraExt::ModelEvaluator::OutArgs outargs = model_->createOutArgs();

  inargs_.set_x( Teuchos::rcp(&x, false) );

  Teuchos::RCP<Epetra_Vector> f = Teuchos::rcp(&F, false);

  if (fillFlag == NOX::Epetra::Interface::Required::Residual)
    eval_f_.reset(f, EpetraExt::ModelEvaluator::EVAL_TYPE_EXACT);
  else if (fillFlag == NOX::Epetra::Interface::Required::Jac)
    eval_f_.reset(f, EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV);
  else
    eval_f_.reset(f, EpetraExt::ModelEvaluator::EVAL_TYPE_VERY_APPROX_DERIV);

  outargs.set_f(eval_f_);

  model_->evalModel(inargs_, outargs);

  inargs_.set_x(Teuchos::null);

  return true;
}

// *****************************************************************
// *****************************************************************
bool NOX::Epetra::ModelEvaluatorInterface::
computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jaco)
{
  EpetraExt::ModelEvaluator::OutArgs outargs = model_->createOutArgs();

  inargs_.set_x( Teuchos::rcp(&x, false) );
  outargs.set_W( Teuchos::rcp(&Jaco, false) );

  model_->evalModel(inargs_, outargs);

  inargs_.set_x(Teuchos::null);

  return true;
}

// *****************************************************************
// *****************************************************************
bool NOX::Epetra::ModelEvaluatorInterface::
computePreconditioner(const Epetra_Vector& x,
              Epetra_Operator& M,
              Teuchos::ParameterList* /* precParams */)
{
  EpetraExt::ModelEvaluator::OutArgs outargs = model_->createOutArgs();

  inargs_.set_x( Teuchos::rcp(&x, false) );
  outargs.set_WPrec( Teuchos::rcp(&M, false) );

  model_->evalModel(inargs_, outargs);

  inargs_.set_x(Teuchos::null);

  return true;
}

// *****************************************************************
// *****************************************************************
bool NOX::Epetra::ModelEvaluatorInterface::
inargs_set_p(const Teuchos::RCP<const Epetra_Vector> p_, const int l)
{
   inargs_.set_p(l, p_);
   return true;
}

// *****************************************************************
// *****************************************************************
bool NOX::Epetra::ModelEvaluatorInterface::
set_inargs(const EpetraExt::ModelEvaluator::InArgs& inargs_in)
{
   inargs_ = inargs_in;
   return true;
}

// *****************************************************************
// *****************************************************************
#endif //HAVE_NOX_EPETRAEXT
