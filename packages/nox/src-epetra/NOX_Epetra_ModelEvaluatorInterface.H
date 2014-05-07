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

#ifndef NOX_EPETRA_MODEL_EVALUATOR_INTERFACE_H
#define NOX_EPETRA_MODEL_EVALUATOR_INTERFACE_H

#include "NOX_Common.H"

#ifdef HAVE_NOX_EPETRAEXT

#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "Teuchos_RCP.hpp"
#include "EpetraExt_ModelEvaluator.h"

namespace NOX {
namespace Epetra {

  /*!
    \brief Wrapper for an EpetraExt::ModelEvaluator

    If an application interfaces their code to solvers using the EpetraExt::ModelEvaluator, this class provides a wrapper so that the model evaluator can be used instead of having the user write concrete versions of the NOX::Epetra::Interface objects.

  */
  class ModelEvaluatorInterface : public NOX::Epetra::Interface::Required,
                  public NOX::Epetra::Interface::Jacobian,
                                  public NOX::Epetra::Interface::Preconditioner
  {

  public:

    //! Constructor
    ModelEvaluatorInterface(const Teuchos::RCP<EpetraExt::ModelEvaluator>& m);

    //! Destructor
    virtual ~ModelEvaluatorInterface();

    virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F,
              const FillType fillFlag);

    virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

    virtual bool computePreconditioner(const Epetra_Vector& x,
                       Epetra_Operator& M,
                       Teuchos::ParameterList* precParams = 0);

    //! Set the parameter vector in inargs of the underlying model evaluator
    virtual bool inargs_set_p(const Teuchos::RCP<const Epetra_Vector> p,
                              const int l = 0);

    //! Set a general inargs for the underlying model evaluator
    virtual bool set_inargs(const EpetraExt::ModelEvaluator::InArgs& inargs_in);

  protected:

    //! RCP to the Model Evaluator.
    Teuchos::RCP<EpetraExt::ModelEvaluator> model_;

    //! Inargs object.
    EpetraExt::ModelEvaluator::InArgs inargs_;

    //! RCP for the residual vector wrapper.
    EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> eval_f_;


  };
} // namespace Epetra
} // namespace NOX

#endif //HAVE_NOX_EPETRAEXT
#endif
