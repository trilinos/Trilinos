// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
