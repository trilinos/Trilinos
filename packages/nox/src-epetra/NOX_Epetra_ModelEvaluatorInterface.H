// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
    
  protected:
    
    //! RCP to the Model Evaluator.
    Teuchos::RCP<EpetraExt::ModelEvaluator> model_;
    
    //! Inargs object.
    EpetraExt::ModelEvaluator::InArgs inargs_;
    
    //! Outargs object.
    EpetraExt::ModelEvaluator::OutArgs outargs_;
    

    //! RCP for the solution vector.
    Teuchos::RCP<const Epetra_Vector> x_;

    //! RCP for the residual vector.
    Teuchos::RCP<Epetra_Vector> f_;

    //! RCP for the Jacobian operator.
    Teuchos::RCP<Epetra_Operator> jacobian_;

    //! RCP for the residual vector wrapper.
    EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> eval_f_;
    
    
  };
} // namespace Epetra
} // namespace NOX

#endif
