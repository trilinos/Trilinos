// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include <cmath>

#include "Piro_Epetra_InvertMassMatrixDecorator.hpp"

#include "Thyra_VectorBase.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_VerboseObject.hpp"


Piro::Epetra::InvertMassMatrixDecorator::InvertMassMatrixDecorator(
                          Teuchos::RCP<Teuchos::ParameterList> stratParams,
                          Teuchos::RCP<EpetraExt::ModelEvaluator>& model_,
                          bool massMatrixIsConstant_, bool lumpMassMatrix_) :
  model(model_),
  calcMassMatrix(true),
  massMatrixIsConstant(massMatrixIsConstant_),
  lumpMassMatrix(lumpMassMatrix_)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  // Create x_dot vector, fill with 0.0 so implicit fill gives
  // correct results for explicit fill
  x_dot = Teuchos::rcp(new Epetra_Vector(*(model->get_x_map())));
  x_dot->PutScalar(0.0);

  // get allocated space for Mass Matrix
  massMatrix = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix> (model->create_W(), true);
  if (lumpMassMatrix) invDiag = Teuchos::rcp(new Epetra_Vector(*(model->get_x_map())));

  Teuchos::RCP<Teuchos::FancyOStream> out
     = Teuchos::VerboseObjectBase::getDefaultOStream();

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
 
  linearSolverBuilder.setParameterList(stratParams);

  // Create a linear solver factory given information read from the
  // parameter list.
  lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");

  // Setup output stream and the verbosity level
  lowsFactory->setOStream(out);
  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

}

Piro::Epetra::InvertMassMatrixDecorator::~InvertMassMatrixDecorator()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::InvertMassMatrixDecorator::get_x_map() const
{
  return model->get_x_map();
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::InvertMassMatrixDecorator::get_f_map() const
{
  return model->get_f_map();
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::InvertMassMatrixDecorator::get_p_map(int l) const
{
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::InvertMassMatrixDecorator::get_g_map(int j) const
{
  return model->get_g_map(j);
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::InvertMassMatrixDecorator::get_x_init() const
{
  return model->get_x_init();
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::InvertMassMatrixDecorator::get_x_dot_init() const
{
  return model->get_x_dot_init();
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::InvertMassMatrixDecorator::get_p_init(int l) const
{
  return model->get_p_init(l);
}

Teuchos::RCP<Epetra_Operator> Piro::Epetra::InvertMassMatrixDecorator::create_W() const
{
  return model->create_W();
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::InvertMassMatrixDecorator::createInArgs() const
{
  return model->createInArgs();
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::InvertMassMatrixDecorator::createOutArgs() const
{
  return model->createOutArgs();
}

void Piro::Epetra::InvertMassMatrixDecorator::evalModel( const InArgs& inArgs,
                                     const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  if (outArgs.Np()>0) {
   if (outArgs.get_DfDp(0).getMultiVector() != Teuchos::null) 
     cout << "InvertMassMatrixDecorator:: NOT IMPLEMENTED FOR dfdp!! " << endl;
  }

  if (outArgs.get_f() == Teuchos::null) {
    // Probably just getting g -- pass through
    model->evalModel(inArgs, outArgs);
  }
  else {
    // Add massMatrix to outargs, with appropriate alpha and beta
    OutArgs modelOutArgs(outArgs);
    InArgs modelInArgs(inArgs);

    modelInArgs.set_x_dot(x_dot);

    
    if (calcMassMatrix) {
      modelOutArgs.set_W(massMatrix);
      modelInArgs.set_alpha(-1.0); 
      modelInArgs.set_beta(0.0);
    }

    //Evaluate the underlying model
    model->evalModel(modelInArgs, modelOutArgs);

    // Invert the mass matrix:   f_exp = M^{-1} f_imp

    if (!lumpMassMatrix) {
      // Create a linear solver based on the forward operator A
      if (calcMassMatrix) {
        A = Thyra::epetraLinearOp( massMatrix );
        lows = Thyra::linearOpWithSolve(*lowsFactory, A);
      }

      // Solve the linear system for x, given b 
      RCP<Thyra::VectorBase<double> >
        x = Thyra::create_Vector( outArgs.get_f(), A->domain() );
      RCP<const Thyra::VectorBase<double> >
        b = Thyra::create_Vector( modelOutArgs.get_f(), A->range() );

      lows->solve(Thyra::NONCONJ_ELE, *b, x.get());
    }
    else { // Lump matrix into inverse of diagonal
      if (calcMassMatrix) {
        invDiag->PutScalar(1.0);
        massMatrix->Multiply(false,*invDiag, *invDiag);
        invDiag->Reciprocal(*invDiag);
      }
      outArgs.get_f()->Multiply(1.0, *invDiag, *modelOutArgs.get_f(), 0.0);
    }

    // Do not recompute mass matrix in future if it is a constant
    if (massMatrixIsConstant) calcMassMatrix = false;
  }
}
