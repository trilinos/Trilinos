// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "FEApp_ModelEvaluator.hpp"
#include "Teuchos_ScalarTraits.hpp"

FEApp::ModelEvaluator::ModelEvaluator(
		const Teuchos::RefCountPtr<FEApp::Application>& app_) :
  app(app_)
{
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RefCountPtr<const Epetra_Map>
FEApp::ModelEvaluator::get_x_map() const
{
  return app->getMap();
}

Teuchos::RefCountPtr<const Epetra_Map>
FEApp::ModelEvaluator::get_f_map() const
{
  return app->getMap();
}

Teuchos::RefCountPtr<const Epetra_Vector>
FEApp::ModelEvaluator::get_x_init() const
{
  return app->getInitialSolution();
}

Teuchos::RefCountPtr<Epetra_Operator>
FEApp::ModelEvaluator::create_W() const
{
  return 
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy, *(app->getJacobianGraph())));
}

EpetraExt::ModelEvaluator::InArgs
FEApp::ModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
FEApp::ModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_UNKNOWN ,DERIV_RANK_FULL ,true)
    );
  return outArgs;
}

void 
FEApp::ModelEvaluator::evalModel(const InArgs& inArgs, 
				 const OutArgs& outArgs) const
{
  //
  // Get the input arguments
  //
  const Epetra_Vector& x = *inArgs.get_x();

  //
  // Get the output arguments
  //
  Teuchos::RefCountPtr<Epetra_Vector> f_out = outArgs.get_f();
  Teuchos::RefCountPtr<Epetra_Operator> W_out = outArgs.get_W();

  // Cast W to a CrsMatrix, throw an exception if this fails
  Teuchos::RefCountPtr<Epetra_CrsMatrix> W_out_crs;
  if (W_out != Teuchos::null)
    W_out_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out, true);
  
  //
  // Compute the functions
  //
  if(f_out != Teuchos::null && W_out != Teuchos::null) {
    app->computeGlobalJacobian(x, *f_out, *W_out_crs);
  }
  else if(f_out != Teuchos::null ) {
    app->computeGlobalResidual(x, *f_out);
  }
  else if(W_out != Teuchos::null ) {
    Epetra_Vector f(x.Map());
    app->computeGlobalJacobian(x, f, *W_out_crs);
  }
}
