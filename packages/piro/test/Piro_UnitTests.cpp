/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "MockModelEval_A.hpp"

#include "Piro_ConfigDefs.hpp"
#ifdef Piro_ENABLE_NOX
#include "Piro_Epetra_NOXSolver.hpp"
#endif

namespace {

#ifdef Piro_ENABLE_NOX

void testSensitivities(const std::string& inputFile, 
		       bool use_op, bool use_op_trans,
		       Teuchos::FancyOStream& out, bool& success)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

#ifdef HAVE_MPI
  MPI_Comm appComm = MPI_COMM_WORLD;
#else
  int appComm=0;
#endif
  // Create (1) a Model Evaluator and (2) a parmeter list
  RCP<EpetraExt::ModelEvaluator> Model = rcp(new MockModelEval_A(appComm));

  RCP<Teuchos::ParameterList> piroParams =
    rcp(new Teuchos::ParameterList("Piro Parameters"));
  Teuchos::updateParametersFromXmlFile(inputFile, piroParams.get());
  
  // Use these two objects to construct a Piro solved application 
  //   EpetraExt::ModelEvaluator is  base class of all Piro::Epetra solvers
  RCP<EpetraExt::ModelEvaluator> piro = 
    rcp(new Piro::Epetra::NOXSolver(piroParams, Model));

  // Now the (somewhat cumbersome) setting of inputs and outputs
  EpetraExt::ModelEvaluator::InArgs inArgs = piro->createInArgs();
  RCP<Epetra_Vector> p1 = rcp(new Epetra_Vector(*(piro->get_p_init(0))));
  int numParams = p1->MyLength(); // Number of parameters in p1 vector
  inArgs.set_p(0,p1);

  // Set output arguments to evalModel call
  EpetraExt::ModelEvaluator::OutArgs outArgs = piro->createOutArgs();
  int num_g = outArgs.Ng(); // Number of *vectors* of responses
  RCP<Epetra_Vector> g1 = rcp(new Epetra_Vector(*(piro->get_g_map(0))));
  RCP<Epetra_Vector> gx = rcp(new Epetra_Vector(*(piro->get_g_map(1))));
  RCP<Epetra_MultiVector> dgdp_mv;
  RCP<Epetra_Operator> dgdp_op;
  outArgs.set_g(0,g1);
  outArgs.set_g(1,gx);
  if (use_op || use_op_trans) {
    dgdp_op = piro->create_DgDp_op(0,0);
    outArgs.set_DgDp(0, 0, dgdp_op);
  }
  else{
    dgdp_mv = rcp(new Epetra_MultiVector(g1->Map(), numParams));
    outArgs.set_DgDp(0, 0, dgdp_mv);
  }

  // Now, solve the problem and return the responses
  piro->evalModel(inArgs, outArgs);

  double tol = 1e-6;
  if (use_op) {
    Epetra_Vector v(p1->Map()), w(g1->Map());
    v.PutScalar(1.0);
    dgdp_op->Apply(v, w);
    TEUCHOS_TEST_FLOATING_EQUALITY(w[0], -6.0, tol, out, success);
  }
  else if (use_op_trans) {
    Epetra_Vector v(g1->Map()), w(p1->Map());
    v.PutScalar(1.0);
    dgdp_op->SetUseTranspose(true);
    dgdp_op->Apply(v, w);
    TEUCHOS_TEST_FLOATING_EQUALITY(w[0],  2.0, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(w[1], -8.0, tol, out, success);
  }
  else {
    TEUCHOS_TEST_FLOATING_EQUALITY((*dgdp_mv)[0][0],  2.0, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY((*dgdp_mv)[1][0], -8.0, tol, out, success);
  }
}

TEUCHOS_UNIT_TEST( Piro, ForwardSensitivities )
{
  std::string inputFile="input_Solve_NOX_1.xml";
  testSensitivities(inputFile, false, false, out, success);
}

TEUCHOS_UNIT_TEST( Piro, AdjointSensitivities )
{
  std::string inputFile="input_Solve_NOX_Adjoint.xml";
  testSensitivities(inputFile, false, false, out, success);
}

TEUCHOS_UNIT_TEST( Piro, ForwardOperatorSensitivities )
{
  std::string inputFile="input_Solve_NOX_1.xml";
  testSensitivities(inputFile, true, false, out, success);
}

TEUCHOS_UNIT_TEST( Piro, AdjointOperatorSensitivities )
{
  std::string inputFile="input_Solve_NOX_Adjoint.xml";
  testSensitivities(inputFile, false, true, out, success);
}
#endif

TEUCHOS_UNIT_TEST( Piro, Basic )
{
  int i1 = 5;
  TEST_EQUALITY_CONST( i1, 5 );
}


TEUCHOS_UNIT_TEST( Piro, Assignment )
{
  int i1 = 4;
  int i2 = i1;
  TEST_EQUALITY( i2, i1 );
}

} // namespace
