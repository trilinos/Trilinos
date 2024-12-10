// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "MockModelEval_A.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "Piro_ConfigDefs.hpp"
#ifdef HAVE_PIRO_NOX
#include "Piro_Epetra_NOXSolver.hpp"
#endif
#include "Piro_Epetra_SolverFactory.hpp"


namespace {

#ifdef HAVE_PIRO_NOX

void setOStream(const Teuchos::RCP<Teuchos::FancyOStream>& out,
                Teuchos::ParameterList& params) 
{
  params.sublist("NOX").sublist("Direction").sublist("Newton").sublist("Stratimikos Linear Solver").sublist("NOX Stratimikos Options").set("Output Stream", out);
  params.sublist("NOX").sublist("Printing").set("Output Stream", out->getOStream());
} 

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
  Teuchos::updateParametersFromXmlFile(inputFile, piroParams.ptr());
  setOStream(rcp(&out,false), *piroParams);
  
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
  TEUCHOS_ASSERT(outArgs.Ng() >= 2); // Number of *vectors* of responses
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
