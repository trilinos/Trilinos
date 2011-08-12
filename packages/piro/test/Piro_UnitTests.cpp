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
#ifdef Piro_ENABLE_Stokhos
#include "Stokhos_Epetra.hpp"
#include "Piro_Epetra_StokhosSolverFactory.hpp"
#include "MockModelEval_C.hpp"
#endif
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

#ifdef Piro_ENABLE_Stokhos
TEUCHOS_UNIT_TEST( Piro, SGResponseStatisticsSensitivity )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  // Create a communicator for Epetra objects
  RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
  globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  globalComm = rcp(new Epetra_SerialComm);
#endif

  std::string xml_filename = "input_SGSolve.xml";

  // Set up application parameters
  RCP<ParameterList> appParams = 
    Teuchos::getParametersFromXmlFile(xml_filename);
    
  // Create stochastic Galerkin solver factory
  RCP<ParameterList> piroParams = 
    rcp(&(appParams->sublist("Piro")),false);
  Piro::Epetra::StokhosSolverFactory sg_solver_factory(piroParams, 
						       globalComm);

  // Get comm for spatial problem
  RCP<const Epetra_Comm> app_comm = sg_solver_factory.getSpatialComm();
  
  // Create application model evaluator
  RCP<EpetraExt::ModelEvaluator> model = rcp(new MockModelEval_C(app_comm));
  
  // Setup rest of solver
  RCP<Stokhos::SGModelEvaluator> sg_model = 
    sg_solver_factory.createSGModel(model);
  RCP<EpetraExt::ModelEvaluator> sg_solver =
    sg_solver_factory.createSGSolver(sg_model);
  RCP<EpetraExt::ModelEvaluator> rs_model =
    sg_solver_factory.createRSModel(sg_solver);
  
  // Evaluate SG responses at SG parameters
  EpetraExt::ModelEvaluator::InArgs sg_inArgs = rs_model->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs sg_outArgs = 
    rs_model->createOutArgs();
  int p_index = 1; // PC expansion coefficients of params
  int g_index = 0;
  int num_g = 2;
  int x_index = num_g-1;
  int g_mean_index = g_index + num_g;
  int g_var_index = g_index + 2*num_g;
  RCP<const Epetra_Vector> p_init = rs_model->get_p_init(p_index);
  RCP<Epetra_Vector> g = 
    rcp(new Epetra_Vector(*(rs_model->get_g_map(g_index))));
  RCP<Epetra_Vector> x = 
    rcp(new Epetra_Vector(*(rs_model->get_g_map(x_index))));
  RCP<Epetra_Vector> g_mean = 
    rcp(new Epetra_Vector(*(rs_model->get_g_map(g_mean_index))));
  RCP<Epetra_Vector> g_var = 
    rcp(new Epetra_Vector(*(rs_model->get_g_map(g_var_index))));
  RCP<Epetra_MultiVector> dgdp_mean = 
    rcp(new Epetra_MultiVector(
	  *(rs_model->get_p_map(p_index)),
	  rs_model->get_g_map(g_mean_index)->NumMyElements()));
  RCP<Epetra_MultiVector> dgdp_var = 
    rcp(new Epetra_MultiVector(
	  *(rs_model->get_p_map(p_index)),
	  rs_model->get_g_map(g_var_index)->NumMyElements()));
  
  sg_outArgs.set_g(g_index, g);
  sg_outArgs.set_g(x_index, x);
  sg_outArgs.set_g(g_mean_index, g_mean);
  sg_outArgs.set_g(g_var_index, g_var);
  sg_outArgs.set_DgDp(
    g_mean_index, p_index, 
    EpetraExt::ModelEvaluator::Derivative(
      dgdp_mean, 
      EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW
      )
    );
  sg_outArgs.set_DgDp(
    g_var_index, p_index, 
    EpetraExt::ModelEvaluator::Derivative(
      dgdp_var, 
      EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW
      )
    );
  
  rs_model->evalModel(sg_inArgs, sg_outArgs);
  
  // Test derivatives with finite differences
  double delta = 1.0e-6;
  int num_p = rs_model->get_p_map(p_index)->NumMyElements();
  int num_resp = model->get_g_map(g_index)->NumMyElements();
  Teuchos::RCP<Epetra_Vector> p_pert = 
    Teuchos::rcp(new Epetra_Vector((*rs_model->get_p_map(p_index))));
  Teuchos::RCP<Epetra_Vector> g_mean_pert = 
    Teuchos::rcp(new Epetra_Vector(*(rs_model->get_g_map(g_mean_index))));
  Teuchos::RCP<Epetra_Vector> g_var_pert = 
    Teuchos::rcp(new Epetra_Vector(*(rs_model->get_g_map(g_var_index))));
  Teuchos::RCP<Epetra_MultiVector> dgdp_mean_fd = 
    Teuchos::rcp(new Epetra_MultiVector(*(rs_model->get_p_map(p_index)),
					num_resp));
  Teuchos::RCP<Epetra_MultiVector> dgdp_var_fd = 
    Teuchos::rcp(new Epetra_MultiVector(*(rs_model->get_p_map(p_index)),
					num_resp));
  EpetraExt::ModelEvaluator::InArgs sg_inArgs_pert = 
    rs_model->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs sg_outArgs_pert = 
    rs_model->createOutArgs();
  sg_inArgs_pert.set_p(p_index, p_pert);
  sg_outArgs_pert.set_g(g_mean_index, g_mean_pert);
  sg_outArgs_pert.set_g(g_var_index, g_var_pert);
  for (int i=0; i<num_p; i++) {
    
    // Perturb p
    //double h = delta*(std::abs((*p_init)[i])+delta);
    double h = delta;
    *p_pert = *p_init;
    (*p_pert)[i] += h;
    
    // Init perturbed g
    g_mean_pert->PutScalar(0.0);
    g_var_pert->PutScalar(0.0);
    
    // Compute perturbed g
    rs_model->evalModel(sg_inArgs_pert, sg_outArgs_pert);
    
    // Compute FD derivatives
    for (int j=0; j<num_resp; j++) {
      (*dgdp_mean_fd)[j][i] = ((*g_mean_pert)[j] - (*g_mean)[j])/h;
      (*dgdp_var_fd)[j][i] = ((*g_var_pert)[j] - (*g_var)[j])/h;
    }
  }

  out.precision(12);
  out << "Solution   = " << std::endl << *x << std::endl;
  out << "Response   = " << std::endl << *g << std::endl;
  out << "Mean       = " << (*g_mean)[0] << std::endl;
  out << "Variance   = " << (*g_var)[0] << std::endl;
  out << "d(Mean)/dp = " << std::endl << *dgdp_mean << std::endl;
  out << "d(Mean)/dp FD = " << std::endl << *dgdp_mean_fd << std::endl;
  out << "d(Var)/dp = " << std::endl << *dgdp_var << std::endl;
  out << "d(Var)/dp FD = " << std::endl << *dgdp_var_fd << std::endl;

  // Check analytic and FD sensitivities agree
  Epetra_MultiVector dgdp_mean_diff(*(rs_model->get_p_map(p_index)), num_resp);
  Epetra_MultiVector dgdp_var_diff(*(rs_model->get_p_map(p_index)), num_resp);
  dgdp_mean_diff.Update(1.0, *dgdp_mean, -1.0, *dgdp_mean_fd, 0.0);
  dgdp_var_diff.Update(1.0, *dgdp_var, -1.0, *dgdp_var_fd, 0.0);
  double nrm_mean, nrm_var;
  dgdp_mean_diff.NormInf(&nrm_mean);
  dgdp_var_diff.NormInf(&nrm_var);
  out << "Infinity norm of d(Mean)/dp error = " << nrm_mean << std::endl;
  out << "Infinity norm of d(Var)/dp error = " << nrm_mean << std::endl;
  double tol = 10*delta;
  success = (nrm_mean < tol) && (nrm_var < tol);
}
#endif
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
