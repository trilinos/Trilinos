/*
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
*/

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "MockModelEval_A.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "Piro_ConfigDefs.hpp"
#ifdef Piro_ENABLE_NOX
#include "Piro_Epetra_NOXSolver.hpp"
#ifdef Piro_ENABLE_Stokhos
#include "Stokhos_Epetra.hpp"
#include "Piro_Epetra_StokhosSolverFactory.hpp"
#include "MockModelEval_C.hpp"

#include "Piro_Epetra_StokhosSolver.hpp"
#include "Piro_Epetra_NECoupledModelEvaluator.hpp"
#include "MockModelEval_D.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Piro_PerformAnalysis.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_DetachedVectorView.hpp"
#endif
#endif
#include "Piro_Epetra_SolverFactory.hpp"


namespace {

#ifdef Piro_ENABLE_NOX

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

int testResponses(const Epetra_Vector& g, 
		  const Teuchos::Array<double> testValues,
		  double absTol, double relTol,
		  const std::string& tag,
                  Teuchos::FancyOStream& out)
{
  int failures = 0;
  TEUCHOS_TEST_FOR_EXCEPTION(g.MyLength() != testValues.size(),
		     std::logic_error,
		     tag << " Test Values array has size " << 
		     testValues.size() << "but expected size " <<
		     g.MyLength());
  for (int i=0; i<testValues.size(); i++) {
    bool success = 
      std::abs(g[i]-testValues[i]) <= relTol*std::abs(testValues[i])+absTol;
    if (!success) 
      ++failures;
    out << tag << " test " << i;
    if (success)
      out << " passed";
    else
      out << " failed";
    out << ":  Expected:  " << testValues[i] << ", got:  " << g[i]
	  << ", abs tol = " << absTol << ", rel tol = " << relTol
	  << "." << std::endl;
  }

  return failures;
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

  RCP<Teuchos::FancyOStream> default_out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::VerboseObjectBase::setDefaultOStream(rcp(&out,false));

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
  setOStream(rcp(&out,false), *piroParams);
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

  Teuchos::VerboseObjectBase::setDefaultOStream(default_out);
}

TEUCHOS_UNIT_TEST( Piro, Coupled )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  RCP<Teuchos::FancyOStream> default_out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::VerboseObjectBase::setDefaultOStream(rcp(&out,false));

  // Create a communicator for Epetra objects
  RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
  globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  globalComm = rcp(new Epetra_SerialComm);
#endif

  std::string problem1_filename = "input_problem1.xml";
  std::string problem2_filename = "input_problem2.xml";
  std::string coupled_filename = "input_coupled.xml";

  // Setup problem 1
  RCP<ParameterList> piroParams1 = 
    Teuchos::getParametersFromXmlFile(problem1_filename);
  setOStream(rcp(&out,false), *piroParams1);
  RCP<EpetraExt::ModelEvaluator> model1 = rcp(new MockModelEval_D(globalComm));
  
  // Setup problem 2
  RCP<ParameterList> piroParams2 = 
    Teuchos::getParametersFromXmlFile(problem2_filename);
  setOStream(rcp(&out,false), *piroParams2);
  RCP<EpetraExt::ModelEvaluator> model2 = rcp(new MockModelEval_D(globalComm));
  
  // Setup coupled model
  RCP<ParameterList> coupledParams = 
    Teuchos::getParametersFromXmlFile(coupled_filename);
  setOStream(rcp(&out,false), *coupledParams);
  Teuchos::Array< RCP<EpetraExt::ModelEvaluator> > models(2);
  models[0] = model1; models[1] = model2;
  Teuchos::Array< RCP<ParameterList> > piroParams(2);
  piroParams[0] = piroParams1; piroParams[1] = piroParams2;
  RCP<Piro::Epetra::AbstractNetworkModel> network_model =
    rcp(new Piro::Epetra::ParamToResponseNetworkModel);
  RCP<Piro::Epetra::NECoupledModelEvaluator> coupledModel =
    rcp(new Piro::Epetra::NECoupledModelEvaluator(models, piroParams,
						  network_model,
						  coupledParams, globalComm));
  coupledModel->setOStream(rcp(&out,false));

  // Setup solver
  Piro::Epetra::SolverFactory solverFactory;
  RCP<EpetraExt::ModelEvaluator> coupledSolver =
    solverFactory.createSolver(coupledParams, coupledModel);
    
  // Solve coupled system
  EpetraExt::ModelEvaluator::InArgs inArgs = coupledSolver->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs outArgs = coupledSolver->createOutArgs();
  for (int i=0; i<inArgs.Np(); i++)
    inArgs.set_p(i, coupledSolver->get_p_init(i));
  for (int i=0; i<outArgs.Ng(); i++) {
    RCP<Epetra_Vector> g = 
      rcp(new Epetra_Vector(*(coupledSolver->get_g_map(i))));
    outArgs.set_g(i, g);
  }
  coupledSolver->evalModel(inArgs, outArgs);

  // Regression tests
  int failures = 0;
  Teuchos::ParameterList& testParams = 
    coupledParams->sublist("Regression Tests");
  double relTol = testParams.get("Relative Tolerance", 1.0e-3);
  double absTol = testParams.get("Absolute Tolerance", 1.0e-8);
  
  // Print results
  for (int i=0; i<outArgs.Ng(); i++) {
    RCP<Epetra_Vector> g = outArgs.get_g(i);
    if (g != Teuchos::null) {
      out << "Response vector " << i << ":" << std::endl;
      g->Print(out);
      
      // Test response
      std::stringstream ss1;
      ss1 << "Response " << i << " Test Values";
      bool testResponse = 
	testParams.isType< Teuchos::Array<double> >(ss1.str());
      if (testResponse) { 
	Teuchos::Array<double> testValues =
	  testParams.get<Teuchos::Array<double> >(ss1.str());
	failures += testResponses(*g, testValues, absTol, relTol, "Response", 
				  out);
      }

    }
  }

  success = failures == 0;
  Teuchos::VerboseObjectBase::setDefaultOStream(default_out);
}

TEUCHOS_UNIT_TEST( Piro, SGCoupled )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  RCP<Teuchos::FancyOStream> default_out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::VerboseObjectBase::setDefaultOStream(rcp(&out,false));

  // Create a communicator for Epetra objects
  RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
  globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  globalComm = rcp(new Epetra_SerialComm);
#endif

  std::string problem1_filename = "input_problem1_sg.xml";
  std::string problem2_filename = "input_problem2_sg.xml";
  std::string coupled_filename = "input_coupled_sg.xml";

  // Setup stochastic coupled problem to get spatial comm's
  RCP<ParameterList> coupledParams = 
    Teuchos::getParametersFromXmlFile(coupled_filename);
  setOStream(rcp(&out,false), *coupledParams);
  RCP<Piro::Epetra::StokhosSolver> coupledSolver =
    rcp(new Piro::Epetra::StokhosSolver(coupledParams, globalComm));
  RCP<const Epetra_Comm> app_comm = coupledSolver->getSpatialComm();

  // Setup problem 1
  RCP<ParameterList> piroParams1 = 
    Teuchos::getParametersFromXmlFile(problem1_filename);
  setOStream(rcp(&out,false), *piroParams1);
  RCP<EpetraExt::ModelEvaluator> model1 = rcp(new MockModelEval_D(app_comm));
  
  // Setup problem 2
  RCP<ParameterList> piroParams2 = 
    Teuchos::getParametersFromXmlFile(problem2_filename);
  setOStream(rcp(&out,false), *piroParams2);
  RCP<EpetraExt::ModelEvaluator> model2 = rcp(new MockModelEval_D(app_comm));
  
  // Setup coupled model
  Teuchos::Array< RCP<EpetraExt::ModelEvaluator> > models(2);
  models[0] = model1; models[1] = model2;
  Teuchos::Array< RCP<ParameterList> > piroParams(2);
  piroParams[0] = piroParams1; piroParams[1] = piroParams2;
  RCP<Piro::Epetra::AbstractNetworkModel> network_model =
    rcp(new Piro::Epetra::ParamToResponseNetworkModel);
  RCP<Piro::Epetra::NECoupledModelEvaluator> coupledModel =
    rcp(new Piro::Epetra::NECoupledModelEvaluator(models, piroParams,
						  network_model,
						  coupledParams, globalComm));
  coupledModel->setOStream(rcp(&out,false));

  // Setup solver
  coupledSolver->setup(coupledModel);
  
  Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly> x_sg_init =
    coupledSolver->get_x_sg_init();
  Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> x_sg_init_new =
    Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(*x_sg_init));
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis =
    coupledSolver->getBasis();
  for (int i=0; i<basis->dimension(); i++)
    (*x_sg_init_new)[i+1].PutScalar(1.0);
  coupledSolver->set_x_sg_init(*x_sg_init_new);
    
  // Solve coupled system
  EpetraExt::ModelEvaluator::InArgs inArgs = coupledSolver->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs outArgs = coupledSolver->createOutArgs();
  for (int i=0; i<inArgs.Np(); i++)
    if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_p_sg, i))
      inArgs.set_p_sg(i, coupledSolver->get_p_sg_init(i));
  for (int i=0; i<outArgs.Ng(); i++) 
    if (outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_g_sg, i)) {
      RCP<Stokhos::EpetraVectorOrthogPoly> g_sg = 
	coupledSolver->create_g_sg(i);
      outArgs.set_g_sg(i, g_sg);
    }
  coupledSolver->evalModel(inArgs, outArgs);

  // Regression tests
  int failures = 0;
  Teuchos::ParameterList& testParams = 
    coupledParams->sublist("Regression Tests");
  double relTol = testParams.get("Relative Tolerance", 1.0e-3);
  double absTol = testParams.get("Absolute Tolerance", 1.0e-8);

  
  // Print results
  for (int i=0; i<outArgs.Ng(); i++) {
    if (outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_g_sg, i)) {
      Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> g_sg = 
	outArgs.get_g_sg(i);
      if (g_sg != Teuchos::null) {
	Epetra_Vector g_mean(*(coupledSolver->get_g_map(i)));
	Epetra_Vector g_std_dev(*(coupledSolver->get_g_map(i)));
	g_sg->computeMean(g_mean);
	g_sg->computeStandardDeviation(g_std_dev);
	out.precision(12);
	out << "Response " << i << " Mean =      " << std::endl 
	    << g_mean << std::endl;
	out << "Response " << i << " Std. Dev. = " << std::endl 
	    << g_std_dev << std::endl;
	out << "Response vector " << i << ":" << std::endl
	    << *(outArgs.get_g_sg(i)) << std::endl;

       	// Test mean
	std::stringstream ss1;
	ss1 << "Response " << i << " Mean Test Values";
	bool testMean = 
	  testParams.isType< Teuchos::Array<double> >(ss1.str());
	if (testMean) { 
	  Teuchos::Array<double> testValues =
	    testParams.get<Teuchos::Array<double> >(ss1.str());
	  failures += testResponses(g_mean, testValues, absTol, relTol, "Mean", 
				    out);
	}

	// Test std. dev.
	std::stringstream ss2;
	ss2 << "Response " << i << " Standard Deviation Test Values";
	bool testSD = 
	  testParams.isType< Teuchos::Array<double> >(ss2.str());
	if (testSD) { 
	  Teuchos::Array<double> testValues =
	    testParams.get<Teuchos::Array<double> >(ss2.str());
	  failures += testResponses(g_std_dev, testValues, absTol, relTol, 
				    "Standard Deviation", out);
	}

      }
    }
  }

  success = failures == 0;
  Teuchos::VerboseObjectBase::setDefaultOStream(default_out);
}

#ifdef Piro_ENABLE_TriKota
TEUCHOS_UNIT_TEST( Piro, SGAnalysis )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  RCP<Teuchos::FancyOStream> default_out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::VerboseObjectBase::setDefaultOStream(rcp(&out,false));

  // Create a communicator for Epetra objects
  RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
    globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    globalComm = rcp(new Epetra_SerialComm);
#endif

  std::string xml_filename = "input_SGAnalysis.xml";

  // Set up application parameters
  RCP<ParameterList> appParams = 
    Teuchos::getParametersFromXmlFile(xml_filename);
    
  // Create stochastic Galerkin solver factory
  RCP<ParameterList> piroParams = 
    rcp(&(appParams->sublist("Piro")),false);
  setOStream(rcp(&out,false), *piroParams);
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

  Thyra::EpetraModelEvaluator rs_model_thyra;
  rs_model_thyra.initialize(rs_model, Teuchos::null);

  RCP< ::Thyra::VectorBase<double> > p;
  ParameterList& analysisParams = piroParams->sublist("Analysis");
  int status = Piro::PerformAnalysis(rs_model_thyra, analysisParams, p); 
  int failures = 0;
  if (status != 0)
    failures += 1000;
  
  out << "Analysis results = " << std::endl << *p;

  // Regression tests
  Teuchos::ParameterList& testParams = 
    appParams->sublist("Regression Tests");
  double relTol = testParams.get("Relative Tolerance", 1.0e-3);
  double absTol = testParams.get("Absolute Tolerance", 1.0e-8);
  if (testParams.isType< Teuchos::Array<double> >("Test Values")) { 
    Teuchos::Array<double> testValues =
      testParams.get<Teuchos::Array<double> >("Test Values");
    Thyra::DetachedVectorView<double> my_p(p);
    TEUCHOS_TEST_FOR_EXCEPTION(my_p.subDim() != testValues.size(),
		       std::logic_error,
		       "Test Values array has size " << 
		       testValues.size() << "but expected size " <<
		       my_p.subDim());
    for (int i=0; i<testValues.size(); i++) {
      bool success = 
	std::abs(my_p[i]-testValues[i]) <= relTol*std::abs(testValues[i])+absTol;
      if (!success) 
	++failures;
      out << " test " << i;
      if (success)
	out << " passed";
      else
	out << " failed";
      out << ":  Expected:  " << testValues[i] << ", got:  " << my_p[i]
	  << ", abs tol = " << absTol << ", rel tol = " << relTol
	  << "." << std::endl;
    }
  }
    
  success = failures == 0;
  Teuchos::VerboseObjectBase::setDefaultOStream(default_out);
}
#endif
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
