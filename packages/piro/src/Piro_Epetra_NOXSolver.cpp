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

#include "Piro_Epetra_NOXSolver.hpp"
#include "Piro_ValidPiroParameters.hpp"
#include "Piro_Epetra_MatrixFreeDecorator.hpp"
#include "NOX_Epetra_LinearSystem_Stratimikos.H"

#include "LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H"
#include "LOCA_Epetra_TransposeLinearSystem_Factory.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Epetra_Factory.H"

Piro::Epetra::NOXSolver::NOXSolver(
  Teuchos::RCP<Teuchos::ParameterList> piroParams_,
  Teuchos::RCP<EpetraExt::ModelEvaluator> model_,
  Teuchos::RCP<NOX::Epetra::Observer> observer_,
  Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> custom_interface,
  Teuchos::RCP<NOX::Epetra::LinearSystem> custom_linsys
) :
  piroParams(piroParams_),
  model(model_),
  observer(observer_),
  utils(piroParams->sublist("NOX").sublist("Printing"))
{
  //piroParams->validateParameters(*Piro::getValidPiroParameters(),0);

  Teuchos::RCP<Teuchos::ParameterList> noxParams =
	Teuchos::rcp(&(piroParams->sublist("NOX")),false);
  Teuchos::ParameterList& printParams = noxParams->sublist("Printing");

  string jacobianSource = piroParams->get("Jacobian Operator", "Have Jacobian");
  bool leanMatrixFree = piroParams->get("Lean Matrix Free",false);

  Teuchos::ParameterList& noxstratlsParams = noxParams->
        sublist("Direction").sublist("Newton").sublist("Stratimikos Linear Solver");

  // Inexact Newton must be set in a second sublist when using 
  // Stratimikos: This code snippet sets it automatically
  bool inexact = (noxParams->sublist("Direction").sublist("Newton").
                  get("Forcing Term Method", "Constant") != "Constant");
  noxstratlsParams.sublist("NOX Stratimikos Options").
                   set("Use Linear Solve Tolerance From NOX", inexact);


  if (jacobianSource == "Matrix-Free") {
    if (piroParams->isParameter("Matrix-Free Perturbation")) {
      model = Teuchos::rcp(new Piro::Epetra::MatrixFreeDecorator(model,
                           piroParams->get<double>("Matrix-Free Perturbation")));
    }
    else model = Teuchos::rcp(new Piro::Epetra::MatrixFreeDecorator(model));
  }

  // Grab some modelEval stuff from underlying model
  EpetraExt::ModelEvaluator::InArgs inargs = model->createInArgs();
  num_p = inargs.Np();
  EpetraExt::ModelEvaluator::OutArgs outargs = model->createOutArgs();
  num_g = outargs.Ng();

  // Create initial guess
  Teuchos::RCP<const Epetra_Vector> u = model->get_x_init();
  currentSolution = Teuchos::rcp(new NOX::Epetra::Vector(*u));

  // Create NOX interface from model evaluator
  if (custom_interface != Teuchos::null)
    interface = custom_interface;
  else
    interface = Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(model));
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;

  // Create the Jacobian matrix (unless flag is set to do it numerically)
  Teuchos::RCP<Epetra_Operator> A;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac;

  if (jacobianSource == "Have Jacobian" || jacobianSource == "Matrix-Free") {
    A = model->create_W();
    iJac = interface;
  }
  else if (jacobianSource == "Finite Difference") {
    A = Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams,
                                            iReq, *currentSolution));
    iJac = Teuchos::rcp_dynamic_cast<NOX::Epetra::Interface::Jacobian>(A);
  }
  else
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::Epetra::NOXSolver " <<
                 "Invalid value for parameter \" Jacobian Operator\"= " <<
                 jacobianSource << std::endl);

  // Create separate preconditioner if the model supports it
  /* NOTE: How do we want to decide between using an
   * available preconditioner: (1) If the model supports
   * it, then we use it, or (2) if a parameter list says
   * User_Defined ?  [Below, logic is ooption (1).]
   */
  Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner> WPrec;
  if (outargs.supports(EpetraExt::ModelEvaluator::OUT_ARG_WPrec))
    WPrec = model->create_WPrec(); 

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;
  if (custom_linsys != Teuchos::null)
    linsys = custom_linsys;
  else {
    if (WPrec != Teuchos::null) {
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;
      linsys = Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
			      printParams,
			      noxstratlsParams, iJac, A, iPrec, WPrec->PrecOp,
			      *currentSolution, WPrec->isAlreadyInverted));
    }
    else {
      linsys = Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
			      printParams,
			      noxstratlsParams, iJac, A, *currentSolution));
    }
  }

  // Build NOX group
  grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq,
                                           *currentSolution, linsys));

  // Saves one resid calculation per solve, but not as safe
  if (leanMatrixFree) grp->disableLinearResidualComputation(true);

  // Create the Solver convergence test
  Teuchos::ParameterList& statusParams = noxParams->sublist("Status Tests");
  Teuchos::RCP<NOX::StatusTest::Generic> statusTests =
    NOX::StatusTest::buildStatusTests(statusParams, utils);

  // Create the solver
  solver = NOX::Solver::buildSolver(grp, statusTests, noxParams);
}

Piro::Epetra::NOXSolver::~NOXSolver()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NOXSolver::get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NOXSolver::get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NOXSolver::get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::NOXSolver::get_p_map():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NOXSolver::get_g_map(int j) const
{
  TEST_FOR_EXCEPTION(j > num_g || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::NOXSolver::get_g_map():  " <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if      (j < num_g) return model->get_g_map(j);
  else if (j == num_g) return model->get_x_map();
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::NOXSolver::get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::NOXSolver::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::NOXSolver::get_p_init():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::NOXSolver::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::NOXSolver::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  // Ng is 1 bigger then model-Ng so that the solution vector can be an outarg
  outArgs.set_Np_Ng(num_p, num_g+1);

  EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
  for (int i=0; i<num_g; i++)
    for (int j=0; j<num_p; j++)
      outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			  model_outargs.supports(OUT_ARG_DgDp, i, j));

  return outArgs;
}

void Piro::Epetra::NOXSolver::evalModel(const InArgs& inArgs,
				const OutArgs& outArgs ) const
{
  // Parse input parameters
  for (int i=0; i<num_p; i++) {
    Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(i);
    if (p_in != Teuchos::null)
      interface->inargs_set_p(p_in, i); // Pass "p_in" through to inargs 
  }

  // Reset initial guess, if the user requests
  if(piroParams->sublist("NOX").get("Reset Initial Guess",false)==true)
    *currentSolution=*model->get_x_init();

  // Solve
  solver->reset(*currentSolution);
  NOX::StatusTest::StatusType status = solver->solve();

  // Print status
  if (status == NOX::StatusTest::Converged) 
    utils.out() << "Step Converged" << endl;
  else {
    utils.out() << "Nonlinear solver failed to converge!" << endl;
    outArgs.setFailed();
  }

  // Get the NOX and Epetra_Vector with the final solution from the solver
  (*currentSolution)=grp->getX();
  Teuchos::RCP<const Epetra_Vector> finalSolution = 
    Teuchos::rcp(&(currentSolution->getEpetraVector()), false);

  // Print solution
  if (utils.isPrintType(NOX::Utils::Details)) {
    utils.out() << endl << "Final Solution" << endl
		<< "****************" << endl;
    finalSolution->Print(std::cout);
  }

  // Output the parameter list
  if (utils.isPrintType(NOX::Utils::Parameters)) {
    utils.out() << endl << "Final Parameters" << endl
		<< "****************" << endl;
    piroParams->print(utils.out());
    utils.out() << endl;
  }

  // Print stats
  {
    static int totalNewtonIters=0;
    static int totalKrylovIters=0;
    static int stepNum=0;
    int NewtonIters = piroParams->sublist("NOX").
      sublist("Output").get("Nonlinear Iterations", -1000);
    int KrylovIters = piroParams->sublist("NOX").sublist("Direction").sublist("Newton").
      sublist("Linear Solver").sublist("Output").
      get("Total Number of Linear Iterations", -1000);
    totalNewtonIters += NewtonIters;
    totalKrylovIters += KrylovIters;
    stepNum++;
    
    utils.out() << "Convergence Stats: for step  #" << stepNum << " : Newton, Krylov, Kr/Ne: " 
	 << NewtonIters << "  " << KrylovIters << "  " 
	 << (double) KrylovIters / (double) NewtonIters << endl;
    if (stepNum > 1)
     utils.out() << "Convergence Stats: running total: Newton, Krylov, Kr/Ne, Kr/Step: " 
           << totalNewtonIters << "  " << totalKrylovIters << "  " 
           << (double) totalKrylovIters / (double) totalNewtonIters 
           << "  " << (double) totalKrylovIters / (double) stepNum << endl;
    
  }
    
  //
  // Do Sensitivity Calc, if requested. See 3 main steps 
  //

  // Set inargs and outargs
  EpetraExt::ModelEvaluator::InArgs model_inargs = model->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
  model_inargs.set_x(finalSolution);

  // We make different choices for layouts of df/dp, dg/dx depending on
  // whether we are doing forward or adjoint sensitivities
  std::string sensitivity_method = piroParams->get("Sensitivity Method",
						   "Forward");

  for (int i=0; i<num_p; i++) {
    // p
    model_inargs.set_p(i, inArgs.get_p(i));
    
    // df/dp
    bool do_sens = false;
    for (int j=0; j<num_g; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp, j, i).none() && 
	  !outArgs.get_DgDp(j,i).isEmpty()) {
	do_sens = true;
        
        // This code does not work with non-empty p_indexes.  The reason is
        // each p_indexes could theoretically be different for each g.
        // We would then need to make one df/dp for all the chosen p's
        // and then index into them properly below.  Note that the number of
        // columns in df/dp should be the number of chosen p's, not the total
        // number of p's.
	if (outArgs.get_DgDp(j,i).getMultiVector() != Teuchos::null) {
	  Teuchos::Array<int> p_indexes = 
	    outArgs.get_DgDp(j,i).getDerivativeMultiVector().getParamIndexes();
	  TEST_FOR_EXCEPTION(p_indexes.size() > 0, 
			     Teuchos::Exceptions::InvalidParameter,
			     std::endl <<
			     "Piro::Epetra::NOXSolver::evalModel():  " <<
			     "Non-empty paramIndexes for dg/dp(" << i << "," <<
			     j << ") is not currently supported." << std::endl);
	}
      }
    }
    if (do_sens) {
      Teuchos::RCP<const Epetra_Map> p_map = model->get_p_map(i);
      Teuchos::RCP<const Epetra_Map> f_map = model->get_f_map();
      int num_params = p_map->NumGlobalElements();
      int num_resids = f_map->NumGlobalElements();
      DerivativeSupport ds =  model_outargs.supports(OUT_ARG_DfDp,i);
      // Determine which layout to use for df/dp.  Ideally one would look
      // at num_params, num_resids, what is supported by the underlying
      // model evaluator, and the sensitivity method, and make the best 
      // choice to minimze the number of solves.  However this choice depends 
      // also on what layout of dg/dx is supported (e.g., if only the operator 
      // form is supported for forward sensitivities, then df/dp must be
      // DERIV_MV_BY_COL).  For simplicity, we order the conditional tests
      // to get the right layout in most situations.
      DerivativeLayout dfdp_layout;
      if (sensitivity_method == "Forward") {
        if (ds.supports(DERIV_MV_BY_COL))
          dfdp_layout = COL;
        else if (ds.supports(DERIV_TRANS_MV_BY_ROW))
          dfdp_layout = ROW;
	else if (ds.supports(DERIV_LINEAR_OP))
          dfdp_layout = OP;
        else
	  TEST_FOR_EXCEPTION(true, std::logic_error,
			     std::endl <<
			     "Piro::Epetra::NOXSolver::evalModel():  " <<
                             "Underlying ModelEvaluator must support " <<
                             "df/dp " << i << " with forward sensitivities."
                             << std::endl);
      }
      else if (sensitivity_method == "Adjoint") {
        if (ds.supports(DERIV_LINEAR_OP))
          dfdp_layout = OP;
	else if (ds.supports(DERIV_TRANS_MV_BY_ROW))
          dfdp_layout = ROW;
	else if (ds.supports(DERIV_MV_BY_COL))
          dfdp_layout = COL;
        else
	  TEST_FOR_EXCEPTION(true, std::logic_error,
			     std::endl <<
			     "Piro::Epetra::NOXSolver::evalModel():  " <<
                             "Underlying ModelEvaluator must support " <<
                             "df/dp " << i << " with adjoint sensitivities."
                             << std::endl);
      }
      else
        TEST_FOR_EXCEPTION(true, 
       		           Teuchos::Exceptions::InvalidParameter,
		           std::endl <<
		           "Piro::Epetra::NOXSolver::evalModel():  " <<
		           "Unknown sensitivity method" << sensitivity_method <<
		           ".  Valid choices are \"Forward\" and \"Adjoint\"." 
                           << std::endl);

      if (dfdp_layout == COL) {
	Teuchos::RCP<Epetra_MultiVector> dfdp =
	  Teuchos::rcp(new Epetra_MultiVector(*f_map, num_params));
	// Teuchos::Array<int> p_indexes = 
	// 	outArgs.get_DgDp(i,0).getDerivativeMultiVector().getParamIndexes();
	// EpetraExt::ModelEvaluator::DerivativeMultiVector 
	// 	dmv_dfdp(dfdp, DERIV_MV_BY_COL, p_indexes);
	EpetraExt::ModelEvaluator::DerivativeMultiVector 
	  dmv_dfdp(dfdp, DERIV_MV_BY_COL);
	model_outargs.set_DfDp(i,dmv_dfdp);
      }
      else if (dfdp_layout == ROW) {
	Teuchos::RCP<Epetra_MultiVector> dfdp =
	  Teuchos::rcp(new Epetra_MultiVector(*p_map, num_resids));
	// Teuchos::Array<int> p_indexes = 
	// 	outArgs.get_DgDp(i,0).getDerivativeMultiVector().getParamIndexes();
	// EpetraExt::ModelEvaluator::DerivativeMultiVector 
	// 	dmv_dfdp(dfdp, DERIV_TRANS_MV_BY_ROW, p_indexes);
	EpetraExt::ModelEvaluator::DerivativeMultiVector 
	  dmv_dfdp(dfdp, DERIV_TRANS_MV_BY_ROW);
	model_outargs.set_DfDp(i,dmv_dfdp);
      }
      else if (dfdp_layout == OP) {
	Teuchos::RCP<Epetra_Operator> dfdp_op = model->create_DfDp_op(i);
	model_outargs.set_DfDp(i,dfdp_op);
      }
    }
  }

  for (int j=0; j<num_g; j++) {
    // g
    Teuchos::RCP<Epetra_Vector> g_out = outArgs.get_g(j);
    if (g_out != Teuchos::null) {
      g_out->PutScalar(0.0);
      model_outargs.set_g(j, g_out);
    }

    // dg/dx
    bool do_sens = false;
    for (int i=0; i<num_p; i++) {
      Teuchos::RCP<Epetra_MultiVector> dgdp_out;
      if (!outArgs.supports(OUT_ARG_DgDp, j, i).none() &&
	  !outArgs.get_DgDp(j,i).isEmpty()) {
	do_sens = true;
      }
    }
    if (do_sens) {
      Teuchos::RCP<const Epetra_Map> g_map = model->get_g_map(j);
      Teuchos::RCP<const Epetra_Map> x_map = model->get_x_map();
      int num_responses = g_map->NumGlobalElements();
      int num_solution = x_map->NumGlobalElements();
      DerivativeSupport ds =  model_outargs.supports(OUT_ARG_DgDx,j);
      DerivativeLayout dgdx_layout;
      if (sensitivity_method == "Forward") {
	if (ds.supports(DERIV_LINEAR_OP))
          dgdx_layout = OP;
        if (ds.supports(DERIV_MV_BY_COL))
          dgdx_layout = COL;
        else if (ds.supports(DERIV_TRANS_MV_BY_ROW))
          dgdx_layout = ROW;
	else
	  TEST_FOR_EXCEPTION(true, std::logic_error,
			     std::endl <<
			     "Piro::Epetra::NOXSolver::evalModel():  " <<
                             "Underlying ModelEvaluator must support " <<
                             "dg/dx " << j << " with forward sensitivities."
                             << std::endl);
      }
      else if (sensitivity_method == "Adjoint") {
	if (ds.supports(DERIV_TRANS_MV_BY_ROW))
          dgdx_layout = ROW;
	else if (ds.supports(DERIV_MV_BY_COL))
          dgdx_layout = COL;
	else if (ds.supports(DERIV_LINEAR_OP))
          dgdx_layout = OP;
        else
	  TEST_FOR_EXCEPTION(true, std::logic_error,
			     std::endl <<
			     "Piro::Epetra::NOXSolver::evalModel():  " <<
                             "Underlying ModelEvaluator must support " <<
                             "dg/dx " << j << " with forward sensitivities."
                             << std::endl);
      }
      else
        TEST_FOR_EXCEPTION(true, 
       		           Teuchos::Exceptions::InvalidParameter,
		           std::endl <<
		           "Piro::Epetra::NOXSolver::evalModel():  " <<
		           "Unknown sensitivity method" << sensitivity_method <<
		           ".  Valid choices are \"Forward\" and \"Adjoint\"." 
                           << std::endl);

      if (dgdx_layout == OP) {
	Teuchos::RCP<Epetra_Operator> dgdx_op = model->create_DgDx_op(j);
	model_outargs.set_DgDx(j,dgdx_op);
      }
      else if (dgdx_layout == ROW) {
	Teuchos::RCP<Epetra_MultiVector> dgdx = 
	  Teuchos::rcp(new Epetra_MultiVector(*x_map, num_responses));
	model_outargs.set_DgDx(j,dgdx);
      }
      else if (dgdx_layout == COL) {
	Teuchos::RCP<Epetra_MultiVector> dgdx = 
	  Teuchos::rcp(new Epetra_MultiVector(*g_map, num_solution));
	model_outargs.set_DgDx(j,dgdx);
      }

      // dg/dp
      for (int i=0; i<num_p; i++) {
	if (!outArgs.supports(OUT_ARG_DgDp,j,i).none()) {
	  Derivative dgdp = outArgs.get_DgDp(j,i);
	  if (dgdp.getLinearOp() != Teuchos::null)
	    model_outargs.set_DgDp(j,i,model->create_DgDp_op(j,i));
	  else
	    model_outargs.set_DgDp(j,i,outArgs.get_DgDp(j,i));
	}
      }
    }
  }
    
  // (1) Calculate g, df/dp, dg/dp, dg/dx
  model->evalModel(model_inargs, model_outargs);
    
  if (sensitivity_method == "Forward") {
    for (int i=0; i<num_p; i++) {
      if (!model_outargs.supports(OUT_ARG_DfDp, i).none()) {
	TEST_FOR_EXCEPTION(
	  model_outargs.get_DfDp(i).getLinearOp()!=Teuchos::null,
	  std::logic_error,
	  std::endl <<
	  "Piro::Epetra::NOXSolver::evalModel():  " <<
	  "Can\'t use df/dp operator " << i << " with forward sensitivities." <<
	  std::endl);
	Teuchos::RCP<Epetra_MultiVector> dfdp  = 
	  model_outargs.get_DfDp(i).getMultiVector();
	if (dfdp != Teuchos::null) {
	  int num_cols = dfdp->NumVectors();
	
	  // (2) Calculate dx/dp multivector from -(J^{-1}*df/dp)
	  Teuchos::RCP<Epetra_MultiVector> dxdp = 
	    Teuchos::rcp(new Epetra_MultiVector(dfdp->Map(), num_cols));
	  NOX::Epetra::MultiVector dfdp_nox(
	    dfdp, NOX::DeepCopy,  
	    NOX::Epetra::MultiVector::CreateView);
	  NOX::Epetra::MultiVector dxdp_nox(
	    dxdp, NOX::DeepCopy,  
	    NOX::Epetra::MultiVector::CreateView);
	  
	  grp->computeJacobian();
	  grp->applyJacobianInverseMultiVector(*piroParams, dfdp_nox, dxdp_nox);
	  dxdp_nox.scale(-1.0);
	
	  // (3) Calculate dg/dp = dg/dx*dx/dp + dg/dp
	  for (int j=0; j<num_g; j++) {
	    if (!outArgs.supports(OUT_ARG_DgDp, j, i).none()) {
	      if (outArgs.get_DgDp(j,i).getLinearOp() != Teuchos::null) {
		// Handle operator case
		TEST_FOR_EXCEPTION(true, std::logic_error,
				   "Handle operator case." << std::endl);
	      }
	      else {
		Teuchos::RCP<Epetra_MultiVector> dgdp_out = 
		  outArgs.get_DgDp(j,i).getMultiVector();
		if (dgdp_out != Teuchos::null) {
		  Derivative dgdx_dv = model_outargs.get_DgDx(j);
		  Teuchos::RCP<Epetra_Operator> dgdx_op = dgdx_dv.getLinearOp();
		  Teuchos::RCP<Epetra_MultiVector> dgdx = 
		    dgdx_dv.getMultiVector();
		  Derivative dfdp_dv = model_outargs.get_DfDp(i);
		  EDerivativeMultiVectorOrientation dgdp_orient =
		    outArgs.get_DgDp(j,i).getMultiVectorOrientation();
		  if (dgdx_op != Teuchos::null) {
		    bool transpose = false;
		    if (dgdp_orient == DERIV_TRANS_MV_BY_ROW)
		      transpose = true;
		    Epetra_MultiVector tmp(dgdx_op->OperatorRangeMap(),
					   dxdp->NumVectors());
		    dgdx_op->Apply(*dxdp, tmp);
		    if (transpose) {
		      TEST_FOR_EXCEPTION(
			dgdp_out->Map().DistributedGlobal(), 
			std::logic_error,
			std::endl <<
			"Piro::Epetra::NOXSolver::evalModel():  " <<
			"Can\'t handle special case:  " << 
			" dg/dx operator, " <<
			" transposed, distributed dg/dp. " << std::endl);
		      for (int j=0; j<dgdp_out->NumVectors(); j++)
			for (int i=0; i<dgdp_out->MyLength(); i++)
			  (*dgdp_out)[j][i] += tmp[i][j];
		    }
		    else
		      dgdp_out->Update(1.0, tmp, 1.0);
		  }
		  else {
		    Teuchos::RCP<Epetra_MultiVector> arg1, arg2;
		    EDerivativeMultiVectorOrientation dfdp_orient =
		      outArgs.get_DfDp(i).getMultiVectorOrientation();
		    EDerivativeMultiVectorOrientation dgdx_orient =
		      outArgs.get_DgDx(j).getMultiVectorOrientation();
		    char flag1, flag2;
		    if (dgdp_orient == DERIV_MV_BY_COL) {
		      arg1 = dgdx;
		      arg2 = dxdp;
		      if (dgdx_orient == DERIV_MV_BY_COL)
			flag1 = 'N';
		      else
			flag1 = 'T';
		      if (dfdp_orient == DERIV_MV_BY_COL)
			flag2 = 'N';
		      else
			flag2 = 'T';
		    }
		    else {
		      arg1 = dxdp;
		      arg2 = dgdx;
		      if (dfdp_orient == DERIV_MV_BY_COL)
			flag1 = 'T';
		      else
			flag1 = 'N';
		      if (dgdx_orient == DERIV_MV_BY_COL)
			flag2 = 'T';
		      else
			flag2 = 'N';
		    }
		    dgdp_out->Multiply(flag1, flag2, 1.0, *arg1, *arg2, 1.0);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  else if (sensitivity_method == "Adjoint") {

    // Create transpose solver
    Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = grp->getLinearSystem();
    Teuchos::RCP<Epetra_Operator> jac = linSys->getJacobianOperator();
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
      Teuchos::rcp(new LOCA::Epetra::Factory);
    Teuchos::RCP<LOCA::GlobalData> globalData = 
      LOCA::createGlobalData(piroParams, epetraFactory);
    LOCA::Epetra::TransposeLinearSystem::Factory tls_factory(globalData);
    Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy = tls_factory.create(piroParams, linSys);
    grp->computeJacobian();
    tls_strategy->createJacobianTranspose();
    NOX::Epetra::Vector x_nox(
      Teuchos::rcp_const_cast<Epetra_Vector>(finalSolution), 
      NOX::Epetra::Vector::CreateView,
      NOX::DeepCopy);
    tls_strategy->createTransposePreconditioner(x_nox, *piroParams);

    for (int j=0; j<num_g; j++) {
      if (!model_outargs.supports(OUT_ARG_DgDx, j).none()) {
	TEST_FOR_EXCEPTION(
	  model_outargs.get_DgDx(j).getLinearOp()!=Teuchos::null,
	  std::logic_error,
	  std::endl <<
	  "Piro::Epetra::NOXSolver::evalModel():  " <<
	  "Can\'t use dg/dx operator " << j << " with adjoint sensitivities." <<
	  std::endl);
	Teuchos::RCP<Epetra_MultiVector> dgdx  = 
	  model_outargs.get_DgDx(j).getMultiVector();
	if (dgdx != Teuchos::null) {
	  int num_cols = dgdx->NumVectors();
	
	  // (2) Calculate xbar multivector from -(J^{-T}*dg/dx)
	  Teuchos::RCP<Epetra_MultiVector> xbar = 
	    Teuchos::rcp(new Epetra_MultiVector(dgdx->Map(), num_cols));
	  for (int col=0; col<num_cols; col++) {
	    Teuchos::RCP<Epetra_Vector> gx =
	      Teuchos::rcp((*dgdx)(col),false);
	    Teuchos::RCP<Epetra_Vector> xb =
	      Teuchos::rcp((*xbar)(col),false);
	    NOX::Epetra::Vector dgdx_nox(
	      gx, NOX::Epetra::Vector::CreateView, NOX::DeepCopy);
	    NOX::Epetra::Vector xbar_nox(
	      xb, NOX::Epetra::Vector::CreateView, NOX::DeepCopy);
	  
	    // Solve
	    tls_strategy->applyJacobianTransposeInverse(
	      *piroParams, dgdx_nox, xbar_nox);
	  }
	  xbar->Scale(-1.0);
	
	  // (3) Calculate dg/dp^T = df/dp^T*xbar + dg/dp^T
	  for (int i=0; i<num_p; i++) {
	    if (!outArgs.supports(OUT_ARG_DgDp, j, i).none()) {
	      if (outArgs.get_DgDp(j,i).getLinearOp() != Teuchos::null) {
		// Handle operator case
		TEST_FOR_EXCEPTION(true, std::logic_error,
				   "Handle operator case." << std::endl);
	      }
	      else {
		Teuchos::RCP<Epetra_MultiVector> dgdp_out = 
		  outArgs.get_DgDp(j,i).getMultiVector();
		if (dgdp_out != Teuchos::null) {
		  Derivative dfdp_dv = model_outargs.get_DfDp(i);
		  Teuchos::RCP<Epetra_Operator> dfdp_op = dfdp_dv.getLinearOp();
		  Teuchos::RCP<Epetra_MultiVector> dfdp = 
		    dfdp_dv.getMultiVector();
		  Derivative dgdx_dv = model_outargs.get_DgDx(j);
		  EDerivativeMultiVectorOrientation dgdp_orient =
		    outArgs.get_DgDp(j,i).getMultiVectorOrientation();
		  if (dfdp_op != Teuchos::null) {
		    bool transpose = false;
		    if (dgdp_orient == DERIV_MV_BY_COL)
		      transpose = true;
		    Epetra_MultiVector tmp(dfdp_op->OperatorDomainMap(),
					   xbar->NumVectors());
		    dfdp_op->SetUseTranspose(true);
		    dfdp_op->Apply(*xbar, tmp);
		    dfdp_op->SetUseTranspose(false);
		    if (transpose) {
		      TEST_FOR_EXCEPTION(
			dgdp_out->Map().DistributedGlobal(), 
			std::logic_error,
			std::endl <<
			"Piro::Epetra::NOXSolver::evalModel():  " <<
			"Can\'t handle special case:  " << 
			" df/dp operator, " <<
			" transposed, distributed dg/dp. " << std::endl);
		      for (int j=0; j<dgdp_out->NumVectors(); j++)
			for (int i=0; i<dgdp_out->MyLength(); i++)
			  (*dgdp_out)[j][i] += tmp[i][j];
		    }
		    else
		      dgdp_out->Update(1.0, tmp, 1.0);
		  }
		  else {
		    Teuchos::RCP<Epetra_MultiVector> arg1, arg2;
		    EDerivativeMultiVectorOrientation dgdp_orient =
		      outArgs.get_DgDp(j,i).getMultiVectorOrientation();
		    EDerivativeMultiVectorOrientation dfdp_orient =
		      outArgs.get_DfDp(i).getMultiVectorOrientation();
		    EDerivativeMultiVectorOrientation dgdx_orient =
		      outArgs.get_DgDx(j).getMultiVectorOrientation();
		    char flag1, flag2;
		    if (dgdp_orient == DERIV_TRANS_MV_BY_ROW) {
		      arg1 = dfdp;
		      arg2 = xbar;
		      if (dfdp_orient == DERIV_TRANS_MV_BY_ROW)
			flag1 = 'N';
		      else
			flag1 = 'T';
		      if (dgdx_orient == DERIV_TRANS_MV_BY_ROW)
			flag2 = 'N';
		      else
			flag2 = 'T';
		    }
		    else {
		      arg1 = xbar;
		      arg2 = dfdp;
		      if (dgdx_orient == DERIV_TRANS_MV_BY_ROW)
			flag1 = 'T';
		      else
			flag1 = 'N';
		      if (dfdp_orient == DERIV_TRANS_MV_BY_ROW)
			flag2 = 'T';
		      else
			flag2 = 'N';
		      
		    }
		    dgdp_out->Multiply(flag1, flag2, 1.0, *arg1, *arg2, 1.0);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  

  if (status == NOX::StatusTest::Converged) 
    if (observer != Teuchos::null)
      observer->observeSolution(*finalSolution);

  // return the final solution as an additional g-vector, if requested
  Teuchos::RCP<Epetra_Vector> gx_out = outArgs.get_g(num_g); 
  if (gx_out != Teuchos::null)  *gx_out = *finalSolution; 

  // Clear RCPs to parameter vectors
  for (int i=0; i<num_p; i++)
    interface->inargs_set_p(Teuchos::null, i); 
 
}
