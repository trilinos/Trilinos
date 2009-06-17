#include "FEApp_NOXSolver.hpp"
#include "NOX_Epetra_MultiVector.H"

FEApp::NOXSolver::
NOXSolver(Teuchos::RCP<Teuchos::ParameterList> appParams_,
	  Teuchos::RCP<EpetraExt::ModelEvaluator> model_,
	  Teuchos::RCP<Epetra_Operator> M) :
  appParams(appParams_),
  model(model_),
  p_init(model_->get_p_init(0)),
  g_map(model_->get_g_map(0))
{
  Teuchos::RCP<Teuchos::ParameterList> noxParams =
    Teuchos::rcp(&(appParams->sublist("NOX")),false);
  Teuchos::ParameterList& printParams = noxParams->sublist("Printing");
  Teuchos::ParameterList& lsParams = noxParams->
    sublist("Direction").sublist("Newton").sublist("Linear Solver");
  
  // Create initial guess
  Teuchos::RCP<const Epetra_Vector> u = model->get_x_init();
  currentSolution = Teuchos::rcp(new NOX::Epetra::Vector(*u));

  // Create Epetra factory
  Teuchos::RefCountPtr<LOCA::Abstract::Factory> epetraFactory =
    Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create global data object
  globalData = LOCA::createGlobalData(appParams, epetraFactory);

  // Create LOCA interface
  interface = 
    Teuchos::rcp(new LOCA::Epetra::ModelEvaluatorInterface(globalData, model));

  // Get LOCA parameter vector
  //pVector = interface->getLOCAParameterVector();
  pVector = Teuchos::rcp(new LOCA::ParameterVector(interface->getLOCAParameterVector()));

  // Create the Jacobian matrix
  Teuchos::RCP<Epetra_Operator> A = model->create_W(); 

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> itmp = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = itmp;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys;
  if (M != Teuchos::null) {
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;
    lsParams.set("Preconditioner", "User Defined");
    linsys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
							lsParams,
							iJac, A, 
							iPrec, M,
							*currentSolution));
  }
  else
    linsys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, 
							iReq, iJac, A, 
							*currentSolution));

  // Build NOX group
  grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, 
					    *currentSolution, linsys));
  
  // Create the Solver convergence test
  Teuchos::ParameterList& statusParams = noxParams->sublist("Status Tests");
  Teuchos::RCP<NOX::StatusTest::Generic> statusTests =
    NOX::StatusTest::buildStatusTests(statusParams, *(globalData->locaUtils));

  // Create the solver
  solver = NOX::Solver::buildSolver(grp, statusTests, noxParams);
}

FEApp::NOXSolver::
~NOXSolver()
{
  LOCA::destroyGlobalData(globalData);
}

Teuchos::RCP<const Epetra_Map> FEApp::NOXSolver::
get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> FEApp::NOXSolver::
get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> FEApp::NOXSolver::
get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  FEApp::NOXSolver::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> FEApp::NOXSolver::
get_g_map(int j) const
{
  TEST_FOR_EXCEPTION(j != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  FEApp::NOXSolver::get_g_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     j << std::endl);
  return g_map;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
FEApp::NOXSolver::get_p_names(int l) const
{
  return model->get_p_names(l);
}

Teuchos::RCP<const Epetra_Vector> FEApp::NOXSolver::
get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> FEApp::NOXSolver::
get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  FEApp::NOXSolver::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_init;
}

EpetraExt::ModelEvaluator::InArgs FEApp::NOXSolver::
createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(1);

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs FEApp::NOXSolver::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 1);
  outArgs.setSupports(OUT_ARG_DgDp, 0, 0, 
		      DerivativeSupport(DERIV_MV_BY_COL));

  return outArgs;
}

void FEApp::NOXSolver::
evalModel(const InArgs& inArgs, const OutArgs& outArgs ) const
{

  // Parse InArgs
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  TEST_FOR_EXCEPTION(p_in == Teuchos::null, std::logic_error,
                     std::endl <<
                     "Error!  FEApp::NOXSolver::evalModel " <<
                     "requirs p as inargs! " << std::endl);
  int numParameters = p_in->GlobalLength();
  for (int i=0; i< numParameters; i++) 
    pVector->setValue(i, (*p_in)[i]);
  interface->setParameters(*pVector);

  // Parse OutArgs
  Teuchos::RCP<Epetra_Vector> g_out = outArgs.get_g(0); 

  // Solve
  solver->reset(*currentSolution);
  status = solver->solve();

  // Create printing utilities
  NOX::Utils utils(appParams->sublist("NOX").sublist("Printing"));

  if (status == NOX::StatusTest::Converged) 
    utils.out() << "Step Converged" << endl;
  else {
    utils.out() << "Nonlinear solver failed to converge!" << endl;
  }

  // Get the NOX and Epetra_Vector with the final solution from the solver
   (*currentSolution)=grp->getX();
   finalSolution = 
     Teuchos::rcp(&(currentSolution->getEpetraVector()), false);
  //finalSolution = 
  //  Teuchos::rcp(&(dynamic_cast<const NOX::Epetra::Vector&>(grp->getX()).getEpetraVector()), false);

  // Output the parameter list
  if (utils.isPrintType(NOX::Utils::Parameters)) {
    utils.out() << endl << "Final Parameters" << endl
		<< "****************" << endl;
    appParams->print(utils.out());
    utils.out() << endl;
  }

  {
    static int totalNewtonIters=0;
    static int totalKrylovIters=0;
    static int stepNum=0;
    int NewtonIters = appParams->sublist("NOX").
      sublist("Output").get("Nonlinear Iterations", -1000);
    int KrylovIters = appParams->sublist("NOX").sublist("Direction").sublist("Newton").
      sublist("Linear Solver").sublist("Output").
      get("Total Number of Linear Iterations", -1000);
    totalNewtonIters += NewtonIters;
    totalKrylovIters += KrylovIters;
    stepNum++;
    
    utils.out() << "Convergence Stats: for step  #" << stepNum 
		<< " : Newton, Krylov, Kr/Ne: " 
		<< NewtonIters << "  " << KrylovIters << "  " 
		<< (double) KrylovIters / (double) NewtonIters << std::endl;
    if (stepNum > 1)
      utils.out() << "Convergence Stats: running total: "
		  << "Newton, Krylov, Kr/Ne, Kr/Step: " 
		  << totalNewtonIters << "  " << totalKrylovIters << "  " 
		  << (double) totalKrylovIters / (double) totalNewtonIters 
		  << "  " << (double) totalKrylovIters / (double) stepNum 
		  << std::endl;  
  }

  // Do Sensitivity Calc, if requested. See 3 main steps 
  Teuchos::RCP<Epetra_MultiVector> dgdp_out;
  dgdp_out = outArgs.get_DgDp(0,0).getMultiVector();

  EpetraExt::ModelEvaluator::InArgs model_inargs = model->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
  model_inargs.set_x(finalSolution);
  model_inargs.set_p(0, p_in);

  if (g_out != Teuchos::null) {
    g_out->PutScalar(0.0);
    model_outargs.set_g(0, g_out);
  }

  Teuchos::RCP<Epetra_MultiVector> dfdp;
  Teuchos::RCP<Epetra_MultiVector> dgdx;
  if (dgdp_out != Teuchos::null) {
    dgdp_out->PutScalar(0.0);
    dfdp = Teuchos::rcp(new Epetra_MultiVector(*(model->get_f_map()), 
					       numParameters) );
    dgdx = Teuchos::rcp(new Epetra_MultiVector(finalSolution->Map(),
					       /*AGS:: Double Check this is right length: */ dgdp_out->GlobalLength()));
    Teuchos::Array<int> p_indexes = 
      outArgs.get_DgDp(0,0).getDerivativeMultiVector().getParamIndexes();
    EpetraExt::ModelEvaluator::DerivativeMultiVector dmv_dfdp(dfdp, 
							      DERIV_MV_BY_COL,
							      p_indexes);
    EpetraExt::ModelEvaluator::DerivativeMultiVector dmv_dgdp(dgdp_out, 
							      DERIV_MV_BY_COL,
							      p_indexes);
    model_outargs.set_DfDp(0,dmv_dfdp);
    model_outargs.set_DgDp(0,0,dmv_dgdp);
    model_outargs.set_DgDx(0,dgdx);
  }

  // (1) Calculate g, df/dp, dg/dp, dg/dx
  model->evalModel(model_inargs, model_outargs);
     
  if (dgdp_out != Teuchos::null) {
    // (2) Calculated dx/dp multivector from -(J^{-1}*df/dp)
    Teuchos::RCP<Epetra_MultiVector> dxdp = 
      Teuchos::rcp(new Epetra_MultiVector(*(model->get_x_map()), 
					  numParameters) );
    NOX::Epetra::MultiVector dfdp_nox(dfdp, NOX::DeepCopy,  
				      NOX::Epetra::MultiVector::CreateView);
    NOX::Epetra::MultiVector dxdp_nox(dxdp, NOX::DeepCopy,  
				      NOX::Epetra::MultiVector::CreateView);
    grp->computeJacobian();
    grp->applyJacobianInverseMultiVector(*appParams, dfdp_nox, dxdp_nox);
    dxdp_nox.scale(-1.0);

    // (3) Calculate dg/dp = dg/dx*dx/dp + dg/dp
    // This may be the transpose of what we want since we specified
    // we want dg/dp by column in createOutArgs(). 
    // In this case just interchange the order of dgdx and dxdp
    // We should really probably check what the underlying ME does
    dgdp_out->Multiply('T', 'N', 1.0, *dgdx, *dxdp, 1.0);
  }
}

Teuchos::RCP<const Epetra_Vector>
FEApp::NOXSolver::
getFinalSolution() const
{
  return finalSolution;
}

NOX::StatusTest::StatusType
FEApp::NOXSolver::
getSolverStatus() const
{
  return status;
}
