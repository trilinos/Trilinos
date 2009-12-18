#include "Piro_Epetra_LOCASolver.hpp"
#include "Piro_ValidPiroParameters.hpp"


Piro::Epetra::LOCASolver::LOCASolver(Teuchos::RCP<Teuchos::ParameterList> appParams_,
                          Teuchos::RCP<EpetraExt::ModelEvaluator> model_,
                          Teuchos::RCP<Piro::Epetra::NOXObserver> observer_
) :
  appParams(appParams_),
  model(model_),
  observer(observer_),
  utils(appParams->sublist("NOX").sublist("Printing"))
{
  appParams->validateParameters(*Piro::getValidPiroParameters(),0);

  Teuchos::RCP<Teuchos::ParameterList> noxParams =
	Teuchos::rcp(&(appParams->sublist("NOX")),false);
  Teuchos::ParameterList& printParams = noxParams->sublist("Printing");
  Teuchos::ParameterList& lsParams = noxParams->
	sublist("Direction").sublist("Newton").sublist("Linear Solver");

  // Grab some modelEval stuff from underlying model
  EpetraExt::ModelEvaluator::InArgs inargs = model->createInArgs();
  num_p = inargs.Np();
  EpetraExt::ModelEvaluator::OutArgs outargs = model->createOutArgs();
  num_g = outargs.Ng();

  // Create initial guess
  Teuchos::RCP<const Epetra_Vector> u = model->get_x_init();
  currentSolution = Teuchos::rcp(new NOX::Epetra::Vector(*u));

  // Create Epetra factory
  Teuchos::RefCountPtr<LOCA::Abstract::Factory> epetraFactory =
	Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create global data object
  Teuchos::RefCountPtr<LOCA::GlobalData> globalData =
	LOCA::createGlobalData(appParams, epetraFactory);

  // Create LOCA interface
  interface = Teuchos::rcp(
	       new LOCA::Epetra::ModelEvaluatorInterface(globalData, model));

  // Get LOCA parameter vector
  //pVector = interface->getLOCAParameterVector();
  pVector = Teuchos::rcp(new LOCA::ParameterVector(interface->getLOCAParameterVector()));

  // Create the Jacobian matrix
  Teuchos::RCP<Epetra_Operator> A;
  Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> itmp = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = itmp;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac;
  string jacobianSource = lsParams.get("Jacobian Operator", "Have Jacobian");

  if (jacobianSource == "Have Jacobian") {
    A = model->create_W();
    iJac = interface;
  }
  else if (jacobianSource == "Matrix-Free") {
    A = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams,
                                     iReq, *currentSolution));
    iJac = Teuchos::rcp_dynamic_cast<NOX::Epetra::Interface::Jacobian>(A);
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
  Teuchos::RCP<Epetra_Operator> M;
  if (outargs.supports(EpetraExt::ModelEvaluator::OUT_ARG_M))
    M = model->create_M();

  // Create the linear system
  // also Build shifted linear system for eigensolver
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> shiftedLinSys;
  if (M != Teuchos::null) {
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;
    lsParams.set("Preconditioner", "User Defined");
    linsys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                  lsParams, iJac, A, iPrec, M, *currentSolution));
    shiftedLinSys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                                lsParams, iReq, iJac, A, *currentSolution));
  }
  else {
     linsys =
        Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                          lsParams, iReq, iJac, A, *currentSolution));
     shiftedLinSys =
        Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                          lsParams, iReq, iJac, A, *currentSolution));
  }

  Teuchos::RefCountPtr<LOCA::Epetra::Interface::TimeDependent> iTime =
    interface;

  // Create the LOCA Group
  grp = Teuchos::rcp(new LOCA::Epetra::Group(globalData, printParams, iTime,
                                       *currentSolution, linsys,
                                        shiftedLinSys, *pVector));
  grp->setDerivUtils(interface);
  
  // Create the Solver convergence test
  Teuchos::ParameterList& statusParams = noxParams->sublist("Status Tests");
  Teuchos::RCP<NOX::StatusTest::Generic> statusTests =
    NOX::StatusTest::buildStatusTests(statusParams,
                                      *(globalData->locaUtils));

  // Create the solver
  stepper = Teuchos::rcp(new LOCA::Stepper(globalData, grp, statusTests, appParams));
}

Piro::Epetra::LOCASolver::~LOCASolver()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::LOCASolver::get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::LOCASolver::get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::LOCASolver::get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::LOCASolver::get_g_map(int j) const
{
  TEST_FOR_EXCEPTION( (j>1 || j<0), Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  Piro::Epetra::NOXSolver::get_g_map() only " <<
                     " supports 2 response vectors.  Supplied index l = " <<
                     j << std::endl);

  if      (j < num_g) return model->get_g_map(j);
  else if (j == num_g) return model->get_x_map();
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::LOCASolver::get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::LOCASolver::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  return model->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::LOCASolver::createInArgs() const
{
  //return underlyingME->createInArgs();
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::LOCASolver::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(num_p, num_g+1);

  EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
  for (int i=0; i<num_g; i++)
    for (int j=0; j<num_p; j++)
      if (!model_outargs.supports(OUT_ARG_DgDp, i, j).none())
        outArgs.setSupports(OUT_ARG_DgDp, i, j,
                            DerivativeSupport(DERIV_MV_BY_COL));

  return outArgs;
}

void Piro::Epetra::LOCASolver::evalModel( const InArgs& inArgs,
                              const OutArgs& outArgs ) const
{
  // Parse InArgs
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  if (!p_in.get()) cout << "ERROR: Piro::Epetra::LOCASolver requires p as inargs" << endl;
  int numParameters = p_in->GlobalLength();

  for (int i=0; i< numParameters; i++) pVector->setValue(i, (*p_in)[i]);
  utils.out() << "eval pVector   " << setprecision(10) << *pVector << endl;
  interface->setParameters(*pVector);

  // Solve
  // AGS:: Need to reset solution?? solver->reset(*currentSolution);
  LOCA::Abstract::Iterator::IteratorStatus status = stepper->run();

  if (status ==  LOCA::Abstract::Iterator::Finished) 
    utils.out() << "Continuation Stepper Finished" << endl;
  else if (status ==  LOCA::Abstract::Iterator::NotFinished) 
    utils.out() << "Continuation Stepper did not reach final value." << endl;
  else {
    utils.out() << "Nonlinear solver failed to converge!" << endl;
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
    appParams->print(utils.out());
    utils.out() << endl;
  }

  if (observer != Teuchos::null)
    observer->observeSolution(*finalSolution);

  // Print stats
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

  for (int i=0; i<num_p; i++) {
    // p
    model_inargs.set_p(i, inArgs.get_p(i));
    
    // df/dp
    bool do_sens = false;
    for (int j=0; j<num_g; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none() && 
	  outArgs.get_DgDp(i,j).getMultiVector() != Teuchos::null) {
	do_sens = true;
	Teuchos::Array<int> p_indexes = 
	  outArgs.get_DgDp(i,j).getDerivativeMultiVector().getParamIndexes();
	TEST_FOR_EXCEPTION(p_indexes.size() > 0, 
			   Teuchos::Exceptions::InvalidParameter,
			   std::endl <<
			   "Piro::Epetra::LOCASolver::evalModel():  " <<
			   "Non-empty paramIndexes for dg/dp(" << i << "," <<
			   j << ") is not currently supported." << std::endl);
      }
    }
    if (do_sens) {
      // This code does not work with non-empty p_indexes.  The reason is
      // each p_indexes could theoretically be different for each g.
      // We would then need to make one df/dp for all the chosen p's
      // and then index into them properly below.  Note that the number of
      // columns in df/dp should be the number of chosen p's, not the total
      // number of p's.
      Teuchos::RCP<const Epetra_Map> p_map = model->get_p_map(i);
      Teuchos::RCP<Epetra_MultiVector> dfdp = 
	Teuchos::rcp(new Epetra_MultiVector(*(model->get_f_map()), 
					    p_map->NumGlobalElements()));
      // Teuchos::Array<int> p_indexes = 
      // 	outArgs.get_DgDp(i,0).getDerivativeMultiVector().getParamIndexes();
      // EpetraExt::ModelEvaluator::DerivativeMultiVector 
      // 	dmv_dfdp(dfdp, DERIV_MV_BY_COL, p_indexes);
      EpetraExt::ModelEvaluator::DerivativeMultiVector 
	dmv_dfdp(dfdp, DERIV_MV_BY_COL);
      model_outargs.set_DfDp(i,dmv_dfdp);
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
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none()) {
	dgdp_out = outArgs.get_DgDp(i,j).getMultiVector();
	if (dgdp_out != Teuchos::null)
	  do_sens = true;
      }
    }
    if (do_sens) {
      Teuchos::RCP<const Epetra_Map> g_map = model->get_g_map(j);
      Teuchos::RCP<Epetra_MultiVector> dgdx = 
	Teuchos::rcp(new Epetra_MultiVector(finalSolution->Map(),
					    g_map->NumGlobalElements()));
      model_outargs.set_DgDx(j,dgdx);

      for (int i=0; i<num_p; i++) {
	// dg/dp
	if (!outArgs.supports(OUT_ARG_DgDp, i, j).none()) {
	  Teuchos::RCP<Epetra_MultiVector> dgdp_out = 
	    outArgs.get_DgDp(i,j).getMultiVector();
	  if (dgdp_out != Teuchos::null) {
	    dgdp_out->PutScalar(0.0);
	    Teuchos::Array<int> p_indexes = 
	      outArgs.get_DgDp(i,j).getDerivativeMultiVector().getParamIndexes();
	    EpetraExt::ModelEvaluator::DerivativeMultiVector 
	      dmv_dgdp(dgdp_out, DERIV_MV_BY_COL,p_indexes);
	    model_outargs.set_DgDp(i,j,dmv_dgdp);
	  }
	}
      }
    }
  }
    
  // (1) Calculate g, df/dp, dg/dp, dg/dx
  model->evalModel(model_inargs, model_outargs);
    
  for (int i=0; i<num_p; i++) {
    if (!model_outargs.supports(OUT_ARG_DfDp, i).none()) {
      Teuchos::RCP<Epetra_MultiVector> dfdp  = 
	model_outargs.get_DfDp(i).getMultiVector();
      if (dfdp != Teuchos::null) {
	int numParameters = dfdp->NumVectors();
	
	// (2) Calculate dx/dp multivector from -(J^{-1}*df/dp)
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
	
	if (observer != Teuchos::null &&
	    appParams->sublist("VTK").get("Visualize Sensitivities", false) ==
	    true) {
	  for (int k=0; k<numParameters; k++) {
	    Epetra_Vector* sv = dxdp_nox.getEpetraMultiVector()(k);
	    observer->observeSolution(*sv);
	  }
	}
	
	// (3) Calculate dg/dp = dg/dx*dx/dp + dg/dp
	// This may be the transpose of what we want since we specified
	// we want dg/dp by column in createOutArgs(). 
	// In this case just interchange the order of dgdx and dxdp
	// We should really probably check what the underlying ME does
	for (int j=0; j<num_g; j++) {
	  if (!outArgs.supports(OUT_ARG_DgDp, i, j).none()) {
	    Teuchos::RCP<Epetra_MultiVector> dgdp_out = 
	      outArgs.get_DgDp(i,j).getMultiVector();
	    if (dgdp_out != Teuchos::null) {
	      Teuchos::RCP<Epetra_MultiVector> dgdx = 
		model_outargs.get_DgDx(j).getMultiVector();
	      dgdp_out->Multiply('T', 'N', 1.0, *dgdx, *dxdp, 1.0);
	    }
	  }
	}
      }
    }
  }

  // return the final solution as an additional g-vector, if requested
  Teuchos::RCP<Epetra_Vector> gx_out = outArgs.get_g(num_g); 
  if (gx_out != Teuchos::null)  *gx_out = *finalSolution; 
}

