#include "ENAT_NOXSolver.hpp"
#include "NOX_Epetra_MultiVector.H"

ENAT::NOXSolver::NOXSolver(Teuchos::RCP<Teuchos::ParameterList> appParams_,
			   Teuchos::RCP<EpetraExt::ModelEvaluator> model_,
			   Teuchos::RCP<Epetra_Operator> M
) :
  appParams(appParams_),
  model(model_)
{

      Teuchos::RCP<Teuchos::ParameterList> noxParams =
	Teuchos::rcp(&(appParams->sublist("NOX")),false);
      Teuchos::ParameterList& printParams = noxParams->sublist("Printing");
      Teuchos::ParameterList& lsParams = noxParams->
	sublist("Direction").sublist("Newton").sublist("Linear Solver");

      // Create initial guess
      Teuchos::RCP<const Epetra_Vector> u = model->get_x_init();
      currentSolution = Teuchos::rcp(new NOX::Epetra::Vector(*u));

      // Create NOX interface from model evaluator
      interface = Teuchos::rcp(
	       new NOX::Epetra::ModelEvaluatorInterface(model));

      // Create the Jacobian matrix (unless flag is set to do it numerically)
      Teuchos::RCP<Epetra_Operator> A;
      if ( !(lsParams.isParameter("Jacobian Operator")) )
        A = model->create_W(); 

      // Create the linear system
      Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> itmp = interface;
      Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = itmp;
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
      Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys;
      if (A != Teuchos::null) {
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
        else {
	  linsys = 
	    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
							      lsParams, 
							      iReq, iJac, A, 
							      *currentSolution));
        }
      }
      else { //Matrix-Free Option
	linsys = 
	  Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
				        lsParams, iReq, *currentSolution));
      }

      // Build NOX group
      grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, *currentSolution, linsys));

      // Create the Solver convergence test
      NOX::Utils outputUtils(printParams);
      Teuchos::ParameterList& statusParams = noxParams->sublist("Status Tests");
      Teuchos::RCP<NOX::StatusTest::Generic> statusTests =
        NOX::StatusTest::buildStatusTests(statusParams, outputUtils);

      // Create the solver
      solver = NOX::Solver::buildSolver(grp, statusTests, noxParams);

      EpetraExt::ModelEvaluator::InArgs inargs = model->createInArgs();
      num_p = inargs.Np();
      EpetraExt::ModelEvaluator::OutArgs outargs = model->createOutArgs();
      num_g = outargs.Ng();
}

ENAT::NOXSolver::~NOXSolver()
{
}

Teuchos::RCP<const Epetra_Map> ENAT::NOXSolver::get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> ENAT::NOXSolver::get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> ENAT::NOXSolver::get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in ENAT::NOXSolver::get_p_map():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> ENAT::NOXSolver::get_g_map(int j) const
{
  TEST_FOR_EXCEPTION(j > num_g || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in ENAT::NOXSolver::get_g_map():  " <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if      (j < num_g) return model->get_g_map(j);
  else if (j == num_g) return model->get_x_map();
}

Teuchos::RCP<const Epetra_Vector> ENAT::NOXSolver::get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> ENAT::NOXSolver::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in ENAT::NOXSolver::get_p_init():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs ENAT::NOXSolver::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs ENAT::NOXSolver::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  // Ng is 1 bigger then model-Ng so that the solution vector can be an outarg
  outArgs.set_Np_Ng(num_p, num_g+1);

  EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
  for (int i=0; i<num_g; i++)
    for (int j=0; j<num_p; j++)
      if (!model_outargs.supports(OUT_ARG_DgDp, i, j).none())
	outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			    DerivativeSupport(DERIV_MV_BY_COL));

  return outArgs;
}

void ENAT::NOXSolver::evalModel(const InArgs& inArgs,
				const OutArgs& outArgs ) const
{
  // Parse input parameters
  for (int i=0; i<num_p; i++) {
    Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(i);
    if (p_in != Teuchos::null)
      interface->inargs_set_p(p_in, i); // Pass "p_in" through to inargs 
  }

  // Solve
  solver->reset(*currentSolution);
  NOX::StatusTest::StatusType status = solver->solve();

  // Print status
  NOX::Utils utils(appParams->sublist("NOX").sublist("Printing"));
  if (status == NOX::StatusTest::Converged) 
    utils.out() << "Step Converged" << endl;
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
    
    cout << "Convergence Stats: for step  #" << stepNum << " : Newton, Krylov, Kr/Ne: " 
	 << NewtonIters << "  " << KrylovIters << "  " 
	 << (double) KrylovIters / (double) NewtonIters << endl;
    if (stepNum > 1)
      cout << "Convergence Stats: running total: Newton, Krylov, Kr/Ne, Kr/Step: " 
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
			   "ENAT::NOXSolver::evalModel():  " <<
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

