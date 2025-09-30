// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_Epetra_LOCAAdaptiveSolver.hpp"
#include "Piro_Epetra_MatrixFreeDecorator.hpp"
#include "Piro_ValidPiroParameters.hpp"
#include "LOCA_Epetra_ModelEvaluatorInterface.H"


Piro::Epetra::LOCAAdaptiveSolver::LOCAAdaptiveSolver(const Teuchos::RCP<Teuchos::ParameterList>& piroParams_,
                          const Teuchos::RCP<EpetraExt::ModelEvaluator>& model_,
                          const Teuchos::RCP<Piro::Epetra::AdaptiveSolutionManager>& solnManager_,
                          Teuchos::RCP<NOX::Epetra::Observer> observer_,
                          Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData_
) :
  piroParams(piroParams_),
  model(model_),
  observer(observer_),
  saveEigData(saveEigData_),
  solnManager(solnManager_),
  utils(piroParams->sublist("NOX").sublist("Printing"))
{

  Teuchos::ParameterList& noxParams = piroParams->sublist("NOX");

  std::string jacobianSource = piroParams->get("Jacobian Operator", "Have Jacobian");

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

  // Create Epetra factory
  Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
	  Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create global data object
  globalData = LOCA::createGlobalData(piroParams, epetraFactory);

  // Create LOCA interface
  interface = Teuchos::rcp(
	       new LOCA::Epetra::ModelEvaluatorInterface(globalData, model));

  // explicitly set the observer for LOCA
  interface->setObserver( observer );


  // Get LOCA parameter vector
  pVector = Teuchos::rcp(new LOCA::ParameterVector(interface->getLOCAParameterVector()));

  // Create the Solver convergence test
  Teuchos::ParameterList& statusParams = noxParams.sublist("Status Tests");
  Teuchos::RCP<NOX::StatusTest::Generic> noxStatusTests =
    NOX::StatusTest::buildStatusTests(statusParams,
                                      *(globalData->locaUtils));

  // Register SaveEigenData strategy if given

  if (saveEigData != Teuchos::null) {
    Teuchos::ParameterList& eigParams =
      piroParams->sublist("LOCA").sublist("Stepper").sublist("Eigensolver");
    eigParams.set("Save Eigen Data Method", "User-Defined");
    eigParams.set("User-Defined Save Eigen Data Name", "Piro Strategy");
    eigParams.set( eigParams.get
      ("User-Defined Save Eigen Data Name", "Piro Strategy"), saveEigData);
  }

  solnManager->initialize(model, interface, pVector, globalData,
    outargs.supports(EpetraExt::ModelEvaluator::OUT_ARG_WPrec));


  // Create the stepper

  stepper = Teuchos::rcp(new LOCA::Epetra::AdaptiveStepper(piroParams, solnManager, globalData, noxStatusTests));

}

Piro::Epetra::LOCAAdaptiveSolver::~LOCAAdaptiveSolver()
{

  // Release solnManager RCPs
  solnManager->destroySolutionGroup();

  LOCA::destroyGlobalData(globalData);

  // Release saveEigenData RCP
  Teuchos::ParameterList& eigParams =
      piroParams->sublist("LOCA").sublist("Stepper").sublist("Eigensolver");
  eigParams.set( eigParams.get
      ("User-Defined Save Eigen Data Name", "Piro Strategy"), Teuchos::null);

}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::LOCAAdaptiveSolver::get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::LOCAAdaptiveSolver::get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::LOCAAdaptiveSolver::get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::LOCAAdaptiveSolver::get_g_map(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( (j>1 || j<0), Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  Piro::Epetra::NOXSolver::get_g_map() only " <<
                     " supports 2 response vectors.  Supplied index l = " <<
                     j << std::endl);

  if      (j < num_g) return model->get_g_map(j);
  else if (j == num_g) return model->get_x_map();
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::LOCAAdaptiveSolver::get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::LOCAAdaptiveSolver::get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  return model->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::LOCAAdaptiveSolver::createInArgs() const
{
  //return underlyingME->createInArgs();
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::LOCAAdaptiveSolver::createOutArgs() const
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

void Piro::Epetra::LOCAAdaptiveSolver::evalModel( const InArgs& inArgs,
                              const OutArgs& outArgs ) const
{
  // Parse InArgs
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  if (!p_in.get()) std::cout << "ERROR: Piro::Epetra::LOCASolver requires p as inargs" << std::endl;
  int numParameters = p_in->GlobalLength();

  for (int i=0; i< numParameters; i++) pVector->setValue(i, (*p_in)[i]);
  utils.out() << "eval pVector   " << std::setprecision(10) << *pVector << std::endl;
  interface->setParameters(*pVector);

  // Solve

  LOCA::Abstract::Iterator::IteratorStatus status = stepper->run();

  // Two acceptable outcomes: LOCA::Abstract::Iterator::Finished || LOCA::Abstract::Iterator::NotFinished

  // Checking the convergence of a nonlinear step
  if (status ==  LOCA::Abstract::Iterator::Finished)
    utils.out() << "Continuation Stepper Finished" << std::endl;
  else if (status ==  LOCA::Abstract::Iterator::NotFinished)
    utils.out() << "Continuation Stepper did not reach final value." << std::endl;
  else {
    utils.out() << "Nonlinear solver failed to converge!" << std::endl;
    outArgs.setFailed();
  }


  // Output the parameter list
/*
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
    globalData->locaUtils->out() << std::endl <<
      "### Final Parameters ############" << std::endl;
    stepper->getList()->print(globalData->locaUtils->out());
    globalData->locaUtils->out() << std::endl;
  }
*/

  // The time spent
  if (p_in->Map().Comm().MyPID() == 0)
    std::cout << std::endl << "#### Statistics ########" << std::endl;

  // Check number of steps
  int numSteps = stepper->getStepNumber();
  if (p_in->Map().Comm().MyPID() == 0)
    std::cout << " Number of continuation Steps = " << numSteps << std::endl;

  // Check number of failed steps
  int numFailedSteps = stepper->getNumFailedSteps();
  if (p_in->Map().Comm().MyPID() == 0)
    std::cout << " Number of failed continuation Steps = " << numFailedSteps << std::endl;

  std::cout << std::endl;

  // Get current solution from solution manager
  Teuchos::RCP<const Epetra_Vector> finalSolution = solnManager->updateSolution();

  // Print solution
  if (utils.isPrintType(NOX::Utils::Details)) {
    utils.out() << std::endl << "Final Solution" << std::endl
		<< "****************" << std::endl;
    finalSolution->Print(std::cout);
  }

  // Output the parameter list
  if (utils.isPrintType(NOX::Utils::Parameters)) {
    utils.out() << std::endl << "Final Parameters" << std::endl
		<< "****************" << std::endl;
    piroParams->print(utils.out());
    utils.out() << std::endl;
  }

  // Don't explicitly observe finalSolution:
  // This is already taken care of by the stepper which observes the solution after each
  // continuation step by default.

  // Print stats
  bool print_stats = piroParams->get("Print Convergence Stats", true);
  if (print_stats) {
    static int totalNewtonIters=0;
    static int totalKrylovIters=0;
    static int totalLinSolves = 0;
    static int stepNum=0;
    int NewtonIters = piroParams->sublist("NOX").
      sublist("Output").get("Nonlinear Iterations", -1000);

    int lastSolveKrylovIters;
    double achtol;
    int curKIters;
    int curLSolves;

    solnManager->getConvergenceData(curKIters, lastSolveKrylovIters, curLSolves, achtol);
    int KrylovIters = curKIters - totalKrylovIters;
    int linSolves = curLSolves - totalLinSolves;

    totalNewtonIters += NewtonIters;
    totalKrylovIters += KrylovIters;
    totalLinSolves += linSolves;
    stepNum++;

    utils.out() << "Convergence Stats: for step  #"
     << stepNum << " : NumLinSolves, Krylov, Kr/Solve; LastKrylov, LastTol: "
	 << linSolves << "  " << KrylovIters << "  "
	 << (double) KrylovIters / (double) linSolves << "  "
         << lastSolveKrylovIters << " " <<  achtol << std::endl;

    if (stepNum > 1)
     utils.out() << "Convergence Stats: running total: NumLinSolves, Krylov, Kr/Solve, Kr/Step: "
           << totalLinSolves << "  " << totalKrylovIters << "  "
           << (double) totalKrylovIters / (double) totalLinSolves
           << "  " << (double) totalKrylovIters / (double) stepNum << std::endl;

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
	TEUCHOS_TEST_FOR_EXCEPTION(p_indexes.size() > 0,
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
	int dfdp_numParameters = dfdp->NumVectors();

	// (2) Calculate dx/dp multivector from -(J^{-1}*df/dp)
	Teuchos::RCP<Epetra_MultiVector> dxdp =
	  Teuchos::rcp(new Epetra_MultiVector(*(model->get_x_map()),
					      dfdp_numParameters) );
	NOX::Epetra::MultiVector dfdp_nox(dfdp, NOX::DeepCopy,
					  NOX::Epetra::MultiVector::CreateView);
	NOX::Epetra::MultiVector dxdp_nox(dxdp, NOX::DeepCopy,
					  NOX::Epetra::MultiVector::CreateView);

	solnManager->applyJacobianInverseMultiVector(dfdp_nox, dxdp_nox);
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

  if (gx_out != Teuchos::null){

    // Resize gx if problem size has changed
    if(gx_out->MyLength() != finalSolution->MyLength()){

      std::cout << "Mismatch in storage sizes between final solution and response vector - not returning final solution."
                << std::endl;

    }
    else

      *gx_out = *finalSolution;

  }


}

