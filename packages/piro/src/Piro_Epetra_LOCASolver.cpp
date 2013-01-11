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

#include "Piro_Epetra_LOCASolver.hpp"
#include "Piro_Epetra_MatrixFreeDecorator.hpp"
#include "Piro_ValidPiroParameters.hpp"
#include "NOX_Epetra_LinearSystem_Stratimikos.H"


Piro::Epetra::LOCASolver::LOCASolver(Teuchos::RCP<Teuchos::ParameterList> piroParams_,
                          Teuchos::RCP<EpetraExt::ModelEvaluator> model_,
                          Teuchos::RCP<NOX::Epetra::Observer> observer_,
                          Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigData_,
                          Teuchos::RCP<LOCA::StatusTest::Abstract> locaStatusTest_
) :
  piroParams(piroParams_),
  model(model_),
  observer(observer_),
  saveEigData(saveEigData_),
  locaStatusTest(locaStatusTest_),
  utils(piroParams->sublist("NOX").sublist("Printing"))
{
  //piroParams->validateParameters(*Piro::getValidPiroParameters(),0);

  Teuchos::ParameterList& noxParams = piroParams->sublist("NOX");
  Teuchos::ParameterList& printParams = noxParams.sublist("Printing");

  string jacobianSource = piroParams->get("Jacobian Operator", "Have Jacobian");
  bool leanMatrixFree = piroParams->get("Lean Matrix Free",false);

  Teuchos::ParameterList& noxstratlsParams = noxParams.
        sublist("Direction").sublist("Newton").sublist("Stratimikos Linear Solver");

  // Inexact Newton must be set in a second sublist when using 
  // Stratimikos: This code snippet sets it automatically
  bool inexact = (noxParams.sublist("Direction").sublist("Newton").
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
  //pVector = interface->getLOCAParameterVector();
  pVector = Teuchos::rcp(new LOCA::ParameterVector(interface->getLOCAParameterVector()));

  // Create the Jacobian matrix
  Teuchos::RCP<Epetra_Operator> A;
  Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> itmp = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = itmp;
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
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::Epetra::NOXSolver " <<
                 "Invalid value for parameter \" Jacobian Operator\"= " <<
                  jacobianSource << std::endl);


  /* Check flag if separate memory is requested for shifted matrix.
   * This allows 2 Matrix versions of eigensolvers, which save lots
   * of matrix recomputations.  */
  Teuchos::RCP<Epetra_Operator> Ashift=A;
  bool separateMatrixMem = piroParams->get("LOCASolver: Create Second Matrix",false);
  if (separateMatrixMem) {
    Teuchos::RCP<Epetra_CrsMatrix> Acrs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A); 
    if (Acrs != Teuchos::null)
     // Ashift = Teuchos::rcp(new Epetra_CrsMatrix(*Acrs));
     Ashift = model->create_W();
    else TEUCHOS_TEST_FOR_EXCEPTION(Acrs == Teuchos::null, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::Epetra::LOCASolver " <<
                 "LOCASolver: Create Second Matrix was requested, but only implemented for CrsMatrix\n");
  }

  // Create separate preconditioner if the model supports it
  /* NOTE: How do we want to decide between using an
   * available preconditioner: (1) If the model supports
   * it, then we use it, or (2) if a parameter list says
   * User_Defined ?  [Below, logic is option (1).]
   */

  Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner> WPrec;
  if (outargs.supports(EpetraExt::ModelEvaluator::OUT_ARG_WPrec))
    WPrec = model->create_WPrec();

  // Create the linear system
  // also Build shifted linear system for eigensolver
  Teuchos::RCP<NOX::Epetra::LinearSystemStratimikos> shiftedLinSys;
  if (WPrec != Teuchos::null) {
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;

    linsys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(printParams,
                  noxstratlsParams, iJac, A, iPrec, WPrec->PrecOp,
                  *currentSolution, WPrec->isAlreadyInverted));
    shiftedLinSys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(printParams,
  		  noxstratlsParams, iJac, Ashift, iPrec, WPrec->PrecOp,
                  *currentSolution, WPrec->isAlreadyInverted));
  }
  else {
     linsys =
        Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(printParams,
                          noxstratlsParams, iJac, A, *currentSolution));
     shiftedLinSys =
        Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(printParams,
                          noxstratlsParams, iJac, Ashift, *currentSolution));
  }

  Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = interface;

  // Create the LOCA Group
  grp = Teuchos::rcp(new LOCA::Epetra::Group(globalData, printParams, iTime,
                                       *currentSolution, linsys,
                                        shiftedLinSys, *pVector));

  grp->setDerivUtils(interface);
  if (separateMatrixMem) grp->declareSeparateMatrixMemory();

  // Saves one resid calculation per solve, but not as safe
  if (leanMatrixFree) grp->disableLinearResidualComputation(true);
  
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
   
  // Create the solver
  stepper = Teuchos::rcp(new LOCA::Stepper(globalData, grp, locaStatusTest, noxStatusTests, piroParams));
}

Piro::Epetra::LOCASolver::~LOCASolver()
{
  LOCA::destroyGlobalData(globalData);

  // Release saveEigenData RCP
  Teuchos::ParameterList& eigParams =
      piroParams->sublist("LOCA").sublist("Stepper").sublist("Eigensolver");
  eigParams.set( eigParams.get
      ("User-Defined Save Eigen Data Name", "Piro Strategy"), Teuchos::null);
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
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::LOCASolver::get_g_map(int j) const
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

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::LOCASolver::get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::LOCASolver::get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
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
  utils.out() << "eval pVector   " << std::setprecision(10) << *pVector << endl;
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

    int KrylovIters = linsys->getLinearItersTotal() - totalKrylovIters;
    int lastSolveKrylovIters = linsys->getLinearItersLastSolve();
    int linSolves = linsys->getNumLinearSolves() - totalLinSolves;

    totalNewtonIters += NewtonIters;
    totalKrylovIters += KrylovIters;
    totalLinSolves += linSolves;
    stepNum++;

    utils.out() << "Convergence Stats: for step  #" << stepNum << " : NumLinSolves, Krylov, Kr/Solve; LastKrylov, LastTol: " 
	 << linSolves << "  " << KrylovIters << "  " 
	 << (double) KrylovIters / (double) linSolves << "  " 
         << lastSolveKrylovIters << " " <<  linsys->getAchievedTol() << endl;

    if (stepNum > 1)
     utils.out() << "Convergence Stats: running total: NumLinSolves, Krylov, Kr/Solve, Kr/Step: " 
           << totalLinSolves << "  " << totalKrylovIters << "  " 
           << (double) totalKrylovIters / (double) totalLinSolves 
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
	grp->applyJacobianInverseMultiVector(*piroParams, dfdp_nox, dxdp_nox);
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

Teuchos::RCP<const LOCA::Stepper>
Piro::Epetra::LOCASolver::
getLOCAStepper() const
{
  return stepper;
}
 
Teuchos::RCP<LOCA::Stepper>
Piro::Epetra::LOCASolver::
getLOCAStepperNonConst()
{
  return stepper;
}

Teuchos::RCP<const LOCA::GlobalData>
Piro::Epetra::LOCASolver::
getGlobalData() const
{
  return globalData; 
}

Teuchos::RCP<LOCA::GlobalData>
Piro::Epetra::LOCASolver::
getGlobalDataNonConst()
{
  return globalData; 
}
