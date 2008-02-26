#include "AppModelEval.H"

AppModelEval::AppModelEval(const int nelem) 

{
#ifdef HAVE_MPI
    Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

    int MyPID = Comm->MyPID();
    
    // Create mesh
    vector<double> x(nelem+1);
    double h = 1.0/nelem;
    for (unsigned int i=0; i<=nelem; i++)
      x[i] = h*i;

    //set up parameters, responses
    int numParameters= 2;
    double left_bc = 2.0;
    double right_bc = 0.1;

    Epetra_LocalMap p_map(numParameters, 0, *Comm);
    p_init = Teuchos::rcp(new Epetra_Vector(p_map), false);
    p_init->operator[](0)=left_bc;
    p_init->operator[](1)=right_bc;

    int numResponses = 3;
    g_map = Teuchos::rcp(new Epetra_LocalMap(numResponses, 0, *Comm));

    // Set up application parameters
    appParams = Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& problemParams = 
      appParams->sublist("Problem");
    //problemParams.set("Name", "Brusselator");
    //problemParams.set("alpha", p_init->operator[](0));
    //problemParams.set("beta",  p_init->operator[](1));
    //problemParams.set("D1",    p_init->operator[](2));
    //problemParams.set("D2",    p_init->operator[](3));

    problemParams.set("Name", "Heat Nonlinear Source");

    // Create application
    Teuchos::RCP<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(x, Comm, appParams, false));

    // Create initial guess
    Teuchos::RCP<const Epetra_Vector> u = app->getInitialSolution();
    currentSolution = Teuchos::rcp(new NOX::Epetra::Vector(*u));

    Teuchos::RefCountPtr< Teuchos::Array<std::string> > free_param_names =
      Teuchos::rcp(new Teuchos::Array<std::string>);
    free_param_names->push_back("Constant Node BC 1");
    free_param_names->push_back("Constant Node BC 2");
//    free_param_names->push_back("Brusselator Alpha");
//    free_param_names->push_back("Brusselator Beta");

    // Create model evaluator
    Teuchos::RCP<FEApp::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app, free_param_names));

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
    pVector = interface->getLOCAParameterVector();
    cout << "constr pVector   " << pVector << endl;

    // Set up NOX parameters
    Teuchos::RCP<Teuchos::ParameterList> noxParams =
      Teuchos::rcp(&(appParams->sublist("NOX")),false);

    // Set the nonlinear solver method
    noxParams->set("Nonlinear Solver", "Line Search Based");

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = noxParams->sublist("Printing");
    printParams.set("MyPID", MyPID); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    printParams.set("Output Information", 
		    NOX::Utils::OuterIteration + 
		    NOX::Utils::OuterIterationStatusTest + 
		    NOX::Utils::InnerIteration +
		   // NOX::Utils::Parameters + 
		    NOX::Utils::Details + 
		    NOX::Utils::LinearSolverDetails +
		    NOX::Utils::Warning + 
		    NOX::Utils::Error);

    // Sublist for line search 
    Teuchos::ParameterList& searchParams = noxParams->sublist("Line Search");
    searchParams.set("Method", "Full Step");

    // Sublist for direction
    Teuchos::ParameterList& dirParams = noxParams->sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");  
    lsParams.set("Max Iterations", 800);  
    lsParams.set("Tolerance", 1e-4); 
    lsParams.set("Output Frequency", 50);
    lsParams.set("Preconditioner", "Ifpack");

    // Create the Jacobian matrix
    Teuchos::RCP<Epetra_Operator> A = model->create_W(); 

    // Create the linear system
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> itmp = interface;
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = itmp;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
			lsParams, iReq, iJac, A, *currentSolution));

    // Create the Group
    Teuchos::RCP<NOX::Epetra::Group> grp =
      Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, *currentSolution, linsys)); 

    // Create the Solver convergence test
    Teuchos::RCP<NOX::StatusTest::NormF> wrms = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8, 
					   NOX::StatusTest::NormF::Unscaled));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(10));
    Teuchos::RCP<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(wrms);
    combo->addStatusTest(maxiters);

    // Create the solver
    solver = NOX::Solver::buildSolver(grp, combo, noxParams);
}

AppModelEval::~AppModelEval()
{
}

Teuchos::RCP<const Epetra_Map> AppModelEval::get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> AppModelEval::get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> AppModelEval::get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> AppModelEval::get_g_map(int j) const
{
  TEST_FOR_EXCEPTION(j != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  AppModelEval::get_g_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     j << std::endl);

  return g_map;
}

Teuchos::RCP<const Epetra_Vector> AppModelEval::get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> AppModelEval::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  return p_init;
}

EpetraExt::ModelEvaluator::InArgs AppModelEval::createInArgs() const
{
  //return underlyingME->createInArgs();
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(1);
//  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs AppModelEval::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 1);
  return outArgs;
}

void AppModelEval::evalModel( const InArgs& inArgs,
                                            const OutArgs& outArgs ) const
{

  // Parse InArgs
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  if (!p_in.get()) cout << "ERROR: AppModelEval requires p as inargs" << endl;

  pVector.setValue(0, p_in->operator[](0));
  pVector.setValue(1, p_in->operator[](1));
  cout << "eval pVector   " << pVector << endl;
  interface->setParameters(pVector);

//COPY THIS AS INITIAL GUESS IF EXISTS : need to set support for x
//  Teuchos::RCP<const Epetra_Vector> x_in = inArgs.get_x();

  // Parse OutArgs

  Teuchos::RCP<Epetra_Vector> g_out = outArgs.get_g(0); 
  g_out->PutScalar(0.0);

    // Solve
    solver->reset(*currentSolution);
    NOX::StatusTest::StatusType status = solver->solve();

    // Create printing utilities
    NOX::Utils utils(appParams->sublist("NOX").sublist("Printing"));

    if (status == NOX::StatusTest::Converged) 
      utils.out() << "Test Passed!" << endl;
    else {
	utils.out() << "Nonlinear solver failed to converge!" << endl;
    }

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    currentSolution->operator=(finalGroup.getX());

    const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    // Output the parameter list
    if (utils.isPrintType(NOX::Utils::Parameters)) {
      utils.out() << endl << "Final Parameters" << endl
		  << "****************" << endl;
      appParams->print(utils.out());
      utils.out() << endl;
    }

   //  finalSolution.Print(utils.out());

    // Compute responses: average, left flux, right flux,
    int endIndex = finalSolution.MyLength() - 1;
    int numUnks = finalSolution.GlobalLength();
    int MyPID = Comm->MyPID();
    double tmp;
    finalSolution.Norm1(&tmp);
    g_out->operator[](0) = tmp / numUnks;

    // Final flux at x=0: grid spacing = 1.0/(numUnks-1.0)
    if (MyPID == 0) tmp = (finalSolution[1] - finalSolution[0]) * (numUnks-1.0);
    else tmp = -1.0e8;
    Comm->MaxAll(&tmp, &(g_out->operator[](1)), 1);

    // Final flux at x=1: grid spacing = 1.0/(numUnks-1.0)
    if (MyPID == Comm->NumProc() - 1) tmp =
       (finalSolution[endIndex] - finalSolution[endIndex - 1]) * (numUnks - 1.0);
    else tmp = -1.0e8;
    Comm->MaxAll(&tmp, &(g_out->operator[](2)), 1);

    cout << finalSolution << endl;

}
