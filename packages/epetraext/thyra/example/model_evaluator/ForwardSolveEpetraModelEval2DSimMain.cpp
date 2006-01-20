#include "EpetraModelEval2DSim.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_DampenedNewtonNonlinearSolver.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main( int argc, char* argv[] )
{

	using Teuchos::CommandLineProcessor;
  typedef Teuchos::RefCountPtr<Thyra::VectorBase<double> > VectorPtr;

  Teuchos::GlobalMPISession  mpiSession(&argc,&argv);

  bool success = true;

	try {
	
		//
		// Get options from the command line
		//
		
		double       d           = 10.0;
    double       p0          = 2.0;
    double       p1          = 0.0;
    double       x00         = 0.0;
    double       x01         = 1.0;

    const int numVerbLevels = 6;
    const Teuchos::EVerbosityLevel verbLevelValues[numVerbLevels]
      = { Teuchos::VERB_DEFAULT, Teuchos::VERB_NONE, Teuchos::VERB_LOW
          ,Teuchos::VERB_MEDIUM, Teuchos::VERB_HIGH, Teuchos::VERB_EXTREME };
    const char* verbLevelNames[numVerbLevels]
      = {"default","none","low","medium","high","extreme"};
    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;

    double       tol         = 1e-10;
    int          maxIters    = 100;

    bool showSetInvalidArg = false;
    bool showGetInvalidArg = false;

		CommandLineProcessor  clp(false); // Don't throw exceptions

		clp.setOption( "d", &d );
		clp.setOption( "p0", &p0 );
		clp.setOption( "p1", &p1 );
		clp.setOption( "x00", &x00 );
		clp.setOption( "x01", &x01 );
    clp.setOption( "verb-level", &verbLevel, numVerbLevels, verbLevelValues, verbLevelNames );
    clp.setOption( "tol", &tol, "Nonlinear solve tolerance" );
    clp.setOption( "max-iters", &maxIters, "Maximum number of nonlinear iterations" );
    clp.setOption( "show-set-invalid-arg", "no-show-set-invalid-arg", &showSetInvalidArg );
    clp.setOption( "show-get-invalid-arg", "no-show-get-invalid-arg", &showGetInvalidArg );
	
		CommandLineProcessor::EParseCommandLineReturn
			parse_return = clp.parse(argc,argv,&std::cerr);

		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
			return parse_return;

    Teuchos::RefCountPtr<Teuchos::FancyOStream>
      out = Teuchos::VerboseObjectBase::getDefaultOStream();

    *out << "\nCreating the nonlinear equations object ...\n";
		
    EpetraModelEval2DSim epetraModel(d,p0,p1,x00,x01,showGetInvalidArg);

    Thyra::EpetraModelEvaluator thyraModel; // Sets default options!
    thyraModel.initialize(
      Teuchos::rcp(&epetraModel,false)
      ,Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory())
      );

    if( showSetInvalidArg ) {
      *out << "\nAttempting to set an invalid input argument that throws an exception ...\n\n";
      Thyra::ModelEvaluatorBase::InArgs<double> inArgs = thyraModel.createInArgs();
      inArgs.set_x_dot(createMember(thyraModel.get_x_space()));
    }
    
    *out << "\nCreating the nonlinear solver and solving the equations ...\n\n";

    Thyra::DampenedNewtonNonlinearSolver<double> newtonSolver; // Set defaults
    newtonSolver.setVerbLevel(verbLevel);
    
    VectorPtr x = createMember(thyraModel.get_x_space());
    V_V( &*x, *thyraModel.get_x_init() );

    Thyra::SolveCriteria<double> solveCriteria; // Sets defaults
    solveCriteria.solveTolType = Thyra::SOLVE_TOL_REL_RESIDUAL_NORM;
    solveCriteria.requestedTol = tol;
    solveCriteria.maxIterations = maxIters;

    Thyra::SolveStatus<double>
      solveStatus = newtonSolver.solve(thyraModel,&*x,&solveCriteria);

    *out << "\nNonlinear solver return status:\n";
    if(1) {
      Teuchos::OSTab tab(out);
      *out << solveStatus;
    }
    *out << "\nFinal solution for x=\n" << *x;

	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)
    
  return  success ? 0 : 1;
}
