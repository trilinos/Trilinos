#include "EpetraModelEval2DSim.hpp"
#include "EpetraModelEval4DOpt.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_DampenedNewtonNonlinearSolver.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main( int argc, char* argv[] )
{

  using Teuchos::CommandLineProcessor;
  typedef Teuchos::RefCountPtr<Thyra::VectorBase<double> > VectorPtr;

  bool success = true;

  try {
  
    //
    // Get options from the command line
    //

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    clp.setDocString(
      "This example program solves a simple 2 x 2 set of nonlinear equations using a simple\n"
      "dampened Newton method.\n\n"

      "The equations that are solved are:\n\n"

      "  f[0] =       x[0]      + x[1]*x[1] - p[0];\n"
      "  f[1] = d * ( x[0]*x[0] - x[1]      - p[1] );\n\n"

      "The Jacobian for these equations is nonsingular for every point except x=(-0.5,0.5)\n"
      "and x=(0.5,-0.5)  You can cause the Jacobian to be singular at the solution by setting\n"
      "p[0]=x[0]+x[1]*x[1] and p[1] = x[0]*x[0]-x[1] for these values of x.\n\n"

      "The equations are solved using a simple dampended Newton method that uses a Armijo\n"
      "line search which is implemented in the general class Thyra::DampenedNewtonNonlinearsolver\n"
      "You can get different levels of detail about the Newton method by adjustingthe command-line\n"
      "option \"verb-level\" (see above)\n"
      );
    
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

    bool         use4DOpt    = false;
    bool         externalFactory = false;

    bool showSetInvalidArg = false;
    bool showGetInvalidArg = false;

    clp.setOption( "d", &d, "Model constant d" );
    clp.setOption( "p0", &p0, "Model constant p[0]" );
    clp.setOption( "p1", &p1, "Model constant p[1]" );
    clp.setOption( "x00", &x00, "Initial guess for x[0]" );
    clp.setOption( "x01", &x01, "Initial guess for x[1]" );
    clp.setOption( "verb-level", &verbLevel, numVerbLevels, verbLevelValues, verbLevelNames, "Verbosity level" );
    clp.setOption( "tol", &tol, "Nonlinear solve tolerance" );
    clp.setOption( "max-iters", &maxIters, "Maximum number of nonlinear iterations" );
    clp.setOption( "use-4D-opt", "use-2D-sim", &use4DOpt
                   ,"Determines if the EpetraModelEval4DOpt or EpetraModelEval2DSim subclasses are used"  );
    clp.setOption( "external-lowsf", "internal-lowsf", &externalFactory
                   ,"Determines of the Thyra::LinearOpWithSolveFactory is used externally or internally to the Thyra::EpetraModelEvaluator object"  );
    clp.setOption( "show-set-invalid-arg", "no-show-set-invalid-arg", &showSetInvalidArg
                   ,"Determines if an attempt is made to set an invalid/unsupported ModelEvaluator input argument"  );
    clp.setOption( "show-get-invalid-arg", "no-show-get-invalid-arg", &showGetInvalidArg
                   ,"Determines if an attempt is made to get an invalid/unsupported ModelEvaluator output argument (2DSim only)"  );
  
    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    Teuchos::RefCountPtr<Teuchos::FancyOStream>
      out = Teuchos::VerboseObjectBase::getDefaultOStream();

    *out << "\nCreating the nonlinear equations object ...\n";
    
    Teuchos::RefCountPtr<EpetraExt::ModelEvaluator> epetraModel;
    if(use4DOpt) {
      epetraModel = rcp(new EpetraModelEval4DOpt(0.0,0.0,p0,p1,d,x00,x01,p0,p1));
    }
    else {
      epetraModel = rcp(new EpetraModelEval2DSim(d,p0,p1,x00,x01,showGetInvalidArg));
    }

    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = rcp(new Thyra::AmesosLinearOpWithSolveFactory());

    Teuchos::RefCountPtr<Thyra::EpetraModelEvaluator>
      epetraThyraModel = rcp(new Thyra::EpetraModelEvaluator());
    
    Teuchos::RefCountPtr<Thyra::ModelEvaluator<double> > thyraModel;
    if(externalFactory) {
      epetraThyraModel->initialize(epetraModel,Teuchos::null);
      thyraModel = Teuchos::rcp(
        new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(
          epetraThyraModel
          ,lowsFactory
          )
        );
    }
    else {
      epetraThyraModel->initialize(epetraModel,lowsFactory);
      thyraModel = epetraThyraModel;
    }
    
    if( showSetInvalidArg ) {
      *out << "\nAttempting to set an invalid input argument that throws an exception ...\n\n";
      Thyra::ModelEvaluatorBase::InArgs<double> inArgs = thyraModel->createInArgs();
      inArgs.set_x_dot(createMember(thyraModel->get_x_space()));
    }
    
    *out << "\nCreating the nonlinear solver and solving the equations ...\n\n";

    Thyra::DampenedNewtonNonlinearSolver<double> newtonSolver; // Set defaults
    newtonSolver.setVerbLevel(verbLevel);
    
    VectorPtr x = createMember(thyraModel->get_x_space());
    V_V( &*x, *thyraModel->getNominalValues().get_x() );

    Thyra::SolveCriteria<double> solveCriteria; // Sets defaults
    solveCriteria.solveMeasureType.set(Thyra::SOLVE_MEASURE_NORM_RESIDUAL,Thyra::SOLVE_MEASURE_NORM_RHS);
    solveCriteria.requestedTol = tol;
    solveCriteria.extraParameters = Teuchos::rcp(new Teuchos::ParameterList("Nonlinear Solve"));
    solveCriteria.extraParameters->set("Max Iters",int(maxIters));

    newtonSolver.setModel(thyraModel);
    Thyra::SolveStatus<double>
      solveStatus = newtonSolver.solve(&*x,&solveCriteria);
    
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
