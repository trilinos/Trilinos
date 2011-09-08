/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "EpetraModelEval2DSim.hpp"
#include "EpetraModelEval4DOpt.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_DampenedNewtonNonlinearSolver.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"


namespace {
  

const Teuchos::RCP<const Epetra_Vector>
createScalingVec(const double &scale, const Epetra_Map &map)
{
  Teuchos::RCP<Epetra_Vector> scalingVec = Teuchos::rcp(new Epetra_Vector(map));
  scalingVec->PutScalar(scale);
  return scalingVec;
}


void scaleEpetraModelEvaluator( const double &s_x, const double &s_f,
  const Teuchos::Ptr<Thyra::EpetraModelEvaluator> &model
  )
{
  if (s_x != 1.0) {
    model->setStateVariableScalingVec(
      createScalingVec(s_x, *model->getEpetraModel()->get_x_map())
      );
  }
  if (s_f != 1.0) {
    model->setStateFunctionScalingVec(
      createScalingVec(s_f, *model->getEpetraModel()->get_f_map())
      );
  }
}


} // namespace


int main( int argc, char* argv[] )
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::CommandLineProcessor;
  using Teuchos::outArg;
  typedef RCP<Thyra::VectorBase<double> > VectorPtr;

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

    double d = 10.0;
    clp.setOption( "d", &d, "Model constant d" );
    double p0 = 2.0;
    clp.setOption( "p0", &p0, "Model constant p[0]" );
    double p1 = 0.0;
    clp.setOption( "p1", &p1, "Model constant p[1]" );
    double x00 = 0.0;
    clp.setOption( "x00", &x00, "Initial guess for x[0]" );
    double x01 = 1.0;
    clp.setOption( "x01", &x01, "Initial guess for x[1]" );
    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
    setVerbosityLevelOption( "verb-level", &verbLevel, "Verbosity level.", &clp );
    double tol = 1e-10;
    clp.setOption( "tol", &tol, "Nonlinear solve tolerance" );
    int maxIters = 100;
    clp.setOption( "max-iters", &maxIters, "Maximum number of nonlinear iterations" );
    bool use4DOpt = false;
    clp.setOption( "use-4D-opt", "use-2D-sim", &use4DOpt,
      "Determines if the EpetraModelEval4DOpt or EpetraModelEval2DSim subclasses are used" );
    bool externalFactory = false;
    clp.setOption( "external-lowsf", "internal-lowsf", &externalFactory,
      "Determines of the Thyra::LinearOpWithSolveFactory is used externally or internally to the Thyra::EpetraModelEvaluator object" );
    bool showSetInvalidArg = false;
    clp.setOption( "show-set-invalid-arg", "no-show-set-invalid-arg", &showSetInvalidArg,
      "Determines if an attempt is made to set an invalid/unsupported ModelEvaluator input argument" );
    bool showGetInvalidArg = false;
    clp.setOption( "show-get-invalid-arg", "no-show-get-invalid-arg", &showGetInvalidArg,
      "Determines if an attempt is made to get an invalid/unsupported ModelEvaluator output argument (2DSim only)" );
    double s_x = 1.0;
    clp.setOption( "x-scale", &s_x, "State variables scaling." );
    double s_f = 1.0;
    clp.setOption( "f-scale", &s_f, "State function scaling." );
  
    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    RCP<Teuchos::FancyOStream>
      out = Teuchos::VerboseObjectBase::getDefaultOStream();

    *out << "\nCreating the nonlinear equations object ...\n";
    
    RCP<EpetraExt::ModelEvaluator> epetraModel;
    if(use4DOpt) {
      epetraModel = rcp(new EpetraModelEval4DOpt(0.0,0.0,p0,p1,d,x00,x01,p0,p1));
    }
    else {
      epetraModel = rcp(new EpetraModelEval2DSim(d,p0,p1,x00,x01,showGetInvalidArg));
    }

    *out << "\nCreating linear solver strategy ...\n";

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(Teuchos::parameterList());
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = linearSolverBuilder.createLinearSolveStrategy("Amesos");
    // Above, we are just using the stratimkikos class
    // DefaultLinearSolverBuilder to create a default Amesos solver
    // (typically KLU) with a default set of options.  By setting a parameter
    // list on linearSolverBuilder, you build from a number of solvers.

    RCP<Thyra::EpetraModelEvaluator>
      epetraThyraModel = rcp(new Thyra::EpetraModelEvaluator());
    
    RCP<Thyra::ModelEvaluator<double> > thyraModel;
    if(externalFactory) {
      epetraThyraModel->initialize(epetraModel, Teuchos::null);
      thyraModel = rcp(
        new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(
          epetraThyraModel, lowsFactory
          )
        );
    }
    else {
      epetraThyraModel->initialize(epetraModel, lowsFactory);
      thyraModel = epetraThyraModel;
    }

    scaleEpetraModelEvaluator( s_x, s_f, epetraThyraModel.ptr() );
    
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
    solveCriteria.extraParameters = Teuchos::parameterList("Nonlinear Solve");
    solveCriteria.extraParameters->set("Max Iters",int(maxIters));

    newtonSolver.setModel(thyraModel);
    Thyra::SolveStatus<double>
      solveStatus = Thyra::solve( newtonSolver, &*x, &solveCriteria );
    
    *out << "\nNonlinear solver return status:\n";
    {
      Teuchos::OSTab tab(out);
      *out << solveStatus;
    }
    *out << "\nFinal solution for x=\n" << *x;
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success)
    
  return  success ? 0 : 1;
}
