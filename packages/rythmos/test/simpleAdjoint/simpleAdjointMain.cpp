//@HEADER

// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER 

#include "EpetraExt_DiagonalTransientModel.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_AdjointModelEvaluator.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearNonlinearSolver.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"
#include "Teuchos_as.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_MPI


/** \file simpleAdjointMain.cpp

\brief Simple example/test problem for basic adjoint solver
    
This simple test program provides the most basic test for the transient
adjoint solver capability in Rythmos.  This example uses the notation of
the SAND report "A derivation of forward and adjoint sensitivities for ODEs
and DAEs" by Roscoe A. Bartlett.

The linear ODE is:

\verbatim

  x_dot(t,i) - gamma(i) * x(t,i) = 0, for i=0...n-1, t in [0,t_final]
  x(0) = p

\verbatim

here <tt>gamma</tt> is a vector of scalar constants, and <tt>p</tt> are
the only parameters which are used to define the initial condition.

The global response function is of the terminal type and takes the form:

\verbatim

  d_hat(p) = h(x(p,t_final)) = 1/2 * x(p,t_final)^T * x(p,t_final)

\endverbatim

where <tt>x(p,t_final)</tt> is the solved-for state at time <tt>t_final</tt>
which is an implicit function of the parameters <tt>p</tt> through the initial
condition to the forward ODE.
  
the reverse-time adjoint ODE for this problem is:

\verbatim
    
  lambda_rev_dot(t_bar,i) - gamma(i) * lambda(t_bar,i) = 0,
     for i=0...n-1, t_bar in [0,t_final-t_0]

  lambda(t_final,i) = x(t_final), for i=0...n-1

\endverbatim

The reduced response function derivative in terms of the adjoint is:
 
\verbatim

  d(d_hat)/d(p) = lambda(t=t_0)

\endverbatim

Using a fixed time-step backward Euler method, the state solution is:

\verbatim

  x(i,k) = beta(i)^k * p(i), for i=0...n-1, k=0...N

\endverbatim

where <tt>k=0...N</tt> is the timestep index (with <tt>k=0</tt> being at
<tt>t_0</tt> and <tt>k=N</tt> being at <tt>t_final</tt>) and with

\vetbatim

  beta(i) = 1 / ( 1 - delta_t * gamma(i) )

\endverbatim

The exact solution for the reduced response function is:

\verbatim

  d_hat(p) = 0.5 * sum( x(i,N)^2, i=0...n-1 )
           = 0.5 * sum( beta(i)^(2*N) * p(i)^2, i=0...n-1 )

\endverbatim

which gives the exact reduced gradient:

\verbatim

  d(d_hat)/d(p(j)) = beta(j)^(2*N) * p(j), for j=0...n-1

\endverbatim

The exact discrete solution to the reverse adjoint equation for backward Euler
is:

\verbatim

  lambda(i,r) = beta(i)^r * x(i,N), for i=0...n-1, k=0...N

\endverbatim

The reduced gradient in terms of the adjoint is:

\verbatim

  d(d_hat)/d(p(i)) = lambda(i,N), for i=0...n-1

\endverbatim

This gives me a simple way to test the forward and adjoint solutions as well
as the redueced gradient computed from the adjoint.

*/

namespace {


const std::string TimeStepNonlinearSolver_name = "TimeStepNonlinearSolver";

const std::string Stratimikos_name = "Stratimikos";

const std::string DiagonalTransientModel_name = "DiagonalTransientModel";

const std::string RythmosStepper_name = "Rythmos Stepper";

const std::string RythmosIntegrationControl_name = "Rythmos Integration Control";

const std::string RythmosIntegrator_name = "Rythmos Integrator";


Teuchos::RCP<const Teuchos::ParameterList>
getValidParameters()
{
  using Teuchos::RCP; using Teuchos::ParameterList;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist(TimeStepNonlinearSolver_name).disableRecursiveValidation();
    pl->sublist(Stratimikos_name).disableRecursiveValidation();
    pl->sublist(DiagonalTransientModel_name).disableRecursiveValidation();
    pl->sublist(RythmosStepper_name).disableRecursiveValidation();
    pl->sublist(RythmosIntegrationControl_name).disableRecursiveValidation();
    pl->sublist(RythmosIntegrator_name).disableRecursiveValidation();
    validPL = pl;
  }
  return validPL;
}


template<class Scalar>
Scalar integralPow( const Scalar &x, const int n )
{
  Scalar rtn = x;
  for ( int i = 0; i < n-1; ++i )
    rtn *= x;
  return rtn;
}



} // namespace


int main(int argc, char *argv[])
{

  using std::endl;
  typedef double Scalar;
  typedef double ScalarMag;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::describe;
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::as;
  using Teuchos::ParameterList;
  using Teuchos::CommandLineProcessor;
  typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;
  typedef Thyra::ModelEvaluatorBase MEB;
  
  bool result, success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  RCP<Epetra_Comm> epetra_comm;
#ifdef HAVE_MPI
  epetra_comm = rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  epetra_comm = rcp( new Epetra_SerialComm );
#endif // HAVE_MPI

  RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read commandline options
    //

    CommandLineProcessor clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    std::string paramsFileName = "";
    clp.setOption( "params-file", &paramsFileName,
      "File name for XML parameters" );

    double t_final = 1e-3;
    clp.setOption( "final-time", &t_final,
      "Final integration time (initial time is 0.0)" );

    int numTimeSteps = 10;
    clp.setOption( "num-time-steps", &numTimeSteps,
      "Number of (fixed) time steps.  If <= 0.0, then variable time steps are taken" );

    double maxStateError = 1e-14;
    clp.setOption( "max-state-error", &maxStateError,
      "Maximum relative error in the integrated state allowed" );

    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
    setVerbosityLevelOption( "verb-level", &verbLevel,
      "Top-level verbosity level.  By default, this gets deincremented as you go deeper into numerical objects.",
      &clp );

    Teuchos::EVerbosityLevel solnVerbLevel = Teuchos::VERB_DEFAULT;
    setVerbosityLevelOption( "soln-verb-level", &solnVerbLevel,
      "Solution verbosity level",
      &clp );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    //
    *out << "\nA) Get the base parameter list ...\n";
    //
    
    RCP<ParameterList>
      paramList = Teuchos::parameterList();
    if (paramsFileName.length())
      updateParametersFromXmlFile( paramsFileName, paramList.ptr() );

    paramList->validateParameters(*getValidParameters());

    const Scalar t_init = 0.0;

    const Rythmos::TimeRange<Scalar> fwdTimeRange(t_init, t_final);
    const Scalar delta_t = t_final / numTimeSteps;
    *out << "\ndelta_t = " << delta_t;

    //
    *out << "\nB) Create the Stratimikos linear solver factory ...\n";
    //
    // This is the linear solve strategy that will be used to solve for the
    // linear system with the W.
    //

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(sublist(paramList,Stratimikos_name));
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      W_factory = createLinearSolveStrategy(linearSolverBuilder);

    //
    *out << "\nC) Create and initalize the forward model ...\n";
    //
    
    // C.1) Create the underlying EpetraExt::ModelEvaluator

    RCP<EpetraExt::DiagonalTransientModel> epetraStateModel =
      EpetraExt::diagonalTransientModel(
        epetra_comm,
        sublist(paramList,DiagonalTransientModel_name)
        );

    *out <<"\nepetraStateModel valid options:\n";
    epetraStateModel->getValidParameters()->print(
      *out, PLPrintOptions().indent(2).showTypes(true).showDoc(true)
      );
    
    // C.2) Create the Thyra-wrapped ModelEvaluator
    
    RCP<Thyra::ModelEvaluator<double> > fwdStateModel =
      epetraModelEvaluator(epetraStateModel, W_factory);

    const RCP<const Thyra::VectorSpaceBase<Scalar> >
      x_space = fwdStateModel->get_x_space();

    const RCP<const Thyra::VectorBase<Scalar> >
      gamma = Thyra::create_Vector(epetraStateModel->get_gamma(), x_space);
    *out << "\ngamma = " << describe(*gamma, solnVerbLevel);

    //
    *out << "\nD) Create the stepper and integrator for the forward problem ...\n";
    //

    RCP<Rythmos::TimeStepNonlinearSolver<double> > fwdTimeStepSolver =
      Rythmos::timeStepNonlinearSolver<double>();
    RCP<Rythmos::StepperBase<Scalar> > fwdStateStepper =
      Rythmos::backwardEulerStepper<double>(fwdStateModel, fwdTimeStepSolver);
    fwdStateStepper->setParameterList(sublist(paramList, RythmosStepper_name));
    RCP<Rythmos::IntegratorBase<Scalar> > fwdStateIntegrator;
    {
      RCP<ParameterList>
        integrationControlPL = sublist(paramList, RythmosIntegrationControl_name);
      integrationControlPL->set( "Take Variable Steps", false );
      integrationControlPL->set( "Fixed dt", as<double>(delta_t) );
      RCP<Rythmos::IntegratorBase<Scalar> >
        defaultIntegrator = Rythmos::controlledDefaultIntegrator<Scalar>(
          Rythmos::simpleIntegrationControlStrategy<Scalar>(integrationControlPL)
          );
      fwdStateIntegrator = defaultIntegrator;
    }
    fwdStateIntegrator->setParameterList(sublist(paramList, RythmosIntegrator_name));

    //
    *out << "\nE) Solve the forward problem ...\n";
    //

    const MEB::InArgs<Scalar>
      state_ic = fwdStateModel->getNominalValues();
    *out << "\nstate_ic:\n" << describe(state_ic,solnVerbLevel);

    fwdStateStepper->setInitialCondition(state_ic);
    fwdStateIntegrator->setStepper(fwdStateStepper, t_final);

    Array<RCP<const Thyra::VectorBase<Scalar> > > x_final_array;
    fwdStateIntegrator->getFwdPoints(
      Teuchos::tuple<Scalar>(t_final), &x_final_array, NULL, NULL
      );
    const RCP<const Thyra::VectorBase<Scalar> > x_final = x_final_array[0];

    *out << "\nx_final:\n" << describe(*x_final, solnVerbLevel);

    //
    *out << "\nF) Check the solution to the forward problem ...\n";
    //

    const RCP<Thyra::VectorBase<Scalar> >
      x_beta = createMember(x_space),
      x_final_be_exact = createMember(x_space);

    {
      Thyra::ConstDetachedVectorView<Scalar> d_gamma(*gamma);
      Thyra::ConstDetachedVectorView<Scalar> d_x_ic(*state_ic.get_x());
      Thyra::DetachedVectorView<Scalar> d_x_beta(*x_beta);
      Thyra::DetachedVectorView<Scalar> d_x_final_be_exact(*x_final_be_exact);
      const int n = d_gamma.subDim();
      for ( int i = 0; i < n; ++i ) {
        d_x_beta(i) = 1.0 / ( 1.0 - delta_t * d_gamma(i) );
        d_x_final_be_exact(i) = integralPow(d_x_beta(i), numTimeSteps) * d_x_ic(i);
      }
    }

    *out << "\nx_final_be_exact:\n" << describe(*x_final_be_exact, solnVerbLevel);

    result = Thyra::testRelNormDiffErr<Scalar>(
      "x_final", *x_final,
      "x_final_be_exact", *x_final_be_exact,
      "maxStateError", maxStateError,
      "warningTol", 1.0, // Don't warn
      &*out, solnVerbLevel
      );
    if (!result) success = false;

    //
    *out << "\nG) Create the Adjoint ME wrapper object ...\n";
    //

    RCP<Thyra::ModelEvaluator<double> > adjModel =
      Rythmos::adjointModelEvaluator<double>(
        fwdStateModel, fwdTimeRange
        );

    //
    *out << "\nH) Create a stepper and integrator for the adjoint ...\n";
    //

    RCP<Thyra::LinearNonlinearSolver<double> > adjTimeStepSolver =
      Thyra::linearNonlinearSolver<double>();
    RCP<Rythmos::StepperBase<Scalar> > adjStepper =
      Rythmos::backwardEulerStepper<double>(adjModel, adjTimeStepSolver);
    adjStepper->setParameterList(sublist(paramList, RythmosStepper_name));
    RCP<Rythmos::IntegratorBase<Scalar> > adjIntegrator =
      fwdStateIntegrator->cloneIntegrator();

    //
    *out << "\nI) Set up the initial condition for the adjoint at the final time ...\n";
    //

    const RCP<const Thyra::VectorSpaceBase<Scalar> >
      f_space = fwdStateModel->get_f_space();

    // lambda(t_final) = x_final
    const RCP<Thyra::VectorBase<Scalar> > lambda_ic = createMember(f_space);
    V_V( outArg(*lambda_ic), *x_final_be_exact );
    
    // lambda_dot(t_final,i) = - gamma(i) * lambda(t_final,i)
    const RCP<Thyra::VectorBase<Scalar> > lambda_dot_ic = createMember(f_space);
    Thyra::V_S<Scalar>( outArg(*lambda_dot_ic), ST::zero() );
    Thyra::ele_wise_prod<Scalar>( -ST::one(), *gamma, *lambda_ic,
      outArg(*lambda_dot_ic) );
    
    MEB::InArgs<Scalar> adj_ic = adjModel->getNominalValues();
    adj_ic.set_x(lambda_ic);
    adj_ic.set_x_dot(lambda_dot_ic);
    *out << "\nadj_ic:\n" << describe(adj_ic,solnVerbLevel);

    //
    *out << "\nJ) Integrate the adjoint backwards in time (using backward time) ...\n";
    //

    adjStepper->setInitialCondition(adj_ic);
    adjIntegrator->setStepper(adjStepper, fwdTimeRange.length());

    Array<RCP<const Thyra::VectorBase<Scalar> > > lambda_final_array;
    adjIntegrator->getFwdPoints(
      Teuchos::tuple<Scalar>(fwdTimeRange.length()), &lambda_final_array, NULL, NULL
      );
    const RCP<const Thyra::VectorBase<Scalar> > lambda_final = lambda_final_array[0];

    *out << "\nlambda_final:\n" << describe(*lambda_final, solnVerbLevel);

    //
    *out << "\nK) Test the final adjoint againt exact discrete solution ...\n";
    //

    {

      const RCP<Thyra::VectorBase<Scalar> >
        lambda_final_be_exact = createMember(x_space);
      
      {
        Thyra::ConstDetachedVectorView<Scalar> d_gamma(*gamma);
        Thyra::ConstDetachedVectorView<Scalar> d_x_final(*x_final);
        Thyra::DetachedVectorView<Scalar> d_x_beta(*x_beta);
        Thyra::DetachedVectorView<Scalar> d_lambda_final_be_exact(*lambda_final_be_exact);
        const int n = d_gamma.subDim();
        for ( int i = 0; i < n; ++i ) {
          d_lambda_final_be_exact(i) = integralPow(d_x_beta(i), numTimeSteps) * d_x_final(i);
        }
      }
      
      *out << "\nlambda_final_be_exact:\n" << describe(*lambda_final_be_exact, solnVerbLevel);
      
      result = Thyra::testRelNormDiffErr<Scalar>(
        "lambda_final", *lambda_final,
        "lambda_final_be_exact", *lambda_final_be_exact,
        "maxStateError", maxStateError,
        "warningTol", 1.0, // Don't warn
        &*out, solnVerbLevel
        );
      if (!result) success = false;

    }

    //
    *out << "\nL) Test the reduced gradient from the adjoint against the discrete forward reduced gradient ...\n";
    //

    {

      const RCP<const Thyra::VectorBase<Scalar> >
        d_d_hat_d_p_from_lambda = lambda_final; // See above

      const RCP<Thyra::VectorBase<Scalar> >
        d_d_hat_d_p_be_exact = createMember(x_space);
      
      {
        Thyra::ConstDetachedVectorView<Scalar> d_x_ic(*state_ic.get_x());
        Thyra::DetachedVectorView<Scalar> d_x_beta(*x_beta);
        Thyra::DetachedVectorView<Scalar> d_d_d_hat_d_p_be_exact(*d_d_hat_d_p_be_exact);
        const int n = d_x_ic.subDim();
        for ( int i = 0; i < n; ++i ) {
          d_d_d_hat_d_p_be_exact(i) = integralPow(d_x_beta(i), 2*numTimeSteps) * d_x_ic(i);
        }
      }
      
      *out << "\nd_d_hat_d_p_be_exact:\n" << describe(*d_d_hat_d_p_be_exact, solnVerbLevel);
      
      result = Thyra::testRelNormDiffErr<Scalar>(
        "d_d_hat_d_p_from_lambda", *d_d_hat_d_p_from_lambda,
        "d_d_hat_d_p_be_exact", *d_d_hat_d_p_be_exact,
        "maxStateError", maxStateError,
        "warningTol", 1.0, // Don't warn
        &*out, solnVerbLevel
        );
      if (!result) success = false;

    }
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success);

  if(success)
    *out << "\nEnd Result: TEST PASSED" << endl;
  else
    *out << "\nEnd Result: TEST FAILED" << endl;
  
  return ( success ? 0 : 1 );

} // end main() [Doxygen looks for this!]

