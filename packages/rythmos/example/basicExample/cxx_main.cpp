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

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"

#include "ExampleApplication.hpp"

// Includes for Rythmos:
#include "Rythmos_ConfigDefs.h"
//#include "ExampleApplicationRythmosInterface.hpp"
#include "Rythmos_Stepper_ForwardEuler.hpp"
#include "Rythmos_Stepper_BackwardEuler.hpp"
#include "Rythmos_Stepper_ExplicitRK.hpp"
#include "Rythmos_Stepper_ImplicitBDF.hpp"

// Includes for Thyra:
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_LinearNonlinearSolver.hpp"
#include "Thyra_TimeStepNewtonNonlinearSolver.hpp"
#include "Thyra_DiagonalEpetraLinearOpWithSolveFactory.hpp"
#include "Thyra_TestingTools.hpp"

// Includes for Amesos:
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"

// Includes for AztecOO:
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"

#include <string>

// Includes for Teuchos:
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"

enum EMethod { METHOD_FE, METHOD_BE, METHOD_ERK, METHOD_BDF };
enum EStepMethod { FIXED_STEP, VARIABLE_STEP };

int main(int argc, char *argv[])
{
  bool verbose = false; // verbosity level.
  bool result, success = true; // determine if the run was successfull

  try { // catch exceptions

#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm mpiComm = MPI_COMM_WORLD;
    int procRank = 0;
    int numProc;
    MPI_Comm_size( mpiComm, &numProc );
    MPI_Comm_rank( mpiComm, &procRank );
#endif // HAVE_MPI

    double lambda_min = -0.9;   // min ODE coefficient
    double lambda_max = -0.01;  // max ODE coefficient
    std::string lambda_fit = "linear"; // Lambda model
    int numElements = 1; // number of elements in vector
    double x0 = 10.0; // ODE initial condition
    double finalTime = 1.0; // ODE final time
    int N = 10;  // number of steps to take
    const int num_methods = 4;
    const EMethod method_values[] = { METHOD_FE, METHOD_BE, METHOD_ERK, METHOD_BDF };
    const char * method_names[] = { "FE", "BE", "ERK", "BDF" };
    EMethod method_val = METHOD_ERK;
    const int num_step_methods = 2;
    const EStepMethod step_method_values[] = { FIXED_STEP, VARIABLE_STEP };
    const char * step_method_names[] = { "fixed", "variable" };
    EStepMethod step_method_val = FIXED_STEP;
    double maxError = 1e-6;
    bool version = false;  // display version information 
    double reltol = 1.0e-2;
    double abstol = 1.0e-4;
    int maxOrder = 5;
#ifdef Rythmos_DEBUG
    int debugLevel = 2; // debugLevel is used when Rythmos_DEBUG ifdef is set.
#endif // Rythmos_DEBUG

    // Parse the command-line options:
    Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "x0", &x0, "Constant ODE initial condition." );
    clp.setOption( "T", &finalTime, "Final time for simulation." );
    clp.setOption( "lambda_min", &lambda_min, "Lower bound for ODE coefficient");
    clp.setOption( "lambda_max", &lambda_max, "Upper bound for ODE coefficient");
    clp.setOption( "lambda_fit", &lambda_fit, "Lambda model:  random, linear");
    clp.setOption( "numelements", &numElements, "Problem size");
    clp.setOption( "method", &method_val, num_methods, method_values, method_names, "Integration method" );
    clp.setOption( "stepmethod", &step_method_val, num_step_methods, step_method_values, step_method_names, "Stepping method" );
    clp.setOption( "numsteps", &N, "Number of integration steps to take" );
    clp.setOption( "maxerror", &maxError, "Maximum error" );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not" );
    clp.setOption( "version", "run", &version, "Version of this code" );
    clp.setOption( "reltol", &reltol, "Relative Error Tolerance" );
    clp.setOption( "abstol", &abstol, "Absolute Error Tolerance" );
    clp.setOption( "maxorder", &maxOrder, "Maximum Implicit BDF order" );
#ifdef Rythmos_DEBUG
    clp.setOption( "debuglevel", &debugLevel, "Debug Level for Rythmos" );
#endif // Rythmos_DEBUG

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    if (version) // Display version information and exit.
    {
      std::cout << Rythmos::Rythmos_Version() << std::endl; 
      std::cout << "basicExample Version 0.1 - 06/23/05" << std::endl;
      return(0);
    }

    if (lambda_min > lambda_max)
    {
      std::cerr << "lamba_min must be less than lambda_max" << std::endl;
      return(1);
    }

    if (finalTime <= 0.0)
    {
      std::cerr << "Final simulation time must be > 0.0." << std::endl;
      return(1);
    }

#ifdef Rythmos_DEBUG
    std::cout.precision(15);
#endif // Rythmos_DEBUG
    
    // Set up the parameter list for the application:
    Teuchos::ParameterList params;
    bool implicitFlag = ((method_val==METHOD_BE) | (method_val==METHOD_BDF));
    //std::cout << "implicitFlag = " << implicitFlag << std::endl;
    params.set( "implicit", implicitFlag );
    params.set( "Lambda_min", lambda_min );
    params.set( "Lambda_max", lambda_max );
    params.set( "Lambda_fit", lambda_fit );
    params.set( "NumElements", numElements );
    params.set( "x0", x0 );
#ifdef HAVE_MPI
    params.set( "MPIComm", mpiComm );
#endif // HAVE_MPI

    // Create the factory for the LinearOpWithSolveBase object
    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >
      W_factory;
    if((method_val == METHOD_BE) | (method_val == METHOD_BDF))
      //W_factory = Teuchos::rcp(new Thyra::DiagonalEpetraLinearOpWithSolveFactory());
      W_factory = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory());
      //W_factory = Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory());

    // create interface to problem
    Teuchos::RefCountPtr<ExampleApplication>
      epetraModel = Teuchos::rcp(new ExampleApplication(params));
    Teuchos::RefCountPtr<Thyra::ModelEvaluator<double> >
      model = Teuchos::rcp(new Thyra::EpetraModelEvaluator(epetraModel,W_factory));

    // Create Stepper object depending on command-line input
    std::string method;
    Teuchos::RefCountPtr<Rythmos::Stepper<double> > stepper_ptr;
    if ( method_val == METHOD_ERK ) {
      stepper_ptr = Teuchos::rcp(new Rythmos::ExplicitRKStepper<double>(model));
      method = "Explicit Runge-Kutta of order 4";
      step_method_val = FIXED_STEP;
    } else if (method_val == METHOD_FE) {
      stepper_ptr = Teuchos::rcp(new Rythmos::ForwardEulerStepper<double>(model));
      method = "Forward Euler";
      step_method_val = FIXED_STEP;
    } else if ((method_val == METHOD_BE) | (method_val == METHOD_BDF)) {
      Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<double> >
        nonlinearSolver;
//    nonlinearSolver = Teuchos::rcp(new Thyra::LinearNonlinearSolver<double>());
      Teuchos::RefCountPtr<Thyra::TimeStepNewtonNonlinearSolver<double> >
        _nonlinearSolver = Teuchos::rcp(new Thyra::TimeStepNewtonNonlinearSolver<double>());
      _nonlinearSolver->defaultTol(1e-3*maxError);
      nonlinearSolver = _nonlinearSolver;
      if (method_val == METHOD_BE)
      {
        stepper_ptr = Teuchos::rcp(new Rythmos::BackwardEulerStepper<double>(model,nonlinearSolver));
        method = "Backward Euler";
        step_method_val = FIXED_STEP;
      } 
      else 
      {
        //TEST_FOR_EXCEPT(true); // RAB: Above is commented out due to lack of other commits
        Teuchos::ParameterList BDFparams;
        BDFparams.set( "stopTime", finalTime );
        BDFparams.set( "maxOrder", maxOrder );
        BDFparams.set( "relErrTol", reltol );
        BDFparams.set( "absErrTol", abstol );
#ifdef Rythmos_DEBUG
        BDFparams.set( "debugLevel", debugLevel );
#endif // Rythmos_DEBUG

        stepper_ptr = Teuchos::rcp(new Rythmos::ImplicitBDFStepper<double>(model,nonlinearSolver,BDFparams));
        method = "Implicit BDF";
        // step_method_val setting is left alone in this case
      }
    } else {
      TEST_FOR_EXCEPT(true);
    }
    Rythmos::Stepper<double> &stepper = *stepper_ptr;
#ifdef Rythmos_DEBUG
    Teuchos::RefCountPtr<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    if (debugLevel > 1)
    {
      stepper.describe(*out,Teuchos::VERB_EXTREME);
    }
#endif // Rythmos_DEBUG

    int numSteps = 0;
    double t0 = 0.0;
    double dt = (finalTime-t0)/N;
    double time = t0;

    if (step_method_val == FIXED_STEP)
    {
      // Integrate forward with fixed step sizes:
      for (int i=1 ; i<=N ; ++i)
      {
        double dt_taken = stepper.TakeStep(dt);
        numSteps++;
#ifdef Rythmos_DEBUG
        if (debugLevel > 1)
        {
          stepper.describe(*out,Teuchos::VERB_EXTREME);
        }
#endif // Rythmos_DEBUG
        if (dt_taken != dt)
        {
          cerr << "Error, stepper took step of dt = " << dt_taken << " when asked to take step of dt = " << dt << std::endl;
          break;
        }
      }
    }
    else // (step_method_val == VARIABLE_STEP)
    {
      while (time < finalTime)
      {
        double dt_taken = stepper.TakeStep();
        numSteps++;
#ifdef Rythmos_DEBUG
        if (debugLevel > 1)
        {
          stepper.describe(*out,Teuchos::VERB_EXTREME);
        }
#endif // Rythmos_DEBUG
        if (dt_taken < 0)
        {
          cerr << "Error, stepper failed for some reason with step taken = " << dt_taken << endl;
          break;
        }
        time += dt_taken;
        cout << "Took stepsize of: " << dt_taken << " time = " << time << endl;
      }
    }
    cout << "Integrated to time = " << time << endl;

    // Get solution out of stepper:
    Teuchos::RefCountPtr<const Thyra::VectorBase<double> > x_computed_thyra_ptr = stepper.get_solution();
    // Convert Thyra::VectorBase to Epetra_Vector
    Teuchos::RefCountPtr<const Epetra_Vector>
      x_computed_ptr = Thyra::get_Epetra_Vector(*(epetraModel->get_x_map()),x_computed_thyra_ptr);
    const Epetra_Vector &x_computed = *x_computed_ptr;

    // compute exact answer
    Teuchos::RefCountPtr<const Epetra_Vector> x_star_ptr = epetraModel->get_exact_solution(finalTime);
    const Epetra_Vector& x_star = *x_star_ptr;
    
    // get lambda from the problem:
    Teuchos::RefCountPtr<const Epetra_Vector> lambda_ptr = epetraModel->get_coeff();
    const Epetra_Vector &lambda = *lambda_ptr;

    // compute numerical exact answer (for FE and BE)
    Teuchos::RefCountPtr<const Epetra_Vector> x_numerical_exact_ptr; 
    if (method_val == METHOD_FE) 
    {
      Teuchos::RefCountPtr<Epetra_Vector> x_exact_ptr = Teuchos::rcp(new Epetra_Vector(x_star.Map()));
      Epetra_Vector& x_exact = *x_exact_ptr;
      int myN = x_exact.MyLength();
      for ( int i=0 ; i<myN ; ++i)
      {
        x_exact[i] = x0*pow(1+lambda[i]*dt,N);
      }
      x_numerical_exact_ptr = x_exact_ptr;
    } 
    else if (method_val == METHOD_BE) 
    {
      Teuchos::RefCountPtr<Epetra_Vector> x_exact_ptr = Teuchos::rcp(new Epetra_Vector(x_star.Map()));
      Epetra_Vector& x_exact = *x_exact_ptr;
      int myN = x_exact.MyLength();
      for ( int i=0 ; i<myN ; ++i)
      {
        x_exact[i] = x0*pow(1/(1-lambda[i]*dt),N);
      }
      x_numerical_exact_ptr = x_exact_ptr;
    }
    else if (method_val == METHOD_BDF)
    {
      // exact bdf solution here?
    }


    // 06/03/05 tscoffe to get an Epetra_Map associated with an Epetra_Vector:
    // x.Map()
    // to get an Epetra_Comm associated with an Epetra_Vector:
    // x.Comm()
    
  //  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm = (*epetraModel).get_epetra_comm();
    //int MyPID = epetraModel->get_Epetra_Map()->Comm()->MyPID();
    int MyPID = x_computed.Comm().MyPID();
  //  int MyPID = epetra_comm->MyPID();
    if (MyPID == 0)
    {
      std::cout << "Integrating \\dot{x}=\\lambda x from t = " << t0 
                << " to t = " << finalTime << std::endl;
      std::cout << "using " << method << std::endl;
      std::cout << "with initial x_0 = " << x0
                << ", \\Delta t = " << dt  << "." << std::endl;
      std::cout << "Took " << numSteps << " steps." << std::endl;
    }
    int MyLength = x_computed.MyLength();
    if (verbose)
    {
      for (int i=0 ; i<MyLength ; ++i)
      {
        std::cout.precision(15);
        std::cout << "lambda[" << MyPID*MyLength+i << "] = " << lambda[i] << std::endl;
      }
      // Print out computed and exact solutions:
      for (int i=0 ; i<MyLength ; ++i)
      {
        std::cout.precision(15);
        std::cout << "Computed: x[" << MyPID*MyLength+i << "] = ";
        std::cout.width(20); std::cout << x_computed[i] << "\t";
        std::cout << "Exact: x[" << MyPID*MyLength+i << "] = ";
        std::cout.width(20); std::cout << x_star[i] << std::endl;
      }
    }
    
    // Check numerics against exact numerical method for FE and BE case:
    double numerical_error = 0;
    double numerical_error_mag = 0;
    if (x_numerical_exact_ptr.get())
    {
      const Epetra_Vector& x_numerical_exact = *x_numerical_exact_ptr;
      for ( int i=0 ; i<MyLength ; ++i)
      {
        if (verbose) 
        {
          std::cout.precision(15);
          std::cout << "Computed: x[" << MyPID*MyLength+i << "] = ";
          std::cout.width(20); std::cout << x_computed[i] << "\t";
          std::cout << "Numerical Exact: x[" << MyPID*MyLength+i << "] = ";
          std::cout.width(20); std::cout << x_numerical_exact[i] << std::endl;
        }
        const double thisError = x_numerical_exact[i]-x_computed[i];
        numerical_error += thisError*thisError;
        numerical_error_mag += x_numerical_exact[i]*x_numerical_exact[i];
      }
      numerical_error = sqrt(numerical_error)/sqrt(numerical_error_mag);
      result = Thyra::testMaxErr(
        "Exact numerical error",numerical_error
        ,"maxError",1.0e-12
        ,"maxWarning",1.0e-11
        ,&std::cerr,""
        );
      if(!result) success = false;
    }

    // Check numerics against exact DE solution:
    double error = 0;
    double errorMag = 0;
    for (int i=0 ; i<MyLength ; ++i)
    {
      const double thisError = x_computed[i]-x_star[i];
      error += thisError*thisError;
      errorMag += x_star[i]*x_star[i];
    }
    error = sqrt(error)/sqrt(errorMag);
    result = Thyra::testMaxErr(
      "Exact DE solution error",error
      ,"maxError",maxError
      ,"maxWarning",10.0*maxError
      ,&std::cerr,""
      );
    if(!result) success = false;

    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif // HAVE_MPI

   } // end try
   catch( const std::exception &excpt ) {
    if(verbose)
      std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  catch( ... ) {
    if(verbose)
      std::cerr << "*** Caught an unknown exception!\n";
    success = false;
  }

  return success ? 0 : 1;
} // end main() [Doxygen looks for this!]

