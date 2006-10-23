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

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif // HAVE_MPI

#include "ExampleApplication.hpp"

// Includes for Rythmos:
#include "Rythmos_ConfigDefs.h"
//#include "ExampleApplicationRythmosInterface.hpp"
#include "Rythmos_ForwardEulerStepper.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
// 10/9/06 tscoffe:  InterpolationBufferAsStepper includes: 
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_InterpolationBufferAsStepper.hpp"

// Includes for Thyra:
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_LinearNonlinearSolver.hpp"
#include "Thyra_TimeStepNewtonNonlinearSolver.hpp"
#include "Thyra_TestingTools.hpp"

// Includes for Stratimikos:
#ifdef HAVE_RYTHMOS_STRATIMIKOS
#  include "Thyra_DefaultRealLinearSolverBuilder.hpp"
#endif

#include <string>

// Includes for Teuchos:
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

enum EMethod { METHOD_FE, METHOD_BE, METHOD_ERK, METHOD_BDF };
enum EStepMethod { FIXED_STEP, VARIABLE_STEP };

int main(int argc, char *argv[])
{
  bool verbose = false; // verbosity level.
  bool result, success = true; // determine if the run was successfull

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

#ifdef HAVE_MPI
  MPI_Comm mpiComm = MPI_COMM_WORLD;
#endif // HAVE_MPI

  try { // catch exceptions

    double lambda_min = -0.9;   // min ODE coefficient
    double lambda_max = -0.01;  // max ODE coefficient
    double coeff_s = 0.0;  // Default is no forcing term
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
    bool useIntegrator = false;
#ifdef Rythmos_DEBUG
    int debugLevel = 2; // debugLevel is used when Rythmos_DEBUG ifdef is set.
#endif // Rythmos_DEBUG
    int outputLevel = -1; // outputLevel determines the level of output / verbosity

    // Parse the command-line options:
    Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.addOutputSetupOptions(true);
#ifdef HAVE_RYTHMOS_STRATIMIKOS
    Thyra::DefaultRealLinearSolverBuilder lowsfCreator;
    lowsfCreator.setupCLP(&clp);
#endif
    clp.setOption( "x0", &x0, "Constant ODE initial condition." );
    clp.setOption( "T", &finalTime, "Final time for simulation." );
    clp.setOption( "lambda_min", &lambda_min, "Lower bound for ODE coefficient");
    clp.setOption( "lambda_max", &lambda_max, "Upper bound for ODE coefficient");
    clp.setOption( "lambda_fit", &lambda_fit, "Lambda model:  random, linear");
    clp.setOption( "force_coeff", &coeff_s, "Forcing term coefficient");
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
    clp.setOption( "useintegrator", "normal", &useIntegrator, "Use InterpolationBufferAsStepper as integrator" );
#ifdef Rythmos_DEBUG
    clp.setOption( "debuglevel", &debugLevel, "Debug Level for Rythmos" );
#endif // Rythmos_DEBUG
    clp.setOption( "outputLevel", &outputLevel, "Verbosity level for Rythmos" );

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

#ifdef HAVE_RYTHMOS_STRATIMIKOS
    lowsfCreator.readParameters(out.get());
    *out << "\nThe parameter list after being read in:\n";
    lowsfCreator.getParameterList()->print(*out,2,true,false);
#endif

    if (version) // Display version information and exit.
    {
      *out << Rythmos::Rythmos_Version() << std::endl; 
      *out << "basicExample Version 0.1 - 06/23/05" << std::endl;
      return(0);
    }

    if (lambda_min > lambda_max)
    {
      *out << "lamba_min must be less than lambda_max" << std::endl;
      return(1);
    }

    if (finalTime <= 0.0)
    {
      *out << "Final simulation time must be > 0.0." << std::endl;
      return(1);
    }

#ifdef Rythmos_DEBUG
    *out << std::setprecision(15);
#endif // Rythmos_DEBUG
    
    // Set up the parameter list for the application:
    Teuchos::ParameterList params;
    bool implicitFlag = ((method_val==METHOD_BE) | (method_val==METHOD_BDF));
    //*out << "implicitFlag = " << implicitFlag << std::endl;
    params.set( "implicit", implicitFlag );
    params.set( "Lambda_min", lambda_min );
    params.set( "Lambda_max", lambda_max );
    params.set( "Lambda_fit", lambda_fit );
    params.set( "NumElements", numElements );
    params.set( "x0", x0 );
    params.set( "Coeff_s", coeff_s );
#ifdef HAVE_MPI
    Teuchos::RefCountPtr<Epetra_Comm> epetra_comm_ptr_ = Teuchos::rcp( new Epetra_MpiComm(mpiComm) );
#else
    Teuchos::RefCountPtr<Epetra_Comm> epetra_comm_ptr_ = Teuchos::rcp( new Epetra_SerialComm  );
#endif // HAVE_MPI

    // Create the factory for the LinearOpWithSolveBase object
    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >
      W_factory;
    if((method_val == METHOD_BE) | (method_val == METHOD_BDF)) {
#ifdef HAVE_RYTHMOS_STRATIMIKOS
      W_factory = lowsfCreator.createLinearSolveStrategy("");
      *out
        << "\nCreated a LinearOpWithSolveFactory described as:\n"
        << Teuchos::describe(*W_factory,Teuchos::VERB_MEDIUM);
#endif
    }

    // create interface to problem
    Teuchos::RefCountPtr<ExampleApplication>
      epetraModel = Teuchos::rcp(new ExampleApplication(epetra_comm_ptr_, params));
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
        Teuchos::RefCountPtr<Teuchos::ParameterList> BDFparams = Teuchos::rcp(new Teuchos::ParameterList);
        BDFparams->set( "stopTime", finalTime );
        BDFparams->set( "maxOrder", maxOrder );
        BDFparams->set( "relErrTol", reltol );
        BDFparams->set( "absErrTol", abstol );
        BDFparams->set( "outputLevel", outputLevel );

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

    Teuchos::RefCountPtr<const Thyra::VectorBase<double> > x_computed_thyra_ptr;
    if (step_method_val == FIXED_STEP)
    {
      if (useIntegrator)
      {
        // Set up fixed-step-size integration:
        Teuchos::RefCountPtr<Teuchos::ParameterList> 
          integratorParams = Teuchos::rcp(new Teuchos::ParameterList);
        integratorParams->set( "fixed_dt", dt );
        // Create integrator using stepper and linear interpolation buffer:
        Teuchos::RefCountPtr<Rythmos::Interpolator<double> > 
          linearInterpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
        Teuchos::RefCountPtr<Rythmos::InterpolationBuffer<double> > 
          IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(linearInterpolator,1000));
        Rythmos::InterpolationBufferAsStepper<double> integrator(stepper_ptr,IB,integratorParams);
        // Ask for desired time value:
        std::vector<double> time_vals;
        time_vals.push_back(finalTime);
        std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<double> > > x_vec;
        std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<double> > > xdot_vec;
        std::vector<double> accuracy_vec;
        bool status = integrator.GetPoints(time_vals,&x_vec,&xdot_vec,&accuracy_vec);
        if (!status) 
        {
          std::cout << "ERROR:  Integrator.GetPoints returned failure" << std::endl;
          return(-1);
        }
        // Get solution out of stepper:
        x_computed_thyra_ptr = x_vec[0];
      }
      else
      {
        // Integrate forward with fixed step sizes:
        for (int i=1 ; i<=N ; ++i)
        {
          double dt_taken = stepper.TakeStep(dt);
          time += dt_taken;
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
        // Get solution out of stepper:
        x_computed_thyra_ptr = stepper.get_solution();
      }
    }
    else // (step_method_val == VARIABLE_STEP)
    {
      if (useIntegrator)
      {
        // Set up fixed-step-size integration:
        Teuchos::RefCountPtr<Teuchos::ParameterList> 
          integratorParams = Teuchos::rcp(new Teuchos::ParameterList);
        //integratorParams->set( "fixed_dt", dt );
        // Create integrator using stepper and linear interpolation buffer:
        Teuchos::RefCountPtr<Rythmos::Interpolator<double> > 
          linearInterpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
        Teuchos::RefCountPtr<Rythmos::InterpolationBuffer<double> > 
          IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(linearInterpolator,1000));
        Rythmos::InterpolationBufferAsStepper<double> integrator(stepper_ptr,IB,integratorParams);
        // Ask for desired time value:
        std::vector<double> time_vals;
        time_vals.push_back(finalTime);
        std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<double> > > x_vec;
        std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<double> > > xdot_vec;
        std::vector<double> accuracy_vec;
        bool status = integrator.GetPoints(time_vals,&x_vec,&xdot_vec,&accuracy_vec);
        if (!status) 
        {
          std::cout << "ERROR:  Integrator.GetPoints returned failure" << std::endl;
          return(-1);
        }
        // Get solution out of stepper:
        x_computed_thyra_ptr = x_vec[0];
      }
      else
      {
#ifdef Rythmos_DEBUG
        // Create a place to store the computed solutions
        x_computed_thyra_ptr = stepper.get_solution();
        // Convert Thyra::VectorBase to Epetra_Vector
        Teuchos::RefCountPtr<const Epetra_Vector> x_computed_ptr = Thyra::get_Epetra_Vector(*(epetraModel->get_x_map()),x_computed_thyra_ptr);
        // Create a place to store the exact numerical solution
        Teuchos::RefCountPtr<Epetra_Vector> x_numerical_exact_ptr = Teuchos::rcp(new Epetra_Vector(x_computed_ptr->Map()));
        Epetra_Vector& x_numerical_exact = *x_numerical_exact_ptr;
        // Create a place to store the relative difference:
        Teuchos::RefCountPtr<Epetra_Vector> x_rel_diff_ptr = Teuchos::rcp(new Epetra_Vector(x_computed_ptr->Map()));
        Epetra_Vector& x_rel_diff = *x_rel_diff_ptr;
        // get lambda from the problem:
        Teuchos::RefCountPtr<const Epetra_Vector> lambda_ptr = epetraModel->get_coeff();
        const Epetra_Vector &lambda = *lambda_ptr;
#endif // Rythmos_DEBUG
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
            *out << "Error, stepper failed for some reason with step taken = " << dt_taken << endl;
            break;
          }
#ifdef Rythmos_DEBUG
          // Get solution out of stepper:
          x_computed_thyra_ptr = stepper.get_solution();
          // Convert Thyra::VectorBase to Epetra_Vector
          x_computed_ptr = Thyra::get_Epetra_Vector(*(epetraModel->get_x_map()),x_computed_thyra_ptr);
          if ((method_val == METHOD_BDF) && (maxOrder == 1))
          {
            int myN = x_numerical_exact.MyLength();
            if (numSteps == 1) // First step
            {
              for (int i=0 ; i<myN ; ++i)
                x_numerical_exact[i] = x0;
            }
            for (int i=0 ; i<myN ; ++i)
              x_numerical_exact[i] = ( x_numerical_exact[i]
                  +dt_taken*epetraModel->evalR(time+dt_taken,lambda[i],coeff_s))
                  /(1-lambda[i]*dt_taken);
            for (int i=0 ; i<myN ; ++i)
              x_rel_diff[i] = (x_numerical_exact[i]-(*x_computed_ptr)[i])/x_numerical_exact[i];
            if (myN == 1)
              *out << "Computed x(" << time+dt_taken << ") = " << (*x_computed_ptr)[0] 
                  << "  Numerical Exact = " << x_numerical_exact[0] 
                  << "  Rel Diff = " << x_rel_diff[0] << std::endl;
            else
            {
              for (int i=0 ; i<myN ; ++i)
                *out << "Computed x_" << i << "(" << time+dt_taken << ") = " << (*x_computed_ptr)[i] 
                    << "  Numerical Exact = " << x_numerical_exact[i] 
                    << "  Rel Diff = " << x_rel_diff[i] <<  std::endl;
            }
          }
          else
          {
            // compute exact answer
            Teuchos::RefCountPtr<const Epetra_Vector> x_star_ptr = epetraModel->get_exact_solution(time);
            const Epetra_Vector& x_star = *x_star_ptr;
            int myN = x_computed_ptr->MyLength();
            for (int i=0 ; i<myN ; ++i)
              x_rel_diff[i] = (x_star[i]-(*x_computed_ptr)[i])/x_star[i];
            if (myN == 1)
              *out << "Computed x(" << time+dt_taken << ") = " << (*x_computed_ptr)[0] 
                  << "  Exact = " << x_star[0] 
                  << "  Rel Diff = " << x_rel_diff[0] << std::endl;
            else
            {
              for (int i=0 ; i<myN ; ++i)
                *out << "Computed x_" << i << "(" << time+dt_taken << ") = " << (*x_computed_ptr)[i] 
                    << "  Exact = " << x_star[i] 
                    << "  Rel Diff = " << x_rel_diff[i] << std::endl;
            }
          }

#endif // Rythmos_DEBUG
          time += dt_taken;
          *out << "Took stepsize of: " << dt_taken << " time = " << time << endl;
        }
        // Get solution out of stepper:
        x_computed_thyra_ptr = stepper.get_solution();
      }
    }
    *out << "Integrated to time = " << time << endl;

    // Convert solution from Thyra::VectorBase to Epetra_Vector
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
        x_exact[i] = x0;
        for (int j=1 ; j<=N ; ++j)
        {
          x_exact[i] = (1+lambda[i]*dt)*x_exact[i]+dt*epetraModel->evalR(0+j*dt,lambda[i],coeff_s);
        }
        //x_exact[i] = x0*pow(1+lambda[i]*dt,N);
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
        x_exact[i] = x0;
        for (int j=1 ; j<=N ; ++j)
        {
          x_exact[i] = (x_exact[i]+dt*epetraModel->evalR(0+j*dt,lambda[i],coeff_s))/(1-lambda[i]*dt);
        }
        //x_exact[i] = x0*pow(1/(1-lambda[i]*dt),N);
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
    
    int MyPID = x_computed.Comm().MyPID();
    if (MyPID == 0)
    {
      *out << "Integrating \\dot{x}=\\lambda x from t = " << t0 
                << " to t = " << finalTime << std::endl;
      *out << "using " << method << std::endl;
      *out << "with initial x_0 = " << x0
                << ", \\Delta t = " << dt  << "." << std::endl;
      *out << "Took " << numSteps << " steps." << std::endl;
    }
    int MyLength = x_computed.MyLength();
    if (verbose)
    {
      for (int i=0 ; i<MyLength ; ++i)
      {
        *out << std::setprecision(15);
        *out << "lambda[" << MyPID*MyLength+i << "] = " << lambda[i] << std::endl;
      }
      // Print out computed and exact solutions:
      for (int i=0 ; i<MyLength ; ++i)
      {
        *out << std::setprecision(15);
        *out << "Computed: x[" << MyPID*MyLength+i << "] = ";
        *out << std::setw(20) << x_computed[i] << "\t";
        *out << "Exact: x[" << MyPID*MyLength+i << "] = ";
        *out << std::setw(20) << x_star[i] << std::endl;
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
          *out << std::setprecision(15);
          *out << "Computed: x[" << MyPID*MyLength+i << "] = ";
          *out << std::setw(20) << x_computed[i] << "\t";
          *out << "Numerical Exact: x[" << MyPID*MyLength+i << "] = ";
          *out << std::setw(20) << x_numerical_exact[i] << std::endl;
        }
        const double thisError = x_numerical_exact[i]-x_computed[i];
        numerical_error += thisError*thisError;
        numerical_error_mag += x_numerical_exact[i]*x_numerical_exact[i];
      }
      numerical_error = sqrt(numerical_error)/sqrt(numerical_error_mag);
      result = Thyra::testMaxErr(
        "Exact numerical error",numerical_error
        ,"maxError",1.0e-10
        ,"maxWarning",1.0e-9
        ,&*out,""
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
      ,&*out,""
      );
    if(!result) success = false;

#ifdef HAVE_RYTHMOS_STRATIMIKOS
    // Write the final parameters to file
    if(W_factory.get())
      lowsfCreator.writeParamsFile(*W_factory);
#endif

   } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

  return success ? 0 : 1;

} // end main() [Doxygen looks for this!]
