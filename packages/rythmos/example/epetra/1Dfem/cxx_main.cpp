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
#endif // HAVE_MPI

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"

#include "ExampleApplication1Dfem.hpp"

// Includes for Rythmos:
#include "Rythmos_ConfigDefs.h"
//#include "ExampleApplicationRythmosInterface.hpp"
#include "Rythmos_ForwardEulerStepper.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"

// Includes for Thyra:
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_TimeStepNewtonNonlinearSolver.hpp"
#include "Thyra_DiagonalEpetraLinearOpWithSolveFactory.hpp"
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
enum STEP_METHOD { STEP_METHOD_FIXED, STEP_METHOD_VARIABLE };

int main(int argc, char *argv[])
{
  bool verbose = true; // verbosity level.
  bool result, success = true; // determine if the run was successfull

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try { // catch exceptions

#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm mpiComm = MPI_COMM_WORLD;
    int procRank = 0;
    int numProc;
    MPI_Comm_size( mpiComm, &numProc );
    MPI_Comm_rank( mpiComm, &procRank );
#endif // HAVE_MPI

    STEP_METHOD step_method = STEP_METHOD_FIXED;
    int numElements = 201; // number of elements in vector
    double finalTime = 1.0; // ODE final time
    int N = 100;  // number of steps to take
    const int num_methods = 4;
    const EMethod method_values[] = { METHOD_FE, METHOD_BE, METHOD_ERK, METHOD_BDF };
    const char * method_names[] = { "FE", "BE", "ERK", "BDF" };
    EMethod method_val = METHOD_BE;
    double maxError = 0.01;
    bool version = false;  // display version information 
    double reltol = 1.0e-2;
    double abstol = 1.0e-4;
    int maxOrder = 5;
    int outputLevel = 2; // outputLevel is used to control Rythmos verbosity

    // Parse the command-line options:
    Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.addOutputSetupOptions(true);
#ifdef HAVE_RYTHMOS_STRATIMIKOS
    Thyra::DefaultRealLinearSolverBuilder lowsfCreator;
    lowsfCreator.setupCLP(&clp);
#endif

    clp.setOption( "T", &finalTime, "Final time for simulation." );
    clp.setOption( "numelements", &numElements, "Problem size");
    clp.setOption( "method", &method_val, num_methods, method_values, method_names, "Integration method" );
    clp.setOption( "numsteps", &N, "Number of integration steps to take" );
    clp.setOption( "maxerror", &maxError, "Maximum error" );
    clp.setOption( "reltol", &reltol, "Relative Error Tolerance" );
    clp.setOption( "abstol", &abstol, "Absolute Error Tolerance" );
    clp.setOption( "maxorder", &maxOrder, "Maximum Implicit BDF order" );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not" );
    clp.setOption( "version", "run", &version, "Version of this code" );
    clp.setOption( "outputLevel", &outputLevel, "Debug Level for Rythmos" );


    Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    // 10/23/06 tscoffe:  bounds on Teuchos::EVerbosityLevel:
    outputLevel = min(max(outputLevel,-1),4);

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

    if (finalTime <= 0.0)
    {
      std::cerr << "Final simulation time must be > 0.0." << std::endl;
      return(1);
    }

    
    // Set up the parameter list for the application:
    Teuchos::ParameterList params;
    params.set( "NumElements", numElements );
#ifdef HAVE_MPI
    Teuchos::RefCountPtr<Epetra_Comm> epetra_comm_ptr_ = Teuchos::rcp( new Epetra_MpiComm(mpiComm) );
#else
    Teuchos::RefCountPtr<Epetra_Comm> epetra_comm_ptr_ = Teuchos::rcp( new Epetra_SerialComm  );
#endif // HAVE_MPI

    // Create the factory for the LinearOpWithSolveBase object
    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >
      W_factory;
    if((method_val == METHOD_BE) or (method_val == METHOD_BDF))
    {
      //W_factory = Teuchos::rcp(new Thyra::DiagonalEpetraLinearOpWithSolveFactory());
      //W_factory = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory());
#ifdef HAVE_RYTHMOS_STRATIMIKOS
      W_factory = lowsfCreator.createLinearSolveStrategy("");
      *out
        << "\nCreated a LinearOpWithSolveFactory described as:\n"
        << Teuchos::describe(*W_factory,Teuchos::VERB_MEDIUM);
#endif

    }

    // create interface to problem
    Teuchos::RefCountPtr<ExampleApplication1Dfem>
      epetraModel = Teuchos::rcp(new ExampleApplication1Dfem(epetra_comm_ptr_,params));
    Teuchos::RefCountPtr<Thyra::ModelEvaluator<double> >
      model = Teuchos::rcp(new Thyra::EpetraModelEvaluator(epetraModel,W_factory));

    // Create Stepper object depending on command-line input
    std::string method;
    Teuchos::RefCountPtr<Rythmos::StepperBase<double> > stepper_ptr;
    if ( method_val == METHOD_ERK ) {
      stepper_ptr = Teuchos::rcp(new Rythmos::ExplicitRKStepper<double>(model));
      method = "Explicit Runge-Kutta of order 4";
    } else if (method_val == METHOD_FE) {
      stepper_ptr = Teuchos::rcp(new Rythmos::ForwardEulerStepper<double>(model));
      method = "Forward Euler";
    } else if (method_val == METHOD_BE) {
      Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<double> >
        nonlinearSolver;
      Teuchos::RefCountPtr<Thyra::TimeStepNewtonNonlinearSolver<double> >
        _nonlinearSolver = Teuchos::rcp(new Thyra::TimeStepNewtonNonlinearSolver<double>());
      _nonlinearSolver->defaultTol(1e-3*maxError);
      nonlinearSolver = _nonlinearSolver;
      stepper_ptr = Teuchos::rcp(new Rythmos::BackwardEulerStepper<double>(model,nonlinearSolver));
      method = "Backward Euler";
    } else if (method_val == METHOD_BDF) {
      Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<double> >
        nonlinearSolver;
      Teuchos::RefCountPtr<Thyra::TimeStepNewtonNonlinearSolver<double> >
        _nonlinearSolver = Teuchos::rcp(new Thyra::TimeStepNewtonNonlinearSolver<double>());
      _nonlinearSolver->defaultTol(1e-3*maxError);
      nonlinearSolver = _nonlinearSolver;
      Teuchos::RefCountPtr<Teuchos::ParameterList> BDFparams = Teuchos::rcp(new Teuchos::ParameterList);
      BDFparams->set( "stopTime", finalTime );
      BDFparams->set( "maxOrder", maxOrder );
      BDFparams->set( "relErrTol", reltol );
      BDFparams->set( "absErrTol", abstol );
      BDFparams->set( "outputLevel", outputLevel );
      stepper_ptr = Teuchos::rcp(new Rythmos::ImplicitBDFStepper<double>(model,nonlinearSolver,BDFparams));
      step_method = STEP_METHOD_VARIABLE;
      method = "Implicit BDF";
    } else {
      TEST_FOR_EXCEPT(true);
    }
    Rythmos::StepperBase<double> &stepper = *stepper_ptr;

    Teuchos::RefCountPtr<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    if (outputLevel >= 3)
    {
      stepper.describe(*out,static_cast<Teuchos::EVerbosityLevel>(outputLevel));
    }

    double t0 = 0.0;
    double t1 = finalTime;
    double dt = (t1-t0)/N;
    double time = t0;
    int numSteps = 0;

    if (step_method == STEP_METHOD_FIXED)
    {
      // Integrate forward with fixed step sizes:
      for (int i=1 ; i<=N ; ++i)
      {
        double dt_taken = stepper.TakeStep(dt,Rythmos::FIXED_STEP);
        numSteps++;
        if (dt_taken != dt)
        {
          cerr << "Error, stepper took step of dt = " << dt_taken << " when asked to take step of dt = " << dt << std::endl;
          break;
        }
        time += dt_taken;
      }
    }
    else // step_method == STEP_METHOD_VARIABLE
    {
      while (time < finalTime)
      {
        double dt_taken = stepper.TakeStep(0.0,Rythmos::VARIABLE_STEP);
        numSteps++;
        if (outputLevel >= 3)
        {
          stepper.describe(*out,static_cast<Teuchos::EVerbosityLevel>(outputLevel));
        }
        if (dt_taken < 0)
        {
          cerr << "Error, stepper failed for some reason with step taken = " << dt_taken << endl;
          break;
        }
        time += dt_taken;
        *out << "Took stepsize of: " << dt_taken << " time = " << time << endl;
      }
    }
    *out << "Integrated to time = " << time << endl;
    // Get solution out of stepper:
    Teuchos::RefCountPtr<const Thyra::VectorBase<double> > x_computed_thyra_ptr = stepper.get_solution();
    // Convert Thyra::VectorBase to Epetra_Vector
    Teuchos::RefCountPtr<const Epetra_Vector>
      x_computed_ptr = Thyra::get_Epetra_Vector(*(epetraModel->get_x_map()),x_computed_thyra_ptr);
    const Epetra_Vector &x_computed = *x_computed_ptr;

    // compute exact answer
    Teuchos::RefCountPtr<Epetra_Vector> x_star_ptr = epetraModel->get_exact_solution(t1);
    Epetra_Vector& x_star = *x_star_ptr;

    int MyPID = x_computed.Comm().MyPID();
    if (MyPID == 0)
    {
      *out << "Integrating 1DfemTransient t = " << t0 
                << " to t = " << t1 << std::endl;
      *out << "using " << method << "." << std::endl;
      *out << "Took " << numSteps << " steps." << std::endl;
    }
    int MyLength = x_computed.MyLength();
    double error = 0;
    double errorMag = 0;
    for (int i=0 ; i<MyLength ; ++i)
    {
      if(verbose)
      {
        *out << std::setprecision(15);
        *out << "Computed: x[" << MyPID*MyLength+i << "] = ";
        *out << std::setw(20); *out << x_computed[i] << "\t";
        *out << "Exact: x[" << MyPID*MyLength+i << "] = ";
        *out << std::setw(20); *out << x_star[i] << std::endl;
      }
      //const double thisError = Thyra::relErr(x_computed[i],x_star[i]);
      const double thisError = x_computed[i]-x_star[i];
      error += thisError*thisError;
      errorMag += x_star[i]*x_star[i];
    }
    error = sqrt(error)/sqrt(errorMag);
    result = Thyra::testMaxErr(
      "error",error
      ,"maxError",maxError
      ,"maxWarning",10.0*maxError
      ,&std::cerr,""
      );
    if(!result) success = false;

#ifdef HAVE_RYTHMOS_STRATIMIKOS
    // Write the final parameters to file
    if(W_factory.get())
      lowsfCreator.writeParamsFile(*W_factory);
#endif
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif // HAVE_MPI

   } // end try
   TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

  return success ? 0 : 1;
} // end main() [Doxygen looks for this!]

