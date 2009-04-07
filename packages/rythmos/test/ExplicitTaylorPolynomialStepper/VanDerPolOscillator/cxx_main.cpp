//@HEADER
// ************************************************************************
// 
//                          Rythmos Package 
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "VanDerPolOscillator.hpp"

// Includes for Rythmos:
#include "Rythmos_ConfigDefs.h"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ExplicitTaylorPolynomialStepper.hpp"

// Includes for Thyra:
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_TestingTools.hpp"

#include <string>

// Includes for Teuchos:
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"

enum EMethod { METHOD_FE, METHOD_BE, METHOD_ERK, METHOD_BDF, METHOD_ETI };
enum EStepMethod { STEP_TYPE_FIXED, STEP_TYPE_VARIABLE };

int main(int argc, char *argv[])
{
  bool verbose = false; // verbosity level.
  bool success = true; // determine if the run was successfull

   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

#ifdef HAVE_MPI
  MPI_Comm mpiComm = MPI_COMM_WORLD;
#endif // HAVE_MPI

  try { // catch exceptions

    double omega = 1.0e-2;   // Van Der Pol coefficient
    double x0_1 = 2.0; // ODE initial condition
    double x0_2 = -0.66; // ODE initial condition
    double initialTime = 0.0; // ODE initial time
    double finalTime = 2.0; // ODE final time
    double minStep = 1.0e-10; // minimum step size
    double maxStep = 1.0; // maximum step size
    int degree = 40; // degree of taylor polynomial expansion
    double tol = 1.0e-10; // local error tolerance
    int N = 1000;  // number of steps to take
    string outfile_name = "vdp.out";
    const int num_methods = 2;
    const EMethod method_values[] = { METHOD_ERK, METHOD_ETI };
    const char * method_names[] = { "ERK", "ETI" };
    EMethod method_val = METHOD_ETI;
    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_NONE;

    // Parse the command-line options:
    Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.addOutputSetupOptions(true);
    clp.setOption( "x0_1", &x0_1, "First coordinate of initial condition." );
    clp.setOption( "x0_2", &x0_2, "Second coordinate of initial condition." );
    clp.setOption( "omega", &omega, "VDP Coefficient" );
    clp.setOption( "T_initial", &initialTime, "Initial time for simulation." );
    clp.setOption( "T_final", &finalTime, "Final time for simulation." );
    clp.setOption( "minStep", &minStep, "Minimum step size." );
    clp.setOption( "maxStep", &maxStep, "Maximum step size." );
    clp.setOption( "degree", &degree, 
		   "Degree of taylor polynomial expansion." );
    clp.setOption( "tol", &tol, "Local error tolerance." );
    clp.setOption( "numsteps", &N, "Number of integration steps to take" );
    clp.setOption( "method", &method_val, num_methods, method_values, 
		   method_names, "Integration method" );
    clp.setOption( "verbose", "quiet", &verbose, 
		   "Set if output is printed or not" );
    clp.setOption( "output", &outfile_name, "Output file name." );
    setVerbosityLevelOption( "verb-level", &verbLevel, "Overall verbosity level.", &clp );

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = 
      clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) 
      return parse_return;

    if (minStep > maxStep)
    {
      std::cerr << "min step must be smaller than max step" << std::endl;
      return(1);
    }

    if (finalTime <= initialTime)
    {
      std::cerr << "Final simulation time must be > initial time." 
		<< std::endl;
      return(1);
    }
    
    // Set up the parameter list for the application:
    Teuchos::ParameterList params;
    Teuchos::ParameterList& vdpoParams = params.sublist("Van Der Pol Oscillator Settings");
      vdpoParams.set( "implicit", false );
      vdpoParams.set( "x0_1", x0_1 );
      vdpoParams.set( "x0_2", x0_2 );
      vdpoParams.set( "omega", omega );
      vdpoParams.set( "Output File Name", outfile_name );

    Teuchos::ParameterList& etpParams = params.sublist("Explicit Taylor Polynomial Settings");
      etpParams.set( "Initial Time", initialTime );
      etpParams.set( "Final Time", finalTime );
      etpParams.set( "Local Error Tolerance", tol );
      etpParams.set( "Minimum Step Size", minStep );
      etpParams.set( "Maximum Step Size", maxStep );
      etpParams.set( "Taylor Polynomial Degree", Teuchos::as<unsigned int>(degree) );
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_Comm> epetra_comm_ptr_ = Teuchos::rcp( new Epetra_MpiComm(mpiComm) );
#else
    Teuchos::RCP<Epetra_Comm> epetra_comm_ptr_ = Teuchos::rcp( new Epetra_SerialComm  );
#endif // HAVE_MPI

    // Create the factory for the LinearOpWithSolveBase object
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
      W_factory;

    // create interface to problem
    Teuchos::RCP<VanDerPolOscillator>
      epetraModel = Teuchos::rcp(new VanDerPolOscillator(epetra_comm_ptr_,
							 vdpoParams));
    Teuchos::RCP<Thyra::ModelEvaluator<double> >
      model = Teuchos::rcp(new Thyra::EpetraModelEvaluator(epetraModel,
							   W_factory));

    // Create Stepper object depending on command-line input
    std::string method;
    Teuchos::RCP<Rythmos::StepperBase<double> > stepper_ptr;
    if ( method_val == METHOD_ERK ) {
      stepper_ptr = Rythmos::explicitRKStepper<double>(model);
      Teuchos::RCP<Teuchos::ParameterList> ERKparams = Teuchos::rcp(new Teuchos::ParameterList);
      ERKparams->sublist("VerboseObject").set(
        "Verbosity Level",
        Teuchos::getVerbosityLevelParameterValueName(verbLevel)
        );
      stepper_ptr->setParameterList(ERKparams);
      method = "Explicit Runge-Kutta of order 4";
    } 
    else if (method_val == METHOD_ETI) {
      method = "Explicit Taylor";
      stepper_ptr = Teuchos::rcp(new Rythmos::ExplicitTaylorPolynomialStepper<double>);
      stepper_ptr->setModel(model);
      stepper_ptr->setParameterList(Teuchos::rcp(&etpParams,false));
    }
    else {
      TEST_FOR_EXCEPT(true);
    }
    Rythmos::StepperBase<double> &stepper = *stepper_ptr;

    double t = initialTime;
    double dt;
    int step = 0;
    if (method_val == METHOD_ETI) {
      while ( fabs(finalTime-t) > 1.0e-14 ) {
	if (verbose)
	  cout << "t = " << t << endl;

	dt = stepper.takeStep(0.0, Rythmos::STEP_TYPE_VARIABLE);
	t += dt;
	step++;

	epetraModel->saveSolution(*get_Epetra_Vector(*(epetraModel->get_x_map()), stepper.getStepStatus().solution),t);
      }
      
    }
    else {
      dt = (finalTime-initialTime)/N;

      // Integrate forward with fixed step sizes:
      for (int i=1 ; i<=N ; ++i, step++) {
	if (verbose)
	  cout << "t = " << t << endl;
	double dt_taken = stepper.takeStep(dt, Rythmos::STEP_TYPE_FIXED);
	if (dt_taken != dt) {
	  cerr << "Error, stepper took step of dt = " << dt_taken 
	       << " when asked to take step of dt = " << dt << std::endl;
	  break;
	}
	t += dt_taken;
	epetraModel->saveSolution(*get_Epetra_Vector(*(epetraModel->get_x_map()), stepper.getStepStatus().solution),t);
      }
    }
    if (verbose)
      cout << "num steps = " << step << endl;

    // Get final solution
    Teuchos::RCP<const Epetra_Vector> final_x = 
      get_Epetra_Vector(*(epetraModel->get_x_map()), stepper.getStepStatus().solution);
    double final_tol = 1.0e-2;
    if (std::abs((*final_x)[0]-1.93704) < final_tol &&
      std::abs((*final_x)[1]+0.70225) < final_tol)
      success = true;
    else 
      success = false;

   } // end try
   catch( const std::exception &excpt ) {
     std::cerr << "*** Caught a standard exception : " 
	       << excpt.what() << std::endl;
     success = false;
   }
   catch( ... ) {
     std::cerr << "*** Caught an unknown exception!\n";
     success = false;
   }

  if(success)
    *out << "\nEnd Result: TEST PASSED" << endl;
  else
    *out << "\nEnd Result: TEST FAILED" << endl;

  return success ? 0 : 1;
} // end main() [Doxygen looks for this!]
