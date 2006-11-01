//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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

// Includes for Rythmos:
#include "Rythmos_InterpolationBufferBaseTester.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ForwardEulerStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_InterpolationBufferAsStepper.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_HermiteInterpolator.hpp"
#include "../../../example/basicExample/ExampleApplication.hpp"
#include "../../../example/epetra/1Dfem/ExampleApplication1Dfem.hpp"
// Includes for Teuchos:
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RefCountPtr.hpp"
// Includes for Thyra:
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_TimeStepNewtonNonlinearSolver.hpp"
// Includes for Stratimikos:
#ifdef HAVE_RYTHMOS_STRATIMIKOS
#  include "Thyra_DefaultRealLinearSolverBuilder.hpp"
#endif
// Includes for Epetra:
#include "Epetra_SerialComm.h"
// Includes for stl:
#include <iostream>

int main(int argc, char *argv[])
{
  bool status = true; 
  try { // catch exceptions

    // Command line Processor directives:
    Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.addOutputSetupOptions(true);
#ifdef HAVE_RYTHMOS_STRATIMIKOS
    Thyra::DefaultRealLinearSolverBuilder lowsfCreator;
    lowsfCreator.setupCLP(&clp);
#endif // HAVE_RYTHMOS_STRATIMIKOS
    int outputLevel = -1;
    clp.setOption( "outputlevel", &outputLevel, "Verbosity level [0-4]" );
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;
    outputLevel = min(max(outputLevel,-1),4);

    Teuchos::RefCountPtr<Teuchos::FancyOStream> 
      out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
    Teuchos::RefCountPtr<Epetra_Comm> epetra_comm_ptr_ = Teuchos::rcp( new Epetra_SerialComm  );
    Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> > W_factory;

    // Create a list of explicit models:
    std::vector<Teuchos::RefCountPtr<Thyra::ModelEvaluator<double> > > explicit_model_vec;
    Teuchos::ParameterList explicitParams;
    explicitParams.set<bool>( "implicit", false );
    explicitParams.set<double>( "Lambda_min", -0.9 );
    explicitParams.set<double>( "Lambda_max", -0.01 );
    std::string lambda_fit = "linear";
    explicitParams.set( "Lambda_fit", lambda_fit );
    explicitParams.set<int>( "NumElements", 100 );
    explicitParams.set<double>( "x0", 10.0 );
    explicitParams.set<double>( "Coeff_s", 0.0 );
    Teuchos::RefCountPtr<ExampleApplication>
      basicExample = Teuchos::rcp(new ExampleApplication(epetra_comm_ptr_, explicitParams));
    Teuchos::RefCountPtr<Thyra::ModelEvaluator<double> >
      basicExamplemodel = Teuchos::rcp(new Thyra::EpetraModelEvaluator(basicExample,W_factory));
    Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<double> > vs = basicExamplemodel->get_x_space();
    explicit_model_vec.push_back(basicExamplemodel);

    /*
    // Create a list of implicit models:
    std::vector<Teuchos::RefCountPtr<Thyra::ModelEvaluator<double> > > implicit_model_vec;
    Teuchos::ParameterList implicitParams;
    implicitParams.set<int>( "NumElements", 201 );
    Teuchos::RefCountPtr<ExampleApplication1Dfem>
      femTransient = Teuchos::rcp(new ExampleApplication1Dfem(epetra_comm_ptr_,implicitParams));
#ifdef HAVE_RYTHMOS_STRATIMIKOS
    lowsfCreator.readParameters(out.get());
    lowsfCreator.getParameterList()->print(*out,2,true,false);
    W_factory = lowsfCreator.createLinearSolveStrategy("");
#endif // HAVE_RYTHMOS_STRATIMIKOS
    Teuchos::RefCountPtr<Thyra::ModelEvaluator<double> >
      femTransientmodel = Teuchos::rcp(new Thyra::EpetraModelEvaluator(femTransient,W_factory));
    //Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<double> > vs = model->get_x_space();
    implicit_model_vec.push_back(femTransientmodel);
    // Create a nonlinear solver for the implicit methods:
    Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<double> > nonlinearSolver;
    Teuchos::RefCountPtr<Thyra::TimeStepNewtonNonlinearSolver<double> >
      _nonlinearSolver = Teuchos::rcp(new Thyra::TimeStepNewtonNonlinearSolver<double>());
    double maxError = 0.01;
    _nonlinearSolver->defaultTol(1e-3*maxError);
    nonlinearSolver = _nonlinearSolver;
    */


    Teuchos::RefCountPtr<Rythmos::StepperBase<double> > stepper_ptr;

    Teuchos::RefCountPtr<Teuchos::ParameterList> fixedStepIntegratorParams = Teuchos::rcp(new Teuchos::ParameterList);
    double dt = 0.1;
    fixedStepIntegratorParams->set( "fixed_dt", dt );
    fixedStepIntegratorParams->set( "outputLevel", outputLevel );

    Teuchos::RefCountPtr<Teuchos::ParameterList> variableStepIntegratorParams = Teuchos::rcp(new Teuchos::ParameterList);
    variableStepIntegratorParams->set( "outputLevel", outputLevel );

    Teuchos::RefCountPtr<Rythmos::InterpolatorBase<double> > interpolator;
    Teuchos::RefCountPtr<Rythmos::InterpolationBuffer<double> > IB;
    Teuchos::RefCountPtr<Rythmos::InterpolationBufferBase<double> > integrator; 
    int buffersize = 10;

    // Create a list of concrete objects derived from InterpolationBufferBase:
    std::vector<Teuchos::RefCountPtr<Rythmos::InterpolationBufferBase<double> > > IB_vec;
    // Explicit models:
    for (unsigned int i=0 ; i < explicit_model_vec.size() ; ++i) 
    {
      // Forward Euler Stepper
      stepper_ptr = Teuchos::rcp( new Rythmos::ForwardEulerStepper<double>(explicit_model_vec[i]) );
      IB_vec.push_back(stepper_ptr);
      // InterpolationBufferAsStepper with ForwardEulerStepper and Linear InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(fixedStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,fixedStepIntegratorParams));
      IB_vec.push_back(integrator);

      // Explicit RK Stepepr
      stepper_ptr = Teuchos::rcp( new Rythmos::ExplicitRKStepper<double>(explicit_model_vec[i]) );
      IB_vec.push_back(stepper_ptr);
      // InterpolationBufferAsStepper with ExplicitRKStepper and Linear InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(fixedStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,fixedStepIntegratorParams));
      IB_vec.push_back(integrator);

      /*
      // Backward Euler Stepper
      stepper_ptr = Teuchos::rcp( new Rythmos::BackwardEulerStepper<double>(explicit_model_vec[i],nonlinearSolver) );
      IB_vec.push_back(stepper_ptr);
      // InterpolationBufferAsStepper with BackwardEulerStepper and Linear InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(fixedStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,fixedStepIntegratorParams));
      IB_vec.push_back(integrator);

      // Implicit BDF Stepper
      stepper_ptr = Teuchos::rcp( new Rythmos::ImplicitBDFStepper<double>(explicit_model_vec[i],nonlinearSolver) );
      IB_vec.push_back(stepper_ptr);
      // InterpolationBufferAsStepper with fixed step ImplicitBDFStepper and Linear InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(fixedStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,fixedStepIntegratorParams));
      IB_vec.push_back(integrator);
      // InterpolationBufferAsStepper with variable step ImplicitBDFStepper and Linear InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(variableStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,variableStepIntegratorParams));
      IB_vec.push_back(integrator);
      // InterpolationBufferAsStepper with fixed step ImplicitBDFStepper and Hermite InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::HermiteInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(fixedStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,fixedStepIntegratorParams));
      IB_vec.push_back(integrator);
      // InterpolationBufferAsStepper with variable step ImplicitBDFStepper and Hermite InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::HermiteInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(variableStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,variableStepIntegratorParams));
      IB_vec.push_back(integrator);
      */
    }
    /*
    // Implicit models:
    for (unsigned int i=0 ; i < implicit_model_vec.size() ; ++i)
    {
      // Backward Euler Stepper
      stepper_ptr = Teuchos::rcp( new Rythmos::BackwardEulerStepper<double>(explicit_model_vec[i],nonlinearSolver) );
      IB_vec.push_back(stepper_ptr);
      // InterpolationBufferAsStepper with BackwardEulerStepper and Linear InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(fixedStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,fixedStepIntegratorParams));
      IB_vec.push_back(integrator);

      // Implicit BDF Stepper
      stepper_ptr = Teuchos::rcp( new Rythmos::ImplicitBDFStepper<double>(explicit_model_vec[i],nonlinearSolver) );
      IB_vec.push_back(stepper_ptr);
      // InterpolationBufferAsStepper with fixed step ImplicitBDFStepper and Linear InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(fixedStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,fixedStepIntegratorParams));
      IB_vec.push_back(integrator);
      // InterpolationBufferAsStepper with variable step ImplicitBDFStepper and Linear InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(variableStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,variableStepIntegratorParams));
      IB_vec.push_back(integrator);
      // InterpolationBufferAsStepper with fixed step ImplicitBDFStepper and Hermite InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::HermiteInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(fixedStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,fixedStepIntegratorParams));
      IB_vec.push_back(integrator);
      // InterpolationBufferAsStepper with variable step ImplicitBDFStepper and Hermite InterpolationBuffer
      interpolator = Teuchos::rcp(new Rythmos::HermiteInterpolator<double>());
      IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
      IB->setParameterList(variableStepIntegratorParams);
      integrator = Teuchos::rcp(new Rythmos::InterpolationBufferAsStepper<double>(stepper_ptr,IB,variableStepIntegratorParams));
      IB_vec.push_back(integrator);
    }
    */
    // Interpolation Buffers:
    // linear:
    interpolator = Teuchos::rcp(new Rythmos::LinearInterpolator<double>());
    IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
    IB_vec.push_back(IB);
    interpolator = Teuchos::rcp(new Rythmos::HermiteInterpolator<double>());
    IB = Teuchos::rcp(new Rythmos::InterpolationBuffer<double>(interpolator,buffersize));
    IB_vec.push_back(IB);
    
    // Actually test the InterpolationBufferBase member functions:
    bool localstatus = false;
    Rythmos::InterpolationBufferBaseTester<double> IBtester;
    for (unsigned int i=0 ; i < IB_vec.size() ; ++i)
    {
      localstatus = IBtester.check(*(IB_vec[i]),*vs);
      if (!localstatus)
        status = false;
    }

   } // end try
    catch( const std::exception &excpt ) {
    std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
    status = false;
  }
  catch( ... ) {
    std::cerr << "*** Caught an unknown exception!\n";
    status = false;
  }

  if (status)
    std::cout << "Test passed." << std::endl;
  else
    std::cout << "Test failed." << std::endl;

  return status ? 0 : 1;
} // end main() [Doxygen looks for this!]
