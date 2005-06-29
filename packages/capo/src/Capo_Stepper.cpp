//
// @HEADER
// ***********************************************************************
// 
//                           Capo Package
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
// ***********************************************************************
// @HEADER

/************************************************************ 
File:      Capo_Stepper.cpp
Purpose:   The main capo driver program.
Date:      6-10-05
Author:    Joseph Simonis
**************************************************************/

/**** Includes ****/
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Capo_Integrator.hpp"
#include "Capo_Parameter_List.hpp"
#include "Capo_Solver.hpp"
#include "Capo_Stepper.hpp"

using namespace CAPO;

//-----------------------------------------------------------------
// Function      : Stepper::Stepper
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
Stepper::Stepper(Teuchos::RefCountPtr<Parameter_List> PL, \
		 Teuchos::RefCountPtr<Solver> App_Solver)
{
  StepSize = PL->get_lambda_stepsize();
  PrevStepSize = StepSize;
  StepNumber = 0;
  PrintProc = PL->get_printproc();
  MaxSteps = PL->get_MaxOuterIts();

  //Problem_Integrator = App_Int; included in solver...
  Problem_Parameters = PL;
  iteration_method = App_Solver;
}

//-----------------------------------------------------------------
// Function      : Stepper::Done
// Purpose       : Check to see if continuation problem is finished
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
bool Stepper::Done() const
{
  if (StepNumber>MaxSteps)
    return true;
  else 
    return false;
}

//-----------------------------------------------------------------
// Function      : Stepper::PrintStart
// Purpose       : Indicate a new continuation step is beginning
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
void Stepper::PrintStart() const
{
  if (PrintProc>0)
    {
      cout << endl <<"---------- Start of Continuation step " << StepNumber << "----------" << endl;
      //cout << "Param = " << iteration_method->Get_lambdafinal() << ", StepSize = " << StepSize << " ~~~~~~" << endl;
    }
}
//-----------------------------------------------------------------
// Function      : Stepper::PrintIter
// Purpose       : Print the results of an inner iteration call.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
void Stepper::PrintIter(const bool converged) const
{
  if (PrintProc>0)
    {
      cout <<"~~~~~~ Step "<< StepNumber << " (Param = "<<iteration_method->Get_lambdafinal() <<"): ";
      
      if (converged) cout<<"*Converged* ";
      else           cout<<"##Failed## to Converge ";
      cout << endl;
      
      /*
      int itr = Problem_Parameters->get_MaxInnerIts();
      cout <<"in "<< itr <<" Iteration";
      if (itr != 1) cout << "s";
      cout << "  ~~~~~~" << endl;
      */
    }
}

//-----------------------------------------------------------------
// Function      : Stepper::StepParamAndPredict
// Purpose       : Step the parameter and make a guess of the solution
//                 at the new parameter value.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
void Stepper::StepParamAndPredict()
{
  StepNumber++;
  
  // Adjust Step Size Here
  StepSize *= 1.0;
  
  iteration_method->Predictor(StepSize, PrevStepSize);
  PrevStepSize = StepSize;
}

//-----------------------------------------------------------------
// Function      : Stepper::Run
// Purpose       : Solve the continuation problem
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
void Stepper::Run()
{
  if (PrintProc>0) cout << "\n:::Starting CAPO Run:::" << endl;

  iteration_method->Initialize();

  do {
    bool converged;

    PrintStart();

    // Perform inner iteration
    converged = iteration_method->InnerIteration();

    PrintIter(converged);

    if (converged) {

      // Save converged solution

      //! J.Simonis June10 Must do something with this output function!!
      //(*Solver).get_xfinal().Output(Step, (*Solver).get_lambdafinal());
      cout << "Period for parameter " << iteration_method->Get_lambdafinal() << " is " << iteration_method->Get_Tfinal() << endl;


      // Some Solvers require work in between inner solves...
      iteration_method->InnerFunctions();

      // Check for basis decrease and maintain accuracy with subspace iter
      //if ((Step+1) % updateBasisFreq == 0)  rpm->UpdateBasis();

      // Step in parameter call predictor 
      StepParamAndPredict();
    }
    else {
      // The inner iterations have failed.  Each Algorithm should
      // have a function too deal with an inner iteration failure.
      cout << "Failed to solve for parameter value "  << iteration_method->Get_lambdafinal() << endl;
      iteration_method->IterationFailure();


      //Increase m this is for Rpm method
      //rpm->Increase_m();
      // If Newton-Picard fails, for now, just output that it failed and stop.
      //cout << "Failure of Newton-Picard Algorithm. " << endl;
      //Step = maxStep+2; // Hopefully this will throw us out when a failure
      // is encountered.
    }

  } while (!Done());

  iteration_method->Finish();

  if (PrintProc>0) cout << "\n:::Ending CAPO Run:::\n" << endl;

}
