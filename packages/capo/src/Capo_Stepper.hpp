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

#ifndef Capo_STEPPER_H
#define Capo_STEPPER_H

/************************************************************ 
File:      Capo_Stepper.hpp
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

namespace CAPO {
  
  class Integrator;
  class Parameter_List;

  /** \brief
     This class provides the main driver of the Capo program.  Once 
     all of the parameters are set, the integrator is built, and
     a solver is chosen.  This class runs the problem by calling
     the predictor function of the solver, and then the corrector
     alternatively until the maximum number of steps has been taken.
  */
  class Stepper {

  public:
    //! Constructor
    Stepper(Teuchos::RefCountPtr<Parameter_List> PL, \
	    Teuchos::RefCountPtr<Solver> App_Solver);

    //! Destructor
    ~Stepper() {};

    /*!
      This function executes the continuation solve.  It fixes the parameter,
      calls the solver, increments the parameter and repeats until the maximum
      number of steps has been achieved.
    */
    void Run();
    
    /*!
      Increment the continuation parameter and predict solution values of
      the state (and period if periodic problem) at the new parameter.

      \note The predicted values are obtained via the Solver algorithms
      Predictor function.
    */
    void StepParamAndPredict();
    
    /*!
      Check to see if the maximum number of continuation steps has been
      reached.
    */
    bool Done() const;
    
    /*!
      Indicate the beginning of each continuation step.
    */
    void PrintStart() const;
    
    /*!
      Print the results of the solve at a particular continuation
      step.
    */
    void PrintIter(const bool converged) const;
    

  private:
    //! Counter to keep track of what continuation step we are on.
    int StepNumber;

    //! The number of continuation steps to take.
    int MaxSteps;

    //! If the PrintProc flag is set >0, then additional information 
    //! will be printed throughout the run.
    int PrintProc;

    //! A list of parameter values.
    Teuchos::RefCountPtr<Parameter_List> Problem_Parameters;

    //! The solver to be used at each continuation step.
    Teuchos::RefCountPtr<Solver> iteration_method;

  };

}

#endif
