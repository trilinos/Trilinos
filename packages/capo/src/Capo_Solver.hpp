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

#ifndef Capo_SOLVER_H
#define Capo_SOLVER_H

/************************************************************ 
File:      Capo_Solver.hpp
Purpose:   The Solver object contains all the necessary functionality
           to be used by the stepper (reword this at somepoint)
Date:      6-13-05
Author:    Joseph Simonis
**************************************************************/

/**** Includes ****/
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Capo_Integrator.hpp"
#include "Capo_Parameter_List.hpp"


namespace CAPO {
  
  typedef double Scalar;

  class Integrator;
  class Parameter_List;

  /** \brief 
      This is the virtual Solver class.  Any future additional
      solvers added to Capo must be implementations of this 
      class.
  */
  class Solver {
  public:

    //! Constructor
    /*!
      \note A Parameter list contains values needed within the
      algorithm.  The Integrator is the only access to the problem
      being solved. x0, T0, and lambda0 are the initial values
      for the state, time, and parameter.
     */    
    Solver(Teuchos::RefCountPtr<Parameter_List> ParamList, \
	   Teuchos::RefCountPtr<Integrator> App_Int, \
	   Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x0, \
	   double lambda0, double T0);

    //! Default Constructor
    Solver() {}

    //! Destructor
    virtual ~Solver() {}

    //! Member Functions

    /*!
      Setup the solver.
    */    
    virtual void Initialize() = 0;

    /*!
      Find a fixed point of the system using the current
      values of x, T, and lambda as initial guesses.
    */    
    virtual bool InnerIteration() = 0 ;

    /*!
      If any computations must be done between solves
      along the continuation branch (e.g. updating
      subspace basis) this function should perform
      them.  Currently in both NPGS and RPM these
      are empty functions.
    */    
    virtual void InnerFunctions() = 0;

    /*!
      In the case where an InnerIteration call fails
      to find a fixed point, the problem must be reset
      or the code must exit cleanly.
    */    
    virtual void IterationFailure() = 0;

    /*!
      This function should perform any tasks desired to 
      occur upon finishing the entire continuation problem.
      For example printing out a final solution or printing
      the final basis.
    */    
    virtual void Finish() = 0;

    /*!
      The predictor increments the parameter and determines
      a suitable initial guess of x and T for the next
      solve.
    */    
    virtual void Predictor() = 0;

    /*!
      Returns the value of xfinal.  
    */    
    virtual Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& Get_xfinal() = 0;

    /*!
      Returns the value of lambdafinal.  
    */    
    virtual double& Get_lambdafinal() = 0;

    /*!
      Returns the value of Tfinal.  
    */    
    virtual double& Get_Tfinal() = 0;

    
  protected:
    //! The starting guess for x at the current continuation step.
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xcurrent;

    //! The final solution x at the current continuation step.
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xfinal;

    //! The final solution x at the previous continuation step.
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xprevious;

    //! The starting guess for lambda at the current continuation step.
    double lambdacurrent;

    //! The final lambda at the current continuation step.
    double lambdafinal;

    //! The final lambda at the previous continuation step.
    double lambdaprevious;

    //! The starting guess for T at the current continuation step.
    double Tcurrent;

    //! The final T at the current continuation step.
    double Tfinal;

    //! The final T at the previous continuation step.
    double Tprevious;

    //! A counter for the number of inner iterations performed
    //! for a particular continuation step.
    int iter;

    //! The system integrator.
    Teuchos::RefCountPtr<Integrator> App_Integrator;

    //! The list of solver parameters.
    Teuchos::RefCountPtr<Parameter_List> SolveParameters;

    //! A bool to indicate if we are on the first Continutaion step.
    bool First_Continuation_Step;

  };

}

#endif
