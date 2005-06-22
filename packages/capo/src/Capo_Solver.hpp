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

  /**** Solver Class ****/
  class Solver {
  public:
    //! Constructor
    
    Solver(Teuchos::RefCountPtr<Parameter_List> ParamList, \
	   Teuchos::RefCountPtr<Integrator> App_Int, \
	   Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x0, \
	   double lambda0, double T0);
    
    Solver() {}

    //! Destructor
    virtual ~Solver() {}

    //! Member Functions
    
    virtual void Initialize() = 0;
    
    virtual bool InnerIteration() = 0 ;

    virtual void InnerFunctions() = 0;

    virtual void IterationFailure() = 0;

    virtual void Finish() = 0;

    virtual void Predictor(double& StepSize, double& PrevStepSize) = 0;

    virtual Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& Get_xfinal() = 0;

    virtual double& Get_lambdafinal() = 0;

    virtual double& Get_Tfinal() = 0;
    
  protected:
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xcurrent;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xfinal;

    double lambdacurrent;
    double lambdafinal;

    double Tcurrent;
    double Tfinal;

    int iter;
    Teuchos::RefCountPtr<Integrator> App_Integrator;
    Teuchos::RefCountPtr<Parameter_List> SolveParameters;
  };
};

#endif
