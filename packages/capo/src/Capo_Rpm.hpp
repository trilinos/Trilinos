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

#ifndef Capo_SOLVER_RPM_H
#define Capo_SOLVER_RPM_H

/************************************************************ 
File:      Capo_Rpm.hpp
Purpose:   The Recursive Projection Method
Date:      6-13-05
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
  
  class Rpm : public CAPO::Solver
  {
  public:
    //! Constructor
    Rpm(Teuchos::RefCountPtr<Parameter_List> ParamList, \
	Teuchos::RefCountPtr<Integrator> App_Int, \
	Teuchos::RefCountPtr<Thyra::VectorBase(Scalar)> x0,\
	double lambda0, double T0);

    //! Destructor
    ~Rpm();

    //! Member Functions
    void Initialize();
    
    void InnerIteration();
    
    void InnerFunctions();
    
    void IterationFailure();
    
    void Finish();
    
    void Predictor();
    
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& Get_xfinal();
    
    double& Get_lambdafinal();
    
    double& Get_Tfinal();

  private:
    
