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

#ifndef Capo_SOLVER_NPGS_H
#define Capo_SOLVER_NPGS_H

/************************************************************ 
File:      Capo_Npgs.hpp
Purpose:   The Newton--Picard Gauss--Seidel Solver method
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

  class Npgs : public virtual Solver
  {
  public:
    //! Constructor
    Npgs(Teuchos::RefCountPtr<Parameter_List> ParamList, \
	 Teuchos::RefCountPtr<Integrator> App_Int, \
	 Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x0,\
	 double lambda0, double T0);
    
    Npgs();
    //! Destructor
    ~Npgs();

    //! Member Functions
    void Initialize();
    
    bool InnerIteration();
    
    void InnerFunctions();
    
    void IterationFailure();
    
    void Finish();
    
    void Predictor();
    
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& Get_xfinal();
    
    double& Get_lambdafinal();
    
    double& Get_Tfinal();

  private:
    Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Ve;

    void Orthonormalize(Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Basis, int Number_of_Columns);

    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > MatVec(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > y);

    Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > MatVecs(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Y);

    void SchurDecomp(Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Se, Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Re);

    void Print(Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Printme);
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
