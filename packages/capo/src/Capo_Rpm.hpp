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
    Rpm(const Teuchos::RefCountPtr<Parameter_List>& ParamList, \
	const Teuchos::RefCountPtr<Integrator>& App_Int, \
	const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x0,\
	double lambda0, double T0);

    Rpm() {}

    //! Destructor
    virtual ~Rpm() {}

    //! Member Functions
    virtual void Initialize();
    
    virtual bool InnerIteration();
    
    virtual void InnerFunctions();
    
    virtual void IterationFailure();
    
    virtual void Finish();
    
    virtual void Predictor(double& StepSize, double& PrevStepSize);
    
    virtual Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& Get_xfinal();
    
    virtual double& Get_lambdafinal();
    
    virtual double& Get_Tfinal();

  private:

    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > finit;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xinit;

    Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Ve;

    int Unstable_Basis_Size;

    void Orthonormalize(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Basis, int Number_of_Columns);

    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > MatVec(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& y);

    Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > MatVecs(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Y);


    bool Converged(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x,
		   const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& y);

    void SubspaceIterations(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se, Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& We, const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re);

    void Calculatedq(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp,
		     const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dq, 
		     const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& r);

    bool ComputeVp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se,
		   const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp);
    bool UpdateVe(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& We,
		  const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se);

    bool Calculatedp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp,
		     const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dq,
		     const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dp,
		     const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re,
		     const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& v,
		     const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& finit,
		     const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& r,
		     double& deltaT);
    bool dphi_dt(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& f);

    bool Solve_Linear(double *Mat,double *rhs, bool resolve, int m, int nrhs);
    void SchurDecomp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se, 
		     const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re);

    void Print(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Printme);
    bool dphi_dlambda(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& f);
    bool ShermanMorrison(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp,
			 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dq,
			 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dp,
			 const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re,
			 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& v,
			 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& finit,
			 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& r,
			 double& deltaT, double& deltalambda);
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xstep;
    double lambdastep;
    double Tstep;
    
  };
}

#endif
