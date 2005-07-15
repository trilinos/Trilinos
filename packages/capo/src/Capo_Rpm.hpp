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

    /*!
      Setup the solver.
    */    
    virtual void Initialize();
    
    /*!
      Find a fixed point of the system using the current
      values of x, T, and lambda as initial guesses.
    */    
    virtual bool InnerIteration();
    
    /*!
      If any computations must be done between solves
      along the continuation branch (e.g. updating
      subspace basis) this function should perform
      them.  Currently this is an empty function.
    */    
    virtual void InnerFunctions();
    
    /*!
      In the case where an InnerIteration call fails
      to find a fixed point, the problem must be reset
      or the code must exit cleanly.

      \note The function resets the problem to start
      the current continuation step again, but halves
      the parameter step and calculates a different 
      starting guess based on the new parameter.  The 
      hope is that if the previous continuation step
      succeeded, we should be able to converge if 
      our step in lambda isn't too large.
    */    
    virtual void IterationFailure();
    
    /*!
      Prints out that the continuation has ended
      successfully.
    */    
    virtual void Finish();
    
    /*!
      The predictor increments the parameter and determines
      a suitable initial guess of x and T for the next
      solve.
    */    
    virtual void Predictor();
    
    /*!
      Returns the value of xfinal.  
    */    
    virtual Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& Get_xfinal();
    
    /*!
      Returns the value of lambdafinal.  
    */    
    virtual double& Get_lambdafinal();
    
    /*!
      Returns the value of Tfinal.  
    */    
    virtual double& Get_Tfinal();

  private:

    /*!
      The current (hard coded) phase condition is given by:
      f(x(0)^{(0)},lambda^{(0)})^T(x(0)-x(0)^{(0)})=0
      which requires that the initial f be stored throughout
      the algorithm.
    */
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > finit;

    /*!
      The current (hard coded) phase condition is given by:
      f(x(0)^{(0)},lambda^{(0)})^T(x(0)-x(0)^{(0)})=0
      which requires that the initial x(0)^{(0)} be stored 
      throughout the algorithm.
    */
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xinit;

    /*!
      Ve stores the schur vectors associated with eigenvalues
      of phi(x,T,lambda) with magnitude greater than a 
      specified value rho.  It also keeps a few additional
      vectors.  This basis can be of size up to 30.  
    */
    Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > Ve;

    /*!
      The number of eigenvalues with magnitude greater 
      than rho.
    */
    int Unstable_Basis_Size;

    /*!
      Orthonormalize the first Number_of_Columns of
      the multivector Basis.
    */
    void Orthonormalize(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Basis, 
			int Number_of_Columns);

    /*!
      Calculate a Matrix Vector product using finite differences:
      Mv=(1/eps)*(phi(x+eps*y,T,lambda)-phi(x,T,lambda))
    */
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > MatVec(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& y);

    /*!
      Perform Matrix Vector products on multiple vectors.
      This function just calls MatVec on each column of the matrix Y.
    */
    Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> > MatVecs(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Y);

    /*!
      Returns true if (||x-y||_2)/n<tol
    */
    bool Converged(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x,
		   const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& y);

    /*!
      Algorithm 5.3 Subspace Iteration with Projection from 
      "Numerical Methods for Large Eigenvalue Problems" by
      Youcef Saad.  The algorithm finds approximations to
      the dominant Schur vectors of a matrix.

      \note Se returns holding the dominant Schur vectors, Re
      contains the multiplication factors from the QR
      decomposition.
    */
    void SubspaceIterations(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se, Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& We, const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re);

    /*!
      Calculate dq using a fixed point iteration.
    */
    void Calculatedq(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp,
		     const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dq, 
		     const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& r);

    /*!
      Compute Vp by multiplying the first "Unstable_Basis_Size" columns of Ve
      by the upper left corner of Se.  See Algorithm 3.2 NPGS(2) in 
      "An Adaptive Newton Picard Algorithm with Subspace Iteration for 
      Computing Periodic Solutions" by Lust, Roose, Spence, Champneys.
    */
    bool ComputeVp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se,
		   const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp);
    /*!
      Ve = We*Se
    */
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
    /*!
      Use finite differences to calculate \partial \phi / \partial t.
      (1/eps)*(\phi(x,T+eps,lambda)-\phi(x,T,lambda))
    */
    bool dphi_dt(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& f);

    /*!
      Solve the linear system:
      Ax=b.
      Mat is the matrix A, rhs is b upon entrence, x upon exit.
      resolve = false upon entering, m is the length of rhs.
      nrhs is the number of right hand sides to be solved.
    */
    bool Solve_Linear(double *Mat,double *rhs, bool resolve, int m, int nrhs);

    /*!
      Perform a Schur Decomposition.  This uses the lapack routine dgees.
    */
    void SchurDecomp(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Se, 
		     const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re);

    /*!
      Print a MultiVector.
    */
    void Print(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Printme);

    /*!
      Use finite differences to calculate \partial \phi / \partial lambda.
      (1/eps)*(\phi(x,T,lambda+eps)-\phi(x,T,lambda))
    */
    bool dphi_dlambda(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& f);
    bool ShermanMorrison(const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Vp,
			 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dq,
			 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& dp,
			 const Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >& Re,
			 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& v,
			 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& finit,
			 const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& r,
			 double& deltaT, double& deltalambda);

    //! The step taken from the previous solution to the guess at this
    //! continuation step.
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xstep;

    //! The increment in lambda from the previous lambda to the guess
    //! at this continuation step.
    double lambdastep;

    //! The increment in T from the previous T to the guess
    //! at this continuation step.
    double Tstep;
    
  };
}

#endif
