// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef SOLVER_NONLINEARCG_H
#define SOLVER_NONLINEARCG_H

#include "NOX_Direction_Generic.H"	    // base class

namespace NOX {
namespace Direction {

//! Nonlinear Conjugate Gradient method with Linesearch
/*!
  
The following parameters are valid for this solver:

- "Direction" -  Allows descent direction to be computed
  as either the negative of the residual using "Richardson" \f${\bf
  d}=-{\bf r} \f$ or as the Steepest Descent direction (default) using
  "Steepest Descent" \f${\bf d}= -{\bf J}^T{\bf r}\f$.  The latter assumes
  the objective function to be minimized is \f$f\equiv\frac{1}{2}{\bf
  r}^T {\bf r}\f$.

- "Diagonal Precondition" - Either "On" or "Off" (default), determines
  whether or not to scale the Descent Vector \f${\bf d}\f$ by the inverse of 
  the diagonal of the current Jacobian.  When "Off", \f${\bf z} = {\bf d}\f$, 
  and when "On" \f${\bf z} = \left({\rm diag} {\bf J}\right)^{-1} 
  {\bf r}\f$.  This means of preconditioning is similar to that described in
  Sandia Report, SAND81-0084J, J. H. Biffle, "The Solution of Quasi-Static
  Nonlinear Mechanics Problems by the Nonlinear Conjugate Gradient Method"

- "Orthogonalize" -  Determines which expression to use for the 
  orthogonalization parameter \f$\beta\f$ which is used to obtain a new search
  direction by \f${\bf s}_{k+1} = {\bf z}_{k+1}+\beta{\bf s}_k\f$.  Currently,
  two choices are available: "Fletcher-Reeves" (default) for which 
  \f$\beta = \frac{{{\bf d}_{k+1}}^T{\bf d}_{k+1}}{{{\bf d}_k}^T{\bf d}_k}\f$
  or "Polak-Ribiere" for which 
  \f$\beta = \frac{{{\bf d}_{k+1}}^T({\bf d}_{k+1}-{\bf d}_k)}{{{\bf d}_k}^T
  {\bf d}_k}\f$.  For more details, see chapter 5 of Nocedal & Wright, 
  "Numerical Optimization", Springer-Verlag, New York, 1999.

-  "Restart Frequency" - Number of nonlinear iterations between restarts
   (default = 10).  Restart corresponds to setting \f$\beta = 0\f$.

- "Output Frequency" - Number of nonlinear iterations between outputs
  (default = 10)

- "Line Search" - Sublist of the line search parameters, passed to the
  NOX::Linesearch::Manager constructor. Defaults to an empty list.

*/

class NonlinearCG : public Generic {

public:

  //! %NonlinearCG-method constructor. 
  NonlinearCG(Parameter::List& params);

  //! Destructor
  virtual ~NonlinearCG();

  virtual bool reset(Parameter::List& p);
//  virtual StatusTest::StatusType getStatus();
  virtual bool compute(Abstract::Vector& dir, Abstract::Group& grp, 
                       const Solver::Generic& solver);
//  virtual StatusTest::StatusType solve();
//  virtual const Abstract::Group& getSolutionGroup() const;
//  virtual const Abstract::Group& getPreviousSolutionGroup() const;
//  virtual int getNumIterations() const;
//  virtual const Parameter::List& getOutputParameters() const;

protected:
  
  //! Print out initialization information and calcuation the RHS.
  //virtual void init();

  //! Prints the current iteration information.
  //virtual void printUpdate();

protected:
  
  //! Current solution.
  Abstract::Group* solnptr;		

  Abstract::Vector* tmpVecPtr;

  //! Previous solution pointer. 
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  const Abstract::Group* oldsolnptr;	
  //! Previous solution reference.
//  Abstract::Group& oldSoln;	

  //! Current search direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
//  Abstract::Vector* dirptr;
  //! Current search direction.reference.
//  Abstract::Vector& dir;

  //! Previous search direction.pointer.
  /*! The current search direction is found by computing the steepest
      descent direction and then orthogonalizing using this previous
      search direction.  We have both a pointer and a reference because
      we need to create a DERIVED object and then want to have a reference
      to it. */
  Abstract::Vector* oldDirPtr;

  //! Pointer to previous descent direction.
  /*! This pointer points into the previous solution Abstract::Group to the 
      descent direction based on the choice of "Richardson" or "Steepest
      Descent" determined by Orthogonalization in the Parameter::List. */
  Abstract::Vector* oldDescentDirPtr;  

  //! Pointer to current preconditioned direction.
  /*! Currently, preconditioning is done using only Jacobi (diagonal)
      scaling as described above.  If preconditioning is turned off, this
      vector is equal to the direction vector. */
  Abstract::Vector* preconditionedDirptr;  

  //! Pointer to previous descent direction preconditioned with current Jacobian.
  /*! This pointer is to the vector needed for the preconditioned Polak-Ribiere
      calculation of /f$\beta/f$ as described above.  Currently, 
      preconditioning is done using only Jacobi (diagonal) scaling.  
      If preconditioning is turned off, this vector is equal to the 
      descent direction from the previous nonlinear iteration. */
  Abstract::Vector* preconditionedOldDirptr;  

  //! Vector pointer to difference between current and previous descent directions.
  /*! This pointer is to the vector needed for the preconditioned Polak-Ribiere
      calculation of /f$\beta/f$ as described above.  
      If preconditioning is turned off, this vector is equal to the 
      difference between the current and previous descent directions. */
  Abstract::Vector* diffVecPtr;  

  //! Stopping test.
  //StatusTest::Generic* testptr;		

  //! Input parameters.
  Parameter::List& iparams;	

  //! Output paramters.
  //mutable Parameter::List oparams; 

  //! Orthogonalization parameter
  double beta;			

  //! Linesearch. 
  //LineSearch::Manager linesearch; 

  //! Current step.
  //double step;			

  //! Number of nonlinear iterations.
  int niter;                    

  //! Restart frequency --> use Steepest descent direction by setting beta = 0
  int restartFrequency;                    

  //! Output frequency 
  //int outputFrequency;                    

  //! %Status of nonlinear solver.
  //StatusTest::StatusType status;

};
} // namespace Direction
} // namespace NOX

#endif

