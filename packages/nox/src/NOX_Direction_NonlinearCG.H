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

//! Calculates a search direction using the Nonlinear Conjugate Gradient method
/*!

Calculates the direction
\f[
s = - M^{-1}(x) F(x) + \beta s_{prev}
\f]

  where \f$ M \f$ is a preconditioner and \f$ \beta \f$ is an orthogonalization
  parameter which can be computed in various ways (see below).

  To use this direction, specify that the "Method" is "NonlinearCG" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).

  The following options may be specified in the "Direction" sublist of
  the solver parameters.

<ul>
<li> "Orthogonalize" can be either of:

    <ul>
    <li>"Fletcher-Reeves" [default] - \f$ \beta =  \frac{F(x)^T M^{-1}(x)
    F(x)}{F(x_{prev})^T M^{-1}(x_{prev}) F(x_{prev})}\f$ 

    <li>"Polak-Ribiere" - \f$ \beta =  \frac{F(x)^T \left[M^{-1}(x) F(x)
    - M^{-1}(x_{prev}) F(x_{prev})\right]}{F(x_{prev})^T M^{-1}(x_{prev})
    F(x_{prev})}\f$

    </ul>

  For more details, see chapter 5 of Nocedal & Wright, "Numerical
  Optimization", Springer-Verlag, New York, 1999.

<li> "Precondition" can be either "On" or "Off" [default]:
  determines whether or not to compute and apply preconditioner \f$ M \f$.
  If "Off" is selected, no preconditioner is computed and the behavior is
  equivalent to \f$ M = I \f$ where \f$ I \f$ is the identity matrix.
  If "On", \f$ M \f$ is computed and applied as determined by the
  underlying implementation.  Currently, NOX::Epetra supports the inverse
  of the current Jacobian, ie \f$ M = {\bf J}^{-1} \f$.  NOX::Petsc
  supports the full suite provided in the user's Petsc installation and
  can be set by proper specification in the .petscrc file.


<li> "Restart Frequency" - An integer specification of the number of nonlinear 
  iterations between restarts [default = 10].  Restart corresponds to
  setting \f$\beta = 0\f$.  A good heuristic is to limit this value to
  the number of problem degrees of freedom.  Setting this value to 1
  forces \f$ \beta = 0 \f$ for every nonlinear iteration and corresponds
  to a variant of NOX::Direction::SteepestDescent in which \f$ \nabla f(x)
  \f$ is equivalent to \f$ F(x) \f$.

</ul>

*/

class NonlinearCG : public Generic {

public:

  //! Constructor. 
  NonlinearCG(Parameter::List& params);

  //! Destructor.
  virtual ~NonlinearCG();

  //! derived
  virtual bool reset(Parameter::List& p);

  //! derived
  virtual bool compute(Abstract::Vector& dir, Abstract::Group& grp, 
                       const Solver::Generic& solver);

private:
  
  //! Previous solution pointer. 
  const Abstract::Group* oldSolnPtr;	

  //! Temporary vector used in applying preconditioner.
  Abstract::Vector* tmpVecPtr;

  //! Previous search direction pointer.
  /*! The current search direction is found by orthogonalizing against the
      previous search direction.  This vector provides storage for the 
      old search direction. */
  Abstract::Vector* oldDirPtr;

  //! Pointer to previous descent direction.
  /*! This vector stores the old descent direction which was orthogonalized
      to determine the old search direction.  It is needed in the numerator
      of the expression for Polak-Ribiere \f$ \beta \f$. */
  Abstract::Vector* oldDescentDirPtr;  

  //! Temporary vector used in computing the numerator of Polak-Ribiere \f$ \beta \f$.
  Abstract::Vector* diffVecPtr;  

  //! Input parameters.
  Parameter::List* paramsPtr;	

  //! Orthogonalization parameter \f$ \beta \f$.
  double beta;			

  //! Number of nonlinear iterations.
  int niter;                    

  //! Restart frequency.
  int restartFrequency;                    

};
} // namespace Direction
} // namespace NOX

#endif

