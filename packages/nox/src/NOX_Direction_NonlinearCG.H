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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef SOLVER_NONLINEARCG_H
#define SOLVER_NONLINEARCG_H

#include "NOX_Direction_Generic.H"        // base class

// Forward declaration
namespace NOX {
  class Utils;
  class GlobalData;
}

namespace NOX {
namespace Direction {

//! Calculates a search direction using the Nonlinear Conjugate Gradient method
/*!

Calculates the direction
\f[
d = - M^{-1}(x) F(x) + \beta d_{prev}
\f]

  where \f$ M \f$ is a preconditioner and \f$ \beta \f$ is an
  orthogonalization parameter which can be computed in various ways
  (see below), and \f$ d_{prev} \f$ is the search direction from the
  previous nonlinear iteration.

  This method provides a generalization of %Linear CG to nonlinear
  problems.  It does this by computing a search direction using
  an expression analogous to that of %Linear CG.  The negative of
  the current residual vector, \f$ F(x) \f$ is taken, allowed to be
  preconditioned, and then orthogonalized against the previous search
  direction.  This direction can sometimes be used successfully with the
  various choices provided in NOX::Linesearch but is intended to be used
  with NOX::Linesearch::NonlinearCG.  In fact, the expected convergence
  behavior of linear problems can only be achieved in this way.

  To use this direction, specify that the "Method" is "NonlinearCG" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).

  The following options may be specified in the "Nonlinear CG" sublist
  of the "Direction" sublist of the solver parameters.

<ul>
<li> "Orthogonalize" can be either of:

    <ul>
    <li>"Fletcher-Reeves" [default] - \f$ \beta =  \frac{F(x)^T M^{-1}(x)
    F(x)}{F(x_{prev})^T M^{-1}(x_{prev}) F(x_{prev})}\f$

    <li>"Polak-Ribiere" - \f$ \beta =  max \left\{ \beta^{PR}, 0 \right\} \f$ ,

    where
    \f$ \beta^{PR} = \frac{F(x)^T \left[M^{-1}(x) F(x)
    - M^{-1}(x_{prev}) F(x_{prev})\right]}{F(x_{prev})^T M^{-1}(x_{prev})
    F(x_{prev})}\f$

    </ul>

  These comprise the two most popular choices for orthogonalization.
  Both reduce to the linear result for linear problems.  "Polak-Ribiere"
  provides an implied restart capability by setting \f$ \beta = 0 \f$
  anytime the computed value is less than zero.

<li> "Precondition" can be either "On" or "Off" [default]:
  determines whether or not to compute and apply preconditioner \f$ M \f$.
  If "Off" is selected, no preconditioner is computed and the behavior is
  equivalent to \f$ M = I \f$ where \f$ I \f$ is the identity matrix.
  If "On", \f$ M \f$ is computed and applied as determined by the
  underlying implementation of the "applyRightPreconditioning" method
  in the Group.


<li> "Restart Frequency" - An integer specification of the number of nonlinear
  iterations between restarts [default = 10].  Restart corresponds to
  setting \f$\beta = 0\f$.  A good heuristic is to limit this value to
  the number of problem degrees of freedom.  Setting this value to 1
  forces \f$ \beta = 0 \f$ for every nonlinear iteration which corresponds
  to suppressing orthogonalization against the previous search direction.

</ul>

<b>References</b>

  <ul>

  Additional information about both linear and nonlinear conjugate gradient
  methods can be found in Chapter 5 of:

  <li>Nocedal & Wright, "Numerical Optimization", Springer-Verlag,
  New York, 1999.</li>
  Though presented within the context of nonlinear optimization, the connection
  to nonlinear systems of equations is made by the correspondence
  \f$ \nabla f(x) \leftrightarrow F(x) \f$ (cf Algorithm 5.4).

  Another useful useful reference is:

  <li>Jonathan Richard Shewchuk,
  <A HREF="http://www-2.cs.cmu.edu/~jrs/jrspapers.html"/>
  "An Introduction to the Conjugate Gradient Method Without the Agonizing
  Pain</A>," 1994.</li> Chapter 14 provides a summary of issues in generalizing
  linear CG to the nonlinear case.  Correspondence to NOX notation
  is made by the equivalence \f$ r \leftrightarrow f' \leftrightarrow F(x) \f$
  (cd Section 14.1).

  </ul>


*/

class NonlinearCG : public Generic {

public:

  //! Constructor.
  NonlinearCG(const Teuchos::RCP<NOX::GlobalData>& gd,
          Teuchos::ParameterList& params);

  //! Destructor.
  virtual ~NonlinearCG();

  //! derived
  virtual bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
             Teuchos::ParameterList& p);

  //! derived
  virtual bool compute(Abstract::Vector& dir, Abstract::Group& grp,
                       const Solver::Generic& solver);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::LineSearchBased& solver);

private:

  //! Global data pointer.  Keep this so the parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Printing Utils
  Teuchos::RCP<NOX::Utils> utils;

  //! Previous solution pointer.
  const Abstract::Group* oldSolnPtr;

  //! Temporary vector used in applying preconditioner.
  Teuchos::RCP<Abstract::Vector> tmpVecPtr;

  //! Previous search direction pointer.
  /*! The current search direction is found by orthogonalizing against the
      previous search direction.  This vector provides storage for the
      old search direction. */
  Teuchos::RCP<Abstract::Vector> oldDirPtr;

  //! Pointer to previous descent direction.
  /*! This vector stores the old descent direction which was orthogonalized
      to determine the old search direction.  It is needed in the numerator
      of the expression for Polak-Ribiere \f$ \beta \f$. */
  Teuchos::RCP<Abstract::Vector> oldDescentDirPtr;

  //! Temporary vector used in computing the numerator of Polak-Ribiere \f$ \beta \f$.
  Teuchos::RCP<Abstract::Vector> diffVecPtr;

  //! Input parameters.
  Teuchos::ParameterList* paramsPtr;

  //! Orthogonalization parameter \f$ \beta \f$.
  double beta;

  //! Number of nonlinear iterations.
  int niter;

  //! Restart frequency.
  int restartFrequency;

  //! Boolean flag indicating whether or not to use preconditioning
  bool doPrecondition;

  //! Boolean flag indicating if beta is computed using Polak-Ribiere
  bool usePRbeta;

};
} // namespace Direction
} // namespace NOX

#endif
