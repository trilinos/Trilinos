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

#ifndef NOX_LINESEARCH_COMMON_H
#define NOX_LINESEARCH_COMMON_H

#include "NOX_LineSearch_Generic.H"
#include "NOX_Common.H"

// Forward declaration
namespace NOX {
class Utils;
}

namespace NOX {

namespace LineSearch {

//! %Common line search utilites
/*!  
  This is an abstract class that other line searches can derive
  from to get access to some common functionality.


  <B> Common Line Search Data </B> <br>

Line searches have a common set of data that should be tracked in saved in the parameter list for the users.  This class provides a set of common data declarations and a routine to print them.  A sublist for output parameters will be created called "Output" in the parameter list used to instantiate or reset the class.  Valid output keys for the parameter list are:

  - "Total Number of Line Search Calls" - Total number of calls to the compute() method of this line search.

  - "Total Number of Non-trivial Line Searches" - The total number of steps that could not directly take a full step and meet the required "Convergence Criteria" (i.e. The line search had to reduce the step length using inner iteration calculations over iterate \f$ k \f$). 

  - "Total Number of Failed Line Searches" - total number of line searches that failed and used a recovery step.

  - "Total Number of Line Search Inner Iterations" - total number of inner iterations \f$ k \f$ performed by this object.

  <B> Flexible step acceptance capability </B> <br>

Line searches typically allow steps to be accepted only if they meet sufficient decrease and curvature criteria, eg Armijo-Goldstein conditions.  It can sometimes happen that such a requirement leads to very many nonlinear iterations being taken before the solution falls into a much steeper basin of attraction.  The following two parameters allow the user to relax the curvature and decrease conditions.  !! WARNING !! Convergence proofs and other guarantees associated with the various Line Searches no longer pertain when using these options.  

  Given the following information at %solver step \f$ n \f$:
  -# Old solution group: \f$ x_{\rm old} \f$
  -# Direction vector: \f$ s \f$
  -# step length \f$ \lambda \f$ 

  - "Allowed Relative Increase" - the allowed value \f$ k > 0 \f$ for which 
  \f$
  \frac{||F (x_{\rm old} + \lambda s)||_2^2}{||F (x_{\rm old})||_2^2} \le k
  \f$

  - "Maximum Increase Steps" - maximum number of outer iterations for which the preceding acceptability criterion is allowed.  After this number of outer iterations, the original curvature and decrease conditions must be met for a step to be accepted.
*/

class Common : public Generic {

public:

  //! Default constructor
  Common(const NOX::Utils& u, NOX::Parameter::List& lineSearchParams);

  //! Destructor
  virtual ~Common();

  //! Reset the Common data components.
  virtual bool reset(NOX::Parameter::List& lineSearchParams);

protected:

  //! Print out step information
  /*!
    NEED TO INSERT EXMAPLE OUTPUT HERE
  */
  void printStep(int n, double step, double oldf, double newf, const string s = "") const;

  //! Compute the inner product of the given direction and the gradient associated with the given group.
  /*!
    Calculates and returns 
    \f[
    \zeta = d^T \nabla f(x).
    \f]

    Here \f$d\f$ represents the input parameter \c dir and \f$\nabla
    f(x)\f$ is the gradient associated with the given group.
  */
  double computeSlope(const NOX::Abstract::Vector& dir, const NOX::Abstract::Group& grp);

  //! This is a variant of the computeSlope() method above optimized to work with out having to compute an explicit Jacobian.  
  /*!
    Calculates and returns 
    \f[
    \zeta = d^T \nabla f(x) = d^TJ^TF
    \f]

    Here \f$d\f$ represents the input parameter \c dir \f$\nabla
    f(x)\f$ is the gradient associated with the given group (for nonlinear solves this equates to \f$ J^TF \f$ where \f$ J \f$ is the Jacobian and \f$ F \f$ is the original nonlinear function).

    We can rewrite this equation as:

    \f[ d^TJ^TF = F^TJd \f]
    
    which allows us to use directional derivatives to estimate \f$ J^TF \f$:

    \f[ F^TJd = F^T \frac{F(x + \eta d) - F(x)}{\eta} \f]

    This may allow for faster computations of the slope if the Jacobian is expensive to evaluate.

where \f$\eta\f$ is a scalar perturbation calculated by:

\f[ \eta = \lambda * (\lambda + \frac{\| x\|}{\| d\|} ) \f]

\f$ \lambda \f$ is a constant fixed at 1.0e-6.

  */
  double computeSlopeWithOutJac(const NOX::Abstract::Vector& dir, const NOX::Abstract::Group& grp);

  //! Reset the common line search information 
  void resetInfo();

  //! Sets the common line search data in an "Output" sublist of the "Line Search" sublist. 
  bool setCommonDataValues();

  //! Provides a check to see if a step leading to an increase in the residual is allowed
  bool isIncreaseAllowed(double, double, int) const;

protected:

  //! Printing utilities
  const NOX::Utils& utils;

  /** @name Common Line Search Data 
  
  All line searches in general should report the following information.  We save a pointer to the parameter list associated with the line search and set an output sublist with each of the following parameters.
  */ 
  //@{

  //! Pointer to the input parameter list.  We need this to create an "Output" sublist to store output parameters.  
  NOX::Parameter::List* paramsPtr;

  //! True is we should allow residual increases constrained by relIncrease
  bool allowIncrease;

  //! Allowable relative increase for residual norm
  double relIncrease;

  //! Maximum number relative increase steps allowed
  int numAllowed;

  //! Counter for the total number of times the compute() method was called for this object. 
  int totalNumLineSearchCalls;

  //! Counter for number of times the line search was actually needed.  This is the number of times the compute() method was called for this object where a step length of 1.0 could not be used.
  int totalNumNonTrivialLineSearches;

  //! Counter for the number of times the compute() method was called where the line search algorithm failed.
  int totalNumFailedLineSearches;

  //! Counter for the total number of line search inner iterations performed.
  int totalNumIterations;

  //@}

private:

  //! A vector that *may* be used in computeSlope()
  NOX::Abstract::Vector* vecPtr;

  //! An extra group that will only be allocated if computeSlopeWithOutJac() is called.
  NOX::Abstract::Group* grpPtr;

};
} // namespace LineSearch
} // namespace NOX

#endif
