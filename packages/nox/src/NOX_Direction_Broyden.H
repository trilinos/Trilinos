// $Id$ 
// $Source$ 

#ifdef WITH_PRERELEASE

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_DIRECTION_BROYDEN_H
#define NOX_DIRECTION_BROYDEN_H

#include "NOX_Common.H"		   

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Direction_Utils_InexactNewton.H" // class data element
#include "NOX_Parameter_List.H"	   // class data element 

// Forward declaration
namespace NOX {
  class Utils;
  namespace Solver {
    class LineSearchBased;
  }
}

namespace NOX {
namespace Direction {

//! %Broyden direction 
/*!  
  
  We will calculate a limited-memory %Broyden direction of the form

  &nbsp;&nbsp; 
  \f$d_k = -B_k^{-1} F_k.\f$ 

  Here \f$B_k\f$ is a limited-memory %Broyden approximation to the
  Jacobian of \f$F\f$ at \f$x_k\f$, and \f$F_k = F(x_k)\f$. It is
  based on apply Broyden updates to the Jacobian from some previous
  step.

  \note
  The %Broyden direction can only be used with NOX::Solver::LineSearchBased.
  It cannot be used with any other solver, include NOX::Solver::TrustRegionBased.

  <b>References</b>

  <ul>
  <li> C. T. Kelley, <em>Iterative Methods for %Linear and Nonlinear Equations</em>, SIAM, 1995.
  </ul>
  

  <b>Parameters</b>

  To use this direction, specify that the "Method" is "Broyden" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).

  In "Direction"/"Broyden":

  <ul>

  <li>"Restart Frequency" - How often the Jacobian should be
      refreshed. A value of 5, for example, means that the Jacobian
      should be updated every 5 iterations. Defaults to 10.

  <li>"Max Convergence Rate" - Maximum convergence rate allowed when
  reusing the Jacobian.  The Jacobian will be refreshed if the
  convergence rate, \f$ \alpha \f$, is larger than this value.  The
  convergence rate is calculated by \f$ \alpha = \frac{\| F_k \| }{\|
  F_{k-1} \|} \f$ where F is the nonlinear residual and \f$ k \f$ is
  the nonlinear iteration.  Defaults to 1.0.

  <li>"Memory" - The maximum number of past updates that can be saved
      in memory.  Defaults to the value of "Restart Frequency".

  <li>"Linear Solver" - optional SUBLIST of linear solver parameters.

  <li>"Linear Solver"/"Tolerance" - Desired tolerance for linear
      solve. Defaults to 1.0e-4. The tolerance can be computed using
      adaptive forcing terms. See NOX::Direction::Utils::InexactNewton 
      for additional options.

  </ul>

 */  

class Broyden : public Generic {

public:
  
  //----------------------------------------------------------------------
  // MemoryUnit helper class
  //----------------------------------------------------------------------
  
  /*! 
    \brief Utility class for NOX::Direction::Broyden::BroydenMemory
    
    Stores an \f$s\f$-vector and associated information for the limited-memory Broyden update.

    Throughout the docuementation of this class, we make the following associations.
    - #sptr stores the vector \f$s\f$
    - #lambda stores the value of \f$\lambda\f$
    - #snormsqr stores the values of \f$\|s\|_2^2\f$
  */
  class BroydenMemoryUnit {
    
  public:
    
    /*! 
      \brief Constructor
      
      Sets #sptr to NULL, and #lambda and #snormsqr to zero.
    */
    BroydenMemoryUnit();
    
    /*! 
      \brief Destuctor
      
      Deletes #sptr.
    */
    ~BroydenMemoryUnit();
    

    /*! 
      \brief Reset this memory unit with a new update vector.

      Let the vector \f$d\f$ represent the input argument.
      Then we set \f$s=d\f$.
      Also recalculates \f$\|s\|_2^2\f$ and sets \f$\lambda = 0\f$.

      \note
      If #sptr is NULL, space is allocated by cloning the input vector 
      (see NOX::Abstract::Vector::clone).
    */
    void reset(const NOX::Abstract::Vector& d);
    
     /*! 
       \brief Update the step length

       Set \f$ \lambda \f$ to the input argument. 
       If \f$ \lambda \neq 1 \f$, then
       reset \f$ s = \lambda s \f$ and \f$ \|s\|_2^2 = \lambda^2 \|s\|_2^2 \f$.
     */
    void setStep(double step);

    //! Get pointer to \f$s\f$
    const NOX::Abstract::Vector* sPtr() const;
    
    //! Get the step, \f$\lambda\f$
    double step() const;
    
    //! Get \f$\|s\|_2^2 \f$
    double sNormSqr() const;
    
  private:
    
    //! A pointer to the \f$s\f$-vector
    NOX::Abstract::Vector* sptr;
    
    /*! 
        \brief The associated step length, i.e., \f$\lambda\f$

        The step is zero if the direction has not yet been scaled via
        setStep().
    */
    double lambda;
    
    //! The norm of \f$s\f$ squared, i.e., \f$s^Ts\f$
    double snormsqr;
    
  };
  
  //----------------------------------------------------------------------
  // Memory helper class
  //----------------------------------------------------------------------
  
  /*! 
    \brief 
    Utility class for NOX::Direction::Broyden method
    to manage the information stored in "limited" memory.
    
    Store up to \f$m\f$ MemoryUnit objects where \f$m\f$ is passed to
    reset().  Every time push() is called, a new
    MemoryUnit is added. If there are already \f$m\f$ MemoryUnit's, the
    oldest is bumped off the list. The zero-th entry is always the oldest.

    \note In order to avoid allocating and deallocating memory, we
    reuse the BroydenMemoryUnit objects rather than destroying and
    re-constructing them. However, this detail should be transparent
    to users of this class.
  */
  class BroydenMemory {
    
  public:
    
    //! Constructor. 
    /*! 
      Does nothing.
    */
    BroydenMemory();
    
    //! Destructor. 
    /*! 
      Does nothing.
    */
    ~BroydenMemory();
    
    //! Reset the memory
    /*! 
      Sets #mMax to the value of \c m.
      Sets the size of the #memory vector to be at least mMax.
      Sets the capacity of the #index vector to be at least mMax.
      Sets the size of the #index vector to be zero.
    */
    void reset(int m);
    
    //! Reset the memory
    /*! 
      Sets the size of the #index vector to be zero.
    */
    void reset();
    
    //! Add new information to the memory
    /*!  
      We need to calculate where the new udpate should be stored in
      #memory and update the information in #index.

      Let k denote the index of where the new update should be stored.
      If there are current m items stored in memory and m < mMax, then
      we set k = m.
      Otherwise, we set k equal to the location of the oldest
      update. The oldest update is deleted to make room for the new
      update.
      In both cases, #index must be updated appropriately so that the
      first (zero) entry points to the oldest update and the last
      entry points to the newest update.
     */
    void push(const NOX::Abstract::Vector& d);
    
    //! Returns true if the memory is empty
    bool empty() const;
    
    //! Number of items in memory
    int size() const;
    
    //! Return the ith MemoryUnit (unchecked access) 
    /*! 
      The zero entry is the oldest memory.
      The m-1 entry is the newest entry (where m denotes the memory size).
    */
    BroydenMemoryUnit& operator[](int i);
    
  private:
    
    //! Maximum number of updates that can be stored
    int mMax;

    //! Index Into memory vector
    /*!
      Keeps track of the appropriate index's into the memory vector.
      <ul>
      <li><tt>index.size()</tt> is equal to the number of updates stored.
      <li><tt>index[0]</tt> is the index to the oldest update
      <li><tt>index.back() = index[index.size()-1]</tt> is the index to the newest update
      </ul>
     */
    vector<int> index;
    
    //! Memory
    /*! 
      The size of this vector is greater than or equal to #mMax. We
      never decrease the size of memory on a reset() call because we
      want to reserve the potential to reuse memory that has already
      been allocated. 
    */
    
    vector<BroydenMemoryUnit> memory;
  };
  
//----------------------------------------------------------------------
// END HELPER CLASSES
//----------------------------------------------------------------------

public:

  //! Constructor
  Broyden(const NOX::Utils& u, NOX::Parameter::List& params);

  //! Destructor
  virtual ~Broyden();

  // derived
  virtual bool reset(NOX::Parameter::List& params);

  //! Not supported for this direction - only works for line search based solver
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::Generic& solver);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::LineSearchBased& solver);

private:

  //! Print an error message and throw an error
  /*!  
    If error printing is enabled for this process (see NOX::Utils),
    prints the following where \<functionName\> and \<errorMsg\> are
    replaced by the strings that are passed in.

\verbatim
NOX::Direction::Broyden::<functionName> - <errorMsg>
\endverbatim

    Then throws the string "NOX Error".
  */
  void throwError(const string& functionName, const string& errorMsg);

  //! Determine if it's time to restart
  /*!
    A restart should be performed if any of the following conditions hold:
    <ul>
    <li> It's the very first iteration
    <li> The current Jacobian has been used #cntMax times
    <li> \f$ \frac{ \|F_{new}\|}{\|F_{old}\|} > \alpha \f$ where 
         \f$ \alpha \f$ is the maximum convergence rate that is 
         specified in the parameter list (defaults to 1.0).
    </ul>
  */
  bool doRestart(NOX::Abstract::Group& soln, 
		       const NOX::Solver::LineSearchBased& solver);

private:

  //! Printing Utils
  const NOX::Utils& utils;

  //! "Linear Solver" sublist from "Direction"/"Broyden"
  /*! 
    \note This is pointer rather than a reference to allow for the
    reset function. It should not be deleted in the destructor.
   */
  NOX::Parameter::List* lsParamsPtr;

  //! Group to hold the "old" Jacobian
  NOX::Abstract::Group* oldJacobianGrpPtr;

  /*! 
    \brief Counter of how many times the Jacobian associated with
    #OldJacobianGrpPtr has been used since its last update 
  */
  int cnt;

  /*! 
    \brief Maximum number of times the Jacboian should be used between updates. 
    \note Should only be changed by reset() 
  */ 
  int cntMax;

  //! Last computed convergence rate.
  double convRate;

  //! Greatest allowable convergence rate before the Jacobian is recomputed.
  double maxConvRate;

  //! Maximum memory size
  /*! \note Should only be changed by reset() */
  int memorySizeMax;

  //! Storage of update vectors and related attributes
  BroydenMemory memory;

  //! Inexact Newton solve utilities
  NOX::Direction::Utils::InexactNewton inexactNewtonUtils;

};
} // namespace Direction
} // namespace NOX

#endif
#endif
