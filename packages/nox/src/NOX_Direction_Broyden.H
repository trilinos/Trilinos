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

#ifndef NOX_DIRECTION_BROYDEN_H
#define NOX_DIRECTION_BROYDEN_H

#include "NOX_Common.H"		   

#include "NOX_Direction_Generic.H" // base class
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
  
  We will calculate a limited-memory Broyden direction 

  &nbsp;&nbsp; 
  \f$
  F(x)=0.
  \f$

  We compute the direction 

  &nbsp;&nbsp; 
  \f$d_k = -B_k^{-1} F_k.\f$ 

  Here \f$B_k\f$ is a limited-memory Broyden approximation to
  the Jacobian of \f$F\f$ at \f$x_k\f$, and \f$F_k =
  F(x_k)\f$.

  The \f$B_k\f$ is only stored implicitly using information
  from the \f$m\f$ most recent changes in \f$x\f$.

  \note
  The Broyden direction can only be used with NOX::Solver::LineSearchBased.
  It cannot be used with any other solver, include NOX::Solver::TrustRegionBased.

  <b>References</b>

  <ul>
  <li> C. T. Kelley, <em>Iterative Methods for Linear and Nonlinear Equations</em>, SIAM, 1995.
  </ul>
  

  <b>Parameters</b>

  To use this direction, specify that the "Method" is "Broyden" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).

  The following options are valid in the "Broyden" sublist of the
  "Direction" sublist of the parameters that are passed to the
  nonlinear solver:

  - "Memory" - Number of saved updates, defaults to 5.

  - "Memory Restrart" - True is we should scrap the old updates every
    time memory gets full. Defaults to false.

  - "Apply Right Linear Preconditioning" - Defaults to false.

  - "Preconditioning" - parameter list to be passed to
    NOX::Abstract::Group::ApplyRightPreconditioning if
    preconditioning is used.

  - "Compute Jacobian" - True if the Broyden method should go ahead
    and compute the Jacobian matrix. This option is here only to
    facilitate certain line searches that require the Jacobian. It is
    not recommended. Defaults to false.

 */  

class Broyden : public Generic {

public:
  
  //----------------------------------------------------------------------
  // MemoryUnit helper class
  //----------------------------------------------------------------------
  
  /*! 
    \brief Utility class for NOX::Direction::Broyden method
    
    Stores an s-vector and associated information for the limited-memory Broyden update 
  */
  class BroydenMemoryUnit {
    
  public:
    
    //! Constructor
    BroydenMemoryUnit();
    
    //! Destuctor
    ~BroydenMemoryUnit();
    

    /*! 
      \brief Reset this memory unit with the current direction. 

      Also calculates \f$\|s\|_2^2\f$ and sets \f$\lambda\f$ to zero. 
    */
    void reset(const NOX::Abstract::Vector& d);
    
     /*! 
       \brief Step the step length

       Set \f$ \lambda \f$. Reset \f$ s = \lambda s \f$ and \f$ \|s\|_2^2 = \lambda^2 \|s\|_2^2 \f$.
     */
    void setStep(double step);

    //! Get pointer to s
    const NOX::Abstract::Vector* sPtr() const;
    
    //! Get the step, \f$\lambda\f$
    double step() const;
    
    //! Get \f$\|s\|_2^2 \f$
    double sNormSqr() const;
    
  private:
    
    //! A pointer to the s-vector
    NOX::Abstract::Vector* sptr;
    
    /*! 
        \brief The associated step length, i.e., \f$\lambda\f$

        The step is zero if the direction has not yet been scaled via
        setStep().
    */
    double lambda;
    
    //! The norm of s squared, i.e., \f$s^Ts\f$
    double snormsqr;
    
  };
  
  //----------------------------------------------------------------------
  // Memory helper class
  //----------------------------------------------------------------------
  
  /*! 
    \brief 
    Utility class for NOX::Direction::Broyden method
    to manage the information stored in "limited" memory.
    
    Store up to \f$m\f$ MemoryUnit objects where \f$m\f$ was passed to
    the constructor or to reset().  Every time push() is called, a new
    MemoryUnit is added. If there are already \f$m\f$ MemoryUnit's, the
    oldest is bumped off the list. The zero entry is always the oldest.
  */
  class BroydenMemory {
    
  public:
    
    //! Constructor. Maximum number of BroydenMemoryUnit's is \c m.
    BroydenMemory();
    
    //! Destructor
    ~BroydenMemory();
    
    //! Reset the memory
    void reset(NOX::Parameter::List& p);
    
    //! Add new information to the memory
    void push(const NOX::Abstract::Vector& d);
    
    //! Returns true if the memory is empty
    bool empty() const;
    
    //! Number of items in memory
    int size() const;
    
    /*! 
      \brief Return the ith MemoryUnit (unchecked access) 
      
      The zero entry is the oldest memory.
      
      The m-1 entry is the newest entry (where m denotes the memory size).
    */
    BroydenMemoryUnit& operator[](int i);
    
  private:
    
    //! True is we should scrap the old updates everytime the memory gets full.
    bool doRestart;

    /*!
      \brief Index Into memory vector

      - <tt>index.size()</tt> is equal to the number of updates stored.
      - <tt>index[0]</tt> is the index to the oldest update
      - <tt>index.back() = index[index.size()-1]</tt> is the index to the newest update

     */
    vector<int> index;
    
    //! Memory
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

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::Generic& solver);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::LineSearchBased& solver);

private:

  //! Print an error message and throw an error
  void throwError(const string& functionName, const string& errorMsg);

private:

  //! Printing Utils
  const NOX::Utils& utils;

  //! "Direction" sublist with parameters for the direction vector
  /*! 
    \note This is pointer rather than a reference to allow for the
    reset function. It should not be deleted.
   */
  NOX::Parameter::List* paramsPtr;

  //! Storage of \f$(s,y)\f$ pairs and related attributes
  BroydenMemory memory;

  //! True if we should compute the Jacobian (even though the Broyden update doesn't use it)
  bool doComputeJacobian;

  //! True if right preconditioning should be used
  bool doRightPreconditioning;

  //! "Parameters" sublist with parameters for preconditioning
  NOX::Parameter::List* precParamsPtr;
  

};
} // namespace Direction
} // namespace NOX

#endif
#endif
