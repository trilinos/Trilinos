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

#ifndef NOX_DIRECTION_QUASINEWTON_H
#define NOX_DIRECTION_QUASINEWTON_H

#include "NOX_Common.H"		   

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Parameter_List.H"	   // class data element 

namespace NOX {

namespace Direction {

//! %Quasi-Newton direction computation
/*!  
  Applied limited-memory BFGS (LM-BFGS) to the problem

  \f[
  f(x) = \frac{1}{2} \|F(x)\|^2
  \f]

  The Quasi-Newton direction is computed by updating a low-rank
  approximation to the inverse Hessian, \f$H\f$, and then computing

  \f[ 
  s = -H_k g_k
  \f]

  Here \f$g_k = J_k^T F_k\f$ represents the gradient of \f$f\f$. 
  
  Suppose at some iteration $k$, we have $m_k$ updates. Then $H_k$ is
  calculated as follows.

  [TO BE FILLED IN LATER]

  To use this direction, specify that the "Method" is "QuasiNewton" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).

  The following options are valid in the "Direction" sublist of the
  parameters that are passed to the nonlinear solver:

  - "Memory" - Number of saved updates, defaults to 5.

 */  

class QuasiNewton : public Generic {

public:

  //! Constructor
  QuasiNewton(NOX::Parameter::List& params);

  //! Destructor
  virtual ~QuasiNewton();

  // derived
  virtual bool reset(NOX::Parameter::List& params);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::Generic& solver);

private:

  //! Print an error message and throw and error
  void throwError(const string& functionName, const string& errorMsg);

protected:

  //! Private auxillary class for NOX::Direction::QuasiNewton method
  /*! Stores an (s,y)-pair for the BFGS update */
  class Update {
    
  public:
    
    //! Constructor
    /*! See reset(). */
    Update(const NOX::Abstract::Vector& newX, const NOX::Abstract::Vector& oldX, 
	   const NOX::Abstract::Vector& newG, const NOX::Abstract::Vector& oldG);
    
    //! Destuctor
    ~Update();

    //! Reset s and y with new vectors
    /* \f[
       s = x_{\rm new} - x_{\rm old} 
       \quad
       y = g_{\rm new} - g_{\rm old} 
       \f]
    */
    void reset(const NOX::Abstract::Vector& newX, const NOX::Abstract::Vector& oldX, 
	       const NOX::Abstract::Vector& newG, const NOX::Abstract::Vector& oldG);
    
    //! Get s
    const NOX::Abstract::Vector& s() const;

    //! Get y
    const NOX::Abstract::Vector& y() const;

    //! Get \f$s^T y \f$
    double sdoty() const;

  private:
    
    //! s-vector pointer
    NOX::Abstract::Vector* sptr;

    //! y-vector pointer
    NOX::Abstract::Vector* yptr;

    //! \f$\rho = s^Ty\f$
    double rho;
    
  };
  
  //! Const iterator for deque of Update pointers
  typedef deque<Update*>::const_iterator UpdateConstIterator;

  //! Const reverse iterator for deque of Update pointers
  typedef deque<Update*>::const_reverse_iterator UpdateConstReverseIterator;

  //! Private auxillary class for NOX::Direction::QuasiNewton method
  //! to manage the limited-memory update vectors.

  /* Acts generally like a deque except the following. 

  <ul> 
  <li>Elements can only be inserted at the end.
  <li>If the deque if full, the first element is automatically pop'd to make room for the insertion.
  <\ul>
  */

  class Updates {

  public:
    
    //! Constructor. Maximum number of Updates is \c m.
    Updates(int m = 0);

    //! Destructor
    ~Updates();
    
    //! Clear the deque and update the maximum number of updates
    void reset(int m = 0);

    //! Add an element at the end of the deque
    void push_back(const NOX::Abstract::Vector& newX, const NOX::Abstract::Vector& oldX, 
		   const NOX::Abstract::Vector& newG, const NOX::Abstract::Vector& oldG);

    //! Returns the last update
    const Update* back() const;

    //! Returns a const iterator that can be used to begin traversing through the updates
    UpdateConstIterator begin() const;

    //! Returns a const iterator that can be used in a comparison for ending traversal through the updates
    UpdateConstIterator end() const;

    //! Returns a const iterator that can be used to begin traversing through the updates in \b reverse order
    UpdateConstReverseIterator rbegin() const;

    //! Returns a const iterator that can be used in a comparison for ending traversal through the 
    //! updates in \b reverse order
    UpdateConstReverseIterator rend() const;

    //! Returns true if there are no updates
    bool empty() const;

  private:

    //! Maximum size of deque
    int maxUpdates;

    //! Update deque
    deque<NOX::Direction::QuasiNewton::Update*> updateDeque;

    //! Deque for storing allocated but unused updates
    deque<NOX::Direction::QuasiNewton::Update*> recycleDeque;
  
  };

  //! "Direction" sublist with parameters for the direction vector
  /*! 
    \note This is pointer rather than a reference to allow for the
    reset function. It should not be deleted.
   */
  NOX::Parameter::List* paramsPtr;

  //! Updates 
  Updates updates;

};
} // namespace Direction
} // namespace NOX

#endif
