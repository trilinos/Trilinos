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

#ifndef NOX_ABSTRACT_GROUP_H
#define NOX_ABSTRACT_GROUP_H

#include "NOX_Abstract_Vector.H" // for NOX::CopyType
#include "NOX_Common.H" // for string 

namespace NOX {
namespace Parameter {
class List;
}
}

namespace NOX { 
namespace Abstract { 

/*!  \brief %NOX pure abstract interface to a "group"; i.e., a
  solution vector and the corresponding F-vector, Jacobian matrix,
  gradient vector, and Newton vector.
*/
/*! 
  This class groups together the objects relating to a particular
  solution vector. Each time the solution vector is updated via the
  "computeX" objects, all the previously computed values become
  invalid.

  Recall that our goal is to solve
  \f[
  F(x) = \left[ \begin{array}{c} F_1(x) \\ F_2(x) \\ \vdots \\ F_n(x) \\ \end{array} \right] = 0.
  \f]
  We may also consider the equivalent optimization problem
  \f[
  \min f(x) \equiv \frac{1}{2} F(x)^T F(x).
  \f]

  The \b Jacobian
  is denoted by \f$J\f$, and defined by
  \f[
  J_{ij} = \frac{\partial F_i}{\partial x_j} (x).
  \f]
  The \b Newton \b direction is the solution of
  \f[
  Js = -F.
  \f]
  The \b gradient (of \f$f\f$) is defined as \f$g \equiv J^T F\f$.
  
*/

class Group {

public:
  
  //! Constructor.
  Group() {};

  //! Destructor.
  virtual ~Group() {};
  
  /*! 
    \brief Copies the values of all vectors and any other data in
    the source group to this group.  (May invalidate shared data for
    source group.)
  */
  virtual Group& operator=(const Group& source) = 0;
      
  //@{ \name "Compute" functions.

  //! Set the solution vector, x, to y.
  virtual bool setX(const Vector& y) = 0;

  //! Compute and return solution vector, x, where this.x = grp.x + step * d.
  /*! Let \f$x\f$ denote this group's solution and \f$\hat x\f$ denote
    the solution associated with group \c grp. Then we compute 
    \f[ x = \hat x + \mbox{step} \; d. \f]
  */
  virtual bool computeX(const Group& grp, const Vector& d, double step) = 0;

  //! Compute and return F(x). 
  /*! In other words, return \f$F \equiv F(x)\f$. It's generally
    useful to also compute and stores norm of F(x) at this point for
    later access by the getNormF function. */
  virtual bool computeF() = 0;

  //! Compute Jacobian.
  /*! In other words, compute \f$  J_{ij} = \frac{\partial F_i}{\partial x_j} (x)\f$.
  */
  virtual bool computeJacobian() = 0;

  //! Compute and return gradient.
  //! Returns \c false if either F or the Jacobian have not been computed.
  /*! In other words, compute \f$g = J^T F\f$. */
  virtual bool computeGradient() = 0;

  //! Compute and return Newton direction, using parameters for nonlinear solve.
  //! Returns \c false if F and/or the Jacobian have not been computed.
  /*! Solve \f$Js = -F\f$ for \f$s\f$.  

    \note The parameter "Tolerance" may be included in the list of
    parameters - this is the ideal solution tolerance. 
  */
  virtual bool computeNewton(Parameter::List& params) = 0;

  //! Compute the Preconditioning matrix, \f$ M \f$,  based on the current solution vector.
  //! Returns true if successfully computed. 
  virtual bool computePrecMatrix() = 0;

  //@}

  /** @name Jacobian operations.
   
    Operations using the Jacobian matrix. These may not be defined in
    matrix-free scenarios.
  */

  //@{
  
  /*! 
    \brief If supported, returns \c true and calculates 
    result = Jacobian * input.  Otherwise, returns \c false.  
    Returns \c false if any errors occur, such as the Jacobian not being
    computed.
  */
  /*! In other words, compute 
    \f[ \mbox{result} = J \; \mbox{input}. \f]
  */
  virtual bool applyJacobian(const Vector& input, Vector& result) const = 0;

  /*! 
    \brief If supported, returns \c true and calculates result =
    Jacobian^T * input.  Otherwise, returns \c false.  Throws an error if
    the Jacobian has not been computed. 
  */
  /*! In other words, compute 
    \f[ \mbox{result} = J^T \; \mbox{input}. \f]
  */
  virtual bool applyJacobianTranspose(const Vector& input, Vector& result) const = 0;
    
  /*!
    \brief Applies the inverse of the preconditioning matrix to the given input vector.
  */
  /*! In other words, compute:
    \f[ \mbox{result} = M^{-1} \mbox{input} \f]
  */
  virtual bool applyPrecMatrixInverse(const Vector& input, Vector& result) const = 0;

  //@}

  /** @name "Is" functions.
   
    Checks to see if various objects have been computed. Returns \c true
    if the corresponding "compute" function has been called since the
    last update to the solution vector (via instantiation or
    computeX).
  */

  //@{
  
  //! Return \c true if F is valid.
  virtual bool isF() const = 0;
  //! Return \c true if the Jacobian is valid.
  virtual bool isJacobian() const = 0;
  //! Return \c true if the gradient is valid.
  virtual bool isGradient() const = 0;
  //! Return \c true if the Newton direction is valid.
  virtual bool isNewton() const = 0;
  //! Return \c true if the Preconditioning Matrix is valid.
  virtual bool isPrecMatrix() const = 0;
  //@}

  /** @name "Get" functions.
   
    Note that these function do not check whether or not the vectors
    are valid. Must use the "Is" functions for that purpose.
  */
  //@{ 

  //! Return solution vector.  
  virtual const Vector& getX() const = 0;

  //! Return F(x)
  virtual const Vector& getF() const = 0;

  //! Return 2-norm of F(x). 
  /*! In other words, \f[ \sqrt{\sum_{i=1}^n F_i^2} \f] */
  virtual double getNormF() const = 0;

  //! Return gradient.
  virtual const Vector& getGradient() const = 0;

  //! Return Newton direction.
  virtual const Vector& getNewton() const = 0;

  //@}


  //@{ \name Creating new Groups.

  /*! 
    \brief Create a new %Group of the same derived type as this one by
    cloning this one, and return a pointer to the new group.  

    If type is \c DeepCopy, then we need to create an exact replica of
    "this". Otherwise, if type is \c ShapeCopy, we need only replicate
    the shape of "this". Returns \c NULL if clone is not supported.
  */
  virtual Group* clone(CopyType type = DeepCopy) const = 0;

  //@}

};
} // namespace Abstract
} // namespace NOX

#endif
