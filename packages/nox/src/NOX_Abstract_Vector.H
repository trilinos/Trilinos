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

#ifndef NOX_ABSTRACT_VECTOR_H
#define NOX_ABSTRACT_VECTOR_H

//! Nonlinear solvers package namespace
namespace NOX {

  //! Specify whether to copy using deep copy or just copy by shape.
  enum CopyType { 
    //! Copy object including all data
    DeepCopy, 
    //! Copy the shape of the object only
    ShapeCopy 
  };

//! %NOX abstract interface for vector and group
namespace Abstract {

/*! 
  \brief %NOX's pure abstract vector interface for vectors that are
  used by the nonlinear solver. We require only limited functionality for
  the vectors.

  The documentation throughout assumes that \f$x\f$ denotes "this"
  vector, \f$x_i\f$ is the i-th component (or entry) of \f$x\f$, and
  \f$n\f$ is the length of \f$x\f$.

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

class Vector {

public:
  
  //! Norm types used in norm() calculations
  enum NormType {
    //! Use the 2-norm
    TwoNorm, 
    //! Use the 1-norm
    OneNorm, 
    //! Use the max-norm, a.k.a. the \f$\infty\f$-norm
    MaxNorm
  };

  //! %Abstract %Vector constructor (does nothing)
  Vector() {};

  //! %Abstract %Vector destructor (does nothing)
  virtual ~Vector() {};

  //@{ \name Initialization methods.

  //! Initialize every element of this vector with \c gamma.
  /*! 
    \f[ x_i = \gamma \quad \mbox{for } i=1,\dots,n \f] 
    \return Reference to x
  */
  virtual NOX::Abstract::Vector& init(double gamma) = 0;

  //! Initialize each element of this vector with a random value
  /*
    If \c useSeed is true, uses the value of \c seed to seed the
    random number generator before filling the entries of x. So, if
    two calls are made where \c useSeed is true and \c seed is the
    same, then the vectors returned should be the same.

    Default implementation throw an error. Only used by a subset of methods.

    \return Reference to x
   */
  virtual NOX::Abstract::Vector& random(bool useSeed = false, double seed = 0.0);

  //! Put element-wise absolute values of source vector \c y into this vector.
  /*! 
    \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to x
  */
  virtual NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y) = 0;

  //! Copy source vector \c y into this vector.
  /*! 
    \f[ x_i = y_i \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to x
  */
  virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y) = 0;

  //! Put element-wise reciprocal of source vector \c y into this vector.
  /*! 
    \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f] 

    \return Reference to x
  */
  virtual NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y) = 0;

  //@}

  //@{ \name Update methods.

  //! Scale each element of this vector by \c gamma.
  /*! 
    \f[ x_i = \gamma x_i \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to x
  */
  virtual NOX::Abstract::Vector& scale(double gamma) = 0;

  //! Scale this vector <em>element-by-element</em> by the vector a.
  /*! 
    \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to x
  */
  virtual NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a) = 0;

  //! Compute x = (alpha * a) + (x * this) where x is this vector.
  /*! 
    \f[ x_i = \alpha \; a_i + \gamma \; x_i \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to x
  */
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, double gamma = 0.0) = 0;

  //! Compute x = (alpha * a) + (beta * b) + (x * this) where x is this vector.
  /*! 
    \f[ x_i = \alpha \; a_i + \beta \; b_i + \gamma \; x_i \quad \mbox{for } i=1,\dots,n \f] 

    \return Reference to x
  */
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, 
			 double beta, const NOX::Abstract::Vector& b,
			 double gamma = 0.0) = 0;

  //@}

  //@{ \name Creating new Vectors.

  /*! 
    \brief Create a new %Vector of the same underlying type by
    cloning "this", and return a pointer to the new vector.  

    If type is NOX::DeepCopy, then we need to create an exact replica
    of "this". Otherwise, if type is NOX::ShapeCopy, we need only
    replicate the shape of "this". 

    \return Pointer to newly created vector or NULL if clone is not supported. 
  */
  virtual NOX::Abstract::Vector* clone(NOX::CopyType type = NOX::DeepCopy) const = 0;

  //@}

  //@{ \name Norms.

  //! Norm.
  /*! Defined as follows for each NOX::Abstract::Vector::NormType:
    <ul>
    <li>\c TWO : 2-norm \f[ \|x\| = \sqrt{\sum_{i=1}^{n} x_i^2} \f] 
    <li>\c ONE : 1-norm \f[ \|x\| = \sum_{i=1}^{n} |x_i| \f]
    <li>\c INF : infinity- or max-norm \f[ \|x\| = \max_{i} |x_i| \f]
    </uL>

    \return \f$\|x\|\f$
  */
  virtual double norm(NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TwoNorm) const = 0;

  //! Weighted 2-Norm.
  /*! Defined as 
    \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f] 
    \return \f$ \|x\|_w \f$
  */
  virtual double norm(const NOX::Abstract::Vector& weights) const = 0;

  //@}

  //@{ \name Dot products.

  //! Dot product with \c y.
  /*! 
    Defined as
    \f[ <x,y> = \sum_{i=1}^n x_i y_i \f] 
    \return \f$<x,y>\f$
   */
  virtual double dot(const NOX::Abstract::Vector& y) const = 0;
  
  //@}

  //! Return the length of vector.
  /*! 
    \return n
    \note Even if the vector is distributed across processors, this
    should return the \e global length of the vector.
   */
  virtual int length() const = 0;

  //! Print the vector.  This is meant for debugging purposes only and does not need to be defined.
  virtual void print() const;

}; // class Vector
} // namespace Abstract
} // namespace NOX

#endif
