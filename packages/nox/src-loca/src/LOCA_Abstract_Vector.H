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

#ifndef LOCA_ABSTRACT_VECTOR_H
#define LOCA_ABSTRACT_VECTOR_H

//! Nonlinear solvers package namespace
namespace LOCA {

  //! Specify whether to copy using deep copy or just copy by shape.
  enum CopyType { DeepCopy, ShapeCopy };

//! %NOX abstract interface for vector and group
namespace Abstract {

/*! 
  \brief %LOCA's pure abstract vector interface for vectors that are
  used by the NLS solver. We require only limited functionality for
  the vectors.

  The documentation throughout assumes that \f$x\f$ denotes "this"
  vector, \f$x_i\f$ is the i-th component (or entry) of \f$x\f$, and
  \f$n\f$ is the length of \f$x\f$.

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

class Vector {

public:
  
  //! Norm types.
  enum NormType {TwoNorm, OneNorm, MaxNorm};

  //! %Vector constructor.
  Vector() {};

  //! %Vector destructor.
  virtual ~Vector() {};

  //@{ \name Initialization methods.

  //! Initialize every element of "this" with the specified value, gamma.
  /*! \f[ x_i = \alpha \quad \mbox{for } i=1,\dots,n \f] */
  virtual Vector& init(double gamma) = 0;

  //! Puts element-wise absolute values of source vector y into "this".
  /*! \f[ x_i = | y_i | \quad \mbox{for } i=1,\dots,n \f] */
  virtual Vector& abs(const Vector& y) = 0;

  //! Fill the vector with random scalars between zero and one.
  virtual Vector& random() = 0;

  //! Copies source vector y into "this".
  /*! \f[ x_i = y_i \quad \mbox{for } i=1,\dots,n \f] */
  virtual Vector& operator=(const Vector& y) = 0;

  //! Puts element-wise reciprocal of source vector into "this".
  /*! \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1,\dots,n  \f] */
  virtual Vector& reciprocal(const Vector& y) = 0;

  //@}

  //@{ \name Update methods.

  //! this = gamma * this.
  /*! \f[ x_i = \gamma x_i \quad \mbox{for } i=1,\dots,n \f] */
  virtual Vector& scale(double gamma) = 0;

  //! Scale "this" by another vector, element-by-element.
  /*! \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1,\dots,n \f] */
  virtual Vector& scale(const Vector& a) = 0;

  //! this = (alpha * a) + (gamma * this).
  /*! \f[ x_i = \alpha \; a_i + \gamma \; x_i \quad \mbox{for } i=1,\dots,n \f] */
  virtual Vector& update(double alpha, const Vector& a, double gamma = 0.0) = 0;

  //! this = (alpha * a) + (beta * b) + (gamma * this).
  /*! \f[ x_i = \alpha \; a_i + \beta \; b_i + \gamma \; x_i \quad \mbox{for } i=1,\dots,n \f] */
  virtual Vector& update(double alpha, const Vector& a, 
			 double beta, const Vector& b,
			 double gamma = 0.0) = 0;

  //@}

  //@{ \name Creating new Vectors.

  /*! \brief Create a new %Vector of the same underlying type by
    cloning "this", and return a pointer to the new vector.  */
  /*! If type is \c DeepCopy, then we need to create an exact replica
    of "this". Otherwise, if type is \c ShapeCopy, we need only
    replicate the shape of "this". Returns \c NULL if clone is not
    supported. */
  virtual Vector* clone(CopyType type = DeepCopy) const = 0;

  //@}

  //@{ \name Norms.

  //! Norm.
  /*! Defined as follows for each NormType:
    <ul>
    <li>\c TWO : 2-norm \f[ \sqrt{\sum_{i=1}^{n} x_i^2} \f] 
    <li>\c ONE : 1-norm \f[ \sum_{i=1}^{n} |x_i| \f]
    <li>\c INF : infinity- or max-norm \f[ \max_{i} |x_i| \f]
    </uL>
  */
  virtual double norm(NormType type = TwoNorm) const = 0;

  //! Weighted 2-Norm.
  /*! Defined as 
    \f[ \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f] 
  */
  virtual double norm(const Vector& weights) const = 0;

  //@}

  //@{ \name Dot products.

  //! Dot product with y-vector.
  /*! \f[\sum_{i=1}^n x_i y_i\f] */
  virtual double dot(const Vector& y) const = 0;
  
  //@}

  //! Return length of vector.
  virtual int length() const = 0;

}; // class Vector
} // namespace Abstract
} // namespace NOX

#endif
