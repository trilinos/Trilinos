// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_ABSTRACT_VECTOR_H
#define _NOX_ABSTRACT_VECTOR_H
 
//! Nonlinear solvers package namespace
namespace NOX {

  //! Specify whether to copy using deep copy or just copy by shape.
  enum CopyType { DeepCopy, CopyShape };

//! NOX abstract interface for vector and group
namespace Abstract {

  /*! \brief NOX's pure abstract vector interface for vectors that are
   * used by the NLS solver. We require only limited functionality for
   * the vectors.*/

class Vector {

public:
  
  //! Norm types.
  enum NormType {TWO, ONE, INF};

  //! Vector constructor.
  Vector() {};

  //! Vector destructor.
  virtual ~Vector() {};

  //@{ \name Initialization methods.

  //! Initialize every element of "this" with specified value.
  virtual Vector& init(double value) = 0;

  //! Puts element-wise absolute values of source vector into "this".
  virtual Vector& abs(const Vector& source) = 0;

  //! Copies source vector into "this".
  virtual Vector& operator=(const Vector& source) = 0;

  //! Puts element-wise reciprocal of source vector into "this".
  virtual Vector& reciprocal(const Vector& source) = 0;

  //@}

  //@{ \name Update methods.

  //! this = gamma * this.
  virtual Vector& scale(double gamma) = 0;

  //! this = (alpha * a) + (gamma * this).
  virtual Vector& update(double alpha, const Vector& a, double gamma = 0.0) = 0;

  //! this = (alpha * a) + (beta * b) + (gamma * this).
  virtual Vector& update(double alpha, const Vector& a, 
			 double beta, const Vector& b,
			 double gamma = 0.0) = 0;
  //@}

  //@{ \name Creating new Vectors.

  /*! \brief Create a new Vector of the same underlying type by
   * cloning "this", and return a pointer to the new vector.  If type
   * is "DeepCopy", then we need to create an exact replica of
   * "this". Otherwise, if type is "CopyShape", we need only replicate
   * the shape of "this". Returns NULL if clone is not supported. */
  virtual Vector* clone(CopyType type = DeepCopy) const = 0;

  //@}

  //@{ \name Norms.

  //! Norm.
  virtual double norm(NormType type = TWO) const = 0;

  //! Weighted Norm.
  virtual double norm(const Vector& weights, NormType type = TWO) const = 0;

  //@}

  //@{ \name Dot products.

  //! Dot product with y-vector.
  virtual double dot(const Vector& y) const = 0;
  
  //@}

  //! Return length of vector.
  virtual int length() const = 0;

};
}
}

#endif
