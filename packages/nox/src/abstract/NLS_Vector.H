// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_VECTOR_H
#define _NLS_VECTOR_H
 
//! NLSPACK pure abstract vector interface for vectors that are used
//! by the NLS solver. We require only limited functionality for the
//! vectors.

class NLS_Vector {

public:

  //! Norm types
  enum NormType {TWO, ONE, INF};

  //! NLS_Vector constructor.
  NLS_Vector() {};

  //! NLS_Vector deconstructor
  virtual ~NLS_Vector() {};

  //@{ \name  Value initialization methods  

  //! Initialize all values in a vector with constant value
  virtual NLS_Vector& init(double value) = 0;

  //@}

  //@{ \name Vector initialization methods

  //! Puts element-wise absolute values of input vector in "this".
  virtual NLS_Vector& abs(const NLS_Vector& y) = 0;

  //! Copies y into "this".
  virtual NLS_Vector& operator=(const NLS_Vector& y) = 0;

  //! this = element-wise reciprocal of a
  virtual NLS_Vector& reciprocal(const NLS_Vector& a) = 0;

  //@}

  //@{ \name Update methods

  //! this = gamma * this
  virtual NLS_Vector& scale(double gamma) = 0;

  //! this = (alpha * a) + (gamma * this) 
  virtual NLS_Vector& update(double alpha, const NLS_Vector& a, double gamma = 0.0) = 0;

  //! this = (alpha * a) + (beta * b) + (gamma * this)
  virtual NLS_Vector& update(double alpha, const NLS_Vector& a, 
			     double beta, const NLS_Vector& b,
			     double gamma = 0.0) = 0;
  //@}

  //@{ \name Creating new NLS_Vectors 

  //! Create a \b new NLS_Vector by copying this one. 
  //! Return a pointer to the new vector. 
  virtual NLS_Vector* clone() const = 0;

  //@}

  //@{ \name Norms

  //! Norm
  virtual double norm(NormType type = TWO) const = 0;

  //! Weighted Norm
  virtual double norm(const NLS_Vector& weights, NormType type = TWO) const = 0;

  //@}

  //@{ \name Dot products

  //! Dot product with y
  virtual double dot(const NLS_Vector& y) const = 0;
  
  //@}

  //! Return length of vector
  virtual int length() const = 0;


};

#endif
