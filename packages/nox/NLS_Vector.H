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

private:

public:

  //! NLS_Vector constructor
  NLS_Vector() {};

  //! NLS_Vector deconstructor
  virtual ~NLS_Vector() {};

  //@{ \name  Value initialization methods  

  //! Initialize all values in a vector with constant value
  virtual NLS_Vector& init(double value) = 0;

  //@}

  //@{ \name Vector initialization methods

  //! Puts element-wise absolute values of input vector in target.
  virtual NLS_Vector& abs(const NLS_Vector& y) = 0;

  //! Copies scaled version of y into "this"
  virtual NLS_Vector& copy(const NLS_Vector& y, double scale = 1.0) = 0;

  //! Operator= (should define or compiler will do it for us!)
  virtual NLS_Vector& operator=(const NLS_Vector& y) = 0;

  //@}

  //@{ \name Update methods

  //! Compute this = alpha * this + beta * y
  virtual NLS_Vector& update(double alpha, const NLS_Vector& y, double beta) = 0;

  //! this = alpha * this
  virtual NLS_Vector& scale(double alpha) = 0;

  //@}

  //@{ \name Creating new NLS_Vectors 

  //! Create a \b new NLS_Vector by copying this one. 
  //! Return a pointer to the new vector. 
  virtual NLS_Vector* newcopy() = 0;

  //@}

  //@{ \name Norms

  //! Infinity-norm. Compute max(|this|).
  virtual double infnorm() const = 0;

  //! 1-norm
  virtual double onenorm() const = 0;

  //! 2-norm
  virtual double norm() const = 0;

  //! Weighted 2-norm
  virtual double norm(const NLS_Vector& weights) const = 0;

  //@}

  //@{ \name Dot products

  //! Dot product with y
  virtual double dot(const NLS_Vector& y) const = 0;
  
  //@}

};

#endif
