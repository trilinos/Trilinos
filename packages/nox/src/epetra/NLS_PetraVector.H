// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_PETRAVECTOR_H
#define _NLS_PETRAVECTOR_H


#include "NLS_Vector.H"	// Base Class
#include "Epetra_Vector.h"

//! Implementation of NLS_Vector pure virtual class for Epetra_Vectors.

class NLS_PetraVector : public NLS_Vector {

 public:			

  //! Constructor type 
  enum ConstructorType { DeepCopy, CopyShape };

  //! Construct an NLS_PetraVector by copying map and/or elements of
  //! a Epetra_Vector.  Creates a deep copy of the vector.
  NLS_PetraVector(const Epetra_Vector& copyFrom, ConstructorType type = DeepCopy);

  //! Deconstruct NLS_PetraVector
  ~NLS_PetraVector();

  //@{ \name Get underlying Petra vector

  //! Returns the underlying vector that NLS wraps (EPetra_Vector).  
  virtual Epetra_Vector& getPetraVector();

  //! Returns the underlying vector that NLS wraps (EPetra_Vector).  
  virtual const Epetra_Vector& getPetraVector() const;
 
  //@}

  //@{ \name  Value initialization methods  

  virtual NLS_Vector& init(double value);

  //@}

  //@{ \name Vector initialization methods

  virtual NLS_Vector& operator=(const NLS_Vector& copyFrom);

  //! See above
  virtual NLS_PetraVector& operator=(const NLS_PetraVector& copyFrom);

  virtual NLS_Vector& abs(const NLS_Vector& y);

  //! See above 
  virtual NLS_Vector& abs(const NLS_PetraVector& y);

  virtual NLS_Vector& reciprocal(const NLS_Vector& a);

  //! See above
  virtual NLS_Vector& reciprocal(const NLS_PetraVector& a);


  //@}

  //@{ \name Update methods

  virtual NLS_Vector& scale(double gamma);

  virtual NLS_Vector& update(double alpha, const NLS_Vector& a, 
			     double gamma = 0.0);

  virtual NLS_Vector& update(double alpha, const NLS_PetraVector& a, 
			     double gamma = 0.0);

  virtual NLS_Vector& update(double alpha, const NLS_Vector& a, 
			     double beta, const NLS_Vector& b,
			     double gamma = 0.0);

  virtual NLS_Vector& update(double alpha, const NLS_PetraVector& a, 
			     double beta, const NLS_PetraVector& b,
			     double gamma = 0.0);


  //@}

  //@{ \name Creating new NLS_Vectors 

  virtual NLS_Vector* clone() const;

  //@}

  //@{ \name Norms

  virtual double norm(NormType type = TWO) const;

  virtual double norm(const NLS_Vector& weights, NormType type = TWO) const;

  //! See above
  virtual double norm(const NLS_PetraVector& weights, NormType type = TWO) const;

  //@}

  //@{ \name Dot products

  virtual double dot(const NLS_Vector& y) const;

  //! See above for definition.
  virtual double dot(const NLS_PetraVector& y) const;

  //@}

  //! Return length of vector
  virtual int length() const;


  //! Print Function
  void print() const;
  
  friend class NLS_PetraGroup;

 private:
  
  //! Pointer to petra vector owned by this object
  Epetra_Vector* petraVec;		
};



#endif
