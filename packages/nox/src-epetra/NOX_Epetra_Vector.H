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

#ifndef NOX_EPETRA_VECTOR_H
#define NOX_EPETRA_VECTOR_H

#include "NOX_Abstract_Vector.H" // base class

// Forward declaration
class Epetra_Vector;

namespace NOX {

//! %NOX %Epetra support.
namespace Epetra {

//! Implementation of NOX::Abstract::Vector for %Epetra vectors.
class Vector : public Abstract::Vector {

 public:			

  //! Construct by copying map and/or elements of an Epetra_Vector.
  Vector(const Epetra_Vector& source, CopyType type = DeepCopy);

  //! Destruct Vector.
  ~Vector();

  //@{ \name Access to underlying Petra vector.

  //! Get reference to underlying Epetra vector.
  virtual Epetra_Vector& getEpetraVector();
  //! Get const reference to underlying Epetra vector.
  virtual const Epetra_Vector& getEpetraVector() const;
 
  //@}

  //@{ \name Initialization methods.

  virtual Abstract::Vector& init(double value);

  //! Copies source vector into "this".
  virtual Abstract::Vector& operator=(const Epetra_Vector& y);

  virtual Abstract::Vector& operator=(const Vector& y);
  //! See above.
  virtual Abstract::Vector& operator=(const Abstract::Vector& y);
  
  virtual Abstract::Vector& abs(const Vector& y);
  //! See above.
  virtual Abstract::Vector& abs(const Abstract::Vector& y);

  virtual Abstract::Vector& reciprocal(const Vector& y);
  //! See above.
  virtual Abstract::Vector& reciprocal(const Abstract::Vector& y);

  //@}

  //@{ \name Update methods.

  virtual Abstract::Vector& scale(double gamma);

  virtual Abstract::Vector& scale(const Vector& a);

  //! See above
  virtual Abstract::Vector& scale(const Abstract::Vector& a);

  virtual Abstract::Vector& update(double alpha, const Vector& a, 
			     double gamma = 0.0);
  //! See above.
  virtual Abstract::Vector& update(double alpha, const Abstract::Vector& a, 
			     double gamma = 0.0);

  virtual Abstract::Vector& update(double alpha, const Vector& a, 
			     double beta, const Vector& b,
			     double gamma = 0.0);
  //! See above.
  virtual Abstract::Vector& update(double alpha, const Abstract::Vector& a, 
			     double beta, const Abstract::Vector& b,
			     double gamma = 0.0);


  //@}

  //@{ \name Creating new Vectors. 

  virtual Abstract::Vector* clone(CopyType type = DeepCopy) const;

  //@}

  //@{ \name Norms.

  virtual double norm(NormType type = TWO) const;

  virtual double norm(const Vector& weights) const;
  //! See above.
  virtual double norm(const Abstract::Vector& weights) const;

  //@}

  //@{ \name Dot products

  virtual double dot(const Vector& y) const;
  //! See above.
  virtual double dot(const Abstract::Vector& y) const;

  //@}

  virtual int length() const;

 protected:
  
  //! Pointer to petra vector owned by this object
  Epetra_Vector* epetraVec;

};
} // namespace Epetra
} // namespace NOX

#endif
