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
class Vector : public virtual NOX::Abstract::Vector {

 public:			

  //! Construct by copying map and/or elements of an Epetra_Vector.
  Vector(const Epetra_Vector& source, NOX::CopyType type = NOX::DeepCopy);

  //! Copy constructor.
  Vector(const NOX::Epetra::Vector& source, NOX::CopyType type = NOX::DeepCopy);

  //! Destruct Vector.
  ~Vector();

  //@{ \name Access to underlying Petra vector.

  //! Get reference to underlying Epetra vector.
  virtual Epetra_Vector& getEpetraVector();

  //! Get const reference to underlying Epetra vector.
  virtual const Epetra_Vector& getEpetraVector() const;
 
  //@}

  //@{ \name Initialization methods.

  // derived
  virtual NOX::Abstract::Vector& init(double value);

  // derived
  virtual NOX::Abstract::Vector& random(bool useSeed = false, double seed = 0.0);

  //! Copies source vector into "this".
  virtual NOX::Abstract::Vector& operator=(const Epetra_Vector& y);

  // derived
  virtual NOX::Abstract::Vector& operator=(const NOX::Epetra::Vector& y);
  virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);
  
  // derived
  virtual NOX::Abstract::Vector& abs(const NOX::Epetra::Vector& y);
  virtual NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);

  // derived
  virtual NOX::Abstract::Vector& reciprocal(const NOX::Epetra::Vector& y);
  virtual NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);

  //@}

  //@{ \name Update methods.

  // derived
  virtual NOX::Abstract::Vector& scale(double gamma);

  // derived
  virtual NOX::Abstract::Vector& scale(const NOX::Epetra::Vector& a);
  virtual NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);

  // derived
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Epetra::Vector& a, 
			     double gamma = 0.0);
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, 
			     double gamma = 0.0);

  // derived
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Epetra::Vector& a, 
			     double beta, const NOX::Epetra::Vector& b,
			     double gamma = 0.0);
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, 
			     double beta, const NOX::Abstract::Vector& b,
			     double gamma = 0.0);


  //@}

  //@{ \name Creating new Vectors. 

  // derived
  virtual NOX::Abstract::Vector* clone(CopyType type = DeepCopy) const;

  //@}

  //@{ \name Norms.

  // derived
  virtual double norm(NOX::Abstract::Vector::NormType type = TwoNorm) const;

  // derived
  virtual double norm(const NOX::Epetra::Vector& weights) const;
  virtual double norm(const NOX::Abstract::Vector& weights) const;

  //@}

  //@{ \name Dot products

  // derived
  virtual double dot(const NOX::Epetra::Vector& y) const;
  virtual double dot(const NOX::Abstract::Vector& y) const;

  //@}

  // derived
  virtual int length() const;

  // derived
  virtual void print() const;

 protected:
  
  //! Pointer to petra vector owned by this object
  Epetra_Vector* epetraVec;

};
} // namespace Epetra
} // namespace NOX

#endif
