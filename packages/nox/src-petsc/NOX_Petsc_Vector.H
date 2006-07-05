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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_PETSC_VECTOR_H
#define NOX_PETSC_VECTOR_H

#include "petscvec.h" 
#include "NOX_Abstract_Vector.H" // base class

namespace NOX {

//! %NOX %Petsc support.
/*!
%NOX provides support for PETSc data structures through concrete 
implementations of the <b>%NOX::Abstract::Group</b> and
<b>%NOX::Abstract::Vector</b>.  A <b>SharedJacobian</b> class
exists which for now simply wraps a PETSc matrix used as the Jacobian.  Future
enhancements may include generalizing this to an operator supporting the
finite-difference and matrix-free capabilities provided in PETSc.

An <b>Interface</b> class provides callbacks to the application needed
for function (residual) evaluations and possibly matrix (Jacobian) fills.

*/

namespace Petsc {

//! Implementation of NOX::Abstract::Vector for %Petsc vectors.
class Vector : public Abstract::Vector {

 public:			

  //! Construct by copying map and/or elements of a Petsc Vec.
  Vector(const Vec& source, CopyType type = DeepCopy);

  //! Construct by copying map and/or elements of a Petsc Vec.
  Vector(const Vec& source, string name="Unnamed", CopyType type = DeepCopy);

  //! Construct by copying a NOX::Petsc::Vector
  Vector(const Vector& source, CopyType type = DeepCopy);

  //! Destruct Vector.
  ~Vector();

  //@{ \name Access to underlying Petsc vector.

  //! Get reference to underlying Petsc vector.
  virtual Vec& getPetscVector();
  //! Get const reference to underlying Petsc vector.
  virtual const Vec& getPetscVector() const;
 
  //@}

  //@{ \name Initialization methods.

  virtual Abstract::Vector& init(double value);

  //! Copies source vector into "this".
  virtual Abstract::Vector& operator=(const Vec& source);

  virtual Abstract::Vector& operator=(const Vector& source);
  //! See above.
  virtual Abstract::Vector& operator=(const Abstract::Vector& source);
  
  virtual Abstract::Vector& abs(const Vector& source);
  //! See above.
  virtual Abstract::Vector& abs(const Abstract::Vector& source);

  virtual Abstract::Vector& reciprocal(const Vector& source);
  //! See above.
  virtual Abstract::Vector& reciprocal(const Abstract::Vector& source);

  //@}

  //@{ \name Update methods.

  virtual Abstract::Vector& scale(double gammaval);

  virtual Abstract::Vector& scale(const Vector& a);

  //! See above
  virtual Abstract::Vector& scale(const Abstract::Vector& a);

  virtual Abstract::Vector& update(double alpha, const Vector& a, 
			     double gammaval = 0.0);
  //! See above.
  virtual Abstract::Vector& update(double alpha, const Abstract::Vector& a, 
			     double gammaval = 0.0);

  virtual Abstract::Vector& update(double alpha, const Vector& a, 
			     double beta, const Vector& b,
			     double gammaval = 0.0);
  //! See above.
  virtual Abstract::Vector& update(double alpha, const Abstract::Vector& a, 
			     double beta, const Abstract::Vector& b,
			     double gammaval = 0.0);

  //@}

  //@{ \name Creating new Vectors. 

  virtual Teuchos::RefCountPtr<NOX::Abstract::Vector> 
  clone(CopyType type = DeepCopy) const;

  //@}

  //@{ \name Norms.

  virtual double norm(Abstract::Vector::NormType type = TwoNorm) const;

  virtual double norm(const Vector& weights) const;
  //! See above.
  virtual double norm(const Abstract::Vector& weights) const;

  //@}

  //@{ \name Inner products

  virtual double innerProduct(const Vector& y) const;
  //! See above.
  virtual double innerProduct(const Abstract::Vector& y) const;

  //@}

  virtual int length() const;
 
  string Name() { return name;};

 protected:

  //@{ \name Protected initialization method used by constructors.
  virtual int allocate(const Vec& source, CopyType type = DeepCopy);
  
  //! Pointer to Petsc vector owned by this object
  Vec petscVec;

  //! Flag to indicate if the vector has been allocated
  bool isAlloc;
  
  string name;

};
} // namespace Petsc
} // namespace NOX

#endif
