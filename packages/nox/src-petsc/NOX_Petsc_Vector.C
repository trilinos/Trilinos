// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "petscvec.h" // Petsc Vec header 

#include "NOX_Common.H"
#include "NOX_Petsc_Vector.H"

using namespace NOX;
using namespace NOX::Petsc;

Vector::Vector(const Vec& source, CopyType type)
{
  allocate(source, type);
}

Vector::Vector(const Vector& source, CopyType type)
{
  allocate(source.getPetscVector(), type);
}

Vector::Vector(const Vec& source, string Name, CopyType type) :
  name(Name)
{
  allocate(source, type);
}

Vector::~Vector()
{
  if(isAlloc)
  {
    //cout << "\n\t\tVecDestroy called ....  " << name << "  " << this << "   " 
    //     << &petscVec << endl;
    //VecDestroy(petscVec);
    isAlloc = false;
  }
}

Abstract::Vector& 
Vector::operator=(const Vec& source)
{
  VecCopy(source, petscVec);
  name = "Unnamed";
  isAlloc = true;
  return *this;
}

Abstract::Vector& 
Vector::operator=(const Abstract::Vector& source)
{
  return operator=(dynamic_cast<const Vector&>(source));
}

Abstract::Vector& 
Vector::operator=(const Vector& source)
{
  VecCopy(source.getPetscVector(), petscVec);
  isAlloc = source.isAlloc;
  name = source.name;
  return *this;
}

int 
Vector::allocate(const Vec& source, CopyType type)
{
 
  int ierr = VecDuplicate(source, &petscVec);
  isAlloc = true;

  switch (type) {

    case DeepCopy:                // default behavior

    ierr += VecCopy(source, petscVec);
    break;

  case ShapeCopy:

    break;
  }

  if(ierr)
    cout << "ERROR: value " << ierr << " returned during "
         << "NOX::Petsc::Vector allocation !!" << endl;

  return ierr;
}

Vec& 
Vector::getPetscVector()
{
  return petscVec;
}

const Vec& 
Vector::getPetscVector() const
{
  return petscVec;
}

Abstract::Vector& 
Vector::init(double value)
{
  VecSet(petscVec, value);
  return *this;
}

Abstract::Vector& 
Vector::abs(const Abstract::Vector& base)
{
  return abs(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& 
Vector::abs(const Vector& base)
{
  VecCopy(base.getPetscVector(), petscVec);
  VecAbs(petscVec);
  return *this;
}

Abstract::Vector& 
Vector::reciprocal(const Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& 
Vector::reciprocal(const Vector& base)
{
  VecCopy(base.getPetscVector(), petscVec);
  VecReciprocal(petscVec);
  return *this;
}

Abstract::Vector& 
Vector::scale(double alpha)
{
  VecScale(petscVec, alpha);
  return *this;
}

Abstract::Vector& 
Vector::scale(const Abstract::Vector& a)
{
  return scale(dynamic_cast<const Petsc::Vector&>(a));
}
  
Abstract::Vector& 
Vector::scale(const Vector& a)
{
  VecPointwiseMult(petscVec, a.getPetscVector(), petscVec);
  return *this;
} 
  
Abstract::Vector& 
Vector::update(double alpha, const Abstract::Vector& a, 
				 double gammaval)
{
  return update(alpha, dynamic_cast<const Vector&>(a), gammaval);
}

Abstract::Vector& 
Vector::update(double alpha, const Vector& a, 
				 double gammaval)
{
  VecAXPBY(petscVec, alpha, gammaval, a.getPetscVector());
  return *this;
}

Abstract::Vector& 
Vector::update(double alpha, const Abstract::Vector& a, 
				 double beta, const Abstract::Vector& b,
				 double gammaval)
{
  return update(alpha, dynamic_cast<const Vector&>(a), 
		beta, dynamic_cast<const Vector&>(b), gammaval);
}

Abstract::Vector& 
Vector::update(double alpha, const Vector& a, 
				 double beta, const Vector& b,
				 double gammaval)
{
  VecAXPBY(petscVec, alpha, gammaval, a.getPetscVector());
  VecAXPY(petscVec, beta, b.getPetscVector());
  return *this;
}


Teuchos::RefCountPtr<NOX::Abstract::Vector> 
Vector::clone(CopyType type) const
{
  Teuchos::RefCountPtr<NOX::Abstract::Vector> newVec = 
    Teuchos::rcp(new NOX::Petsc::Vector(petscVec, type));
  return newVec;
}

double 
Vector::norm(Abstract::Vector::NormType type) const
{
  double n;
  switch (type) {
  case MaxNorm:
    VecNorm(petscVec, NORM_INFINITY, &n);
    break;
  case OneNorm:
    VecNorm(petscVec, NORM_1, &n);
    break;
  case TwoNorm:
  default:
   VecNorm(petscVec, NORM_2, &n);
   break;
  }
  return n;
}

double 
Vector::norm(const Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const Vector&>(weights));
}

double 
Vector::norm(const Vector& weights) const
{
  double n = 0.0;
  cerr << "Norm type not supported for weighted norm" << endl;
  throw;
  return n;
}

double 
Vector::innerProduct(const Abstract::Vector& y) const
{
  return innerProduct(dynamic_cast<const Vector&>(y));
}

double 
Vector::innerProduct(const Vector& y) const
{
  double dotprod;
  VecDot(y.getPetscVector(), petscVec, &dotprod);
  return dotprod;
}

int 
Vector::length() const
{
  int size;
  VecGetSize(petscVec, &size);
  return size;
}
