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

#include "petscvec.h" // Petsc Vec header 

#include "NOX_Common.H"
#include "NOX_Petsc_Vector.H"

using namespace NOX;
using namespace NOX::Petsc;

Vector::Vector(const Vec& source, CopyType type)
{
  switch (type) {

  case DeepCopy:		// default behavior

    petscVec = new Vec;
    VecDuplicate(source, petscVec);
    VecCopy(source, *petscVec);
    break;

  case ShapeCopy:

    petscVec = new Vec;
    VecDuplicate(source, petscVec);
    break;  

  }
}

Vector::~Vector()
{
  delete petscVec;
}

Abstract::Vector& Vector::operator=(const Vec& source)
{
  VecCopy(source, *petscVec);
  return *this;
}

Abstract::Vector& Vector::operator=(const Abstract::Vector& source)
{
  return operator=(dynamic_cast<const Vector&>(source));
}

Abstract::Vector& Vector::operator=(const Vector& source)
{
  VecCopy(source.getPetscVector(), *petscVec);
  return *this;
}

Vec& Vector::getPetscVector()
{
  return *petscVec;
}

const Vec& Vector::getPetscVector() const
{
  return *petscVec;
}

Abstract::Vector& Vector::init(double value)
{
  VecSet(&value, *petscVec);
  return *this;
}

Abstract::Vector& Vector::abs(const Abstract::Vector& base)
{
  return abs(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::abs(const Vector& base)
{
  VecCopy(base.getPetscVector(), *petscVec);
  VecAbs(*petscVec);
  return *this;
}

Abstract::Vector& Vector::reciprocal(const Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::reciprocal(const Vector& base)
{
  VecCopy(base.getPetscVector(), *petscVec);
  VecReciprocal(*petscVec);
  return *this;
}

Abstract::Vector& Vector::scale(double alpha)
{
  VecScale(&alpha, *petscVec);
  return *this;
}

Abstract::Vector& Vector::scale(const Abstract::Vector& a)
{
  return scale(dynamic_cast<const Petsc::Vector&>(a));
}
  
Abstract::Vector& Vector::scale(const Vector& a)
{
  VecPointwiseMult(*petscVec, a.getPetscVector(), *petscVec);
  return *this;
} 
  
Abstract::Vector& Vector::update(double alpha, const Abstract::Vector& a, 
				 double gammaval)
{
  return update(alpha, dynamic_cast<const Vector&>(a), gammaval);
}

Abstract::Vector& Vector::update(double alpha, const Vector& a, 
				 double gammaval)
{
  VecAXPBY(&alpha, &gammaval, a.getPetscVector(), *petscVec);
  return *this;
}

Abstract::Vector& Vector::update(double alpha, const Abstract::Vector& a, 
				 double beta, const Abstract::Vector& b,
				 double gammaval)
{
  return update(alpha, dynamic_cast<const Vector&>(a), 
		beta, dynamic_cast<const Vector&>(b), gammaval);
}

Abstract::Vector& Vector::update(double alpha, const Vector& a, 
				 double beta, const Vector& b,
				 double gammaval)
{
  VecAXPBY(&alpha, &gammaval, a.getPetscVector(), *petscVec);
  VecAXPY(&beta, b.getPetscVector(), *petscVec);
  return *this;
}


Abstract::Vector* Vector::clone(CopyType type) const
{
  Vector* newVec = new Vector(*petscVec, type);
  return newVec;
}

double Vector::norm(Abstract::Vector::NormType type) const
{
  double n;
  switch (type) {
  case MaxNorm:
    VecNorm(*petscVec, NORM_INFINITY, &n);
    break;
  case OneNorm:
    VecNorm(*petscVec, NORM_1, &n);
    break;
  case TwoNorm:
  default:
   VecNorm(*petscVec, NORM_2, &n);
   break;
  }
  return n;
}

double Vector::norm(const Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const Vector&>(weights));
}

double Vector::norm(const Vector& weights) const
{
  double n = 0.0;
  cerr << "Norm type not supported for weighted norm" << endl;
  throw;
  return n;
}

double Vector::dot(const Abstract::Vector& y) const
{
  return dot(dynamic_cast<const Vector&>(y));
}

double Vector::dot(const Vector& y) const
{
  double dotprod;
  VecDot(y.getPetscVector(), *petscVec, &dotprod);
  return dotprod;
}

int Vector::length() const
{
  int size;
  VecGetSize(*petscVec, &size);
  return size;
}
