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

#include "LOCA_Epetra_Vector.H"
#include "NOX_Epetra_Vector.H"
#include "Epetra_Vector.h"

using namespace LOCA;
using namespace LOCA::Epetra;

Vector::Vector(const Epetra_Vector& source, CopyType type)
{
  switch (type) {

  case DeepCopy:		// default behavior

    noxVec = new NOX::Epetra::Vector(source, NOX::DeepCopy); 
    break;

  case ShapeCopy:

    noxVec = new NOX::Epetra::Vector(source, NOX::ShapeCopy); 
    break;  

  }
}

Vector::~Vector()
{
  delete noxVec;
}

Abstract::Vector& Vector::operator=(const Epetra_Vector& source)
{
  *noxVec = source;
  return *this;
}

Abstract::Vector& Vector::operator=(const Abstract::Vector& source)
{
  return operator=(dynamic_cast<const Vector&>(source));
}

Abstract::Vector& Vector::operator=(const Vector& source)
{
  *noxVec=source.getEpetraVector();
  return *this;
}

NOX::Epetra::Vector& Vector::getNoxVector()
{
  return *noxVec;
}

const NOX::Epetra::Vector& Vector::getNoxVector() const
{
  return *noxVec;
}

Epetra_Vector& Vector::getEpetraVector()
{
  return noxVec->getEpetraVector();
}

const Epetra_Vector& Vector::getEpetraVector() const
{
  return noxVec->getEpetraVector();
}

Abstract::Vector& Vector::init(double value)
{
  noxVec->init(value);
  return *this;
}

Abstract::Vector& Vector::abs(const Abstract::Vector& base)
{
  return abs(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::abs(const Vector& base)
{
  noxVec->abs(base.getNoxVector());
  return *this;
}

Abstract::Vector& Vector::reciprocal(const Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::reciprocal(const Vector& base)
{
  noxVec->reciprocal(base.getNoxVector());
  return *this;
}

Abstract::Vector& Vector::scale(double alpha)
{
  noxVec->scale(alpha);
  return *this;
}

Abstract::Vector& Vector::update(double alpha, const Abstract::Vector& a, 
				 double gamma)
{
  return update(alpha, dynamic_cast<const Vector&>(a), gamma);
}

Abstract::Vector& Vector::update(double alpha, const Vector& a, 
				 double gamma)
{
  noxVec->update(alpha, a.getNoxVector(), gamma);
  return *this;
}

Abstract::Vector& Vector::update(double alpha, const Abstract::Vector& a, 
				 double beta, const Abstract::Vector& b,
				 double gamma)
{
  return update(alpha, dynamic_cast<const Vector&>(a), 
		beta, dynamic_cast<const Vector&>(b), gamma);
}

Abstract::Vector& Vector::update(double alpha, const Vector& a, 
				 double beta, const Vector& b,
				 double gamma)
{
  noxVec->update(alpha, a.getNoxVector(), beta, b.getNoxVector(), gamma);
  return *this;
}

Abstract::Vector& Vector::scale(const Abstract::Vector& a)
{  
  return scale(dynamic_cast<const Epetra::Vector&>(a));
}

Abstract::Vector& Vector::scale(const Vector& a)
{  
  noxVec->scale(a.getNoxVector());
  return *this;
}

Abstract::Vector& Vector::random()
{
  (noxVec->getEpetraVector()).Random();
  return *this;
}

Abstract::Vector* Vector::clone(CopyType type) const
{
  Vector* newVec = new Vector(noxVec->getEpetraVector(), type);
  return newVec;
}

double Vector::norm(Abstract::Vector::NormType type) const
{
  double n;
  switch (type) {
  case MaxNorm:
    noxVec->norm(NOX::Abstract::Vector::MaxNorm);
    break;
  case OneNorm:
    noxVec->norm(NOX::Abstract::Vector::OneNorm);
    break;
  case TwoNorm:
  default:
    noxVec->norm(NOX::Abstract::Vector::TwoNorm);
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
  return noxVec->norm(weights.getNoxVector());
}

double Vector::dot(const Abstract::Vector& y) const
{
  return dot(dynamic_cast<const Vector&>(y));
}

double Vector::dot(const Vector& y) const
{
  return noxVec->dot(y.getNoxVector());
}

int Vector::length() const
{
  return noxVec->length();
}
