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

#include "LOCA_Epetra_Vector.H"    // Class definition

#include "NOX_Epetra_Vector.H"     // Composition
#include "Epetra_Vector.h"         // Epetra_vector in Trilinos

using namespace LOCA;
using namespace LOCA::Epetra;

Vector::Vector(const Epetra_Vector& source, NOX::CopyType type) :
  NOX::Epetra::Vector::Vector(source, type)
{
  
}

Vector::~Vector()
{

}

Abstract::Vector& Vector::operator=(const Vector& source)
{
  *epetraVec=source.getEpetraVector();
  return *this;
}

Abstract::Vector& Vector::operator=(const Abstract::Vector& source)
{
  return operator=(dynamic_cast<const Vector&>(source));
}

Abstract::Vector& Vector::operator=(const NOX::Abstract::Vector& source)
{
  return operator=(dynamic_cast<const Vector&>(source));
}

Abstract::Vector& Vector::operator=(const Epetra_Vector& source)
{
  *epetraVec = source;
  return *this;
}

Epetra_Vector& Vector::getEpetraVector()
{
  return *epetraVec;
}

const Epetra_Vector& Vector::getEpetraVector() const
{
  return *epetraVec;
}

Abstract::Vector& Vector::init(double value)
{
  NOX::Epetra::Vector::init(value);
  return *this;
}

Abstract::Vector& Vector::abs(const NOX::Abstract::Vector& base)
{
  return abs(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::abs(const Vector& base)
{
  NOX::Epetra::Vector::abs(base);
  return *this;
}

Abstract::Vector& Vector::reciprocal(const NOX::Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::reciprocal(const Vector& base)
{
  NOX::Epetra::Vector::reciprocal(base);
  return *this;
}

Abstract::Vector& Vector::scale(double alpha)
{
  NOX::Epetra::Vector::scale(alpha);
  return *this;
}

Abstract::Vector& Vector::update(double alpha, const NOX::Abstract::Vector& a, 
				 double gamma)
{
  return update(alpha, dynamic_cast<const Vector&>(a), gamma);
}

Abstract::Vector& Vector::update(double alpha, const Vector& a, 
				 double gamma)
{
  NOX::Epetra::Vector::update(alpha, a, gamma);
  return *this;
}

Abstract::Vector& Vector::update(double alpha, const NOX::Abstract::Vector& a, 
				 double beta, const NOX::Abstract::Vector& b,
				 double gamma)
{
  return update(alpha, dynamic_cast<const Vector&>(a), 
		beta, dynamic_cast<const Vector&>(b), gamma);
}

Abstract::Vector& Vector::update(double alpha, const Vector& a, 
				 double beta, const Vector& b,
				 double gamma)
{
  NOX::Epetra::Vector::update(alpha, a, beta, b, gamma);
  return *this;
}

Abstract::Vector& Vector::scale(const NOX::Abstract::Vector& a)
{  
  return scale(dynamic_cast<const Epetra::Vector&>(a));
}

Abstract::Vector& Vector::scale(const Vector& a)
{  
  NOX::Epetra::Vector::scale(a);
  return *this;
}

Abstract::Vector& Vector::random()
{
  epetraVec->Random();
  return *this;
}

Abstract::Vector* Vector::clone(NOX::CopyType type) const
{
  Vector* newVec = new Vector(*epetraVec, type);
  return newVec;
}

double Vector::norm(NOX::Abstract::Vector::NormType type) const
{
  return NOX::Epetra::Vector::norm(type);
}

double Vector::norm(const NOX::Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const Vector&>(weights));
}

double Vector::norm(const Vector& weights) const
{
  return NOX::Epetra::Vector::norm(weights);
}

double Vector::dot(const NOX::Abstract::Vector& y) const
{
  return dot(dynamic_cast<const Vector&>(y));
}

double Vector::dot(const Vector& y) const
{
  return NOX::Epetra::Vector::dot(y);
}

int Vector::length() const
{
  return NOX::Epetra::Vector::length();
}

bool Vector::print() const
{
  return NOX::Epetra::Vector::print();
}
