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

#include "NOX_Epetra_Vector.H"
#include "Epetra_Vector.h"

using namespace NOX;
using namespace NOX::Epetra;

Vector::Vector(const Epetra_Vector& source, CopyType type)
{
  switch (type) {

  case DeepCopy:		// default behavior

    epetraVec = new Epetra_Vector(source); 
    break;

  case ShapeCopy:

    epetraVec = new Epetra_Vector(source.Map()); 
    break;  

  }
}

Vector::Vector(const Vector& source, CopyType type)
{

  switch (type) {

  case DeepCopy:		// default behavior

    epetraVec = new Epetra_Vector(source.getEpetraVector()); 
    break;

  case ShapeCopy:

    epetraVec = new Epetra_Vector(source.getEpetraVector().Map()); 
    break;  

  }
}

Vector::~Vector()
{
  delete epetraVec;
}

Abstract::Vector& Vector::operator=(const Epetra_Vector& source)
{
  epetraVec->Scale(1.0, source);
  return *this;
}

Abstract::Vector& Vector::operator=(const Abstract::Vector& source)
{
  return operator=(dynamic_cast<const Vector&>(source));
}

Abstract::Vector& Vector::operator=(const Vector& source)
{
  epetraVec->Scale(1.0, source.getEpetraVector());
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
  epetraVec->PutScalar(value);
  return *this;
}

Abstract::Vector& Vector::abs(const Abstract::Vector& base)
{
  return abs(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::abs(const Vector& base)
{
  epetraVec->Abs(base.getEpetraVector());
  return *this;
}

Abstract::Vector& Vector::reciprocal(const Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::reciprocal(const Vector& base)
{
  epetraVec->Reciprocal(base.getEpetraVector());
  return *this;
}

Abstract::Vector& Vector::scale(double alpha)
{
  epetraVec->Scale(alpha);
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
  epetraVec->Update(alpha, a.getEpetraVector(), gamma);
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
  epetraVec->Update(alpha, a.getEpetraVector(), beta, b.getEpetraVector(), gamma);
  return *this;
}

Abstract::Vector& Vector::scale(const Abstract::Vector& a)
{  
  return scale(dynamic_cast<const Epetra::Vector&>(a));
}

Abstract::Vector& Vector::scale(const Vector& a)
{  
  epetraVec->Multiply(1.0, *epetraVec, a.getEpetraVector(), 0.0);
  return *this;
}

Abstract::Vector* Vector::clone(CopyType type) const
{
  Vector* newVec = new Vector(*epetraVec, type);
  return newVec;
}

double Vector::norm(Abstract::Vector::NormType type) const
{
  double n;
  switch (type) {
  case MaxNorm:
    epetraVec->NormInf(&n);
    break;
  case OneNorm:
    epetraVec->Norm1(&n);
    break;
  case TwoNorm:
  default:
   epetraVec->Norm2(&n);
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
    cerr << "NOX::Epetra::Vector - Weighted norm not supported" << endl;
    throw "NOX-Epetra Error";
}

double Vector::dot(const Abstract::Vector& y) const
{
  return dot(dynamic_cast<const Vector&>(y));
}

double Vector::dot(const Vector& y) const
{
  double dot;
  epetraVec->Dot(y.getEpetraVector(), &dot);
  return dot;
}

int Vector::length() const
{
  return epetraVec->GlobalLength();
}

void Vector::print() const
{
  epetraVec->Print(cout);
  return;
}
