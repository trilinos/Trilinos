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

#include "NOX_Epetra_Vector.H"

// Included multivector declarations if requested
#ifdef HAVE_NOX_MULTIVECS
#include "NOX_Epetra_MultiVector.H"
#endif

#include "Epetra_Vector.h"

NOX::Epetra::Vector::Vector(Epetra_Vector& source, NOX::CopyType type,
			    bool createView)
{
  if (createView) {
    ownsEpetraVector = false;
    epetraVec = &source;
  }
  else {
    ownsEpetraVector = true;

    switch (type) {
      
    case DeepCopy:		// default behavior
      
      epetraVec = new Epetra_Vector(source); 
      break;
      
    case ShapeCopy:
      
      epetraVec = new Epetra_Vector(source.Map()); 
      break;  
    }
  }
}

NOX::Epetra::Vector::Vector(const Epetra_Vector& source, NOX::CopyType type) :
  ownsEpetraVector(true)
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

NOX::Epetra::Vector::Vector(const NOX::Epetra::Vector& source, 
			    NOX::CopyType type) :
  ownsEpetraVector(true)
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

NOX::Epetra::Vector::~Vector()
{
  if (ownsEpetraVector)
    delete epetraVec;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::operator=(const Epetra_Vector& source)
{
  epetraVec->Scale(1.0, source);
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::operator=(const NOX::Abstract::Vector& source)
{
  return operator=(dynamic_cast<const NOX::Epetra::Vector&>(source));
}

NOX::Abstract::Vector& NOX::Epetra::Vector::operator=(const NOX::Epetra::Vector& source)
{
  epetraVec->Scale(1.0, source.getEpetraVector());
  return *this;
}

Epetra_Vector& NOX::Epetra::Vector::getEpetraVector()
{
  return *epetraVec;
}

const Epetra_Vector& NOX::Epetra::Vector::getEpetraVector() const
{
  return *epetraVec;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::init(double value)
{
  epetraVec->PutScalar(value);
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::random(bool useSeed, int seed)
{
  if (useSeed)
    epetraVec->SetSeed(seed);
  epetraVec->Random();
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::abs(const NOX::Abstract::Vector& base)
{
  return abs(dynamic_cast<const NOX::Epetra::Vector&>(base));
}

NOX::Abstract::Vector& NOX::Epetra::Vector::abs(const NOX::Epetra::Vector& base)
{
  epetraVec->Abs(base.getEpetraVector());
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::reciprocal(const NOX::Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const NOX::Epetra::Vector&>(base));
}

NOX::Abstract::Vector& NOX::Epetra::Vector::reciprocal(const NOX::Epetra::Vector& base)
{
  epetraVec->Reciprocal(base.getEpetraVector());
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::scale(double alpha)
{
  epetraVec->Scale(alpha);
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::update(double alpha, const NOX::Abstract::Vector& a, 
				 double gamma)
{
  return update(alpha, dynamic_cast<const NOX::Epetra::Vector&>(a), gamma);
}

NOX::Abstract::Vector& NOX::Epetra::Vector::update(double alpha, const NOX::Epetra::Vector& a, 
				 double gamma)
{
  epetraVec->Update(alpha, a.getEpetraVector(), gamma);
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::update(double alpha, const NOX::Abstract::Vector& a, 
				 double beta, const NOX::Abstract::Vector& b,
				 double gamma)
{
  return update(alpha, dynamic_cast<const NOX::Epetra::Vector&>(a), 
		beta, dynamic_cast<const NOX::Epetra::Vector&>(b), gamma);
}

NOX::Abstract::Vector& NOX::Epetra::Vector::update(double alpha, const NOX::Epetra::Vector& a, 
				 double beta, const NOX::Epetra::Vector& b,
				 double gamma)
{
  epetraVec->Update(alpha, a.getEpetraVector(), beta, b.getEpetraVector(), gamma);
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::scale(const NOX::Abstract::Vector& a)
{  
  return scale(dynamic_cast<const Epetra::Vector&>(a));
}

NOX::Abstract::Vector& NOX::Epetra::Vector::scale(const NOX::Epetra::Vector& a)
{  
  epetraVec->Multiply(1.0, *epetraVec, a.getEpetraVector(), 0.0);
  return *this;
}

NOX::Abstract::Vector* NOX::Epetra::Vector::clone(CopyType type) const
{
  NOX::Abstract::Vector* newVec = new NOX::Epetra::Vector(*epetraVec, type);
  return newVec;
}

#ifdef HAVE_NOX_MULTIVECS

NOX::Abstract::MultiVector* 
NOX::Epetra::Vector::createMultiVector(
				    const NOX::Abstract::Vector* const* vecs,
				    int numVecs, NOX::CopyType type) const
{
  if (numVecs < 0) {
    cerr << "NOX::Epetra::Vector::createMultiVector:  Error!  Multivector" 
	 << " must have postive number of columns!" << endl;
    throw "NOX Error";
  }

  double** v = new double*[numVecs+1];
  const Epetra_BlockMap& map = epetraVec->Map();
  const NOX::Epetra::Vector* noxEpetraVecPtr;
  const Epetra_Vector* vecPtr;
 
  epetraVec->ExtractView(&(v[0]));
  for (int i=0; i<numVecs; i++) {
    noxEpetraVecPtr = dynamic_cast<const NOX::Epetra::Vector*>(vecs[i]);
    vecPtr = &(noxEpetraVecPtr->getEpetraVector());
    vecPtr->ExtractView(&(v[i+1]));
  }

  Epetra_MultiVector epetra_mv(View, map, v, numVecs+1);

  NOX::Epetra::MultiVector* mv = 
    new NOX::Epetra::MultiVector(epetra_mv, type);

  delete [] v;

  return mv;
}

NOX::Abstract::MultiVector* 
NOX::Epetra::Vector::createMultiVector(int numVecs, NOX::CopyType type) const
{
  if (numVecs <= 0) {
    cerr << "NOX::Epetra::Vector::createMultiVector:  Error!  Multivector" 
	 << " must have postive number of columns!" << endl;
    throw "NOX Error";
  }

  const Epetra_BlockMap& map = epetraVec->Map();
  Epetra_MultiVector *epetra_mv;
  
  if (type == NOX::ShapeCopy)
    epetra_mv = new Epetra_MultiVector(map, numVecs, true);
  else {
    epetra_mv = new Epetra_MultiVector(map, numVecs, false);
    Epetra_Vector* v;
    for (int i=0; i<numVecs; i++) {
      v = (*epetra_mv)(i);
      *v = *epetraVec;
    }
  }

  NOX::Epetra::MultiVector* mv = 
    new NOX::Epetra::MultiVector(*epetra_mv, type, true);

  delete epetra_mv;

  return mv;
}

#endif

double NOX::Epetra::Vector::norm(NOX::Abstract::Vector::NormType type) const
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

double NOX::Epetra::Vector::norm(const NOX::Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const NOX::Epetra::Vector&>(weights));
}

double NOX::Epetra::Vector::norm(const NOX::Epetra::Vector& weights) const
{
    cerr << "NOX::Epetra::Vector - Weighted norm not supported" << endl;
    throw "NOX-Epetra Error";
}

double NOX::Epetra::Vector::dot(const NOX::Abstract::Vector& y) const
{
  return dot(dynamic_cast<const NOX::Epetra::Vector&>(y));
}

double NOX::Epetra::Vector::dot(const NOX::Epetra::Vector& y) const
{
  double dot;
  epetraVec->Dot(y.getEpetraVector(), &dot);
  return dot;
}

int NOX::Epetra::Vector::length() const
{
  return epetraVec->GlobalLength();
}

void NOX::Epetra::Vector::print() const
{
  epetraVec->Print(cout);
  return;
}
