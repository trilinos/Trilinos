// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include <iostream>
#include "NLS_PetraVector.H"

NLS_PetraVector::NLS_PetraVector(const Epetra_Vector& copyFrom, ConstructorType type)
{
  switch (type) {

  case DeepCopy:		// default behavior

    petraVec = new Epetra_Vector(copyFrom); 
    break;

  case CopyShape:

    petraVec = new Epetra_Vector(copyFrom.Map()); 
    break;  

  }
}

NLS_PetraVector::~NLS_PetraVector()
{
  delete petraVec;
}

NLS_Vector& NLS_PetraVector::operator=(const NLS_Vector& copyFrom)
{
  return operator=(dynamic_cast<const NLS_PetraVector&>(copyFrom));
}

NLS_PetraVector& NLS_PetraVector::operator=(const NLS_PetraVector& copyFrom)
{
  petraVec->Scale(1.0, copyFrom.getPetraVector());
  return *this;
}

Epetra_Vector& NLS_PetraVector::getPetraVector()
{
  return *petraVec;
}

const Epetra_Vector& NLS_PetraVector::getPetraVector() const
{
  return *petraVec;
}

NLS_Vector& NLS_PetraVector::init(double value)
{
  petraVec->PutScalar(value);
  return *this;
}

NLS_Vector& NLS_PetraVector::abs(const NLS_Vector& base)
{
  return abs(dynamic_cast<const NLS_PetraVector&>(base));
}

NLS_Vector& NLS_PetraVector::abs(const NLS_PetraVector& base)
{
  petraVec->Abs(base.getPetraVector());
  return *this;
}

NLS_Vector& NLS_PetraVector::scale(double alpha)
{
  petraVec->Scale(alpha);
  return *this;
}

NLS_Vector& NLS_PetraVector::update(double alpha, const NLS_Vector& a, 
				    double gamma)
{
  return update(alpha, dynamic_cast<const NLS_PetraVector&>(a), gamma);
}

NLS_Vector& NLS_PetraVector::update(double alpha, const NLS_PetraVector& a, 
				    double gamma)
{
  petraVec->Update(alpha, a.getPetraVector(), gamma);
  return *this;
}

NLS_Vector& NLS_PetraVector::update(double alpha, const NLS_Vector& a, 
				    double beta, const NLS_Vector& b,
				    double gamma)
{
  return update(alpha, dynamic_cast<const NLS_PetraVector&>(a), 
		beta, dynamic_cast<const NLS_PetraVector&>(b), gamma);
}

NLS_Vector& NLS_PetraVector::update(double alpha, const NLS_PetraVector& a, 
			     double beta, const NLS_PetraVector& b,
				    double gamma)
{
  petraVec->Update(alpha, a.getPetraVector(), beta, b.getPetraVector(), gamma);
  return *this;
}


NLS_Vector* NLS_PetraVector::newcopy() const
{
  NLS_PetraVector* newVec = new NLS_PetraVector(*petraVec);
  return newVec;
}

double NLS_PetraVector::norm(NLS_Vector::NormType type) const
{
  double n;
  switch (type) {
  case INF:
    petraVec->NormInf(&n);
    break;
  case ONE:
    petraVec->Norm1(&n);
    break;
  case TWO:
  default:
   petraVec->Norm2(&n);
   break;
  }
  return n;
}

double NLS_PetraVector::norm(const NLS_Vector& weights, NLS_Vector::NormType type) const
{
  return norm(dynamic_cast<const NLS_PetraVector&>(weights), type);
}

double  NLS_PetraVector::norm(const NLS_PetraVector& weights, NLS_Vector::NormType type) const
{
  double n;
  switch (type) {
  case INF:
  case ONE:
    cerr << "Norm type not supported for weighted norm" << endl;
    throw;
    break;
  case TWO:
  default:
    petraVec->NormWeighted(weights.getPetraVector(), &n);
    break;
  }
  return n;
}

double NLS_PetraVector::dot(const NLS_Vector& y) const
{
  return dot(dynamic_cast<const NLS_PetraVector&>(y));
}

double NLS_PetraVector::dot(const NLS_PetraVector& y) const
{
  double dot;
  petraVec->Dot(y.getPetraVector(), &dot);
  return dot;
}

void NLS_PetraVector::print() const
{
  cout << *petraVec << endl;
}

