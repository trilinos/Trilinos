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

NLS_PetraVector::NLS_PetraVector(const Epetra_Vector& copyFrom, bool doCopyEntries)
{
  if (doCopyEntries) {
    // deep copy
    petraVec = new Epetra_Vector(copyFrom); 
  }
  else {
    // copy map and fill with zeros
    petraVec = new Epetra_Vector(copyFrom.Map()); 
  }

  doDeletePetraVec = true;
}

NLS_PetraVector::NLS_PetraVector(Epetra_Vector& pointTo)
{
  petraVec = &pointTo;		// copy pointer only
  doDeletePetraVec = false;	// do not delete when this is deleted
}

NLS_PetraVector::~NLS_PetraVector()
{
  if (doDeletePetraVec)
    delete petraVec;
  petraVec = NULL;
}

NLS_Vector& NLS_PetraVector::operator=(const NLS_Vector& copyFrom)
{
  throw;
}

NLS_PetraVector& NLS_PetraVector::operator=(const NLS_PetraVector& copyFrom)
{
  if (petraVec == NULL) {
    // If petraVec is empty, fill it...
    petraVec = new Epetra_Vector(*(copyFrom.petraVec)); // deep copy
    doDeletePetraVec = true;
  }

  else {
    // Otherwise, copy into existing petraVec
    int errcode = petraVec->Update(1.0, *(copyFrom.petraVec), 0.0);
    if (errcode != 0) 
      cerr << "Error in NLS_Epetra_Vec::operator=!" << endl;
  }
  
  return *this;
}

NLS_Vector& NLS_PetraVector::init(double value)
{
  petraVec->PutScalar(value);
  return *this;
}

NLS_Vector& NLS_PetraVector::abs(const NLS_Vector& base)
{
  throw; 
}

NLS_Vector& NLS_PetraVector::abs(const NLS_PetraVector& base)
{
  petraVec->Abs(*(base.petraVec));
  return *this;
}

NLS_Vector& NLS_PetraVector::copy(const NLS_Vector& y, double scale = 1.0) 
{
  throw;
}

NLS_Vector& NLS_PetraVector::copy(const NLS_PetraVector& y, double scale = 1.0)
{
  petraVec->Update(scale, *(y.petraVec), 0.0);
  return *this;
}

NLS_Vector& NLS_PetraVector::update(double alpha, 
			       const NLS_Vector& y, 
			       double beta)
{
  throw;
}

NLS_Vector& NLS_PetraVector::update(double alpha, 
				    const NLS_PetraVector& y, 
				    double beta)
{
  petraVec->Update(beta, *(y.petraVec), alpha);
  return *this;
}

NLS_Vector& NLS_PetraVector::scale(double alpha)
{
  petraVec->Scale(alpha);
  return *this;
}

NLS_Vector* NLS_PetraVector::newcopy() 
{
  NLS_PetraVector *newVec;
  newVec = new NLS_PetraVector(*this->petraVec,true);
  newVec->doDeletePetraVec = true;
  return newVec;
}

double NLS_PetraVector::infnorm() const
{
  double norm;
  petraVec->NormInf(&norm);
  return norm;
}

double NLS_PetraVector::onenorm() const
{
  double norm;
  petraVec->Norm1(&norm);
  return norm;
}

double NLS_PetraVector::norm() const
{
  double norm;
  petraVec->Norm2(&norm);
  return norm;
}

double NLS_PetraVector::norm(const NLS_Vector& weights) const
{
  throw;
}

double NLS_PetraVector::norm(const NLS_PetraVector& weights) const
{
  double norm;
  petraVec->NormWeighted(*(weights.petraVec), &norm);
  return norm;
}

double NLS_PetraVector::dot(const NLS_Vector& y) const
{
  throw;
}

double NLS_PetraVector::dot(const NLS_PetraVector& y) const
{
  double dot;
  petraVec->Dot(*(y.petraVec), &dot);
  return dot;
}

void NLS_PetraVector::print() const
{
  cout << petraVec << endl;
}

