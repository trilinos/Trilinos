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
    doDeletePetraVec = true;	
    break;

  case CopyShape:

    petraVec = new Epetra_Vector(copyFrom.Map()); 
    doDeletePetraVec = true;
    break;  

  }
}

NLS_PetraVector::~NLS_PetraVector()
{
  // delete this if it was created by a deep copy
  if (doDeletePetraVec)
    delete petraVec;

  // regardless, set the pointer to NULL
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
    petraVec = new Epetra_Vector(copyFrom.getPetraVector()); // deep copy
    doDeletePetraVec = true;
  }

  else {
    // Otherwise, copy into existing petraVec
    petraVec->Update(1.0, copyFrom.getPetraVector(), 0.0);
  }
  
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
  throw; 
}

NLS_Vector& NLS_PetraVector::abs(const NLS_PetraVector& base)
{
  petraVec->Abs(base.getPetraVector());
  return *this;
}

NLS_Vector& NLS_PetraVector::copy(const NLS_Vector& y, double scale) 
{
  throw;
}

NLS_Vector& NLS_PetraVector::copy(const NLS_PetraVector& y, double scale)
{
  petraVec->Update(scale, y.getPetraVector(), 0.0);
  return *this;
}

NLS_Vector& NLS_PetraVector::update(const NLS_Vector& x, double alpha, 
				    const NLS_Vector& y, double beta,
				    double gamma)
{
  cout << "ERROR: NLS_PetraVector::update - must pass specific NLS_Vector class!" << endl;
  throw;
}

NLS_Vector& NLS_PetraVector::update(const NLS_PetraVector& x, double alpha, 
				    const NLS_PetraVector& y, double beta,
				    double gamma)
{
  petraVec->Update(alpha, x.getPetraVector(), beta, y.getPetraVector(), gamma);
  return *this;
}

NLS_Vector& NLS_PetraVector::update(const NLS_Vector& x, const NLS_Vector& d, double step)
{
  petraVec->Update(1.0, 
		   dynamic_cast<const NLS_PetraVector&>(x).getPetraVector(), 
		   step, 
		   dynamic_cast<const NLS_PetraVector&>(d).getPetraVector(), 
		   0.0);
  return *this;
}

NLS_Vector& NLS_PetraVector::update(const NLS_PetraVector& x, const NLS_PetraVector& d, double step)
{
  petraVec->Update(1.0, x.getPetraVector(), step, d.getPetraVector(), 0.0);
  return *this;
}

NLS_Vector& NLS_PetraVector::scale(double alpha)
{
  petraVec->Scale(alpha);
  return *this;
}

NLS_Vector* NLS_PetraVector::newcopy() const
{
  NLS_PetraVector* newVec = new NLS_PetraVector(*petraVec);
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
  petraVec->NormWeighted(weights.getPetraVector(), &norm);
  return norm;
}

double NLS_PetraVector::dot(const NLS_Vector& y) const
{
  throw;
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

