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

#include "Epetra_config.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"

#include "LOCA_ErrorCheck.H"
#include "LOCA_Epetra_HouseholderJacOp.H"

LOCA::Epetra::HouseholderJacOp::HouseholderJacOp(
					 Epetra_Operator& jac, 
					 const Epetra_MultiVector& dfdpVec, 
					 const Epetra_MultiVector& houseVec_x,
					 double houseVec_p, double b) :
  label("LOCA::Epetra::HouseholderJacOp"),
  jacOperator(jac),
  dfdp(dfdpVec),
  house_x(houseVec_x),
  house_p(houseVec_p),
  beta(b),
  tmp(NULL)
{
}

LOCA::Epetra::HouseholderJacOp::~HouseholderJacOp()
{
  if (tmp != NULL)
    delete tmp;
}

int 
LOCA::Epetra::HouseholderJacOp::SetUseTranspose(bool UseTranspose) 
{
  if (UseTranspose == false)
    return 0;
  else {
    LOCA::ErrorCheck::throwError(
	  "LOCA::Epetra::HouseholderJacOp::SetUseTranspose",
	  "Operator does not support transpose");
    return -1;
  }
}

int 
LOCA::Epetra::HouseholderJacOp::Apply(const Epetra_MultiVector& Input, 
				      Epetra_MultiVector& Result) const
{
  double tmp_p;

  if (tmp == NULL)
    cerr << "tmp is NULL!" << endl;
  
  // Apply Householder transformation using temporary vector
  applyHouse(Input, 0.0, *tmp, tmp_p);
  
  // Compute J*tmp
  jacOperator.Apply(*tmp, Result);

  // Compute J*tmp + tmp_p*df/dp
  Result.Update(tmp_p, dfdp, 1.0);

  return 0;
}

int 
LOCA::Epetra::HouseholderJacOp::ApplyInverse(const Epetra_MultiVector& cInput, 
			       Epetra_MultiVector& Result) const
{
  LOCA::ErrorCheck::throwError(
	  "LOCA::Epetra::HouseholderJacOp::SetUseTranspose",
	  "Operator does not support ApplyInverse");
    return -1;
}

double 
LOCA::Epetra::HouseholderJacOp::NormInf() const
{
  double Jn, an;

  Jn = jacOperator.NormInf();
  dfdp.NormInf(&an);

  return Jn + an;
}


char* 
LOCA::Epetra::HouseholderJacOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
LOCA::Epetra::HouseholderJacOp::UseTranspose() const
{
  return false;
}

bool 
LOCA::Epetra::HouseholderJacOp::HasNormInf() const
{
  return jacOperator.HasNormInf();
}

const Epetra_Comm & 
LOCA::Epetra::HouseholderJacOp::Comm() const
{
  return jacOperator.Comm();
}
const Epetra_Map& 
LOCA::Epetra::HouseholderJacOp::OperatorDomainMap() const
{
  return jacOperator.OperatorDomainMap();
}

const Epetra_Map& 
LOCA::Epetra::HouseholderJacOp::OperatorRangeMap() const
{
  return jacOperator.OperatorRangeMap();
}

void
LOCA::Epetra::HouseholderJacOp::init(const Epetra_MultiVector& x)
{
  if (tmp != NULL)
    delete tmp;
  tmp = new Epetra_MultiVector(x.Map(), x.NumVectors(), true);
}

void
LOCA::Epetra::HouseholderJacOp::finish()
{
  if (tmp != NULL)
    delete tmp;
  tmp = NULL;
}

void
LOCA::Epetra::HouseholderJacOp::applyHouse(Epetra_MultiVector& x, 
					   double& p) const
{
  // Compute -beta * (u_x^T*x + u_p*p)
  double t;
  house_x.Dot(x, &t); 
  t = -beta*(t + house_p*p);

  // Compute [x; p] = [x; p] + t*[u_x; u_p]
  x.Update(t, house_x, 1.0);
  p += t*house_p;
}

void
LOCA::Epetra::HouseholderJacOp::applyHouse(const Epetra_MultiVector& x, 
					   double p,
					   Epetra_MultiVector& result_x, 
					   double& result_p) const
{
  // Compute -beta * (u_x^T*x + u_p*p)
  double t;
  house_x.Dot(x, &t); 
  t = -beta*(t + house_p*p);

  // Compute [r_x; r_p] = [x; p] + t*[u_x; u_p]
  result_x.Update(1.0, x, t, house_x, 0.0);
  result_p = p + t*house_p;
}
