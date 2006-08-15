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

#include "Epetra_config.h"
#include "LOCA_Epetra_LowRankUpdateOp.H"

#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Epetra::LowRankUpdateOp::LowRankUpdateOp(
        const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	const Teuchos::RefCountPtr<Epetra_Operator>& jacOperator, 
	const Teuchos::RefCountPtr<const Epetra_MultiVector>& U_multiVec, 
	const Teuchos::RefCountPtr<const Epetra_MultiVector>& V_multiVec) :
  globalData(global_data),
  label("LOCA::Epetra::LowRankUpdateOp"),
  localMap(V_multiVec->NumVectors(), 0, jacOperator->Comm()),
  J(jacOperator),
  U(U_multiVec),
  V(V_multiVec),
  useTranspose(false),
  tmpMat()
{
}

LOCA::Epetra::LowRankUpdateOp::~LowRankUpdateOp()
{
}

int 
LOCA::Epetra::LowRankUpdateOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  return J->SetUseTranspose(UseTranspose);
}

int 
LOCA::Epetra::LowRankUpdateOp::Apply(const Epetra_MultiVector& Input, 
				     Epetra_MultiVector& Result) const
{
  // Number of input vectors
  int m = Input.NumVectors();

  // Compute J*Input or J^T*input
  int res = J->Apply(Input, Result);

  // Create temporary matrix to store V^T*input or U^T*input
  if (tmpMat == Teuchos::null || tmpMat->NumVectors() != m) {
    tmpMat = Teuchos::rcp(new Epetra_MultiVector(localMap, m, false));
  }

  if (!useTranspose) {

    // Compute V^T*Input
    tmpMat->Multiply('T', 'N', 1.0, *V, Input, 0.0);

    // Compute J*Input + U*(V^T*input)
    Result.Multiply('N', 'N', 1.0, *U, *tmpMat, 1.0);

  }
  else {

    // Compute U^T*Input
    tmpMat->Multiply('T', 'N', 1.0, *U, Input, 0.0);

    // Compute J^T*Input + V*(U^T*input)
    Result.Multiply('N', 'N', 1.0, *V, *tmpMat, 1.0);

  }

  return res;
}

int 
LOCA::Epetra::LowRankUpdateOp::ApplyInverse(const Epetra_MultiVector& cInput, 
					Epetra_MultiVector& Result) const
{
  globalData->locaErrorCheck->throwError(
	  "LOCA::Epetra::LowRankUpdateOp::ApplyInverse",
	  "Operator does not support ApplyInverse");
    return -1;
}

double 
LOCA::Epetra::LowRankUpdateOp::NormInf() const
{
  double Jn;
  vector<double> un(U->NumVectors());
  vector<double> vn(V->NumVectors());

  Jn = J->NormInf();
  U->NormInf(&un[0]);
  V->NormInf(&vn[0]);

  for (unsigned int i=0; i<un.size(); i++)
    Jn += un[i];
  for (unsigned int i=0; i<vn.size(); i++)
    Jn += vn[i];

  return Jn;
}


const char* 
LOCA::Epetra::LowRankUpdateOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
LOCA::Epetra::LowRankUpdateOp::UseTranspose() const
{
  return useTranspose;
}

bool 
LOCA::Epetra::LowRankUpdateOp::HasNormInf() const
{
  return J->HasNormInf();
}

const Epetra_Comm & 
LOCA::Epetra::LowRankUpdateOp::Comm() const
{
  return J->Comm();
}
const Epetra_Map& 
LOCA::Epetra::LowRankUpdateOp::OperatorDomainMap() const
{
  return J->OperatorDomainMap();
}

const Epetra_Map& 
LOCA::Epetra::LowRankUpdateOp::OperatorRangeMap() const
{
  return J->OperatorRangeMap();
}
