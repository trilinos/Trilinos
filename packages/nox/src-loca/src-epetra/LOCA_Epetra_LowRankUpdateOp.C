// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<Epetra_Operator>& jacOperator, 
	const Teuchos::RCP<const Epetra_MultiVector>& U_multiVec, 
	const Teuchos::RCP<const Epetra_MultiVector>& V_multiVec,
	bool setup_for_solve) :
  globalData(global_data),
  label("LOCA::Epetra::LowRankUpdateOp"),
  localMap(V_multiVec->NumVectors(), 0, jacOperator->Comm()),
  J(jacOperator),
  U(U_multiVec),
  V(V_multiVec),
  useTranspose(false),
  tmpMat(),
  JinvU(),
  lu(),
  ipiv(),
  lapack()
{
  if (setup_for_solve) {
    int m = U->NumVectors();

    // Compute J^{-1}*U
    JinvU = Teuchos::rcp(new Epetra_MultiVector(U->Map(), m));
    J->ApplyInverse(*U, *JinvU);

    // Compute I + V^T * J^{-1} * U
    lu = Teuchos::rcp(new Epetra_MultiVector(localMap, m));
    lu->Multiply('T', 'N', 1.0, *V, *JinvU, 0.0);
    for (int i=0; i<m; i++)
      (*lu)[i][i] += 1.0;

    // Compute LU factorization of I + V^T * J^{-1} * U
    ipiv.resize(m);
    int info;
    lapack.GETRF(m, m, lu->Values(), m, &ipiv[0], &info);
  }
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
LOCA::Epetra::LowRankUpdateOp::ApplyInverse(const Epetra_MultiVector& Input, 
					    Epetra_MultiVector& Result) const
{
  // Number of input vectors
  int k = Input.NumVectors();

  // Size of update
  int m = U->NumVectors();

  // Compute J^{-1}*Input or J^{-T}*Input
  int res = J->ApplyInverse(Input, Result);

  // Create temporary matrix to store V^T*input or U^T*input
  if (tmpMat == Teuchos::null || tmpMat->NumVectors() != k) {
    tmpMat = Teuchos::rcp(new Epetra_MultiVector(localMap, k, false));
  }

  if (!useTranspose) {
    
    // Compute V^T*Result
    tmpMat->Multiply('T', 'N', 1.0, *V, Result, 0.0);

    // Backsolve LU factorization against tmpMat
    int info;
    lapack.GETRS('N', m, k, lu->Values(), m, &ipiv[0], tmpMat->Values(), m,
		 &info);

    // Compute Result - JinvU*tmpMat
    Result.Multiply('N', 'N', -1.0, *JinvU, *tmpMat, 1.0);
  }
  else {
    globalData->locaErrorCheck->throwError(
		      "LOCA::Epetra::LowRankUpdateOp::ApplyInverse",
		      "Operator does not support transpose");
    return -1;
  }
  
  return res;
}

double 
LOCA::Epetra::LowRankUpdateOp::NormInf() const
{
  double Jn;
  std::vector<double> un(U->NumVectors());
  std::vector<double> vn(V->NumVectors());

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
