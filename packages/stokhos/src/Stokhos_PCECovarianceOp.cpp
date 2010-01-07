// $Id: Stokhos_MatrixFreeEpetraOp.cpp,v 1.7 2009/09/14 18:35:48 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/stokhos/src/Stokhos_MatrixFreeEpetraOp.cpp,v $ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_PCECovarianceOp.hpp"

#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

Stokhos::PCECovarianceOp::
PCECovarianceOp(const Stokhos::VectorOrthogPoly<Epetra_Vector>& X_poly) 
  : label("Stokhos::PCECovarianceOp"),
    X(),
    s(X_poly.basis()->norm_squared()),
    useTranspose(false),
    tmp_map(),
    tmp()
{
  const Epetra_BlockMap& base_map = X_poly[0].Map();
  int sz = X_poly.basis()->size();
  Teuchos::RCP<Epetra_MultiVector> XX = 
    Teuchos::rcp(new Epetra_MultiVector(base_map, sz-1));
  for (int i=0; i<sz-1; i++)
    (*XX)(i)->Scale(1.0, X_poly[i+1]);
  X = XX;
  tmp_map = 
    Teuchos::rcp(new Epetra_LocalMap(X->NumVectors(), 0, X->Map().Comm()));
}

Stokhos::PCECovarianceOp::
PCECovarianceOp(const Teuchos::RCP<const EpetraExt::BlockVector>& X_bv,
		const Stokhos::OrthogPolyBasis<int,double>& basis) 
  : label("Stokhos::PCECovarianceOp"),
    X(),
    s(basis.norm_squared()),
    useTranspose(false),
    tmp_map(),
    tmp()
{
  const Epetra_BlockMap& base_map = X_bv->GetBaseMap();
  int N = base_map.NumMyElements();
  int sz = basis.size();
  X = Teuchos::rcp(new Epetra_MultiVector(View, base_map, X_bv->Values()+N, 
					  N, sz-1));
  tmp_map = 
    Teuchos::rcp(new Epetra_LocalMap(X->NumVectors(), 0, X->Map().Comm()));
}

Stokhos::PCECovarianceOp::
PCECovarianceOp(const Teuchos::RCP<const Epetra_MultiVector>& X_,
		const Stokhos::OrthogPolyBasis<int,double>& basis) 
  : label("Stokhos::PCECovarianceOp"),
    X(X_),
    s(basis.norm_squared()),
    useTranspose(false),
    tmp_map(),
    tmp()
{
  tmp_map = 
    Teuchos::rcp(new Epetra_LocalMap(X->NumVectors(), 0, X->Map().Comm()));
}

Stokhos::PCECovarianceOp::~PCECovarianceOp()
{
}

int 
Stokhos::PCECovarianceOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  return 0;
}

int 
Stokhos::PCECovarianceOp::Apply(const Epetra_MultiVector& Input, 
				Epetra_MultiVector& Result) const
{
  // Allocate temporary storage
  int m = Input.NumVectors();
  if (tmp == Teuchos::null || tmp->NumVectors() != m)
    tmp = Teuchos::rcp(new Epetra_MultiVector(*tmp_map, m));

  // Compute X^T*Input
  tmp->Multiply('T', 'N', 1.0, *X, Input, 0.0);

  // Compute S*tmp
  for (int j=0; j<m; j++)
    for (int i=0; i<X->NumVectors(); i++)
      (*tmp)[j][i] *= s[i+1];

  // Compute X*tmp
  Result.Multiply('N', 'N', 1.0, *X, *tmp, 0.0);

  return 0;
}

int 
Stokhos::PCECovarianceOp::ApplyInverse(const Epetra_MultiVector& Input, 
					  Epetra_MultiVector& Result) const
{
  throw "PCECovarianceOp::ApplyInverse not defined!";
  return -1;
}

double 
Stokhos::PCECovarianceOp::NormInf() const
{
  return 1.0;
}


const char* 
Stokhos::PCECovarianceOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::PCECovarianceOp::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::PCECovarianceOp::HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
Stokhos::PCECovarianceOp::Comm() const
{
  return X->Map().Comm();
}
const Epetra_Map& 
Stokhos::PCECovarianceOp::OperatorDomainMap() const
{
  return dynamic_cast<const Epetra_Map&>(X->Map());
}

const Epetra_Map& 
Stokhos::PCECovarianceOp::OperatorRangeMap() const
{
  return dynamic_cast<const Epetra_Map&>(X->Map());
}

const Epetra_BlockMap&
Stokhos::PCECovarianceOp::CoeffMap() const
{
  return X->Map();
}
