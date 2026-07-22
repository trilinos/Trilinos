// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  int sz = X_poly.size();
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
Stokhos::PCECovarianceOp::SetUseTranspose(bool UseTheTranspose)
{
  useTranspose = UseTheTranspose;
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
