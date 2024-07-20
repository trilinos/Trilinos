// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_EpetraMultiVectorOperator.hpp"
#include "Epetra_Vector.h"
#include "Epetra_LocalMap.h"

Stokhos::EpetraMultiVectorOperator::
EpetraMultiVectorOperator(
  const Teuchos::RCP<const Epetra_MultiVector>& multi_vec_,
  bool is_multi_vec_transposed_) 
  : label("Epetra MultiVector Operator"),
    multi_vec(multi_vec_),
    nonconst_multi_vec(),
    is_multi_vec_transposed(is_multi_vec_transposed_),
    useTranspose(is_multi_vec_transposed),
    domain_map()
{
  domain_map = Teuchos::rcp(new Epetra_LocalMap(multi_vec->NumVectors(), 0,
						multi_vec->Map().Comm()));
}

Stokhos::EpetraMultiVectorOperator::
EpetraMultiVectorOperator(
  const Teuchos::RCP<Epetra_MultiVector>& multi_vec_,
  bool is_multi_vec_transposed_) 
  : label("Epetra MultiVector Operator"),
    multi_vec(multi_vec_),
    nonconst_multi_vec(multi_vec_),
    is_multi_vec_transposed(is_multi_vec_transposed_),
    useTranspose(is_multi_vec_transposed),
    domain_map()
{
  domain_map = Teuchos::rcp(new Epetra_LocalMap(multi_vec->NumVectors(), 0,
						multi_vec->Map().Comm()));
}

Stokhos::EpetraMultiVectorOperator::
~EpetraMultiVectorOperator()
{
}

int 
Stokhos::EpetraMultiVectorOperator::
SetUseTranspose(bool UseTranspose) 
{
  if (is_multi_vec_transposed)
    useTranspose = !UseTranspose;
  else
    useTranspose = UseTranspose;

  return 0;
}

int 
Stokhos::EpetraMultiVectorOperator::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  char trans = 'N';
  if (useTranspose)
    trans = 'T';

  int ret = Result.Multiply(trans, 'N', 1.0, *multi_vec, Input, 0.0);
  TEUCHOS_TEST_FOR_EXCEPTION(ret != 0, std::logic_error,
		     "Error!  Stokhos::EpetraMultiVectorOperator:  " <<
		     "Result.Multiply() returned " << ret << "!");
  
  return ret;
}

int 
Stokhos::EpetraMultiVectorOperator::
ApplyInverse(const Epetra_MultiVector& Input, 
	     Epetra_MultiVector& Result) const
{
  throw "EpetraMultiVectorOperator::ApplyInverse not defined!";
  return -1;
}

double 
Stokhos::EpetraMultiVectorOperator::
NormInf() const
{
  // ||A||_inf = || (|A_1| + ... + |A_n|) ||_inf where A_i is the i-th column
  // of the multi-vector A
  Epetra_Vector tmp1(multi_vec->Map());
  Epetra_Vector tmp2(multi_vec->Map());
  for (int j=0; j<multi_vec->NumVectors(); j++) {
    tmp1.Abs(*((*multi_vec)(j)));
    tmp2.Update(1.0, tmp1, 1.0);
  }
  double nrm;
  tmp2.NormInf(&nrm);

  return nrm;
}


const char* 
Stokhos::EpetraMultiVectorOperator::
Label() const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::EpetraMultiVectorOperator::
UseTranspose() const
{
  if (is_multi_vec_transposed)
    return !useTranspose;
  return useTranspose;
}

bool 
Stokhos::EpetraMultiVectorOperator::
HasNormInf() const
{
  return true;
}

const Epetra_Comm & 
Stokhos::EpetraMultiVectorOperator::
Comm() const
{
  return domain_map->Comm();
}
const Epetra_Map& 
Stokhos::EpetraMultiVectorOperator::
OperatorDomainMap() const
{
  if (useTranspose)
    return dynamic_cast<const Epetra_Map&>(multi_vec->Map());
  return *domain_map;
}

const Epetra_Map& 
Stokhos::EpetraMultiVectorOperator::
OperatorRangeMap() const
{
  if (useTranspose)
    return *domain_map;
  return dynamic_cast<const Epetra_Map&>(multi_vec->Map());
}
