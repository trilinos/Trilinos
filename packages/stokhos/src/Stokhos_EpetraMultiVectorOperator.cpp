// $Id$ 
// $Source$ 
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
  TEST_FOR_EXCEPTION(ret != 0, std::logic_error,
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
