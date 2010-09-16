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
// Questions? Contact Christopher W. Miller (cmiller@math.umd.edu).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_KLMatrixFreeOperator.hpp"
#include "Epetra_config.h"
#include "EpetraExt_BlockMultiVector.h"

Stokhos::KLMatrixFreeOperator::
KLMatrixFreeOperator(const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
		     const Teuchos::RCP<const Epetra_Map>& range_base_map_,
		     const Teuchos::RCP<const Epetra_Map>& domain_sg_map_,
		     const Teuchos::RCP<const Epetra_Map>& range_sg_map_,
		     const Teuchos::RCP<Teuchos::ParameterList>& params) : 
  label("Stokhos KL MatrixFree Operator"),
  domain_base_map(domain_base_map_),
  range_base_map(range_base_map_),
  domain_sg_map(domain_sg_map_),
  range_sg_map(range_sg_map_),
  Cijk(),
  block_ops(),
  scale_op(true),
  include_mean(true),
  useTranspose(false),
  expansion_size(0),
  num_blocks(0),
  result_block(),
  block_products(),
  block_products_trans()
{
  if (params != Teuchos::null) {
    scale_op = params->get("Scale Operator by Inverse Basis Norms", true);
    include_mean = params->get("Include Mean", true);
  }
}

Stokhos::KLMatrixFreeOperator::~KLMatrixFreeOperator()
{
}

void 
Stokhos::KLMatrixFreeOperator::
setupOperator(
   const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> >& ops,
   const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_)
{
  block_ops = ops;
  Cijk = Cijk_;
  sg_basis = block_ops->basis();
  expansion_size = sg_basis->size();
  num_blocks = sg_basis->dimension() + 1;
  if (result_block.size() != expansion_size)
    result_block.resize(expansion_size);
}

Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > 
Stokhos::KLMatrixFreeOperator::
getSGPolynomial()
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::VectorOrthogPoly<Epetra_Operator> > 
Stokhos::KLMatrixFreeOperator::
getSGPolynomial() const
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > 
Stokhos::KLMatrixFreeOperator::
getTripleProduct() const
{
  return Cijk;
}

int 
Stokhos::KLMatrixFreeOperator::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  for (int i=0; i<num_blocks; i++)
    (*block_ops)[i].SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::KLMatrixFreeOperator::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (Input.Values() == Result.Values()) {
    input = new Epetra_MultiVector(Input);
    made_copy = true;
  }

  // Initialize
  Result.PutScalar(0.0);

  const Epetra_Map* input_base_map = domain_base_map.get();
  const Epetra_Map* result_base_map = range_base_map.get();
  if (useTranspose == true) {
    input_base_map = range_base_map.get();
    result_base_map = domain_base_map.get();
  }

  // Allocate temporary storage
  int m = Input.NumVectors();
  if (useTranspose == false && 
      (block_products == Teuchos::null || 
       block_products->NumVectors() != m*expansion_size))
    block_products = 
      Teuchos::rcp(new Epetra_MultiVector(*result_base_map, 
					  m*expansion_size));
  else if (useTranspose == true && 
	   (block_products_trans == Teuchos::null || 
	    block_products_trans->NumVectors() != m*expansion_size))
    block_products_trans = 
      Teuchos::rcp(new Epetra_MultiVector(*result_base_map, 
					  m*expansion_size));
  Epetra_MultiVector *bp;
  if (useTranspose == false)
    bp = block_products.get();
  else
    bp = block_products_trans.get();

  // Extract blocks
  int N = input_base_map->NumMyElements();
  Epetra_MultiVector input_tmp(View, *input_base_map, input->Values(), N,
			       expansion_size*m);
  EpetraExt::BlockMultiVector sg_result(View, *result_base_map, Result);
  for (int i=0; i<expansion_size; i++)
    result_block[i] = sg_result.GetBlock(i);

  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  int d = sg_basis->dimension();
  Teuchos::Array<double> zero(d), one(d);
  for(int j = 0; j<d; j++) {
    zero[j] = 0.0;
    one[j] = 1.0;
  }
  Teuchos::Array< double > phi_0(expansion_size), phi_1(expansion_size);
  sg_basis->evaluateBases(zero, phi_0);
  sg_basis->evaluateBases(one, phi_1);

  int k_begin = 0;
  if (!include_mean)
    k_begin = 1;
  for (int k=k_begin; k<=d; k++) {
    (*block_ops)[k].Apply(input_tmp, *bp); 
    int nj = Cijk->num_j(k);
    const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
    for (int jj=0; jj<nj; jj++) {
      int j = j_indices[jj];
      const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
      const Teuchos::Array<double>& c_values = Cijk->values(k,jj);
      for (int ii=0; ii<i_indices.size(); ii++) {
  	int i = i_indices[ii];
	double c = c_values[ii];
	if (k == 0)
	  c /= phi_0[0];
	else {
	  c /= phi_1[k];
	  if (i == j)
	    c -= phi_0[k]/(phi_1[k]*phi_0[0])*norms[i];
	}
	if (scale_op)
	  c /= norms[i];
	for (int mm=0; mm<m; mm++)
	  (*result_block[i])(mm)->Update(c, *((*bp)(j*m+mm)), 1.0);
      }
    }
  }

  // Destroy blocks
  for (int i=0; i<num_blocks; i++)
    result_block[i] = Teuchos::null;

  if (made_copy)
    delete input;

  return 0;
}

int 
Stokhos::KLMatrixFreeOperator::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  throw "KLMatrixFreeOperator::ApplyInverse not defined!";
  return -1;
}

double 
Stokhos::KLMatrixFreeOperator::
NormInf() const
{
  return 1.0;
}


const char* 
Stokhos::KLMatrixFreeOperator::
Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::KLMatrixFreeOperator::
UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::KLMatrixFreeOperator::
HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
Stokhos::KLMatrixFreeOperator::
Comm() const
{
  return domain_base_map->Comm();
}
const Epetra_Map& 
Stokhos::KLMatrixFreeOperator::
OperatorDomainMap() const
{
  return *domain_sg_map;
}

const Epetra_Map& 
Stokhos::KLMatrixFreeOperator::
OperatorRangeMap() const
{
  return *range_sg_map;
}
