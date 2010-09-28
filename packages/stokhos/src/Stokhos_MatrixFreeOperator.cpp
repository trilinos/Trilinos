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

#include "Stokhos_MatrixFreeOperator.hpp"
#include "Epetra_config.h"
#include "EpetraExt_BlockMultiVector.h"

Stokhos::MatrixFreeOperator::
MatrixFreeOperator(const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
		   const Teuchos::RCP<const Epetra_Map>& range_base_map_,
		   const Teuchos::RCP<const Epetra_Map>& domain_sg_map_,
		   const Teuchos::RCP<const Epetra_Map>& range_sg_map_,
		   const Teuchos::RCP<Teuchos::ParameterList>& params) : 
  label("Stokhos Matrix Free Operator"),
  domain_base_map(domain_base_map_),
  range_base_map(range_base_map_),
  domain_sg_map(domain_sg_map_),
  range_sg_map(range_sg_map_),
  Cijk(),
  block_ops(),
  scale_op(true),
  include_mean(true),
  only_use_linear(false),
  useTranspose(false),
  expansion_size(0),
  num_blocks(0),
  input_block(),
  result_block(),
  tmp(),
  tmp_trans()
{
  if (params != Teuchos::null) {
    scale_op = params->get("Scale Operator by Inverse Basis Norms", true);
    include_mean = params->get("Include Mean", true);
    only_use_linear = params->get("Only Use Linear Terms", false);
  }
}

Stokhos::MatrixFreeOperator::
~MatrixFreeOperator()
{
}

void 
Stokhos::MatrixFreeOperator::
setupOperator(
   const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> >& ops,
   const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_)
{
  block_ops = ops;
  Cijk = Cijk_;
  sg_basis = block_ops->basis();
  expansion_size = sg_basis->size();
  num_blocks = block_ops->size();
  if (input_block.size() != expansion_size)
    input_block.resize(expansion_size);
  if (result_block.size() != expansion_size)
    result_block.resize(expansion_size);
}

Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > 
Stokhos::MatrixFreeOperator::
getSGPolynomial()
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::VectorOrthogPoly<Epetra_Operator> > 
Stokhos::MatrixFreeOperator::
getSGPolynomial() const
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > 
Stokhos::MatrixFreeOperator::
getTripleProduct() const
{
  return Cijk;
}

int 
Stokhos::MatrixFreeOperator::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  for (int i=0; i<num_blocks; i++)
    (*block_ops)[i].SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::MatrixFreeOperator::
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
      (tmp == Teuchos::null || tmp->NumVectors() != m*expansion_size))
    tmp = Teuchos::rcp(new Epetra_MultiVector(*result_base_map, 
					      m*expansion_size));
  else if (useTranspose == true && 
	   (tmp_trans == Teuchos::null || 
	    tmp_trans->NumVectors() != m*expansion_size))
    tmp_trans = Teuchos::rcp(new Epetra_MultiVector(*result_base_map, 
						    m*expansion_size));
  Epetra_MultiVector *tmp_result;
  if (useTranspose == false)
    tmp_result = tmp.get();
  else
    tmp_result = tmp_trans.get();

  // Extract blocks
  EpetraExt::BlockMultiVector sg_input(View, *input_base_map, *input);
  EpetraExt::BlockMultiVector sg_result(View, *result_base_map, Result);
  for (int i=0; i<expansion_size; i++) {
    input_block[i] = sg_input.GetBlock(i);
    result_block[i] = sg_result.GetBlock(i);
  }
  int N = input_block[0]->MyLength();

  // Apply block SG operator via
  // w_i = 
  //    \sum_{j=0}^P \sum_{k=0}^L J_k v_j < \psi_i \psi_j \psi_k > / <\psi_i^2>
  // for i=0,...,P where P = expansion_size, L = num_blocks, w_j is the jth 
  // input block, w_i is the ith result block, and J_k is the kth block operator
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  int k_begin = 0;
  if (!include_mean)
    k_begin = 1;
  int k_end = num_blocks;
  int dim = sg_basis->dimension();
  if (only_use_linear && num_blocks > dim+1)
    k_end = dim + 1;
  for (int k=k_begin; k<k_end; ++k) {
    Cijk_type::kj_iterator j_begin = Cijk->j_begin(k);
    Cijk_type::kj_iterator j_end = Cijk->j_end(k);
    int nj = Cijk->num_j(k);
    if (nj > 0) {
      Teuchos::Array<double*> j_ptr(nj*m);
      Teuchos::Array<int> mj_indices(nj*m);
      int l = 0;
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	int j = index(j_it);
	for (int mm=0; mm<m; mm++) {
	  j_ptr[l*m+mm] = input_block[j]->Values()+mm*N;
	  mj_indices[l*m+mm] = j*m+mm;
	}
	l++;
      }
      Epetra_MultiVector input_tmp(View, *input_base_map, &j_ptr[0], nj*m);
      Epetra_MultiVector result_tmp(View, *tmp_result, &mj_indices[0], nj*m);
      (*block_ops)[k].Apply(input_tmp, result_tmp);
      l = 0;
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	     i_it != Cijk->i_end(j_it); ++i_it) {
	  int i = index(i_it);
	  double c = value(i_it);
	  if (scale_op)
	    c /= norms[i];
	  for (int mm=0; mm<m; mm++)
	    (*result_block[i])(mm)->Update(c, *result_tmp(l*m+mm), 1.0);
	}
	l++;
      }
    }
  }

  // Destroy blocks
  for (int i=0; i<expansion_size; i++) {
    input_block[i] = Teuchos::null;
    result_block[i] = Teuchos::null;
  }

  if (made_copy)
    delete input;

  return 0;
}

int 
Stokhos::MatrixFreeOperator::ApplyInverse(const Epetra_MultiVector& Input, 
					  Epetra_MultiVector& Result) const
{
  throw "MatrixFreeOperator::ApplyInverse not defined!";
  return -1;
}

double 
Stokhos::MatrixFreeOperator::NormInf() const
{
  return 1.0;
}


const char* 
Stokhos::MatrixFreeOperator::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::MatrixFreeOperator::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::MatrixFreeOperator::HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
Stokhos::MatrixFreeOperator::Comm() const
{
  return domain_base_map->Comm();
}
const Epetra_Map& 
Stokhos::MatrixFreeOperator::OperatorDomainMap() const
{
  return *domain_sg_map;
}

const Epetra_Map& 
Stokhos::MatrixFreeOperator::OperatorRangeMap() const
{
  return *range_sg_map;
}
