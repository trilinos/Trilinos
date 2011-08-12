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
#include "EpetraExt_BlockMultiVector.h"
#include "EpetraExt_BlockUtility.h"

Stokhos::MatrixFreeOperator::
MatrixFreeOperator(
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk_,
  const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
  const Teuchos::RCP<const Epetra_Map>& range_base_map_,
  const Teuchos::RCP<const Epetra_Map>& domain_sg_map_,
  const Teuchos::RCP<const Epetra_Map>& range_sg_map_,
  const Teuchos::RCP<Teuchos::ParameterList>& params) : 
  label("Stokhos Matrix Free Operator"),
  sg_comm(sg_comm_),
  sg_basis(sg_basis_),
  epetraCijk(epetraCijk_),
  domain_base_map(domain_base_map_),
  range_base_map(range_base_map_),
  domain_sg_map(domain_sg_map_),
  range_sg_map(range_sg_map_),
  is_stoch_parallel(epetraCijk->isStochasticParallel()),
  global_col_map(),
  global_col_map_trans(),
  stoch_col_map(epetraCijk->getStochasticColMap()),
  col_importer(),
  col_importer_trans(),
  Cijk(epetraCijk->getParallelCijk()),
  block_ops(),
  scale_op(true),
  include_mean(true),
  only_use_linear(false),
  useTranspose(false),
  expansion_size(sg_basis->size()),
  num_blocks(0),
  input_col(),
  input_col_trans(),
  input_block(),
  result_block(),
  tmp(),
  tmp_trans(),
  k_begin(Cijk->k_begin()),
  k_end(Cijk->k_end())
{
  scale_op = params->get("Scale Operator by Inverse Basis Norms", true);
  include_mean = params->get("Include Mean", true);
  only_use_linear = params->get("Only Use Linear Terms", false);

  // Compute maximum number of mat-vec's needed
  if (!include_mean && index(k_begin) == 0)
    ++k_begin;
  if (only_use_linear) {
    int dim = sg_basis->dimension();
    k_end = Cijk->find_k(dim+1);
  }
  max_num_mat_vec = 0;
  for (Cijk_type::k_iterator k=k_begin; k!=k_end; ++k) {
    int nj = Cijk->num_j(k);
    if (max_num_mat_vec < nj)
      max_num_mat_vec = nj;
  }

  // Build up column map of SG operator
  int num_col_blocks = expansion_size;
  int num_row_blocks = expansion_size;
  if (is_stoch_parallel) {

    // Build column map from base domain map.  This will communicate 
    // stochastic components to column map, but not deterministic.  It would
    // be more efficient to do both, but the Epetra_Operator interface
    // doesn't have the concept of a column map.
    global_col_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*domain_base_map,
							     *stoch_col_map,
							     *sg_comm));
    global_col_map_trans = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*range_base_map,
							     *stoch_col_map,
							     *sg_comm));

    // Build importer from Domain Map to Column Map
    col_importer = 
      Teuchos::rcp(new Epetra_Import(*global_col_map, *domain_sg_map));
    col_importer_trans = 
      Teuchos::rcp(new Epetra_Import(*global_col_map_trans, *range_sg_map));

    num_col_blocks = epetraCijk->numMyCols();
    num_row_blocks = epetraCijk->numMyRows();
  } 

  input_block.resize(num_col_blocks);
  result_block.resize(num_row_blocks);
}

Stokhos::MatrixFreeOperator::
~MatrixFreeOperator()
{
}

void 
Stokhos::MatrixFreeOperator::
setupOperator(
   const Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >& ops)
{
  block_ops = ops;
  num_blocks = block_ops->size();
  if (num_blocks < Cijk->num_k())
    k_end = Cijk->find_k(num_blocks);
}

Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly > 
Stokhos::MatrixFreeOperator::
getSGPolynomial()
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::EpetraOperatorOrthogPoly > 
Stokhos::MatrixFreeOperator::
getSGPolynomial() const
{
  return block_ops;
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
  // Note for transpose:
  // The stochastic matrix is symmetric, however the matrix blocks may not
  // be.  So the algorithm here is the same whether we are using the transpose
  // or not.  We just apply the transpose of the blocks in the case of
  // applying the global transpose, and make sure the imported Input
  // vectors use the right map.

  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (Input.Values() == Result.Values() && !is_stoch_parallel) {
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
      (tmp == Teuchos::null || tmp->NumVectors() != m*max_num_mat_vec))
    tmp = Teuchos::rcp(new Epetra_MultiVector(*result_base_map, 
					      m*max_num_mat_vec));
  else if (useTranspose == true && 
	   (tmp_trans == Teuchos::null || 
	    tmp_trans->NumVectors() != m*max_num_mat_vec))
    tmp_trans = Teuchos::rcp(new Epetra_MultiVector(*result_base_map, 
						    m*max_num_mat_vec));
  Epetra_MultiVector *tmp_result;
  if (useTranspose == false)
    tmp_result = tmp.get();
  else
    tmp_result = tmp_trans.get();

  // Map input into column map
  const Epetra_MultiVector *tmp_col;
  if (!is_stoch_parallel)
    tmp_col = input;
  else {
    if (useTranspose == false) {
      if (input_col == Teuchos::null || input_col->NumVectors() != m)
	input_col = Teuchos::rcp(new Epetra_MultiVector(*global_col_map, m));
      input_col->Import(*input, *col_importer, Insert);
      tmp_col = input_col.get();
    }
    else {
      if (input_col_trans == Teuchos::null || 
	  input_col_trans->NumVectors() != m)
	input_col_trans = 
	  Teuchos::rcp(new Epetra_MultiVector(*global_col_map_trans, m));
      input_col_trans->Import(*input, *col_importer_trans, Insert);
      tmp_col = input_col_trans.get();
    }
  }

  // Extract blocks
  EpetraExt::BlockMultiVector sg_input(View, *input_base_map, *tmp_col);
  EpetraExt::BlockMultiVector sg_result(View, *result_base_map, Result);
  for (int i=0; i<input_block.size(); i++)
    input_block[i] = sg_input.GetBlock(i);
  for (int i=0; i<result_block.size(); i++)
    result_block[i] = sg_result.GetBlock(i);
  int N = result_block[0]->MyLength();

  // Apply block SG operator via
  // w_i = 
  //    \sum_{j=0}^P \sum_{k=0}^L J_k v_j < \psi_i \psi_j \psi_k > / <\psi_i^2>
  // for i=0,...,P where P = expansion_size, L = num_blocks, w_j is the jth 
  // input block, w_i is the ith result block, and J_k is the kth block operator
  
  // k_begin and k_end are initialized in the constructor
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
    int k = index(k_it);
    Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
    Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
    int nj = Cijk->num_j(k_it);
    if (nj > 0) {
      Teuchos::Array<double*> j_ptr(nj*m);
      Teuchos::Array<int> mj_indices(nj*m);
      int l = 0;
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	int j = index(j_it);
	for (int mm=0; mm<m; mm++) {
	  j_ptr[l*m+mm] = input_block[j]->Values()+mm*N;
	  mj_indices[l*m+mm] = l*m+mm;
	}
	l++;
      }
      Epetra_MultiVector input_tmp(View, *input_base_map, &j_ptr[0], nj*m);
      Epetra_MultiVector result_tmp(View, *tmp_result, &mj_indices[0], nj*m);
      (*block_ops)[k].Apply(input_tmp, result_tmp);
      l = 0;
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	int j = index(j_it);
	for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	     i_it != Cijk->i_end(j_it); ++i_it) {
	  int i = index(i_it);
	  double c = value(i_it);
	  if (scale_op) {
	    int i_gid;
	    if (useTranspose)
	      i_gid = epetraCijk->GCID(j);
	    else
	      i_gid = epetraCijk->GRID(i);
	    c /= norms[i_gid];
	  }
	  for (int mm=0; mm<m; mm++)
	    (*result_block[i])(mm)->Update(c, *result_tmp(l*m+mm), 1.0);
	}
	l++;
      }
    }
  }

  // Destroy blocks
  for (int i=0; i<input_block.size(); i++)
    input_block[i] = Teuchos::null;
  for (int i=0; i<result_block.size(); i++)
    result_block[i] = Teuchos::null;

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
  return *sg_comm;
}
const Epetra_Map& 
Stokhos::MatrixFreeOperator::OperatorDomainMap() const
{
  if (useTranspose)
    return *range_sg_map;
  return *domain_sg_map;
}

const Epetra_Map& 
Stokhos::MatrixFreeOperator::OperatorRangeMap() const
{
  if (useTranspose)
    return *domain_sg_map;
  return *range_sg_map;
}
