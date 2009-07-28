// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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

#include "Epetra_config.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_MatrixFreeEpetraOp.hpp"

Stokhos::MatrixFreeEpetraOp::MatrixFreeEpetraOp(
 const Teuchos::RCP<const Epetra_Map>& base_map_,
 const Teuchos::RCP<const Epetra_Map>& sg_map_,
 const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
 const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
 const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> >& ops_) 
  : label("Stokhos::MatrixFreeEpetraOp"),
    base_map(base_map_),
    sg_map(sg_map_),
    sg_basis(sg_basis_),
    Cijk(Cijk_),
    block_ops(ops_),
    useTranspose(false),
    num_blocks(block_ops->size()),
    input_block(num_blocks),
    result_block(num_blocks),
    tmp(),
    tmp2()
{
  tmp2 = Teuchos::rcp(new Epetra_MultiVector(*base_map, num_blocks));
}

Stokhos::MatrixFreeEpetraOp::~MatrixFreeEpetraOp()
{
}

void 
Stokhos::MatrixFreeEpetraOp::reset(
   const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> >& ops)
{
  block_ops = ops;
}

const Stokhos::VectorOrthogPoly<Epetra_Operator>&
Stokhos::MatrixFreeEpetraOp::getOperatorBlocks()
{
  return *block_ops;
}

int 
Stokhos::MatrixFreeEpetraOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  for (unsigned int i=0; i<num_blocks; i++)
    (*block_ops)[i].SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::MatrixFreeEpetraOp::Apply(const Epetra_MultiVector& Input, 
				   Epetra_MultiVector& Result) const
{
  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (&Input == &Result) {
    input = new Epetra_MultiVector(Input);
    made_copy = true;
  }

  // Initialize
  Result.PutScalar(0.0);

  // Allocate temporary storage
  int m = Input.NumVectors();
  if (tmp == Teuchos::null || tmp->NumVectors() != m)
    tmp = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));

  // Extract blocks
  EpetraExt::BlockMultiVector sg_input(View, *base_map, *input);
  EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);
  for (unsigned int i=0; i<num_blocks; i++) {
    input_block[i] = sg_input.GetBlock(i);
    result_block[i] = sg_result.GetBlock(i);
  }

  // Apply block SG operator via
  // w_i = 
  //    \sum_{j=0}^P \sum_{k=0}^P J_k v_j < \psi_i \psi_j \psi_k > / <\psi_i^2>
  // for i=0,...,P where P = num_blocks w_j is the jth input block, w_i
  // is the ith result block, and J_k is the kth block operator
  const std::vector<double>& norms = sg_basis->norm_squared();
  for (unsigned int k=0; k<num_blocks; k++) {
    unsigned int nj = Cijk->num_j(k);
    const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
    std::vector<double*> j_ptr(nj);
    for (unsigned int l=0; l<nj; l++)
      j_ptr[l] = input_block[j_indices[l]]->Values();
    Epetra_MultiVector input_tmp(View, *base_map, &j_ptr[0], nj);
    Epetra_MultiVector result_tmp(View, *tmp2, const_cast<int*>(&j_indices[0]), nj);
    (*block_ops)[k].Apply(input_tmp, result_tmp);
    for (unsigned int l=0; l<nj; l++) {
      const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,l);
      const Teuchos::Array<double>& c_values = Cijk->values(k,l);
      for (unsigned int i=0; i<i_indices.size(); i++) {
  	int ii = i_indices[i];
  	result_block[ii]->Update(c_values[i]/norms[ii], *result_tmp(l), 1.0);
      }
    }
  }

  // Destroy blocks
  for (unsigned int i=0; i<num_blocks; i++) {
    input_block[i] = Teuchos::null;
    result_block[i] = Teuchos::null;
  }

  if (made_copy)
    delete input;

  return 0;
}

int 
Stokhos::MatrixFreeEpetraOp::ApplyInverse(const Epetra_MultiVector& Input, 
					  Epetra_MultiVector& Result) const
{
  throw "MatrixFreeEpetraOp::ApplyInverse not defined!";
  return -1;
}

double 
Stokhos::MatrixFreeEpetraOp::NormInf() const
{
  return 1.0;
}


const char* 
Stokhos::MatrixFreeEpetraOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::MatrixFreeEpetraOp::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::MatrixFreeEpetraOp::HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
Stokhos::MatrixFreeEpetraOp::Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
Stokhos::MatrixFreeEpetraOp::OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::MatrixFreeEpetraOp::OperatorRangeMap() const
{
  return *sg_map;
}
