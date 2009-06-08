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
    sg_input(),
    sg_result(),
    input_block(num_blocks),
    result_block(num_blocks),
    tmp()
{
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
  int m = Input.NumVectors();
  if (sg_input == Teuchos::null || sg_input->NumVectors() != m) {
    sg_input = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
    sg_result = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
    for (unsigned int i=0; i<num_blocks; i++) {
      input_block[i] = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));
      result_block[i] = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));
    }
    tmp = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));
  }

  // Fill input blocks
  sg_input->Scale(1.0, Input);
  for (unsigned int i=0; i<num_blocks; i++) {
    sg_input->ExtractBlockValues(*input_block[i], i);
    result_block[i]->PutScalar(0.0);
  }

  // Apply block SG operator via
  // w_i = 
  //    \sum_{j=0}^P \sum_{k=0}^P J_k v_j < \psi_i \psi_j \psi_k > / <\psi_i^2>
  // for i=0,...,P where P = num_blocks w_j is the jth input block, w_i
  // is the ith result block, and J_k is the kth block operator
  double cijk;
  int i, j;
  for (unsigned int k=0; k<num_blocks; k++) {
    unsigned int nl = Cijk->num_values(k);
    for (unsigned int l=0; l<nl; l++) {
      Cijk->value(k, l, i, j, cijk);
      cijk /= sg_basis->norm_squared(i);
      (*block_ops)[k].Apply(*input_block[j], *tmp);
      result_block[i]->Update(cijk, *tmp, 1.0);
    }
  }

  // Get result from blocks
  for (unsigned int i=0; i<num_blocks; i++)
    sg_result->LoadBlockValues(*result_block[i], i);
  Result.Scale(1.0, *sg_result);

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
