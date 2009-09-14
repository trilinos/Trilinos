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

#include "Epetra_config.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_KLMatrixFreeEpetraOp.hpp"

Stokhos::KLMatrixFreeEpetraOp::KLMatrixFreeEpetraOp(
 const Teuchos::RCP<const Epetra_Map>& base_map_,
 const Teuchos::RCP<const Epetra_Map>& sg_map_,
 const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
 const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
 const Teuchos::Array<Teuchos::RCP<Epetra_CrsMatrix> >& ops_) 
  : label("Stokhos::MatrixFreeEpetraOp"),
    base_map(base_map_),
    sg_map(sg_map_),
    sg_basis(sg_basis_),
    Cijk(Cijk_),
    block_ops(ops_),
    useTranspose(false),
    num_blocks(sg_basis->size()),
    input_block(num_blocks),
    result_block(num_blocks),
    tmp()
{
}

Stokhos::KLMatrixFreeEpetraOp::~KLMatrixFreeEpetraOp()
{
}

void 
Stokhos::KLMatrixFreeEpetraOp::reset(
   const Teuchos::Array<Teuchos::RCP<Epetra_CrsMatrix> >& ops)
{
  block_ops = ops;
}

const Teuchos::Array<Teuchos::RCP<Epetra_CrsMatrix> >&
Stokhos::KLMatrixFreeEpetraOp::getOperatorBlocks()
{
  return block_ops;
}

int 
Stokhos::KLMatrixFreeEpetraOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  for (unsigned int i=0; i<num_blocks; i++)
    (block_ops)[i]->SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::KLMatrixFreeEpetraOp::Apply(const Epetra_MultiVector& Input, 
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
  
  //Compute K_i*x_k for all i and k.
  int d = sg_basis->dimension();
  Teuchos::Array<Teuchos::RCP<Epetra_MultiVector> > blockProducts(d+1);
  for(int k = 0; k<=d; k++){
    blockProducts[k] = Teuchos::rcp(new Epetra_MultiVector(*base_map,num_blocks)); 
    for(unsigned int i = 0; i<num_blocks; i++){
     (block_ops)[k]->Apply(*input_block[i],*(*blockProducts[k])(i));
    }
  }
  
  double cijk, gamma;
  int i, k;
  Teuchos::Array< double > norms = sg_basis->norm_squared();
  Teuchos::Array<double> one(d);
  for(int j = 0; j<d; j++)one[j] = 1;
  Teuchos::Array< double > values(num_blocks);
  sg_basis->evaluateBases(one, values);
    
  for (int j=0; j<=d; j++) {
    gamma = values[j]/(1-sg_basis->evaluateZero(j));
    unsigned int nl = Cijk->num_values(j);
    for (unsigned int l=0; l<nl; l++) {
      Cijk->value(j, l, i, k, cijk);
      if( j!=0 ){
        cijk /= gamma;
        if( i == k) cijk = cijk -sg_basis->evaluateZero(j)*norms[i];
      }else{
        if( i == k ) cijk = norms[i];
      }
      result_block[i]->Update(cijk, *(*blockProducts[j])(k), 1.0);
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
Stokhos::KLMatrixFreeEpetraOp::ApplyInverse(const Epetra_MultiVector& Input, 
					  Epetra_MultiVector& Result) const
{
  throw "MatrixFreeEpetraOp::ApplyInverse not defined!";
  return -1;
}

double 
Stokhos::KLMatrixFreeEpetraOp::NormInf() const
{
  return 1.0;
}


const char* 
Stokhos::KLMatrixFreeEpetraOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::KLMatrixFreeEpetraOp::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::KLMatrixFreeEpetraOp::HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
Stokhos::KLMatrixFreeEpetraOp::Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
Stokhos::KLMatrixFreeEpetraOp::OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::KLMatrixFreeEpetraOp::OperatorRangeMap() const
{
  return *sg_map;
}
