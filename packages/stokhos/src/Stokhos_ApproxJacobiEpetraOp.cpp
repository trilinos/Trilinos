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

#include "Epetra_config.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_ApproxJacobiEpetraOp.hpp"

Stokhos::ApproxJacobiEpetraOp::ApproxJacobiEpetraOp(
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  unsigned int num_blocks_,
  const Teuchos::RCP<Epetra_Operator>& mean_prec_op_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
  const Teuchos::RCP<Epetra_Operator>& J,
  bool symmetric_):
  label("Stokhos Approximate Jacobi Preconditioner"),
  base_map(base_map_),
  sg_map(sg_map_),
  useTranspose(false),
  num_blocks(num_blocks_),
  mean_prec_op(mean_prec_op_),
  Cijk(Cijk_),
  symmetric(symmetric_),
  mat_vec_tmp()
{
  stokhos_op = Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(J, true);
  sg_J_poly = stokhos_op->getOperatorBlocks();

  if (mean_prec_op != Teuchos::null)
    label = std::string("Stokhos Approximate Jacobi Preconditioner:\n") + 
      std::string("		***** ") + 
      std::string(mean_prec_op->Label());
}

Stokhos::ApproxJacobiEpetraOp::~ApproxJacobiEpetraOp()
{
}

void
Stokhos::ApproxJacobiEpetraOp::setOperator(
  const Teuchos::RCP<Epetra_Operator>& J) 
{
  stokhos_op = Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(J, true);
  sg_J_poly = stokhos_op->getOperatorBlocks();
}

void
Stokhos::ApproxJacobiEpetraOp::setPrecOperator(
  const Teuchos::RCP<Epetra_Operator>& M) 
{
  mean_prec_op = M;
  label = std::string("Stokhos Approximate Jacobi Preconditioner:\n") + 
      std::string("		***** ") + 
      std::string(mean_prec_op->Label());
}

int 
Stokhos::ApproxJacobiEpetraOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;

  return 0;
}

int 
Stokhos::ApproxJacobiEpetraOp::Apply(const Epetra_MultiVector& Input, 
					  Epetra_MultiVector& Result) const
{
  return stokhos_op->Apply(Input,Result);
}

int 
Stokhos::ApproxJacobiEpetraOp::ApplyInverse(
  const Epetra_MultiVector& Input, 
  Epetra_MultiVector& Result) const
{
  TEUCHOS_FUNC_TIME_MONITOR("Total Approximate Jacobi Time");

  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (Input.Values() == Result.Values()) {
    input = new Epetra_MultiVector(Input);
    made_copy = true;
  } 

  int m = input->NumVectors();
  if (mat_vec_tmp == Teuchos::null || mat_vec_tmp->NumVectors() != m)
    mat_vec_tmp = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));

  // Extract blocks
  EpetraExt::BlockMultiVector input_block(View, *base_map, *input);
  EpetraExt::BlockMultiVector result_block(View, *base_map, Result);
  EpetraExt::BlockMultiVector rhs_block(*base_map, *sg_map, m);
//  EpetraExt::BlockMultiVector kx_block(Copy, *base_map, *input);
  std::vector< Teuchos::RCP< Epetra_MultiVector> > sg_kx_vec_all ;

  int sz = sg_J_poly->basis()->size();
  for (int i=0;i<sz;i++) 
    sg_kx_vec_all.push_back(Teuchos::rcp(new Epetra_MultiVector(*base_map,1)));

  result_block.PutScalar(0.0);

  int numKL = sg_J_poly->basis()->dimension();
  const Teuchos::Array<double>& norms = sg_J_poly->basis()->norm_squared();
//  int i,j,nl;
//  double c;
  
  for (int iter=0; iter<2; iter++) {
    rhs_block.Update(1.0, input_block, 0.0);
    if (iter !=0 ) {
  for (int k=1; k<sg_J_poly->size(); k++) {
   if (k<= numKL+1) {
    int nj = Cijk->num_j(k);
    const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
    for (int jj=0; jj<nj; jj++) {
      int j = j_indices[jj];
      //(*sg_J_poly)[k].Apply(*(result_block.GetBlock(j)),*(kx_block.GetBlock(j)));
      (*sg_J_poly)[k].Apply(*(result_block.GetBlock(j)),*(sg_kx_vec_all[j]));
//      (*sg_J_poly)[k].Apply(*(result_block.GetBlock(j)),*mat_vec_tmp);
    }
    for (int jj=0; jj<nj; jj++) {
      int j = j_indices[jj];
      const Teuchos::Array<double>& cijk_values = Cijk->values(k,jj);
      const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
      int ni = i_indices.size();
      for (int ii=0; ii<ni; ii++) {
        int i = i_indices[ii];
        double c = cijk_values[ii];  // C(i,j,k)
        //rhs_block.GetBlock(i)->Update(-1.0*c/norms[i],*(kx_block.GetBlock(j)),1.0);
        rhs_block.GetBlock(i)->Update(-1.0*c/norms[i],*(sg_kx_vec_all[j]),1.0);
       // rhs_block.GetBlock(i)->Update(-1.0*c/norms[i],*mat_vec_tmp,1.0);
      }
    }
   } //if(k<=numKL+1) loop
  } //End of k loop
  }// if(iter!=0) loop 
  for(int i=0; i<sz; i++) {
    // Apply deterministic preconditioner
    TEUCHOS_FUNC_TIME_MONITOR("Total Deterministic Solve Time");
     mean_prec_op->ApplyInverse(*(rhs_block.GetBlock(i)),
                                *(result_block.GetBlock(i)));
  }
 }
/*  rhs_block.Update(1.0, input_block, 0.0);
  for (int k=0; k<sz; k++) {
    nl = Cijk->num_values(k);
    for (int l=0; l<nl; l++) {
      Cijk->value(k,l,i,j,c); 
      if (i!=0) {
	(*sg_J_poly)[i].Apply(*(result_block.GetBlock(j)), *mat_vec_tmp);
	rhs_block.GetBlock(k)->Update(-1.0*c/norms[k], *mat_vec_tmp, 1.0);
      }
    }
      
    result_block.GetBlock(k)->PutScalar(0.0);
    {
      // Apply deterministic preconditioner
      TEUCHOS_FUNC_TIME_MONITOR("Total Deterministic Solve Time");
      mean_prec_op->ApplyInverse(*(rhs_block.GetBlock(k)), 
				 *(result_block.GetBlock(k)));
    }
    
  }
*/
/*  // For symmetric Gauss-Seidel
  if (symmetric) {
    
    rhs_block.Update(1.0, input_block, 0.0);
    for (int k=sz-1; k>=0; k--) {
      nl = Cijk->num_values(k);
      for (int l=0; l<nl; l++) {
	Cijk->value(k,l,i,j,c); 
	if (i!=0) {
	  (*sg_J_poly)[i].Apply(*(result_block.GetBlock(j)), *mat_vec_tmp);
	  rhs_block.GetBlock(k)->Update(-1.0*c/norms[k], *mat_vec_tmp, 1.0);
	}
      }
      
      {
	// Apply deterministic preconditioner
	TEUCHOS_FUNC_TIME_MONITOR("Total Deterministic Solve Time");
	mean_prec_op->ApplyInverse(*(rhs_block.GetBlock(k)), 
				   *(result_block.GetBlock(k)));
      }
    
    }
  }
*/
  if (made_copy)
    delete input;

  return 0; 
}

double 
Stokhos::ApproxJacobiEpetraOp::NormInf() const
{
  return stokhos_op->NormInf();
}


const char* 
Stokhos::ApproxJacobiEpetraOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::ApproxJacobiEpetraOp::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::ApproxJacobiEpetraOp::HasNormInf() const
{
  return stokhos_op->HasNormInf();
}

const Epetra_Comm & 
Stokhos::ApproxJacobiEpetraOp::Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
Stokhos::ApproxJacobiEpetraOp::OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::ApproxJacobiEpetraOp::OperatorRangeMap() const
{
  return *sg_map;
}
