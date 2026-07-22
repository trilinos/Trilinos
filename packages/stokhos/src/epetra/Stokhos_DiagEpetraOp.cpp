// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Epetra_config.h"
#include "EpetraExt_BlockMultiVector.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Stokhos_DiagEpetraOp.hpp"

Stokhos::DiagEpetraOp::DiagEpetraOp(
 const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
 const Teuchos::RCP<const Epetra_Map>& range_base_map_,
 const Teuchos::RCP<const Epetra_Map>& domain_sg_map_,
 const Teuchos::RCP<const Epetra_Map>& range_sg_map_,
 const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
 const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
 const Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >& ops_) 
  : label("Stokhos Diagonal Operator"),
    domain_base_map(domain_base_map_),
    range_base_map(range_base_map_),
    domain_sg_map(domain_sg_map_),
    range_sg_map(range_sg_map_),
    sg_basis(sg_basis_),
    Cijk(Cijk_),
    block_ops(ops_),
    useTranspose(false),
    expansion_size(sg_basis->size()),
    num_blocks(block_ops->size()),
    input_block(expansion_size),
    result_block(expansion_size),
    tmp(),
    tmp_trans()
{
}

Stokhos::DiagEpetraOp::~DiagEpetraOp()
{
}

void 
Stokhos::DiagEpetraOp::reset(
   const Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >& ops)
{
  block_ops = ops;
}

Teuchos::RCP<const Stokhos::EpetraOperatorOrthogPoly >
Stokhos::DiagEpetraOp::getOperatorBlocks() const
{
  return block_ops;
}

Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >
Stokhos::DiagEpetraOp::getOperatorBlocks()
{
  return block_ops;
}

int 
Stokhos::DiagEpetraOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  for (int i=0; i<num_blocks; i++)
    (*block_ops)[i].SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::DiagEpetraOp::Apply(std::vector< Teuchos::RCP< const Epetra_CrsMatrix> >& sg_J_all, std::vector< Teuchos::RCP< Epetra_CrsMatrix> >& sg_Kkk_all) const
{
  
//  std::vector< Teuchos::RCP< Epetra_CrsMatrix> > sg_Kkk_all ;
//  std::vector< Teuchos::RCP< const Epetra_CrsMatrix> > sg_J_all; 

  for (int i=0;i<num_blocks+1;i++) {
      Teuchos::RCP<Epetra_CrsMatrix> sg_J_poly_Crs =
      	Teuchos::rcp_dynamic_cast< Epetra_CrsMatrix>((*block_ops).getCoeffPtr(i),true);
      sg_J_all.push_back(sg_J_poly_Crs); 
  }

/*  Teuchos::RCP<Epetra_CrsMatrix> sg_J_poly_Crs =
    Teuchos::rcp_dynamic_cast< Epetra_CrsMatrix>((*sg_J_poly).getCoeffPtr(0),true);

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(sg_J_poly_Crs), std::runtime_error,
		     "Dynamic cast of sg_J_poly failed!");
*/

  for(int k=0;k<expansion_size;k++) {
     Teuchos::RCP<Epetra_CrsMatrix> Kkk = 
	Teuchos::rcp(new Epetra_CrsMatrix(*(sg_J_all[0])));  
     sg_Kkk_all.push_back(Kkk);
     sg_Kkk_all[k]->PutScalar(0.0); 
  }

  //Compute diagonal blocks of SG matrix
  for (int k=0; k<expansion_size; k++) {
    int nj = Cijk->num_j(k);
    const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
    for (int jj=0; jj<nj; jj++) {
      int j = j_indices[jj];
      if (j==k) {
      const Teuchos::Array<double>& cijk_values = Cijk->values(k,jj);
      const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
      int ni = i_indices.size();
      for (int ii=0; ii<ni; ii++) {
        int i = i_indices[ii];
        if (i<num_blocks+1) {
         double cikk = cijk_values[ii];  // C(i,j,k)
         EpetraExt::MatrixMatrix::Add((*sg_J_all[i]), false, cikk, *(sg_Kkk_all[k]), 1.0);
        }          
      }
      }
    }
  }

    /* // Apply block SG operator via
  // w_i = 
  //    \sum_{j=0}^P \sum_{k=0}^L J_k v_j < \psi_i \psi_j \psi_k > / <\psi_i^2>
  // for i=0,...,P where P = expansion_size, L = num_blocks, w_j is the jth 
  // input block, w_i is the ith result block, and J_k is the kth block operator
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  for (int k=0; k<num_blocks; k++) {
    int nj = Cijk->num_j(k);
    const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
    Teuchos::Array<double*> j_ptr(nj*m);
    Teuchos::Array<int> mj_indices(nj*m);
    for (int l=0; l<nj; l++) {
      for (int mm=0; mm<m; mm++) {
	j_ptr[l*m+mm] = input_block[j_indices[l]]->Values()+mm*N;
	mj_indices[l*m+mm] = j_indices[l]*m+mm;
      }
    }
    Epetra_MultiVector input_tmp(View, *input_base_map, &j_ptr[0], nj*m);
    Epetra_MultiVector result_tmp(View, *tmp_result, &mj_indices[0], nj*m);
    (*block_ops)[k].Apply(input_tmp, result_tmp);
    for (int l=0; l<nj; l++) {
      const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,l);
      const Teuchos::Array<double>& c_values = Cijk->values(k,l);
      for (int i=0; i<i_indices.size(); i++) {
  	int ii = i_indices[i];
	for (int mm=0; mm<m; mm++)
	  (*result_block[ii])(mm)->Update(c_values[i]/norms[ii], 
					  *result_tmp(l*m+mm), 1.0);
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
    */
  return 0;
}

int 
Stokhos::DiagEpetraOp::ApplyInverse(const Epetra_MultiVector& Input, 
					  Epetra_MultiVector& Result) const
{
  throw "DiagEpetraOp::ApplyInverse not defined!";
  return -1;
}

double 
Stokhos::DiagEpetraOp::NormInf() const
{
  return 1.0;
}


const char* 
Stokhos::DiagEpetraOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::DiagEpetraOp::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::DiagEpetraOp::HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
Stokhos::DiagEpetraOp::Comm() const
{
  return domain_base_map->Comm();
}
const Epetra_Map& 
Stokhos::DiagEpetraOp::OperatorDomainMap() const
{
  return *domain_sg_map;
}

const Epetra_Map& 
Stokhos::DiagEpetraOp::OperatorRangeMap() const
{
  return *range_sg_map;
}
