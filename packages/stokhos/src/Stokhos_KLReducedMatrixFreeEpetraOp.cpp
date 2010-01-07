// $Id: Stokhos_KLReducedMatrixFreeEpetraOp.cpp,v 1.7 2009/09/14 18:35:48 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/stokhos/src/Stokhos_KLReducedMatrixFreeEpetraOp.cpp,v $ 
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
#include "Stokhos_KLReducedMatrixFreeEpetraOp.hpp"
#include "Stokhos_PCEAnasaziKL.hpp"
#include "Teuchos_TestForException.hpp"

#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_MatrixFreeEpetraOp.hpp"

Stokhos::KLReducedMatrixFreeEpetraOp::KLReducedMatrixFreeEpetraOp(
 const Teuchos::RCP<const Epetra_Map>& base_map_,
 const Teuchos::RCP<const Epetra_Map>& sg_map_,
 const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
 const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
 const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> >& ops_,
 int num_KL_) 
  : label("Stokhos KL Reduced Matrix Free Operator"),
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
    tmp2(),
    num_KL(num_KL_),
    num_KL_computed(0),
    mean(),
    dot_products(),
    kl_coeffs(),
    kl_blocks()
{
  tmp2 = Teuchos::rcp(new Epetra_MultiVector(*base_map, num_blocks));
 
  // Setup KL blockss
  setup();
}

Stokhos::KLReducedMatrixFreeEpetraOp::~KLReducedMatrixFreeEpetraOp()
{
}

void 
Stokhos::KLReducedMatrixFreeEpetraOp::reset(
   const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> >& ops)
{
  block_ops = ops;
  setup();
}

const Stokhos::VectorOrthogPoly<Epetra_Operator>&
Stokhos::KLReducedMatrixFreeEpetraOp::getOperatorBlocks()
{
  return *block_ops;
}

int 
Stokhos::KLReducedMatrixFreeEpetraOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  for (int i=0; i<num_blocks; i++)
    (*block_ops)[i].SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::KLReducedMatrixFreeEpetraOp::Apply(const Epetra_MultiVector& Input, 
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
  if (tmp == Teuchos::null || tmp->NumVectors() != m*num_blocks)
    tmp = Teuchos::rcp(new Epetra_MultiVector(*base_map, num_blocks*m));

  Epetra_MultiVector input_tmp(View, *base_map, input->Values(), 
			       base_map->NumMyElements(), num_blocks*m);
  Epetra_MultiVector result_tmp(View, *base_map, Result.Values(),
				base_map->NumMyElements(), num_blocks*m);

  EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);
  for (int i=0; i<num_blocks; i++)
    result_block[i] = sg_result.GetBlock(i);

  // Apply block SG operator via
  // w_i = 
  //    \sum_{j=0}^P \sum_{k=0}^P J_k v_j < \psi_i \psi_j \psi_k > / <\psi_i^2>
  // for i=0,...,P where P = num_blocks w_j is the jth input block, w_i
  // is the ith result block, and J_k is the kth block operator
  mean->Apply(input_tmp, result_tmp);
  for (int k=0; k<num_KL_computed; k++) {
    kl_blocks[k]->Apply(input_tmp, *tmp2);
    for (int i=0; i<num_blocks; i++)
      for (int j=0; j<num_blocks; j++)
	result_block[i]->Update(kl_coeffs[k][i][j], *(*tmp2)(j), 1.0);
  }

  // Destroy blocks
  for (int i=0; i<num_blocks; i++)
    result_block[i] = Teuchos::null;

  if (made_copy)
    delete input;

  return 0;
}

int 
Stokhos::KLReducedMatrixFreeEpetraOp::ApplyInverse(const Epetra_MultiVector& Input, 
					  Epetra_MultiVector& Result) const
{
  throw "KLReducedMatrixFreeEpetraOp::ApplyInverse not defined!";
  return -1;
}

double 
Stokhos::KLReducedMatrixFreeEpetraOp::NormInf() const
{
  return 1.0;
}


const char* 
Stokhos::KLReducedMatrixFreeEpetraOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::KLReducedMatrixFreeEpetraOp::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::KLReducedMatrixFreeEpetraOp::HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
Stokhos::KLReducedMatrixFreeEpetraOp::Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
Stokhos::KLReducedMatrixFreeEpetraOp::OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::KLReducedMatrixFreeEpetraOp::OperatorRangeMap() const
{
  return *sg_map;
}

void
Stokhos::KLReducedMatrixFreeEpetraOp::setup()
{
#ifdef HAVE_STOKHOS_ANASAZI
  // Build a vector polynomial out of matrix nonzeros
  mean = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(
    block_ops->getCoeffPtr(0));
  Epetra_Map block_vec_map(-1, mean->NumMyNonzeros(), 0, base_map->Comm());
  Stokhos::VectorOrthogPoly<Epetra_Vector> block_vec_poly(
    sg_basis, Stokhos::EpetraVectorCloner(block_vec_map));
  for (int coeff=0; coeff<num_blocks; coeff++) {
    Teuchos::RCP<const Epetra_CrsMatrix> block_coeff = 
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>
        (block_ops->getCoeffPtr(coeff));
    int row = 0;
    for (int i=0; i<mean->NumMyRows(); i++) {
      int num_col;
      mean->NumMyRowEntries(i, num_col);
      for (int j=0; j<num_col; j++)
	block_vec_poly[coeff][row++] = (*block_coeff)[i][j];
    }
  }

  // Compute KL expansion of solution sg_J_vec_poly
  Stokhos::PCEAnasaziKL pceKL(block_vec_poly, num_KL);
  Teuchos::ParameterList anasazi_params = pceKL.getDefaultParams();
  //anasazi_params.set("Num Blocks", 10);
  bool result = pceKL.computeKL(anasazi_params);
  if (!result)
    std::cout << "KL Eigensolver did not converge!" << std::endl;
  //Teuchos::Array<double> evals = pceKL.getEigenvalues();
  Teuchos::RCP<Epetra_MultiVector> evecs = pceKL.getEigenvectors();
  num_KL_computed = evecs->NumVectors();
  std::cout << "num computed evecs = " << num_KL_computed << std::endl;

  // Compute dot products of Jacobian blocks and KL eigenvectors
  dot_products.resize(num_KL_computed);
  for (int rv=0; rv < num_KL_computed; rv++) {
    dot_products[rv].resize(num_blocks-1);
    for (int coeff=1; coeff < num_blocks; coeff++) {
      double dot;
      block_vec_poly[coeff].Dot(*((*evecs)(rv)), &dot);
      dot_products[rv][coeff-1] = dot;
    }
  }

  // Compute KL coefficients
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  kl_coeffs.resize(num_KL_computed);
  for (int k=0; k<num_KL_computed; k++) {
    kl_coeffs[k].resize(num_blocks);
    for (int i=0; i<num_blocks; i++) {
      kl_coeffs[k][i].resize(num_blocks);
      for (int j=0; j<num_blocks; j++)
	kl_coeffs[k][i][j] = 0.0;
    }
  }
  for (int l=1; l<num_blocks; l++) {
    int nj = Cijk->num_j(l);
    const Teuchos::Array<int>& j_indices = Cijk->Jindices(l);
    for (int p=0; p<nj; p++) {
      int j = j_indices[p];
      const Teuchos::Array<int>& i_indices = Cijk->Iindices(l,p);
      const Teuchos::Array<double>& c_values = Cijk->values(l,p);
      for (unsigned int ii=0; ii<i_indices.size(); ii++) {
	int i = i_indices[ii];
	double nrm = norms[i];
	double c  = c_values[ii];
	for (int k=0; k<num_KL_computed; k++) {
	  double dp = dot_products[k][l-1];
	  kl_coeffs[k][i][j] += dp*c/nrm;
	}
      }
    }
  }
    

  // Transform eigenvectors back to matrices
  kl_blocks.resize(num_KL_computed);
  for (int rv=0; rv<num_KL_computed; rv++) {
    kl_blocks[rv] = Teuchos::rcp(new Epetra_CrsMatrix(*mean));
    int row = 0;
    for (int i=0; i<mean->NumMyRows(); i++) {
      int num_col;
      mean->NumMyRowEntries(i, num_col);
      for (int j=0; j<num_col; j++)
	(*kl_blocks[rv])[i][j] = (*evecs)[rv][row++];
    }
  }

  // Check accuracy of KL expansion
  Teuchos::Array<double> point(sg_basis->dimension());
  for (unsigned int i=0; i<sg_basis->dimension(); i++)
    point[i] = 0.5;
  Teuchos::Array<double> basis_vals(num_blocks);
  sg_basis->evaluateBases(point, basis_vals);
  Epetra_Vector val(block_vec_map);
  Epetra_Vector val_kl(block_vec_map);
  block_vec_poly.evaluate(basis_vals, val);
  val_kl.Update(1.0, block_vec_poly[0], 0.0);
  Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > rvs(num_KL_computed);
  Teuchos::Array<double> val_rvs(num_KL_computed);
  for (int rv=0; rv<num_KL_computed; rv++) {
    rvs[rv].reset(sg_basis);
    rvs[rv][0] = 0.0;
    for (int coeff=1; coeff<num_blocks; coeff++)
      rvs[rv][coeff] = dot_products[rv][coeff-1];
    val_rvs[rv] = rvs[rv].evaluate(point, basis_vals);
    val_kl.Update(val_rvs[rv], *((*evecs)(rv)), 1.0);
  } 
  val.Update(-1.0, val_kl, 1.0);
  double diff;
  val.NormInf(&diff);
  std::cout << "Infinity norm of random field difference = " << diff 
	    << std::endl;
  
  // Check accuracy of operator
  Epetra_Vector op_input(*sg_map), op_result(*sg_map), op_kl_result(*sg_map);
  op_input.PutScalar(1.0);
  Stokhos::MatrixFreeEpetraOp op(base_map, sg_map, sg_basis, Cijk, block_ops);
  op.Apply(op_input, op_result);
  this->Apply(op_input, op_kl_result);
  op_result.Update(-1.0, op_kl_result, 1.0);
  op_result.NormInf(&diff);
  std::cout << "Infinity norm of operator difference = " << diff 
	    << std::endl;
#else
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::KLReducedMatrixFreeEpetraOp is available " <<
		     "only when configured with Anasazi support!")
#endif
}
