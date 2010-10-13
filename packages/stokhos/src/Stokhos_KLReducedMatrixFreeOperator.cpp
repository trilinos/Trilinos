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
#include "Stokhos_KLReducedMatrixFreeOperator.hpp"
#include "Stokhos_PCEAnasaziKL.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_MatrixFreeOperator.hpp"
#include "Stokhos_Sparse3TensorUtilities.hpp"
#include <sstream>

Stokhos::KLReducedMatrixFreeOperator::
KLReducedMatrixFreeOperator(
 const Teuchos::RCP<const Epetra_Map>& domain_base_map_,
 const Teuchos::RCP<const Epetra_Map>& range_base_map_,
 const Teuchos::RCP<const Epetra_Map>& domain_sg_map_,
 const Teuchos::RCP<const Epetra_Map>& range_sg_map_,
 const Teuchos::RCP<Teuchos::ParameterList>& params_) : 
  label("Stokhos KL Reduced Matrix Free Operator"),
  domain_base_map(domain_base_map_),
  range_base_map(range_base_map_),
  domain_sg_map(domain_sg_map_),
  range_sg_map(range_sg_map_),
  sg_basis(),
  Cijk(),
  block_ops(),
  params(params_),
  useTranspose(false),
  expansion_size(0),
  num_blocks(0),
  num_KL(0),
  num_KL_computed(0),
  mean(),
  block_vec_map(),
  block_vec_poly(),
  dot_products(),
  sparse_kl_coeffs(),
  kl_blocks()
{
  num_KL = params->get("Number of KL Terms", 5);
  drop_tolerance = params->get("Sparse 3 Tensor Drop Tolerance", 1e-6);
  do_error_tests = params->get("Do Error Tests", false);
}

void 
Stokhos::KLReducedMatrixFreeOperator::
setupOperator(
   const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> >& ops,
   const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_)
{
  block_ops = ops;
  Cijk = Cijk_;
  sg_basis = block_ops->basis();
  expansion_size = sg_basis->size();
  num_blocks = block_ops->size();

  // Build a vector polynomial out of matrix nonzeros
  mean = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(
    block_ops->getCoeffPtr(0));
  block_vec_map = 
    Teuchos::rcp(new Epetra_Map(-1, mean->NumMyNonzeros(), 0, 
				domain_base_map->Comm()));
  block_vec_poly = Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_Vector>(
				  sg_basis, 
				  Stokhos::EpetraVectorCloner(*block_vec_map)));
  
  // Setup KL blocks
  setup();
}

Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > 
Stokhos::KLReducedMatrixFreeOperator::
getSGPolynomial()
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::VectorOrthogPoly<Epetra_Operator> > 
Stokhos::KLReducedMatrixFreeOperator::
getSGPolynomial() const
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > 
Stokhos::KLReducedMatrixFreeOperator::
getTripleProduct() const
{
  return Cijk;
}

Stokhos::KLReducedMatrixFreeOperator::
~KLReducedMatrixFreeOperator()
{
}

int 
Stokhos::KLReducedMatrixFreeOperator::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  for (int i=0; i<num_blocks; i++)
    (*block_ops)[i].SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::KLReducedMatrixFreeOperator::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  return kl_mat_free_op->Apply(Input, Result);
}

int 
Stokhos::KLReducedMatrixFreeOperator::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  throw "KLReducedMatrixFreeOperator::ApplyInverse not defined!";
  return -1;
}

double 
Stokhos::KLReducedMatrixFreeOperator::
NormInf() const
{
  return 1.0;
}


const char* 
Stokhos::KLReducedMatrixFreeOperator::
Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::KLReducedMatrixFreeOperator::
UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::KLReducedMatrixFreeOperator::
HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
Stokhos::KLReducedMatrixFreeOperator::
Comm() const
{
  return domain_base_map->Comm();
}
const Epetra_Map& 
Stokhos::KLReducedMatrixFreeOperator::
OperatorDomainMap() const
{
  return *domain_sg_map;
}

const Epetra_Map& 
Stokhos::KLReducedMatrixFreeOperator::
OperatorRangeMap() const
{
  return *range_sg_map;
}

void
Stokhos::KLReducedMatrixFreeOperator::
setup()
{
#ifdef HAVE_STOKHOS_ANASAZI
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::KLReducedMatrixFreeOperator -- Calculation/setup of KL opeator");
#endif

  mean = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(
    block_ops->getCoeffPtr(0));

  // Copy matrix coefficients into vectors
  for (int coeff=0; coeff<num_blocks; coeff++) {
    Teuchos::RCP<const Epetra_CrsMatrix> block_coeff = 
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>
        (block_ops->getCoeffPtr(coeff));
    int row = 0;
    for (int i=0; i<mean->NumMyRows(); i++) {
      int num_col;
      mean->NumMyRowEntries(i, num_col);
      for (int j=0; j<num_col; j++)
	(*block_vec_poly)[coeff][row++] = (*block_coeff)[i][j];
    }
  }

  // Compute KL expansion of solution sg_J_vec_poly
  Stokhos::PCEAnasaziKL pceKL(*block_vec_poly, num_KL);
  Teuchos::ParameterList anasazi_params = pceKL.getDefaultParams();
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
      (*block_vec_poly)[coeff].Dot(*((*evecs)(rv)), &dot);
      dot_products[rv][coeff-1] = dot;
    }
  }

  // Compute KL coefficients
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  sparse_kl_coeffs = 
    //Teuchos::rcp(new Stokhos::Sparse3Tensor<int,double>(num_KL_computed+1));
    Teuchos::rcp(new Stokhos::Sparse3Tensor<int,double>(expansion_size));
  for (int i=0; i<expansion_size; i++)
    sparse_kl_coeffs->sum_term(i, i, 0, norms[i]);
  for (int l=1; l<num_blocks; l++) {
    for (Cijk_type::kj_iterator j_it = Cijk->j_begin(l); 
	 j_it != Cijk->j_end(l); ++j_it) {
      int j = index(j_it);
      for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	   i_it != Cijk->i_end(j_it); ++i_it) {
	int i = index(i_it);
	double c  = value(i_it);
	for (int k=1; k<num_KL_computed+1; k++) {
	  double dp = dot_products[k-1][l-1];
	  double v = dp*c;
	  if (std::abs(v) > drop_tolerance)
	    sparse_kl_coeffs->sum_term(i,j,k,v);
	}
      }
    }
  }
   
  bool save_tensor = params->get("Save Sparse 3 Tensor To File", false);
  if (save_tensor) {
    static int idx = 0;
    std::string basename = params->get("Sparse 3 Tensor Base Filename", 
				       "sparse_KL_coeffs");
    std::stringstream ss;
    ss << basename << "_" << idx++ << ".mm";
    sparse3Tensor2MatrixMarket(*sg_basis, *sparse_kl_coeffs, 
			       domain_base_map->Comm(), ss.str());
  }

  // Transform eigenvectors back to matrices
  kl_blocks.resize(num_KL_computed);
  kl_ops = 
    Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_Operator>(
		   sg_basis, num_KL_computed+1));
  kl_ops->setCoeffPtr(0, mean);
  for (int rv=0; rv<num_KL_computed; rv++) {
    if (kl_blocks[rv] == Teuchos::null) 
      kl_blocks[rv] = Teuchos::rcp(new Epetra_CrsMatrix(*mean));
    int row = 0;
    for (int i=0; i<mean->NumMyRows(); i++) {
      int num_col;
      mean->NumMyRowEntries(i, num_col);
      for (int j=0; j<num_col; j++)
	(*kl_blocks[rv])[i][j] = (*evecs)[rv][row++];
    }
    kl_ops->setCoeffPtr(rv+1, kl_blocks[rv]);
  }
  if (kl_mat_free_op == Teuchos::null)
    kl_mat_free_op = Teuchos::rcp(new Stokhos::MatrixFreeOperator(domain_base_map, range_base_map, domain_sg_map, range_sg_map));
  kl_mat_free_op->setupOperator(kl_ops, sparse_kl_coeffs);

  // Check accuracy of KL expansion
  if (do_error_tests) {
    Teuchos::Array<double> point(sg_basis->dimension());
    for (int i=0; i<sg_basis->dimension(); i++)
      point[i] = 0.5;
    Teuchos::Array<double> basis_vals(num_blocks);
    sg_basis->evaluateBases(point, basis_vals);
    Epetra_Vector val(*block_vec_map);
    Epetra_Vector val_kl(*block_vec_map);
    block_vec_poly->evaluate(basis_vals, val);
    val_kl.Update(1.0, (*block_vec_poly)[0], 0.0);
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
    Epetra_Vector op_input(*domain_sg_map), op_result(*range_sg_map), op_kl_result(*range_sg_map);
    op_input.PutScalar(1.0);
    Stokhos::MatrixFreeOperator op(domain_base_map, range_base_map, domain_sg_map, range_sg_map, params);
    op.setupOperator(block_ops, Cijk);
    op.Apply(op_input, op_result);
    this->Apply(op_input, op_kl_result);
    op_result.Update(-1.0, op_kl_result, 1.0);
    op_result.NormInf(&diff);
    std::cout << "Infinity norm of operator difference = " << diff 
	      << std::endl;
  }

#else
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::KLReducedMatrixFreeOperator is available " <<
		     "only when configured with Anasazi support!")
#endif
}
