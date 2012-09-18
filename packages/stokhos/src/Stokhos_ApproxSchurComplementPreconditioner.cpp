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

#include "Stokhos_ApproxSchurComplementPreconditioner.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include <algorithm>
#include <iostream>
#include <iterator>

Stokhos::ApproxSchurComplementPreconditioner::
ApproxSchurComplementPreconditioner(
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk_,
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& prec_factory_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  label("Stokhos Approximate Schur Complement Preconditioner"),
  sg_comm(sg_comm_),
  sg_basis(sg_basis_),
  epetraCijk(epetraCijk_),
  base_map(base_map_),
  sg_map(sg_map_),
  prec_factory(prec_factory_),
  mean_prec(),
  useTranspose(false),
  sg_op(),
  sg_poly(),
  Cijk(epetraCijk->getParallelCijk()),
  P(sg_basis->order()),
  block_indices(P+2),
  upper_block_Cijk(P+1),
  lower_block_Cijk(P+1),
  scale_op(true),
  symmetric(false),
  only_use_linear(true),
  rhs_block()
{
  // Check if parallel, which we don't support
  TEUCHOS_TEST_FOR_EXCEPTION(
    epetraCijk->isStochasticParallel(), std::logic_error, 
    "Stokhos::ApproxSchurComplementPreconditioner does not support " << 
    "a parallel stochastic distribution.");

  scale_op = params_->get("Scale Operator by Inverse Basis Norms", true);
  symmetric = params_->get("Symmetric Gauss-Seidel", false);
  only_use_linear = params_->get("Only Use Linear Terms", true);

  Cijk_type::k_iterator k_begin = Cijk->k_begin();
  Cijk_type::k_iterator k_end = Cijk->k_end();
  if (only_use_linear)
    k_end = Cijk->find_k(sg_basis()->dimension() + 1);

  max_num_mat_vec = 0;
  for (Cijk_type::k_iterator k=k_begin; k!=k_end; ++k) {
    int nj = Cijk->num_j(k);
    if (max_num_mat_vec < nj)
      max_num_mat_vec = nj;
  }

  // Get indices for each block
  Teuchos::RCP<const Stokhos::ProductBasis<int,double> > prod_basis =
    Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int,double> >(sg_basis, true);
  int d = prod_basis->dimension();
  Teuchos::Array<int> term(d);
  for (int p=0; p<=P; p++) {
    term[0] = p;
    block_indices[p] = prod_basis->getIndex(term);
    upper_block_Cijk[p] = Teuchos::rcp(new Cijk_type);
    lower_block_Cijk[p] = Teuchos::rcp(new Cijk_type);
  }
  block_indices[P+1] = sg_basis->size();

  // std::cout << "block_indices = [";
  // std::copy(block_indices.begin(), block_indices.end(), 
  // 	    std::ostream_iterator<int>(std::cout, " "));
  // std::cout << "]" << std::endl;

  // Build Cijk tensors for each order block
  for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
    int k = index(k_it);
    Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
    Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
    for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
      int j = index(j_it);
      Teuchos::Array<int>::iterator col_it = 
	std::upper_bound(block_indices.begin(), block_indices.end(), j);
      int p_col = col_it - block_indices.begin() - 1;
      Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
      Cijk_type::kji_iterator i_end = Cijk->i_end(j_it);
      for (Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
	int i = index(i_it);
	double c = value(i_it);
	Teuchos::Array<int>::iterator row_it = 
	  std::upper_bound(block_indices.begin(), block_indices.end(), i);
	int p_row = row_it - block_indices.begin() - 1;
	//std::cout << "i = " << i << ", p_row = " << p_row << ", j = " << j << ", p_col = " << p_col;
	if (p_col > p_row) {
	  upper_block_Cijk[p_col]->add_term(i,j,k,c);
	  //std::cout << " upper" << std::endl;
	}
	else if (p_col < p_row) {
	  lower_block_Cijk[p_row]->add_term(i,j,k,c);
	  //std::cout << " lower" << std::endl;
	}
      }
    }
  }
  for (int p=0; p<=P; p++) {
    upper_block_Cijk[p]->fillComplete();
    lower_block_Cijk[p]->fillComplete();
  }
}

Stokhos::ApproxSchurComplementPreconditioner::
~ApproxSchurComplementPreconditioner()
{
}

void
Stokhos::ApproxSchurComplementPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op_, 
		    const Epetra_Vector& x)
{
  sg_op = sg_op_;
  sg_poly = sg_op->getSGPolynomial();
  mean_prec = prec_factory->compute(sg_poly->getCoeffPtr(0));
  label = std::string("Stokhos Approximate Schur Complement Preconditioner:\n") 
    + std::string("		***** ") + std::string(mean_prec->Label());
}

int 
Stokhos::ApproxSchurComplementPreconditioner::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  sg_op->SetUseTranspose(UseTranspose);
  mean_prec->SetUseTranspose(UseTranspose);

  return 0;
}

int 
Stokhos::ApproxSchurComplementPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  return sg_op->Apply(Input, Result);
}

int 
Stokhos::ApproxSchurComplementPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Approximate Schur Complement Time");
#endif

  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (Input.Values() == Result.Values()) {
    input = new Epetra_MultiVector(Input);
    made_copy = true;
  } 

  // Allocate temporary storage
  int m = input->NumVectors();
  if (rhs_block == Teuchos::null || rhs_block->NumVectors() != m)
    rhs_block = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
  if (tmp == Teuchos::null || tmp->NumVectors() != m*max_num_mat_vec)
    tmp = Teuchos::rcp(new Epetra_MultiVector(*base_map, 
					      m*max_num_mat_vec));
  j_ptr.resize(m*max_num_mat_vec);
  mj_indices.resize(m*max_num_mat_vec);
  
  // Extract blocks
  EpetraExt::BlockMultiVector input_block(View, *base_map, *input);
  EpetraExt::BlockMultiVector result_block(View, *base_map, Result);

  result_block.PutScalar(0.0);

  // Set right-hand-side to input_block
  rhs_block->Update(1.0, input_block, 0.0);

  // At level l, linear system has the structure
  // [ A_{l-1} B_l ][ u_l^{l-1} ] = [ r_l^{l-1} ]
  // [ C_l     D_l ][ u_l^l     ]   [ r_l^l     ]

  for (int l=P; l>=1; l--) {
    // Compute D_l^{-1} r_l^l
    divide_diagonal_block(block_indices[l], block_indices[l+1], 
			  *rhs_block, result_block);

    // Compute r_l^{l-1} = r_l^{l-1} - B_l D_l^{-1} r_l^l
    multiply_block(upper_block_Cijk[l], -1.0, result_block, *rhs_block);
  }

  // Solve A_0 u_0 = r_0
  divide_diagonal_block(0, 1, *rhs_block, result_block);

  for (int l=1; l<=P; l++) {
    // Compute r_l^l - C_l*u_l^{l-1}
    multiply_block(lower_block_Cijk[l], -1.0, result_block, *rhs_block);

    // Compute D_l^{-1} (r_l^l - C_l*u_l^{l-1})
    divide_diagonal_block(block_indices[l], block_indices[l+1], 
			  *rhs_block, result_block);
  }

  if (made_copy)
    delete input;

  return 0; 
}

double 
Stokhos::ApproxSchurComplementPreconditioner::
NormInf() const
{
  return sg_op->NormInf();
}


const char* 
Stokhos::ApproxSchurComplementPreconditioner::
Label() const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::ApproxSchurComplementPreconditioner::
UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::ApproxSchurComplementPreconditioner::
HasNormInf() const
{
  return sg_op->HasNormInf();
}

const Epetra_Comm & 
Stokhos::ApproxSchurComplementPreconditioner::
Comm() const
{
  return *sg_comm;
}
const Epetra_Map& 
Stokhos::ApproxSchurComplementPreconditioner::
OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::ApproxSchurComplementPreconditioner::
OperatorRangeMap() const
{
  return *sg_map;
}

void
Stokhos::ApproxSchurComplementPreconditioner::
multiply_block(
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& cijk,
  double alpha,
  const EpetraExt::BlockMultiVector& Input, 
  EpetraExt::BlockMultiVector& Result) const
{
  // Input and Result are the whole vector/multi-vector, not just the portion
  // needed for the particular sub-block
  int m = Input.NumVectors();
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  Cijk_type::k_iterator k_begin = cijk->k_begin();
  Cijk_type::k_iterator k_end = cijk->k_end();
  if (only_use_linear)
    k_end = cijk->find_k(sg_basis()->dimension() + 1);
  for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
    int k = index(k_it);
    Cijk_type::kj_iterator j_begin = cijk->j_begin(k_it);
    Cijk_type::kj_iterator j_end = cijk->j_end(k_it);
    int nj = cijk->num_j(k_it);
    if (nj > 0) {
      int l = 0;
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
        int j = index(j_it);
	for (int mm=0; mm<m; mm++) {
	  j_ptr[l*m+mm] = (*(Input.GetBlock(j)))[mm];
	  mj_indices[l*m+mm] = l*m+mm;
	}
	l++;
      }
      Epetra_MultiVector input_tmp(View, *base_map, &j_ptr[0], nj*m);
      Epetra_MultiVector result_tmp(View, *tmp, &mj_indices[0], nj*m);
      (*sg_poly)[k].Apply(input_tmp, result_tmp);
      l = 0;
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	Cijk_type::kji_iterator i_begin = cijk->i_begin(j_it);
	Cijk_type::kji_iterator i_end = cijk->i_end(j_it);
	for (Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
	  int i = index(i_it);
	  double c = value(i_it);
	  if (scale_op)
	    c /= norms[i];
	  for (int mm=0; mm<m; mm++)
	    (*Result.GetBlock(i))(mm)->Update(alpha*c, *result_tmp(l*m+mm), 1.0);
	}
	l++;
      }
    }
  }
}

void
Stokhos::ApproxSchurComplementPreconditioner::
divide_diagonal_block(int row_begin, int row_end, 
		      const EpetraExt::BlockMultiVector& Input, 
		      EpetraExt::BlockMultiVector& Result) const
{
  for (int i=row_begin; i<row_end; i++) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
      TEUCHOS_FUNC_TIME_MONITOR(
	"Stokhos: ASC Deterministic Preconditioner Time");
#endif
      mean_prec->ApplyInverse(*(Input.GetBlock(i)), *(Result.GetBlock(i)));
  }
}
