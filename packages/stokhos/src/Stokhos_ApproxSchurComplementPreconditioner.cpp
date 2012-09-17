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
  scale_op(true),
  symmetric(false),
  only_use_linear(true),
  rhs_block()
{
  scale_op = params_->get("Scale Operator by Inverse Basis Norms", true);
  symmetric = params_->get("Symmetric Gauss-Seidel", false);
  only_use_linear = params_->get("Only Use Linear Terms", true);

  k_begin = Cijk->k_begin();
  k_end = Cijk->k_end();
  if (only_use_linear)
    k_end = Cijk->find_k(sg_basis()->dimension() + 1);

  max_num_mat_vec = 0;
  for (Cijk_type::k_iterator k=k_begin; k!=k_end; ++k) {
    int nj = Cijk->num_j(k);
    if (max_num_mat_vec < nj)
      max_num_mat_vec = nj;
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
  label = std::string("Stokhos Approximate Gauss-Seidel Preconditioner:\n") + 
    std::string("		***** ") + 
    std::string(mean_prec->Label());
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

/*
int 
Stokhos::ApproxSchurComplementPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Approximate Gauss-Seidel Time");
#endif

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
  if (rhs_block == Teuchos::null || rhs_block->NumVectors() != m)
    rhs_block = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
  
  // Extract blocks
  EpetraExt::BlockMultiVector input_block(View, *base_map, *input);
  EpetraExt::BlockMultiVector result_block(View, *base_map, Result);

  result_block.PutScalar(0.0);

  int k_limit = sg_poly->size();
  if (only_use_linear)
    k_limit = sg_poly->basis()->dimension() + 1;
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();

  // Set right-hand-side to input_block
  rhs_block->Update(1.0, input_block, 0.0);

  // Preconditioner is recursive on the total polynomial order
  int P = sg_basis->order();
  int n_begin, n_end;
  Teuchos::RCP<Epetra_MultiVector> res_i;
  for (int l=P; l>=1; l--) {

    // At level l, linear system has the structure
    // [ A_{l-1} B_l ][ u_l^{l-1} ] = [ r_l^{l-1} ]
    // [ C_l     D_l ][ u_l^l     ]   [ r_l^l     ]
    
    // Overwrite r_l^{l-1} = r_l^{l-1} - B_l D_l^{-1} r_l^l
    // i here refers to the column of B_l
    compute_index(l, n_begin, n_end);
    //std::cout << "n_begin = " << n_begin << " n_end = " << n_end << std::endl;
    Cijk_type::i_iterator i_begin = Cijk->find_i(n_begin);
    Cijk_type::i_iterator i_end = Cijk->find_i(n_end);
    for (Cijk_type::i_iterator i_it = i_begin; i_it != i_end; ++i_it) {
      int i = index(i_it);
      res_i = result_block.GetBlock(i);

      //std::cout << "i = " << i << std::endl;

      // Compute ith row of D_l^{-1} r_l^l
      {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
	TEUCHOS_FUNC_TIME_MONITOR("Stokhos: ASC Forward Deterministic Preconditioner Time");
#endif
	mean_prec->ApplyInverse(*(rhs_block->GetBlock(i)), *res_i);
      }

      for (Cijk_type::ik_iterator k_it = Cijk->k_begin(i_it);
	   k_it != Cijk->k_end(i_it); ++k_it) {
	int k = index(k_it);
	if (k < k_limit) {

	  // Compute J[k]*res_i
	  (*sg_poly)[k].Apply(*res_i, *mat_vec_tmp);

	  for (Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	       j_it != Cijk->j_end(k_it); ++j_it) {
	    int j = index(j_it);
	    if (j < n_begin) {
	      double c = value(j_it);
	      if (scale_op) {
		if (useTranspose)
		  c /= norms[i];
		else
		  c /= norms[j];
	      }
	      rhs_block->GetBlock(j)->Update(-c, *mat_vec_tmp, 1.0);
	    }
	  }
	}
      }
      
    }
  }

  // Solve A_0 u_0 = r_0
  {
    // Apply deterministic preconditioner
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
    TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total ASC Deterministic Preconditioner Time");
#endif
    mean_prec->ApplyInverse(*(rhs_block->GetBlock(0)), 
			    *(result_block.GetBlock(0)));
  }

  for (int l=1; l<=P; l++) {

    // Overwrite r_l^{l-1} = D_l^{-1} r_l^l - D_l^{-1}*C_l*u_l^{l-1})
    // Note the first diagonal solve has been done above and is stored in 
    // result_block
    // Here i refers to the row of C_l/D_l
    compute_index(l, n_begin, n_end);
    Cijk_type::i_iterator i_begin = Cijk->find_i(n_begin);
    Cijk_type::i_iterator i_end = Cijk->find_i(n_end);
    for (Cijk_type::i_iterator i_it = i_begin; i_it != i_end; ++i_it) {
      int i = index(i_it);
      Teuchos::RCP<Epetra_MultiVector> res_i = result_block.GetBlock(i);

      // compute row i of C_l*u_l^{l-1}
      for (Cijk_type::ik_iterator k_it = Cijk->k_begin(i_it);
	   k_it != Cijk->k_end(i_it); ++k_it) {
	int k = index(k_it);
	if (k < k_limit) {

	  // Compute J[k]*res_i
	  (*sg_poly)[k].Apply(*res_i, *mat_vec_tmp);

	  for (Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	       j_it != Cijk->j_end(k_it); ++j_it) {
	    int j = index(j_it);
	    if (j < n_begin) {
	      double c = value(j_it);
	      if (scale_op) {
		if (useTranspose)
		  c /= norms[i];
		else
		  c /= norms[j];
	      }
	      rhs_block->GetBlock(j)->Update(-c, *mat_vec_tmp, 1.0);
	    }
	  }
	}
      }

      // Compute row i of D_l^{-1}*C_l*u_l^{l-1}
      {
	// Apply deterministic preconditioner
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
	TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total ASC Deterministic Preconditioner Time");
#endif
	mean_prec->ApplyInverse(*(rhs_block->GetBlock(i)), *mat_vec_tmp);
      }

      // Update u_l^l = u_l^l - D_l^{-1}*C_l*u_l^{l-1}
      res_i->Update(-1.0, *mat_vec_tmp, 1.0);
    }
  }

  if (made_copy)
    delete input;

  return 0; 
}

*/

int 
Stokhos::ApproxSchurComplementPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Approximate Gauss-Seidel Time");
#endif

  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (Input.Values() == Result.Values()) {
    input = new Epetra_MultiVector(Input);
    made_copy = true;
  } 

  int m = input->NumVectors();
  if (rhs_block == Teuchos::null || rhs_block->NumVectors() != m)
    rhs_block = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
  
  // Extract blocks
  EpetraExt::BlockMultiVector input_block(View, *base_map, *input);
  EpetraExt::BlockMultiVector result_block(View, *base_map, Result);

  result_block.PutScalar(0.0);

  // Set right-hand-side to input_block
  rhs_block->Update(1.0, input_block, 0.0);

  // At level l, linear system has the structure
  // [ A_{l-1} B_l ][ u_l^{l-1} ] = [ r_l^{l-1} ]
  // [ C_l     D_l ][ u_l^l     ]   [ r_l^l     ]

  int P = sg_basis->order();
  int n_begin, n_end;
  for (int l=P; l>=1; l--) {
    compute_index(l, n_begin, n_end);

    // Compute D_l^{-1} r_l^l
    divide_diagonal_block(n_begin, n_end, *rhs_block, result_block);

    // Compute r_l^{l-1} = r_l^{l-1} - B_l D_l^{-1} r_l^l
    multiply_block(0, n_begin, n_begin, n_end, -1.0, result_block, *rhs_block);
  }

  // Solve A_0 u_0 = r_0
  compute_index(0, n_begin, n_end);
  divide_diagonal_block(n_begin, n_end, *rhs_block, result_block);

  for (int l=1; l<=P; l++) {
    compute_index(l, n_begin, n_end);

    // Compute r_l^l - C_l*u_l^{l-1}
    multiply_block(n_begin, n_end, 0, n_begin, -1.0, result_block, *rhs_block);

    // Compute D_l^{-1} (r_l^l - C_l*u_l^{l-1})
    divide_diagonal_block(n_begin, n_end, *rhs_block, result_block);
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
compute_index(int l, int& n_begin, int& n_end) const
{
  Teuchos::RCP<const Stokhos::ProductBasis<int,double> > prod_basis =
    Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int,double> >(sg_basis, true);
  int d = prod_basis->dimension();
  Teuchos::Array<int> term(d);
  term[0] = l;
  n_begin = prod_basis->getIndex(term);
  term[0] = 0;
  term[d-1] = l;
  n_end = prod_basis->getIndex(term) + 1;
}

void
Stokhos::ApproxSchurComplementPreconditioner::
multiply_block(int row_begin, int row_end, int col_begin, int col_end,
	       double alpha,
	       const EpetraExt::BlockMultiVector& Input, 
	       EpetraExt::BlockMultiVector& Result) const
{
  // Input and Result are the whole vector/multi-vector, not just the portion
  // needed for the particular sub-block

  // Allocate temporary storage
  int m = Input.NumVectors();
  if (tmp == Teuchos::null || tmp->NumVectors() != m*max_num_mat_vec)
    tmp = Teuchos::rcp(new Epetra_MultiVector(*base_map, 
					      m*max_num_mat_vec));

  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
    int k = index(k_it);
    Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
    Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
    int nj = 0;
    for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
      int j = index(j_it);
      if (j >= col_begin && j < col_end) 
	++nj;
    }
    if (nj > 0) {
      Teuchos::Array<double*> j_ptr(nj*m);
      Teuchos::Array<int> mj_indices(nj*m);
      int l = 0;
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
        int j = index(j_it);
	if (j >= col_begin && j < col_end) {
	  for (int mm=0; mm<m; mm++) {
	    j_ptr[l*m+mm] = (*(Input.GetBlock(j)))[mm];
	    mj_indices[l*m+mm] = l*m+mm;
	  }
	  l++;
	}
      }
      Epetra_MultiVector input_tmp(View, *base_map, &j_ptr[0], nj*m);
      Epetra_MultiVector result_tmp(View, *tmp, &mj_indices[0], nj*m);
      (*sg_poly)[k].Apply(input_tmp, result_tmp);
      l = 0;
      for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	int j = index(j_it);
	if (j >= col_begin && j < col_end) {
	  Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
	  Cijk_type::kji_iterator i_end = Cijk->i_end(j_it);
	  for (Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
	    int i = index(i_it);
	    if (i >= row_begin && i < row_end) {
	      double c = value(i_it);
	      if (scale_op)
		c /= norms[i];
	      for (int mm=0; mm<m; mm++)
		(*Result.GetBlock(i))(mm)->Update(alpha*c, *result_tmp(l*m+mm), 1.0);
	    }
	  }
	  l++;
	}
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
