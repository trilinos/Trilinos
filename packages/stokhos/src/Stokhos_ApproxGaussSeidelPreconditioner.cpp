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

#include "Stokhos_ApproxGaussSeidelPreconditioner.hpp"
#include "Epetra_config.h"
#include "Teuchos_TimeMonitor.hpp"

Stokhos::ApproxGaussSeidelPreconditioner::
ApproxGaussSeidelPreconditioner(
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<Stokhos::PreconditionerFactory>& prec_factory_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  label("Stokhos Approximate Gauss-Seidel Preconditioner"),
  base_map(base_map_),
  sg_map(sg_map_),
  prec_factory(prec_factory_),
  mean_prec(),
  useTranspose(false),
  sg_op(),
  sg_poly(),
  Cijk(),
  symmetric(false),
  only_use_linear(true),
  mat_vec_tmp(),
  rhs_block()
{
  symmetric = params_->get("Symmetric Gauss-Seidel", false);
  only_use_linear = params_->get("Only Use Linear Terms", true);
}

Stokhos::ApproxGaussSeidelPreconditioner::
~ApproxGaussSeidelPreconditioner()
{
}

void
Stokhos::ApproxGaussSeidelPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op_, 
		    const Epetra_Vector& x)
{
  sg_op = sg_op_;
  sg_poly = sg_op->getSGPolynomial();
  mean_prec = prec_factory->compute(sg_poly->getCoeffPtr(0));
  label = std::string("Stokhos Approximate Gauss-Seidel Preconditioner:\n") + 
    std::string("		***** ") + 
    std::string(mean_prec->Label());
  Cijk = sg_op->getTripleProduct();
}

int 
Stokhos::ApproxGaussSeidelPreconditioner::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;

  return 0;
}

int 
Stokhos::ApproxGaussSeidelPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  return sg_op->Apply(Input, Result);
}

int 
Stokhos::ApproxGaussSeidelPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  TEUCHOS_FUNC_TIME_MONITOR("Total Approximate Gauss-Seidel Time");

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

  int sz = sg_poly->basis()->size();
  int i_limit = sg_poly->size();
  if (only_use_linear)
    i_limit = sg_poly->basis()->dimension() + 1;
  const Teuchos::Array<double>& norms = sg_poly->basis()->norm_squared();
  int i,j,nl;
  double c;

  rhs_block->Update(1.0, input_block, 0.0);
  for (int k=0; k<sz; k++) {
    nl = Cijk->num_values(k);
    for (int l=0; l<nl; l++) {
      Cijk->value(k,l,i,j,c); 
      if (i!=0 && i<i_limit) {
	(*sg_poly)[i].Apply(*(result_block.GetBlock(j)), *mat_vec_tmp);
	rhs_block->GetBlock(k)->Update(-1.0*c/norms[k], *mat_vec_tmp, 1.0);
      }
    }
      
    {
      // Apply deterministic preconditioner
      TEUCHOS_FUNC_TIME_MONITOR("Total AGS Deterministic Preconditioner Time");
      mean_prec->ApplyInverse(*(rhs_block->GetBlock(k)), 
			      *(result_block.GetBlock(k)));
    }
  }

  // For symmetric Gauss-Seidel
  if (symmetric) {
    
    rhs_block->Update(1.0, input_block, 0.0);
    for (int k=sz-1; k>=0; k--) {
      nl = Cijk->num_values(k);
      for (int l=0; l<nl; l++) {
	Cijk->value(k,l,i,j,c); 
	if (i!=0 && i<i_limit) {
	  (*sg_poly)[i].Apply(*(result_block.GetBlock(j)), *mat_vec_tmp);
	  rhs_block->GetBlock(k)->Update(-1.0*c/norms[k], *mat_vec_tmp, 1.0);
	}
      }
      
      {
	// Apply deterministic preconditioner
	TEUCHOS_FUNC_TIME_MONITOR("Total AGS Deterministic Preconditioner Time");
	mean_prec->ApplyInverse(*(rhs_block->GetBlock(k)), 
				*(result_block.GetBlock(k)));
      }
    
    }
  }

  if (made_copy)
    delete input;

  return 0; 
}

double 
Stokhos::ApproxGaussSeidelPreconditioner::
NormInf() const
{
  return sg_op->NormInf();
}


const char* 
Stokhos::ApproxGaussSeidelPreconditioner::
Label() const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::ApproxGaussSeidelPreconditioner::
UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::ApproxGaussSeidelPreconditioner::
HasNormInf() const
{
  return sg_op->HasNormInf();
}

const Epetra_Comm & 
Stokhos::ApproxGaussSeidelPreconditioner::
Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
Stokhos::ApproxGaussSeidelPreconditioner::
OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::ApproxGaussSeidelPreconditioner::
OperatorRangeMap() const
{
  return *sg_map;
}
