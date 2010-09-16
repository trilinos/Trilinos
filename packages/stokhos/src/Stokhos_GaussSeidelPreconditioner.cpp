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

#include "Stokhos_GaussSeidelPreconditioner.hpp"
#include "Epetra_config.h"
#include "Teuchos_TimeMonitor.hpp"

Stokhos::GaussSeidelPreconditioner::
GaussSeidelPreconditioner(
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& det_solver_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_):
  label("Stokhos Gauss-Seidel Preconditioner"),
  base_map(base_map_),
  sg_map(sg_map_),
  det_solver(det_solver_),
  params(params_),
  useTranspose(false),
  sg_op(),
  sg_poly(),
  Cijk()
{
}

Stokhos::GaussSeidelPreconditioner::
~GaussSeidelPreconditioner()
{
}

void
Stokhos::GaussSeidelPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op_, 
		    const Epetra_Vector& x)
{
  sg_op = sg_op_;
  sg_poly = sg_op->getSGPolynomial();
  EpetraExt::BlockVector sg_x_block(View, *base_map, x);
  det_solver->setJacobianOperatorForSolve(sg_poly->getCoeffPtr(0));
  det_solver->createPreconditioner(*(sg_x_block.GetBlock(0)),
				   params->sublist("Deterministic Solver Parameters"),
				   false);
  Cijk = sg_op->getTripleProduct();
}

int 
Stokhos::GaussSeidelPreconditioner::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;

  return 0;
}

int 
Stokhos::GaussSeidelPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  return sg_op->Apply(Input,Result);
}

int 
Stokhos::GaussSeidelPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  int max_iter = params->get("Max Iterations",100 );
  double sg_tol = params->get("Tolerance", 1e-12);
  bool MatVecTable = params->get("Save MatVec Table", true);
  bool only_use_linear = params->get("Only Use Linear Terms", true);
  
  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (Input.Values() == Result.Values()) {
    input = new Epetra_MultiVector(Input);
    made_copy = true;
  }

  int sz = sg_poly->basis()->size();
  int K_limit = sg_poly->size();
  int dim = sg_poly->basis()->dimension();
  if (only_use_linear && sg_poly->size() > dim+1)
    K_limit = dim + 1;
  const Teuchos::Array<double>& norms = sg_poly->basis()->norm_squared();

  int m = input->NumVectors();
  if (sg_df_block == Teuchos::null || sg_df_block->NumVectors() != m) {
    sg_df_block = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
    kx = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));
    if (MatVecTable) {
      Kx_table.resize(sz);
      for (int i=0;i<sz;i++) {
	Kx_table[i].resize(K_limit);
	for (int j=0;j<K_limit;j++) {
	  Kx_table[i][j] = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));
	}
      }
    }
  }
  else {
    if (MatVecTable) {
      for (int i=0;i<sz;i++) {
  	for (int j=0;j<K_limit;j++) {
  	  Kx_table[i][j]->PutScalar(0.0);
  	}
      }
    }
  }
  
  // Extract blocks
  EpetraExt::BlockMultiVector sg_dx_block(View, *base_map, Result);
  EpetraExt::BlockMultiVector sg_f_block(View, *base_map, *input);
  sg_dx_block.PutScalar(0.0);
   
  // Compute initial residual norm
  double norm_f,norm_df;
  sg_f_block.Norm2(&norm_f);
  sg_op->Apply(sg_dx_block, *sg_df_block);
  sg_df_block->Update(-1.0, sg_f_block, 1.0);
  sg_df_block->Norm2(&norm_df);
  
  Teuchos::RCP<Epetra_MultiVector> df, dx;

  int iter = 0;
  while (((norm_df/norm_f)>sg_tol) && (iter<max_iter)) {
    TEUCHOS_FUNC_TIME_MONITOR("Total global solve Time");
    iter++;
    
    // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
    //  ordinal_type Cijk_size = Cijk.size();
    for (int k=0; k<sz; k++) {
      df = sg_df_block->GetBlock(k);
      dx = sg_dx_block.GetBlock(k);
      df->Update(1.0, *sg_f_block.GetBlock(k), 0.0);

      int nl = Cijk->num_values(k);
      for (int l=0; l<nl; l++) {
	int i,j;
	double c;
	Cijk->value(k,l,i,j,c); 
	if (i!=0 && i<K_limit) {
	  if (MatVecTable) {
	    df->Update(-1.0*c/norms[k],*(Kx_table[j][i]),1.0);      
	  }
	  else {
	    (*sg_poly)[i].Apply(*(sg_dx_block.GetBlock(j)),*kx);
	    df->Update(-1.0*c/norms[k],*kx,1.0);      
	  }  
	}
      }
      
      (*sg_poly)[0].Apply(*dx, *kx);

      dx->PutScalar(0.0);
      Teuchos::ParameterList& det_solver_params = 
	params->sublist("Deterministic Solver Parameters");
      for (int col=0; col<m; col++) {
	NOX::Epetra::Vector nox_df(Teuchos::rcp((*df)(col),false), 
				   NOX::Epetra::Vector::CreateView);
	NOX::Epetra::Vector nox_dx(Teuchos::rcp((*dx)(col),false), 
				   NOX::Epetra::Vector::CreateView);
	// Solve linear system
	{
	  TEUCHOS_FUNC_TIME_MONITOR("Total deterministic solve Time");
	  det_solver->applyJacobianInverse(det_solver_params, nox_df, nox_dx);
	}
      }

      df->Update(-1.0, *kx, 1.0);
      
      if (MatVecTable) {
	for(int i=1;i<K_limit;i++) {
	  (*sg_poly)[i].Apply(*dx, *(Kx_table[k][i]));
	}
      }
    } //End of k loop
    
    sg_df_block->Norm2(&norm_df);
  } //End of iter loop

  if (made_copy)
    delete input;

  return 0; 
}

double 
Stokhos::GaussSeidelPreconditioner::
NormInf() const
{
  return sg_op->NormInf();
}


const char* 
Stokhos::GaussSeidelPreconditioner::
Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::GaussSeidelPreconditioner::
UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::GaussSeidelPreconditioner::
HasNormInf() const
{
  return sg_op->HasNormInf();
}

const Epetra_Comm & 
Stokhos::GaussSeidelPreconditioner::
Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
Stokhos::GaussSeidelPreconditioner::
OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::GaussSeidelPreconditioner::
OperatorRangeMap() const
{
  return *sg_map;
}
