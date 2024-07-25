// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_GaussSeidelPreconditioner.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "EpetraExt_BlockUtility.h"
#include "Teuchos_Assert.hpp"

Stokhos::GaussSeidelPreconditioner::
GaussSeidelPreconditioner(
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk_,
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& det_solver_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_):
  label("Stokhos Gauss-Seidel Preconditioner"),
  sg_comm(sg_comm_),
  sg_basis(sg_basis_),
  epetraCijk(epetraCijk_),
  base_map(base_map_),
  sg_map(sg_map_),
  is_stoch_parallel(epetraCijk->isStochasticParallel()),
  stoch_row_map(epetraCijk->getStochasticRowMap()),
  det_solver(det_solver_),
  params(params_),
  useTranspose(false),
  sg_op(),
  sg_poly(),
  Cijk(epetraCijk->getParallelCijk()),
  is_parallel(epetraCijk->isStochasticParallel())
{
  if (is_parallel) {
    Teuchos::RCP<const Epetra_BlockMap> stoch_col_map =
      epetraCijk->getStochasticColMap();
    sg_col_map =
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*base_map,
							     *stoch_col_map,
							     sg_map->Comm()));
    col_exporter = Teuchos::rcp(new Epetra_Export(*sg_col_map, *sg_map));
  }
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
}

int 
Stokhos::GaussSeidelPreconditioner::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  TEUCHOS_TEST_FOR_EXCEPTION(
      UseTranspose == true, std::logic_error,
      "Stokhos::GaussSeidelPreconditioner::SetUseTranspose():  " <<
      "Preconditioner does not support transpose!" << std::endl);

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
  bool scale_op = params->get("Scale Operator by Inverse Basis Norms", true);
  bool only_use_linear = params->get("Only Use Linear Terms", true);
  
  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (Input.Values() == Result.Values()) {
    input = new Epetra_MultiVector(Input);
    made_copy = true;
  }

  int k_limit = sg_poly->size();
  int dim = sg_poly->basis()->dimension();
  if (only_use_linear && sg_poly->size() > dim+1)
    k_limit = dim + 1;
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  
  int m = input->NumVectors();
  if (sg_df_block == Teuchos::null || sg_df_block->NumVectors() != m) {
    sg_df_block = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
    sg_y_block = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
    kx = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));
    if (is_parallel) {
      sg_df_col = Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, 
							       *sg_col_map, m));
      sg_df_tmp = Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, 
							       *sg_map, m));
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
  
  Teuchos::RCP<Epetra_MultiVector> f, df, dx;

  sg_df_block->Update(1.0, sg_f_block, 0.0);

  int iter = 0;
  while (((norm_df/norm_f)>sg_tol) && (iter<max_iter)) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
    TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total global solve Time");
#endif
    iter++;

    sg_y_block->Update(1.0, sg_f_block, 0.0);
    if (is_parallel)
      sg_df_col->PutScalar(0.0);
    
    for (Cijk_type::i_iterator i_it=Cijk->i_begin(); 
	 i_it!=Cijk->i_end(); ++i_it) {
      int i = index(i_it);

      f = sg_f_block.GetBlock(i);
      df = sg_df_block->GetBlock(i);
      dx = sg_dx_block.GetBlock(i);

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
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
	  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total deterministic solve Time");
#endif
	  det_solver->applyJacobianInverse(det_solver_params, nox_df, nox_dx);
	}
      }

      df->Update(1.0, *f, 0.0);
      
      for (Cijk_type::ik_iterator k_it = Cijk->k_begin(i_it);
	   k_it != Cijk->k_end(i_it); ++k_it) {
	int k = index(k_it);
	if (k!=0 && k<k_limit) {
	  (*sg_poly)[k].Apply(*dx, *kx);
	  for (Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	       j_it != Cijk->j_end(k_it); ++j_it) {
	    int j = index(j_it);
	    int j_gid = epetraCijk->GCID(j);
	    double c = value(j_it);
	    if (scale_op)
	      c /= norms[j_gid];
	    bool owned = epetraCijk->myGRID(j_gid);
	    if (!is_parallel || owned) {
	      sg_df_block->GetBlock(j)->Update(-c, *kx, 1.0);
	      sg_y_block->GetBlock(j)->Update(-c, *kx, 1.0);
	    }
	    else
	      sg_df_col->GetBlock(j)->Update(-c, *kx, 1.0);
	  }
	}
      }
      
      (*sg_poly)[0].Apply(*dx, *kx);
      sg_y_block->GetBlock(i)->Update(-1.0, *kx, 1.0);
      
    } //End of k loop

    // Add in contributions from off-process
    if (is_parallel) {
      sg_df_tmp->Export(*sg_df_col, *col_exporter, InsertAdd);
      sg_df_block->Update(1.0, *sg_df_tmp, 1.0);
      sg_y_block->Update(1.0, *sg_df_tmp, 1.0);
    }
    
    sg_y_block->Norm2(&norm_df);
    //std::cout << "norm_df = " << norm_df << std::endl;
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
  return *sg_comm;
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
