// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_ConfigDefs.h"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Stokhos_Sparse3TensorUtilities.hpp"
#include "Teuchos_Assert.hpp"
#ifdef HAVE_STOKHOS_ISORROPIA
#include "Isorropia_Epetra.hpp"
#endif
#include "Epetra_Map.h"

Stokhos::EpetraSparse3Tensor::
EpetraSparse3Tensor(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& globalMultiComm_,
  int k_begin_, int k_end_) :
  basis(basis_),
  Cijk(Cijk_),
  globalMultiComm(globalMultiComm_),
  k_begin(k_begin_),
  k_end(k_end_)
{
  // Get stochastic distribution
  num_global_stoch_blocks = basis->size();
  stoch_comm = Teuchos::rcp(&(globalMultiComm->TimeDomainComm()), false);
  is_parallel = (stoch_comm->NumProc() > 1);
  if (k_end == -1)
    k_end = Cijk->num_k();

  // Build stochastic row map -- stochastic blocks evenly distributed
  // across stochastic procs
  stoch_row_map = Teuchos::rcp(new Epetra_Map(num_global_stoch_blocks, 0,
					      *stoch_comm));
  
  // Build Cijk tensor parallel over i
  if (!is_parallel)
    Cijk_parallel = Cijk;
  else
    Cijk_parallel = buildParallelCijk();

  // Build stochastic graph
  stoch_graph = Stokhos::sparse3Tensor2CrsGraph(*Cijk_parallel, *stoch_row_map);
  
  // Build stochastic column map
  stoch_col_map = Teuchos::rcp(&(stoch_graph->ColMap()), false);
}

Stokhos::EpetraSparse3Tensor::
EpetraSparse3Tensor(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& globalMultiComm_,
  const Teuchos::RCP<const Epetra_BlockMap>& stoch_row_map_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_parallel_,
  int k_begin_, int k_end_) :
  basis(basis_),
  Cijk(Cijk_),
  globalMultiComm(globalMultiComm_),
  k_begin(k_begin_),
  k_end(k_end_),
  stoch_row_map(stoch_row_map_),
  Cijk_parallel(Cijk_parallel_)
{
  // Get stochastic distribution
  num_global_stoch_blocks = basis->size();
  stoch_comm = Teuchos::rcp(&(globalMultiComm->TimeDomainComm()), false);
  is_parallel = (stoch_comm->NumProc() > 1);
  if (k_end == -1)
    k_end = Cijk->num_k();

  // Build Cijk tensor parallel over i
  if (Cijk_parallel == Teuchos::null) {
    if (!is_parallel)
      Cijk_parallel = Cijk;
    else
      Cijk_parallel = buildParallelCijk();
  }

  // Build stochastic graph
  stoch_graph = Stokhos::sparse3Tensor2CrsGraph(*Cijk_parallel, *stoch_row_map);
  
  // Build stochastic column map
  stoch_col_map = Teuchos::rcp(&(stoch_graph->ColMap()), false);
}

Stokhos::EpetraSparse3Tensor::
EpetraSparse3Tensor(const EpetraSparse3Tensor& epetraCijk,
		    int k_begin_, int k_end_) :
  basis(epetraCijk.basis),
  Cijk(epetraCijk.Cijk),
  globalMultiComm(epetraCijk.globalMultiComm),
  num_global_stoch_blocks(epetraCijk.num_global_stoch_blocks),
  k_begin(k_begin_),
  k_end(k_end_),
  stoch_comm(epetraCijk.stoch_comm),
  is_parallel(epetraCijk.is_parallel),
  stoch_row_map(epetraCijk.stoch_row_map)
{
  if (k_end == -1)
    k_end = Cijk->num_k();

  // Build Cijk tensor parallel over i
  if (!is_parallel)
    Cijk_parallel = Cijk;
  else
    Cijk_parallel = buildParallelCijk();
  
  // Build stochastic graph
  stoch_graph = Stokhos::sparse3Tensor2CrsGraph(*Cijk_parallel, *stoch_row_map);
  
  // Build stochastic column map
  stoch_col_map = Teuchos::rcp(&(stoch_graph->ColMap()), false);
}

void
Stokhos::EpetraSparse3Tensor::
rebalance(Teuchos::ParameterList& isorropia_params)
{
#ifdef HAVE_STOKHOS_ISORROPIA
  // Reblance with Isorropia
  Teuchos::RCP<const Epetra_CrsGraph> rebalanced_graph = 
    Teuchos::rcp(Isorropia::Epetra::createBalancedCopy(
		   *stoch_graph, isorropia_params));
  stoch_graph = rebalanced_graph;
  
  // Get new maps
  stoch_row_map = Teuchos::rcp(&(stoch_graph->RowMap()), false);
  stoch_col_map = Teuchos::rcp(&(stoch_graph->ColMap()), false);
  
  // Build new parallel Cijk
  Cijk_parallel = buildParallelCijk();
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Error!  Stokhos must be configured with " <<
		     "isorropia support to rebalance the stochastic " <<
		     "graph.");
#endif
}

Teuchos::RCP<Stokhos::EpetraSparse3Tensor::Cijk_type> 
Stokhos::EpetraSparse3Tensor::
buildParallelCijk() const
{
  Teuchos::RCP<Cijk_type> Cijk_tmp = Teuchos::rcp(new Cijk_type);
  int *rowIndices = stoch_row_map->MyGlobalElements();
  int myBlockRows = stoch_row_map->NumMyElements();
  for (int i=0; i<myBlockRows; i++) {
    int myRow = rowIndices[i];
    Cijk_type::i_iterator i_it = Cijk->find_i(myRow);
    for (Cijk_type::ik_iterator k_it = Cijk->k_begin(i_it);
	 k_it != Cijk->k_end(i_it); ++k_it) {
      int k = index(k_it);
      if (k >= k_begin && k < k_end) {
	for (Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	     j_it != Cijk->j_end(k_it); ++j_it) {
	  int j = index(j_it);
	  double c = value(j_it);
	  Cijk_tmp->add_term(myRow, j, k, c);
	}
      }
    }
  }
  Cijk_tmp->fillComplete();

  return Cijk_tmp;
}

void
Stokhos::EpetraSparse3Tensor::
transformToLocal()
{
  Teuchos::RCP<Cijk_type> Cijk_tmp = Teuchos::rcp(new Cijk_type);
  for (Cijk_type::i_iterator i_it = Cijk_parallel->i_begin();
       i_it != Cijk_parallel->i_end(); ++i_it) {
    int i = stoch_row_map->LID(index(i_it));
    TEUCHOS_TEST_FOR_EXCEPTION(i == -1, std::logic_error,
		       "Error!  Global row index " << index(i_it) <<
		       " is not on processor " << globalMultiComm->MyPID());
    for (Cijk_type::ik_iterator k_it = Cijk_parallel->k_begin(i_it);
	 k_it != Cijk_parallel->k_end(i_it); ++k_it) {
      int k = index(k_it);
      for (Cijk_type::ikj_iterator j_it = Cijk_parallel->j_begin(k_it);
	   j_it != Cijk_parallel->j_end(k_it); ++j_it) {
	int j= stoch_col_map->LID(index(j_it));
	TEUCHOS_TEST_FOR_EXCEPTION(j == -1, std::logic_error,
			   "Error!  Global col index " << index(j_it) <<
			   " is not on processor " << globalMultiComm->MyPID());
	double c = value(j_it);
	Cijk_tmp->add_term(i, j, k, c);
      }
    }
  }
  Cijk_tmp->fillComplete();
  
  Cijk_parallel = Cijk_tmp;
}
