// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_ConfigDefs.h"
#include "Stokhos_ParallelData.hpp"
#ifdef HAVE_MPI
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif

Stokhos::ParallelData::
ParallelData(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk,
  const Teuchos::RCP<const Epetra_Comm>& globalComm,
  Teuchos::ParameterList& params)
{
  int num_global_stoch_blocks = basis->size();

  int num_spatial_procs = params.get("Number of Spatial Processors", -1);

  // Build multi-comm
  globalMultiComm = 
    Stokhos::buildMultiComm(*globalComm, num_global_stoch_blocks,
			    num_spatial_procs);

  // Get stochastic and spatial comm's
  stoch_comm = Stokhos::getStochasticComm(globalMultiComm);
  spatial_comm = Stokhos::getSpatialComm(globalMultiComm);

  if (Cijk != Teuchos::null) {
    // Build Epetra Cijk
    epetraCijk = 
      Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(basis, Cijk, 
						    globalMultiComm));
    
    // Rebalance graphs
    bool use_isorropia = params.get("Rebalance Stochastic Graph", false);
    if (use_isorropia)
    epetraCijk->rebalance(params.sublist("Isorropia"));
    
    // Transform to local indices
    epetraCijk->transformToLocal();
  }
}

Stokhos::ParallelData::
ParallelData(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk,
  const Teuchos::RCP<const EpetraExt::MultiComm>& globalMultiComm_,
  Teuchos::ParameterList& params) :
  globalMultiComm(globalMultiComm_)
{
  // Get stochastic and spatial comm's
  stoch_comm = Stokhos::getStochasticComm(globalMultiComm);
  spatial_comm = Stokhos::getSpatialComm(globalMultiComm);

  if (Cijk != Teuchos::null) {
    // Build Epetra Cijk
    epetraCijk = 
      Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(basis, Cijk, 
						    globalMultiComm));

    // Rebalance graphs
    bool use_isorropia = params.get("Rebalance Stochastic Graph", false);
    if (use_isorropia)
      epetraCijk->rebalance(params.sublist("Isorropia"));
    
    // Transform to local indices
    epetraCijk->transformToLocal();
  }
}
 
Teuchos::RCP<const EpetraExt::MultiComm> 
Stokhos::buildMultiComm(const Epetra_Comm& globalComm,
			int num_global_stochastic_blocks,
			int num_spatial_procs)
{
  Teuchos::RCP<const EpetraExt::MultiComm> globalMultiComm;

#ifdef HAVE_MPI
  if (num_spatial_procs == -1) {
    // By default, use all procs for spatial parallelism
    //MPI_Comm_size(MPI_COMM_WORLD, &num_spatial_procs);
    num_spatial_procs = globalComm.NumProc();
  }
  const Epetra_MpiComm& globalMpiComm = 
    dynamic_cast<const Epetra_MpiComm&>(globalComm);
  globalMultiComm = 
    Teuchos::rcp(new EpetraExt::MultiMpiComm(globalMpiComm.Comm(), 
					     num_spatial_procs, 
					     num_global_stochastic_blocks,
					     Teuchos::VERB_NONE));
#else
  globalMultiComm = 
    Teuchos::rcp(new EpetraExt::MultiSerialComm(num_global_stochastic_blocks));
#endif

  return globalMultiComm;
}

Teuchos::RCP<const Epetra_Comm> 
Stokhos::getSpatialComm(
  const Teuchos::RCP<const EpetraExt::MultiComm>& globalMultiComm)
{
  return Teuchos::rcp(&(globalMultiComm->SubDomainComm()), false);
}

Teuchos::RCP<const Epetra_Comm> 
Stokhos::getStochasticComm(
  const Teuchos::RCP<const EpetraExt::MultiComm>& globalMultiComm)
{
  return Teuchos::rcp(&(globalMultiComm->TimeDomainComm()), false);
}
