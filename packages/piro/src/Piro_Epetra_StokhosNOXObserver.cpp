// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_Epetra_StokhosNOXObserver.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"

Piro::Epetra::StokhosNOXObserver::StokhosNOXObserver (
     const Teuchos::RCP<NOX::Epetra::Observer>& noxObserver_,
     const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis_,        const Teuchos::RCP<const Epetra_BlockMap>& stoch_map_,
     const Teuchos::RCP<const Epetra_BlockMap>& spatial_map_,
     const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
     const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
     const Teuchos::RCP<const Epetra_Import>& importer_,
     int save_moments_) :
  noxObserver(noxObserver_),
  basis(basis_),  
  stoch_map(stoch_map_),
  spatial_map(spatial_map_),
  product_map(product_map_),
  product_comm(product_comm_),
  importer(importer_),
  numSGBlocks(basis->size()),
  save_moments(save_moments_)
{
 //if (noxObserver == Teuchos::null) cout << "XXX1" << endl;
  if (save_moments > 0)
    moment = Teuchos::rcp(new Epetra_Vector(*spatial_map));
  if (product_map != Teuchos::null)
    overlap_vec = Teuchos::rcp(new Epetra_Vector(*product_map));
}

void Piro::Epetra::StokhosNOXObserver::observeSolution(
    const Epetra_Vector& solution)
{

  if (noxObserver == Teuchos::null)
    return;

  // Copy into block vector, so Block access is available
  overlap_vec->Import(solution, *importer, Insert);
  Stokhos::EpetraVectorOrthogPoly vec_poly(
    basis, stoch_map, spatial_map, product_map, product_comm, View, 
    *overlap_vec);
  if (save_moments <= 0) {
    for (int i=0; i< numSGBlocks; i++) {
      noxObserver->observeSolution(vec_poly[i], i);
    }
  }
  else {
    // Always write out first moment
    vec_poly.computeMean(*moment);
    noxObserver->observeSolution(*moment, 1);
    if (save_moments >= 2) {
      vec_poly.computeStandardDeviation(*moment);
      noxObserver->observeSolution(*moment, 2);
    }
  }
  

}
