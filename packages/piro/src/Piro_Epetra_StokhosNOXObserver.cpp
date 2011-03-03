// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
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
  save_moments(save_moments_),
  numSGBlocks(basis->size())
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
      noxObserver->observeSolution(vec_poly[i]);
    }
  }
  else {
    // Always write out first moment
    vec_poly.computeMean(*moment);
    noxObserver->observeSolution(*moment);
    if (save_moments >= 2) {
      vec_poly.computeStandardDeviation(*moment);
      noxObserver->observeSolution(*moment);
    }
  }
  

}
