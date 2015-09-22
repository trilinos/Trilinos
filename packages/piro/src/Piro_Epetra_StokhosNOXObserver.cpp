// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
