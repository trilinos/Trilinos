// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_STOKHOSNOXOBSERVER
#define PIRO_EPETRA_STOKHOSNOXOBSERVER

#include "NOX_Epetra_Observer.H"
#include "EpetraExt_BlockVector.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "EpetraExt_MultiComm.h"
#include "Epetra_Import.h"

namespace Piro {
namespace Epetra {

class StokhosNOXObserver : public NOX::Epetra::Observer
{
public:
  StokhosNOXObserver (
    const Teuchos::RCP<NOX::Epetra::Observer>& noxObserver_,
    const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis_,
    const Teuchos::RCP<const Epetra_BlockMap>& stoch_map_,
    const Teuchos::RCP<const Epetra_BlockMap>& spatial_map_,
    const Teuchos::RCP<const Epetra_BlockMap>& product_map_,
    const Teuchos::RCP<const EpetraExt::MultiComm>& product_comm_,
    const Teuchos::RCP<const Epetra_Import>& importer_,
    int save_moments_ = -1);

  void observeSolution(const Epetra_Vector& soln);

private:

   Teuchos::RCP<NOX::Epetra::Observer> noxObserver;
   Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;
   Teuchos::RCP<const Epetra_BlockMap> stoch_map;
   Teuchos::RCP<const Epetra_BlockMap> spatial_map;
   Teuchos::RCP<const Epetra_BlockMap> product_map;
   Teuchos::RCP<const EpetraExt::MultiComm> product_comm;
   Teuchos::RCP<const Epetra_Import> importer;
   const int numSGBlocks;
   int save_moments;
   Teuchos::RCP<Epetra_Vector> moment;
   Teuchos::RCP<Epetra_Vector> overlap_vec;
};

}
}

#endif //PIRO_EPETRA_STOKHOSNOXOBSERVER
