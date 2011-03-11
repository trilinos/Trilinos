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
