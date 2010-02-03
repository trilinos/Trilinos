#ifndef PIRO_EPETRA_STOKHOSNOXOBSERVER
#define PIRO_EPETRA_STOKHOSNOXOBSERVER

#include "Piro_Epetra_NOXObserver.hpp"
#include "EpetraExt_BlockVector.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"

namespace Piro {
namespace Epetra {

class StokhosNOXObserver : public NOXObserver
{
public:
  StokhosNOXObserver (
         Teuchos::RCP<NOXObserver> noxObserver_,
         const Epetra_Map& map_,
         const int sz_);

  void observeSolution(const Epetra_Vector& soln);

private:

   Teuchos::RCP<NOXObserver> noxObserver;
   Epetra_Map map;
   const int numSGBlocks;
};

}
}

#endif //PIRO_EPETRA_STOKHOSNOXOBSERVER
