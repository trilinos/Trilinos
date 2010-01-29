#ifndef PIRO_EPETRA_NOXOBSERVER
#define PIRO_EPETRA_NOXOBSERVER

#include "Teuchos_RCP.hpp"
#include "EpetraExt_BlockVector.h"
#include "Epetra_Map.h"


namespace Piro {
namespace Epetra {

class NOXObserver 
{
public:
  NOXObserver () {};
  virtual ~NOXObserver () {};

  virtual void observeSolution(const Epetra_Vector& soln) = 0;

private:

};

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

#endif //PIRO_EPETRA_NOXOBSERVER
