#include "Piro_Epetra_StokhosNOXObserver.hpp"

Piro::Epetra::StokhosNOXObserver::StokhosNOXObserver (
     Teuchos::RCP<Piro::Epetra::NOXObserver> noxObserver_,
     const Epetra_Map& map_,
     const int sz_) :
  noxObserver(noxObserver_),
  map(map_),
  numSGBlocks(sz_)
{
 //if (noxObserver == Teuchos::null) cout << "XXX1" << endl;
}

void Piro::Epetra::StokhosNOXObserver::observeSolution(
    const Epetra_Vector& solution)
{

cout << solution << "\nendsolution" << endl;
  // Copy into block vector, so Block access is available
  EpetraExt::BlockVector blockVec(View, map, solution);
cout << blockVec << "\nendblockvec" << endl;

  Teuchos::RCP<Epetra_Vector> deterministicSizedSolnVec;

  for (int i=0; i< numSGBlocks; i++) {
    deterministicSizedSolnVec = blockVec.GetBlock(i);

cout << *deterministicSizedSolnVec << "\nenddeterministicSizedSolnVec  " << i << endl;

    if (noxObserver != Teuchos::null)
      noxObserver->observeSolution(*deterministicSizedSolnVec);
  }

}
