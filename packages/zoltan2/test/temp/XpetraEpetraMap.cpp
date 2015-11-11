#include <Zoltan2_Standards.hpp>

#ifdef HAVE_ZOLTAN2_EPETRA

#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>

#include <Tpetra_Map.hpp>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_EpetraVector.hpp>
#include <Xpetra_EpetraUtils.hpp>

int main(int narg, char **arg)
{
  Teuchos::GlobalMPISession session(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > tcomm = 
               Teuchos::DefaultComm<int>::getComm();
  Teuchos::RCP<const Epetra_Comm> ecomm = Xpetra::toEpetra(tcomm);

  const int nGlobRows = 50;
  Epetra_Map emap(nGlobRows, 0, *ecomm);

  const int nRows = emap.NumMyElements();
  int *myGlobRows = emap.MyGlobalElements();

  Epetra_CrsMatrix A(Copy, emap, 3);
  double tmpVal[3] = {-1., 2., -1.};
  int tmpCol[3];
  for (int i = 0; i < nRows; ++i) {
    tmpCol[0] = myGlobRows[i] - 1;
    tmpCol[1] = myGlobRows[i];
    tmpCol[2] = myGlobRows[i] + 1;

    if (myGlobRows[i] == 0)
      A.InsertGlobalValues(myGlobRows[i], 2, &tmpVal[1], &tmpCol[1]);
    else if (myGlobRows[i] == nGlobRows - 1)
      A.InsertGlobalValues (myGlobRows[i], 2, tmpVal, tmpCol);
    else 
      A.InsertGlobalValues (myGlobRows[i], 3, tmpVal, tmpCol);
  }
  A.FillComplete();

  typedef Tpetra::Map<>::node_type znode_t;
  typedef Xpetra::EpetraMapT<int, znode_t> xemap_t;

  Teuchos::RCP<const Epetra_BlockMap> ebmap = Teuchos::rcpFromRef(A.RowMap());
  Teuchos::RCP<const xemap_t> xmap(new xemap_t(ebmap));
  
  const Teuchos::RCP<const Teuchos::Comm<int> > &xcomm = xmap->getComm();

  std::cout << "Teuchos:  Hello from " 
            << tcomm->getRank() << " of " 
            << tcomm->getSize() << std::endl;
  std::cout << "Epetra:   Hello from " 
            << ecomm->MyPID() << " of " 
            << ecomm->NumProc() << std::endl;
  std::cout << "Xpetra:   Hello from " 
            << xcomm->getRank() << " of " 
            << xcomm->getSize() << std::endl;
}
#endif
