// Program to debug segfaults being reported in CDASH when
//  -D Kokkos_ENABLE_Pthreadi:BOOL=ON 
//  -D Tpetra_INST_PTHREAD:BOOL=ON 
// Problem appears to be in creation of Xpetra::EpetraCrsMatrixT

#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Epetra_Map.h>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_EpetraUtils.hpp>


int main(int narg, char **arg)
{
  using Teuchos::rcp;

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > tcomm = Tpetra::getDefaultComm();
  Teuchos::RCP<const Epetra_Comm> ecomm = Xpetra::toEpetra(tcomm);

  ////////////////////////
  // Build a boring matrix

  const int nGlobRows = 50;
  const Epetra_Map emap(nGlobRows, 0, *ecomm);
  Epetra_CrsMatrix emat(Copy, emap, 1, true);
  const double one = 1.;
  for (int i = 0; i < emat.NumMyRows(); i++) {
    int gid = emat.GCID(i);
    emat.InsertGlobalValues(gid, 1, &one, &gid);
  }
  emat.FillComplete();
  
  ////////////////////////////////////////////////////////////////////////////
  // Test whether conversion from Epetra_CrsMatrix to Xpetra::EpetraCrsMatrixT
  // gives a valid resulting matrix.

  std::cout << "Building Xpetra::EpetraCrsMatrixT from Epetra_CrsMatrix: "
            << std::endl;

  Teuchos::RCP<Epetra_CrsMatrix> ematrcp = Teuchos::rcpFromRef(emat);
  typedef Xpetra::EpetraCrsMatrixT<int, Tpetra::Map<>::node_type> xemat_t;
  Teuchos::RCP<const xemat_t> xmat;

  bool aok_mat = true;
  try {
    xmat = rcp(new xemat_t(ematrcp));
  }
  catch (std::exception &e) {
    std::cout << "Xpetra::EpetraCrsMatrixT threw an error " 
              << e.what() << std::endl;
    aok_mat = false;
  }
  
  if (aok_mat)
    std::cout << "Building Xpetra::EpetraCrsMatrixT from Epetra_CrsMatrix: "
              << "DONE with no errors caught " << std::endl;

  ////////////////////////////////////////////////////////////
  // Try the same thing with Epetra_Map and Xpetra::EpetraMapT

  std::cout << "Building Xpetra::EpetraMapT from Epetra_Map: "
            << std::endl;

  Teuchos::RCP<const Epetra_BlockMap> emaprcp = Teuchos::rcpFromRef(emap);
  typedef Xpetra::EpetraMapT<int, Tpetra::Map<>::node_type> xemap_t;
  Teuchos::RCP<const xemap_t> xmap;

  bool aok_map = true;
  try {
    xmap = rcp(new xemap_t(emaprcp));
  }
  catch (std::exception &e) {
    std::cout << "Xpetra::EpetraMapT threw an error " 
              << e.what() << std::endl;
    aok_map = false;
  }

  if (aok_map)
    std::cout << "Building Xpetra::EpetraMapT from Epetra_Map: "
              << "DONE with no errors caught " << std::endl;
  
  ///////////////////////////////////
  // Print some info from the classes

  std::cout << "Teuchos:                    Hello from " 
            << tcomm->getRank() << " of " 
            << tcomm->getSize() << std::endl;
  std::cout << "Epetra_CrsMatrix:           Hello from " 
            << ematrcp->Comm().MyPID() << " of " 
            << ematrcp->Comm().NumProc() << std::endl;
  std::cout << "Epetra_Map:                 Hello from " 
            << emaprcp->Comm().MyPID() << " of " 
            << emaprcp->Comm().NumProc() << std::endl;
  if (aok_mat)
    std::cout << "Xpetra::EpetraCrsMatrixT:   Hello from " 
              << xmat->getRowMap()->getComm()->getRank() << " of " 
              << xmat->getRowMap()->getComm()->getSize() << std::endl;
  if (aok_map)
    std::cout << "Xpetra::EpetraMapT:         Hello from " 
              << xmap->getComm()->getRank() << " of " 
              << xmap->getComm()->getSize() << std::endl;

  return 0;
}
