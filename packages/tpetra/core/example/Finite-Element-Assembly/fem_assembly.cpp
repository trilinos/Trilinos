
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>

#include "typedefs.hpp"
#include "MeshDatabase.hpp"



int main (int argc, char *argv[]) {
  using Teuchos::RCP;

  // MPI boilerplate
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

  // Generate the mesh
  int nex = 10;
  int ney = 10;
  MeshDatabase mesh(comm,nex,ney);


  // Build Tpetra Maps


  // Build graphs multiple ways


  // Build matrices 



  // Build RHS vectors



  return 0;
}

