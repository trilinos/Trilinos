#include <iostream>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

#define CTHULHU_ENABLED //TODO

#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"

#include "Cthulhu_Map.hpp"
#include "Cthulhu_CrsMatrix.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraMap.hpp"
#include "Cthulhu_TpetraCrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#endif

#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraCrsMatrix.hpp"
#endif

// Define data types
typedef double SC;
typedef int    LO;
typedef int    GO;

using Teuchos::RCP;

// This is a test to see if the gallery is working with allowed matrices types.

int main(int argc, char* argv[]) 
{
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  GO nx=4;
  GO ny=4;
  GO nz=4;
  std::string matrixType("Laplace1D");

  Teuchos::ParameterList pl;
  GO numGlobalElements = nx*ny;
  LO indexBase = 0;

  Teuchos::ParameterList matrixList;
  matrixList.set("nx",nx);
  matrixList.set("ny",ny);
  matrixList.set("nz",nz);

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

#ifdef HAVE_CTHULHU_TPETRA
  // Tpetra::CrsMatrix
  {
    RCP<const Tpetra::Map<LO,GO> > map = rcp( new Tpetra::Map<LO,GO> (numGlobalElements, indexBase, comm) ); 
    RCP<Tpetra::CrsMatrix<SC,LO,GO> > A = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Tpetra::Map<LO,GO>, Tpetra::CrsMatrix<SC,LO,GO> > (matrixType,map,matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
#endif

#ifdef HAVE_CTHULHU_TPETRA
  // Cthulhu::TpetraCrsMatrix
  {
    RCP<const Cthulhu::TpetraMap<LO,GO> > map = rcp( new Cthulhu::TpetraMap<LO,GO> (numGlobalElements, indexBase, comm) );
    RCP<Cthulhu::TpetraCrsMatrix<SC,LO,GO> > A = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Cthulhu::TpetraMap<LO,GO>, Cthulhu::TpetraCrsMatrix<SC,LO,GO> > (matrixType,map,matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
#endif

#ifdef HAVE_CTHULHU_EPETRA
  // Cthulhu::EpetraCrsMatrix
  { 
    RCP<const Cthulhu::EpetraMap > map = rcp( new Cthulhu::EpetraMap (numGlobalElements, indexBase, comm) );
    RCP<Cthulhu::EpetraCrsMatrix> A = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Cthulhu::EpetraMap, Cthulhu::EpetraCrsMatrix> (matrixType,map,matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  } 
#endif

#ifdef HAVE_CTHULHU_TPETRA
  // Cthulhu::CrsMatrix (Tpetra)
  {
    RCP<const Cthulhu::Map<LO,GO> > map = rcp( new Cthulhu::TpetraMap<LO,GO> (numGlobalElements, indexBase, comm) ); 
    RCP<Cthulhu::CrsMatrix<SC,LO,GO> > A = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Cthulhu::Map<LO,GO>, Cthulhu::CrsMatrix<SC,LO,GO> >  (matrixType,map,matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  } 
#endif

#ifdef HAVE_CTHULHU_EPETRA
  // Cthulhu::CrsMatrix (Epetra)
  {
    RCP<const Cthulhu::Map<LO,GO> > map = rcp( new Cthulhu::EpetraMap (numGlobalElements, indexBase, comm) );
    RCP<Cthulhu::CrsMatrix<SC,LO,GO> > A = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Cthulhu::Map<LO,GO>, Cthulhu::CrsMatrix<SC,LO,GO> >  (matrixType,map,matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  } 
#endif

#ifdef HAVE_CTHULHU_TPETRA
  // Cthulhu::Operator (Tpetra)
  {
    RCP<const Cthulhu::Map<LO,GO> > map = rcp( new Cthulhu::TpetraMap<LO,GO> (numGlobalElements, indexBase, comm) );
    RCP<Cthulhu::Operator<SC,LO,GO> > A = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Cthulhu::Map<LO,GO>, Cthulhu::Operator<SC,LO,GO> >  (matrixType,map,matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  } 
#endif

#ifdef HAVE_CTHULHU_EPETRA
  // Cthulhu::Operator (Epetra)
  {
    RCP<const Cthulhu::Map<LO,GO> > map = rcp( new Cthulhu::EpetraMap (numGlobalElements, indexBase, comm) );
    RCP<Cthulhu::Operator<SC,LO,GO> > A = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Cthulhu::Map<LO,GO>, Cthulhu::Operator<SC,LO,GO> >  (matrixType,map,matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  } 
#endif
  
 return EXIT_SUCCESS;
}
