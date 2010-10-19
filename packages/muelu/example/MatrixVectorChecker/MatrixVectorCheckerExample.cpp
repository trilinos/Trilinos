#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Using Galeri:
// #include <Galeri_Maps.h>
// #include <Galeri_CrsMatrices.h>

// Using MueLu gallery:
#include <Cthulhu_Map.hpp>
#include <Cthulhu_CrsMatrix.hpp>
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"

#include "MatrixVectorChecker.hpp"

// Define template types
typedef double ScalarType;
typedef int    LocalOrdinal;
typedef int    GlobalOrdinal;
typedef Kokkos::DefaultNode::DefaultNodeType Node;

// Convenient shortcut for template types
typedef ScalarType    SC;
typedef LocalOrdinal  LO;
typedef GlobalOrdinal GO;
typedef Node          NO;

// Get ride of template parameters for some classes
typedef Cthulhu::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
typedef Cthulhu::CrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
typedef Cthulhu::MultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVector;

// Use Tpetra
#include <Cthulhu_TpetraMap.hpp>
#include <Cthulhu_TpetraCrsMatrix.hpp>
typedef Cthulhu::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> MyMap;
typedef Cthulhu::TpetraCrsMatrix<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MyCrsMatrix;
typedef Cthulhu::TpetraMultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MyMultiVector;

int main(int argc, char *argv[])
{
  
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::ParameterList paramList;
  paramList.set("nx", 10 * comm->getSize());
  paramList.set("ny", 10);
  
  paramList.set("mx", comm->getSize());
  paramList.set("my", 1);

  // Using Galeri:
  //
  // RCP<Tpetra_Map> map = rcp( Galeri::CreateMap("Cartesian2D", comm, paramList) );
  // RCP<Tpetra_CrsMatrix> matrix = rcp( Galeri::CreateCrsMatrix("Laplace2D", map.get(), paramList) );

  // Using MueLu gallery:
  GO nx = paramList.get<GO>("nx");
  GO ny = paramList.get<GO>("ny");

  RCP<const Map> map = rcp( new MyMap(nx*ny, 0, comm) );
  RCP<const CrsMatrix> matrix = CreateCrsMatrix<SC,LO,GO,NO>("Laplace2D", map, paramList);

  MatrixVectorChecker<SC,LO,GO,NO>(matrix);

  return EXIT_SUCCESS;

}

// JG TODO:
// - add a method CreateMap for cthulhu/gallery (as Galeri)
// - wrap galeri matrix for the new MatrixVectorChecker
