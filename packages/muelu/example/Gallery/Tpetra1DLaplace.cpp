#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_Array.hpp>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include <algorithm>
#include <functional>

#include <string>
#include <sstream>
#include <fstream>
std::string PrintMemoryUsage() {
  std::ostringstream mem;
  std::ifstream proc("/proc/self/status");
  std::string s;
  while(getline(proc, s), !proc.fail()) {
    if(s.substr(0, 6) == "VmSize") {
      mem << s;
      return mem.str();
    }
  }
  return mem.str();
}

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  typedef double                                                  Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType            Magnitude;
  typedef int                                                     Ordinal;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType            Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
  typedef Tpetra::Map<Ordinal,Ordinal,Node>                       Map;
  typedef Tpetra::CrsMatrix<Scalar,Ordinal,Ordinal,Node>          CrsMatrix;
  using Teuchos::RCP;
  using Teuchos::tuple;

  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  RCP<Node>                      node = platform.getNode();
  const int myRank = comm->getRank();

  int numGlobalElements = 10000000;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("numGlobalElements",&numGlobalElements,"Global problem size.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  RCP<const Map> map = Tpetra::createUniformContigMap<Ordinal,Ordinal>(numGlobalElements, comm);
  const size_t numMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const Ordinal> myGlobalElements = map->getNodeElementList();

  std::cout << myRank << ": " << "Inital memory usage: " << PrintMemoryUsage() << std::endl;

  RCP<CrsMatrix> A = Tpetra::createCrsMatrix<Scalar>(map,3);

  std::cout << myRank << ": " << "Memory after CrsMatrix constructor: " << PrintMemoryUsage() << std::endl;

  for (size_t i=0; i<numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues( myGlobalElements[i],
                             tuple<Ordinal>( myGlobalElements[i], myGlobalElements[i]+1 ),
                             tuple<Scalar> ( 2.0, -1.0 ) );
    }
    else if (myGlobalElements[i] == numGlobalElements-1) {
      A->insertGlobalValues( myGlobalElements[i],
                             tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i] ),
                             tuple<Scalar> ( -1.0, 2.0 ) );
    }
    else {
      A->insertGlobalValues( myGlobalElements[i],
                             tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1 ),
                             tuple<Scalar> ( -1.0, 2.0, -1.0 ) );
    }
  }

  std::cout << myRank << ": " << "Memory after InsertGlobalValues(): " << PrintMemoryUsage() << std::endl;

  A->fillComplete(Tpetra::DoOptimizeStorage);

  std::cout << myRank << ": " << "Memory after FillComplete(): " << PrintMemoryUsage() << std::endl;

  return 0;
}
