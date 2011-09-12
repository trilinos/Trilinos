// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayView.hpp>

// Xpetra
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsOperatorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

// MueLu
#include "MueLu_Hierarchy.hpp"
#include "MueLu_Utilities.hpp"

// Define default template types
typedef double Scalar;
typedef int    LocalOrdinal;
typedef int    GlobalOrdinal;

int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  //
  // MPI initialization
  //

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Construct the problem
  //

  // Construct a Map that puts approximately the same number of equations on each processor.
  RCP<const Map<Ordinal> > map = Xpetra::createUniformContigMap<Ordinal,Ordinal>(numGlobalElements, comm);

  // Get update list and number of local equations from newly created map.
  const size_t numMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const Ordinal> myGlobalElements = map->getNodeElementList();

  // Create a CrsMatrix using the map, with a dynamic allocation of 3 entries per row
  RCP<CrsMatrix> A = Xpetra::CrsOperatorFactory(map,3);

  // Add rows one-at-a-time
  for (size_t i = 0; i < numMyElements; i++) {
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

  // Complete the fill, ask that storage be reallocated and optimized
  A->fillComplete();

  //
  // Construct AMG preconditioner
  //

  // AMG Hierarchy
  Hierarchy H(A);
  H.setDefaultVerbLevel(Teuchos::VERB_HIGH);

  // AMG Setup phase
  H->Setup();

  //
  // Solve Ax = b
  //

  RCP<Vector> X = VectorFactory::Build(map,1);
  RCP<Vector> B = VectorFactory::Build(map,1);
  
  X->putScalar((SC) 0.0);
  B->setSeed(846930886); B->randomize();

  // Use AMG directly as an iterative method
  int nIts = 9;

  H.Iterate(*B, its, *X);

  // Print relative residual norm
  ST::magnitudeType residualNorms = MueLu::Utils::ResidualNorm(A, X, B);
  std::out << "||Residual|| = " << std::setiosflags(ios::fixed) << std::setprecision(20) << residualNorms << std::endl;

  return EXIT_SUCCESS;
}
