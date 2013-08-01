// This will be a driver to test the development and implementation of
// the Vanka Smoother in Ifpack2. 
//

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include "Teuchos_GlobalMPISession.hpp"
#include <Teuchos_DefaultComm.hpp>
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#else
#include <Teuchos_DefaultSerialComm.hpp>
#endif
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <Teuchos_FancyOStream.hpp>

#include <Ifpack2_Preconditioner.hpp>
#include <Ifpack2_Partitioner.hpp>
#include <Ifpack2_OverlappingPartitioner.hpp>
#include <Ifpack2_Details_UserPartitioner_def.hpp>
#include <Ifpack2_LinearPartitioner_def.hpp>
#include <Ifpack2_BlockRelaxation_def.hpp>
#include <Ifpack2_ILUT.hpp>
#include <Ifpack2_SparseContainer.hpp>
#include <Ifpack2_AdditiveSchwarz_def.hpp>

#include "Ifpack2_Parameters.hpp"

//#include <Kokkos_ConfigDefs.hpp>
//#include <Kokkos_SerialNode.hpp>

#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>


int main(int argc, char *argv[]){

  Teuchos::GlobalMPISession mpisess(&argc, &argv);\

  //Define communicator:
  Tpetra::DefaultPlatform::DefaultPlatformType& platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();

  
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef double Scalar;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

  
  //typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
  //typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crsgraph;
  //typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;

  // typedef our matrix/local inverse type to make life a little easier
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CRS;
  typedef Ifpack2::BlockRelaxation<CRS, Ifpack2::SparseContainer<CRS,Ifpack2::ILUT<CRS> > > BlockRelax;
  
  //typedef KokkosClassic::SerialNode node_type;
  using Teuchos::RCP;
  


  // Initialize a "FancyOStream" to output to standard out (cout)
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);




  // **************************************** //
  // 1D Poisson Test                          //
  // **************************************** //

  // ************************************************************************************************** //
  // Create our matrix!

  // Parameters
  GlobalOrdinal numGlobalDOFs = 32;
  
  // Create a map
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > myMap = Tpetra::createUniformContigMap<LocalOrdinal,GlobalOrdinal>(numGlobalDOFs,comm);

  // Get update list and number of local equations from newly created map
  const size_t numMyDOFs = myMap->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> myGlobalDOFs = myMap->getNodeElementList();

  // Create a CrsMatrix using the map
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(myMap,3));

  // Add rows
  for (size_t ii = 0; ii < numMyDOFs; ii++){
    if (myGlobalDOFs[ii] == 0) { //left boundary
      A->insertGlobalValues(myGlobalDOFs[ii],
                            Teuchos::tuple<GlobalOrdinal>(myGlobalDOFs[ii],myGlobalDOFs[ii]+1),
                            Teuchos::tuple<Scalar>(2.0,-1.0));
    }
    else if (myGlobalDOFs[ii] == numGlobalDOFs - 1) { //right boundary
      A->insertGlobalValues(myGlobalDOFs[ii],
                            Teuchos::tuple<GlobalOrdinal>(myGlobalDOFs[ii]-1,myGlobalDOFs[ii]),
                            Teuchos::tuple<Scalar>(-1.0,2.0));
    }
    else { //interior
      A->insertGlobalValues(myGlobalDOFs[ii],
                            Teuchos::tuple<GlobalOrdinal>(myGlobalDOFs[ii]-1,myGlobalDOFs[ii],myGlobalDOFs[ii]+1),
                            Teuchos::tuple<Scalar>(-1.0,2.0,-1.0));
    }
  }

  // Complete the fill
  A->fillComplete();

  // ************************************************************************************************** //
  // Set up ifpack2 parameters

  Teuchos::ParameterList MyList;
  // Distribute the local dofs linearly across nonoverlapping partitions
  MyList.set("partitioner: type"       ,"linear");
  // Distribute the local dofs over how many partitions?
  MyList.set("partitioner: local parts",(LocalOrdinal) A->getNodeNumRows());
  // How much overlap in the partitions
  MyList.set("partitioner: overlap"    ,(int) 1);

  // What type of block relaxation should we use?
  MyList.set("relaxation: type"        ,"Gauss-Seidel");
  // How many sweeps?
  MyList.set("relaxation: sweeps"      ,(int) 1);
  // How much damping?
  MyList.set("relaxation: damping factor",1.0);
  // Do we have a zero initial guess
  MyList.set("relaxation: zero starting solution", false);

  // Solving local blocks with ILUT, so we need some parameters
  MyList.set("fact: absolute threshold", 0.0);
  MyList.set("fact: relative threshold", 1.0);
  MyList.set("fact: ilut level-of-fill", 3.0);

  // Schwarz for the parallel
  MyList.set("schwarz: overlap level", (int) 1);



  // Create our preconditioner
  //  Ifpack2::BlockRelaxation<CRS, Ifpack2::SparseContainer<CRS,Ifpack2::ILUT<CRS> > > prec(A);
  Ifpack2::AdditiveSchwarz<CRS,BlockRelax> prec(A);
  

  // Set the parameters
  prec.setParameters(MyList);
  // Initialize everything
  prec.initialize();
  // Compute - this will form the blocks and factor the local matrices
  prec.compute();

  // Define RHS / Initial Guess
  Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > X = Tpetra::createVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(myMap);
  Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > B = Tpetra::createVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(myMap);
  
  Teuchos::ScalarTraits<Scalar>::seedrandom(846930883);
  B->putScalar((Scalar) 2.0);  
  //B->randomize();
  X->putScalar((Scalar) 0.0);
  //X->randomize();

  
  // Apply the preconditioner
  prec.apply(*B,*X);

  // This writes B = B - A*x;
  A->apply(*X,*B,Teuchos::NO_TRANS,-1.0,1.0);

  // Print X and B (multi)vectors
  //X->describe(*out,Teuchos::VERB_EXTREME);
  //B->describe(*out,Teuchos::VERB_EXTREME);
  
  // Print final residual norm
  std::cout << B->norm2() << std::endl;

  // Print stuff about our precondtioner (timings and whatnot)
  prec.describe(*out,Teuchos::VERB_EXTREME);

  

}//int main
