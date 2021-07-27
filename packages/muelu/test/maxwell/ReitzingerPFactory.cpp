// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>
#include <map>
#include <Teuchos_UnitTestRepository.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_IO.hpp>

// MueLu
#include <MueLu_config.hpp>
#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>
#include <MueLu_ReitzingerPFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_Exceptions.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_RAPFactory.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

// Copy & paste from MueLu_TestHelpers.cpp 
namespace MueLuTests {

  // static members initialization of the class TestHelpers::Parameters
  Xpetra::Parameters TestHelpers::Parameters::xpetraParameters = Xpetra::Parameters(Teuchos::UnitTestRepository::getCLP());

}

template<typename Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void read_matrix(Xpetra::UnderlyingLib & lib,RCP<const Teuchos::Comm<int> > & comm, 
                 RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & SM_Matrix, 
                 RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & D0_Matrix,
                 RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & Kn_Matrix,
                 RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> > & coords) {
#include <MueLu_UseShortNames.hpp>
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));


  std::string S_file, M1_file, D0_file;
  if (!Teuchos::ScalarTraits<Scalar>::isComplex) {
    S_file  = "S.mat";
    D0_file = "D0.mat";
    M1_file = "M1.mat";
  } else {
    S_file  = "S_complex.mat";
    D0_file = "D0_complex.mat";
    M1_file = "M1_complex.mat";
  }
  std::string coords_file = "coords.mat";

  // maps for nodal and edge matrices
  RCP<const Map> node_map;
  RCP<const Map> edge_map;

  // gradient matrix
  try {
    std::string base = D0_file.substr(0, D0_file.find_last_of('/')+1);
    std::string D0_filename = D0_file.substr(D0_file.find_last_of('/')+1, std::string::npos);
    std::string edgeMap_file = base + "rowmap_"    + D0_filename;
    std::string nodeMap_file = base + "domainmap_" + D0_filename;
    std::string colMap_file  = base + "colmap_"    + D0_filename;
    node_map = Xpetra::IO<SC, LO, GO, NO>::ReadMap(nodeMap_file, lib, comm);
    edge_map = Xpetra::IO<SC, LO, GO, NO>::ReadMap(edgeMap_file, lib, comm);
    RCP<const Map> colMap;
    if (comm->getSize() > 1)
      colMap = Xpetra::IO<SC, LO, GO, NO>::ReadMap(colMap_file, lib, comm);
    D0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(D0_file, edge_map, colMap, node_map, edge_map);
  } catch (const std::exception& e) {
    // *out << "Skipping D0 maps, because: " << e.what() << std::endl;
    D0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(D0_file, lib, comm);
    node_map = D0_Matrix->getDomainMap();
    edge_map = D0_Matrix->getRangeMap();
  }

  // build stiffness plus mass matrix (SM_Matrix)
  // edge stiffness matrix
  RCP<Matrix> S_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(S_file, edge_map);
  RCP<Matrix> M1_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M1_file, edge_map);
  Scalar one = Teuchos::ScalarTraits<SC>::one();
  Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*S_Matrix,false,one,*M1_Matrix,false,one,SM_Matrix,*out);
  SM_Matrix->fillComplete();


  // coordinates
  coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::ReadMultiVector(coords_file, node_map);

  // Generate Kn matrix
  {
    RCP<Matrix> temp = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(SM_Matrix->getRangeMap());
    Xpetra::MatrixMatrix<SC,LO,GO,NO>::Multiply(*SM_Matrix,false,*D0_Matrix,false,*temp,true,true);
    Kn_Matrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(D0_Matrix->getDomainMap());
    Xpetra::MatrixMatrix<SC,LO,GO,NO>::Multiply(*D0_Matrix,true,*temp,false,*Kn_Matrix,true,true);
  }

}

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ReitzingerPFactory, Setup2Level, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#   include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  
  RCP<Matrix> SM_Matrix, D0_Matrix, Kn_Matrix;
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coords;
  read_matrix<SC,LO,GO,NO>(lib,comm,SM_Matrix,D0_Matrix,Kn_Matrix,coords);

  printf("CMS: Aedge = %d x %d, Anode = %d x %d\n, D0 = %d x %d\n",
    (int)SM_Matrix->getGlobalNumRows(),
    (int)SM_Matrix->getGlobalNumCols(),
    (int)Kn_Matrix->getGlobalNumRows(),    
    (int)Kn_Matrix->getGlobalNumCols(),
    (int)D0_Matrix->getGlobalNumRows(),
    (int)D0_Matrix->getGlobalNumCols());
  

  int NumLevels = 2;
  // This guy works by generating a nodal hierarchy, copying the relevant data to the edge hierarchy
  // and then generating the edge hierarchy

  Hierarchy NodeH, EdgeH;
  NodeH.EnableGraphDumping("node_graph_0",0);
  EdgeH.EnableGraphDumping("edge_graph_0",0);
  NodeH.SetMaxCoarseSize(10);

  // Generate Node Hierarchy
  out<<"*** Setting Up Node Hierarchy *** "<<std::endl;
  {
    RCP<MueLu::Level> Finest = NodeH.GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);  
    Finest->Set("A", Kn_Matrix);

    // Level 0
    FactoryManager M0; // how to build aggregates and smoother of the first level
    M0.SetKokkosRefactor(false);
    
    // Level 1 (Plain aggregation)
    FactoryManager M1; 
    M1.SetKokkosRefactor(false);
    RCP<FactoryBase>        PnodalFact = rcp(new TentativePFactory());
    M1.SetFactory("P",PnodalFact);


    bool r;
    r = NodeH.Setup(0, Teuchos::null,  rcpFromRef(M0), rcpFromRef(M1));  TEST_EQUALITY(r, false);
    r = NodeH.Setup(1, rcpFromRef(M0), rcpFromRef(M1), Teuchos::null); TEST_EQUALITY(r, true);
    RCP<Level> l0 = NodeH.GetLevel(0);
    RCP<Level> l1 = NodeH.GetLevel(1);
  }  

  out<<"*** Copy Node->Edge Data *** "<<std::endl;
  // Copy Data to Edge Hierarchy
  for(int i=0; i<NumLevels; i++) {  
    EdgeH.AddNewLevel();
    RCP<Level> NodeL = NodeH.GetLevel(i);
    RCP<Level> EdgeL = EdgeH.GetLevel(i);

    EdgeL->Set("NodeMatrix",NodeL->Get<RCP<Matrix> >("A"));
    if(i!=0) 
      EdgeL->Set("Pnodal",NodeL->Get<RCP<Matrix> >("P"));    


    if(i==0) {
      EdgeL->Set("Coordinates", coords);
      EdgeL->Set("A", SM_Matrix);
      EdgeL->Set("D0", D0_Matrix);
    }
  }

  // Generate the Edge Hierarchy
  out<<"*** Setting Up Egde Hierarchy *** "<<std::endl;
  {
    FactoryManager M0; // how to build aggregates and smoother of the first level
    M0.SetKokkosRefactor(false);
    
    // Level 1 (Plain aggregation)
    FactoryManager M1; 
    M1.SetKokkosRefactor(false);
    RCP<ReitzingerPFactory> PedgeFact  = rcp(new ReitzingerPFactory());
    M1.SetFactory("P",PedgeFact);

    // Do the setup
    bool r;
    r = EdgeH.Setup(0, Teuchos::null,  rcpFromRef(M0), rcpFromRef(M1));  TEST_EQUALITY(r, false);
    r = EdgeH.Setup(1, rcpFromRef(M0), rcpFromRef(M1), Teuchos::null); TEST_EQUALITY(r, true);
  }

  
  RCP<Level> l0 = EdgeH.GetLevel(0);
  RCP<Level> l1 = EdgeH.GetLevel(1);

  TEST_EQUALITY(0,0);
}

# define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ReitzingerPFactory, Setup2Level, Scalar, LO, GO, Node) \



# include <MueLu_ETI_4arg.hpp>

}//end namespace




int main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

#ifdef HAVE_MUELU_KOKKOSCORE
  Kokkos::initialize(argc , argv);
#endif

  bool success = false;
  bool verbose = true;
  int ierr = -1;
  try {
    // Note: the command line parameter --linAlgebra= is take into account.
    // Xpetra parameters are added to the Teuchos::CommandLineProcessor of Teuchos::UnitTestRepository in MueLu_TestHelpers.cpp

#ifdef ParallelDebug
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    int mypid = comm->getRank();

    if (mypid  == 0) std::cout << "Host and Process Ids for tasks" << std::endl;
    for (int i = 0; i <comm->getSize(); i++) {
      if (i == mypid ) {
        char buf[80];
        char hostname[80];
        gethostname(hostname, sizeof(hostname));
        int pid = getpid();
        sprintf(buf, "Host: %s\tMPI rank: %d,\tPID: %d\n\tattach %d\n\tcontinue\n",
            hostname, mypid, pid, pid);
        printf("%s\n",buf);
        fflush(stdout);
        sleep(1);
      }
    }

    if (mypid == 0) {
      printf( "** Enter a character to continue > "); fflush(stdout);
      char go = ' ';
      scanf("%c",&go);
    }
    comm->barrier();
#endif

    ierr = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MUELU_KOKKOSCORE
  Kokkos::finalize();
#endif

  return (success ? ierr : EXIT_FAILURE);
}
