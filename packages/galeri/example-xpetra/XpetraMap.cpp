// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_Parameters.hpp>

#include "Galeri_XpetraMaps.hpp"

#ifndef XPETRA_TEST_USE_LONGLONG_GO
#define GO int
#else
#define GO long long
#endif

int main(int argc, char *argv[]) {

  using Teuchos::RCP;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Read the --linAlgebra option
  ::Xpetra::UnderlyingLib lib;
  {
    // Note: use --help to list available options.
    Teuchos::CommandLineProcessor clp(false);
    ::Xpetra::Parameters xpetraParameters(clp);

    switch (clp.parse(argc,argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }
    if (comm->getRank() == 0) std::cout << xpetraParameters;

    lib = xpetraParameters.GetLib();
  }

  {
    // Creates an Xpetra::Map corresponding to a 2D Cartesian grid
    // on the unit square. For parallel runs, the nodes are divided into
    // strips, so that the total number of subdomains is comm->getSize() x 1.
    
    // Type of the object
    std::string mapType = "Cartesian2D";
    
    // Container for parameters
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", (GO) 2 * comm->getSize()); 
    galeriList.set("ny", (GO) 2);
    galeriList.set("mx", (GO) comm->getSize());
    galeriList.set("my", (GO) 1);
    
    // Creation of the map
    RCP< ::Xpetra::Map<int, GO, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType> > map = Galeri::Xpetra::CreateMap<int, GO, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>(lib, "Cartesian2D", comm, galeriList);

    // Print out the parameters
    cout << galeriList;
  
    // Print out the map
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    map->describe(*fos, Teuchos::VERB_EXTREME);
  }

  {
    // Creates an Xpetra::Map corresponding to a 3D Cartesian grid
    
    // Type of the object
    std::string mapType = "Cartesian3D";
    
    // Container for parameters
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", (GO) 2 * comm->getSize()); 
    galeriList.set("ny", (GO) 2);
    galeriList.set("nz", (GO) 2);
    
    // Creation of the map
    RCP< ::Xpetra::Map<int, GO, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType> > map = Galeri::Xpetra::CreateMap<int, GO, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>(lib, "Cartesian3D", comm, galeriList);

    // Print out the parameters
    cout << galeriList;
  
    // Print out the map
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    map->describe(*fos, Teuchos::VERB_EXTREME);
  }

  return EXIT_SUCCESS;
}
