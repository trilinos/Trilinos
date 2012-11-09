// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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

using namespace Galeri;

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
    string mapType = "Cartesian2D";
    
    // Container for parameters
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", 2 * comm->getSize()); 
    galeriList.set("ny", 2);
    galeriList.set("mx", comm->getSize());
    galeriList.set("my", 1);
    
    // Creation of the map
    RCP< ::Xpetra::Map<int, int, Kokkos::DefaultNode::DefaultNodeType> > map = Galeri::Xpetra::CreateMap<int, int, Kokkos::DefaultNode::DefaultNodeType>(lib, "Cartesian2D", comm, galeriList);

    // Print out the parameters
    cout << galeriList;
  
    // Print out the map
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    map->describe(*fos, Teuchos::VERB_EXTREME);
  }

  {
    // Creates an Xpetra::Map corresponding to a 3D Cartesian grid
    
    // Type of the object
    string mapType = "Cartesian3D";
    
    // Container for parameters
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", 2 * comm->getSize()); 
    galeriList.set("ny", 2);
    galeriList.set("nz", 2);
    
    // Creation of the map
    RCP< ::Xpetra::Map<int, int, Kokkos::DefaultNode::DefaultNodeType> > map = Galeri::Xpetra::CreateMap<int, int, Kokkos::DefaultNode::DefaultNodeType>(lib, "Cartesian3D", comm, galeriList);

    // Print out the parameters
    cout << galeriList;
  
    // Print out the map
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    map->describe(*fos, Teuchos::VERB_EXTREME);
  }

  return EXIT_SUCCESS;
}
