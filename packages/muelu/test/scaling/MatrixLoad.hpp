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
#ifndef MATRIXLOAD_HPP
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <unistd.h>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_IO.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>



// This is a standard Galeri-or-MatrixFile loading routine designed to be shared between the various scaling tests

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixLoad(Teuchos::RCP<const Teuchos::Comm<int> > &comm,  Xpetra::UnderlyingLib& lib,
                bool binaryFormat,const std::string & matrixFile, const std::string & rhsFile, 
                const std::string & rowMapFile, 
                const std::string & colMapFile, 
                const std::string & domainMapFile, 
                const std::string & rangeMapFile, 
                const std::string & coordFile, const std::string &nullFile,
                Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >          & map, 
                Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >      & A, 
                Teuchos::RCP<Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> > & coordinates,
                Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & nullspace,
                Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & X, 
                Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & B,
                Galeri::Xpetra::Parameters<GlobalOrdinal> & galeriParameters,  Xpetra::Parameters & xpetraParameters, 
                std::ostringstream & galeriStream) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();


  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
  galeriStream << "========================================================\n" << xpetraParameters << galeriParameters;
  if (matrixFile.empty()) {

    // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
    //                                 d1  d2  d3
    //                                 d4  d5  d6
    //                                 d7  d8  d9
    //                                 d10 d11 d12
    // A perfect distribution is only possible when the #processors is a perfect square.
    // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
    // size. For example, np=14 will give a 7-by-2 distribution.
    // If you don't want Galeri to do this, specify mx or my on the galeriList.
    std::string matrixType = galeriParameters.GetMatrixType();

    // Create map and coordinates
    // In the future, we hope to be able to first create a Galeri problem, and then request map and coordinates from it
    // At the moment, however, things are fragile as we hope that the Problem uses same map and coordinates inside
    if (matrixType == "Laplace1D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double,LO,GO,Map,RealValuedMultiVector>("1D", map, galeriList);

    } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
               matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double,LO,GO,Map,RealValuedMultiVector>("2D", map, galeriList);

    } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double,LO,GO,Map,RealValuedMultiVector>("3D", map, galeriList);
    }

    // Expand map to do multiple DOF per node for block problems
    if (matrixType == "Elasticity2D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 2);
    if (matrixType == "Elasticity3D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 3);

    galeriStream << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
                 << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
                 << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
                 << "========================================================" << std::endl;

    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
      // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
      galeriList.set("right boundary" , "Neumann");
      galeriList.set("bottom boundary", "Neumann");
      galeriList.set("top boundary"   , "Neumann");
      galeriList.set("front boundary" , "Neumann");
      galeriList.set("back boundary"  , "Neumann");
    }

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);
    A = Pr->BuildMatrix();

    if (matrixType == "Elasticity2D" ||
        matrixType == "Elasticity3D") {
      nullspace = Pr->BuildNullspace();
      A->SetFixedBlockSize((galeriParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);
    }

  } else {
    if (!rowMapFile.empty())
      map = Xpetra::IO<SC,LO,GO,Node>::ReadMap(rowMapFile, lib, comm);
    comm->barrier();

    if (!binaryFormat && !map.is_null()) {
      RCP<const Map> colMap    = (!colMapFile.empty()    ? Xpetra::IO<SC,LO,GO,Node>::ReadMap(colMapFile,    lib, comm) : Teuchos::null);
      RCP<const Map> domainMap = (!domainMapFile.empty() ? Xpetra::IO<SC,LO,GO,Node>::ReadMap(domainMapFile, lib, comm) : Teuchos::null);
      RCP<const Map> rangeMap  = (!rangeMapFile.empty()  ? Xpetra::IO<SC,LO,GO,Node>::ReadMap(rangeMapFile,  lib, comm) : Teuchos::null);
      A = Xpetra::IO<SC,LO,GO,Node>::Read(matrixFile, map, colMap, domainMap, rangeMap);

    } else {
      A = Xpetra::IO<SC,LO,GO,Node>::Read(matrixFile, lib, comm, binaryFormat);

      if (!map.is_null()) {
        RCP<Matrix> newMatrix = MatrixFactory::Build(map, 1);
        RCP<Import> importer  = ImportFactory::Build(A->getRowMap(), map);
        newMatrix->doImport(*A, *importer, Xpetra::INSERT);
        newMatrix->fillComplete();

        A.swap(newMatrix);
      }
    }
    map = A->getMap();

    comm->barrier();

    if (!coordFile.empty()) {
      // NOTE: currently we only allow reading scalar matrices, thus coordinate
      // map is same as matrix map
      coordinates = Xpetra::IO<double,LO,GO,Node>::ReadMultiVector(coordFile, map);
    }

    if (!nullFile.empty())
      nullspace = Xpetra::IO<SC,LO,GO,Node>::ReadMultiVector(nullFile, map);
  }

  X = VectorFactory::Build(map);
  B = VectorFactory::Build(map);

  if (rhsFile.empty()) {
    // we set seed for reproducibility
    Utilities::SetRandomSeed(*comm);
    X->randomize();
    A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

    Teuchos::Array<typename STS::magnitudeType> norms(1);
    B->norm2(norms);
    B->scale(one/norms[0]);

  } else {
    // read in B
    B = Xpetra::IO<SC,LO,GO,Node>::ReadMultiVector(rhsFile, map);
  }
  galeriStream << "Galeri complete.\n========================================================" << std::endl;
}

#endif
