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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixLoad(Teuchos::RCP<const Teuchos::Comm<int> >& comm, Xpetra::UnderlyingLib& lib,
                bool binaryFormat, const std::string& matrixFile, const std::string& rhsFile,
                const std::string& rowMapFile,
                const std::string& colMapFile,
                const std::string& domainMapFile,
                const std::string& rangeMapFile,
                const std::string& coordFile,
                const std::string& coordMapFile, const std::string& nullFile, const std::string& materialFile,
                Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& map,
                Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A,
                Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> >& coordinates,
                Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& nullspace,
                Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& material,
                Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& X,
                Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& B,
                const int numVectors,
                Galeri::Xpetra::Parameters<GlobalOrdinal>& galeriParameters, Xpetra::Parameters& xpetraParameters,
                std::ostringstream& galeriStream) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
  galeriStream << "========================================================\n"
               << xpetraParameters;
  if (matrixFile.empty()) {
    galeriStream << galeriParameters;

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
    if (matrixType == "Laplace1D" || matrixType == "Identity") {
      map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real_type, LO, GO, Map, RealValuedMultiVector>("1D", map, galeriList);

    } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
               matrixType == "BigStar2D" || matrixType == "AnisotropicDiffusion" || matrixType == "Elasticity2D") {
      map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real_type, LO, GO, Map, RealValuedMultiVector>("2D", map, galeriList);

    } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
      map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real_type, LO, GO, Map, RealValuedMultiVector>("3D", map, galeriList);
    }

    // Expand map to do multiple DOF per node for block problems
    if (matrixType == "Elasticity2D")
      map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 2);
    if (matrixType == "Elasticity3D")
      map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 3);

    galeriStream << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
                 << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
                 << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
                 << "========================================================" << std::endl;

    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
      // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
      galeriList.set("right boundary", "Neumann");
      galeriList.set("bottom boundary", "Neumann");
      galeriList.set("top boundary", "Neumann");
      galeriList.set("front boundary", "Neumann");
      galeriList.set("back boundary", "Neumann");
    }

    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);
    A         = Pr->BuildMatrix();
    nullspace = Pr->BuildNullspace();
    //  The coordinates used by Galeri here might not match the coordinates that come into this function.
    //  In particular, they might correspond to different stretch factors. To fix this, we overwrite the
    //  coordinate array with those that Galeri now provides.

    if (!coordinates.is_null() && (matrixType == "Elasticity2D" || matrixType == "Elasticity3D")) {
      Teuchos::RCP<Xpetra::MultiVector<real_type, LocalOrdinal, GlobalOrdinal, Node> > newcoordinates;
      newcoordinates = Pr->BuildCoords();

      // Galeri makes multiple copies of coordinates to deal with
      // some issues when Ndofs != Nmeshnodes

      for (size_t kkk = 0; kkk < coordinates->getNumVectors(); kkk++) {
        Teuchos::ArrayRCP<real_type> old     = coordinates->getDataNonConst(kkk);
        Teuchos::ArrayRCP<real_type> newvals = newcoordinates->getDataNonConst(kkk);
        int numCopies                        = newvals.size() / old.size();
        for (int jj = 0; jj < old.size(); jj++) old[jj] = newvals[numCopies * jj];
      }
    }

    if (matrixType == "Elasticity2D" ||
        matrixType == "Elasticity3D") {
      A->SetFixedBlockSize((galeriParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);
    }

  } else {
    if (!rowMapFile.empty())
      map = Xpetra::IO<SC, LO, GO, Node>::ReadMap(rowMapFile, lib, comm);
    comm->barrier();

    {
      RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1a - Matrix read")));
      if (!binaryFormat && !map.is_null()) {
        RCP<const Map> colMap    = (!colMapFile.empty() ? Xpetra::IO<SC, LO, GO, Node>::ReadMap(colMapFile, lib, comm) : Teuchos::null);
        RCP<const Map> domainMap = (!domainMapFile.empty() ? Xpetra::IO<SC, LO, GO, Node>::ReadMap(domainMapFile, lib, comm) : Teuchos::null);
        RCP<const Map> rangeMap  = (!rangeMapFile.empty() ? Xpetra::IO<SC, LO, GO, Node>::ReadMap(rangeMapFile, lib, comm) : Teuchos::null);
        A                        = Xpetra::IO<SC, LO, GO, Node>::Read(matrixFile, map, colMap, domainMap, rangeMap);

      } else {
        A = Xpetra::IO<SC, LO, GO, Node>::Read(matrixFile, lib, comm, binaryFormat);
      }
      comm->barrier();
    }  // timer scope

    // If no rowmap file has been provided and the driver is being run in parallel,
    // create a uniformly distributed map and use it as A's row map.
    if (map.is_null() && comm->getSize() > 1) {
      RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1b - Matrix redistribute")));
      if (comm->getRank() == 0)
        std::cout << "No rowmap file specified, redistributing matrix using a uniformly distributed rowmap." << std::endl;
      map = Xpetra::MapFactory<LO, GO, Node>::Build(lib, A->getRowMap()->getGlobalNumElements(), (int)0, comm);

      RCP<Matrix> newMatrix = MatrixFactory::Build(map, 1);
      RCP<Import> importer  = ImportFactory::Build(A->getRowMap(), map);
      newMatrix->doImport(*A, *importer, Xpetra::INSERT);
      newMatrix->fillComplete();

      A.swap(newMatrix);
    } else {
      map = A->getRowMap();
    }

    comm->barrier();

    if (!coordFile.empty()) {
      RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1c - Read coordinates")));
      RCP<const Map> coordMap;
      if (!coordMapFile.empty())
        coordMap = Xpetra::IO<SC, LO, GO, Node>::ReadMap(coordMapFile, lib, comm);
      else
        coordMap = map;
      coordinates = Xpetra::IO<real_type, LO, GO, Node>::ReadMultiVector(coordFile, coordMap);
      comm->barrier();
    }

    if (!nullFile.empty())
      nullspace = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector(nullFile, map);

    if (!materialFile.empty())
      material = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector(materialFile, map);
  }

  X = MultiVectorFactory::Build(map, numVectors);
  B = MultiVectorFactory::Build(map, numVectors);

  if (rhsFile.empty()) {
    // we set seed for reproducibility
    Utilities::SetRandomSeed(*comm);
    X->randomize();
    A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

    Teuchos::Array<real_type> norms(numVectors);
    B->norm2(norms);
    B->scale(one / norms[0]);

  } else {
    // read in B
    RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1d - Read RHS")));
    B                   = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector(rhsFile, map);
    comm->barrier();
  }
  galeriStream << "Galeri complete.\n========================================================" << std::endl;
}

// Read in block information when doing user based blocks smoothing via "partitioner: global ID parts"
// Block information is read from file userBlkFileName. This file is assumed to be 'MatrixMarket
// matrix coordinate real general' where each row corresponds defines one
// block. Specifically, each point dof j  residing in the ith block is
// indicated by having a nonzero (i,j) in userBlkFileName. The file reader
// makes the following additional assumptions
//     1) The first line in the file is a comment
//     2) all file entries defining blk k occur before
//        file entries defining blk k+1.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void readUserBlks(const std::string& userBlkFileName, const std::string& smootherOrCoarse, Teuchos::ParameterList& mueluList,
                  const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A) {
  //     userBlkFileName          MatrixMarket file defining blocks
  //     smootherOrCoarse         Whether blocks are associated with "smoother: params" or "coarse: params"
  //     mueluList                Parameter list that is modified to reflect block informaiton found in file
  //     A                        Matrix that block smoother will be applied to
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();

  TEUCHOS_TEST_FOR_EXCEPTION((smootherOrCoarse != "smoother") && (smootherOrCoarse != "coarse"),
                             std::runtime_error, "2nd argument to readUserBlks must be either \"smoother\" or \"coarse\" and not \"" + smootherOrCoarse + "\"");
  if ((!userBlkFileName.empty()) && (mueluList.isSublist(smootherOrCoarse + ": params"))) {
    bool hasSubdomSolver = mueluList.sublist(smootherOrCoarse + ": params").isSublist("subdomain solver parameters");
    if (hasSubdomSolver) {
      if (mueluList.sublist(smootherOrCoarse + ": params").sublist("subdomain solver parameters").isParameter("partitioner: type")) {
        if (mueluList.sublist(smootherOrCoarse + ": params").sublist("subdomain solver parameters").get<std::string>("partitioner: type") == "user") {
          FILE* fp;
          int retVal;
          int nBlks, nRows, nnzs, ch, row, col;
          int procId = comm->getRank();
          double val;

          /* read block information from file only to count the # of rows */
          /* we own within each block. We'll reopen the file later to fill*/
          /* the Teuchos::Array for "partitioner: global ID parts"        */

          fp = fopen(&userBlkFileName[0], "r");
          TEUCHOS_TEST_FOR_EXCEPTION(fp == NULL, std::runtime_error, userBlkFileName + " file not found");

          while ((ch = getc(fp) != '\n'))
            ;  // read first line

          retVal = fscanf(fp, "%d %d %d\n", &nBlks, &nRows, &nnzs);
          TEUCHOS_TEST_FOR_EXCEPTION(retVal != 3,
                                     std::runtime_error, "unable to parse nBlks, nRows, nnzs in user file " + userBlkFileName);

          TEUCHOS_TEST_FOR_EXCEPTION(nRows != (int)A->getRowMap()->getGlobalNumElements(),
                                     std::runtime_error, "number of global rows in " + userBlkFileName + " does not match those in A");

          Teuchos::ArrayRCP<int> myOwnedRowsPerBlock(nBlks, 0);

          for (int i = 0; i < nnzs; i++) {
            retVal = fscanf(fp, "%d %d %lf", &row, &col, &val);
            row--;
            col--;
            if (A->getRowMap()->getLocalElement((GlobalOrdinal)col) != Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) (myOwnedRowsPerBlock[row])++;
          }
          fclose(fp);

          // Assign ownership of each block to one mpirank corresponding
          // to which mpirank owns the most rows within a block. When
          // ties occur, highest procId wins.

          // first, set to procId+1 if we own the largest # of rows
          // in block.  Otherwise, set to 0.

          Teuchos::ArrayRCP<int> maxOwnedRowsPerBlock(nBlks, 0);
          Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, (LocalOrdinal)nBlks, myOwnedRowsPerBlock.getRawPtr(), maxOwnedRowsPerBlock.getRawPtr());

          Teuchos::ArrayRCP<int> haveMaxOwned(nBlks, 0);
          for (int i = 0; i < nBlks; i++) {
            haveMaxOwned[i] = 0;
            if ((myOwnedRowsPerBlock[i] == maxOwnedRowsPerBlock[i]) && (myOwnedRowsPerBlock[i] > 0)) haveMaxOwned[i] = procId + 1;
          }
          Teuchos::ArrayRCP<int> maxProcIdPlusOne(nBlks, 0);
          Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, (LocalOrdinal)nBlks, haveMaxOwned.getRawPtr(), maxProcIdPlusOne.getRawPtr());

          Teuchos::ArrayRCP<bool> ownBlk(nBlks, false);
          int nOwned            = 0;
          int maxDofsInAnyBlock = 0;
          for (int i = 0; i < nBlks; i++) {
            TEUCHOS_TEST_FOR_EXCEPTION(maxProcIdPlusOne[i] == 0, std::runtime_error,
                                       "some blocks are not owned by any processor?");
            if (maxProcIdPlusOne[i] == procId + 1) {
              ownBlk[i] = true;
              nOwned++;
              maxDofsInAnyBlock += myOwnedRowsPerBlock[i];
            }
          }
          maxDofsInAnyBlock *= 2;  // We only have an estimate as to the largest number of dofs in any block (as some dofs
          maxDofsInAnyBlock++;     // might reside on other processors. So we multiple our estimate by 2 hoping to be large
                                   // enough. Later, we check to see if things are not large enough and re-allocate space
                                   // in this case.

          Teuchos::Array<Teuchos::ArrayRCP<GlobalOrdinal> > blockLists(nOwned, Teuchos::null);

          Teuchos::Array<int> buffer(maxDofsInAnyBlock, 0);

          int curRow, currentOwnedRow = 0;

          // reopen userBlkFilename and record block information in blockLists

          fp = fopen(&userBlkFileName[0], "r");
          TEUCHOS_TEST_FOR_EXCEPTION(fp == NULL, std::runtime_error, userBlkFileName + " file not found");

          while ((ch = getc(fp) != '\n'))
            ;

          retVal = fscanf(fp, "%d %d %d\n", &nBlks, &nRows, &nnzs);
          TEUCHOS_TEST_FOR_EXCEPTION(retVal != 3,
                                     std::runtime_error, "unable to parse nBlks, nRows, nnzs in user file " + userBlkFileName);

          retVal = fscanf(fp, "%d %d %lf", &row, &col, &val);
          row    = row - 1;
          col    = col - 1;

          int jj = 1;
          while (jj <= nnzs) {
            if (ownBlk[row] == true) {
              int i  = 0;
              curRow = row;
              while ((row == curRow) && (jj <= nnzs)) {
                if (i == maxDofsInAnyBlock) {
                  maxDofsInAnyBlock *= 2;
                  buffer.resize(maxDofsInAnyBlock);
                }
                buffer[i] = col;
                i++;
                retVal = fscanf(fp, "%d %d %lf", &row, &col, &val);
                jj++;
                row = row - 1;
                col = col - 1;
              }
              blockLists[currentOwnedRow] = Teuchos::arcp<GlobalOrdinal>(i);
              for (int k = 0; k < i; k++) {
                blockLists[currentOwnedRow][k] = (GlobalOrdinal)buffer[k];
              }
              currentOwnedRow++;
            } else {
              retVal = fscanf(fp, "%d %d %lf", &row, &col, &val);
              jj++;
              row = row - 1;
              col = col - 1;
            }
          }
          fclose(fp);
          mueluList.sublist(smootherOrCoarse + ": params").sublist("subdomain solver parameters").set<Teuchos::Array<Teuchos::ArrayRCP<GlobalOrdinal> > >("partitioner: global ID parts", blockLists);
        }
      }
    }
  }
}

#endif
