// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

int main (int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using Teuchos::Array;
  using Teuchos::tuple;

  // Tpetra::ScopeGuard initializes MPI and Kokkos, if they haven't
  // been initialized already.  Its destructor will finalize MPI and
  // Kokkos.
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    // Put all Tpetra objects in an inner scope.  They must not be
    // allowed to persist past the end of the Tpetra::ScopeGuard
    // object's lifetime.
    auto comm = Tpetra::getDefaultComm ();

    int num_halo_expansions = 1;
    Teuchos::CommandLineProcessor clp;
    clp.addOutputSetupOptions(true);
    clp.setOption("halo-expansions", &num_halo_expansions ,
                  "Number of times to expand the halo for Additive Schwarz");
    switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    using LO = Tpetra::Map<>::local_ordinal_type;
    using GO = Tpetra::Map<>::global_ordinal_type;
    using Node = Tpetra::Map<>::node_type;
    using MAT = Tpetra::CrsMatrix<>;
    using Scalar = MAT::scalar_type;
    using ST = Teuchos::ScalarTraits<Scalar>;
    //using MV = Tpetra::MultiVector<>;
    using IMP = Tpetra::Import<>;

    const GO ONE = Teuchos::OrdinalTraits<GO>::one();
    const GO GO_INVALID = Teuchos::OrdinalTraits<GO>::invalid();

    RCP<MAT> A;
    {
      TimeMonitor tm (*TimeMonitor::getNewTimer("1) Matrix Fill"));

      const size_t numImages = comm->getSize();
      const size_t myImageID = comm->getRank();

      if (numImages < 2) return -1;
      // create a Map
      RCP<const Tpetra::Map<> > map =
        Tpetra::createContigMapWithNode<LO,GO,Node>(GO_INVALID, ONE, comm);
      /* create the following matrix:
         [2 1           ]
         [1 1 1         ]
         [  1 1 1       ]
         [   . . .      ]
         [     . . .    ]
         [       . . .  ]
         [         1 1 1]
         [           1 2]
         this matrix has an eigenvalue lambda=3, with eigenvector v = [1 ... 1]
      */
      A = rcp (new MAT(map, 3));
      if (myImageID == 0) {
        Array<Scalar> vals(tuple<Scalar>(static_cast<Scalar>(2)*ST::one(), ST::one()));
        Array<GO> cols(tuple<GO>(myImageID, myImageID+1));
        A->insertGlobalValues(myImageID,cols(),vals());
      }
      else if (myImageID == numImages-1) {
        Array<Scalar> vals(tuple<Scalar>(ST::one(), static_cast<Scalar>(2)*ST::one()));
        Array<GO> cols(tuple<GO>(myImageID-1,myImageID));
        A->insertGlobalValues(myImageID,cols(),vals());
      }
      else {
        Array<Scalar> vals(3,ST::one());
        Array<GO> cols(tuple<GO>(myImageID-1, myImageID, myImageID+1));
        A->insertGlobalValues(myImageID,cols(),vals());
      }
      A->fillComplete();
    }

    RCP<MAT> Mold, Mnew;
    {
      TimeMonitor tm (*TimeMonitor::getNewTimer("2) Halo Generation"));

      Mold = A;
      Mnew = Mold;

      for (int i=0; i<num_halo_expansions; ++i) {
        RCP<const IMP> rowImporter = Mold->getGraph()->getImporter();
        Mnew = Tpetra::importAndFillCompleteCrsMatrix<MAT>(Mold,*rowImporter);
        Mold = Mnew;
      }
    }

#if DEBUG_OUTPUT
    int rank = comm->getRank();
    {
      std::ofstream ofs1( std::string("orig.") + std::to_string(rank) + std::string(".dat") );
      auto out1 = Teuchos::getFancyOStream(Teuchos::rcpFromRef(ofs1));
      A.describe(*out1,Teuchos::VERB_EXTREME);
    }
    {
      std::ofstream ofs1( std::string("new.") + std::to_string(rank) + std::string(".dat") );
      auto out1 = Teuchos::getFancyOStream(Teuchos::rcpFromRef(ofs1));
      Mnew->describe(*out1,Teuchos::VERB_EXTREME);
    }
#endif

    if(comm->getRank() == 0) {
      // Tell the Trilinos test framework that the test passed.
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
  }
  return EXIT_SUCCESS;
}  // END main()
