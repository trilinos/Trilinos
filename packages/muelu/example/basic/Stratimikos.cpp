// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

/*
   Call MueLu via the Stratimikos interface.

Usage:
./MueLu_Stratimikos.exe : use xml configuration file stratimikos_ParameterList.xml

Note:
The source code is not MueLu specific and can be used with any Stratimikos strategy.
*/

// Teuchos includes
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_YamlParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>

// Stratimikos includes
#include <Stratimikos_LinearSolverBuilder.hpp>
#include <Stratimikos_MueLuHelpers.hpp>

// Xpetra include
#include <Xpetra_Parameters.hpp>

// MueLu includes
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include <MatrixLoad.hpp>

// Galeri includes
#include <Galeri_XpetraParameters.hpp>

template <typename Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  TEUCHOS_ASSERT(lib == Xpetra::UseTpetra);
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::coordinateType real_type;
  typedef Xpetra::MultiVector<real_type, LocalOrdinal, GlobalOrdinal, Node> RealValuedMultiVector;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  try {
    //
    // MPI initialization
    //
    RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

    //
    // Parameters
    //
    // manage parameters of the test case
    Galeri::Xpetra::Parameters<GlobalOrdinal> matrixParameters(clp, 100, 100, 100, "Laplace2D");
    // manage parameters of Xpetra
    Xpetra::Parameters xpetraParameters(clp);

    // command line parameters
    std::string xmlFileName = "stratimikos_ParameterList.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from an xml file");
    std::string yamlFileName = "";
    clp.setOption("yaml", &yamlFileName, "read parameters from a yaml file");
    bool printTimings = false;
    clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
    bool use_stacked_timer = false;
    clp.setOption("stacked-timer", "no-stacked-timer", &use_stacked_timer, "Run with or without stacked timer output");
    std::string timingsFormat = "table-fixed";
    clp.setOption("time-format", &timingsFormat, "timings format (table-fixed | table-scientific | yaml)");
    bool binaryFormat = false;
    clp.setOption("binary", "ascii", &binaryFormat, "read matrices in binary format");
    std::string rowMapFile;
    clp.setOption("rowmap", &rowMapFile, "map data file");
    std::string colMapFile;
    clp.setOption("colmap", &colMapFile, "colmap data file");
    std::string domainMapFile;
    clp.setOption("domainmap", &domainMapFile, "domainmap data file");
    std::string rangeMapFile;
    clp.setOption("rangemap", &rangeMapFile, "rangemap data file");
    std::string matrixFile;
    clp.setOption("matrix", &matrixFile, "matrix data file");
    std::string rhsFile;
    clp.setOption("rhs", &rhsFile, "rhs data file");
    std::string coordFile;
    clp.setOption("coords", &coordFile, "coordinates data file");
    std::string coordMapFile;
    clp.setOption("coordsmap", &coordMapFile, "coordinates map data file");
    std::string nullFile;
    clp.setOption("nullspace", &nullFile, "nullspace data file");
    std::string blockNumberFile;
    clp.setOption("blocknumber", &blockNumberFile, "block number data file");
    std::string materialFile;
    clp.setOption("material", &materialFile, "material data file");
    int numVectors = 1;
    clp.setOption("multivector", &numVectors, "number of rhs to solve simultaneously");
    int numSolves = 1;
    clp.setOption("numSolves", &numSolves, "number of times the system should be solved");
    bool useGaleriForMatrixConstruction = true;
    clp.setOption("useGaleriForMatrixConstruction", "useDirectMatrixConstruction", &useGaleriForMatrixConstruction, "Use matrices from Galeri or from direct inline assembly");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream &out       = *fancy;
    out.setOutputToRootOnly(0);

    // Set up timers
    Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
    if (use_stacked_timer)
      stacked_timer = rcp(new Teuchos::StackedTimer("Main"));
    TimeMonitor::setStackedTimer(stacked_timer);

    // Read in parameter list
    TEUCHOS_TEST_FOR_EXCEPTION(xmlFileName == "" && yamlFileName == "", std::runtime_error,
                               "Need to provide xml or yaml input file");
    RCP<ParameterList> paramList = rcp(new ParameterList("params"));
    if (yamlFileName != "")
      Teuchos::updateParametersFromYamlFileAndBroadcast(yamlFileName, paramList.ptr(), *comm);
    else
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, paramList.ptr(), *comm);

    //
    // Construct the problem
    //

    using TpMap                   = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
    using TpMatrix                = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using TpMultiVector           = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using TpRealValuedMultiVector = Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>;
    using TpLOVector              = Tpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;

    RCP<const TpMatrix> A;
    RCP<const TpMap> map;
    RCP<TpRealValuedMultiVector> coordinates;
    RCP<TpMultiVector> nullspace;
    RCP<TpMultiVector> material;
    RCP<TpLOVector> blocknumber;
    RCP<TpMultiVector> X;
    RCP<const TpMultiVector> B;

    if (useGaleriForMatrixConstruction) {
      /////////////////////////////////////////////////////////////////////////
      // Use Galeri package for assembly of matrix and vectors.

      RCP<Matrix> xpetraA;
      RCP<const Map> xpetraMap;
      RCP<RealValuedMultiVector> xpetraCoordinates;
      RCP<MultiVector> xpetraNullspace;
      RCP<MultiVector> xpetraMaterial;
      RCP<LOVector> xpetraBlocknumber;
      RCP<MultiVector> xpetraX;
      RCP<MultiVector> xpetraB;

      std::ostringstream galeriStream;
      MatrixLoad<SC, LocalOrdinal, GlobalOrdinal, Node>(comm, lib, binaryFormat, matrixFile, rhsFile, rowMapFile, colMapFile, domainMapFile, rangeMapFile, coordFile, coordMapFile, nullFile, materialFile, blockNumberFile, xpetraMap, xpetraA, xpetraCoordinates, xpetraNullspace, xpetraMaterial, xpetraBlocknumber, xpetraX, xpetraB, numVectors, matrixParameters, xpetraParameters, galeriStream);
      out << galeriStream.str();

      A   = toTpetra(xpetraA);
      map = toTpetra(xpetraMap);
      if (!xpetraCoordinates.is_null())
        coordinates = toTpetra(xpetraCoordinates);
      if (!xpetraNullspace.is_null())
        nullspace = toTpetra(xpetraNullspace);
      if (!xpetraMaterial.is_null())
        material = toTpetra(xpetraMaterial);
      if (!xpetraBlocknumber.is_null())
        blocknumber = toTpetra(xpetraBlocknumber);
      X = toTpetra(xpetraX);
      B = toTpetra(xpetraB);

      X->putScalar(0);

    } else {
      /////////////////////////////////////////////////////////////////////////
      // We demonstrate how to set up a sparse matrix using 3 arrays.
      // This can serve as a starting point for using Trilinos solvers with a
      // matrix that has been assembled outside of Trilinos.

      // We need to provide the following inputs:
      //
      // - three arrays describing the local part of the sparse CRS matrix
      //   (row offsets, column indices, values),
      // - mappings from local indices to global indices for rows and columns.

      // Kokkos memory space
      using memory_space = typename Node::memory_space;

      // Types of the three data arrays.
      using rowptr_type  = typename TpMatrix::row_ptrs_device_view_type::non_const_type;
      using indices_type = typename TpMatrix::local_inds_device_view_type::non_const_type;
      using values_type  = typename TpMatrix::values_device_view_type::non_const_type;

      // The three arrays describing the sparse matrix.
      // These could just be unmanaged views for the memory locations of user provided data, see:
      // https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/View.html#unmanaged-views
      rowptr_type rowptr;
      indices_type indices;
      values_type values;

      // The mappings between local indices and global indices for columns of the matrix.
      RCP<const TpMap> colmap;

      // In this example, we set up a 1D Laplacian with Dirichlet conditions eliminated from
      // the system. This is a tri-diagonal matrix:
      //
      //  |  2 -1           |
      //  | -1  2 -1        |
      //  |    -1  2 -1     |
      //  |        ...      |
      //  |        -1  2 -1 |
      //  |           -1  2 |

      // This matrix is distributed row-wise over the ranks of the MPI
      // communicator, i.e. like this
      //
      //  |  2 -1             |  rank 0
      //  | -1  2 -1          |
      //  |    -1  2 -1       |
      //  |        ...        |
      // -----------------------------------------
      //  |           ...     |
      //  |         -1  2 -1  |  rank 1
      //  |             ...   |
      // -----------------------------------------
      //  ...
      // -----------------------------------------
      //  |        ...        |  rank P-1
      //  |        -1  2 -1   |
      //  |           -1  2   |

      // Each rank will own the same number of rows, so we split the global
      // index space into equal size parts. This information will get encode in the
      // row map of the matrix.
      //
      // Each rank (expect rank 0 and P-1) also has two off-rank column indices.
      // This information gets encode in the column map of the operator. For more
      // details regarding the concept of maps, have a look at the Tpetra documentation.

      // number of rows on each rank
      LocalOrdinal localNumRows = matrixParameters.GetNumGlobalElements() / comm->getSize();

      // Each row has 3 entries:
      LocalOrdinal localNNZ = 3 * localNumRows;
      // Except for the first and last row.
      if (comm->getRank() == 0)
        --localNNZ;
      if (comm->getRank() + 1 == comm->getSize())
        --localNNZ;

      // Allocate memory for local matrix
      rowptr  = rowptr_type(Kokkos::ViewAllocateWithoutInitializing("rowptr"), localNumRows + 1);
      indices = indices_type(Kokkos::ViewAllocateWithoutInitializing("indices"), localNNZ);
      values  = values_type(Kokkos::ViewAllocateWithoutInitializing("values"), localNNZ);

      {
        // Set up the map for the row indices.
        //
        // Simple map construction:
        // The data is distributed in contiguous fashion, i.e.
        // that the global ids are
        // on rank 0: [0, 1, ..., numLocalElements-1]
        // on rank 1: [numLocalElements, numLocalElements+1, ..., 2*numLocalElements-1]
        // ....
        map = rcp(new TpMap(
            /*numGlobalElements=*/Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
            localNumRows,
            /*indexBase=*/0, comm));
      }

      {
        // Set up the map for the column indices.
        //
        // This map is slightly more tricky because we need to keep track of the off-rank indices.

        // Each rank has 2 off-rank entries, except for first and last rank.
        LocalOrdinal localNumColumnMapEntries = localNumRows + 2;
        if (comm->getRank() == 0)
          --localNumColumnMapEntries;
        if (comm->getRank() + 1 == comm->getSize())
          --localNumColumnMapEntries;
        Kokkos::View<GlobalOrdinal *, memory_space> columnMapEntries("columnMap_entries",
                                                                     localNumColumnMapEntries);
        // We use the same local indices for the on-rank entries of the row map and the column map.
        // Copy global indices over for on-rank entries.
        Kokkos::deep_copy(
            Kokkos::subview(columnMapEntries, Kokkos::make_pair(0, localNumRows)),
            map->getMyGlobalIndicesDevice());

        // For the off-rank entries, we use local indices localNumRows and localNumRows+1.
        auto offset = localNumRows;
        if (comm->getRank() > 0) {
          Kokkos::deep_copy(Kokkos::subview(columnMapEntries, offset), map->getGlobalElement(0) - 1);
          ++offset;
        }
        if (comm->getRank() + 1 < comm->getSize()) {
          Kokkos::deep_copy(Kokkos::subview(columnMapEntries, offset), map->getGlobalElement(localNumRows - 1) + 1);
          ++offset;
        }

        // We now have constructed our local mapping between local indices and global indices
        // and call the Tpetra::Map constructor.
        colmap = rcp(new TpMap(
            /*numGlobalElements=*/Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
            columnMapEntries,
            /*indexBase=*/0, comm));
      }

      {
        // We fill the arrays describing the local part of the matrix.
        // This step can obviously be skipped if the user has pre-assembled arrays.

#if KOKKOS_VERSION >= 40799
        using ATS = KokkosKernels::ArithTraits<Scalar>;
#else
        using ATS      = Kokkos::ArithTraits<Scalar>;
#endif
        using impl_SC = typename ATS::val_type;
#if KOKKOS_VERSION >= 40799
        using impl_ATS = KokkosKernels::ArithTraits<impl_SC>;
#else
        using impl_ATS = Kokkos::ArithTraits<impl_SC>;
#endif

        LocalOrdinal offset         = (comm->getRank() == 0) ? 2 : 3;
        GlobalOrdinal numGlobalRows = map->getGlobalNumElements();
        auto lclRowMap              = map->getLocalMap();
        LocalOrdinal comm_rank      = comm->getRank();
        Kokkos::parallel_for(
            "matrix_construction",
            Kokkos::RangePolicy(0, localNumRows + 1),
            KOKKOS_LAMBDA(const LocalOrdinal &row) {
              if (row < localNumRows) {
                GlobalOrdinal globalRow = lclRowMap.getGlobalElement(row);

                LocalOrdinal k = (offset + (row - 1) * 3 >= 0) ? (offset + (row - 1) * 3) : 0;

                rowptr(row) = k;

                // left neighbor
                if (globalRow - 1 >= 0) {
                  LocalOrdinal col = (row - 1 >= 0) ? row - 1 : localNumRows;
                  indices(k)       = col;
                  values(k)        = -impl_ATS::one();
                  ++k;
                }

                // diagonal entry
                {
                  indices(k) = row;
                  values(k)  = 2 * impl_ATS::one();
                  ++k;
                }

                // right neighbor
                if (globalRow + 1 < numGlobalRows) {
                  LocalOrdinal col;
                  if (row + 1 < localNumRows)
                    col = row + 1;
                  else {
                    col = (comm_rank == 0) ? localNumRows : localNumRows + 1;
                  }
                  indices(k) = col;
                  values(k)  = -impl_ATS::one();
                }
              } else {
                // Set last entry of the rowptr.
                rowptr(row) = localNNZ;
              }
            });
      }

      // Construct the Tpetra::CrsMatrix
      {
        auto A_nonconst = rcp(new TpMatrix(map, colmap,
                                           rowptr, indices, values));
        A_nonconst->fillComplete();
        A = A_nonconst;
      }

      // Print out the matrix that we just constructed:
      // A->describe(*Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)), Teuchos::VERB_EXTREME);

      // Set up left-hand side and right-hand side of the linear system.
      X = Teuchos::rcp(new TpMultiVector(map, 1));
      {
        auto B_nonconst = Teuchos::rcp(new TpMultiVector(map, 1));
        B_nonconst->putScalar(1.);
        B = B_nonconst;
      }

      // Set up a near-nullspace for the linear system. (Only needed for multigrid)
      nullspace = Teuchos::rcp(new TpMultiVector(map, 1));
      nullspace->putScalar(1.);

      // Set up coordinates for the linear system. (Only needed for multigrid)
      coordinates = Teuchos::rcp(new TpMultiVector(map, 1));
      {
        GlobalOrdinal numGlobalRows = map->getGlobalNumElements();
        auto lclCoordinates         = coordinates->getLocalViewDevice(Tpetra::Access::OverwriteAll);
        auto lclMap                 = map->getLocalMap();
        Kokkos::parallel_for(
            "fill_coords",
            Kokkos::RangePolicy(0, map->getLocalNumElements()),
            KOKKOS_LAMBDA(const LocalOrdinal &i) {
              lclCoordinates(i, 0) = ((double)lclMap.getGlobalElement(i)) / ((double)numGlobalRows);
            });
      }
    }

    //
    // Build Stratimikos solver
    //
    RCP<Thyra::LinearOpWithSolveBase<Scalar>> solver;
    RCP<Thyra::PreconditionerBase<Scalar>> prec;

    // This is the Stratimikos main class (= factory of solver factory).
    Stratimikos::LinearSolverBuilder<Scalar> linearSolverBuilder;
    // Register MueLu as a Stratimikos preconditioner strategy.
    Stratimikos::enableMueLu<Scalar, LocalOrdinal, GlobalOrdinal, Node>(linearSolverBuilder);
    // add coordinates and nullspace to parameter list for multigrid preconditioner
    if (paramList->isSublist("Preconditioner Types") &&
        paramList->sublist("Preconditioner Types").isSublist("MueLu")) {
      ParameterList &userParamList = paramList->sublist("Preconditioner Types").sublist("MueLu").sublist("user data");
      if (!coordinates.is_null())
        userParamList.set("Coordinates", coordinates);
      if (!nullspace.is_null())
        userParamList.set("Nullspace", nullspace);
      if (!material.is_null())
        userParamList.set("Material", material);
      if (!blocknumber.is_null())
        userParamList.set("BlockNumber", blocknumber);
    }
    // Setup solver parameters using a Stratimikos parameter list.
    linearSolverBuilder.setParameterList(paramList);

    // Build a new "solver factory" according to the previously specified parameter list.
    auto solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
    auto precFactory   = solverFactory->getPreconditionerFactory();

    // Build Thyra solver
    if (!precFactory.is_null()) {
      prec   = Thyra::initializePrec(*precFactory, A);
      solver = Thyra::initializePreconditionedOp(*solverFactory, A, prec);
    } else {
      solver = Thyra::linearOpWithSolve(*solverFactory, A);
    }

    // Solve Ax = b.
    auto status = Thyra::solve(*solver, Thyra::NOTRANS, B, X);
    success     = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);

    for (int solveno = 1; solveno < numSolves; solveno++) {
      if (!precFactory.is_null())
        Thyra::initializePrec<Scalar>(*precFactory, A, prec);
      X->putScalar(0.);

      status  = Thyra::solve(*solver, Thyra::NOTRANS, B, X);
      success = success && (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);
    }

    // print timings
    if (printTimings) {
      if (use_stacked_timer) {
        stacked_timer->stop("Main");
        Teuchos::StackedTimer::OutputOptions options;
        options.output_fraction = options.output_histogram = options.output_minmax = true;
        stacked_timer->report(out, comm, options);
      } else {
        RCP<ParameterList> reportParams = rcp(new ParameterList);
        if (timingsFormat == "yaml") {
          reportParams->set("Report format", "YAML");  // "Table" or "YAML"
          reportParams->set("YAML style", "compact");  // "spacious" or "compact"
        }
        reportParams->set("How to merge timer sets", "Union");
        reportParams->set("alwaysWriteLocal", false);
        reportParams->set("writeGlobalStats", true);
        reportParams->set("writeZeroTimers", false);

        const std::string filter = "";

        std::ios_base::fmtflags ff(out.flags());
        if (timingsFormat == "table-fixed")
          out << std::fixed;
        else
          out << std::scientific;
        TimeMonitor::report(comm.ptr(), out, filter, reportParams);
        out << std::setiosflags(ff);
      }
    }

    TimeMonitor::clearCounters();
    out << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
