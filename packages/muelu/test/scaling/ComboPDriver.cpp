// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Teuchos
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

// MueLu
#include <MueLu_Exceptions.hpp>
#include <MueLu_MultiPhys.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_TestHelpers.hpp>
#include "MueLu_SemiCoarsenPFactory.hpp"

// Xpetra
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosXpetraAdapter.hpp>
#endif

#define FOURBYFOUR

namespace {

std::string defaultSuffix() {
#ifdef BIGVERSION
  return "big";
#else
  return "";
#endif
}

int defaultBlockCount() {
#ifdef FOURBYFOUR
  return 4;
#else
  return 2;
#endif
}

std::string defaultMultiPhysFile(const std::string& suffix) {
#ifdef FOURBYFOUR
  return "multiPhys4x4" + suffix + ".mat";
#else
  return "multiPhys2x2" + suffix + ".mat";
#endif
}

std::string trim(std::string s) {
  auto notSpace = [](unsigned char ch) { return !std::isspace(ch); };

  s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
  s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());

  return s;
}

std::string normalizeOptionalFilename(std::string value) {
  value = trim(std::move(value));

  if (value == "none" || value == "null" ||
      value == "NONE" || value == "NULL") {
    return "";
  }

  return value;
}

std::vector<std::string> splitCommaSeparated(const std::string& input) {
  if (input.empty()) {
    return {};
  }

  std::vector<std::string> result;
  std::size_t start = 0;

  while (true) {
    const auto pos = input.find(',', start);
    if (pos == std::string::npos) {
      result.emplace_back(trim(input.substr(start)));
      break;
    }

    result.emplace_back(trim(input.substr(start, pos - start)));
    start = pos + 1;
  }

  return result;
}

int readDimensionFromMatFile(const std::string& filename) {
  std::ifstream in(filename);
  if (!in.is_open()) {
    throw std::runtime_error("File not found: " + filename);
  }

  std::string line;
  while (std::getline(in, line)) {
    line = trim(line);

    if (line.empty()) {
      continue;
    }

    if (line[0] == '%') {
      continue;
    }

    std::istringstream iss(line);
    int nRows = 0;
    int nCols = 0;
    int nnz   = 0;

    if (!(iss >> nRows >> nCols)) {
      throw std::runtime_error("Failed to parse matrix dimension line from file: " + filename);
    }

    (void)nnz;
    return nRows;
  }

  throw std::runtime_error("Failed to locate matrix dimension line in file: " + filename);
}

struct BlockSpec {
  std::string auxFile;
  std::string coordsFile;
  std::string nullFile;
  std::string materialFile;
};

struct DriverOptions {
  std::string suffix        = defaultSuffix();
  std::string xmlFile       = "comboP.xml";
  std::string multiPhysFile = defaultMultiPhysFile(defaultSuffix());
  int nBlks                 = defaultBlockCount();
  int expectedIters         = 50;
  std::vector<BlockSpec> blocks;
  double tolerance  = 1e-6;
  int maxIterations = 100;
};

DriverOptions makeDefaultOptions() {
  DriverOptions opts;
  const auto s = opts.suffix;

#ifdef FOURBYFOUR
  opts.blocks = {
      {"aux1" + s + ".mat", "coords1" + s + ".mat", "null1" + s + ".mat", ""},
      {"aux2" + s + ".mat", "coords2" + s + ".mat", "null2" + s + ".mat", ""},
      {"aux1" + s + ".mat", "coords1" + s + ".mat", "null1" + s + ".mat", ""},
      {"aux2" + s + ".mat", "coords2" + s + ".mat", "null2" + s + ".mat", ""}};
#else
  opts.blocks = {
      {"aux1" + s + ".mat", "coords1" + s + ".mat", "null1" + s + ".mat", ""},
      {"aux2" + s + ".mat", "coords2" + s + ".mat", "null2" + s + ".mat", ""}};
#endif

  return opts;
}

void parseCommandLine(Teuchos::CommandLineProcessor& clp, int argc, char* argv[], DriverOptions& opts) {
  std::string auxFilesArg;
  std::string coordFilesArg;
  std::string nullFilesArg;
  std::string materialFilesArg;

  clp.setOption("xml", &opts.xmlFile, "XML parameter file");
  clp.setOption("matrix", &opts.multiPhysFile, "Multiphysics matrix file");
  clp.setOption("nblks", &opts.nBlks, "Number of multiphysics blocks");
  clp.setOption("its", &opts.expectedIters, "Expected maximum iteration count for success");
  clp.setOption("tol", &opts.tolerance, "Target residual reduction for convergence");
  clp.setOption("max-its", &opts.maxIterations, "Maximum number of iterations");
  clp.setOption("auxFiles", &auxFilesArg, "Comma-separated auxiliary matrix files");
  clp.setOption("coordFiles", &coordFilesArg, "Comma-separated coordinate multivector files (empty entry allowed)");
  clp.setOption("nullFiles", &nullFilesArg, "Comma-separated nullspace multivector files (empty entry allowed)");
  clp.setOption("materialFiles", &materialFilesArg, "Comma-separated material multivector files (empty entry allowed)");

  const auto parseResult = clp.parse(argc, argv);
  if (parseResult != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    throw std::runtime_error("Command line parsing failed");
  }

  if (opts.nBlks <= 0) {
    throw std::runtime_error("nBlks must be positive");
  }

  std::vector<BlockSpec> blocks;
  if (static_cast<int>(opts.blocks.size()) == opts.nBlks) {
    blocks = opts.blocks;
  } else {
    blocks.assign(opts.nBlks, BlockSpec{});
  }

  const auto assignList = [&](const std::string& arg,
                              const std::string& fieldName,
                              auto setter) {
    if (trim(arg).empty()) {
      return;
    }

    auto values = splitCommaSeparated(arg);
    if (static_cast<int>(values.size()) != opts.nBlks) {
      throw std::runtime_error("Expected " + fieldName + " list to have length " +
                               std::to_string(opts.nBlks) + ", but got " +
                               std::to_string(values.size()));
    }

    for (int i = 0; i < opts.nBlks; ++i) {
      setter(blocks[i], normalizeOptionalFilename(values[i]));
    }
  };

  assignList(auxFilesArg, "auxFiles", [](BlockSpec& b, const std::string& v) { b.auxFile = v; });
  assignList(coordFilesArg, "coordFiles", [](BlockSpec& b, const std::string& v) { b.coordsFile = v; });
  assignList(nullFilesArg, "nullFiles", [](BlockSpec& b, const std::string& v) { b.nullFile = v; });
  assignList(materialFilesArg, "materialFiles", [](BlockSpec& b, const std::string& v) { b.materialFile = v; });

  for (int i = 0; i < opts.nBlks; ++i) {
    if (blocks[i].auxFile.empty()) {
      throw std::runtime_error("Auxiliary matrix file for block " + std::to_string(i) +
                               " must not be empty");
    }
  }

  opts.blocks = std::move(blocks);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
readOptionalMultiVector(
    const std::string& filename,
    const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map) {
  if (filename.empty()) {
    return Teuchos::null;
  }

  return Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(filename, map);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,
                                 LocalOrdinal, GlobalOrdinal, Node>>
readOptionalCoordinateMultiVector(
    const std::string& filename,
    const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map) {
  using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  if (filename.empty()) {
    return Teuchos::null;
  }

  return Xpetra::IO<magnitude_type, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(filename, map);
}

}  // namespace

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>

  using SC = Scalar;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using NO = Node;

  using STS                   = Teuchos::ScalarTraits<SC>;
  using Magnitude             = typename STS::magnitudeType;
  using Coordinate            = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<Coordinate, LO, GO, NO>;

  bool success = true;
  bool verbose = true;

  try {
    const auto comm = Teuchos::DefaultComm<int>::getComm();

    auto opts = makeDefaultOptions();
    parseCommandLine(clp, argc, argv, opts);

    std::vector<GO> blockDims(opts.nBlks, 0);
    for (int i = 0; i < opts.nBlks; ++i) {
      blockDims[i] = static_cast<GO>(readDimensionFromMatFile(opts.blocks[i].auxFile));
    }

    const GO expectedMultiPhysDim =
        std::accumulate(blockDims.begin(), blockDims.end(), static_cast<GO>(0));

    const GO actualMultiPhysDim =
        static_cast<GO>(readDimensionFromMatFile(opts.multiPhysFile));

    if (actualMultiPhysDim != expectedMultiPhysDim) {
      if (comm->getRank() == 0) {
        std::cerr << "ERROR: multiphysics matrix dimension is "
                  << actualMultiPhysDim
                  << ", expected " << expectedMultiPhysDim
                  << std::endl;
      }
      return EXIT_FAILURE;
    }

    std::vector<Teuchos::RCP<const Map>> nodeMaps(opts.nBlks);
    std::vector<Teuchos::RCP<const Map>> dofMaps(opts.nBlks);
    std::vector<GO> blockOffsets(opts.nBlks, 0);

    for (int i = 1; i < opts.nBlks; ++i) {
      blockOffsets[i] = blockOffsets[i - 1] + blockDims[i - 1];
    }

    for (int i = 0; i < opts.nBlks; ++i) {
      const GO blockDim = blockDims[i];

      GO nodalCount = blockDim;

      if (!opts.blocks[i].coordsFile.empty()) {
        const GO probeSize  = blockDim;
        const auto probeMap = MapFactory::Build(lib, probeSize, 0, comm);
        auto coordsProbe =
            Xpetra::IO<Magnitude, LO, GO, NO>::ReadMultiVector(opts.blocks[i].coordsFile, probeMap);

        nodalCount = static_cast<GO>(coordsProbe->getMap()->getGlobalNumElements());

        if (nodalCount <= 0) {
          throw std::runtime_error("Invalid nodal count inferred from coordinates for block " +
                                   std::to_string(i));
        }

        if (blockDim % nodalCount != 0) {
          throw std::runtime_error(
              "Cannot infer integer dofs-per-node for block " + std::to_string(i) +
              ": block dimension " + std::to_string(blockDim) +
              " is not divisible by coordinate map size " + std::to_string(nodalCount));
        }
      }

      nodeMaps[i] = MapFactory::Build(lib, nodalCount, 0, comm);

      const GO dofsPerNode = blockDim / nodalCount;
      if (dofsPerNode == 1) {
        dofMaps[i] = nodeMaps[i];
      } else {
        dofMaps[i] = Xpetra::MapFactory<LO, GO, NO>::Build(
            nodeMaps[i], static_cast<size_t>(dofsPerNode));
      }
    }

    Teuchos::Array<GO> gidsCombo;
    {
      std::size_t localTotal = 0;
      for (int i = 0; i < opts.nBlks; ++i) {
        localTotal += dofMaps[i]->getLocalNumElements();
      }

      gidsCombo.resize(localTotal);

      std::size_t cursor = 0;
      for (int i = 0; i < opts.nBlks; ++i) {
        const auto localGids = dofMaps[i]->getLocalElementList();
        for (int j = 0; j < static_cast<int>(localGids.size()); ++j) {
          gidsCombo[cursor++] = localGids[j] + blockOffsets[i];
        }
      }
    }

    const auto mapMultiPhysA =
        Xpetra::MapFactory<LO, GO, NO>::Build(
            lib,
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            gidsCombo,
            0,
            comm);

    auto multiPhysA = Xpetra::IO<SC, LO, GO, NO>::Read(opts.multiPhysFile, mapMultiPhysA);

    Teuchos::ArrayRCP<Teuchos::RCP<Matrix>> arrayOfAuxMatrices(opts.nBlks);
    Teuchos::ArrayRCP<Teuchos::RCP<RealValuedMultiVector>> arrayOfCoords(opts.nBlks);
    Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfNullspaces(opts.nBlks);
    Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfMaterials(opts.nBlks);

    for (int i = 0; i < opts.nBlks; ++i) {
      arrayOfAuxMatrices[i] = Xpetra::IO<SC, LO, GO, NO>::Read(opts.blocks[i].auxFile, dofMaps[i]);

      arrayOfCoords[i] = readOptionalCoordinateMultiVector<SC, LO, GO, NO>(
          opts.blocks[i].coordsFile, nodeMaps[i]);

      arrayOfNullspaces[i] = readOptionalMultiVector<SC, LO, GO, NO>(
          opts.blocks[i].nullFile, dofMaps[i]);

      arrayOfMaterials[i] = readOptionalMultiVector<SC, LO, GO, NO>(
          opts.blocks[i].materialFile, dofMaps[i]);
    }

    Teuchos::ParameterList comboList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        opts.xmlFile,
        Teuchos::Ptr<Teuchos::ParameterList>(&comboList),
        *comm);

    auto preconditioner = Teuchos::rcp(
        new MueLu::MultiPhys<SC, LO, GO, NO>(
            multiPhysA,
            arrayOfAuxMatrices,
            arrayOfNullspaces,
            arrayOfCoords,
            opts.nBlks,
            comboList,
            true,
            arrayOfMaterials,
            true));

    auto globalTimeMonitor =
        Teuchos::rcp(new Teuchos::TimeMonitor(
            *Teuchos::TimeMonitor::getNewTimer("Timings: Global Time")));

#ifdef HAVE_MUELU_BELOS
    const SC zero = STS::zero();
    const SC one  = STS::one();

    auto X = VectorFactory::Build(multiPhysA->getRowMap());
    auto B = VectorFactory::Build(multiPhysA->getRowMap());

    {
      Utilities::SetRandomSeed(*comm);
      X->randomize();
      multiPhysA->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

      Teuchos::Array<typename STS::magnitudeType> norms(1);
      B->norm2(norms);
      B->scale(one / norms[0]);
      X->putScalar(zero);
    }

    using MV = MultiVector;
    using OP = Belos::OperatorT<MV>;

    auto belosOp     = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(multiPhysA));
    auto belosPrecOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(preconditioner));

    auto belosProblem =
        Teuchos::rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
    belosProblem->setRightPrec(belosPrecOp);

    const bool set = belosProblem->setProblem();
    if (!set) {
      if (comm->getRank() == 0) {
        std::cerr << "ERROR: Belos::LinearProblem failed to set up correctly!" << std::endl;
      }
      return EXIT_FAILURE;
    }

    auto belosList = Teuchos::parameterList();
    belosList->set("Maximum Iterations", opts.maxIterations);
    belosList->set("Convergence Tolerance", opts.tolerance);
    belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList->set("Output Frequency", 1);
    belosList->set("Output Style", Belos::Brief);
    belosList->set("Implicit Residual Scaling", "None");

    Belos::SolverFactory<SC, MV, OP> solverFactory;
    auto solver = solverFactory.create("gmres", belosList);
    solver->setProblem(belosProblem);

    const auto retStatus = solver->solve();
    const int iters      = solver->getNumIters();

    success = (iters <= opts.expectedIters && retStatus == Belos::Converged);

    if (comm->getRank() == 0) {
      if (success) {
        std::cout << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
      } else if (retStatus == Belos::Converged) {
        std::cout << "FAILURE! Belos converged, but required "
                  << iters << " iterations which exceeds the allowed "
                  << opts.expectedIters << "." << std::endl;
      } else {
        std::cout << "FAILURE! Belos did not converge." << std::endl;
      }
    }
#else
    if (comm->getRank() == 0) {
      std::cout << "Belos not enabled; preconditioner was constructed, but no solve was performed." << std::endl;
    }
#endif

    globalTimeMonitor = Teuchos::null;
    Teuchos::TimeMonitor::summarize();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char* argv[]) {
  return Automatic_Test_ETI(argc, argv);
}