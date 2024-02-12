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

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"

#include "MatrixLoad.hpp"

class Measurement {
 public:
  virtual void measure() = 0;
};

class Reporter {
 public:
  Reporter(std::string n)
    : name(n)
    , comm(Teuchos::DefaultComm<int>::getComm())
    , out(Teuchos::rcpFromRef(std::cout))
    , timer(Teuchos::rcp(new Teuchos::StackedTimer(name.c_str()))) {
    out.setOutputToRootOnly(0);
    Teuchos::TimeMonitor::setStackedTimer(timer);
  }

  void record(Measurement& measurement) const {
    comm->barrier();
    measurement.measure();
  }

  void finalReport() {
    timer->stopBaseTimer();

    reportCommandLine();
    reportXML();

    out << "Run Complete!" << std::endl;
  }

 private:
  void reportCommandLine() {
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction  = true;
    options.output_histogram = true;
    options.output_minmax    = true;

    timer->report(out, comm, options);
  }

  void reportXML() {
    std::string xmlName = name + " " + std::to_string(comm->getSize()) + " ranks";
    auto xmlOut         = timer->reportWatchrXML(xmlName, comm);

    if (xmlOut.length()) {
      out << std::endl
          << "Created Watchr performance report " << xmlOut << std::endl;
    } else {
      out << std::endl
          << "Did not create Watchr performance report" << std::endl;
    }
  }

  const std::string name;
  Teuchos::RCP<const Teuchos::Comm<int>> comm;
  Teuchos::FancyOStream out;
  Teuchos::RCP<Teuchos::StackedTimer> timer;
};

class System {
 public:
  virtual void apply()             = 0;
  virtual std::string name() const = 0;
};

class SpMVMeasurement : public Measurement {
 public:
  SpMVMeasurement(System& sys, int n)
    : numRuns(n)
    , system(sys)
    , name(sys.name() + " SpMV") {}

  void measure() override {
    for (int i = 0; i < numRuns; i++) {
      Teuchos::TimeMonitor subTimer(*Teuchos::TimeMonitor::getNewTimer(name));
      system.apply();
    }
  }

 private:
  const int numRuns;
  System& system;
  const std::string name;
};

std::string systemName(const Xpetra::UnderlyingLib& lib) {
  if (lib == Xpetra::UseEpetra) return "Epetra";
  if (lib == Xpetra::UseTpetra) return "Tpetra";
  return "Unknown";
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class LinearSystem : public System {
  using MultiVector = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Matrix      = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

 public:
  LinearSystem(Xpetra::UnderlyingLib l, int n)
    : lib(l)
    , numVectors(n) {}

  void apply() override {
    A->apply(*X, *B);
  }

  std::string name() const override {
    return systemName(lib);
  }

  Xpetra::UnderlyingLib lib;
  int numVectors;

  Teuchos::RCP<Matrix> A;
  Teuchos::RCP<MultiVector> X;
  Teuchos::RCP<MultiVector> B;
};

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class SystemLoader {
  using Map                   = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using MultiVector           = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using CoordScalar           = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<CoordScalar, LocalOrdinal, GlobalOrdinal, Node>;

 public:
  SystemLoader(Teuchos::CommandLineProcessor& c)
    : clp(c)
    , comm(Teuchos::DefaultComm<int>::getComm())
    , galeriParameters(clp, 100, 100, 100, "Laplace2D")
    , xpetraParameters(clp) {}

  void fill(LinearSystem<Scalar, LocalOrdinal, GlobalOrdinal, Node>& system) {
    MatrixLoad<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
        comm, system.lib,
        binaryFormat, matrixFile, rhsFile, rowMapFile, colMapFile,
        domainMapFile, rangeMapFile, coordFile, coordMapFile, nullFile, materialFile,
        map, system.A,
        coordinates, nullspace, material,
        system.X, system.B, system.numVectors,
        galeriParameters, xpetraParameters,
        galeriStream);
  }

 private:
  Teuchos::CommandLineProcessor& clp;
  Teuchos::RCP<const Teuchos::Comm<int>> comm;
  Galeri::Xpetra::Parameters<GlobalOrdinal> galeriParameters;
  Xpetra::Parameters xpetraParameters;

  bool binaryFormat;
  std::string matrixFile    = "";
  std::string rhsFile       = "";
  std::string rowMapFile    = "";
  std::string colMapFile    = "";
  std::string domainMapFile = "";
  std::string rangeMapFile  = "";
  std::string coordFile     = "";
  std::string coordMapFile  = "";
  std::string nullFile      = "";
  std::string materialFile  = "";

  Teuchos::RCP<const Map> map;
  Teuchos::RCP<RealValuedMultiVector> coordinates;
  Teuchos::RCP<MultiVector> nullspace;
  Teuchos::RCP<MultiVector> material;
  std::ostringstream galeriStream;
};

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
int main_ETI(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib& lib, int argc, char* argv[]) {
  Reporter reporter("SpMV Performance " + systemName(lib));
  SystemLoader<Scalar, LocalOrdinal, GlobalOrdinal, Node> systemLoader(clp);

  int numRuns = 1;
  clp.setOption("num-runs", &numRuns, "number of times to run operation");

  int numVectors = 1;
  clp.setOption("multivector", &numVectors, "number of rhs to solve simultaneously");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  LinearSystem<Scalar, LocalOrdinal, GlobalOrdinal, Node> system(lib, numVectors);
  systemLoader.fill(system);

  SpMVMeasurement SpMV(system, numRuns);
  reporter.record(SpMV);

  reporter.finalReport();
  return 0;
}

#define MUELU_AUTOMATIC_TEST_ETI_NAME main_ETI
#include "MueLu_Test_ETI.hpp"

int main(int argc, char* argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
