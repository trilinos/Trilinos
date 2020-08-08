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

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_Matrix.hpp"

#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraParameters.hpp"

#include "MatrixLoad.hpp"

#include "MueLu.hpp"

class Measurement {
public:
  virtual void measure() = 0;
};

class Reporter {
public:
  Reporter(std::string n)
    : name(n),
      comm(Teuchos::DefaultComm<int>::getComm()),
      out(Teuchos::rcpFromRef(std::cout)),
      timer(Teuchos::rcp(new Teuchos::StackedTimer(name.c_str())))
  {
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
    options.output_fraction = true;
    options.output_histogram = true;
    options.output_minmax = true;

    timer->report(out, comm, options);
  }

  void reportXML() {
    std::string xmlName = name + " " + std::to_string(comm->getSize()) + " ranks";
    auto xmlOut = timer->reportWatchrXML(xmlName, comm);

    if(xmlOut.length()) {
      out << std::endl << "Created Watchr performance report " << xmlOut << std::endl;
    }
    else {
      out << std::endl << "Did not create Watchr performance report" << std::endl;
    }
  }

  const std::string name;
  Teuchos::RCP<const Teuchos::Comm<int>> comm;
  Teuchos::FancyOStream out;
  Teuchos::RCP<Teuchos::StackedTimer> timer;
};

class System {
public:
  virtual void apply() = 0;
  virtual std::string name() const = 0;
};

class SpMVMeasurement : public Measurement {
public:
  SpMVMeasurement(System& sys, int n)
    : numRuns(n),
      system(sys),
      name(sys.name() + " SpMV")
  { }

  void measure() override {
    Teuchos::TimeMonitor subTimer(*Teuchos::TimeMonitor::getNewTimer(name));

    for (int i=0; i<numRuns; i++) {
      system.apply();
    }
  }

private:
  const int numRuns;
  System& system;
  const std::string name;
};

template<typename S, typename LO, typename GO, typename N>
class LinearSystem : public System {
public:
  using Scalar = S;
  using LocalOrdinal = LO;
  using GlobalOrdinal = GO;
  using Node = N;

  using MultiVector = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  LinearSystem(Xpetra::UnderlyingLib l)
    : lib(l)
  { }

  void apply() override {
    A->apply(*X, *B);
  }

  std::string name() const override {
    if (lib == Xpetra::UseEpetra) return "Epetra";
    if (lib == Xpetra::UseTpetra) return "Tpetra";
    return "Unknown";
  }

  Xpetra::UnderlyingLib lib;

  Teuchos::RCP<Matrix> A;
  Teuchos::RCP<MultiVector> X;
  Teuchos::RCP<MultiVector> B;
};

class DirectSystemLoader {
public:
  template<typename LinearSystem>
  void fill(LinearSystem& data) {
    using Scalar = typename LinearSystem::Scalar;
    using LocalOrdinal = typename LinearSystem::LocalOrdinal;
    using GlobalOrdinal = typename LinearSystem::GlobalOrdinal;
    using Node = typename LinearSystem::Node;

    using MultiVector = typename LinearSystem::MultiVector;
    using Matrix = typename LinearSystem::Matrix;

    using Map = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
    using MultiVectorFactory = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using Utilities = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using CrsMatrixWrap = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr;
    Teuchos::RCP<Map> map;

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", 100);
    galeriList.set("ny", 100);
    map = Galeri::Xpetra::CreateMap<int, int, Node>(data.lib, "Cartesian2D", comm, galeriList);

    Pr = Galeri::Xpetra::BuildProblem<double, int, int, Map, CrsMatrixWrap, MultiVector>("Laplace2D", map, galeriList);
    data.A = Pr->BuildMatrix();

    int numVecs = 4;
    Utilities::SetRandomSeed(*comm);
    data.X = MultiVectorFactory::Build(map, numVecs);
    data.X->randomize();
    data.B = MultiVectorFactory::Build(map, numVecs);
    data.B->randomize();
  }
};

class ParameterSystemLoader {
public:
  ParameterSystemLoader(Teuchos::CommandLineProcessor& c)
    : clp(c)
  { }

  template<typename LinearSystem>
  void fill(LinearSystem& data) {
    using Scalar = typename LinearSystem::Scalar;
    using LocalOrdinal = typename LinearSystem::LocalOrdinal;
    using GlobalOrdinal = typename LinearSystem::GlobalOrdinal;
    using Node = typename LinearSystem::Node;

    using Map = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
    using MultiVector = typename LinearSystem::MultiVector;
    using CoordScalar = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
    using RealValuedMultiVector = Xpetra::MultiVector<CoordScalar, LocalOrdinal, GlobalOrdinal, Node>;

    bool binaryFormat = "";
    std::string matrixFile = "";
    std::string rhsFile = "";
    std::string rowMapFile = "";
    std::string colMapFile = "";
    std::string domainMapFile = "";
    std::string rangeMapFile = "";
    std::string coordFile = "";
    std::string nullFile = "";
    std::string materialFile = "";
    Teuchos::RCP<RealValuedMultiVector> coordinates;
    Teuchos::RCP<MultiVector> nullspace, material;
    std::ostringstream galeriStream;

    Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<const Map> map;
    int numVectors = 4;
    GlobalOrdinal nx = 100, ny = 100, nz = 100;
    Galeri::Xpetra::Parameters<GlobalOrdinal> galeriParameters(clp, nx, ny, nz, "Laplace2D");
    Xpetra::Parameters xpetraParameters(clp);

    MatrixLoad<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
        comm, data.lib,
        binaryFormat, matrixFile, rhsFile, rowMapFile, colMapFile,
        domainMapFile, rangeMapFile, coordFile, nullFile, materialFile,
        map, data.A,
        coordinates, nullspace, material,
        data.X, data.B, numVectors,
        galeriParameters, xpetraParameters,
        galeriStream);
  }

private:
  Teuchos::CommandLineProcessor& clp;
};

using EpetraNode = Xpetra::EpetraNode;
using EpetraSystem = LinearSystem<double, int, int, EpetraNode>;

using TpetraNode = KokkosClassic::DefaultNode::DefaultNodeType;
using TpetraSystem = LinearSystem<double, int, int, TpetraNode>;

int main(int argc, char* argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  Teuchos::CommandLineProcessor clp(false);
  clp.recogniseAllOptions(false);
  switch (clp.parse(argc, argv, NULL)) {
    case Teuchos::CommandLineProcessor::PARSE_ERROR:                return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:         break;
  }

  Reporter reporter("SpMV Performance");
  DirectSystemLoader systemLoader;
  //ParameterSystemLoader systemLoader(clp);

  int numRuns = 4;

  EpetraSystem epetraSystem(Xpetra::UseEpetra);
  systemLoader.fill(epetraSystem);
  SpMVMeasurement epetraSpMV(epetraSystem, numRuns);
  reporter.record(epetraSpMV);

  TpetraSystem tpetraSystem(Xpetra::UseTpetra);
  systemLoader.fill(tpetraSystem);
  SpMVMeasurement tpetraSpMV(tpetraSystem, numRuns);
  reporter.record(tpetraSpMV);

  reporter.finalReport();
  return 0;
}
