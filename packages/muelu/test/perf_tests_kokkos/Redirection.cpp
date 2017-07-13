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
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <unistd.h>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra
#include <Xpetra_Matrix.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>

#include <MueLu_BaseClass.hpp>
#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif

template<class Matrix, class LO, class GO>
class RowFunctor {
  public:
    RowFunctor(Matrix localMatrix_) :
      localMatrix(localMatrix_)
    { }

    KOKKOS_INLINE_FUNCTION
    void operator()(const LO row, GO& r) const {
      auto rowView = localMatrix.row (row);
      auto length  = rowView.length;

      r += length;
    }

  private:
    Matrix localMatrix;
};

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::TimeMonitor;
  using Teuchos::ParameterList;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace3D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  int loops = 10; clp.setOption("l", &loops, "Number of loops");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }
  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  // Retrieve matrix parameters (they may have been changed on the command line)
  // [for instance, if we changed matrix type from 2D to 3D we need to update nz]
  ParameterList galeriList = galeriParameters.GetParameterList();

  // =========================================================================
  // Problem construction
  // =========================================================================
  std::string matrixType = galeriParameters.GetMatrixType();

  RCP<const Map> map;
  if (matrixType == "Laplace1D")
    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
  else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
           matrixType == "BigStar2D" || matrixType == "Elasticity2D")
    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
  else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D")
    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);

  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixType, map, galeriList);
  RCP<Matrix> A = Pr->BuildMatrix();

  LO numRows = A->getNodeNumRows();

  RCP<TimeMonitor> tm;

  // Loop 1
  {
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Loop #1: Xpetra")));

    GO validation = 0;

    ArrayView<const LocalOrdinal> indices;
    ArrayView<const Scalar> vals;
    for (int i = 0; i < loops; i++) {
      for (LocalOrdinal row = 0; row < numRows; row++) {

        A->getLocalRowView(row, indices, vals);

        validation += indices.size();
      }
    }
    std::cout << "validation = " << validation << std::endl;

    tm = Teuchos::null;
  }

  // Loop 2
  {
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Loop #2: Tpetra/Epetra")));

    GO validation = 0;

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<SC,LO,GO,NO> tCrsMatrix;
      RCP<const tCrsMatrix> tA = Utilities::Op2TpetraCrs(A);
      TEUCHOS_TEST_FOR_EXCEPTION(tA.is_null(), MueLu::Exceptions::RuntimeError,
                                 "A is not a Tpetra CrsMatrix");

      ArrayView<const LocalOrdinal> indices;
      ArrayView<const Scalar> vals;
      for (int i = 0; i < loops; i++)  {
        for (LocalOrdinal row = 0; row < numRows; row++) {
          tA->getLocalRowView(row, indices, vals);

          validation += indices.size();
        }
      }
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
                                 "Tpetra is not available");
#endif
    }
    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      typedef Epetra_CrsMatrix eCrsMatrix;
      RCP<const eCrsMatrix> eA = Utilities::Op2EpetraCrs(A);
      TEUCHOS_TEST_FOR_EXCEPTION(eA.is_null(), MueLu::Exceptions::RuntimeError,
                                 "A is not a Epetra CrsMatrix");

      for (int i = 0; i < loops; i++) {
        for (LocalOrdinal row = 0; row < numRows; row++) {
          int      numEntries;
          double * eValues;
          int    * eIndices;

          eA->ExtractMyRowView(row, numEntries, eValues, eIndices);

          validation += numEntries;
        }
      }
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
                                 "Epetra is not available");
#endif
    }
      std::cout << "validation = " << validation << std::endl;

    tm = Teuchos::null;
  }


  // Loop 3
  {
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Loop #3: Kokkos-serial")));

    typedef Kokkos::ArithTraits<SC> ATS;

    auto localMatrix = A->getLocalMatrix();

    GO validation = 0;
    Kokkos::View<bool*, typename NO::device_type> boundaryNodes("boundaryNodes", numRows);
    for (int i = 0; i < loops; i++) {
      for (LocalOrdinal row = 0; row < numRows; row++) {
        auto rowView = localMatrix.row (row);
        auto length  = rowView.length;

        validation += length;
      }
    }
    std::cout << "validation = " << validation << std::endl;

    tm = Teuchos::null;
  }

  // Loop 4
  {
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Loop #4: Kokkos-parallel (node)")));

    typedef Kokkos::ArithTraits<SC> ATS;

    auto localMatrix = A->getLocalMatrix();

    GO validation = 0;

    typedef Kokkos::RangePolicy<typename NO::execution_space> RangePolicy;

    Kokkos::View<bool*, typename NO::device_type> boundaryNodes("boundaryNodes", numRows);
    RowFunctor<decltype(localMatrix), LO, GO> functor(localMatrix);
    for (int i = 0; i < loops; i++)
      Kokkos::parallel_reduce("Utils::DetectDirichletRows", RangePolicy(0, numRows), functor, validation);
    std::cout << "validation = " << validation << std::endl;

    tm = Teuchos::null;
  }

  {
    const bool alwaysWriteLocal = false;
    const bool writeGlobalStats = true;
    const bool writeZeroTimers  = false;
    const bool ignoreZeroTimers = true;
    const std::string filter    = "";
    TimeMonitor::summarize(comm.ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                           writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
  }

  return EXIT_SUCCESS;
}

int main(int argc, char* argv[]) {
  bool success = false;
  bool verbose = true;

  Kokkos::initialize(argc, argv);

  try {
    const bool throwExceptions = false;

    Teuchos::CommandLineProcessor clp(throwExceptions);
    Xpetra::Parameters xpetraParameters(clp);

    std::string node = "";  clp.setOption("node", &node, "node type (serial | openmp | cuda)");

    clp.recogniseAllOptions(false);
    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR:               return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      return main_<double,int,int,Xpetra::EpetraNode>(clp, argc, argv);
#else
      throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      if (node == "") {
        typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, argc, argv);
#else
#  if defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        return main_<double,int,int,Node> (clp, argc, argv);
#  elif defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        return main_<double,int,long,Node>(clp, argc, argv);
#  elif defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        return main_<double,int,long long,Node>(clp, argc, argv);
#  else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#  endif
#endif
      } else if (node == "serial") {
#ifdef KOKKOS_HAVE_SERIAL
        typedef Kokkos::Compat::KokkosSerialWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, argc, argv);
#  else
#    if   defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        return main_<double,int,int,Node> (clp, argc, argv);
#    elif defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        return main_<double,int,long,Node>(clp, argc, argv);
#    elif defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        return main_<double,int,long long,Node>(clp, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("Serial node type is disabled");
#endif
      } else if (node == "openmp") {
#ifdef KOKKOS_HAVE_OPENMP
        typedef Kokkos::Compat::KokkosOpenMPWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, argc, argv);
#  else
#    if   defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        return main_<double,int,int,Node> (clp, argc, argv);
#    elif defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        return main_<double,int,long,Node>(clp, argc, argv);
#    elif defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        return main_<double,int,long long,Node>(clp, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("OpenMP node type is disabled");
#endif
      } else if (node == "cuda") {
#ifdef KOKKOS_HAVE_CUDA
        typedef Kokkos::Compat::KokkosCudaWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, argc, argv);
#  else
#    if   defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
        return main_<double,int,int,Node> (clp, argc, argv);
#    elif defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
        return main_<double,int,long,Node>(clp, argc, argv);
#    elif defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
        return main_<double,int,long long,Node>(clp, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("CUDA node type is disabled");
#endif
      } else {
        throw MueLu::Exceptions::RuntimeError("Unrecognized node type");
      }
#else
      throw MueLu::Exceptions::RuntimeError("Tpetra is not available");
#endif
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  Kokkos::finalize();

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
