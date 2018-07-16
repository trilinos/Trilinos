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
#include <unistd.h>
#include <iostream>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include "MueLu.hpp"

#ifdef HAVE_MUELU_TPETRA
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp>
#endif

#ifdef HAVE_MUELU_STRATIMIKOS
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_MueLuHelpers.hpp>
#endif

#ifdef HAVE_MUELU_TEKO
#include <Teko_EpetraBlockPreconditioner.hpp>
#include <Teko_InverseFactory.hpp>
#include <Teko_InverseLibrary.hpp>
#include <Teko_InvLSCStrategy.hpp>
#include <Teko_LSCPreconditionerFactory.hpp>
#include <Teko_SIMPLEPreconditionerFactory.hpp>
#include <Teko_StridedEpetraOperator.hpp>
#include <Teko_Utilities.hpp>
#endif

#include <Thyra_DefaultPreconditioner.hpp>
#include <Thyra_DefaultScaledAdjointLinearOp.hpp>
#include <Thyra_Ifpack2PreconditionerFactory.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include <Thyra_MueLuTpetraQ2Q1PreconditionerFactory.hpp>
#include <Thyra_MultiVectorStdOps.hpp>
#include <Thyra_PreconditionerFactoryHelpers.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_TpetraVectorSpace.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorStdOps.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_Utilities.hpp"

// FIXME
#include "MueLu_UseDefaultTypes.hpp"

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
ReadBinary(const std::string& fileName, const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  typedef Scalar        SC;
  typedef LocalOrdinal  LO;
  typedef GlobalOrdinal GO;
  typedef Node          NO;
  TEUCHOS_TEST_FOR_EXCEPTION(comm->getSize() != 1, MueLu::Exceptions::RuntimeError, "Serial read only");

  std::ifstream ifs(fileName.c_str(), std::ios::binary);
  TEUCHOS_TEST_FOR_EXCEPTION(!ifs.good(), MueLu::Exceptions::RuntimeError, "Can not read \"" << fileName << "\"");

  int m, n, nnz;
  ifs.read(reinterpret_cast<char*>(&m),   sizeof(m));
  ifs.read(reinterpret_cast<char*>(&n),   sizeof(n));
  ifs.read(reinterpret_cast<char*>(&nnz), sizeof(nnz));

  typedef Tpetra::Map      <LO,GO,NO>       tMap;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NO>    tCrsMatrix;

  GO indexBase = 0;
  Teuchos::RCP<tMap>       rowMap = rcp(new tMap(m, indexBase, comm)), rangeMap  = rowMap;
  Teuchos::RCP<tMap>       colMap = rcp(new tMap(n, indexBase, comm)), domainMap = colMap;;
  Teuchos::RCP<tCrsMatrix> A      = rcp(new tCrsMatrix(rowMap, colMap, 9));

  TEUCHOS_TEST_FOR_EXCEPTION(sizeof(int) != sizeof(GO), MueLu::Exceptions::RuntimeError, "Incompatible sizes");

  Teuchos::Array<GO> inds;
  Teuchos::Array<SC> vals;
  for (int i = 0; i < m; i++) {
    int row, rownnz;
    ifs.read(reinterpret_cast<char*>(&row),    sizeof(row));
    ifs.read(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));
    inds.resize(rownnz);
    vals.resize(rownnz);
    for (int j = 0; j < rownnz; j++) {
      int index;
      ifs.read(reinterpret_cast<char*>(&index), sizeof(index));
      inds[j] = Teuchos::as<GO>(index);
    }
    for (int j = 0; j < rownnz; j++) {
      double value;
      ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
      vals[j] = Teuchos::as<SC>(value);
    }
    A->insertGlobalValues(row, inds, vals);
  }

  A->fillComplete(domainMap, rangeMap);

  return A;
}

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::null;
  using Teuchos::as;
  using Teuchos::TimeMonitor;
  using Tpetra::MatrixMarket::Reader;
  using Thyra::tpetraVectorSpace;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = true;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    TEUCHOS_TEST_FOR_EXCEPTION(comm->getSize() > 1, MueLu::Exceptions::RuntimeError,
                               "For now, Q2Q1 only works in serial. We are working on the parallel implementation.");

    // =========================================================================
    // Parameter initialization
    // =========================================================================
    Teuchos::CommandLineProcessor clp(false);

    Xpetra::Parameters xpetraParameters(clp);

    // configure problem
    std::string prefix = "./Q2Q1_9x9_";      clp.setOption("prefix",     &prefix,        "prefix for data files");
    std::string rhs    = "";                    clp.setOption("rhs",        &rhs,           "rhs");

    // configure run
    std::string xmlFileName  = "driver.xml";    clp.setOption("xml",        &xmlFileName,   "read parameters from a file [default = 'driver.xml']");
    double      tol          = 1e-8;            clp.setOption("tol",        &tol,           "solver convergence tolerance");
    std::string type         = "unstructured";  clp.setOption("type",       &type,          "structured/unstructured");
    int         use9ptPatA   = 1;               clp.setOption("use9pt",     &use9ptPatA,    "use 9-point stencil matrix for velocity prolongator construction");
    int         useFilters   = 1;               clp.setOption("usefilters", &useFilters,    "use filters on A and BB^T");

    double      tau_1        = 0.06;            clp.setOption("tau_1",      &tau_1,         "tau_1 parameter from paper (for matrix filtering)");
    double      tau_2        = sqrt(0.0015);    clp.setOption("tau_2",      &tau_2,         "tau_2 parameter from paper (for mid point detection)");

    int         binary       = 0;               clp.setOption("binary",     &binary,        "read matrix in binary format");

    // configure misc
    int         printTimings = 0;               clp.setOption("timings",    &printTimings,  "print timings to screen");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    typedef Tpetra::CrsMatrix<SC,LO,GO,NO>        tCrsMatrix;
    typedef Tpetra::Operator<SC,LO,GO,NO>         tOperator;
    typedef Tpetra::MultiVector<SC,LO,GO,NO>      tMultiVector;
    typedef Tpetra::Map<LO,GO,NO>                 tMap;
    typedef Thyra::TpetraVectorSpace<SC,LO,GO,NO> THTP_Vs;
    RCP<NO> node = Tpetra::DefaultPlatform::getDefaultPlatform().getNode();

    // Read data from files
    RCP<tOperator> A11, A119Pt, A21, A12;
    RCP<tMultiVector> Vcoords, Pcoords;
    ArrayRCP<LO> p2vMap;

    std::string filename;
    try {
      if (!binary) {
        filename = prefix + "A.mm";
        A11     = Reader<tCrsMatrix>::readSparseFile(filename.c_str(), comm, node);
        A119Pt  = A11;
        if (use9ptPatA) {
          filename = prefix + "AForPat.mm";
          A119Pt = Reader<tCrsMatrix>::readSparseFile(filename.c_str(), comm, node);
        }
      } else {
        filename = prefix + "A.dat";
        A11     = ReadBinary<SC,LO,GO,NO>(filename.c_str(), comm);
        A119Pt  = A11;
        if (use9ptPatA) {
          filename = prefix + "AForPat.dat";
          A119Pt = ReadBinary<SC,LO,GO,NO>(filename.c_str(), comm);
        }
      }
      filename = prefix + "B.mm";
      A21 = Reader<tCrsMatrix>::readSparseFile(filename.c_str(), comm, node);
      filename = prefix + "Bt.mm";
      A12 = Reader<tCrsMatrix>::readSparseFile(filename.c_str(), comm, node);

      RCP<const tMap> cmap1 = A11->getDomainMap(), cmap2 = A12->getDomainMap();
      filename = prefix + "VelCoords.mm";
      Vcoords = Reader<tCrsMatrix>::readDenseFile(filename.c_str(),  comm, node, cmap1);
      filename = prefix + "PresCoords.mm";
      Pcoords = Reader<tCrsMatrix>::readDenseFile(filename.c_str(), comm, node, cmap2);

      // For now, we assume that p2v maps local pressure DOF to a local x-velocity DOF
      filename = prefix + "p2vMap.mm";
      ArrayRCP<const SC> slop = Xpetra::IO<SC,LO,GO,NO>::ReadMultiVector(filename.c_str(),
                                                        Xpetra::toXpetra(A21->getRangeMap()))->getData(0);
      p2vMap.resize(slop.size());
      for (int i = 0; i < slop.size(); i++)
        p2vMap[i] = as<LO>(slop[i]);
    } catch (...) {
      throw MueLu::Exceptions::RuntimeError("Error reading file: \"" + filename + "\"");
    }

    // Convert matrices to Teko/Thyra operators
    RCP<const THTP_Vs> domain11 = tpetraVectorSpace<SC>(A11->getDomainMap());
    RCP<const THTP_Vs> domain12 = tpetraVectorSpace<SC>(A12->getDomainMap());
    RCP<const THTP_Vs> domain21 = tpetraVectorSpace<SC>(A21->getDomainMap());

    RCP<const THTP_Vs> range11  = tpetraVectorSpace<SC>(A11->getRangeMap());
    RCP<const THTP_Vs> range12  = tpetraVectorSpace<SC>(A12->getRangeMap());
    RCP<const THTP_Vs> range21  = tpetraVectorSpace<SC>(A21->getRangeMap());

    Teko::LinearOp thA11        = Thyra::tpetraLinearOp<double>(range11, domain11, A11);
    Teko::LinearOp thA12        = Thyra::tpetraLinearOp<double>(range12, domain12, A12);
    Teko::LinearOp thA21        = Thyra::tpetraLinearOp<double>(range21, domain21, A21);
    Teko::LinearOp thA11_9Pt    = Thyra::tpetraLinearOp<double>(range11, domain11, A119Pt);

    // Bang together the parameter list. Right now, all the MueLu details is
    // hardwired in the MueLu-TpetraQ2Q1 classes. We probably want to switch
    // things so that several of these hardwired features could be modified
    // via parameter lists.

    RCP<ParameterList> stratimikosList = rcp(new ParameterList);
    stratimikosList->set("Linear Solver Type",  "Belos");
    stratimikosList->set("Preconditioner Type", "MueLu-TpetraQ2Q1");

    ParameterList& BelosList = stratimikosList->sublist("Linear Solver Types").sublist("Belos");
    BelosList.set("Solver Type", "Block GMRES"); // FIXME: should it be "Pseudo Block GMRES"?
    BelosList.sublist("VerboseObject").set("Verbosity Level", "low"); // this is needed, as otherwise Stratimikos ignores Belos output

    ParameterList& GmresDetails = BelosList.sublist("Solver Types").sublist("Block GMRES");
    GmresDetails.set("Maximum Iterations",      100);
    GmresDetails.set("Convergence Tolerance",   tol);
    GmresDetails.set("Verbosity",               Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    GmresDetails.set("Output Frequency",        1);
    GmresDetails.set("Output Style",            Belos::Brief);

    ParameterList& Q2Q1List = stratimikosList->sublist("Preconditioner Types").sublist("MueLu-TpetraQ2Q1");
    Q2Q1List.set("useFilters",  useFilters);
    Q2Q1List.set("tau_1",       tau_1);
    Q2Q1List.set("tau_2",       tau_2);
    Q2Q1List.set("Velcoords",   Vcoords);
    Q2Q1List.set("Prescoords",  Pcoords);
    Q2Q1List.set("p2vMap",      p2vMap);
    Q2Q1List.set("A11",         thA11);
    Q2Q1List.set("A12",         thA12);
    Q2Q1List.set("A21",         thA21);
    Q2Q1List.set("A11_9Pt",     thA11_9Pt);

    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&Q2Q1List), *comm);

    std::cout << "Input parameters: " << *stratimikosList << std::endl;

    // Stratimikos vodou
    typedef Thyra::PreconditionerFactoryBase<SC>             Base;
    typedef Thyra::Ifpack2PreconditionerFactory<tCrsMatrix > Impl;
    typedef Thyra::LinearOpWithSolveFactoryBase<SC>          LOWSFB;
    typedef Thyra::LinearOpWithSolveBase<SC>                 LOWSB;
    typedef Thyra::MultiVectorBase<SC>                       TH_Mvb;

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

    //Thyra::addMueLuToStratimikosBuilder(linearSolverBuilder);

    Stratimikos::enableMueLu(linearSolverBuilder);                                          // Epetra
    Stratimikos::enableMueLuTpetraQ2Q1<LO,GO,NO>(linearSolverBuilder, "MueLu-TpetraQ2Q1");  // Tpetra

    linearSolverBuilder.setParameterList(stratimikosList);
    RCP<const LOWSFB> lowsFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
    RCP<LOWSB       > nsA = lowsFactory->createOp();

    // I've hacked together a big matrix that does not use strided maps by
    // simply reading the data again. Normally, this would be supplied by Eric
    // Cyr and would be Teko operators.

    int numElem = A12->getRangeMap()->getNodeNumElements() + A21->getRangeMap()->getNodeNumElements();
    RCP<const tMap> fullMap = Utilities::Map2TpetraMap(*(MapFactory::createUniformContigMap(Xpetra::UseTpetra, numElem, comm)));

    RCP<tOperator> A;
    if (!binary)
      A = Reader<tCrsMatrix>::readSparseFile((prefix + "BigA.mm").c_str(), fullMap, fullMap, fullMap, fullMap, true, true, false);
    else
      A = ReadBinary<SC,LO,GO,NO>((prefix + "BigA.dat").c_str(), comm);

    const RCP<Thyra::LinearOpBase<SC> > thA = Thyra::createLinearOp(A);
    Thyra::initializeOp<SC>(*lowsFactory, thA, nsA.ptr());

    RCP<tMultiVector> tX = Tpetra::createVector<SC>(fullMap), tB;
    if (rhs == "") {
      tX->randomize();

      tB = Tpetra::createVector<SC>(fullMap);
      A->apply(*tX, *tB);
    } else {
      typedef Tpetra::MatrixMarket::Reader<tCrsMatrix> reader_type;
      tB = reader_type::readDenseFile(rhs.c_str(), fullMap->getComm(), fullMap->getNode(), fullMap);
    }

    tX->putScalar(0.0);

    // Set the initial guess Dirichlet points to the proper value.
    // This step is pretty important as the preconditioner may return zero at Dirichlet points
    ArrayRCP<const bool> dirBCs = Utilities::DetectDirichletRows(*MueLu::TpetraCrs_To_XpetraMatrix(rcp_dynamic_cast<tCrsMatrix>(A)));
    ArrayRCP<SC>         tXdata = tX->getDataNonConst(0);
    ArrayRCP<const SC>   tBdata = tB->getData(0);
    for (LO i = 0; i < tXdata.size(); i++)
      if (dirBCs[i])
        tXdata[i] = tBdata[i];

    RCP<TH_Mvb> sX = Thyra::createMultiVector(tX);
    RCP<TH_Mvb> sB = Thyra::createMultiVector(tB);

    Thyra::SolveStatus<SC> solveStatus = Thyra::solve(*nsA, Thyra::NOTRANS, *sB, sX.ptr());

    if (printTimings) {
      const bool alwaysWriteLocal = false;
      const bool writeGlobalStats = true;
      const bool writeZeroTimers  = false;
      const bool ignoreZeroTimers = true;
      const std::string filter    = "";
      TimeMonitor::summarize(comm.ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                             writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
