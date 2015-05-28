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
#include <Stratimikos_MueluTpetraHelpers.hpp>
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

#include <Xpetra_MapFactory.hpp>

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_Utilities.hpp"


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

    std::string xmlFileName  = "driver.xml";   clp.setOption("xml",      &xmlFileName,   "read parameters from a file [default = 'driver.xml']");
    double      tol          = 1e-12;          clp.setOption("tol",      &tol,           "solver convergence tolerance");
    int         n            = 17;             clp.setOption("n",        &n,             "problem size (1D)");
    int         maxLevels    = 4;              clp.setOption("nlevels",  &maxLevels,     "max num levels");
    std::string type         = "structured";   clp.setOption("type",     &type,          "structured/unstructured");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    typedef Tpetra::CrsMatrix<SC,LO,GO>           TP_Crs;
    typedef Tpetra::Operator<SC,LO,GO>            TP_Op;
    typedef Tpetra::MultiVector<SC,LO,GO,NO>      TP_Mv;
    typedef Tpetra::Map<LO,GO,NO>                 TP_Map;
    typedef Thyra::TpetraVectorSpace<SC,LO,GO,NO> THTP_Vs;
    RCP<NO> node = Tpetra::DefaultPlatform::getDefaultPlatform().getNode();

    // Read data from files
    std::string prefix = "Q2Q1_" + MueLu::toString(n) + "x" + MueLu::toString(n) + "_";
    RCP<TP_Op> A11     = Reader<TP_Crs>::readSparseFile((prefix + "A.mm").c_str(),          comm, node);
    RCP<TP_Op> A21     = Reader<TP_Crs>::readSparseFile((prefix + "B.mm").c_str(),          comm, node);
    RCP<TP_Op> A12     = Reader<TP_Crs>::readSparseFile((prefix + "Bt.mm").c_str(),         comm, node);
    RCP<TP_Op> A119Pt  = Reader<TP_Crs>::readSparseFile((prefix + "AForPat.mm").c_str(),    comm, node);
    RCP<const TP_Map> cmap1 = A11->getDomainMap(), cmap2 = A12->getDomainMap();
    RCP<TP_Mv> Vcoords = Reader<TP_Crs>::readDenseFile ((prefix + "VelCoords.mm").c_str(),  comm, node, cmap1);
    RCP<TP_Mv> Pcoords = Reader<TP_Crs>::readDenseFile ((prefix + "PresCoords.mm").c_str(), comm, node, cmap2);

    // For now, we assume that p2v maps local pressure DOF to a local x-velocity DOF
    ArrayRCP<const SC> slop = Utils2::ReadMultiVector((prefix + "p2vMap.mm").c_str(),
                                                      Xpetra::toXpetra(A21->getRangeMap()))->getData(0);
    ArrayRCP<LO> p2vMap(n*n);
    for (int i = 0; i < n*n; i++)
      p2vMap[i] = as<LO>(slop[i]);

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
    GmresDetails.set("Maximum Iterations",      20);
    GmresDetails.set("Convergence Tolerance",   1e-12);
    GmresDetails.set("Verbosity",               Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    GmresDetails.set("Output Frequency",        1);
    GmresDetails.set("Output Style",            Belos::Brief);

    ParameterList& Q2Q1List = stratimikosList->sublist("Preconditioner Types").sublist("MueLu-TpetraQ2Q1");
    Q2Q1List.set("Velcoords",   Vcoords);
    Q2Q1List.set("Prescoords",  Pcoords);
    Q2Q1List.set("p2vMap",      p2vMap);
    Q2Q1List.set("A11",         thA11);
    Q2Q1List.set("A12",         thA12);
    Q2Q1List.set("A21",         thA21);
    Q2Q1List.set("A11_9Pt",     thA11_9Pt);

    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&Q2Q1List), *comm);

    // Stratimikos vodou
    typedef Thyra::PreconditionerFactoryBase<SC>         Base;
    typedef Thyra::Ifpack2PreconditionerFactory<TP_Crs > Impl;
    typedef Thyra::LinearOpWithSolveFactoryBase<SC>      LOWSFB;
    typedef Thyra::LinearOpWithSolveBase<SC>             LOWSB;
    typedef Thyra::MultiVectorBase<SC>                   TH_Mvb;

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

    Thyra::addMueLuToStratimikosBuilder(linearSolverBuilder);
    Stratimikos::enableMueLuTpetraQ2Q1<LO,GO,NO>(linearSolverBuilder, "MueLu-TpetraQ2Q1");

    linearSolverBuilder.setParameterList(stratimikosList);
    RCP<const LOWSFB> lowsFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
    RCP<LOWSB       > nsA = lowsFactory->createOp();

    // I've hacked together a big matrix that does not use strided maps by
    // simply reading the data again. Normally, this would be supplied by Eric
    // Cyr and would be Teko operators.

    int numElem = A12->getRangeMap()->getNodeNumElements() + A21->getRangeMap()->getNodeNumElements();
    RCP<const TP_Map> fullMap = Utils::Map2TpetraMap(*(MapFactory::createUniformContigMap(Xpetra::UseTpetra, numElem, comm)));

    RCP<TP_Op> A = Tpetra::MatrixMarket::Reader<TP_Crs>::readSparseFile((prefix + "BigA.mm").c_str(), fullMap, fullMap, fullMap, fullMap, true, true, false);

    const RCP<Thyra::LinearOpBase<SC> > thA = Thyra::createLinearOp(A);
    Thyra::initializeOp<SC>(*lowsFactory, thA, nsA.ptr());

    RCP<TP_Mv> tX = Tpetra::createVector<SC>(fullMap);
#if 1
    tX->randomize();

    RCP<TP_Mv> tB = Tpetra::createVector<SC>(fullMap);
    A->apply(*tX, *tB);
#else
    typedef Tpetra::MatrixMarket::Reader<TP_Crs> reader_type;
    RCP<TP_Mv> tB = reader_type::readDenseFile((prefix + "rhs.mm").c_str(), fullMap->getComm(), fullMap->getNode(), fullMap);
#endif

    tX->putScalar(0.0);

    RCP<TH_Mvb> sX = Thyra::createMultiVector(tX);
    RCP<TH_Mvb> sB = Thyra::createMultiVector(tB);

    Thyra::SolveStatus<SC> solveStatus = Thyra::solve(*nsA, Thyra::NOTRANS, *sB, sX.ptr());
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
