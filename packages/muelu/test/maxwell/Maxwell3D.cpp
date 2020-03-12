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

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StackedTimer.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_IO.hpp>

// MueLu
#include <MueLu_RefMaxwell.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_TestHelpers_Common.hpp>
#include <MueLu_Exceptions.hpp>

#ifdef HAVE_MUELU_TPETRA
#include <MueLu_TpetraOperator.hpp>
#endif

// Belos
#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
//#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

// Stratimikos
#ifdef HAVE_MUELU_STRATIMIKOS
// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
// Stratimikos includes
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_MueLuHelpers.hpp>
#endif

// Support for ML interface
#if defined(HAVE_MUELU_ML) and defined(HAVE_MUELU_EPETRA)
#include <Xpetra_EpetraOperator.hpp>
#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "ml_RefMaxwell.h"

// Helper functions for compilation purposes
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct ML_Wrapper{
  static void Generate_ML_MaxwellPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& SM,
                                                Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& D0,
                                                Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Kn,
                                                Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& nullspace,
                                                Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node> >& coords,
                                                Teuchos::ParameterList & mueluList,
                                                Teuchos::RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& mlopX) {
    throw std::runtime_error("Template parameter mismatch");
  }

  static void Generate_ML_RefMaxwellPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& SM,
                                                   Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& D0,
                                                   Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Ms,
                                                   Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& M0inv,
                                                   Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& M1,
                                                   Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& nullspace,
                                                   Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& node_material,
                                                   Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node> >& coords,
                                                   Teuchos::ParameterList & mueluList,
                                                   Teuchos::RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & mlopX) {
    throw std::runtime_error("Template parameter mismatch");
  }
};


template<class GlobalOrdinal>
struct ML_Wrapper<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> {
  static void Generate_ML_MaxwellPreconditioner(Teuchos::RCP<Xpetra::Matrix<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& SM,
                                                Teuchos::RCP<Xpetra::Matrix<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& D0,
                                                Teuchos::RCP<Xpetra::Matrix<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& Kn,
                                                Teuchos::RCP<Xpetra::MultiVector<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& nullspace,
                                                Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<double>::coordinateType,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& coords,
                                                Teuchos::ParameterList & mueluList,
                                                Teuchos::RCP<Xpetra::Operator<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& mlopX) {
    typedef double SC;
    typedef int LO;
    typedef GlobalOrdinal GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
    typedef typename Teuchos::ScalarTraits<SC>::coordinateType coordinate_type;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef typename Xpetra::Matrix<SC,LO,GO,NO> Matrix;

    RCP<const Epetra_CrsMatrix> epetraSM = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(SM);
    RCP<const Epetra_CrsMatrix> epetraD0 = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(D0);
    if(!coords.is_null()) {
      RCP<const Epetra_MultiVector> epetraCoord =  MueLu::Utilities<coordinate_type,LO,GO,NO>::MV2EpetraMV(coords);
      if(epetraCoord->NumVectors() > 0)  mueluList.set("x-coordinates",(*epetraCoord)[0]);
      if(epetraCoord->NumVectors() > 1)  mueluList.set("y-coordinates",(*epetraCoord)[1]);
      if(epetraCoord->NumVectors() > 2)  mueluList.set("z-coordinates",(*epetraCoord)[2]);
    }
    if(!nullspace.is_null()) {
      RCP<const Epetra_MultiVector> epetraNullspace =  MueLu::Utilities<SC,LO,GO,NO>::MV2EpetraMV(nullspace);
      mueluList.set("null space: dimension",epetraNullspace->NumVectors());
      mueluList.set("null space: vectors",(*epetraNullspace)[0]);
      mueluList.set("null space: type","pre-computed");
    }
    RCP<const Epetra_CrsMatrix> epetraKn;
    if (Kn.is_null()) {
      RCP<Matrix> temp = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(SM->getRangeMap());
      Xpetra::MatrixMatrix<SC,LO,GO,NO>::Multiply(*SM,false,*D0,false,*temp,true,true);
      RCP<Matrix> Kn2 = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(D0->getDomainMap());
      Xpetra::MatrixMatrix<SC,LO,GO,NO>::Multiply(*D0,true,*temp,false,*Kn2,true,true);
      epetraKn = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(Kn2);
    } else
      epetraKn = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(Kn);

    RCP<Epetra_Operator> mlop  = rcp<Epetra_Operator>(new ML_Epetra::MultiLevelPreconditioner(*epetraSM,*epetraD0,*epetraKn,mueluList,true));
#if defined(HAVE_MUELU_BELOS)
    // NOTE: Belos needs the Apply() and AppleInverse() routines of ML swapped.  So...
    mlop = rcp<Belos::EpetraPrecOp>(new Belos::EpetraPrecOp(mlop));
#endif

    mlopX = rcp(new Xpetra::EpetraOperator<GO,NO>(mlop));
  }

  static void Generate_ML_RefMaxwellPreconditioner(Teuchos::RCP<Xpetra::Matrix<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& SM,
                                                   Teuchos::RCP<Xpetra::Matrix<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& D0,
                                                   Teuchos::RCP<Xpetra::Matrix<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& Ms,
                                                   Teuchos::RCP<Xpetra::Matrix<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& M0inv,
                                                   Teuchos::RCP<Xpetra::Matrix<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& M1,
                                                   Teuchos::RCP<Xpetra::MultiVector<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& nullspace,
                                                   Teuchos::RCP<Xpetra::MultiVector<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& node_material,
                                                   Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<double>::coordinateType,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& coords,
                                                   Teuchos::ParameterList & mueluList,
                                                   Teuchos::RCP<Xpetra::Operator<double,int,GlobalOrdinal,Kokkos::Compat::KokkosSerialWrapperNode> >& mlopX) {
    typedef double SC;
    typedef int LO;
    typedef GlobalOrdinal GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
    typedef typename Teuchos::ScalarTraits<SC>::coordinateType coordinate_type;
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<const Epetra_CrsMatrix> epetraSM    = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(SM);
    RCP<const Epetra_CrsMatrix> epetraD0    = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(D0);
    RCP<const Epetra_CrsMatrix> epetraM0inv = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(M0inv);
    RCP<const Epetra_CrsMatrix> epetraMs;
    RCP<const Epetra_CrsMatrix> epetraM1    = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(M1);
    if (!Ms.is_null())
      epetraMs    = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(Ms);
    else
      epetraMs = epetraM1;
    mueluList.set("D0",epetraD0);
    mueluList.set("Ms",epetraMs);
    mueluList.set("M0inv",epetraM0inv);
    mueluList.set("M1",epetraM1);
    if(!coords.is_null()) {
      RCP<const Epetra_MultiVector> epetraCoord =  MueLu::Utilities<coordinate_type,LO,GO,NO>::MV2EpetraMV(coords);
      if(epetraCoord->NumVectors() > 0)  mueluList.sublist("refmaxwell: 11list").set("x-coordinates",(*epetraCoord)[0]);
      if(epetraCoord->NumVectors() > 1)  mueluList.sublist("refmaxwell: 11list").set("y-coordinates",(*epetraCoord)[1]);
      if(epetraCoord->NumVectors() > 2)  mueluList.sublist("refmaxwell: 11list").set("z-coordinates",(*epetraCoord)[2]);
    }
    if(!node_material.is_null()) {
      RCP<const Epetra_MultiVector> epetraMaterial =  MueLu::Utilities<coordinate_type,LO,GO,NO>::MV2EpetraMV(node_material);
      mueluList.sublist("refmaxwell: 11list").set("material coordinates",(*epetraMaterial)[0]);
    }
    if(!nullspace.is_null()) {
      RCP<const Epetra_MultiVector> epetraNullspace =  MueLu::Utilities<SC,LO,GO,NO>::MV2EpetraMV(nullspace);
      mueluList.sublist("refmaxwell: 11list").set("null space: dimension",epetraNullspace->NumVectors());
      mueluList.sublist("refmaxwell: 11list").set("null space: vectors",(*epetraNullspace)[0]);
      mueluList.sublist("refmaxwell: 11list").set("null space: type","pre-computed");
    }

    RCP<Epetra_Operator> mlop  = rcp<Epetra_Operator>(new ML_Epetra::RefMaxwellPreconditioner(*epetraSM,mueluList,true));
#if defined(HAVE_MUELU_BELOS)
    // NOTE: Belos needs the Apply() and AppleInverse() routines of ML swapped.  So...
    mlop = rcp<Belos::EpetraPrecOp>(new Belos::EpetraPrecOp(mlop));
#endif

    mlopX = rcp(new Xpetra::EpetraOperator<GO,NO>(mlop));
  }
};
#endif

// Main wrappers struct
// Because C++ doesn't support partial template specialization of functions.
// By default, do not try to run Stratimikos, since that only works for Scalar=double.
template<typename Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
struct MainWrappers {
  static int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]);
};


// Partial template specialization on SC=double
// This code branch gives the option to run with Stratimikos.
template<class LocalOrdinal, class GlobalOrdinal, class Node>
struct MainWrappers<double,LocalOrdinal,GlobalOrdinal,Node> {
  static int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]);
};


// By default, do not try to run Stratimikos, since that only works for Scalar=double.
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int MainWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)

#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP; using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);

    bool        printTimings      = true;               clp.setOption("timings", "notimings",  &printTimings,      "print timings to screen");
    std::string timingsFormat     = "table-fixed";      clp.setOption("time-format",           &timingsFormat,     "timings format (table-fixed | table-scientific | yaml)");
    double scaling                = 1.0;                clp.setOption("scaling",               &scaling,           "scale mass term");
    std::string solverName = "Belos";                   clp.setOption("solverName",            &solverName,        "Name of iterative linear solver "
                                                                                                                  "to use for solving the linear system. "
                                                                                                                  "(\"Belos\")");
    std::string belosSolverType   = "Block CG";         clp.setOption("belosSolverType",       &belosSolverType,   "Name of the Belos linear solver");
    std::string precType          = "MueLu-RefMaxwell"; clp.setOption("precType",              &precType,          "preconditioner to use (MueLu-RefMaxwell|ML-RefMaxwell|none)");
    std::string xml;
    if (!TYPE_EQUAL(SC, std::complex<double>) && !TYPE_EQUAL(SC, std::complex<float>))
      xml = "Maxwell.xml";
    else
      xml = "Maxwell_complex.xml";
                                                        clp.setOption("xml",                   &xml,               "xml file with solver parameters");
    double      tol               = 1e-10;              clp.setOption("tol",                   &tol,               "solver convergence tolerance");
    int         maxIts            = 200;                clp.setOption("its",                   &maxIts,            "maximum number of solver iterations");
    bool        use_stacked_timer = false;              clp.setOption("stacked-timer", "no-stacked-timer", &use_stacked_timer, "use stacked timer");

    std::string S_file, SM_file, M1_file, M0_file, M0inv_file, D0_file, coords_file, rhs_file="", nullspace_file="", material_file = "", Ms_file="";

    if (!TYPE_EQUAL(SC, std::complex<double>)) {
      S_file = "S.mat";
      SM_file = "";
      M1_file = "M1.mat";
      M0_file = "M0.mat";
      M0inv_file = "";
      D0_file = "D0.mat";
    } else {
      S_file = "S_complex.mat";
      SM_file = "";
      M1_file = "M1_complex.mat";
      M0_file = "M0_complex.mat";
      M0inv_file = "";
      D0_file = "D0_complex.mat";
    }
    coords_file = "coords.mat";

    clp.setOption("S", &S_file);
    clp.setOption("SM", &SM_file);
    clp.setOption("Ms", &Ms_file);
    clp.setOption("M1", &M1_file);
    clp.setOption("M0", &M0_file);
    clp.setOption("M0inv", &M0inv_file);
    clp.setOption("D0", &D0_file);
    clp.setOption("coords", &coords_file);
    clp.setOption("material", &material_file);
    clp.setOption("nullspace", &nullspace_file);
    clp.setOption("rhs", &rhs_file);

    clp.recogniseAllOptions(true);
    switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    comm->barrier();

    Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
    if (use_stacked_timer)
      stacked_timer = rcp(new Teuchos::StackedTimer("Maxwell Driver"));
    Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

    auto globalTimeMonitor = TimeMonitor::getNewTimer("Maxwell: S - Global Time");
    auto tm                = TimeMonitor::getNewTimer("Maxwell: 1 - Read and Build Matrices");

    // Read matrices in from files
    // gradient matrix
    RCP<Matrix> D0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(D0_file, lib, comm);
    // maps for nodal and edge matrices
    RCP<const Map> node_map = D0_Matrix->getDomainMap();
    RCP<const Map> edge_map = D0_Matrix->getRangeMap();
    // edge mass matrix
    RCP<Matrix> M1_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M1_file, edge_map);
    RCP<Matrix> Ms_Matrix;
    if (Ms_file != "")
      Ms_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(Ms_file, edge_map);
    else
      Ms_Matrix = M1_Matrix;
    // build stiffness plus mass matrix (SM_Matrix)
    RCP<Matrix> SM_Matrix;
    if (SM_file == "") {
      // edge stiffness matrix
      RCP<Matrix> S_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(S_file, edge_map);
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*S_Matrix,false,(SC)1.0,*M1_Matrix,false,scaling,SM_Matrix,*out);
      SM_Matrix->fillComplete();
    } else {
      SM_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(SM_file, edge_map);
    }
    RCP<Matrix> M0inv_Matrix;
    if (M0inv_file == "") {
      // nodal mass matrix
      RCP<Matrix> M0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M0_file, node_map);
      // build lumped mass matrix inverse (M0inv_Matrix)
      RCP<Vector> diag = Utilities::GetLumpedMatrixDiagonal(M0_Matrix);
      RCP<CrsMatrixWrap> M0inv_MatrixWrap = Teuchos::rcp(new CrsMatrixWrap(node_map, node_map, 0));
      RCP<CrsMatrix> M0inv_CrsMatrix = M0inv_MatrixWrap->getCrsMatrix();
      Teuchos::ArrayRCP<size_t> rowPtr;
      Teuchos::ArrayRCP<LO> colInd;
      Teuchos::ArrayRCP<SC> values;
      Teuchos::ArrayRCP<const SC> diags = diag->getData(0);
      size_t nodeNumElements = node_map->getNodeNumElements();
      M0inv_CrsMatrix->allocateAllValues(nodeNumElements, rowPtr, colInd, values);
      SC ONE = (SC)1.0;
      for (size_t i = 0; i < nodeNumElements; i++) {
        rowPtr[i] = i;  colInd[i] = i;  values[i] = ONE / diags[i];
      }
      rowPtr[nodeNumElements] = nodeNumElements;
      M0inv_CrsMatrix->setAllValues(rowPtr, colInd, values);
      M0inv_CrsMatrix->expertStaticFillComplete(node_map, node_map);
      M0inv_Matrix = Teuchos::rcp_dynamic_cast<Matrix>(M0inv_MatrixWrap);
    } else if (M0inv_file == "none") {
      // pass
    } else {
      M0inv_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M0inv_file, node_map);
    }
    // coordinates
    RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO> > coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO>::ReadMultiVector(coords_file, node_map);

    RCP<MultiVector> nullspace = Teuchos::null;
    if (nullspace_file != "")
      nullspace = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(nullspace_file, edge_map);

    RCP<MultiVector> material = Teuchos::null;
    if (material_file != "") {
      material  = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(material_file, node_map);
    }

    // set parameters
    Teuchos::ParameterList params;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xml,Teuchos::Ptr<Teuchos::ParameterList>(&params),*comm);

    // setup LHS, RHS
    RCP<MultiVector> B;
    if (rhs_file == "") {
      B = MultiVectorFactory::Build(edge_map,1);
      RCP<MultiVector> vec = MultiVectorFactory::Build(edge_map,1);
      vec -> putScalar((SC)1.0);
      SM_Matrix->apply(*vec,*B);
    } else
      B = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(rhs_file, edge_map);

    RCP<MultiVector> X = MultiVectorFactory::Build(edge_map,1);
    X -> putScalar((SC)0.0);

    comm->barrier();
    tm = Teuchos::null;

    if (solverName == "Belos") {
      auto tm2  = TimeMonitor::getNewTimer("Maxwell: 2 - Build Belos solver etc");

      // construct preconditioner
      RCP<MueLu::RefMaxwell<SC,LO,GO,NO> > preconditioner;
      if (precType=="MueLu-RefMaxwell")
        preconditioner
          = rcp( new MueLu::RefMaxwell<SC,LO,GO,NO>(SM_Matrix,D0_Matrix,Ms_Matrix,M0inv_Matrix,
                                                    M1_Matrix,nullspace,coords,params) );


#ifdef HAVE_MUELU_TPETRA
      {
        // A test to make sure we can wrap this guy as a MueLu::TpetraOperator
        RCP<Operator> precOp = Teuchos::rcp_dynamic_cast<Operator>(preconditioner);
        MueLu::TpetraOperator<SC,LO,GO,NO> OpT(precOp);
      }
#endif

      // Belos linear problem
#ifdef HAVE_MUELU_BELOS
      typedef MultiVector          MV;
      typedef Belos::OperatorT<MV> OP;
      Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(SM_Matrix)); // Turns a Xpetra::Matrix object into a Belos operator

      RCP<Belos::LinearProblem<SC, MV, OP> > problem = rcp( new Belos::LinearProblem<SC, MV, OP>() );
      problem -> setOperator( belosOp );
      Teuchos::RCP<OP> belosPrecOp;
      if (precType!="none") {
        belosPrecOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(preconditioner)); // Turns a Xpetra::Matrix object into a Belos operator
        problem -> setRightPrec( belosPrecOp );
      }

      problem -> setProblem( X, B );

      bool set = problem->setProblem();
      if (set == false) {
        *out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        return EXIT_FAILURE;
      }

      // Belos solver
      RCP< Belos::SolverManager<SC, MV, OP> > solver;
      RCP< Belos::SolverFactory<SC, MV,OP> > factory = rcp( new  Belos::SolverFactory<SC,MV,OP>() );
      RCP<Teuchos::ParameterList> belosParams
        = rcp( new Teuchos::ParameterList() );
      belosParams->set("Maximum Iterations", maxIts);
      belosParams->set("Convergence Tolerance",tol);
      belosParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosParams->set("Output Frequency",1);
      belosParams->set("Output Style",Belos::Brief);
      belosParams->set("Implicit Residual Scaling","None");
      belosParams->set("Explicit Residual Test",true);
      belosParams->set("Explicit Residual Scaling","None");
      solver = factory->create(belosSolverType,belosParams);

      comm->barrier();
      tm2=Teuchos::null;
      auto tm3  = TimeMonitor::getNewTimer("Maxwell: 3 - Solve");

      // set problem and solve
      solver -> setProblem( problem );
      Belos::ReturnType status = solver -> solve();
      int iters = solver -> getNumIters();
      success = (iters<50 && status == Belos::Converged);
      if (success)
        *out << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
      else
        *out << "FAILURE! Belos did not converge fast enough." << std::endl;
      tm3 = Teuchos::null;
    }
    comm->barrier();
    globalTimeMonitor = Teuchos::null;

    if (printTimings) {
      RCP<Teuchos::ParameterList> reportParams = rcp(new Teuchos::ParameterList);
      if (timingsFormat == "yaml") {
        reportParams->set("Report format",             "YAML");            // "Table" or "YAML"
        reportParams->set("YAML style",                "compact");         // "spacious" or "compact"
      }
      reportParams->set("How to merge timer sets",   "Union");
      reportParams->set("alwaysWriteLocal",          false);
      reportParams->set("writeGlobalStats",          true);
      reportParams->set("writeZeroTimers",           false);
      // FIXME: no "ignoreZeroTimers"

      const std::string filter = "";

      std::ios_base::fmtflags ff(out->flags());
      if (timingsFormat == "table-fixed") *out << std::fixed;
      else * out << std::scientific;
      if (use_stacked_timer) {
        stacked_timer->stop("Maxwell Driver");
        Teuchos::StackedTimer::OutputOptions options;
        options.output_fraction = options.output_histogram = options.output_minmax = true;
        stacked_timer->report(*out, comm, options);
      } else
        TimeMonitor::report(comm.ptr(), *out, filter, reportParams);
       *out << std::setiosflags(ff);
    }

    TimeMonitor::clearCounters();
#endif // #ifdef HAVE_MUELU_BELOS
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
#else
  return EXIT_SUCCESS;
#endif
} // main


// This code branch gives the option to run with Stratimikos.
template<class LocalOrdinal, class GlobalOrdinal, class Node>
int MainWrappers<double,LocalOrdinal,GlobalOrdinal,Node>::main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
  typedef double Scalar;
#include <MueLu_UseShortNames.hpp>

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)

#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP; using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);

    bool        printTimings      = true;               clp.setOption("timings", "notimings",  &printTimings,      "print timings to screen");
    std::string timingsFormat     = "table-fixed";      clp.setOption("time-format",           &timingsFormat,     "timings format (table-fixed | table-scientific | yaml)");
    double scaling                = 1.0;                clp.setOption("scaling",               &scaling,           "scale mass term");
    std::string solverName        = "Belos";            clp.setOption("solverName",            &solverName,        "Name of iterative linear solver "
                                                                                                                  "to use for solving the linear system. "
                                                                                                                  "(\"Belos\" or \"Stratimikos\")");
    std::string belosSolverType   = "Block CG";         clp.setOption("belosSolverType",       &belosSolverType,   "Name of the Belos linear solver");
    std::string precType          = "MueLu-RefMaxwell"; clp.setOption("precType",              &precType,          "preconditioner to use (MueLu-RefMaxwell|ML-RefMaxwell|none)");
    std::string xml               = "";                 clp.setOption("xml",                   &xml,               "xml file with solver parameters (default: \"Maxwell.xml\")");
    double      tol               = 1e-10;              clp.setOption("tol",                   &tol,               "solver convergence tolerance");
    int         maxIts            = 200;                clp.setOption("its",                   &maxIts,            "maximum number of solver iterations");
    bool        use_stacked_timer = false;              clp.setOption("stacked-timer", "no-stacked-timer", &use_stacked_timer, "use stacked timer");


    std::string S_file, SM_file, M1_file, M0_file, M0inv_file, D0_file, coords_file, rhs_file="", nullspace_file="", material_file = "", Ms_file="", sol_file="", Kn_file="";

    S_file = "S.mat";
    SM_file = "";
    M1_file = "M1.mat";
    M0_file = "M0.mat";
    M0inv_file = "";
    D0_file = "D0.mat";
    coords_file = "coords.mat";

    clp.setOption("S", &S_file);
    clp.setOption("SM", &SM_file);
    clp.setOption("Ms", &Ms_file);
    clp.setOption("M1", &M1_file);
    clp.setOption("M0", &M0_file);
    clp.setOption("Kn", &Kn_file);
    clp.setOption("M0inv", &M0inv_file);
    clp.setOption("D0", &D0_file);
    clp.setOption("coords", &coords_file);
    clp.setOption("rhs", &rhs_file);
    clp.setOption("nullspace", &nullspace_file);
    clp.setOption("material", &material_file);

    clp.recogniseAllOptions(true);
    switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    if (xml == ""){
      if (precType == "MueLu-RefMaxwell")
        xml = "Maxwell.xml";
      else if (precType == "ML-RefMaxwell")
        xml = "Maxwell_ML.xml";
      else if (precType == "ML-Maxwell")
        xml = "Maxwell_ML1.xml";
    }

    comm->barrier();

    Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
    if (use_stacked_timer)
      stacked_timer = rcp(new Teuchos::StackedTimer("Maxwell Driver"));
    Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

    auto globalTimeMonitor = TimeMonitor::getNewTimer("Maxwell: S - Global Time");
    auto tm                = TimeMonitor::getNewTimer("Maxwell: 1 - Read and Build Matrices");

    // Read matrices in from files
    // gradient matrix
    RCP<Matrix> D0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(D0_file, lib, comm);
    // maps for nodal and edge matrices
    RCP<const Map> node_map = D0_Matrix->getDomainMap();
    RCP<const Map> edge_map = D0_Matrix->getRangeMap();
    // edge mass matrix
    RCP<Matrix> M1_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M1_file, edge_map);
    RCP<Matrix> Ms_Matrix;
    if (Ms_file != "")
      Ms_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(Ms_file, edge_map);
    else
      Ms_Matrix = M1_Matrix;
    // build stiffness plus mass matrix (SM_Matrix)
    RCP<Matrix> SM_Matrix, S_Matrix;
    if (SM_file == "") {
      // edge stiffness matrix
      S_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(S_file, edge_map);
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*S_Matrix,false,(SC)1.0,*M1_Matrix,false,scaling,SM_Matrix,*out);
      SM_Matrix->fillComplete();
    } else {
      SM_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(SM_file, edge_map);
    }
    RCP<Matrix> M0inv_Matrix, M0_Matrix;
    if (M0inv_file == "") {
      // nodal mass matrix
      M0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M0_file, node_map);
      // build lumped mass matrix inverse (M0inv_Matrix)
      RCP<Vector> diag = Utilities::GetLumpedMatrixDiagonal(M0_Matrix);
      RCP<CrsMatrixWrap> M0inv_MatrixWrap = Teuchos::rcp(new CrsMatrixWrap(node_map, node_map, 0));
      RCP<CrsMatrix> M0inv_CrsMatrix = M0inv_MatrixWrap->getCrsMatrix();
      Teuchos::ArrayRCP<size_t> rowPtr;
      Teuchos::ArrayRCP<LO> colInd;
      Teuchos::ArrayRCP<SC> values;
      Teuchos::ArrayRCP<const SC> diags = diag->getData(0);
      size_t nodeNumElements = node_map->getNodeNumElements();
      M0inv_CrsMatrix->allocateAllValues(nodeNumElements, rowPtr, colInd, values);
      SC ONE = (SC)1.0;
      for (size_t i = 0; i < nodeNumElements; i++) {
        rowPtr[i] = i;  colInd[i] = i;  values[i] = ONE / diags[i];
      }
      rowPtr[nodeNumElements] = nodeNumElements;
      M0inv_CrsMatrix->setAllValues(rowPtr, colInd, values);
      M0inv_CrsMatrix->expertStaticFillComplete(node_map, node_map);
      M0inv_Matrix = Teuchos::rcp_dynamic_cast<Matrix>(M0inv_MatrixWrap);
    } else if (M0inv_file == "none") {
      // pass
    } else {
      M0inv_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(M0inv_file, node_map);
    }
    // coordinates
    RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::ReadMultiVector(coords_file, node_map);

    RCP<MultiVector> nullspace = Teuchos::null;
    if (nullspace_file != "")
      nullspace = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(nullspace_file, edge_map);
    
    RCP<MultiVector> material = Teuchos::null;
    if (material_file != "") {
      material  = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(material_file, node_map);
    }

    // set parameters
    Teuchos::ParameterList params;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xml,Teuchos::Ptr<Teuchos::ParameterList>(&params),*comm);

    // setup LHS, RHS
    RCP<MultiVector> B;
    if (rhs_file == "") {
      B = MultiVectorFactory::Build(edge_map,1);
      RCP<MultiVector> vec = MultiVectorFactory::Build(edge_map,1);
      vec -> putScalar((SC)1.0);
      SM_Matrix->apply(*vec,*B);
    } else
      B = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(rhs_file, edge_map);
    RCP<MultiVector> X = MultiVectorFactory::Build(edge_map,1);
    X -> putScalar((SC)0.0);

    comm->barrier();
    tm = Teuchos::null;

    if (solverName == "Belos") {
      auto tm2 = TimeMonitor::getNewTimer("Maxwell: 2 - Build Belos solver etc");

      // construct preconditioner
      RCP<Operator> preconditioner;
      if (precType=="MueLu-RefMaxwell")
        preconditioner = rcp( new MueLu::RefMaxwell<SC,LO,GO,NO>(SM_Matrix,D0_Matrix,Ms_Matrix,M0inv_Matrix,
                                                                 M1_Matrix,nullspace,coords,params) );
      else if (precType=="ML-RefMaxwell") {
#if defined(HAVE_MUELU_ML) and defined(HAVE_MUELU_EPETRA)
        TEUCHOS_ASSERT(lib==Xpetra::UseEpetra);
        ML_Wrapper<SC, LO, GO, NO>::Generate_ML_RefMaxwellPreconditioner(SM_Matrix,D0_Matrix,Ms_Matrix,M0inv_Matrix,
                                                                         M1_Matrix,nullspace,material,coords,params,preconditioner);
#else
        throw std::runtime_error("ML not available.");
#endif
      } else if (precType=="ML-Maxwell") {
#if defined(HAVE_MUELU_ML) and defined(HAVE_MUELU_EPETRA)
        TEUCHOS_ASSERT(lib==Xpetra::UseEpetra);
        RCP<Matrix> Kn_Matrix;
        if (Kn_file!="")
          Kn_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read(Kn_file, node_map);
        ML_Wrapper<SC, LO, GO, NO>::Generate_ML_MaxwellPreconditioner(SM_Matrix,D0_Matrix,Kn_Matrix,
                                                                      nullspace,coords,params,preconditioner);
#else
        throw std::runtime_error("ML not available.");
#endif
      }

#ifdef HAVE_MUELU_TPETRA
      {
        // A test to make sure we can wrap this guy as a MueLu::TpetraOperator
        RCP<Operator> precOp = Teuchos::rcp_dynamic_cast<Operator>(preconditioner);
        MueLu::TpetraOperator<SC,LO,GO,NO> OpT(precOp);
      }
#endif


      // Belos linear problem
#ifdef HAVE_MUELU_BELOS
      typedef MultiVector          MV;
      typedef Belos::OperatorT<MV> OP;
      Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(SM_Matrix)); // Turns a Xpetra::Matrix object into a Belos operator

      RCP<Belos::LinearProblem<SC, MV, OP> > problem = rcp( new Belos::LinearProblem<SC, MV, OP>() );
      problem -> setOperator( belosOp );
      Teuchos::RCP<OP> belosPrecOp;
      if (precType != "none") {
        belosPrecOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(preconditioner)); // Turns a Xpetra::Matrix object into a Belos operator
        problem -> setRightPrec( belosPrecOp );
      }

      problem -> setProblem( X, B );

      bool set = problem->setProblem();
      if (set == false) {
        *out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        return EXIT_FAILURE;
      }

      // Belos solver
      RCP< Belos::SolverManager<SC, MV, OP> > solver;
      RCP< Belos::SolverFactory<SC, MV,OP> > factory = rcp( new  Belos::SolverFactory<SC,MV,OP>() );
      RCP<Teuchos::ParameterList> belosParams
        = rcp( new Teuchos::ParameterList() );
      belosParams->set("Maximum Iterations", maxIts);
      belosParams->set("Convergence Tolerance",tol);
      belosParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosParams->set("Output Frequency",1);
      belosParams->set("Output Style",Belos::Brief);
      belosParams->set("Implicit Residual Scaling","None");
      belosParams->set("Explicit Residual Test",true);
      belosParams->set("Explicit Residual Scaling","None");

      solver = factory->create(belosSolverType,belosParams);

      comm->barrier();

      auto tm3 = TimeMonitor::getNewTimer("Maxwell: 3 - Solve");

      // set problem and solve
      solver -> setProblem( problem );
      Belos::ReturnType status = solver -> solve();
      int iters = solver -> getNumIters();
      success = (iters<50 && status == Belos::Converged);
      if (success)
        *out << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
      else
        *out << "FAILURE! Belos did not converge fast enough." << std::endl;
      tm3 = Teuchos::null;
    }
#ifdef HAVE_MUELU_STRATIMIKOS
    if (solverName == "Stratimikos") {
      auto tm4 = TimeMonitor::getNewTimer("Maxwell: 2 - Build Stratimikos solver");

      // Build the rest of the Stratimikos list
      Teuchos::ParameterList SList;
      SList.set("Linear Solver Type","Belos");
      SList.sublist("Linear Solver Types").sublist("Belos").set("Solver Type", belosSolverType);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Output Frequency",1);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Maximum Iterations",maxIts);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Convergence Tolerance",tol);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Output Style",1);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Verbosity",33);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Implicit Residual Scaling","None");
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Explicit Residual Test",true);
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist(belosSolverType).set("Explicit Residual Scaling","None");
      SList.sublist("Linear Solver Types").sublist("Belos").sublist("VerboseObject").set("Verbosity Level", "medium");
      SList.set("Preconditioner Type","MueLuRefMaxwell");
      params.set("parameterlist: syntax","muelu");
      SList.sublist("Preconditioner Types").set("MueLuRefMaxwell",params);
      // Add matrices to parameterlist
      SList.sublist("Preconditioner Types").sublist("MueLuRefMaxwell").set("D0",D0_Matrix);
      SList.sublist("Preconditioner Types").sublist("MueLuRefMaxwell").set("M0inv",M0inv_Matrix);
      SList.sublist("Preconditioner Types").sublist("MueLuRefMaxwell").set("M1",M1_Matrix);
      SList.sublist("Preconditioner Types").sublist("MueLuRefMaxwell").set("Ms",Ms_Matrix);
      SList.sublist("Preconditioner Types").sublist("MueLuRefMaxwell").set("Coordinates",coords);

      // Build Thyra linear algebra objects
      RCP<const Thyra::LinearOpBase<Scalar> > thyraA = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(SM_Matrix)->getCrsMatrix());
      RCP<      Thyra::VectorBase<Scalar> >thyraX = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraVector(X->getVectorNonConst(0)));
      // TODO: Why do we loose a reference when running this with Epetra?
      RCP<const Thyra::VectorBase<Scalar> >thyraB = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraVector(B->getVector(0));

      // Build Stratimikos solver
      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;  // This is the Stratimikos main class (= factory of solver factory).
      Stratimikos::enableMueLuRefMaxwell<LocalOrdinal,GlobalOrdinal,Node>(linearSolverBuilder);                // Register MueLu as a Stratimikos preconditioner strategy.
      linearSolverBuilder.setParameterList(rcp(&SList,false));              // Setup solver parameters using a Stratimikos parameter list.

      // Build a new "solver factory" according to the previously specified parameter list.
      RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);

      // Build a Thyra operator corresponding to A^{-1} computed using the Stratimikos solver.
      Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > thyraInverseA = Thyra::linearOpWithSolve(*solverFactory, thyraA);

      comm->barrier();

      tm4 = Teuchos::null;
      auto tm5 = TimeMonitor::getNewTimer("Maxwell: 3 - Solve");

      // Solve Ax = b.
      Thyra::SolveStatus<Scalar> status = Thyra::solve<Scalar>(*thyraInverseA, Thyra::NOTRANS, *thyraB, thyraX.ptr());
      std::cout << status << std::endl;

      success = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);
      tm5 = Teuchos::null;
    }
#endif
    comm->barrier();
    globalTimeMonitor = Teuchos::null;

    if (printTimings) {
      RCP<Teuchos::ParameterList> reportParams = rcp(new Teuchos::ParameterList);
      if (timingsFormat == "yaml") {
        reportParams->set("Report format",             "YAML");            // "Table" or "YAML"
        reportParams->set("YAML style",                "compact");         // "spacious" or "compact"
      }
      reportParams->set("How to merge timer sets",   "Union");
      reportParams->set("alwaysWriteLocal",          false);
      reportParams->set("writeGlobalStats",          true);
      reportParams->set("writeZeroTimers",           false);
      // FIXME: no "ignoreZeroTimers"

      const std::string filter = "";

      std::ios_base::fmtflags ff(out->flags());
      if (timingsFormat == "table-fixed") *out << std::fixed;
      else * out << std::scientific;
      if (use_stacked_timer) {
        stacked_timer->stop("Maxwell Driver");
        Teuchos::StackedTimer::OutputOptions options;
        options.output_fraction = options.output_histogram = options.output_minmax = true;
        stacked_timer->report(*out, comm, options);
      } else
        TimeMonitor::report(comm.ptr(), *out, filter, reportParams);
       *out << std::setiosflags(ff);
    }

    TimeMonitor::clearCounters();
#endif // #ifdef HAVE_MUELU_BELOS
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
#else
  return EXIT_SUCCESS;
#endif
} // main


template<typename Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
  return MainWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Node>::main_(clp, lib, argc, argv);
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
