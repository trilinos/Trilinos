// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <MueLu_Exceptions.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

#include "MueLu_SemiCoarsenPFactory.hpp"  // for semi-coarsening constants
#include <MueLu_TestHelpers.hpp>
#include <MueLu_MultiPhys.hpp>

/**********************************************************************************/
/* CREATE INITAL MATRIX                                                           */
/**********************************************************************************/
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#ifdef HAVE_MUELU_TPETRA
#include <BelosTpetraAdapter.hpp>
#endif
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp
#endif
#define FOURBYFOUR

#include <unistd.h>
/**********************************************************************************/

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

#ifdef BIGVERSION
  std::string suffix("big");
#else
  std::string suffix("");
#endif
  Teuchos::oblackholestream blackhole;

  bool success = true;
  bool verbose = true;
  try {
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

    // Four different versions of essentially the same problem can be run.
    // depending on whether or not FOURBYFOUR is defined and whether or not
    // BIGVERSION is defined  BIGVERSION just uses a finer mesh than the
    // standard version. Data for the big and the 2x2 versions is kept
    // in https://gitlab-ex.sandia.gov/muelu/mueludata in multiphys directory
    // In the description below we use terms like
    // filename[big].mat.  When BIGVERSION is defined, this refers to
    // filenamebig.mat while it refers to filename.mat when BIGVERSION is not
    // defined.  FOURBYFOUR decides whether or not a 4x4 or a 2x2 block system
    // is solved.
    //
    // The 2 x 2 block system
    //
    // [ A  B ; C  D]
    //
    // Presently, the file multiPhys2x2[big].mat is read to define the entries of A, B, C, D.
    // A is supposed to be a 2D elasticity problem while D is a Laplace operator. In particular,
    // the map for D has one DoF per node while the map for A has 2 DoFs per node. Currently,
    // the B and C associated with multiPhys2x2[big].mat are both zero. The larger version
    // of this problem (when FOURBYFOUR is defined) just replicates this system
    //
    //             [ A  B 0  0 ; C  D 0 0 ;  0 0 A  B ; 0  0 C  D ]
    //
    // and the file read is multiPhys4x4[big].mat.
    //
    // AMG is first applied twice for the regular problem. Specifically, the
    // file aux1[big].mat defines an auxiliary system that AMG is applied to.
    // This auxiliary matrix should have the same dimensions as A. The file
    // aux2[big].mat defines a 2nd auxiliary system that AMG is applied to.
    // This auxiliary matrix should have the same dimensions as D. The grid
    // transfer for the 2x2 matrix is a block diagonal matrix
    //
    //       [  P1  0  ; 0  P2 ]
    //
    // where P1 is pulled out of the hierarchy for aux1[big].mat and P2 is pulled
    // out of the hieararchy for aux2[big].mat (on each level). For the 4x4
    // version, 4 hierarchies are built and a 4x4 block diagonal prolongator
    // is constructed.
    //
    // Additionally, the files coords*[big].mat and null*[big].mat are ready to
    // define the auxiliary coordinates and null spaces needed by the AMG
    // algorithm. For our specific problem, the meshes for the elasticity
    // and the Poisson operator partially overlap (as this is similar
    // to what happens on some ARIA ablation problems).

    int nElas, nPois;

    // read in sizes of elasticity and Poisson problem

    FILE *fp = fopen(std::string("aux1" + suffix + ".mat").c_str(), "r");
    if (fp == NULL) {
      std::cout << "\nERROR:  File aux1" + suffix + ".mat not found" << std::endl;
      return EXIT_FAILURE;
    }

    int ch;
    int ret_val = 0;
    (void)ret_val;  // suppress fscanf return value and unused variable warnings
    while ((ch = getc(fp)) != '\n')
      ;
    ret_val = fscanf(fp, "%d", &nElas);

    fp = fopen(std::string("aux2" + suffix + ".mat").c_str(), "r");
    if (fp == NULL) {
      std::cout << "\nERROR:  File aux2" + suffix + ".mat not found" << std::endl;
      return EXIT_FAILURE;
    }
    while ((ch = getc(fp)) != '\n')
      ;
    ret_val = fscanf(fp, "%d", &nPois);
    fclose(fp);

    // check that multiphysics problem and aux matrices are consistent with respect to sizes

    int itemp;

#ifdef FOURBYFOUR
    fp = fopen(std::string("multiPhys4x4" + suffix + ".mat").c_str(), "r");
    if (fp == NULL) {
      std::cout << "\nERROR:  File multiPhys4x4" + suffix + ".mat not found" << std::endl;
      return EXIT_FAILURE;
    }
    while ((ch = getc(fp)) != '\n')
      ;
    ret_val = fscanf(fp, "%d", &itemp);
    if (itemp != nElas + nPois + nElas + nPois) {
      std::cout << "\nERROR:  multiPhys4x4" + suffix + ".mat dimension is" << itemp << "was expecting it to be " << 2 * (nElas + nPois) << std::endl;
      return EXIT_FAILURE;
    }
#else
    fp = fopen(std::string("multiPhys2x2" + suffix + ".mat").c_str(), "r");
    if (fp == NULL) {
      std::cout << "\nERROR:  File multiPhys2x2" + suffix + ".mat not found" << std::endl;
      return EXIT_FAILURE;
    }
    while ((ch = getc(fp)) != '\n')
      ;
    ret_val = fscanf(fp, "%d", &itemp);
    if (itemp != nElas + nPois) {
      std::cout << "\nERROR:  multiPhys2x2" + suffix + ".mat dimension is" << itemp << "was expecting it to be " << (nElas + nPois) << std::endl;
      return EXIT_FAILURE;
    }
#endif

    Teuchos::RCP<const Map> mapPois            = MapFactory::Build(lib, nPois, 0, comm);
    Teuchos::RCP<const Map> mapElas1DofPerNode = MapFactory::Build(lib, nElas / 2, 0, comm);
    Teuchos::RCP<const Map> mapElas2DofPerNode = Xpetra::MapFactory<LO, GO, Node>::Build(mapElas1DofPerNode, 2);

    Teuchos::ArrayView<const GO> GIDsElas2DofPerNode = mapElas2DofPerNode->getLocalElementList();
    Teuchos::ArrayView<const GO> GIDsPois            = mapPois->getLocalElementList();

#ifdef FOURBYFOUR
    Teuchos::Array<GO> GIDsCombo(2 * (GIDsElas2DofPerNode.size() + GIDsPois.size()));
#else
    Teuchos::Array<GO> GIDsCombo(GIDsElas2DofPerNode.size() + GIDsPois.size());
#endif

    for (int ii = 0; ii < GIDsElas2DofPerNode.size(); ii++) GIDsCombo[ii] = GIDsElas2DofPerNode[ii];
    for (int ii = 0; ii < GIDsPois.size(); ii++) GIDsCombo[ii + GIDsElas2DofPerNode.size()] = GIDsPois[ii] + nElas;
#ifdef FOURBYFOUR
    for (int ii = 0; ii < GIDsElas2DofPerNode.size(); ii++) GIDsCombo[ii + GIDsElas2DofPerNode.size() + GIDsPois.size()] = GIDsElas2DofPerNode[ii] + nElas + nPois;
    for (int ii = 0; ii < GIDsPois.size(); ii++) GIDsCombo[ii + 2 * GIDsElas2DofPerNode.size() + GIDsPois.size()] = GIDsPois[ii] + 2 * nElas + nPois;
#endif

    Teuchos::RCP<const Map> mapMultiPhysA = Xpetra::MapFactory<LO, GO, Node>::Build(lib, Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), GIDsCombo, mapElas2DofPerNode->getIndexBase(), comm);
#ifdef FOURBYFOUR
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> multiPhysA = Xpetra::IO<SC, LO, GO, Node>::Read("multiPhys4x4" + suffix + ".mat", mapMultiPhysA);
#else
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> multiPhysA = Xpetra::IO<SC, LO, GO, Node>::Read("multiPhys2x2" + suffix + ".mat", mapMultiPhysA);
#endif

    // To solve the block system a block diagonal auxiliary operator such as
    //
    //     [ aux1    0     0      0 ]
    //     [ 0    aux1     0      0 ]
    //     [ 0       0  aux1      0 ]
    //     [ 0       0     0   aux2 ]
    //
    //
    // is used in conjunction with multiple invocations of Hierarchy Setup

    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> aux1 = Xpetra::IO<SC, LO, GO, Node>::Read("aux1" + suffix + ".mat", mapElas2DofPerNode);
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> aux2 = Xpetra::IO<SC, LO, GO, Node>::Read("aux2" + suffix + ".mat", mapPois);

    typedef Teuchos::ScalarTraits<SC> STS;
    typedef typename STS::magnitudeType real_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coordinate_type;
    typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;
    typedef Xpetra::MultiVectorFactory<coordinate_type, LO, GO, NO> RealValuedMultiVectorFactory;

    Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>> elasCoords = Xpetra::IO<real_type, LO, GO, Node>::ReadMultiVector("coords1" + suffix + ".mat", mapElas1DofPerNode);
    Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>> poisCoords = Xpetra::IO<real_type, LO, GO, Node>::ReadMultiVector("coords2" + suffix + ".mat", mapPois);

    Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> nullspace1 = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector("null1" + suffix + ".mat", mapElas2DofPerNode);
    Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> nullspace2 = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector("null2" + suffix + ".mat", mapPois);

    // Note: Entries of arrayOfNullspaces[] may be set to null.
#ifdef FOURBYFOUR
    Teuchos::ArrayRCP<Teuchos::RCP<Matrix>> arrayOfAuxMatrices(4);
    Teuchos::ArrayRCP<Teuchos::RCP<RealValuedMultiVector>> arrayOfCoords(4);
    Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfNullspaces(4);

    arrayOfAuxMatrices[2] = aux1;
    arrayOfAuxMatrices[3] = aux2;
    arrayOfCoords[2]      = elasCoords;
    arrayOfCoords[3]      = poisCoords;
    arrayOfNullspaces[2]  = nullspace1;
    arrayOfNullspaces[3]  = nullspace2;
#else
    Teuchos::ArrayRCP<Teuchos::RCP<Matrix>> arrayOfAuxMatrices(2);
    Teuchos::ArrayRCP<Teuchos::RCP<RealValuedMultiVector>> arrayOfCoords(2);
    Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>> arrayOfNullspaces(2);
#endif

    arrayOfAuxMatrices[0] = aux1;
    arrayOfAuxMatrices[1] = aux2;
    arrayOfCoords[0]      = elasCoords;
    arrayOfCoords[1]      = poisCoords;
    arrayOfNullspaces[0]  = nullspace1;
    arrayOfNullspaces[1]  = nullspace2;

    Teuchos::ParameterList comboList;
    Teuchos::updateParametersFromXmlFileAndBroadcast("comboP.xml", Teuchos::Ptr<Teuchos::ParameterList>(&comboList), *comm);
#ifdef FOURBYFOUR
    Teuchos::RCP<Operator> preconditioner = rcp(new MueLu::MultiPhys<SC, LO, GO, NO>(multiPhysA, arrayOfAuxMatrices, arrayOfNullspaces, arrayOfCoords, 4, comboList, true));
#else
    Teuchos::RCP<Operator> preconditioner = rcp(new MueLu::MultiPhys<SC, LO, GO, NO>(multiPhysA, arrayOfAuxMatrices, arrayOfNullspaces, arrayOfCoords, 2, comboList, true));
#endif

    Teuchos::RCP<Teuchos::TimeMonitor> globalTimeMonitor = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Timings: Global Time")));

    // =========================================================================
    // System solution (Ax = b)
    // =========================================================================
    typedef Teuchos::ScalarTraits<SC> STS;
    SC zero = STS::zero(), one = STS::one();

    Teuchos::RCP<Vector> X = VectorFactory::Build(multiPhysA->getRowMap());
    Teuchos::RCP<Vector> B = VectorFactory::Build(multiPhysA->getRowMap());

    {
      // set seed for reproducibility
      Utilities::SetRandomSeed(*comm);
      X->randomize();
      multiPhysA->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

      Teuchos::Array<typename STS::magnitudeType> norms(1);
      B->norm2(norms);
      B->scale(one / norms[0]);
      X->putScalar(zero);
    }

    // Belos linear problem
    typedef MultiVector MV;
    typedef Belos::OperatorT<MV> OP;

    Teuchos::RCP<OP> belosOp                                    = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(multiPhysA));
    Teuchos::RCP<OP> belosPrecOp                                = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(preconditioner));  // Turns a Xpetra::Matrix object into a Belos operator
    Teuchos::RCP<Belos::LinearProblem<SC, MV, OP>> belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
    belosProblem->setRightPrec(belosPrecOp);
    bool set = belosProblem->setProblem();
    if (set == false) {
      std::cout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }
    Teuchos::RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
    belosList->set("Maximum Iterations", 100);
    belosList->set("Convergence Tolerance", 1e-6);
    belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList->set("Output Frequency", 1);
    belosList->set("Output Style", Belos::Brief);
    belosList->set("Implicit Residual Scaling", "None");

    Belos::SolverFactory<SC, MV, OP> solverFactory;
    Teuchos::RCP<Belos::SolverManager<SC, MV, OP>> solver = solverFactory.create("gmres", belosList);
    solver->setProblem(belosProblem);
    Belos::ReturnType retStatus = Belos::Unconverged;
    retStatus                   = solver->solve();

    int iters = solver->getNumIters();
    success   = (iters < 50 && retStatus == Belos::Converged);
    if (comm->getRank() == 0) {
      if (success)
        std::cout << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
      else
        std::cout << "FAILURE! Belos did not converge fast enough." << std::endl;
    }

    // Timer final summaries
    globalTimeMonitor = Teuchos::null;  // stop this timer before summary
    Teuchos::TimeMonitor::summarize();
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
