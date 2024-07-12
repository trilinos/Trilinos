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

#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_MatrixMatrix.h>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Parameters.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_PermutationFactory.hpp"

#include "MueLu_Exceptions.hpp"

#include <unistd.h>
/**********************************************************************************/

namespace MueLuTests {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Epetra_CrsMatrix> GetEpetraMatrix(std::string name, const Teuchos::RCP<MueLu::Level> level, const Teuchos::RCP<MueLu::Factory>& fct) {
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > result        = level->Get<Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >(name, fct.get());
  Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > crsres = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(result);
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > crsmat     = crsres->getCrsMatrix();
  Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node> > epcrsmat                  = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node> >(crsmat);
  Teuchos::RCP<const Epetra_CrsMatrix> epres                                             = epcrsmat->getEpetra_CrsMatrix();
  return epres;
}

// run tests with "Algebraic" permutation strategy and nDofsPerNode = 1
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool runPermutationTest(const std::string input_filename, const std::string expected_filename, const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
#include <MueLu_UseShortNames.hpp>

#ifndef HAVE_MUELU_INST_COMPLEX_INT_INT
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

  Epetra_CrsMatrix* ptrA        = NULL;
  Epetra_CrsMatrix* ptrExpected = NULL;
  int ret                       = EpetraExt::MatlabFileToCrsMatrix(input_filename.c_str(),
                                                                   *Xpetra::toEpetra(comm),
                                                                   ptrA);

  if (ret != 0)
    std::cout << "failed to read matrix from file" << std::endl;

  if (expected_filename.size() > 0) {
    int ret2 = EpetraExt::MatlabFileToCrsMatrix(expected_filename.c_str(),
                                                *Xpetra::toEpetra(comm),
                                                ptrExpected);

    if (ret2 != 0)
      std::cout << "failed to read matrix from file" << std::endl;
  }
  Teuchos::RCP<Epetra_CrsMatrix> epA        = Teuchos::rcp(ptrA);
  Teuchos::RCP<Epetra_CrsMatrix> epExpected = Teuchos::rcp(ptrExpected);

  // Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<CrsMatrix> exA       = Teuchos::rcp(new EpetraCrsMatrix(epA));
  Teuchos::RCP<CrsMatrixWrap> crsOp = Teuchos::rcp(new CrsMatrixWrap(exA));
  Teuchos::RCP<Matrix> A            = Teuchos::rcp_dynamic_cast<Matrix>(crsOp);
  A->SetFixedBlockSize(1);

  Teuchos::RCP<Level> Finest = Teuchos::rcp(new Level());
  Finest->SetLevelID(0);  // must be level 0 for NullspaceFactory
  Finest->Set("A", A);

  // permute full matrix
  Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new MueLu::PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
  PermFact->SetParameter("PermutationStrategy", Teuchos::ParameterEntry(std::string("Algebraic")));
  // PermFact->SetParameter("PermutationStrategy",Teuchos::ParameterEntry(std::string("Local")));
  PermFact->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
  PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);

  // setup main factory manager
  Teuchos::RCP<FactoryManager> M = Teuchos::rcp(new FactoryManager());
  M->SetFactory("permQT", PermFact);
  M->SetFactory("A", MueLu::NoFactory::getRCP());  // this is the input matrix
  MueLu::SetFactoryManager SFMFinest(Finest, M);

  // prepare building process for permutation operators
  Finest->Request("A", PermFact.get());
  Finest->Request("permA", PermFact.get());
  Finest->Request("permP", PermFact.get());
  Finest->Request("permQT", PermFact.get());
  Finest->Request("permScaling", PermFact.get());
  Finest->Request("#RowPermutations", PermFact.get());
  Finest->Request("#ColPermutations", PermFact.get());
  Finest->Request("#WideRangeRowPermutations", PermFact.get());
  Finest->Request("#WideRangeColPermutations", PermFact.get());

  // build permutation operators
  PermFact->Build(*Finest);

  // std::cout << "P" <<  *GetEpetraMatrix("permP", Finest, PermFact) << std::endl;
  // std::cout << "Q^T" << *GetEpetraMatrix("permQT", Finest, PermFact) << std::endl;
  // std::cout << "permA" <<  *GetEpetraMatrix("A", Finest, PermFact) << std::endl;

  Teuchos::RCP<const Epetra_CrsMatrix> epResult = GetEpetraMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>("A", Finest, PermFact);
  // std::cout << *epResult << std::endl;

  if (epExpected != Teuchos::null) {
    Epetra_CrsMatrix* comparison = NULL;
    EpetraExt::MatrixMatrix::Add(*epResult, false, -1.0, *epExpected, false, 1.0, comparison);
    comparison->FillComplete();
    double norm = comparison->NormInf();
    delete comparison;
    comparison = NULL;

    if (norm < 1.0e-14) {
      *out << "** PASSED **: " << input_filename << std::endl;
      return true;
    } else {
      *out << "-- FAILED --: " << input_filename << std::endl;
      return false;
    }
  }
#endif
  *out << "-- FAILED --: " << input_filename << " no result file found" << std::endl;
  return false;  // no result for comparison available
}

// run tests with "Local" permutation strategy and nDofsPerNode = 3
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool runPermutationTest2(const std::string input_filename, const std::string expected_filename, const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
#include <MueLu_UseShortNames.hpp>

#ifndef HAVE_MUELU_INST_COMPLEX_INT_INT
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

  Epetra_CrsMatrix* ptrA        = NULL;
  Epetra_CrsMatrix* ptrExpected = NULL;
  int ret                       = EpetraExt::MatlabFileToCrsMatrix(input_filename.c_str(),
                                                                   *Xpetra::toEpetra(comm),
                                                                   ptrA);

  if (ret != 0)
    std::cout << "failed to read matrix from file" << std::endl;

  if (expected_filename.size() > 0) {
    int ret2 = EpetraExt::MatlabFileToCrsMatrix(expected_filename.c_str(),
                                                *Xpetra::toEpetra(comm),
                                                ptrExpected);

    if (ret2 != 0)
      std::cout << "failed to read matrix from file" << std::endl;
  }
  Teuchos::RCP<Epetra_CrsMatrix> epA        = Teuchos::rcp(ptrA);
  Teuchos::RCP<Epetra_CrsMatrix> epExpected = Teuchos::rcp(ptrExpected);

  // Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<CrsMatrix> exA       = Teuchos::rcp(new EpetraCrsMatrix(epA));
  Teuchos::RCP<CrsMatrixWrap> crsOp = Teuchos::rcp(new CrsMatrixWrap(exA));
  Teuchos::RCP<Matrix> A            = Teuchos::rcp_dynamic_cast<Matrix>(crsOp);
  A->SetFixedBlockSize(3);

  Teuchos::RCP<Level> Finest = Teuchos::rcp(new Level());
  Finest->SetLevelID(0);  // must be level 0 for NullspaceFactory
  Finest->Set("A", A);

  // permute full matrix
  Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new MueLu::PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
  PermFact->SetParameter("PermutationStrategy", Teuchos::ParameterEntry(std::string("Local")));
  PermFact->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
  PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);

  // setup main factory manager
  Teuchos::RCP<FactoryManager> M = Teuchos::rcp(new FactoryManager());
  M->SetFactory("permQT", PermFact);
  M->SetFactory("A", MueLu::NoFactory::getRCP());  // this is the input matrix
  MueLu::SetFactoryManager SFMFinest(Finest, M);

  // prepare building process for permutation operators
  Finest->Request("A", PermFact.get());
  Finest->Request("permA", PermFact.get());
  Finest->Request("permP", PermFact.get());
  Finest->Request("permQT", PermFact.get());
  Finest->Request("permScaling", PermFact.get());
  Finest->Request("#RowPermutations", PermFact.get());
  Finest->Request("#ColPermutations", PermFact.get());
  Finest->Request("#WideRangeRowPermutations", PermFact.get());
  Finest->Request("#WideRangeColPermutations", PermFact.get());

  // build permutation operators
  PermFact->Build(*Finest);

  // std::cout << "P" <<  *GetEpetraMatrix("permP", Finest, PermFact) << std::endl;
  // std::cout << "Q^T" << *GetEpetraMatrix("permQT", Finest, PermFact) << std::endl;
  // std::cout << "permA" <<  *GetEpetraMatrix("A", Finest, PermFact) << std::endl;

  Teuchos::RCP<const Epetra_CrsMatrix> epResult = GetEpetraMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>("A", Finest, PermFact);
  // std::cout << *epResult << std::endl;

  if (epExpected != Teuchos::null) {
    Epetra_CrsMatrix* comparison = NULL;
    EpetraExt::MatrixMatrix::Add(*epResult, false, -1.0, *epExpected, false, 1.0, comparison);
    comparison->FillComplete();
    // std::cout << *comparison << std::endl;
    double norm = comparison->NormInf();
    delete comparison;
    comparison = NULL;

    if (norm < 1.0e-14) {
      *out << "** PASSED **: " << input_filename << std::endl;
      return true;
    } else {
      *out << "-- FAILED --: " << input_filename << std::endl;
      return false;
    }
  }
#endif
  *out << "-- FAILED --: " << input_filename << " no result file found" << std::endl;
  return false;  // no result for comparison available
}

}  // namespace MueLuTests

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using namespace MueLuTests;

  Teuchos::oblackholestream blackhole;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int numProcs                        = comm->getSize();

  bool bSuccess = true;
#ifndef HAVE_MUELU_INST_COMPLEX_INT_INT
  // runPermutationTest(MatrixFileName, ExpectedFileName, comm)
  if (runPermutationTest<Scalar, LocalOrdinal, GlobalOrdinal, Node>("test1.txt", "exp1.txt", comm) == false) bSuccess = false;
  if (runPermutationTest<Scalar, LocalOrdinal, GlobalOrdinal, Node>("test2.txt", "exp2.txt", comm) == false) bSuccess = false;
  if (runPermutationTest<Scalar, LocalOrdinal, GlobalOrdinal, Node>("test3.txt", "exp3.txt", comm) == false) bSuccess = false;

  // the following tests work only on 1 or 2 processors
  if (numProcs == 1 || numProcs == 2) {
    if (runPermutationTest2<Scalar, LocalOrdinal, GlobalOrdinal, Node>("test4.txt", "exp4.txt", comm) == false) bSuccess = false;
    // test seems to be ok, but matrix addition is not working
    // has wrong entries on the diagonal on proc1... -> wrong handling of index base?
    // if(runPermutationTest2("test5.txt", "exp5.txt" /*"exp5.txt"*/, comm) == false) bSuccess = false;
  }
#endif
  if (bSuccess == false)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

int main(int argc, char* argv[]) {
  bool success = false;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  const bool throwExceptions     = false;
  const bool recogniseAllOptions = false;

  Teuchos::CommandLineProcessor clp(throwExceptions, recogniseAllOptions);
  Xpetra::Parameters xpetraParameters(clp);

  std::string node = "";
  clp.setOption("node", &node, "node type (serial | openmp | cuda | hip)");

  switch (clp.parse(argc, argv, NULL)) {
    case Teuchos::CommandLineProcessor::PARSE_ERROR: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
    return main_<double, int, int, Xpetra::EpetraNode>(clp, lib, argc, argv);
#else
    throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
  }

  if (lib == Xpetra::UseTpetra) {
    std::cout << "Skip permutation tests for Tpetra." << std::endl;
  }

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
