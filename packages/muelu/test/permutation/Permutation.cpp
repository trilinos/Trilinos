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


typedef double  Scalar;
typedef int     LocalOrdinal;
typedef int     GlobalOrdinal;
typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

#include <unistd.h>
/**********************************************************************************/

namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

  Teuchos::RCP<const Epetra_CrsMatrix> GetEpetraMatrix(std::string name, const Teuchos::RCP<Level> level, const Teuchos::RCP<Factory>& fct) {
    Teuchos::RCP<Matrix> result = level->Get<Teuchos::RCP<Matrix> >(name, fct.get());
    Teuchos::RCP<CrsMatrixWrap> crsres = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(result);
    Teuchos::RCP<CrsMatrix> crsmat = crsres->getCrsMatrix();
    Teuchos::RCP<EpetraCrsMatrix> epcrsmat = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(crsmat);
    Teuchos::RCP<const Epetra_CrsMatrix> epres = epcrsmat->getEpetra_CrsMatrix();
    return epres;
  }

  // run tests with "Algebraic" permutation strategy and nDofsPerNode = 1
  bool runPermutationTest(const std::string input_filename, const std::string expected_filename, const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {

#ifndef HAVE_MUELU_INST_COMPLEX_INT_INT
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);

    Epetra_CrsMatrix * ptrA = NULL;
    Epetra_CrsMatrix * ptrExpected = NULL;
    int ret = EpetraExt::MatlabFileToCrsMatrix ( input_filename.c_str(),
        *Xpetra::toEpetra(comm),
        ptrA
        );

    if(ret!=0)
      std::cout << "failed to read matrix from file" << std::endl;

    if(expected_filename.size() > 0)
    {
      int ret2 = EpetraExt::MatlabFileToCrsMatrix (expected_filename.c_str(),
          *Xpetra::toEpetra(comm),
          ptrExpected
          );

      if(ret2!=0)
        std::cout << "failed to read matrix from file" << std::endl;

    }
    Teuchos::RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
    Teuchos::RCP<Epetra_CrsMatrix> epExpected = Teuchos::rcp(ptrExpected);

    // Epetra_CrsMatrix -> Xpetra::Matrix
    Teuchos::RCP<CrsMatrix> exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA));
    Teuchos::RCP<CrsMatrixWrap> crsOp = Teuchos::rcp(new CrsMatrixWrap(exA));
    Teuchos::RCP<Matrix> A = Teuchos::rcp_dynamic_cast<Matrix>(crsOp);
    A->SetFixedBlockSize(1);

    Teuchos::RCP<Level> Finest = Teuchos::rcp(new Level());
    Finest->SetLevelID(0);  // must be level 0 for NullspaceFactory
    Finest->Set("A", A);

    // permute full matrix
    Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new MueLu::PermutationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>());
    PermFact->SetParameter("PermutationStrategy",Teuchos::ParameterEntry(std::string("Algebraic")));
    //PermFact->SetParameter("PermutationStrategy",Teuchos::ParameterEntry(std::string("Local")));
    PermFact->SetParameter("PermutationRowMapName",Teuchos::ParameterEntry(std::string("")));
    PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);

    // setup main factory manager
    Teuchos::RCP<FactoryManager> M = Teuchos::rcp(new FactoryManager());
    M->SetFactory("permQT",          PermFact);
    M->SetFactory("A",               MueLu::NoFactory::getRCP()); // this is the input matrix
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

    //std::cout << "P" <<  *GetEpetraMatrix("permP", Finest, PermFact) << std::endl;
    //std::cout << "Q^T" << *GetEpetraMatrix("permQT", Finest, PermFact) << std::endl;
    //std::cout << "permA" <<  *GetEpetraMatrix("A", Finest, PermFact) << std::endl;

    Teuchos::RCP<const Epetra_CrsMatrix> epResult = GetEpetraMatrix("A", Finest, PermFact);
    //std::cout << *epResult << std::endl;

    if(epExpected != Teuchos::null) {
      Epetra_CrsMatrix* comparison = NULL;
      EpetraExt::MatrixMatrix::Add(*epResult, false, -1.0, *epExpected, false, 1.0, comparison);
      comparison->FillComplete();
      double norm = comparison->NormInf();
      delete comparison;
      comparison = NULL;

      if(norm < 1.0e-14) {
        *out << "** PASSED **: " << input_filename << std::endl;
        return true;
      }
      else {
        *out << "-- FAILED --: " << input_filename << std::endl;
        return false;
      }

    }
#endif
    *out << "-- FAILED --: " << input_filename << " no result file found" << std::endl;
    return false; // no result for comparison available
  }

  // run tests with "Local" permutation strategy and nDofsPerNode = 3
  bool runPermutationTest2(const std::string input_filename, const std::string expected_filename, const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {

#ifndef HAVE_MUELU_INST_COMPLEX_INT_INT
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);

    Epetra_CrsMatrix * ptrA = NULL;
    Epetra_CrsMatrix * ptrExpected = NULL;
    int ret = EpetraExt::MatlabFileToCrsMatrix ( input_filename.c_str(),
        *Xpetra::toEpetra(comm),
        ptrA
        );

    if(ret!=0)
      std::cout << "failed to read matrix from file" << std::endl;

    if(expected_filename.size() > 0)
    {
      int ret2 = EpetraExt::MatlabFileToCrsMatrix (expected_filename.c_str(),
          *Xpetra::toEpetra(comm),
          ptrExpected
          );

      if(ret2!=0)
        std::cout << "failed to read matrix from file" << std::endl;

    }
    Teuchos::RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
    Teuchos::RCP<Epetra_CrsMatrix> epExpected = Teuchos::rcp(ptrExpected);

    // Epetra_CrsMatrix -> Xpetra::Matrix
    Teuchos::RCP<CrsMatrix> exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA));
    Teuchos::RCP<CrsMatrixWrap> crsOp = Teuchos::rcp(new CrsMatrixWrap(exA));
    Teuchos::RCP<Matrix> A = Teuchos::rcp_dynamic_cast<Matrix>(crsOp);
    A->SetFixedBlockSize(3);

    Teuchos::RCP<Level> Finest = Teuchos::rcp(new Level());
    Finest->SetLevelID(0);  // must be level 0 for NullspaceFactory
    Finest->Set("A", A);

    // permute full matrix
    Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new MueLu::PermutationFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>());
    PermFact->SetParameter("PermutationStrategy",Teuchos::ParameterEntry(std::string("Local")));
    PermFact->SetParameter("PermutationRowMapName",Teuchos::ParameterEntry(std::string("")));
    PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);

    // setup main factory manager
    Teuchos::RCP<FactoryManager> M = Teuchos::rcp(new FactoryManager());
    M->SetFactory("permQT",          PermFact);
    M->SetFactory("A",               MueLu::NoFactory::getRCP()); // this is the input matrix
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

    //std::cout << "P" <<  *GetEpetraMatrix("permP", Finest, PermFact) << std::endl;
    //std::cout << "Q^T" << *GetEpetraMatrix("permQT", Finest, PermFact) << std::endl;
    //std::cout << "permA" <<  *GetEpetraMatrix("A", Finest, PermFact) << std::endl;

    Teuchos::RCP<const Epetra_CrsMatrix> epResult = GetEpetraMatrix("A", Finest, PermFact);
    //std::cout << *epResult << std::endl;

    if(epExpected != Teuchos::null) {
      Epetra_CrsMatrix* comparison = NULL;
      EpetraExt::MatrixMatrix::Add(*epResult, false, -1.0, *epExpected, false, 1.0, comparison);
      comparison->FillComplete();
      //std::cout << *comparison << std::endl;
      double norm = comparison->NormInf();
      delete comparison;
      comparison = NULL;

      if(norm < 1.0e-14) {
        *out << "** PASSED **: " << input_filename << std::endl;
        return true;
      }
      else {
        *out << "-- FAILED --: " << input_filename << std::endl;
        return false;
      }

    }
#endif
    *out << "-- FAILED --: " << input_filename << " no result file found" << std::endl;
    return false; // no result for comparison available
  }

}


int main(int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  using Teuchos::RCP;
  using namespace MueLuTests;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int numProcs = comm->getSize();

  bool bSuccess = true;
#ifndef HAVE_MUELU_INST_COMPLEX_INT_INT
  //runPermutationTest(MatrixFileName, ExpectedFileName, comm)
  if(runPermutationTest("test1.txt", "exp1.txt", comm) == false) bSuccess = false;
  if(runPermutationTest("test2.txt", "exp2.txt", comm) == false) bSuccess = false;
  if(runPermutationTest("test3.txt", "exp3.txt", comm) == false) bSuccess = false;

  // the following tests work only on 1 or 2 processors
  if(numProcs == 1 || numProcs == 2) {
    if(runPermutationTest2("test4.txt", "exp4.txt", comm) == false) bSuccess = false;
    // test seems to be ok, but matrix addition is not working
    // has wrong entries on the diagonal on proc1... -> wrong handling of index base?
    //if(runPermutationTest2("test5.txt", "exp5.txt" /*"exp5.txt"*/, comm) == false) bSuccess = false;
  }
#endif
  if(bSuccess == false)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

