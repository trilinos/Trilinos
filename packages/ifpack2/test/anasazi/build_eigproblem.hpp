/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER
*/

#ifndef _build_eigproblem_hpp_
#define _build_eigproblem_hpp_

#include <string>
#include <sstream>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_Comm.hpp>

#include <Ifpack2_Diagonal.hpp>

#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziTpetraAdapter.hpp>

#include <Tpetra_MatrixIO.hpp>
#include "read_matrix.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Anasazi::Eigenproblem<Scalar,Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > > 
build_eigproblem_mm(Teuchos::ParameterList& test_params, 
                 const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A, 
                 const Teuchos::RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& dvec)
{
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>    TOP;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TCRS;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                TMap;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>       TV;
  typedef Anasazi::OperatorTraits<Scalar,TMV,TOP>                    AOPT;
  typedef Anasazi::MultiVecTraits<Scalar,TMV>                        AMVT;
  typedef Anasazi::BasicEigenproblem<Scalar,TMV,TOP>             AEigProb;

  Teuchos::RCP<const TMap> rowmap = A->getRowMap();

  Teuchos::ParameterList aparams;
  if (test_params.isSublist("Anasazi")) {
    aparams = test_params.sublist("Anasazi");
  }
  int blockSize = aparams.get<int>("Block Size",1);
  int NEV       = aparams.get<int>("NEV",1);
  Teuchos::RCP<TMV> X0 = Teuchos::rcp(new TMV(rowmap, blockSize));
  X0->randomize();

  Teuchos::RCP<AEigProb> problem = Teuchos::rcp(new AEigProb(A,X0));

  std::string ifpack2_precond("not specified");
  Ifpack2::getParameter(test_params, "Ifpack2::Preconditioner", ifpack2_precond);
  if (ifpack2_precond == "DIAGONAL") {
    typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
    Teuchos::RCP<Tprec> prec;
    prec = Teuchos::rcp(new Ifpack2::Diagonal<TCRS>(dvec));
    prec->compute();
    problem->setPrec(prec);
  }
  else if (ifpack2_precond != "not specified") {
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  problem->setHermitian(false);
  problem->setNEV(NEV);
  TEUCHOS_TEST_FOR_EXCEPTION( problem->setProblem() != true, std::runtime_error, "build_eigproblem_mm(): error returned finalizing BasicEigenproblem.");

  return problem;
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<
    Anasazi::Eigenproblem<
        Scalar,
        Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,
        Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
   > build_eigproblem(Teuchos::ParameterList& test_params,
                   const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                   Teuchos::RCP<Node> node,
                   int myWeight)
{
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>      TV;
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>    TOP;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TCRS;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                TMap;
  typedef Anasazi::Eigenproblem<Scalar,TMV,TOP>                   AEigProb;

  Teuchos::RCP<TCRS> A;
  Teuchos::RCP<TV> dvec;

  std::string mm_file("not specified");
  std::string gen_mat("not specified");
  // std::string rhs_mm_file("not specified");
  // std::string hb_file("not specified");
  Ifpack2::getParameter(test_params, "mm_file", mm_file);
  int mm_file_maxnnz = test_params.get<int>("mm_file_maxnnz", -1);
  // Ifpack2::getParameter(test_params, "hb_file", hb_file);
  // Ifpack2::getParameter(test_params, "rhs_mm_file", rhs_mm_file);
  Ifpack2::getParameter(test_params, "gen_mat", gen_mat);

  const bool IAmRoot = (comm->getRank() == 0);

  if (mm_file != "not specified") {
    if (IAmRoot) {
      std::cout << "Matrix-Market file: " << mm_file << std::endl;
    }
    read_matrix_mm<Scalar,LocalOrdinal,GlobalOrdinal,Node>(mm_file, comm, node, myWeight, mm_file_maxnnz, A, dvec);
  }
  // else if (hb_file != "not specified") {
  //   if (IAmRoot) {
  //     std::cout << "Harwell-Boeing file: " << hb_file << std::endl;
  //   }
  //   A = read_matrix_hb<Scalar,LocalOrdinal,GlobalOrdinal,Node>(hb_file, comm, node, myWeight);
  // }
  else if (gen_mat != "not specified") {
    Teuchos::RCP<Teuchos::ParameterList> genparams = sublist(rcpFromRef(test_params),"GenMat");
    if (genparams == Teuchos::null) {
      throw std::runtime_error("Missing parameter list \"GenMat\".");
    }
    Teuchos::RCP<TCRS> ncA;
    Tpetra::Utils::generateMatrix(genparams,comm,node,ncA);
    A = ncA;
    dvec = rcp(new TV(A->getRowMap(), false));
    dvec->putScalar(6);
  }
  else {
    throw std::runtime_error("No matrix specified.");
  }

  // print node details
  {
    Teuchos::RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer = A->getGraph()->getImporter();
    for (int i=0; i<comm->getSize(); ++i) {
      if (comm->getRank() == i) {
        if (importer != Teuchos::null) {
          if (IAmRoot) {
            std::cout << "\nPartitioning details" << std::endl;
          }
          std::cout << "node: " << i 
            << "   same: " << importer->getNumSameIDs() 
            << "   permute: " << importer->getNumPermuteIDs()
            << "   remote: " << importer->getNumRemoteIDs()
            << "   export: " << importer->getNumExportIDs() << std::endl;
        }
      }
      comm->barrier();
      comm->barrier();
      comm->barrier();
    }
    if (IAmRoot) {
      std::cout << std::endl;
    }
  }

  Teuchos::RCP<AEigProb> problem = build_eigproblem_mm<Scalar,LocalOrdinal,GlobalOrdinal,Node>(test_params, A, dvec);

  return problem;
}


#endif

