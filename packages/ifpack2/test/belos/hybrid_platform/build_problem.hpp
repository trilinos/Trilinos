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

#ifndef _build_problem_hpp_
#define _build_problem_hpp_

#include <string>
#include <sstream>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_Comm.hpp>

#include <Ifpack2_Diagonal.hpp>

#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>

#include <Tpetra_MatrixIO.hpp>

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Belos::LinearProblem<Scalar,Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > > 
build_problem(Teuchos::ParameterList& test_params, 
              const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A, 
              const Teuchos::RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& b_in,
              const Teuchos::RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& dvec)
{
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>    TOP;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TCRS;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                TMap;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Belos::OperatorTraits<Scalar,TMV,TOP>                       BOPT;
  typedef Belos::MultiVecTraits<Scalar,TMV>                           BMVT;
  typedef Belos::LinearProblem<Scalar,TMV,TOP>                        BLinProb;

  Teuchos::RCP<const TMap> rowmap = A->getRowMap();

  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(rowmap, 1));
  Teuchos::RCP<const TMV> b(b_in);

  if (b == Teuchos::null) {
    Teuchos::RCP<TMV> b_nc;
    b_nc = Teuchos::rcp(new TMV(rowmap, 1, false));
    x->randomize();

    BOPT::Apply(*A, *x, *b_nc);
    BMVT::MvInit(*x, 0);
    b = b_nc;
  }
  else {
    x->putScalar(0);
  }

  Teuchos::RCP<BLinProb> problem = Teuchos::rcp(new BLinProb(A,x,b));

  std::string ifpack2_precond("not specified");
  Ifpack2::getParameter(test_params, "Ifpack2::Preconditioner", ifpack2_precond);
  if (ifpack2_precond == "DIAGONAL") {
    typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
    Teuchos::RCP<Tprec> prec;
    prec = Teuchos::rcp(new Ifpack2::Diagonal<TCRS>(dvec));
    prec->compute();
    problem->setLeftPrec(prec);
  }
  else if (ifpack2_precond != "not specified") {
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  problem->setProblem();

  return problem;
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<
    Belos::LinearProblem<
        Scalar,
        Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,
        Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
   > build_problem(Teuchos::ParameterList& test_params,
                   const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                   Teuchos::RCP<Node> node,
                   int myWeight)
{
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>      TV;
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>    TOP;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TCRS;
  typedef Belos::LinearProblem<Scalar,TMV,TOP>                        BLinProb;

  Teuchos::RCP<TCRS> A;
  Teuchos::RCP<TMV> b = Teuchos::null;
  Teuchos::RCP<TV> dvec;

  std::string mm_file("not specified");
  //std::string rhs_mm_file("not specified");
  std::string gen_mat("not specified");
  std::string hb_file("not specified");
  Ifpack2::getParameter(test_params, "mm_file", mm_file);
  Ifpack2::getParameter(test_params, "hb_file", hb_file);
  //Ifpack2::getParameter(test_params, "rhs_mm_file", rhs_mm_file);
  Ifpack2::getParameter(test_params, "gen_mat", gen_mat);

  const bool IAmRoot = (comm->getRank() == 0);

  if (hb_file != "not specified") {
    if (IAmRoot) {
      std::cout << "Harwell-Boeing file: " << hb_file << std::endl;
    }
    Tpetra::Utils::readHBMatrix(hb_file,comm,node,A);
  }
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

  Teuchos::RCP<BLinProb> problem = build_problem<Scalar,LocalOrdinal,GlobalOrdinal,Node>(test_params, A, b, dvec);

  return problem;
}


#endif

