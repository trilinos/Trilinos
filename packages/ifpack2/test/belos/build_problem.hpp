/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef _build_problem_hpp_
#define _build_problem_hpp_

#include <string>
#include <sstream>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Comm.hpp"

#include "Ifpack2_BorderedOperator.hpp"
#include "Ifpack2_Preconditioner.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"

#include "read_matrix.hpp"
#include "build_precond.hpp"

template< class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node >
Teuchos::RCP<Belos::LinearProblem<Scalar,
                                  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,
                                  Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > >
build_problem_mm (Teuchos::ParameterList& test_params,
                  const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,
                  Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& b,
                  Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& nullVec)
{
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>             TOP;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>          TMV;
  typedef Belos::OperatorTraits<Scalar,TMV,TOP>                                BOPT;
  typedef Belos::MultiVecTraits<Scalar,TMV>                                    BMVT;
  typedef Belos::LinearProblem<Scalar,TMV,TOP>                                 BLinProb;
  typedef Ifpack2::BorderedOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>    IBOP;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                         TMap;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  Teuchos::RCP<const TMap> rowmap = A->getRowMap();

  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(rowmap, 1));

  if (b == Teuchos::null) {
    b = Teuchos::rcp (new TMV (rowmap, 1));
    x->randomize ();
    BOPT::Apply (*A, *x, *b);
    BMVT::MvInit (*x, STS::zero ());
  }
  else {
    x->putScalar (STS::zero ());
  }

  Teuchos::RCP< BLinProb > problem;
  Teuchos::RCP<IBOP> borderedA;
  if (nullVec == Teuchos::null) {
     problem = Teuchos::rcp (new BLinProb (A, x, b));
  }
  else {
    borderedA = Teuchos::rcp (new IBOP (A));
    problem = Teuchos::rcp (new BLinProb (borderedA, x, b));
  }

  std::string tifpack_precond("not specified");
  Ifpack2::getParameter (test_params, "Ifpack2::Preconditioner", tifpack_precond);
  if (tifpack_precond != "not specified") {
    Teuchos::RCP<TOP> precond = build_precond<Scalar,LocalOrdinal,GlobalOrdinal,Node> (test_params, A);
    problem->setLeftPrec (precond);
  }

  problem->setProblem ();
  return problem;
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Belos::LinearProblem<
               Scalar,
               Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,
               Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > >
build_problem (Teuchos::ParameterList& test_params,
               const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
               const Teuchos::RCP<Node>& node = Teuchos::null)
{
  using Teuchos::ArrayView;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node>       crs_matrix_type;
  typedef Tpetra::Map<LO, GO, Node>                     map_type;
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node>     TMV;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  //typedef Belos::LinearProblem<Scalar, TMV, TOP>        BLinProb; // unused

  RCP<const crs_matrix_type> A;
  RCP<TMV> b = Teuchos::null;
  RCP<TMV> nullVec = Teuchos::null;

  std::string mm_file("not specified");
  std::string rhs_mm_file("not specified");
  std::string nullMvec_mm_file("not specified");
  Ifpack2::getParameter(test_params, "mm_file", mm_file);
  Ifpack2::getParameter(test_params, "rhs_mm_file", rhs_mm_file);
  std::string hb_file("not specified");
  Ifpack2::getParameter(test_params, "hb_file", hb_file);
  bool useMatrixWithConstGraph = false;
  Ifpack2::getParameter(test_params, "Use matrix with const graph", useMatrixWithConstGraph);

  if (mm_file != "not specified") {
    if (comm->getRank() == 0) {
      std::cout << "Matrix Market file for sparse matrix A: " << mm_file << std::endl;
    }
    RCP<ParameterList> constructorParams = parameterList ("CrsMatrix");
    RCP<ParameterList> fillCompleteParams = parameterList ("fillComplete");
    if (useMatrixWithConstGraph) {
      // We need to keep the local graph so that we can create a new
      // matrix with a const graph, using the graph of the original
      // matrix read in here.
      // fillCompleteParams->set ("Optimize Storage", false);
      fillCompleteParams->set ("Preserve Local Graph", true);
    }
    A = reader_type::readSparseFile (mm_file, comm, constructorParams,
                                     fillCompleteParams);

    RCP<const map_type> domainMap = A->getDomainMap ();
    RCP<const map_type> rangeMap = A->getRangeMap ();

    if (rhs_mm_file != "not specified") {
      if (comm->getRank() == 0) {
        std::cout << "Matrix Market file for right-hand-side(s) B: " << rhs_mm_file << std::endl;
      }
      b = reader_type::readDenseFile (rhs_mm_file, comm, rangeMap);
    }

    if (nullMvec_mm_file != "not specified") {
      if (comm->getRank() == 0) {
        std::cout << "Matrix Market file for null multivector: " << nullMvec_mm_file << std::endl;
      }
      // mfh 31 Jan 2013: I'm not sure what a "null multivector" means
      // in this context, so I'm only guessing that it's a domain Map
      // multivector.
      nullVec = reader_type::readDenseFile (nullMvec_mm_file, comm, domainMap);
    }

  }
  else if (hb_file != "not specified") {
    if (comm->getRank() == 0) {
      std::cout << "Harwell-Boeing file: " << hb_file << std::endl;
    }
    A = read_matrix_hb<Scalar,LO,GO,Node> (hb_file, comm);
  }
  else {
    throw std::runtime_error("No matrix file specified.");
  }

  if (useMatrixWithConstGraph) {
    // Some Ifpack2 preconditioners that extract diagonal entries have
    // a special path for doing so more efficiently when the matrix
    // has a const graph (isStaticGraph()).  In order to test this, we
    // specifically create a matrix with a const graph, by extracting
    // the original matrix's graph and copying all the values into the
    // new matrix.
    RCP<crs_matrix_type> A_constGraph (new crs_matrix_type (A->getCrsGraph ()));
    // Copy the values row by row from A into A_constGraph.
    ArrayView<const LO> ind;
    ArrayView<const Scalar> val;
    const LO numLocalRows = static_cast<LO> (A->getNodeNumRows ());
    for (LO localRow = 0; localRow < numLocalRows; ++localRow) {
      A->getLocalRowView (localRow, ind, val);
      A_constGraph->replaceLocalValues (localRow, ind, val);
    }
    A_constGraph->fillComplete (A->getDomainMap (), A->getRangeMap ());
    A = A_constGraph; // Replace A with A_constGraph.
  }

  return build_problem_mm<Scalar,LO,GO,Node> (test_params, A, b, nullVec);
}


#endif

