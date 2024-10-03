// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
# include "Zoltan2_PartitioningProblem.hpp"
# include "Zoltan2_XpetraCrsMatrixAdapter.hpp"
# include "Zoltan2_XpetraMultiVectorAdapter.hpp"
#endif

#include "read_matrix.hpp"
#include "build_precond.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Belos::LinearProblem<
               Scalar,
               Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,
               Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > >
build_problem (Teuchos::ParameterList& test_params,
               const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
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
  typedef Tpetra::Operator<Scalar,LO,GO,Node>           TOP;
  typedef Belos::OperatorTraits<Scalar,TMV,TOP>         BOPT;
  typedef Belos::MultiVecTraits<Scalar,TMV>             BMVT;
  typedef Belos::LinearProblem<Scalar,TMV,TOP>          BLinProb;
  typedef Ifpack2::BorderedOperator<Scalar,LO,GO,Node>  IBOP;
  typedef Teuchos::ScalarTraits<Scalar>                 STS;

  RCP<const crs_matrix_type> A;
  RCP<TMV> b = Teuchos::null;
  RCP<TMV> nullVec = Teuchos::null;

  std::string mm_file("not specified");
  std::string map_mm_file("not specified");
  std::string rhs_mm_file("not specified");
  std::string nullMvec_mm_file("not specified");
  Ifpack2::getParameter(test_params, "mm_file", mm_file);
  Ifpack2::getParameter(test_params, "map_mm_file", map_mm_file);
  Ifpack2::getParameter(test_params, "rhs_mm_file", rhs_mm_file);
  std::string hb_file("not specified");
  Ifpack2::getParameter(test_params, "hb_file", hb_file);
  bool useMatrixWithConstGraph = false;
  bool useZoltan2 = false;
  bool useParMETIS = false;
  Ifpack2::getParameter(test_params, "Use matrix with const graph", useMatrixWithConstGraph);
  Ifpack2::getParameter(test_params, "Use Zoltan2",  useZoltan2);
  Ifpack2::getParameter(test_params, "Use ParMetis", useParMETIS);

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
    RCP<const map_type> rowMap;
    RCP<const map_type> colMap = Teuchos::null;
    if (map_mm_file != "not specified") {
      if (comm->getRank() == 0) {
        std::cout << "Matrix Market file for row Map of the sparse matrix A: " << map_mm_file << std::endl;
      }
      rowMap = reader_type::readMapFile(map_mm_file, comm);
      A = reader_type::readSparseFile (mm_file, rowMap, colMap, rowMap, rowMap);
    }
    else {
      A = reader_type::readSparseFile (mm_file, comm, constructorParams,
                                     fillCompleteParams);
    }

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
    using lids_type = typename crs_matrix_type::local_inds_host_view_type;
    using vals_type = typename crs_matrix_type::values_host_view_type;
    lids_type ind;
    vals_type val;
    const LO numLocalRows = static_cast<LO> (A->getLocalNumRows ());
    for (LO localRow = 0; localRow < numLocalRows; ++localRow) {
      A->getLocalRowView (localRow, ind, val);
      A_constGraph->replaceLocalValues (localRow, ind, val);
    }
    A_constGraph->fillComplete (A->getDomainMap (), A->getRangeMap ());
    A = A_constGraph; // Replace A with A_constGraph.
  }

  Teuchos::RCP<const map_type> rowmap = A->getRowMap();

  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(rowmap, 1));

  if (b == Teuchos::null) {
    bool rhs_unit = false;
    int rhs_option = 2;
    b = Teuchos::rcp (new TMV (rowmap, 1));
    if (rhs_option == 0) {
      // random B
      b->randomize ();
    } else if (rhs_option == 1) {
      // b = ones
      b->putScalar (STS::one ());
    } else {
      if (rhs_option == 2) {
        // b = A * random
        x->randomize ();
      } else {
        // b = A * ones
        x->putScalar (STS::one ());
      }
      BOPT::Apply (*A, *x, *b);
    }
    if (rhs_unit) {
      // scale B to unit-norm
      Teuchos::Array<typename STS::magnitudeType> normsB(b->getNumVectors());
      b->norm2(normsB);
      for (size_t j = 0; j < b->getNumVectors(); j++) {
        b->getVectorNonConst(j)->scale(STS::one() / normsB[j]);
      }
    }
    // X = zero
    BMVT::MvInit (*x, STS::zero ());
  }
  else {
    x->putScalar (STS::zero ());
  }

  if (useZoltan2) {
#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
    // Create an input adapter for the Tpetra matrix.
    Zoltan2::XpetraCrsMatrixAdapter<crs_matrix_type>
      zoltan_matrix(A);

    // Specify partitioning parameters
    Teuchos::ParameterList zoltan_params;
    zoltan_params.set("partitioning_approach", "partition");
    //
    if (useParMETIS) {
      if (comm->getRank() == 0) {
        std::cout << "Using Zoltan2(ParMETIS)" << std::endl;
      }
      zoltan_params.set("algorithm", "parmetis");
      zoltan_params.set("symmetrize_input", "transpose");
      zoltan_params.set("partitioning_objective", "minimize_cut_edge_weight");
    } else {
      if (comm->getRank() == 0) {
        std::cout << "Using Zoltan2(HyperGraph)" << std::endl;
      }
      zoltan_params.set("algorithm", "phg");
    }

    // Create and solve partitioning problem
    Zoltan2::PartitioningProblem<Zoltan2::XpetraCrsMatrixAdapter<crs_matrix_type>>
      problem(&zoltan_matrix, &zoltan_params);
    problem.solve();

    // Redistribute matrix
    RCP<crs_matrix_type> zoltan_A;
    zoltan_matrix.applyPartitioningSolution (*A, zoltan_A, problem.getSolution());
    // Set it as coefficient matrix
    A = zoltan_A;

    // Redistribute RHS
    RCP<TMV> zoltan_b;
    Zoltan2::XpetraMultiVectorAdapter<TMV> adapterRHS(rcpFromRef (*b));
    adapterRHS.applyPartitioningSolution (*b, zoltan_b, problem.getSolution());
    // Set it as RHS
    b = zoltan_b;

    // Redistribute Sol
    RCP<TMV> zoltan_x;
    Zoltan2::XpetraMultiVectorAdapter<TMV> adapterSol(rcpFromRef (*x));
    adapterSol.applyPartitioningSolution (*x, zoltan_x, problem.getSolution());
    // Set it as Sol
    x = zoltan_x;
#else
    TEUCHOS_TEST_FOR_EXCEPTION(
      useZoltan2, std::invalid_argument,
      "Both Xpetra and Zoltan2 are needed to use Zoltan2.");
#endif
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

  return problem;
}


#endif

