// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
 * \file   Superlu_Solve_Tests.cpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Thu May 19 08:56:28 2011
 *
 * \brief This is meant to test all the SuperLU solver against a
 * battery of different matrices.  We try to test with matrices of
 * it's supported ordinal and scalar types.
 */

#include <string>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <MatrixMarket_Tpetra.hpp> // for loading matrices from file

#include "Amesos2.hpp"

namespace {

  using std::string;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ScalarTraits;
  using Teuchos::OrdinalTraits;
  using Teuchos::ETransp;
  using Teuchos::CONJ_TRANS;
  using Teuchos::TRANS;
  using Teuchos::NO_TRANS;


  using Tpetra::global_size_t;
  using Tpetra::CrsMatrix;
  using Tpetra::MultiVector;
  using Tpetra::Map;

  bool testMpi = false;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-mpi", "test-serial", &testMpi,
                  "Test Serial by default (for now) or force MPI test.  In a serial build,"
                  " this option is ignored and a serial comm is always used." );
  }

  // The Process:
  //
  // 1. Create random multi-vector X
  // 2. Perform AMat * X to get rhs multi-vector B
  // 3. Use Superlu to solve AMat * Xhat = B
  // 4. Check that Xhat ~= X
  //
  // Small note on the final comparison: The Teuchos utilities to
  // compare two floating arrays don't work that well when one of the
  // arrays contains only zeros, since the relative error will always
  // turn out to be somewhere right around 1, no matter what the
  // number in the other array is.  This is the way that Tpetra
  // compares the output from its solve routine in its unit tests, but
  // I believe this method could return a false-positive.  What I have
  // done here is to instead compare the norm2 of X and Xhat directly,
  // with a reasonably small tolerance.  While this could also produce
  // a false-positive, it is less likely.
#define TEST_WITH_MATRIX(MATNAME, transpose)                            \
  typedef CrsMatrix<SCALAR,LO,GO> MAT;                                  \
  typedef MultiVector<SCALAR,LO,GO> MV;                                 \
  typedef ScalarTraits<SCALAR> ST;                                      \
  typedef typename ST::magnitudeType Mag;                               \
  typedef ScalarTraits<Mag> MT;                                         \
  const size_t numVecs = 5;                                             \
  ETransp trans = ((transpose) ? CONJ_TRANS : NO_TRANS);                \
                                                                        \
  RCP<const Comm<int> > comm = Tpetra.getDefaultComm();                 \
                                                                        \
  string path = string("../matrices/") + (MATNAME);                     \
  RCP<MAT> AMat =                                                       \
    Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(path,comm);       \
                                                                        \
  RCP<const Map<LO,GO> > dmnmap = AMat->getDomainMap();                 \
  RCP<const Map<LO,GO> > rngmap = AMat->getRangeMap();                  \
                                                                        \
  RCP<MV> X, B, Xhat;                                                   \
  if( transpose ){                                                      \
    X = rcp(new MV(dmnmap,numVecs));                                    \
    B = rcp(new MV(rngmap,numVecs));                                    \
    Xhat = rcp(new MV(dmnmap,numVecs));                                 \
  } else {                                                              \
    X = rcp(new MV(rngmap,numVecs));                                    \
    B = rcp(new MV(dmnmap,numVecs));                                    \
    Xhat = rcp(new MV(rngmap,numVecs));                                 \
  }                                                                     \
  X->setObjectLabel("X");                                               \
  B->setObjectLabel("B");                                               \
  Xhat->setObjectLabel("Xhat");                                         \
  X->randomize();                                                       \
  AMat->apply(*X,*B,trans);                                             \
                                                                        \
  RCP<Amesos2::Solver<MAT,MV> > solver                                  \
  = Amesos2::create<MAT,MV>("Superlu", AMat, Xhat, B );                 \
                                                                        \
  Teuchos::ParameterList amesos2_params("Amesos2");                     \
  if( transpose ){                                                      \
    amesos2_params.sublist("SuperLU").set("Trans","CONJ","Solve with transpose"); \
  } else {                                                              \
    amesos2_params.sublist("SuperLU"). set("Trans","NOTRANS","Do not solve with transpose"); \
  }                                                                     \
                                                                        \
  solver->setParameters( rcpFromRef(amesos2_params) );                  \
  solver->symbolicFactorization().numericFactorization().solve();       \
                                                                        \
  Array<Mag> xhatnorms(numVecs), xnorms(numVecs);                       \
  Xhat->norm2(xhatnorms());                                             \
  X->norm2(xnorms());                                                   \
  TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 )

  /**************
   * UNIT TESTS *
   **************/

#define SUPERLU_MATRIX_TEST_DECL(MATNAME)                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Superlu, MATNAME, LO, GO, SCALAR)   \
  {                                                                     \
    string matfile = #MATNAME + string(".mtx");                         \
    TEST_WITH_MATRIX(matfile, false);                                   \
    Xhat->describe(out, Teuchos::VERB_EXTREME);                         \
    X->describe(out, Teuchos::VERB_EXTREME);                            \
  }                                                                     \
                                                                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Superlu, MATNAME##_trans, LO, GO, SCALAR) \
  {                                                                     \
    string matfile = #MATNAME + string(".mtx");                         \
    TEST_WITH_MATRIX(matfile, true);                                    \
    Xhat->describe(out, Teuchos::VERB_EXTREME);                         \
    X->describe(out, Teuchos::VERB_EXTREME);                            \
  }

  /*************************************
   * Declarations for all the matrices *
   *************************************/

  SUPERLU_MATRIX_TEST_DECL(arc130)
  SUPERLU_MATRIX_TEST_DECL(bcsstk01)
  SUPERLU_MATRIX_TEST_DECL(bcsstk18)
  SUPERLU_MATRIX_TEST_DECL(bcsstm01)
  SUPERLU_MATRIX_TEST_DECL(beacxc)
  SUPERLU_MATRIX_TEST_DECL(gemat12)
  SUPERLU_MATRIX_TEST_DECL(sherman3)
  SUPERLU_MATRIX_TEST_DECL(young1c)


#define SUPERLU_MATRIX_TEST(MATNAME, L, G, S)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(Superlu, MATNAME, L, G, S)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(Superlu, MATNAME##_trans, L, G, S)


#define ARC130_SOLVE(LO, GO, SCALAR)            \
  SUPERLU_MATRIX_TEST(arc130, LO, GO, SCALAR)

  // MatrixMarket read error? ::
  //
  // Cannot add entry A(37,18) = -2.08333e+06 to matrix; already have expected number of entries 224.
#define BCSSTK01_SOLVE(LO, GO, SCALAR)                  \
  // SUPERLU_MATRIX_TEST(bcsstk01, LO, GO, SCALAR)

  // MatrixMarket read error? ::
  //
  // Cannot add entry A(6510,7637) = 328839 to matrix; already have expected number of entries 80519.
#define BCSSTK18_SOLVE(LO, GO, SCALAR)                  \
  // SUPERLU_MATRIX_TEST(bcsstk18, LO, GO, SCALAR)

  // Integer matrices not yet supported
#define BCSSTM01_SOLVE(LO, GO, SCALAR)                  \
  // SUPERLU_MATRIX_TEST(bcsstm01, LO, GO, SCALAR)

  // This is a rectangular matrix
  //
  //   Throw test that evaluated to true: *A.getMap() != *importer.getSourceMap()
  //   Source Maps don't match.
#define BEACXC_SOLVE(LO, GO, SCALAR)                    \
  // SUPERLU_MATRIX_TEST(beacxc, LO, GO, SCALAR)

#define GEMAT12_SOLVE(LO, GO, SCALAR)           \
  SUPERLU_MATRIX_TEST(gemat12, LO, GO, SCALAR)

#define SHERMAN3_SOLVE(LO, GO, SCALAR)          \
  SUPERLU_MATRIX_TEST(sherman3, LO, GO, SCALAR)

  // A complex valued matrix
  //
  // Currently, the transpose solve for complex problems is not
  // functioning, so we do just the standard solve
#define YOUNG1C_SOLVE(LO, GO, COMPLEX)                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(Superlu, young1c, LO, GO, COMPLEX)



  /*****************************
   * Instantiations with types *
   *****************************/

#define UNIT_TEST_GROUP_ORDINALS_SCALAR(LO, GO, SCALAR) \
  ARC130_SOLVE(LO, GO, SCALAR)                          \
  BCSSTK01_SOLVE(LO, GO, SCALAR)                        \
  BCSSTK18_SOLVE(LO, GO, SCALAR)                        \
  BCSSTM01_SOLVE(LO, GO, SCALAR)                        \
  BEACXC_SOLVE(LO, GO, SCALAR)                          \
  GEMAT12_SOLVE(LO, GO, SCALAR)                         \
  SHERMAN3_SOLVE(LO, GO, SCALAR)

#ifdef HAVE_TPETRA_INST_FLOAT
#  define UNIT_TEST_GROUP_ORDINALS_FLOAT(LO, GO)        \
  UNIT_TEST_GROUP_ORDINALS_SCALAR(LO, GO, float)
#else
#  define UNIT_TEST_GROUP_ORDINALS_FLOAT(LO, GO)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
#  define UNIT_TEST_GROUP_ORDINALS_DOUBLE(LO, GO)       \
  UNIT_TEST_GROUP_ORDINALS_SCALAR(LO, GO, double)
#else
#  define UNIT_TEST_GROUP_ORDINALS_DOUBLE(LO, GO)
#endif

#define UNIT_TEST_GROUP_ORDINALS_REALS(LO, GO)          \
  UNIT_TEST_GROUP_ORDINALS_FLOAT(LO, GO)                \
  UNIT_TEST_GROUP_ORDINALS_DOUBLE(LO, GO)

#ifdef HAVE_TEUCHOS_COMPLEX
#  ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
#    define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO) \
  typedef std::complex<float> ComplexFloat;               \
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexFloat)
#  else
#    define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  endif
#  ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
#    define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)      \
  typedef std::complex<double> ComplexDouble;                   \
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexDouble)
#  else
#    define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#  endif
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX(LO, GO)       \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)         \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#  define UNIT_TEST_GROUP_ORDINALS_COMPLEX(LO, GO)
#endif

#define UNIT_TEST_GROUP_ORDINALS(LO, GO)        \
  UNIT_TEST_GROUP_ORDINALS_REALS(LO, GO)        \
  UNIT_TEST_GROUP_ORDINALS_COMPLEX(LO, GO)

  // Entry:
  UNIT_TEST_GROUP_ORDINALS(int, int)

#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
  typedef long int LongInt;
  UNIT_TEST_GROUP_ORDINALS(int, LongInt)
#  ifdef HAVE_TPETRA_INT_LONG_LONG
  typedef long long int LongLongInt;
  UNIT_TEST_GROUP_ORDINALS(int, LongLongInt)
#  endif
#endif  // EXPL-INST

} // end anonymous namespace
