// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2010 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

/**
 * @file   Superlu_Solve_Tests.cpp
 * @author Eric Bavier <etbavie@sandia.gov>
 * @date   Thu May 19 08:56:28 2011
 *
 * @brief This is meant to test all the SuperLU solver against a
 * battery of different matrices.  We try to test with matrices of
 * it's supported ordinal and scalar types.
 */

#include <string>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp> // for loading matrices from file

#include "Amesos2_Factory.hpp"
#include "Amesos2_Util_is_same.hpp"
#include "Amesos2_Superlu.hpp"

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

  using Amesos::Util::is_same;

  typedef Tpetra::DefaultPlatform::DefaultPlatformType           Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

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
  // 4. Check that Xhat == X
#define TEST_WITH_MATRIX(MATNAME, transpose, TOL)			\
  typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;				\
  typedef MultiVector<SCALAR,LO,GO,Node> MV;				\
  const size_t numVecs = 5;						\
  ETransp trans = ((transpose) ? CONJ_TRANS : NO_TRANS);		\
									\
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();	\
  RCP<const Comm<int> > comm = platform.getComm();			\
  RCP<Node>             node = platform.getNode();			\
									\
  string path = string("../matrices/") + (MATNAME);			\
  RCP<MAT> AMat =							\
    Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(path,comm,node);	\
									\
  RCP<const Map<LO,GO,Node> > rowmap = AMat->getRowMap();		\
  RCP<const Map<LO,GO,Node> > colmap = AMat->getColMap();		\
									\
  RCP<MV> X = rcp(new MV(rowmap,numVecs));				\
  RCP<MV> B = rcp(new MV(colmap,numVecs));				\
  RCP<MV> Xhat = rcp(new MV(rowmap,numVecs));				\
  X->setObjectLabel("X");						\
  B->setObjectLabel("B");						\
  Xhat->setObjectLabel("Xhat");						\
  X->randomize();							\
  AMat->apply(*X,*B,trans);						\
									\
  RCP<Amesos::SolverBase> solver					\
  = Amesos::Factory<MAT,MV>::create("Superlu", AMat, Xhat, B );		\
									\
  Teuchos::ParameterList params;					\
  if( transpose ){							\
    params.set("Trans","CONJ","Solve with transpose");			\
  } else {								\
    params.set("Trans","NOTRANS","Do not solve with transpose");	\
  }									\
									\
  solver->setParameters( rcpFromRef(params) );				\
  solver->symbolicFactorization().numericFactorization().solve();	\
									\
  TEST_COMPARE_FLOATING_ARRAYS(Xhat->get1dView(), X->get1dView(), TOL)


  /**************
   * UNIT TESTS *
   **************/

#define SUPERLU_MATRIX_TEST_WITH_TOL_DECL(MATNAME, TI, TF)		\
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Superlu, MATNAME##_##TI##TF, LO, GO, SCALAR) \
  {									\
    string matfile = #MATNAME + string(".mtx");				\
    double tol = TI##.##TF;						\
    TEST_WITH_MATRIX(matfile, false, tol);				\
  }

#define SUPERLU_MATRIX_TEST(MATNAME, L, G, S, TI, TF)			\
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(Superlu, MATNAME##_##TI##TF, L, G, S)

  /*************************************
   * Declarations for all the matrices *
   *************************************/

#define ALL_MATS_WITH_TOL(TI, TF)			\
  SUPERLU_MATRIX_TEST_WITH_TOL_DECL(arc130, TI, TF)	\
  SUPERLU_MATRIX_TEST_WITH_TOL_DECL(bcsstk01, TI, TF)	\
  SUPERLU_MATRIX_TEST_WITH_TOL_DECL(bcsstk18, TI, TF)	\
  SUPERLU_MATRIX_TEST_WITH_TOL_DECL(bcsstm01, TI, TF)	\
  SUPERLU_MATRIX_TEST_WITH_TOL_DECL(beacxc, TI, TF)	\
  SUPERLU_MATRIX_TEST_WITH_TOL_DECL(gemat12, TI, TF)	\
  SUPERLU_MATRIX_TEST_WITH_TOL_DECL(sherman3, TI, TF)	\
  SUPERLU_MATRIX_TEST_WITH_TOL_DECL(young1c, TI, TF)

  /*** Predeclarations of all tolerance classes ***/
  ALL_MATS_WITH_TOL(5, 0)	// i.e. tol := 5.0
  ALL_MATS_WITH_TOL(0, 2)	// i.e. tol := 0.2
  //  ALL_MATS_WITH_TOL(1, 0)	// i.e. tol := 1.0

#define ARC130_SOLVE(LO, GO, SCALAR, TI, TF)		\
  SUPERLU_MATRIX_TEST(arc130, LO, GO, SCALAR, TI, TF)

  // MatrixMarket read error? ::
  //
  // Cannot add entry A(37,18) = -2.08333e+06 to matrix; already have expected number of entries 224.
#define BCSSTK01_SOLVE(LO, GO, SCALAR, TI, TF)		\
  // SUPERLU_MATRIX_TEST(bcsstk01, LO, GO, SCALAR, TI, TF)

  // MatrixMarket read error? ::
  //
  // Cannot add entry A(6510,7637) = 328839 to matrix; already have expected number of entries 80519.
#define BCSSTK18_SOLVE(LO, GO, SCALAR, TI, TF)		\
  // SUPERLU_MATRIX_TEST(bcsstk18, LO, GO, SCALAR, TI, TF)

  // Integer matrices not yet supported
#define BCSSTM01_SOLVE(LO, GO, SCALAR, TI, TF)		\
  // SUPERLU_MATRIX_TEST(bcsstm01, LO, GO, SCALAR, TI, TF)

  // This is a rectangular matrix
  //
  //   Throw test that evaluated to true: *A.getMap() != *importer.getSourceMap()
  //   Source Maps don't match.
#define BEACXC_SOLVE(LO, GO, SCALAR, TI, TF)		\
  SUPERLU_MATRIX_TEST(beacxc, LO, GO, SCALAR, TI, TF)

#define GEMAT12_SOLVE(LO, GO, SCALAR, TI, TF)		\
  SUPERLU_MATRIX_TEST(gemat12, LO, GO, SCALAR, TI, TF)

#define SHERMAN3_SOLVE(LO, GO, SCALAR, TI, TF)		\
  SUPERLU_MATRIX_TEST(sherman3, LO, GO, SCALAR, TI, TF)

  // A complex valued matrix
#define YOUNG1C_SOLVE(LO, GO, COMPLEX, TI, TF)		\
  SUPERLU_MATRIX_TEST(young1c, LO, GO, COMPLEX, TI, TF)


  /*****************************
   * Instantiations with types *
   *****************************/

		    // Note: The tolerances used here much be predeclared above

#define UNIT_TEST_GROUP_ORDINALS_SCALAR(LO, GO, SCALAR, TI, TF)	\
  ARC130_SOLVE(LO, GO, SCALAR, TI, TF)				\
  BCSSTK01_SOLVE(LO, GO, SCALAR, TI, TF)			\
  BCSSTK18_SOLVE(LO, GO, SCALAR, TI, TF)			\
  BCSSTM01_SOLVE(LO, GO, SCALAR, TI, TF)			\
  BEACXC_SOLVE(LO, GO, SCALAR, TI, TF)				\
  GEMAT12_SOLVE(LO, GO, SCALAR, TI, TF)				\
  SHERMAN3_SOLVE(LO, GO, SCALAR, TI, TF)			

#define UNIT_TEST_GROUP_ORDINALS_REALS(LO, GO)          \
  UNIT_TEST_GROUP_ORDINALS_SCALAR(LO, GO, float, 5, 0)	\
  UNIT_TEST_GROUP_ORDINALS_SCALAR(LO, GO, double, 0, 2)

#define UNIT_TEST_GROUP_ORDINALS_COMPLEX(LO, GO)        \
  typedef std::complex<float>  ComplexFloat;            \
  YOUNG1C_SOLVE(LO, GO, ComplexFloat, 5, 0)		\
  typedef std::complex<double> ComplexDouble;           \
  YOUNG1C_SOLVE(LO, GO, ComplexDouble, 0, 2)

#define UNIT_TEST_GROUP_ORDINALS(LO, GO)        \
  UNIT_TEST_GROUP_ORDINALS_REALS(LO, GO)        \
  UNIT_TEST_GROUP_ORDINALS_COMPLEX(LO, GO)

  // Entry:
  UNIT_TEST_GROUP_ORDINALS(int, int)

  typedef long int LongInt;
  UNIT_TEST_GROUP_ORDINALS(int, LongInt)
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
  typedef long long int LongLongInt;
  UNIT_TEST_GROUP_ORDINALS(int, LongLongInt)
#    endif

} // end anonymous namespace
