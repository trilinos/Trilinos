// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHelpers.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

#include "Stokhos_Sacado_Kokkos.hpp"

// Tpetra
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Kokkos device types
#include "Kokkos_Serial.hpp"

#ifdef HAVE_STOKHOS_BELOS
// Belos solver
#include "Belos_TpetraAdapter_MP_Vector.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#endif

template <typename scalar, typename ordinal>
inline
scalar generate_vector_coefficient( const ordinal nFEM,
                                    const ordinal nStoch,
                                    const ordinal iColFEM,
                                    const ordinal iStoch )
{
  const scalar X_fem = 100.0 + scalar(iColFEM) / scalar(nFEM);
  const scalar X_stoch =  1.0 + scalar(iStoch) / scalar(nStoch);
  return X_fem + X_stoch;
  //return 1.0;
}

// Currently only using Kokkos::Serial device for every node
template <typename Node>
struct DeviceForNode {
  typedef Kokkos::Serial type;
};

//
// Tests
//

//
// Test matrix-vector multiplication for a simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, Basic, BaseScalar, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  const LocalOrdinal VectorSize = 3;
  typedef typename DeviceForNode<Node>::type Device;
  typedef Stokhos::StaticFixedStorage<LocalOrdinal,BaseScalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if (!Device::is_initialized())
    Device::initialize();

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  RCP<Node> node = rcp(new Node);
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal>(
      nrow, comm, node);
  RCP<Tpetra_CrsGraph> graph = Tpetra::createCrsGraph(map, size_t(2));
  Array<GlobalOrdinal> columnIndices(2);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
  const size_t num_my_row = myGIDs.size();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    size_t ncol = 1;
    if (row != nrow-1) {
      columnIndices[1] = row+1;
      ncol = 2;
    }
    graph->insertGlobalIndices(row, columnIndices(0,ncol));
  }
  graph->fillComplete();
  RCP<Tpetra_CrsMatrix> matrix = rcp(new Tpetra_CrsMatrix(graph));

  // Set values in matrix
  Array<Scalar> vals(2);
  Scalar val(VectorSize, BaseScalar(0.0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    for (LocalOrdinal j=0; j<VectorSize; ++j)
      val.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, VectorSize, row, j);
    vals[0] = val;
    size_t ncol = 1;

    if (row != nrow-1) {
      columnIndices[1] = row+1;
      for (LocalOrdinal j=0; j<VectorSize; ++j)
        val.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(
          nrow, VectorSize, row+1, j);
      vals[1] = val;
      ncol = 2;
    }
    matrix->replaceGlobalValues(row, columnIndices(0,ncol), vals(0,ncol));
  }
  matrix->fillComplete();

  // Fill vector
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<VectorSize; ++j)
      val.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, VectorSize, row, j);
    x_view[i] = val;
  }

  // Multiply
  RCP<Tpetra_Vector> y = Tpetra::createVector<Scalar>(map);
  matrix->apply(*x, *y);

  // y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check
  ArrayRCP<Scalar> y_view = y->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      BaseScalar v = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, VectorSize, row, j);
      val.fastAccessCoeff(j) = v*v;
    }
    if (row != nrow-1) {
      for (LocalOrdinal j=0; j<VectorSize; ++j) {
        BaseScalar v = generate_vector_coefficient<BaseScalar,size_t>(
          nrow, VectorSize, row+1, j);
        val.fastAccessCoeff(j) += v*v;
      }
    }
    TEST_EQUALITY( y_view[i].size(), VectorSize );
    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_EQUALITY( y_view[i].fastAccessCoeff(j), val.fastAccessCoeff(j) );
  }
}

#ifdef HAVE_STOKHOS_BELOS

//
// Test Belos GMRES solve for a simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES, BaseScalar, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  const LocalOrdinal VectorSize = 3;
  typedef typename DeviceForNode<Node>::type Device;
  typedef Stokhos::StaticFixedStorage<LocalOrdinal,BaseScalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if (!Device::is_initialized())
    Device::initialize();

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  RCP<Node> node = rcp(new Node);
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal>(
      nrow, comm, node);
  RCP<Tpetra_CrsGraph> graph = Tpetra::createCrsGraph(map, size_t(2));
  Array<GlobalOrdinal> columnIndices(2);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
  const size_t num_my_row = myGIDs.size();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    size_t ncol = 1;
    if (row != nrow-1) {
      columnIndices[1] = row+1;
      ncol = 2;
    }
    graph->insertGlobalIndices(row, columnIndices(0,ncol));
  }
  graph->fillComplete();
  RCP<Tpetra_CrsMatrix> matrix = rcp(new Tpetra_CrsMatrix(graph));

  // Set values in matrix
  Array<Scalar> vals(2);
  Scalar val(VectorSize, BaseScalar(0.0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    for (LocalOrdinal j=0; j<VectorSize; ++j)
      val.fastAccessCoeff(j) = j+1;
    vals[0] = val;
    size_t ncol = 1;

    if (row != nrow-1) {
      columnIndices[1] = row+1;
      for (LocalOrdinal j=0; j<VectorSize; ++j)
        val.fastAccessCoeff(j) = j+1;
      vals[1] = val;
      ncol = 2;
    }
    matrix->replaceGlobalValues(row, columnIndices(0,ncol), vals(0,ncol));
  }
  matrix->fillComplete();

  // Fill RHS vector
  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> b_view = b->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    b_view[i] = Scalar(1.0);
  }

  // Solve
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
  typedef BaseScalar BelosScalar;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
  typedef Belos::OperatorTraits<BelosScalar,MV,OP> BOPT;
  typedef Belos::MultiVecTraits<BelosScalar,MV> BMVT;
  typedef Belos::LinearProblem<BelosScalar,MV,OP> BLinProb;
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  RCP< BLinProb > problem = rcp(new BLinProb(matrix, x, b));
  RCP<ParameterList> belosParams = rcp(new ParameterList);
  typename ST::magnitudeType tol = 1e-12;
  belosParams->set("Flexible Gmres", false);
  belosParams->set("Num Blocks", 100);
  belosParams->set("Convergence Tolerance", BelosScalar(tol));
  belosParams->set("Maximum Iterations", 100);
  belosParams->set("Verbosity", 33);
  belosParams->set("Output Style", 1);
  belosParams->set("Output Frequency", 1);
  belosParams->set("Output Stream", out.getOStream());
  RCP<Belos::SolverManager<BelosScalar,MV,OP> > solver =
    rcp(new Belos::PseudoBlockGmresSolMgr<BelosScalar,MV,OP>(problem, belosParams));
  problem->setProblem();
  Belos::ReturnType ret = solver->solve();
  TEST_EQUALITY_CONST( ret, Belos::Converged );

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check -- Correct answer is:
  //     [ 0, 0,   ..., 0            ]
  //     [ 1, 1/2, ..., 1/VectorSize ]
  //     [ 0, 0,   ..., 0            ]
  //     [ 1, 1/2, ..., 1/VectorSize ]
  //     ....
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row % 2) {
      for (LocalOrdinal j=0; j<VectorSize; ++j) {
        val.fastAccessCoeff(j) = BaseScalar(1.0) / BaseScalar(j+1);
      }
    }
    else
      val = Scalar(0.0);
    TEST_EQUALITY( x_view[i].size(), VectorSize );

    // Set small values to zero
    Scalar v = x_view[i];
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      if (ST::magnitude(v.coeff(j)) < tol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), val.coeff(j), tol);
  }
}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES, BaseScalar, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

#define CRSMATRIX_MP_VECTOR_TESTS_SLGN(S, LO, GO, N)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Basic, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, BelosGMRES, S, LO, GO, N )
