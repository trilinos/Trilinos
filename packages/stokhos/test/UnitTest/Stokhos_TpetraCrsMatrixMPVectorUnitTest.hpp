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
#include "Stokhos_Tpetra_Utilities_MP_Vector.hpp"

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

// Belos solver
#ifdef HAVE_STOKHOS_BELOS
#include "Belos_TpetraAdapter_MP_Vector.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#endif

// Ifpack2 preconditioner
#ifdef HAVE_STOKHOS_IFPACK2
#include "Stokhos_Ifpack2_MP_Vector.hpp"
#endif

// MueLu preconditioner
#ifdef HAVE_STOKHOS_MUELU
#include "MueLu_CreateTpetraPreconditioner.hpp"
#endif

// Amesos2 solver
#ifdef HAVE_STOKHOS_AMESOS2
#include "Stokhos_Amesos2_MP_Vector.hpp"
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

//
// Test flattening MP::Vector matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, Flatten, BaseScalar, LocalOrdinal, GlobalOrdinal, Node )
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

  typedef Tpetra::Vector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_Vector;
  typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_CrsMatrix;

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

  // Flatten matrix
  RCP<const Tpetra_Map> flat_x_map, flat_y_map;
  RCP<const Tpetra_CrsGraph> flat_graph =
    Stokhos::create_flat_mp_graph(*graph, flat_x_map, flat_y_map, VectorSize);
  RCP<Flat_Tpetra_CrsMatrix> flat_matrix =
    Stokhos::create_flat_matrix(*matrix, flat_graph);

  // Multiply with flattened matix
  RCP<Tpetra_Vector> y2 = Tpetra::createVector<Scalar>(map);
  RCP<Flat_Tpetra_Vector> flat_x =
    Stokhos::create_flat_vector_view(*x, flat_x_map);
  RCP<Flat_Tpetra_Vector> flat_y =
    Stokhos::create_flat_vector_view(*y2, flat_y_map);
  flat_matrix->apply(*flat_x, *flat_y);

  // flat_y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //                  Teuchos::VERB_EXTREME);

  // Check
  ArrayRCP<Scalar> y_view = y->get1dViewNonConst();
  ArrayRCP<Scalar> y2_view = y2->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    TEST_EQUALITY( y_view[i].size(), VectorSize );
    TEST_EQUALITY( y2_view[i].size(), VectorSize );
    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_EQUALITY( y_view[i].fastAccessCoeff(j),
                     y2_view[i].fastAccessCoeff(j) );
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

#if defined(HAVE_STOKHOS_BELOS) && defined(HAVE_STOKHOS_IFPACK2)

//
// Test Belos GMRES solve with Ifpack2 RILUK preconditioning for a
// simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES_RILUK, BaseScalar, LocalOrdinal, GlobalOrdinal, Node )
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

  // Create preconditioner
  typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Prec;
  Ifpack2::Factory factory;
  RCP<Prec> M = factory.create<Tpetra_CrsMatrix>("RILUK", matrix);
  M->initialize();
  M->compute();

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
  problem->setRightPrec(M);
  problem->setProblem();
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
  //belosParams->set("Orthogonalization", "TSQR");
  RCP<Belos::SolverManager<BelosScalar,MV,OP> > solver =
    rcp(new Belos::PseudoBlockGmresSolMgr<BelosScalar,MV,OP>(problem, belosParams));
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
  Tpetra_CrsMatrix_MP, BelosGMRES_RILUK, BaseScalar, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

#if defined(HAVE_STOKHOS_BELOS) && defined(HAVE_STOKHOS_IFPACK2) && defined(HAVE_STOKHOS_MUELU)

//
// Test Belos CG solve with MueLu preconditioning for a 1-D Laplacian matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosCG_Muelu, BaseScalar, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::getParametersFromXmlFile;

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

  // 1-D Laplacian matrix
  GlobalOrdinal nrow = 50;
  BaseScalar h = 1.0 / static_cast<BaseScalar>(nrow-1);
  RCP<const Tpetra_Comm> comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  RCP<Node> node = rcp(new Node);
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal>(
      nrow, comm, node);
  RCP<Tpetra_CrsGraph> graph = Tpetra::createCrsGraph(map, size_t(3));
  Array<GlobalOrdinal> columnIndices(3);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
  const size_t num_my_row = myGIDs.size();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row == 0 || row == nrow-1) { // Boundary nodes
      columnIndices[0] = row;
      graph->insertGlobalIndices(row, columnIndices(0,1));
    }
    else { // Interior nodes
      columnIndices[0] = row-1;
      columnIndices[1] = row;
      columnIndices[2] = row+1;
      graph->insertGlobalIndices(row, columnIndices(0,3));
    }
  }
  graph->fillComplete();
  RCP<Tpetra_CrsMatrix> matrix = rcp(new Tpetra_CrsMatrix(graph));

  // Set values in matrix
  Array<Scalar> vals(3);
  Scalar a_val(VectorSize, BaseScalar(0.0));
  for (LocalOrdinal j=0; j<VectorSize; ++j) {
    a_val.fastAccessCoeff(j) =
      BaseScalar(1.0) + BaseScalar(j) / BaseScalar(VectorSize);
  }
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row == 0 || row == nrow-1) { // Boundary nodes
      columnIndices[0] = row;
      vals[0] = Scalar(1.0);
      matrix->replaceGlobalValues(row, columnIndices(0,1), vals(0,1));
    }
    else {
      columnIndices[0] = row-1;
      columnIndices[1] = row;
      columnIndices[2] = row+1;
      vals[0] = Scalar(1.0) * a_val;
      vals[1] = Scalar(-2.0) * a_val;
      vals[2] = Scalar(1.0) * a_val;
      matrix->replaceGlobalValues(row, columnIndices(0,3), vals(0,3));
    }
  }
  matrix->fillComplete();

  // Fill RHS vector
  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> b_view = b->get1dViewNonConst();
  Scalar b_val(VectorSize, BaseScalar(0.0));
  for (LocalOrdinal j=0; j<VectorSize; ++j) {
    b_val.fastAccessCoeff(j) =
      BaseScalar(-1.0) + BaseScalar(j) / BaseScalar(VectorSize);
  }
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row == 0 || row == nrow-1)
      b_view[i] = Scalar(0.0);
    else
      b_view[i] = Scalar(b_val * h * h);
  }

  // Create preconditioner
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
  RCP<ParameterList> muelu_params =
    getParametersFromXmlFile("muelu_cheby.xml");
  RCP<OP> M =
    MueLu::CreateTpetraPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>(matrix);
  // typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Prec;
  // Ifpack2::Factory factory;
  // RCP<Prec> M = factory.create<Tpetra_CrsMatrix>("RILUK", matrix);
  // M->initialize();
  // M->compute();

  // Solve
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
  typedef BaseScalar BelosScalar;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Belos::OperatorTraits<BelosScalar,MV,OP> BOPT;
  typedef Belos::MultiVecTraits<BelosScalar,MV> BMVT;
  typedef Belos::LinearProblem<BelosScalar,MV,OP> BLinProb;
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  RCP< BLinProb > problem = rcp(new BLinProb(matrix, x, b));
  problem->setRightPrec(M);
  problem->setProblem();
  RCP<ParameterList> belosParams = rcp(new ParameterList);
  typename ST::magnitudeType tol = 1e-12;
  belosParams->set("Num Blocks", 100);
  belosParams->set("Convergence Tolerance", BelosScalar(tol));
  belosParams->set("Maximum Iterations", 100);
  belosParams->set("Verbosity", 33);
  belosParams->set("Output Style", 1);
  belosParams->set("Output Frequency", 1);
  belosParams->set("Output Stream", out.getOStream());
  //belosParams->set("Orthogonalization", "TSQR");
  RCP<Belos::SolverManager<BelosScalar,MV,OP> > solver =
    rcp(new Belos::PseudoBlockGmresSolMgr<BelosScalar,MV,OP>(problem, belosParams));
  Belos::ReturnType ret = solver->solve();
  TEST_EQUALITY_CONST( ret, Belos::Converged );

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check -- For a*y'' = b, correct answer is y = 0.5 *(b/a) * x * (x-1)
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  Scalar val(VectorSize, BaseScalar(0.0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    BaseScalar xx = row * h;
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      val.fastAccessCoeff(j) =
        BaseScalar(0.5) * (b_val.coeff(j)/a_val.coeff(j)) * xx * (xx - BaseScalar(1.0));
    }
    TEST_EQUALITY( x_view[i].size(), VectorSize );

    // Set small values to zero
    Scalar v = x_view[i];
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      if (ST::magnitude(v.coeff(j)) < tol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
      if (ST::magnitude(val.coeff(j)) < tol)
        val.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), val.coeff(j), tol);
  }

}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosCG_Muelu, BaseScalar, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

#if defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_AMESOS2_SUPERLU)

//
// Test Amesos2 solve for a 1-D Laplacian matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, Amesos2, BaseScalar, LocalOrdinal, GlobalOrdinal, Node )
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
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if (!Device::is_initialized())
    Device::initialize();

  // 1-D Laplacian matrix
  GlobalOrdinal nrow = 50;
  BaseScalar h = 1.0 / static_cast<BaseScalar>(nrow-1);
  RCP<const Tpetra_Comm> comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  RCP<Node> node = rcp(new Node);
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal>(
      nrow, comm, node);
  RCP<Tpetra_CrsGraph> graph = Tpetra::createCrsGraph(map, size_t(3));
  Array<GlobalOrdinal> columnIndices(3);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
  const size_t num_my_row = myGIDs.size();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row == 0 || row == nrow-1) { // Boundary nodes
      columnIndices[0] = row;
      graph->insertGlobalIndices(row, columnIndices(0,1));
    }
    else { // Interior nodes
      columnIndices[0] = row-1;
      columnIndices[1] = row;
      columnIndices[2] = row+1;
      graph->insertGlobalIndices(row, columnIndices(0,3));
    }
  }
  graph->fillComplete();
  RCP<Tpetra_CrsMatrix> matrix = rcp(new Tpetra_CrsMatrix(graph));

  // Set values in matrix
  Array<Scalar> vals(3);
  Scalar a_val(VectorSize, BaseScalar(0.0));
  for (LocalOrdinal j=0; j<VectorSize; ++j) {
    a_val.fastAccessCoeff(j) =
      BaseScalar(1.0) + BaseScalar(j) / BaseScalar(VectorSize);
  }
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row == 0 || row == nrow-1) { // Boundary nodes
      columnIndices[0] = row;
      vals[0] = Scalar(1.0);
      matrix->replaceGlobalValues(row, columnIndices(0,1), vals(0,1));
    }
    else {
      columnIndices[0] = row-1;
      columnIndices[1] = row;
      columnIndices[2] = row+1;
      vals[0] = Scalar(1.0) * a_val;
      vals[1] = Scalar(-2.0) * a_val;
      vals[2] = Scalar(1.0) * a_val;
      matrix->replaceGlobalValues(row, columnIndices(0,3), vals(0,3));
    }
  }
  matrix->fillComplete();

  // Fill RHS vector
  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> b_view = b->get1dViewNonConst();
  Scalar b_val(VectorSize, BaseScalar(0.0));
  for (LocalOrdinal j=0; j<VectorSize; ++j) {
    b_val.fastAccessCoeff(j) =
      BaseScalar(-1.0) + BaseScalar(j) / BaseScalar(VectorSize);
  }
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row == 0 || row == nrow-1)
      b_view[i] = Scalar(0.0);
    else
      b_view[i] = Scalar(b_val * h * h);
  }

  // Solve
  typedef Amesos2::Solver<Tpetra_CrsMatrix,Tpetra_MultiVector> Solver;
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  std::string solver_name;
#if defined(HAVE_AMESOS2_KLU2)
  solver_name = "klu2";
#elif defined(HAVE_AMESOS2_SUPERLUDIST)
  solver_name = "superlu_dist";
#elif defined(HAVE_AMESOS2_SUPERLUMT)
  solver_name = "superlu_mt";
#elif defined(HAVE_AMESOS2_SUPERLU)
  solver_name = "superlu";
#elif defined(HAVE_AMESOS2_PARDISO_MKL)
  solver_name = "pardisomkl";
#elif defined(HAVE_AMESOS2_LAPACK)
  solver_name = "lapack";
#elif defined(HAVE_AMESOS2_CHOLMOD) && defined (HAVE_AMESOS2_EXPERIMENTAL)
  solver_name = "lapack";
#else
  // if there are no solvers, we just return as a successful test
  success = true;
  return;
#endif
  RCP<Solver> solver = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>(
    solver_name, matrix, x, b);
  solver->solve();

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check -- For a*y'' = b, correct answer is y = 0.5 *(b/a) * x * (x-1)
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
  typename ST::magnitudeType tol = 1e-12;
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  Scalar val(VectorSize, BaseScalar(0.0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    BaseScalar xx = row * h;
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      val.fastAccessCoeff(j) =
        BaseScalar(0.5) * (b_val.coeff(j)/a_val.coeff(j)) * xx * (xx - BaseScalar(1.0));
    }
    TEST_EQUALITY( x_view[i].size(), VectorSize );

    // Set small values to zero
    Scalar v = x_view[i];
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      if (ST::magnitude(v.coeff(j)) < tol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
      if (ST::magnitude(val.coeff(j)) < tol)
        val.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), val.coeff(j), tol);
  }
}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, Amesos2, BaseScalar, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

#define CRSMATRIX_MP_VECTOR_TESTS_SLGN(S, LO, GO, N)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Basic, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Flatten, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, BelosGMRES, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, BelosGMRES_RILUK, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, BelosCG_Muelu, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Amesos2, S, LO, GO, N )
