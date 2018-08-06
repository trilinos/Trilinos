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

// Teuchos
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

// Tpetra
#include "Stokhos_Tpetra_UQ_PCE.hpp"
#include "Stokhos_Tpetra_Utilities.hpp"
#include "Stokhos_Tpetra_Utilities_UQ_PCE.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Stokhos_Tpetra_CG.hpp"

// Belos solver
#ifdef HAVE_STOKHOS_BELOS
#include "Belos_TpetraAdapter_UQ_PCE.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#endif

// Ifpack2 preconditioner
#ifdef HAVE_STOKHOS_IFPACK2
#include "Stokhos_Ifpack2_UQ_PCE.hpp"
#include "Ifpack2_Factory.hpp"
#endif

// MueLu preconditioner
#ifdef HAVE_STOKHOS_MUELU
#include "Stokhos_MueLu_UQ_PCE.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#endif

// Amesos2 solver
#ifdef HAVE_STOKHOS_AMESOS2
#include "Stokhos_Amesos2_UQ_PCE.hpp"
#include "Amesos2_Factory.hpp"
#endif

// Stokhos
#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

// Use "scalar" version of mean-based preconditioner (i.e., a preconditioner
// with double as the scalar type).  This is currently necessary to get the
// MueLu tests to pass on OpenMP and Cuda due to various kernels that don't
// work with the PCE scalar type.
#define USE_SCALAR_MEAN_BASED_PREC 1

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

template <typename scalar, typename ordinal>
inline
scalar generate_multi_vector_coefficient( const ordinal nFEM,
                                          const ordinal nVec,
                                          const ordinal nStoch,
                                          const ordinal iColFEM,
                                          const ordinal iVec,
                                          const ordinal iStoch)
{
  const scalar X_fem  = 100.0  + scalar(iColFEM) / scalar(nFEM);
  const scalar X_stoch =  1.0  + scalar(iStoch)  / scalar(nStoch);
  const scalar X_col    = 0.01 + scalar(iVec)    / scalar(nVec);
  return X_fem + X_stoch + X_col;
  //return 1.0;
}

template <typename scalar, typename ordinal>
inline
scalar generate_matrix_coefficient( const ordinal nFEM,
                                    const ordinal nStoch,
                                    const ordinal iRowFEM,
                                    const ordinal iColFEM,
                                    const ordinal iStoch )
{
  const scalar A_fem = ( 10.0 + scalar(iRowFEM) / scalar(nFEM) ) +
    (  5.0 + scalar(iColFEM) / scalar(nFEM) );

  const scalar A_stoch = ( 1.0 + scalar(iStoch) / scalar(nStoch) );

  return A_fem + A_stoch;
  //return 1.0;
}

template <typename kokkos_cijk_type, typename ordinal_type>
kokkos_cijk_type build_cijk(ordinal_type stoch_dim,
                            ordinal_type poly_ord)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;

  typedef typename kokkos_cijk_type::value_type value_type;
  typedef typename kokkos_cijk_type::execution_space execution_space;
  typedef Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> one_d_basis;
  typedef Stokhos::LegendreBasis<ordinal_type,value_type> legendre_basis;
  typedef Stokhos::CompletePolynomialBasis<ordinal_type,value_type> product_basis;
  typedef Stokhos::Sparse3Tensor<ordinal_type,value_type> Cijk;

  // Create product basis
  Array< RCP<const one_d_basis> > bases(stoch_dim);
  for (ordinal_type i=0; i<stoch_dim; i++)
    bases[i] = rcp(new legendre_basis(poly_ord, true));
  RCP<const product_basis> basis = rcp(new product_basis(bases));

  // Triple product tensor
  RCP<Cijk> cijk = basis->computeTripleProductTensor();

  // Kokkos triple product tensor
  kokkos_cijk_type kokkos_cijk =
    Stokhos::create_product_tensor<execution_space>(*basis, *cijk);

  return kokkos_cijk;
}

//
// Tests
//

// Stochastic discretizaiton used in the tests
const int stoch_dim = 2;
const int poly_ord = 3;

//
// Test vector addition
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, VectorAdd, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  Kokkos::setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // Comm
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  RCP<Tpetra_Vector> x1 = Tpetra::createVector<Scalar>(map);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> x1_view = x1->get1dViewNonConst();
  ArrayRCP<Scalar> x2_view = x2->get1dViewNonConst();
  Scalar val1(cijk), val2(cijk);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<pce_size; ++j) {
      val1.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(nrow, pce_size, row, j);
      val2.fastAccessCoeff(j) = 0.12345 * generate_vector_coefficient<BaseScalar,size_t>(nrow, pce_size, row, j);
    }
    x1_view[i] = val1;
    x2_view[i] = val2;
  }
  x1_view = Teuchos::null;
  x2_view = Teuchos::null;

  // Add
  Scalar alpha = 2.1;
  Scalar beta = 3.7;
  RCP<Tpetra_Vector> y = Tpetra::createVector<Scalar>(map);
  y->update(alpha, *x1, beta, *x2, Scalar(0.0));

  // y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check
  ArrayRCP<Scalar> y_view = y->get1dViewNonConst();
  Scalar val(cijk);
  BaseScalar tol = 1.0e-14;
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<pce_size; ++j) {
      BaseScalar v = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, pce_size, row, j);
      val.fastAccessCoeff(j) = alpha.coeff(0)*v + 0.12345*beta.coeff(0)*v;
    }
    TEST_EQUALITY( y_view[i].size(), pce_size );
    for (LocalOrdinal j=0; j<pce_size; ++j)
      TEST_FLOATING_EQUALITY( y_view[i].fastAccessCoeff(j), val.fastAccessCoeff(j), tol );
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

//
// Test vector dot product
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, VectorDot, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef typename Tpetra_Vector::dot_type dot_type;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // Comm
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  RCP<Tpetra_Vector> x1 = Tpetra::createVector<Scalar>(map);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> x1_view = x1->get1dViewNonConst();
  ArrayRCP<Scalar> x2_view = x2->get1dViewNonConst();
  Scalar val1(cijk), val2(cijk);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<pce_size; ++j) {
      val1.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(nrow, pce_size, row, j);
      val2.fastAccessCoeff(j) = 0.12345 * generate_vector_coefficient<BaseScalar,size_t>(nrow, pce_size, row, j);
    }
    x1_view[i] = val1;
    x2_view[i] = val2;
  }
  x1_view = Teuchos::null;
  x2_view = Teuchos::null;

  // Dot product
  dot_type dot = x1->dot(*x2);

  // Check

  // Local contribution
  dot_type local_val(0);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<pce_size; ++j) {
      BaseScalar v = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, pce_size, row, j);
      local_val += 0.12345 * v * v;
    }
  }

  // Global reduction
  dot_type val(0);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, local_val,
                     Teuchos::outArg(val));

  out << "dot = " << dot << " expected = " << val << std::endl;

  BaseScalar tol = 1.0e-14;
  TEST_FLOATING_EQUALITY( dot, val, tol );

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

//
// Test multi-vector addition
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, MultiVectorAdd, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // Comm
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  size_t ncol = 5;
  RCP<Tpetra_MultiVector> x1 = Tpetra::createMultiVector<Scalar>(map, ncol);
  RCP<Tpetra_MultiVector> x2 = Tpetra::createMultiVector<Scalar>(map, ncol);
  ArrayRCP< ArrayRCP<Scalar> > x1_view = x1->get2dViewNonConst();
  ArrayRCP< ArrayRCP<Scalar> > x2_view = x2->get2dViewNonConst();
  Scalar val1(cijk), val2(cijk);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol; ++j) {
      for (LocalOrdinal k=0; k<pce_size; ++k) {
        BaseScalar v =
          generate_multi_vector_coefficient<BaseScalar,size_t>(
            nrow, ncol, pce_size, row, j, k);
        val1.fastAccessCoeff(k) = v;
        val2.fastAccessCoeff(k) = 0.12345 * v;
      }
      x1_view[j][i] = val1;
      x2_view[j][i] = val2;
    }
  }
  x1_view = Teuchos::null;
  x2_view = Teuchos::null;

  // Add
  Scalar alpha = 2.1;
  Scalar beta = 3.7;
  RCP<Tpetra_MultiVector> y = Tpetra::createMultiVector<Scalar>(map, ncol);
  y->update(alpha, *x1, beta, *x2, Scalar(0.0));

  // y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check
  ArrayRCP< ArrayRCP<Scalar> > y_view = y->get2dViewNonConst();
  Scalar val(cijk);
  BaseScalar tol = 1.0e-14;
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol; ++j) {
      for (LocalOrdinal k=0; k<pce_size; ++k) {
        BaseScalar v = generate_multi_vector_coefficient<BaseScalar,size_t>(
          nrow, ncol, pce_size, row, j, k);
        val.fastAccessCoeff(k) = alpha.coeff(0)*v + 0.12345*beta.coeff(0)*v;
      }
      TEST_EQUALITY( y_view[j][i].size(), pce_size );
      for (LocalOrdinal k=0; k<pce_size; ++k)
        TEST_FLOATING_EQUALITY( y_view[j][i].fastAccessCoeff(k),
                                val.fastAccessCoeff(k), tol );
    }
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

//
// Test multi-vector dot product
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, MultiVectorDot, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector;
  typedef typename Tpetra_MultiVector::dot_type dot_type;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // Comm
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  size_t ncol = 5;
  RCP<Tpetra_MultiVector> x1 = Tpetra::createMultiVector<Scalar>(map, ncol);
  RCP<Tpetra_MultiVector> x2 = Tpetra::createMultiVector<Scalar>(map, ncol);
  ArrayRCP< ArrayRCP<Scalar> > x1_view = x1->get2dViewNonConst();
  ArrayRCP< ArrayRCP<Scalar> > x2_view = x2->get2dViewNonConst();
  Scalar val1(cijk), val2(cijk);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol; ++j) {
      for (LocalOrdinal k=0; k<pce_size; ++k) {
        BaseScalar v =
          generate_multi_vector_coefficient<BaseScalar,size_t>(
            nrow, ncol, pce_size, row, j, k);
        val1.fastAccessCoeff(k) = v;
        val2.fastAccessCoeff(k) = 0.12345 * v;
      }
      x1_view[j][i] = val1;
      x2_view[j][i] = val2;
    }
  }
  x1_view = Teuchos::null;
  x2_view = Teuchos::null;

  // Dot product
  Array<dot_type> dots(ncol);
  x1->dot(*x2, dots());

  // Check

  // Local contribution
  Array<dot_type> local_vals(ncol, dot_type(0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol; ++j) {
      for (LocalOrdinal k=0; k<pce_size; ++k) {
        BaseScalar v = generate_multi_vector_coefficient<BaseScalar,size_t>(
          nrow, ncol, pce_size, row, j, k);
        local_vals[j] += 0.12345 * v * v;
      }
    }
  }

  // Global reduction
  Array<dot_type> vals(ncol, dot_type(0));
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<int>(ncol),
                     local_vals.getRawPtr(), vals.getRawPtr());

  BaseScalar tol = 1.0e-14;
  for (size_t j=0; j<ncol; ++j) {
    out << "dots(" << j << ") = " << dots[j]
        << " expected(" << j << ") = " << vals[j] << std::endl;
    TEST_FLOATING_EQUALITY( dots[j], vals[j], tol );
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

//
// Test multi-vector dot product using subviews
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, MultiVectorDotSub, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector;
  typedef typename Tpetra_MultiVector::dot_type dot_type;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // Comm
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  ArrayView<const GlobalOrdinal> myGIDs = map->getNodeElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  size_t ncol = 5;
  RCP<Tpetra_MultiVector> x1 = Tpetra::createMultiVector<Scalar>(map, ncol);
  RCP<Tpetra_MultiVector> x2 = Tpetra::createMultiVector<Scalar>(map, ncol);
  ArrayRCP< ArrayRCP<Scalar> > x1_view = x1->get2dViewNonConst();
  ArrayRCP< ArrayRCP<Scalar> > x2_view = x2->get2dViewNonConst();
  Scalar val1(cijk), val2(cijk);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol; ++j) {
      for (LocalOrdinal k=0; k<pce_size; ++k) {
        BaseScalar v =
          generate_multi_vector_coefficient<BaseScalar,size_t>(
            nrow, ncol, pce_size, row, j, k);
        val1.fastAccessCoeff(k) = v;
        val2.fastAccessCoeff(k) = 0.12345 * v;
      }
      x1_view[j][i] = val1;
      x2_view[j][i] = val2;
    }
  }
  x1_view = Teuchos::null;
  x2_view = Teuchos::null;

  // Get subviews
  size_t ncol_sub = 2;
  Teuchos::Array<size_t> cols(ncol_sub);
  cols[0] = 4; cols[1] = 2;
  RCP<const Tpetra_MultiVector> x1_sub = x1->subView(cols());
  RCP<const Tpetra_MultiVector> x2_sub = x2->subView(cols());

  // Dot product
  Array<dot_type> dots(ncol_sub);
  x1_sub->dot(*x2_sub, dots());

  // Check

  // Local contribution
  Array<dot_type> local_vals(ncol_sub, dot_type(0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol_sub; ++j) {
      for (LocalOrdinal k=0; k<pce_size; ++k) {
        BaseScalar v = generate_multi_vector_coefficient<BaseScalar,size_t>(
          nrow, ncol, pce_size, row, cols[j], k);
        local_vals[j] += 0.12345 * v * v;
      }
    }
  }

  // Global reduction
  Array<dot_type> vals(ncol_sub, dot_type(0));
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM,
                     Teuchos::as<int>(ncol_sub), local_vals.getRawPtr(),
                     vals.getRawPtr());

  BaseScalar tol = 1.0e-14;
  for (size_t j=0; j<ncol_sub; ++j) {
    out << "dots(" << j << ") = " << dots[j]
        << " expected(" << j << ") = " << vals[j] << std::endl;
    TEST_FLOATING_EQUALITY( dots[j], vals[j], tol );
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

//
// Test matrix-vector multiplication for a simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, MatrixVectorMultiply, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // Build banded matrix
  GlobalOrdinal nrow = 13;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2), Tpetra::StaticProfile));
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
  Scalar val(cijk);
  for (size_t local_row=0; local_row<num_my_row; ++local_row) {
    const GlobalOrdinal row = myGIDs[local_row];
    const size_t num_col = row == nrow - 1 ? 1 : 2;
    for (size_t local_col=0; local_col<num_col; ++local_col) {
      const GlobalOrdinal col = row + local_col;
      columnIndices[local_col] = col;
      for (LocalOrdinal k=0; k<pce_size; ++k)
        val.fastAccessCoeff(k) =
          generate_matrix_coefficient<BaseScalar,size_t>(
            nrow, pce_size, row, col, k);
      vals[local_col] = val;
    }
    matrix->replaceGlobalValues(row, columnIndices(0,num_col), vals(0,num_col));
  }
  matrix->fillComplete();

  // Fill vector
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  for (size_t local_row=0; local_row<num_my_row; ++local_row) {
    const GlobalOrdinal row = myGIDs[local_row];
    for (LocalOrdinal j=0; j<pce_size; ++j)
      val.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, pce_size, row, j);
    x_view[local_row] = val;
  }
  x_view = Teuchos::null;

  // matrix->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //                  Teuchos::VERB_EXTREME);

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Multiply
  RCP<Tpetra_Vector> y = Tpetra::createVector<Scalar>(map);
  matrix->apply(*x, *y);

  // y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check
  ArrayRCP<Scalar> y_view = y->get1dViewNonConst();
  BaseScalar tol = 1.0e-14;
  typename Cijk::HostMirror host_cijk =
    Kokkos::create_mirror_view(cijk);
  Kokkos::deep_copy(host_cijk, cijk);
  for (size_t local_row=0; local_row<num_my_row; ++local_row) {
    const GlobalOrdinal row = myGIDs[local_row];
    const size_t num_col = row == nrow - 1 ? 1 : 2;
    val = 0.0;
    for (size_t local_col=0; local_col<num_col; ++local_col) {
      const GlobalOrdinal col = row + local_col;
      for (LocalOrdinal i=0; i<pce_size; ++i) {
        const LocalOrdinal num_entry = host_cijk.num_entry(i);
        const LocalOrdinal entry_beg = host_cijk.entry_begin(i);
        const LocalOrdinal entry_end = entry_beg + num_entry;
        BaseScalar tmp = 0;
        for (LocalOrdinal entry = entry_beg; entry < entry_end; ++entry) {
          const LocalOrdinal j = host_cijk.coord(entry,0);
          const LocalOrdinal k = host_cijk.coord(entry,1);
          const BaseScalar a_j =
            generate_matrix_coefficient<BaseScalar,size_t>(
              nrow, pce_size, row, col, j);
          const BaseScalar a_k =
            generate_matrix_coefficient<BaseScalar,size_t>(
              nrow, pce_size, row, col, k);
          const BaseScalar x_j =
            generate_vector_coefficient<BaseScalar,size_t>(
              nrow, pce_size, col, j);
          const BaseScalar x_k =
            generate_vector_coefficient<BaseScalar,size_t>(
              nrow, pce_size, col, k);
          tmp += host_cijk.value(entry) * ( a_j * x_k + a_k * x_j );
        }
        val.fastAccessCoeff(i) += tmp;
      }
    }
    TEST_EQUALITY( y_view[local_row].size(), pce_size );
    for (LocalOrdinal i=0; i<pce_size; ++i)
      TEST_FLOATING_EQUALITY( y_view[local_row].fastAccessCoeff(i),
                              val.fastAccessCoeff(i), tol );
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

//
// Test matrix-multi-vector multiplication for a simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, MatrixMultiVectorMultiply, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2), Tpetra::StaticProfile));
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
  Scalar val(cijk);
  for (size_t local_row=0; local_row<num_my_row; ++local_row) {
    const GlobalOrdinal row = myGIDs[local_row];
    const size_t num_col = row == nrow - 1 ? 1 : 2;
    for (size_t local_col=0; local_col<num_col; ++local_col) {
      const GlobalOrdinal col = row + local_col;
      columnIndices[local_col] = col;
      for (LocalOrdinal k=0; k<pce_size; ++k)
        val.fastAccessCoeff(k) =
          generate_matrix_coefficient<BaseScalar,size_t>(
            nrow, pce_size, row, col, k);
      vals[local_col] = val;
    }
    matrix->replaceGlobalValues(row, columnIndices(0,num_col), vals(0,num_col));
  }
  matrix->fillComplete();

  // Fill multi-vector
  size_t ncol = 5;
  RCP<Tpetra_MultiVector> x = Tpetra::createMultiVector<Scalar>(map, ncol);
  ArrayRCP< ArrayRCP<Scalar> > x_view = x->get2dViewNonConst();
  for (size_t local_row=0; local_row<num_my_row; ++local_row) {
    const GlobalOrdinal row = myGIDs[local_row];
    for (size_t col=0; col<ncol; ++col) {
      for (LocalOrdinal k=0; k<pce_size; ++k) {
        BaseScalar v =
          generate_multi_vector_coefficient<BaseScalar,size_t>(
            nrow, ncol, pce_size, row, col, k);
        val.fastAccessCoeff(k) = v;
      }
      x_view[col][local_row] = val;
    }
  }
  x_view = Teuchos::null;

  // matrix->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //                  Teuchos::VERB_EXTREME);

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Multiply
  RCP<Tpetra_MultiVector> y = Tpetra::createMultiVector<Scalar>(map, ncol);
  matrix->apply(*x, *y);

  // y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check
  ArrayRCP< ArrayRCP<Scalar> > y_view = y->get2dViewNonConst();
  BaseScalar tol = 1.0e-14;
  typename Cijk::HostMirror host_cijk =
    Kokkos::create_mirror_view(cijk);
  Kokkos::deep_copy(host_cijk, cijk);
  for (size_t local_row=0; local_row<num_my_row; ++local_row) {
    const GlobalOrdinal row = myGIDs[local_row];
    for (size_t xcol=0; xcol<ncol; ++xcol) {
      const size_t num_col = row == nrow - 1 ? 1 : 2;
      val = 0.0;
      for (size_t local_col=0; local_col<num_col; ++local_col) {
        const GlobalOrdinal col = row + local_col;
        for (LocalOrdinal i=0; i<pce_size; ++i) {
          const LocalOrdinal num_entry = host_cijk.num_entry(i);
          const LocalOrdinal entry_beg = host_cijk.entry_begin(i);
          const LocalOrdinal entry_end = entry_beg + num_entry;
          BaseScalar tmp = 0;
          for (LocalOrdinal entry = entry_beg; entry < entry_end; ++entry) {
            const LocalOrdinal j = host_cijk.coord(entry,0);
            const LocalOrdinal k = host_cijk.coord(entry,1);
            const BaseScalar a_j =
              generate_matrix_coefficient<BaseScalar,size_t>(
                nrow, pce_size, row, col, j);
            const BaseScalar a_k =
              generate_matrix_coefficient<BaseScalar,size_t>(
                nrow, pce_size, row, col, k);
            const BaseScalar x_j =
              generate_multi_vector_coefficient<BaseScalar,size_t>(
                nrow, ncol, pce_size, col, xcol, j);
            const BaseScalar x_k =
              generate_multi_vector_coefficient<BaseScalar,size_t>(
                nrow, ncol, pce_size, col, xcol, k);
            tmp += host_cijk.value(entry) * ( a_j * x_k + a_k * x_j );
          }
          val.fastAccessCoeff(i) += tmp;
        }
      }
      TEST_EQUALITY( y_view[xcol][local_row].size(), pce_size );
      for (LocalOrdinal i=0; i<pce_size; ++i)
        TEST_FLOATING_EQUALITY( y_view[xcol][local_row].fastAccessCoeff(i),
                                val.fastAccessCoeff(i), tol );
    }
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

//
// Test flattening UQ::PCE matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, Flatten, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  typedef Tpetra::Vector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_Vector;
  typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_CrsMatrix;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2), Tpetra::StaticProfile));
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
  Scalar val(cijk);
  for (size_t local_row=0; local_row<num_my_row; ++local_row) {
    const GlobalOrdinal row = myGIDs[local_row];
    const size_t num_col = row == nrow - 1 ? 1 : 2;
    for (size_t local_col=0; local_col<num_col; ++local_col) {
      const GlobalOrdinal col = row + local_col;
      columnIndices[local_col] = col;
      for (LocalOrdinal k=0; k<pce_size; ++k)
        val.fastAccessCoeff(k) =
          generate_matrix_coefficient<BaseScalar,size_t>(
            nrow, pce_size, row, col, k);
      vals[local_col] = val;
    }
    matrix->replaceGlobalValues(row, columnIndices(0,num_col), vals(0,num_col));
  }
  matrix->fillComplete();

  // Fill vector
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<pce_size; ++j)
      val.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, pce_size, row, j);
    x_view[i] = val;
  }

  // Multiply
  RCP<Tpetra_Vector> y = Tpetra::createVector<Scalar>(map);
  matrix->apply(*x, *y);

  /*
  graph->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
                  Teuchos::VERB_EXTREME);

  matrix->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
                   Teuchos::VERB_EXTREME);

  x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
              Teuchos::VERB_EXTREME);

  y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
              Teuchos::VERB_EXTREME);
  */

  // Flatten matrix
  RCP<const Tpetra_Map> flat_x_map, flat_y_map;
  RCP<const Tpetra_CrsGraph> flat_graph, cijk_graph;
  flat_graph =
    Stokhos::create_flat_pce_graph(*graph, cijk, flat_x_map, flat_y_map,
                                   cijk_graph, pce_size);
  RCP<Flat_Tpetra_CrsMatrix> flat_matrix =
    Stokhos::create_flat_matrix(*matrix, flat_graph, cijk_graph, cijk);

  // Multiply with flattened matix
  RCP<Tpetra_Vector> y2 = Tpetra::createVector<Scalar>(map);
  RCP<Flat_Tpetra_Vector> flat_x =
    Stokhos::create_flat_vector_view(*x, flat_x_map);
  RCP<Flat_Tpetra_Vector> flat_y =
    Stokhos::create_flat_vector_view(*y2, flat_y_map);
  flat_matrix->apply(*flat_x, *flat_y);

  /*
  cijk_graph->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
                       Teuchos::VERB_EXTREME);

  flat_graph->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
                       Teuchos::VERB_EXTREME);

  flat_matrix->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
                        Teuchos::VERB_EXTREME);

  flat_x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
                   Teuchos::VERB_EXTREME);

  flat_y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
                   Teuchos::VERB_EXTREME);
  */

  // Check
  BaseScalar tol = 1.0e-14;
  ArrayRCP<Scalar> y_view = y->get1dViewNonConst();
  ArrayRCP<Scalar> y2_view = y2->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    TEST_EQUALITY( y_view[i].size(), pce_size );
    TEST_EQUALITY( y2_view[i].size(), pce_size );
    for (LocalOrdinal j=0; j<pce_size; ++j)
      TEST_FLOATING_EQUALITY( y_view[i].fastAccessCoeff(j),
                              y2_view[i].fastAccessCoeff(j), tol );
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

//
// Test simple CG solve without preconditioning for a 1-D Laplacian matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, SimpleCG, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // 1-D Laplacian matrix
  GlobalOrdinal nrow = 10;
  BaseScalar h = 1.0 / static_cast<BaseScalar>(nrow-1);
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(3), Tpetra::StaticProfile));
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
  Scalar a_val(cijk);
  for (LocalOrdinal j=0; j<pce_size; ++j) {
    a_val.fastAccessCoeff(j) =
      BaseScalar(1.0) + BaseScalar(1.0) / BaseScalar(j+1);
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

  // matrix->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //                  Teuchos::VERB_EXTREME);

  // Fill RHS vector
  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> b_view = b->get1dViewNonConst();
  Scalar b_val(cijk);
  for (LocalOrdinal j=0; j<pce_size; ++j) {
    b_val.fastAccessCoeff(j) =
      BaseScalar(2.0) - BaseScalar(1.0) / BaseScalar(j+1);
  }
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row == 0 || row == nrow-1)
      b_view[i] = Scalar(0.0);
    else
      b_view[i] = b_val * (h*h);
  }

  // b->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Solve
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  typedef Kokkos::Details::ArithTraits<BaseScalar> BST;
  typedef typename BST::mag_type base_mag_type;
  typedef typename Tpetra_Vector::mag_type mag_type;
  base_mag_type btol = 1e-9;
  mag_type tol = btol;
  int max_its = 1000;
  bool solved = Stokhos::CG_Solve(*matrix, *x, *b, tol, max_its,
                                  out.getOStream().get());
  TEST_EQUALITY_CONST( solved, true );

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check by solving flattened system
  typedef Tpetra::Vector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_Vector;
  typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_CrsMatrix;
  RCP<const Tpetra_Map> flat_x_map, flat_b_map;
  RCP<const Tpetra_CrsGraph> flat_graph, cijk_graph;
  flat_graph =
    Stokhos::create_flat_pce_graph(*graph, cijk, flat_x_map, flat_b_map,
                                   cijk_graph, pce_size);
  RCP<Flat_Tpetra_CrsMatrix> flat_matrix =
    Stokhos::create_flat_matrix(*matrix, flat_graph, cijk_graph, cijk);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  RCP<Flat_Tpetra_Vector> flat_x =
    Stokhos::create_flat_vector_view(*x2, flat_x_map);
  RCP<Flat_Tpetra_Vector> flat_b =
    Stokhos::create_flat_vector_view(*b, flat_b_map);
  bool solved_flat = Stokhos::CG_Solve(*flat_matrix, *flat_x, *flat_b,
                                       tol, max_its, out.getOStream().get());
  TEST_EQUALITY_CONST( solved_flat, true );

  btol = 500*btol;
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  ArrayRCP<Scalar> x2_view = x2->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    TEST_EQUALITY( x_view[i].size(),  pce_size );
    TEST_EQUALITY( x2_view[i].size(), pce_size );

    // Set small values to zero
    Scalar v = x_view[i];
    Scalar v2 = x2_view[i];
    for (LocalOrdinal j=0; j<pce_size; ++j) {
      if (j < v.size() && BST::abs(v.coeff(j)) < btol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
      if (j < v2.size() && BST::abs(v2.coeff(j)) < btol)
        v2.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<pce_size; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), v2.coeff(j), btol);
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());

}

#if defined(HAVE_STOKHOS_MUELU) && defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_STOKHOS_IFPACK2) && defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)

//
// Test simple CG solve with MueLu preconditioning for a 1-D Laplacian matrix
//
// Currently requires ETI since the specializations needed for mean-based
// are only brought in with ETI
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, SimplePCG_Muelu, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::getParametersFromXmlFile;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // 1-D Laplacian matrix
  GlobalOrdinal nrow = 10;
  BaseScalar h = 1.0 / static_cast<BaseScalar>(nrow-1);
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(3), Tpetra::StaticProfile));
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
  Scalar a_val(cijk);
  for (LocalOrdinal j=0; j<pce_size; ++j) {
    a_val.fastAccessCoeff(j) =
      BaseScalar(1.0) + BaseScalar(1.0) / BaseScalar(j+1);
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
  Scalar b_val(cijk);
  for (LocalOrdinal j=0; j<pce_size; ++j) {
    b_val.fastAccessCoeff(j) =
      BaseScalar(2.0) - BaseScalar(1.0) / BaseScalar(j+1);
  }
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row == 0 || row == nrow-1)
      b_view[i] = Scalar(0.0);
    else
      b_view[i] = b_val * (h*h);
  }

  // Create preconditioner
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
  RCP<ParameterList> muelu_params =
    getParametersFromXmlFile("muelu_cheby.xml");
#if USE_SCALAR_MEAN_BASED_PREC
  typedef Tpetra::Operator<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Scalar_OP;
  typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Scalar_Tpetra_CrsMatrix;
  RCP<Scalar_Tpetra_CrsMatrix> mean_matrix =
    Stokhos::build_mean_scalar_matrix(*matrix);
  RCP<Scalar_OP> mean_matrix_op = mean_matrix;
  RCP<Scalar_OP> M_s =
    MueLu::CreateTpetraPreconditioner<BaseScalar,LocalOrdinal,GlobalOrdinal,Node>(mean_matrix_op, *muelu_params);
  RCP<OP> M = rcp(new Stokhos::MeanBasedTpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>(M_s));
#else
  Cijk mean_cijk =
    Stokhos::create_mean_based_product_tensor<typename Storage::execution_space,typename Storage::ordinal_type,BaseScalar>();
  Kokkos::setGlobalCijkTensor(mean_cijk);
  RCP<Tpetra_CrsMatrix> mean_matrix = Stokhos::build_mean_matrix(*matrix);
  RCP<OP> mean_matrix_op = mean_matrix;
  RCP<OP> M =
    MueLu::CreateTpetraPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>(mean_matrix_op, *muelu_params);
  Kokkos::setGlobalCijkTensor(cijk);
#endif

  // Solve
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  typedef Kokkos::Details::ArithTraits<BaseScalar> BST;
  typedef typename BST::mag_type base_mag_type;
  typedef typename Tpetra_Vector::mag_type mag_type;
  base_mag_type btol = 1e-9;
  mag_type tol = btol;
  int max_its = 1000;
  bool solved = Stokhos::PCG_Solve(*matrix, *x, *b, *M, tol, max_its,
                                   out.getOStream().get());
  TEST_EQUALITY_CONST( solved, true );

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check by solving flattened system
   typedef Tpetra::Vector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_Vector;
  typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_CrsMatrix;
  RCP<const Tpetra_Map> flat_x_map, flat_b_map;
  RCP<const Tpetra_CrsGraph> flat_graph, cijk_graph;
  flat_graph =
    Stokhos::create_flat_pce_graph(*graph, cijk, flat_x_map, flat_b_map,
                                   cijk_graph, pce_size);
  RCP<Flat_Tpetra_CrsMatrix> flat_matrix =
    Stokhos::create_flat_matrix(*matrix, flat_graph, cijk_graph, cijk);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  RCP<Flat_Tpetra_Vector> flat_x =
    Stokhos::create_flat_vector_view(*x2, flat_x_map);
  RCP<Flat_Tpetra_Vector> flat_b =
    Stokhos::create_flat_vector_view(*b, flat_b_map);
  // typedef Tpetra::Operator<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FlatPrec;
  // RCP<FlatPrec> flat_M =
  //   MueLu::CreateTpetraPreconditioner<BaseScalar,LocalOrdinal,GlobalOrdinal,Node>(flat_matrix, *muelu_params);
  // bool solved_flat = Stokhos::PCG_Solve(*flat_matrix, *flat_x, *flat_b, *flat_M,
  //                                      tol, max_its, out.getOStream().get());
  bool solved_flat = Stokhos::CG_Solve(*flat_matrix, *flat_x, *flat_b,
                                       tol, max_its, out.getOStream().get());
  TEST_EQUALITY_CONST( solved_flat, true );

  btol = 500*btol;
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  ArrayRCP<Scalar> x2_view = x2->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    TEST_EQUALITY( x_view[i].size(),  pce_size );
    TEST_EQUALITY( x2_view[i].size(), pce_size );

    // Set small values to zero
    Scalar v = x_view[i];
    Scalar v2 = x2_view[i];
    for (LocalOrdinal j=0; j<pce_size; ++j) {
      if (j < v.size() && BST::abs(v.coeff(j)) < btol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
      if (j < v2.size() && BST::abs(v2.coeff(j)) < btol)
        v2.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<pce_size; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), v2.coeff(j), btol);
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());

}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, SimplePCG_Muelu, Storage, LocalOrdinal, GlobalOrdinal, Node ) {}

#endif

#if defined(HAVE_STOKHOS_BELOS)

//
// Test Belos GMRES solve for a simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, BelosGMRES, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2), Tpetra::StaticProfile));
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
  Scalar val(cijk);
  for (size_t local_row=0; local_row<num_my_row; ++local_row) {
    const GlobalOrdinal row = myGIDs[local_row];
    const size_t num_col = row == nrow - 1 ? 1 : 2;
    for (size_t local_col=0; local_col<num_col; ++local_col) {
      const GlobalOrdinal col = row + local_col;
      columnIndices[local_col] = col;
      for (LocalOrdinal k=0; k<pce_size; ++k)
        val.fastAccessCoeff(k) =
          BaseScalar(1.0) + BaseScalar(1.0) / BaseScalar(k+1);
      vals[local_col] = val;
    }
    matrix->replaceGlobalValues(row, columnIndices(0,num_col), vals(0,num_col));
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
  typedef Belos::LinearProblem<BelosScalar,MV,OP> BLinProb;
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  RCP< BLinProb > problem = rcp(new BLinProb(matrix, x, b));
  RCP<ParameterList> belosParams = rcp(new ParameterList);
  typename ST::magnitudeType tol = 1e-9;
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

  // Check by solving flattened system
  typedef Tpetra::Vector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_Vector;
  typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_CrsMatrix;
  RCP<const Tpetra_Map> flat_x_map, flat_b_map;
  RCP<const Tpetra_CrsGraph> flat_graph, cijk_graph;
  flat_graph =
    Stokhos::create_flat_pce_graph(*graph, cijk, flat_x_map, flat_b_map,
                                   cijk_graph, pce_size);
  RCP<Flat_Tpetra_CrsMatrix> flat_matrix =
    Stokhos::create_flat_matrix(*matrix, flat_graph, cijk_graph, cijk);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  RCP<Flat_Tpetra_Vector> flat_x =
    Stokhos::create_flat_vector_view(*x2, flat_x_map);
  RCP<Flat_Tpetra_Vector> flat_b =
    Stokhos::create_flat_vector_view(*b, flat_b_map);
  typedef Tpetra::MultiVector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FMV;
  typedef Tpetra::Operator<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FOP;
  typedef Belos::LinearProblem<BelosScalar,FMV,FOP> FBLinProb;
  RCP< FBLinProb > flat_problem =
    rcp(new FBLinProb(flat_matrix, flat_x, flat_b));
  RCP<Belos::SolverManager<BelosScalar,FMV,FOP> > flat_solver =
    rcp(new Belos::PseudoBlockGmresSolMgr<BelosScalar,FMV,FOP>(flat_problem,
                                                               belosParams));
  flat_problem->setProblem();
  Belos::ReturnType flat_ret = flat_solver->solve();
  TEST_EQUALITY_CONST( flat_ret, Belos::Converged );

  typename ST::magnitudeType btol = 100*tol;
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  ArrayRCP<Scalar> x2_view = x2->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    TEST_EQUALITY( x_view[i].size(),  pce_size );
    TEST_EQUALITY( x2_view[i].size(), pce_size );

    // Set small values to zero
    Scalar v = x_view[i];
    Scalar v2 = x2_view[i];
    for (LocalOrdinal j=0; j<pce_size; ++j) {
      if (j < v.size() && ST::magnitude(v.coeff(j)) < btol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
      if (j < v2.size() && ST::magnitude(v2.coeff(j)) < btol)
        v2.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<pce_size; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), v2.coeff(j), btol);
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, BelosGMRES, Storage, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

// Test currently doesn't work (in serial) because of our bad division strategy

#if defined(HAVE_STOKHOS_BELOS) && defined(HAVE_STOKHOS_IFPACK2)

//
// Test Belos GMRES solve with Ifpack2 RILUK preconditioning for a
// simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, BelosGMRES_RILUK, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2), Tpetra::StaticProfile));
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
  Scalar val(cijk);
  for (size_t local_row=0; local_row<num_my_row; ++local_row) {
    const GlobalOrdinal row = myGIDs[local_row];
    const size_t num_col = row == nrow - 1 ? 1 : 2;
    for (size_t local_col=0; local_col<num_col; ++local_col) {
      const GlobalOrdinal col = row + local_col;
      columnIndices[local_col] = col;
      for (LocalOrdinal k=0; k<pce_size; ++k)
        val.fastAccessCoeff(k) =
          BaseScalar(1.0) + BaseScalar(1.0) / BaseScalar(k+1);
      vals[local_col] = val;
    }
    matrix->replaceGlobalValues(row, columnIndices(0,num_col), vals(0,num_col));
  }
  matrix->fillComplete();

  // Create mean matrix for preconditioning
  RCP<Tpetra_CrsMatrix> mean_matrix = Stokhos::build_mean_matrix(*matrix);

  // Fill RHS vector
  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  ArrayRCP<Scalar> b_view = b->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    b_view[i] = Scalar(1.0);
  }

  // Create preconditioner
  typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Prec;
  Ifpack2::Factory factory;
  RCP<Prec> M = factory.create<Tpetra_CrsMatrix>("RILUK", mean_matrix);
  M->initialize();
  M->compute();

  // Solve
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
  typedef BaseScalar BelosScalar;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
  typedef Belos::LinearProblem<BelosScalar,MV,OP> BLinProb;
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  RCP< BLinProb > problem = rcp(new BLinProb(matrix, x, b));
  problem->setRightPrec(M);
  problem->setProblem();
  RCP<ParameterList> belosParams = rcp(new ParameterList);
  typename ST::magnitudeType tol = 1e-9;
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

  // Check by solving flattened system
  typedef Tpetra::Vector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_Vector;
  typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_CrsMatrix;
  RCP<const Tpetra_Map> flat_x_map, flat_b_map;
  RCP<const Tpetra_CrsGraph> flat_graph, cijk_graph;
  flat_graph =
    Stokhos::create_flat_pce_graph(*graph, cijk, flat_x_map, flat_b_map,
                                   cijk_graph, pce_size);
  RCP<Flat_Tpetra_CrsMatrix> flat_matrix =
    Stokhos::create_flat_matrix(*matrix, flat_graph, cijk_graph, cijk);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  RCP<Flat_Tpetra_Vector> flat_x =
    Stokhos::create_flat_vector_view(*x2, flat_x_map);
  RCP<Flat_Tpetra_Vector> flat_b =
    Stokhos::create_flat_vector_view(*b, flat_b_map);
  typedef Tpetra::MultiVector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FMV;
  typedef Tpetra::Operator<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FOP;
  typedef Belos::LinearProblem<BelosScalar,FMV,FOP> FBLinProb;
  RCP< FBLinProb > flat_problem =
    rcp(new FBLinProb(flat_matrix, flat_x, flat_b));
  RCP<Belos::SolverManager<BelosScalar,FMV,FOP> > flat_solver =
    rcp(new Belos::PseudoBlockGmresSolMgr<BelosScalar,FMV,FOP>(flat_problem,
                                                               belosParams));
  flat_problem->setProblem();
  Belos::ReturnType flat_ret = flat_solver->solve();
  TEST_EQUALITY_CONST( flat_ret, Belos::Converged );

  typename ST::magnitudeType btol = 100*tol;
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  ArrayRCP<Scalar> x2_view = x2->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    TEST_EQUALITY( x_view[i].size(),  pce_size );
    TEST_EQUALITY( x2_view[i].size(), pce_size );

    // Set small values to zero
    Scalar v = x_view[i];
    Scalar v2 = x2_view[i];
    for (LocalOrdinal j=0; j<pce_size; ++j) {
      if (j < v.size() && ST::magnitude(v.coeff(j)) < btol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
      if (j < v2.size() && ST::magnitude(v2.coeff(j)) < btol)
        v2.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<pce_size; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), v2.coeff(j), btol);
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, BelosGMRES_RILUK, Storage, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

#if defined(HAVE_STOKHOS_BELOS) && defined(HAVE_STOKHOS_IFPACK2) && defined(HAVE_STOKHOS_MUELU) && defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)

//
// Test Belos CG solve with MueLu preconditioning for a 1-D Laplacian matrix
//
// Currently requires ETI since the specializations needed for mean-based
// are only brought in with ETI
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, BelosCG_Muelu, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::getParametersFromXmlFile;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // 1-D Laplacian matrix
  GlobalOrdinal nrow = 10;
  BaseScalar h = 1.0 / static_cast<BaseScalar>(nrow-1);
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(3), Tpetra::StaticProfile));
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
  Scalar a_val(cijk);
  for (LocalOrdinal j=0; j<pce_size; ++j) {
    a_val.fastAccessCoeff(j) =
      BaseScalar(1.0) + BaseScalar(1.0) / BaseScalar(j+1);
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
  Scalar b_val(cijk);
  for (LocalOrdinal j=0; j<pce_size; ++j) {
    b_val.fastAccessCoeff(j) =
      BaseScalar(2.0) - BaseScalar(1.0) / BaseScalar(j+1);
  }
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row == 0 || row == nrow-1)
      b_view[i] = Scalar(0.0);
    else
      b_view[i] = b_val * (h*h);
  }

  // Create preconditioner
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
  RCP<ParameterList> muelu_params =
    getParametersFromXmlFile("muelu_cheby.xml");
#if USE_SCALAR_MEAN_BASED_PREC
  typedef Tpetra::Operator<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Scalar_OP;
  typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Scalar_Tpetra_CrsMatrix;
  RCP<Scalar_Tpetra_CrsMatrix> mean_matrix =
    Stokhos::build_mean_scalar_matrix(*matrix);
  RCP<Scalar_OP> mean_matrix_op = mean_matrix;
  RCP<Scalar_OP> M_s =
    MueLu::CreateTpetraPreconditioner<BaseScalar,LocalOrdinal,GlobalOrdinal,Node>(mean_matrix_op, *muelu_params);
  RCP<OP> M = rcp(new Stokhos::MeanBasedTpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>(M_s));
#else
  Cijk mean_cijk =
    Stokhos::create_mean_based_product_tensor<typename Storage::execution_space,typename Storage::ordinal_type,BaseScalar>();
  Kokkos::setGlobalCijkTensor(mean_cijk);
  RCP<Tpetra_CrsMatrix> mean_matrix = Stokhos::build_mean_matrix(*matrix);
  RCP<OP> mean_matrix_op = mean_matrix;
  RCP<OP> M =
    MueLu::CreateTpetraPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>(mean_matrix_op, *muelu_params);
  Kokkos::setGlobalCijkTensor(cijk);
#endif

  // Solve
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
  typedef BaseScalar BelosScalar;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Belos::LinearProblem<BelosScalar,MV,OP> BLinProb;
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  RCP< BLinProb > problem = rcp(new BLinProb(matrix, x, b));
  problem->setRightPrec(M);
  problem->setProblem();
  RCP<ParameterList> belosParams = rcp(new ParameterList);
  typename ST::magnitudeType tol = 1e-9;
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
  // RCP<Belos::SolverManager<BelosScalar,MV,OP> > solver =
  //   rcp(new Belos::PseudoBlockCGSolMgr<BelosScalar,MV,OP>(problem, belosParams));
  Belos::ReturnType ret = solver->solve();
  TEST_EQUALITY_CONST( ret, Belos::Converged );

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check by solving flattened system
  typedef Tpetra::Vector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_Vector;
  typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_CrsMatrix;
  RCP<const Tpetra_Map> flat_x_map, flat_b_map;
  RCP<const Tpetra_CrsGraph> flat_graph, cijk_graph;
  flat_graph =
    Stokhos::create_flat_pce_graph(*graph, cijk, flat_x_map, flat_b_map,
                                   cijk_graph, pce_size);
  RCP<Flat_Tpetra_CrsMatrix> flat_matrix =
    Stokhos::create_flat_matrix(*matrix, flat_graph, cijk_graph, cijk);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  RCP<Flat_Tpetra_Vector> flat_x =
    Stokhos::create_flat_vector_view(*x2, flat_x_map);
  RCP<Flat_Tpetra_Vector> flat_b =
    Stokhos::create_flat_vector_view(*b, flat_b_map);
  typedef Tpetra::MultiVector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FMV;
  typedef Tpetra::Operator<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FOP;
  typedef Belos::LinearProblem<BelosScalar,FMV,FOP> FBLinProb;
  RCP< FBLinProb > flat_problem =
    rcp(new FBLinProb(flat_matrix, flat_x, flat_b));
  RCP<Belos::SolverManager<BelosScalar,FMV,FOP> > flat_solver =
    rcp(new Belos::PseudoBlockGmresSolMgr<BelosScalar,FMV,FOP>(flat_problem,
                                                               belosParams));
  // RCP<Belos::SolverManager<BelosScalar,FMV,FOP> > flat_solver =
  //   rcp(new Belos::PseudoBlockCGSolMgr<BelosScalar,FMV,FOP>(flat_problem,
  //                                                           belosParams));
  flat_problem->setProblem();
  Belos::ReturnType flat_ret = flat_solver->solve();
  TEST_EQUALITY_CONST( flat_ret, Belos::Converged );

  typename ST::magnitudeType btol = 100*tol;
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  ArrayRCP<Scalar> x2_view = x2->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    TEST_EQUALITY( x_view[i].size(),  pce_size );
    TEST_EQUALITY( x2_view[i].size(), pce_size );

    // Set small values to zero
    Scalar v = x_view[i];
    Scalar v2 = x2_view[i];
    for (LocalOrdinal j=0; j<pce_size; ++j) {
      if (j < v.size() && ST::magnitude(v.coeff(j)) < btol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
      if (j < v2.size() && ST::magnitude(v2.coeff(j)) < btol)
        v2.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<pce_size; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), v2.coeff(j), btol);
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());

}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, BelosCG_Muelu, Storage, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

#if defined(HAVE_STOKHOS_AMESOS2)

//
// Test Amesos2 solve for a 1-D Laplacian matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, Amesos2, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typedef typename Scalar::cijk_type Cijk;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Cijk
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);
  setGlobalCijkTensor(cijk);
  LocalOrdinal pce_size = cijk.dimension();

  // 1-D Laplacian matrix
  GlobalOrdinal nrow = 10;
  BaseScalar h = 1.0 / static_cast<BaseScalar>(nrow-1);
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
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
  Scalar a_val(cijk);
  for (LocalOrdinal j=0; j<pce_size; ++j) {
    a_val.fastAccessCoeff(j) =
      BaseScalar(1.0) + BaseScalar(1.0) / BaseScalar(j+1);
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
  Scalar b_val(cijk);
  for (LocalOrdinal j=0; j<pce_size; ++j) {
    b_val.fastAccessCoeff(j) =
      BaseScalar(2.0) - BaseScalar(1.0) / BaseScalar(j+1);
  }
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row == 0 || row == nrow-1)
      b_view[i] = Scalar(0.0);
    else
      b_view[i] = b_val * (h*h);
  }

  // Solve
  typedef Amesos2::Solver<Tpetra_CrsMatrix,Tpetra_MultiVector> Solver;
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  std::string solver_name;
#if defined(HAVE_AMESOS2_BASKER)
  solver_name = "basker";
#elif defined(HAVE_AMESOS2_KLU2)
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
  out << "Solving linear system with " << solver_name << std::endl;
  RCP<Solver> solver = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>(
    solver_name, matrix, x, b);
  solver->solve();

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check by solving flattened system
  typedef Tpetra::Vector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_Vector;
  typedef Tpetra::MultiVector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_MultiVector;
  typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> Flat_Tpetra_CrsMatrix;
  RCP<const Tpetra_Map> flat_x_map, flat_b_map;
  RCP<const Tpetra_CrsGraph> flat_graph, cijk_graph;
  flat_graph =
    Stokhos::create_flat_pce_graph(*graph, cijk, flat_x_map, flat_b_map,
                                   cijk_graph, pce_size);
  RCP<Flat_Tpetra_CrsMatrix> flat_matrix =
    Stokhos::create_flat_matrix(*matrix, flat_graph, cijk_graph, cijk);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  RCP<Flat_Tpetra_Vector> flat_x =
    Stokhos::create_flat_vector_view(*x2, flat_x_map);
  RCP<Flat_Tpetra_Vector> flat_b =
    Stokhos::create_flat_vector_view(*b, flat_b_map);
  typedef Amesos2::Solver<Flat_Tpetra_CrsMatrix,Flat_Tpetra_MultiVector> Flat_Solver;
  RCP<Flat_Solver> flat_solver =
    Amesos2::create<Flat_Tpetra_CrsMatrix,Flat_Tpetra_MultiVector>(
      solver_name, flat_matrix, flat_x, flat_b);
  flat_solver->solve();

  typedef Kokkos::Details::ArithTraits<BaseScalar> ST;
  typename ST::mag_type btol = 1e-12;
  ArrayRCP<Scalar> x_view = x->get1dViewNonConst();
  ArrayRCP<Scalar> x2_view = x2->get1dViewNonConst();
  for (size_t i=0; i<num_my_row; ++i) {
    TEST_EQUALITY( x_view[i].size(),  pce_size );
    TEST_EQUALITY( x2_view[i].size(), pce_size );

    // Set small values to zero
    Scalar v = x_view[i];
    Scalar v2 = x2_view[i];
    for (LocalOrdinal j=0; j<pce_size; ++j) {
      if (j < v.size() && ST::magnitude(v.coeff(j)) < btol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
      if (j < v2.size() && ST::magnitude(v2.coeff(j)) < btol)
        v2.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<pce_size; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), v2.coeff(j), btol);
  }

  // Clear global tensor
  Kokkos::setGlobalCijkTensor(Cijk());
}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_PCE, Amesos2, Storage, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

#define CRSMATRIX_UQ_PCE_TESTS_SLGN(S, LO, GO, N)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, VectorAdd, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, VectorDot, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, MultiVectorAdd, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, MultiVectorDot, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, MultiVectorDotSub, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, MatrixVectorMultiply, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, MatrixMultiVectorMultiply, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, Flatten, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, SimpleCG, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, SimplePCG_Muelu, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, BelosGMRES, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, BelosGMRES_RILUK, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, BelosCG_Muelu, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_PCE, Amesos2, S, LO, GO, N )

#define CRSMATRIX_UQ_PCE_TESTS_N(N)                                     \
  typedef Stokhos::DeviceForNode2<N>::type Device;                      \
  typedef Stokhos::DynamicStorage<int,double,Device::execution_space> DS; \
  CRSMATRIX_UQ_PCE_TESTS_SLGN(DS, int, int, N)
