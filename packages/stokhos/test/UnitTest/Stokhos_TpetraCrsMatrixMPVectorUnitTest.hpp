// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHelpers.hpp"
#include "Stokhos_UnitTestHelpers.hpp"
#include "Stokhos_Ensemble_Sizes.hpp"

// Teuchos
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

// Tpetra
#include "Stokhos_Tpetra_MP_Vector.hpp"
#include "Stokhos_Tpetra_Utilities_MP_Vector.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Details_WrappedDualView.hpp"
#include "Stokhos_Tpetra_CG.hpp"

// Belos solver
#ifdef HAVE_STOKHOS_BELOS
#include "Belos_Tpetra_MP_Vector.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#endif

// Ifpack2 preconditioner
#ifdef HAVE_STOKHOS_IFPACK2
#include "Stokhos_Ifpack2_MP_Vector.hpp"
#include "Ifpack2_Factory.hpp"
#endif

// MueLu preconditioner
#ifdef HAVE_STOKHOS_MUELU
#include "Stokhos_MueLu_MP_Vector.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#endif

// Amesos2 solver
#ifdef HAVE_STOKHOS_AMESOS2
#include "Stokhos_Amesos2_MP_Vector.hpp"
#include "Amesos2_Factory.hpp"
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

//
// Tests
//

const int VectorSize = STOKHOS_DEFAULT_ENSEMBLE_SIZE;

//
// Test vector addition
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, VectorAdd, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Comm
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  RCP<Tpetra_Vector> x1 = Tpetra::createVector<Scalar>(map);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  {
    auto x1_view = x1->getLocalViewHost(Tpetra::Access::OverwriteAll);
    auto x2_view = x2->getLocalViewHost(Tpetra::Access::OverwriteAll);
    Scalar val1(VectorSize, BaseScalar(0.0)), val2(VectorSize, BaseScalar(0.0));
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      for (LocalOrdinal j=0; j<VectorSize; ++j) {
        val1.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(nrow, VectorSize, row, j);
        val2.fastAccessCoeff(j) = 0.12345 * generate_vector_coefficient<BaseScalar,size_t>(nrow, VectorSize, row, j);
      }
      x1_view(i,0) = val1;
      x2_view(i,0) = val2;
    }
  }

  // Add
  Scalar alpha = 2.1;
  Scalar beta = 3.7;
  RCP<Tpetra_Vector> y = Tpetra::createVector<Scalar>(map);
  y->update(alpha, *x1, beta, *x2, Scalar(0.0));

  // y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check
  auto y_view = y->getLocalViewHost(Tpetra::Access::ReadOnly);
  Scalar val(VectorSize, BaseScalar(0.0));
  BaseScalar tol = 1.0e-14;
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      BaseScalar v = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, VectorSize, row, j);
      val.fastAccessCoeff(j) = alpha.coeff(j)*v + 0.12345*beta.coeff(j)*v;
    }
    TEST_EQUALITY( y_view(i,0).size(), VectorSize );
    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY( y_view(i,0).fastAccessCoeff(j), val.fastAccessCoeff(j), tol );
  }
}

//
// Test vector dot product
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, VectorDot, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef typename Tpetra_Vector::dot_type dot_type;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Comm
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  RCP<Tpetra_Vector> x1 = Tpetra::createVector<Scalar>(map);
  RCP<Tpetra_Vector> x2 = Tpetra::createVector<Scalar>(map);
  {
    auto x1_view = x1->getLocalViewHost(Tpetra::Access::OverwriteAll);
    auto x2_view = x2->getLocalViewHost(Tpetra::Access::OverwriteAll);
    Scalar val1(VectorSize, BaseScalar(0.0)), val2(VectorSize, BaseScalar(0.0));
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      for (LocalOrdinal j=0; j<VectorSize; ++j) {
        val1.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(nrow, VectorSize, row, j);
        val2.fastAccessCoeff(j) = 0.12345 * generate_vector_coefficient<BaseScalar,size_t>(nrow, VectorSize, row, j);
      }
      x1_view(i,0) = val1;
      x2_view(i,0) = val2;
    }
  }

  // Dot product
  dot_type dot = x1->dot(*x2);

  // Check

#ifdef HAVE_STOKHOS_ENSEMBLE_REDUCT

  // Local contribution
  dot_type local_val(0);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      BaseScalar v = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, VectorSize, row, j);
      local_val += 0.12345 * v * v;
    }
  }

  // Global reduction
  dot_type val(0);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, local_val,
                     Teuchos::outArg(val));

#else

  // Local contribution
  dot_type local_val(VectorSize, 0.0);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      BaseScalar v = generate_vector_coefficient<BaseScalar,size_t>(
        nrow, VectorSize, row, j);
      local_val.fastAccessCoeff(j) += 0.12345 * v * v;
    }
  }

  // Global reduction
  dot_type val(VectorSize, 0.0);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, local_val,
                     Teuchos::outArg(val));

#endif

  out << "dot = " << dot << " expected = " << val << std::endl;

  BaseScalar tol = 1.0e-14;
  TEST_FLOATING_EQUALITY( dot, val, tol );
}

//
// Test multi-vector addition
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, MultiVectorAdd, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Comm
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  size_t ncol = 5;
  RCP<Tpetra_MultiVector> x1 = Tpetra::createMultiVector<Scalar>(map, ncol);
  RCP<Tpetra_MultiVector> x2 = Tpetra::createMultiVector<Scalar>(map, ncol);
  {
    auto x1_view = x1->getLocalViewHost(Tpetra::Access::OverwriteAll);
    auto x2_view = x2->getLocalViewHost(Tpetra::Access::OverwriteAll);
    Scalar val1(VectorSize, BaseScalar(0.0)), val2(VectorSize, BaseScalar(0.0));
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      for (size_t j=0; j<ncol; ++j) {
        for (LocalOrdinal k=0; k<VectorSize; ++k) {
          BaseScalar v =
            generate_multi_vector_coefficient<BaseScalar,size_t>(
              nrow, ncol, VectorSize, row, j, k);
          val1.fastAccessCoeff(k) = v;
          val2.fastAccessCoeff(k) = 0.12345 * v;
        }
        x1_view(i,j) = val1;
        x2_view(i,j) = val2;
      }
    }
  }

  // Add
  Scalar alpha = 2.1;
  Scalar beta = 3.7;
  RCP<Tpetra_MultiVector> y = Tpetra::createMultiVector<Scalar>(map, ncol);
  y->update(alpha, *x1, beta, *x2, Scalar(0.0));

  // y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check
  auto y_view = y->getLocalViewHost(Tpetra::Access::ReadOnly);
  Scalar val(VectorSize, BaseScalar(0.0));
  BaseScalar tol = 1.0e-14;
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol; ++j) {
      for (LocalOrdinal k=0; k<VectorSize; ++k) {
        BaseScalar v = generate_multi_vector_coefficient<BaseScalar,size_t>(
          nrow, ncol, VectorSize, row, j, k);
        val.fastAccessCoeff(k) = alpha.coeff(k)*v + 0.12345*beta.coeff(k)*v;
      }
      TEST_EQUALITY( y_view(i,j).size(), VectorSize );
      for (LocalOrdinal k=0; k<VectorSize; ++k)
        TEST_FLOATING_EQUALITY( y_view(i,j).fastAccessCoeff(k),
                                val.fastAccessCoeff(k), tol );
    }
  }
}

//
// Test multi-vector dot product
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, MultiVectorDot, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector;
  typedef typename Tpetra_MultiVector::dot_type dot_type;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Comm
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  size_t ncol = 5;
  RCP<Tpetra_MultiVector> x1 = Tpetra::createMultiVector<Scalar>(map, ncol);
  RCP<Tpetra_MultiVector> x2 = Tpetra::createMultiVector<Scalar>(map, ncol);
  {
    auto x1_view = x1->getLocalViewHost(Tpetra::Access::OverwriteAll);
    auto x2_view = x2->getLocalViewHost(Tpetra::Access::OverwriteAll);
    Scalar val1(VectorSize, BaseScalar(0.0)), val2(VectorSize, BaseScalar(0.0));
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      for (size_t j=0; j<ncol; ++j) {
        for (LocalOrdinal k=0; k<VectorSize; ++k) {
          BaseScalar v =
            generate_multi_vector_coefficient<BaseScalar,size_t>(
              nrow, ncol, VectorSize, row, j, k);
          val1.fastAccessCoeff(k) = v;
          val2.fastAccessCoeff(k) = 0.12345 * v;
        }
        x1_view(i,j) = val1;
        x2_view(i,j) = val2;
      }
    }
  }

  // Dot product
  Array<dot_type> dots(ncol);
  x1->dot(*x2, dots());

  // Check

#ifdef HAVE_STOKHOS_ENSEMBLE_REDUCT

  // Local contribution
  Array<dot_type> local_vals(ncol, dot_type(0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol; ++j) {
      for (LocalOrdinal k=0; k<VectorSize; ++k) {
        BaseScalar v = generate_multi_vector_coefficient<BaseScalar,size_t>(
          nrow, ncol, VectorSize, row, j, k);
        local_vals[j] += 0.12345 * v * v;
      }
    }
  }

  // Global reduction
  Array<dot_type> vals(ncol, dot_type(0));
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<int>(ncol),
                     local_vals.getRawPtr(), vals.getRawPtr());

#else

  // Local contribution
  Array<dot_type> local_vals(ncol, dot_type(VectorSize, 0.0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol; ++j) {
      for (LocalOrdinal k=0; k<VectorSize; ++k) {
        BaseScalar v = generate_multi_vector_coefficient<BaseScalar,size_t>(
          nrow, ncol, VectorSize, row, j, k);
        local_vals[j].fastAccessCoeff(k) += 0.12345 * v * v;
      }
    }
  }

  // Global reduction
  Array<dot_type> vals(ncol, dot_type(VectorSize, 0.0));
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<int>(ncol),
                     local_vals.getRawPtr(), vals.getRawPtr());

#endif

  BaseScalar tol = 1.0e-14;
  for (size_t j=0; j<ncol; ++j) {
    out << "dots(" << j << ") = " << dots[j]
        << " expected(" << j << ") = " << vals[j] << std::endl;
    TEST_FLOATING_EQUALITY( dots[j], vals[j], tol );
  }
}

//
// Test multi-vector dot product using subviews
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, MultiVectorDotSub, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector;
  typedef typename Tpetra_MultiVector::dot_type dot_type;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Comm
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();

  // Map
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
  const size_t num_my_row = myGIDs.size();

  // Fill vectors
  size_t ncol = 5;
  RCP<Tpetra_MultiVector> x1 = Tpetra::createMultiVector<Scalar>(map, ncol);
  RCP<Tpetra_MultiVector> x2 = Tpetra::createMultiVector<Scalar>(map, ncol);
  {
    auto x1_view = x1->getLocalViewHost(Tpetra::Access::OverwriteAll);
    auto x2_view = x2->getLocalViewHost(Tpetra::Access::OverwriteAll);
    Scalar val1(VectorSize, BaseScalar(0.0)), val2(VectorSize, BaseScalar(0.0));
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      for (size_t j=0; j<ncol; ++j) {
        for (LocalOrdinal k=0; k<VectorSize; ++k) {
          BaseScalar v =
            generate_multi_vector_coefficient<BaseScalar,size_t>(
              nrow, ncol, VectorSize, row, j, k);
          val1.fastAccessCoeff(k) = v;
          val2.fastAccessCoeff(k) = 0.12345 * v;
        }
        x1_view(i,j) = val1;
        x2_view(i,j) = val2;
      }
    }
  }

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

#ifdef HAVE_STOKHOS_ENSEMBLE_REDUCT

  // Local contribution
  Array<dot_type> local_vals(ncol_sub, dot_type(0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol_sub; ++j) {
      for (LocalOrdinal k=0; k<VectorSize; ++k) {
        BaseScalar v = generate_multi_vector_coefficient<BaseScalar,size_t>(
          nrow, ncol, VectorSize, row, cols[j], k);
        local_vals[j] += 0.12345 * v * v;
      }
    }
  }

  // Global reduction
  Array<dot_type> vals(ncol_sub, dot_type(0));
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM,
                     Teuchos::as<int>(ncol_sub), local_vals.getRawPtr(),
                     vals.getRawPtr());

#else

  // Local contribution
  Array<dot_type> local_vals(ncol_sub, dot_type(VectorSize, 0.0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol_sub; ++j) {
      for (LocalOrdinal k=0; k<VectorSize; ++k) {
        BaseScalar v = generate_multi_vector_coefficient<BaseScalar,size_t>(
          nrow, ncol, VectorSize, row, cols[j], k);
        local_vals[j].fastAccessCoeff(k) += 0.12345 * v * v;
      }
    }
  }

  // Global reduction
  Array<dot_type> vals(ncol_sub, dot_type(VectorSize, 0.0));
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM,
                     Teuchos::as<int>(ncol_sub), local_vals.getRawPtr(),
                     vals.getRawPtr());

#endif

  BaseScalar tol = 1.0e-14;
  for (size_t j=0; j<ncol_sub; ++j) {
    out << "dots(" << j << ") = " << dots[j]
        << " expected(" << j << ") = " << vals[j] << std::endl;
    TEST_FLOATING_EQUALITY( dots[j], vals[j], tol );
  }
}

//
// Test matrix-vector multiplication for a simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, MatrixVectorMultiply, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2)));
  Array<GlobalOrdinal> columnIndices(2);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
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
  {
    auto x_view = x->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      for (LocalOrdinal j=0; j<VectorSize; ++j)
        val.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(
          nrow, VectorSize, row, j);
      x_view(i,0) = val;
    }
  }

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
  auto y_view = y->getLocalViewHost(Tpetra::Access::ReadOnly);
  BaseScalar tol = 1.0e-14;
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
    TEST_EQUALITY( y_view(i,0).size(), VectorSize );
    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY( y_view(i,0).fastAccessCoeff(j), val.fastAccessCoeff(j), tol );
  }
}

//
// Test matrix-multi-vector multiplication for a simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, MatrixMultiVectorMultiply, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2)));
  Array<GlobalOrdinal> columnIndices(2);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
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

  // Fill multi-vector
  size_t ncol = 5;
  RCP<Tpetra_MultiVector> x = Tpetra::createMultiVector<Scalar>(map, ncol);
  {
    auto x_view = x->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      for (size_t j=0; j<ncol; ++j) {
        for (LocalOrdinal k=0; k<VectorSize; ++k) {
          BaseScalar v =
            generate_multi_vector_coefficient<BaseScalar,size_t>(
              nrow, ncol, VectorSize, row, j, k);
          val.fastAccessCoeff(k) = v;
        }
        x_view(i,j) = val;
      }
    }
  }

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
  auto y_view = y->getLocalViewHost(Tpetra::Access::ReadOnly);
  BaseScalar tol = 1.0e-14;
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    for (size_t j=0; j<ncol; ++j) {
      for (LocalOrdinal k=0; k<VectorSize; ++k) {
        BaseScalar v1 = generate_vector_coefficient<BaseScalar,size_t>(
          nrow, VectorSize, row, k);
        BaseScalar v2 = generate_multi_vector_coefficient<BaseScalar,size_t>(
          nrow, ncol, VectorSize, row, j, k);
        val.fastAccessCoeff(k) = v1*v2;
      }
      if (row != nrow-1) {
        for (LocalOrdinal k=0; k<VectorSize; ++k) {
          BaseScalar v1 = generate_vector_coefficient<BaseScalar,size_t>(
            nrow, VectorSize, row+1, k);
          BaseScalar v2 = generate_multi_vector_coefficient<BaseScalar,size_t>(
            nrow, ncol, VectorSize, row+1, j, k);
          val.fastAccessCoeff(k) += v1*v2;
        }
      }
      TEST_EQUALITY( y_view(i,j).size(), VectorSize );
      for (LocalOrdinal k=0; k<VectorSize; ++k)
        TEST_FLOATING_EQUALITY( y_view(i,j).fastAccessCoeff(k),
                                val.fastAccessCoeff(k), tol );
    }
  }
}

//
// Test flattening MP::Vector matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, Flatten, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

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

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2)));
  Array<GlobalOrdinal> columnIndices(2);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
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
  {
    auto x_view = x->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      for (LocalOrdinal j=0; j<VectorSize; ++j)
        val.fastAccessCoeff(j) = generate_vector_coefficient<BaseScalar,size_t>(
          nrow, VectorSize, row, j);
        x_view(i,0) = val;
    }
  }

  // Multiply
  RCP<Tpetra_Vector> y = Tpetra::createVector<Scalar>(map);
  matrix->apply(*x, *y);

  // Flatten matrix
  RCP<const Tpetra_Map> flat_x_map, flat_y_map;
  RCP<const Tpetra_CrsGraph> flat_graph =
    Stokhos::create_flat_mp_graph(*graph, flat_x_map, flat_y_map, VectorSize);
  RCP<Flat_Tpetra_CrsMatrix> flat_matrix =
    Stokhos::create_flat_matrix(*matrix, flat_graph, VectorSize);

  // Multiply with flattened matix
  RCP<Tpetra_Vector> y2 = Tpetra::createVector<Scalar>(map);
  {
    RCP<Flat_Tpetra_Vector> flat_x =
      Stokhos::create_flat_vector_view(*x, flat_x_map);
    RCP<Flat_Tpetra_Vector> flat_y =
      Stokhos::create_flat_vector_view(*y2, flat_y_map);
    flat_matrix->apply(*flat_x, *flat_y);
  }

  // flat_y->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //                  Teuchos::VERB_EXTREME);

  // Check
  BaseScalar tol = 1.0e-14;
  auto y_view  = y-> getLocalViewHost(Tpetra::Access::ReadOnly);
  auto y2_view = y2->getLocalViewHost(Tpetra::Access::ReadOnly);
  for (size_t i=0; i<num_my_row; ++i) {
    TEST_EQUALITY( y_view( i,0).size(), VectorSize );
    TEST_EQUALITY( y2_view(i,0).size(), VectorSize );
    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY( y_view( i,0).fastAccessCoeff(j),
                              y2_view(i,0).fastAccessCoeff(j), tol );
  }
}

//
// Test interaction between Tpetra WrappedDualView and MP::Vector
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, WrappedDualView, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  //BMK 6-2021: This test is required because a View of MP::Vector has slightly different behavior than a typical Kokkos::View.
  //If you construct a Kokkos::View with a label and 0 extent, it gets a non-null allocation.
  //But for View<MP::Vector>, the same constructor produces a null data pointer but
  //an active reference counting node (use_count() > 0).
  //This test makes sure that Tpetra WrappedDualView works correctly with a View where data() == nullptr but use_count() > 0.
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  //typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  using DualViewType = Kokkos::DualView<Scalar*, typename Node::device_type>;
  using WDV = Tpetra::Details::WrappedDualView<DualViewType>;
  using values_view = typename DualViewType::t_dev;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  WDV wdv;
  {
    values_view myView("emptyTestView", 0);
    wdv = WDV(myView);
  }
  size_t use_h = wdv.getHostView(Tpetra::Access::ReadOnly).use_count();
  size_t use_d = wdv.getDeviceView(Tpetra::Access::ReadOnly).use_count();
  //The WrappedDualView is now the only object holding references to the host and device views,
  //so they should have identical use counts.
  TEST_EQUALITY(use_h, use_d);
}

//
// Test simple CG solve without preconditioning for a 1-D Laplacian matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, SimpleCG, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // 1-D Laplacian matrix
  GlobalOrdinal nrow = 50;
  BaseScalar h = 1.0 / static_cast<BaseScalar>(nrow-1);
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(3)));
  Array<GlobalOrdinal> columnIndices(3);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
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
      vals[0] = Scalar(-1.0) * a_val;
      vals[1] = Scalar(2.0) * a_val;
      vals[2] = Scalar(-1.0) * a_val;
      matrix->replaceGlobalValues(row, columnIndices(0,3), vals(0,3));
    }
  }
  matrix->fillComplete();

  matrix->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
                   Teuchos::VERB_EXTREME);

  // Fill RHS vector
  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  Scalar b_val;
  {
    auto b_view = b->getLocalViewHost(Tpetra::Access::OverwriteAll);
    b_val = Scalar(VectorSize, BaseScalar(0.0));
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      b_val.fastAccessCoeff(j) =
        BaseScalar(-1.0) + BaseScalar(j) / BaseScalar(VectorSize);
    }
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      if (row == 0 || row == nrow-1)
        b_view(i,0) = Scalar(0.0);
      else
        b_view(i,0) = -Scalar(b_val * h * h);
    }
  }

  // Solve
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  typedef Kokkos::ArithTraits<BaseScalar> BST;
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

  // Check -- For a*y'' = b, correct answer is y = 0.5 *(b/a) * x * (x-1)
  btol = 1000*btol;
  auto x_view = x->getLocalViewHost(Tpetra::Access::ReadOnly);
  Scalar val(VectorSize, BaseScalar(0.0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    BaseScalar xx = row * h;
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      val.fastAccessCoeff(j) =
        BaseScalar(0.5) * (b_val.coeff(j)/a_val.coeff(j)) * xx * (xx - BaseScalar(1.0));
    }
    TEST_EQUALITY( x_view(i,0).size(), VectorSize );

    // Set small values to zero
    Scalar v = x_view(i,0);
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      if (BST::abs(v.coeff(j)) < btol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
      if (BST::abs(val.coeff(j)) < btol)
        val.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), val.coeff(j), btol);
  }

}

#if defined(HAVE_STOKHOS_MUELU) && defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_STOKHOS_IFPACK2)

//
// Test simple CG solve with MueLu preconditioning for a 1-D Laplacian matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, SimplePCG_Muelu, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::getParametersFromXmlFile;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // 1-D Laplacian matrix
  GlobalOrdinal nrow = 50;
  BaseScalar h = 1.0 / static_cast<BaseScalar>(nrow-1);
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(3)));
  Array<GlobalOrdinal> columnIndices(3);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
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
      vals[0] = Scalar(-1.0) * a_val;
      vals[1] = Scalar(2.0) * a_val;
      vals[2] = Scalar(-1.0) * a_val;
      matrix->replaceGlobalValues(row, columnIndices(0,3), vals(0,3));
    }
  }
  matrix->fillComplete();

  // Fill RHS vector
  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  Scalar b_val;
  {
    auto b_view = b->getLocalViewHost(Tpetra::Access::OverwriteAll);
    b_val = Scalar(VectorSize, BaseScalar(0.0));
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      b_val.fastAccessCoeff(j) =
        BaseScalar(-1.0) + BaseScalar(j) / BaseScalar(VectorSize);
    }
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      if (row == 0 || row == nrow-1)
        b_view(i,0) = Scalar(0.0);
      else
        b_view(i,0) = -Scalar(b_val * h * h);
    }
  }

  // Create preconditioner
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
  RCP<OP> matrix_op = matrix;
  RCP<ParameterList> muelu_params =
    getParametersFromXmlFile("muelu_cheby.xml");
  RCP<OP> M =
    MueLu::CreateTpetraPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>(matrix_op, *muelu_params);

  // Solve
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map);
  typedef Kokkos::ArithTraits<BaseScalar> BST;
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

  // Check -- For a*y'' = b, correct answer is y = 0.5 *(b/a) * x * (x-1)
  btol = 1000*btol;
  auto x_view = x->getLocalViewHost(Tpetra::Access::ReadOnly);
  Scalar val(VectorSize, BaseScalar(0.0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    BaseScalar xx = row * h;
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      val.fastAccessCoeff(j) =
        BaseScalar(0.5) * (b_val.coeff(j)/a_val.coeff(j)) * xx * (xx - BaseScalar(1.0));
    }
    TEST_EQUALITY( x_view(i,0).size(), VectorSize );

    // Set small values to zero
    Scalar v = x_view(i,0);
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      if (BST::magnitude(v.coeff(j)) < btol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
      if (BST::magnitude(val.coeff(j)) < btol)
        val.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), val.coeff(j), btol);
  }

}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, SimplePCG_Muelu, Storage, LocalOrdinal, GlobalOrdinal, Node ) {}

#endif

#if defined(HAVE_STOKHOS_BELOS)

//
// Test Belos GMRES solve for a simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2)));
  Array<GlobalOrdinal> columnIndices(2);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
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
  {
    auto b_view = b->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (size_t i=0; i<num_my_row; ++i) {
      b_view(i,0) = Scalar(1.0);
    }
  }

  // Solve
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
#ifdef HAVE_STOKHOS_ENSEMBLE_REDUCT
  typedef BaseScalar BelosScalar;
#else
  typedef Scalar BelosScalar;
#endif
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

  // Check -- Correct answer is:
  //     [ 0, 0,   ..., 0            ]
  //     [ 1, 1/2, ..., 1/VectorSize ]
  //     [ 0, 0,   ..., 0            ]
  //     [ 1, 1/2, ..., 1/VectorSize ]
  //     ....
  tol = 1000*tol;
  auto x_view = x->getLocalViewHost(Tpetra::Access::ReadOnly);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row % 2) {
      for (LocalOrdinal j=0; j<VectorSize; ++j) {
        val.fastAccessCoeff(j) = BaseScalar(1.0) / BaseScalar(j+1);
      }
    }
    else
      val = Scalar(VectorSize, BaseScalar(0.0));
    TEST_EQUALITY( x_view(i,0).size(), VectorSize );

    // Set small values to zero
    Scalar v = x_view(i,0);
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      if (ST::magnitude(v.coeff(j)) < tol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), val.coeff(j), tol);
  }
}

//
// Test Belos GMRES solve for a simple lower-triangular matrix with lucky breakdown with DGKS orthogonalization
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES_DGKS, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Build diagonal matrix
  GlobalOrdinal nrow = 20;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(1)));
  Array<GlobalOrdinal> columnIndices(1);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
  const size_t num_my_row = myGIDs.size();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    size_t ncol = 1;
    graph->insertGlobalIndices(row, columnIndices(0,ncol));
  }
  graph->fillComplete();
  RCP<Tpetra_CrsMatrix> matrix = rcp(new Tpetra_CrsMatrix(graph));

  Array<Scalar> vals(1);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    vals[0] = Scalar(row+1);
    matrix->replaceGlobalValues(row, columnIndices(0,1), vals(0,1));
  }
  matrix->fillComplete();

  // Fill RHS vector:
  //     [ 0, 0, ..., 0, 0, 1]
  //     [ 0, 0, ..., 0, 2, 2]
  //     [ 0, 0, ..., 3, 3, 3]
  //     [ 0, 0, ..., 4, 4, 4]
  //     ...

  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  {
    auto b_view = b->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      b_view(i,0) = Scalar(0.0);
      for (LocalOrdinal j=0; j<VectorSize; ++j)
        if (int(j+2+row-VectorSize) > 0)
          b_view(i,0).fastAccessCoeff(j) = BaseScalar(row+1);
    }
  }

  // Solve
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
#ifdef HAVE_STOKHOS_ENSEMBLE_REDUCT
  typedef BaseScalar BelosScalar;
#else
  typedef Scalar BelosScalar;
#endif
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
  belosParams->set("Orthogonalization","DGKS");
  RCP<Belos::PseudoBlockGmresSolMgr<BelosScalar,MV,OP> > solver =
    rcp(new Belos::PseudoBlockGmresSolMgr<BelosScalar,MV,OP>(problem, belosParams));
  problem->setProblem();
  Belos::ReturnType ret = solver->solve();
  TEST_EQUALITY_CONST( ret, Belos::Converged );

#ifndef HAVE_STOKHOS_ENSEMBLE_REDUCT
  int numItersExpected = nrow;
  int numIters = solver->getNumIters();
  out << "numIters = " << numIters << std::endl;
  TEST_EQUALITY( numIters, numItersExpected);

  // Get and print number of ensemble iterations
  std::vector<int> ensemble_iterations =
    static_cast<const Belos::StatusTestImpResNorm<BelosScalar, MV, OP> *>(solver->getResidualStatusTest())->getEnsembleIterations();
  out << "Ensemble iterations = ";
  for (auto ensemble_iteration : ensemble_iterations)
    out << ensemble_iteration << " ";
  out << std::endl;

  for (LocalOrdinal j=0; j<VectorSize; ++j) {
    if (int(j+1+nrow-VectorSize) > 0) {
      TEST_EQUALITY(int(j+1+nrow-VectorSize), ensemble_iterations[j]);
    }
    else {
      TEST_EQUALITY(int(0), ensemble_iterations[j]);
    }
  }
#endif

  // Check -- Correct answer is:
  //     [ 0, 0, ..., 0, 0, 1]
  //     [ 0, 0, ..., 0, 1, 1]
  //     [ 0, 0, ..., 1, 1, 1]
  //     [ 0, 0, ..., 1, 1, 1]
  //     ...
  tol = 1000*tol;
  auto x_view = x->getLocalViewHost(Tpetra::Access::ReadOnly);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    Scalar v = x_view(i,0);

    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      if (ST::magnitude(v.coeff(j)) < tol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    Scalar val = Scalar(0.0);

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      if (j+2+row-VectorSize > 0)
        val.fastAccessCoeff(j) = BaseScalar(1.0);

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), val.coeff(j), tol);
  }
}

//
// Test Belos GMRES solve for a simple lower-triangular matrix with lucky breakdown with ICGS orthogonalization
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES_ICGS, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Build diagonal matrix
  GlobalOrdinal nrow = 20;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(1)));
  Array<GlobalOrdinal> columnIndices(1);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
  const size_t num_my_row = myGIDs.size();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    size_t ncol = 1;
    graph->insertGlobalIndices(row, columnIndices(0,ncol));
  }
  graph->fillComplete();
  RCP<Tpetra_CrsMatrix> matrix = rcp(new Tpetra_CrsMatrix(graph));

  Array<Scalar> vals(1);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    vals[0] = Scalar(row+1);
    matrix->replaceGlobalValues(row, columnIndices(0,1), vals(0,1));
  }
  matrix->fillComplete();

  // Fill RHS vector:
  //     [ 0, 0, ..., 0, 0, 1]
  //     [ 0, 0, ..., 0, 2, 2]
  //     [ 0, 0, ..., 3, 3, 3]
  //     [ 0, 0, ..., 4, 4, 4]
  //     ...

  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  {
    auto b_view = b->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      b_view(i,0) = Scalar(0.0);
      for (LocalOrdinal j=0; j<VectorSize; ++j)
        if (int(j+2+row-VectorSize) > 0)
          b_view(i,0).fastAccessCoeff(j) = BaseScalar(row+1);
    }
  }

  // Solve
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
#ifdef HAVE_STOKHOS_ENSEMBLE_REDUCT
  typedef BaseScalar BelosScalar;
#else
  typedef Scalar BelosScalar;
#endif
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
  belosParams->set("Orthogonalization","ICGS");
  RCP<Belos::PseudoBlockGmresSolMgr<BelosScalar,MV,OP> > solver =
    rcp(new Belos::PseudoBlockGmresSolMgr<BelosScalar,MV,OP>(problem, belosParams));
  problem->setProblem();
  Belos::ReturnType ret = solver->solve();
  TEST_EQUALITY_CONST( ret, Belos::Converged );

#ifndef HAVE_STOKHOS_ENSEMBLE_REDUCT
  int numItersExpected = nrow;
  int numIters = solver->getNumIters();
  out << "numIters = " << numIters << std::endl;
  TEST_EQUALITY( numIters, numItersExpected);

  // Get and print number of ensemble iterations
  std::vector<int> ensemble_iterations =
    static_cast<const Belos::StatusTestImpResNorm<BelosScalar, MV, OP> *>(solver->getResidualStatusTest())->getEnsembleIterations();
  out << "Ensemble iterations = ";
  for (auto ensemble_iteration : ensemble_iterations)
    out << ensemble_iteration << " ";
  out << std::endl;

  for (LocalOrdinal j=0; j<VectorSize; ++j) {
    if (int(j+1+nrow-VectorSize) > 0) {
      TEST_EQUALITY(int(j+1+nrow-VectorSize), ensemble_iterations[j]);
    }
    else {
      TEST_EQUALITY(int(0), ensemble_iterations[j]);
    }
  }
#endif

  // Check -- Correct answer is:
  //     [ 0, 0, ..., 0, 0, 1]
  //     [ 0, 0, ..., 0, 1, 1]
  //     [ 0, 0, ..., 1, 1, 1]
  //     [ 0, 0, ..., 1, 1, 1]
  //     ...
  tol = 1000*tol;
  auto x_view = x->getLocalViewHost(Tpetra::Access::ReadOnly);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    Scalar v = x_view(i,0);

    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      if (ST::magnitude(v.coeff(j)) < tol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    Scalar val = Scalar(0.0);

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      if (j+2+row-VectorSize > 0)
        val.fastAccessCoeff(j) = BaseScalar(1.0);

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), val.coeff(j), tol);
  }
}

//
// Test Belos GMRES solve for a simple lower-triangular matrix with lucky breakdown with IMGS orthogonalization
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES_IMGS, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Build diagonal matrix
  GlobalOrdinal nrow = 20;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(1)));
  Array<GlobalOrdinal> columnIndices(1);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
  const size_t num_my_row = myGIDs.size();
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    size_t ncol = 1;
    graph->insertGlobalIndices(row, columnIndices(0,ncol));
  }
  graph->fillComplete();
  RCP<Tpetra_CrsMatrix> matrix = rcp(new Tpetra_CrsMatrix(graph));

  Array<Scalar> vals(1);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    columnIndices[0] = row;
    vals[0] = Scalar(row+1);
    matrix->replaceGlobalValues(row, columnIndices(0,1), vals(0,1));
  }
  matrix->fillComplete();

  // Fill RHS vector:
  //     [ 0, 0, ..., 0, 0, 1]
  //     [ 0, 0, ..., 0, 2, 2]
  //     [ 0, 0, ..., 3, 3, 3]
  //     [ 0, 0, ..., 4, 4, 4]
  //     ...

  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  {
    auto b_view = b->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      b_view(i,0) = Scalar(0.0);
      for (LocalOrdinal j=0; j<VectorSize; ++j)
        if (int(j+2+row-VectorSize) > 0)
          b_view(i,0).fastAccessCoeff(j) = BaseScalar(row+1);
    }
  }

  // Solve
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
#ifdef HAVE_STOKHOS_ENSEMBLE_REDUCT
  typedef BaseScalar BelosScalar;
#else
  typedef Scalar BelosScalar;
#endif
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
  belosParams->set("Orthogonalization","IMGS");
  RCP<Belos::PseudoBlockGmresSolMgr<BelosScalar,MV,OP> > solver =
    rcp(new Belos::PseudoBlockGmresSolMgr<BelosScalar,MV,OP>(problem, belosParams));
  problem->setProblem();
  Belos::ReturnType ret = solver->solve();
  TEST_EQUALITY_CONST( ret, Belos::Converged );

#ifndef HAVE_STOKHOS_ENSEMBLE_REDUCT
  int numItersExpected = nrow;
  int numIters = solver->getNumIters();
  out << "numIters = " << numIters << std::endl;
  TEST_EQUALITY( numIters, numItersExpected);

  // Get and print number of ensemble iterations
  std::vector<int> ensemble_iterations =
    static_cast<const Belos::StatusTestImpResNorm<BelosScalar, MV, OP> *>(solver->getResidualStatusTest())->getEnsembleIterations();
  out << "Ensemble iterations = ";
  for (auto ensemble_iteration : ensemble_iterations)
    out << ensemble_iteration << " ";
  out << std::endl;

  for (LocalOrdinal j=0; j<VectorSize; ++j) {
    if (int(j+1+nrow-VectorSize) > 0) {
      TEST_EQUALITY(int(j+1+nrow-VectorSize), ensemble_iterations[j]);
    }
    else {
      TEST_EQUALITY(int(0), ensemble_iterations[j]);
    }
  }
#endif

  // Check -- Correct answer is:
  //     [ 0, 0, ..., 0, 0, 1]
  //     [ 0, 0, ..., 0, 1, 1]
  //     [ 0, 0, ..., 1, 1, 1]
  //     [ 0, 0, ..., 1, 1, 1]
  //     ...
  tol = 1000*tol;
  auto x_view = x->getLocalViewHost(Tpetra::Access::ReadOnly);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    Scalar v = x_view(i,0);

    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      if (ST::magnitude(v.coeff(j)) < tol)
        v.fastAccessCoeff(j) = BaseScalar(0.0);
    }

    Scalar val = Scalar(0.0);

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      if (j+2+row-VectorSize > 0)
        val.fastAccessCoeff(j) = BaseScalar(1.0);

    for (LocalOrdinal j=0; j<VectorSize; ++j)
      TEST_FLOATING_EQUALITY(v.coeff(j), val.coeff(j), tol);
  }
}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES, Storage, LocalOrdinal, GlobalOrdinal, Node )
{}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES_DGKS, Storage, LocalOrdinal, GlobalOrdinal, Node )
{}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES_ICGS, Storage, LocalOrdinal, GlobalOrdinal, Node )
{}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES_IMGS, Storage, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

#if defined(HAVE_STOKHOS_BELOS) && defined(HAVE_STOKHOS_IFPACK2)

//
// Test Belos GMRES solve with Ifpack2 RILUK preconditioning for a
// simple banded upper-triangular matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosGMRES_RILUK, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // Build banded matrix
  GlobalOrdinal nrow = 10;
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(2)));
  Array<GlobalOrdinal> columnIndices(2);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
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
  {
    auto b_view = b->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (size_t i=0; i<num_my_row; ++i) {
      b_view(i,0) = Scalar(1.0);
    }
  }

  // Create preconditioner
  typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Prec;
  Ifpack2::Factory factory;
  RCP<Prec> M = factory.create<Tpetra_CrsMatrix>("RILUK", matrix);
  M->initialize();
  M->compute();

  // Solve
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
#ifdef HAVE_STOKHOS_ENSEMBLE_REDUCT
  typedef BaseScalar BelosScalar;
#else
  typedef Scalar BelosScalar;
#endif
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

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check -- Correct answer is:
  //     [ 0, 0,   ..., 0            ]
  //     [ 1, 1/2, ..., 1/VectorSize ]
  //     [ 0, 0,   ..., 0            ]
  //     [ 1, 1/2, ..., 1/VectorSize ]
  //     ....
  tol = 1000*tol;
  auto x_view = x->getLocalViewHost(Tpetra::Access::ReadOnly);
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    if (row % 2) {
      for (LocalOrdinal j=0; j<VectorSize; ++j) {
        val.fastAccessCoeff(j) = BaseScalar(1.0) / BaseScalar(j+1);
      }
    }
    else
      val = Scalar(VectorSize, BaseScalar(0.0));
    TEST_EQUALITY( x_view(i,0).size(), VectorSize );

    // Set small values to zero
    Scalar v = x_view(i,0);
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
  Tpetra_CrsMatrix_MP, BelosGMRES_RILUK, Storage, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

#if defined(HAVE_STOKHOS_BELOS) && defined(HAVE_STOKHOS_IFPACK2) && defined(HAVE_STOKHOS_MUELU) && defined(HAVE_STOKHOS_AMESOS2)

//
// Test Belos CG solve with MueLu preconditioning for a 1-D Laplacian matrix
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(
  Tpetra_CrsMatrix_MP, BelosCG_Muelu, Storage, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::getParametersFromXmlFile;

  typedef typename Storage::value_type BaseScalar;
  typedef Sacado::MP::Vector<Storage> Scalar;

  typedef Teuchos::Comm<int> Tpetra_Comm;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;

  // Ensure device is initialized
  if ( !Kokkos::is_initialized() )
    Kokkos::initialize();

  // 1-D Laplacian matrix
  GlobalOrdinal nrow = 50;
  BaseScalar h = 1.0 / static_cast<BaseScalar>(nrow-1);
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm();
  RCP<const Tpetra_Map> map =
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
      nrow, comm);
  RCP<Tpetra_CrsGraph> graph =
    rcp(new Tpetra_CrsGraph(map, size_t(3)));
  Array<GlobalOrdinal> columnIndices(3);
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList();
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
      vals[0] = Scalar(-1.0) * a_val;
      vals[1] = Scalar(2.0) * a_val;
      vals[2] = Scalar(-1.0) * a_val;
      matrix->replaceGlobalValues(row, columnIndices(0,3), vals(0,3));
    }
  }
  matrix->fillComplete();

  // Fill RHS vector
  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map);
  Scalar b_val;
  {
    auto b_view = b->getLocalViewHost(Tpetra::Access::OverwriteAll);
    b_val = Scalar(VectorSize, BaseScalar(0.0));
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      b_val.fastAccessCoeff(j) =
        BaseScalar(-1.0) + BaseScalar(j) / BaseScalar(VectorSize);
    }
    for (size_t i=0; i<num_my_row; ++i) {
      const GlobalOrdinal row = myGIDs[i];
      if (row == 0 || row == nrow-1)
        b_view(i,0) = Scalar(0.0);
      else
        b_view(i,0) = -Scalar(b_val * h * h);
    }
  }

  // Create preconditioner
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
  RCP<ParameterList> muelu_params =
    getParametersFromXmlFile("muelu_cheby.xml");
  RCP<OP> matrix_op = matrix;
  RCP<OP> M =
    MueLu::CreateTpetraPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>(matrix_op, *muelu_params);

  // Solve
  typedef Teuchos::ScalarTraits<BaseScalar> ST;
#ifdef HAVE_STOKHOS_ENSEMBLE_REDUCT
  typedef BaseScalar BelosScalar;
#else
  typedef Scalar BelosScalar;
#endif
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
  // Turn off residual scaling so we can see some variation in the number
  // of iterations across the ensemble when not doing ensemble reductions
  belosParams->set("Implicit Residual Scaling", "None");

  RCP<Belos::PseudoBlockCGSolMgr<BelosScalar,MV,OP,true> > solver =
    rcp(new Belos::PseudoBlockCGSolMgr<BelosScalar,MV,OP,true>(problem, belosParams));
  Belos::ReturnType ret = solver->solve();
  TEST_EQUALITY_CONST( ret, Belos::Converged );

#ifndef HAVE_STOKHOS_ENSEMBLE_REDUCT
  // Get and print number of ensemble iterations
  std::vector<int> ensemble_iterations =
    solver->getResidualStatusTest()->getEnsembleIterations();
  out << "Ensemble iterations = ";
  for (int i=0; i<VectorSize; ++i)
    out << ensemble_iterations[i] << " ";
  out << std::endl;
#endif

  // x->describe(*(Teuchos::fancyOStream(rcp(&std::cout,false))),
  //             Teuchos::VERB_EXTREME);

  // Check -- For a*y'' = b, correct answer is y = 0.5 *(b/a) * x * (x-1)
  tol = 1000*tol;
  auto x_view = x->getLocalViewHost(Tpetra::Access::ReadOnly);
  Scalar val(VectorSize, BaseScalar(0.0));
  for (size_t i=0; i<num_my_row; ++i) {
    const GlobalOrdinal row = myGIDs[i];
    BaseScalar xx = row * h;
    for (LocalOrdinal j=0; j<VectorSize; ++j) {
      val.fastAccessCoeff(j) =
        BaseScalar(0.5) * (b_val.coeff(j)/a_val.coeff(j)) * xx * (xx - BaseScalar(1.0));
    }
    TEST_EQUALITY( x_view(i,0).size(), VectorSize );

    // Set small values to zero
    Scalar v = x_view(i,0);
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
  Tpetra_CrsMatrix_MP, BelosCG_Muelu, Storage, LocalOrdinal, GlobalOrdinal, Node )
{}

#endif

#if defined(HAVE_STOKHOS_AMESOS2)

//
// Test Amesos2 solve for a 1-D Laplacian matrix
//
#define TEST_AMESOS2_SOLVER(SolverName) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_CrsMatrix_MP, Amesos2_##SolverName, Storage, LocalOrdinal, GlobalOrdinal, Node) \
{ \
  using Teuchos::RCP; \
  using Teuchos::rcp; \
  using Teuchos::ArrayView; \
  using Teuchos::Array; \
  using Teuchos::ArrayRCP; \
  using Teuchos::ParameterList; \
  \
  typedef typename Storage::value_type BaseScalar; \
  typedef Sacado::MP::Vector<Storage> Scalar; \
  \
  typedef Teuchos::Comm<int> Tpetra_Comm; \
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map; \
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector; \
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_MultiVector; \
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix; \
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph; \
  \
  /* Ensure device is initialized */ \
  if ( !Kokkos::is_initialized() ) \
    Kokkos::initialize(); \
  \
  /* 1-D Laplacian matrix */ \
  GlobalOrdinal nrow = 50; \
  BaseScalar h = 1.0 / static_cast<BaseScalar>(nrow-1); \
  RCP<const Tpetra_Comm> comm = Tpetra::getDefaultComm(); \
  RCP<const Tpetra_Map> map = \
    Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>( \
      nrow, comm); \
  RCP<Tpetra_CrsGraph> graph = Tpetra::createCrsGraph(map, size_t(3)); \
  Array<GlobalOrdinal> columnIndices(3); \
  ArrayView<const GlobalOrdinal> myGIDs = map->getLocalElementList(); \
  const size_t num_my_row = myGIDs.size(); \
  for (size_t i=0; i<num_my_row; ++i) { \
    const GlobalOrdinal row = myGIDs[i]; \
    if (row == 0 || row == nrow-1) { /* Boundary nodes */ \
      columnIndices[0] = row; \
      graph->insertGlobalIndices(row, columnIndices(0,1)); \
    } \
    else { /* Interior nodes */ \
      columnIndices[0] = row-1; \
      columnIndices[1] = row; \
      columnIndices[2] = row+1; \
      graph->insertGlobalIndices(row, columnIndices(0,3)); \
    } \
  } \
  graph->fillComplete(); \
  RCP<Tpetra_CrsMatrix> matrix = rcp(new Tpetra_CrsMatrix(graph)); \
  \
  /* Set values in matrix */ \
  Array<Scalar> vals(3); \
  Scalar a_val(VectorSize, BaseScalar(0.0)); \
  for (LocalOrdinal j=0; j<VectorSize; ++j) { \
    a_val.fastAccessCoeff(j) = \
      BaseScalar(1.0) + BaseScalar(j) / BaseScalar(VectorSize); \
  } \
  for (size_t i=0; i<num_my_row; ++i) { \
    const GlobalOrdinal row = myGIDs[i]; \
    if (row == 0 || row == nrow-1) { /* Boundary nodes */ \
      columnIndices[0] = row; \
      vals[0] = Scalar(1.0); \
      matrix->replaceGlobalValues(row, columnIndices(0,1), vals(0,1)); \
    } \
    else { \
      columnIndices[0] = row-1; \
      columnIndices[1] = row; \
      columnIndices[2] = row+1; \
      vals[0] = Scalar(-1.0) * a_val; \
      vals[1] = Scalar(2.0) * a_val; \
      vals[2] = Scalar(-1.0) * a_val; \
      matrix->replaceGlobalValues(row, columnIndices(0,3), vals(0,3)); \
    } \
  } \
  matrix->fillComplete();\
  \
  /* Fill RHS vector */ \
  RCP<Tpetra_Vector> b = Tpetra::createVector<Scalar>(map); \
  Scalar b_val; \
  { \
    auto b_view = b->getLocalViewHost(Tpetra::Access::OverwriteAll); \
    b_val = Scalar(VectorSize, BaseScalar(0.0)); \
    for (LocalOrdinal j=0; j<VectorSize; ++j) { \
      b_val.fastAccessCoeff(j) = \
        BaseScalar(-1.0) + BaseScalar(j) / BaseScalar(VectorSize); \
    } \
    for (size_t i=0; i<num_my_row; ++i) { \
      const GlobalOrdinal row = myGIDs[i]; \
      if (row == 0 || row == nrow-1) \
        b_view(i,0) = Scalar(0.0); \
      else \
        b_view(i,0) = -Scalar(b_val * h * h); \
    } \
  } \
  \
  /* Solve */ \
  typedef Amesos2::Solver<Tpetra_CrsMatrix,Tpetra_MultiVector> Solver; \
  RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(map); \
  std::string solver_name(#SolverName); \
  out << "Solving linear system with " << solver_name << std::endl; \
  RCP<Solver> solver = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>( \
    solver_name, matrix, x, b); \
  solver->solve(); \
  \
  /* Check -- For a*y'' = b, correct answer is y = 0.5 *(b/a) * x * (x-1) */ \
  solver = Teuchos::null; /* Delete solver to eliminate live device views of x */ \
  typedef Teuchos::ScalarTraits<BaseScalar> ST; \
  typename ST::magnitudeType tol = 1e-9; \
  auto x_view = x->getLocalViewHost(Tpetra::Access::ReadOnly); \
  Scalar val(VectorSize, BaseScalar(0.0)); \
  for (size_t i=0; i<num_my_row; ++i) { \
    const GlobalOrdinal row = myGIDs[i]; \
    BaseScalar xx = row * h; \
    for (LocalOrdinal j=0; j<VectorSize; ++j) { \
      val.fastAccessCoeff(j) = \
        BaseScalar(0.5) * (b_val.coeff(j)/a_val.coeff(j)) * xx * (xx - BaseScalar(1.0)); \
    } \
    TEST_EQUALITY( x_view(i,0).size(), VectorSize ); \
  \
    /* Set small values to zero */ \
    Scalar v = x_view(i,0); \
    for (LocalOrdinal j=0; j<VectorSize; ++j) { \
      if (ST::magnitude(v.coeff(j)) < tol) \
        v.fastAccessCoeff(j) = BaseScalar(0.0); \
      if (ST::magnitude(val.coeff(j)) < tol) \
        val.fastAccessCoeff(j) = BaseScalar(0.0); \
    } \
  \
    for (LocalOrdinal j=0; j<VectorSize; ++j) \
      TEST_FLOATING_EQUALITY(v.coeff(j), val.coeff(j), tol); \
  } \
}

#endif // defined(HAVE_STOKHOS_AMESOS2)

#if defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_AMESOS2_BASKER)
  TEST_AMESOS2_SOLVER(basker)
#else
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_CrsMatrix_MP, Amesos2_basker, Storage, LocalOrdinal, GlobalOrdinal, Node) {}
#endif

#if defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_AMESOS2_KLU2)
  TEST_AMESOS2_SOLVER(klu2)
#else
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_CrsMatrix_MP, Amesos2_klu2, Storage, LocalOrdinal, GlobalOrdinal, Node) {}
#endif

#if defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_AMESOS2_SUPERLUDIST)
  TEST_AMESOS2_SOLVER(superlu_dist)
#else
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_CrsMatrix_MP, Amesos2_superlu_dist, Storage, LocalOrdinal, GlobalOrdinal, Node) {}
#endif

#if defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_AMESOS2_SUPERLUMT)
  TEST_AMESOS2_SOLVER(superlu_mt)
#else
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_CrsMatrix_MP, Amesos2_superlu_mt, Storage, LocalOrdinal, GlobalOrdinal, Node) {}
#endif

#if defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_AMESOS2_SUPERLU)
  TEST_AMESOS2_SOLVER(superlu)
#else
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_CrsMatrix_MP, Amesos2_superlu, Storage, LocalOrdinal, GlobalOrdinal, Node) {}
#endif

#if defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_AMESOS2_PARDISO_MKL)
  TEST_AMESOS2_SOLVER(pardisomkl)
#else
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_CrsMatrix_MP, Amesos2_pardisomkl, Storage, LocalOrdinal, GlobalOrdinal, Node) {}
#endif

#if defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_AMESOS2_LAPACK)
  TEST_AMESOS2_SOLVER(lapack)
#else
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_CrsMatrix_MP, Amesos2_lapack, Storage, LocalOrdinal, GlobalOrdinal, Node) {}
#endif

#if defined(HAVE_STOKHOS_AMESOS2) && defined(HAVE_AMESOS2_CHOLMOD) && defined (HAVE_AMESOS2_EXPERIMENTAL)
  TEST_AMESOS2_SOLVER(cholmod)
#else
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_CrsMatrix_MP, Amesos2_cholmod, Storage, LocalOrdinal, GlobalOrdinal, Node) {}
#endif

#define CRSMATRIX_MP_VECTOR_TESTS_SLGN(S, LO, GO, N)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, VectorAdd, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, VectorDot, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, MultiVectorAdd, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, MultiVectorDot, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, MultiVectorDotSub, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, MatrixVectorMultiply, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, MatrixMultiVectorMultiply, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, WrappedDualView, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Flatten, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, SimpleCG, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, SimplePCG_Muelu, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, BelosGMRES, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, BelosGMRES_DGKS, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, BelosGMRES_ICGS, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, BelosGMRES_IMGS, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, BelosGMRES_RILUK, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, BelosCG_Muelu, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Amesos2_basker, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Amesos2_klu2, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Amesos2_superlu_dist, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Amesos2_superlu_mt, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Amesos2_superlu, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Amesos2_pardisomkl, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Amesos2_lapack, S, LO, GO, N ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_CrsMatrix_MP, Amesos2_cholmod, S, LO, GO, N )

#define CRSMATRIX_MP_VECTOR_TESTS_N_SFS(N)                              \
  typedef Stokhos::DeviceForNode<N>::type Device;              \
  typedef Stokhos::StaticFixedStorage<int,double,VectorSize,Device::execution_space> SFS; \
  using default_global_ordinal_type = ::Tpetra::Map<>::global_ordinal_type; \
  using default_local_ordinal_type = ::Tpetra::Map<>::local_ordinal_type; \
  CRSMATRIX_MP_VECTOR_TESTS_SLGN(SFS, default_local_ordinal_type, default_global_ordinal_type, N)

#define CRSMATRIX_MP_VECTOR_TESTS_N(N)                                  \
  CRSMATRIX_MP_VECTOR_TESTS_N_SFS(N)

// Disabling testing of dynamic storage -- we don't really need it
  // typedef Stokhos::DynamicStorage<int,double,Device> DS;
  // using default_global_ordinal_type = ::Tpetra::Map<>::global_ordinal_type;
  // using default_local_ordinal_type = ::Tpetra::Map<>::local_ordinal_type;
  // CRSMATRIX_MP_VECTOR_TESTS_SLGN(DS, default_global_ordinal_type, default_local_ordinal_type, N)
