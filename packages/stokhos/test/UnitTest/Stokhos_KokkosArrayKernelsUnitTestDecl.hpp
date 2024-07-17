// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CrsMatrixFree, Scalar, Device ) {
  typedef Stokhos::DefaultMultiply SparseMatOps;
  success = test_crs_matrix_free<Scalar,Device,SparseMatOps>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CrsMatrixFreeView, Scalar, Device ) {
  typedef Stokhos::DefaultMultiply SparseMatOps;
  success = test_crs_matrix_free_view<Scalar,Device,SparseMatOps>(setup, out);
}

#ifdef HAVE_STOKHOS_KOKKOSLINALG
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CrsMatrixFreeKokkos, Scalar, Device ) {
  success = test_crs_matrix_free_kokkos<Scalar,Device>(setup, out);
}
#else
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CrsMatrixFreeKokkos, Scalar, Device ) {}
#endif

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CrsMatrixFreeSingleCol, Scalar, Device ) {
  typedef Stokhos::SingleColumnMultivectorMultiply SparseMatOps;
  success = test_crs_matrix_free<Scalar,Device,SparseMatOps>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CrsDenseBlock, Scalar, Device) {
  success = test_crs_dense_block<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CrsFlatCommuted, Scalar, Device ) {
  success = test_crs_flat_commuted<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CrsFlatOriginal, Scalar, Device ) {
  success = test_crs_flat_original<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CrsProductTensor, Scalar, Device ) {
  typedef Stokhos::CrsProductTensor<Scalar,Device> Tensor;
  success = test_crs_product_tensor<Scalar,Tensor,Device>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, TiledCrsProductTensor, Scalar, Device ) {
  typedef Stokhos::TiledCrsProductTensor<Scalar,Device> Tensor;
  Teuchos::ParameterList params;
  params.set("Tile Size", 10);
  params.set("Max Tiles", 10000);
  success = test_crs_product_tensor<Scalar,Tensor,Device>(setup, out, params);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, SimpleTiledCrsProductTensor, Scalar, Device ){
  typedef Stokhos::SimpleTiledCrsProductTensor<Scalar,Device> Tensor;
  Teuchos::ParameterList params;
  params.set("Tile Size", 10);
  success = test_crs_product_tensor<Scalar,Tensor,Device>(setup, out, params);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CooProductTensorPacked, Scalar, Device ) {
  typedef Stokhos::CooProductTensor<Scalar,Device,true> Tensor;
  success = test_crs_product_tensor<Scalar,Tensor,Device>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CooProductTensorUnpacked, Scalar, Device ) {
  typedef Stokhos::CooProductTensor<Scalar,Device,false> Tensor;
  success = test_crs_product_tensor<Scalar,Tensor,Device>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, FlatSparse3Tensor, Scalar, Device ) {
  typedef Stokhos::FlatSparse3Tensor<Scalar,Device> Tensor;
  success = test_crs_product_tensor<Scalar,Tensor,Device>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, FlatSparse3Tensor_kji, Scalar, Device ) {
  typedef Stokhos::FlatSparse3Tensor_kji<Scalar,Device> Tensor;
  success = test_crs_product_tensor<Scalar,Tensor,Device>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, LinearTensorSymmetric, Scalar, Device ) {
  const bool symmetric = true;
  UnitTestSetup<Device> s;
  s.setup(1, 10);
  success = test_linear_tensor<Scalar,Device,4>(s, out, symmetric);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, LinearTensorAsymmetric, Scalar, Device ) {
  const bool symmetric = false;
  UnitTestSetup<Device> s;
  s.setup(1, 10);
  success = test_linear_tensor<Scalar,Device,4>(s, out, symmetric);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, LexoBlockTensor, Scalar, Device ) {
  success = test_lexo_block_tensor<Scalar,Device>(setup, out);
}

// ETP 6/23/14:  CooProductTensor tests are failing with Intel compiler
// (optimized with AVX).  Probably alignment issue

// ETP 5/12/16:  LinearProductTensor tests are failing with Intel compiler too.

#define UNIT_TEST_GROUP_SCALAR_DEVICE( SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsMatrixFree, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsMatrixFreeView, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsMatrixFreeKokkos, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsMatrixFreeSingleCol, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsDenseBlock, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsFlatCommuted, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsFlatOriginal, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsProductTensor, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, TiledCrsProductTensor, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, SimpleTiledCrsProductTensor, SCALAR, DEVICE ) \
  /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CooProductTensorPacked, SCALAR, DEVICE )*/ \
  /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CooProductTensorUnpacked, SCALAR, DEVICE )*/ \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, FlatSparse3Tensor, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, FlatSparse3Tensor_kji, SCALAR, DEVICE ) \
  /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, LinearTensorSymmetric, SCALAR, DEVICE )*/ \
  /*TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, LinearTensorAsymmetric, SCALAR, DEVICE )*/ \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, LexoBlockTensor, SCALAR, DEVICE )
