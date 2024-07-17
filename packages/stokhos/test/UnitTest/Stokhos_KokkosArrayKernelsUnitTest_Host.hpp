// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Host-specific tests
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CrsProductTensorCijk, Scalar, Device ) {
  success = true;

  typedef Scalar value_type;
  typedef Stokhos::CrsProductTensor< value_type , Device > tensor_type ;

  tensor_type tensor =
   Stokhos::create_product_tensor<Device>( *setup.basis, *setup.Cijk );

  for (int i=0; i<setup.stoch_length; ++i) {
    const int iEntryBeg = tensor.entry_begin(i);
    const int iEntryEnd = tensor.entry_end(i);
    for (int iEntry = iEntryBeg ; iEntry < iEntryEnd ; ++iEntry ) {
      const int kj = tensor.coord( iEntry );
      const int j  = kj & 0x0ffff;
      const int k  = kj >> 16;
      // const int j = tensor.coord(iEntry,0);
      // const int k = tensor.coord(iEntry,1);
      value_type c2 = tensor.value(iEntry);
      if (j == k) c2 *= 2.0;

      int ii = setup.inv_perm[i];
      int jj = setup.inv_perm[j];
      int kk = setup.inv_perm[k];
      value_type c = setup.Cijk->getValue(ii,jj,kk);

      if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
        out << "(" << ii << "," << jj << "," << kk << "):  " << c
            << " == " << c2 << " failed!" << std::endl;
        success = false;
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, TiledCrsProductTensorCijk, Scalar, Device ) {
  success = true;

  typedef Scalar value_type;
  typedef Stokhos::TiledCrsProductTensor< value_type , Device > tensor_type ;

  Teuchos::ParameterList params;
  params.set("Tile Size",10);
  params.set("Max Tiles",10000);

  tensor_type tensor =
    Stokhos::create_tiled_product_tensor<Device>( *setup.basis, *setup.Cijk,
                                                  params );

  // This is a valid test only with no symmetry
  // TEUCHOS_TEST_EQUALITY( tensor.entry_count(), setup.Cijk->num_entries(),
  //                        out, success );

  const size_t n_tile = tensor.num_tiles();
  for ( size_t tile = 0 ; tile < n_tile ; ++tile ) {
    const size_t i_offset = tensor.offset(tile, 0);
    const size_t j_offset = tensor.offset(tile, 1);
    const size_t k_offset = tensor.offset(tile, 2);
    const size_t n_row = tensor.num_rows(tile);

    for (size_t i=0; i<n_row; ++i) {
      const size_t iEntryBeg = tensor.entry_begin(tile,i);
      const size_t iEntryEnd = tensor.entry_end(tile,i);
      for (size_t iEntry = iEntryBeg ; iEntry < iEntryEnd ; ++iEntry ) {
        const size_t j = tensor.coord(iEntry,0);
        const size_t k = tensor.coord(iEntry,1);
        value_type c2 = tensor.value(iEntry);
        int ii = i + i_offset;
        int jj = j + j_offset;
        int kk = k + k_offset;
        if (jj == kk)
          c2 *= 2.0;
        value_type c = setup.Cijk->getValue(ii,jj,kk);

        if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
          out << "(" << ii << "," << jj << "," << kk << "):  " << c
              << " == " << c2 << " failed!" << std::endl;
          success = false;
        }
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, SimpleTiledCrsProductTensorCijk, Scalar, Device ) {
  success = true;

  typedef Scalar value_type;
  typedef Stokhos::SimpleTiledCrsProductTensor< value_type , Device > tensor_type ;

  Teuchos::ParameterList params;
  params.set("Tile Size",10);

  tensor_type tensor =
    Stokhos::create_simple_tiled_product_tensor<Device>(
      *setup.basis, *setup.Cijk, params);

  int num_entry = 0;
  const size_t n_i_tile = tensor.num_i_tiles();
  for (size_t i_tile = 0; i_tile<n_i_tile; ++i_tile) {
    const size_t i_begin = tensor.i_begin(i_tile);
    const size_t i_size  = tensor.i_size(i_tile);

    const size_t n_j_tile = tensor.num_j_tiles(i_tile);
    for (size_t j_tile = 0; j_tile<n_j_tile; ++j_tile) {
      const size_t j_begin = tensor.j_begin(i_tile, j_tile);
      //const size_t j_size  = tensor.j_size(i_tile, j_tile);

      const size_t n_k_tile = tensor.num_k_tiles(i_tile, j_tile);
      for (size_t k_tile = 0; k_tile<n_k_tile; ++k_tile) {
        const size_t k_begin = tensor.k_begin(i_tile, j_tile, k_tile);
        //const size_t k_size  = tensor.k_size(i_tile, j_tile, k_tile);

        for (size_t i=0; i<i_size; ++i) {
          const size_t iEntryBeg = tensor.entry_begin(i_tile,j_tile,k_tile,i);
          const size_t iEntryEnd = tensor.entry_end(i_tile,j_tile,k_tile,i);
          for (size_t iEntry = iEntryBeg ; iEntry < iEntryEnd ; ++iEntry ) {
            const size_t j = tensor.coord(iEntry,0);
            const size_t k = tensor.coord(iEntry,1);
            value_type c2 = tensor.value(iEntry);
            int ii = i + i_begin;
            int jj = j + j_begin;
            int kk = k + k_begin;
            ++num_entry;
            if (jj == kk)
              c2 *= 2.0;
            else
              ++num_entry;
            value_type c = setup.Cijk->getValue(ii,jj,kk);

            if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
              out << "(" << ii << "," << jj << "," << kk << "):  " << c
                  << " == " << c2 << " failed!" << std::endl;
              success = false;
            }
          }
        }
      }
    }
  }
  TEUCHOS_TEST_EQUALITY( num_entry, setup.Cijk->num_entries(), out, success );
}

template <typename Scalar, typename Device, bool Pack>
bool test_coo_product_tensor_cijk(
  const KokkosKernelsUnitTest::UnitTestSetup<Device>& setup,
  Teuchos::FancyOStream& out) {
  bool success = true;

  typedef Scalar value_type;
  typedef Stokhos::CooProductTensor< value_type , Device , Pack > tensor_type ;

  tensor_type tensor =
    Stokhos::create_coo_product_tensor<Device, Pack>(
      *setup.basis, *setup.Cijk );

  const size_t nEntry = tensor.entry_count();
  size_t i, j, k;
  for ( size_t entry = 0 ; entry < nEntry ; ++entry ) {
    tensor.coord(entry, i, j, k);
    value_type c2 = tensor.value(entry);
    if (j == k) c2 *= 2.0;
    value_type c = setup.Cijk->getValue(i,j,k);

    if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
      out << "(" << i << "," << j << "," << k << "):  " << c
          << " == " << c2 << " failed!" << std::endl;
      success = false;
    }
  }

  return success;
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CooProductTensorCijk_Packed, Scalar, Device ) {
  success = test_coo_product_tensor_cijk<Scalar,Device,true>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, CooProductTensorCijk_Unpacked, Scalar, Device ) {
  success = test_coo_product_tensor_cijk<Scalar,Device,false>(setup, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, FlatSparseCijk, Scalar, Device ) {
  success = true;

  typedef Scalar value_type;
  typedef Stokhos::FlatSparse3Tensor< value_type , Device > tensor_type ;
  typedef size_t size_type;

  tensor_type tensor =
   Stokhos::create_flat_sparse_3_tensor<Device>( *setup.basis, *setup.Cijk );

  for (int i=0; i<setup.stoch_length; ++i) {
    const size_type nk = tensor.num_k(i);
    const size_type kBeg = tensor.k_begin(i);
    const size_type kEnd = kBeg + nk;
    for (size_type kEntry = kBeg; kEntry < kEnd; ++kEntry) {
      const size_type k = tensor.k_coord(kEntry);
      const size_type nj = tensor.num_j(kEntry);
      const size_type jBeg = tensor.j_begin(kEntry);
      const size_type jEnd = jBeg + nj;
      for (size_type jEntry = jBeg; jEntry < jEnd; ++jEntry) {
        const size_type j = tensor.j_coord(jEntry);
        value_type c2 = tensor.value(jEntry);
        if (j == k) c2 *= 2.0;
        value_type c = setup.Cijk->getValue(i,j,k);
        if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
          out << "(" << i << "," << j << "," << k << "):  " << c
              << " == " << c2 << " failed!" << std::endl;
          success = false;
        }
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_SG_SpMv, FlatSparseCijk_kji, Scalar, Device ) {
  success = true;

  typedef Scalar value_type;
  typedef Stokhos::FlatSparse3Tensor_kji< value_type , Device > tensor_type ;
  typedef size_t size_type;

  tensor_type tensor =
   Stokhos::create_flat_sparse_3_tensor_kji<Device>(*setup.basis, *setup.Cijk);
  const size_type nk = tensor.num_k();

  for ( size_type k = 0; k < nk; ++k) {
    const size_type nj = tensor.num_j(k);
    const size_type jBeg = tensor.j_begin(k);
    const size_type jEnd = jBeg + nj;
    for (size_type jEntry = jBeg; jEntry < jEnd; ++jEntry) {
      const size_type j = tensor.j_coord(jEntry);
      const size_type ni = tensor.num_i(jEntry);
      const size_type iBeg = tensor.i_begin(jEntry);
      const size_type iEnd = iBeg + ni;
      for (size_type iEntry = iBeg; iEntry < iEnd; ++iEntry) {
        const size_type i = tensor.i_coord(iEntry);
        value_type c2 = tensor.value(iEntry);
        if (j == k) c2 *= 2.0;
        value_type c = setup.Cijk->getValue(i,j,k);
        if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
          out << "(" << i << "," << j << "," << k << "):  " << c
              << " == " << c2 << " failed!" << std::endl;
          success = false;
        }
      }
    }
  }
}

#define UNIT_TEST_GROUP_SCALAR_HOST_DEVICE( SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsProductTensorCijk, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, TiledCrsProductTensorCijk, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, SimpleTiledCrsProductTensorCijk, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CooProductTensorCijk_Packed, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CooProductTensorCijk_Unpacked, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, FlatSparseCijk, SCALAR, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, FlatSparseCijk_kji, SCALAR, DEVICE )
