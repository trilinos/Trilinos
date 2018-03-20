#ifndef __TACHO_TEST_NUMERIC_HPP__
#define __TACHO_TEST_NUMERIC_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_MatrixMarket.hpp"

#include "Tacho_Graph.hpp"
#include "Tacho_SymbolicTools.hpp"

#if defined(TACHO_HAVE_SCOTCH)
#include "Tacho_GraphTools_Scotch.hpp"
#endif

#if defined(TACHO_HAVE_METIS)
#include "Tacho_GraphTools_Metis.hpp"
#endif

#include "Tacho_GraphTools_CAMD.hpp"

#include "Tacho_NumericTools.hpp"

using namespace Tacho;

typedef CrsMatrixBase<ValueType,HostSpaceType> CrsMatrixBaseHostType;
typedef CrsMatrixBase<ValueType,DeviceSpaceType> CrsMatrixBaseDeviceType;

TEST( Numeric, constructor ) {
  TEST_BEGIN;
  const ordinal_type
    m = 4,
    n = 4,
    nnz = 16;

  CrsMatrixBaseHostType A(m, n, nnz);
  ordinal_type cnt = 0;
  for (ordinal_type i=0;i<m;++i) {
    A.RowPtrBegin(i) = cnt;
    for (ordinal_type j=0;j<n;++j,++cnt) {
      A.Col(cnt) = j;
      A.Value(cnt) = i*n+j;
    }
    A.RowPtrEnd(i) = cnt;
  }

  typedef Kokkos::View<ordinal_type*,HostSpaceType> ordinal_type_array;

  ordinal_type_array idx("idx", m);
  for (ordinal_type i=0;i<m;++i) idx(i) = i;

  SymbolicTools S(m, A.RowPtr(), A.Cols(), idx, idx);
  S.symbolicFactorize();

  typedef typename DeviceSpaceType::memory_space device_memory_space; 

  auto a_row_ptr              = Kokkos::create_mirror_view(device_memory_space(), A.RowPtr());
  auto a_cols                 = Kokkos::create_mirror_view(device_memory_space(), A.Cols());

  auto d_idx                  = Kokkos::create_mirror_view(device_memory_space(), idx);
  auto s_supernodes           = Kokkos::create_mirror_view(device_memory_space(), S.Supernodes());
  auto s_gid_spanel_ptr       = Kokkos::create_mirror_view(device_memory_space(), S.gidSuperPanelPtr());
  auto s_gid_spanel_colidx    = Kokkos::create_mirror_view(device_memory_space(), S.gidSuperPanelColIdx());
  auto s_sid_spanel_ptr       = Kokkos::create_mirror_view(device_memory_space(), S.sidSuperPanelPtr());
  auto s_sid_spanel_colidx    = Kokkos::create_mirror_view(device_memory_space(), S.sidSuperPanelColIdx());
  auto s_blk_spanel_colidx    = Kokkos::create_mirror_view(device_memory_space(), S.blkSuperPanelColIdx());
  auto s_snodes_tree_parent   = Kokkos::create_mirror_view(device_memory_space(), S.SupernodesTreeParent());
  auto s_snodes_tree_ptr      = Kokkos::create_mirror_view(device_memory_space(), S.SupernodesTreePtr());
  auto s_snodes_tree_children = Kokkos::create_mirror_view(device_memory_space(), S.SupernodesTreeChildren());

  Kokkos::deep_copy(a_row_ptr              , A.RowPtr());
  Kokkos::deep_copy(a_cols                 , A.Cols());

  Kokkos::deep_copy(s_supernodes           , S.Supernodes());
  Kokkos::deep_copy(s_gid_spanel_ptr       , S.gidSuperPanelPtr());
  Kokkos::deep_copy(s_gid_spanel_colidx    , S.gidSuperPanelColIdx());
  Kokkos::deep_copy(s_sid_spanel_ptr       , S.sidSuperPanelPtr());
  Kokkos::deep_copy(s_sid_spanel_colidx    , S.sidSuperPanelColIdx());
  Kokkos::deep_copy(s_blk_spanel_colidx    , S.blkSuperPanelColIdx());
  Kokkos::deep_copy(s_snodes_tree_parent   , S.SupernodesTreeParent());
  Kokkos::deep_copy(s_snodes_tree_ptr      , S.SupernodesTreePtr());
  Kokkos::deep_copy(s_snodes_tree_children , S.SupernodesTreeChildren());

  NumericTools<ValueType,DeviceSpaceType> N(m, a_row_ptr, a_cols,
                                            d_idx, d_idx,
                                            S.NumSupernodes(), s_supernodes,
                                            s_gid_spanel_ptr, s_gid_spanel_colidx,
                                            s_sid_spanel_ptr, s_sid_spanel_colidx, s_blk_spanel_colidx,
                                            s_snodes_tree_parent, s_snodes_tree_ptr, s_snodes_tree_children,
                                            S.SupernodesTreeRoots());
  TEST_END;
}

#if !defined (KOKKOS_ENABLE_CUDA)
TEST( Numeric, Cholesky_Serial ) {
  TEST_BEGIN;
  std::string inputfilename = MM_TEST_FILE + ".mtx";
  CrsMatrixBaseDeviceType A;
  MatrixMarket<ValueType>::read(inputfilename, A);

  Graph G(A);

#if   defined(TACHO_HAVE_METIS)
  GraphTools_Metis T(G);
#elif defined(TACHO_HAVE_SCOTCH)
  GraphTools_Scotch T(G);
#else
  GraphTools_CAMD T(G);
#endif
  T.reorder();

  SymbolicTools S(A, T);
  S.symbolicFactorize();
  
  NumericTools<ValueType,DeviceSpaceType> N(A.NumRows(), A.RowPtr(), A.Cols(), 
                                            T.PermVector(), T.InvPermVector(),
                                            S.NumSupernodes(), S.Supernodes(),
                                            S.gidSuperPanelPtr(), S.gidSuperPanelColIdx(),
                                            S.sidSuperPanelPtr(), S.sidSuperPanelColIdx(), S.blkSuperPanelColIdx(),
                                            S.SupernodesTreeParent(), S.SupernodesTreePtr(), S.SupernodesTreeChildren(), S.SupernodesTreeRoots());
  N.factorizeCholesky_Serial(A.Values());
  
  CrsMatrixBaseDeviceType F;
  N.exportFactorsToCrsMatrix(F);
  
  std::ofstream out("test_numeric_factorize_serial.mtx");
  MatrixMarket<ValueType>::write(out, F);
  
  const ordinal_type m = A.NumRows(), n = 2;
  Kokkos::View<ValueType**,Kokkos::LayoutLeft,DeviceSpaceType> 
    x("x", m, n), b("b", m, n), t("t", m, n);
  
  Random<ValueType> random;
  for (ordinal_type j=0;j<n;++j)
    for (ordinal_type i=0;i<m;++i) 
      b(i,j) = random.value();
  
  N.solveCholesky_Serial(x, b, t);
  
  const double eps = std::numeric_limits<double>::epsilon()*100;
  EXPECT_TRUE(N.computeRelativeResidual(x,b) < eps);
  TEST_END;
}
#endif

TEST( Numeric, factorizeCholesky_Parallel ) {
  TEST_BEGIN;
  std::string inputfilename = MM_TEST_FILE + ".mtx";
  CrsMatrixBaseHostType A;
  MatrixMarket<ValueType>::read(inputfilename, A);

  Graph G(A);

#if   defined(TACHO_HAVE_METIS)
  GraphTools_Metis T(G);
#elif defined(TACHO_HAVE_SCOTCH)
  GraphTools_Scotch T(G);
#else
  GraphTools_CAMD T(G);
#endif
  T.reorder();

  SymbolicTools S(A, T);
  S.symbolicFactorize();

  typedef typename DeviceSpaceType::memory_space device_memory_space; 

  auto a_row_ptr              = Kokkos::create_mirror_view(device_memory_space(), A.RowPtr());
  auto a_cols                 = Kokkos::create_mirror_view(device_memory_space(), A.Cols());
  auto a_values               = Kokkos::create_mirror_view(device_memory_space(), A.Values());

  auto t_perm                 = Kokkos::create_mirror_view(device_memory_space(), T.PermVector());
  auto t_peri                 = Kokkos::create_mirror_view(device_memory_space(), T.InvPermVector());
  auto s_supernodes           = Kokkos::create_mirror_view(device_memory_space(), S.Supernodes());
  auto s_gid_spanel_ptr       = Kokkos::create_mirror_view(device_memory_space(), S.gidSuperPanelPtr());
  auto s_gid_spanel_colidx    = Kokkos::create_mirror_view(device_memory_space(), S.gidSuperPanelColIdx());
  auto s_sid_spanel_ptr       = Kokkos::create_mirror_view(device_memory_space(), S.sidSuperPanelPtr());
  auto s_sid_spanel_colidx    = Kokkos::create_mirror_view(device_memory_space(), S.sidSuperPanelColIdx());
  auto s_blk_spanel_colidx    = Kokkos::create_mirror_view(device_memory_space(), S.blkSuperPanelColIdx());
  auto s_snodes_tree_parent   = Kokkos::create_mirror_view(device_memory_space(), S.SupernodesTreeParent());
  auto s_snodes_tree_ptr      = Kokkos::create_mirror_view(device_memory_space(), S.SupernodesTreePtr());
  auto s_snodes_tree_children = Kokkos::create_mirror_view(device_memory_space(), S.SupernodesTreeChildren());

  Kokkos::deep_copy(a_row_ptr              , A.RowPtr());
  Kokkos::deep_copy(a_cols                 , A.Cols());
  Kokkos::deep_copy(a_values               , A.Values());

  Kokkos::deep_copy(t_perm                 , T.PermVector());
  Kokkos::deep_copy(t_peri                 , T.InvPermVector());
  Kokkos::deep_copy(s_supernodes           , S.Supernodes());
  Kokkos::deep_copy(s_gid_spanel_ptr       , S.gidSuperPanelPtr());
  Kokkos::deep_copy(s_gid_spanel_colidx    , S.gidSuperPanelColIdx());
  Kokkos::deep_copy(s_sid_spanel_ptr       , S.sidSuperPanelPtr());
  Kokkos::deep_copy(s_sid_spanel_colidx    , S.sidSuperPanelColIdx());
  Kokkos::deep_copy(s_blk_spanel_colidx    , S.blkSuperPanelColIdx());
  Kokkos::deep_copy(s_snodes_tree_parent   , S.SupernodesTreeParent());
  Kokkos::deep_copy(s_snodes_tree_ptr      , S.SupernodesTreePtr());
  Kokkos::deep_copy(s_snodes_tree_children , S.SupernodesTreeChildren());

  NumericTools<ValueType,DeviceSpaceType> N(A.NumRows(), a_row_ptr, a_cols,
                                            t_perm, t_peri,
                                            S.NumSupernodes(), s_supernodes,
                                            s_gid_spanel_ptr, s_gid_spanel_colidx,
                                            s_sid_spanel_ptr, s_sid_spanel_colidx, s_blk_spanel_colidx,
                                            s_snodes_tree_parent, s_snodes_tree_ptr, s_snodes_tree_children,
                                            S.SupernodesTreeRoots());
  
  N.factorizeCholesky_Parallel(a_values);
  TEST_END;
}

#endif



// debug print

  // const auto m = A.NumRows();
  // const auto ap = A.RowPtr();
  // const auto aj = A.Cols();
  // const auto ax = A.Values();
  
  // const auto perm = T.PermVector();
  // const auto peri = T.InvPermVector();
  
  // const size_type ns = S.NumSupernodes();
  // const auto supernodes = S.Supernodes();
  // const auto gid_super_panel_ptr =  S.gidSuperPanelPtr();
  // const auto gid_super_panel_colidx = S.gidSuperPanelColIdx();
  // const auto sid_super_panel_ptr = S.sidSuperPanelPtr();
  // const auto sid_super_panel_colidx = S.sidSuperPanelColIdx();
  // const auto blk_super_panel_colidx = S.blkSuperPanelColIdx();

  // printf("supernodes = %d\n", ns);
  // for (size_type i=0;i<ns;++i) {
  //   printf("sid=  %d\n", supernodes(i));
  //   printf("-- gid = ");
  //   for (size_type j=gid_super_panel_ptr(i);j<gid_super_panel_ptr(i+1);++j)
  //     printf("  %d", gid_super_panel_colidx(j));
  //   printf("\n");

  //   printf("-- supernodes id connected = ");
  //   for (size_type j=sid_super_panel_ptr(i);j<(sid_super_panel_ptr(i+1)-1);++j)
  //     printf("  %d", sid_super_panel_colidx(j));
  //   printf("\n");
  //   printf("-- supernodes blocks connected = ");
  //   for (size_type j=sid_super_panel_ptr(i);j<sid_super_panel_ptr(i+1);++j)
  //     printf("  %d", blk_super_panel_colidx(j));
  //   printf("\n");
  // }

