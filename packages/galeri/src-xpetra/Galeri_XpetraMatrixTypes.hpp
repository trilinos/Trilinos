// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/*
  Direct translation of parts of Galeri to use Tpetra or Xpetra rather than Epetra. Epetra also supported.
*/

// TODO: rename variables (camelCase)

#ifndef GALERI_XPETRAMATRIXTYPES_HPP
#define GALERI_XPETRAMATRIXTYPES_HPP
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Galeri_MapTraits.hpp"
#include "Galeri_config.h"
#include "Galeri_MatrixTraits.hpp"
#include "Galeri_Problem.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_ScalarTraits.hpp"
#if defined(HAVE_GALERI_KOKKOS) && defined(HAVE_GALERI_KOKKOSKERNELS)
#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Tpetra_Distributor.hpp"
#endif

namespace Galeri {

namespace Xpetra {

/* prototypes */
template <typename GlobalOrdinal>
bool IsBoundary2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny);
template <typename GlobalOrdinal>
bool IsBoundary3d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz);

template <typename GlobalOrdinal>
void GetNeighboursCartesian2d(const GlobalOrdinal i,
                              const GlobalOrdinal nx, const GlobalOrdinal ny,
                              GlobalOrdinal& left, GlobalOrdinal& right,
                              GlobalOrdinal& lower, GlobalOrdinal& upper);

template <typename GlobalOrdinal>
void GetNeighboursCartesian2d(GlobalOrdinal i, GlobalOrdinal nx, const GlobalOrdinal ny,
                              GlobalOrdinal& left, GlobalOrdinal& right, GlobalOrdinal& lower, GlobalOrdinal& upper,
                              GlobalOrdinal& left2, GlobalOrdinal& right2, GlobalOrdinal& lower2, GlobalOrdinal& upper2);

template <typename GlobalOrdinal>
void GetNeighboursCartesian3d(const GlobalOrdinal i,
                              const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                              GlobalOrdinal& left, GlobalOrdinal& right,
                              GlobalOrdinal& front, GlobalOrdinal& back,
                              GlobalOrdinal& bottom, GlobalOrdinal& top);

#if defined(HAVE_GALERI_KOKKOS) && defined(HAVE_GALERI_KOKKOSKERNELS)

template <typename GlobalOrdinal>
KOKKOS_FORCEINLINE_FUNCTION void
GetNeighboursCartesian1dKokkos(const GlobalOrdinal i,
                               const GlobalOrdinal nx, GlobalOrdinal& left,
                               GlobalOrdinal& right,
                               const GlobalOrdinal INVALID) {
  if (i == 0)
    left = INVALID;
  else
    left = i - 1;
  if (i == nx - 1)
    right = INVALID;
  else
    right = i + 1;
}

template <typename GlobalOrdinal>
KOKKOS_FORCEINLINE_FUNCTION void GetNeighboursCartesian2dKokkos(const GlobalOrdinal i,
                                                                const GlobalOrdinal nx, const GlobalOrdinal ny,
                                                                GlobalOrdinal& left, GlobalOrdinal& right,
                                                                GlobalOrdinal& lower, GlobalOrdinal& upper,
                                                                const GlobalOrdinal INVALID) {
  GlobalOrdinal ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  if (ix == 0)
    left = INVALID;
  else
    left = i - 1;
  if (ix == nx - 1)
    right = INVALID;
  else
    right = i + 1;
  if (iy == 0)
    lower = INVALID;
  else
    lower = i - nx;
  if (iy == ny - 1)
    upper = INVALID;
  else
    upper = i + nx;
}

template <typename GlobalOrdinal>
KOKKOS_FORCEINLINE_FUNCTION void GetNeighboursCartesian3dKokkos(const GlobalOrdinal i,
                                                                const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                                                                GlobalOrdinal& left, GlobalOrdinal& right,
                                                                GlobalOrdinal& front, GlobalOrdinal& back,
                                                                GlobalOrdinal& bottom, GlobalOrdinal& top,
                                                                const GlobalOrdinal INVALID) {
  GlobalOrdinal ixy, iz;
  ixy = i % (nx * ny);

  iz = (i - ixy) / (nx * ny);

  if (iz == 0)
    bottom = INVALID;
  else
    bottom = i - nx * ny;
  if (iz == nz - 1)
    top = INVALID;
  else
    top = i + nx * ny;

  GetNeighboursCartesian2dKokkos(ixy, nx, ny, left, right, front, back, INVALID);

  if (left != INVALID) left += iz * (nx * ny);
  if (right != INVALID) right += iz * (nx * ny);
  if (front != INVALID) front += iz * (nx * ny);
  if (back != INVALID) back += iz * (nx * ny);
}

// This macro is used when we enter count up the number of nonzero entries per row.
#define Galeri_processEntry(entry)                      \
  {                                                     \
    if (entry != INVALID) {                             \
      ++partial_nnz;                                    \
      if (is_final) {                                   \
        if (lclMap.getLocalElement(entry) == INVALID) { \
          off_rank_indices.insert(entry);               \
        }                                               \
      }                                                 \
    }                                                   \
  }

// This macro is used when we enter values in the matrix.
// We insert in a way that guarantees for rows to be sorted.
#define Galeri_enterValue(rowStart, column_index, value)      \
  {                                                           \
    if (column_index != INVALID) {                            \
      auto clid = lclColMap.getLocalElement(column_index);    \
      auto K    = entryPtr;                                   \
      while ((rowStart + 1 <= K) && (clid < colidx(K - 1))) { \
        colidx(K) = colidx(K - 1);                            \
        values(K) = values(K - 1);                            \
        --K;                                                  \
      }                                                       \
      colidx(K) = clid;                                       \
      values(K) = value;                                      \
      if (i != column_index)                                  \
        offDiagonalSum += value;                              \
      ++entryPtr;                                             \
    }                                                         \
  }

template <class Scalar, class Map>
class ScaledIdentityStencil {
  // a

#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<Scalar>;
#else
  using ATS                   = Kokkos::ArithTraits<Scalar>;
#endif
  using impl_scalar_type  = typename ATS::val_type;
  using LocalOrdinal      = typename Map::local_ordinal_type;
  using GlobalOrdinal     = typename Map::global_ordinal_type;
  using Node              = typename Map::node_type;
  using local_map_type    = typename Map::local_map_type;
  using local_matrix_type = typename ::Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;
  using rowptr_type       = typename local_matrix_type::row_map_type::non_const_type;
  using colidx_type       = typename local_matrix_type::index_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using exec_space        = typename Node::execution_space;
  using memory_space      = typename Node::memory_space;
  using hashmap_type      = typename Kokkos::UnorderedMap<GlobalOrdinal, void, exec_space>;

  GlobalOrdinal nx;
  impl_scalar_type a;
  local_map_type lclMap;

#if KOKKOS_VERSION >= 40799
  const impl_scalar_type one = KokkosKernels::ArithTraits<impl_scalar_type>::one();
#else
  const impl_scalar_type one  = Kokkos::ArithTraits<impl_scalar_type>::one();
#endif
#if KOKKOS_VERSION >= 40799
  const impl_scalar_type zero = KokkosKernels::ArithTraits<impl_scalar_type>::zero();
#else
  const impl_scalar_type zero = Kokkos::ArithTraits<impl_scalar_type>::zero();
#endif
  const GlobalOrdinal INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();

 public:
  hashmap_type off_rank_indices;
  local_map_type lclColMap;
  colidx_type colidx;
  values_type values;

  ScaledIdentityStencil(const Map& map, GlobalOrdinal nx_, Scalar a_)
    : nx(nx_)
    , a(a_) {
    lclMap = map.getLocalMap();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void CountRowNNZ(const GlobalOrdinal i, LocalOrdinal& partial_nnz, const bool is_final) const {
    GlobalOrdinal center = lclMap.getGlobalElement(i);
    Galeri_processEntry(center);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void EnterValues(const LocalOrdinal i, typename rowptr_type::value_type& entryPtr) const {
    GlobalOrdinal center                      = lclMap.getGlobalElement(i);
    impl_scalar_type offDiagonalSum           = zero;
    typename rowptr_type::value_type rowStart = entryPtr;
    Galeri_enterValue(rowStart, center, a);
  }
};

template <class Scalar, class Map, bool keepBCs>
class TriDiagStencil {
  //  b a c

#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<Scalar>;
#else
  using ATS                   = Kokkos::ArithTraits<Scalar>;
#endif
  using impl_scalar_type  = typename ATS::val_type;
  using LocalOrdinal      = typename Map::local_ordinal_type;
  using GlobalOrdinal     = typename Map::global_ordinal_type;
  using Node              = typename Map::node_type;
  using local_map_type    = typename Map::local_map_type;
  using local_matrix_type = typename ::Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;
  using rowptr_type       = typename local_matrix_type::row_map_type::non_const_type;
  using colidx_type       = typename local_matrix_type::index_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using exec_space        = typename Node::execution_space;
  using memory_space      = typename Node::memory_space;
  using hashmap_type      = typename Kokkos::UnorderedMap<GlobalOrdinal, void, exec_space>;

  GlobalOrdinal nx;
  impl_scalar_type a, b, c;
  DirBC DirichletBC;
  local_map_type lclMap;

#if KOKKOS_VERSION >= 40799
  const impl_scalar_type one = KokkosKernels::ArithTraits<impl_scalar_type>::one();
#else
  const impl_scalar_type one  = Kokkos::ArithTraits<impl_scalar_type>::one();
#endif
#if KOKKOS_VERSION >= 40799
  const impl_scalar_type zero = KokkosKernels::ArithTraits<impl_scalar_type>::zero();
#else
  const impl_scalar_type zero = Kokkos::ArithTraits<impl_scalar_type>::zero();
#endif
  const GlobalOrdinal INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();

 public:
  hashmap_type off_rank_indices;
  local_map_type lclColMap;
  colidx_type colidx;
  values_type values;

  TriDiagStencil(const Map& map, GlobalOrdinal nx_, Scalar a_, Scalar b_, Scalar c_, DirBC DirichletBC_)
    : nx(nx_)
    , a(a_)
    , b(b_)
    , c(c_)
    , DirichletBC(DirichletBC_) {
    lclMap = map.getLocalMap();
    TEUCHOS_ASSERT((GlobalOrdinal)map.getGlobalNumElements() == nx);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void GetNeighbours(const GlobalOrdinal i,
                     GlobalOrdinal& left, GlobalOrdinal& right, bool& isDirichlet) const {
    GetNeighboursCartesian1dKokkos(i, nx, left, right, INVALID);
    isDirichlet = (left == INVALID && (DirichletBC & DIR_LEFT)) ||
                  (right == INVALID && (DirichletBC & DIR_RIGHT));
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool IsBoundary(const GlobalOrdinal i) const {
    return (i == 0 || i == nx - 1);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void CountRowNNZ(const GlobalOrdinal i, LocalOrdinal& partial_nnz, const bool is_final) const {
    GlobalOrdinal center, left, right;
    bool isDirichlet;

    center = lclMap.getGlobalElement(i);
    GetNeighbours(center, left, right, isDirichlet);

    if (isDirichlet && keepBCs) {
      // Dirichlet unknown we want to keep
      Galeri_processEntry(center);
    } else {
      Galeri_processEntry(left);
      Galeri_processEntry(right);
      Galeri_processEntry(center);
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void EnterValues(const LocalOrdinal i, typename rowptr_type::value_type& entryPtr) const {
    GlobalOrdinal center, left, right;
    bool isDirichlet;

    center = lclMap.getGlobalElement(i);
    GetNeighbours(center, left, right, isDirichlet);

    impl_scalar_type offDiagonalSum           = zero;
    typename rowptr_type::value_type rowStart = entryPtr;
    if (isDirichlet && keepBCs) {
      // Dirichlet unknown we want to keep
      Galeri_enterValue(rowStart, center, one);
    } else {
      Galeri_enterValue(rowStart, left, b);
      Galeri_enterValue(rowStart, right, c);
      Galeri_enterValue(rowStart, center, (IsBoundary(center) && !isDirichlet) ? -offDiagonalSum : a);
    }
  }
};

template <class Scalar, class Map, bool keepBCs>
class Cross2DStencil {
  //    e
  //  b a c
  //    d

#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<Scalar>;
#else
  using ATS                   = Kokkos::ArithTraits<Scalar>;
#endif
  using impl_scalar_type  = typename ATS::val_type;
  using LocalOrdinal      = typename Map::local_ordinal_type;
  using GlobalOrdinal     = typename Map::global_ordinal_type;
  using Node              = typename Map::node_type;
  using local_map_type    = typename Map::local_map_type;
  using local_matrix_type = typename ::Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;
  using rowptr_type       = typename local_matrix_type::row_map_type::non_const_type;
  using colidx_type       = typename local_matrix_type::index_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using exec_space        = typename Node::execution_space;
  using memory_space      = typename Node::memory_space;
  using hashmap_type      = typename Kokkos::UnorderedMap<GlobalOrdinal, void, exec_space>;

  GlobalOrdinal nx, ny;
  impl_scalar_type a, b, c, d, e;
  DirBC DirichletBC;
  local_map_type lclMap;

#if KOKKOS_VERSION >= 40799
  const impl_scalar_type one = KokkosKernels::ArithTraits<impl_scalar_type>::one();
#else
  const impl_scalar_type one  = Kokkos::ArithTraits<impl_scalar_type>::one();
#endif
#if KOKKOS_VERSION >= 40799
  const impl_scalar_type zero = KokkosKernels::ArithTraits<impl_scalar_type>::zero();
#else
  const impl_scalar_type zero = Kokkos::ArithTraits<impl_scalar_type>::zero();
#endif
  const GlobalOrdinal INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();

 public:
  hashmap_type off_rank_indices;
  local_map_type lclColMap;
  colidx_type colidx;
  values_type values;

  Cross2DStencil(const Map& map, GlobalOrdinal nx_, GlobalOrdinal ny_, Scalar a_, Scalar b_, Scalar c_, Scalar d_, Scalar e_, DirBC DirichletBC_)
    : nx(nx_)
    , ny(ny_)
    , a(a_)
    , b(b_)
    , c(c_)
    , d(d_)
    , e(e_)
    , DirichletBC(DirichletBC_) {
    lclMap = map.getLocalMap();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void GetNeighbours(const GlobalOrdinal i,
                     GlobalOrdinal& left, GlobalOrdinal& right,
                     GlobalOrdinal& lower, GlobalOrdinal& upper, bool& isDirichlet) const {
    GetNeighboursCartesian2dKokkos(i, nx, ny, left, right, lower, upper, INVALID);
    isDirichlet = (left == INVALID && (DirichletBC & DIR_LEFT)) ||
                  (right == INVALID && (DirichletBC & DIR_RIGHT)) ||
                  (lower == INVALID && (DirichletBC & DIR_BOTTOM)) ||
                  (upper == INVALID && (DirichletBC & DIR_TOP));
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool IsBoundary(const GlobalOrdinal i) const {
    GlobalOrdinal ix = i % nx;
    GlobalOrdinal iy = (i - ix) / nx;
    return (ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void CountRowNNZ(const GlobalOrdinal i, LocalOrdinal& partial_nnz, const bool is_final) const {
    GlobalOrdinal center, left, right, bottom, top;
    bool isDirichlet;

    center = lclMap.getGlobalElement(i);
    GetNeighbours(center, left, right, bottom, top, isDirichlet);

    if (isDirichlet && keepBCs) {
      // Dirichlet unknown we want to keep
      Galeri_processEntry(center);
    } else {
      Galeri_processEntry(left);
      Galeri_processEntry(right);
      Galeri_processEntry(bottom);
      Galeri_processEntry(top);
      Galeri_processEntry(center);
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void EnterValues(const LocalOrdinal i, typename rowptr_type::value_type& entryPtr) const {
    GlobalOrdinal center, left, right, bottom, top;
    bool isDirichlet;

    center = lclMap.getGlobalElement(i);
    GetNeighbours(center, left, right, bottom, top, isDirichlet);

    impl_scalar_type offDiagonalSum           = zero;
    typename rowptr_type::value_type rowStart = entryPtr;
    if (isDirichlet && keepBCs) {
      // Dirichlet unknown we want to keep
      Galeri_enterValue(rowStart, center, one);
    } else {
      Galeri_enterValue(rowStart, left, b);
      Galeri_enterValue(rowStart, right, c);
      Galeri_enterValue(rowStart, bottom, d);
      Galeri_enterValue(rowStart, top, e);
      Galeri_enterValue(rowStart, center, (IsBoundary(center) && !isDirichlet) ? -offDiagonalSum : a);
    }
  }
};

template <class Scalar, class Map, bool keepBCs>
class Cross3DStencil {
  //    e
  //  b a c
  //    d
  // + f bottom and g top

#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<Scalar>;
#else
  using ATS                   = Kokkos::ArithTraits<Scalar>;
#endif
  using impl_scalar_type  = typename ATS::val_type;
  using LocalOrdinal      = typename Map::local_ordinal_type;
  using GlobalOrdinal     = typename Map::global_ordinal_type;
  using Node              = typename Map::node_type;
  using local_map_type    = typename Map::local_map_type;
  using local_matrix_type = typename ::Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;
  using rowptr_type       = typename local_matrix_type::row_map_type::non_const_type;
  using colidx_type       = typename local_matrix_type::index_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using exec_space        = typename Node::execution_space;
  using memory_space      = typename Node::memory_space;
  using hashmap_type      = typename Kokkos::UnorderedMap<GlobalOrdinal, void, exec_space>;

  GlobalOrdinal nx, ny, nz;
  impl_scalar_type a, b, c, d, e, f, g;
  DirBC DirichletBC;
  local_map_type lclMap;

#if KOKKOS_VERSION >= 40799
  const impl_scalar_type one = KokkosKernels::ArithTraits<impl_scalar_type>::one();
#else
  const impl_scalar_type one  = Kokkos::ArithTraits<impl_scalar_type>::one();
#endif
#if KOKKOS_VERSION >= 40799
  const impl_scalar_type zero = KokkosKernels::ArithTraits<impl_scalar_type>::zero();
#else
  const impl_scalar_type zero = Kokkos::ArithTraits<impl_scalar_type>::zero();
#endif
  const GlobalOrdinal INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();

 public:
  hashmap_type off_rank_indices;
  local_map_type lclColMap;
  colidx_type colidx;
  values_type values;

  Cross3DStencil(const Map& map, GlobalOrdinal nx_, GlobalOrdinal ny_, GlobalOrdinal nz_, Scalar a_, Scalar b_, Scalar c_, Scalar d_, Scalar e_, Scalar f_, Scalar g_, DirBC DirichletBC_)
    : nx(nx_)
    , ny(ny_)
    , nz(nz_)
    , a(a_)
    , b(b_)
    , c(c_)
    , d(d_)
    , e(e_)
    , f(f_)
    , g(g_)
    , DirichletBC(DirichletBC_) {
    lclMap = map.getLocalMap();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void GetNeighbours(const GlobalOrdinal i,
                     GlobalOrdinal& left, GlobalOrdinal& right,
                     GlobalOrdinal& front, GlobalOrdinal& back,
                     GlobalOrdinal& bottom, GlobalOrdinal& top, bool& isDirichlet) const {
    GetNeighboursCartesian3dKokkos(i, nx, ny, nz, left, right, front, back, bottom, top, INVALID);
    isDirichlet = (left == INVALID && (DirichletBC & DIR_LEFT)) ||
                  (right == INVALID && (DirichletBC & DIR_RIGHT)) ||
                  (bottom == INVALID && (DirichletBC & DIR_BOTTOM)) ||
                  (top == INVALID && (DirichletBC & DIR_TOP)) ||
                  (front == INVALID && (DirichletBC & DIR_FRONT)) ||
                  (back == INVALID && (DirichletBC & DIR_BACK));
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool IsBoundary(const GlobalOrdinal i) const {
    GlobalOrdinal ix  = i % nx;
    GlobalOrdinal ixy = i % (nx * ny);
    GlobalOrdinal iy  = (ixy - ix) / nx;
    GlobalOrdinal iz  = (i - ixy) / (nx * ny);
    return (ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 || iz == 0 || iz == nz - 1);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void CountRowNNZ(const GlobalOrdinal i, LocalOrdinal& partial_nnz, const bool is_final) const {
    GlobalOrdinal center, left, right, front, back, bottom, top;
    bool isDirichlet;

    center = lclMap.getGlobalElement(i);
    GetNeighbours(center, left, right, front, back, bottom, top, isDirichlet);

    if (isDirichlet && keepBCs) {
      // Dirichlet unknown we want to keep
      Galeri_processEntry(center);
    } else {
      Galeri_processEntry(left);
      Galeri_processEntry(right);
      Galeri_processEntry(front);
      Galeri_processEntry(back);
      Galeri_processEntry(bottom);
      Galeri_processEntry(top);
      Galeri_processEntry(center);
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void EnterValues(const LocalOrdinal i, typename rowptr_type::value_type& entryPtr) const {
    GlobalOrdinal center, left, right, front, back, bottom, top;
    bool isDirichlet;

    center = lclMap.getGlobalElement(i);
    GetNeighbours(center, left, right, front, back, bottom, top, isDirichlet);

    impl_scalar_type offDiagonalSum           = zero;
    typename rowptr_type::value_type rowStart = entryPtr;
    if (isDirichlet && keepBCs) {
      // Dirichlet unknown we want to keep
      Galeri_enterValue(rowStart, center, one);
    } else {
      Galeri_enterValue(rowStart, left, b);
      Galeri_enterValue(rowStart, right, c);
      Galeri_enterValue(rowStart, front, d);
      Galeri_enterValue(rowStart, back, e);
      Galeri_enterValue(rowStart, bottom, f);
      Galeri_enterValue(rowStart, top, g);
      Galeri_enterValue(rowStart, center, (IsBoundary(center) && !isDirichlet) ? -offDiagonalSum : a);
    }
  }
};

template <class Scalar, class Map, bool keepBCs>
class Brick3DStencil {
  // lower plane
  //   e  d  e
  //   d  b  d
  //   e  d  e

  // middle plane
  //   c  b  c
  //   b  a  b
  //   c  b  c

  // upper plane
  //   e  d  e
  //   d  b  d
  //   e  d  e

#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<Scalar>;
#else
  using ATS                   = Kokkos::ArithTraits<Scalar>;
#endif
  using impl_scalar_type  = typename ATS::val_type;
  using LocalOrdinal      = typename Map::local_ordinal_type;
  using GlobalOrdinal     = typename Map::global_ordinal_type;
  using Node              = typename Map::node_type;
  using local_map_type    = typename Map::local_map_type;
  using local_matrix_type = typename ::Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;
  using rowptr_type       = typename local_matrix_type::row_map_type::non_const_type;
  using colidx_type       = typename local_matrix_type::index_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using exec_space        = typename Node::execution_space;
  using memory_space      = typename Node::memory_space;
  using hashmap_type      = typename Kokkos::UnorderedMap<GlobalOrdinal, void, exec_space>;
  using StencilIndices    = GlobalOrdinal[3][3][3];

  GlobalOrdinal nx, ny, nz;
  impl_scalar_type a, b, c, d, e, f, g;
  DirBC DirichletBC;
  local_map_type lclMap;

#if KOKKOS_VERSION >= 40799
  const impl_scalar_type one = KokkosKernels::ArithTraits<impl_scalar_type>::one();
#else
  const impl_scalar_type one  = Kokkos::ArithTraits<impl_scalar_type>::one();
#endif
#if KOKKOS_VERSION >= 40799
  const impl_scalar_type zero = KokkosKernels::ArithTraits<impl_scalar_type>::zero();
#else
  const impl_scalar_type zero = Kokkos::ArithTraits<impl_scalar_type>::zero();
#endif
  GlobalOrdinal INVALID;

  Kokkos::View<impl_scalar_type[3][3][3], memory_space> stencil_entries;

  const size_t LOW  = 0;
  const size_t MID  = 1;
  const size_t HIGH = 2;

 public:
  hashmap_type off_rank_indices;
  local_map_type lclColMap;
  colidx_type colidx;
  values_type values;

  Brick3DStencil(const Map& map, GlobalOrdinal nx_, GlobalOrdinal ny_, GlobalOrdinal nz_, Scalar a_, Scalar b_, Scalar c_, Scalar d_, Scalar e_, Scalar f_, Scalar g_, DirBC DirichletBC_)
    : nx(nx_)
    , ny(ny_)
    , nz(nz_)
    , a(a_)
    , b(b_)
    , c(c_)
    , d(d_)
    , e(e_)
    , f(f_)
    , g(g_)
    , DirichletBC(DirichletBC_) {
    lclMap  = map.getLocalMap();
    INVALID = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();

    stencil_entries        = Kokkos::View<impl_scalar_type[3][3][3], memory_space>("stencil_entries");
    auto stencil_entries_h = Kokkos::create_mirror_view(stencil_entries);

    // lower plane
    stencil_entries_h(LOW, LOW, LOW)   = e;
    stencil_entries_h(LOW, LOW, MID)   = d;
    stencil_entries_h(LOW, LOW, HIGH)  = e;
    stencil_entries_h(LOW, MID, LOW)   = d;
    stencil_entries_h(LOW, MID, MID)   = b;
    stencil_entries_h(LOW, MID, HIGH)  = d;
    stencil_entries_h(LOW, HIGH, LOW)  = e;
    stencil_entries_h(LOW, HIGH, MID)  = d;
    stencil_entries_h(LOW, HIGH, HIGH) = e;

    // middle plane
    stencil_entries_h(MID, LOW, LOW)   = c;
    stencil_entries_h(MID, LOW, MID)   = b;
    stencil_entries_h(MID, LOW, HIGH)  = c;
    stencil_entries_h(MID, MID, LOW)   = b;
    stencil_entries_h(MID, MID, MID)   = a;
    stencil_entries_h(MID, MID, HIGH)  = b;
    stencil_entries_h(MID, HIGH, LOW)  = c;
    stencil_entries_h(MID, HIGH, MID)  = b;
    stencil_entries_h(MID, HIGH, HIGH) = c;

    // upper plane
    stencil_entries_h(HIGH, LOW, LOW)   = e;
    stencil_entries_h(HIGH, LOW, MID)   = d;
    stencil_entries_h(HIGH, LOW, HIGH)  = e;
    stencil_entries_h(HIGH, MID, LOW)   = d;
    stencil_entries_h(HIGH, MID, MID)   = b;
    stencil_entries_h(HIGH, MID, HIGH)  = d;
    stencil_entries_h(HIGH, HIGH, LOW)  = e;
    stencil_entries_h(HIGH, HIGH, MID)  = d;
    stencil_entries_h(HIGH, HIGH, HIGH) = e;

    Kokkos::deep_copy(stencil_entries, stencil_entries_h);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void GetNeighbours(const GlobalOrdinal i,
                     StencilIndices stencil_indices,
                     bool& isDirichlet) const {
    GlobalOrdinal& below = stencil_indices[LOW][MID][MID];
    GlobalOrdinal& front = stencil_indices[MID][LOW][MID];
    GlobalOrdinal& left  = stencil_indices[MID][MID][LOW];
    GlobalOrdinal& right = stencil_indices[MID][MID][HIGH];
    GlobalOrdinal& back  = stencil_indices[MID][HIGH][MID];
    GlobalOrdinal& above = stencil_indices[HIGH][MID][MID];

    GetNeighboursCartesian3dKokkos(i, nx, ny, nz, left, right, front, back, below, above, INVALID);

    stencil_indices[LOW][LOW][LOW]   = below - nx - 1;
    stencil_indices[LOW][LOW][MID]   = below - nx;
    stencil_indices[LOW][LOW][HIGH]  = below - nx + 1;
    stencil_indices[LOW][MID][LOW]   = below - 1;
    stencil_indices[LOW][MID][MID]   = below;
    stencil_indices[LOW][MID][HIGH]  = below + 1;
    stencil_indices[LOW][HIGH][LOW]  = below + nx - 1;
    stencil_indices[LOW][HIGH][MID]  = below + nx;
    stencil_indices[LOW][HIGH][HIGH] = below + nx + 1;

    stencil_indices[MID][LOW][LOW]   = front - 1;
    stencil_indices[MID][LOW][MID]   = front;
    stencil_indices[MID][LOW][HIGH]  = front + 1;
    stencil_indices[MID][MID][LOW]   = left;
    stencil_indices[MID][MID][MID]   = i;
    stencil_indices[MID][MID][HIGH]  = right;
    stencil_indices[MID][HIGH][LOW]  = back - 1;
    stencil_indices[MID][HIGH][MID]  = back;
    stencil_indices[MID][HIGH][HIGH] = back + 1;

    stencil_indices[HIGH][LOW][LOW]   = above - nx - 1;
    stencil_indices[HIGH][LOW][MID]   = above - nx;
    stencil_indices[HIGH][LOW][HIGH]  = above - nx + 1;
    stencil_indices[HIGH][MID][LOW]   = above - 1;
    stencil_indices[HIGH][MID][MID]   = above;
    stencil_indices[HIGH][MID][HIGH]  = above + 1;
    stencil_indices[HIGH][HIGH][LOW]  = above + nx - 1;
    stencil_indices[HIGH][HIGH][MID]  = above + nx;
    stencil_indices[HIGH][HIGH][HIGH] = above + nx + 1;

    if (left == INVALID) {
      for (size_t k0 = 0; k0 < 3; ++k0)
        for (size_t k1 = 0; k1 < 3; ++k1)
          stencil_indices[k0][k1][LOW] = INVALID;
    }
    if (right == INVALID) {
      for (size_t k0 = 0; k0 < 3; ++k0)
        for (size_t k1 = 0; k1 < 3; ++k1)
          stencil_indices[k0][k1][HIGH] = INVALID;
    }
    if (front == INVALID) {
      for (size_t k0 = 0; k0 < 3; ++k0)
        for (size_t k2 = 0; k2 < 3; ++k2)
          stencil_indices[k0][LOW][k2] = INVALID;
    }
    if (back == INVALID) {
      for (size_t k0 = 0; k0 < 3; ++k0)
        for (size_t k2 = 0; k2 < 3; ++k2)
          stencil_indices[k0][HIGH][k2] = INVALID;
    }
    if (below == INVALID) {
      for (size_t k1 = 0; k1 < 3; ++k1)
        for (size_t k2 = 0; k2 < 3; ++k2)
          stencil_indices[LOW][k1][k2] = INVALID;
    }
    if (above == INVALID) {
      for (size_t k1 = 0; k1 < 3; ++k1)
        for (size_t k2 = 0; k2 < 3; ++k2)
          stencil_indices[HIGH][k1][k2] = INVALID;
    }

    isDirichlet = (left == INVALID && (DirichletBC & DIR_LEFT)) ||
                  (right == INVALID && (DirichletBC & DIR_RIGHT)) ||
                  (below == INVALID && (DirichletBC & DIR_BOTTOM)) ||
                  (above == INVALID && (DirichletBC & DIR_TOP)) ||
                  (front == INVALID && (DirichletBC & DIR_FRONT)) ||
                  (back == INVALID && (DirichletBC & DIR_BACK));
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool IsBoundary(const GlobalOrdinal i) const {
    GlobalOrdinal ix  = i % nx;
    GlobalOrdinal ixy = i % (nx * ny);
    GlobalOrdinal iy  = (ixy - ix) / nx;
    GlobalOrdinal iz  = (i - ixy) / (nx * ny);
    return (ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 || iz == 0 || iz == nz - 1);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void CountRowNNZ(const GlobalOrdinal i, LocalOrdinal& partial_nnz, const bool is_final) const {
    GlobalOrdinal center;
    StencilIndices stencil_indices;
    bool isDirichlet;

    center = lclMap.getGlobalElement(i);
    GetNeighbours(center, stencil_indices, isDirichlet);

    if (isDirichlet && keepBCs) {
      // Dirichlet unknown we want to keep
      Galeri_processEntry(center);
    } else {
      for (size_t k0 = 0; k0 < 3; ++k0)
        for (size_t k1 = 0; k1 < 3; ++k1)
          for (size_t k2 = 0; k2 < 3; ++k2)
            Galeri_processEntry(stencil_indices[k0][k1][k2]);
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void EnterValues(const LocalOrdinal i, typename rowptr_type::value_type& entryPtr) const {
    GlobalOrdinal center;
    StencilIndices stencil_indices;
    bool isDirichlet;

    center = lclMap.getGlobalElement(i);
    GetNeighbours(center, stencil_indices, isDirichlet);

    impl_scalar_type offDiagonalSum           = zero;
    typename rowptr_type::value_type rowStart = entryPtr;
    if (isDirichlet && keepBCs) {
      // Dirichlet unknown we want to keep
      Galeri_enterValue(rowStart, center, one);
    } else {
      for (size_t k0 = 0; k0 < 3; ++k0)
        for (size_t k1 = 0; k1 < 3; ++k1)
          for (size_t k2 = 0; k2 < 3; ++k2)
            if ((k0 != MID) || (k1 != MID) || (k2 != MID))
              Galeri_enterValue(rowStart, stencil_indices[k0][k1][k2], stencil_entries(k0, k1, k2));
      Galeri_enterValue(rowStart, center, (IsBoundary(center) && !isDirichlet) ? -offDiagonalSum : stencil_entries(MID, MID, MID));
    }
  }
};

#undef Galeri_processEntry
#undef Galeri_enterValue
#endif

template <typename GlobalOrdinal, typename Scalar>
void Fill9PointStencil(const GlobalOrdinal center,
                       std::vector<Scalar>& Values, std::vector<GlobalOrdinal>& Indices, size_t& numEntries,
                       const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                       const Scalar b, const Scalar c, const Scalar d, const Scalar e,
                       const Scalar z1, const Scalar z2, const Scalar z3, const Scalar z4,
                       GlobalOrdinal left = -2, GlobalOrdinal right = -2,
                       GlobalOrdinal lower = -2, GlobalOrdinal upper = -2);
/* end of prototypes */

/* ****************************************************************************************************** *
 *    (Scaled) Identity
 * ****************************************************************************************************** */
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
Identity(const Teuchos::RCP<const Map>& map, const Scalar a) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Teuchos::RCP<Matrix> mtx = MatrixTraits<Map, Matrix>::Build(map, 1);

  LocalOrdinal NumMyElements                               = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();

  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Scaled Identity Generation")));

    for (LocalOrdinal i = 0; i < NumMyElements; i++)
      mtx->insertGlobalValues(MyGlobalElements[i],
                              Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                              Teuchos::tuple<Scalar>(a));
  }
  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Scaled Identity fillComplete")));
    mtx->fillComplete();
  }
  return mtx;
}

/* ****************************************************************************************************** *
 *    Laplace 1D
 * ****************************************************************************************************** */
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
TriDiag(const Teuchos::RCP<const Map>& map,
        const GlobalOrdinal nx,  // note: nx unused
        const Scalar a, const Scalar b, const Scalar c) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Teuchos::RCP<Matrix> mtx = MatrixTraits<Map, Matrix>::Build(map, 3);

  LocalOrdinal NumMyElements                               = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
  GlobalOrdinal indexBase                                  = map->getIndexBase();

  Teuchos::RCP<const Teuchos::Comm<int>> comm = map->getComm();

  GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();

  GlobalOrdinal NumEntries;
  LocalOrdinal nnz = 2;
  std::vector<Scalar> Values(nnz);
  std::vector<GlobalOrdinal> Indices(nnz);

  comm->barrier();

  // c a b
  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 1D Generation")));

    for (LocalOrdinal i = 0; i < NumMyElements; i++) {
      if (MyGlobalElements[i] == indexBase) {
        // off-diagonal for first row
        Indices[0] = 1 + indexBase;
        Values[0]  = c;
        NumEntries = 1;

      } else if (MyGlobalElements[i] == NumGlobalElements + indexBase - 1) {
        // off-diagonal for last row
        Indices[0] = NumGlobalElements - 2 + indexBase;
        Values[0]  = b;
        NumEntries = 1;

      } else {
        // off-diagonal for internal row
        Indices[0] = MyGlobalElements[i] - 1;
        Values[0]  = b;
        Indices[1] = MyGlobalElements[i] + 1;
        Values[1]  = c;
        NumEntries = 2;
      }

      // put the off-diagonal entries
      // Xpetra wants ArrayViews (sigh)
      Teuchos::ArrayView<Scalar> av(&Values[0], NumEntries);
      Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], NumEntries);
      mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

      // Put in the diagonal entry
      mtx->insertGlobalValues(MyGlobalElements[i],
                              Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                              Teuchos::tuple<Scalar>(a));
    }
  }

  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 1D fillComplete")));
    mtx->fillComplete();
  }

  return mtx;
}

/* ****************************************************************************************************** *
 *    Laplace 2D
 * ****************************************************************************************************** */
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
Cross2D(const Teuchos::RCP<const Map>& map,
        const GlobalOrdinal nx, const GlobalOrdinal ny,
        const Scalar a, const Scalar b, const Scalar c,
        const Scalar d, const Scalar e,
        const DirBC DirichletBC = 0, const bool keepBCs = false) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  LocalOrdinal nnz = 5;

  RCP<Matrix> mtx = MatrixTraits<Map, Matrix>::Build(map, nnz);

  LocalOrdinal numMyElements = map->getLocalNumElements();
  GlobalOrdinal indexBase    = map->getIndexBase();

  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

  GlobalOrdinal center, left, right, lower, upper;
  std::vector<Scalar> vals(nnz);
  std::vector<GlobalOrdinal> inds(nnz);

  //    e
  //  b a c
  //    d
  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 2D Generation")));
    for (LocalOrdinal i = 0; i < numMyElements; ++i) {
      size_t n = 0;

      center = myGlobalElements[i] - indexBase;
      GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper);

      bool isDirichlet = (left == -1 && (DirichletBC & DIR_LEFT)) ||
                         (right == -1 && (DirichletBC & DIR_RIGHT)) ||
                         (lower == -1 && (DirichletBC & DIR_BOTTOM)) ||
                         (upper == -1 && (DirichletBC & DIR_TOP));

      if (isDirichlet && keepBCs) {
        // Dirichlet unknown we want to keep
        inds[n]   = center;
        vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

      } else {
        // The Neumann b.c. are treated in a sane way. The Dirichlet b.c., however, are treated
        // insane when the option keepBCs=false. Speicifically, in this case we don't want to keep
        // Dirichlet b.c., but that would result in inconsistency between the map and the number of
        // degrees of freedom, plus the problem with GIDs. Therefore, we virtually expand domain by
        // one node in the direction of the Dirichlet b.c., and then assume that that node was
        // not kept. But we use an old GIDs. So yes, that's weird.

        if (left != -1) {
          inds[n]   = left;
          vals[n++] = b;
        }
        if (right != -1) {
          inds[n]   = right;
          vals[n++] = c;
        }
        if (lower != -1) {
          inds[n]   = lower;
          vals[n++] = d;
        }
        if (upper != -1) {
          inds[n]   = upper;
          vals[n++] = e;
        }

        // diagonal
        Scalar z = a;
        if (IsBoundary2d(center, nx, ny) && !isDirichlet) {
          // Neumann boundary unknown (diagonal = sum of all offdiagonal)
          z = Teuchos::ScalarTraits<Scalar>::zero();
          for (size_t j = 0; j < n; j++)
            z -= vals[j];
        }
        inds[n]   = center;
        vals[n++] = z;
      }

      for (size_t j = 0; j < n; j++)
        inds[j] += indexBase;

      Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
      Teuchos::ArrayView<Scalar> av(&vals[0], n);
      mtx->insertGlobalValues(myGlobalElements[i], iv, av);
    }
  }

  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 2D FillComplete")));
    mtx->fillComplete();
  }

  return mtx;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/* ****************************************************************************************************** *
 *    Recirc 2D
 * ****************************************************************************************************** */
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename Vector>
Teuchos::RCP<Matrix>
Cross2D(const Teuchos::RCP<const Map>& map, const GlobalOrdinal nx, const GlobalOrdinal ny,
        const Vector& A, const Vector& B, const Vector& C,
        const Vector& D, const Vector& E) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  LocalOrdinal nnz = 5;
  RCP<Matrix> mtx  = MatrixTraits<Map, Matrix>::Build(map, nnz);

  LocalOrdinal numMyElements                               = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

  GlobalOrdinal left, right, lower, upper;
  std::vector<Scalar> Values(nnz);
  std::vector<GlobalOrdinal> Indices(nnz);

  auto Adata = A->getData(0);
  auto Bdata = B->getData(0);
  auto Cdata = C->getData(0);
  auto Ddata = D->getData(0);
  auto Edata = E->getData(0);

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Recirc 2D Generation")));
  for (int i = 0; i < numMyElements; i++) {
    int NumEntries = 0;
    GetNeighboursCartesian2d(myGlobalElements[i], nx, ny, left, right, lower, upper);

    if (left != -1) {
      Indices[NumEntries] = left;
      Values[NumEntries]  = Bdata[i];
      ++NumEntries;
    }
    if (right != -1) {
      Indices[NumEntries] = right;
      Values[NumEntries]  = Cdata[i];
      ++NumEntries;
    }
    if (lower != -1) {
      Indices[NumEntries] = lower;
      Values[NumEntries]  = Ddata[i];
      ++NumEntries;
    }
    if (upper != -1) {
      Indices[NumEntries] = upper;
      Values[NumEntries]  = Edata[i];
      ++NumEntries;
    }

    Indices[NumEntries] = myGlobalElements[i];
    Values[NumEntries]  = Adata[i];
    ++NumEntries;

    Teuchos::ArrayView<GlobalOrdinal> inds(Indices.data(), NumEntries);
    Teuchos::ArrayView<Scalar> vals(Values.data(), NumEntries);
    mtx->insertGlobalValues(myGlobalElements[i], inds, vals);
  }
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Recirc2D FillComplete")));
  mtx->fillComplete();
  tm = Teuchos::null;

  return (mtx);
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

/* ****************************************************************************************************** *
 *    Star2D
 * ****************************************************************************************************** */
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
Star2D(const Teuchos::RCP<const Map>& map,
       const GlobalOrdinal nx, const GlobalOrdinal ny,
       const Scalar a, const Scalar b, const Scalar c,
       const Scalar d, const Scalar e,
       const Scalar z1, const Scalar z2,
       const Scalar z3, const Scalar z4,
       const DirBC DirichletBC = 0, const bool keepBCs = false) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  LocalOrdinal nnz = 9;

  Teuchos::RCP<Matrix> mtx = MatrixTraits<Map, Matrix>::Build(map, nnz);

  LocalOrdinal numMyElements = map->getLocalNumElements();
  GlobalOrdinal indexBase    = map->getIndexBase();

  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

  GlobalOrdinal center, left, right, lower, upper;
  std::vector<Scalar> vals(nnz);
  std::vector<GlobalOrdinal> inds(nnz);

  //  z3  e  z4
  //   b  a  c
  //  z1  d  z2
  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Star2D generation")));

    for (LocalOrdinal i = 0; i < numMyElements; ++i) {
      size_t n = 0;

      center = myGlobalElements[i] - indexBase;
      GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper);

      bool isDirichlet = (left == -1 && (DirichletBC & DIR_LEFT)) ||
                         (right == -1 && (DirichletBC & DIR_RIGHT)) ||
                         (lower == -1 && (DirichletBC & DIR_BOTTOM)) ||
                         (upper == -1 && (DirichletBC & DIR_TOP));

      if (isDirichlet && keepBCs) {
        // Dirichlet unknown we want to keep
        inds[n]   = center;
        vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

      } else {
        // See comments about weird in Cross2D
        if (left != -1) {
          inds[n]   = left;
          vals[n++] = b;
        }
        if (right != -1) {
          inds[n]   = right;
          vals[n++] = c;
        }
        if (lower != -1) {
          inds[n]   = lower;
          vals[n++] = d;
        }
        if (upper != -1) {
          inds[n]   = upper;
          vals[n++] = e;
        }
        if (left != -1 && lower != -1) {
          inds[n]   = lower - 1;
          vals[n++] = z1;
        }
        if (right != -1 && lower != -1) {
          inds[n]   = lower + 1;
          vals[n++] = z2;
        }
        if (left != -1 && upper != -1) {
          inds[n]   = upper - 1;
          vals[n++] = z3;
        }
        if (right != -1 && upper != -1) {
          inds[n]   = upper + 1;
          vals[n++] = z4;
        }

        // diagonal
        Scalar z = a;
        if (IsBoundary2d(center, nx, ny) && !isDirichlet) {
          // Neumann boundary unknown (diagonal = sum of all offdiagonal)
          z = Teuchos::ScalarTraits<Scalar>::zero();
          for (size_t j = 0; j < n; j++)
            z -= vals[j];
        }
        inds[n]   = center + indexBase;
        vals[n++] = z;
      }

      for (size_t j = 0; j < n; j++)
        inds[j] += indexBase;

      Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
      Teuchos::ArrayView<Scalar> av(&vals[0], n);
      mtx->insertGlobalValues(myGlobalElements[i], iv, av);
    }
  }

  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: star2D fillComplete")));
    mtx->fillComplete();
  }

  return mtx;
}

/* ****************************************************************************************************** *
 *    BigStar2D (2D Biharmonic operator, for example)
 * ****************************************************************************************************** */
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
BigStar2D(const Teuchos::RCP<const Map>& map,
          const GlobalOrdinal nx, const GlobalOrdinal ny,
          const Scalar a, const Scalar b, const Scalar c,
          const Scalar d, const Scalar e,
          const Scalar z1, const Scalar z2,
          const Scalar z3, const Scalar z4,
          const Scalar bb, const Scalar cc, const Scalar dd, const Scalar ee,
          const DirBC DirichletBC = 0, const bool keepBCs = false) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  LocalOrdinal nnz = 13;

  Teuchos::RCP<Matrix> mtx = MatrixTraits<Map, Matrix>::Build(map, nnz);

  LocalOrdinal numMyElements = map->getLocalNumElements();
  GlobalOrdinal indexBase    = map->getIndexBase();

  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

  GlobalOrdinal center, left, right, lower, upper, left2, right2, lower2, upper2;
  std::vector<Scalar> vals(nnz);
  std::vector<GlobalOrdinal> inds(nnz);

  //        ee
  //    z3  e  z4
  // bb  b  a  c  cc
  //    z1  d  z2
  //        dd
  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: bigStar2D Generation")));

    for (LocalOrdinal i = 0; i < numMyElements; ++i) {
      size_t n = 0;

      center = myGlobalElements[i] - indexBase;
      GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper, left2, right2, lower2, upper2);

      bool isDirichlet = (left == -1 && (DirichletBC & DIR_LEFT)) ||
                         (right == -1 && (DirichletBC & DIR_RIGHT)) ||
                         (lower == -1 && (DirichletBC & DIR_BOTTOM)) ||
                         (upper == -1 && (DirichletBC & DIR_TOP));

      if (isDirichlet && keepBCs) {
        // Dirichlet unknown we want to keep
        inds[n]   = center;
        vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

      } else {
        // See comments about weird in Cross2D
        if (left != -1) {
          inds[n]   = left;
          vals[n++] = b;
        }
        if (right != -1) {
          inds[n]   = right;
          vals[n++] = c;
        }
        if (lower != -1) {
          inds[n]   = lower;
          vals[n++] = d;
        }
        if (upper != -1) {
          inds[n]   = upper;
          vals[n++] = e;
        }
        if (left != -1 && lower != -1) {
          inds[n]   = lower - 1;
          vals[n++] = z1;
        }
        if (right != -1 && lower != -1) {
          inds[n]   = lower + 1;
          vals[n++] = z2;
        }
        if (left != -1 && upper != -1) {
          inds[n]   = upper - 1;
          vals[n++] = z3;
        }
        if (right != -1 && upper != -1) {
          inds[n]   = upper + 1;
          vals[n++] = z4;
        }
        if (left2 != -1) {
          inds[n]   = left2;
          vals[n++] = bb;
        }
        if (right2 != -1) {
          inds[n]   = right2;
          vals[n++] = cc;
        }
        if (lower2 != -1) {
          inds[n]   = lower2;
          vals[n++] = dd;
        }
        if (upper2 != -1) {
          inds[n]   = upper2;
          vals[n++] = ee;
        }

        // diagonal
        Scalar z = a;
        if (IsBoundary2d(center, nx, ny) && !isDirichlet) {
          // Neumann boundary unknown (diagonal = sum of all offdiagonal)
          z = Teuchos::ScalarTraits<Scalar>::zero();
          for (size_t j = 0; j < n; j++)
            z -= vals[j];
        }
        inds[n]   = center + indexBase;
        vals[n++] = z;
      }

      for (size_t j = 0; j < n; j++)
        inds[j] += indexBase;

      Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
      Teuchos::ArrayView<Scalar> av(&vals[0], n);
      mtx->insertGlobalValues(myGlobalElements[i], iv, av);
    }
  }

  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: bigStar2D fillComplete")));

    mtx->fillComplete();
  }

  return mtx;
}

/* ****************************************************************************************************** *
 *    Laplace 3D
 * ****************************************************************************************************** */
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
Cross3D(const Teuchos::RCP<const Map>& map,
        const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
        const Scalar a, const Scalar b, const Scalar c,
        const Scalar d, const Scalar e,
        const Scalar f, const Scalar g,
        const DirBC DirichletBC = 0, const bool keepBCs = false) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  LocalOrdinal nnz = 7;

  Teuchos::RCP<Matrix> mtx = MatrixTraits<Map, Matrix>::Build(map, nnz);

  LocalOrdinal numMyElements = map->getLocalNumElements();
  GlobalOrdinal indexBase    = map->getIndexBase();

  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

  GlobalOrdinal center, left, right, bottom, top, front, back;
  std::vector<GlobalOrdinal> inds(nnz);
  std::vector<Scalar> vals(nnz);

  //    e
  //  b a c
  //    d
  // + f bottom and g top
  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 3D generation")));
    for (LocalOrdinal i = 0; i < numMyElements; ++i) {
      size_t n = 0;

      center = myGlobalElements[i] - indexBase;
      GetNeighboursCartesian3d(center, nx, ny, nz,
                               left, right, front, back, bottom, top);

      bool isDirichlet = (left == -1 && (DirichletBC & DIR_LEFT)) ||
                         (right == -1 && (DirichletBC & DIR_RIGHT)) ||
                         (bottom == -1 && (DirichletBC & DIR_BOTTOM)) ||
                         (top == -1 && (DirichletBC & DIR_TOP)) ||
                         (front == -1 && (DirichletBC & DIR_FRONT)) ||
                         (back == -1 && (DirichletBC & DIR_BACK));

      if (isDirichlet && keepBCs) {
        // Dirichlet unknown we want to keep
        inds[n]   = center;
        vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

      } else {
        // See comments about weird in Cross2D
        if (left != -1) {
          inds[n]   = left;
          vals[n++] = b;
        }
        if (right != -1) {
          inds[n]   = right;
          vals[n++] = c;
        }
        if (front != -1) {
          inds[n]   = front;
          vals[n++] = d;
        }
        if (back != -1) {
          inds[n]   = back;
          vals[n++] = e;
        }
        if (bottom != -1) {
          inds[n]   = bottom;
          vals[n++] = f;
        }
        if (top != -1) {
          inds[n]   = top;
          vals[n++] = g;
        }

        // diagonal
        Scalar z = a;
        if (IsBoundary3d(center, nx, ny, nz) && !isDirichlet) {
          // Neumann boundary unknown (diagonal = sum of all offdiagonal)
          z = Teuchos::ScalarTraits<Scalar>::zero();
          for (size_t j = 0; j < n; j++)
            z -= vals[j];
        }
        inds[n]   = center;
        vals[n++] = z;
      }

      for (size_t j = 0; j < n; j++)
        inds[j] += indexBase;

      Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
      Teuchos::ArrayView<Scalar> av(&vals[0], n);
      mtx->insertGlobalValues(myGlobalElements[i], iv, av);
    }
  }

  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 3D fillComplete")));
    mtx->fillComplete();
  }

  return mtx;
}

#if defined(HAVE_GALERI_KOKKOS) && defined(HAVE_GALERI_KOKKOSKERNELS)

template <typename Matrix, typename Map, typename Stencil>
Teuchos::RCP<Matrix>
StencilMatrixKokkos(const Teuchos::RCP<const Map>& map,
                    Stencil& stencil,
                    const std::string& matLabel) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using Scalar            = typename Matrix::scalar_type;
  using LocalOrdinal      = typename Map::local_ordinal_type;
  using GlobalOrdinal     = typename Map::global_ordinal_type;
  using Node              = typename Map::node_type;
  using local_matrix_type = typename ::Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;
  using rowptr_type       = typename local_matrix_type::row_map_type::non_const_type;
  using colidx_type       = typename local_matrix_type::index_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using tpetra_map        = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using exec_space        = typename Node::execution_space;
  using memory_space      = typename Node::memory_space;
  using range_type        = Kokkos::RangePolicy<LocalOrdinal, exec_space>;

  // We perform the assembly in 3 steps:
  //
  // 1) We count the number of entries per row and collect the off-rank entries in a hashmap.
  // 2) We construct the column map from the row map and the hashmap.
  // 3)_We allocate indices and values and enter the values.

  Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: " + matLabel + " generation")));

  ///////////////////////////
  // Step 1 - rowptr

  LocalOrdinal numMyElements = map->getLocalNumElements();
  // The hash map that collects all off-rank entries that we encounter.
  // We need to provide sufficient capacity so that no insertions fail.
  stencil.off_rank_indices = Kokkos::UnorderedMap<GlobalOrdinal, void, exec_space>(std::max(numMyElements, 100));

  auto lclMap = map->getLocalMap();
  // The row offsets of the matrix, initialized to zero.
  auto rowptr = rowptr_type("rowptr", numMyElements + 1);

  // We loop over the rows and over the stencil.
  // We count entries of row i using rowptr(i+2) and accumulate the results.
  //
  // E.g if we have a matrix with three rows with 2, 3 and 4 entries respectively, rowptr will look like
  // [0 0 2 5]
  // and myNNZ = 2+3+4 = 9.
  //
  // Note that this is not yet the correct rowptr for the matrix which looks like
  // [0 2 5 9]
  //
  // The reason why we are not constructing the correct rowptr right away is that we will use
  // rowptr(i+1) as the offset for filling the colidx and values views.
  LocalOrdinal myNNZ = 0;  // The number of nonzero entries on this rank.
  Kokkos::parallel_scan(
      "Galeri::" + matLabel + "::Graph", range_type(0, numMyElements),
      KOKKOS_LAMBDA(const LocalOrdinal i, LocalOrdinal& partial_nnz,
                    const bool is_final) {
        stencil.CountRowNNZ(i, partial_nnz, is_final);
        if (is_final) {
          if (i + 2 < numMyElements + 1)
            rowptr(i + 2) = partial_nnz;
        }
      },
      myNNZ);
  // We check that no hashmap insertions failed.
  // If this turns out to not be enough we can bump up the capacity.
  TEUCHOS_ASSERT(!stencil.off_rank_indices.failed_insert());

  ///////////////////////////
  // Step 2 - Column map

  auto columnMap_entries = Kokkos::View<GlobalOrdinal*, memory_space>("columnMap_entries", numMyElements + stencil.off_rank_indices.size());
  // Copy over on-rank entries from rowmap
  Kokkos::deep_copy(Kokkos::subview(columnMap_entries, Kokkos::make_pair(0, numMyElements)),
                    map->getMyGlobalIndicesDevice());
  // Set off-rank entries from hashmap.
  // This is the Kokkos recommended way of iterating, see
  // https://kokkos.org/kokkos-core-wiki/API/containers/Unordered-Map.html#iteration
  Kokkos::parallel_scan(
      "Galeri::" + matLabel + "::Graph", range_type(0, stencil.off_rank_indices.capacity()), KOKKOS_LAMBDA(const LocalOrdinal i, size_t& pos, const bool is_final) {
        if (stencil.off_rank_indices.valid_at(i)) {
          if (is_final) {
            auto key                               = stencil.off_rank_indices.key_at(i);
            columnMap_entries(numMyElements + pos) = key;
          }
          ++pos;
        }
      });
  // We no longer need the hashmap.
  stencil.off_rank_indices.clear();
  stencil.off_rank_indices.rehash(0);

  {
    // Sort the off-rank entries by PID
    auto remoteGIDs = Kokkos::subview(columnMap_entries, Kokkos::make_pair(numMyElements, columnMap_entries.extent_int(0)));
    auto remotePIDs = Kokkos::View<int*, memory_space>("remotePIDs", remoteGIDs.extent(0));
    {
      auto remoteGIDs_h = Kokkos::create_mirror_view(remoteGIDs);
      auto remotePIDs_h = Kokkos::create_mirror_view(remotePIDs);
      Kokkos::deep_copy(remoteGIDs_h, remoteGIDs);
      Teuchos::ArrayView<GlobalOrdinal> remoteGIDs_av(remoteGIDs_h.data(), remoteGIDs_h.extent(0));
      Teuchos::ArrayView<int> remotePIDs_av(remotePIDs_h.data(), remotePIDs_h.extent(0));
      auto ret = map->getRemoteIndexList(remoteGIDs_av, remotePIDs_av);
      if constexpr (std::is_same_v<Map, tpetra_map>) {
        TEUCHOS_ASSERT(ret == Tpetra::AllIDsPresent);
      } else {
        TEUCHOS_ASSERT(ret == ::Xpetra::AllIDsPresent);
      }
      Kokkos::deep_copy(remotePIDs, remotePIDs_h);
    }
    {
      typename decltype(remoteGIDs)::execution_space exec;
      Kokkos::Experimental::sort_by_key(exec, remotePIDs, remoteGIDs);
    }
  }

  RCP<Map> ghosted_map;
  if constexpr (std::is_same_v<Map, tpetra_map>) {
    ghosted_map = rcp(new Map(Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                              columnMap_entries,
                              map->getIndexBase(),
                              map->getComm()));
  } else {
    ghosted_map = ::Xpetra::MapFactory<typename Matrix::local_ordinal_type, typename Matrix::global_ordinal_type, typename Matrix::node_type>::Build(map->lib(),
                                                                                                                                                     Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                                                                                                                                     columnMap_entries,
                                                                                                                                                     map->getIndexBase(),
                                                                                                                                                     map->getComm());
  }

  ///////////////////////////
  // Step 3 - Fill

  stencil.lclColMap = ghosted_map->getLocalMap();
  stencil.colidx    = colidx_type("colidx", myNNZ);
  stencil.values    = values_type("values", myNNZ);

  // Loop over rows and stencil and fill the matrix.
  // To enter values in row i we use rowptr(i+1) as offset.
  //
  // Using the same example as above, after this step rowptr is
  // [0 2 5 9]
  Kokkos::parallel_for(
      "Galeri::" + matLabel + "::fill", range_type(0, numMyElements),
      KOKKOS_LAMBDA(const LocalOrdinal i) {
        stencil.EnterValues(i, rowptr(i + 1));
      });

  auto lclA = local_matrix_type(matLabel, numMyElements, ghosted_map->getLocalNumElements(), myNNZ, stencil.values, rowptr, stencil.colidx);
  auto mtx  = MatrixTraits<Map, Matrix>::Build(lclA, map, ghosted_map, map, map);

#ifdef HAVE_GALERI_DEBUG
  // Checks to run in debug mode:
  // - local CRS graph is sorted by row
  // - imports can be aliased to target MV
  // - remote entries of column map are sorted by PID so that transfers are faster

  TEUCHOS_ASSERT(::KokkosSparse::isCrsGraphSorted(lclA.graph.row_map,
                                                  lclA.graph.entries));

  if constexpr (std::is_same_v<Map, tpetra_map>) {
    auto import = mtx->getCrsGraph()->getImporter();
    if (!import.is_null()) {
      TEUCHOS_ASSERT(import->areRemoteLIDsContiguous());

      auto distor = import->getDistributor();
      TEUCHOS_ASSERT(distor.getPlan().getIndicesTo().is_null());
    }
  } else if constexpr (std::is_same_v<Matrix, ::Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>) {
    auto import = mtx->getCrsGraph()->getImporter();
    if (!import.is_null()) {
      auto xp_tp_import = Teuchos::rcp_dynamic_cast<const ::Xpetra::TpetraImport<LocalOrdinal, GlobalOrdinal, Node>>(import, true);
      auto tp_import    = xp_tp_import->getTpetra_Import();
      TEUCHOS_ASSERT(tp_import->areRemoteLIDsContiguous());

      auto distor = tp_import->getDistributor();
      TEUCHOS_ASSERT(distor.getPlan().getIndicesTo().is_null());
    }
  } else {
    auto import = mtx->getCrsMatrix()->getCrsGraph()->getImporter();
    if (!import.is_null()) {
      auto xp_tp_import                  = Teuchos::rcp_dynamic_cast<const ::Xpetra::TpetraImport<LocalOrdinal, GlobalOrdinal, Node>>(import, true);
      auto tp_import                     = xp_tp_import->getTpetra_Import();
      const bool areRemoteLIDsContiguous = tp_import->areRemoteLIDsContiguous();
      TEUCHOS_ASSERT(areRemoteLIDsContiguous);

      auto distor         = tp_import->getDistributor();
      const bool fastPath = distor.getPlan().getIndicesTo().is_null();
      TEUCHOS_ASSERT(fastPath);
    }
  }
#endif  // HAVE_GALERI_DEBUG

  return mtx;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
ScaledIdentityKokkos(const Teuchos::RCP<const Map>& map,
                     const GlobalOrdinal nx,
                     const Scalar a,
                     const std::string& label = "ScaledIdentity") {
  ScaledIdentityStencil<Scalar, Map> stencil(*map,
                                             nx,
                                             a);
  return StencilMatrixKokkos<Matrix>(map, stencil, label);
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
TriDiagKokkos(const Teuchos::RCP<const Map>& map,
              const GlobalOrdinal nx,
              const Scalar a, const Scalar b, const Scalar c,
              const DirBC DirichletBC  = 0,
              const bool keepBCs       = false,
              const std::string& label = "TriDiag") {
  if (keepBCs) {
    TriDiagStencil<Scalar, Map, true> stencil(*map,
                                              nx,
                                              a, b, c,
                                              DirichletBC);
    return StencilMatrixKokkos<Matrix>(map, stencil, label);
  } else {
    TriDiagStencil<Scalar, Map, false> stencil(*map,
                                               nx,
                                               a, b, c,
                                               DirichletBC);
    return StencilMatrixKokkos<Matrix>(map, stencil, label);
  }
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
Cross2DKokkos(const Teuchos::RCP<const Map>& map,
              const GlobalOrdinal nx, const GlobalOrdinal ny,
              const Scalar a, const Scalar b, const Scalar c,
              const Scalar d, const Scalar e,
              const DirBC DirichletBC  = 0,
              const bool keepBCs       = false,
              const std::string& label = "Cross2D") {
  if (keepBCs) {
    Cross2DStencil<Scalar, Map, true> stencil(*map,
                                              nx, ny,
                                              a, b, c, d, e,
                                              DirichletBC);
    return StencilMatrixKokkos<Matrix>(map, stencil, label);
  } else {
    Cross2DStencil<Scalar, Map, false> stencil(*map,
                                               nx, ny,
                                               a, b, c, d, e,
                                               DirichletBC);
    return StencilMatrixKokkos<Matrix>(map, stencil, label);
  }
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
Cross3DKokkos(const Teuchos::RCP<const Map>& map,
              const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
              const Scalar a, const Scalar b, const Scalar c,
              const Scalar d, const Scalar e,
              const Scalar f, const Scalar g,
              const DirBC DirichletBC  = 0,
              const bool keepBCs       = false,
              const std::string& label = "Cross3D") {
  if (keepBCs) {
    Cross3DStencil<Scalar, Map, true> stencil(*map,
                                              nx, ny, nz,
                                              a, b, c, d, e, f, g,
                                              DirichletBC);
    return StencilMatrixKokkos<Matrix>(map, stencil, label);
  } else {
    Cross3DStencil<Scalar, Map, false> stencil(*map,
                                               nx, ny, nz,
                                               a, b, c, d, e, f, g,
                                               DirichletBC);
    return StencilMatrixKokkos<Matrix>(map, stencil, label);
  }
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
Brick3DKokkos(const Teuchos::RCP<const Map>& map,
              const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
              const Scalar a, const Scalar b, const Scalar c,
              const Scalar d, const Scalar e,
              const Scalar f, const Scalar g,
              const DirBC DirichletBC  = 0,
              const bool keepBCs       = false,
              const std::string& label = "Brick3D") {
  if (keepBCs) {
    Brick3DStencil<Scalar, Map, true> stencil(*map,
                                              nx, ny, nz,
                                              a, b, c, d, e, f, g,
                                              DirichletBC);
    return StencilMatrixKokkos<Matrix>(map, stencil, label);
  } else {
    Brick3DStencil<Scalar, Map, false> stencil(*map,
                                               nx, ny, nz,
                                               a, b, c, d, e, f, g,
                                               DirichletBC);
    return StencilMatrixKokkos<Matrix>(map, stencil, label);
  }
}
#endif

/* ****************************************************************************************************** *
 *    3D 27-point stencil
 * ****************************************************************************************************** */
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
Teuchos::RCP<Matrix>
Brick3D(const Teuchos::RCP<const Map>& map,
        const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
        const Scalar a, const Scalar b, const Scalar c,
        const Scalar d, const Scalar e,
        const Scalar f, const Scalar g,
        const DirBC DirichletBC = 0, const bool keepBCs = false) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  LocalOrdinal nnz = 27;

  Teuchos::RCP<Matrix> mtx = MatrixTraits<Map, Matrix>::Build(map, nnz);

  LocalOrdinal numMyElements = map->getLocalNumElements();
  GlobalOrdinal indexBase    = map->getIndexBase();

  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

  GlobalOrdinal center, left, right, front, back, below, above;
  std::vector<Scalar> vals(nnz);
  std::vector<GlobalOrdinal> inds(nnz);

  // upper plane
  //   e  d  e
  //   d  b  d
  //   e  d  e

  // middle plane
  //   c  b  c
  //   b  a  b
  //   c  b  c

  // lower plane
  //   e  d  e
  //   d  b  d
  //   e  d  e

  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: 3D 27 point stencil generation")));
    for (LocalOrdinal i = 0; i < numMyElements; ++i) {
      size_t n = 0;

      center = myGlobalElements[i] - indexBase;
      GetNeighboursCartesian3d(center, nx, ny, nz,
                               left, right, front, back, below, above);

      bool isDirichlet = (left == -1 && (DirichletBC & DIR_LEFT)) ||
                         (right == -1 && (DirichletBC & DIR_RIGHT)) ||
                         (below == -1 && (DirichletBC & DIR_BOTTOM)) ||
                         (above == -1 && (DirichletBC & DIR_TOP)) ||
                         (front == -1 && (DirichletBC & DIR_FRONT)) ||
                         (back == -1 && (DirichletBC & DIR_BACK));

      if (isDirichlet && keepBCs) {
        // Dirichlet unknown we want to keep
        inds[n]   = center;
        vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

      } else {
        // See comments about weird in Cross2D

        // center plance (centered on center)
        Fill9PointStencil(center, vals, inds,
                          n, nx, ny, nz,
                          b, b, b, b, c, c, c, c,
                          left, right, front, back);
        // lower plane (centered on "below")
        if (below != -1) {
          inds[n]   = below;
          vals[n++] = b;
          Fill9PointStencil(below, vals, inds, n, nx, ny, nz,
                            d, d, d, d, e, e, e, e);
        }
        // upper plane (centered on "upper")
        if (above != -1) {
          inds[n]   = above;
          vals[n++] = b;
          Fill9PointStencil(above, vals, inds, n, nx, ny, nz,
                            d, d, d, d, e, e, e, e);
        }

        // diagonal
        Scalar z = a;
        if (IsBoundary3d(center, nx, ny, nz) && !isDirichlet) {
          // Neumann boundary unknown (diagonal = sum of all offdiagonal)
          z = Teuchos::ScalarTraits<Scalar>::zero();
          for (size_t j = 0; j < n; j++)
            z -= vals[j];
        }
        inds[n]   = center;
        vals[n++] = z;
      }

      for (size_t j = 0; j < n; j++)
        inds[j] += indexBase;

      Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
      Teuchos::ArrayView<Scalar> av(&vals[0], n);
      mtx->insertGlobalValues(myGlobalElements[i], iv, av);
    }
  }

  {
    Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: 3D 27 point stencil fillComplete")));
    mtx->fillComplete();
  }
  return mtx;
}

/* ****************************************************************************************************** *
 *    Utilities
 * ****************************************************************************************************** */

/* IsBoundary2d */
template <typename GlobalOrdinal>
bool IsBoundary2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny) {
  GlobalOrdinal ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  return (ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1);
}

/* IsBoundary3d */
template <typename GlobalOrdinal>
bool IsBoundary3d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz) {
  GlobalOrdinal ix, iy, ixy, iz;
  ix  = i % nx;
  ixy = i % (nx * ny);
  iy  = (ixy - ix) / nx;
  iz  = (i - ixy) / (nx * ny);

  return (ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 || iz == 0 || iz == nz - 1);
}

/* GetNeighboursCartesian2d */
template <typename GlobalOrdinal>
void GetNeighboursCartesian2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny,
                              GlobalOrdinal& left, GlobalOrdinal& right,
                              GlobalOrdinal& lower, GlobalOrdinal& upper) {
  GlobalOrdinal ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  if (ix == 0)
    left = -1;
  else
    left = i - 1;
  if (ix == nx - 1)
    right = -1;
  else
    right = i + 1;
  if (iy == 0)
    lower = -1;
  else
    lower = i - nx;
  if (iy == ny - 1)
    upper = -1;
  else
    upper = i + nx;
}

/* GetNeighboursCartesian2d */
template <typename GlobalOrdinal>
void GetNeighboursCartesian2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny,
                              GlobalOrdinal& left, GlobalOrdinal& right, GlobalOrdinal& lower, GlobalOrdinal& upper,
                              GlobalOrdinal& left2, GlobalOrdinal& right2, GlobalOrdinal& lower2, GlobalOrdinal& upper2) {
  GlobalOrdinal ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  if (ix == 0)
    left = -1;
  else
    left = i - 1;
  if (ix == nx - 1)
    right = -1;
  else
    right = i + 1;
  if (iy == 0)
    lower = -1;
  else
    lower = i - nx;
  if (iy == ny - 1)
    upper = -1;
  else
    upper = i + nx;

  if (ix <= 1)
    left2 = -1;
  else
    left2 = i - 2;
  if (ix >= nx - 2)
    right2 = -1;
  else
    right2 = i + 2;
  if (iy <= 1)
    lower2 = -1;
  else
    lower2 = i - 2 * nx;
  if (iy >= ny - 2)
    upper2 = -1;
  else
    upper2 = i + 2 * nx;
}

/* GetNeighboursCartesian3d */
template <typename GlobalOrdinal>
void GetNeighboursCartesian3d(const GlobalOrdinal i,
                              const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                              GlobalOrdinal& left, GlobalOrdinal& right,
                              GlobalOrdinal& front, GlobalOrdinal& back,
                              GlobalOrdinal& bottom, GlobalOrdinal& top) {
  GlobalOrdinal ixy, iz;
  ixy = i % (nx * ny);

  iz = (i - ixy) / (nx * ny);

  if (iz == 0)
    bottom = -1;
  else
    bottom = i - nx * ny;
  if (iz == nz - 1)
    top = -1;
  else
    top = i + nx * ny;

  GetNeighboursCartesian2d(ixy, nx, ny, left, right, front, back);

  if (left != -1) left += iz * (nx * ny);
  if (right != -1) right += iz * (nx * ny);
  if (front != -1) front += iz * (nx * ny);
  if (back != -1) back += iz * (nx * ny);
}

/* Fill9PointStencil */
template <typename GlobalOrdinal, typename Scalar>
void Fill9PointStencil(const GlobalOrdinal center,
                       std::vector<Scalar>& vals, std::vector<GlobalOrdinal>& inds, size_t& n,
                       const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                       const Scalar b, const Scalar c, const Scalar d, const Scalar e,
                       const Scalar z1, const Scalar z2, const Scalar z3, const Scalar z4,
                       GlobalOrdinal left, GlobalOrdinal right,
                       GlobalOrdinal lower, GlobalOrdinal upper) {
  //  z3  e  z4
  //   b  .  c
  //  z1  d  z2
  GlobalOrdinal below, above;
  if (left == -2)
    GetNeighboursCartesian3d(center, nx, ny, nz, left, right, lower, upper, below, above);

  if (left != -1) {
    inds[n]   = left;
    vals[n++] = b;
  }
  if (right != -1) {
    inds[n]   = right;
    vals[n++] = c;
  }
  if (lower != -1) {
    inds[n]   = lower;
    vals[n++] = d;
  }
  if (upper != -1) {
    inds[n]   = upper;
    vals[n++] = e;
  }
  if (left != -1 && lower != -1) {
    inds[n]   = lower - 1;
    vals[n++] = z1;
  }
  if (right != -1 && lower != -1) {
    inds[n]   = lower + 1;
    vals[n++] = z2;
  }
  if (left != -1 && upper != -1) {
    inds[n]   = upper - 1;
    vals[n++] = z3;
  }
  if (right != -1 && upper != -1) {
    inds[n]   = upper + 1;
    vals[n++] = z4;
  }
}

}  // namespace Xpetra

}  // namespace Galeri

#endif  // ifndef GALERI_XPETRAMATRIXTYPES_HPP
