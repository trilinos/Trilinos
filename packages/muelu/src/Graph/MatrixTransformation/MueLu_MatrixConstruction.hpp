// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MATRIXCONSTRUCTION_HPP
#define MUELU_MATRIXCONSTRUCTION_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"

#include "MueLu_DroppingCommon.hpp"

#ifdef MUELU_COALESCE_DROP_DEBUG
// For demangling function names
#include <cxxabi.h>
#endif

namespace MueLu::MatrixConstruction {

template <class local_matrix_type,
          class functor_type_0 = int,
          class functor_type_1 = int,
          class functor_type_2 = int,
          class functor_type_3 = int,
          class functor_type_4 = int,
          class functor_type_5 = int,
          class functor_type_6 = int,
          class functor_type_7 = int,
          class functor_type_8 = int,
          class functor_type_9 = int>
class PointwiseCountingFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  using rowptr_type = typename local_matrix_type::row_map_type::non_const_type;

  local_matrix_type A;
  results_view results;
  rowptr_type rowptr;
  functor_type_0 functor0;
  functor_type_1 functor1;
  functor_type_2 functor2;
  functor_type_3 functor3;
  functor_type_4 functor4;
  functor_type_5 functor5;
  functor_type_6 functor6;
  functor_type_7 functor7;
  functor_type_8 functor8;
  functor_type_9 functor9;

  std::vector<std::string> functorNames;

  void init() {
#ifdef MUELU_COALESCE_DROP_DEBUG
    functorNames = std::vector<std::string>();

#define DEBUG_PREP(functor, functor_number)                                                 \
  {                                                                                         \
    if constexpr (!std::is_same_v<decltype(functor), int>) {                                \
      std::string functorName    = typeid(decltype(functor)).name();                        \
      int status                 = 0;                                                       \
      char* demangledFunctorName = 0;                                                       \
      demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status); \
      functorNames.push_back(demangledFunctorName);                                         \
    }                                                                                       \
  }
    DEBUG_PREP(functor0, 0);
    DEBUG_PREP(functor1, 1);
    DEBUG_PREP(functor2, 2);
    DEBUG_PREP(functor3, 3);
    DEBUG_PREP(functor4, 4);
    DEBUG_PREP(functor5, 5);
    DEBUG_PREP(functor6, 6);
    DEBUG_PREP(functor7, 7);
    DEBUG_PREP(functor8, 8);
    DEBUG_PREP(functor9, 9);
#undef DEBUG_PREP
#endif
  }

 public:
  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(0)
    , functor1(0)
    , functor2(0)
    , functor3(0)
    , functor4(0)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type_0& functor0_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(functor0_)
    , functor1(0)
    , functor2(0)
    , functor3(0)
    , functor4(0)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type_0& functor0_, functor_type_1& functor1_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(0)
    , functor3(0)
    , functor4(0)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(0)
    , functor4(0)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(0)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                           functor_type_4& functor4_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                           functor_type_4& functor4_, functor_type_5& functor5_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(functor5_)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                           functor_type_4& functor4_, functor_type_5& functor5_, functor_type_6& functor6_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(functor5_)
    , functor6(functor6_)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                           functor_type_4& functor4_, functor_type_5& functor5_, functor_type_6& functor6_, functor_type_7& functor7_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(functor5_)
    , functor6(functor6_)
    , functor7(functor7_)
    , functor8(0)
    , functor9(0) {
    init();
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                           functor_type_4& functor4_, functor_type_5& functor5_, functor_type_6& functor6_, functor_type_7& functor7_, functor_type_8& functor8_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(functor5_)
    , functor6(functor6_)
    , functor7(functor7_)
    , functor8(functor8_)
    , functor9(0) {
    init();
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                           functor_type_4& functor4_, functor_type_5& functor5_, functor_type_6& functor6_, functor_type_7& functor7_, functor_type_8& functor8_, functor_type_9& functor9_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(functor5_)
    , functor6(functor6_)
    , functor7(functor7_)
    , functor8(functor8_)
    , functor9(functor9_) {
    init();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid, local_ordinal_type& nnz, const bool& final) const {
#ifdef MUELU_COALESCE_DROP_DEBUG
    {
      Kokkos::printf("\nStarting on row %d\n", rlid);

      auto row = A.rowConst(rlid);

      Kokkos::printf("indices:    ");
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        auto clid = row.colidx(k);
        Kokkos::printf("%5d ", clid);
      }
      Kokkos::printf("\n");

      Kokkos::printf("values:     ");
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        auto val = row.value(k);
        Kokkos::printf("%5f ", val);
      }
      Kokkos::printf("\n");
    }
#endif

#ifdef MUELU_COALESCE_DROP_DEBUG
#define APPLY_FUNCTOR(functor, functor_number)                        \
  {                                                                   \
    if constexpr (!std::is_same_v<decltype(functor), int>) {          \
      functor(rlid);                                                  \
      {                                                               \
        Kokkos::printf("%s\n", functorNames[functor_number].c_str()); \
                                                                      \
        auto row            = A.rowConst(rlid);                       \
        const size_t offset = A.graph.row_map(rlid);                  \
                                                                      \
        Kokkos::printf("decisions:  ");                               \
        for (local_ordinal_type k = 0; k < row.length; ++k) {         \
          Kokkos::printf("%5d ", results(offset + k));                \
        }                                                             \
        Kokkos::printf("\n");                                         \
      }                                                               \
    }                                                                 \
  }
#else
#define APPLY_FUNCTOR(functor, functor_number)               \
  {                                                          \
    if constexpr (!std::is_same_v<decltype(functor), int>) { \
      functor(rlid);                                         \
    }                                                        \
  }
#endif
    APPLY_FUNCTOR(functor0, 0);
    APPLY_FUNCTOR(functor1, 1);
    APPLY_FUNCTOR(functor2, 2);
    APPLY_FUNCTOR(functor3, 3);
    APPLY_FUNCTOR(functor4, 4);
    APPLY_FUNCTOR(functor5, 5);
    APPLY_FUNCTOR(functor6, 6);
    APPLY_FUNCTOR(functor7, 7);
    APPLY_FUNCTOR(functor8, 8);
    APPLY_FUNCTOR(functor9, 9);

#undef APPLY_FUNCTOR

#ifdef MUELU_COALESCE_DROP_DEBUG
    {
      Kokkos::printf("Done with row %d\n", rlid);
    }
#endif

    size_t start = A.graph.row_map(rlid);
    size_t end   = A.graph.row_map(rlid + 1);
    for (size_t i = start; i < end; ++i) {
      if (results(i) == KEEP) {
        ++nnz;
      }
    }
    if (final)
      rowptr(rlid + 1) = nnz;
  }
};

template <class local_matrix_type, class local_graph_type, bool lumping>
class PointwiseFillReuseFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;
  using ATS                = Kokkos::ArithTraits<scalar_type>;

  local_matrix_type A;
  results_view results;
  local_matrix_type filteredA;
  local_graph_type graph;
  const scalar_type zero = ATS::zero();

 public:
  PointwiseFillReuseFunctor(local_matrix_type& A_, results_view& results_, local_matrix_type& filteredA_, local_graph_type& graph_)
    : A(A_)
    , results(results_)
    , filteredA(filteredA_)
    , graph(graph_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto rowA                       = A.row(rlid);
    size_t K                        = A.graph.row_map(rlid);
    auto rowFilteredA               = filteredA.row(rlid);
    local_ordinal_type j            = 0;
    local_ordinal_type jj           = 0;
    local_ordinal_type graph_offset = graph.row_map(rlid);
    scalar_type diagCorrection      = zero;
    local_ordinal_type diagOffset   = -1;
    for (local_ordinal_type k = 0; k < rowA.length; ++k) {
      if constexpr (lumping) {
        local_ordinal_type clid = rowA.colidx(k);
        if (rlid == clid) {
          diagOffset = j;
        }
      }
      if (results(K + k) == KEEP) {
        rowFilteredA.colidx(j) = rowA.colidx(k);
        rowFilteredA.value(j)  = rowA.value(k);
        ++j;
        graph.entries(graph_offset + jj) = rowA.colidx(k);
        ++jj;
      } else if constexpr (lumping) {
        diagCorrection += rowA.value(k);
        rowFilteredA.colidx(j) = rowA.colidx(k);
        rowFilteredA.value(j)  = zero;
        ++j;
      } else {
        rowFilteredA.colidx(j) = rowA.colidx(k);
        rowFilteredA.value(j)  = zero;
        ++j;
      }
    }
    if constexpr (lumping) {
      rowFilteredA.value(diagOffset) += diagCorrection;
    }
  }
};

template <class local_matrix_type, bool lumping>
class PointwiseFillNoReuseFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;
  using ATS                = Kokkos::ArithTraits<scalar_type>;

  local_matrix_type A;
  results_view results;
  local_matrix_type filteredA;
  const scalar_type zero = ATS::zero();

 public:
  PointwiseFillNoReuseFunctor(local_matrix_type& A_, results_view& results_, local_matrix_type& filteredA_)
    : A(A_)
    , results(results_)
    , filteredA(filteredA_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto rowA                     = A.row(rlid);
    size_t K                      = A.graph.row_map(rlid);
    auto rowFilteredA             = filteredA.row(rlid);
    local_ordinal_type j          = 0;
    scalar_type diagCorrection    = zero;
    local_ordinal_type diagOffset = -1;
    for (local_ordinal_type k = 0; k < rowA.length; ++k) {
      if constexpr (lumping) {
        local_ordinal_type clid = rowA.colidx(k);
        if (rlid == clid) {
          diagOffset = j;
        }
      }
      if (results(K + k) == KEEP) {
        rowFilteredA.colidx(j) = rowA.colidx(k);
        rowFilteredA.value(j)  = rowA.value(k);
        ++j;
      } else if constexpr (lumping) {
        diagCorrection += rowA.value(k);
      }
    }
    if constexpr (lumping) {
      rowFilteredA.value(diagOffset) += diagCorrection;
    }
  }
};

template <class local_matrix_type,
          class functor_type_0 = int,
          class functor_type_1 = int,
          class functor_type_2 = int,
          class functor_type_3 = int,
          class functor_type_4 = int,
          class functor_type_5 = int,
          class functor_type_6 = int,
          class functor_type_7 = int,
          class functor_type_8 = int,
          class functor_type_9 = int>
class VectorCountingFunctor {
 private:
  using scalar_type             = typename local_matrix_type::value_type;
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;

  using rowptr_type = typename local_matrix_type::row_map_type::non_const_type;
  using ATS         = Kokkos::ArithTraits<local_ordinal_type>;

  local_matrix_type A;
  local_ordinal_type blockSize;
  block_indices_view_type ghosted_point_to_block;
  results_view results;
  rowptr_type filtered_rowptr;
  rowptr_type graph_rowptr;

  functor_type_0 functor0;
  functor_type_1 functor1;
  functor_type_2 functor2;
  functor_type_3 functor3;
  functor_type_4 functor4;
  functor_type_5 functor5;
  functor_type_6 functor6;
  functor_type_7 functor7;
  functor_type_8 functor8;
  functor_type_9 functor9;

  std::vector<std::string> functorNames;

  void init() {
#ifdef MUELU_COALESCE_DROP_DEBUG
    functorNames = std::vector<std::string>();

#define DEBUG_PREP(functor, functor_number)                                                 \
  {                                                                                         \
    if constexpr (!std::is_same_v<decltype(functor), int>) {                                \
      std::string functorName    = typeid(decltype(functor)).name();                        \
      int status                 = 0;                                                       \
      char* demangledFunctorName = 0;                                                       \
      demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status); \
      functorNames.push_back(demangledFunctorName);                                         \
    }                                                                                       \
  }
    DEBUG_PREP(functor0, 0);
    DEBUG_PREP(functor1, 1);
    DEBUG_PREP(functor2, 2);
    DEBUG_PREP(functor3, 3);
    DEBUG_PREP(functor4, 4);
    DEBUG_PREP(functor5, 5);
    DEBUG_PREP(functor6, 6);
    DEBUG_PREP(functor7, 7);
    DEBUG_PREP(functor8, 8);
    DEBUG_PREP(functor9, 9);
#undef DEBUG_PREP
#endif
  }

 public:
  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(0)
    , functor1(0)
    , functor2(0)
    , functor3(0)
    , functor4(0)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_,
                        functor_type_0& functor0_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(functor0_)
    , functor1(0)
    , functor2(0)
    , functor3(0)
    , functor4(0)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_,
                        functor_type_0& functor0_, functor_type_1& functor1_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(0)
    , functor3(0)
    , functor4(0)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_,
                        functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(0)
    , functor4(0)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_,
                        functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(0)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_,
                        functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                        functor_type_4& functor4_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(0)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_,
                        functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                        functor_type_4& functor4_, functor_type_5& functor5_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(functor5_)
    , functor6(0)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_,
                        functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                        functor_type_4& functor4_, functor_type_5& functor5_, functor_type_6& functor6_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(functor5_)
    , functor6(functor6_)
    , functor7(0)
    , functor8(0)
    , functor9(0) {
    init();
  }

  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_,
                        functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                        functor_type_4& functor4_, functor_type_5& functor5_, functor_type_6& functor6_, functor_type_7& functor7_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(functor5_)
    , functor6(functor6_)
    , functor7(functor7_)
    , functor8(0)
    , functor9(0) {
    init();
  }

  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_,
                        functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                        functor_type_4& functor4_, functor_type_5& functor5_, functor_type_6& functor6_, functor_type_7& functor7_, functor_type_8& functor8_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(functor5_)
    , functor6(functor6_)
    , functor7(functor7_)
    , functor8(functor8_)
    , functor9(0) {
    init();
  }

  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_,
                        functor_type_0& functor0_, functor_type_1& functor1_, functor_type_2& functor2_, functor_type_3& functor3_,
                        functor_type_4& functor4_, functor_type_5& functor5_, functor_type_6& functor6_, functor_type_7& functor7_, functor_type_8& functor8_, functor_type_9& functor9_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor0(functor0_)
    , functor1(functor1_)
    , functor2(functor2_)
    , functor3(functor3_)
    , functor4(functor4_)
    , functor5(functor5_)
    , functor6(functor6_)
    , functor7(functor7_)
    , functor8(functor8_)
    , functor9(functor9_) {
    init();
  }

  KOKKOS_INLINE_FUNCTION
  void join(Kokkos::pair<local_ordinal_type, local_ordinal_type>& dest, const Kokkos::pair<local_ordinal_type, local_ordinal_type>& src) const {
    dest.first += src.first;
    dest.second += src.second;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type brlid, Kokkos::pair<local_ordinal_type, local_ordinal_type>& nnz, const bool& final) const {
    auto nnz_filtered = &nnz.first;
    auto nnz_graph    = &nnz.second;
#ifdef MUELU_COALESCE_DROP_DEBUG
    Kokkos::printf("\nStarting on block row %d\n", brlid);
#endif
    for (local_ordinal_type rlid = blockSize * brlid; rlid < blockSize * (brlid + 1); ++rlid) {
#ifdef MUELU_COALESCE_DROP_DEBUG
      {
        Kokkos::printf("\nStarting on row %d\n", rlid);

        auto row = A.rowConst(rlid);

        Kokkos::printf("indices:    ");
        for (local_ordinal_type k = 0; k < row.length; ++k) {
          auto clid = row.colidx(k);
          Kokkos::printf("%5d ", clid);
        }
        Kokkos::printf("\n");

        Kokkos::printf("values:     ");
        for (local_ordinal_type k = 0; k < row.length; ++k) {
          auto val = row.value(k);
          Kokkos::printf("%5f ", val);
        }
        Kokkos::printf("\n");
      }
#endif

#ifdef MUELU_COALESCE_DROP_DEBUG
#define APPLY_FUNCTOR(functor, functor_number)                        \
  {                                                                   \
    if constexpr (!std::is_same_v<decltype(functor), int>) {          \
      functor(rlid);                                                  \
      {                                                               \
        Kokkos::printf("%s\n", functorNames[functor_number].c_str()); \
                                                                      \
        auto row            = A.rowConst(rlid);                       \
        const size_t offset = A.graph.row_map(rlid);                  \
                                                                      \
        Kokkos::printf("decisions:  ");                               \
        for (local_ordinal_type k = 0; k < row.length; ++k) {         \
          Kokkos::printf("%5d ", results(offset + k));                \
        }                                                             \
        Kokkos::printf("\n");                                         \
      }                                                               \
    }                                                                 \
  }
#else
#define APPLY_FUNCTOR(functor, functor_number)               \
  {                                                          \
    if constexpr (!std::is_same_v<decltype(functor), int>) { \
      functor(rlid);                                         \
    }                                                        \
  }
#endif

      APPLY_FUNCTOR(functor0, 0);
      APPLY_FUNCTOR(functor1, 1);
      APPLY_FUNCTOR(functor2, 2);
      APPLY_FUNCTOR(functor3, 3);
      APPLY_FUNCTOR(functor4, 4);
      APPLY_FUNCTOR(functor5, 5);
      APPLY_FUNCTOR(functor6, 6);
      APPLY_FUNCTOR(functor7, 7);
      APPLY_FUNCTOR(functor8, 8);
      APPLY_FUNCTOR(functor9, 9);

#undef APPLY_FUNCTOR

#ifdef MUELU_COALESCE_DROP_DEBUG
      Kokkos::printf("Done with row %d\n", rlid);
#endif

      size_t start = A.graph.row_map(rlid);
      size_t end   = A.graph.row_map(rlid + 1);
      for (size_t i = start; i < end; ++i) {
        if (results(i) == KEEP) {
          ++(*nnz_filtered);
        }
      }
      if (final)
        filtered_rowptr(rlid + 1) = *nnz_filtered;
    }

#ifdef MUELU_COALESCE_DROP_DEBUG
    Kokkos::printf("Done with block row %d\nGraph indices ", brlid);
#endif

    local_ordinal_type* nextIndices = new local_ordinal_type[blockSize];
    for (local_ordinal_type block_index = 0; block_index < blockSize; ++block_index) {
      nextIndices[block_index] = 0;
    }
    local_ordinal_type prev_bclid = -1;
    while (true) {
      local_ordinal_type min_block_index = -1;
      local_ordinal_type min_clid        = ATS::max();
      local_ordinal_type min_offset      = -1;
      for (local_ordinal_type block_index = 0; block_index < blockSize; ++block_index) {
        auto rlid   = blockSize * brlid + block_index;
        auto offset = A.graph.row_map(rlid) + nextIndices[block_index];
        if (offset == A.graph.row_map(rlid + 1))
          continue;
        auto clid = A.graph.entries(offset);
        if (clid < min_clid) {
          min_block_index = block_index;
          min_clid        = clid;
          min_offset      = offset;
        }
      }
      if (min_block_index == -1)
        break;
      ++nextIndices[min_block_index];
      auto bclid = ghosted_point_to_block(min_clid);
      if (prev_bclid < bclid) {
        if (results(min_offset) == KEEP) {
          ++(*nnz_graph);
#ifdef MUELU_COALESCE_DROP_DEBUG
          Kokkos::printf("%5d ", bclid);
#endif
          prev_bclid = bclid;
        }
      }
    }
#ifdef MUELU_COALESCE_DROP_DEBUG
    Kokkos::printf("\n");
#endif
    if (final)
      graph_rowptr(brlid + 1) = *nnz_graph;
  }
};

template <class local_matrix_type, bool lumping, bool reuse>
class VectorFillFunctor {
 private:
  using scalar_type             = typename local_matrix_type::value_type;
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using local_graph_type        = typename local_matrix_type::staticcrsgraph_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using ATS                     = Kokkos::ArithTraits<scalar_type>;
  using OTS                     = Kokkos::ArithTraits<local_ordinal_type>;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;

  local_matrix_type A;
  local_ordinal_type blockSize;
  block_indices_view_type ghosted_point_to_block;
  results_view results;
  local_matrix_type filteredA;
  local_graph_type graph;
  const scalar_type zero = ATS::zero();

 public:
  VectorFillFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, local_matrix_type& filteredA_, local_graph_type& graph_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filteredA(filteredA_)
    , graph(graph_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type brlid) const {
    for (local_ordinal_type rlid = blockSize * brlid; rlid < blockSize * (brlid + 1); ++rlid) {
      auto rowA                     = A.row(rlid);
      size_t K                      = A.graph.row_map(rlid);
      auto rowFilteredA             = filteredA.row(rlid);
      local_ordinal_type j          = 0;
      scalar_type diagCorrection    = zero;
      local_ordinal_type diagOffset = -1;
      for (local_ordinal_type k = 0; k < rowA.length; ++k) {
        if constexpr (lumping) {
          local_ordinal_type clid = rowA.colidx(k);
          if (rlid == clid) {
            diagOffset = j;
          }
        }
        if (results(K + k) == KEEP) {
          rowFilteredA.colidx(j) = rowA.colidx(k);
          rowFilteredA.value(j)  = rowA.value(k);
          ++j;
        } else if constexpr (lumping) {
          diagCorrection += rowA.value(k);
          if constexpr (reuse) {
            rowFilteredA.colidx(j) = rowA.colidx(k);
            rowFilteredA.value(j)  = zero;
            ++j;
          }
        } else if constexpr (reuse) {
          rowFilteredA.colidx(j) = rowA.colidx(k);
          rowFilteredA.value(j)  = zero;
          ++j;
        }
      }
      if constexpr (lumping) {
        rowFilteredA.value(diagOffset) += diagCorrection;
      }
    }

    local_ordinal_type* nextIndices = new local_ordinal_type[blockSize];
    for (local_ordinal_type block_index = 0; block_index < blockSize; ++block_index) {
      nextIndices[block_index] = 0;
    }
    local_ordinal_type prev_bclid = -1;

    local_ordinal_type j = graph.row_map(brlid);
    while (true) {
      local_ordinal_type min_block_index = -1;
      local_ordinal_type min_clid        = OTS::max();
      local_ordinal_type min_offset      = -1;
      for (local_ordinal_type block_index = 0; block_index < blockSize; ++block_index) {
        auto rlid   = blockSize * brlid + block_index;
        auto offset = A.graph.row_map(rlid) + nextIndices[block_index];
        if (offset == A.graph.row_map(rlid + 1))
          continue;
        auto clid = A.graph.entries(offset);
        if (clid < min_clid) {
          min_block_index = block_index;
          min_clid        = clid;
          min_offset      = offset;
        }
      }
      if (min_block_index == -1)
        break;
      ++nextIndices[min_block_index];
      auto bclid = ghosted_point_to_block(min_clid);
      if (prev_bclid < bclid) {
        if (results(min_offset) == KEEP) {
          graph.entries(j) = bclid;
          ++j;
          prev_bclid = bclid;
        }
      }
    }
  }
};

}  // namespace MueLu::MatrixConstruction

#endif
