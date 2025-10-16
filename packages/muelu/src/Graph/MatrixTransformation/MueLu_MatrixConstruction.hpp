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
#if KOKKOS_VERSION >= 40799
#include "KokkosKernels_ArithTraits.hpp"
#else
#include "Kokkos_ArithTraits.hpp"
#endif

#include "MueLu_DroppingCommon.hpp"

#include "Xpetra_MatrixFactory.hpp"

#ifdef MUELU_COALESCE_DROP_DEBUG
// For demangling function names
#include <cxxabi.h>
#endif

namespace MueLu::MatrixConstruction {
/*!
  @class PointwiseCountingFunctor
  @brief Functor that executes a sequence of sub-functors on each row for a problem with blockSize == 1.

  The functor applies a series of functors to each row of the matrix.
  Each sub-functor can modify the decision to drop or keep any matrix entry in the given row.
  These decisions are applied to the results_view.
  Once a row has been processed by all sub-functors, the number of entries in the row after dropping is determined.
  The result is saved as offsets in rowptr.
*/
template <class local_matrix_type, class functor_type, class... remaining_functor_types>
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
  functor_type functor;
  PointwiseCountingFunctor<local_matrix_type, remaining_functor_types...> remainingFunctors;
  bool firstFunctor;

#ifdef MUELU_COALESCE_DROP_DEBUG
  std::string functorName;
#endif

 public:
  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type& functor_, remaining_functor_types&... remainingFunctors_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor(functor_)
    , remainingFunctors(A_, results_, rowptr_, false, remainingFunctors_...)
    , firstFunctor(true) {
#ifdef MUELU_COALESCE_DROP_DEBUG
    std::string mangledFunctorName = typeid(decltype(functor)).name();
    int status                     = 0;
    char* demangledFunctorName     = 0;
    demangledFunctorName           = abi::__cxa_demangle(mangledFunctorName.c_str(), 0, 0, &status);
    functorName                    = demangledFunctorName;
#endif
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, bool firstFunctor_, functor_type& functor_, remaining_functor_types&... remainingFunctors_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor(functor_)
    , remainingFunctors(A_, results_, rowptr_, false, remainingFunctors_...)
    , firstFunctor(firstFunctor_) {
#ifdef MUELU_COALESCE_DROP_DEBUG
    std::string mangledFunctorName = typeid(decltype(functor)).name();
    int status                     = 0;
    char* demangledFunctorName     = 0;
    demangledFunctorName           = abi::__cxa_demangle(mangledFunctorName.c_str(), 0, 0, &status);
    functorName                    = demangledFunctorName;
#endif
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid, local_ordinal_type& nnz, const bool& final) const {
#ifdef MUELU_COALESCE_DROP_DEBUG
    if (firstFunctor) {
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

    functor(rlid);

#ifdef MUELU_COALESCE_DROP_DEBUG
    {
      Kokkos::printf("%s\n", functorName.c_str());

      auto row            = A.rowConst(rlid);
      const size_t offset = A.graph.row_map(rlid);

      Kokkos::printf("decisions:  ");
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        Kokkos::printf("%5d ", results(offset + k));
      }
      Kokkos::printf("\n");
    }
#endif

    remainingFunctors(rlid, nnz, final);
  }
};

template <class local_matrix_type, class functor_type>
class PointwiseCountingFunctor<local_matrix_type, functor_type> {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  using rowptr_type = typename local_matrix_type::row_map_type::non_const_type;

  local_matrix_type A;
  results_view results;
  rowptr_type rowptr;
  functor_type functor;
  bool firstFunctor;

#ifdef MUELU_COALESCE_DROP_DEBUG
  std::string functorName;
#endif

 public:
  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, functor_type& functor_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor(functor_)
    , firstFunctor(true) {
#ifdef MUELU_COALESCE_DROP_DEBUG
    std::string mangledFunctorName = typeid(decltype(functor)).name();
    int status                     = 0;
    char* demangledFunctorName     = 0;
    demangledFunctorName           = abi::__cxa_demangle(mangledFunctorName.c_str(), 0, 0, &status);
    functorName                    = demangledFunctorName;
#endif
  }

  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, bool firstFunctor_, functor_type& functor_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functor(functor_)
    , firstFunctor(firstFunctor_) {
#ifdef MUELU_COALESCE_DROP_DEBUG
    std::string mangledFunctorName = typeid(decltype(functor)).name();
    int status                     = 0;
    char* demangledFunctorName     = 0;
    demangledFunctorName           = abi::__cxa_demangle(mangledFunctorName.c_str(), 0, 0, &status);
    functorName                    = demangledFunctorName;
#endif
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid, local_ordinal_type& nnz, const bool& final) const {
#ifdef MUELU_COALESCE_DROP_DEBUG
    if (firstFunctor) {
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

    functor(rlid);

#ifdef MUELU_COALESCE_DROP_DEBUG
    Kokkos::printf("%s\n", functorName.c_str());

    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);

    Kokkos::printf("decisions:  ");
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      Kokkos::printf("%5d ", results(offset + k));
    }

    Kokkos::printf("\n");
    Kokkos::printf("Done with row %d\n", rlid);
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

/*!
  @class PointwiseFillReuseFunctor
  @brief Functor that fills the filtered matrix while reusing the graph of the matrix before dropping, blockSize == 1.

  The dropped graph is built from scratch.
  The filtered matrix reuses the graph of the matrix before dropping.
  Lumps dropped entries to the diagonal if lumping==true.
*/
template <class local_matrix_type, class local_graph_type, bool lumping>
class PointwiseFillReuseFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
  using magnitudeType = typename ATS::magnitudeType;

  local_matrix_type A;
  results_view results;
  local_matrix_type filteredA;
  local_graph_type graph;
  magnitudeType dirichletThreshold;
  const scalar_type zero = ATS::zero();
  const scalar_type one  = ATS::one();

 public:
  PointwiseFillReuseFunctor(local_matrix_type& A_, results_view& results_, local_matrix_type& filteredA_, local_graph_type& graph_, magnitudeType dirichletThreshold_)
    : A(A_)
    , results(results_)
    , filteredA(filteredA_)
    , graph(graph_)
    , dirichletThreshold(dirichletThreshold_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto rowA                       = A.row(rlid);
    size_t row_start                = A.graph.row_map(rlid);
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
      if (results(row_start + k) == KEEP) {
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
      if ((dirichletThreshold >= 0.0) && (ATS::real(rowFilteredA.value(diagOffset)) <= dirichletThreshold))
        rowFilteredA.value(diagOffset) = one;
    }
  }
};
enum lumpingType : int { no_lumping,
                         diag_lumping,
                         distributed_lumping };

/*!
  @class PointwiseFillNoReuseFunctor
  @brief Functor does not reuse the graph of the matrix for a problem with blockSize == 1.

  The dropped graph and the filtered matrix are built from scratch.
  Lumps dropped entries to the diagonal if lumpingChoice==diag_lumping.
  Lumps dropped entries across all kept entries (proportional to their magnitude) if lumpingChoice==ddistributed_lumping.
*/
template <class local_matrix_type, lumpingType lumpingChoice>
class PointwiseFillNoReuseFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
  using magnitudeType = typename ATS::magnitudeType;

  local_matrix_type A;
  results_view results;
  local_matrix_type filteredA;
  magnitudeType dirichletThreshold;
  const scalar_type zero = ATS::zero();
  const scalar_type one  = ATS::one();

 public:
  PointwiseFillNoReuseFunctor(local_matrix_type& A_, results_view& results_, local_matrix_type& filteredA_, magnitudeType dirichletThreshold_)
    : A(A_)
    , results(results_)
    , filteredA(filteredA_)
    , dirichletThreshold(dirichletThreshold_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto rowA                     = A.row(rlid);
    size_t K                      = A.graph.row_map(rlid);
    auto rowFilteredA             = filteredA.row(rlid);
    local_ordinal_type j          = 0;
    scalar_type droppedSum        = zero;
    scalar_type keptRowSumAbs     = zero;
    local_ordinal_type diagOffset = -1;
    for (local_ordinal_type k = 0; k < rowA.length; ++k) {
      if constexpr (lumpingChoice != no_lumping) {
        local_ordinal_type clid = rowA.colidx(k);
        if (rlid == clid) {
          diagOffset = j;
        }
      }
      if (results(K + k) == KEEP) {
        rowFilteredA.colidx(j) = rowA.colidx(k);
        rowFilteredA.value(j)  = rowA.value(k);
        if constexpr (lumpingChoice == distributed_lumping) {
          keptRowSumAbs += ATS::magnitude(rowFilteredA.value(j));
        }
        ++j;
      } else if constexpr (lumpingChoice != no_lumping) {
        droppedSum += rowA.value(k);
      }
    }
    if constexpr (lumpingChoice == diag_lumping) {
      rowFilteredA.value(diagOffset) += droppedSum;
      if ((dirichletThreshold >= 0.0) && (ATS::real(rowFilteredA.value(diagOffset)) <= dirichletThreshold))
        rowFilteredA.value(diagOffset) = one;
    } else if constexpr (lumpingChoice == distributed_lumping) {
      if (ATS::real(droppedSum) >= ATS::real(zero)) {
        rowFilteredA.value(diagOffset) += droppedSum;

      } else {
        for (local_ordinal_type k = 0; k < j; ++k) {
          rowFilteredA.value(k) += droppedSum * ATS::magnitude(rowFilteredA.value(k)) / keptRowSumAbs;
        }
      }
    }
  }
};

template <class local_matrix_type>
class BlockRowComparison {
 public:
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;

  local_matrix_type A;
  local_ordinal_type bsize;
  block_indices_view_type ghosted_point_to_block;

 public:
  BlockRowComparison(local_matrix_type& A_, local_ordinal_type bsize_, block_indices_view_type ghosted_point_to_block_)
    : A(A_)
    , bsize(bsize_)
    , ghosted_point_to_block(ghosted_point_to_block_) {}

  template <class local_matrix_type2>
  struct Comparator {
   private:
    using local_ordinal_type      = typename local_matrix_type2::ordinal_type;
    using memory_space            = typename local_matrix_type2::memory_space;
    using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;

    const local_matrix_type2 A;
    const local_ordinal_type offset;
    const block_indices_view_type ghosted_point_to_block;

   public:
    KOKKOS_INLINE_FUNCTION
    Comparator(const local_matrix_type2& A_, local_ordinal_type bsize_, local_ordinal_type brlid_, block_indices_view_type ghosted_point_to_block_)
      : A(A_)
      , offset(A_.graph.row_map(bsize_ * brlid_))
      , ghosted_point_to_block(ghosted_point_to_block_) {}

    KOKKOS_INLINE_FUNCTION
    bool operator()(size_t x, size_t y) const {
      return ghosted_point_to_block(A.graph.entries(offset + x)) < ghosted_point_to_block(A.graph.entries(offset + y));
    }
  };

  using comparator_type = Comparator<local_matrix_type>;

  KOKKOS_INLINE_FUNCTION
  comparator_type getComparator(local_ordinal_type brlid) const {
    return comparator_type(A, bsize, brlid, ghosted_point_to_block);
  }
};

/*!
  @class VectorCountingFunctor
  @brief Functor that executes a sequence of sub-functors on each block of rows.

  The functor applies a series of functors to each row of the matrix.
  Each sub-functor can modify the decision to drop or keep any matrix entry in the given row.
  These decisions are applied to the results_view.
  Once a row has been processed by all sub-functors, the number of entries in the row after dropping is determined.
  The result is saved as offsets in rowptr.
*/
template <class local_matrix_type,
          class functor_type,
          class... remaining_functor_types>
class VectorCountingFunctor {
 private:
  using scalar_type             = typename local_matrix_type::value_type;
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;
  using permutation_type        = Kokkos::View<local_ordinal_type*, memory_space>;

  using rowptr_type = typename local_matrix_type::row_map_type::non_const_type;
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<local_ordinal_type>;
#else
  using ATS = Kokkos::ArithTraits<local_ordinal_type>;
#endif

  local_matrix_type A;
  local_ordinal_type blockSize;
  block_indices_view_type ghosted_point_to_block;
  results_view results;
  rowptr_type filtered_rowptr;
  rowptr_type graph_rowptr;

  functor_type functor;

  BlockRowComparison<local_matrix_type> comparison;
  permutation_type permutation;

  VectorCountingFunctor<local_matrix_type, remaining_functor_types...> remainingFunctors;

#ifdef MUELU_COALESCE_DROP_DEBUG
  std::string functorName;
#endif

 public:
  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_, functor_type& functor_, remaining_functor_types&... remainingFunctors_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor(functor_)
    , comparison(BlockRowComparison(A, blockSize_, ghosted_point_to_block))
    , remainingFunctors(A_, blockSize_, ghosted_point_to_block_, results_, filtered_rowptr_, graph_rowptr_, remainingFunctors_...) {
    permutation = permutation_type("permutation", A.nnz());
#ifdef MUELU_COALESCE_DROP_DEBUG
    std::string mangledFunctorName = typeid(decltype(functor)).name();
    int status                     = 0;
    char* demangledFunctorName     = 0;
    demangledFunctorName           = abi::__cxa_demangle(mangledFunctorName.c_str(), 0, 0, &status);
    functorName                    = demangledFunctorName;
#endif
  }

  KOKKOS_INLINE_FUNCTION
  void join(Kokkos::pair<local_ordinal_type, local_ordinal_type>& dest, const Kokkos::pair<local_ordinal_type, local_ordinal_type>& src) const {
    dest.first += src.first;
    dest.second += src.second;
  }

  KOKKOS_INLINE_FUNCTION
  void operatorRow(const local_ordinal_type rlid) const {
    functor(rlid);
    remainingFunctors.operatorRow(rlid);
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

      functor(rlid);
      remainingFunctors.operatorRow(rlid);

#ifdef MUELU_COALESCE_DROP_DEBUG
      {
        Kokkos::printf("%s\n", functorName.c_str());

        auto row            = A.rowConst(rlid);
        const size_t offset = A.graph.row_map(rlid);

        Kokkos::printf("decisions:  ");
        for (local_ordinal_type k = 0; k < row.length; ++k) {
          Kokkos::printf("%5d ", results(offset + k));
        }
        Kokkos::printf("\n");
      }
#endif

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

    // column lids for all rows in the block
    auto block_clids = Kokkos::subview(A.graph.entries, Kokkos::make_pair(A.graph.row_map(blockSize * brlid),
                                                                          A.graph.row_map(blockSize * (brlid + 1))));
    // set up a permutatation index
    auto block_permutation = Kokkos::subview(permutation, Kokkos::make_pair(A.graph.row_map(blockSize * brlid),
                                                                            A.graph.row_map(blockSize * (brlid + 1))));
    for (size_t i = 0; i < block_permutation.extent(0); ++i)
      block_permutation(i) = i;
    // get permuatation for sorted column indices of the entire block
    auto comparator = comparison.getComparator(brlid);
    Misc::serialHeapSort(block_permutation, comparator);

    local_ordinal_type prev_bclid = -1;
    bool alreadyAdded             = false;

    // loop over all sorted entries in block
    auto offset = A.graph.row_map(blockSize * brlid);
    for (size_t i = 0; i < block_permutation.extent(0); ++i) {
      auto idx   = offset + block_permutation(i);
      auto clid  = A.graph.entries(idx);
      auto bclid = ghosted_point_to_block(clid);

      // unseen block column index
      if (bclid > prev_bclid)
        alreadyAdded = false;

      // add entry to graph
      if (!alreadyAdded && (results(idx) == KEEP)) {
        ++(*nnz_graph);
        alreadyAdded = true;
#ifdef MUELU_COALESCE_DROP_DEBUG
        Kokkos::printf("%5d ", bclid);
#endif
      }
      prev_bclid = bclid;
    }
#ifdef MUELU_COALESCE_DROP_DEBUG
    Kokkos::printf("\n");
#endif
    if (final)
      graph_rowptr(brlid + 1) = *nnz_graph;
  }
};

template <class local_matrix_type,
          class functor_type>
class VectorCountingFunctor<local_matrix_type, functor_type> {
 private:
  using scalar_type             = typename local_matrix_type::value_type;
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;
  using permutation_type        = Kokkos::View<local_ordinal_type*, memory_space>;

  using rowptr_type = typename local_matrix_type::row_map_type::non_const_type;
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<local_ordinal_type>;
#else
  using ATS = Kokkos::ArithTraits<local_ordinal_type>;
#endif

  local_matrix_type A;
  local_ordinal_type blockSize;
  block_indices_view_type ghosted_point_to_block;
  results_view results;
  rowptr_type filtered_rowptr;
  rowptr_type graph_rowptr;

  bool firstFunctor;
  functor_type functor;

#ifdef MUELU_COALESCE_DROP_DEBUG
  std::string functorName;
#endif

  BlockRowComparison<local_matrix_type> comparison;
  permutation_type permutation;

 public:
  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_, functor_type& functor_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functor(functor_)
    , comparison(BlockRowComparison(A, blockSize_, ghosted_point_to_block)) {
    permutation = permutation_type("permutation", A.nnz());
#ifdef MUELU_COALESCE_DROP_DEBUG
    std::string mangledFunctorName = typeid(decltype(functor)).name();
    int status                     = 0;
    char* demangledFunctorName     = 0;
    demangledFunctorName           = abi::__cxa_demangle(mangledFunctorName.c_str(), 0, 0, &status);
    functorName                    = demangledFunctorName;
#endif
  }

  KOKKOS_INLINE_FUNCTION
  void join(Kokkos::pair<local_ordinal_type, local_ordinal_type>& dest, const Kokkos::pair<local_ordinal_type, local_ordinal_type>& src) const {
    dest.first += src.first;
    dest.second += src.second;
  }

  KOKKOS_INLINE_FUNCTION
  void operatorRow(const local_ordinal_type rlid) const {
    functor(rlid);
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

      functor(rlid);

#ifdef MUELU_COALESCE_DROP_DEBUG
      {
        Kokkos::printf("%s\n", functorName.c_str());

        auto row            = A.rowConst(rlid);
        const size_t offset = A.graph.row_map(rlid);

        Kokkos::printf("decisions:  ");
        for (local_ordinal_type k = 0; k < row.length; ++k) {
          Kokkos::printf("%5d ", results(offset + k));
        }
        Kokkos::printf("\n");
      }
#endif

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

    // column lids for all rows in the block
    auto block_clids = Kokkos::subview(A.graph.entries, Kokkos::make_pair(A.graph.row_map(blockSize * brlid),
                                                                          A.graph.row_map(blockSize * (brlid + 1))));
    // set up a permutation index
    auto block_permutation = Kokkos::subview(permutation, Kokkos::make_pair(A.graph.row_map(blockSize * brlid),
                                                                            A.graph.row_map(blockSize * (brlid + 1))));
    for (size_t i = 0; i < block_permutation.extent(0); ++i)
      block_permutation(i) = i;
    // get permutation for sorted column indices of the entire block
    auto comparator = comparison.getComparator(brlid);
    Misc::serialHeapSort(block_permutation, comparator);

    local_ordinal_type prev_bclid = -1;
    bool alreadyAdded             = false;

    // loop over all sorted entries in block
    auto offset = A.graph.row_map(blockSize * brlid);
    for (size_t i = 0; i < block_permutation.extent(0); ++i) {
      auto idx   = offset + block_permutation(i);
      auto clid  = A.graph.entries(idx);
      auto bclid = ghosted_point_to_block(clid);

      // unseen block column index
      if (bclid > prev_bclid)
        alreadyAdded = false;

      // add entry to graph
      if (!alreadyAdded && (results(idx) == KEEP)) {
        ++(*nnz_graph);
        alreadyAdded = true;
#ifdef MUELU_COALESCE_DROP_DEBUG
        Kokkos::printf("%5d ", bclid);
#endif
      }
      prev_bclid = bclid;
    }
#ifdef MUELU_COALESCE_DROP_DEBUG
    Kokkos::printf("\n");
#endif
    if (final)
      graph_rowptr(brlid + 1) = *nnz_graph;
  }
};

/*!
  @class VectorFillNoReuseFunctor
  @brief Functor does not reuse the graph of the matrix for a problem with blockSize>1.

  The dropped graph and the filtered matrix are built from scratch.
  Lumps dropped entries to the diagonal if lumping==true.
*/
template <class local_matrix_type, bool lumping, bool reuse>
class VectorFillFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using local_graph_type   = typename local_matrix_type::staticcrsgraph_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
#if KOKKOS_VERSION >= 40799
  using OTS = KokkosKernels::ArithTraits<local_ordinal_type>;
#else
  using OTS = Kokkos::ArithTraits<local_ordinal_type>;
#endif
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;
  using permutation_type        = Kokkos::View<local_ordinal_type*, memory_space>;
  using magnitudeType           = typename ATS::magnitudeType;

  local_matrix_type A;
  local_ordinal_type blockSize;
  block_indices_view_type ghosted_point_to_block;
  results_view results;
  local_matrix_type filteredA;
  local_graph_type graph;
  magnitudeType dirichletThreshold;
  const scalar_type zero = ATS::zero();
  const scalar_type one  = ATS::one();

  BlockRowComparison<local_matrix_type> comparison;
  permutation_type permutation;

 public:
  VectorFillFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, local_matrix_type& filteredA_, local_graph_type& graph_, magnitudeType dirichletThreshold_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filteredA(filteredA_)
    , graph(graph_)
    , dirichletThreshold(dirichletThreshold_)
    , comparison(BlockRowComparison(A, blockSize_, ghosted_point_to_block)) {
    permutation = permutation_type("permutation", A.nnz());
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type brlid) const {
    for (local_ordinal_type rlid = blockSize * brlid; rlid < blockSize * (brlid + 1); ++rlid) {
      auto rowA                     = A.row(rlid);
      size_t row_start              = A.graph.row_map(rlid);
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
        if (results(row_start + k) == KEEP) {
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
        if ((dirichletThreshold >= 0.0) && (ATS::real(rowFilteredA.value(diagOffset)) <= dirichletThreshold))
          rowFilteredA.value(diagOffset) = one;
      }
    }

    // column lids for all rows in the block
    auto block_clids = Kokkos::subview(A.graph.entries, Kokkos::make_pair(A.graph.row_map(blockSize * brlid),
                                                                          A.graph.row_map(blockSize * (brlid + 1))));
    // set up a permuatation index
    auto block_permutation = Kokkos::subview(permutation, Kokkos::make_pair(A.graph.row_map(blockSize * brlid),
                                                                            A.graph.row_map(blockSize * (brlid + 1))));
    for (size_t i = 0; i < block_permutation.extent(0); ++i)
      block_permutation(i) = i;
    // get permutation for sorted column indices of the entire block
    auto comparator = comparison.getComparator(brlid);
    Misc::serialHeapSort(block_permutation, comparator);

    local_ordinal_type prev_bclid = -1;
    bool alreadyAdded             = false;
    local_ordinal_type j          = graph.row_map(brlid);

    // loop over all sorted entries in block
    auto offset = A.graph.row_map(blockSize * brlid);
    for (size_t i = 0; i < block_permutation.extent(0); ++i) {
      auto idx   = offset + block_permutation(i);
      auto clid  = A.graph.entries(idx);
      auto bclid = ghosted_point_to_block(clid);

      // unseen block column index
      if (bclid > prev_bclid)
        alreadyAdded = false;

      // add entry to graph
      if (!alreadyAdded && (results(idx) == KEEP)) {
        graph.entries(j) = bclid;
        ++j;
        alreadyAdded = true;
      }
      prev_bclid = bclid;
    }
  }
};

template <class local_matrix_type>
class MergeCountFunctor {
 private:
  using scalar_type             = typename local_matrix_type::value_type;
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;
  using permutation_type        = Kokkos::View<local_ordinal_type*, memory_space>;

  using rowptr_type = typename local_matrix_type::row_map_type::non_const_type;
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<local_ordinal_type>;
#else
  using ATS = Kokkos::ArithTraits<local_ordinal_type>;
#endif

  local_matrix_type A;
  local_ordinal_type blockSize;
  block_indices_view_type ghosted_point_to_block;
  rowptr_type merged_rowptr;

  BlockRowComparison<local_matrix_type> comparison;
  permutation_type permutation;

 public:
  MergeCountFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, rowptr_type& merged_rowptr_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , merged_rowptr(merged_rowptr_)
    , comparison(BlockRowComparison(A, blockSize_, ghosted_point_to_block)) {
    permutation = permutation_type("permutation", A.nnz());
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type brlid, local_ordinal_type& nnz_graph, const bool& final) const {
    // column lids for all rows in the block
    auto block_clids = Kokkos::subview(A.graph.entries, Kokkos::make_pair(A.graph.row_map(blockSize * brlid),
                                                                          A.graph.row_map(blockSize * (brlid + 1))));
    // set up a permutation index
    auto block_permutation = Kokkos::subview(permutation, Kokkos::make_pair(A.graph.row_map(blockSize * brlid),
                                                                            A.graph.row_map(blockSize * (brlid + 1))));
    for (size_t i = 0; i < block_permutation.extent(0); ++i)
      block_permutation(i) = i;
    // get permutation for sorted column indices of the entire block
    auto comparator = comparison.getComparator(brlid);
    Misc::serialHeapSort(block_permutation, comparator);

    local_ordinal_type prev_bclid = -1;
    bool alreadyAdded             = false;

    // loop over all sorted entries in block
    auto offset = A.graph.row_map(blockSize * brlid);
    for (size_t i = 0; i < block_permutation.extent(0); ++i) {
      auto idx   = offset + block_permutation(i);
      auto clid  = A.graph.entries(idx);
      auto bclid = ghosted_point_to_block(clid);

      // unseen block column index
      if (bclid > prev_bclid)
        alreadyAdded = false;

      // add entry to graph
      if (!alreadyAdded) {
        ++nnz_graph;
        alreadyAdded = true;
#ifdef MUELU_COALESCE_DROP_DEBUG
        Kokkos::printf("%5d ", bclid);
#endif
      }
      prev_bclid = bclid;
    }
#ifdef MUELU_COALESCE_DROP_DEBUG
    Kokkos::printf("\n");
#endif
    if (final)
      merged_rowptr(brlid + 1) = nnz_graph;
  }
};

template <class local_matrix_type>
class MergeFillFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using local_graph_type   = typename local_matrix_type::staticcrsgraph_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
#if KOKKOS_VERSION >= 40799
  using OTS = KokkosKernels::ArithTraits<local_ordinal_type>;
#else
  using OTS = Kokkos::ArithTraits<local_ordinal_type>;
#endif
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;
  using permutation_type        = Kokkos::View<local_ordinal_type*, memory_space>;
  using magnitudeType           = typename ATS::magnitudeType;

  local_matrix_type A;
  local_ordinal_type blockSize;
  block_indices_view_type ghosted_point_to_block;
  local_matrix_type mergedA;
  const scalar_type zero = ATS::zero();
  const scalar_type one  = ATS::one();

  BlockRowComparison<local_matrix_type> comparison;
  permutation_type permutation;

 public:
  MergeFillFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, local_matrix_type& mergedA_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , mergedA(mergedA_)
    , comparison(BlockRowComparison(A, blockSize_, ghosted_point_to_block)) {
    permutation = permutation_type("permutation", A.nnz());
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type brlid) const {
    // column lids for all rows in the block
    auto block_clids = Kokkos::subview(A.graph.entries, Kokkos::make_pair(A.graph.row_map(blockSize * brlid),
                                                                          A.graph.row_map(blockSize * (brlid + 1))));
    // set up a permuatation index
    auto block_permutation = Kokkos::subview(permutation, Kokkos::make_pair(A.graph.row_map(blockSize * brlid),
                                                                            A.graph.row_map(blockSize * (brlid + 1))));
    for (size_t i = 0; i < block_permutation.extent(0); ++i)
      block_permutation(i) = i;
    // get permutation for sorted column indices of the entire block
    auto comparator = comparison.getComparator(brlid);
    Misc::serialHeapSort(block_permutation, comparator);

    local_ordinal_type prev_bclid = -1;
    bool alreadyAdded             = false;
    local_ordinal_type j          = mergedA.graph.row_map(brlid);

    // loop over all sorted entries in block
    auto offset = A.graph.row_map(blockSize * brlid);
    for (size_t i = 0; i < block_permutation.extent(0); ++i) {
      auto idx   = offset + block_permutation(i);
      auto clid  = A.graph.entries(idx);
      auto bclid = ghosted_point_to_block(clid);

      // unseen block column index
      if (bclid > prev_bclid)
        alreadyAdded = false;

      // add entry to graph
      if (!alreadyAdded) {
        mergedA.graph.entries(j) = bclid;
        mergedA.values(j)        = one;
        ++j;
        alreadyAdded = true;
      }
      prev_bclid = bclid;
    }
  }
};

template <class local_matrix_type, class local_graph_type>
class GraphConstruction {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;
  local_graph_type graph;

 public:
  GraphConstruction(local_matrix_type& A_, results_view& results_, local_graph_type& graph_)
    : A(A_)
    , results(results_)
    , graph(graph_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto rowA                       = A.row(rlid);
    size_t row_start                = A.graph.row_map(rlid);
    local_ordinal_type jj           = 0;
    local_ordinal_type graph_offset = graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < rowA.length; ++k) {
      if (results(row_start + k) == KEEP) {
        graph.entries(graph_offset + jj) = rowA.colidx(k);
        ++jj;
      }
    }
  }
};

}  // namespace MueLu::MatrixConstruction

#endif
