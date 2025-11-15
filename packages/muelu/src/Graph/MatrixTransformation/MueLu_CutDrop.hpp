// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_CUTDROP_HPP
#define MUELU_CUTDROP_HPP

#include "Kokkos_Core.hpp"
#if KOKKOS_VERSION >= 40799
#include "KokkosKernels_ArithTraits.hpp"
#else
#include "Kokkos_ArithTraits.hpp"
#endif
#include "MueLu_DroppingCommon.hpp"
#include "MueLu_Utilities.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MultiVector.hpp"
#include "MueLu_DistanceLaplacianDropping.hpp"

namespace MueLu::CutDrop {

/*! Cut drop algorithm options*/
enum decisionAlgoType { defaultAlgo,
                        unscaled_cut,
                        scaled_cut,
                        scaled_cut_symmetric };

/*!
  @class UnscaledComparison
  @brief Orders entries of row \f$i\f$ by \f$|A_{ij}|^2\f$.
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class UnscaledComparison {
 public:
  using matrix_type = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  using local_matrix_type  = typename matrix_type::local_matrix_device_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 private:
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
  using magnitudeType = typename ATS::magnitudeType;
  using values_view   = Kokkos::View<magnitudeType*, memory_space>;
  mutable values_view values;

 public:
  UnscaledComparison(matrix_type& A_, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , results(results_)
    , values("UnscaledComparison::values", A.nnz()) {}

  template <class local_matrix_type2>
  struct Comparator {
   private:
    using scalar_type        = typename local_matrix_type2::value_type;
    using local_ordinal_type = typename local_matrix_type2::ordinal_type;
    using memory_space       = typename local_matrix_type2::memory_space;
    using results_view       = Kokkos::View<DecisionType*, memory_space>;

#if KOKKOS_VERSION >= 40799
    using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
    using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
    using magnitudeType = typename ATS::magnitudeType;
    using values_view   = Kokkos::View<magnitudeType*, memory_space>;

    const local_matrix_type2 A;
    const local_ordinal_type offset;
    const results_view results;
    mutable values_view values;

   public:
    KOKKOS_INLINE_FUNCTION
    Comparator(const local_matrix_type2& A_, local_ordinal_type rlid_, const results_view& results_, values_view& values_)
      : A(A_)
      , offset(A_.graph.row_map(rlid_))
      , results(results_)
      , values(Kokkos::subview(values_, Kokkos::make_pair(A.graph.row_map(rlid_), A.graph.row_map(rlid_ + 1)))) {
      for (auto i = 0U; i < values.extent(0); ++i) {
        values(i) = get_value_impl(i);
      }
    }

    KOKKOS_FORCEINLINE_FUNCTION
    magnitudeType get_value(size_t x) const {
      return values(x);
    }

    KOKKOS_INLINE_FUNCTION
    bool operator()(size_t x, size_t y) const {
      if (results(offset + x) != UNDECIDED) {
        if (results(offset + y) != UNDECIDED) {
          // does not matter
          return (x < y);
        } else {
          // sort undecided to the right
          return true;
        }
      } else {
        if (results(offset + y) != UNDECIDED) {
          // sort undecided to the right
          return false;
        } else {
          return get_value(x) > get_value(y);
        }
      }
    }

   private:
    KOKKOS_INLINE_FUNCTION
    magnitudeType get_value_impl(size_t x) const {
      return ATS::magnitude(A.values(offset + x) * A.values(offset + x));
    }
  };

  using comparator_type = Comparator<local_matrix_type>;

  KOKKOS_INLINE_FUNCTION
  comparator_type getComparator(local_ordinal_type rlid) const {
    return comparator_type(A, rlid, results, values);
  }
};

/*!
  @class ScaledComparison
  @brief Orders entries of row \f$i\f$ by \f$\frac{|A_{ij}|^2}{|A_{ii}| |A_{jj}|}\f$.
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, Misc::StrengthMeasure measure>
class ScaledComparison {
 public:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_device_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using diag_vec_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_view_type     = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 private:
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
  using magnitudeType = typename ATS::magnitudeType;

  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;
  using values_view = Kokkos::View<magnitudeType*, memory_space>;
  mutable values_view values;

 public:
  ScaledComparison(matrix_type& A_, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , results(results_)
    , values("ScaledComparison::values", A.nnz()) {
    if constexpr ((measure == Misc::SmoothedAggregationMeasure) || (measure == Misc::SignedSmoothedAggregationMeasure)) {
      diagVec        = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixOverlappedDiagonal(A_);
      auto lclDiag2d = diagVec->getLocalViewDevice(Tpetra::Access::ReadOnly);
      diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
    } else if constexpr (measure == Misc::SignedRugeStuebenMeasure) {
      diagVec    = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixMaxMinusOffDiagonal(A_);
      auto lcl2d = diagVec->getLocalViewDevice(Tpetra::Access::ReadOnly);
      diag       = Kokkos::subview(lcl2d, Kokkos::ALL(), 0);
    }
  }

  template <class local_matrix_type2, class diag_view_type2>
  struct Comparator {
   private:
    using scalar_type        = typename local_matrix_type2::value_type;
    using local_ordinal_type = typename local_matrix_type2::ordinal_type;
    using memory_space       = typename local_matrix_type2::memory_space;
    using results_view       = Kokkos::View<DecisionType*, memory_space>;

#if KOKKOS_VERSION >= 40799
    using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
    using ATS  = Kokkos::ArithTraits<scalar_type>;
#endif
    using magnitudeType = typename ATS::magnitudeType;
#if KOKKOS_VERSION >= 40799
    using mATS = KokkosKernels::ArithTraits<magnitudeType>;
#else
    using mATS = Kokkos::ArithTraits<magnitudeType>;
#endif
    using values_view = Kokkos::View<magnitudeType*, memory_space>;

    const local_matrix_type2 A;
    const diag_view_type2 diag;
    const local_ordinal_type rlid;
    const local_ordinal_type offset;
    const results_view results;
    mutable values_view values;

   public:
    KOKKOS_INLINE_FUNCTION
    Comparator(const local_matrix_type2& A_, const diag_view_type2& diag_, const local_ordinal_type rlid_, const results_view& results_, values_view& values_)
      : A(A_)
      , diag(diag_)
      , rlid(rlid_)
      , offset(A_.graph.row_map(rlid_))
      , results(results_)
      , values(Kokkos::subview(values_, Kokkos::make_pair(A.graph.row_map(rlid_), A.graph.row_map(rlid_ + 1)))) {
      for (auto i = 0U; i < values.extent(0); ++i) {
        values(i) = get_value_impl(i);
      }
    }

    KOKKOS_INLINE_FUNCTION
    magnitudeType get_value(size_t x) const {
      return values(x);
    }

    KOKKOS_INLINE_FUNCTION
    bool operator()(size_t x, size_t y) const {
      if (results(offset + x) != UNDECIDED) {
        if (results(offset + y) != UNDECIDED) {
          // does not matter
          return (x < y);
        } else {
          // sort undecided to the right
          return true;
        }
      } else {
        if (results(offset + y) != UNDECIDED) {
          // sort undecided to the right
          return false;
        } else {
          return get_value(x) > get_value(y);
        }
      }
    }

   private:
    KOKKOS_INLINE_FUNCTION
    magnitudeType get_value_impl(size_t x) const {
      if constexpr (measure == Misc::SmoothedAggregationMeasure) {
        auto x_aij    = ATS::magnitude(A.values(offset + x) * A.values(offset + x));
        auto x_aiiajj = ATS::magnitude(diag(rlid) * diag(A.graph.entries(offset + x)));
        return (x_aij / x_aiiajj);
      } else if constexpr (measure == Misc::SignedRugeStuebenMeasure) {
        auto val         = A.values(offset + x);
        auto neg_aij     = -ATS::real(val);
        auto max_neg_aik = ATS::real(diag(rlid));
        auto v           = neg_aij / max_neg_aik;
        if (ATS::real(v) <= mATS::zero()) {
          return -ATS::magnitude(v * v);
        } else {
          return ATS::magnitude(v * v);
        }
      } else if constexpr (measure == Misc::SignedSmoothedAggregationMeasure) {
        auto val                  = A.values(offset + x);
        auto x_aiiajj             = ATS::magnitude(diag(rlid) * diag(A.graph.entries(offset + x)));
        const bool is_nonpositive = ATS::real(val) <= mATS::zero();
        magnitudeType aij2        = ATS::magnitude(val) * ATS::magnitude(val);
        if (!is_nonpositive)
          aij2 = -aij2;
        return (aij2 / x_aiiajj);
      }
    }
  };

  using comparator_type = Comparator<local_matrix_type, diag_view_type>;

  KOKKOS_INLINE_FUNCTION
  comparator_type getComparator(local_ordinal_type rlid) const {
    return comparator_type(A, diag, rlid, results, values);
  }
};

// helper function to allow partial template deduction
template <Misc::StrengthMeasure measure, class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
auto make_comparison_functor(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A_,
                             typename ScaledComparison<Scalar, LocalOrdinal, GlobalOrdinal, Node, measure>::results_view& results_) {
  if constexpr (measure == Misc::UnscaledMeasure) {
    auto functor = UnscaledComparison<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A_, results_);
    return functor;
  } else {
    auto functor = ScaledComparison<Scalar, LocalOrdinal, GlobalOrdinal, Node, measure>(A_, results_);
    return functor;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class DistanceFunctorType>
class UnscaledDistanceLaplacianComparison {
 public:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_device_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using diag_vec_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_view_type     = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 private:
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
  using magnitudeType = typename ATS::magnitudeType;
  using values_view   = Kokkos::View<magnitudeType*, memory_space>;

  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;
  DistanceFunctorType dist2;
  mutable values_view values;

 public:
  UnscaledDistanceLaplacianComparison(matrix_type& A_, DistanceFunctorType& dist2_, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , results(results_)
    , dist2(dist2_)
    , values("UnscaledDistanceLaplacianComparison::values", A.nnz()) {
    // Construct ghosted distance Laplacian diagonal
    diagVec        = DistanceLaplacian::getDiagonal(A_, dist2);
    auto lclDiag2d = diagVec->getLocalViewDevice(Tpetra::Access::ReadOnly);
    diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
  }

  template <class local_matrix_type2, class DistanceFunctorType2, class diag_view_type2>
  struct Comparator {
   private:
    using scalar_type        = typename local_matrix_type2::value_type;
    using local_ordinal_type = typename local_matrix_type2::ordinal_type;
    using memory_space       = typename local_matrix_type2::memory_space;
    using results_view       = Kokkos::View<DecisionType*, memory_space>;

#if KOKKOS_VERSION >= 40799
    using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
    using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
    using magnitudeType = typename ATS::magnitudeType;
    using values_view   = Kokkos::View<magnitudeType*, memory_space>;

    const local_matrix_type2 A;
    const diag_view_type2 diag;
    const DistanceFunctorType2* dist2;
    const local_ordinal_type rlid;
    const local_ordinal_type offset;
    const results_view results;
    mutable values_view values;

    const scalar_type one = ATS::one();

   public:
    KOKKOS_INLINE_FUNCTION
    Comparator(const local_matrix_type2& A_, const diag_view_type2& diag_, const DistanceFunctorType2* dist2_, local_ordinal_type rlid_, const results_view& results_, values_view& values_)
      : A(A_)
      , diag(diag_)
      , dist2(dist2_)
      , rlid(rlid_)
      , offset(A_.graph.row_map(rlid_))
      , results(results_)
      , values(Kokkos::subview(values_, Kokkos::make_pair(A.graph.row_map(rlid_), A.graph.row_map(rlid_ + 1)))) {
      for (auto i = 0U; i < values.extent(0); ++i) {
        values(i) = get_value_impl(i);
      }
    }

    KOKKOS_FORCEINLINE_FUNCTION
    magnitudeType get_value(size_t x) const {
      return values(x);
    }

    KOKKOS_INLINE_FUNCTION
    bool operator()(size_t x, size_t y) const {
      if (results(offset + x) != UNDECIDED) {
        if (results(offset + y) != UNDECIDED) {
          // does not matter
          return (x < y);
        } else {
          // sort undecided to the right
          return true;
        }
      } else {
        if (results(offset + y) != UNDECIDED) {
          // sort undecided to the right
          return false;
        } else {
          return get_value(x) > get_value(y);
        }
      }
    }

   private:
    KOKKOS_INLINE_FUNCTION
    magnitudeType get_value_impl(size_t x) const {
      auto clid = A.graph.entries(offset + x);
      scalar_type val;
      if (rlid != clid) {
        val = -one / dist2->distance2(rlid, clid);
      } else {
        val = diag(rlid);
      }
      auto aij2 = ATS::magnitude(val) * ATS::magnitude(val);  // |a_ij|^2
      return aij2;
    }
  };

  using comparator_type = Comparator<local_matrix_type, DistanceFunctorType, diag_view_type>;

  KOKKOS_INLINE_FUNCTION
  comparator_type getComparator(local_ordinal_type rlid) const {
    return comparator_type(A, diag, &dist2, rlid, results, values);
  }
};

/*!
  @class ScaledDistanceLaplacianComparison
  @brief Orders entries of row \f$i\f$ by \f$\frac{|d_{ij}|^2}{|d_{ii}| |d_{jj}|}\f$ where \f$d_ij\f$ is the distance Laplacian.
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class DistanceFunctorType, Misc::StrengthMeasure measure>
class ScaledDistanceLaplacianComparison {
 public:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_device_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using diag_vec_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_view_type     = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 private:
#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
  using magnitudeType = typename ATS::magnitudeType;

  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;
  DistanceFunctorType dist2;

  using values_view = Kokkos::View<magnitudeType*, memory_space>;
  mutable values_view values;

 public:
  ScaledDistanceLaplacianComparison(matrix_type& A_, DistanceFunctorType& dist2_, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , results(results_)
    , dist2(dist2_)
    , values("ScaledDistanceLaplacianComparison::values", A.nnz()) {
    // Construct ghosted distance Laplacian diagonal
    if constexpr ((measure == Misc::SmoothedAggregationMeasure) || (measure == Misc::SignedSmoothedAggregationMeasure)) {
      diagVec        = DistanceLaplacian::getDiagonal(A_, dist2);
      auto lclDiag2d = diagVec->getLocalViewDevice(Tpetra::Access::ReadOnly);
      diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
    } else if constexpr (measure == Misc::SignedRugeStuebenMeasure) {
      diagVec        = DistanceLaplacian::getMaxMinusOffDiagonal(A_, dist2);
      auto lclDiag2d = diagVec->getLocalViewDevice(Tpetra::Access::ReadOnly);
      diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
    }
  }

  template <class local_matrix_type2, class DistanceFunctorType2, class diag_view_type2>
  struct Comparator {
   private:
    using scalar_type        = typename local_matrix_type2::value_type;
    using local_ordinal_type = typename local_matrix_type2::ordinal_type;
    using memory_space       = typename local_matrix_type2::memory_space;
    using results_view       = Kokkos::View<DecisionType*, memory_space>;

#if KOKKOS_VERSION >= 40799
    using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
    using ATS  = Kokkos::ArithTraits<scalar_type>;
#endif
    using magnitudeType = typename ATS::magnitudeType;
#if KOKKOS_VERSION >= 40799
    using mATS = KokkosKernels::ArithTraits<magnitudeType>;
#else
    using mATS = Kokkos::ArithTraits<magnitudeType>;
#endif

    using values_view = Kokkos::View<magnitudeType*, memory_space>;

    const local_matrix_type2 A;
    const diag_view_type2 diag;
    const DistanceFunctorType2* dist2;
    const local_ordinal_type rlid;
    const local_ordinal_type offset;
    const results_view results;
    mutable values_view values;

    const scalar_type one     = ATS::one();
    const scalar_type zero    = ATS::zero();
    const magnitudeType mzero = mATS::zero();

   public:
    KOKKOS_INLINE_FUNCTION
    Comparator(const local_matrix_type2& A_, const diag_view_type2& diag_, const DistanceFunctorType2* dist2_, local_ordinal_type rlid_, const results_view& results_, values_view& values_)
      : A(A_)
      , diag(diag_)
      , dist2(dist2_)
      , rlid(rlid_)
      , offset(A_.graph.row_map(rlid_))
      , results(results_)
      , values(Kokkos::subview(values_, Kokkos::make_pair(A.graph.row_map(rlid), A.graph.row_map(rlid + 1)))) {
      for (auto i = 0U; i < values.extent(0); ++i) {
        values(i) = get_value_impl(i);
      }
    }

    KOKKOS_FORCEINLINE_FUNCTION
    magnitudeType get_value(size_t x) const {
      return values(x);
    }

    KOKKOS_INLINE_FUNCTION
    bool operator()(size_t x, size_t y) const {
      if (results(offset + x) != UNDECIDED) {
        if (results(offset + y) != UNDECIDED) {
          // does not matter
          return (x < y);
        } else {
          // sort undecided to the right
          return true;
        }
      } else {
        if (results(offset + y) != UNDECIDED) {
          // sort undecided to the right
          return false;
        } else {
          return get_value(x) > get_value(y);
        }
      }
    }

   private:
    KOKKOS_INLINE_FUNCTION
    magnitudeType get_value_impl(size_t x) const {
      auto clid = A.graph.entries(offset + x);
      scalar_type val;
      if (rlid != clid) {
        val = -one / dist2->distance2(rlid, clid);
      } else {
        val = diag(rlid);
      }

      if constexpr (measure == Misc::SmoothedAggregationMeasure) {
        auto aiiajj = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
        auto aij2   = ATS::magnitude(val) * ATS::magnitude(val);                // |a_ij|^2
        return (aij2 / aiiajj);
      } else if constexpr (measure == Misc::SignedRugeStuebenMeasure) {
        auto neg_aij     = -ATS::real(val);
        auto max_neg_aik = ATS::real(diag(rlid));
        auto v           = ATS::magnitude(neg_aij / max_neg_aik);
        if (ATS::real(neg_aij) >= mzero)
          return v * v;
        else
          return -v * v;
      } else if constexpr (measure == Misc::SignedSmoothedAggregationMeasure) {
        auto aiiajj               = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
        const bool is_nonpositive = ATS::real(val) <= mATS::zero();
        magnitudeType aij2        = ATS::magnitude(val) * ATS::magnitude(val);  // |a_ij|^2
        // + |a_ij|^2, if a_ij < 0, - |a_ij|^2 if a_ij >=0
        if (!is_nonpositive)
          aij2 = -aij2;
        return aij2 / aiiajj;
      }
    }
  };

  using comparator_type = Comparator<local_matrix_type, DistanceFunctorType, diag_view_type>;

  KOKKOS_INLINE_FUNCTION
  comparator_type getComparator(local_ordinal_type rlid) const {
    return comparator_type(A, diag, &dist2, rlid, results, values);
  }
};

// helper function to allow partial template deduction
template <Misc::StrengthMeasure measure, class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class DistanceFunctorType>
auto make_dlap_comparison_functor(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A_,
                                  DistanceFunctorType& dist2_,
                                  typename ScaledDistanceLaplacianComparison<Scalar, LocalOrdinal, GlobalOrdinal, Node, DistanceFunctorType, measure>::results_view& results_) {
  if constexpr (measure == Misc::UnscaledMeasure) {
    auto functor = UnscaledDistanceLaplacianComparison<Scalar, LocalOrdinal, GlobalOrdinal, Node, DistanceFunctorType>(A_, dist2_, results_);
    return functor;
  } else {
    auto functor = ScaledDistanceLaplacianComparison<Scalar, LocalOrdinal, GlobalOrdinal, Node, DistanceFunctorType, measure>(A_, dist2_, results_);
    return functor;
  }
}

/*!
  @class CutDropFunctor
  @brief Order each row by a criterion, compare the ratio of values and drop all entries once the ratio is below the threshold.
*/
template <class comparison_type>
class CutDropFunctor {
 private:
  using local_matrix_type  = typename comparison_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS = Kokkos::ArithTraits<scalar_type>;
#endif
  using magnitudeType       = typename ATS::magnitudeType;
  using boundary_nodes_view = Kokkos::View<const bool*, memory_space>;

  local_matrix_type A;
  comparison_type comparison;
  magnitudeType eps;
  results_view results;
  Kokkos::View<local_ordinal_type*, memory_space> index;

 public:
  CutDropFunctor(comparison_type& comparison_, magnitudeType threshold)
    : A(comparison_.A)
    , comparison(comparison_)
    , eps(threshold)
    , results(comparison_.results) {
    index = Kokkos::View<local_ordinal_type*, memory_space>("indices", A.nnz());
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(const local_ordinal_type& rlid) const {
    auto row   = A.rowConst(rlid);
    size_t nnz = row.length;

    auto drop_view       = Kokkos::subview(results, Kokkos::make_pair(A.graph.row_map(rlid), A.graph.row_map(rlid + 1)));
    auto row_permutation = Kokkos::subview(index, Kokkos::make_pair(A.graph.row_map(rlid), A.graph.row_map(rlid + 1)));

    auto comparator = comparison.getComparator(rlid);

#ifdef MUELU_COALESCE_DROP_DEBUG
    {
      Kokkos::printf("SoC:        ");
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        Kokkos::printf("%5f ", comparator.get_value(k));
      }
      Kokkos::printf("\n");
    }
#endif

    for (size_t i = 0; i < nnz; ++i) {
      row_permutation(i) = i;
    }
    Misc::serialHeapSort(row_permutation, comparator);

    size_t keepStart = 0;
    size_t dropStart = nnz;
    // find index where dropping starts
    for (size_t i = 1; i < nnz; ++i) {
      auto const& x = row_permutation(i - 1);
      auto const& y = row_permutation(i);
      if ((drop_view(x) != UNDECIDED) && (drop_view(y) == UNDECIDED))
        keepStart = i;
      if ((drop_view(x) != UNDECIDED) || (drop_view(y) != UNDECIDED))
        continue;
      magnitudeType x_aij = comparator.get_value(x);
      magnitudeType y_aij = comparator.get_value(y);
      if (eps * eps * x_aij > y_aij) {
        if (i < dropStart) {
          dropStart = i;
        }
      }
    }

    // drop everything to the right of where values stop passing threshold
    for (size_t i = keepStart; i < nnz; ++i) {
      drop_view(row_permutation(i)) = Kokkos::max(dropStart <= i ? DROP : KEEP, drop_view(row_permutation(i)));
    }
  }
};

}  // namespace MueLu::CutDrop

#endif
