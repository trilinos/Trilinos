//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef _KOKKOSKERNELS_SIMPLEUTILS_HPP
#define _KOKKOSKERNELS_SIMPLEUTILS_HPP
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <type_traits>

#define KOKKOSKERNELS_MACRO_MIN(x, y) ((x) < (y) ? (x) : (y))
#define KOKKOSKERNELS_MACRO_MAX(x, y) ((x) < (y) ? (y) : (x))
#define KOKKOSKERNELS_MACRO_ABS(x) Kokkos::ArithTraits<typename std::decay<decltype(x)>::type>::abs(x)

namespace KokkosKernels {

namespace Impl {

template <class ViewType>
class SquareRootFunctor {
 public:
  typedef typename ViewType::execution_space execution_space;
  typedef typename ViewType::size_type size_type;

  SquareRootFunctor(const ViewType &theView) : theView_(theView) {}

  KOKKOS_INLINE_FUNCTION void operator()(const size_type i) const {
    typedef typename ViewType::value_type value_type;
    theView_(i) = Kokkos::ArithTraits<value_type>::sqrt(theView_(i));
  }

 private:
  ViewType theView_;
};

template <typename view_t>
struct ExclusiveParallelPrefixSum {
  typedef typename view_t::value_type value_type;
  view_t array_sum;
  ExclusiveParallelPrefixSum(view_t arr_) : array_sum(arr_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, value_type &update, const bool final) const {
    value_type val = (ii == array_sum.extent(0) - 1) ? value_type(0) : array_sum(ii);
    if (final) {
      array_sum(ii) = value_type(update);
    }
    update += val;
  }
};

template <typename array_type>
struct InclusiveParallelPrefixSum {
  typedef typename array_type::value_type idx;
  array_type array_sum;
  InclusiveParallelPrefixSum(array_type arr_) : array_sum(arr_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, size_t &update, const bool final) const {
    update += array_sum(ii);
    if (final) {
      array_sum(ii) = idx(update);
    }
  }
};

/***
 * \brief Function performs the exclusive parallel prefix sum. That is each
 * entry holds the sum until itself.
 * \param exec: the execution space instance on which to run
 * \param num_elements: size of the array
 * \param arr: the array for which the prefix sum will be performed.
 */
template <typename MyExecSpace, typename view_t>
inline void kk_exclusive_parallel_prefix_sum(const MyExecSpace &exec, typename view_t::value_type num_elements,
                                             view_t arr) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_scan("KokkosKernels::Common::PrefixSum", my_exec_space(exec, 0, num_elements),
                        ExclusiveParallelPrefixSum<view_t>(arr));
}

/***
 * \brief Function performs the exclusive parallel prefix sum. That is each
 * entry holds the sum until itself.
 * \param num_elements: size of the array
 * \param arr: the array for which the prefix sum will be performed.
 */
template <typename MyExecSpace, typename view_t>
inline void kk_exclusive_parallel_prefix_sum(typename view_t::value_type num_elements, view_t arr) {
  kk_exclusive_parallel_prefix_sum(MyExecSpace(), num_elements, arr);
}

/***
 * \brief Function performs the exclusive parallel prefix sum. That is each
 * entry holds the sum until itself. This version also returns the final sum
 * equivalent to the sum-reduction of arr before doing the scan.
 * \param exec: the execution space instance on which to run
 * \param num_elements: size of the array
 * \param arr: the array for which the prefix sum will be performed.
 * \param finalSum: will be set to arr[num_elements - 1] after computing the
 * prefix sum.
 */
template <typename MyExecSpace, typename view_t>
inline void kk_exclusive_parallel_prefix_sum(const MyExecSpace &exec, typename view_t::value_type num_elements,
                                             view_t arr, typename view_t::non_const_value_type &finalSum) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_scan("KokkosKernels::Common::PrefixSum", my_exec_space(exec, 0, num_elements),
                        ExclusiveParallelPrefixSum<view_t>(arr), finalSum);
}

/***
 * \brief Function performs the exclusive parallel prefix sum. That is each
 * entry holds the sum until itself. This version also returns the final sum
 * equivalent to the sum-reduction of arr before doing the scan.
 * \param num_elements: size of the array
 * \param arr: the array for which the prefix sum will be performed.
 * \param finalSum: will be set to arr[num_elements - 1] after computing the
 * prefix sum.
 */
template <typename MyExecSpace, typename view_t>
inline void kk_exclusive_parallel_prefix_sum(typename view_t::value_type num_elements, view_t arr,
                                             typename view_t::non_const_value_type &finalSum) {
  kk_exclusive_parallel_prefix_sum(MyExecSpace(), num_elements, arr, finalSum);
}

///
/// \brief Function performs the inclusive parallel prefix sum. That is each
///        entry holds the sum until itself including itself.
/// \param my_exec_space: The execution space instance
/// \param num_elements: size of the array
/// \param arr: the array for which the prefix sum will be performed.
///
template <typename MyExecSpace, typename forward_array_type>
void kk_inclusive_parallel_prefix_sum(MyExecSpace my_exec_space, typename forward_array_type::value_type num_elements,
                                      forward_array_type arr) {
  typedef Kokkos::RangePolicy<MyExecSpace> range_policy_t;
  Kokkos::parallel_scan("KokkosKernels::Common::PrefixSum", range_policy_t(my_exec_space, 0, num_elements),
                        InclusiveParallelPrefixSum<forward_array_type>(arr));
}

///
/// \brief Function performs the inclusive parallel prefix sum. That is each
///        entry holds the sum until itself including itself.
/// \param num_elements: size of the array
/// \param arr: the array for which the prefix sum will be performed.
///
template <typename MyExecSpace, typename forward_array_type>
void kk_inclusive_parallel_prefix_sum(typename forward_array_type::value_type num_elements, forward_array_type arr) {
  MyExecSpace my_exec_space;
  return kk_inclusive_parallel_prefix_sum(my_exec_space, num_elements, arr);
}

template <typename view_t>
struct ReductionFunctor {
  view_t array_sum;
  ReductionFunctor(view_t arr_) : array_sum(arr_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, typename view_t::value_type &update) const { update += array_sum(ii); }
};

template <typename view_t>
struct ReductionFunctor2 {
  view_t array_sum;
  ReductionFunctor2(view_t arr_) : array_sum(arr_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, size_t &update) const { update += array_sum(ii); }
};

template <typename view_t, typename view2_t>
struct DiffReductionFunctor {
  view_t array_begins;
  view2_t array_ends;
  DiffReductionFunctor(view_t begins, view2_t ends) : array_begins(begins), array_ends(ends) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, typename view_t::non_const_value_type &update) const {
    update += (array_ends(ii) - array_begins(ii));
  }
};

template <typename view_t, typename view2_t, typename MyExecSpace>
inline void kk_reduce_diff_view(size_t num_elements, view_t smaller, view2_t bigger,
                                typename view_t::non_const_value_type &reduction) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_reduce("KokkosKernels::Common::ReduceDiffView", my_exec_space(0, num_elements),
                          DiffReductionFunctor<view_t, view2_t>(smaller, bigger), reduction);
}

template <typename it>
struct DiffReductionFunctorP {
  const it *array_begins;
  const it *array_ends;
  DiffReductionFunctorP(const it *begins, const it *ends) : array_begins(begins), array_ends(ends) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, it &update) const { update += (array_ends[ii] - array_begins[ii]); }
};

template <typename it, typename MyExecSpace>
inline void kkp_reduce_diff_view(const size_t num_elements, const it *smaller, const it *bigger, it &reduction) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_reduce("KokkosKernels::Common::ReduceDiffView", my_exec_space(0, num_elements),
                          DiffReductionFunctorP<it>(smaller, bigger), reduction);
}

/***
 * \brief Function performs the a reduction
 * until itself.
 * \param num_elements: size of the array
 * \param arr: the array for which the prefix sum will be performed.
 */
template <typename view_t, typename MyExecSpace>
inline void kk_reduce_view(size_t num_elements, view_t arr, typename view_t::value_type &reduction) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_reduce("KokkosKernels::Common::ReduceView", my_exec_space(0, num_elements),
                          ReductionFunctor<view_t>(arr), reduction);
}

template <typename view_t, typename MyExecSpace>
inline void kk_reduce_view2(size_t num_elements, view_t arr, size_t &reduction) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_reduce("KokkosKernels::Common::ReduceView2", my_exec_space(0, num_elements),
                          ReductionFunctor2<view_t>(arr), reduction);
}

template <typename view_type1, typename view_type2,
          typename eps_type = typename Kokkos::ArithTraits<typename view_type2::non_const_value_type>::mag_type>
struct IsIdenticalFunctor {
  view_type1 view1;
  view_type2 view2;
  eps_type eps;

  IsIdenticalFunctor(view_type1 view1_, view_type2 view2_, eps_type eps_) : view1(view1_), view2(view2_), eps(eps_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, size_t &is_equal) const {
    typedef typename view_type2::non_const_value_type val_type;
    typedef Kokkos::ArithTraits<val_type> KAT;
    typedef typename KAT::mag_type mag_type;
    const mag_type val_diff = KAT::abs(view1(i) - view2(i));

    if (val_diff > eps) {
      is_equal += 1;
    }
  }
};

template <typename view_type1, typename view_type2, typename eps_type, typename MyExecSpace>
bool kk_is_identical_view(view_type1 view1, view_type2 view2, eps_type eps) {
  if (view1.extent(0) != view2.extent(0)) {
    return false;
  }

  size_t num_elements = view1.extent(0);

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  size_t issame = 0;
  Kokkos::parallel_reduce("KokkosKernels::Common::IsIdenticalView", my_exec_space(0, num_elements),
                          IsIdenticalFunctor<view_type1, view_type2, eps_type>(view1, view2, eps), issame);
  MyExecSpace().fence();
  if (issame > 0) {
    return false;
  } else {
    return true;
  }
}

template <typename view_type1, typename view_type2,
          typename eps_type = typename Kokkos::ArithTraits<typename view_type2::non_const_value_type>::mag_type>
struct IsRelativelyIdenticalFunctor {
  view_type1 view1;
  view_type2 view2;
  eps_type eps;

  IsRelativelyIdenticalFunctor(view_type1 view1_, view_type2 view2_, eps_type eps_)
      : view1(view1_), view2(view2_), eps(eps_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, size_t &num_diffs) const {
    typedef typename view_type2::non_const_value_type val_type;
    typedef Kokkos::ArithTraits<val_type> KAT;
    typedef typename KAT::mag_type mag_type;
    typedef Kokkos::ArithTraits<mag_type> KATM;

    mag_type val_diff = KATM::zero();
    if (KAT::abs(view1(i)) > mag_type(eps) || KAT::abs(view2(i)) > mag_type(eps)) {
      val_diff = KAT::abs(view1(i) - view2(i)) / (KAT::abs(view1(i)) + KAT::abs(view2(i)));
    }

    if (val_diff > mag_type(eps)) {
      Kokkos::printf(
          "Values at index %d, %.6f + %.6fi and %.6f + %.6fi, differ too much "
          "(eps = %e, rel err = %e)\n",
          (int)i, KAT::real(view1(i)), KAT::imag(view1(i)), KAT::real(view2(i)), KAT::imag(view2(i)), eps, val_diff);
      num_diffs++;
    }
  }
};

template <typename view_type1, typename view_type2, typename eps_type, typename MyExecSpace>
bool kk_is_relatively_identical_view(view_type1 view1, view_type2 view2, eps_type eps) {
  if (view1.extent(0) != view2.extent(0)) {
    return false;
  }

  size_t num_elements = view1.extent(0);

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  size_t numDifferences = 0;
  Kokkos::parallel_reduce("KokkosKernels::Common::IsRelativelyIdenticalView", my_exec_space(0, num_elements),
                          IsRelativelyIdenticalFunctor<view_type1, view_type2, eps_type>(view1, view2, eps),
                          numDifferences);
  return numDifferences == 0;
}

template <typename view_type>
struct ReduceMaxFunctor {
  view_type view_to_reduce;
  typedef typename view_type::non_const_value_type value_type;
  const value_type min_val;
  ReduceMaxFunctor(view_type view_to_reduce_)
      : view_to_reduce(view_to_reduce_), min_val((std::numeric_limits<value_type>::lowest())) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, value_type &max_reduction) const {
    value_type val = view_to_reduce(i);
    if (max_reduction < val) {
      max_reduction = val;
    }
  }
  KOKKOS_INLINE_FUNCTION
  void join(value_type &dst, const value_type &src) const {
    if (dst < src) {
      dst = src;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type &dst) const {
    // The identity under max is -Inf.
    // Kokkos does not come with a portable way to access
    // floating -point Inf and NaN. Trilinos does , however;
    // see Kokkos :: ArithTraits in the Tpetra package.
    dst = min_val;
  }
};

template <typename view_type, typename MyExecSpace>
void kk_view_reduce_max(const MyExecSpace &exec, size_t num_elements, view_type view_to_reduce,
                        typename view_type::non_const_value_type &max_reduction) {
  typedef Kokkos::RangePolicy<MyExecSpace> policy_t;
  Kokkos::parallel_reduce("KokkosKernels::Common::ReduceMax", policy_t(exec, 0, num_elements),
                          ReduceMaxFunctor<view_type>(view_to_reduce), max_reduction);
}

template <typename view_type, typename MyExecSpace>
void kk_view_reduce_max(size_t num_elements, view_type view_to_reduce,
                        typename view_type::non_const_value_type &max_reduction) {
  kk_view_reduce_max(MyExecSpace(), num_elements, view_to_reduce, max_reduction);
}

// xorshift hash/pseudorandom function (supported for 32- and 64-bit integer
// types only)
template <typename Value>
KOKKOS_FORCEINLINE_FUNCTION Value xorshiftHash(Value v) {
  static_assert(std::is_unsigned<Value>::value, "xorshiftHash: value must be an unsigned integer type");
  uint64_t x = v;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  return std::is_same<Value, uint32_t>::value ? static_cast<Value>((x * 2685821657736338717ULL - 1) >> 16)
                                              : static_cast<Value>(x * 2685821657736338717ULL - 1);
}

struct ViewHashFunctor {
  ViewHashFunctor(const uint8_t *data_) : data(data_) {}

  KOKKOS_INLINE_FUNCTION void operator()(size_t i, uint32_t &lhash) const {
    // Compute a hash/digest of both the index i, and data[i]. Then add that to
    // overall hash.
    uint32_t x = uint32_t(i);
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    x ^= uint32_t(data[i]);
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    lhash += x;
  }

  const uint8_t *data;
};

/// \brief Compute a hash of a view.
/// \param v: the view to hash. Must be contiguous, and its element type must
/// not contain any padding bytes.
template <typename View>
uint32_t hashView(const View &v) {
  assert(v.span_is_contiguous());
  // Note: This type trait is supposed to be part of C++17,
  // but it's not defined on Intel 19 (with GCC 7.2.0 standard library).
  // So just check if it's available before using.
#ifdef __cpp_lib_has_unique_object_representations
  static_assert(std::has_unique_object_representations<typename View::non_const_value_type>::value,
                "KokkosKernels::Impl::hashView: the view's element type must "
                "not have any padding bytes.");
#endif
  size_t nbytes = v.span() * sizeof(typename View::value_type);
  uint32_t h;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<typename View::execution_space, size_t>(0, nbytes),
                          ViewHashFunctor(reinterpret_cast<const uint8_t *>(v.data())), h);
  return h;
}

template <typename V>
struct SequentialFillFunctor {
  using size_type = typename V::size_type;
  using val_type  = typename V::non_const_value_type;
  SequentialFillFunctor(const V &v_, val_type start_) : v(v_), start(start_) {}
  KOKKOS_INLINE_FUNCTION void operator()(size_type i) const { v(i) = start + (val_type)i; }
  V v;
  val_type start;
};

template <typename ExecSpace, typename V>
void sequential_fill(const ExecSpace &exec, const V &v, typename V::non_const_value_type start = 0) {
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(exec, 0, v.extent(0)), SequentialFillFunctor<V>(v, start));
}

template <typename V>
void sequential_fill(const V &v, typename V::non_const_value_type start = 0) {
  sequential_fill(typename V::execution_space(), v, start);
}

}  // namespace Impl
}  // namespace KokkosKernels

#endif
