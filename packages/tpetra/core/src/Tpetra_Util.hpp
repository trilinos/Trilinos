// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_UTIL_HPP
#define TPETRA_UTIL_HPP

/*!
  \file Tpetra_Util.hpp
  \brief Stand-alone utility functions and macros.

  Tpetra_Util contains utility functions macros that are used
  throughout Tpetra, by many classes and functions. They are placed
  here so that they can be updated and maintained in a single spot.
*/

#include "Tpetra_ConfigDefs.hpp"
#include "Kokkos_DualView.hpp"
#include "KokkosCompat_View.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_Utils.hpp"
#include <numeric>
#include <algorithm>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream>

#if defined(HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS)
/// \brief Print or throw an efficency warning.
///
/// This macro is only for use by Tpetra developers.  It is only to be
/// used in the implementation of a Tpetra class' instance method.
///
/// If HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS is defined, print the
/// given message to std::cerr, along with other useful information.
///
/// This macro must be called in an instance method of a class.
/// (Alas, C++ gives us no way to search the scope for a list of
/// defined names.  Otherwise, the macro could search for "this" to
/// determine whether it is in a class' instance method.)
///
/// Macro arguments:
///
/// throw_exception_test: Boolean expression to evaluate.  If true,
/// this macro will trigger the efficiency warning.  The test will be
/// evaluated at most once.  Nevertheless, the test should not have
/// side effects, since it will not be evaluated at all if
/// HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS is defined.
///
/// msg: The message to include in the warning.  The warning also
/// includes the name of the class, and other useful information.
///
#define TPETRA_EFFICIENCY_WARNING(throw_exception_test,msg)  \
{ \
  const bool tpetraEfficiencyWarningTest = (throw_exception_test); \
  if (tpetraEfficiencyWarningTest) { \
    std::ostringstream errStream; \
    errStream << Teuchos::typeName(*this) << ":" << std::endl; \
    errStream << "Efficiency warning: " << #throw_exception_test << std::endl; \
    errStream << msg; \
    std::string err = errStream.str(); \
    if (TPETRA_PRINTS_EFFICIENCY_WARNINGS && tpetraEfficiencyWarningTest) { \
      std::cerr << err << std::endl; \
    } \
  } \
}
#else
/// \brief Print or throw an efficency warning.
///
/// This macro is only for use by Tpetra developers.  It is only to be
/// used in the implementation of a Tpetra class' instance method.
///
/// If HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS is defined, print the
/// given message to std::cerr, along with other useful information.
///
/// This macro must be called in an instance method of a class.
/// (Alas, C++ gives us no way to search the scope for a list of
/// defined names.  Otherwise, the macro could search for "this" to
/// determine whether it is in a class' instance method.)
///
/// Macro arguments:
///
/// throw_exception_test: Boolean expression to evaluate.  If true,
/// this macro will trigger the efficiency warning.  The test will be
/// evaluated at most once.  Nevertheless, the test should not have
/// side effects, since it will not be evaluated at all if neither
/// HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS is defined.
///
/// msg: The message to include in the warning.  The warning also
/// includes the name of the class, and other useful information.
///
#define TPETRA_EFFICIENCY_WARNING(throw_exception_test,msg)
#endif

// handle an abuse warning, according to HAVE_TPETRA_THROW_ABUSE_WARNINGS and HAVE_TPETRA_PRINT_ABUSE_WARNINGS
#if defined(HAVE_TPETRA_THROW_ABUSE_WARNINGS) || defined(HAVE_TPETRA_PRINT_ABUSE_WARNINGS)
//! Handle an abuse warning, according to HAVE_TPETRA_THROW_ABUSE_WARNINGS and HAVE_TPETRA_PRINT_ABUSE_WARNINGS
#define TPETRA_ABUSE_WARNING(throw_exception_test,Exception,msg)                               \
{                                                                                              \
  std::ostringstream errStream;                                                                \
  errStream << Teuchos::typeName(*this) << msg;                                                \
  std::string err = errStream.str();                                                           \
  if (TPETRA_PRINTS_ABUSE_WARNINGS && (throw_exception_test)) {                                \
    std::cerr << err << std::endl;                                                             \
  }                                                                                            \
  TEUCHOS_TEST_FOR_EXCEPTION(TPETRA_THROWS_ABUSE_WARNINGS && (throw_exception_test), Exception, err);  \
}
#else
//! Handle an abuse warning, according to HAVE_TPETRA_THROW_ABUSE_WARNINGS and HAVE_TPETRA_PRINT_ABUSE_WARNINGS
#define TPETRA_ABUSE_WARNING(throw_exception_test,Exception,msg)
#endif

/// \brief Test for exception, with reduction over the given communicator.
///
/// This is like Teuchos' TEUCHOS_TEST_FOR_EXCEPTION macro, except
/// that it performs an all-reduce over the given communicator to
/// ensure that all processes throw the exception if the test is true
/// on at least one process.
///
/// Macro arguments:
///
/// throw_exception_test: Boolean expression to evaluate.  If true on
/// at least one calling process in the given communicator, this macro
/// will throw an exception of the given type on all processes in the
/// communicator.  The exception message may differ on different
/// processes.  The expression will only be evaluated once on each
/// process.
///
/// Exception: The type of exception to throw, if throw_exception_test
/// evaluates to true on at least one process in the given
/// communicator.  The Exception should be a subclass of
/// std::exception.
///
/// msg: The message to include in the warning.  The warning also
/// includes the name of the class, and the maximum process rank on
/// which the test evaluated to true.  The message may be different on
/// different processes.
///
/// comm: The communicator (instance of a subclass of
/// Teuchos::Comm<int>) over which to test.  This must evaluate to a
/// class instance or reference, not a Teuchos::RCP of an instance.
#define SHARED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm) \
{ \
    using Teuchos::outArg; \
    const int lcl_throw_exception = (throw_exception_test) ? Teuchos::rank(comm)+1 : 0; \
    int gbl_throw; \
    Teuchos::reduceAll(comm,Teuchos::REDUCE_MAX,lcl_throw_exception,outArg(gbl_throw)); \
    TEUCHOS_TEST_FOR_EXCEPTION(gbl_throw,Exception,  \
      msg << " Failure on at least one process, including process " << gbl_throw-1 << "." << std::endl); \
}

namespace Tpetra {

  /// \brief Efficiently insert or replace an entry in an std::map.
  ///
  /// \tparam MapType Specialization of std::map
  /// \tparam KeyArgType Type of keys of the std::map
  /// \tparam ValueArgType Type of values of the std::map
  ///
  /// This function is taken from Scott Meyers' "Effective STL", Item
  /// 24.  If the given std::map m already contains an entry with key
  /// k, replace its value with the given value v.  Otherwise, insert
  /// (k,v) into the std::map.  In both cases, return an iterator that
  /// points to the inserted or updated entry.
  template<typename MapType, typename KeyArgType, typename ValueArgType>
  typename MapType::iterator efficientAddOrUpdate(MapType& m,
                          const KeyArgType & k,
                          const ValueArgType & v)
  {
    typename MapType::iterator lb = m.lower_bound(k);
    if(lb != m.end() && !(m.key_comp()(k, lb->first))) {
      lb->second = v;
      return(lb);
    }
    else {
      typedef typename MapType::value_type MVT;
      return(m.insert(lb, MVT(k, v)));
    }
  }

  /// \namespace SortDetails
  /// \brief Implementation details of sort routines used by Tpetra.
  ///
  /// \warning The contents of this namespace are an implementation
  ///   detail of Tpetra.  Tpetra users should not use this namespace
  ///   or any of its contents, as they may change at any time without
  ///   deprecation warnings.
  namespace SortDetails {

  /**
   * \brief Determines whether or not a random access sequence is already sorted.
   *
   * @param first An iterator pointing to the beginning of the sequence.
   * @param last An iterator pointing to the end of the sequence.
   *
   * @return True if the sequence is already sorted, else false.
   */
  template<class IT1>
  bool isAlreadySorted(const IT1& first, const IT1& last){
    typedef typename std::iterator_traits<IT1>::difference_type DT;
    DT myit = Teuchos::OrdinalTraits<DT>::one();
    const DT sz  = last - first;
    for(;myit < sz; ++myit){
      if(first[myit] < first[myit-1]){
        return false;
      }
    }
    return true;
  }

  /**
   * \brief Determines the pivot point as part of the quicksort routine.
   *
   * @param first An iterator pointing to the beginning of the array segment
   * we are trying to find a pivot for.
   * @param last An iterator pointing to the ending of the array segment
   * we are trying to find a pivot for.
   * @return The pivot point for the given array segment.
   */
  template<class IT>
  IT getPivot(const IT& first, const IT& last){
    IT pivot(first+(last-first)/2);
    if(*first<=*pivot && *(last-1)<=*first) pivot=first;
    else if(*(last-1)<=*pivot && *first<= *(last-1)) pivot = last-1;
    return pivot;
  }

  /**
   * \brief Partition operation for \c quicksort2().
   *
   * This is a helper routine for \c quicksort2().  It performs the
   * partition step in Quicksort.  All of the input argument iterators
   * are random access iterators.
   *
   * @param first1 An iterator pointing to the beginning of the first array.
   * @param last1 An iterator pointing to the end of the first array.
   * @param first2 An iterator pointing to the beginning of the second array.
   * @param last2 An iterator pointing to the end of the second array.
   * @param pivot A pivot point calculated by the pivot function.
   * @return An iterator pointing to where the partition should be made.
   */
  template<class IT1, class IT2>
  IT1 partition2(
    const IT1& first1,
    const IT1& last1,
    const IT2& first2,
    const IT2& last2,
    const IT1& pivot)
  {
    typename std::iterator_traits<IT1>::value_type piv(*pivot);
    std::swap(*pivot, *(last1-1));
    std::swap(first2[(pivot-first1)], *(last2-1));
    IT1 store1=first1;
    for(IT1 it=first1; it!=last1-1; ++it){
      if(*it<=piv){
        std::swap(*store1, *it);
        std::swap(first2[(store1-first1)], first2[(it-first1)]);
        ++store1;
      }
    }
    std::swap(*(last1-1), *store1);
    std::swap(*(last2-1), first2[store1-first1]);
    return store1;
  }

  /**
   * \brief Partition operation for \c quicksort3().
   *
   * This is a helper routine for \c quicksort3().  It performs the
   * partition step in Quicksort.  All of the input argument iterators
   * are random access iterators.
   *
   * @param first1 An iterator pointing to the beginning of the first array.
   * @param last1 An iterator pointing to the end of the first array.
   * @param first2 An iterator pointing to the beginning of the second array.
   * @param last2 An iterator pointing to the end of the second array.
   * @param first3 An iterator pointing to the beginning of the third array.
   * @param last3 An iterator pointing to the end of the third array.
   * @param pivot A pivot point calculated by the pivot function.
   * @return An iterator pointing to where the partition should be made.
   */
  template<class IT1, class IT2, class IT3>
  IT1 partition3(
    const IT1& first1,
    const IT1& last1,
    const IT2& first2,
    const IT2& last2,
    const IT3& first3,
    const IT3& last3,
    const IT1& pivot)
  {
    typename std::iterator_traits<IT1>::value_type piv(*pivot);
    std::swap(*pivot, *(last1-1));
    std::swap(first2[(pivot-first1)], *(last2-1));
    std::swap(first3[(pivot-first1)], *(last3-1));
    IT1 store1=first1;
    for(IT1 it=first1; it!=last1-1; ++it){
      if(*it<=piv){
        std::swap(*store1, *it);
        std::swap(first2[(store1-first1)], first2[(it-first1)]);
        std::swap(first3[(store1-first1)], first3[(it-first1)]);
        ++store1;
      }
    }
    std::swap(*(last1-1), *store1);
    std::swap(*(last2-1), first2[store1-first1]);
    std::swap(*(last3-1), first3[store1-first1]);
    return store1;
  }

  /**
   * \brief Sort the first array using Quicksort, and apply the resulting permutation to the second array.
   *
   * All of the input argument iterators are random access iterators.
   *
   * @param first1 An iterator pointing to the beginning of the first array.
   * @param last1 An iterator pointing to the end (exclusive) of the first array.
   * @param first2 An iterator pointing to the beginning of the second array.
   * @param last2 An iterator pointing to the end of the second array.
   */
  template<class IT1, class IT2>
  void quicksort2(
    const IT1& first1,
    const IT1& last1,
    const IT2& first2,
    const IT2& last2)
  {
    typedef typename std::iterator_traits<IT1>::difference_type DT;
    DT DT1 = Teuchos::OrdinalTraits<DT>::one();
    if(last1-first1 > DT1){
      IT1 pivot = getPivot(first1, last1);
      pivot = partition2(first1, last1, first2, last2, pivot);
      quicksort2(first1, pivot, first2, first2+(pivot-first1));
      quicksort2(pivot+1, last1, first2+(pivot-first1)+1, last2);
    }
  }

  /**
   * \brief Sort the first array using Quicksort, and apply the resulting permutation to the second and third arrays.
   *
   * All of the input argument iterators are random access iterators.
   *
   * @param first1 An iterator pointing to the beginning of the first array.
   * @param last1 An iterator pointing to the end of the first array.
   * @param first2 An iterator pointing to the beginning of the second array.
   * @param last2 An iterator pointing to the end of the second array.
   * @param first3 An iterator pointing to the beginning of the third array.
   * @param last3 An iterator pointing to the end of the third array.
   */
  template<class IT1, class IT2, class IT3>
  void quicksort3(
    const IT1& first1,
    const IT1& last1,
    const IT2& first2,
    const IT2& last2,
    const IT3& first3,
    const IT3& last3)
  {
    typedef typename std::iterator_traits<IT1>::difference_type DT;
    DT DT1 = Teuchos::OrdinalTraits<DT>::one();
    if(last1-first1 > DT1){
      IT1 pivot = getPivot(first1, last1);
      pivot = partition3(first1, last1, first2, last2, first3, last3, pivot);
      quicksort3(first1, pivot, first2, first2+(pivot-first1), first3, first3+(pivot-first1));
      quicksort3(pivot+1, last1, first2+(pivot-first1)+1, last2, first3+(pivot-first1)+1, last3);
    }
  }

  /**
   * \brief Sort the first array using shell sort, and apply the resulting
   * permutation to the second and third arrays.
   *
   * All of the input argument iterators are random access iterators.
   */
  template<class IT1, class IT2, class IT3>
  void sh_sort3(
    const IT1& first1,
    const IT1& last1,
    const IT2& first2,
    const IT2& /* last2 */,
    const IT3& first3,
    const IT3& /* last3 */)
   {
        typedef typename std::iterator_traits<IT1>::difference_type DT;
        DT n = last1 - first1;
        DT m = n / 2;
        DT z = Teuchos::OrdinalTraits<DT>::zero();
        while (m > z)
        {
            DT max = n - m;
            for (DT j = 0; j < max; j++)
            {
                for (DT k = j; k >= 0; k-=m)
                {
                    if (first1[k+m] >= first1[k])
                        break;
                    std::swap(first1[k+m], first1[k]);
                    std::swap(first2[k+m], first2[k]);
                    std::swap(first3[k+m], first3[k]);
                }
            }
            m = m/2;
        }
   }

  /**
   * \brief Sort the first array using shell sort, and apply the resulting
   * permutation to the second array.
   *
   * All of the input argument iterators are random access iterators.
   */
  template<class IT1, class IT2>
  void sh_sort2(
    const IT1& first1,
    const IT1& last1,
    const IT2& first2,
    const IT2& /* last2 */)
   {
        typedef typename std::iterator_traits<IT1>::difference_type DT;
        DT n = last1 - first1;
        DT m = n / 2;
        DT z = Teuchos::OrdinalTraits<DT>::zero();
        while (m > z)
        {
            DT max = n - m;
            for (DT j = 0; j < max; j++)
            {
                for (DT k = j; k >= 0; k-=m)
                {
                    if (first1[k+m] >= first1[k])
                        break;
                    std::swap(first1[k+m], first1[k]);
                    std::swap(first2[k+m], first2[k]);
                }
            }
            m = m/2;
        }
   }

  /**
   * \brief Compute the permutation of the ordering that 
   * sorts the array using a stable sort.
   */
  template<typename IT1>
  std::vector<typename std::iterator_traits<IT1>::difference_type> sort_indexes(IT1 first, IT1 last) {

    typedef typename std::iterator_traits<IT1>::difference_type DT;

    DT length = last - first;
    std::vector<DT> idx(length);
    std::iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(),
        [&first](size_t i1, size_t i2) {return first[i1] < first[i2];});

    return idx;
  }

  /**
   * \brief Apply a permutation of the ordering of an array in place.
   */
  template<typename IT1, typename IT2>
  void
  apply_permutation(
      IT1 first,
      IT1 last,
      IT2 indices)
  {
    typedef typename std::iterator_traits<IT1>::difference_type DT;
    typedef typename std::iterator_traits<IT1>::value_type T;

    DT n = last - first;

    Kokkos::View<T*, Kokkos::HostSpace> tmp(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp"), n);

    Kokkos::parallel_for("apply_permutation_1",
                        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
                        [=](const int j) {
                            tmp(j) = first[indices[j]];
    });
    Kokkos::parallel_for("apply_permutation_2",
                        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
                        [=](const int j) {
                            first[j] = tmp(j);
    });
  }

  /**
   * \brief Sort the first array using std sort, and apply the resulting
   * permutation to the second array.
   *
   * All of the input argument iterators are random access iterators.
   */
  template<class IT1, class IT2>
  void std_sort2(const IT1& first1,
    const IT1& last1,
    const IT2& first2,
    const IT2& last2)
   {
    auto indices = sort_indexes(first1, last1);
    apply_permutation(first1, last1, indices);
    apply_permutation(first2, last2, indices);
  }

  /**
   * \brief Sort the first array using std sort, and apply the resulting
   * permutation to the second and third arrays.
   *
   * All of the input argument iterators are random access iterators.
   */
  template<class IT1, class IT2, class IT3>
  void std_sort3(
    const IT1& first1,
    const IT1& last1,
    const IT2& first2,
    const IT2& last2,
    const IT3& first3,
    const IT3& last3)
   {
    auto indices = sort_indexes(first1, last1);
    apply_permutation(first1, last1, indices);
    apply_permutation(first2, last2, indices);
    apply_permutation(first3, last3, indices);
   }

  } //end namespace SortDetails


  /**
   * \brief Sort the first array, and apply the resulting permutation to the second array.
   *
   * Sort the values in the first array (represented by the exclusive
   * iterator range first1,last1) in ascending order.  Apply the
   * permutation resulting from the sort to the second array
   * (represented by a starting iterator first2).
   *
   * @param first1 A random access iterator pointing to the beginning
   *   of the first array.
   * @param last1 A random access iterator pointing to the end
   *   (exclusive) of the first array.
   * @param first2 A random access iterator pointing to the beginning
   *   of the second array.  The second array must have no fewer
   *   elements than the first array.  If the first array has N
   *   elements, then the permutation will only be applied to the
   *   first N elements of the second array.
   */
  template<class IT1, class IT2>
  void sort2(const IT1 &first1, const IT1 &last1, const IT2 &first2, const bool stableSort=false) {
    // Quicksort uses best-case N log N time whether or not the input
    // data is sorted.  However, the common case in Tpetra is that the
    // input data are sorted, so we first check whether this is the
    // case before reverting to Quicksort.
    if(SortDetails::isAlreadySorted(first1, last1)){
      return;
    }
    if(stableSort)
      SortDetails::std_sort2(first1, last1, first2, first2+(last1-first1));
    else
      SortDetails::sh_sort2(first1, last1, first2, first2+(last1-first1));
#ifdef HAVE_TPETRA_DEBUG
    if(!SortDetails::isAlreadySorted(first1, last1)){
      std::cout << "Trouble: sort() did not sort !!" << std::endl;
      return;
    }
#endif
  }


/**
   * \brief Sort the first array, and apply the resulting permutation to the second array.
   *
   * Sort the values in the first array (of length size)
   * in ascending order.  Apply the
   * permutation resulting from the sort to the second array
   *
   * @param view1 A host-accessible 1D Kokkos::View
   *   of the first array.
   * @param size Length of the first array.
   * @param view2 A host-accessible 1D Kokkos::View
   *   of the second array.  The second array must have no fewer
   *   elements than the first array.  If the first array has N
   *   elements, then the permutation will only be applied to the
   *   first N elements of the second array.
   */
  template<class View1, class View2>
  void sort2(View1 &view1, const size_t &size, View2 &view2) {
    // NOTE: This assumes the view is host-accessible.

    // Wrap the views as rcps (this happens to preserve the reference counting, but that doesn't really matter here)
    Teuchos::ArrayRCP<typename View1::non_const_value_type> view1_rcp =  Kokkos::Compat::persistingView(view1, 0, size);
    Teuchos::ArrayRCP<typename View2::non_const_value_type> view2_rcp =  Kokkos::Compat::persistingView(view2, 0, size);

    sort2(view1_rcp.begin(),view1_rcp.end(),view2_rcp.begin());    
  }

/**
   * \brief Convenience wrapper for std::sort for host-accessible views
   *
   * Sort the values in the array (of length size) in ascending order.

   * @param view A host-accessible 1D Kokkos::View.
   * @param size Length of the first array (or portion of which to sort).
   */
  template<class View>
  void sort(View &view, const size_t &size) {
    // NOTE: This assumes the view is host-accessible.

    // Wrap the view as rcps (this happens to preserve the reference counting, but that doesn't really matter here)
    Teuchos::ArrayRCP<typename View::non_const_value_type> view_rcp =  Kokkos::Compat::persistingView(view, 0, size);

    std::sort(view_rcp.begin(),view_rcp.end());    
  }

  /**
   * \brief Convenience wrapper for a reversed std::sort for host-accessible views
   *
   * Reverse Sort the values in the array (of length size) in ascending order.

   * @param view A host-accessible 1D Kokkos::View.
   * @param size Length of the array (or portion of which to sort, from the *end*)
   */
  template<class View>
  void reverse_sort(View &view, const size_t &size) {
    // NOTE: This assumes the view is host-accessible.
    // Wrap the view as rcps (this happens to preserve the reference counting, but that doesn't really matter here)
    Teuchos::ArrayRCP<typename View::non_const_value_type> view_rcp =  Kokkos::Compat::persistingView(view, 0, size);

    std::sort(view_rcp.rbegin(),view_rcp.rend());    
  }
  



  /**
   * \brief Sort the first array, and apply the same permutation to the second
   * and third arrays.
   *
   * Sort the values in the first array (represented by the exclusive
   * iterator range first1,last1) in ascending order.  Apply the
   * permutation resulting from the sort to the second array
   * (represented by a starting iterator first2) and third array
   * (represented by a starting iterator first3).
   *
   * @param first1 A random access iterator pointing to the beginning of the first array.
   * @param last1 A random access iterator pointing to the end (exclusive) of the first array.
   * @param first2 A random access iterator pointing to the beginning of the second array.
   * @param first3 A random access iterator pointing to the beginning of the third array.
   */
  template<class IT1, class IT2, class IT3>
  void sort3(const IT1 &first1, const IT1 &last1, const IT2 &first2,
    const IT3 &first3, const bool stableSort=false)
  {
    // Quicksort uses best-case N log N time whether or not the input
    // data is sorted.  However, the common case in Tpetra is that the
    // input data are sorted, so we first check whether this is the
    // case before reverting to Quicksort.
    if(SortDetails::isAlreadySorted(first1, last1)){
      return;
    }
    if(stableSort)
      SortDetails::std_sort3(first1, last1, first2, first2+(last1-first1), first3,
                      first3+(last1-first1));
    else
      SortDetails::sh_sort3(first1, last1, first2, first2+(last1-first1), first3,
                      first3+(last1-first1));
#ifdef HAVE_TPETRA_DEBUG
    if(!SortDetails::isAlreadySorted(first1, last1)){
        std::cout << " Trouble sort did not actually sort... !!!!!!" <<
                        std::endl;
      return;
    }
#endif
  }

  /// \brief Merge values in place, additively, with the same index.
  ///
  /// \tparam IT1 Iterator type for the range of indices
  /// \tparam IT2 Iterator type for the range of values
  ///
  /// indBeg, indEnd defines a half-exclusive (does not include the
  /// end) range of indices, and valBeg, valEnd its corresponding
  /// range of values.  The range of values must have the same number
  /// of entries as the range of indices.  In every nondecreasing
  /// subsequence of indices, this method will merge values that have
  /// the same index, by adding the values together.  When done, it
  /// assigns the new end (exclusive) of the index range to
  /// indResultOut, and the new end (exclusive) of the value range to
  /// valResultOut.  (It is legal for the index range not to be
  /// sorted, but then only nondecreasing subsequences will get
  /// merged.)
  ///
  /// For example, if the indices on input are {0, 1, 1, 3, -1, -1,
  /// -1, 0}, and their corresponding values on input are {42.0, -4.0,
  /// -3.0, 1.5, 1.0, 2.0, 3.0}, then on exit from this function, the
  /// indices are {0, 1, 3, -1, 0}, and the values are {42.0, -7.0,
  /// 1.5, 6.0, 100.0}.
  ///
  /// On entry to the function, indResultOut may alias indEnd, and
  /// valResultOut may alias valEnd.  For example, the following code
  /// is legal:
  /// \code
  /// std::vector<int> ind (...);
  /// std::vector<double> val (...);
  /// // ... fill ind and val ...
  /// std::vector<int>::iterator indEnd = ind.end ();
  /// std::vector<int>::iterator valEnd = val.end ();
  /// merge2 (indEnd, valEnd, ind.begin (), indEnd, val.begin (), valEnd);
  /// \endcode
  /// However, the following code is <i>not</i> legal, because the
  /// return value of <tt>std::vector::end()</tt> cannot be modified:
  /// \code
  /// std::vector<int> ind (...);
  /// std::vector<double> val (...);
  /// // ... fill ind and val ...
  /// merge2 (ind.end (), val.end (), ind.begin (), ind.end (),
  ///         val.begin (), val.end ());
  /// \endcode
  template<class IT1, class IT2>
  void
  merge2 (IT1& indResultOut, IT2& valResultOut,
          IT1 indBeg, IT1 indEnd,
          IT2 valBeg, IT2 /* valEnd */)
  {
    if (indBeg == indEnd) {
      indResultOut = indBeg; // It's allowed for indResultOut to alias indEnd
      valResultOut = valBeg; // It's allowed for valResultOut to alias valEnd
    }
    else {
      IT1 indResult = indBeg;
      IT2 valResult = valBeg;
      if (indBeg != indEnd) {
        ++indBeg;
        ++valBeg;
        while (indBeg != indEnd) {
          if (*indResult == *indBeg) { // adjacent column indices equal
            *valResult += *valBeg; // merge entries by adding their values together
          } else { // adjacent column indices not equal
            *(++indResult) = *indBeg; // shift over the index
            *(++valResult) = *valBeg; // shift over the value
          }
          ++indBeg;
          ++valBeg;
        }
        ++indResult; // exclusive end of merged result
        ++valResult; // exclusive end of merged result

        // mfh 24 Feb 2014: Setting these is technically correct, but
        // since the resulting variables aren't used after this point,
        // it may result in a build error.
        //
        // indEnd = indResult;
        // valEnd = valResult;
      }
      indResultOut = indResult;
      valResultOut = valResult;
    }
  }

  /// \brief Merge values in place with the same index, using any
  ///   associative binary function.
  ///
  /// \tparam IT1 Iterator type for the range of indices
  /// \tparam IT2 Iterator type for the range of values
  /// \tparam BinaryFunction The type of a function that takes two
  ///   values and returns another value.
  ///
  /// indBeg, indEnd defines a half-exclusive (does not include the
  /// end) range of indices, and valBeg, valEnd its corresponding
  /// range of values.  The range of values must have the same number
  /// of entries as the range of indices.  In every nondecreasing
  /// subsequence of indices, this method will merge values that have
  /// the same index, by using the given binary function.  When done,
  /// it assigns the new end (exclusive) of the index range to
  /// indResultOut, and the new end (exclusive) of the value range to
  /// valResultOut.  (It is legal for the index range not to be
  /// sorted, but then only nondecreasing subsequences will get
  /// merged.)
  ///
  /// For example, if the indices on input are {0, 1, 1, 3, -1, -1,
  /// -1, 0}, their corresponding values on input are {42.0, -4.0,
  /// -3.0, 1.5, 1.0, 2.0, 3.0}, and the binary function is an
  /// instance of <tt>std::plus<double></tt>, then on exit from this
  /// function, the indices are {0, 1, 3, -1, 0}, and the values are
  /// {42.0, -7.0, 1.5, 6.0, 100.0}.
  ///
  /// On entry to the function, indResultOut may alias indEnd, and
  /// valResultOut may alias valEnd.  For example, the following code
  /// is legal:
  /// \code
  /// std::vector<int> ind (...);
  /// std::vector<double> val (...);
  /// // ... fill ind and val ...
  /// std::vector<int>::iterator indEnd = ind.end ();
  /// std::vector<int>::iterator valEnd = val.end ();
  /// merge2 (indEnd, valEnd, ind.begin (), indEnd,
  ///         val.begin (), valEnd, std::plus<double> ());
  /// \endcode
  /// However, the following code is <i>not</i> legal, because the
  /// return value of <tt>std::vector::end()</tt> cannot be modified:
  /// \code
  /// std::vector<int> ind (...);
  /// std::vector<double> val (...);
  /// // ... fill ind and val ...
  /// merge2 (ind.end (), val.end (), ind.begin (), ind.end (),
  ///         val.begin (), val.end (), std::plus<double> ());
  /// \endcode
  template<class IT1, class IT2, class BinaryFunction>
  void
  merge2 (IT1& indResultOut, IT2& valResultOut,
          IT1 indBeg, IT1 indEnd,
          IT2 valBeg, IT2 valEnd,
          BinaryFunction f)
  {
    if (indBeg == indEnd) {
      indResultOut = indBeg; // It's allowed for indResultOut to alias indEnd
      valResultOut = valBeg; // It's allowed for valResultOut to alias valEnd
    }
    else {
      IT1 indResult = indBeg;
      IT2 valResult = valBeg;
      if (indBeg != indEnd) {
        ++indBeg;
        ++valBeg;
        while (indBeg != indEnd) {
          if (*indResult == *indBeg) { // adjacent column indices equal
            *valResult = f (*valResult, *valBeg); // merge entries via values
          } else { // adjacent column indices not equal
            *(++indResult) = *indBeg; // shift over the index
            *(++valResult) = *valBeg; // shift over the value
          }
          ++indBeg;
          ++valBeg;
        }
        ++indResult; // exclusive end of merged result
        ++valResult; // exclusive end of merged result
        indEnd = indResult;
        valEnd = valResult;
      }
      indResultOut = indResult;
      valResultOut = valResult;
    }
  }


  /// \brief Merge two sorted (by keys) sequences of unique
  ///   (key,value) pairs by combining pairs with equal keys.
  ///
  /// \param keyBeg1 [in] Start of first sequence of keys.
  /// \param keyEnd1 [in] End (exclusive) of first sequence of keys.
  /// \param valBeg1 [in] Start of first sequence of values.
  /// \param valEnd1 [in] End (exclusive) of first sequence of values.
  ///
  /// \param keyBeg2 [in] Start of second sequence of keys.
  /// \param keyEnd2 [in] End (exclusive) of second sequence of keys.
  /// \param valBeg2 [in] Start of second sequence of values.
  /// \param valEnd2 [in] End (exclusive) of second sequence of values.
  ///
  /// \param keyOut [in/out] Output sequence of keys.
  /// \param valOut [in/out] Output sequence of values.
  ///
  /// \param f [in] Binary associative function to use to combine
  ///   values whose keys are equal.  For example, for simple
  ///   replacement, use a function like std::project1st (in the SGI
  ///   extensions to the STL) with both template parameters equal to
  ///   the value type.  For addition, use std::plus with template
  ///   parameter equal to the value type.
  ///
  /// \return Number of (key,value) pairs in the merged sequence.
  ///
  /// \warning For now, this function requires that the two input
  ///   sequences be made unique by key.  Later, we plan to relax that
  ///   requirement.
  template<class KeyInputIterType, class ValueInputIterType,
           class KeyOutputIterType, class ValueOutputIterType,
           class BinaryFunction>
  void
  keyValueMerge (KeyInputIterType keyBeg1, KeyInputIterType keyEnd1,
                 ValueInputIterType valBeg1, ValueInputIterType valEnd1,
                 KeyInputIterType keyBeg2, KeyInputIterType keyEnd2,
                 ValueInputIterType valBeg2, ValueInputIterType valEnd2,
                 KeyOutputIterType keyOut, ValueOutputIterType valOut,
                 BinaryFunction f)
  {
    while (keyBeg1 != keyEnd1 && keyBeg2 != keyEnd2) {
      if (*keyBeg1 < *keyBeg2) {
        *keyOut++ = *keyBeg1++;
        *valOut++ = *valBeg1++;
      } else if (*keyBeg1 > *keyBeg2) {
        *keyOut++ = *keyBeg2++;
        *valOut++ = *valBeg2++;
      } else { // *keyBeg1 == *keyBeg2
        *keyOut++ = *keyBeg1;
        *valOut++ = f (*valBeg1, *valBeg2);
        ++keyBeg1;
        ++valBeg1;
        ++keyBeg2;
        ++valBeg2;
      }
    }
    // At most one of the two sequences will be nonempty.
    std::copy (keyBeg1, keyEnd1, keyOut);
    std::copy (valBeg1, valEnd1, valOut);
    std::copy (keyBeg2, keyEnd2, keyOut);
    std::copy (valBeg2, valEnd2, valOut);
  }

  template<class KeyInputIterType>
  size_t
  keyMergeCount (KeyInputIterType keyBeg1, KeyInputIterType keyEnd1,
                 KeyInputIterType keyBeg2, KeyInputIterType keyEnd2)
  {
    size_t count = 0;
    while (keyBeg1 != keyEnd1 && keyBeg2 != keyEnd2) {
      if (*keyBeg1 < *keyBeg2) {
        ++keyBeg1;
      } else if (*keyBeg1 > *keyBeg2) {
        ++keyBeg2;
      } else { // *keyBeg1 == *keyBeg2
        ++keyBeg1;
        ++keyBeg2;
      }
      ++count;
    }
    // At most one of the two sequences will be nonempty.
    count += static_cast<size_t> (keyEnd1 - keyBeg1) +
      static_cast<size_t> (keyEnd2 - keyBeg2);
    return count;
  }

  namespace Details {
    /// \brief Whether the two communicators are congruent.
    ///
    /// Two communicators are <i>congruent</i> when they have the same
    /// number of processes, and those processes occur in the same
    /// rank order.
    ///
    /// If both communicators are Teuchos::MpiComm instances, this
    /// function returns <tt>true</tt> exactly when
    /// <tt>MPI_Comm_compare</tt> returns <tt>MPI_IDENT</tt> (the
    /// communicators are handles for the same object) or
    /// <tt>MPI_CONGRUENT</tt> on their MPI_Comm handles.  Any two
    /// Teuchos::SerialComm instances are always congruent.  An
    /// MpiComm instance is congruent to a SerialComm instance if and
    /// only if the MpiComm has one process.  This function is
    /// symmetric in its arguments.
    ///
    /// If either Teuchos::Comm instance is neither an MpiComm nor a
    /// SerialComm, this method cannot do any better than to compare
    /// their process counts.
    bool
    congruent (const Teuchos::Comm<int>& comm1,
               const Teuchos::Comm<int>& comm2);

    /// \brief Get a Teuchos::ArrayView which views the host
    ///   Kokkos::View of the input 1-D Kokkos::DualView.
    ///
    /// \pre The input DualView must be sync'd to host.
    ///
    /// \param x [in] A specialization of Kokkos::DualView.
    ///
    /// \return Teuchos::ArrayView that views the host version of the
    ///   DualView's data.
    template<class DualViewType>
    Teuchos::ArrayView<typename DualViewType::t_dev::value_type>
    getArrayViewFromDualView (const DualViewType& x)
    {
      static_assert (static_cast<int> (DualViewType::t_dev::rank) == 1,
                     "The input DualView must have rank 1.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (x.need_sync_host (), std::logic_error, "The "
         "input Kokkos::DualView was most recently modified on device, but this "
         "function needs the host view of the data to be the most recently "
         "modified.");

      auto x_host = x.view_host ();
      typedef typename DualViewType::t_dev::value_type value_type;
      // mfh 15 Jan 2019: In debug mode, Teuchos::ArrayView's
      // constructor throws if the pointer is nonnull but the length
      // is nonpositive.
      const auto len = x_host.extent (0);
      return Teuchos::ArrayView<value_type> (len != 0 ? x_host.data () : nullptr,
                                             len);
    }

    /// \brief Get a 1-D Kokkos::DualView which is a deep copy of the
    ///   input Teuchos::ArrayView (which views host memory).
    ///
    /// \tparam T The type of the entries of the input Teuchos::ArrayView.
    /// \tparam DT The Kokkos Device type.
    ///
    /// \param x_av [in] The Teuchos::ArrayView to copy.
    /// \param label [in] String label for the Kokkos::DualView.
    /// \param leaveOnHost [in] If true, the host version of the
    ///   returned Kokkos::DualView is most recently updated (and the
    ///   DualView may need a sync to device).  If false, the device
    ///   version is most recently updated (and the DualView may need
    ///   a sync to host).
    ///
    /// \return Kokkos::DualView that is a deep copy of the input
    ///   Teuchos::ArrayView.
    template<class T, class DT>
    Kokkos::DualView<T*, DT>
    getDualViewCopyFromArrayView (const Teuchos::ArrayView<const T>& x_av,
                                  const char label[],
                                  const bool leaveOnHost)
    {
      using Kokkos::MemoryUnmanaged;
      typedef typename DT::memory_space DMS;
      typedef typename DT::execution_space execution_space;
      typedef Kokkos::HostSpace HMS;

      const size_t len = static_cast<size_t> (x_av.size ());
      Kokkos::View<const T*, HMS, MemoryUnmanaged> x_in (x_av.getRawPtr (), len);
      Kokkos::DualView<T*, DT> x_out (label, len);
      if (leaveOnHost) {
        x_out.modify_host ();
        // DEEP_COPY REVIEW - NOT TESTED FOR CUDA BUILD
        Kokkos::deep_copy (x_out.view_host (), x_in);
      }
      else {
        x_out.template modify<DMS> ();
        // DEEP_COPY REVIEW - HOST-TO-DEVICE
        Kokkos::deep_copy (execution_space(), x_out.template view<DMS> (), x_in);
      }
      return x_out;
    }

    /// \brief Return the status of the given Kokkos::DualView, as a
    ///   human-readable string.
    ///
    /// This is meant for Tpetra developers as a debugging aid.
    ///
    /// \param dv [in] Kokkos::DualView
    /// \param name [in] Human-readable name of the Kokkos::DualView
    template<class DualViewType>
    std::string dualViewStatusToString (const DualViewType& dv, const char name[])
    {
      const auto host = dv.need_sync_device();
      const auto dev = dv.need_sync_host();

      std::ostringstream os;
      os << name << ": {size: " << dv.extent (0)
         << ", sync: {host: " << host << ", dev: " << dev << "}";
      return os.str ();
    }

    /// \brief Print min(x.size(), maxNumToPrint) entries of x.
    ///
    /// \return void, because returning std::ostream& won't work
    ///   if \c out is an std::ostringstream.
    template<class ArrayType>
    void
    verbosePrintArray(std::ostream& out,
                      const ArrayType& x,
                      const char name[],
                      const size_t maxNumToPrint)
    {
      out << name << ": [";

      const size_t numEnt(x.size());
      if (maxNumToPrint == 0) {
        if (numEnt != 0) {
          out << "...";
        }
      }
      else {
        const size_t numToPrint = numEnt > maxNumToPrint ?
          maxNumToPrint : numEnt;
        size_t k = 0;
        for ( ; k < numToPrint; ++k) {
          out << x[k];
          if (k + size_t(1) < numToPrint) {
            out << ", ";
          }
        }
        if (k < numEnt) {
          out << ", ...";
        }
      }
      out << "]";
    }

    /// \brief Create string prefix for each line of verbose output.
    ///
    /// \return "Proc ${myRank}: ${prefix}: " (using Python notation).
    std::unique_ptr<std::string>
    createPrefix(const int myRank,
                 const char prefix[]);

    /// \brief Create string prefix for each line of verbose output,
    ///   for a Tpetra function (not a class or instance method).
    ///
    /// \param comm [in] May be null; if not, the communicator from
    ///   which to draw the (MPI) process rank.
    ///
    /// \param functionName [in] Name of the function.
    std::unique_ptr<std::string>
    createPrefix(const Teuchos::Comm<int>* comm,
                 const char functionName[]);

    /// \brief Create string prefix for each line of verbose output,
    ///   for a method of a Tpetra class.
    ///
    /// \param className [in] Name of the class.
    ///
    /// \param methodName [in] Name of the (class or instance) method.
    std::unique_ptr<std::string>
    createPrefix(const Teuchos::Comm<int>*,
                 const char className[],
                 const char methodName[]);

  } // namespace Details
} // namespace Tpetra

#endif // TPETRA_UTIL_HPP
