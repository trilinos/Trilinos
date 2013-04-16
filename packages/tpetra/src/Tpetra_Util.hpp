// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER

#ifndef TPETRA_UTIL_HPP
#define TPETRA_UTIL_HPP

/*! 
  \file Tpetra_Util.hpp
  \brief Stand-alone utility functions and macros.

  Tpetra_Util contains utility functions macros that are used
  throughout Tpetra, by many classes and functions. They are placed
  here so that they can be updated and maintained in a single spot.

  Here are some of the utility functions found in this file:
  <ul>
  <li>An efficientAddOrUpdate for inserting data into an std::map.
    
  <li>Functions for converting Ordinals to Scalars and for converting
  Scalars to Ordinals.

  <li>A templated toString function, which is mainly used to easily
  output the contents of STL containers.

  <li>A multiple-array sort function, similar to the one found in
  Epetra_Util.  

  <li>Macros for reporting efficiency warnings and synchronizing 
  tests for exceptions over a given communicator.
  </ul>
*/

#include "Tpetra_ConfigDefs.hpp" // for map, vector, string, and iostream 
#include <iterator>
#include <algorithm>
#include <Teuchos_Utils.hpp>
#include <Teuchos_Assert.hpp>
#include <sstream>

#if defined(HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS) || defined(HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS)
/// \brief Print or throw an efficency warning.
///
/// This macro is only for use by Tpetra developers.  It is only to be
/// used in the implementation of a Tpetra class' instance method.
///
/// If HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS is defined, throw an
/// exception of type Exception, whose exception message will include
/// \c msg (along with other useful information).  
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
/// HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS nor
/// HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS are defined.
///
/// Exception: The type of exception to throw, if throw_exception_test
/// evaluates to true and TPETRA_THROWS_EFFICIENCY_WARNINGS is
/// defined.  The Exception should be a subclass of std::exception.
///
/// msg: The message to include in the warning.  The warning also
/// includes the name of the class, and other useful information.
///
#define TPETRA_EFFICIENCY_WARNING(throw_exception_test,Exception,msg)  \
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
    TEUCHOS_TEST_FOR_EXCEPTION(TPETRA_THROWS_EFFICIENCY_WARNINGS && tpetraEfficiencyWarningTest, Exception, err); \
  } \
}
#else
/// \brief Print or throw an efficency warning.
///
/// This macro is only for use by Tpetra developers.  It is only to be
/// used in the implementation of a Tpetra class' instance method.
///
/// If HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS is defined, throw an
/// exception of type Exception, whose exception message will include
/// \c msg (along with other useful information).  
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
/// HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS nor
/// HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS are defined.
///
/// Exception: The type of exception to throw, if throw_exception_test
/// evaluates to true and TPETRA_THROWS_EFFICIENCY_WARNINGS is
/// defined.  The Exception should be a subclass of std::exception.
///
/// msg: The message to include in the warning.  The warning also
/// includes the name of the class, and other useful information.
///
#define TPETRA_EFFICIENCY_WARNING(throw_exception_test,Exception,msg)
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

#ifdef HAVE_TEUCHOS_DEBUG
//! If TEUCHOS_DEBUG is defined, then it calls SHARED_TEST_FOR_EXCEPTION. Otherwise, it calls TEUCHOS_TEST_FOR_EXCEPTION
#define SWITCHED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm) \
{ \
    SHARED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm); \
}
#else 
//! If TEUCHOS_DEBUG is defined, then it calls SHARED_TEST_FOR_EXCEPTION. Otherwise, it calls TEUCHOS_TEST_FOR_EXCEPTION
#define SWITCHED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm) \
{ \
    TEUCHOS_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg); \
}
#endif

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
    DT myit =OrdinalTraits<DT>::one();
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
    DT DT1 = OrdinalTraits<DT>::one();
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
    DT DT1 = OrdinalTraits<DT>::one();
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
    const IT2& last2,
    const IT3& first3,
    const IT3& last3)
   {
        typedef typename std::iterator_traits<IT1>::difference_type DT;
        DT n = last1 - first1;
        DT m = n / 2;
        DT z = OrdinalTraits<DT>::zero();
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
    const IT2& last2)
   {
        typedef typename std::iterator_traits<IT1>::difference_type DT;
        DT n = last1 - first1;
        DT m = n / 2;
        DT z = OrdinalTraits<DT>::zero();
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
  void sort2(const IT1 &first1, const IT1 &last1, const IT2 &first2) {
    // Quicksort uses best-case N log N time whether or not the input
    // data is sorted.  However, the common case in Tpetra is that the
    // input data are sorted, so we first check whether this is the
    // case before reverting to Quicksort.
    if(SortDetails::isAlreadySorted(first1, last1)){
      return;
    }
    SortDetails::sh_sort2(first1, last1, first2, first2+(last1-first1));
#ifdef HAVE_TPETRA_DEBUG
    if(!SortDetails::isAlreadySorted(first1, last1)){
      std::cout << "Trouble: sort() did not sort !!" << std::endl;
      return;
    }
#endif
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
    const IT3 &first3)
  {
    // Quicksort uses best-case N log N time whether or not the input
    // data is sorted.  However, the common case in Tpetra is that the
    // input data are sorted, so we first check whether this is the
    // case before reverting to Quicksort.
    if(SortDetails::isAlreadySorted(first1, last1)){
      return;
    }
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

} // namespace Tpetra


#endif // TPETRA_UTIL_HPP
