// HAVE_@HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_UTIL_HPP
#define TPETRA_UTIL_HPP

#include "Tpetra_ConfigDefs.hpp" // for map, vector, string, and iostream 
#include <iterator>
#include <algorithm>
#include <Teuchos_Utils.hpp>
#include <Teuchos_TestForException.hpp>
#include <sstream>

#if defined(HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS) || defined(HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS)
//! Handle an efficiency warning, according to HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS and HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS
#define TPETRA_EFFICIENCY_WARNING(throw_exception_test,Exception,msg)                                 \
{                                                                                                     \
  std::ostringstream errStream;                                                                       \
  errStream << Teuchos::typeName(*this) << msg;                                                       \
  std::string err = errStream.str();                                                                  \
  if (TPETRA_PRINTS_EFFICIENCY_WARNINGS && (throw_exception_test)) {                                  \
    std::cerr << err << std::endl;                                                                    \
  }                                                                                                   \
  TEST_FOR_EXCEPTION(TPETRA_THROWS_EFFICIENCY_WARNINGS && (throw_exception_test), Exception, err);    \
}
#else
//! Handle an efficiency warning, according to HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS and HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS
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
  TEST_FOR_EXCEPTION(TPETRA_THROWS_ABUSE_WARNINGS && (throw_exception_test), Exception, err);  \
}
#else
//! Handle an abuse warning, according to HAVE_TPETRA_THROW_ABUSE_WARNINGS and HAVE_TPETRA_PRINT_ABUSE_WARNINGS
#define TPETRA_ABUSE_WARNING(throw_exception_test,Exception,msg)
#endif


/** Shared test for exception
   Just like Teuchos TEST_FOR_EXCEPTION, but with the assurance that all nodes test and throw the exception together.
 */
#define SHARED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm) \
{ \
    using Teuchos::outArg; \
    const int lcl_throw_exception = (throw_exception_test) ? Teuchos::rank(comm)+1 : 0; \
    int gbl_throw; \
    Teuchos::reduceAll(comm,Teuchos::REDUCE_MAX,lcl_throw_exception,outArg(gbl_throw)); \
    TEST_FOR_EXCEPTION(gbl_throw,Exception,  \
                       msg << " Failure on node " << gbl_throw-1 << "." << std::endl); \
}

#ifdef HAVE_TEUCHOS_DEBUG
//! If TEUCHOS_DEBUG is defined, then it calls SHARED_TEST_FOR_EXCEPTION. Otherwise, it calls TEST_FOR_EXCEPTION
#define SWITCHED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm) \
{ \
    SHARED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm); \
}
#else 
//! If TEUCHOS_DEBUG is defined, then it calls SHARED_TEST_FOR_EXCEPTION. Otherwise, it calls TEST_FOR_EXCEPTION
#define SWITCHED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm) \
{ \
    TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg); \
}
#endif

namespace Tpetra {

  /*! 
    \file Tpetra_Util.hpp
    \brief Stand-alone utility functions.

    Tpetra_Util contains utility functions that are used throughout
    Tpetra, by many classes and functions. They are placed here
    so that they can be updated and maintained in a single spot.

    Utility functions housed here:
    <ul>
    <li>An efficientAddOrUpdate for inserting data into a STL map.
    
    <li>Functions for converting Ordinals to Scalars and for converting
    Scalars to Ordinals.

    <li>A templated toString function, which is mainly used to easily
    output the contents of STL containers.

    <li>A multiple-array sort function, similar to the one found in Epetra_Util.
    </ul>
  */

  /** efficientAddOrUpdate is taken from Scott Meyers' "Effective STL", Item 24.
     if m already contains an entry with key k, use operator [].
     if it doesn't, insert is used.
   */
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


  namespace SortDetails{


  /**
   * \brief Determines whether or not a data structure is already sorted.
   *
   * @param first An iterator pointing to the beginning of the 
   * data structure.
   * @param last An iterator pointing to the end of the data structure.
   * @return True if the datasctructure is already sorted, false otherwise.
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
   * \brief Preforms the partition operation that is part of the quicksort routine
   * for two arrays
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
   * \brief Preforms the partition operation that is part of the quicksort routine
   * for three arrays
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
   * \brief Prefoms a quick sort on the two arrays using the permutations
   * that are done on the first array on the second array as well.
   *
   * @param first1 An iterator pointing to the beginning of the first array.
   * @param last1 An iterator pointing to the end of the first array.
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
   * \brief Prefoms a quick sort on the three arrays using the permutations
   * that are done on the first array on the second and third arrays as well.
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


  } //end namespace SortDetails


  /** 
   * \brief Sort function for two arrays.
   *
   * The values in sortVals will be sorted in ascending order.
   * The same permutation required to sort sortVals will be applied
   * to otherVals.
   * 
   * @param first1 An iterator pointing to the beginning of the first array.
   * @param last1 An iterator pointing to the end of the first array.
   * @param first2 An iterator pointing to the beginning of the second array.
   */
  template<class IT1, class IT2>
  void sort2(const IT1 &first1, const IT1 &last1, const IT2 &first2) {
    if(SortDetails::isAlreadySorted(first1, last1)){
      return;
    }
    SortDetails::quicksort2(first1, last1, first2, first2+(last1-first1));
  }

  
  /** 
   * \brief Sort function for three arrays
   *
   * The values in sortVals will be sorted in ascending order.
   * The same permutation required to sort sortVals will be applied
   * to otherVals.
   *
   * @param first1 An iterator pointing to the beginning of the first array.
   * @param last1 An iterator pointing to the end of the first array.
   * @param first2 An iterator pointing to the beginning of the second array.
   * @param first3 An iterator pointing to the beginning of the third array.
   */
  template<class IT1, class IT2, class IT3>
  void sort3(const IT1 &first1, const IT1 &last1, const IT2 &first2, const IT3 &first3)
  {
    if(SortDetails::isAlreadySorted(first1, last1)){
      return;
    }
    SortDetails::quicksort3(first1, last1, first2, first2+(last1-first1), first3, first3+(last1-first1));
  }

  /**
   * \brief Preforms a binary search.
   *
   * Preforms a binary search for value amount the entries between
   * the first iterator and the last iterator.
   */
  template<class IT1, class T>
  IT1 binary_serach(IT1 first, IT1 last, const T& value){
	first = std::lower_bound(first,last,value);
	if(first!=last && !(value<*first)){
		return first;
	}
	else{
		return last;
	}
  }
  
} // namespace Tpetra


#endif // TPETRA_UTIL_HPP
