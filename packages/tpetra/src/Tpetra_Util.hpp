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


  /** sort function for two arrays
     The values in sortVals will be sorted in ascending order.
     The same permutation required to sort sortVals will be applied
     to otherVals.
   */
  template<class IT1, class IT2>
  void sort2(const IT1 &first1, const IT1 &last1, const IT2 &first2) {
    sort2Shell(first1,last1,first2);
  }

  /** sort function for three arrays
     The values in sortVals will be sorted in ascending order.
     The same permutation required to sort sortVals will be applied
     to otherVals.
   */
  template<class IT1, class IT2, class IT3>
  void sort3(const IT1 &first1, const IT1 &last1, const IT2 &first2, const IT3 &first3)
  {
    sort3Shell(first1, last1, first2, first3);
  }


  template<class IT1, class IT2>
  void sort2Shell(const IT1 &first1, const IT1 &last1, const IT2 &first2)
  {
    typedef typename std::iterator_traits<IT1>::difference_type DT;
    const DT sz  = last1 - first1;
    const DT DT0 = OrdinalTraits<DT>::zero();
    for (DT i = DT0; i < sz; ++i){
      DT m = sz/2;
      while (m > DT0){
        DT max = sz - m;
        for (DT j = DT0; j < max; ++j){
          for (DT k = j; k >= DT0; k-=m){
            if (first1[k+m] >= first1[k]){
              break;
            }
            std::swap( first1[k+m], first1[k] );
            std::swap( first2[k+m], first2[k] );
          }
        }
        m = m/2;
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    for (DT i=DT0; i < sz-1; ++i) {
      TEST_FOR_EXCEPTION( first1[i] > first1[i+1], std::logic_error,
          "Tpetra::sort2Shell(): internal Tpetra error. Please contact Tpetra team.");
    }
#endif
  }


  template<class IT1, class IT2, class IT3>
  void sort3Shell(const IT1 &first1, const IT1 &last1, const IT2 &first2, const IT3 &first3)
  {
    typedef typename std::iterator_traits<IT1>::difference_type DT;
    const DT sz  = last1 - first1;
    const DT DT0 = OrdinalTraits<DT>::zero();
    for (DT i = DT0; i < sz; ++i){
      DT m = sz/2;
      while (m > DT0){
        DT max = sz - m;
        for (DT j = DT0; j < max; ++j){
          for (DT k = j; k >= DT0; k-=m){
            if (first1[k+m] >= first1[k]){
              break;
            }
            std::swap( first1[k+m], first1[k] );
            std::swap( first2[k+m], first2[k] );
            std::swap( first3[k+m], first3[k] );
          }
        }
        m = m/2;
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    for (DT i=DT0; i < sz-1; ++i) {
      TEST_FOR_EXCEPTION( first1[i] > first1[i+1], std::logic_error,
          "Tpetra::sort3Shell(): internal Tpetra error. Please contact Tpetra team.");
    }
#endif
  }


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
