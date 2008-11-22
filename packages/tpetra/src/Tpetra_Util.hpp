// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
#include <Teuchos_Utils.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Tpetra {

// shared test for exception
// just like Teuchos TEST_FOR_EXCEPTION, but with the assurance 
// that all nodes test and throw the exception together
#define SHARED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm) \
{ \
    const Ordinal lcl_throw_exception = (throw_exception_test) ? Teuchos::rank(comm)+1 : 0; \
    Ordinal gbl_throw; \
    Teuchos::reduceAll(comm,Teuchos::REDUCE_MAX,lcl_throw_exception,&gbl_throw); \
    TEST_FOR_EXCEPTION(gbl_throw,Exception,  \
                       msg << " Failure on node " << gbl_throw-1 << "." << std::endl); \
}

// if TEUCHOS_DEBUG is defined, then it calls SHARED_TEST_FOR_EXCEPTION
// otherwise, it calls TEST_FOR_EXCEPTION
#ifdef TEUCHOS_DEBUG
#define SWITCHED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm) \
{ \
    SHARED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm); \
}
#else 
#define SWITCHED_TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg,comm) \
{ \
    TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg); \
}
#endif


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

  // efficientAddOrUpdate is taken from Scott Meyers' "Effective STL", Item 24.
  // if m already contains an entry with key k, use operator [].
  // if it doesn't, insert is used.
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

  // this function works much the way Teuchos::Array::toString works.
  // it allows std::vector to be used with an ostream.
  // The contents of the vector are printed in the following format:
  // "{4, 7, 18, 23, 6, 2}"
  template <typename T>
  std::string toString(const std::vector<T> & x) {
    std::ostringstream os;
    os << "{";
    typename std::vector<T>::const_iterator i = x.begin();
    if (i != x.end()) {
      os << *i;
      i++;
      for(; i != x.end(); i++) {
        os << "," << *i;
      }
    }
    os << "}";
    return(os.str());
  }

  template <typename T>
  std::string toString(const std::complex<T> & x) {
    return("(" + Teuchos::toString(x.real()) + "," + Teuchos::toString(x.imag()) + ")");
  }
 
  // sort function for multiple arrays
  // The values in sortVals will be sorted in ascending order.
  // The same permutation required to sort sortVals will be applied
  // to otherVals.

  // for two arrays
  template<typename T1, typename T2>
  void sortArrays(const Teuchos::ArrayView<T1> &sortVals, const Teuchos::ArrayView<T2> &otherVals) {
    // if sortVals and otherVals are not the same length, throw exception
    TEST_FOR_EXCEPTION(sortVals.size() != otherVals.size(), std::runtime_error,
        "Error in Tpetra_Util::sortArrays(sortVals,otherVals): sortVals and otherVals are not equally sized");
    
    // copy sortVals and otherVals into a multimap
    // (using a multimap instead of a map because values in sortVals may be duplicated)
    std::multimap<T1,T2> tempMap;
    typename Teuchos::ArrayView<T1>::iterator keyIter   = sortVals.begin();
    typename Teuchos::ArrayView<T2>::iterator valueIter = otherVals.begin();
    for(; keyIter != sortVals.end(); ++keyIter, ++valueIter) {
      tempMap.insert(std::pair<T1,T2>(*keyIter, *valueIter));
    }

    // multimap will automatically sort them, we just need to pull them out in order
    // and write them back to the original arrays
    keyIter   = sortVals.begin();
    valueIter = otherVals.begin();
    for(typename std::multimap<T1,T2>::iterator i = tempMap.begin(); 
        i != tempMap.end(); ++i, ++keyIter, ++valueIter) 
    {
      *keyIter   = i->first;
      *valueIter = i->second;
    }
  }

  // for three arrays
  template<typename T1, typename T2, typename T3>
  void sortArrays(const Teuchos::ArrayView<T1> &sortVals, 
                  const Teuchos::ArrayView<T2> &otherVals1, 
                  const Teuchos::ArrayView<T3> &otherVals2) 
  {
    // if sortVals and otherVals are not the same length, throw exception
    TEST_FOR_EXCEPTION((sortVals.size() != otherVals1.size()) || (sortVals.size() != otherVals2.size()), std::runtime_error,
        "Error in Tpetra_Util::sortArrays(sortVals,otherVals1,otherVals2): sortVals and otherVals are not equally sized");
    
    // copy sortVals and otherVals into a multimap
    // (using a multimap instead of a map because values in sortVals may be duplicated)
    typedef typename std::pair<T2, T3> ValuePair;
    std::multimap<T1, ValuePair> tempMap;
    typename Teuchos::ArrayView<T1>::iterator keyIter    = sortVals.begin();
    typename Teuchos::ArrayView<T2>::iterator valueIter1 = otherVals1.begin();
    typename Teuchos::ArrayView<T3>::iterator valueIter2 = otherVals2.begin();
    for(; keyIter != sortVals.end(); ++keyIter, ++valueIter1, ++valueIter2) {
      tempMap.insert(std::pair<T1, ValuePair>(*keyIter, ValuePair(*valueIter1, *valueIter2)));
    }

    // multimap will automatically sort them, we just need to pull them out in order
    // and write them back to the original arrays
    keyIter    = sortVals.begin();
    valueIter1 = otherVals1.begin();
    valueIter2 = otherVals2.begin();
    for(typename std::multimap<T1,ValuePair>::iterator i = tempMap.begin(); 
        i != tempMap.end(); i++, keyIter++, valueIter1++, valueIter2++) 
    {
      *keyIter    = i->first;
      *valueIter1 = i->second.first;
      *valueIter2 = i->second.second;
    }
  }

} // namespace Tpetra

#endif // TPETRA_UTIL_HPP
