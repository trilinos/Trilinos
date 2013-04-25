/*
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
*/

#ifndef TPETRA_DETAILS_FIXEDHASHTABLE_DEF_HPP
#define TPETRA_DETAILS_FIXEDHASHTABLE_DEF_HPP

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_Details_FixedHashTable_decl.hpp"
#endif

#include <Teuchos_as.hpp>
#include "MurmurHash3.hpp"

namespace Mine {
  template<class T>
  std::string toString (const T& x);

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  template<>
  std::string
  toString<Teuchos::ArrayView<const std::pair<long long, int> > > (const Teuchos::ArrayView<const std::pair<long long, int> >& x) {
    typedef Teuchos::ArrayView<int>::size_type size_type;
    std::ostringstream os;
    os << "{";
    for (size_type k = 0; k < x.size (); ++k) {
      const long long first = x[k].first;
      const int second = x[k].second;
      os << "(" << first << "," << second << ")";
      if (k + 1 < x.size ()) {
	os << ", ";
      }
    }
    os << "}";
    return os.str ();
  }

  template<>
  std::string
  toString<Teuchos::ArrayView<std::pair<long long, int> > > (const Teuchos::ArrayView<std::pair<long long, int> >& x) {
    return Mine::toString<Teuchos::ArrayView<const std::pair<long long, int> > > (x.getConst ());
  }
#endif // HAVE_TEUCHOS_LONG_LONG_INT

  template<>
  std::string
  toString<Teuchos::ArrayView<const std::pair<unsigned int, int> > > (const Teuchos::ArrayView<const std::pair<unsigned int, int> >& x) {
    typedef Teuchos::ArrayView<int>::size_type size_type;
    std::ostringstream os;
    os << "{";
    for (size_type k = 0; k < x.size (); ++k) {
      const unsigned int first = x[k].first;
      const int second = x[k].second;
      os << "(" << first << "," << second << ")";
      if (k + 1 < x.size ()) {
	os << ", ";
      }
    }
    os << "}";
    return os.str ();
  }

  template<>
  std::string
  toString<Teuchos::ArrayView<std::pair<unsigned int, int> > > (const Teuchos::ArrayView<std::pair<unsigned int, int> >& x) {
    return Mine::toString<Teuchos::ArrayView<const std::pair<unsigned int, int> > > (x.getConst ());
  }

  template<>
  std::string
  toString<Teuchos::ArrayView<const std::pair<long, int> > > (const Teuchos::ArrayView<const std::pair<long, int> >& x) {
    typedef Teuchos::ArrayView<int>::size_type size_type;
    std::ostringstream os;
    os << "{";
    for (size_type k = 0; k < x.size (); ++k) {
      const long first = x[k].first;
      const int second = x[k].second;
      os << "(" << first << "," << second << ")";
      if (k + 1 < x.size ()) {
	os << ", ";
      }
    }
    os << "}";
    return os.str ();
  }

  template<>
  std::string
  toString<Teuchos::ArrayView<std::pair<long, int> > > (const Teuchos::ArrayView<std::pair<long, int> >& x) {
    return Mine::toString<Teuchos::ArrayView<const std::pair<long, int> > > (x.getConst ());
  }

  template<>
  std::string
  toString<Teuchos::ArrayView<const std::pair<int, int> > > (const Teuchos::ArrayView<const std::pair<int, int> >& x) {
    typedef Teuchos::ArrayView<int>::size_type size_type;
    std::ostringstream os;
    os << "{";
    for (size_type k = 0; k < x.size (); ++k) {
      const int first = x[k].first;
      const int second = x[k].second;
      os << "(" << first << "," << second << ")";
      if (k + 1 < x.size ()) {
	os << ", ";
      }
    }
    os << "}";
    return os.str ();
  }

  template<>
  std::string
  toString<Teuchos::ArrayView<std::pair<int, int> > > (const Teuchos::ArrayView<std::pair<int, int> >& x) {
    return Mine::toString<Teuchos::ArrayView<const std::pair<int, int> > > (x.getConst ());
  }
} // namespace Mine


namespace Tpetra {
namespace Details {

template<typename KeyType, typename ValueType>
int FixedHashTable<KeyType, ValueType>::hashFunc (const KeyType key) const {
#ifdef TPETRA_USE_MURMUR_HASH
  uint32_t k;
  MurmurHash3_x86_32 ((void *) &key, sizeof (KeyType), 1, (void *) &k);
  return (int) (k % size_);
#else
  // We are using Epetra's hash function by default, as we have
  // observed that it is much faster than the Murmur hash
  // function. However, this is not a good hash function for general
  // sets of keys.  For our typical use case, this is good.  Use
  // Murmur hash if the maps are sparse.
  const unsigned int seed = (2654435761U);

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    size_ == 0, std::logic_error, "Tpetra::Details::FixedHashTable::hashFunc: "
    "size_ == 0.  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  const int intkey = (int) ((key & 0x000000007fffffffLL) +
			    ((key & 0x7fffffff80000000LL) >> 31));
  return (int) ((seed ^ intkey) % size_);
#endif
}

template<typename KeyType, typename ValueType>
int
FixedHashTable<KeyType, ValueType>::getRecommendedSize (const int size)
{
  // A large list of prime numbers.
  // Based on a recommendation by Andres Valloud in hash forums.
  //  There are only enough primes here so that between any number N and 2*N,
  //  there will be at least about 8 to choose from (except the first 10).
  //  This is a balance between a small list of primes, and getting a
  //  collection size that doesn't waste too much space.  In addition,
  //  the primes in this table were chosen so that they do not divide
  //  256^k +- a, for 1<=k<=8, and -32<=a<=32.  This is so that
  //  using them as a modulo will not have a tendency to just throw away
  //  the most significant bits of the object's hash.  The primes (except the
  //  first ten) or not close to any power of two to avoid aliasing
  //  between hash functions based on bit manipulation and the moduli.
  int primes [ ] = {
    3, 7, 13, 23, 53, 97, 193, 389, 769, 1543,
    2237, 2423, 2617, 2797, 2999, 3167, 3359, 3539,
    3727, 3911, 4441 , 4787 , 5119 , 5471 , 5801 , 6143 , 6521 , 6827
    , 7177 , 7517 , 7853 , 8887 , 9587 , 10243 , 10937 , 11617 , 12289
    , 12967 , 13649 , 14341 , 15013 , 15727
    , 17749 , 19121 , 20479 , 21859 , 23209 , 24593 , 25939 , 27329
    , 28669 , 30047 , 31469 , 35507 , 38231 , 40961 , 43711 , 46439
    , 49157 , 51893 , 54617 , 57347 , 60077 , 62801 , 70583 , 75619
    , 80669 , 85703 , 90749 , 95783 , 100823 , 105871 , 110909 , 115963
    , 120997 , 126031 , 141157 , 151237 , 161323 , 171401 , 181499 , 191579
    , 201653 , 211741 , 221813 , 231893 , 241979 , 252079
    , 282311 , 302483 , 322649 , 342803 , 362969 , 383143 , 403301 , 423457
    , 443629 , 463787 , 483953 , 504121 , 564617 , 604949 , 645313 , 685609
    , 725939 , 766273 , 806609 , 846931 , 887261 , 927587 , 967919 , 1008239
    , 1123477 , 1198397 , 1273289 , 1348177 , 1423067 , 1497983 , 1572869
    , 1647761 , 1722667 , 1797581 , 1872461 , 1947359 , 2022253
    , 2246953 , 2396759 , 2546543 , 2696363 , 2846161 , 2995973 , 3145739
    , 3295541 , 3445357 , 3595117 , 3744941 , 3894707 , 4044503
    , 4493921 , 4793501 , 5093089 , 5392679 , 5692279 , 5991883 , 6291469
    , 6591059 , 6890641 , 7190243 , 7489829 , 7789447 , 8089033
    , 8987807 , 9586981 , 10186177 , 10785371 , 11384539 , 11983729
    , 12582917 , 13182109 , 13781291 , 14380469 , 14979667 , 15578861
    , 16178053 , 17895707 , 19014187 , 20132683 , 21251141 , 22369661
    , 23488103 , 24606583 , 25725083 , 26843549 , 27962027 , 29080529
    , 30198989 , 31317469 , 32435981 , 35791397 , 38028379 , 40265327
    , 42502283 , 44739259 , 46976221 , 49213237 , 51450131 , 53687099
    , 55924061 , 58161041 , 60397993 , 62634959 , 64871921
    , 71582857 , 76056727 , 80530643 , 85004567 , 89478503 , 93952427
    , 98426347 , 102900263 , 107374217 , 111848111 , 116322053 , 120795971
    , 125269877 , 129743807 , 143165587 , 152113427 , 161061283 , 170009141
    , 178956983 , 187904819 , 196852693 , 205800547 , 214748383 , 223696237
    , 232644089 , 241591943 , 250539763 , 259487603 , 268435399 };

  int hsize = primes[220] ;
  for (int i = 0; i < 221; ++i) {
    if (size <= primes[i]) {
      hsize = primes[i];
      break;
    }
  }
  return hsize;
}

template<typename KeyType, typename ValueType>
FixedHashTable<KeyType, ValueType>::
FixedHashTable (const ArrayView<const KeyType>& keys) :
  size_ (0),
  rawPtr_ (NULL),
  rawVal_ (NULL)
{
  init (keys, Teuchos::as<ValueType> (0));
}

template<typename KeyType, typename ValueType>
FixedHashTable<KeyType, ValueType>::
FixedHashTable (const ArrayView<const KeyType>& keys,
		const ValueType startingValue) :
  size_ (0),
  rawPtr_ (NULL),
  rawVal_ (NULL)
{
  init (keys, startingValue);
}

template<typename KeyType, typename ValueType>
void
FixedHashTable<KeyType, ValueType>::
init (const ArrayView<const KeyType>& keys,
      const ValueType startingValue)
{
  using Teuchos::arcp;
  using Teuchos::arcp_const_cast;
  using Teuchos::ArrayRCP;
  using Teuchos::as;

  const size_type numKeys = keys.size ();
  const size_type size = getRecommendedSize (as<int> (numKeys));
#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    size == 0 && numKeys != 0, std::logic_error, 
    "Tpetra::Details::FixedHashTable constructor: "
    "getRecommendedSize(" << numKeys << ") returned zero, "
    "even though the number of keys " << numKeys << " is nonzero.  "
    "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  
  // We have to set the size_ internal state before calling the hash
  // function, since the hash function uses it.
  size_ = as<KeyType> (size);

  ArrayRCP<size_type> ptr (size + 1, 0);
  // The constructor that takes just a size argument automatically
  // fills the data.  We don't need to waste time filling it here
  // because we will do so below.  The try/catch block isn't strictly
  // necessary; we could just give "new std::pair<...> [numKeys]" to
  // the ArrayRCP constructor, since no other arguments of the
  // constructor might throw before entering the constructor.
  ArrayRCP<std::pair<KeyType, ValueType> > val;
  std::pair<KeyType, ValueType>* rawVal = NULL;
  try {
    rawVal = new std::pair<KeyType, ValueType> [numKeys];
    val = arcp<std::pair<KeyType, ValueType> > (rawVal, 0, numKeys, true);
  } catch (...) {
    if (rawVal != NULL) {
      delete [] rawVal;
    }
    throw;
  }

  // Compute number of entries in each hash table position.
  for (size_type k = 0; k < numKeys; ++k) {
    const int hashVal = hashFunc (keys[k]);
    // Shift over one, so that counts[j] = ptr[j+1].  See below.
    ++ptr[hashVal+1];
  }

  // Compute row offsets via prefix sum:
  //
  // ptr[i+1] = \sum_{j=0}^{i} counts[j].
  //
  // Thus, ptr[i+1] - ptr[i] = counts[i], so that ptr[i+1] = ptr[i] +
  // counts[i].  If we stored counts[i] in ptr[i+1] on input, then the
  // formula is ptr[i+1] += ptr[i].
  for (size_type i = 0; i < size; ++i) {
    ptr[i+1] += ptr[i];
  }
  //ptr[0] = 0; // We've already done this when initializing ptr above.

  // curRowStart[i] is the offset of the next element in row i.
  ArrayRCP<size_type> curRowStart (size, 0);

  // Fill in the hash table.
  for (size_type k = 0; k < numKeys; ++k) {
    const KeyType key = keys[k];
    const ValueType theVal = startingValue + as<ValueType> (k);
    const int hashVal = hashFunc (key);

    const size_type offset = curRowStart[hashVal];
    const size_type curPos = ptr[hashVal] + offset;

    val[curPos].first = key;
    val[curPos].second = theVal;
    ++curRowStart[hashVal];
  }

  // "Commit" the computed arrays.
  ptr_ = arcp_const_cast<const size_type> (ptr);

  // FIXME (mfh 25 Apr 2013) arcp_const_cast on val_ and
  // val_.getRawPtr() both cause a hang with MPI for some reason.  Not
  // sure what's going on.  Anyway, calling val.release(), recreating
  // val_ by hand, and using the released raw pointer as rawVal_,
  // seems to fix the problem.

  //val_ = arcp_const_cast<const std::pair<KeyType, ValueType> > (val);
  const std::pair<KeyType, ValueType>* valRaw = val.release ();
  val_ = arcp<const std::pair<KeyType, ValueType> > (valRaw, 0, numKeys, true);
  rawPtr_ = ptr_.getRawPtr ();
  //  rawVal_ = val_.getRawPtr ();
  rawVal_ = valRaw;
}

template<typename KeyType, typename ValueType>
FixedHashTable<KeyType, ValueType>::
FixedHashTable (const FixedHashTable & obj) :
  size_ (obj.size_),
  ptr_ (obj.ptr_),
  val_ (obj.val_),
  rawPtr_ (obj.rawPtr_),
  rawVal_ (obj.rawVal_)
{}

template<typename KeyType, typename ValueType>
ValueType
FixedHashTable<KeyType, ValueType>::
get (const KeyType key) const
{
  const int hashVal = hashFunc (key);
#ifdef HAVE_TPETRA_DEBUG

  const size_type start = ptr_[hashVal];
  const size_type end = ptr_[hashVal+1];
  for (size_type k = start; k < end; ++k) {
    if (val_[k].first == key) {
      return val_[k].second;
    }
  }
#else
  const size_type start = rawPtr_[hashVal];
  const size_type end = rawPtr_[hashVal+1];
  for (size_type k = start; k < end; ++k) {
    if (rawVal_[k].first == key) {
      return rawVal_[k].second;
    }
  }
#endif // HAVE_TPETRA_DEBUG
  return Teuchos::OrdinalTraits<ValueType>::invalid ();
}

template <typename KeyType, typename ValueType>
std::string FixedHashTable<KeyType, ValueType>::description() const {
  std::ostringstream oss;
  oss << "FixedHashTable<"
      << Teuchos::TypeNameTraits<KeyType>::name () << ","
      << Teuchos::TypeNameTraits<ValueType>::name () << ">: "
      << "{ numKeys: " << val_.size ()
      << ", tableSize: " << ptr_.size () << " }";
  return oss.str();
}

template <typename KeyType, typename ValueType>
void FixedHashTable<KeyType, ValueType>::describe(
 Teuchos::FancyOStream &out,
 const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  using std::setw;
  using Teuchos::OSTab;
  using Teuchos::rcpFromRef;
  using Teuchos::TypeNameTraits;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_EXTREME;

  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;

  if (vl == VERB_NONE) {
    // do nothing
  }
  else if (vl == VERB_LOW) {
    out << this->description() << endl;
  }
  else {  // MEDIUM, HIGH or EXTREME
    out << "FixedHashTable: {" << endl;
    {
      OSTab tab1 (rcpFromRef (out));

      const std::string label = this->getObjectLabel ();
      if (label != "") {
        out << "label: " << label << endl;
      }
      out << "Template parameters: {" << endl;
      {
        OSTab tab2 (rcpFromRef (out));
        out << "KeyType: " << TypeNameTraits<KeyType>::name () << endl
            << "ValueType" << TypeNameTraits<ValueType>::name () << endl;
      }

      const size_type tableSize = ptr_.size ();
      const size_type numKeys = val_.size ();

      out << "}" << endl << "Table parameters: {" << endl;
      {
        OSTab tab2 (rcpFromRef (out));
        out << "numKeys: " << numKeys << endl
            << "tableSize: " << tableSize << endl;
      }
      out << "}" << endl;

      if (vl >= VERB_EXTREME) {
        out << "Contents: ";
        if (tableSize == 0 || numKeys == 0) {
          out << "[]" << endl;
        } else {
          out << "[ " << endl;
          {
            OSTab tab2 (rcpFromRef (out));
            for (size_type i = 0; i < tableSize; ++i) {
              OSTab tab3 (rcpFromRef (out));
              for (size_type k = ptr_[i]; k < ptr_[i+1]; ++k) {
                out << "(" << val_[k].first << "," << val_[k].second << ")";
                if (k + 1 < ptr_[i+1]) {
                  out << ", ";
                }
              }
              out << endl;
            } // for each table position i
          }
          out << "]" << endl;
        } // The table contains entries
      } // vl >= VERB_EXTREME
    }
    out << "}" << endl;
  } // if vl > VERB_LOW
}

} // namespace Details
} // namespace Tpetra

// Macro that explicitly instantiates FixedHashTable for the given local
// ordinal (LO) and global ordinal (GO) types.  Note that FixedHashTable's
// template parameters occur in the opposite order of most Tpetra
// classes.  This is because FixedHashTable performs global-to-local
// lookup, and the convention in templated C++ lookup tables (such as
// std::map) is <KeyType, ValueType>.
//
// This macro must be explanded within the Tpetra::Details namespace.
#define TPETRA_DETAILS_FIXEDHASHTABLE_INSTANT_DEFAULTNODE(LO,GO) \
  template class FixedHashTable< GO , LO >;                      \

#endif // TPETRA_DETAILS_FIXEDHASHTABLE_DEF_HPP
