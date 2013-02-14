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

#ifndef TPETRA_HASHTABLE_DEF_HPP_
#define TPETRA_HASHTABLE_DEF_HPP_

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_HashTable_decl.hpp"
#endif


#include "MurmurHash3.hpp"


namespace Tpetra
{

namespace Details
{

template<typename KeyType, typename ValueType>
int
Tpetra_HashTable<KeyType, ValueType>:: hashFunc( const KeyType key ) {
    uint32_t k;
    MurmurHash3_x86_32((void *)&key, sizeof(KeyType),
                          1, (void *)&k);
    return (int) (k%Size_);
    // Epetra's version of hash, keeping it here for testing temporarily
    //int intkey = (int) ((key & 0x000000007fffffffLL) +
       //((key & 0x7fffffff80000000LL) >> 31));
    //return (int) ((Seed_ ^ intkey)%Size_);
}

template<typename KeyType, typename ValueType>
int
Tpetra_HashTable<KeyType, ValueType>:: getRecommendedSize( const int size ) {
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
    for (int i = 0 ; i < 221 ; i++)
    {
        if (size <= primes[i])
        {
            hsize = primes[i] ;
            break ;
        }
    }

    return hsize ;
}

template<typename KeyType, typename ValueType>
Tpetra_HashTable<KeyType, ValueType>::
  Tpetra_HashTable( const int size, const unsigned int seed )
  : Container_(NULL),
    Seed_(seed)
  {
  TEUCHOS_TEST_FOR_EXCEPTION(size < 0, std::runtime_error,
    "Tpetra_HashTable : ERROR, size cannot be less than zero");

    Size_ = getRecommendedSize(size);
    Container_ = new Node * [Size_];
    for( int i = 0; i < Size_; ++i ) Container_[i] = NULL;
#ifdef HAVE_TEUCHOS_DEBUG
    maxc_ = 0;
    nc_ = 0;
#endif
  }

template<typename KeyType, typename ValueType>
Tpetra_HashTable<KeyType, ValueType>::
  Tpetra_HashTable( const Tpetra_HashTable & obj )
  : Container_(NULL),
    Size_(obj.Size_),
    Seed_(obj.Seed_)
  {
#ifdef HAVE_TEUCHOS_DEBUG
    maxc_ = 0;
    nc_ = 0;
#endif
    Container_ = new Node * [Size_];
    for( int i = 0; i < Size_; ++i ) Container_[i] = NULL;
    for( int i = 0; i < Size_; ++i ) {
      Node * ptr = obj.Container_[i];
      while( ptr ) { add( ptr->Key, ptr->Value ); ptr = ptr->Ptr; }
    }
  }

template<typename KeyType, typename ValueType>
Tpetra_HashTable<KeyType, ValueType>::
  ~Tpetra_HashTable() {
    Node * ptr1;
    Node * ptr2;
    for( int i = 0; i < Size_; ++i ) {
      ptr1 = Container_[i];
      while( ptr1 ) { ptr2 = ptr1; ptr1 = ptr1->Ptr; delete ptr2; }
    }

    delete [] Container_;
  }

template<typename KeyType, typename ValueType>
void
Tpetra_HashTable<KeyType, ValueType>::
  add( const KeyType key, const ValueType value ) {
    int v = hashFunc(key);
    Node * n1 = Container_[v];
    Container_[v] = new Node(key,value,n1);
  }

template<typename KeyType, typename ValueType>
ValueType
Tpetra_HashTable<KeyType, ValueType>::
  get( const KeyType key ) {
    Node * n = Container_[ hashFunc(key) ];

#ifdef HAVE_TEUCHOS_DEBUG
    int k = 0;
#endif

    while( n && (n->Key != key) ){
         n = n->Ptr;
#ifdef HAVE_TEUCHOS_DEBUG
         ((k+1 > maxc_) ? maxc_ = k+1 : 0) ;
         k++;
#endif
    }

#ifdef HAVE_TEUCHOS_DEBUG
    if (k != 0) nc_++;
#endif
    if( n ) return n->Value;
    else    return -1;
  }

template <typename KeyType, typename ValueType>
std::string Tpetra_HashTable<KeyType, ValueType>::description() const {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "Tpetra_HashTable size=" << Size_ << " "; 
    return oss.str();
  }

template <typename KeyType, typename ValueType>
void Tpetra_HashTable<KeyType, ValueType>::describe(
 Teuchos::FancyOStream &out,
 const Teuchos::EVerbosityLevel verbLevel) const {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;

    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;

    if (vl == VERB_NONE) {
      // do nothing
    }
    else if (vl == VERB_LOW) {
      out << this->description() << endl;
    }
    else {  // MEDIUM, HIGH or EXTREME
      out << this->description() << endl;
#ifdef HAVE_TEUCHOS_DEBUG
    out << endl
        << "Maximum number of collisions in Tpetra_HashTable (for a key) "
        << maxc_ << std::endl;
    out << "Number of collisions in Tpetra_HashTable " << nc_
        << std::endl;
#endif
    }
  }

} // namespace Details

} // namespace Tpetra

//! Explicit instantiation macro supporting the HashTable class,
//on the default node for specified ordinals.
// Must be explanded within the Tpetra::Details namespace
#define TPETRA_HASHTABLE_INSTANT_DEFAULTNODE(LO,GO) \
  template class Tpetra_HashTable< GO , LO >;                         \

#endif
