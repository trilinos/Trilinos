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

#ifndef TPETRA_HASHTABLE_DECL_HPP
#define TPETRA_HASHTABLE_DECL_HPP

#include <Teuchos_Describable.hpp>
#include "Tpetra_ConfigDefs.hpp"

namespace Tpetra
{

namespace Details
{


template<typename KeyType, typename ValueType>
class Tpetra_HashTable : public Teuchos::Describable {
  struct Node {
     KeyType Key;
     ValueType Value;
     Node * Ptr;

     Node( const KeyType key = 0, const ValueType value = 0, Node * ptr = NULL )
     : Key(key), Value(value), Ptr(ptr) {}

    private:
     Node(const Node& src)
       : Key(src.Key), Value(src.Value), Ptr(src.Ptr) {}

    Node& operator=(const Node& src) {
        Key = src.Key;
        Value = src.Value;
        Ptr = src.Ptr;
        return(*this);
    }
  };

  Node ** Container_;
  KeyType Size_;
  unsigned int Seed_;
#ifdef HAVE_TEUCHOS_DEBUG
  int maxc_; // max size of the list among all entriesw/ collisons. debug only.
  int nc_;   // number of entries with collisions, should use only in debug
#endif

  // Hash Function
  int hashFunc( const KeyType key );

  int getRecommendedSize( const int size );

public:

  Tpetra_HashTable( const int size, const unsigned int seed = (2654435761U) );

  Tpetra_HashTable( const Tpetra_HashTable & obj );

  ~Tpetra_HashTable();

  void add( const KeyType key, const ValueType value );

  ValueType get( const KeyType key );

  //! Implementation of \c Teuchos::Describable
  //@{
  //! Return a simple one-line description of this object.
  std::string description() const;

  //! Print this object with the given verbosity level to the given
  //\c FancyOStream.
  void
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel=
            Teuchos::Describable::verbLevel_default) const;
  //@}

private:

  //! Assignment operator (declared but not defined; do not use).
  Tpetra_HashTable<KeyType,ValueType>&
  operator= (const Tpetra_HashTable<KeyType,ValueType> & source);

};

} // Details namespace

} // Tpetra namespace

#endif
