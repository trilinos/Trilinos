/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
//@HEADER
*/

#ifndef Epetra_HashTable_H_
#define Epetra_HashTable_H_

#include "Epetra_Object.h"

// TODO the code needs to be changed for long long keys, only
// the types have been changed right now.

class Epetra_HashTable : public Epetra_Object
{
  struct Node
  {
     int Key;
     int Value;
     Node * Ptr;

     Node( const int key = 0, const int value = 0, Node * ptr = 0 )
     : Key(key), Value(value), Ptr(ptr) {}

    private:
     Node(const Node& src)
       : Key(src.Key), Value(src.Value), Ptr(src.Ptr) {}

    Node& operator=(const Node& src)
    { Key = src.Key; Value = src.Value; Ptr = src.Ptr; return(*this); }
  };

  Node ** Container_;
  int Size_;
  unsigned int Seed_;

  int Func( const int key ) { return (Seed_ ^ key)%Size_; }
     
 public:

  Epetra_HashTable( const int size, const unsigned int seed = (2654435761U) )
  : Container_(NULL),
    Size_(size),
    Seed_(seed)
  {
    if (size<=0)
      throw ReportError( "Bad Hash Table Size: " + toString(size), -1 );

    Container_ = new Node * [size];
    for( int i = 0; i < size; ++i ) Container_[i] = 0;
  }

  Epetra_HashTable( const Epetra_HashTable & obj )
  : Container_(NULL),
    Size_(obj.Size_),
    Seed_(obj.Seed_)
  {
    Container_ = new Node * [Size_];
    for( int i = 0; i < Size_; ++i ) Container_[i] = 0;
    for( int i = 0; i < Size_; ++i )
    {
      Node * ptr = obj.Container_[i];
      while( ptr ) { Add( ptr->Key, ptr->Value ); ptr = ptr->Ptr; }
    }
  }

  ~Epetra_HashTable()
  {
    Node * ptr1;
    Node * ptr2;
    for( int i = 0; i < Size_; ++i )
    {
      ptr1 = Container_[i];
      while( ptr1 ) { ptr2 = ptr1; ptr1 = ptr1->Ptr; delete ptr2; }
    }

    delete [] Container_;
  }

  void Add( const int key, const int value )
  {
    int v = Func(key);
    Node * n1 = Container_[v];
    Container_[v] = new Node(key,value,n1);
  }

  int Get( const int key )
  {
    Node * n = Container_[ Func(key) ];
    while( n && (n->Key != key) ) n = n->Ptr;
    if( n ) return n->Value;
    else    return -1;
  }

 private:
  Epetra_HashTable& operator=(const Epetra_HashTable& src)
    {
      (void)src;
      //not currently supported
      bool throw_error = true;
      if (throw_error) {
	throw ReportError("Epetra_HashTable::operator= not supported.",-1);
      }
      return(*this);
    }

};

#endif
