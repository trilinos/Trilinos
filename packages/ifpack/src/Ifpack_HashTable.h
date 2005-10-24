//@HEADER
/*
************************************************************************

              IFPACK: Robust Algebraic Preconditioning Package
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef IFPACK_HASHTABLE_H
#define IFPACK_HASHTABLE_H

#include "Ifpack_Utils.h"

// ============================================================================ 
// Contains a single entry of the hash table
struct Node
{
  int Key;
  double Value;
  Node * Ptr;
};

// ============================================================================ 
// Contains an array of nodes with a linked list
struct NodeArray
{
  int pos;
  int size;
  Node* Array;
  NodeArray* Next;
};

// ============================================================================
// Copied from ~/Trilinos/packages/epetra/src/Epetra_HashTable.h;
// then modified to make it faster for IFPACK's needs. The idea is the 
// following: each (Key, Value) is stored in a Node; to avoid too many calls
// to new/delete, method GetNewNode() returns a pointer from the array
// NodeArray, of size ChunkSize.
//
// Usage:
//
// 1) before using any Ifpack_HashTable object, allocate memory:
//
//    Ifpack_HashTable::Init(ChunkSize);
//
// 2) allocate one or more objects:
//
//    Ifpack_HashTable Hash(size);
//
// 3) use it, then delete it:
//
//    Hash.Add(...)
//    Hash.Get(...)
//    Hash.Arrayify(...)
//
// 4) clean memory:
//
//    Ifpack_HashTable::Finalize()
//
//
// \note This class is "delicate", there is no copy constructor, several
// checks are skipped. It is NOT intended to be a general purpose hash table,
// for this please refer to
// Trilinos/packages/teuchos/src/Teuchos_HashTable.hpp. 
//
// \author Marzio Sala, ETHZ/COLAB
//
// \date 21-Oct-05
// ============================================================================ 
class Ifpack_HashTable
{
 public:

  Ifpack_HashTable (const int size, const unsigned int seed = (2654435761U) )
  {
    Size_ = size;
    Seed_ = seed;
    start_ = Size_;
    end_ = 0;

    if (size<=0)
    {
      cout << "Bad table size (" << size << ")" << endl;
      throw(-1);
    }

    Container_ = new Node * [size];

    for (int i = 0; i < size; ++i) Container_[i] = 0;

    NumValues_ = 0;
    current_ = Container_[0];
    Data_ = FirstData_;
  }

  static void Init(const int size)
  {
    FirstData_ = new NodeArray;
    FirstData_->size = size;
    FirstData_->Array = new Node[FirstData_->size];
    FirstData_->pos = 0;
    FirstData_->Next = 0;
  }

  static void Finalize()
  {
    NodeArray* Next;
    do 
    {
      Next = FirstData_->Next;
      delete[] FirstData_->Array;
      delete FirstData_;
      FirstData_ = Next;
    }
    while (Next != 0);
  }

  ~Ifpack_HashTable()
  {
    delete[] Container_;
  }

  void Replace(const int key, const double value)
  {
    Node * n = Container_[ Func(key) ];
    while( n && (n->Key != key) ) n = n->Ptr;
    if (n) n->Value = value;
    else    
      Add(key, value);
  }

  inline Node* GetNewNode()
  {
    int& pos = Data_->pos;
    Node* Next; // node to be returned back ("allocated")
    if (pos != Data_->size - 1)
    {
      // If I have space in this NodeArray, use it
      Next = &(Data_->Array[pos]);
      ++pos;
    }
    else
    {
      // Need to go to the next array
      if (Data_->Next != 0)
      {
        // array has already been allocated, use it
        Data_ = Data_->Next;
        Data_->pos = 0;
      }
      else
      {
        // need to allocate one new array
        NodeArray* NextArray = new NodeArray;
        NextArray->size = Data_->size;
        NextArray->Array = new Node[Data_->size];
        NextArray->pos = 0;
        NextArray->Next = 0;
        Data_->Next = NextArray;
        Data_ = NextArray;
      }
      Next = &(Data_->Array[Data_->pos]);
      ++(Data_->pos);
    }
    return (Next);
  }

  void Add(const int key, const double value)
  {
    int v = Func(key);
    Node * n1 = Container_[v];
    Container_[v] = GetNewNode(); 
    Container_[v]->Key = key;
    Container_[v]->Value = value;
    Container_[v]->Ptr = n1;
    ++NumValues_;

    if (v < start_) start_ = v;
    if (v > end_)   end_ = v;
  }

  inline int Size()
  {
    return(NumValues_);
  }

  double Get(const int key)
  {
    Node * n = Container_[ Func(key) ];
    while( n && (n->Key != key) ) n = n->Ptr;
    if( n ) return n->Value;
    else    return 0.0;
  }

  // \warning: we don't check for bounds, the arrays are supposed to be long
  // enough (i.e., NumValues()) to store all entries.
  void Arrayify(int* keys, double* values)
  {
    int count = 0;
    for (int i = start_ ; i <= end_ ; ++i)
    {
      current_ = Container_[i];
      while (current_ != 0)
      {
        values[count] = current_->Value;
        keys[count] = current_->Key;
        ++count;
        current_ = current_->Ptr;
      }
    }
  }

 private:
  Ifpack_HashTable& operator=(const Ifpack_HashTable& src)
    {
      //not currently supported
      abort();
      return(*this);
    }

  int NumValues_;
  Node* current_;
  int start_; // first non-empty container
  int end_; // last non-empty container
  static NodeArray* Data_; // used to store data
  static NodeArray* FirstData_; // pointer to the first NodeArray.
  Node** Container_;

  int Size_;
  unsigned int Seed_;
  inline int Func( const int key ) { return (Seed_ ^ key)%Size_; }
};
#endif
