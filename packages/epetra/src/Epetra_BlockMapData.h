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

#ifndef EPETRA_BLOCKMAPDATA_H
#define EPETRA_BLOCKMAPDATA_H

#include "Epetra_Data.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_LongLongSerialDenseVector.h"

class Epetra_Comm;
class Epetra_Directory;
template<typename value_type> class Epetra_HashTable;

//! Epetra_BlockMapData:  The Epetra BlockMap Data Class.
/*! The Epetra_BlockMapData class is an implementation detail of Epetra_BlockMap.
    It is reference-counted, and can be shared by multiple Epetra_BlockMap instances. 
    It derives from Epetra_Data, and inherits reference-counting from it.
*/

class Epetra_BlockMapData : public Epetra_Data {
  friend class Epetra_BlockMap;

 private:

  //! @name Constructor/Destructor Methods
  //@{ 

  //! Epetra_BlockMapData Default Constructor.
  Epetra_BlockMapData(long long NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm & Comm, bool IsLongLong);

  //! Epetra_BlockMapData Destructor.
  ~Epetra_BlockMapData();

  //@}

  const Epetra_Comm * Comm_;

  mutable Epetra_Directory* Directory_;

  Epetra_IntSerialDenseVector LID_;
  Epetra_IntSerialDenseVector MyGlobalElements_int_;
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  Epetra_LongLongSerialDenseVector MyGlobalElements_LL_;
#endif
  Epetra_IntSerialDenseVector FirstPointInElementList_;
  Epetra_IntSerialDenseVector ElementSizeList_;
  Epetra_IntSerialDenseVector PointToElementList_;
  
  long long NumGlobalElements_;
  int NumMyElements_;
  int IndexBase_;
  int ElementSize_;
  int MinMyElementSize_;
  int MaxMyElementSize_;
  int MinElementSize_;
  int MaxElementSize_;
  long long MinAllGID_;
  long long MaxAllGID_;
  long long MinMyGID_;
  long long MaxMyGID_;
  int MinLID_;
  int MaxLID_;
  long long NumGlobalPoints_;
  int NumMyPoints_;
  
  bool ConstantElementSize_;
  bool LinearMap_;
  bool DistributedGlobal_;
  mutable bool OneToOneIsDetermined_;
  mutable bool OneToOne_;
  bool GlobalIndicesInt_;
  bool GlobalIndicesLongLong_;

  long long LastContiguousGID_;
  int LastContiguousGIDLoc_;
  Epetra_HashTable<int> * LIDHash_;

  // these are intentionally declared but not defined. See Epetra Developer's Guide for details.
  Epetra_BlockMapData(const Epetra_BlockMapData & BlockMapData);
  Epetra_BlockMapData& operator=(const Epetra_BlockMapData & BlockMapData);

};
#endif /* EPETRA_BLOCKMAPDATA_H */
