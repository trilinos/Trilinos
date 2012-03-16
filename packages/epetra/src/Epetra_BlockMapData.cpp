
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

#include "Epetra_BlockMapData.h"
#include "Epetra_HashTable.h"
#include "Epetra_Comm.h"
#include "Epetra_Directory.h"
//#include "Epetra_ConfigDefs.h" //DATA_DEBUG
// Use the new LID hash table approach by default
#define EPETRA_BLOCKMAP_NEW_LID

//=============================================================================
Epetra_BlockMapData::Epetra_BlockMapData(long long NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm & Comm) 
  : Comm_(Comm.Clone()),
    Directory_(0),
    LID_(0),
    MyGlobalElements_int_(0),
    MyGlobalElements_LL_(0),
    FirstPointInElementList_(0),
    ElementSizeList_(0),
    PointToElementList_(0),
    NumGlobalElements_(NumGlobalElements),
    NumMyElements_(0),
    IndexBase_(IndexBase),
    ElementSize_(ElementSize),
    MinMyElementSize_(0),
    MaxMyElementSize_(0),
    MinElementSize_(0),
    MaxElementSize_(0),
    MinAllGID_(0),
    MaxAllGID_(0),
    MinMyGID_(0),
    MaxMyGID_(-1),
    MinLID_(0),
    MaxLID_(0),
    NumGlobalPoints_(0),
    NumMyPoints_(0),
    ConstantElementSize_(false),
    LinearMap_(false),
    DistributedGlobal_(false),
    OneToOneIsDetermined_(false),
    OneToOne_(false),
    GlobalIndicesInt_(false),
    GlobalIndicesLongLong_(false),
    LastContiguousGID_(0),
    LastContiguousGIDLoc_(0),
    LIDHash_(0)
{
  //cout << "--BMD created, addr: " << this << endl; //DATA_DEBUG
}

//=============================================================================
Epetra_BlockMapData::~Epetra_BlockMapData()
{
  if(LIDHash_ != 0) {
    delete LIDHash_;
    LIDHash_ = 0;
  }

  if (Directory_ != 0) {
    delete Directory_;
    Directory_ = 0;
  }

  if(Comm_ != 0) {
    delete Comm_;
    Comm_ = 0;
  }
  //cout << "--BMD destroyed, addr: " << this << endl; //DATA_DEBUG
}
