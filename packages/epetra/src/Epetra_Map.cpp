
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

#include "Epetra_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_HashTable.h"

//==============================================================================
// Epetra_Map constructor for a Epetra-defined uniform linear distribution of elements.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map::Epetra_Map(int numGlobalElements, int indexBase, const Epetra_Comm& comm)
  : Epetra_BlockMap(numGlobalElements, 1, indexBase, comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map::Epetra_Map(long long numGlobalElements, int indexBase, const Epetra_Comm& comm)
  : Epetra_BlockMap(numGlobalElements, 1, indexBase, comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}

Epetra_Map::Epetra_Map(long long numGlobalElements, long long indexBase, const Epetra_Comm& comm)
  : Epetra_BlockMap(numGlobalElements, 1, indexBase, comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
#endif
//==============================================================================
// Epetra_Map constructor for a user-defined linear distribution of constant block size elements.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map::Epetra_Map(int numGlobalElements, int numMyElements, int indexBase, const Epetra_Comm& comm)
  : Epetra_BlockMap(numGlobalElements, numMyElements, 1, indexBase, comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map::Epetra_Map(long long numGlobalElements, int numMyElements, int indexBase, const Epetra_Comm& comm)
  : Epetra_BlockMap(numGlobalElements, numMyElements, 1, indexBase, comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}

Epetra_Map::Epetra_Map(long long numGlobalElements, int numMyElements, long long indexBase, const Epetra_Comm& comm)
  : Epetra_BlockMap(numGlobalElements, numMyElements, 1, indexBase, comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
#endif
//==============================================================================
// Epetra_Map constructor for a user-defined arbitrary distribution of constant block size elements.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map::Epetra_Map(int numGlobalElements, int numMyElements,
                       const int * myGlobalElements,
                       int indexBase, const Epetra_Comm& comm)
  : Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, 1, indexBase, comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map::Epetra_Map(long long numGlobalElements, int numMyElements,
                       const long long * myGlobalElements,
                       int indexBase, const Epetra_Comm& comm)
  : Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, 1, indexBase, comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}

Epetra_Map::Epetra_Map(long long numGlobalElements, int numMyElements,
                       const long long * myGlobalElements,
                       long long indexBase, const Epetra_Comm& comm)
  : Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, 1, indexBase, comm) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
#endif

//==============================================================================
// Epetra_Map constructor for a user-defined arbitrary distribution of constant block size elements w/ user provided globals.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map::Epetra_Map(int numGlobalElements, int numMyElements,
                       const int * myGlobalElements,
                       int indexBase, const Epetra_Comm& comm,
		       bool UserIsDistributedGlobal,
		       int UserMinAllGID, int UserMaxAllGID)
  : Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, 1, indexBase, comm, UserIsDistributedGlobal, UserMinAllGID, UserMaxAllGID) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map::Epetra_Map(long long numGlobalElements, int numMyElements,
                       const long long * myGlobalElements,
                       int indexBase, const Epetra_Comm& comm,
		       bool UserIsDistributedGlobal,
		       long long UserMinAllGID, long long UserMaxAllGID)
  : Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, 1, indexBase, comm, UserIsDistributedGlobal, UserMinAllGID, UserMaxAllGID) // Map is just a special case of BlockMap
{
  SetLabel("Epetra::Map");
}
#endif

//==============================================================================
Epetra_Map::Epetra_Map(const Epetra_Map& map)
  : Epetra_BlockMap(map) // Map is just a special case of BlockMap
{
}

//==============================================================================
Epetra_Map::~Epetra_Map(void)
{
}

//=============================================================================
Epetra_Map & Epetra_Map::operator= (const Epetra_Map& map)
{
  if(this != &map) {
    Epetra_BlockMap::operator=(map); // call this->Epetra_BlockMap::operator=
  }
  return(*this);
}

//=============================================================================
Epetra_Map * Epetra_Map::RemoveEmptyProcesses() const 
{  
#ifdef HAVE_MPI
  const Epetra_MpiComm * MpiComm = dynamic_cast<const Epetra_MpiComm*>(&Comm());

  // If the Comm isn't MPI, just treat this as a copy constructor
  if(!MpiComm) return new Epetra_Map(*this);      

  MPI_Comm NewComm,MyMPIComm = MpiComm->Comm();

  // Create the new communicator.  MPI_Comm_split returns a valid
  // communicator on all processes.  On processes where color == MPI_UNDEFINED,
  // ignore the result.  Passing key == 0 tells MPI to order the
  // processes in the new communicator by their rank in the old
  // communicator.
  const int color = (NumMyElements() == 0) ? MPI_UNDEFINED : 1;

  // MPI_Comm_split must be called collectively over the original
  // communicator.  We can't just call it on processes with color
  // one, even though we will ignore its result on processes with
  // color zero.
  int rv = MPI_Comm_split(MyMPIComm,color,0,&NewComm);
  if(rv!=MPI_SUCCESS) throw ReportError("Epetra_Map::RemoveEmptyProcesses: MPI_Comm_split failed.",-1);

  if(color == MPI_UNDEFINED)
    return 0; // We're not in the new map
  else {
    Epetra_MpiComm * NewEpetraComm = new Epetra_MpiComm(NewComm);

    // Use the copy constructor for a new map, but basically because it does nothing useful
    Epetra_Map * NewMap = new Epetra_Map(*this);

    // Get rid of the old BlockMapData, now make a new one from scratch...
    NewMap->CleanupData();
    if(GlobalIndicesInt()) 
      NewMap->BlockMapData_ = new Epetra_BlockMapData(NumGlobalElements(),0,IndexBase(),*NewEpetraComm,false);
    else
      NewMap->BlockMapData_ = new Epetra_BlockMapData(NumGlobalElements64(),0,IndexBase64(),*NewEpetraComm,true);

    // Now copy all of the relevent bits of BlockMapData...
    //    NewMap->BlockMapData_->Comm_                    = NewEpetraComm;
    NewMap->BlockMapData_->LID_                     = BlockMapData_->LID_;
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    NewMap->BlockMapData_->MyGlobalElements_int_    = BlockMapData_->MyGlobalElements_int_;
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    NewMap->BlockMapData_->MyGlobalElements_LL_     = BlockMapData_->MyGlobalElements_LL_;
#endif
    NewMap->BlockMapData_->FirstPointInElementList_ = BlockMapData_->FirstPointInElementList_;
    NewMap->BlockMapData_->ElementSizeList_         = BlockMapData_->ElementSizeList_;
    NewMap->BlockMapData_->PointToElementList_      = BlockMapData_->PointToElementList_;

    NewMap->BlockMapData_->NumGlobalElements_       = BlockMapData_->NumGlobalElements_;
    NewMap->BlockMapData_->NumMyElements_           = BlockMapData_->NumMyElements_;
    NewMap->BlockMapData_->IndexBase_               = BlockMapData_->IndexBase_;
    NewMap->BlockMapData_->ElementSize_             = BlockMapData_->ElementSize_;
    NewMap->BlockMapData_->MinMyElementSize_        = BlockMapData_->MinMyElementSize_;
    NewMap->BlockMapData_->MaxMyElementSize_        = BlockMapData_->MaxMyElementSize_;
    NewMap->BlockMapData_->MinElementSize_          = BlockMapData_->MinElementSize_;
    NewMap->BlockMapData_->MaxElementSize_          = BlockMapData_->MaxElementSize_;
    NewMap->BlockMapData_->MinAllGID_               = BlockMapData_->MinAllGID_;
    NewMap->BlockMapData_->MaxAllGID_               = BlockMapData_->MaxAllGID_;
    NewMap->BlockMapData_->MinMyGID_                = BlockMapData_->MinMyGID_;
    NewMap->BlockMapData_->MaxMyGID_                = BlockMapData_->MaxMyGID_;
    NewMap->BlockMapData_->MinLID_                  = BlockMapData_->MinLID_;
    NewMap->BlockMapData_->MaxLID_                  = BlockMapData_->MaxLID_;
    NewMap->BlockMapData_->NumGlobalPoints_         = BlockMapData_->NumGlobalPoints_;
    NewMap->BlockMapData_->NumMyPoints_             = BlockMapData_->NumMyPoints_;
    NewMap->BlockMapData_->ConstantElementSize_     = BlockMapData_->ConstantElementSize_;
    NewMap->BlockMapData_->LinearMap_               = BlockMapData_->LinearMap_;
    NewMap->BlockMapData_->DistributedGlobal_       = NewEpetraComm->NumProc()==1 ? false : BlockMapData_->DistributedGlobal_;
    NewMap->BlockMapData_->OneToOneIsDetermined_    = BlockMapData_->OneToOneIsDetermined_;
    NewMap->BlockMapData_->OneToOne_                = BlockMapData_->OneToOne_;
    NewMap->BlockMapData_->GlobalIndicesInt_        = BlockMapData_->GlobalIndicesInt_;
    NewMap->BlockMapData_->GlobalIndicesLongLong_   = BlockMapData_->GlobalIndicesLongLong_;
    NewMap->BlockMapData_->LastContiguousGID_       = BlockMapData_->LastContiguousGID_;
    NewMap->BlockMapData_->LastContiguousGIDLoc_    = BlockMapData_->LastContiguousGIDLoc_;
    NewMap->BlockMapData_->LIDHash_                 = BlockMapData_->LIDHash_ ? new Epetra_HashTable<int>(*BlockMapData_->LIDHash_) : 0;

    // Delay directory construction
    NewMap->BlockMapData_->Directory_               = 0;

    // Cleanup
    delete NewEpetraComm;

    return NewMap;
  }
#else
    // MPI isn't compiled, so just treat this as a copy constructor
    return new Epetra_Map(*this);      
#endif
}

//=============================================================================
Epetra_Map* Epetra_Map::ReplaceCommWithSubset(const Epetra_Comm * Comm) const
{
  // mfh 26 Mar 2013: The lazy way to do this is simply to recreate
  // the Map by calling its ordinary public constructor, using the
  // original Map's data.  This only involves O(1) all-reduces over
  // the new communicator, which in the common case only includes a
  // small number of processes.
  Epetra_Map * NewMap=0;
  
  // Create the Map to return (unless Comm is zero, in which case we return zero).
  if(Comm) {
    // Map requires that the index base equal the global min GID.
    // Figuring out the global min GID requires a reduction over all
    // processes in the new communicator.  It could be that some (or
    // even all) of these processes contain zero entries.  (Recall
    // that this method, unlike removeEmptyProcesses(), may remove
    // an arbitrary subset of processes.)  We deal with this by
    // doing a min over the min GID on each process if the process
    // has more than zero entries, or the global max GID, if that
    // process has zero entries.  If no processes have any entries,
    // then the index base doesn't matter anyway.

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    if(GlobalIndicesInt()) {
      int MyMin, IndexBase;
      MyMin  = NumMyElements() > 0 ? MinMyGID() : MaxAllGID();
      Comm->MinAll(&MyMin,&IndexBase,1);
      NewMap = new Epetra_Map(-1,NumMyElements(),MyGlobalElements(),IndexBase,*Comm);
    }
    else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    if(GlobalIndicesLongLong()) {
      long long MyMin, IndexBase;
      MyMin = NumMyElements() > 0 ? MinMyGID64() : MaxAllGID64();
      Comm->MinAll(&MyMin,&IndexBase,1);
      NewMap = new Epetra_Map(-1,NumMyElements(),MyGlobalElements64(),IndexBase,*Comm);
    }
    else
#endif
    throw ReportError("Epetra_Map::ReplaceCommWithSubset ERROR, GlobalIndices type unknown.",-1);
  }
  return NewMap;
}

