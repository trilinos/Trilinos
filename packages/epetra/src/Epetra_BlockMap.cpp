
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Epetra_BlockMap.h"
#include "Epetra_Comm.h"
#include "Epetra_Directory.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_HashTable.h"

// Use the new LID hash table approach by default
#define EPETRA_BLOCKMAP_NEW_LID

//==============================================================================
// Epetra_BlockMap constructor for a Epetra-defined uniform linear distribution of constant size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_Object("Epetra::BlockMap"),
    BlockMapData_(new Epetra_BlockMapData(NumGlobalElements, ElementSize, IndexBase, Comm))
{
  BlockMapData_->ConstantElementSize_ = true;
  BlockMapData_->LinearMap_ = true;
  
  // Each processor gets roughly numGlobalPoints/p points
  // This routine automatically defines a linear partitioning of a
  // map with numGlobalPoints across the processors
  // specified in the given Epetra_Comm
  
  if (BlockMapData_->NumGlobalElements_ < 0) 
    throw ReportError("NumGlobalElements = " + toString(BlockMapData_->NumGlobalElements_) + ".  Should be >= 0.", -1);
  if (BlockMapData_->ElementSize_ <= 0) 
    throw ReportError("ElementSize = " + toString(BlockMapData_->ElementSize_) + ".  Should be > 0.", -2);
  
  int NumProc = Comm.NumProc();
  int MyPID = Comm.MyPID();
  BlockMapData_->NumMyElements_ = BlockMapData_->NumGlobalElements_ / NumProc;
  int remainder = BlockMapData_->NumGlobalElements_ % NumProc;
  int start_index = MyPID * (BlockMapData_->NumMyElements_ + 1);

  if (MyPID < remainder) 
    BlockMapData_->NumMyElements_++;
  else 
    start_index -= (MyPID - remainder);

  BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumGlobalElements_ * BlockMapData_->ElementSize_;
  BlockMapData_->NumMyPoints_ = BlockMapData_->NumMyElements_ * BlockMapData_->ElementSize_;

  BlockMapData_->MinMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MinElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxElementSize_ = BlockMapData_->ElementSize_;

  BlockMapData_->MinAllGID_ = BlockMapData_->IndexBase_;
  BlockMapData_->MaxAllGID_ = BlockMapData_->MinAllGID_ + BlockMapData_->NumGlobalElements_ - 1;
  BlockMapData_->MinMyGID_ = start_index + BlockMapData_->IndexBase_;
  BlockMapData_->MaxMyGID_ = BlockMapData_->MinMyGID_ + BlockMapData_->NumMyElements_ - 1;
  BlockMapData_->DistributedGlobal_ = IsDistributedGlobal(BlockMapData_->NumGlobalElements_, BlockMapData_->NumMyElements_);

  EndOfConstructorOps();
}

//==============================================================================
// Epetra_BlockMap constructor for a user-defined linear distribution of constant size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int NumMyElements, 
				 int ElementSize, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_Object("Epetra::BlockMap"),
    BlockMapData_(new Epetra_BlockMapData(NumGlobalElements, ElementSize, IndexBase, Comm))
{
  BlockMapData_->NumMyElements_ = NumMyElements;
  BlockMapData_->MinMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MinElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->ConstantElementSize_ = true;
  BlockMapData_->LinearMap_ = true;

  // Each processor gets NumMyElements points
  
  if (BlockMapData_->NumGlobalElements_ < -1) 
    throw ReportError("NumGlobalElements = " + toString(BlockMapData_->NumGlobalElements_) + ".  Should be >= -1.", -1);
  if (BlockMapData_->NumMyElements_ < 0) 
    throw ReportError("NumMyElements = " + toString(BlockMapData_->NumMyElements_) + ".  Should be >= 0.", -2);
  if (BlockMapData_->ElementSize_ <= 0) 
    throw ReportError("ElementSize = " + toString(BlockMapData_->ElementSize_) + ". Should be > 0.", -3);

  // Get processor information

  int NumProc = Comm.NumProc();

  BlockMapData_->DistributedGlobal_ = IsDistributedGlobal(NumGlobalElements, NumMyElements);

  // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!BlockMapData_->DistributedGlobal_ || NumProc==1) {
    BlockMapData_->NumGlobalElements_ = BlockMapData_->NumMyElements_;
    CheckValidNGE(NumGlobalElements);
    
    BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumGlobalElements_ * BlockMapData_->ElementSize_;
    BlockMapData_->NumMyPoints_ = BlockMapData_->NumMyElements_ * BlockMapData_->ElementSize_;
    
    BlockMapData_->	MinAllGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MaxAllGID_ = BlockMapData_->MinAllGID_ + BlockMapData_->NumGlobalElements_ - 1;
    BlockMapData_->MinMyGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MaxMyGID_ = BlockMapData_->MinMyGID_ + BlockMapData_->NumMyElements_ - 1;
  }
  else if (NumProc > 1) {
    // Sum up all local element counts to get global count
    BlockMapData_->Comm_->SumAll(&BlockMapData_->NumMyElements_, &BlockMapData_->NumGlobalElements_, 1);
    
    CheckValidNGE(NumGlobalElements);
    
    BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumGlobalElements_ * BlockMapData_->ElementSize_;
    BlockMapData_->NumMyPoints_ = BlockMapData_->NumMyElements_ * BlockMapData_->ElementSize_;
    
    BlockMapData_->MinAllGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MaxAllGID_ = BlockMapData_->MinAllGID_ + BlockMapData_->NumGlobalElements_ - 1;
    
    // Use the ScanSum function to compute a prefix sum of the number of points
    BlockMapData_->Comm_->ScanSum(&BlockMapData_->NumMyElements_, &BlockMapData_->MaxMyGID_, 1);
    
    int start_index = BlockMapData_->MaxMyGID_ - BlockMapData_->NumMyElements_;
    BlockMapData_->MinMyGID_ = start_index + BlockMapData_->IndexBase_;
    BlockMapData_->MaxMyGID_ = BlockMapData_->MinMyGID_ + BlockMapData_->NumMyElements_ - 1;
  }
  else
    throw ReportError("Internal Error.  Report to Epetra developer", -99);
  
  EndOfConstructorOps();
}

//==============================================================================
// Epetra_BlockMap constructor for a user-defined arbitrary distribution of constant size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int * MyGlobalElements, 
				 int ElementSize, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_Object("Epetra::BlockMap"),
    BlockMapData_(new Epetra_BlockMapData(NumGlobalElements, ElementSize, IndexBase, Comm))
{
  BlockMapData_->NumMyElements_ = NumMyElements;
  BlockMapData_->MinMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxMyElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MinElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->MaxElementSize_ = BlockMapData_->ElementSize_;
  BlockMapData_->ConstantElementSize_ = true;
  BlockMapData_->LinearMap_ = false;

  int i;
  // Each processor gets NumMyElements points

  if (BlockMapData_->NumGlobalElements_ < -1) 
    throw ReportError("NumGlobalElements = " + toString(BlockMapData_->NumGlobalElements_) + ".  Should be >= -1.", -1);
  if (BlockMapData_->NumMyElements_ < 0) 
    throw ReportError("NumMyElements = " + toString(BlockMapData_->NumMyElements_) + ".  Should be >= 0.", -2);
  if (BlockMapData_->ElementSize_ <= 0) 
    throw ReportError("ElementSize = " + toString(BlockMapData_->ElementSize_) + ". Should be > 0.", -3);

  // Allocate storage for global index list information

  if (NumMyElements > 0) {
    int errorcode = BlockMapData_->MyGlobalElements_.Size(NumMyElements);
    if(errorcode != 0)
      throw ReportError("Error with MyGlobalElements allocation.", -99);
  }

  // Get processor information

  int NumProc = Comm.NumProc();
  if (NumMyElements > 0) {
    // Compute min/max GID on this processor
    BlockMapData_->MinMyGID_ = MyGlobalElements[0];
    BlockMapData_->MaxMyGID_ = MyGlobalElements[0];
    for (i = 0; i < NumMyElements; i++) {
      BlockMapData_->MyGlobalElements_[i] = MyGlobalElements[i];
      BlockMapData_->MinMyGID_ = EPETRA_MIN(BlockMapData_->MinMyGID_,MyGlobalElements[i]);
      BlockMapData_->MaxMyGID_ = EPETRA_MAX(BlockMapData_->MaxMyGID_,MyGlobalElements[i]);
    }
  }
  else {
    BlockMapData_->MinMyGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MaxMyGID_ = BlockMapData_->IndexBase_;
  }
	
  BlockMapData_->DistributedGlobal_ = IsDistributedGlobal(NumGlobalElements, NumMyElements);

  // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!BlockMapData_->DistributedGlobal_ || NumProc==1) {
    BlockMapData_->NumGlobalElements_ = BlockMapData_->NumMyElements_;
    CheckValidNGE(NumGlobalElements);
    BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumGlobalElements_ * BlockMapData_->ElementSize_;
    BlockMapData_->NumMyPoints_ = BlockMapData_->NumMyElements_ * BlockMapData_->ElementSize_;
    
    BlockMapData_->MinAllGID_ = BlockMapData_->MinMyGID_;
    BlockMapData_->MaxAllGID_ = BlockMapData_->MaxMyGID_;
  }
  else if (NumProc > 1) {
    // Sum up all local element counts to get global count
    BlockMapData_->Comm_->SumAll(&BlockMapData_->NumMyElements_, &BlockMapData_->NumGlobalElements_, 1);
    CheckValidNGE(NumGlobalElements);
    
    BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumGlobalElements_ * BlockMapData_->ElementSize_;
    BlockMapData_->NumMyPoints_ = BlockMapData_->NumMyElements_ * BlockMapData_->ElementSize_;
    
    // Use the Allreduce function to find min/max GID 
    int *tmp_send = new int[2];
    int *tmp_recv = new int[2];
    tmp_send[0] = - BlockMapData_->MinMyGID_; // Negative sign lets us do one reduction
    tmp_send[1] =   BlockMapData_->MaxMyGID_;
    BlockMapData_->Comm_->MaxAll(tmp_send, tmp_recv, 2);
    BlockMapData_->MinAllGID_ = - tmp_recv[0];
    BlockMapData_->MaxAllGID_ =   tmp_recv[1];
    delete [] tmp_send;
    delete [] tmp_recv;
    if (BlockMapData_->MinAllGID_ < BlockMapData_->IndexBase_)
      throw ReportError("Minimum global element index = " + toString(BlockMapData_->MinAllGID_) + 
			" is less than index base = " + toString(BlockMapData_->IndexBase_) +".", -5);
  }
  else
    throw ReportError("Internal Error.  Report to Epetra developer", -99);
  
  EndOfConstructorOps();
}

//==============================================================================
// Epetra_BlockMap constructor for a user-defined arbitrary distribution of variable size elements.
Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int * MyGlobalElements, 
				 int *ElementSizeList, int IndexBase, const Epetra_Comm& Comm)
  : Epetra_Object("Epetra::BlockMap"),
    BlockMapData_(new Epetra_BlockMapData(NumGlobalElements, 0, IndexBase, Comm))
{
  BlockMapData_->NumMyElements_ = NumMyElements;
  BlockMapData_->ConstantElementSize_ = false;
  BlockMapData_->LinearMap_ = false;

  int i;
  // Each processor gets NumMyElements points

  if (BlockMapData_->NumGlobalElements_ < -1) 
    throw ReportError("NumGlobalElements = " + toString(BlockMapData_->NumGlobalElements_) + ".  Should be >= -1.", -1);
  if (BlockMapData_->NumMyElements_ < 0) 
    throw ReportError("NumMyElements = " + toString(BlockMapData_->NumMyElements_) + ".  Should be >= 0.", -2);
  for (i = 0; i < BlockMapData_->NumMyElements_; i++)
    if (ElementSizeList[i] <= 0) 
      throw ReportError("ElementSizeList["+toString(i)+"] = " + toString(ElementSizeList[i]) + ". Should be > 0.", -3);
  
  // Allocate storage for global index list and element size information

  if (NumMyElements > 0) {
    int errorcode = BlockMapData_->MyGlobalElements_.Size(NumMyElements);
    if(errorcode != 0)
      throw ReportError("Error with MyGlobalElements allocation.", -99);
    errorcode = BlockMapData_->ElementSizeList_.Size(NumMyElements);
    if(errorcode != 0)
      throw ReportError("Error with ElementSizeList allocation.", -99);
  }
  // Get processor information

  int NumProc = Comm.NumProc();
  
  if (NumMyElements > 0) {
    // Compute min/max GID and element size, number of points on this processor
    BlockMapData_->MinMyGID_ = MyGlobalElements[0];
    BlockMapData_->MaxMyGID_ = MyGlobalElements[0];
    BlockMapData_->MinMyElementSize_ = ElementSizeList[0];
    BlockMapData_->MaxMyElementSize_ = ElementSizeList[0];
    BlockMapData_->NumMyPoints_ = 0;
    for (i = 0; i < NumMyElements; i++) {
      BlockMapData_->MyGlobalElements_[i] = MyGlobalElements[i];
      BlockMapData_->ElementSizeList_[i] = ElementSizeList[i];
      BlockMapData_->MinMyGID_ = EPETRA_MIN(BlockMapData_->MinMyGID_,MyGlobalElements[i]);
      BlockMapData_->MaxMyGID_ = EPETRA_MAX(BlockMapData_->MaxMyGID_,MyGlobalElements[i]);
      BlockMapData_->MinMyElementSize_ = EPETRA_MIN(BlockMapData_->MinMyElementSize_,ElementSizeList[i]);
      BlockMapData_->MaxMyElementSize_ = EPETRA_MAX(BlockMapData_->MaxMyElementSize_,ElementSizeList[i]);
      BlockMapData_->NumMyPoints_ += ElementSizeList[i];
    }
  }
  else {
    BlockMapData_->MinMyGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MaxMyGID_ = BlockMapData_->IndexBase_;
    BlockMapData_->MinMyElementSize_ = 1;
    BlockMapData_->MaxMyElementSize_ = 1;
    BlockMapData_->NumMyPoints_ = 0;
  }

  BlockMapData_->DistributedGlobal_ = IsDistributedGlobal(NumGlobalElements, NumMyElements);  
  
  // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!BlockMapData_->DistributedGlobal_ || NumProc == 1) {
    BlockMapData_->NumGlobalElements_ = BlockMapData_->NumMyElements_;
    CheckValidNGE(NumGlobalElements);
    BlockMapData_->NumGlobalPoints_ = BlockMapData_->NumMyPoints_;
    
    BlockMapData_->MinAllGID_ = BlockMapData_->MinMyGID_;
    BlockMapData_->MaxAllGID_ = BlockMapData_->MaxMyGID_;
    BlockMapData_->MinElementSize_ = BlockMapData_->MinMyElementSize_;
    BlockMapData_->MaxElementSize_ = BlockMapData_->MaxMyElementSize_;
  }
  else if (NumProc > 1) {
    // Sum up all local element and point counts to get global counts
    int *tmp_send = new int[4];
    int *tmp_recv = new int[4];
    tmp_send[0] = BlockMapData_->NumMyElements_;
    tmp_send[1] = BlockMapData_->NumMyPoints_;
    BlockMapData_->Comm_->SumAll(tmp_send, tmp_recv, 2);
    BlockMapData_->NumGlobalElements_ =  tmp_recv[0];
    BlockMapData_->NumGlobalPoints_ = tmp_recv[1];
    
    CheckValidNGE(NumGlobalElements);
    
    // Use the MaxAll function to find min/max GID 
    tmp_send[0] = - BlockMapData_->MinMyGID_; // Negative signs lets us do one reduction
    tmp_send[1] =   BlockMapData_->MaxMyGID_;
    tmp_send[2] = - BlockMapData_->MinMyElementSize_;
    if (BlockMapData_->NumMyElements_ == 0) 
      tmp_send[2] = - BlockMapData_->NumGlobalPoints_; // This processor has no elements, so should not sizes.
    tmp_send[3] =   BlockMapData_->MaxMyElementSize_;
    
    BlockMapData_->Comm_->MaxAll(tmp_send, tmp_recv, 4);
    
    BlockMapData_->MinAllGID_ =      - tmp_recv[0];
    BlockMapData_->MaxAllGID_ =        tmp_recv[1];
    BlockMapData_->MinElementSize_ = - tmp_recv[2];
    BlockMapData_->MaxElementSize_ =   tmp_recv[3];
    
    delete [] tmp_send;
    delete [] tmp_recv;
    
    if (BlockMapData_->MinAllGID_ < BlockMapData_->IndexBase_)
      throw ReportError("Minimum global element index = " + toString(BlockMapData_->MinAllGID_) + 
			" is less than index base = " + toString(BlockMapData_->IndexBase_) +".", -5);
  }
  else
    throw ReportError("Internal Error.  Report to Epetra developer", -99);
  
  EndOfConstructorOps();
}

//==============================================================================
Epetra_BlockMap::Epetra_BlockMap(const Epetra_BlockMap& map)
  : Epetra_Object(map.Label()),
    BlockMapData_(map.BlockMapData_)
{
  BlockMapData_->IncrementReferenceCount();
  
  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.
}

//==============================================================================
bool Epetra_BlockMap::SameAs(const Epetra_BlockMap & Map) const {

  // Quickest test: See if both maps share an inner data class
  if (this->BlockMapData_ == Map.BlockMapData_) 
    return(true);


  // Next check other global properties that are easy global attributes
  if (BlockMapData_->MinAllGID_ != Map.MinAllGID() ||
      BlockMapData_->MaxAllGID_ != Map.MaxAllGID() ||
      BlockMapData_->NumGlobalElements_ != Map.NumGlobalElements() ||
      BlockMapData_->IndexBase_ != Map.IndexBase()) 
    return(false);
  
  // Last possible global check for constant element sizes
  if (BlockMapData_->ConstantElementSize_ && BlockMapData_->ElementSize_!=Map.ElementSize()) 
    return(false);

  // If we get this far, we need to check local properties and then check across
  // all processors to see if local properties are all true
 
  int numMyElements = BlockMapData_->NumMyElements_;

  int MySameMap = 1; // Assume not needed
  
  // First check if number of element is the same in each map
  if (numMyElements != Map.NumMyElements()) MySameMap = 0;
  
  if (MySameMap==1) // If numMyElements is the same, check to see that list of GIDs is the same
    for (int i = 0; i < numMyElements; i++)
      if (GID(i) != Map.GID(i)) MySameMap = 0;

  // If GIDs are the same, check to see element sizes are the same
  if (MySameMap==1 && !BlockMapData_->ConstantElementSize_) {
    int * sizeList1 = ElementSizeList();
    int * sizeList2 = Map.ElementSizeList();
    for (int i = 0; i < numMyElements; i++) if (sizeList1[i] != sizeList2[i]) MySameMap=0;
  }
  // Now get min of MySameMap across all processors

  int GlobalSameMap = 0;
  int err = Comm().MinAll(&MySameMap, &GlobalSameMap, 1);
  assert(err==0);

  return(GlobalSameMap==1);
}

//==============================================================================
bool Epetra_BlockMap::PointSameAs(const Epetra_BlockMap & Map) const
{
  if (this->BlockMapData_ == Map.BlockMapData_) 
    return(true);
  
  if (BlockMapData_->NumGlobalPoints_ != Map.NumGlobalPoints() ) 
    return(false);
  
  // If we get this far, we need to check local properties and then check across
  // all processors to see if local properties are all true

  int MySameMap = 1; // Assume not needed
  if (BlockMapData_->NumMyPoints_ != Map.NumMyPoints())
    MySameMap = 0;

  int GlobalSameMap = 0;
  int err = Comm().MinAll(&MySameMap, &GlobalSameMap, 1);
  assert( err == 0 );

  return(GlobalSameMap==1);
}

//==============================================================================
int Epetra_BlockMap::MyGlobalElements(int * MyGlobalElements) const
{
  // If the global element list is not create, then do so.  This can only happen when
  // a linear distribution has been specified.  Thus we can easily construct the update
  // list in this case.

  int i;
  int numMyElements = BlockMapData_->NumMyElements_;
  
  if (BlockMapData_->MyGlobalElements_.Length() == 0)
    for (i = 0; i < numMyElements; i++)
      MyGlobalElements[i] = BlockMapData_->MinMyGID_ + i;
  else
    for (i = 0; i < numMyElements; i++)
      MyGlobalElements[i] = BlockMapData_->MyGlobalElements_[i];
  return(0);
}

//==============================================================================
int * Epetra_BlockMap::MyGlobalElements() const {
  int numMyElements = BlockMapData_->NumMyElements_;  

  // If ElementSizeList not built, do so
  if(BlockMapData_->MyGlobalElements_.Length() == 0 && numMyElements > 0) {
    int errorcode = BlockMapData_->MyGlobalElements_.Size(numMyElements + 1);
    if(errorcode != 0)
      throw ReportError("Error with MyGlobalElements allocation.", -99);
    
    // Build the array
    for (int i = 0; i < numMyElements; i++)
      BlockMapData_->MyGlobalElements_[i] = BlockMapData_->MinMyGID_ + i;
  }
  return(BlockMapData_->MyGlobalElements_.Values());
}

//==============================================================================
int Epetra_BlockMap::FirstPointInElement(int LID) const
{
  if (!MyLID(LID)) 
    EPETRA_CHK_ERR(-1);
  
  int entry;

  if (ConstantElementSize())
    entry = MaxElementSize() * LID; // convert to vector entry
  else {
    int * entrylist = FirstPointInElementList(); // get entry list
    entry = entrylist[LID];
  }
  return(entry);
}

//==============================================================================
int Epetra_BlockMap::FirstPointInElementList(int * FirstPointInElementList) const
{
  // If the first element entry list is not create, then do so.  

  // Note: This array is of length NumMyElement+1

  int i;
  int numMyElements = BlockMapData_->NumMyElements_;

  if (BlockMapData_->FirstPointInElementList_.Length() == 0) {
    FirstPointInElementList[0] = 0; // First element of first entry is always zero
    
    if (BlockMapData_->ConstantElementSize_)
      for (i = 0; i < numMyElements; i++)
	FirstPointInElementList[i+1] = FirstPointInElementList[i] + BlockMapData_->ElementSize_;
    else
      for (i = 0; i < numMyElements; i++)
	FirstPointInElementList[i+1] = FirstPointInElementList[i] + BlockMapData_->ElementSizeList_[i];
  }
  else 
    for (i = 0; i <= numMyElements; i++)
      FirstPointInElementList[i] = BlockMapData_->FirstPointInElementList_[i];
  return(0);
}

//==============================================================================
int * Epetra_BlockMap::FirstPointInElementList() const {
  int numMyElements = BlockMapData_->NumMyElements_;

  // If ElementSizeList not built, do so
  if ((BlockMapData_->FirstPointInElementList_.Length() == 0) && (numMyElements > 0)) {
    BlockMapData_->FirstPointInElementList_.Size(BlockMapData_->NumMyElements_ + 1);
    BlockMapData_->FirstPointInElementList_[0] = 0; // First element of first entry is always zero
    if (BlockMapData_->ConstantElementSize_)
      for (int i = 0; i < numMyElements; i++)
	BlockMapData_->FirstPointInElementList_[i+1] = BlockMapData_->FirstPointInElementList_[i] + BlockMapData_->ElementSize_;
    else
      for (int i = 0; i < numMyElements; i++)
	BlockMapData_->FirstPointInElementList_[i+1] = BlockMapData_->FirstPointInElementList_[i] + BlockMapData_->ElementSizeList_[i];
  }
  return(BlockMapData_->FirstPointInElementList_.Values());
}

//==============================================================================
int Epetra_BlockMap::ElementSizeList(int * ElementSizeList) const
{
  // If the element size list is not create, then do so.  This can only happen when
  // a constant element size has been specified.  Thus we can easily construct the element size
  // list in this case.

  int i;
  int numMyElements = BlockMapData_->NumMyElements_;

  if (BlockMapData_->ElementSizeList_.Length() == 0)
    for (i = 0; i < numMyElements; i++)
      ElementSizeList[i] = BlockMapData_->ElementSize_;
  else
    for (i = 0; i < numMyElements; i++)
      ElementSizeList[i] = BlockMapData_->ElementSizeList_[i];
  
  return(0);
}

//==============================================================================
int * Epetra_BlockMap::ElementSizeList() const {
  int numMyElements = BlockMapData_->NumMyElements_;

  // If ElementSizeList not built, do so
  if ((BlockMapData_->ElementSizeList_.Length() == 0) && (numMyElements > 0)) {
    BlockMapData_->ElementSizeList_.Size(numMyElements);
    for (int i = 0; i < numMyElements; i++)
      BlockMapData_->ElementSizeList_[i] = BlockMapData_->ElementSize_;
  }
  return(BlockMapData_->ElementSizeList_.Values());
}

//==============================================================================
int Epetra_BlockMap::PointToElementList(int * PointToElementList) const {
  // Build an array such that the local element ID is stored for each point

  int i;
  if (BlockMapData_->PointToElementList_.Length() == 0) {
    int numMyElements = BlockMapData_->NumMyElements_;
    int * ptr = PointToElementList;
    for (i = 0; i < numMyElements; i++) {
      int Size = ElementSize(i);
      for (int j = 0; j < Size; j++) 
	*ptr++ = i;
    }
  }
  else {
    int numMyPoints = BlockMapData_->NumMyPoints_;
    for (i = 0; i < numMyPoints; i++)
      PointToElementList[i] = BlockMapData_->PointToElementList_[i];
  }
  return(0);
}

//==============================================================================
int * Epetra_BlockMap::PointToElementList() const {

  // If PointToElementList not built, do so
  if ((BlockMapData_->PointToElementList_.Length() == 0) && (BlockMapData_->NumMyPoints_ > 0)) {
    BlockMapData_->PointToElementList_.Size(BlockMapData_->NumMyPoints_);
    int numMyElements = BlockMapData_->NumMyElements_;
    int * ptr = BlockMapData_->PointToElementList_.Values();
    for (int i = 0; i < numMyElements; i++) {
      int Size = ElementSize(i);
      for (int j = 0; j < Size; j++) 
	*ptr++ = i;
    }
  }
  return(BlockMapData_->PointToElementList_.Values());
}

//==============================================================================
int Epetra_BlockMap::ElementSize(int LID) const {

  if (ConstantElementSize()) 
    return(BlockMapData_->ElementSize_);
  else
    return(BlockMapData_->ElementSizeList_[LID]);
}

//==============================================================================
void Epetra_BlockMap::GlobalToLocalSetup()
{
  int i;
  int numMyElements = BlockMapData_->NumMyElements_;

  if (BlockMapData_->NumGlobalElements_ == 0) {
    return; // Nothing to do
  }

  if (LinearMap() || numMyElements == 0) {
    return; // Nothing else to do
  }

  // Build LID_ vector to make look up of local index values fast

#ifdef EPETRA_BLOCKMAP_NEW_LID

  //check for initial contiguous block
  int val = BlockMapData_->MyGlobalElements_[0];
  for( i = 0 ; i < numMyElements; ++i ) {
    if (val != BlockMapData_->MyGlobalElements_[i]) break;
    ++val;
  }
  BlockMapData_->LastContiguousGIDLoc_ = i - 1;
  if (BlockMapData_->LastContiguousGIDLoc_ < 0) {
    BlockMapData_->LastContiguousGID_ = BlockMapData_->MyGlobalElements_[0];
  }
  else {
    BlockMapData_->LastContiguousGID_ =
      BlockMapData_->MyGlobalElements_[BlockMapData_->LastContiguousGIDLoc_];
  }

  //Hash everything else
  if(i < numMyElements) {
    if (BlockMapData_->LIDHash_ != NULL) {
      delete BlockMapData_->LIDHash_;
    }

    BlockMapData_->LIDHash_ = new Epetra_HashTable(numMyElements - i + 1 );
    for(; i < numMyElements; ++i )
      BlockMapData_->LIDHash_->Add( BlockMapData_->MyGlobalElements_[i], i );
  }
    
#else
    
  int SpanGID = BlockMapData_->MaxMyGID_ - BlockMapData_->MinMyGID_ + 1;
  BlockMapData_->LID_.Size(SpanGID);
    
  for (i = 0; i < SpanGID; i++) 
    BlockMapData_->LID_[i] = -1; // Fill all locations with -1
    
  for (i = 0; i < numMyElements; i++) {
    int tmp = BlockMapData_->MyGlobalElements_[i] - BlockMapData_->MinMyGID_;
    assert(tmp >= 0); 
    assert(tmp < SpanGID);
    BlockMapData_->LID_[BlockMapData_->MyGlobalElements_[i] - BlockMapData_->MinMyGID_] = i; // Spread local indices
  }

#endif

}

//==============================================================================
int Epetra_BlockMap::LID(int GID) const
{
  if ((GID < BlockMapData_->MinMyGID_) || (GID > BlockMapData_->MaxMyGID_)) {
    return(-1); // Out of range
  }

  if (BlockMapData_->LinearMap_) {
    return(GID - BlockMapData_->MinMyGID_); // Can compute with an offset
  }

  if( GID >= BlockMapData_->MyGlobalElements_[0] &&
      GID <= BlockMapData_->LastContiguousGID_ ) {
    return( GID - BlockMapData_->MyGlobalElements_[0] );
  }

#ifdef EPETRA_BLOCKMAP_NEW_LID
  return BlockMapData_->LIDHash_->Get( GID );
#else
  return(BlockMapData_->LID_[GID - BlockMapData_->MinMyGID_]); // Find it in LID array  
#endif
}

//==============================================================================
int Epetra_BlockMap::GID(int LID) const
{
  if ((LID < BlockMapData_->MinLID_) || (LID > BlockMapData_->MaxLID_) ||
      (BlockMapData_->NumMyElements_ == 0)) {
    return(BlockMapData_->IndexBase_ - 1); // Out of range
  }

  if (LinearMap()) {
    return(LID + BlockMapData_->MinMyGID_); // Can compute with an offset
  }

  return(BlockMapData_->MyGlobalElements_[LID]); // Find it in MyGlobalElements array
}

//==============================================================================
int Epetra_BlockMap::FindLocalElementID(int PointID, int & ElementID, int & ElementOffset) const {

  if (PointID >= BlockMapData_->NumMyPoints_)
    return(-1); // Point is out of range
  
  if (ConstantElementSize()) {
    ElementID = PointID / BlockMapData_->MaxElementSize_;
    ElementOffset = PointID % BlockMapData_->MaxElementSize_;
    return(0);
  }
  else {
    int * tmpPointToElementList = PointToElementList();
    int * tmpFirstPointInElementList = FirstPointInElementList();
    ElementID = tmpPointToElementList[PointID];
    ElementOffset = PointID - tmpFirstPointInElementList[ElementID];
    return(0);
  }
}

//==============================================================================
int Epetra_BlockMap::RemoteIDList(int NumIDs, const int * GIDList,
				  int * PIDList, int * LIDList,
				  int * SizeList) const
{
  if (BlockMapData_->Directory_ == NULL) {
    BlockMapData_->Directory_ = Comm().CreateDirectory(*this);
  }

  Epetra_Directory* directory = BlockMapData_->Directory_;
  if (directory == NULL) {
    return(-1);
  }

  EPETRA_CHK_ERR( directory->GetDirectoryEntries(*this, NumIDs, GIDList,
						 PIDList, LIDList, SizeList) );

  return(0);
}

//==============================================================================
bool Epetra_BlockMap::IsDistributedGlobal(int NumGlobalElements, int NumMyElements) const {

  bool DistributedGlobal = false; // Assume map is not global distributed
  if (BlockMapData_->Comm_->NumProc() > 1) {
    int LocalReplicated = 0;
    int AllLocalReplicated;
    if (NumGlobalElements == NumMyElements) 
      LocalReplicated=1;
    BlockMapData_->Comm_->MinAll(&LocalReplicated, &AllLocalReplicated, 1);
    
    // If any PE has LocalReplicated=0, then map is distributed global
    if (AllLocalReplicated != 1) 
      DistributedGlobal = true;
  }
  return(DistributedGlobal);
}

//==============================================================================
void Epetra_BlockMap::CheckValidNGE(int NumGlobalElements) const {
  // Check to see if user's value for NumGlobalElements is either -1 
  // (in which case we use our computed value) or matches ours.
  if ((NumGlobalElements != -1) && (NumGlobalElements != BlockMapData_->NumGlobalElements_))
    throw ReportError("Invalid NumGlobalElements.  NumGlobalElements = " + toString(NumGlobalElements) + 
		      ".  Should equal " + toString(BlockMapData_->NumGlobalElements_) + 
		      ", or be set to -1 to compute automatically", -4);
}

//==============================================================================
void Epetra_BlockMap::EndOfConstructorOps() {
  BlockMapData_->MinLID_ = 0;
  BlockMapData_->MaxLID_ = EPETRA_MAX(BlockMapData_->NumMyElements_ - 1, 0);
  
  GlobalToLocalSetup(); // Setup any information for making global index to local index translation fast.
}

//==============================================================================
void Epetra_BlockMap::Print(ostream & os) const
{
  int * MyGlobalElements1 = MyGlobalElements();
  int * FirstPointInElementList1 = 0;
  int * ElementSizeList1 = 0;
  if (!ConstantElementSize()) {
    FirstPointInElementList1 = FirstPointInElementList();
    ElementSizeList1 = ElementSizeList();
  }
  int MyPID = Comm().MyPID();
  int NumProc = Comm().NumProc();
  
  for (int iproc = 0; iproc < NumProc; iproc++) {
    if (MyPID == iproc) {
      if (MyPID == 0) {
	os <<  "\nNumber of Global Elements  = "; os << NumGlobalElements(); os << endl;
	os <<    "Number of Global Points = "; os << NumGlobalPoints(); os << endl;
	os <<    "Maximum of all GIDs        = "; os << MaxAllGID(); os << endl;
	os <<    "Minimum of all GIDs        = "; os << MinAllGID(); os << endl;
	os <<    "Index Base                 = "; os << IndexBase(); os << endl;
	if (ConstantElementSize())
	  os <<  "Constant Element Size      = "; os << ElementSize(); os << endl;
      }
      os << endl;
      
      os <<    "Number of Local Elements   = "; os << NumMyElements(); os << endl;
      os <<    "Number of Local Points  = "; os << NumMyPoints(); os << endl;
      os <<    "Maximum of my GIDs         = "; os << MaxMyGID(); os << endl;
      os <<    "Minimum of my GIDs         = "; os << MinMyGID(); os << endl;
      os << endl;
			
      os.width(14);
      os <<  "     MyPID"; os << "    ";
      os.width(14);
      os <<  "       Local Index "; os << " ";
      os.width(14);
      os <<  "      Global Index "; os << " ";
      if (!ConstantElementSize()) {
	os.width(14);
	os <<" FirstPointInElement "; os << " ";
	os.width(14);
	os <<"   ElementSize "; os << " ";
      }
      os << endl;
			
      for (int i = 0; i < NumMyElements(); i++) {
	os.width(14);
	os <<  MyPID; os << "    ";
	os.width(14);
	os <<  i; os << "    ";
	os.width(14);
	os <<  MyGlobalElements1[i]; os << "    ";
	if (!ConstantElementSize()) {	  
	  os.width(14);
	  os << FirstPointInElementList1[i]; os << "    ";
	  os.width(14);
	  os << ElementSizeList1[i]; os << "    ";
	}
	os << endl;
      }
      
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }
  return;
}

//==============================================================================
Epetra_BlockMap::~Epetra_BlockMap()  {
  CleanupData();
}

//==============================================================================
void Epetra_BlockMap::CleanupData()
{
  if(BlockMapData_ != 0) {

    BlockMapData_->DecrementReferenceCount();
    if(BlockMapData_->ReferenceCount() == 0) {
      delete BlockMapData_;
      BlockMapData_ = 0;
    }
  }
}

//=============================================================================
Epetra_BlockMap & Epetra_BlockMap::operator= (const Epetra_BlockMap & map)
{
  if((this != &map) && (BlockMapData_ != map.BlockMapData_)) {
    CleanupData();
    BlockMapData_ = map.BlockMapData_;
    BlockMapData_->IncrementReferenceCount();
  }

  return(*this);
}
