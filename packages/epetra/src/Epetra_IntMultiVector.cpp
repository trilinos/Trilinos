
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
#include "Epetra_IntMultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Comm.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Distributor.h"

//=============================================================================

// Epetra_BlockMap Constructor

Epetra_IntMultiVector::Epetra_IntMultiVector(const Epetra_BlockMap& map, int numVectors, bool zeroOut)
  : Epetra_DistObject(map, "Epetra::MultiVector"),
    Epetra_CompObject(),
    Values_(0),
    Pointers_(0),
    MyLength_(map.NumMyPoints()),
    GlobalLength_(map.NumGlobalPoints64()),
    NumVectors_(numVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(map.NumMyPoints()),
    Allocated_(false)
{
  Util_.SetSeed(1);

    AllocateForCopy();

    for (int i = 0; i < NumVectors_; i++) Pointers_[i] = Values_+i*Stride_;

  if(zeroOut) PutScalar(0); // Fill all vectors with zero.
}
//==========================================================================

// Copy Constructor

Epetra_IntMultiVector::Epetra_IntMultiVector(const Epetra_IntMultiVector& Source)
  : Epetra_DistObject(Source),
    Epetra_CompObject(Source),
    Values_(0),
    Pointers_(0),
    MyLength_(Source.MyLength_),
    GlobalLength_(Source.GlobalLength_),
    NumVectors_(Source.NumVectors_),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(0),
    Allocated_(false),
    Util_(Source.Util_)
{
  AllocateForCopy();

  int ** Source_Pointers = Source.Pointers();
  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Source_Pointers[i];

  DoCopy();

}
//==========================================================================

// This constructor copies in or makes view of a standard Fortran array

Epetra_IntMultiVector::Epetra_IntMultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& map,
               int *A, int MyLDA, int numVectors)
  : Epetra_DistObject(map, "Epetra::MultiVector"),
    Epetra_CompObject(),
    Values_(0),
    Pointers_(0),
    MyLength_(map.NumMyPoints()),
    GlobalLength_(map.NumGlobalPoints64()),
    NumVectors_(numVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(map.NumMyPoints()),
    Allocated_(false)
{
  Util_.SetSeed(1);

  if (CV==Copy) AllocateForCopy();
  else AllocateForView();

  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = A + i*MyLDA;

   if (CV==Copy) DoCopy();
   else DoView();

}

//==========================================================================

// This constructor copies in or makes view of a C/C++ array of pointer

Epetra_IntMultiVector::Epetra_IntMultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& map,
                int **ArrayOfPointers, int numVectors)
  : Epetra_DistObject(map, "Epetra::MultiVector"),
    Epetra_CompObject(),
    Values_(0),
    Pointers_(0),
    MyLength_(map.NumMyPoints()),
    GlobalLength_(map.NumGlobalPoints64()),
    NumVectors_(numVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(map.NumMyPoints()),
    Allocated_(false)
{
  Util_.SetSeed(1);

  if (CV==Copy) AllocateForCopy();
  else AllocateForView();

  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = ArrayOfPointers[i];

   if (CV==Copy) DoCopy();
   else DoView();

}

//==========================================================================

// This constructor copies or makes view of selected vectors, specified in Indices,
// from an existing MultiVector

Epetra_IntMultiVector::Epetra_IntMultiVector(Epetra_DataAccess CV, const Epetra_IntMultiVector& Source,
               int *Indices, int numVectors)
  : Epetra_DistObject(Source.Map(), "Epetra::MultiVector"),
    Epetra_CompObject(),
    Values_(0),
    Pointers_(0),
    MyLength_(Source.MyLength_),
    GlobalLength_(Source.GlobalLength_),
    NumVectors_(numVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(0),
    Allocated_(false)
{
  Util_.SetSeed(1);

  if (CV==Copy) AllocateForCopy();
  else AllocateForView();

  int ** Source_Pointers = Source.Pointers();
  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Source_Pointers[Indices[i]];

   if (CV==Copy) DoCopy();
   else DoView();

}

//==========================================================================

// This interface copies or makes view of a range of vectors from an existing IntMultiVector

Epetra_IntMultiVector::Epetra_IntMultiVector(Epetra_DataAccess CV, const Epetra_IntMultiVector& Source,
               int StartIndex, int numVectors)
  : Epetra_DistObject(Source.Map(), "Epetra::MultiVector"),
    Epetra_CompObject(),
    Values_(0),
    Pointers_(0),
    MyLength_(Source.MyLength_),
    GlobalLength_(Source.GlobalLength_),
    NumVectors_(numVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(0),
    Allocated_(false)
{
  Util_.SetSeed(1);

  if (CV==Copy) AllocateForCopy();
  else AllocateForView();

  int ** Source_Pointers = Source.Pointers();
  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Source_Pointers[StartIndex+i];

   if (CV==Copy) DoCopy();
   else DoView();
}

//=========================================================================
Epetra_IntMultiVector::~Epetra_IntMultiVector(){

  if (!Allocated_) return;

  delete [] Pointers_;
  if (!UserAllocated_ && Values_!=0) delete [] Values_;

  if (IntVectors_!=0) {
    for (int i=0; i<NumVectors_; i++) if (IntVectors_[i]!=0) delete IntVectors_[i];
    delete [] IntVectors_;
  }


  if (OrdinalTemp_!=0) delete [] OrdinalTemp_;

}

//=========================================================================
int Epetra_IntMultiVector::AllocateForCopy(void)
{

  if (Allocated_) return(0);

  if (NumVectors_<=0)
    throw ReportError("Number of Vectors = " + toString(NumVectors_) + ", but must be greater than zero", -1);

  Stride_ = Map_.NumMyPoints();
  if (Stride_>0) Values_ = new int[Stride_ * NumVectors_];
  Pointers_ = new int *[NumVectors_];

  OrdinalTemp_ = 0;
  IntVectors_ = 0;

  int randval = rand(); // Use POSIX standard random function
  if (DistributedGlobal())
    Util_.SetSeed(2*Comm_->MyPID() + randval);
  else {
    int locrandval = randval;
    Comm_->MaxAll(&locrandval, &randval, 1);
    Util_.SetSeed(randval); // Replicated Local MultiVectors must have same seeds
  }

  Allocated_ = true;
  UserAllocated_ = false;
  return(0);
}

//=========================================================================
int Epetra_IntMultiVector::DoCopy(void)
{
  // On entry Pointers_ contains pointers to the incoming vectors.  These
  // pointers are the only unique piece of information for each of the
  // constructors.

  // \internal { Optimization of this function can impact performance since it
  //           involves a fair amount of memory traffic.}

  // On exit, Pointers_ is redefined to point to its own MultiVector vectors.

  for (int i = 0; i< NumVectors_; i++)
    {
      int * from = Pointers_[i];
      int * to = Values_+i*Stride_;
      Pointers_[i] = to;
      const int myLength = MyLength_;
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(to,from)
      for (int j=0; j<myLength; j++) to[j] = from[j];
#else
      memcpy(to, from, myLength*sizeof(int));
#endif
    }

  return(0);
}
//=========================================================================
int Epetra_IntMultiVector::AllocateForView(void)
{

  if (NumVectors_<=0)
    throw ReportError("Number of Vectors = " + toString(NumVectors_) + ", but must be greater than zero", -1);

  Pointers_ = new int *[NumVectors_];

  OrdinalTemp_ = 0;
  IntVectors_ = 0;

  int randval = rand(); // Use POSIX standard random function
  if (DistributedGlobal())
    Util_.SetSeed(2*Comm_->MyPID() + randval);
  else {
    int locrandval = randval;
    Comm_->MaxAll(&locrandval, &randval, 1);
    Util_.SetSeed(randval); // Replicated Local MultiVectors must have same seeds
  }

  Allocated_ = true;
  UserAllocated_ = true;

  return(0);
}

//=========================================================================
int Epetra_IntMultiVector::DoView(void)
{
  // On entry Pointers_ contains pointers to the incoming vectors.  These
  // pointers are the only unique piece of information for each of the
  // constructors.


  Values_ = Pointers_[0];

  if (NumVectors_ == 1) {
    Stride_ = Map_.NumMyPoints();
    ConstantStride_ = true;
    return(0);
  }

  // Remainder of code checks if MultiVector has regular stride

  Stride_ = (int)(Pointers_[1] - Pointers_[0]);
  ConstantStride_ = false;

  for (int i = 1; i < NumVectors_-1; i++) {
    if (Pointers_[i+1] - Pointers_[i] != Stride_) return(0);
  }

  ConstantStride_ = true;

  return(0);
}
//=========================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_IntMultiVector::ReplaceGlobalValue(int GlobalRow, int VectorIndex, int OrdinalValue) {

 // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue<int>(GlobalRow, 0, VectorIndex, OrdinalValue, false));
  return(0);
}
#endif
//=========================================================================
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_IntMultiVector::ReplaceGlobalValue(long long GlobalRow, int VectorIndex, int OrdinalValue) {

 // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue<long long>(GlobalRow, 0, VectorIndex, OrdinalValue, false));
  return(0);
}
#endif
//=========================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_IntMultiVector::ReplaceGlobalValue(int GlobalBlockRow, int BlockRowOffset,
                                              int VectorIndex, int OrdinalValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue<int>(GlobalBlockRow, BlockRowOffset, VectorIndex, OrdinalValue, false));
  return(0);
}
#endif
//=========================================================================
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_IntMultiVector::ReplaceGlobalValue(long long GlobalBlockRow, int BlockRowOffset,
                                              int VectorIndex, int OrdinalValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue<long long>(GlobalBlockRow, BlockRowOffset, VectorIndex, OrdinalValue, false));
  return(0);
}
#endif
//=========================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_IntMultiVector::SumIntoGlobalValue(int GlobalRow, int VectorIndex, int OrdinalValue) {

  // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue<int>(GlobalRow, 0, VectorIndex, OrdinalValue, true));
  return(0);
}
#endif
//=========================================================================
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_IntMultiVector::SumIntoGlobalValue(long long GlobalRow, int VectorIndex, int OrdinalValue) {

  // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue<long long>(GlobalRow, 0, VectorIndex, OrdinalValue, true));
  return(0);
}
#endif
//=========================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_IntMultiVector::SumIntoGlobalValue(int GlobalBlockRow, int BlockRowOffset,
             int VectorIndex, int OrdinalValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue<int>(GlobalBlockRow, BlockRowOffset, VectorIndex, OrdinalValue, true));
  return(0);
}
#endif
//=========================================================================
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_IntMultiVector::SumIntoGlobalValue(long long GlobalBlockRow, int BlockRowOffset,
             int VectorIndex, int OrdinalValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue<long long>(GlobalBlockRow, BlockRowOffset, VectorIndex, OrdinalValue, true));
  return(0);
}
#endif
//=========================================================================
int Epetra_IntMultiVector::ReplaceMyValue(int MyRow, int VectorIndex, int OrdinalValue) {

  // Use the more general method below
  EPETRA_CHK_ERR(ChangeMyValue(MyRow, 0, VectorIndex, OrdinalValue, false));
  return(0);
}
//=========================================================================
int Epetra_IntMultiVector::ReplaceMyValue(int MyBlockRow, int BlockRowOffset,
             int VectorIndex, int OrdinalValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeMyValue(MyBlockRow, BlockRowOffset, VectorIndex, OrdinalValue, false));
  return(0);
}
//=========================================================================
int Epetra_IntMultiVector::SumIntoMyValue(int MyRow, int VectorIndex, int OrdinalValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeMyValue(MyRow, 0, VectorIndex, OrdinalValue, true));
  return(0);
}
//=========================================================================
int Epetra_IntMultiVector::SumIntoMyValue(int MyBlockRow, int BlockRowOffset,
             int VectorIndex, int OrdinalValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeMyValue(MyBlockRow, BlockRowOffset, VectorIndex, OrdinalValue, true));
  return(0);
}
//=========================================================================
template<typename int_type>
int Epetra_IntMultiVector::ChangeGlobalValue(int_type GlobalBlockRow, int BlockRowOffset,
             int VectorIndex, int OrdinalValue, bool SumInto) {

  if(!Map().template GlobalIndicesIsType<int_type>())
  throw ReportError("Epetra_IntMultiVector::ChangeGlobalValues mismatch between argument types (int/long long) and map type.", -1);

  // Convert GID to LID and call LID version
  EPETRA_CHK_ERR(ChangeMyValue(Map().LID(GlobalBlockRow), BlockRowOffset, VectorIndex, OrdinalValue, SumInto));
  return(0);
}
//=========================================================================
int Epetra_IntMultiVector::ChangeMyValue(int MyBlockRow, int BlockRowOffset,
             int VectorIndex, int OrdinalValue, bool SumInto) {

  if (!Map().MyLID(MyBlockRow)) EPETRA_CHK_ERR(1); // I don't own this one, return a warning flag
  if (VectorIndex>= NumVectors_) EPETRA_CHK_ERR(-1); // Consider this a real error
  if (BlockRowOffset<0 || BlockRowOffset>=Map().ElementSize(MyBlockRow)) EPETRA_CHK_ERR(-2); // Offset is out-of-range

  int entry = Map().FirstPointInElement(MyBlockRow);

  if (SumInto)
    Pointers_[VectorIndex][entry+BlockRowOffset] += OrdinalValue;
  else
    Pointers_[VectorIndex][entry+BlockRowOffset] = OrdinalValue;

  return(0);
}

//=========================================================================

// Extract a copy of a Epetra_IntMultiVector.  Put in a user's Fortran-style array

int Epetra_IntMultiVector::ExtractCopy(int *A, int MyLDA) const {
  if (NumVectors_>1 && Stride_ > MyLDA) EPETRA_CHK_ERR(-1); // LDA not big enough

  const int myLength = MyLength_;
  for (int i=0; i< NumVectors_; i++)
    {
      int * from = Pointers_[i];
      int * to = A + i*MyLDA;
      for (int j=0; j<myLength; j++) *to++ = *from++;
    }

  return(0);
}

//=========================================================================
int Epetra_IntMultiVector::ReplaceMap(const Epetra_BlockMap& map)
{
  // mfh 28 Mar 2013: We can't check for compatibility across the
  // whole communicator, unless we know that the current and new
  // Maps are nonnull on _all_ participating processes.

  // So, we'll check to make sure that the maps are the same size on this processor and then
  // just go with it.
  if(Map().NumMyElements() == map.NumMyElements() && Map().NumMyPoints() == map.NumMyPoints()) {
    Epetra_DistObject::Map_ = map;
    return(0);
  }

  return(-1);
}

//=========================================================================

// Extract a copy of a Epetra_IntMultiVector.  Put in a user's array of pointers

int Epetra_IntMultiVector::ExtractCopy(int **ArrayOfPointers) const {
  const int myLength = MyLength_;
  for (int i=0; i< NumVectors_; i++)
    {
      int * from = Pointers_[i];
      int * to = ArrayOfPointers[i];
      memcpy(to, from, myLength*sizeof(int));
    }

  return(0);
}



//=========================================================================

// Extract a view of a Epetra_IntMultiVector.  Set up a user's Fortran-style array

int Epetra_IntMultiVector::ExtractView(int **A, int *MyLDA) const {
  if (!ConstantStride_) EPETRA_CHK_ERR(-1);  // Can't make a Fortran-style view if not constant stride
  *MyLDA = Stride_; // Set user's LDA
  *A = Values_; // Set user's value pointer
  return(0);
}



//=========================================================================

// Extract a view of a Epetra_IntMultiVector.  Put in a user's array of pointers

int Epetra_IntMultiVector::ExtractView(int ***ArrayOfPointers) const {
  *ArrayOfPointers = Pointers_;

  return(0);
}


//=========================================================================
int Epetra_IntMultiVector::PutScalar(int ScalarConstant) {

  // Fills MultiVector with the value ScalarConstant **/

  const int myLength = MyLength_;
  for (int i = 0; i < NumVectors_; i++) {
    int * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarConstant)
#endif
    for (int j=0; j<myLength; j++) to[j] = ScalarConstant;
  }
  return(0);
}
//=========================================================================
int Epetra_IntMultiVector::CheckSizes(const Epetra_SrcDistObject& Source) {

  const Epetra_IntMultiVector & A = dynamic_cast<const Epetra_IntMultiVector &>(Source);

  if (NumVectors()!=A.NumVectors()) {EPETRA_CHK_ERR(-3)};
  return(0);
}

//=========================================================================
int Epetra_IntMultiVector::CopyAndPermute(const Epetra_SrcDistObject& Source,
                                          int NumSameIDs,
                                          int NumPermuteIDs,
                                          int * PermuteToLIDs,
                                          int *PermuteFromLIDs,
                                          const Epetra_OffsetIndex * Indexor,
                                          Epetra_CombineMode CombineMode)
{
  (void)Indexor;

  const Epetra_IntMultiVector & A = dynamic_cast<const Epetra_IntMultiVector &>(Source);

  int **From = A.Pointers();
  int **To = Pointers_;
  int numVectors = NumVectors_;

  int * ToFirstPointInElementList = 0;
  int * FromFirstPointInElementList = 0;
  int * FromElementSizeList = 0;
  int MaxElementSize = Map().MaxElementSize();
  bool ConstantElementSize = Map().ConstantElementSize();

  if (!ConstantElementSize) {
    ToFirstPointInElementList =   Map().FirstPointInElementList();
    FromFirstPointInElementList = A.Map().FirstPointInElementList();
    FromElementSizeList = A.Map().ElementSizeList();
  }
  int jj, jjj, k;

  int NumSameEntries;

  bool Case1 = false;
  bool Case2 = false;
  // bool Case3 = false;

  if (MaxElementSize==1) {
    Case1 = true;
    NumSameEntries = NumSameIDs;
  }
  else if (ConstantElementSize) {
    Case2 = true;
    NumSameEntries = NumSameIDs * MaxElementSize;
  }
  else {
    // Case3 = true;
    NumSameEntries = FromFirstPointInElementList[NumSameIDs];
  }

  // Short circuit for the case where the source and target vector is the same.
  if (To==From) NumSameEntries = 0;

  // Do copy first
  if (NumSameIDs>0)
    if (To!=From) {
      for (int i=0; i < numVectors; i++) {
        if (CombineMode==Epetra_AddLocalAlso) for (int j=0; j<NumSameEntries; j++) To[i][j] += From[i][j]; // Add to existing value
        else for (int j=0; j<NumSameEntries; j++) To[i][j] = From[i][j];
    }
    }
  // Do local permutation next
  if (NumPermuteIDs>0) {

    // Point entry case
    if (Case1) {

      if (numVectors==1) {
        if (CombineMode==Epetra_AddLocalAlso) for (int j=0; j<NumPermuteIDs; j++) To[0][PermuteToLIDs[j]] += From[0][PermuteFromLIDs[j]]; // Add to existing value
        else for (int j=0; j<NumPermuteIDs; j++) To[0][PermuteToLIDs[j]] = From[0][PermuteFromLIDs[j]];
  }
      else {
  for (int j=0; j<NumPermuteIDs; j++) {
    jj = PermuteToLIDs[j];
    jjj = PermuteFromLIDs[j];
    if (CombineMode==Epetra_AddLocalAlso) for (int i=0; i<numVectors; i++) To[i][jj] += From[i][jjj]; // Add to existing value
    else for (int i=0; i<numVectors; i++) To[i][jj] = From[i][jjj];
  }
      }
    }
    // constant element size case
    else if (Case2) {

      for (int j=0; j<NumPermuteIDs; j++) {
        jj = MaxElementSize*PermuteToLIDs[j];
        jjj = MaxElementSize*PermuteFromLIDs[j];
        if (CombineMode==Epetra_AddLocalAlso) for (int i=0; i<numVectors; i++) for (k=0; k<MaxElementSize; k++) To[i][jj+k] += From[i][jjj+k]; // Add to existing value
        else for(int i=0; i<numVectors; i++) for (k=0; k<MaxElementSize; k++) To[i][jj+k] = From[i][jjj+k];
      }
    }

    // variable element size case
    else {

      for (int j=0; j<NumPermuteIDs; j++) {
        jj = ToFirstPointInElementList[PermuteToLIDs[j]];
        jjj = FromFirstPointInElementList[PermuteFromLIDs[j]];
        int ElementSize = FromElementSizeList[PermuteFromLIDs[j]];
        if (CombineMode==Epetra_AddLocalAlso) for (int i=0; i<numVectors; i++) for (k=0; k<ElementSize; k++) To[i][jj+k] += From[i][jjj+k]; // Add to existing value
        else for (int i=0; i<numVectors; i++) for (k=0; k<ElementSize; k++) To[i][jj+k] = From[i][jjj+k];
      }
    }
  }
  return(0);
}

//=========================================================================
int Epetra_IntMultiVector::PackAndPrepare(const Epetra_SrcDistObject & Source,
                                          int NumExportIDs,
                                          int * ExportLIDs,
                                          int & LenExports,
                                          char * & Exports,
                                          int & SizeOfPacket,
                                          int * Sizes,
                                          bool & VarSizes,
                                          Epetra_Distributor & Distor)
{
  (void)Sizes;
  (void)VarSizes;
  (void)Distor;

  const Epetra_IntMultiVector & A = dynamic_cast<const Epetra_IntMultiVector &>(Source);
  int jj, k;

  int **From = A.Pointers();
  int MaxElementSize = Map().MaxElementSize();
  int numVectors = NumVectors_;
  bool ConstantElementSize = Map().ConstantElementSize();

  int * FromFirstPointInElementList = 0;
  int * FromElementSizeList = 0;

  if (!ConstantElementSize) {
    FromFirstPointInElementList = A.Map().FirstPointInElementList();
    FromElementSizeList = A.Map().ElementSizeList();
  }

  int * OrdinalExports = 0;

  SizeOfPacket = numVectors*MaxElementSize*(int)sizeof(int);

  if(SizeOfPacket*NumExportIDs>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = SizeOfPacket*NumExportIDs;
    OrdinalExports = new int[numVectors*MaxElementSize*NumExportIDs];
    Exports = (char *) OrdinalExports;
  }

  int * ptr;

  if (NumExportIDs>0) {
    ptr = (int *) Exports;

    // Point entry case
    if (MaxElementSize==1) {

      if (numVectors==1)
        for (int j=0; j<NumExportIDs; j++)
          *ptr++ = From[0][ExportLIDs[j]];

      else {
        for (int j=0; j<NumExportIDs; j++) {
          jj = ExportLIDs[j];
          for (int i=0; i<numVectors; i++)
            *ptr++ = From[i][jj];
        }
      }
    }

    // constant element size case
    else if (ConstantElementSize) {

      for (int j=0; j<NumExportIDs; j++) {
        jj = MaxElementSize*ExportLIDs[j];
        for (int i=0; i<numVectors; i++)
          for (k=0; k<MaxElementSize; k++)
            *ptr++ = From[i][jj+k];
      }
    }

    // variable element size case
    else {

      int thisSizeOfPacket = numVectors*MaxElementSize;
      for (int j=0; j<NumExportIDs; j++) {
        ptr = (int *) Exports + j*thisSizeOfPacket;
        jj = FromFirstPointInElementList[ExportLIDs[j]];
        int ElementSize = FromElementSizeList[ExportLIDs[j]];
        for (int i=0; i<numVectors; i++)
          for (k=0; k<ElementSize; k++)
            *ptr++ = From[i][jj+k];
      }
    }
  }

  return(0);
}

//=========================================================================
int Epetra_IntMultiVector::UnpackAndCombine(const Epetra_SrcDistObject & Source,
                                            int NumImportIDs,
                                            int * ImportLIDs,
                                            int LenImports,
                                            char * Imports,
                                            int & SizeOfPacket,
                                            Epetra_Distributor & Distor,
                                            Epetra_CombineMode CombineMode,
                                            const Epetra_OffsetIndex * Indexor )
{
  (void)Source;
  (void)LenImports;
  (void)SizeOfPacket;
  (void)Distor;
  (void)Indexor;
  int jj, k;

  if(    CombineMode != Add
      && CombineMode != Zero
      && CombineMode != Insert
      && CombineMode != InsertAdd
      && CombineMode != Average
      && CombineMode != Epetra_Max
      && CombineMode != Epetra_Min
      && CombineMode != AbsMin
      && CombineMode != AbsMax )
      EPETRA_CHK_ERR(-1); //Unsupported CombinedMode, will default to Zero

  if (NumImportIDs<=0) return(0);

  int ** To = Pointers_;
  int numVectors = NumVectors_;
  int MaxElementSize = Map().MaxElementSize();
  bool ConstantElementSize = Map().ConstantElementSize();

  int * ToFirstPointInElementList = 0;
  int * ToElementSizeList = 0;

  if (!ConstantElementSize) {
    ToFirstPointInElementList = Map().FirstPointInElementList();
    ToElementSizeList = Map().ElementSizeList();
  }

  int * ptr;
  // Unpack it...

  ptr = (int *) Imports;

  // Point entry case
  if (MaxElementSize==1) {

    if (numVectors==1) {
      if (CombineMode==InsertAdd) for (int j=0; j<NumImportIDs; j++) To[0][ImportLIDs[j]] = 0.0; // Zero out first
      if (CombineMode==Add || CombineMode==InsertAdd) for (int j=0; j<NumImportIDs; j++) To[0][ImportLIDs[j]] += *ptr++; // Add to existing value
      else if(CombineMode==Insert) for (int j=0; j<NumImportIDs; j++) To[0][ImportLIDs[j]] = *ptr++;
      else if(CombineMode==AbsMax) for (int j=0; j<NumImportIDs; j++) {To[0][ImportLIDs[j]] = EPETRA_MAX( To[0][ImportLIDs[j]],std::abs(*ptr)); ptr++;}
      else if(CombineMode==AbsMin) for (int j=0; j<NumImportIDs; j++) {To[0][ImportLIDs[j]] = EPETRA_MIN( To[0][ImportLIDs[j]],std::abs(*ptr)); ptr++;}
      else if(CombineMode==Epetra_Max) for (int j=0; j<NumImportIDs; j++) {To[0][ImportLIDs[j]] = EPETRA_MAX( To[0][ImportLIDs[j]],*ptr); ptr++;}
      else if(CombineMode==Epetra_Min) for (int j=0; j<NumImportIDs; j++) {To[0][ImportLIDs[j]] = EPETRA_MIN( To[0][ImportLIDs[j]],*ptr); ptr++;}
      else if(CombineMode==Average) for (int j=0; j<NumImportIDs; j++) {To[0][ImportLIDs[j]] += *ptr++; To[0][ImportLIDs[j]] *= 0.5;} // Not a true avg if >2 occurance of an ID
      // Note:  The following form of averaging is not a true average if more that one value is combined.
      //        This might be an issue in the future, but we leave this way for now.
/*
      if (CombineMode==Add)
  for (int j=0; j<NumImportIDs; j++) To[0][ImportLIDs[j]] += *ptr++; // Add to existing value
      else if(CombineMode==Insert)
  for (int j=0; j<NumImportIDs; j++) To[0][ImportLIDs[j]] = *ptr++;
      else if(CombineMode==InsertAdd) {
  for (int j=0; j<NumImportIDs; j++) To[0][ImportLIDs[j]] = 0.0;
  for (int j=0; j<NumImportIDs; j++) To[0][ImportLIDs[j]] += *ptr++;
      }
      else if(CombineMode==AbsMax)
        for (int j=0; j<NumImportIDs; j++) {
    To[0][ImportLIDs[j]] = EPETRA_MAX( To[0][ImportLIDs[j]],std::abs(*ptr));
    ptr++;
  }
      // Note:  The following form of averaging is not a true average if more that one value is combined.
      //        This might be an issue in the future, but we leave this way for now.
      else if(CombineMode==Average)
  for (int j=0; j<NumImportIDs; j++) {To[0][ImportLIDs[j]] += *ptr++; To[0][ImportLIDs[j]] *= 0.5;}
*/
    }

    else {  // numVectors>1

      for (int j=0; j<NumImportIDs; j++) {
        jj = ImportLIDs[j];
        for (int i=0; i<numVectors; i++) {
        if (CombineMode==InsertAdd) To[i][jj] = 0.0; // Zero out if needed
        if (CombineMode==Add || CombineMode==InsertAdd) To[i][jj] += *ptr++; // Add to existing value
        else if (CombineMode==Insert) To[i][jj] = *ptr++; // Insert values
        else if (CombineMode==AbsMax) {To[i][jj] = EPETRA_MAX( To[i][jj], std::abs(*ptr)); ptr++; } // max of absolutes
        else if (CombineMode==AbsMin) {To[i][jj] = EPETRA_MIN( To[i][jj], std::abs(*ptr)); ptr++; } // max of absolutes
        else if (CombineMode==Epetra_Max)    {To[i][jj] = EPETRA_MAX( To[i][jj], *ptr); ptr++; } // simple max
        else if (CombineMode==Epetra_Min)    {To[i][jj] = EPETRA_MIN( To[i][jj], *ptr); ptr++; } // simple min
        else if (CombineMode==Average){To[i][jj] += *ptr++; To[i][jj] *= 0.5;}} // Not a true avg if >2 occurance of an ID

      }
/*
      if (CombineMode==Add) {
  for (int j=0; j<NumImportIDs; j++) {
    jj = ImportLIDs[j];
    for (int i=0; i<numVectors; i++)
      To[i][jj] += *ptr++; // Add to existing value
  }
      }
      else if(CombineMode==Insert) {
  for (int j=0; j<NumImportIDs; j++) {
    jj = ImportLIDs[j];
    for (int i=0; i<numVectors; i++)
      To[i][jj] = *ptr++;
  }
      }
      else if(CombineMode==InsertAdd) {
  for (int j=0; j<NumImportIDs; j++) {
    jj = ImportLIDs[j];
    for (int i=0; i<numVectors; i++)
      To[i][jj] = 0.0;
  }
  for (int j=0; j<NumImportIDs; j++) {
    jj = ImportLIDs[j];
    for (int i=0; i<numVectors; i++)
      To[i][jj] += *ptr++;
  }
      }
      else if(CombineMode==AbsMax) {
        for (int j=0; j<NumImportIDs; j++) {
          jj = ImportLIDs[j];
          for (int i=0; i<numVectors; i++) {
            To[i][jj] = EPETRA_MAX( To[i][jj], std::abs(*ptr) );
      ptr++;
    }
        }
      }
      // Note:  The following form of averaging is not a true average if more that one value is combined.
      //        This might be an issue in the future, but we leave this way for now.
      else if(CombineMode==Average) {
  for (int j=0; j<NumImportIDs; j++) {
    jj = ImportLIDs[j];
    for (int i=0; i<numVectors; i++)
      { To[i][jj] += *ptr++;  To[i][jj] *= 0.5;}
  }
      }
*/
    }
  }

  // constant element size case

  else if (ConstantElementSize) {

    for (int j=0; j<NumImportIDs; j++) {
      jj = MaxElementSize*ImportLIDs[j];
      for (int i=0; i<numVectors; i++) {
        if (CombineMode==InsertAdd) for (k=0; k<MaxElementSize; k++) To[i][jj+k] = 0.0; // Zero out if needed
        if (CombineMode==Add || CombineMode==InsertAdd) for (k=0; k<MaxElementSize; k++) To[i][jj+k] += *ptr++; // Add to existing value
        else if (CombineMode==Insert) for (k=0; k<MaxElementSize; k++) To[i][jj+k] = *ptr++; // Insert values
        else if (CombineMode==AbsMax) {for (k=0; k<MaxElementSize; k++) { To[i][jj+k] = EPETRA_MAX( To[i][jj+k], std::abs(*ptr)); ptr++; }} // max of absolutes
        else if (CombineMode==AbsMin) {for (k=0; k<MaxElementSize; k++) { To[i][jj+k] = EPETRA_MIN( To[i][jj+k], std::abs(*ptr)); ptr++; }} // max of absolutes
        else if (CombineMode==Epetra_Max) {for (k=0; k<MaxElementSize; k++) { To[i][jj+k] = EPETRA_MAX( To[i][jj+k], *ptr); ptr++; }} // simple max
        else if (CombineMode==Epetra_Min) {for (k=0; k<MaxElementSize; k++) { To[i][jj+k] = EPETRA_MIN( To[i][jj+k], *ptr); ptr++; }} // simple min
        else if (CombineMode==Average) {for (k=0; k<MaxElementSize; k++) { To[i][jj+k] += *ptr++; To[i][jj+k] *= 0.5;}} // Not a true avg if >2 occurance of an ID
     }
   }
/*
    if (CombineMode==Add) {
      for (int j=0; j<NumImportIDs; j++) {
  jj = MaxElementSize*ImportLIDs[j];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<MaxElementSize; k++)
      To[i][jj+k] += *ptr++; // Add to existing value
      }
    }
    else if(CombineMode==Insert) {
      for (int j=0; j<NumImportIDs; j++) {
  jj = MaxElementSize*ImportLIDs[j];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<MaxElementSize; k++)
      To[i][jj+k] = *ptr++;
      }
    }
    else if(CombineMode==InsertAdd) {
      for (int j=0; j<NumImportIDs; j++) {
  jj = MaxElementSize*ImportLIDs[j];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<MaxElementSize; k++)
      To[i][jj+k] = 0.0;
      }
      for (int j=0; j<NumImportIDs; j++) {
  jj = MaxElementSize*ImportLIDs[j];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<MaxElementSize; k++)
      To[i][jj+k] += *ptr++;
      }
    }
    else if(CombineMode==AbsMax) {
      for (int j=0; j<NumImportIDs; j++) {
  jj = MaxElementSize*ImportLIDs[j];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<MaxElementSize; k++) {
      To[i][jj+k] = EPETRA_MAX( To[i][jj+k], std::abs(*ptr) );
      ptr++;
    }
      }
    }
    // Note:  The following form of averaging is not a true average if more that one value is combined.
    //        This might be an issue in the future, but we leave this way for now.
    else if(CombineMode==Average) {
      for (int j=0; j<NumImportIDs; j++) {
  jj = MaxElementSize*ImportLIDs[j];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<MaxElementSize; k++)
      { To[i][jj+k] += *ptr++; To[i][jj+k] *= 0.5;}
      }
    }
*/
  }

  // variable element size case

  else {
    int thisSizeOfPacket = numVectors*MaxElementSize;

    for (int j=0; j<NumImportIDs; j++) {
      ptr = (int *) Imports + j*thisSizeOfPacket;
      jj = ToFirstPointInElementList[ImportLIDs[j]];
      int ElementSize = ToElementSizeList[ImportLIDs[j]];
      for (int i=0; i<numVectors; i++) {
        if (CombineMode==InsertAdd) for (k=0; k<ElementSize; k++) To[i][jj+k] = 0.0; // Zero out if needed
        if (CombineMode==Add || CombineMode==InsertAdd) for (k=0; k<ElementSize; k++) To[i][jj+k] += *ptr++; // Add to existing value
        else if (CombineMode==Insert) for (k=0; k<ElementSize; k++) To[i][jj+k] = *ptr++; // Insert values
        else if (CombineMode==AbsMax) {for (k=0; k<ElementSize; k++) { To[i][jj+k] = EPETRA_MAX( To[i][jj+k], std::abs(*ptr)); ptr++; }} // max of absolutes
        else if (CombineMode==Epetra_Max) {for (k=0; k<ElementSize; k++) { To[i][jj+k] = EPETRA_MAX( To[i][jj+k], *ptr); ptr++; }} // simple max
        else if (CombineMode==Epetra_Min) {for (k=0; k<ElementSize; k++) { To[i][jj+k] = EPETRA_MIN( To[i][jj+k], *ptr); ptr++; }} // simple min
        else if (CombineMode==Average) {for (k=0; k<ElementSize; k++) { To[i][jj+k] += *ptr++; To[i][jj+k] *= 0.5;}} // Not a true avg if >2 occurance of an ID
     }
   }
/*
    if (CombineMode==Add) {
      for (int j=0; j<NumImportIDs; j++) {
  ptr = (double *) Imports + j*thisSizeOfPacket;
  jj = ToFirstPointInElementList[ImportLIDs[j]];
  int ElementSize = ToElementSizeList[ImportLIDs[j]];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<ElementSize; k++)
      To[i][jj+k] += *ptr++; // Add to existing value
      }
    }
    else  if(CombineMode==Insert){
      for (int j=0; j<NumImportIDs; j++) {
  ptr = (double *) Imports + j*thisSizeOfPacket;
  jj = ToFirstPointInElementList[ImportLIDs[j]];
  int ElementSize = ToElementSizeList[ImportLIDs[j]];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<ElementSize; k++)
      To[i][jj+k] = *ptr++;
      }
    }
    else  if(CombineMode==InsertAdd){
      for (int j=0; j<NumImportIDs; j++) {
  ptr = (double *) Imports + j*thisSizeOfPacket;
  jj = ToFirstPointInElementList[ImportLIDs[j]];
  int ElementSize = ToElementSizeList[ImportLIDs[j]];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<ElementSize; k++)
      To[i][jj+k] = 0.0;
      }
      for (int j=0; j<NumImportIDs; j++) {
  ptr = (double *) Imports + j*thisSizeOfPacket;
  jj = ToFirstPointInElementList[ImportLIDs[j]];
  int ElementSize = ToElementSizeList[ImportLIDs[j]];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<ElementSize; k++)
      To[i][jj+k] += *ptr++;
      }
    }
    else  if(CombineMode==AbsMax){
      for (int j=0; j<NumImportIDs; j++) {
  ptr = (double *) Imports + j*thisSizeOfPacket;
  jj = ToFirstPointInElementList[ImportLIDs[j]];
  int ElementSize = ToElementSizeList[ImportLIDs[j]];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<ElementSize; k++) {
      To[i][jj+k] = EPETRA_MAX( To[i][jj+k], std::abs(*ptr));
      ptr++;
    }
      }
    }
    // Note:  The following form of averaging is not a true average if more that one value is combined.
    //        This might be an issue in the future, but we leave this way for now.
    else if(CombineMode==Average) {
      for (int j=0; j<NumImportIDs; j++) {
  ptr = (double *) Imports + j*thisSizeOfPacket;
  jj = ToFirstPointInElementList[ImportLIDs[j]];
  int ElementSize = ToElementSizeList[ImportLIDs[j]];
  for (int i=0; i<numVectors; i++)
    for (k=0; k<ElementSize; k++)
      { To[i][jj+k] += *ptr++; To[i][jj+k] *= 0.5;}
      }
    }
*/
  }

  return(0);
}

//=========================================================================
int  Epetra_IntMultiVector::MinValue (int* Result) const {

  // Minimum value of each vector in MultiVector

  int ierr = 0;

  const int myLength = MyLength_;
  UpdateOrdinalTemp();

  for (int i=0; i < NumVectors_; i++)
    {
      const int * const from = Pointers_[i];
      int MinVal = 2000000000; // 2 billion is close to largest 32 bit int
      if (myLength>0) MinVal = from[0];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel default(none) shared(MinVal)
{
      int localMinVal = MinVal;
#pragma omp for
      for (int j=0; j< myLength; j++) localMinVal = EPETRA_MIN(localMinVal,from[j]);
#pragma omp critical
      {
      MinVal = EPETRA_MIN(MinVal,localMinVal);
      }
}
      OrdinalTemp_[i] = MinVal;
#else
      for (int j=0; j< myLength; j++) MinVal = EPETRA_MIN(MinVal,from[j]);
      OrdinalTemp_[i] = MinVal;
#endif
    }

  if (myLength > 0) {
    for(int i=0; i<NumVectors_; ++i) {
      Result[i] = OrdinalTemp_[i];
    }
  }

  //If myLength == 0 and Comm_->NumProc() == 1, then Result has
  //not been referenced. Also, if vector contents are uninitialized
  //then Result contents are not well defined...

  if (Comm_->NumProc() == 1 || !DistributedGlobal()) return(ierr);

  //We're going to use MPI_Allgather to gather every proc's local-
  //min values onto every other proc. We'll use the last position
  //of the OrdinalTemp_ array to indicate whether this proc has
  //valid data that should be considered by other procs when forming
  //the global-min results.

  if (myLength == 0) OrdinalTemp_[NumVectors_] = 0;
  else OrdinalTemp_[NumVectors_] = 1;

  //Now proceed to handle the parallel case. We'll gather local-min
  //values from other procs and form a global-min. If any processor
  //has myLength>0, we'll end up with a valid result.

#ifdef EPETRA_MPI
  const Epetra_MpiComm* epetrampicomm =
    dynamic_cast<const Epetra_MpiComm*>(Comm_);
  if (!epetrampicomm) {
    return(-2);
  }

  MPI_Comm mpicomm = epetrampicomm->GetMpiComm();
  int numProcs = epetrampicomm->NumProc();
  int* owork = new int[numProcs*(NumVectors_+1)];

  MPI_Allgather(OrdinalTemp_, NumVectors_+1, MPI_INT,
                owork, NumVectors_+1, MPI_INT, mpicomm);

  //if myLength==0, then our Result array currently contains
  //Epetra_MaxOrdinal from the local-min calculations above. In this
  //case we'll overwrite our Result array with values from the first
  //processor that sent valid data.
  bool overwrite = myLength == 0 ? true : false;

  int myPID = epetrampicomm->MyPID();
  int* owork_ptr = owork;

  for(int j=0; j<numProcs; ++j) {

    //skip data from self, and skip data from
    //procs with OrdinalTemp_[NumVectors_] == 0.
    if (j == myPID || owork_ptr[NumVectors_] == 0) {
      owork_ptr += NumVectors_+1;
      continue;
    }

    for(int i=0; i<NumVectors_; ++i) {
      int val = owork_ptr[i];

      //Set val into our Result array if overwrite is true (see above),
      //or if val is less than our current Result[i].
      if (overwrite || (Result[i] > val)) Result[i] = val;
    }

    //Now set overwrite to false so that we'll do the right thing
    //when processing data from subsequent processors.
    if (overwrite) overwrite = false;

    owork_ptr += NumVectors_+1;
  }

  delete [] owork;
#endif

  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine

  return(ierr);
}

//=========================================================================
int  Epetra_IntMultiVector::MaxValue (int* Result) const {

  // Maximum value of each vector in MultiVector

  int ierr = 0;

  const int myLength = MyLength_;
  UpdateOrdinalTemp();

  for (int i=0; i < NumVectors_; i++)
    {
      const int * const from = Pointers_[i];
      int MaxVal = -2000000000; // Negative 2 billion is close to smallest 32 bit int
      if (myLength>0) MaxVal = from[0];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel default(none) shared(MaxVal)
{
      int localMaxVal = MaxVal;
#pragma omp for
      for (int j=0; j< myLength; j++) localMaxVal = EPETRA_MAX(localMaxVal,from[j]);
#pragma omp critical
      {
      MaxVal = EPETRA_MAX(MaxVal,localMaxVal);
      }
}
      OrdinalTemp_[i] = MaxVal;
#else
      for (int j=0; j< myLength; j++) MaxVal = EPETRA_MAX(MaxVal,from[j]);
      OrdinalTemp_[i] = MaxVal;
#endif
    }

  if (myLength > 0) {
    for(int i=0; i<NumVectors_; ++i) {
      Result[i] = OrdinalTemp_[i];
    }
  }

  //If myLength == 0 and Comm_->NumProc() == 1, then Result has
  //not been referenced. Also, if vector contents are uninitialized
  //then Result contents are not well defined...

  if (Comm_->NumProc() == 1  || !DistributedGlobal()) return(ierr);

  //We're going to use MPI_Allgather to gather every proc's local-
  //max values onto every other proc. We'll use the last position
  //of the OrdinalTemp_ array to indicate whether this proc has
  //valid data that should be considered by other procs when forming
  //the global-max results.

  if (myLength == 0) OrdinalTemp_[NumVectors_] = 0.0;
  else OrdinalTemp_[NumVectors_] = 1.0;

  //Now proceed to handle the parallel case. We'll gather local-max
  //values from other procs and form a global-max. If any processor
  //has myLength>0, we'll end up with a valid result.

#ifdef EPETRA_MPI
  const Epetra_MpiComm* epetrampicomm =
    dynamic_cast<const Epetra_MpiComm*>(Comm_);
  if (!epetrampicomm) {
    return(-2);
  }

  MPI_Comm mpicomm = epetrampicomm->GetMpiComm();
  int numProcs = epetrampicomm->NumProc();
  int* owork = new int[numProcs*(NumVectors_+1)];

  MPI_Allgather(OrdinalTemp_, NumVectors_+1, MPI_INT,
                owork, NumVectors_+1, MPI_INT, mpicomm);

  //if myLength==0, then our Result array currently contains
  //-Epetra_MaxOrdinal from the local-max calculations above. In this
  //case we'll overwrite our Result array with values from the first
  //processor that sent valid data.
  bool overwrite = myLength == 0 ? true : false;

  int myPID = epetrampicomm->MyPID();
  int* owork_ptr = owork;

  for(int j=0; j<numProcs; ++j) {

    //skip data from self, and skip data from
    //procs with OrdinalTemp_[NumVectors_] == 0.
    if (j == myPID || owork_ptr[NumVectors_] == 0) {
      owork_ptr += NumVectors_+1;
      continue;
    }

    for(int i=0; i<NumVectors_; ++i) {
      int val = owork_ptr[i];

      //Set val into our Result array if overwrite is true (see above),
      //or if val is larger than our current Result[i].
      if (overwrite || (Result[i] < val)) Result[i] = val;
    }

    //Now set overwrite to false so that we'll do the right thing
    //when processing data from subsequent processors.
    if (overwrite) overwrite = false;

    owork_ptr += NumVectors_+1;
  }

  delete [] owork;
#endif

  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine

  return(ierr);
}


//=======================================================================
Epetra_IntVector *& Epetra_IntMultiVector::operator () (int index)  {

  //  Epetra_IntMultiVector::operator () --- return non-const reference

  if (index < 0 || index >=NumVectors_)
    throw ReportError("Vector index = " + toString(index) + "is out of range. Number of Vectors = " + toString(NumVectors_), -1);

  UpdateIntVectors();

  // Create a new Epetra_IntVector that is a view of ith vector, if not already present
  if (IntVectors_[index]==0)
    IntVectors_[index] = new Epetra_IntVector(View, Map(), Pointers_[index]);
  return(IntVectors_[index]);
}

//=======================================================================
const Epetra_IntVector *& Epetra_IntMultiVector::operator () (int index) const {

  //  Epetra_IntMultiVector::operator () --- return non-const reference

  if (index < 0 || index >=NumVectors_)
    throw ReportError("Vector index = " + toString(index) + "is out of range. Number of Vectors = " + toString(NumVectors_), -1);

  UpdateIntVectors();

  if (IntVectors_[index]==0)
    IntVectors_[index] = new Epetra_IntVector(View, Map(), Pointers_[index]);

  const Epetra_IntVector * & temp = (const Epetra_IntVector * &) (IntVectors_[index]);
  return(temp);
}

//========================================================================
Epetra_IntMultiVector& Epetra_IntMultiVector::operator = (const Epetra_IntMultiVector& Source) {

  // Check for special case of this=Source
  if (this != &Source) Assign(Source);

  return(*this);
}

//=========================================================================
void Epetra_IntMultiVector::Assign(const Epetra_IntMultiVector& A) {

  const int myLength = MyLength_;
  if (NumVectors_ != A.NumVectors())
    throw ReportError("Number of vectors incompatible in Assign.  The this MultiVector has NumVectors = " + toString(NumVectors_)
          + ".  The A MultiVector has NumVectors = " + toString(A.NumVectors()), -3);
  if (myLength != A.MyLength())
    throw ReportError("Length of MultiVectors incompatible in Assign.  The this MultiVector has MyLength = " + toString(myLength)
          + ".  The A MultiVector has MyLength = " + toString(A.MyLength()), -4);

  int ** A_Pointers = A.Pointers();
  for (int i = 0; i< NumVectors_; i++) {
      int * const to = Pointers_[i];
      const int * const from = A_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none)
#endif
      for (int j=0; j<myLength; j++) to[j] = from[j];
    }
    return;
  }

//=========================================================================
int  Epetra_IntMultiVector::Reduce() {

  // Global reduction on each entry of a Replicated Local MultiVector
  const int myLength = MyLength_;
  int * source = 0;
  if (myLength>0) source = new int[myLength*NumVectors_];
  int * target = 0;
  bool packed = (ConstantStride_ && (Stride_==myLength));
  if (packed) {
    for (int i=0; i<myLength*NumVectors_; i++) source[i] = Values_[i];
    target = Values_;
  }
  else {
    int * tmp1 = source;
    for (int i = 0; i < NumVectors_; i++) {
      int * tmp2 = Pointers_[i];
      for (int j=0; j< myLength; j++) *tmp1++ = *tmp2++;
    }
    if (myLength>0) target = new int[myLength*NumVectors_];
  }

  Comm_->SumAll(source, target, myLength*NumVectors_);
  if (myLength>0) delete [] source;
  if (!packed) {
    int * tmp2 = target;
    for (int i = 0; i < NumVectors_; i++) {
      int * tmp1 = Pointers_[i];
      for (int j=0; j< myLength; j++) *tmp1++ = *tmp2++;
    }
    if (myLength>0) delete [] target;
  }
  // UpdateFlops(0);  No serial Flops in this function
  return(0);
}

//=======================================================================
int Epetra_IntMultiVector::ResetView(int ** ArrayOfPointers) {

  if (!UserAllocated_) {
    EPETRA_CHK_ERR(-1); // Can't reset view if multivector was not allocated as a view
  }

  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = ArrayOfPointers[i];
  DoView();

  return(0);
}

//=======================================================================
void Epetra_IntMultiVector::Print(std::ostream& os) const {
  int MyPID = Map().Comm().MyPID();
  int NumProc = Map().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumVectors1 = NumVectors();
      int NumMyElements1 =Map(). NumMyElements();
      int MaxElementSize1 = Map().MaxElementSize();
      int * FirstPointInElementList1 = NULL;
      if (MaxElementSize1!=1) FirstPointInElementList1 = Map().FirstPointInElementList();
      int ** A_Pointers = Pointers();

      if (MyPID==0) {
  os.width(8);
  os <<  "     MyPID"; os << "    ";
  os.width(12);
  if (MaxElementSize1==1)
    os <<  "GID  ";
  else
    os <<  "     GID/Point";
  for (int j = 0; j < NumVectors1 ; j++)
    {
      os.width(20);
      os <<  "Value  ";
    }
  os << std::endl;
      }
      for (int i=0; i < NumMyElements1; i++) {
  for (int ii=0; ii< Map().ElementSize(i); ii++) {
       int iii;
    os.width(10);
    os <<  MyPID; os << "    ";
    os.width(10);
    if (MaxElementSize1==1) {
      if(Map().GlobalIndicesInt())
      {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      int * MyGlobalElements1 = Map().MyGlobalElements();
      os << MyGlobalElements1[i] << "    ";
#else
            throw ReportError("Epetra_IntMultiVector::Print: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
      }
      else if(Map().GlobalIndicesLongLong())
      {
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
      long long * MyGlobalElements1 = Map().MyGlobalElements64();
      os << MyGlobalElements1[i] << "    ";
#else
            throw ReportError("Epetra_IntMultiVector::Print: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
      }
      else
        throw ReportError("Epetra_IntMultiVector::Print ERROR, Don't know map global index type.",-1);

       iii = i;
       }
    else {
      if(Map().GlobalIndicesInt())
      {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      int * MyGlobalElements1 = Map().MyGlobalElements();
      os <<  MyGlobalElements1[i]<< "/" << ii << "    ";
#else
            throw ReportError("Epetra_IntMultiVector::Print: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
      }
      else if(Map().GlobalIndicesLongLong())
      {
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
      long long * MyGlobalElements1 = Map().MyGlobalElements64();
      os <<  MyGlobalElements1[i]<< "/" << ii << "    ";
#else
            throw ReportError("Epetra_IntMultiVector::Print: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
      }
      else
        throw ReportError("Epetra_IntMultiVector::Print ERROR, Don't know map global index type.",-1);

         iii = FirstPointInElementList1[i]+ii;
       }
    for (int j = 0; j < NumVectors1 ; j++)
      {
        os.width(20);
        os <<  A_Pointers[j][iii];
      }
    os << std::endl;
  }
      }
      os << std::flush;
    }

    // Do a few global ops to give I/O a chance to complete
    Map().Comm().Barrier();
    Map().Comm().Barrier();
    Map().Comm().Barrier();
  }
  return;
}
