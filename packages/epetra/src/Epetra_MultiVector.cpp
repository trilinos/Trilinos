
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

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
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

Epetra_MultiVector::Epetra_MultiVector(const Epetra_BlockMap& map, int numVectors, bool zeroOut)
  : Epetra_DistObject(map, "Epetra::MultiVector"),
    Epetra_CompObject(),
    Values_(0),
    Pointers_(0),
    MyLength_(map.NumMyPoints()),
    GlobalLength_(map.NumGlobalPoints()),
    NumVectors_(numVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(map.NumMyPoints()),
    Allocated_(false)
{
	Util_.SetSeed(1);

    AllocateForCopy();
    
    for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Values_+i*Stride_;

	if(zeroOut) PutScalar(0.0); // Fill all vectors with zero.
}
//==========================================================================

// Copy Constructor

Epetra_MultiVector::Epetra_MultiVector(const Epetra_MultiVector& Source)
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
  
  double ** Source_Pointers = Source.Pointers();
  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Source_Pointers[i];
  
  DoCopy();
  
}
//==========================================================================

// This constructor copies in or makes view of a standard Fortran array

Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& map, 
					     double *A, int MyLDA, int numVectors)
  : Epetra_DistObject(map, "Epetra::MultiVector"),
    Epetra_CompObject(),
    Values_(0),
    Pointers_(0),
    MyLength_(map.NumMyPoints()),
    GlobalLength_(map.NumGlobalPoints()),
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

Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& map, 
					      double **ArrayOfPointers, int numVectors)
  : Epetra_DistObject(map, "Epetra::MultiVector"),
    Epetra_CompObject(),
    Values_(0),
    Pointers_(0),
    MyLength_(map.NumMyPoints()),
    GlobalLength_(map.NumGlobalPoints()),
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

Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, 
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

  double ** Source_Pointers = Source.Pointers();  
  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Source_Pointers[Indices[i]];
  
   if (CV==Copy) DoCopy();
   else DoView();
  
}

//==========================================================================

// This interface copies or makes view of a range of vectors from an existing MultiVector

Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, 
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

  double ** Source_Pointers = Source.Pointers();
  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = Source_Pointers[StartIndex+i];

   if (CV==Copy) DoCopy();
   else DoView();
}

//=========================================================================
Epetra_MultiVector::~Epetra_MultiVector(){

  if (!Allocated_) return;

  delete [] Pointers_;
  if (!UserAllocated_ && Values_!=0) delete [] Values_;

  if (Vectors_!=0) {
    for (int i=0; i<NumVectors_; i++) if (Vectors_[i]!=0) delete Vectors_[i];
    delete [] Vectors_;
  }


  if (DoubleTemp_!=0) delete [] DoubleTemp_;

}

//=========================================================================
int Epetra_MultiVector::AllocateForCopy(void)
{
  
  if (Allocated_) return(0);
    
  if (NumVectors_<=0) 
    throw ReportError("Number of Vectors = " + toString(NumVectors_) + ", but must be greater than zero", -1);

  Stride_ = Map_.NumMyPoints();
  if (Stride_>0) Values_ = new double[Stride_ * NumVectors_];
  Pointers_ = new double *[NumVectors_];

  DoubleTemp_ = 0;
  Vectors_ = 0;
  
  int randval = rand(); // Use POSIX standard random function
  if (DistributedGlobal_)
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
int Epetra_MultiVector::DoCopy(void)
{
  // On entry Pointers_ contains pointers to the incoming vectors.  These
  // pointers are the only unique piece of information for each of the 
  // constructors.

  // \internal { Optimization of this function can impact performance since it 
  //           involves a fair amount of memory traffic.}

  // On exit, Pointers_ is redefined to point to its own MultiVector vectors.

  for (int i = 0; i< NumVectors_; i++)
    {
      double * from = Pointers_[i];
      double * to = Values_+i*Stride_;
      Pointers_[i] = to;
      const int myLength = MyLength_;
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(to,from)
#endif
      for (int j=0; j<myLength; j++) to[j] = from[j];
    }

  return(0);
}
//=========================================================================
int Epetra_MultiVector::AllocateForView(void)
{
  
  if (NumVectors_<=0) 
    throw ReportError("Number of Vectors = " + toString(NumVectors_) + ", but must be greater than zero", -1);
 
  Pointers_ = new double *[NumVectors_];
  
  DoubleTemp_ = 0;
  Vectors_ = 0;
  
  int randval = rand(); // Use POSIX standard random function
  if (DistributedGlobal_)
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
int Epetra_MultiVector::DoView(void)
{
  // On entry Pointers_ contains pointers to the incoming vectors.  These
  // pointers are the only unique piece of information for each of the 
  // constructors.


  Values_ = Pointers_[0];

  if (NumVectors_==1) {
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
int Epetra_MultiVector::ReplaceGlobalValue(long long GlobalRow, int VectorIndex, double ScalarValue) {

 // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue(GlobalRow, 0, VectorIndex, ScalarValue, false));
  return(0);
}
//=========================================================================
int Epetra_MultiVector::ReplaceGlobalValue(long long GlobalBlockRow, int BlockRowOffset, 
					   int VectorIndex, double ScalarValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue(GlobalBlockRow, BlockRowOffset, VectorIndex, ScalarValue, false)); 
  return(0);
}
//=========================================================================
int Epetra_MultiVector::SumIntoGlobalValue(long long GlobalRow, int VectorIndex, double ScalarValue) {

  // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue(GlobalRow, 0, VectorIndex, ScalarValue, true)); 
  return(0);
}
//=========================================================================
int Epetra_MultiVector::SumIntoGlobalValue(int GlobalBlockRow, int BlockRowOffset, 
					   int VectorIndex, double ScalarValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue(GlobalBlockRow, BlockRowOffset, VectorIndex, ScalarValue, true));
  return(0);
}
//=========================================================================
int Epetra_MultiVector::ReplaceMyValue(int MyRow, int VectorIndex, double ScalarValue) {

  // Use the more general method below
  EPETRA_CHK_ERR(ChangeMyValue(MyRow, 0, VectorIndex, ScalarValue, false)); 
  return(0);
}
//=========================================================================
int Epetra_MultiVector::ReplaceMyValue(int MyBlockRow, int BlockRowOffset, 
					   int VectorIndex, double ScalarValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeMyValue(MyBlockRow, BlockRowOffset, VectorIndex, ScalarValue, false)); 
  return(0);
}
//=========================================================================
int Epetra_MultiVector::SumIntoMyValue(int MyRow, int VectorIndex, double ScalarValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeMyValue(MyRow, 0, VectorIndex, ScalarValue, true)); 
  return(0);
}
//=========================================================================
int Epetra_MultiVector::SumIntoMyValue(int MyBlockRow, int BlockRowOffset, 
					   int VectorIndex, double ScalarValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeMyValue(MyBlockRow, BlockRowOffset, VectorIndex, ScalarValue, true));
  return(0);
}
//=========================================================================
int Epetra_MultiVector::ChangeGlobalValue(long long GlobalBlockRow, int BlockRowOffset, 
				     int VectorIndex, double ScalarValue, bool SumInto) {

  // Convert GID to LID and call LID version
  EPETRA_CHK_ERR(ChangeMyValue(Map().LID(GlobalBlockRow), BlockRowOffset, VectorIndex, ScalarValue, SumInto));
  return(0);
}
//=========================================================================
int Epetra_MultiVector::ChangeMyValue(int MyBlockRow, int BlockRowOffset, 
				     int VectorIndex, double ScalarValue, bool SumInto) {
  
  if (!Map().MyLID(MyBlockRow)) EPETRA_CHK_ERR(1); // I don't own this one, return a warning flag
  if (VectorIndex>= NumVectors_) EPETRA_CHK_ERR(-1); // Consider this a real error
  if (BlockRowOffset<0 || BlockRowOffset>=Map().ElementSize(MyBlockRow)) EPETRA_CHK_ERR(-2); // Offset is out-of-range

  int entry = Map().FirstPointInElement(MyBlockRow);

  if (SumInto)
    Pointers_[VectorIndex][entry+BlockRowOffset] += ScalarValue;
  else
    Pointers_[VectorIndex][entry+BlockRowOffset] = ScalarValue;

  return(0);
}
//=========================================================================
int Epetra_MultiVector::Random() {
  // Generate random numbers drawn from a uniform distribution on
  // the interval (-1,1) using a multiplicative congruential generator
  // with modulus 2^31 - 1.
	/*  
  const double a = 16807.0, BigInt=2147483647.0, DbleOne=1.0, DbleTwo=2.0;
  
  for(int i=0; i < NumVectors_; i++)
    for (int j=0; j<MyLength_; j++){
      Seed_ = fmod( a*Seed_, BigInt );
      Pointers_[i][j] = DbleTwo*(Seed_/BigInt)-DbleOne;
    }
	*/

   const int myLength = MyLength_;
	for(int i = 0; i < NumVectors_; i++) {
          double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
#pragma omp parallel for default(none)
#endif
		for(int j = 0; j < myLength; j++)
			to[j] = Util_.RandomDouble();
        }

  return(0);
}
 
//=========================================================================

// Extract a copy of a Epetra_MultiVector.  Put in a user's Fortran-style array

int Epetra_MultiVector::ExtractCopy(double *A, int MyLDA) const {
  if (NumVectors_>1 && Stride_ > MyLDA) EPETRA_CHK_ERR(-1); // LDA not big enough
  
  const int myLength = MyLength_;
  for (int i=0; i< NumVectors_; i++)
    {
      double * from = Pointers_[i];
      double * to = A + i*MyLDA;
      for (int j=0; j<myLength; j++) *to++ = *from++;
    }

  return(0);
}

//=========================================================================
int Epetra_MultiVector::ReplaceMap(const Epetra_BlockMap& map)
{
  if (Map().PointSameAs(map)) {
    Epetra_DistObject::Map_ = map;
    return(0);
  }

  return(-1);  
}

//=========================================================================

// Extract a copy of a Epetra_MultiVector.  Put in a user's array of pointers

int Epetra_MultiVector::ExtractCopy(double **ArrayOfPointers) const {
  const int myLength = MyLength_;
  for (int i=0; i< NumVectors_; i++)
    {
      double * from = Pointers_[i];
      double * to = ArrayOfPointers[i];
      for (int j=0; j<myLength; j++) *to++ = *from++;
    }
  
  return(0);
}
      
      

//=========================================================================

// Extract a view of a Epetra_MultiVector.  Set up a user's Fortran-style array

int Epetra_MultiVector::ExtractView(double **A, int *MyLDA) const {
  if (!ConstantStride_) EPETRA_CHK_ERR(-1);  // Can't make a Fortran-style view if not constant stride
  *MyLDA = Stride_; // Set user's LDA
  *A = Values_; // Set user's value pointer
  return(0);
}
      
      

//=========================================================================

// Extract a view of a Epetra_MultiVector.  Put in a user's array of pointers

int Epetra_MultiVector::ExtractView(double ***ArrayOfPointers) const {
  *ArrayOfPointers = Pointers_;
  
  return(0);
}
      
      
//=========================================================================
int Epetra_MultiVector::PutScalar(double ScalarConstant) {

  // Fills MultiVector with the value ScalarConstant **/

  const int myLength = MyLength_;
  for (int i = 0; i < NumVectors_; i++) {
    double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarConstant)
#endif
    for (int j=0; j<myLength; j++) to[j] = ScalarConstant;
  }
  return(0);
}
//=========================================================================
int Epetra_MultiVector::CheckSizes(const Epetra_SrcDistObject& Source) {

  const Epetra_MultiVector & A = dynamic_cast<const Epetra_MultiVector &>(Source);

  if (NumVectors()!=A.NumVectors()) {EPETRA_CHK_ERR(-3)};
  return(0);
}

//=========================================================================
int Epetra_MultiVector::CopyAndPermute(const Epetra_SrcDistObject& Source,
                                       int NumSameIDs, 
                                       int NumPermuteIDs,
                                       int * PermuteToLIDs, 
                                       int *PermuteFromLIDs,
                                       const Epetra_OffsetIndex * Indexor)
{
  (void)Indexor;

  const Epetra_MultiVector & A = dynamic_cast<const Epetra_MultiVector &>(Source);

  double **From = A.Pointers();
  double **To = Pointers_;
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
      for (int i=0; i < numVectors; i++)
	for (int j=0; j<NumSameEntries; j++)
	  To[i][j] = From[i][j];
    }
  // Do local permutation next
  if (NumPermuteIDs>0) {
  
    // Point entry case
    if (Case1) {
      
      if (numVectors==1)
	for (int j=0; j<NumPermuteIDs; j++) 
	  To[0][PermuteToLIDs[j]] = From[0][PermuteFromLIDs[j]];
      
      else {
	for (int j=0; j<NumPermuteIDs; j++) {
	  jj = PermuteToLIDs[j];
	  jjj = PermuteFromLIDs[j];
	  for (int i=0; i<numVectors; i++)
	    To[i][jj] = From[i][jjj];
	}
      }
    }
    // constant element size case
    else if (Case2) {
      
      for (int j=0; j<NumPermuteIDs; j++) {
	jj = MaxElementSize*PermuteToLIDs[j];
	jjj = MaxElementSize*PermuteFromLIDs[j];
	for (int i=0; i<numVectors; i++)
	  for (k=0; k<MaxElementSize; k++)
	    To[i][jj+k] = From[i][jjj+k];
      }
    }
    
    // variable element size case
    else {
      
      for (int j=0; j<NumPermuteIDs; j++) {
	jj = ToFirstPointInElementList[PermuteToLIDs[j]];
	jjj = FromFirstPointInElementList[PermuteFromLIDs[j]];
	int ElementSize = FromElementSizeList[PermuteFromLIDs[j]];
	for (int i=0; i<numVectors; i++)
	  for (k=0; k<ElementSize; k++)
	    To[i][jj+k] = From[i][jjj+k];
      }
    }
  }
  return(0);
}

//=========================================================================
int Epetra_MultiVector::PackAndPrepare(const Epetra_SrcDistObject & Source,
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

  const Epetra_MultiVector & A = dynamic_cast<const Epetra_MultiVector &>(Source);
  int jj, k;

  double **From = A.Pointers();
  int MaxElementSize = Map().MaxElementSize();
  int numVectors = NumVectors_;
  bool ConstantElementSize = Map().ConstantElementSize();

  int * FromFirstPointInElementList = 0;
  int * FromElementSizeList = 0;

  if (!ConstantElementSize) {
    FromFirstPointInElementList = A.Map().FirstPointInElementList();
    FromElementSizeList = A.Map().ElementSizeList();
  }

  double * DoubleExports = 0;

  SizeOfPacket = numVectors*MaxElementSize*(int)sizeof(double);

  if(SizeOfPacket*NumExportIDs>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = SizeOfPacket*NumExportIDs;
    DoubleExports = new double[numVectors*MaxElementSize*NumExportIDs];
    Exports = (char *) DoubleExports;
  }

  double * ptr;

  if (NumExportIDs>0) {
    ptr = (double *) Exports;
    
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
	ptr = (double *) Exports + j*thisSizeOfPacket;
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
int Epetra_MultiVector::UnpackAndCombine(const Epetra_SrcDistObject & Source,
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

  double ** To = Pointers_;
  int numVectors = NumVectors_;
  int MaxElementSize = Map().MaxElementSize();
  bool ConstantElementSize = Map().ConstantElementSize();

  int * ToFirstPointInElementList = 0;
  int * ToElementSizeList = 0;

  if (!ConstantElementSize) {
    ToFirstPointInElementList = Map().FirstPointInElementList();
    ToElementSizeList = Map().ElementSizeList();
  }
  
  double * ptr;
  // Unpack it...

  ptr = (double *) Imports;
    
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
      ptr = (double *) Imports + j*thisSizeOfPacket;
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
int Epetra_MultiVector::Dot(const Epetra_MultiVector& A, double *Result) const {

  // Dot product of two MultiVectors 

  if (NumVectors_ != A.NumVectors()) EPETRA_CHK_ERR(-1);
  const int myLength = MyLength_;
  if (myLength != A.MyLength()) EPETRA_CHK_ERR(-2);
  UpdateDoubleTemp();
    
  double **A_Pointers = A.Pointers();

#ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
  for (int i=0; i < NumVectors_; i++) 
    {
      const double * const from = Pointers_[i];
      const double * const fromA = A_Pointers[i];
      double sum = 0.0;
#pragma omp parallel for reduction (+:sum) default(none)
      for (int j=0; j < myLength; j++) sum += from[j] * fromA[j];
      DoubleTemp_[i] = sum;
    }
#else
  for (int i=0; i < NumVectors_; i++) DoubleTemp_[i] = DOT(myLength, Pointers_[i], A_Pointers[i]);
#endif

  if (DistributedGlobal())
    Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  else
    for (int i=0; i< NumVectors_; ++i) Result[i] = DoubleTemp_[i];
  
  UpdateFlops(2*GlobalLength_*NumVectors_);

  return(0);
}
//=========================================================================
int Epetra_MultiVector::Abs(const Epetra_MultiVector& A) {

  // this[i][j] = std::abs(A[i][j])

  const int myLength = MyLength_;
  if (NumVectors_ != A.NumVectors()) EPETRA_CHK_ERR(-1);
  if (myLength != A.MyLength()) EPETRA_CHK_ERR(-2);

  double **A_Pointers = A.Pointers();

  for (int i=0; i < NumVectors_; i++) {
    double * const to = Pointers_[i];
    const double * const from = A_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none)
#endif
    for (int j=0; j < myLength; j++) to[j] = std::abs(from[j]);
  }

  return(0);
}
//=========================================================================
int Epetra_MultiVector::Reciprocal(const Epetra_MultiVector& A) {

  // this[i][j] = 1.0/(A[i][j])

  int ierr = 0;
#ifndef EPETRA_HAVE_OMP
  int localierr = 0;
#endif
  const int myLength = MyLength_;
  if (NumVectors_ != A.NumVectors()) EPETRA_CHK_ERR(-1);
  if (myLength != A.MyLength()) EPETRA_CHK_ERR(-2);

  double **A_Pointers = A.Pointers();

  for (int i=0; i < NumVectors_; i++) {
    double * const to = Pointers_[i];
    const double * const from = A_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel default(none) shared(ierr)
{
    int localierr = 0;
#pragma omp for
#endif
    for (int j=0; j < myLength; j++) {
      double value = from[j];
      // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)     
      if (std::abs(value)<Epetra_MinDouble) {
        if (value==0.0) localierr = 1;
        else if (localierr!=1) localierr = 2;
        to[j] = EPETRA_SGN(value) * Epetra_MaxDouble;
      }
      else
        to[j] = 1.0/value;
    }
#ifdef EPETRA_HAVE_OMP
#pragma omp critical 
#endif
{
   if (localierr==1) ierr = 1;
   else if (localierr==2 && ierr!=1) ierr = 2;
}
#ifdef EPETRA_HAVE_OMP
}
#endif
  }
  EPETRA_CHK_ERR(ierr);
  return(0);
}
  //=========================================================================
  int Epetra_MultiVector::Scale (double ScalarValue) {

    // scales a MultiVector in place by a scalar
  

  const int myLength = MyLength_;
#ifdef EPETRA_HAVE_OMP
    for (int i = 0; i < NumVectors_; i++) {
      double * const to = Pointers_[i];
#pragma omp parallel for default(none) shared(ScalarValue)
      for (int j = 0; j < myLength; j++) to[j] = ScalarValue * to[j];
    }
#else
    for (int i = 0; i < NumVectors_; i++)
      SCAL(myLength, ScalarValue, Pointers_[i]);
#endif
    UpdateFlops(GlobalLength_*NumVectors_);

    return(0);
  }

  //=========================================================================
  int Epetra_MultiVector::Scale (double ScalarA, const Epetra_MultiVector& A) {

    // scales a MultiVector by a scalar and put in the this:
    // this = ScalarA * A

  const int myLength = MyLength_;
  if (NumVectors_ != A.NumVectors()) EPETRA_CHK_ERR(-1);
  if (myLength != A.MyLength()) EPETRA_CHK_ERR(-2);

    double **A_Pointers = (double**)A.Pointers();
    
    for (int i = 0; i < NumVectors_; i++) {
      double * const to = Pointers_[i];
      const double * const from = A_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA)
#endif
      for (int j = 0; j < myLength; j++) to[j] = ScalarA * from[j];
    }
    UpdateFlops(GlobalLength_*NumVectors_);

    return(0);
  }

  //=========================================================================
  int Epetra_MultiVector::Update(double ScalarA, const Epetra_MultiVector& A, double ScalarThis) {


    // linear combination of two MultiVectors: this = ScalarThis * this + ScalarA * A

  const int myLength = MyLength_;
  if (NumVectors_ != A.NumVectors()) EPETRA_CHK_ERR(-1);
  if (myLength != A.MyLength()) EPETRA_CHK_ERR(-2);

    double **A_Pointers = (double**)A.Pointers();

    if (ScalarThis==0.0)
      {
	for (int i = 0; i < NumVectors_; i++) {
          double * const to = Pointers_[i];
          const double * const from = A_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA)
#endif
	  for (int j = 0; j < myLength; j++) to[j] = ScalarA * from[j];
        }
	UpdateFlops(GlobalLength_*NumVectors_);
      }
    else if (ScalarThis==1.0)
      {
	for (int i = 0; i < NumVectors_; i++) {
          double * const to = Pointers_[i];
          const double * const from = A_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA)
#endif
	  for (int j = 0; j < myLength; j++) to[j] = to[j] + ScalarA * from[j];
        }
	UpdateFlops(2*GlobalLength_*NumVectors_);
      }
    else if (ScalarA==1.0)
      {
	for (int i = 0; i < NumVectors_; i++) {
          double * const to = Pointers_[i];
          const double * const from = A_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarThis)
#endif
	  for (int j = 0; j < myLength; j++) to[j] = ScalarThis * to[j] + from[j];
        }
	UpdateFlops(2*GlobalLength_*NumVectors_);
      }
    else
      {
	for (int i = 0; i < NumVectors_; i++) {
          double * const to = Pointers_[i];
          const double * const from = A_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA,ScalarThis)
#endif
	  for (int j = 0; j < myLength; j++) to[j] = ScalarThis * to[j] +
					                    ScalarA *  from[j];
        }
	UpdateFlops(3*GlobalLength_*NumVectors_);
      }

    return(0);
  }

//=========================================================================
int Epetra_MultiVector::Update(double ScalarA, const Epetra_MultiVector& A, 
				  double ScalarB, const Epetra_MultiVector& B, double ScalarThis) {
  
  
  // linear combination of three MultiVectors: 
  // this = ScalarThis * this + ScalarA * A + ScalarB * B
  
  if (ScalarA==0.0) {
    EPETRA_CHK_ERR(Update(ScalarB, B, ScalarThis));
    return(0);
  }
  if (ScalarB==0.0) {
    EPETRA_CHK_ERR(Update(ScalarA, A, ScalarThis));
    return(0);
  }
			   
  const int myLength = MyLength_;
  if (NumVectors_ != A.NumVectors() || NumVectors_ != B.NumVectors()) EPETRA_CHK_ERR(-1);
  if (myLength != A.MyLength() || myLength != B.MyLength()) EPETRA_CHK_ERR(-2);
  
    double **A_Pointers = (double**)A.Pointers();
    double **B_Pointers = (double**)B.Pointers();

    if (ScalarThis==0.0)
      {
	if (ScalarA==1.0)
	  {
	    for (int i = 0; i < NumVectors_; i++) {
              double * const to = Pointers_[i];
              const double * const fromA = A_Pointers[i];
              const double * const fromB = B_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarB)
#endif
	      for (int j = 0; j < myLength; j++) to[j] =           fromA[j] + 
						                ScalarB * fromB[j];
            }
	    UpdateFlops(2*GlobalLength_*NumVectors_);
	  }
	else if (ScalarB==1.0)
	  {
	    for (int i = 0; i < NumVectors_; i++) {
              double * const to = Pointers_[i];
              const double * const fromA = A_Pointers[i];
              const double * const fromB = B_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA)
#endif
	      for (int j = 0; j < myLength; j++) to[j] = ScalarA * fromA[j] +
						                          fromB[j];
            }
	    UpdateFlops(2*GlobalLength_*NumVectors_);
	  }
	else
	  {
	    for (int i = 0; i < NumVectors_; i++) {
              double * const to = Pointers_[i];
              const double * const fromA = A_Pointers[i];
              const double * const fromB = B_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA,ScalarB)
#endif
	      for (int j = 0; j < myLength; j++) to[j] = ScalarA * fromA[j] + 
						                ScalarB * fromB[j];
            }
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
      }
    else if (ScalarThis==1.0)
      {
	if (ScalarA==1.0)
	  {
	    for (int i = 0; i < NumVectors_; i++) {
              double * const to = Pointers_[i];
              const double * const fromA = A_Pointers[i];
              const double * const fromB = B_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarB)
#endif
	      for (int j = 0; j < myLength; j++) to[j] +=           fromA[j] + 
						                 ScalarB * fromB[j];
            }
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
	else if (ScalarB==1.0)
	  {
	    for (int i = 0; i < NumVectors_; i++) {
              double * const to = Pointers_[i];
              const double * const fromA = A_Pointers[i];
              const double * const fromB = B_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA)
#endif
	      for (int j = 0; j < myLength; j++) to[j] += ScalarA * fromA[j] +
						                           fromB[j];
            }
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
	else
	  {
	    for (int i = 0; i < NumVectors_; i++) {
              double * const to = Pointers_[i];
              const double * const fromA = A_Pointers[i];
              const double * const fromB = B_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA,ScalarB)
#endif
	      for (int j = 0; j < myLength; j++) to[j] += ScalarA * fromA[j] + 
						                 ScalarB * fromB[j];
            }
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
      }
    else
      {
	if (ScalarA==1.0)
	  {
	    for (int i = 0; i < NumVectors_; i++) {
              double * const to = Pointers_[i];
              const double * const fromA = A_Pointers[i];
              const double * const fromB = B_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA,ScalarB,ScalarThis)
#endif
	      for (int j = 0; j < myLength; j++) to[j] =  ScalarThis *    to[j] +
						                           fromA[j] + 
						                 ScalarB * fromB[j];
            }
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
	else if (ScalarB==1.0)
	  {
	    for (int i = 0; i < NumVectors_; i++) {
              double * const to = Pointers_[i];
              const double * const fromA = A_Pointers[i];
              const double * const fromB = B_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA,ScalarThis)
#endif
	      for (int j = 0; j < myLength; j++) to[j] =  ScalarThis *    to[j] +
						                 ScalarA * fromA[j] +
						                           fromB[j];
            }
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
	else
	  {
	    for (int i = 0; i < NumVectors_; i++) {
              double * const to = Pointers_[i];
              const double * const fromA = A_Pointers[i];
              const double * const fromB = B_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarA,ScalarB,ScalarThis)
#endif
	      for (int j = 0; j < myLength; j++) to[j] =  ScalarThis *    to[j] +
						                 ScalarA * fromA[j] + 
						                 ScalarB * fromB[j];
            }
	    UpdateFlops(5*GlobalLength_*NumVectors_);
	  }
      }


    return(0);
  }

//=========================================================================
int  Epetra_MultiVector::Norm1 (double* Result) const {
  
  // 1-norm of each vector in MultiVector 
    
  if (!Map().UniqueGIDs()) {EPETRA_CHK_ERR(-1);}

  const int myLength = MyLength_;
  UpdateDoubleTemp();
#ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
  for (int i=0; i < NumVectors_; i++) 
    {
      const double * const from = Pointers_[i];
      double asum = 0.0;
#pragma omp parallel default(none) shared(asum)
{
      double localasum = 0.0;
#pragma omp for 
      for (int j=0; j< myLength; j++) localasum += std::abs(from[j]); 
#pragma omp critical
      asum += localasum;
}
      DoubleTemp_[i] = asum;
  }
#else

  for (int i=0; i < NumVectors_; i++) DoubleTemp_[i] = ASUM(myLength, Pointers_[i]);
#endif
  
  if (DistributedGlobal())
    Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  else
    for (int i=0; i< NumVectors_; ++i) Result[i] = DoubleTemp_[i];
  
  UpdateFlops(2*GlobalLength_*NumVectors_);

  return(0);
}

//=========================================================================
int  Epetra_MultiVector::Norm2 (double* Result) const {
  
  // 2-norm of each vector in MultiVector 
  

  if (!Map().UniqueGIDs()) {EPETRA_CHK_ERR(-1);}

  const int myLength = MyLength_;
  UpdateDoubleTemp();

  for (int i=0; i < NumVectors_; i++) 
    {
      const double * const from = Pointers_[i];
      double sum = 0.0;
#ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
#pragma omp parallel for reduction (+:sum) default(none)
#endif
      for (int j=0; j < myLength; j++) sum += from[j] * from[j];
      DoubleTemp_[i] = sum;
    }
  if (DistributedGlobal())
    Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  else
    for (int i=0; i< NumVectors_; ++i) Result[i] = DoubleTemp_[i];

  for (int i=0; i < NumVectors_; i++) Result[i] = std::sqrt(Result[i]);
  
  UpdateFlops(2*GlobalLength_*NumVectors_);
  
  return(0);
}

//=========================================================================
int  Epetra_MultiVector::NormInf (double* Result) const {
  
  // Inf-norm of each vector in MultiVector 
    

  const int myLength = MyLength_;
  UpdateDoubleTemp();

  for (int i=0; i < NumVectors_; i++) 
    {
      DoubleTemp_[i] = 0.0;
      double normval = 0.0;
      const double * const from = Pointers_[i];
      if (myLength>0) normval = from[0];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel default(none) shared(normval)
{
      double localnormval = 0.0;
#pragma omp for
      for (int j=0; j< myLength; j++) {
         localnormval = EPETRA_MAX(localnormval,std::abs(from[j])); 
      }
#pragma omp critical
      {
      normval = EPETRA_MAX(normval,localnormval);
      }
}
      DoubleTemp_[i] = normval;
#else
      int jj = IAMAX(myLength, Pointers_[i]);
      if (jj>-1) DoubleTemp_[i] = std::abs(Pointers_[i][jj]);
#endif
    }
  if (DistributedGlobal())
    Comm_->MaxAll(DoubleTemp_, Result, NumVectors_);
  else
    for (int i=0; i< NumVectors_; ++i) Result[i] = DoubleTemp_[i];
  
  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine  
  return(0);
}

//=========================================================================
int  Epetra_MultiVector::NormWeighted (const Epetra_MultiVector& Weights, double* Result) const {
  
  // Weighted 2-norm of each vector in MultiVector 

  // If only one vector in Weights, we assume it will be used as the weights for all vectors

  const int myLength = MyLength_;
  bool OneW = false;
  if (Weights.NumVectors()==1) OneW = true;
  else if (NumVectors_ != Weights.NumVectors()) EPETRA_CHK_ERR(-1);

  if (myLength != Weights.MyLength()) EPETRA_CHK_ERR(-2);


  UpdateDoubleTemp();

  double *W = Weights.Values();
  double **W_Pointers = Weights.Pointers();
  
  for (int i=0; i < NumVectors_; i++) 
    {
      if (!OneW) W = W_Pointers[i]; // If Weights has the same number of vectors as this, use each weight vector
      double sum = 0.0;
      const double * const from = Pointers_[i];
#ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
#pragma omp parallel for reduction (+:sum) default(none) shared(W)
#endif
      for (int j=0; j < myLength; j++) {
        double tmp = from[j]/W[j];
        sum += tmp * tmp;
      }
      DoubleTemp_[i] = sum;
    }
  double OneOverN;
  if (DistributedGlobal()) {
    Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
    OneOverN = 1.0 / (double) GlobalLength_;
  }
  else {
    for (int i=0; i< NumVectors_; ++i) Result[i] = DoubleTemp_[i];
    OneOverN = 1.0 / (double) myLength;
  }
  for (int i=0; i < NumVectors_; i++) Result[i] = std::sqrt(Result[i]*OneOverN);
  
  UpdateFlops(3*GlobalLength_*NumVectors_);
  
  return(0);
  }

//=========================================================================
int  Epetra_MultiVector::MinValue (double* Result) const {
  
  // Minimum value of each vector in MultiVector 
  
  int ierr = 0;

  const int myLength = MyLength_;
  UpdateDoubleTemp();

  for (int i=0; i < NumVectors_; i++) 
    {
      const double * const from = Pointers_[i];
      double MinVal = Epetra_MaxDouble;
      if (myLength>0) MinVal = from[0];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel default(none) shared(MinVal)
{
      double localMinVal = MinVal;
#pragma omp for
      for (int j=0; j< myLength; j++) localMinVal = EPETRA_MIN(localMinVal,from[j]); 
#pragma omp critical 
      {
      MinVal = EPETRA_MIN(MinVal,localMinVal);
      }
}
      DoubleTemp_[i] = MinVal;
#else
      for (int j=0; j< myLength; j++) MinVal = EPETRA_MIN(MinVal,from[j]); 
      DoubleTemp_[i] = MinVal;
#endif
    }

  if (myLength > 0) {
    for(int i=0; i<NumVectors_; ++i) {
      Result[i] = DoubleTemp_[i];
    }
  }

  //If myLength == 0 and Comm_->NumProc() == 1, then Result has
  //not been referenced. Also, if vector contents are uninitialized
  //then Result contents are not well defined...

  if (Comm_->NumProc() == 1 || !DistributedGlobal()) return(ierr);

  //We're going to use MPI_Allgather to gather every proc's local-
  //min values onto every other proc. We'll use the last position
  //of the DoubleTemp_ array to indicate whether this proc has
  //valid data that should be considered by other procs when forming
  //the global-min results.

  if (myLength == 0) DoubleTemp_[NumVectors_] = 0.0;
  else DoubleTemp_[NumVectors_] = 1.0;

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
  double* dwork = new double[numProcs*(NumVectors_+1)];

  MPI_Allgather(DoubleTemp_, NumVectors_+1, MPI_DOUBLE,
                dwork, NumVectors_+1, MPI_DOUBLE, mpicomm);

  //if myLength==0, then our Result array currently contains
  //Epetra_MaxDouble from the local-min calculations above. In this
  //case we'll overwrite our Result array with values from the first
  //processor that sent valid data.
  bool overwrite = myLength == 0 ? true : false;

  int myPID = epetrampicomm->MyPID();
  double* dwork_ptr = dwork;

  for(int j=0; j<numProcs; ++j) {

    //skip data from self, and skip data from
    //procs with DoubleTemp_[NumVectors_] == 0.0.
    if (j == myPID || dwork_ptr[NumVectors_] < 0.5) {
      dwork_ptr += NumVectors_+1;
      continue;
    }

    for(int i=0; i<NumVectors_; ++i) {
      double val = dwork_ptr[i];

      //Set val into our Result array if overwrite is true (see above),
      //or if val is less than our current Result[i].
      if (overwrite || (Result[i] > val)) Result[i] = val;
    }

    //Now set overwrite to false so that we'll do the right thing
    //when processing data from subsequent processors.
    if (overwrite) overwrite = false;

    dwork_ptr += NumVectors_+1;
  }

  delete [] dwork;
#endif

  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine
  
  return(ierr);
}

//=========================================================================
int  Epetra_MultiVector::MaxValue (double* Result) const {
  
  // Maximum value of each vector in MultiVector 
  
  int ierr = 0;

  const int myLength = MyLength_;
  UpdateDoubleTemp();

  for (int i=0; i < NumVectors_; i++) 
    {
      const double * const from = Pointers_[i];
      double MaxVal = -Epetra_MaxDouble;
      if (myLength>0) MaxVal = from[0];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel default(none) shared(MaxVal)
{
      double localMaxVal = MaxVal;
#pragma omp for
      for (int j=0; j< myLength; j++) localMaxVal = EPETRA_MAX(localMaxVal,from[j]); 
#pragma omp critical 
      {
      MaxVal = EPETRA_MAX(MaxVal,localMaxVal);
      }
}
      DoubleTemp_[i] = MaxVal;
#else
      for (int j=0; j< myLength; j++) MaxVal = EPETRA_MAX(MaxVal,from[j]); 
      DoubleTemp_[i] = MaxVal;
#endif
    }

  if (myLength > 0) {
    for(int i=0; i<NumVectors_; ++i) {
      Result[i] = DoubleTemp_[i];
    }
  }

  //If myLength == 0 and Comm_->NumProc() == 1, then Result has
  //not been referenced. Also, if vector contents are uninitialized
  //then Result contents are not well defined...

  if (Comm_->NumProc() == 1  || !DistributedGlobal()) return(ierr);

  //We're going to use MPI_Allgather to gather every proc's local-
  //max values onto every other proc. We'll use the last position
  //of the DoubleTemp_ array to indicate whether this proc has
  //valid data that should be considered by other procs when forming
  //the global-max results.

  if (myLength == 0) DoubleTemp_[NumVectors_] = 0.0;
  else DoubleTemp_[NumVectors_] = 1.0;

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
  double* dwork = new double[numProcs*(NumVectors_+1)];

  MPI_Allgather(DoubleTemp_, NumVectors_+1, MPI_DOUBLE,
                dwork, NumVectors_+1, MPI_DOUBLE, mpicomm);

  //if myLength==0, then our Result array currently contains
  //-Epetra_MaxDouble from the local-max calculations above. In this
  //case we'll overwrite our Result array with values from the first
  //processor that sent valid data.
  bool overwrite = myLength == 0 ? true : false;

  int myPID = epetrampicomm->MyPID();
  double* dwork_ptr = dwork;

  for(int j=0; j<numProcs; ++j) {

    //skip data from self, and skip data from
    //procs with DoubleTemp_[NumVectors_] == 0.0.
    if (j == myPID || dwork_ptr[NumVectors_] < 0.5) {
      dwork_ptr += NumVectors_+1;
      continue;
    }

    for(int i=0; i<NumVectors_; ++i) {
      double val = dwork_ptr[i];

      //Set val into our Result array if overwrite is true (see above),
      //or if val is larger than our current Result[i].
      if (overwrite || (Result[i] < val)) Result[i] = val; 
    }

    //Now set overwrite to false so that we'll do the right thing
    //when processing data from subsequent processors.
    if (overwrite) overwrite = false;

    dwork_ptr += NumVectors_+1;
  }

  delete [] dwork;
#endif

  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine
  
  return(ierr);
}

//=========================================================================
int  Epetra_MultiVector::MeanValue (double* Result) const {
  
  // Mean value of each vector in MultiVector 
  
  const int myLength = MyLength_;

  if (!Map().UniqueGIDs()) {EPETRA_CHK_ERR(-1);}

  double fGlobalLength = 1.0/EPETRA_MAX((double) GlobalLength_, 1.0);
  

  UpdateDoubleTemp();

  for (int i=0; i < NumVectors_; i++) 
    {
      double sum = 0.0;
      const double * const from = Pointers_[i];
#ifdef EPETRA_HAVE_OMP_NONASSOCIATIVE
#pragma omp parallel for reduction (+:sum) default(none)
#endif
      for (int j=0; j < myLength; j++) sum += from[j];
      DoubleTemp_[i] = sum;
    }
  if (DistributedGlobal())
    Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  else
    for (int i=0; i< NumVectors_; ++i) Result[i] = DoubleTemp_[i];

  for (int i=0; i < NumVectors_; i++) Result[i] = Result[i]*fGlobalLength;
  
  UpdateFlops(GlobalLength_*NumVectors_);

  return(0);
}

  //=========================================================================
  int  Epetra_MultiVector::Multiply (char TransA, char TransB, double ScalarAB, 
					const Epetra_MultiVector& A, 
					const Epetra_MultiVector& B,
					double ScalarThis ) {

    // This routine performs a variety of matrix-matrix multiply operations, interpreting
    // the Epetra_MultiVector (this-aka C , A and B) as 2D matrices.  Variations are due to
    // the fact that A, B and C can be local replicated or global distributed
    // Epetra_MultiVectors and that we may or may not operate with the transpose of 
    // A and B.  Possible cases are:

    //                                       Num
    //     OPERATIONS                        case  Notes
    // 1) C(local) = A^X(local) * B^X(local)  4   (X=Trans or Not, No comm needed) 
    // 2) C(local) = A^T(distr) * B  (distr)  1   (2D dot product, replicate C)
    // 3) C(distr) = A  (distr) * B^X(local)  2   (2D vector update, no comm needed)

    // Note that the following operations are not meaningful for 
    // 1D distributions:

    // 1) C(local) = A^T(distr) * B^T(distr)  1
    // 2) C(local) = A  (distr) * B^X(distr)  2
    // 3) C(distr) = A^X(local) * B^X(local)  4
    // 4) C(distr) = A^X(local) * B^X(distr)  4
    // 5) C(distr) = A^T(distr) * B^X(local)  2
    // 6) C(local) = A^X(distr) * B^X(local)  4
    // 7) C(distr) = A^X(distr) * B^X(local)  4
    // 8) C(local) = A^X(local) * B^X(distr)  4

    // Total of 32 case (2^5).

    
    //if (!ConstantStride_    ||
    //!A.ConstantStride() ||
    //!B.ConstantStride()    ) EPETRA_CHK_ERR(-1); // Return error

    // Check for compatible dimensions

    int A_nrows = (TransA=='T') ? A.NumVectors() : A.MyLength();
    int A_ncols = (TransA=='T') ? A.MyLength() : A.NumVectors();
    int B_nrows = (TransB=='T') ? B.NumVectors() : B.MyLength();
    int B_ncols = (TransB=='T') ? B.MyLength() : B.NumVectors();

    double Scalar_local = ScalarThis; // local copy of Scalar
  const int myLength = MyLength_;

    if( myLength      != A_nrows     ||   // RAB: 2002/01/25: Minor reformat to allow
		A_ncols        != B_nrows     ||   //   setting breakpoint at error return.
		NumVectors_    != B_ncols  )
		EPETRA_CHK_ERR(-2); // Return error

    bool A_is_local = (!A.DistributedGlobal());
    bool B_is_local = (!B.DistributedGlobal());
    bool C_is_local = (!DistributedGlobal_);
    bool Case1 = ( A_is_local &&  B_is_local &&  C_is_local);  // Case 1 above
    bool Case2 = (!A_is_local && !B_is_local &&  C_is_local && TransA=='T' );// Case 2
    bool Case3 = (!A_is_local &&  B_is_local && !C_is_local && TransA=='N');// Case 3

    if (Case2 && (!A.Map().UniqueGIDs() || !B.Map().UniqueGIDs())) {EPETRA_CHK_ERR(-4);}
    if (Case3 && (!A.Map().UniqueGIDs() || !Map().UniqueGIDs())) {EPETRA_CHK_ERR(-5);}
  
    // Test for meaningful cases

    if (Case1 || Case2 || Case3)
      {
	if (ScalarThis!=0.0 && Case2)
	  {
	    const int MyPID = Comm_->MyPID();
	    if (MyPID!=0) Scalar_local = 0.0;
	  }

        // Check if A, B, C have constant stride, if not then make temp copy (strided)

        Epetra_MultiVector * A_tmp, * B_tmp, *C_tmp;
        if (!ConstantStride_) C_tmp = new Epetra_MultiVector(*this);
        else C_tmp = this;
          
        if (!A.ConstantStride()) A_tmp = new Epetra_MultiVector(A);
        else A_tmp = (Epetra_MultiVector *) &A;
    
        if (!B.ConstantStride()) B_tmp = new Epetra_MultiVector(B);
        else B_tmp = (Epetra_MultiVector *) &B;
    	
    
	int m = myLength;
	int n = NumVectors_;
	int k = A_ncols;
	int lda = EPETRA_MAX(A_tmp->Stride(),1); // The reference BLAS implementation requires lda, ldb, ldc > 0 even if m, n, or k = 0
	int ldb = EPETRA_MAX(B_tmp->Stride(),1);
	int ldc = EPETRA_MAX(C_tmp->Stride(),1);
	double *Ap = A_tmp->Values();
	double *Bp = B_tmp->Values();
	double *Cp = C_tmp->Values();
   
	GEMM(TransA, TransB,  m, n, k, ScalarAB,
	     Ap, lda, Bp, ldb, Scalar_local, Cp, ldc);
      
	// FLOP Counts
	//                                       Num
	//     OPERATIONS                        case  Notes 
	// 1) C(local) = A^X(local) * B^X(local)  4   (X=Trans or Not, No comm needed)
	// 2) C(local) = A^T(distr) * B  (distr)  1   (2D dot product, replicate C)      
	// 3) C(distr) = A  (distr) * B^X(local)  2   (2D vector update, no comm needed)

	// For Case 1 we only count the local operations, since we are interested in serial
	// cost.  Computation on other processors is redundant.
	if (Case1)
	  {	    
	    UpdateFlops(2*m*n*k);
	    if (ScalarAB!=1.0) UpdateFlops(m*n);
	    if (ScalarThis==1.0) UpdateFlops(m*n); else if (ScalarThis!=0.0) UpdateFlops(2*m*n);
	  }
	else if (Case2)
	  {	    
	    UpdateFlops(2*m*n*A.GlobalLength());
	    if (ScalarAB!=1.0) UpdateFlops(m*n);
	    if (ScalarThis==1.0) UpdateFlops(m*n); else if (ScalarThis!=0.0) UpdateFlops(2*m*n);
	  }
	else
	  {
	    UpdateFlops(2*GlobalLength_*n*k);
	    if (ScalarAB!=1.0) UpdateFlops(GlobalLength_*n);
	    if (ScalarThis==1.0) UpdateFlops(GlobalLength_*n);
	    else if (ScalarThis!=0.0) UpdateFlops(2*GlobalLength_*n);
	  }

	// If A or B were not strided, dispose of extra copies.
	if (!A.ConstantStride()) delete A_tmp;
	if (!B.ConstantStride()) delete B_tmp;

	// If C was not strided, copy from strided version and delete
	if (!ConstantStride_) 
	  {
	    C_tmp->ExtractCopy(Pointers_);
	    delete C_tmp;
	  }

	// If Case 2 then sum up C and distribute it to all processors.

	if (Case2) {EPETRA_CHK_ERR(Reduce());}

	return(0);

      }
    else {EPETRA_CHK_ERR(-3)}; // Return error: not a supported operation

  return(0);
  }


//=========================================================================
int Epetra_MultiVector::Multiply(double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B,
		       double ScalarThis) {
  
  
  // Hadamard product of two MultiVectors: 
  // this = ScalarThis * this + ScalarAB * A * B (element-wise)
  
  if (ScalarAB==0.0) {
    EPETRA_CHK_ERR(Scale(ScalarThis));
    return(0);
  }
  const int myLength = MyLength_;
			   
  if (A.NumVectors() != 1 && A.NumVectors() != B.NumVectors()) EPETRA_CHK_ERR(-1); // A must have one column or be the same as B.
  if (NumVectors_ != B.NumVectors()) EPETRA_CHK_ERR(-2);
  if (myLength != A.MyLength() || myLength != B.MyLength()) EPETRA_CHK_ERR(-3);
  
  double **A_Pointers = (double**)A.Pointers();
  double **B_Pointers = (double**)B.Pointers();

  int IncA = 1;
  if (A.NumVectors() == 1 ) IncA = 0;

    if (ScalarThis==0.0) {
      if (ScalarAB==1.0)
	{
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] =  Aptr[j] * Bptr[j];
            }
	  }
	  UpdateFlops(GlobalLength_*NumVectors_);
	}
      else
	{
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarAB)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] = ScalarAB * Aptr[j] * Bptr[j];
            }
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
    }
    else if (ScalarThis==1.0) {
      if (ScalarAB==1.0)
	{
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] +=  Aptr[j] * Bptr[j];
            }
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
      else {
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarAB)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] += ScalarAB * Aptr[j] * Bptr[j];
            }
          }
          UpdateFlops(3*GlobalLength_*NumVectors_);
      }
    }
    else { // if (ScalarThis!=1.0 && ScalarThis !=0 ) 
      if (ScalarAB==1.0)
	{
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarThis)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] =  ScalarThis * to[j] +
                  Aptr[j] * Bptr[j];
            }
	  }
	  UpdateFlops(3*GlobalLength_*NumVectors_);
	}
      else
	{
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarThis,ScalarAB)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] = ScalarThis * to[j] +
                     ScalarAB * Aptr[j] * Bptr[j];
            }
	  }
	  UpdateFlops(4*GlobalLength_*NumVectors_);
	}
    }
  return(0);
}
//=========================================================================
int Epetra_MultiVector::ReciprocalMultiply(double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B,
		       double ScalarThis) {
  
  
  // Hadamard product of two MultiVectors: 
  // this = ScalarThis * this + ScalarAB * B / A (element-wise)
  
  if (ScalarAB==0.0) {
    EPETRA_CHK_ERR(Scale(ScalarThis));
    return(0);
  }
  const int myLength = MyLength_;
			   
  if (A.NumVectors() != 1 && A.NumVectors() != B.NumVectors()) EPETRA_CHK_ERR(-1); // A must have one column or be the same as B.
  if (NumVectors_ != B.NumVectors()) EPETRA_CHK_ERR(-2);
  if (myLength != A.MyLength() || myLength != B.MyLength()) EPETRA_CHK_ERR(-3);
  
  double **A_Pointers = (double**)A.Pointers();
  double **B_Pointers = (double**)B.Pointers();

  int IncA = 1;
  if (A.NumVectors() == 1 ) IncA = 0;

    if (ScalarThis==0.0) {
      if (ScalarAB==1.0)
	{
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] = Bptr[j] / Aptr[j];
            }
	  }
	  UpdateFlops(GlobalLength_*NumVectors_);
	}
      else
	{
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarAB)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] = ScalarAB * Bptr[j] / Aptr[j];
            }
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
    }
    else if (ScalarThis==1.0) {
      if (ScalarAB==1.0)
	{
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] +=  Bptr[j] / Aptr[j];
            }
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
      else
	{
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarAB)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] += ScalarAB * Bptr[j] / Aptr[j];
            }
	  }
	  UpdateFlops(3*GlobalLength_*NumVectors_);
        }
    }
    else { // if (ScalarThis!=1.0 && ScalarThis !=0 ) 
      if (ScalarAB==1.0)
	{
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarThis)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] =  ScalarThis * to[j] +
                   Bptr[j] / Aptr[j];
            }
	  }
	  UpdateFlops(3*GlobalLength_*NumVectors_);
	}
      else {
	  for (int i = 0; i < NumVectors_; i++) {
	    const double * const Aptr = A_Pointers[i*IncA];
            const double * const Bptr = B_Pointers[i];
            double * const to = Pointers_[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(ScalarAB,ScalarThis)
#endif
	    for (int j = 0; j < myLength; j++) {
              to[j] = ScalarThis * to[j] + ScalarAB * 
					      Bptr[j] / Aptr[j];
            }
	  }
          UpdateFlops(4*GlobalLength_*NumVectors_);
      }
    }
  return(0);
}


//=======================================================================
Epetra_Vector *& Epetra_MultiVector::operator () (int index)  {
  
  //  Epetra_MultiVector::operator () --- return non-const reference 
  
  if (index < 0 || index >=NumVectors_) 
    throw ReportError("Vector index = " + toString(index) + "is out of range. Number of Vectors = " + toString(NumVectors_), -1);

  UpdateVectors();

  // Create a new Epetra_Vector that is a view of ith vector, if not already present
  if (Vectors_[index]==0)
    Vectors_[index] = new Epetra_Vector(View, Map(), Pointers_[index]);
  return(Vectors_[index]);
}

//=======================================================================
const Epetra_Vector *& Epetra_MultiVector::operator () (int index) const {
  
  //  Epetra_MultiVector::operator () --- return non-const reference 

  if (index < 0 || index >=NumVectors_) 
    throw ReportError("Vector index = " + toString(index) + "is out of range. Number of Vectors = " + toString(NumVectors_), -1);

  UpdateVectors();

  if (Vectors_[index]==0)
    Vectors_[index] = new Epetra_Vector(View, Map(), Pointers_[index]);

  const Epetra_Vector * & temp = (const Epetra_Vector * &) (Vectors_[index]);
  return(temp);
}

//========================================================================
Epetra_MultiVector& Epetra_MultiVector::operator = (const Epetra_MultiVector& Source) {
  
  // Check for special case of this=Source
  if (this != &Source) Assign(Source);
  
  return(*this);
}

//=========================================================================
void Epetra_MultiVector::Assign(const Epetra_MultiVector& A) {
  
  const int myLength = MyLength_;
  if (NumVectors_ != A.NumVectors())
    throw ReportError("Number of vectors incompatible in Assign.  The this MultiVector has NumVectors = " + toString(NumVectors_)
		      + ".  The A MultiVector has NumVectors = " + toString(A.NumVectors()), -3);
  if (myLength != A.MyLength())
    throw ReportError("Length of MultiVectors incompatible in Assign.  The this MultiVector has MyLength = " + toString(myLength)
		      + ".  The A MultiVector has MyLength = " + toString(A.MyLength()), -4);
  
  double ** A_Pointers = A.Pointers();
  for (int i = 0; i< NumVectors_; i++) {
      double * const to = Pointers_[i];
      const double * const from = A_Pointers[i];
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none)
#endif
      for (int j=0; j<myLength; j++) to[j] = from[j];
    }
    return;    
  }

  //=========================================================================
  int  Epetra_MultiVector::Reduce() {

    // Global reduction on each entry of a Replicated Local MultiVector
  const int myLength = MyLength_;
    double * source = 0;
      if (myLength>0) source = new double[myLength*NumVectors_];
    double * target = 0;
    bool packed = (ConstantStride_ && (Stride_==myLength));
    if (packed) {
       for (int i=0; i<myLength*NumVectors_; i++) source[i] = Values_[i];
       target = Values_;
    }
    else {
       double * tmp1 = source;
       for (int i = 0; i < NumVectors_; i++) {
         double * tmp2 = Pointers_[i];
         for (int j=0; j< myLength; j++) *tmp1++ = *tmp2++;
       }
       if (myLength>0) target = new double[myLength*NumVectors_];
    }

    Comm_->SumAll(source, target, myLength*NumVectors_);
    if (myLength>0) delete [] source;
    if (!packed) {
       double * tmp2 = target;
       for (int i = 0; i < NumVectors_; i++) {
         double * tmp1 = Pointers_[i];
         for (int j=0; j< myLength; j++) *tmp1++ = *tmp2++;
       }
       if (myLength>0) delete [] target;
    }
    // UpdateFlops(0);  No serial Flops in this function
    return(0);
  }
//=======================================================================
int Epetra_MultiVector::ResetView(double ** ArrayOfPointers) {
	
	if (!UserAllocated_) {
		EPETRA_CHK_ERR(-1); // Can't reset view if multivector was not allocated as a view
	}

  for (int i = 0; i< NumVectors_; i++) Pointers_[i] = ArrayOfPointers[i];
  DoView();
  
	return(0);
  }
//=======================================================================
void Epetra_MultiVector::Print(ostream& os) const {
  int MyPID = Map().Comm().MyPID();
  int NumProc = Map().Comm().NumProc();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumVectors1 = NumVectors();
      int NumMyElements1 =Map(). NumMyElements();
      int MaxElementSize1 = Map().MaxElementSize();
      int * FirstPointInElementList1 = NULL;
      if (MaxElementSize1!=1) FirstPointInElementList1 = Map().FirstPointInElementList();
      double ** A_Pointers = Pointers();

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
	os << endl;
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
			int * MyGlobalElements1 = Map().MyGlobalElements();
			os << MyGlobalElements1[i] << "    ";
		  }
		  else if(Map().GlobalIndicesLongLong())
		  {
			long long * MyGlobalElements1 = Map().MyGlobalElements_LL();
			os << MyGlobalElements1[i] << "    ";
		  }
		  else
		    throw ReportError("Epetra_MultiVector::Print ERROR, Don't know map global index type.",-1);

       iii = i;
       }
	  else {
		  if(Map().GlobalIndicesInt())
		  {
			int * MyGlobalElements1 = Map().MyGlobalElements();
			os <<  MyGlobalElements1[i]<< "/" << ii << "    ";
		  }
		  else if(Map().GlobalIndicesLongLong())
		  {
			long long * MyGlobalElements1 = Map().MyGlobalElements_LL();
			os <<  MyGlobalElements1[i]<< "/" << ii << "    ";
		  }
		  else
		    throw ReportError("Epetra_MultiVector::Print ERROR, Don't know map global index type.",-1);

         iii = FirstPointInElementList1[i]+ii;
       }
	  for (int j = 0; j < NumVectors1 ; j++)
	    {   
	      os.width(20);
	      os <<  A_Pointers[j][iii];
	    }
	  os << endl;
	}
      }
      os << flush; 
    }

    // Do a few global ops to give I/O a chance to complete
    Map().Comm().Barrier();
    Map().Comm().Barrier();
    Map().Comm().Barrier();
  }
  return;
}
