
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

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Distributor.h"

//=============================================================================

// Epetra_BlockMap Constructor

Epetra_MultiVector::Epetra_MultiVector(const Epetra_BlockMap& Map, int NumVectors, bool zeroOut)
  : Epetra_DistObject(Map, "Epetra::MultiVector"),
    Epetra_CompObject(),
    IndexBase_(Map.IndexBase()),
    Values_(0),
    Pointers_(0),
    MyLength_(Map.NumMyPoints()),
    GlobalLength_(Map.NumGlobalPoints()),
    NumVectors_(NumVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(Map.NumMyPoints()),
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
    IndexBase_(Source.IndexBase_),
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

Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, 
					     double *A, int MyLDA, int NumVectors)
  : Epetra_DistObject(Map, "Epetra::MultiVector"),
    Epetra_CompObject(),
    IndexBase_(Map.IndexBase()),
    Values_(0),
    Pointers_(0),
    MyLength_(Map.NumMyPoints()),
    GlobalLength_(Map.NumGlobalPoints()),
    NumVectors_(NumVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(Map.NumMyPoints()),
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

Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, 
					      double **ArrayOfPointers, int NumVectors)
  : Epetra_DistObject(Map, "Epetra::MultiVector"),
    Epetra_CompObject(),
    IndexBase_(Map.IndexBase()),
    Values_(0),
    Pointers_(0),
    MyLength_(Map.NumMyPoints()),
    GlobalLength_(Map.NumGlobalPoints()),
    NumVectors_(NumVectors),
    UserAllocated_(false),
    ConstantStride_(true),
    Stride_(Map.NumMyPoints()),
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
					     int *Indices, int NumVectors)
  : Epetra_DistObject(Source.Map(), "Epetra::MultiVector"),
    Epetra_CompObject(),
    IndexBase_(Source.IndexBase_),
    Values_(0),
    Pointers_(0),
    MyLength_(Source.MyLength_),
    GlobalLength_(Source.GlobalLength_),
    NumVectors_(NumVectors),
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
					     int StartIndex, int NumVectors)
  : Epetra_DistObject(Source.Map(), "Epetra::MultiVector"),
    Epetra_CompObject(),
    IndexBase_(Source.IndexBase_),
    Values_(0),
    Pointers_(0),
    MyLength_(Source.MyLength_),
    GlobalLength_(Source.GlobalLength_),
    NumVectors_(NumVectors),
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

  for (int i=0; i<NumVectors_; i++) if (Vectors_[i]!=0) delete Vectors_[i];

  delete [] Vectors_;

  if (!UserAllocated_ && Values_!=0) delete [] Values_;

  delete [] Pointers_;

  delete [] DoubleTemp_;
  delete [] IntTemp_;

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

  DoubleTemp_ = new double[NumVectors_];
  IntTemp_ = new int[NumVectors_];

  Vectors_ = new Epetra_Vector *[NumVectors_];
  for (int i=0; i<NumVectors_; i++) Vectors_[i] = 0;
  
  int randval = rand(); // Use POSIX standard random function
  if (DistributedGlobal_)
    Util_.SetSeed(2*Comm_->MyPID() + randval);
  else
    Util_.SetSeed(randval); // Replicated Local MultiVectors get the same seed on all PEs

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

  // \interal { Optimization of this function can impact performance since it 
  //           involves a fair amount of memory traffic.}

  // On exit, Pointers_ is redefined to point to its own MultiVector vectors.

  for (int i = 0; i< NumVectors_; i++)
    {
      double * from = Pointers_[i];
      double * to = Values_+i*Stride_;
      Pointers_[i] = to;
      for (int j=0; j<MyLength_; j++) *to++ = *from++;
    }

  return(0);
}
//=========================================================================
int Epetra_MultiVector::AllocateForView(void)
{
  
  if (NumVectors_<=0) 
    throw ReportError("Number of Vectors = " + toString(NumVectors_) + ", but must be greater than zero", -1);
 
  Pointers_ = new double *[NumVectors_];
  
  DoubleTemp_ = new double[NumVectors_];
  IntTemp_ = new int[NumVectors_];

  Vectors_ = new Epetra_Vector *[NumVectors_];
  for (int i=0; i<NumVectors_; i++) Vectors_[i] = 0;
  
  int randval = rand(); // Use POSIX standard random function
  if (DistributedGlobal_)
    Util_.SetSeed(2*Comm_->MyPID() + randval);
  else
    Util_.SetSeed(randval); // Replicated Local MultiVectors get the same seed on all PEs

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

  Stride_ = Pointers_[1] - Pointers_[0];
  ConstantStride_ = false;

  for (int i = 1; i < NumVectors_-1; i++) if (Pointers_[i+1] - Pointers_[i] != Stride_) return(0);

  ConstantStride_ = true;

  return(0);
}
//=========================================================================
int Epetra_MultiVector::ReplaceGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue) {

 // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue(GlobalRow, 0, VectorIndex, ScalarValue, false));
  return(0);
}
//=========================================================================
int Epetra_MultiVector::ReplaceGlobalValue(int GlobalBlockRow, int BlockRowOffset, 
					   int VectorIndex, double ScalarValue) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeGlobalValue(GlobalBlockRow, BlockRowOffset, VectorIndex, ScalarValue, false)); 
  return(0);
}
//=========================================================================
int Epetra_MultiVector::SumIntoGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue) {

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
int Epetra_MultiVector::ChangeGlobalValue(int GlobalBlockRow, int BlockRowOffset, 
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
  int i,j;
  const double a = 16807.0, BigInt=2147483647.0, DbleOne=1.0, DbleTwo=2.0;
  
  for(i=0; i < NumVectors_; i++)
    for (j=0; j<MyLength_; j++){
      Seed_ = fmod( a*Seed_, BigInt );
      Pointers_[i][j] = DbleTwo*(Seed_/BigInt)-DbleOne;
    }
	*/

	for(int i = 0; i < NumVectors_; i++)
		for(int j = 0; j < MyLength_; j++)
			Pointers_[i][j] = Util_.RandomDouble();

  return(0);
}
 
//=========================================================================

// Extract a copy of a Epetra_MultiVector.  Put in a user's Fortran-style array

int Epetra_MultiVector::ExtractCopy(double *A, int MyLDA) const {
  if (NumVectors_>1 && Stride_ > MyLDA) EPETRA_CHK_ERR(-1); // LDA not big enough
  
  for (int i=0; i< NumVectors_; i++)
    {
      double * from = Pointers_[i];
      double * to = A + i*MyLDA;
      for (int j=0; j<MyLength_; j++) *to++ = *from++;
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
  for (int i=0; i< NumVectors_; i++)
    {
      double * from = Pointers_[i];
      double * to = ArrayOfPointers[i];
      for (int j=0; j<MyLength_; j++) *to++ = *from++;
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

  for (int i = 0; i < NumVectors_; i++)
    for (int j=0; j<MyLength_; j++) Pointers_[i][j] = ScalarConstant;
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
                                       const Epetra_OffsetIndex * Indexor) {

  const Epetra_MultiVector & A = dynamic_cast<const Epetra_MultiVector &>(Source);

  double **From = A.Pointers();
  double **To = Pointers_;
  int NumVectors = NumVectors_;

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
  int i, j, jj, jjj, k;
  
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
      for (i=0; i < NumVectors; i++)
	for (j=0; j<NumSameEntries; j++)
	  To[i][j] = From[i][j];
    }
  // Do local permutation next
  if (NumPermuteIDs>0) {
  
    // Point entry case
    if (Case1) {
      
      if (NumVectors==1)
	for (j=0; j<NumPermuteIDs; j++) 
	  To[0][PermuteToLIDs[j]] = From[0][PermuteFromLIDs[j]];
      
      else {
	for (j=0; j<NumPermuteIDs; j++) {
	  jj = PermuteToLIDs[j];
	  jjj = PermuteFromLIDs[j];
	  for (i=0; i<NumVectors; i++)
	    To[i][jj] = From[i][jjj];
	}
      }
    }
    // constant element size case
    else if (Case2) {
      
      for (j=0; j<NumPermuteIDs; j++) {
	jj = MaxElementSize*PermuteToLIDs[j];
	jjj = MaxElementSize*PermuteFromLIDs[j];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<MaxElementSize; k++)
	    To[i][jj+k] = From[i][jjj+k];
      }
    }
    
    // variable element size case
    else {
      
      for (j=0; j<NumPermuteIDs; j++) {
	jj = ToFirstPointInElementList[PermuteToLIDs[j]];
	jjj = FromFirstPointInElementList[PermuteFromLIDs[j]];
	int ElementSize = FromElementSizeList[PermuteFromLIDs[j]];
	for (i=0; i<NumVectors; i++)
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
                                       Epetra_Distributor & Distor) {



  const Epetra_MultiVector & A = dynamic_cast<const Epetra_MultiVector &>(Source);
  int i, j, jj, k;

  double **From = A.Pointers();
  int MaxElementSize = Map().MaxElementSize();
  int NumVectors = NumVectors_;
  bool ConstantElementSize = Map().ConstantElementSize();

  int * FromFirstPointInElementList = 0;
  int * FromElementSizeList = 0;

  if (!ConstantElementSize) {
    FromFirstPointInElementList = A.Map().FirstPointInElementList();
    FromElementSizeList = A.Map().ElementSizeList();
  }

  double * DoubleExports = 0;
  double * DoubleImports = 0;

  SizeOfPacket = NumVectors*MaxElementSize;

  if (SizeOfPacket*NumExportIDs>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = SizeOfPacket*NumExportIDs;
    DoubleExports = new double[LenExports];
    Exports = (char *) DoubleExports;
  }

  SizeOfPacket *= sizeof(double);

  double * ptr;

  if (NumExportIDs>0) {
    ptr = (double *) Exports;
    
    
    // Point entry case
    if (MaxElementSize==1) {
      
      if (NumVectors==1) for (j=0; j<NumExportIDs; j++) *ptr++ = From[0][ExportLIDs[j]];

      else {
	for (j=0; j<NumExportIDs; j++) {
	  jj = ExportLIDs[j];
	  for (i=0; i<NumVectors; i++)
	    *ptr++ = From[i][jj];
	}
      }
    }

    // constant element size case
    else if (ConstantElementSize) {
      
      for (j=0; j<NumExportIDs; j++) {
	jj = MaxElementSize*ExportLIDs[j];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<MaxElementSize; k++)
	    *ptr++ = From[i][jj+k];
      }
    }
    
    // variable element size case
    else {
      
      int SizeOfPacket = NumVectors*MaxElementSize;
      for (j=0; j<NumExportIDs; j++) {
	ptr = (double *) Exports + j*SizeOfPacket;
	jj = FromFirstPointInElementList[ExportLIDs[j]];
	int ElementSize = FromElementSizeList[ExportLIDs[j]];
	for (i=0; i<NumVectors; i++)
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
                                         const Epetra_OffsetIndex * Indexor ) {
  int i, j, jj, k;
  
  if(    CombineMode != Add
      && CombineMode != Zero
      && CombineMode != Insert
      && CombineMode != Average
      && CombineMode != AbsMax )
    EPETRA_CHK_ERR(-1); //Unsupported CombinedMode, will default to Zero

  if (NumImportIDs<=0) return(0);

  double ** To = Pointers_;
  int NumVectors = NumVectors_;
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
      
    if (NumVectors==1) {
      if (CombineMode==Add)
	for (j=0; j<NumImportIDs; j++) To[0][ImportLIDs[j]] += *ptr++; // Add to existing value
      else if(CombineMode==Insert)
	for (j=0; j<NumImportIDs; j++) To[0][ImportLIDs[j]] = *ptr++;
      else if(CombineMode==AbsMax)
        for (j=0; j<NumImportIDs; j++) To[0][ImportLIDs[j]] = EPETRA_MAX( To[0][ImportLIDs[j]],fabs(*ptr++)); //
      // Note:  The following form of averaging is not a true average if more that one value is combined.
      //        This might be an issue in the future, but we leave this way for now.
      else if(CombineMode==Average)
	for (j=0; j<NumImportIDs; j++) {To[0][ImportLIDs[j]] += *ptr++; To[0][ImportLIDs[j]] *= 0.5;}
    }

    else {  // NumVectors>1

      if (CombineMode==Add) {
	for (j=0; j<NumImportIDs; j++) {
	  jj = ImportLIDs[j];
	  for (i=0; i<NumVectors; i++)
	    To[i][jj] += *ptr++; // Add to existing value
	}
      }
      else if(CombineMode==Insert) {
	for (j=0; j<NumImportIDs; j++) {
	  jj = ImportLIDs[j];
	  for (i=0; i<NumVectors; i++)
	    To[i][jj] = *ptr++;
	}
      }
      else if(CombineMode==AbsMax) {
        for (j=0; j<NumImportIDs; j++) {
          jj = ImportLIDs[j];
          for (i=0; i<NumVectors; i++)
            To[i][jj] = EPETRA_MAX( To[i][jj], fabs(*ptr++) );
        }
      }
      // Note:  The following form of averaging is not a true average if more that one value is combined.
      //        This might be an issue in the future, but we leave this way for now.
      else if(CombineMode==Average) {
	for (j=0; j<NumImportIDs; j++) {
	  jj = ImportLIDs[j];
	  for (i=0; i<NumVectors; i++)
	    { To[i][jj] += *ptr++;  To[i][jj] *= 0.5;}
	}
      }
    }
  }

  // constant element size case

  else if (ConstantElementSize) {
   
    if (CombineMode==Add) {
      for (j=0; j<NumImportIDs; j++) {
	jj = MaxElementSize*ImportLIDs[j];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<MaxElementSize; k++)
	    To[i][jj+k] += *ptr++; // Add to existing value
      }
    }
    else if(CombineMode==Insert) {
      for (j=0; j<NumImportIDs; j++) {
	jj = MaxElementSize*ImportLIDs[j];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<MaxElementSize; k++)
	    To[i][jj+k] = *ptr++;
      }
    }
    else if(CombineMode==AbsMax) {
      for (j=0; j<NumImportIDs; j++) {
	jj = MaxElementSize*ImportLIDs[j];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<MaxElementSize; k++)
	    To[i][jj+k] = EPETRA_MAX( To[i][jj+k], fabs(*ptr++) );
      }
    }
    // Note:  The following form of averaging is not a true average if more that one value is combined.
    //        This might be an issue in the future, but we leave this way for now.
    else if(CombineMode==Average) {
      for (j=0; j<NumImportIDs; j++) {
	jj = MaxElementSize*ImportLIDs[j];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<MaxElementSize; k++)
	    { To[i][jj+k] += *ptr++; To[i][jj+k] *= 0.5;}
      }
    }
  }
    
  // variable element size case

  else {
      
    int SizeOfPacket = NumVectors*MaxElementSize;

    if (CombineMode==Add) {
      for (j=0; j<NumImportIDs; j++) {
	ptr = (double *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<ElementSize; k++)
	    To[i][jj+k] += *ptr++; // Add to existing value
      }
    }
    else  if(CombineMode==Insert){
      for (j=0; j<NumImportIDs; j++) {
	ptr = (double *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<ElementSize; k++)
	    To[i][jj+k] = *ptr++;
      }
    }
    else  if(CombineMode==AbsMax){
      for (j=0; j<NumImportIDs; j++) {
	ptr = (double *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<ElementSize; k++)
	    To[i][jj+k] = EPETRA_MAX( To[i][jj+k], fabs(*ptr++) );
      }
    }
    // Note:  The following form of averaging is not a true average if more that one value is combined.
    //        This might be an issue in the future, but we leave this way for now.
    else if(CombineMode==Average) {
      for (j=0; j<NumImportIDs; j++) {
	ptr = (double *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	for (i=0; i<NumVectors; i++)
	  for (k=0; k<ElementSize; k++)
	    { To[i][jj+k] += *ptr++; To[i][jj+k] *= 0.5;}
      }
    }
  }

  return(0);
}

//=========================================================================
int Epetra_MultiVector::Dot(const Epetra_MultiVector& A, double *Result) const {

  // Dot product of two MultiVectors 

  int i;
  if (NumVectors_ != A.NumVectors()) EPETRA_CHK_ERR(-1);
  if (MyLength_ != A.MyLength()) EPETRA_CHK_ERR(-2);
    
  double **A_Pointers = A.Pointers();

  for (i=0; i < NumVectors_; i++) DoubleTemp_[i] = DOT(MyLength_, Pointers_[i], A_Pointers[i]);
  
  Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  
  UpdateFlops(2*GlobalLength_*NumVectors_);

  return(0);
}
//=========================================================================
int Epetra_MultiVector::Abs(const Epetra_MultiVector& A) {

  // this[i][j] = fabs(A[i][j])

  int i, j;
  if (NumVectors_ != A.NumVectors()) EPETRA_CHK_ERR(-1);
  if (MyLength_ != A.MyLength()) EPETRA_CHK_ERR(-2);

  double **A_Pointers = A.Pointers();

  for (i=0; i < NumVectors_; i++) 
    for (j=0; j < MyLength_; j++)
      Pointers_[i][j] = fabs(A_Pointers[i][j]);

  return(0);
}
//=========================================================================
int Epetra_MultiVector::Reciprocal(const Epetra_MultiVector& A) {

  // this[i][j] = 1.0/(A[i][j])

  int ierr = 0;
  int i, j;
  if (NumVectors_ != A.NumVectors()) EPETRA_CHK_ERR(-1);
  if (MyLength_ != A.MyLength()) EPETRA_CHK_ERR(-2);

  double **A_Pointers = A.Pointers();

  for (i=0; i < NumVectors_; i++) 
    for (j=0; j < MyLength_; j++) {
      double value = A_Pointers[i][j];
      // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)     
      if (fabs(value)<Epetra_MinDouble) {
	if (value==0.0) ierr = 1;
	else if (ierr!=1) ierr = 2;
	Pointers_[i][j] = EPETRA_SGN(value) * Epetra_MaxDouble;
      }
      else
	Pointers_[i][j] = 1.0/value;
    }
  EPETRA_CHK_ERR(ierr);
  return(0);
}
  //=========================================================================
  int Epetra_MultiVector::Scale (double ScalarValue) {

    // scales a MultiVector in place by a scalar
  
    for (int i = 0; i < NumVectors_; i++)
      SCAL(MyLength_, ScalarValue, Pointers_[i]);

    UpdateFlops(GlobalLength_*NumVectors_);

    return(0);
  }

  //=========================================================================
  int Epetra_MultiVector::Scale (double ScalarA, const Epetra_MultiVector& A) {

    // scales a MultiVector by a scalar and put in the this:
    // this = ScalarA * A

  if (NumVectors_ != A.NumVectors()) EPETRA_CHK_ERR(-1);
  if (MyLength_ != A.MyLength()) EPETRA_CHK_ERR(-2);

    double **A_Pointers = (double**)A.Pointers();
    
    for (int i = 0; i < NumVectors_; i++)
      for (int j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarA * A_Pointers[i][j];

    UpdateFlops(GlobalLength_*NumVectors_);

    return(0);
  }

  //=========================================================================
  int Epetra_MultiVector::Update(double ScalarA, const Epetra_MultiVector& A, double ScalarThis) {

    int i, j;

    // linear combination of two MultiVectors: this = ScalarThis * this + ScalarA * A

  if (NumVectors_ != A.NumVectors()) EPETRA_CHK_ERR(-1);
  if (MyLength_ != A.MyLength()) EPETRA_CHK_ERR(-2);

    double **A_Pointers = (double**)A.Pointers();

    if (ScalarThis==0.0)
      {
	for (i = 0; i < NumVectors_; i++)
	  for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarA * A_Pointers[i][j];
	UpdateFlops(GlobalLength_*NumVectors_);
      }
    else if (ScalarThis==1.0)
      {
	for (i = 0; i < NumVectors_; i++)
	  for (j = 0; j < MyLength_; j++) Pointers_[i][j] = Pointers_[i][j] + ScalarA * A_Pointers[i][j];
	UpdateFlops(2*GlobalLength_*NumVectors_);
      }
    else if (ScalarA==1.0)
      {
	for (i = 0; i < NumVectors_; i++)
	  for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarThis * Pointers_[i][j] + A_Pointers[i][j];
	UpdateFlops(2*GlobalLength_*NumVectors_);
      }
    else
      {
	for (i = 0; i < NumVectors_; i++)
	  for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarThis * Pointers_[i][j] +
					                    ScalarA *  A_Pointers[i][j];
	UpdateFlops(3*GlobalLength_*NumVectors_);
      }

    return(0);
  }

//=========================================================================
int Epetra_MultiVector::Update(double ScalarA, const Epetra_MultiVector& A, 
				  double ScalarB, const Epetra_MultiVector& B, double ScalarThis) {
  
  int i, j;
  
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
			   
  if (NumVectors_ != A.NumVectors() || NumVectors_ != B.NumVectors()) EPETRA_CHK_ERR(-1);
  if (MyLength_ != A.MyLength() || MyLength_ != B.MyLength()) EPETRA_CHK_ERR(-2);
  
    double **A_Pointers = (double**)A.Pointers();
    double **B_Pointers = (double**)B.Pointers();

    if (ScalarThis==0.0)
      {
	if (ScalarA==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] =           A_Pointers[i][j] + 
						                ScalarB * B_Pointers[i][j];
	    UpdateFlops(2*GlobalLength_*NumVectors_);
	  }
	else if (ScalarB==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarA * A_Pointers[i][j] +
						                          B_Pointers[i][j];
	    UpdateFlops(2*GlobalLength_*NumVectors_);
	  }
	else
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarA * A_Pointers[i][j] + 
						                ScalarB * B_Pointers[i][j];
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
      }
    else if (ScalarThis==1.0)
      {
	if (ScalarA==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] +=           A_Pointers[i][j] + 
						                 ScalarB * B_Pointers[i][j];
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
	else if (ScalarB==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] += ScalarA * A_Pointers[i][j] +
						                           B_Pointers[i][j];
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
	else
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] += ScalarA * A_Pointers[i][j] + 
						                 ScalarB * B_Pointers[i][j];
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
      }
    else
      {
	if (ScalarA==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  ScalarThis *    Pointers_[i][j] +
						                           A_Pointers[i][j] + 
						                 ScalarB * B_Pointers[i][j];
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
	else if (ScalarB==1.0)
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  ScalarThis *    Pointers_[i][j] +
						                 ScalarA * A_Pointers[i][j] +
						                           B_Pointers[i][j];
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
	else
	  {
	    for (i = 0; i < NumVectors_; i++)
	      for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  ScalarThis *    Pointers_[i][j] +
						                 ScalarA * A_Pointers[i][j] + 
						                 ScalarB * B_Pointers[i][j];
	    UpdateFlops(5*GlobalLength_*NumVectors_);
	  }
      }


    return(0);
  }

//=========================================================================
int  Epetra_MultiVector::Norm1 (double* Result) const {
  
  // 1-norm of each vector in MultiVector 
    
  int i;
  for (i=0; i < NumVectors_; i++) DoubleTemp_[i] = ASUM(MyLength_, Pointers_[i]);
  
  Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  
  UpdateFlops(2*GlobalLength_*NumVectors_);

  return(0);
}

//=========================================================================
int  Epetra_MultiVector::Norm2 (double* Result) const {
  
  // 2-norm of each vector in MultiVector 
  
  int i, j;
  for (i=0; i < NumVectors_; i++) 
    {
      double sum = 0.0;
      for (j=0; j < MyLength_; j++) sum += Pointers_[i][j] * Pointers_[i][j];
      DoubleTemp_[i] = sum;
    }
  Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  for (i=0; i < NumVectors_; i++) Result[i] = sqrt(Result[i]);
  
  UpdateFlops(2*GlobalLength_*NumVectors_);
  
  return(0);
}

//=========================================================================
int  Epetra_MultiVector::NormInf (double* Result) const {
  
  // Inf-norm of each vector in MultiVector 
    
  int i, j;
  for (i=0; i < NumVectors_; i++) 
    {
      j = IAMAX(MyLength_, Pointers_[i]);
      DoubleTemp_[i] = fabs(Pointers_[i][j]);
    }
  Comm_->MaxAll(DoubleTemp_, Result, NumVectors_);
  
  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine  
  return(0);
}

//=========================================================================
int  Epetra_MultiVector::NormWeighted (const Epetra_MultiVector& Weights, double* Result) const {
  
  // Weighted 2-norm of each vector in MultiVector 

  // If only one vector in Weights, we assume it will be used as the weights for all vectors

  int i, j;
  bool OneW = false;
  if (Weights.NumVectors()==1) OneW = true;
  else if (NumVectors_ != Weights.NumVectors()) EPETRA_CHK_ERR(-1);

  if (MyLength_ != Weights.MyLength()) EPETRA_CHK_ERR(-2);

  double *W = Weights.Values();
  double **W_Pointers = Weights.Pointers();
  
  for (i=0; i < NumVectors_; i++) 
    {
      if (!OneW) W = W_Pointers[i]; // If Weights has the same number of vectors as this, use each weight vector
      double sum = 0.0;
      for (j=0; j < MyLength_; j++) {
        double tmp = Pointers_[i][j]/W[j];
        sum += tmp * tmp;
      }
      DoubleTemp_[i] = sum;
    }
  Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  double OneOverN = 1.0 / (double) GlobalLength_;
  for (i=0; i < NumVectors_; i++) Result[i] = sqrt(Result[i]*OneOverN);
  
  UpdateFlops(3*GlobalLength_*NumVectors_);
  
  return(0);
  }

//=========================================================================
int  Epetra_MultiVector::MinValue (double* Result) const {
  
  // Minimum value of each vector in MultiVector 
  
  int i, j, ierr = 0;

  for (i=0; i < NumVectors_; i++) 
    {
      double MinVal = Epetra_MaxDouble;
      if (MyLength_>0) MinVal = Pointers_[i][0];
      for (j=0; j< MyLength_; j++) MinVal = EPETRA_MIN(MinVal,Pointers_[i][j]); 
      DoubleTemp_[i] = MinVal;
    }
  Comm_->MinAll(DoubleTemp_, Result, NumVectors_);

  for (i=0; i < NumVectors_; i++) 
    if (Result[i]==Epetra_MaxDouble) ierr = -1; // Report a problem, numbers are too small
  
  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine
  
  return(ierr);
}

//=========================================================================
int  Epetra_MultiVector::MaxValue (double* Result) const {
  
  // Maximum value of each vector in MultiVector 
  
  int i, j, ierr = 0;
  for (i=0; i < NumVectors_; i++) 
    {
      double MaxVal = - Epetra_MaxDouble;
      if (MyLength_>0) MaxVal = Pointers_[i][0];
      for (j=0; j< MyLength_; j++) MaxVal = EPETRA_MAX(MaxVal,Pointers_[i][j]); 
      DoubleTemp_[i] = MaxVal;
    }
  Comm_->MaxAll(DoubleTemp_, Result, NumVectors_);

  for (i=0; i < NumVectors_; i++) 
    if (Result[i]==Epetra_MinDouble) ierr = -1; // Report a problem, numbers are too small
  
  
  // UpdateFlops(0);  Strictly speaking there are not FLOPS in this routine
  
  return(ierr);
}

//=========================================================================
int  Epetra_MultiVector::MeanValue (double* Result) const {
  
  // Mean value of each vector in MultiVector 
  
  int i, j;
  double fGlobalLength = 1.0/EPETRA_MAX((double) GlobalLength_, 1.0);
  
  for (i=0; i < NumVectors_; i++) 
    {
      double sum = 0.0;
      for (j=0; j < MyLength_; j++) sum += Pointers_[i][j];
      DoubleTemp_[i] = sum;
    }
  Comm_->SumAll(DoubleTemp_, Result, NumVectors_);
  for (i=0; i < NumVectors_; i++) Result[i] = Result[i]*fGlobalLength;
  
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

    if( MyLength_      != A_nrows     ||   // RAB: 2002/01/25: Minor reformat to allow
		A_ncols        != B_nrows     ||   //   setting breakpoint at error return.
		NumVectors_    != B_ncols  )
		EPETRA_CHK_ERR(-2); // Return error

    bool A_is_local = (!A.DistributedGlobal());
    bool B_is_local = (!B.DistributedGlobal());
    bool C_is_local = (!DistributedGlobal_);
    bool Case1 = ( A_is_local &&  B_is_local &&  C_is_local);  // Case 1 above
    bool Case2 = (!A_is_local && !B_is_local &&  C_is_local && TransA=='T' );// Case 2
    bool Case3 = (!A_is_local &&  B_is_local && !C_is_local && TransA=='N');// Case 3
  
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
    	
    
	int m = MyLength_;
	int n = NumVectors_;
	int k = A_ncols;
	int lda = A_tmp->Stride();
	int ldb = B_tmp->Stride();
	int ldc = C_tmp->Stride();
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
  
  int i, j;
  
  // Hadamard product of two MultiVectors: 
  // this = ScalarThis * this + ScalarAB * A * B (element-wise)
  
  if (ScalarAB==0.0) {
    EPETRA_CHK_ERR(Scale(ScalarThis));
    return(0);
  }
			   
  if (A.NumVectors() != 1 && A.NumVectors() != B.NumVectors()) EPETRA_CHK_ERR(-1); // A must have one column or be the same as B.
  if (NumVectors_ != B.NumVectors()) EPETRA_CHK_ERR(-2);
  if (MyLength_ != A.MyLength() || MyLength_ != B.MyLength()) EPETRA_CHK_ERR(-3);
  
  double **A_Pointers = (double**)A.Pointers();
  double **B_Pointers = (double**)B.Pointers();

  int IncA = 1;
  if (A.NumVectors() == 1 ) IncA = 0;

    if (ScalarThis==0.0) {
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  A[j] * B_Pointers[i][j];
	  }
	  UpdateFlops(GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarAB * A[j] *
					      B_Pointers[i][j];
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
    }
    else if (ScalarThis==1.0) {
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] +=  A[j] * B_Pointers[i][j];
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] += ScalarAB * A[j] *
					      B_Pointers[i][j];
	    }
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
    }
    else { // if (ScalarThis!=1.0 && ScalarThis !=0 ) 
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  ScalarThis * Pointers_[i][j] + A[j] * B_Pointers[i][j];
	  }
	  UpdateFlops(3*GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarThis * Pointers_[i][j] + ScalarAB * A[j] *
					      B_Pointers[i][j];
	    }
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
    }
  return(0);
}
//=========================================================================
int Epetra_MultiVector::ReciprocalMultiply(double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B,
		       double ScalarThis) {
  
  int i, j;
  
  // Hadamard product of two MultiVectors: 
  // this = ScalarThis * this + ScalarAB * B / A (element-wise)
  
  if (ScalarAB==0.0) {
    EPETRA_CHK_ERR(Scale(ScalarThis));
    return(0);
  }
			   
  if (A.NumVectors() != 1 && A.NumVectors() != B.NumVectors()) EPETRA_CHK_ERR(-1); // A must have one column or be the same as B.
  if (NumVectors_ != B.NumVectors()) EPETRA_CHK_ERR(-2);
  if (MyLength_ != A.MyLength() || MyLength_ != B.MyLength()) EPETRA_CHK_ERR(-3);
  
  double **A_Pointers = (double**)A.Pointers();
  double **B_Pointers = (double**)B.Pointers();

  int IncA = 1;
  if (A.NumVectors() == 1 ) IncA = 0;

    if (ScalarThis==0.0) {
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] = B_Pointers[i][j] / A[j];
	  }
	  UpdateFlops(GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarAB * 
					      B_Pointers[i][j] / A[j];
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
    }
    else if (ScalarThis==1.0) {
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] +=  B_Pointers[i][j] / A[j];
	  }
	  UpdateFlops(2*GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] += ScalarAB *
					      B_Pointers[i][j] / A[j];
	    }
	    UpdateFlops(3*GlobalLength_*NumVectors_);
	  }
    }
    else { // if (ScalarThis!=1.0 && ScalarThis !=0 ) 
      if (ScalarAB==1.0)
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] =  ScalarThis * Pointers_[i][j] + B_Pointers[i][j] / A[j];
	  }
	  UpdateFlops(3*GlobalLength_*NumVectors_);
	}
      else
	{
	  for (i = 0; i < NumVectors_; i++) {
	    double * A = A_Pointers[i*IncA];
	    for (j = 0; j < MyLength_; j++) Pointers_[i][j] = ScalarThis * Pointers_[i][j] + ScalarAB * 
					      B_Pointers[i][j] / A[j];
	    }
	    UpdateFlops(4*GlobalLength_*NumVectors_);
	  }
    }
  return(0);
}

/*
//=======================================================================
double*& Epetra_MultiVector::operator [] (int index)  {
  
  //  Epetra_MultiVector::operator [] --- return non-const reference 
  
  if (index < 0 || index >=NumVectors_) 
    throw ReportError("Vector index = " + toString(index) + "is out of range. Number of Vectors = " + toString(NumVectors_), -1);
  return(Pointers_[index]);
}

//=======================================================================
const double*& Epetra_MultiVector::operator [] (int index) const {
  
  //  Epetra_MultiVector::operator [] --- return non-const reference 

  if (index < 0 || index >=NumVectors_) 
    throw ReportError("Vector index = " + toString(index) + "is out of range. Number of Vectors = " + toString(NumVectors_), -1);

  const double * & temp = (const double * &) (Pointers_[index]);
  return(temp);
}
*/

//=======================================================================
Epetra_Vector *& Epetra_MultiVector::operator () (int index)  {
  
  //  Epetra_MultiVector::operator [] --- return non-const reference 
  
  if (index < 0 || index >=NumVectors_) 
    throw ReportError("Vector index = " + toString(index) + "is out of range. Number of Vectors = " + toString(NumVectors_), -1);

  // Create a new Epetra_Vector that is a view of ith vector, if not already present
  if (Vectors_[index]==0)
    Vectors_[index] = new Epetra_Vector(View, Map(), Pointers_[index]);
  return(Vectors_[index]);
}

//=======================================================================
const Epetra_Vector *& Epetra_MultiVector::operator () (int index) const {
  
  //  Epetra_MultiVector::operator [] --- return non-const reference 

  if (index < 0 || index >=NumVectors_) 
    throw ReportError("Vector index = " + toString(index) + "is out of range. Number of Vectors = " + toString(NumVectors_), -1);

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
  
  if (NumVectors_ != A.NumVectors())
    throw ReportError("Number of vectors incompatible in Assign.  The this MultiVector has NumVectors = " + toString(NumVectors_)
		      + ".  The A MultiVector has NumVectors = " + toString(A.NumVectors()), -3);
  if (MyLength_ != A.MyLength())
    throw ReportError("Length of MultiVectors incompatible in Assign.  The this MultiVector has MyLength = " + toString(MyLength_)
		      + ".  The A MultiVector has MyLength = " + toString(A.MyLength()), -4);
  
  double ** A_Pointers = A.Pointers();
  for (int i = 0; i< NumVectors_; i++)
      for (int j=0; j<MyLength_; j++) Pointers_[i][j] = A_Pointers[i][j];
    return;    
  }

  //=========================================================================
  int  Epetra_MultiVector::Reduce() {

    // Global reduction on each entry of a Replicated Local MultiVector

    int i, j;
    double * source = new double[MyLength_*NumVectors_];
    double * target = 0;
    bool packed = (ConstantStride_ && (Stride_==MyLength_));
    if (packed) {
       for (i=0; i<MyLength_*NumVectors_; i++) source[i] = Values_[i];
       target = Values_;
    }
    else {
       double * tmp1 = source;
       for (i = 0; i < NumVectors_; i++) {
         double * tmp2 = Pointers_[i];
         for (j=0; j< MyLength_; j++) *tmp1++ = *tmp2++;
       }
       target = new double[MyLength_*NumVectors_];
    }

    Comm_->SumAll(source, target, MyLength_*NumVectors_);
    delete [] source;
    if (!packed) {
       double * tmp2 = target;
       for (i = 0; i < NumVectors_; i++) {
         double * tmp1 = Pointers_[i];
         for (j=0; j< MyLength_; j++) *tmp1++ = *tmp2++;
       }
       delete [] target;
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
      int * MyGlobalElements1 = Map().MyGlobalElements();
      int * FirstPointInElementList1;
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
	for (int ii=0; ii< Map().ElementSize(ii); ii++) {
       int iii;
	  os.width(10);
	  os <<  MyPID; os << "    ";
	  os.width(10);
	  if (MaxElementSize1==1) {
	    os << MyGlobalElements1[i] << "    ";
       iii = i;
       }
	  else {
	    os <<  MyGlobalElements1[i]<< "/" << ii << "    ";
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
