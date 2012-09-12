
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

#include "Epetra_LongLongVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include <limits>

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

//=============================================================================
Epetra_LongLongVector::Epetra_LongLongVector(const Epetra_BlockMap& map, bool zeroOut)
  : Epetra_DistObject(map, "Epetra::LongLongVector"),
    Values_(0),
    UserAllocated_(false),
    Allocated_(false)
{
  if(!map.GlobalIndicesLongLong())
     throw ReportError("Epetra_IntVector::Epetra_IntVector: cannot be called with non long long map index type", -1);

  AllocateForCopy();
  if(zeroOut) PutValue(0); // Zero out values
}
//=============================================================================
Epetra_LongLongVector::Epetra_LongLongVector(const Epetra_LongLongVector& Source)
  : Epetra_DistObject(Source),
    Values_(0),
    UserAllocated_(false),
    Allocated_(false)
{
  if(!Source.Map().GlobalIndicesLongLong())
     throw ReportError("Epetra_IntVector::Epetra_IntVector: cannot be called with non long long map index type", -1);

  AllocateForCopy();
  DoCopy(Source.Values_);
}
//=============================================================================
Epetra_LongLongVector::Epetra_LongLongVector(Epetra_DataAccess CV, const Epetra_BlockMap& map, long long *V)
  : Epetra_DistObject(map, "Epetra::LongLongVector"),
    Values_(0),
    UserAllocated_(false),
    Allocated_(false)
{
  if(!map.GlobalIndicesLongLong())
     throw ReportError("Epetra_IntVector::Epetra_IntVector: cannot be called with non long long map index type", -1);

  if (CV==Copy) {
    AllocateForCopy();
    DoCopy(V);
  }
  else {
    AllocateForView();
    DoView(V);
  }
}
//=========================================================================
Epetra_LongLongVector::~Epetra_LongLongVector(){


 if (Allocated_ && (!UserAllocated_) && Values_!=0) delete [] Values_;
}

//=========================================================================
int Epetra_LongLongVector::AllocateForCopy()
{
  
  if (Allocated_) return(0);
    
  int myLength = MyLength();
  if (myLength>0)
    Values_ = new long long[myLength];
  else
    Values_ = 0;

  Allocated_ = true;
  UserAllocated_ = false;
  return(0);
}

//=========================================================================
int Epetra_LongLongVector::DoCopy(long long * V)
{
  int iend = MyLength();
  for (int i=0; i<iend; i++) Values_[i] = V[i];

  return(0);
}
//=========================================================================
int Epetra_LongLongVector::AllocateForView()
{

  Allocated_ = true;
  UserAllocated_ = true;
  
  return(0);
}

//=========================================================================
int Epetra_LongLongVector::DoView(long long * V)
{
  
  Values_ = V;

  return(0);
}
//=============================================================================
int Epetra_LongLongVector::ExtractCopy(long long *V) const {
  
  int iend = MyLength();
  for (int i=0; i<iend; i++) V[i] = Values_[i];
  return(0);
}

//=============================================================================
int Epetra_LongLongVector::ExtractView(long long **V) const
{
  *V = Values_;
  return(0);
}

//=============================================================================
int Epetra_LongLongVector::PutValue(long long Value) {
  int iend = MyLength();
  for (int i=0; i<iend; i++) Values_[i] = Value;
  return(0);
}
//=============================================================================
long long Epetra_LongLongVector::MaxValue() {
  long long result = std::numeric_limits<long long>::min(); // smallest 64 bit int
  int iend = MyLength();
  if (iend>0) result = Values_[0];
  for (int i=0; i<iend; i++) result = EPETRA_MAX(result, Values_[i]);
  long long globalResult;
  this->Comm().MaxAll(&result, &globalResult, 1);
  return(globalResult);
}
//=============================================================================
long long Epetra_LongLongVector::MinValue() {
  long long result = std::numeric_limits<long long>::max(); // largest 64 bit int
  int iend = MyLength();
  if (iend>0) result = Values_[0];
  for (int i=0; i<iend; i++) result = EPETRA_MIN(result, Values_[i]);
  long long globalResult;
  this->Comm().MinAll(&result, &globalResult, 1);
  return(globalResult);
}
//========================================================================
Epetra_LongLongVector& Epetra_LongLongVector::operator = (const Epetra_LongLongVector& V) {
  

  if (MyLength() != V.MyLength())
    throw ReportError("Length of LongLongVectors incompatible in Assign.  The this LongLongVector has MyLength = " + toString(MyLength())
          + ".  The V LongLongVector has MyLength = " + toString(V.MyLength()), -1);
  
  int iend = MyLength();
  for (int i=0; i<iend; i++) Values_[i] =V[i];
  return(*this);
}

void Epetra_LongLongVector::Print(ostream& os) const {
  int MyPID = Map().Comm().MyPID();
  int NumProc = Map().Comm().NumProc();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumMyElements1 =Map(). NumMyElements();
      int MaxElementSize1 = Map().MaxElementSize();
    long long * MyGlobalElements1 = Map().MyGlobalElements64();
      int * FirstPointInElementList1=0;
      if (MaxElementSize1!=1) FirstPointInElementList1 = Map().FirstPointInElementList();

      if (MyPID==0) {
  os.width(8);
  os <<  "     MyPID"; os << "    ";
  os.width(12);
  if (MaxElementSize1==1)
    os <<  "GID  ";
  else
    os <<  "     GID/Point";
  os.width(20);
  os <<  "Value  ";
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
        os.width(20);
        os <<  Values_[iii];
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
//=========================================================================
int Epetra_LongLongVector::CheckSizes(const Epetra_SrcDistObject& Source)
{
  (void)Source;
  return(0);
}

//=========================================================================
int Epetra_LongLongVector::CopyAndPermute(const Epetra_SrcDistObject& Source,
                                     int NumSameIDs, 
                                     int NumPermuteIDs,
                                     int * PermuteToLIDs, 
                                     int *PermuteFromLIDs,
                                     const Epetra_OffsetIndex * Indexor)
{
  (void)Indexor;
  const Epetra_LongLongVector & A = dynamic_cast<const Epetra_LongLongVector &>(Source);

  long long * From;
  A.ExtractView(&From);
  long long *To = Values_;

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
  int j, jj, jjj, k;
  
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
  for (j=0; j<NumSameEntries; j++)
    To[j] = From[j];
    }
  // Do local permutation next
  if (NumPermuteIDs>0) {
  
    // Point entry case
    if (Case1) {
      
      for (j=0; j<NumPermuteIDs; j++) 
  To[PermuteToLIDs[j]] = From[PermuteFromLIDs[j]];
    }
    // constant element size case
    else if (Case2) {
      
      for (j=0; j<NumPermuteIDs; j++) {
  jj = MaxElementSize*PermuteToLIDs[j];
  jjj = MaxElementSize*PermuteFromLIDs[j];
    for (k=0; k<MaxElementSize; k++)
      To[jj+k] = From[jjj+k];
      }
    }
    
    // variable element size case
    else {
      
      for (j=0; j<NumPermuteIDs; j++) {
  jj = ToFirstPointInElementList[PermuteToLIDs[j]];
  jjj = FromFirstPointInElementList[PermuteFromLIDs[j]];
  int ElementSize = FromElementSizeList[PermuteFromLIDs[j]];
    for (k=0; k<ElementSize; k++)
      To[jj+k] = From[jjj+k];
      }
    }
  }
  return(0);
}

//=========================================================================
int Epetra_LongLongVector::PackAndPrepare(const Epetra_SrcDistObject & Source,
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
  const Epetra_LongLongVector & A = dynamic_cast<const Epetra_LongLongVector &>(Source);

  int j, jj, k;

  long long  * From;
  A.ExtractView(&From);
  int MaxElementSize = Map().MaxElementSize();
  bool ConstantElementSize = Map().ConstantElementSize();

  int * FromFirstPointInElementList = 0;
  int * FromElementSizeList = 0;

  if (!ConstantElementSize) {
    FromFirstPointInElementList = A.Map().FirstPointInElementList();
    FromElementSizeList = A.Map().ElementSizeList();
  }

  SizeOfPacket = MaxElementSize * (int)sizeof(long long); 

  if(NumExportIDs*SizeOfPacket>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = NumExportIDs*SizeOfPacket;
    Exports = new char[LenExports];
  }

  long long * ptr;

  if (NumExportIDs>0) {
    ptr = (long long *) Exports;
    
    // Point entry case
    if (MaxElementSize==1) for (j=0; j<NumExportIDs; j++) *ptr++ = From[ExportLIDs[j]];

    // constant element size case
    else if (ConstantElementSize) {
      
      for (j=0; j<NumExportIDs; j++) {
  jj = MaxElementSize*ExportLIDs[j];
    for (k=0; k<MaxElementSize; k++)
      *ptr++ = From[jj+k];
      }
    }
    
    // variable element size case
    else {
      
      SizeOfPacket = MaxElementSize;
      for (j=0; j<NumExportIDs; j++) {
  ptr = (long long *) Exports + j*SizeOfPacket;
  jj = FromFirstPointInElementList[ExportLIDs[j]];
  int ElementSize = FromElementSizeList[ExportLIDs[j]];
    for (k=0; k<ElementSize; k++)
      *ptr++ = From[jj+k];
      }
    }
  }

  return(0);
}

//=========================================================================
int Epetra_LongLongVector::UnpackAndCombine(const Epetra_SrcDistObject & Source,
                                       int NumImportIDs,
                                       int * ImportLIDs, 
                                       int LenImports, 
                                       char * Imports,
                                       int & SizeOfPacket, 
                                       Epetra_Distributor & Distor, 
                                       Epetra_CombineMode CombineMode,
                                       const Epetra_OffsetIndex * Indexor)
{
  (void)Source;
  (void)LenImports;
  (void)Distor;
  (void)Indexor;
  int j, jj, k;
  
  if(    CombineMode != Add
      && CombineMode != Zero
      && CombineMode != Insert
      && CombineMode != Average
      && CombineMode != AbsMax )
    EPETRA_CHK_ERR(-1); //Unsupported CombinedMode, will default to Zero

  if (NumImportIDs<=0) return(0);

  long long * To = Values_;
  int MaxElementSize = Map().MaxElementSize();
  bool ConstantElementSize = Map().ConstantElementSize();

  int * ToFirstPointInElementList = 0;
  int * ToElementSizeList = 0;

  if (!ConstantElementSize) {
    ToFirstPointInElementList = Map().FirstPointInElementList();
    ToElementSizeList = Map().ElementSizeList();
  }
  
  long long * ptr;
  // Unpack it...

  ptr = (long long *) Imports;
    
  // Point entry case
  if (MaxElementSize==1) {
      
      if (CombineMode==Add)
  for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] += *ptr++; // Add to existing value
      else if(CombineMode==Insert)
  for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] = *ptr++;
      else if(CombineMode==AbsMax)
        for (j=0; j<NumImportIDs; j++) {
    To[ImportLIDs[j]] = EPETRA_MAX( To[ImportLIDs[j]],(*ptr > 0 ? *ptr : -*ptr)); // CJ: No std::abs(long long) in VS2005
    ptr++;
  }
      // Note:  The following form of averaging is not a true average if more that one value is combined.
      //        This might be an issue in the future, but we leave this way for now.
      else if(CombineMode==Average)
  for (j=0; j<NumImportIDs; j++) {To[ImportLIDs[j]] += *ptr++; To[ImportLIDs[j]] /= 2;}
  }

  // constant element size case

  else if (ConstantElementSize) {
   
    if (CombineMode==Add) {
      for (j=0; j<NumImportIDs; j++) {
  jj = MaxElementSize*ImportLIDs[j];
    for (k=0; k<MaxElementSize; k++)
      To[jj+k] += *ptr++; // Add to existing value
      }
    }
    else if(CombineMode==Insert) {
      for (j=0; j<NumImportIDs; j++) {
  jj = MaxElementSize*ImportLIDs[j];
    for (k=0; k<MaxElementSize; k++)
      To[jj+k] = *ptr++;
      }
    }
    else if(CombineMode==AbsMax) {
      for (j=0; j<NumImportIDs; j++) {
  jj = MaxElementSize*ImportLIDs[j];
  for (k=0; k<MaxElementSize; k++) {
      To[jj+k] = EPETRA_MAX( To[jj+k], (*ptr > 0 ? *ptr : -*ptr)); // CJ: No std::abs(long long) in VS2005
      ptr++;
  }
      }
    }
    // Note:  The following form of averaging is not a true average if more that one value is combined.
    //        This might be an issue in the future, but we leave this way for now.
    else if(CombineMode==Average) {
      for (j=0; j<NumImportIDs; j++) {
  jj = MaxElementSize*ImportLIDs[j];
    for (k=0; k<MaxElementSize; k++)
      { To[jj+k] += *ptr++; To[jj+k] /= 2;}
      }
    }
  }
    
  // variable element size case

  else {
      
    SizeOfPacket = MaxElementSize;

    if (CombineMode==Add) {
      for (j=0; j<NumImportIDs; j++) {
  ptr = (long long *) Imports + j*SizeOfPacket;
  jj = ToFirstPointInElementList[ImportLIDs[j]];
  int ElementSize = ToElementSizeList[ImportLIDs[j]];
    for (k=0; k<ElementSize; k++)
      To[jj+k] += *ptr++; // Add to existing value
      }
    }
    else  if(CombineMode==Insert){
      for (j=0; j<NumImportIDs; j++) {
  ptr = (long long *) Imports + j*SizeOfPacket;
  jj = ToFirstPointInElementList[ImportLIDs[j]];
  int ElementSize = ToElementSizeList[ImportLIDs[j]];
    for (k=0; k<ElementSize; k++)
      To[jj+k] = *ptr++;
      }
    }
    else  if(CombineMode==AbsMax){
      for (j=0; j<NumImportIDs; j++) {
  ptr = (long long *) Imports + j*SizeOfPacket;
  jj = ToFirstPointInElementList[ImportLIDs[j]];
  int ElementSize = ToElementSizeList[ImportLIDs[j]];
  for (k=0; k<ElementSize; k++) {
      To[jj+k] = EPETRA_MAX( To[jj+k], (*ptr > 0 ? *ptr : -*ptr)); // CJ: No std::abs(long long) in VS2005
      ptr++;
  }
      }
    }
    // Note:  The following form of averaging is not a true average if more that one value is combined.
    //        This might be an issue in the future, but we leave this way for now.
    else if(CombineMode==Average) {
      for (j=0; j<NumImportIDs; j++) {
  ptr = (long long *) Imports + j*SizeOfPacket;
  jj = ToFirstPointInElementList[ImportLIDs[j]];
  int ElementSize = ToElementSizeList[ImportLIDs[j]];
    for (k=0; k<ElementSize; k++)
      { To[jj+k] += *ptr++; To[jj+k] /= 2;}
      }
    }
  }
  
  return(0);
}

#endif // EPETRA_NO_64BIT_GLOBAL_INDICES
