
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

#include "Epetra_IntVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
//=============================================================================
Epetra_IntVector::Epetra_IntVector(const Epetra_BlockMap& Map, bool zeroOut)
  : Epetra_DistObject(Map, "Epetra::IntVector"),
    UserAllocated_(false),
    Allocated_(false)
{
  AllocateForCopy();
  if(zeroOut) PutValue(0); // Zero out values
}
//=============================================================================
Epetra_IntVector::Epetra_IntVector(const Epetra_IntVector& Source)
  : Epetra_DistObject(Source),
    UserAllocated_(false),
    Allocated_(false)
{
  AllocateForCopy();
  DoCopy(Source.Values_);
}
//=============================================================================
Epetra_IntVector::Epetra_IntVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, int *V)
  : Epetra_DistObject(Map, "Epetra::IntVector"),
    UserAllocated_(false),
    Allocated_(false)
{
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
Epetra_IntVector::~Epetra_IntVector(){


 if (Allocated_ && (!UserAllocated_) && Values_!=0) delete [] Values_;
}

//=========================================================================
int Epetra_IntVector::AllocateForCopy()
{
  
  if (Allocated_) return(0);
    
  int myLength = MyLength();
  if (myLength>0)
    Values_ = new int[myLength];
  else
    Values_ = 0;

  Allocated_ = true;
  UserAllocated_ = false;
  return(0);
}

//=========================================================================
int Epetra_IntVector::DoCopy(int * V)
{
  int iend = MyLength();
  for (int i=0; i<iend; i++) Values_[i] = V[i];

  return(0);
}
//=========================================================================
int Epetra_IntVector::AllocateForView()
{

  Allocated_ = true;
  UserAllocated_ = true;
  
  return(0);
}

//=========================================================================
int Epetra_IntVector::DoView(int * V)
{
  
  Values_ = V;

  return(0);
}
//=============================================================================
int Epetra_IntVector::ExtractCopy(int *V) const {
  
  int iend = MyLength();
  for (int i=0; i<iend; i++) V[i] = Values_[i];
  return(0);
}

//=============================================================================
int Epetra_IntVector::ExtractView(int **V) const
{
  *V = Values_;
  return(0);
}

//=============================================================================
int Epetra_IntVector::PutValue(int Value) {
  int iend = MyLength();
  for (int i=0; i<iend; i++) Values_[i] = Value;
  return(0);
}
//=============================================================================
int Epetra_IntVector::MaxValue() {

  int result = -2000000000; // Negative 2 billion is close to smallest 32 bit int
  int iend = MyLength();
  if (iend>0) result = Values_[0];
  for (int i=0; i<iend; i++) result = EPETRA_MAX(result, Values_[i]);
  int globalResult;
  this->Comm().MaxAll(&result, &globalResult, 1);
  return(globalResult);
}
//=============================================================================
int Epetra_IntVector::MinValue() {

  int result = 2000000000; // 2 billion is close to largest 32 bit int
  int iend = MyLength();
  if (iend>0) result = Values_[0];
  for (int i=0; i<iend; i++) result = EPETRA_MIN(result, Values_[i]);
  int globalResult;
  this->Comm().MinAll(&result, &globalResult, 1);
  return(globalResult);
}
//========================================================================
Epetra_IntVector& Epetra_IntVector::operator = (const Epetra_IntVector& V) {
  

  if (MyLength() != V.MyLength())
    throw ReportError("Length of IntVectors incompatible in Assign.  The this IntVector has MyLength = " + toString(MyLength())
		      + ".  The V IntVector has MyLength = " + toString(V.MyLength()), -1);
  
  int iend = MyLength();
  for (int i=0; i<iend; i++) Values_[i] =V[i];
  return(*this);
}

void Epetra_IntVector::Print(ostream& os) const {
  int MyPID = Map().Comm().MyPID();
  int NumProc = Map().Comm().NumProc();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumMyElements1 =Map(). NumMyElements();
      int MaxElementSize1 = Map().MaxElementSize();
      int * MyGlobalElements1 = Map().MyGlobalElements();
      int * FirstPointInElementList1;
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
int Epetra_IntVector::CheckSizes(const Epetra_SrcDistObject& Source) {
  return(0);
}

//=========================================================================
int Epetra_IntVector::CopyAndPermute(const Epetra_SrcDistObject& Source,
                                     int NumSameIDs, 
                                     int NumPermuteIDs,
                                     int * PermuteToLIDs, 
                                     int *PermuteFromLIDs,
                                     const Epetra_OffsetIndex * Indexor) {

  const Epetra_IntVector & A = dynamic_cast<const Epetra_IntVector &>(Source);

  int * From;
  A.ExtractView(&From);
  int *To = Values_;

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
int Epetra_IntVector::PackAndPrepare(const Epetra_SrcDistObject & Source,
                                     int NumExportIDs,
                                     int * ExportLIDs,
                                     int & LenExports,
                                     char * & Exports,
                                     int & SizeOfPacket,
                                     int * Sizes,
                                     bool & VarSizes,
                                     Epetra_Distributor & Distor) {

  const Epetra_IntVector & A = dynamic_cast<const Epetra_IntVector &>(Source);

  int j, jj, k;

  int  * From;
  A.ExtractView(&From);
  int MaxElementSize = Map().MaxElementSize();
  bool ConstantElementSize = Map().ConstantElementSize();

  int * FromFirstPointInElementList = 0;
  int * FromElementSizeList = 0;

  if (!ConstantElementSize) {
    FromFirstPointInElementList = A.Map().FirstPointInElementList();
    FromElementSizeList = A.Map().ElementSizeList();
  }

  SizeOfPacket = MaxElementSize * sizeof(int); 

  if(NumExportIDs*SizeOfPacket>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = NumExportIDs*SizeOfPacket;
    Exports = new char[LenExports];
  }

  int * ptr;

  if (NumExportIDs>0) {
    ptr = (int *) Exports;
    
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
      
      int SizeOfPacket = MaxElementSize;
      for (j=0; j<NumExportIDs; j++) {
	ptr = (int *) Exports + j*SizeOfPacket;
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
int Epetra_IntVector::UnpackAndCombine(const Epetra_SrcDistObject & Source,
                                       int NumImportIDs,
                                       int * ImportLIDs, 
                                       int LenImports, 
                                       char * Imports,
                                       int & SizeOfPacket, 
                                       Epetra_Distributor & Distor, 
                                       Epetra_CombineMode CombineMode,
                                       const Epetra_OffsetIndex * Indexor) {
  int j, jj, k;
  
  if(    CombineMode != Add
      && CombineMode != Zero
      && CombineMode != Insert
      && CombineMode != Average
      && CombineMode != AbsMax )
    EPETRA_CHK_ERR(-1); //Unsupported CombinedMode, will default to Zero

  if (NumImportIDs<=0) return(0);

  int * To = Values_;
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
      
      if (CombineMode==Add)
	for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] += *ptr++; // Add to existing value
      else if(CombineMode==Insert)
	for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] = *ptr++;
      else if(CombineMode==AbsMax)
        for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] = EPETRA_MAX( To[ImportLIDs[j]],abs(*ptr++)); //
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
	  for (k=0; k<MaxElementSize; k++)
	    To[jj+k] = EPETRA_MAX( To[jj+k], abs(*ptr++) );
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
      
    int SizeOfPacket = MaxElementSize;

    if (CombineMode==Add) {
      for (j=0; j<NumImportIDs; j++) {
	ptr = (int *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	  for (k=0; k<ElementSize; k++)
	    To[jj+k] += *ptr++; // Add to existing value
      }
    }
    else  if(CombineMode==Insert){
      for (j=0; j<NumImportIDs; j++) {
	ptr = (int *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	  for (k=0; k<ElementSize; k++)
	    To[jj+k] = *ptr++;
      }
    }
    else  if(CombineMode==AbsMax){
      for (j=0; j<NumImportIDs; j++) {
	ptr = (int *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	  for (k=0; k<ElementSize; k++)
	    To[jj+k] = EPETRA_MAX( To[jj+k], abs(*ptr++) );
      }
    }
    // Note:  The following form of averaging is not a true average if more that one value is combined.
    //        This might be an issue in the future, but we leave this way for now.
    else if(CombineMode==Average) {
      for (j=0; j<NumImportIDs; j++) {
	ptr = (int *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	  for (k=0; k<ElementSize; k++)
	    { To[jj+k] += *ptr++; To[jj+k] /= 2;}
      }
    }
  }
  
  return(0);
}

