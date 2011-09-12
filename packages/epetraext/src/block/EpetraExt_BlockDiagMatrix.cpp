/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#include "EpetraExt_BlockDiagMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Distributor.h"

#define AM_MULTIPLY 0
#define AM_INVERT   1
#define AM_FACTOR   2 

//=========================================================================
// Constructor
EpetraExt_BlockDiagMatrix::EpetraExt_BlockDiagMatrix(const Epetra_BlockMap& Map,bool zero_out)
  : Epetra_DistObject(Map, "EpetraExt::BlockDiagMatrix"),
    HasComputed_(false),
    ApplyMode_(AM_MULTIPLY),
    DataMap_(0),
    Values_(0),
    Pivots_(0)
{
  Allocate();
  if(zero_out) PutScalar(0.0);
}


//=========================================================================
// Destructor
EpetraExt_BlockDiagMatrix::~EpetraExt_BlockDiagMatrix(){
  if(DataMap_) delete DataMap_;
  if(Pivots_) delete [] Pivots_;
  if(Values_) delete [] Values_;
}


//=========================================================================
// Copy constructor.  
EpetraExt_BlockDiagMatrix::EpetraExt_BlockDiagMatrix(const EpetraExt_BlockDiagMatrix& Source)
  : Epetra_DistObject(Source),
    HasComputed_(false),
    ApplyMode_(AM_MULTIPLY),
    Values_(0),
    Pivots_(0)
{
  DataMap_=new Epetra_BlockMap(*Source.DataMap_);
  Pivots_=new int[NumMyUnknowns()];
  Values_=new double[DataMap_->NumMyPoints()];
  DoCopy(Source);
}

//=========================================================================
// Allocate
void EpetraExt_BlockDiagMatrix::Allocate(){

  int DataSize=0, NumBlocks=NumMyBlocks();
  Pivots_=new int[NumMyUnknowns()];
  int *ElementSize=new int[NumBlocks];
  
  for(int i=0;i<NumBlocks;i++) {
    ElementSize[i]=BlockSize(i)*BlockSize(i);
    DataSize+=ElementSize[i];
  }
  
  DataMap_=new Epetra_BlockMap(-1,Map().NumMyElements(),Map().MyGlobalElements(),ElementSize,0,Map().Comm());
  Values_=new double[DataSize];  
  delete [] ElementSize;
}


//=========================================================================
//! SetParameters
int EpetraExt_BlockDiagMatrix::SetParameters(Teuchos::ParameterList & List){
  List_=List;

  // Inverse storage mode
  string dummy("multiply");
  string ApplyMode=List_.get("apply mode",dummy);
  if(ApplyMode==string("multiply"))    ApplyMode_=AM_MULTIPLY;
  else if(ApplyMode==string("invert")) ApplyMode_=AM_INVERT;
  else if(ApplyMode==string("factor")) ApplyMode_=AM_FACTOR;
  else EPETRA_CHK_ERR(-1);

  return 0;
}


//=========================================================================
void EpetraExt_BlockDiagMatrix::PutScalar(double value){
  int MaxData=NumData();
  for (int i=0;i<MaxData;i++) Values_[i]=value;
}

//=========================================================================
// Assignment operator: Needs the same maps
EpetraExt_BlockDiagMatrix& EpetraExt_BlockDiagMatrix::operator=(const EpetraExt_BlockDiagMatrix& Source){
  DoCopy(Source);
  return(*this);
}
//=========================================================================
int EpetraExt_BlockDiagMatrix::DoCopy(const EpetraExt_BlockDiagMatrix& Source){
  // Need the same map
  if(!Map().SameAs(Source.Map()) || !DataMap_->SameAs(*Source.DataMap_))
    throw ReportError("Maps of DistBlockMatrices incompatible in assignment",-1);

  int MaxData=Source.NumData();

  for(int i=0;i<MaxData;i++)                 Values_[i]=Source.Values_[i];
  for(int i=0;i<Source.NumMyUnknowns();i++)  Pivots_[i]=Source.Pivots_[i];

  List_=Source.List_;
  ApplyMode_=Source.ApplyMode_;
  HasComputed_=Source.HasComputed_;

  return 0;
}


//=========================================================================
//! Computes the inverse / factorization if such is set on the list
int EpetraExt_BlockDiagMatrix::Compute(){
  int info;

  if(ApplyMode_==AM_MULTIPLY)
    // Multiply mode - noop    
    return 0;
  else {
    // Factorization - Needed for both 'factor' and 'invert' modes
    int NumBlocks=NumMyBlocks();
    for(int i=0;i<NumBlocks;i++){
      int Nb=BlockSize(i);
      if(Nb==1) {
        // Optimize for Block Size 1        
        Values_[DataMap_->FirstPointInElement(i)]=1.0/Values_[DataMap_->FirstPointInElement(i)];
      }
      else if(Nb==2) {
        // Optimize for Block Size 2
        double * v=&Values_[DataMap_->FirstPointInElement(i)];          
        double d=1/(v[0]*v[3]-v[1]*v[2]);
        double v0old=v[0];
        v[0]=v[3]*d;
        v[1]=-v[1]*d;
        v[2]=-v[2]*d;
        v[3]=v0old*d;
      }
      else{
        // "Large" Block - Use LAPACK
        LAPACK.GETRF(Nb,Nb,&Values_[DataMap_->FirstPointInElement(i)],Nb,&Pivots_[Map().FirstPointInElement(i)],&info);
        if(info) EPETRA_CHK_ERR(-2);
      }
    }
    
    if(ApplyMode_==AM_INVERT){
      // Invert - Use the factorization and invert the blocks in-place
      int lwork=3*DataMap_->MaxMyElementSize();
      std::vector<double> work(lwork);
      for(int i=0;i<NumBlocks;i++){
        int Nb=BlockSize(i);
        if(Nb==1 || Nb==2){
          // Optimize for Block Size 1 and 2
          // No need to do anything - factorization has already inverted the value
        }
        else{
          // "Large" Block - Use LAPACK
          LAPACK.GETRI(Nb,&Values_[DataMap_->FirstPointInElement(i)],Nb,&Pivots_[Map().FirstPointInElement(i)],&work[0],&lwork,&info);
          if(info) EPETRA_CHK_ERR(-3);
        }
      }
    }      
  }
  HasComputed_=true;
  return 0;
}


//=========================================================================
int EpetraExt_BlockDiagMatrix::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  int info;
  // Sanity Checks
  int NumVectors=X.NumVectors();  
  if(NumVectors!=Y.NumVectors())
    EPETRA_CHK_ERR(-1);
  if(!HasComputed_ && (ApplyMode_==AM_INVERT || ApplyMode_==AM_FACTOR))
    EPETRA_CHK_ERR(-2);
  
  //NTS: MultiVector's MyLength and [] Operators are  "points" level operators
  //not a "block/element" level operators.

  const int *vlist=DataMap_->FirstPointInElementList();
  const int *xlist=Map().FirstPointInElementList();
  const int *blocksize=Map().ElementSizeList();
  
  if(ApplyMode_==AM_MULTIPLY || ApplyMode_==AM_INVERT){
    // Multiply & Invert mode have the same apply
    int NumBlocks=NumMyBlocks();
    for(int i=0;i<NumBlocks;i++){
      int Nb=blocksize[i];
      int vidx0=vlist[i];
      int xidx0=xlist[i];
      for(int j=0;j<NumVectors;j++){	
	if(Nb==1) {
	  // Optimize for size = 1
	  Y[j][xidx0]=Values_[vidx0]*X[j][xidx0];
	}
	else if(Nb==2){
	  // Optimize for size = 2
	  Y[j][xidx0  ]=Values_[vidx0  ]*X[j][xidx0] + Values_[vidx0+2]*X[j][xidx0+1];
	  Y[j][xidx0+1]=Values_[vidx0+1]*X[j][xidx0] + Values_[vidx0+3]*X[j][xidx0+1];
	}
	else{
	  // "Large" Block - Use BLAS
	  //void 	GEMV (const char TRANS, const int M, const int N, const double ALPHA, const double *A, const int LDA, const double *X, const double BETA, double *Y, const int INCX=1, const int INCY=1) const 
	  GEMV('N',Nb,Nb,1.0,&Values_[vidx0],Nb,&X[j][xidx0],0.0,&Y[j][xidx0]);
	}
      }   
    }
  }
  else{
    // Factorization mode has a different apply
    int NumBlocks=NumMyBlocks();
    for(int i=0;i<NumBlocks;i++){
      int Nb=blocksize[i];
      int vidx0=vlist[i];
      int xidx0=xlist[i];      
      for(int j=0;j<NumVectors;j++){
	if(Nb==1) {
	  // Optimize for size = 1 - use the inverse
	  Y[j][xidx0]=Values_[vidx0]*X[j][xidx0];
	}
	else if(Nb==2){
	  // Optimize for size = 2 - use the inverse
	  Y[j][xidx0  ]=Values_[vidx0  ]*X[j][xidx0] + Values_[vidx0+2]*X[j][xidx0+1];
	  Y[j][xidx0+1]=Values_[vidx0+1]*X[j][xidx0] + Values_[vidx0+3]*X[j][xidx0+1];
	}
	else{
	  // "Large" Block - use LAPACK
	  //    void 	GETRS (const char TRANS, const int N, const int NRHS, const double *A, const int LDA, const int *IPIV, double *X, const int LDX, int *INFO) const 
	  for(int k=0;k<Nb;k++) Y[j][xidx0+k]=X[j][xidx0+k];
	  LAPACK.GETRS('N',Nb,1,&Values_[vidx0],Nb,&Pivots_[xidx0],&Y[j][xidx0],Nb,&info);
	  if(info) EPETRA_CHK_ERR(info);
	}
      }
    }    
  }  
  return 0;
}




//=========================================================================
// Print method
void EpetraExt_BlockDiagMatrix::Print(ostream & os) const{
  int MyPID = DataMap_->Comm().MyPID();
  int NumProc = DataMap_->Comm().NumProc();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumMyElements1 =DataMap_->NumMyElements();
      int MaxElementSize1 = DataMap_->MaxElementSize();
      int * MyGlobalElements1 = DataMap_->MyGlobalElements();
      int * FirstPointInElementList1;
      if (MaxElementSize1!=1) FirstPointInElementList1 = DataMap_->FirstPointInElementList();

      if (MyPID==0) {
	os.width(8);
	os <<  "     MyPID"; os << "    ";
	os.width(12);
	if (MaxElementSize1==1)
	  os <<  "GID  ";
	else
	  os <<  "     GID/Point";
	os.width(20);
	os <<  "Values ";
	os << endl;
      }
      for (int i=0; i < NumMyElements1; i++) {
	for (int ii=0; ii< DataMap_->ElementSize(i); ii++) {
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
// Allows the source and target (\e this) objects to be compared for compatibility, return nonzero if not.
int EpetraExt_BlockDiagMatrix::CheckSizes(const Epetra_SrcDistObject& Source){
  return &Map() == &Source.Map();
}


 //=========================================================================
// Perform ID copies and permutations that are on processor.
int EpetraExt_BlockDiagMatrix::CopyAndPermute(const Epetra_SrcDistObject& Source,
                                           int NumSameIDs, 
                                           int NumPermuteIDs,
                                           int * PermuteToLIDs,
                                           int * PermuteFromLIDs,
                                           const Epetra_OffsetIndex * Indexor){
  (void)Indexor;

  const EpetraExt_BlockDiagMatrix & A = dynamic_cast<const EpetraExt_BlockDiagMatrix &>(Source);

  double *From=A.Values();
  double *To = Values_;

  int * ToFirstPointInElementList = 0;
  int * FromFirstPointInElementList = 0;
  int * FromElementSizeList = 0;
  int MaxElementSize = DataMap().MaxElementSize();
  bool ConstantElementSize = DataMap().ConstantElementSize();

  if (!ConstantElementSize) {
    ToFirstPointInElementList =   DataMap().FirstPointInElementList();
    FromFirstPointInElementList = A.DataMap().FirstPointInElementList();
    FromElementSizeList = A.DataMap().ElementSizeList();
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
// Perform any packing or preparation required for call to DoTransfer().
int EpetraExt_BlockDiagMatrix::PackAndPrepare(const Epetra_SrcDistObject& Source,
                                           int NumExportIDs,
                                           int* ExportLIDs,
                                           int& LenExports,
                                           char*& Exports,
                                           int& SizeOfPacket,
                                           int* Sizes,
                                           bool & VarSizes,
                                           Epetra_Distributor& Distor){
  (void)Sizes;
  (void)VarSizes;
  (void)Distor;
  const EpetraExt_BlockDiagMatrix & A = dynamic_cast<const EpetraExt_BlockDiagMatrix &>(Source);

  int j, jj, k;

  double *From=A.Values();
  int MaxElementSize = DataMap().MaxElementSize();
  bool ConstantElementSize = DataMap().ConstantElementSize();

  int * FromFirstPointInElementList = 0;
  int * FromElementSizeList = 0;

  if (!ConstantElementSize) {
    FromFirstPointInElementList = A.DataMap().FirstPointInElementList();
    FromElementSizeList = A.DataMap().ElementSizeList();
  }

  SizeOfPacket = MaxElementSize * (int)sizeof(double); 

  if(NumExportIDs*SizeOfPacket>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = NumExportIDs*SizeOfPacket;
    Exports = new char[LenExports];
  }

  double * ptr;

  if (NumExportIDs>0) {
    ptr = (double *) Exports;
    
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
	ptr = (double *) Exports + j*SizeOfPacket;
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
// Perform any unpacking and combining after call to DoTransfer().
int EpetraExt_BlockDiagMatrix::UnpackAndCombine(const Epetra_SrcDistObject& Source, 
                                             int NumImportIDs,
                                             int* ImportLIDs, 
                                             int LenImports,
                                             char* Imports,
                                             int& SizeOfPacket, 
                                             Epetra_Distributor& Distor,
                                             Epetra_CombineMode CombineMode,
                                             const Epetra_OffsetIndex * Indexor){
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

  double * To = Values_;
  int MaxElementSize = DataMap().MaxElementSize();
  bool ConstantElementSize = DataMap().ConstantElementSize();

  int * ToFirstPointInElementList = 0;
  int * ToElementSizeList = 0;

  if (!ConstantElementSize) {
    ToFirstPointInElementList = DataMap().FirstPointInElementList();
    ToElementSizeList = DataMap().ElementSizeList();
  }
  
  double * ptr;
  // Unpack it...

  ptr = (double *) Imports;
    
  // Point entry case
  if (MaxElementSize==1) {
      
      if (CombineMode==Add)
	for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] += *ptr++; // Add to existing value
      else if(CombineMode==Insert)
	for (j=0; j<NumImportIDs; j++) To[ImportLIDs[j]] = *ptr++;
      else if(CombineMode==AbsMax)
        for (j=0; j<NumImportIDs; j++) {
	  To[ImportLIDs[j]] = EPETRA_MAX( To[ImportLIDs[j]],std::abs(*ptr));
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
	    To[jj+k] = EPETRA_MAX( To[jj+k], std::abs(*ptr));
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
	ptr = (double *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	  for (k=0; k<ElementSize; k++)
	    To[jj+k] += *ptr++; // Add to existing value
      }
    }
    else  if(CombineMode==Insert){
      for (j=0; j<NumImportIDs; j++) {
	ptr = (double *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	  for (k=0; k<ElementSize; k++)
	    To[jj+k] = *ptr++;
      }
    }
    else  if(CombineMode==AbsMax){
      for (j=0; j<NumImportIDs; j++) {
	ptr = (double *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	for (k=0; k<ElementSize; k++) {
	    To[jj+k] = EPETRA_MAX( To[jj+k], std::abs(*ptr));
	    ptr++;
	}
      }
    }
    // Note:  The following form of averaging is not a true average if more that one value is combined.
    //        This might be an issue in the future, but we leave this way for now.
    else if(CombineMode==Average) {
      for (j=0; j<NumImportIDs; j++) {
	ptr = (double *) Imports + j*SizeOfPacket;
	jj = ToFirstPointInElementList[ImportLIDs[j]];
	int ElementSize = ToElementSizeList[ImportLIDs[j]];
	  for (k=0; k<ElementSize; k++)
	    { To[jj+k] += *ptr++; To[jj+k] /= 2;}
      }
    }
  }
  
  return(0);
}
  
