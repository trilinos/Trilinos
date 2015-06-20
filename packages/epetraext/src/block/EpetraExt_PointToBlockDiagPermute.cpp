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

#include "Epetra_ConfigDefs.h"
#include "EpetraExt_PointToBlockDiagPermute.h"
#include "EpetraExt_BlockDiagMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Comm.h"

#include <stdio.h>
#include <fstream>

#define MAX(x,y) ((x)>(y)?(x):(y))

//=========================================================================
EpetraExt_PointToBlockDiagPermute::EpetraExt_PointToBlockDiagPermute(const Epetra_CrsMatrix& MAT)
  :Epetra_DistObject(MAT.RowMap()),
   Matrix_(&MAT),
   PurelyLocalMode_(true),
   ContiguousBlockMode_(false),
   ContiguousBlockSize_(0),
   NumBlocks_(0),
   Blockstart_(0),
   Blockids_int_(0),
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   Blockids_LL_(0),
#endif
   BDMap_(0),
   CompatibleMap_(0),
   BDMat_(0),
   Importer_(0),
   Exporter_(0),
   ImportVector_(0),
   ExportVector_(0)
{

}
  
//=========================================================================  
// Destructor
EpetraExt_PointToBlockDiagPermute::~EpetraExt_PointToBlockDiagPermute()
{
  if(BDMap_) delete BDMap_;
  if(CompatibleMap_) delete CompatibleMap_;
  if(BDMat_) delete BDMat_;
  if(Importer_) delete Importer_;
  if(Exporter_) delete Exporter_;
  if(ImportVector_) delete ImportVector_;
  if(ExportVector_) delete ExportVector_;
}

//=========================================================================  
// Set list
template<typename int_type>
int EpetraExt_PointToBlockDiagPermute::TSetParameters(Teuchos::ParameterList & List){
  List_=List;

  // Check for contiguous blocking first
  ContiguousBlockSize_=List_.get("contiguous block size",0);
  if(ContiguousBlockSize_!=0){
    ContiguousBlockMode_=true;
    PurelyLocalMode_=false;
  }  
  
  // Local vs. global ids & mode
  NumBlocks_=List_.get("number of local blocks",0);  
  Blockstart_=List_.get("block start index",(int*)0);    
  Blockids_ref<int>()=List_.get("block entry lids",(int*)0);
  if(Blockids_const_ptr<int>())
    PurelyLocalMode_=true;
  else{
    Blockids_ref<int_type>()=List_.get("block entry gids",(int_type*)0);
    if(Blockids_const_ptr<int_type>()) {
      PurelyLocalMode_=false;
      if(sizeof(int) != sizeof(int_type)) {
        Blockids_ref<int>() = new int[Matrix_->RowMap().NumMyElements()];
      }
    }
  }
  
  // Sanity checks
  if(ContiguousBlockMode_){
    // Can't use contiguous at the same time as the other modes
    if(NumBlocks_ || Blockstart_ || Blockids_const_ptr<int>() || Blockids_const_ptr<int_type>()) EPETRA_CHK_ERR(-4);
  }
  else {
    if(NumBlocks_ <= 0) EPETRA_CHK_ERR(-1);
    if(!Blockstart_) EPETRA_CHK_ERR(-2);
    if(!Blockids_const_ptr<int>() && !Blockids_const_ptr<int_type>()) EPETRA_CHK_ERR(-3);
  }
  
  return 0;
}

int EpetraExt_PointToBlockDiagPermute::SetParameters(Teuchos::ParameterList & List){
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Matrix_->RowMap().GlobalIndicesInt()) {
	return TSetParameters<int>(List);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Matrix_->RowMap().GlobalIndicesLongLong()) {
	return TSetParameters<long long>(List);
  }
  else
#endif
    throw "EpetraExt_PointToBlockDiagPermute::SetParameters: ERROR, GlobalIndices type unknown.";
}

//=========================================================================  
// Extracts the block-diagonal, builds maps, etc.
int EpetraExt_PointToBlockDiagPermute::Compute(){
  int rv=ExtractBlockDiagonal();
  return rv;
}

//=========================================================================  
// Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
int EpetraExt_PointToBlockDiagPermute::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  return -1;

}

//=========================================================================  
// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
int EpetraExt_PointToBlockDiagPermute::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  // Stuff borrowed from Epetra_CrsMatrix
  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) {
    EPETRA_CHK_ERR(-2); // Need same number of vectors in each MV
  }

  const Epetra_MultiVector *Xp=&X;
  Epetra_MultiVector *Yp=&Y;

  // Allocate temp workspace if X==Y and there are no imports or exports
  Epetra_MultiVector * Xcopy = 0;
  if (&X==&Y && Importer_==0 && Exporter_==0) {
    Xcopy = new Epetra_MultiVector(X);
    Xp=Xcopy;
  }
  
  UpdateImportVector(NumVectors); // Make sure Import and Export Vectors are compatible
  UpdateExportVector(NumVectors);

  // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
  if (Importer_){
    EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer_, Insert));
    Xp=ImportVector_;
  }
  
  // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
  if (Exporter_) {
    Yp=ExportVector_;
  }
  
  // Do the matvec 
  BDMat_->ApplyInverse(*Xp,*Yp);

  // Export if needed
  if (Exporter_) {
    Y.PutScalar(0.0);  // Make sure target is zero
    Y.Export(*ExportVector_, *Exporter_, Add); // Fill Y with Values from export vector
  }
  
  // Cleanup
  if(Xcopy) {
    delete Xcopy;
    EPETRA_CHK_ERR(1); // Return positive code to alert the user about needing extra copy of X
    return 1;
  }

  return 0;
}

//=========================================================================  
// Print method
void EpetraExt_PointToBlockDiagPermute::Print(std::ostream& os) const{
  if(Importer_) std::cout<<*Importer_<<std::endl;
  if(Exporter_) std::cout<<*Exporter_<<std::endl;
  if(BDMat_) std::cout<<*BDMat_<<std::endl;
}

//=========================================================================  
// Pulls the block diagonal of the matrix and then builds the BDMat_
template<typename int_type>
int EpetraExt_PointToBlockDiagPermute::TExtractBlockDiagonal(){  
  int i,j;
  std::vector<int_type> ExtRows;
  int* LocalColIDS=0;
  int ExtSize=0,ExtCols=0,MainCols=0;
  int Nrows=Matrix_->NumMyRows();
  int *l2b,*block_offset;
  int_type *l2blockid;
  const Epetra_Map &RowMap=Matrix_->RowMap();
  int index,col,row_in_block,col_in_block,length,*colind;
  double *values;

  bool verbose=(bool)(List_.get("output",0) > 0);
  
  // Contiguous Setup
  SetupContiguousMode();

  // Compute block size lists
  int *bsize=new int[NumBlocks_]; 
  for(i=0;i<NumBlocks_;i++) 
    bsize[i]=Blockstart_[i+1]-Blockstart_[i];
  
  // Use the ScanSum function to compute a prefix sum of the number of points
  int_type MyMinGID;
  int_type NumBlocks_int_type = NumBlocks_;
  Matrix_->Comm().ScanSum(&NumBlocks_int_type,&MyMinGID, 1);
  MyMinGID-=NumBlocks_;
  int_type *MyBlockGIDs=new int_type[NumBlocks_];
  for(i=0;i<NumBlocks_;i++)
    MyBlockGIDs[i]=MyMinGID+i;

  BDMap_=new Epetra_BlockMap((int_type) -1,NumBlocks_,MyBlockGIDs,bsize,0,Matrix_->Comm());
  
  // Allocate the Epetra_DistBlockMatrix
  BDMat_=new EpetraExt_BlockDiagMatrix(*BDMap_,true);  
  
  // First check to see if we can switch back to PurelyLocalMode_ if we're not
  // in it
  if(!PurelyLocalMode_){
    // Find the non-local row IDs
    ExtSize=MAX(0,Blockstart_[NumBlocks_]-RowMap.NumMyElements());// estimate
    ExtRows.reserve(ExtSize);
    
    for(i=0;i<Blockstart_[NumBlocks_];i++){
      if(RowMap.LID(Blockids_const_ptr<int_type>()[i])==-1)
        ExtRows.push_back(Blockids_const_ptr<int_type>()[i]);
    }

    // Switch to PurelyLocalMode_ if we never need GIDs
    int size_sum;
    ExtSize=ExtRows.size();   
    RowMap.Comm().SumAll(&ExtSize,&size_sum,1);
    if(size_sum==0){
      if(verbose && !Matrix_->Comm().MyPID()) printf("EpetraExt_PointToBlockDiagPermute: Switching to purely local mode\n");
      PurelyLocalMode_=true;
      for(i=0;i<Blockstart_[NumBlocks_];i++){
        Blockids_ref<int>()[i]=RowMap.LID(Blockids_const_ptr<int_type>()[i]);
      }     
    }
  }
  
  if(PurelyLocalMode_){
    /*******************************************************/
    // Allocations
    l2b=new int[Nrows];
    block_offset=new int[Nrows];

    // Build the local-id-to-block-id list, block offset list    
    for(i=0;i<NumBlocks_;i++) {
      for(j=Blockstart_[i];j<Blockstart_[i+1];j++){
        block_offset[Blockids_const_ptr<int>()[j]]=j-Blockstart_[i];
        l2b[Blockids_const_ptr<int>()[j]]=i;
      }
    }
    
    // Copy the data to the EpetraExt_BlockDiagMatrix        
    for(i=0;i<Nrows;i++) {
      int block_num = l2b[i];
      if(block_num>=0 && block_num<NumBlocks_) {
        row_in_block=block_offset[i];
        Matrix_->ExtractMyRowView(i,length,values,colind);
        int Nnz=0;
        for (j = 0; j < length; j++) {
          col = colind[j];
          if(col < Nrows && l2b[col]==block_num){
            Nnz++;
            col_in_block = block_offset[col];
            index = col_in_block * bsize[block_num] + row_in_block;
            (*BDMat_)[block_num][index]=values[j];
          }
        }
        /* Handle the case of a zero row. */
        /* By just putting a 1 on the diagonal. */
        if (Nnz == 0) {
          index = row_in_block * bsize[block_num] + row_in_block;
          (*BDMat_)[block_num][index]=1.0;
        }
      }
    }

    // Build the compatible map for import/export
    l2blockid=new int_type[Nrows];
    for(i=0;i<Nrows;i++) 
      l2blockid[Blockstart_[l2b[i]]+block_offset[i]]= (int_type) Matrix_->RowMap().GID64(i);

    // Build the Compatible Map, Import/Export Objects
    CompatibleMap_=new Epetra_Map(-1,Nrows,l2blockid,0,Matrix_->Comm());
    
  }
  else{
    /*******************************************************/      
    // Do the import to grab matrix entries
    // Allocate temporaries for import
    int_type* ExtRowsPtr = ExtRows.size()>0 ? &ExtRows[0] : NULL;
    Epetra_Map TmpMap((int_type) -1,ExtSize, ExtRowsPtr,0,Matrix_->Comm());; 
    Epetra_CrsMatrix TmpMatrix(Copy,TmpMap,0);
    Epetra_Import TmpImporter(TmpMap,RowMap);

    TmpMatrix.Import(*Matrix_,TmpImporter,Insert);
    TmpMatrix.FillComplete();
      
    ExtCols=TmpMatrix.NumMyCols();
    MainCols=Matrix_->NumMyCols();
        

    // Build the column reidex - main matrix
    LocalColIDS=new int[MainCols+ExtCols];      
    for(i=0;i<MainCols;i++){
      int_type GID= (int_type) Matrix_->GCID64(i);
      int MainLID=RowMap.LID(GID);
      if(MainLID!=-1) LocalColIDS[i]=MainLID;
      else{
        int ExtLID=TmpMatrix.LRID(GID);
        if(ExtLID!=-1) LocalColIDS[i]=Nrows+ExtLID;
        else LocalColIDS[i]=-1;
      }
    }

    // Build the column reidex - ext matrix
    for(i=0;i<ExtCols;i++){
      int_type GID= (int_type) TmpMatrix.GCID64(i);
      int MainLID=RowMap.LID(GID);
      if(MainLID!=-1) LocalColIDS[MainCols+i]=MainLID;
      else{
        int ExtLID=TmpMatrix.LRID(GID);          
        if(ExtLID!=-1) LocalColIDS[MainCols+i]=Nrows+ExtLID;
        else LocalColIDS[MainCols+i]=-1;
      }
    }

    // Allocations
    l2b=new int[Nrows+ExtSize];
    block_offset=new int[Nrows+ExtSize];  
  
    // Build l2b/block_offset with the expanded local index
    //NTS: Uses the same ordering of operation as the above routine, which is why
    //it works.
    for(i=0;i<Nrows+ExtSize;i++) block_offset[i]=l2b[i]=-1;    
    int ext_idx=0;
    for(i=0;i<NumBlocks_;i++) {
      for(j=Blockstart_[i];j<Blockstart_[i+1];j++){
        int LID=RowMap.LID(Blockids_const_ptr<int_type>()[j]);
        if(LID==-1) {LID=Nrows+ext_idx;ext_idx++;}
        block_offset[LID]=j-Blockstart_[i];
        l2b[LID]=i;          
      }
    }
    
    // Copy the data to the EpetraExt_BlockDiagMatrix from Matrix_
    for(i=0;i<Nrows;i++) {
      int block_num = l2b[i];
      if(block_num>=0 && block_num<NumBlocks_) {
        row_in_block=block_offset[i];
        Matrix_->ExtractMyRowView(i,length,values,colind);
        int Nnz=0;
        for (j = 0; j < length; j++) {
          col = LocalColIDS[colind[j]];
          if(col!=-1 && l2b[col]==block_num){
            Nnz++;
            col_in_block = block_offset[col];
            index = col_in_block * bsize[block_num] + row_in_block;
            (*BDMat_)[block_num][index]=values[j];
          }
        }
        /* Handle the case of a zero row. */
        /* By just putting a 1 on the diagonal. */
        if (Nnz == 0) {
          index = row_in_block * bsize[block_num] + row_in_block;
          (*BDMat_)[block_num][index]=values[j];
        }
      }
    }
    
    // Copy the data to the EpetraExt_BlockDiagMatrix from TmpMatrix
    for(i=0;i<ExtSize;i++) {
      int block_num = l2b[Nrows+i];
      if(block_num>=0 && block_num<NumBlocks_) {
        row_in_block=block_offset[Nrows+i];
        TmpMatrix.ExtractMyRowView(i,length,values,colind);
        int Nnz=0;
        for (j = 0; j < length; j++) {
          col = LocalColIDS[MainCols+colind[j]];
          if(col!=-1 && l2b[col]==block_num){
            Nnz++;
            col_in_block = block_offset[col];
            index = col_in_block * bsize[block_num] + row_in_block;
            (*BDMat_)[block_num][index]=values[j];
          }
        }
        /* Handle the case of a zero row. */
        /* By just putting a 1 on the diagonal. */
        if (Nnz == 0) {
          index = row_in_block * bsize[block_num] + row_in_block;
          (*BDMat_)[block_num][index]=1.0;
        }
      }      
    }
    
    // Build the compatible map for import/export
    l2blockid=new int_type[Blockstart_[NumBlocks_]];
    for(i=0;i<Nrows;i++){
      int bkid=l2b[i];
      if(bkid>-1) l2blockid[Blockstart_[bkid]+block_offset[i]]= (int_type) RowMap.GID64(i);
    }
    // NTS: This is easier - if we imported it, we know we need it
    for(i=0;i<ExtSize;i++)
      l2blockid[Blockstart_[l2b[Nrows+i]]+block_offset[Nrows+i]]= (int_type) TmpMatrix.GRID64(i);
    
    // Build the Compatible Map, Import/Export Objects
    CompatibleMap_=new Epetra_Map(-1,Blockstart_[NumBlocks_],l2blockid,0,Matrix_->Comm());
        
  }//end else

  // Set BDMat_'s parameter list and compute!
  Teuchos::ParameterList dummy,inList;
  inList=List_.get("blockdiagmatrix: list",dummy);
  BDMat_->SetParameters(inList);  
  BDMat_->Compute();
  
  // Build importer/exporter
  if(!CompatibleMap_->SameAs(Matrix_->DomainMap()))
    Importer_ = new Epetra_Import(*CompatibleMap_,Matrix_->DomainMap());
  if(!CompatibleMap_->SameAs(Matrix_->RangeMap()))
    Exporter_ = new Epetra_Export(*CompatibleMap_,Matrix_->RangeMap());

  // Cleanup
  delete [] LocalColIDS;
  delete [] block_offset;
  delete [] l2b;
  delete [] l2blockid;
  delete [] bsize;
  delete [] MyBlockGIDs;

  // Contiguous Cleanup
  CleanupContiguousMode();

  return 0;
}

int EpetraExt_PointToBlockDiagPermute::ExtractBlockDiagonal(){
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Matrix_->RowMap().GlobalIndicesInt()) {
	return TExtractBlockDiagonal<int>();
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Matrix_->RowMap().GlobalIndicesLongLong()) {
	return TExtractBlockDiagonal<long long>();
  }
  else
#endif
    throw "EpetraExt_PointToBlockDiagPermute::ExtractBlockDiagonal: ERROR, GlobalIndices type unknown.";
}

//=======================================================================================================
template<typename int_type>
int EpetraExt_PointToBlockDiagPermute::TSetupContiguousMode(){
  if(!ContiguousBlockMode_) return 0;
  // NTS: In case of processor-crossing blocks, the lowest PID always gets the block;
  const Epetra_Map &RowMap=Matrix_->RowMap();
  
  int_type MinMyGID= (int_type) RowMap.MinMyGID64(); 
  int_type MaxMyGID= (int_type) RowMap.MaxMyGID64(); 
  int_type Base= (int_type) Matrix_->IndexBase64();

  // Find the GID that begins my first block
  int_type MyFirstBlockGID=ContiguousBlockSize_*(int_type)ceil(((double)(MinMyGID - Base))/ContiguousBlockSize_)+Base;
  NumBlocks_=(int)ceil((double)((MaxMyGID-MyFirstBlockGID+1.0)) / ContiguousBlockSize_);

  // Allocate memory
  Blockstart_=new int[NumBlocks_+1];
  Blockids_ref<int_type>()=new int_type[NumBlocks_*ContiguousBlockSize_];
  if(sizeof(int) != sizeof(int_type))
    Blockids_ref<int>()=new int[NumBlocks_*ContiguousBlockSize_];
  Blockstart_[NumBlocks_]=NumBlocks_*ContiguousBlockSize_;

  // Fill the arrays
  for(int i=0,ct=0;i<NumBlocks_;i++){
    Blockstart_[i]=ct;
    for(int j=0;j<ContiguousBlockSize_;j++,ct++){
      Blockids_ref<int_type>()[ct]=MyFirstBlockGID+ct;
    }
  }
  
  return 0;
}

int EpetraExt_PointToBlockDiagPermute::SetupContiguousMode(){
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Matrix_->RowMap().GlobalIndicesInt()) {
	return TSetupContiguousMode<int>();
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Matrix_->RowMap().GlobalIndicesLongLong()) {
	return TSetupContiguousMode<long long>();
  }
  else
#endif
    throw "EpetraExt_PointToBlockDiagPermute::SetupContiguousMode: ERROR, GlobalIndices type unknown.";
}

//=======================================================================================================
int EpetraExt_PointToBlockDiagPermute::CleanupContiguousMode(){
  if(!ContiguousBlockMode_) return 0;
  NumBlocks_=0;
  if(Blockstart_) {delete [] Blockstart_; Blockstart_=0;}
  if(Blockids_int_)   {delete [] Blockids_int_; Blockids_int_=0;}
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Blockids_LL_)   {delete [] Blockids_LL_; Blockids_LL_=0;}
#endif
  return 0;
}


//=======================================================================================================
// Creates an Epetra_CrsMatrix from the BlockDiagMatrix.  This is generally only useful if you want to do a matrix-matrix multiply.
template<typename int_type>
Epetra_FECrsMatrix * EpetraExt_PointToBlockDiagPermute::TCreateFECrsMatrix(){
  Epetra_FECrsMatrix * NewMat=new Epetra_FECrsMatrix(Copy,Matrix_->RowMap(),0);
  
  const Epetra_BlockMap &BlockMap=BDMat_->BlockMap();
  const Epetra_BlockMap &DataMap=BDMat_->DataMap();
  const int *vlist=DataMap.FirstPointInElementList();
  const int *xlist=BlockMap.FirstPointInElementList();
  const int *blocksize=BlockMap.ElementSizeList();
  const double *values=BDMat_->Values();
  int NumBlocks=BDMat_->NumMyBlocks();

  // Maximum size vector for stashing GIDs
  std::vector<int_type> GIDs;
  GIDs.resize(BlockMap.MaxMyElementSize());


  for(int i=0;i<NumBlocks;i++){
    int Nb=blocksize[i];
    int vidx0=vlist[i];
    int xidx0=xlist[i];
    // Get global indices
    for(int j=0;j<Nb;j++)
      GIDs[j]= (int_type) CompatibleMap_->GID64(xidx0+j);
    
    // Remember: We're using column-major storage for LAPACK's benefit    
    int ierr=NewMat->InsertGlobalValues(Nb,&GIDs[0],&values[vidx0],Epetra_FECrsMatrix::COLUMN_MAJOR);
    if(ierr < 0) throw "CreateFECrsMatrix: ERROR in InsertGlobalValues";
  }   
  NewMat->GlobalAssemble();
  return NewMat;
}

Epetra_FECrsMatrix * EpetraExt_PointToBlockDiagPermute::CreateFECrsMatrix(){
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Matrix_->RowMap().GlobalIndicesInt()) {
	return TCreateFECrsMatrix<int>();
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Matrix_->RowMap().GlobalIndicesLongLong()) {
	return TCreateFECrsMatrix<long long>();
  }
  else
#endif
    throw "EpetraExt_PointToBlockDiagPermute::CreateFECrsMatrix: ERROR, GlobalIndices type unknown.";
}

//=======================================================================================================
void EpetraExt_PointToBlockDiagPermute::UpdateImportVector(int NumVectors) const {    
  if(Importer_ != 0) {
    if(ImportVector_ != 0) {
      if(ImportVector_->NumVectors() != NumVectors) { 
	delete ImportVector_; 
	ImportVector_= 0;
      }
    }
    if(ImportVector_ == 0) 
      ImportVector_ = new Epetra_MultiVector(*CompatibleMap_,NumVectors); // Create import vector if needed
  }
  return;
}
//=======================================================================================================
void EpetraExt_PointToBlockDiagPermute::UpdateExportVector(int NumVectors) const {    
  if(Exporter_ != 0) {
    if(ExportVector_ != 0) {
      if(ExportVector_->NumVectors() != NumVectors) { 
	delete ExportVector_; 
	ExportVector_= 0;
      }
    }
    if(ExportVector_ == 0) 
      ExportVector_ = new Epetra_MultiVector(*CompatibleMap_,NumVectors); // Create Export vector if needed
  }
  return;
}




//=========================================================================  
int EpetraExt_PointToBlockDiagPermute::Import(const Epetra_SrcDistObject& A, const Epetra_Import& Importer, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex *Indexor){
 return -1;
}

//=========================================================================  
int EpetraExt_PointToBlockDiagPermute::Import(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex *Indexor){
 return -1;
}

//=========================================================================  
int EpetraExt_PointToBlockDiagPermute::Export(const Epetra_SrcDistObject& A, const Epetra_Import & Importer, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex *Indexor){
 return -1;
}

//=========================================================================  
int EpetraExt_PointToBlockDiagPermute::Export(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex *Indexor){
 return -1;
}

//=========================================================================  
// Allows the source and target (\e this) objects to be compared for compatibility, return nonzero if not.
int EpetraExt_PointToBlockDiagPermute::CheckSizes(const Epetra_SrcDistObject& Source){
  return -1;
}

//=========================================================================  
// Perform ID copies and permutations that are on processor.
int EpetraExt_PointToBlockDiagPermute::CopyAndPermute(const Epetra_SrcDistObject& Source,
                   int NumSameIDs, 
                   int NumPermuteIDs,
                   int * PermuteToLIDs,
                   int * PermuteFromLIDs,
                   const Epetra_OffsetIndex * Indexor,
                   Epetra_CombineMode CombineMode){
  return -1;
}

//=========================================================================  
// Perform any packing or preparation required for call to DoTransfer().
int EpetraExt_PointToBlockDiagPermute::PackAndPrepare(const Epetra_SrcDistObject& Source,
                   int NumExportIDs,
                   int* ExportLIDs,
                   int& LenExports,
                   char*& Exports,
                   int& SizeOfPacket,
                   int* Sizes,
                   bool & VarSizes,
                   Epetra_Distributor& Distor){
  return -1;
}

//=========================================================================  
// Perform any unpacking and combining after call to DoTransfer().
int EpetraExt_PointToBlockDiagPermute::UnpackAndCombine(const Epetra_SrcDistObject& Source, 
                     int NumImportIDs,
                     int* ImportLIDs, 
                     int LenImports,
                     char* Imports,
                     int& SizeOfPacket, 
                     Epetra_Distributor& Distor,
                     Epetra_CombineMode CombineMode,
                     const Epetra_OffsetIndex * Indexor){
  return -1;
}


