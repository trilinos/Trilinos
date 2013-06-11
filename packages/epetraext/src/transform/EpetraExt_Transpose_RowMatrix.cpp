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

#include <EpetraExt_Transpose_RowMatrix.h>

#include <Epetra_Export.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>

#include <Teuchos_TimeMonitor.hpp>
#include <vector>

//#define ENABLE_TRANSPOSE_TIMINGS

namespace EpetraExt {

// Provide a "resize" operation for double*'s. 
inline void resize_doubles(int nold,int nnew,double*& d){
  if(nnew > nold){
    double *tmp = new double[nnew];
    for(int i=0; i<nold; i++)
      tmp[i]=d[i];
    delete [] d;
    d=tmp;
  }
}


RowMatrix_Transpose::
~RowMatrix_Transpose()
{
  if( TransposeMatrix_ ) delete TransposeMatrix_;

  if( !OrigMatrixIsCrsMatrix_ )
  {
    delete [] Indices_;
    delete [] Values_;
  }
}

//=========================================================================
  Epetra_CrsMatrix* EpetraExt::RowMatrix_Transpose::BuildTempTrans() 
{
#ifdef ENABLE_TRANSPOSE_TIMINGS
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor MM(myTime);
  Teuchos::RCP<Teuchos::Time> mtime;
  mtime=MM.getNewTimer("Transpose: BuildTempTrans 1");
  mtime->start();
#endif

  int i,j,err;
  const Epetra_RowMatrix & orig    = *origObj_;
  const Epetra_CrsMatrix * OrigCrsMatrix = dynamic_cast<const Epetra_CrsMatrix*>(&orig);
  const Epetra_Map & TransMap      = orig.RowMatrixColMap();
  int TransNnz                     = orig.NumMyNonzeros();
  int NumIndices;

  Epetra_CrsMatrix *TempTransA1 = new Epetra_CrsMatrix(Copy, TransMap,orig.RowMatrixRowMap(),0);
  Epetra_IntSerialDenseVector & TransRowptr = TempTransA1->ExpertExtractIndexOffset();
  Epetra_IntSerialDenseVector & TransColind = TempTransA1->ExpertExtractIndices();
  double *&                     TransVals   = TempTransA1->ExpertExtractValues();  

  TransRowptr.Resize(NumMyCols_+1);
  TransColind.Resize(TransNnz);
  resize_doubles(0,TransNnz,TransVals);
  std::vector<int> CurrentStart(NumMyCols_,0);

  // Count up nnz using the Rowptr to count the number of non-nonzeros.
  if (OrigMatrixIsCrsMatrix_)
  {
    const Epetra_CrsGraph & OrigGraph = OrigCrsMatrix->Graph(); // Get matrix graph

    for (i=0; i<NumMyRows_; i++)
    {
      err = OrigGraph.ExtractMyRowView(i, NumIndices, Indices_); // Get view of ith row
      if (err != 0) throw OrigGraph.ReportError("ExtractMyRowView failed",err);
      for (j=0; j<NumIndices; j++) ++CurrentStart[Indices_[j]];
    }
  }
  else // Original is not a CrsMatrix
  {
    MaxNumEntries_ = 0;
    MaxNumEntries_ = orig.MaxNumEntries();
    delete [] Indices_; delete [] Values_;
    Indices_ = new int[MaxNumEntries_];
    Values_  = new double[MaxNumEntries_];

    for (i=0; i<NumMyRows_; i++)
    {
      err = orig.ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_); 
      if (err != 0) {
        std::cerr << "ExtractMyRowCopy failed."<<std::endl;
        throw err;
      }
      for (j=0; j<NumIndices; j++) ++CurrentStart[Indices_[j]];
    }
  }

  // Scansum the rowptr; reset currentstart
  TransRowptr[0] = 0;
  for (i=1;i<NumMyCols_+1; i++)  TransRowptr[i]   = CurrentStart[i-1] + TransRowptr[i-1];
  for (i=0;i<NumMyCols_;   i++)  CurrentStart[i]  = TransRowptr[i];

  // Now copy values and global indices into newly create transpose storage
  for (i=0; i<NumMyRows_; i++)
  {
    if (OrigMatrixIsCrsMatrix_)
      err = OrigCrsMatrix->ExtractMyRowView(i, NumIndices, Values_, Indices_);
    else
      err = orig.ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_);
    if (err != 0) {
      std::cerr << "ExtractMyRowCopy failed."<<std::endl;
      throw err;
    }

    for (j=0; j<NumIndices; j++)
    {
      int idx = CurrentStart[Indices_[j]];
      TransColind[idx] = i;
      TransVals[idx]    = Values_[j];
      ++CurrentStart[Indices_[j]];// increment the counter into the local row
    }
  }

#ifdef ENABLE_TRANSPOSE_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("Transpose: BuildTempTrans 2");
  mtime->start();
#endif

  // Prebuild the importers and exporters the no-communication way, flipping the importers
  // and exporters around.
  Epetra_Import * myimport = 0;
  Epetra_Export * myexport = 0;
  if(OrigMatrixIsCrsMatrix_ && OrigCrsMatrix->Importer())
    myexport = new Epetra_Export(*OrigCrsMatrix->Importer());
  if(OrigMatrixIsCrsMatrix_ && OrigCrsMatrix->Exporter())
    myimport = new Epetra_Import(*OrigCrsMatrix->Exporter());

#ifdef ENABLE_TRANSPOSE_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("Transpose: BuildTempTrans 3");
  mtime->start();
#endif

  // Call ExpertStaticFillComplete
  err = TempTransA1->ExpertStaticFillComplete(orig.OperatorRangeMap(),*TransposeRowMap_,myimport,myexport);
  if (err != 0) {
    throw TempTransA1->ReportError("ExpertStaticFillComplete failed.",err);
  }

#ifdef ENABLE_TRANSPOSE_TIMINGS
  mtime->stop();
#endif
 
  return TempTransA1;
}


//=========================================================================
RowMatrix_Transpose::NewTypeRef
RowMatrix_Transpose::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  if( !TransposeRowMap_ )
  {
    if( IgnoreNonLocalCols_ )
      TransposeRowMap_ = (Epetra_Map *) &(orig.OperatorRangeMap()); // Should be replaced with refcount =
    else
      TransposeRowMap_ = (Epetra_Map *) &(orig.OperatorDomainMap()); // Should be replaced with refcount =
  }

  NumMyRows_ = orig.NumMyRows();
  NumMyCols_ = orig.NumMyCols();

  // This routine will work for any RowMatrix object, but will attempt cast the matrix to a CrsMatrix if
  // possible (because we can then use a View of the matrix and graph, which is much cheaper).

  // First get the local indices to count how many nonzeros will be in the 
  // transpose graph on each processor
  Epetra_CrsMatrix * TempTransA1 = BuildTempTrans();  

  if(!TempTransA1->Exporter()) {
    // Short Circuit: There is no need to make another matrix since TransposeRowMap_== TransMap
    newObj_ = TransposeMatrix_ = TempTransA1;
    return *newObj_;    
  }

#ifdef ENABLE_TRANSPOSE_TIMINGS
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor MM(myTime);
  Teuchos::RCP<Teuchos::Time> mtime;
  mtime=MM.getNewTimer("Transpose: Final FusedExport");
  mtime->start();
#endif

  // Now that transpose matrix with shared rows is entered, create a new matrix that will
  // get the transpose with uniquely owned rows (using the same row distribution as A).
  TransposeMatrix_ = new Epetra_CrsMatrix(*TempTransA1,*TempTransA1->Exporter(),0,TransposeRowMap_);

#ifdef ENABLE_TRANSPOSE_TIMINGS
  mtime->stop();
#endif

  newObj_ = TransposeMatrix_;
  delete TempTransA1;
  return *newObj_;
}

//=========================================================================
bool EpetraExt::RowMatrix_Transpose::fwd()
{
  Epetra_CrsMatrix * TempTransA1 = BuildTempTrans();  
  const Epetra_Export * TransposeExporter=0;
  bool DeleteExporter = false;
  
  if(TempTransA1->Exporter()) TransposeExporter = TempTransA1->Exporter();
  else {
    DeleteExporter=true;
    TransposeExporter = new Epetra_Export(TransposeMatrix_->DomainMap(),TransposeMatrix_->RowMap());
  }

  TransposeMatrix_->PutScalar(0.0);  // Zero out all values of the matrix

  EPETRA_CHK_ERR(TransposeMatrix_->Export(*TempTransA1, *TransposeExporter, Add));

  if(DeleteExporter) delete TransposeExporter;
  delete TempTransA1;
  return 0;
}

//=========================================================================
bool EpetraExt::RowMatrix_Transpose::rvs()
{
  EPETRA_CHK_ERR(-1); //Not Implemented Yet
  return false;
}

} // namespace EpetraExt
