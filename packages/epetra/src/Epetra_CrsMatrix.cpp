/*
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
*/

#include "Epetra_ConfigDefs.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_BLAS_wrappers.h"

#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_CrsGraphData.h"
#include "Epetra_HashTable.h"
#include "Epetra_Util.h"

#include <cstdlib>

#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
# include "Teuchos_VerboseObject.hpp"
bool Epetra_CrsMatrixTraceDumpMultiply = false;
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY

#ifdef HAVE_EPETRA_TEUCHOS
// Define this macro to see some timers for some of these functions
# define EPETRA_CRSMATRIX_TEUCHOS_TIMERS
#endif

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
#  include "Teuchos_TimeMonitor.hpp"
#endif

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& rowMap, const int* NumEntriesPerRow, bool StaticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(CV, rowMap, NumEntriesPerRow, StaticProfile),
    Allocated_(false),
    StaticGraph_(false),
    UseTranspose_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(rowMap.NumMyPoints()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Allocate();
}

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& rowMap, int NumEntriesPerRow, bool StaticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(CV, rowMap, NumEntriesPerRow, StaticProfile),
    Allocated_(false),
    StaticGraph_(false),
    UseTranspose_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(rowMap.NumMyPoints()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Allocate();
}
//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& rowMap, 
           const Epetra_Map& colMap, const int* NumEntriesPerRow, bool StaticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(CV, rowMap, colMap, NumEntriesPerRow, StaticProfile),
    Allocated_(false),
    StaticGraph_(false),
    UseTranspose_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(rowMap.NumMyPoints()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  Allocate();
}

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& rowMap, 
           const Epetra_Map& colMap, int NumEntriesPerRow, bool StaticProfile) 
  : Epetra_DistObject(rowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(CV, rowMap, colMap,  NumEntriesPerRow, StaticProfile),
    Allocated_(false),
    StaticGraph_(false),
    UseTranspose_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(rowMap.NumMyPoints()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  if(!rowMap.GlobalIndicesTypeMatch(colMap))
    throw ReportError("Epetra_CrsMatrix::Epetra_CrsMatrix: cannot be called with different indices types for rowMap and colMap", -1);

  InitializeDefaults();
  Allocate();
}
//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph& graph) 
  : Epetra_DistObject(graph.Map(), "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(graph),
    Allocated_(false),
    StaticGraph_(true),
    UseTranspose_(false),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(graph.NumMyRows()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(CV),
    squareFillCompleteCalled_(false)
{
  constructedWithFilledGraph_ = graph.Filled();
  InitializeDefaults();
  Allocate();
}

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(const Epetra_CrsMatrix& Matrix) 
  : Epetra_DistObject(Matrix),
    Epetra_CompObject(Matrix),
    Epetra_BLAS(),
    Graph_(Matrix.Graph()),
    Allocated_(false),
    StaticGraph_(true),
    UseTranspose_(Matrix.UseTranspose_),
    constructedWithFilledGraph_(false),
    matrixFillCompleteCalled_(false),
    StorageOptimized_(false),
    Values_(0),
    Values_alloc_lengths_(0),
    All_Values_(0),
    NormInf_(0.0),
    NormOne_(0.0),
    NormFrob_(0.0),
    NumMyRows_(Matrix.NumMyRows()),
    ImportVector_(0),
    ExportVector_(0),
    CV_(Copy),
    squareFillCompleteCalled_(false)
{
  InitializeDefaults();
  operator=(Matrix);
}

//==============================================================================
Epetra_CrsMatrix& Epetra_CrsMatrix::operator=(const Epetra_CrsMatrix& src)
{
  if (this == &src) {
    return( *this );
  }

  if (!src.Filled()) throw ReportError("Copying an Epetra_CrsMatrix requires source matrix to have Filled()==true", -1);

  Graph_ = src.Graph_; // Copy graph

  DeleteMemory();

  StaticGraph_ = true;
  UseTranspose_ = src.UseTranspose_;
  constructedWithFilledGraph_ = src.constructedWithFilledGraph_;
  matrixFillCompleteCalled_ = src.matrixFillCompleteCalled_;
  Values_ = 0;
  Values_alloc_lengths_ = 0;
  All_Values_ = 0;
  NormInf_ = -1.0;
  NormOne_ = -1.0;
  NormFrob_ = -1.0;
  NumMyRows_ = src.NumMyRows_;
  ImportVector_ = 0;
  ExportVector_ = 0;

  CV_ = Copy;

  StorageOptimized_ = src.StorageOptimized_;
  if (src.StorageOptimized()) { // Special copy for case where storage is optimized

    int numMyNonzeros = Graph().NumMyEntries();
    if (numMyNonzeros>0) All_Values_ = new double[numMyNonzeros];
    double * srcValues = src.All_Values();
    for (int i=0; i<numMyNonzeros; ++i) All_Values_[i] = srcValues[i];
    Allocated_ = true;
#ifdef Epetra_ENABLE_CASK
    if( matrixFillCompleteCalled_  ) {
      cask = cask_handler_copy(src.cask);
    }
#endif
  }
  else { // copy for non-optimized storage
    
    Allocate();
    for (int i=0; i<NumMyRows_; i++) {
      int NumEntries = src.NumMyEntries(i);
      double * const srcValues = src.Values(i);
      double * targValues = Values(i);
      for (int j=0; j< NumEntries; j++) targValues[j] = srcValues[j];
    }
  }

  return( *this );
}

//==============================================================================
void Epetra_CrsMatrix::InitializeDefaults() { // Initialize all attributes that have trivial default values

  UseTranspose_ = false;
  Values_ = 0;
  Values_alloc_lengths_ = 0;
  All_Values_ = 0;
  NormInf_ = -1.0;
  NormOne_ = -1.0;
  NormFrob_ = -1.0;
  ImportVector_ = 0;
  ExportVector_ = 0;

  return;
}

//==============================================================================
int Epetra_CrsMatrix::Allocate() {

  int i, j;

  // Allocate Values array
  Values_ = NumMyRows_ > 0 ? new double*[NumMyRows_] : NULL;
  Values_alloc_lengths_ = NumMyRows_ > 0 ? new int[NumMyRows_] : NULL;
  if (NumMyRows_ > 0) {
    for(j=0; j<NumMyRows_; ++j) Values_alloc_lengths_[j] = 0;
  }

  // Allocate and initialize entries if we are copying data
  if (CV_==Copy) {
    if (Graph().StaticProfile() || Graph().StorageOptimized()) {
      int numMyNonzeros = Graph().NumMyEntries();
      if (numMyNonzeros>0) All_Values_ = new double[numMyNonzeros];
      if(Graph().StorageOptimized()){
        StorageOptimized_ = true;
      }
    }
    double * all_values = All_Values_;
    for (i=0; i<NumMyRows_; i++) {
      int NumAllocatedEntries = Graph().NumAllocatedMyIndices(i);

      if (NumAllocatedEntries > 0) {
        if (Graph().StaticProfile() || Graph().StorageOptimized()) {
          Values_[i] = all_values;
          all_values += NumAllocatedEntries;
        }
        else {
          Values_[i] = new double[NumAllocatedEntries];
          Values_alloc_lengths_[i] = NumAllocatedEntries;
        }
      }
      else 
        Values_[i] = 0;

      for(j=0; j< NumAllocatedEntries; j++) 
        Values_[i][j] = 0.0; // Fill values with zero
    }
  }   
  else {
    for (i=0; i<NumMyRows_; i++) {
      Values_[i] = 0;
    }
  }
  SetAllocated(true);
#ifdef Epetra_ENABLE_CASK
  cask=NULL;
#endif
  return(0);
}
//==============================================================================
Epetra_CrsMatrix::~Epetra_CrsMatrix()
{
  DeleteMemory();
}

//==============================================================================
void Epetra_CrsMatrix::DeleteMemory()
{
  int i;

  if (CV_==Copy) {
    if (All_Values_!=0)
      delete [] All_Values_;
    else if (Values_!=0)
      for (i=0; i<NumMyRows_; i++) 
        if (Graph().NumAllocatedMyIndices(i) >0) 
          delete [] Values_[i];
  }

  if (ImportVector_!=0) 
    delete ImportVector_;
  ImportVector_=0;

  if (ExportVector_!=0)
    delete ExportVector_;
  ExportVector_=0;

  delete [] Values_;
  Values_ = 0;

  delete [] Values_alloc_lengths_;
  Values_alloc_lengths_ = 0;

  NumMyRows_ = 0;

#ifdef Epetra_ENABLE_CASK
  if( StorageOptimized_  )
  {
    if( cask != NULL ) 
      cask_handler_destroy(cask);

    cask = NULL;
  }
#endif

  Allocated_ = false;
}

//==============================================================================
int Epetra_CrsMatrix::ReplaceRowMap(const Epetra_BlockMap& newmap)
{
  int err = Graph_.ReplaceRowMap(newmap);
  if (err == 0) {
    //update export vector.

    if (ExportVector_ != 0) {
      delete ExportVector_; 
      ExportVector_= 0;
    }

    ExportVector_ = new Epetra_MultiVector(RowMap(),1);
  }
  return(err);
}

//==============================================================================
int Epetra_CrsMatrix::ReplaceColMap(const Epetra_BlockMap& newmap)
{
  int err = Graph_.ReplaceColMap(newmap);
  if (err == 0) {
    //update import vector.

    if (ImportVector_ != 0) {
      delete ImportVector_;
      ImportVector_= 0;
    }

    ImportVector_ = new Epetra_MultiVector(ColMap(),1);
  }
  return(err);
}

//==============================================================================
int Epetra_CrsMatrix::ReplaceDomainMapAndImporter(const Epetra_Map& NewDomainMap, const Epetra_Import * NewImporter) {
  return Graph_.ReplaceDomainMapAndImporter(NewDomainMap,NewImporter);
}

//==============================================================================
int Epetra_CrsMatrix::PutScalar(double ScalarConstant) 
{
  if (StorageOptimized()) {
    int length = NumMyNonzeros();
    for (int i=0; i<length; ++i) All_Values_[i] = ScalarConstant;
  }
  else {
    for(int i=0; i<NumMyRows_; i++) {
      int NumEntries = Graph().NumMyIndices(i);
      double * targValues = Values(i);
      for(int j=0; j< NumEntries; j++) 
  targValues[j] = ScalarConstant;
    }
  }
  return(0);
}
//==============================================================================
int Epetra_CrsMatrix::Scale(double ScalarConstant) 
{
  if (StorageOptimized()) {
    int length = NumMyNonzeros();
    for (int i=0; i<length; ++i) All_Values_[i] *= ScalarConstant;
  }
  else {
    for(int i=0; i<NumMyRows_; i++) {
      int NumEntries = Graph().NumMyIndices(i);
      double * targValues = Values(i);
      for(int j=0; j< NumEntries; j++) 
  targValues[j] *= ScalarConstant;
    }
  }
  return(0);
}

//==========================================================================
template<typename int_type>
int Epetra_CrsMatrix::TInsertGlobalValues(int_type Row, int NumEntries,
           const double* values,
           const int_type* Indices)
{
  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into local graph
  if(IndicesAreContiguous()) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed
  Graph_.SetIndicesAreGlobal(true);
  int locRow = Graph_.LRID(Row); // Find local row number for this global row index

  EPETRA_CHK_ERR( InsertValues(locRow, NumEntries, values, Indices) );

  return(0);
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::InsertGlobalValues(int Row, int NumEntries,
           const double* values,
           const int* Indices)
{
  if(RowMap().GlobalIndicesInt())
    return TInsertGlobalValues<int>(Row, NumEntries, values, Indices);
  else
    throw ReportError("Epetra_CrsMatrix::InsertGlobalValues int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::InsertGlobalValues(long long Row, int NumEntries,
           const double* values,
           const long long* Indices)
{
  if(RowMap().GlobalIndicesLongLong())
    return TInsertGlobalValues<long long>(Row, NumEntries, values, Indices);
  else
    throw ReportError("Epetra_CrsMatrix::InsertGlobalValues long long version called for a matrix that is not long long.", -1);
}
#endif

//==========================================================================
template<typename int_type>
int Epetra_CrsMatrix::TInsertGlobalValues(int_type Row, int NumEntries,
           double* values,
           int_type* Indices)
{
  if(IndicesAreLocal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into local graph
  if(IndicesAreContiguous()) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed
  Graph_.SetIndicesAreGlobal(true);
  int locRow = Graph_.LRID(Row); // Find local row number for this global row index

  EPETRA_CHK_ERR( InsertValues(locRow, NumEntries, values, Indices) );

  return(0);
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::InsertGlobalValues(int Row, int NumEntries,
           double* values,
           int* Indices)
{
  if(RowMap().GlobalIndicesInt())
    return TInsertGlobalValues<int>(Row, NumEntries, values, Indices);
  else 
    throw ReportError("Epetra_CrsMatrix::InsertGlobalValues int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::InsertGlobalValues(long long Row, int NumEntries,
           double* values,
           long long* Indices)
{
  if(RowMap().GlobalIndicesLongLong()) {
    return TInsertGlobalValues<long long>(Row, NumEntries, values, Indices);
  }
  else
    throw ReportError("Epetra_CrsMatrix::InsertGlobalValues long long version called for a matrix that is not long long.", -1);
}
#endif
//==========================================================================
int Epetra_CrsMatrix::InsertMyValues(int Row, int NumEntries,
             const double* values,
             const int* Indices)
{
  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into filled graph
  if(IndicesAreContiguous() && CV_==Copy) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and new
  Graph_.SetIndicesAreLocal(true);

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  EPETRA_CHK_ERR( InsertValues(Row, NumEntries, values, Indices) );
#else
  throw ReportError("Epetra_CrsMatrix::InsertMyValues: Failure because neither 32 bit nor 64 bit indices insertable.", -1);
#endif

  return(0);

}

//==========================================================================
int Epetra_CrsMatrix::InsertMyValues(int Row, int NumEntries,
             double* values,
             int* Indices)
{
  if(IndicesAreGlobal()) 
    EPETRA_CHK_ERR(-2); // Cannot insert global values into filled graph
  if(IndicesAreContiguous() && CV_==Copy) 
    EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and new
  Graph_.SetIndicesAreLocal(true);

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  EPETRA_CHK_ERR( InsertValues(Row, NumEntries, values, Indices) );
#else
    throw ReportError("Epetra_CrsMatrix::InsertMyValues: Failure because neither 32 bit nor 64 bit indices insertable.", -1);
#endif

  return(0);

}

//==========================================================================
template<typename int_type>
int Epetra_CrsMatrix::InsertValues(int Row, int NumEntries,
           const double* values,
           const int_type* Indices)
{
  if(CV_ == View){
    //cannot allow View mode with const pointers
    EPETRA_CHK_ERR(-4);
  }
  else{
    //to avoid code duplication I am using a cheap tactic of removing the constness of
    //values and Indices. Since this is only called in copy mode the passed in variables
    //will not be modified.
    return(InsertValues(Row, NumEntries, const_cast<double*>(values), const_cast<int_type*>(Indices)));
  }
  return 0;
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::InsertValues(int Row, int NumEntries,
           const double* values,
           const int* Indices)
{
  if(RowMap().GlobalIndicesInt() || (RowMap().GlobalIndicesLongLong() && IndicesAreLocal()))
    return InsertValues<int>(Row, NumEntries, values, Indices);
  else
    throw ReportError("Epetra_CrsMatrix::InsertValues int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::InsertValues(int Row, int NumEntries,
           const double* values,
           const long long* Indices)
{
  if(RowMap().GlobalIndicesLongLong())
    return InsertValues<long long>(Row, NumEntries, values, Indices);
  else
    throw ReportError("Epetra_CrsMatrix::InsertValues long long version called for a matrix that is not long long.", -1);
}
#endif
//==========================================================================
template<typename int_type>
int Epetra_CrsMatrix::InsertValues(int Row, int NumEntries,
           double* values,
           int_type* Indices)
{
  int j;
  int ierr = 0;

  if(Row < 0 || Row >= NumMyRows_) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  if(CV_ == View) {
    //test indices in static graph
    if(StaticGraph()) {
      int testNumEntries;
      int* testIndices;
      int testRow = Row;
      if(IndicesAreGlobal()) 
        testRow = Graph_.LRID( Row );
      EPETRA_CHK_ERR(Graph_.ExtractMyRowView(testRow, testNumEntries, testIndices));

      bool match = true;
      if(NumEntries != testNumEntries) 
        match = false;
      for(int i = 0; i < NumEntries; ++i)
        match = match && (Indices[i]==testIndices[i]);

      if(!match)
        ierr = -3;
    }

    if(Values_[Row] != 0) 
      ierr = 2; // This row has been defined already.  Issue warning.
    Values_[Row] = values;
  }
  else {    
    if(StaticGraph()) 
      EPETRA_CHK_ERR(-2); // If the matrix graph is fully constructed, we cannot insert new values

    int tmpNumEntries = NumEntries;

    if(Graph_.HaveColMap()) { //must insert only valid indices, values
      const double* tmpValues = values;
      values = new double[NumEntries];
      int loc = 0;
      if(IndicesAreLocal()) {
        for(int i = 0; i < NumEntries; ++i)
          if(Graph_.ColMap().MyLID(static_cast<int>(Indices[i])))
            values[loc++] = tmpValues[i];
      }
      else {
        for(int i = 0; i < NumEntries; ++i)
          if(Graph_.ColMap().MyGID(Indices[i])) 
            values[loc++] = tmpValues[i];
      }
      if(NumEntries != loc) 
        ierr = 2; //Some columns excluded
      NumEntries = loc;
    } 

    int start = Graph().NumMyIndices(Row);
    int stop = start + NumEntries;
    int NumAllocatedEntries = Values_alloc_lengths_[Row];
    if(stop > NumAllocatedEntries) {
      if (Graph().StaticProfile() && stop > Graph().NumAllocatedMyIndices(Row)) {
        EPETRA_CHK_ERR(-2); // Cannot expand graph storage if graph created using StaticProfile
      }
      if(NumAllocatedEntries == 0) {
        Values_[Row] = new double[NumEntries]; // Row was never allocated, so do it
        Values_alloc_lengths_[Row] = NumEntries;
      }
      else {
        ierr = 1; // Out of room.  Must delete and allocate more space...
        double* tmp_Values = new double[stop];
        for(j = 0; j < start; j++) 
          tmp_Values[j] = Values_[Row][j]; // Copy existing entries
        delete[] Values_[Row]; // Delete old storage
        Values_[Row] = tmp_Values; // Set pointer to new storage
        Values_alloc_lengths_[Row] = stop;
      }
    }

    for(j = start; j < stop; j++) 
      Values_[Row][j] = values[j-start];


    NumEntries = tmpNumEntries;
    if(Graph_.HaveColMap()) 
      delete[] values;
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  if(!StaticGraph()) {
    EPETRA_CHK_ERR(Graph_.InsertIndices(Row, NumEntries, Indices));
  }

  EPETRA_CHK_ERR(ierr);
  return(0);
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::InsertValues(int Row, int NumEntries,
           double* values,
           int* Indices)
{
  if(RowMap().GlobalIndicesInt() || (RowMap().GlobalIndicesLongLong() && IndicesAreLocal()))
  return InsertValues<int>(Row, NumEntries, values, Indices);
  else
  throw ReportError("Epetra_CrsMatrix::InsertValues int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::InsertValues(int Row, int NumEntries,
           double* values,
           long long* Indices)
{
  if(RowMap().GlobalIndicesLongLong())
  return InsertValues<long long>(Row, NumEntries, values, Indices);
  else
  throw ReportError("Epetra_CrsMatrix::InsertValues long long version called for a matrix that is not long long.", -1);
}
#endif

//==========================================================================
int Epetra_CrsMatrix::InsertOffsetValues(long long Row, int NumEntries,
           double* values,
           int* Indices)
{
  return ReplaceOffsetValues(Row, NumEntries, values, Indices);
}

//==========================================================================
template<typename int_type>
int Epetra_CrsMatrix::TReplaceGlobalValues(int_type Row, int NumEntries, const double * srcValues, const int_type *Indices) {

  int j;
  int ierr = 0;
  int Loc;

  int locRow = Graph_.LRID(Row); // Normalize row range
    
  if (locRow < 0 || locRow >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }
  double * targValues = Values(locRow);
  for (j=0; j<NumEntries; j++) {
    int_type Index = Indices[j];
    if (Graph_.FindGlobalIndexLoc(locRow,Index,j,Loc)) 
      targValues[Loc] = srcValues[j];
    else 
      ierr = 2; // Value Excluded
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);
  return(0);
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::ReplaceGlobalValues(int Row, int NumEntries, const double * srcValues, const int *Indices) {
  if(RowMap().GlobalIndicesInt())
  return TReplaceGlobalValues<int>(Row, NumEntries, srcValues, Indices);
  else
  throw ReportError("Epetra_CrsMatrix::ReplaceGlobalValues int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::ReplaceGlobalValues(long long Row, int NumEntries, const double * srcValues, const long long *Indices) {
  if(RowMap().GlobalIndicesLongLong())
  return TReplaceGlobalValues<long long>(Row, NumEntries, srcValues, Indices);
  else
  throw ReportError("Epetra_CrsMatrix::ReplaceGlobalValues long long version called for a matrix that is not long long.", -1);
}
#endif

//==========================================================================
int Epetra_CrsMatrix::ReplaceMyValues(int Row, int NumEntries, const double * srcValues, const int *Indices) {

  if (!IndicesAreLocal()) 
    EPETRA_CHK_ERR(-4); // Indices must be local.

  int j;
  int ierr = 0;
  int Loc;

  if (Row < 0 || Row >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }

  double* RowValues = Values(Row); 
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_.FindMyIndexLoc(Row,Index,j,Loc)) 
      RowValues[Loc] = srcValues[j];
    else 
      ierr = 2; // Value Excluded
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::ReplaceOffsetValues(long long Row, int NumEntries,
            const double * srcValues, const int *Offsets)
{
  int j;
  int ierr = 0;

  // Normalize row range
#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
  int locRow = LRID((int) Row);
#else
  int locRow = LRID(Row);
#endif
    
  if (locRow < 0 || locRow >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }

  double* RowValues = Values(locRow); 
  for(j=0; j<NumEntries; j++) {
    if( Offsets[j] != -1 )
      RowValues[Offsets[j]] = srcValues[j];
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
template<typename int_type>
int Epetra_CrsMatrix::TSumIntoGlobalValues(int_type Row,
            int NumEntries,
            const double * srcValues,
            const int_type *Indices)
{
  int j;
  int ierr = 0;
  int Loc = 0;

  int locRow = Graph_.LRID(Row); // Normalize row range
    
  if (locRow < 0 || locRow >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }

  if (StaticGraph() && !Graph_.HaveColMap()) {
    EPETRA_CHK_ERR(-1);
  }

  double * RowValues = Values(locRow);

  if (!StaticGraph()) {
    for (j=0; j<NumEntries; j++) {
      int_type Index = Indices[j];
      if (Graph_.FindGlobalIndexLoc(locRow,Index,j,Loc))
        RowValues[Loc] += srcValues[j];
      else
        ierr = 2; // Value Excluded
    }
  }
  else {
    const Epetra_BlockMap& colmap = Graph_.ColMap();
    int NumColIndices = Graph_.NumMyIndices(locRow);
    const int* ColIndices = Graph_.Indices(locRow);

    if (Graph_.Sorted()) {
      int insertPoint;
      for (j=0; j<NumEntries; j++) {
        int Index = colmap.LID(Indices[j]);

        // Check whether the next added element is the subsequent element in
        // the graph indices, then we can skip the binary search
        if (Loc < NumColIndices && Index == ColIndices[Loc])
          RowValues[Loc] += srcValues[j];
        else {
          Loc = Epetra_Util_binary_search(Index, ColIndices, NumColIndices, insertPoint);
          if (Loc > -1)
            RowValues[Loc] += srcValues[j];
          else 
            ierr = 2; // Value Excluded
        }
        ++Loc;
      }
    }
    else
      for (j=0; j<NumEntries; j++) {
        int Index = colmap.LID(Indices[j]);
        if (Graph_.FindMyIndexLoc(NumColIndices,ColIndices,Index,j,Loc)) 
          RowValues[Loc] += srcValues[j];
        else 
          ierr = 2; // Value Excluded
      }
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);

  return(0);
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::SumIntoGlobalValues(int Row,
            int NumEntries,
            const double * srcValues,
            const int *Indices)
{
  if(RowMap().GlobalIndicesInt())
  return TSumIntoGlobalValues<int>(Row, NumEntries, srcValues, Indices);
  else
  throw ReportError("Epetra_CrsMatrix::SumIntoGlobalValues int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::SumIntoGlobalValues(long long Row,
            int NumEntries,
            const double * srcValues,
            const long long *Indices)
{
  if(RowMap().GlobalIndicesLongLong())
  return TSumIntoGlobalValues<long long>(Row, NumEntries, srcValues, Indices);
  else
  throw ReportError("Epetra_CrsMatrix::SumIntoGlobalValues long long version called for a matrix that is not long long.", -1);
}
#endif

//==========================================================================
int Epetra_CrsMatrix::SumIntoMyValues(int Row, int NumEntries, const double * srcValues, const int *Indices) {

  if (!IndicesAreLocal()) 
    EPETRA_CHK_ERR(-4); // Indices must be local.

  int j;
  int ierr = 0;
  int Loc = 0;
  int insertPoint;

  if (Row < 0 || Row >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }

  double* RowValues = Values(Row);
  int NumColIndices = Graph_.NumMyIndices(Row);
  const int* ColIndices = Graph_.Indices(Row);
  if (Graph_.Sorted()) {
    for (j=0; j<NumEntries; j++) {
      int Index = Indices[j];

      // Check whether the next added element is the subsequent element in
      // the graph indices.
      if (Loc < NumColIndices && Index == ColIndices[Loc])
        RowValues[Loc] += srcValues[j];
      else {
        Loc = Epetra_Util_binary_search(Index, ColIndices, NumColIndices, insertPoint);
        if (Loc > -1)
          RowValues[Loc] += srcValues[j];
        else 
          ierr = 2; // Value Excluded
      }
      ++Loc;
    }
  }
  else {
    for (j=0; j<NumEntries; j++) {
      int Index = Indices[j];
      if (Graph_.FindMyIndexLoc(Row,Index,j,Loc)) 
        RowValues[Loc] += srcValues[j];
      else 
        ierr = 2; // Value Excluded
    }
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);

  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::SumIntoOffsetValues(long long Row, int NumEntries, const double * srcValues, const int *Offsets) {

  int j;
  int ierr = 0;

  // Normalize row range
#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
  int locRow = LRID((int) Row);
#else
  int locRow = LRID(Row);
#endif
    
  if (locRow < 0 || locRow >= NumMyRows_) {
    EPETRA_CHK_ERR(-1); // Not in Row range
  }

  double* RowValues = Values(locRow);
  for (j=0; j<NumEntries; j++) {
    if( Offsets[j] != -1 )
      RowValues[Offsets[j]] += srcValues[j];
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);

  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::FillComplete(bool OptimizeDataStorage) {
  squareFillCompleteCalled_ = true;
  EPETRA_CHK_ERR(FillComplete(RowMap(), RowMap(), OptimizeDataStorage));
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::FillComplete(const Epetra_Map& domain_map,
           const Epetra_Map& range_map, bool OptimizeDataStorage)
{
  int returnValue = 0;

  if (Graph_.Filled()) {
    if (!constructedWithFilledGraph_ && !matrixFillCompleteCalled_) {
      returnValue = 2;
    }
  }

  if (!StaticGraph()) {
    if (Graph_.MakeIndicesLocal(domain_map, range_map) < 0) {
      return(-1);
    }
  }
  SortEntries();  // Sort column entries from smallest to largest
  MergeRedundantEntries(); // Get rid of any redundant index values
  if (!StaticGraph()) {
    if (Graph_.FillComplete(domain_map, range_map) < 0) {
      return(-2);
    }
  }

  matrixFillCompleteCalled_ = true;

  if (squareFillCompleteCalled_) {
    if (DomainMap().NumGlobalElements64() != RangeMap().NumGlobalElements64()) {
      returnValue = 3;
    }
    squareFillCompleteCalled_ = false;
    EPETRA_CHK_ERR(returnValue);
  }

  if (OptimizeDataStorage) { EPETRA_CHK_ERR(OptimizeStorage()); }

  return(returnValue);
}

//==========================================================================
int Epetra_CrsMatrix::TransformToLocal() {
  EPETRA_CHK_ERR(FillComplete());
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::TransformToLocal(const Epetra_Map* domainMap, const Epetra_Map* rangeMap) {
  EPETRA_CHK_ERR(FillComplete(*domainMap, *rangeMap));
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::SortEntries() {

  if(!IndicesAreLocal()) 
    EPETRA_CHK_ERR(-1);
  if(Sorted())
    return(0);

  // For each row, sort column entries from smallest to largest.
  // Use shell sort. Stable sort so it is fast if indices are already sorted.

  
  for(int i = 0; i < NumMyRows_; i++){

    double* locValues = Values(i);
    int NumEntries = Graph().NumMyIndices(i);
    int* locIndices = Graph().Indices(i);
    
    int n = NumEntries;
    int m = n/2;
    
    while(m > 0) {
      int max = n - m;
      for(int j = 0; j < max; j++) {
  for(int k = j; k >= 0; k-=m) {
    if(locIndices[k+m] >= locIndices[k])
      break;
    double dtemp = locValues[k+m];
    locValues[k+m] = locValues[k];
    locValues[k] = dtemp;
    int itemp = locIndices[k+m];
    locIndices[k+m] = locIndices[k];
    locIndices[k] = itemp;
  }
      }
      m = m/2;
    }
  }
  Graph_.SetSorted(true); // This also sorted the graph
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::MergeRedundantEntries() {

  int i;

  if(NoRedundancies()) 
    return(0);
  if(!Sorted()) 
    EPETRA_CHK_ERR(-1);  // Must have sorted entries

  // For each row, remove column indices that are repeated.
  // Also, determine if matrix is upper or lower triangular or has no diagonal (Done in graph)
  // Note:  This function assumes that SortEntries was already called.

  for(i = 0; i<NumMyRows_; i++) {
    int NumEntries = Graph().NumMyIndices(i);
    if(NumEntries > 1) {
      double* const locValues = Values(i);
      int* const locIndices = Graph().Indices(i);    
      int curEntry =0;
      double curValue = locValues[0];
      for(int k = 1; k < NumEntries; k++) {
  if(locIndices[k] == locIndices[k-1]) 
    curValue += locValues[k];
  else {
    locValues[curEntry++] = curValue;
    curValue = locValues[k];
  }
      }
      locValues[curEntry] = curValue;
      
    }
  }
  
  EPETRA_CHK_ERR(Graph_.RemoveRedundantIndices()); // Remove redundant indices and then return
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::OptimizeStorage() {


  if (StorageOptimized()) 
    return(0); // Have we been here before?
  if (!Filled()) EPETRA_CHK_ERR(-1); // Cannot optimize storage before calling FillComplete()


  int ierr = Graph_.OptimizeStorage();
  if (ierr!=0) EPETRA_CHK_ERR(ierr);  // In order for OptimizeStorage to make sense for the matrix, it must work on the graph.

  bool Contiguous = true; // Assume contiguous is true
  for (int i=1; i<NumMyRows_; i++){
    int NumEntries = Graph().NumMyIndices(i-1);
    
    // check if end of beginning of current row starts immediately after end of previous row.
    if (Values_[i]!=Values_[i-1]+NumEntries) {
      Contiguous = false;
      break;
    }
  }

  // NOTE:  At the end of the above loop set, there is a possibility that NumEntries and NumAllocatedEntries
  //        for the last row could be different, but I don't think it matters.


  if ((CV_==View) && !Contiguous) 
    EPETRA_CHK_ERR(-1);  // This is user data, it's not contiguous and we can't make it so.

  
  if(!Contiguous) { // Must pack indices if not already contiguous.  Pack values into All_values_
    
    if (All_Values_==0) {
      // Compute Number of Nonzero entries (Done in FillComplete, but we may not have been there yet.)
      int numMyNonzeros = Graph_.NumMyNonzeros();
    
      // Allocate one big array for all values
      All_Values_ = new double[numMyNonzeros];
      if(All_Values_ == 0) throw ReportError("Error with All_Values_ allocation.", -99);

      const int *const  IndexOffset = Graph().IndexOffset();
      const int numMyRows = NumMyRows_;
      double ** Values_s = Values_;
      double * All_Values_s = All_Values_;
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(Values_s,All_Values_s)
#endif
      for (int i=0; i<numMyRows; i++) {
        int NumEntries = Graph().NumMyIndices(i);
        int curOffset = IndexOffset[i];
        double * values = Values_s[i];
        double * newValues = All_Values_s+curOffset;
  for (int j=0; j<NumEntries; j++) newValues[j] = values[j];
      }
    }
    else { // Static Profile, so just pack into existing storage (can't be threaded)
      double * tmp = All_Values_;
      for (int i=0; i<NumMyRows_; i++) {
        int NumEntries = Graph().NumMyIndices(i);
        double * values = Values_[i];
        if (tmp!=values) // Copy values if not pointing to same location
          for (int j=0; j<NumEntries; j++) tmp[j] = values[j];
        tmp += NumEntries;
      }
    }

    // Free Values_ arrays
    for (int i=0;i<NumMyRows_; ++i) {
      if (Values_alloc_lengths_[i] != 0) delete [] Values_[i];
    }
    
    delete [] Values_alloc_lengths_; Values_alloc_lengths_ = 0;
  } // End of !Contiguous section
  else {
    //if already contiguous, we'll simply set All_Values_ to be
    //a copy of Values_[0].
    if (All_Values_!=0 && All_Values_!=Values_[0]) delete [] All_Values_;
    All_Values_ = NumMyRows_ > 0 ? Values_[0] : NULL;
  }
  
  // Delete unneeded storage
  delete [] Values_; Values_=0;

#ifdef Epetra_ENABLE_CASK 
  if (cask == NULL  && Graph().StorageOptimized() )  {
     int * Indices = Graph().All_Indices();
     int * IndexOffset = Graph().IndexOffset();
     int NumMyCols_ = NumMyCols();
     cask_handler_initialize(&cask);
     cask_csr_analysis(NumMyRows_, NumMyCols_, IndexOffset, Indices,
                       NumGlobalNonzeros64(),cask);
  }
#endif


  StorageOptimized_ = true;

  
  return(0);
}

//==========================================================================
template<typename int_type>
int Epetra_CrsMatrix::ExtractGlobalRowCopy(int_type Row, int Length, int & NumEntries, double * values,
             int_type * Indices) const 
{

  int ierr = Graph_.ExtractGlobalRowCopy(Row, Length, NumEntries, Indices);
  if (ierr) 
    EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractGlobalRowCopy(Row, Length, NumEntries, values));
  return(0);
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::ExtractGlobalRowCopy(int Row, int Length, int & NumEntries, double * values,
             int * Indices) const 
{
  if(RowMap().GlobalIndicesInt())
    return ExtractGlobalRowCopy<int>(Row, Length, NumEntries, values, Indices);
  else
    throw ReportError("Epetra_CrsMatrix::ExtractGlobalRowCopy int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::ExtractGlobalRowCopy(long long Row, int Length, int & NumEntries, double * values,
             long long * Indices) const 
{
  if(RowMap().GlobalIndicesLongLong())
    return ExtractGlobalRowCopy<long long>(Row, Length, NumEntries, values, Indices);
  else
    throw ReportError("Epetra_CrsMatrix::ExtractGlobalRowCopy long long version called for a matrix that is not long long.", -1);
}
#endif

//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, double * values,
               int * Indices) const 
{

  int ierr = Graph_.ExtractMyRowCopy(Row, Length, NumEntries, Indices);
  if (ierr) 
    EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractMyRowCopy(Row, Length, NumEntries, values));
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::NumMyRowEntries(int Row, int & NumEntries) const 
{

  if (!MyLRID(Row)) 
    EPETRA_CHK_ERR(-1); // Not in the range of local rows
  NumEntries = NumMyEntries(Row);
  return(0);
}

//==========================================================================
template<typename int_type>
int Epetra_CrsMatrix::ExtractGlobalRowCopy(int_type Row, int Length, int & NumEntries, double * values) const 
{
  int Row0 = Graph_.RowMap().LID(Row); // Normalize row range

  EPETRA_CHK_ERR(ExtractMyRowCopy(Row0, Length, NumEntries, values));
  return(0);
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::ExtractGlobalRowCopy(int Row, int Length, int & NumEntries, double * values) const
{
  if(RowMap().GlobalIndicesInt())
    return ExtractGlobalRowCopy<int>(Row, Length, NumEntries, values);
  else
    throw ReportError("Epetra_CrsMatrix::ExtractGlobalRowCopy int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::ExtractGlobalRowCopy(long long Row, int Length, int & NumEntries, double * values) const
{
  if(RowMap().GlobalIndicesLongLong())
    return ExtractGlobalRowCopy<long long>(Row, Length, NumEntries, values);
  else
    throw ReportError("Epetra_CrsMatrix::ExtractGlobalRowCopy long long version called for a matrix that is not long long.", -1);
}
#endif

//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, double * targValues) const 
{
  int j;

  if (Row < 0 || Row >= NumMyRows_) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  NumEntries = Graph().NumMyIndices(Row);
  if (Length < NumEntries) 
    EPETRA_CHK_ERR(-2); // Not enough space for copy. Needed size is passed back in NumEntries

  double * srcValues = Values(Row);

  for(j=0; j<NumEntries; j++)
    targValues[j] = srcValues[j];
  
  return(0);
}

//==============================================================================
int Epetra_CrsMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const {
  
  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Can't get diagonal unless matrix is filled (and in local index space)
  if(!RowMap().SameAs(Diagonal.Map())) 
    EPETRA_CHK_ERR(-2); // Maps must be the same

  for(int i = 0; i < NumMyRows_; i++) {
    long long ii = GRID64(i);
    int NumEntries = Graph().NumMyIndices(i);
    int* Indices = Graph().Indices(i);
    double * srcValues = Values(i);
    
    Diagonal[i] = 0.0;
    for(int j = 0; j < NumEntries; j++) {
      if(ii == GCID64(Indices[j])) {
  Diagonal[i] = srcValues[j];
  break;
      }
    }
  }
  return(0);
}
//==============================================================================
int Epetra_CrsMatrix::ReplaceDiagonalValues(const Epetra_Vector & Diagonal) {
  
  if(!Filled())
    EPETRA_CHK_ERR(-1); // Can't replace diagonal unless matrix is filled (and in local index space)
  if(!RowMap().SameAs(Diagonal.Map())) 
    EPETRA_CHK_ERR(-2); // Maps must be the same

  int ierr = 0;
  for(int i = 0; i < NumMyRows_; i++) {
    long long ii = GRID64(i);
    int NumEntries = Graph().NumMyIndices(i);
    int* Indices = Graph().Indices(i);
    double * targValues = Values(i);
    bool DiagMissing = true;
    for(int j = 0; j < NumEntries; j++) {
      if(ii == GCID64(Indices[j])) {
  targValues[j] = Diagonal[i];
  DiagMissing = false;
  break;
      }
    }
    if(DiagMissing) 
      ierr = 1; // flag a warning error
  }

  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  EPETRA_CHK_ERR(ierr);

  return(0);
}

//==========================================================================
template<typename int_type>
int Epetra_CrsMatrix::ExtractGlobalRowView(int_type Row, int & NumEntries, double *& values, int_type *& Indices) const 
{
  int ierr = Graph_.ExtractGlobalRowView(Row, NumEntries, Indices);
  if (ierr) 
    EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractGlobalRowView(Row, NumEntries, values));
  return(0);
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::ExtractGlobalRowView(int Row, int & NumEntries, double *& values, int *& Indices) const 
{
  if(RowMap().GlobalIndicesInt())
    return ExtractGlobalRowView<int>(Row, NumEntries, values, Indices);
  else
    throw ReportError("Epetra_CrsMatrix::ExtractGlobalRowView int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::ExtractGlobalRowView(long long Row, int & NumEntries, double *& values, long long *& Indices) const 
{
  if(RowMap().GlobalIndicesLongLong())
    return ExtractGlobalRowView<long long>(Row, NumEntries, values, Indices);
  else
    throw ReportError("Epetra_CrsMatrix::ExtractGlobalRowView long long version called for a matrix that is not long long.", -1);
}
#endif

//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowView(int Row, int & NumEntries, double *& values, int *& Indices) const 
{
  int ierr = Graph_.ExtractMyRowView(Row, NumEntries, Indices);
  if (ierr) 
    EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractMyRowView(Row, NumEntries, values));
  return(0);
}
//==========================================================================
template<typename int_type>
int Epetra_CrsMatrix::ExtractGlobalRowView(int_type Row, int & NumEntries, double *& values) const 
{
  int Row0 = Graph_.RowMap().LID(Row); // Normalize row range

  EPETRA_CHK_ERR(ExtractMyRowView(Row0, NumEntries, values));
  return(0);
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::ExtractGlobalRowView(int Row, int & NumEntries, double *& values) const 
{
  if(RowMap().GlobalIndicesInt())
    return ExtractGlobalRowView<int>(Row, NumEntries, values);
  else
    throw ReportError("Epetra_CrsMatrix::ExtractGlobalRowView int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_CrsMatrix::ExtractGlobalRowView(long long Row, int & NumEntries, double *& values) const 
{
  if(RowMap().GlobalIndicesLongLong())
    return ExtractGlobalRowView<long long>(Row, NumEntries, values);
  else
    throw ReportError("Epetra_CrsMatrix::ExtractGlobalRowView long long version called for a matrix that is not long long.", -1);
}
#endif

//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowView(int Row, int & NumEntries, double *& targValues) const 
{
  if (Row < 0 || Row >= NumMyRows_) 
    EPETRA_CHK_ERR(-1); // Not in Row range

  NumEntries = Graph().NumMyIndices(Row);

  targValues = Values(Row);
  
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal,
          const Epetra_Vector& x, Epetra_Vector& y) const
{

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Solve(Upper,Trans,UnitDiag,x,y)");
#endif

  //
  // This function finds y such that Ly = x or Uy = x or the transpose cases.
  //

  // First short-circuit to the pre-5.0 version if no storage optimization was performed
  if (!StorageOptimized() && !Graph().StorageOptimized()) {
    EPETRA_CHK_ERR(Solve1(Upper, Trans, UnitDiagonal, x, y));
    return(0);
  }

  if (!Filled()) {
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  }

  if ((Upper) && (!UpperTriangular())) 
    EPETRA_CHK_ERR(-2);
  if ((!Upper) && (!LowerTriangular())) 
    EPETRA_CHK_ERR(-3);
  if ((!UnitDiagonal) && (NoDiagonal())) 
    EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) 
    EPETRA_CHK_ERR(-5); // Need each row to have a diagonal

  double *xp = (double*)x.Values();
  double *yp = (double*)y.Values();

      
  GeneralSV(Upper, Trans, UnitDiagonal, xp, yp);

  UpdateFlops(2*NumGlobalNonzeros64());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Solve(Upper,Trans,UnitDiag,X,Y)");
#endif

  //
  // This function find Y such that LY = X or UY = X or the transpose cases.
  //

  // First short-circuit to the pre-5.0 version if no storage optimization was performed
  if (!StorageOptimized() && !Graph().StorageOptimized()) {
    EPETRA_CHK_ERR(Solve1(Upper, Trans, UnitDiagonal, X, Y));
    return(0);
  }
  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  if((Upper) && (!UpperTriangular()))
    EPETRA_CHK_ERR(-2);
  if((!Upper) && (!LowerTriangular()))
    EPETRA_CHK_ERR(-3);
  if((!UnitDiagonal) && (NoDiagonal()))
    EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_))
    EPETRA_CHK_ERR(-5); // Need each row to have a diagonal

  double** Xp = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();
  int LDX = X.ConstantStride() ? X.Stride() : 0;
  int LDY = Y.ConstantStride() ? Y.Stride() : 0;
  int NumVectors = X.NumVectors();


    // Do actual computation
  if (NumVectors==1)
    GeneralSV(Upper, Trans, UnitDiagonal, *Xp, *Yp);
  else
    GeneralSM(Upper, Trans, UnitDiagonal, Xp, LDX, Yp, LDY, NumVectors);

  UpdateFlops(2 * NumVectors * NumGlobalNonzeros64());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::InvRowSums(Epetra_Vector& x) const {
  //
  // Put inverse of the sum of absolute values of the ith row of A in x[i].
  //

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  int ierr = 0;
  int i, j;
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double * xp = (double*)x.Values();
  if (Graph().RangeMap().SameAs(x.Map()) && Exporter() != 0) {
    Epetra_Vector x_tmp(RowMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = NumMyEntries(i);
      double * RowValues  = Values(i);
      for (j=0; j < NumEntries; j++)  x_tmp_p[i] += std::abs(RowValues[j]);
    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *Exporter(), Add)); //Export partial row sums to x.
    int myLength = x.MyLength();
    for (i=0; i<myLength; i++) { 
      if (xp[i]<Epetra_MinDouble) {
        if (xp[i]==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
        else if (ierr!=1) ierr = 2;
        xp[i] = Epetra_MaxDouble;
      }
      else
        xp[i] = 1.0/xp[i];
    }
  }
  else if (Graph().RowMap().SameAs(x.Map())) {
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = NumMyEntries(i);
      double * RowValues  = Values(i);
      double scale = 0.0;
      for (j=0; j < NumEntries; j++) scale += std::abs(RowValues[j]);
      if (scale<Epetra_MinDouble) {
        if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
        else if (ierr!=1) ierr = 2;
        xp[i] = Epetra_MaxDouble;
      }
      else
        xp[i] = 1.0/scale;
    }
  }
  else { // x.Map different than both Graph().RowMap() and Graph().RangeMap()
    EPETRA_CHK_ERR(-2); // The map of x must be the RowMap or RangeMap of A.
  }
  UpdateFlops(NumGlobalNonzeros64());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::InvRowMaxs(Epetra_Vector& x) const {
  //
  // Put inverse of the max of absolute values of the ith row of A in x[i].
  //

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  int ierr = 0;
  int i, j;
  bool needExport = false;
  double * xp = (double*)x.Values();
  Epetra_Vector* x_tmp = 0;
  if (Graph().RangeMap().SameAs(x.Map())) {
    if (Exporter() != 0) {
      needExport = true; //Having this information later avoids a .SameAs
      x_tmp = new Epetra_Vector(RowMap()); // Create import vector if needed
      xp = (double*)x_tmp->Values();
    }
  }
  else if (!Graph().RowMap().SameAs(x.Map())) {
    EPETRA_CHK_ERR(-2); // The map of x must be the RowMap or RangeMap of A.
  }
  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = NumMyEntries(i);
    double * RowValues  = Values(i);
    double scale = 0.0;
    for (j=0; j < NumEntries; j++) scale = EPETRA_MAX(std::abs(RowValues[j]),scale);
    if (scale<Epetra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowmax found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      xp[i] = Epetra_MaxDouble;
    }
    else
      xp[i] = 1.0/scale;
  }
  if(needExport) {
    x.PutScalar(0.0);
    EPETRA_CHK_ERR(x.Export(*x_tmp, *Exporter(), Insert)); // Fill x with values from temp vector
    delete x_tmp;
  }
  UpdateFlops(NumGlobalNonzeros64());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::InvColSums(Epetra_Vector& x) const {
  //
  // Put inverse of the sum of absolute values of the jth column of A in x[j].
  //

  if(!Filled())  EPETRA_CHK_ERR(-1); // Matrix must be filled.
  int ierr = 0;
  int i, j;
  int MapNumMyElements = x.Map().NumMyElements();
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double* xp = (double*)x.Values();
  if(Graph().DomainMap().SameAs(x.Map()) && Importer() != 0) {
    Epetra_Vector x_tmp(ColMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    ColIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for(j = 0; j < NumEntries; j++) 
        x_tmp_p[ColIndices[j]] += std::abs(RowValues[j]);
    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *Importer(), Add)); // Fill x with partial column sums
  }
  else if(Graph().ColMap().SameAs(x.Map())) {
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    ColIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for(j = 0; j < NumEntries; j++) 
        xp[ColIndices[j]] += std::abs(RowValues[j]);
    }
  }
  else { //x.Map different than both Graph().ColMap() and Graph().DomainMap()
    EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A
  }

  // Invert values, don't allow them to get too large
  for(i = 0; i < MapNumMyElements; i++) {
    double scale = xp[i];
    if(scale < Epetra_MinDouble) {
      if(scale == 0.0) 
  ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if(ierr != 1) 
  ierr = 2;
      xp[i] = Epetra_MaxDouble;
    }
    else
      xp[i] = 1.0 / scale;
  }

  UpdateFlops(NumGlobalNonzeros64());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::InvColMaxs(Epetra_Vector& x) const {
  //
  // Put inverse of the max of absolute values of the jth column of A in x[j].
  //

  if(!Filled())  EPETRA_CHK_ERR(-1); // Matrix must be filled.
  int ierr = 0;
  int i, j;
  int MapNumMyElements = x.Map().NumMyElements();
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double* xp = (double*)x.Values();
  if(Graph().DomainMap().SameAs(x.Map()) && Importer() != 0) {
    Epetra_Vector x_tmp(ColMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    ColIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for(j = 0; j < NumEntries; j++) 
        x_tmp_p[ColIndices[j]] = EPETRA_MAX(std::abs(RowValues[j]),x_tmp_p[ColIndices[j]]);
    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *Importer(), AbsMax)); // Fill x with partial column sums
  }
  else if(Graph().ColMap().SameAs(x.Map())) {
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    ColIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for(j = 0; j < NumEntries; j++) 
        xp[ColIndices[j]] = EPETRA_MAX(std::abs(RowValues[j]),xp[ColIndices[j]]);
    }
  }
  else { //x.Map different than both Graph().ColMap() and Graph().DomainMap()
    EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A
  }

  // Invert values, don't allow them to get too large
  for(i = 0; i < MapNumMyElements; i++) {
    double scale = xp[i];
    if(scale < Epetra_MinDouble) {
      if(scale == 0.0) 
  ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if(ierr != 1) 
  ierr = 2;
      xp[i] = Epetra_MaxDouble;
    }
    else
      xp[i] = 1.0 / scale;
  }

  UpdateFlops(NumGlobalNonzeros64());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::LeftScale(const Epetra_Vector& x) {
  //
  // This function scales the ith row of A by x[i].
  //

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  double* xp = 0;
  if(Graph().RangeMap().SameAs(x.Map()))  
    // If we have a non-trivial exporter, we must import elements that are 
    // permuted or are on other processors.  (We will use the exporter to
    // perform the import.)
    if(Exporter() != 0) {
      UpdateExportVector(1);
      EPETRA_CHK_ERR(ExportVector_->Import(x,*Exporter(), Insert));
      xp = (double*) ExportVector_->Values();
    }
    else
      xp = (double*)x.Values();
  else if (Graph().RowMap().SameAs(x.Map()))
    xp = (double*)x.Values();
  else {
    EPETRA_CHK_ERR(-2); // The Map of x must be the RowMap or RangeMap of A.
  }
  int i, j;

  for(i = 0; i < NumMyRows_; i++) {
    int      NumEntries = NumMyEntries(i);
    double* RowValues  = Values(i);
    double scale = xp[i];
    for(j = 0; j < NumEntries; j++)  
      RowValues[j] *= scale;
  }
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  UpdateFlops(NumGlobalNonzeros64());

  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::RightScale(const Epetra_Vector& x) {
  //
  // This function scales the jth column of A by x[j].
  //

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  double* xp = 0;
  if(Graph().DomainMap().SameAs(x.Map())) 
    // If we have a non-trivial exporter, we must import elements that are 
    // permuted or are on other processors.
    if(Importer() != 0) {
      UpdateImportVector(1);
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*) ImportVector_->Values();
    }
    else
      xp = (double*)x.Values();
  else if(Graph().ColMap().SameAs(x.Map()))
    xp = (double*)x.Values(); 
  else
    EPETRA_CHK_ERR(-2); // The Map of x must be the RowMap or RangeMap of A.
  int i, j;

  for(i = 0; i < NumMyRows_; i++) {
    int     NumEntries = NumMyEntries(i);
    int*    ColIndices = Graph().Indices(i);
    double* RowValues  = Values(i);
    for(j = 0; j < NumEntries; j++)  
      RowValues[j] *=  xp[ColIndices[j]];
  }
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  NormFrob_ = -1.0;

  UpdateFlops(NumGlobalNonzeros64());
  return(0);
}

//=============================================================================
double Epetra_CrsMatrix::NormInf() const {

#if 0
  //
  //  Commenting this section out disables caching, ie. 
  //  causes the norm to be computed each time NormInf is called.
  //  See bug #1151 for a full discussion.  
  //
  double MinNorm ; 
  Comm().MinAll( &NormInf_, &MinNorm, 1 ) ; 

  if( MinNorm >= 0.0) 
    return(NormInf_);
#endif

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  Epetra_Vector x(RangeMap()); // Need temp vector for row sums
  double* xp = (double*)x.Values();
  Epetra_MultiVector* x_tmp = 0;

  // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
  if(Exporter() != 0) {
    x_tmp = new Epetra_Vector(RowMap()); // Create temporary import vector if needed
    xp = (double*)x_tmp->Values();
  }
  int i, j;

  for(i = 0; i < NumMyRows_; i++) {
    xp[i] = 0.0;
    int     NumEntries = NumMyEntries(i);
    double* RowValues  = Values(i);
    for(j = 0; j < NumEntries; j++) 
      xp[i] += std::abs(RowValues[j]);
  }
  if(Exporter() != 0) {
    x.PutScalar(0.0);
    EPETRA_CHK_ERR(x.Export(*x_tmp, *Exporter(), Add)); // Fill x with Values from temp vector
  }
  x.MaxValue(&NormInf_); // Find max
  if(x_tmp != 0) 
    delete x_tmp;
  UpdateFlops(NumGlobalNonzeros64());
  return(NormInf_);
}
//=============================================================================
double Epetra_CrsMatrix::NormOne() const {

#if 0
  //
  //  Commenting this section out disables caching, ie. 
  //  causes the norm to be computed each time NormOne is called.  
  //  See bug #1151 for a full discussion.  
  //
  double MinNorm ; 
  Comm().MinAll( &NormOne_, &MinNorm, 1 ) ; 

  if( MinNorm >= 0.0) 
    return(NormOne_);
#endif

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  Epetra_Vector x(DomainMap()); // Need temp vector for column sums
  
  double* xp = (double*)x.Values();
  Epetra_MultiVector* x_tmp = 0;
  int NumCols = NumMyCols();
  

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if(Importer() != 0) {
    x_tmp = new Epetra_Vector(ColMap()); // Create temporary import vector if needed
    xp = (double*)x_tmp->Values();
  }
  int i, j;

  for(i = 0; i < NumCols; i++) 
    xp[i] = 0.0;

  for(i = 0; i < NumMyRows_; i++) {
    int     NumEntries = NumMyEntries(i);
    int*    ColIndices = Graph().Indices(i);
    double* RowValues  = Values(i);
    for(j = 0; j < NumEntries; j++) 
      xp[ColIndices[j]] += std::abs(RowValues[j]);
  }
  if(Importer() != 0) {
    x.PutScalar(0.0);
    EPETRA_CHK_ERR(x.Export(*x_tmp, *Importer(), Add)); // Fill x with Values from temp vector
  }
  x.MaxValue(&NormOne_); // Find max
  if(x_tmp != 0) 
    delete x_tmp;
  UpdateFlops(NumGlobalNonzeros64());
  return(NormOne_);
}
//=============================================================================
double Epetra_CrsMatrix::NormFrobenius() const {

#if 0
  //
  //  Commenting this section out disables caching, ie. 
  //  causes the norm to be computed each time NormFrobenius is called.  
  //  See bug #1151 for a full discussion.  
  //
  double MinNorm ; 
  Comm().MinAll( &NormFrob_, &MinNorm, 1 ) ; 

  if( MinNorm >= 0.0) 
    return(NormFrob_);
#endif

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  double local_sum = 0.0;

  for(int i = 0; i < NumMyRows_; i++) {
    int     NumEntries = NumMyEntries(i);
    double* RowValues  = Values(i);
    for(int j = 0; j < NumEntries; j++) {
      local_sum += RowValues[j]*RowValues[j];
    }
  }

  double global_sum = 0.0;
  Comm().SumAll(&local_sum, &global_sum, 1);

  NormFrob_ = std::sqrt(global_sum);

  UpdateFlops(NumGlobalNonzeros64());

  return(NormFrob_);
}
//=========================================================================
int Epetra_CrsMatrix::CheckSizes(const Epetra_SrcDistObject & Source) {
  try {
    const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Source);
    if (!A.Graph().GlobalConstantsComputed()) EPETRA_CHK_ERR(-1); // Must have global constants to proceed
  }
  catch (...) {
    return(0); // No error at this point, object could be a RowMatrix
  }
  return(0);
}

//=========================================================================
int Epetra_CrsMatrix::CopyAndPermute(const Epetra_SrcDistObject & Source,
             int NumSameIDs, 
             int NumPermuteIDs,
                                     int * PermuteToLIDs,
             int *PermuteFromLIDs,
                                     const Epetra_OffsetIndex * Indexor ) {
 
  if(!Source.Map().GlobalIndicesTypeMatch(RowMap()))
    throw ReportError("Epetra_CrsMatrix::CopyAndPermute: Incoming global index type does not match the one for *this",-1);

  try {
    const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Source);
    EPETRA_CHK_ERR(CopyAndPermuteCrsMatrix(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs,
             PermuteFromLIDs,Indexor));
  }
  catch (...) {
    try {
      const Epetra_RowMatrix & A = dynamic_cast<const Epetra_RowMatrix &>(Source);
      EPETRA_CHK_ERR(CopyAndPermuteRowMatrix(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs,
               PermuteFromLIDs,Indexor));
    }
    catch (...) {
      EPETRA_CHK_ERR(-1); // Incompatible SrcDistObject
    }
  }
  
  return(0);
}

//=========================================================================
template<typename int_type>
int Epetra_CrsMatrix::TCopyAndPermuteCrsMatrix(const Epetra_CrsMatrix & A,
                                              int NumSameIDs, 
                int NumPermuteIDs,
                                              int * PermuteToLIDs,
                int *PermuteFromLIDs,
                                              const Epetra_OffsetIndex * Indexor) {
  
  int i, ierr;
  
  int_type Row;
  int NumEntries;
  int maxNumEntries = A.MaxNumEntries();
  int_type * Indices = 0;
  double * values = 0;

  if (maxNumEntries>0 && A.IndicesAreLocal() ) { //Need Temp Space
    Indices = new int_type[maxNumEntries];
    values = new double[maxNumEntries];
  }
  
  // Do copy first
  if (NumSameIDs>0) {
    if (A.IndicesAreLocal()) {
      if (StaticGraph() || IndicesAreLocal()) {
        if(Indexor) {
          for (i=0; i<NumSameIDs; i++) {
      Row = (int_type) GRID64(i);
      EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, maxNumEntries, NumEntries, values, Indices)); // Set pointers
            ierr = ReplaceOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumSameIDs; i++) {
      Row = (int_type) GRID64(i);
      EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, maxNumEntries, NumEntries, values, Indices)); // Set pointers
            ierr = ReplaceGlobalValues(Row, NumEntries, values, Indices);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
      }
      else {
        if(Indexor) {
          for (i=0; i<NumSameIDs; i++) {
      Row = (int_type) GRID64(i);
      EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, maxNumEntries, NumEntries, values, Indices)); // Set pointers
            ierr = InsertOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumSameIDs; i++) {
      Row = (int_type) GRID64(i);
      EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(Row, maxNumEntries, NumEntries, values, Indices)); // Set pointers
            ierr = InsertGlobalValues(Row, NumEntries, values, Indices);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
      } 
    }
    else { // A.IndicesAreGlobal()
      if (StaticGraph() || IndicesAreLocal()) {
        if(Indexor) {
          for (i=0; i<NumSameIDs; i++) {
      Row = (int_type) GRID64(i);
      EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumEntries, values, Indices)); // Set pointers
            ierr = ReplaceOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumSameIDs; i++) {
      Row = (int_type) GRID64(i);
      EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumEntries, values, Indices)); // Set pointers
            ierr = ReplaceGlobalValues(Row, NumEntries, values, Indices);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
      }
      else {
        if(Indexor) {
          for (i=0; i<NumSameIDs; i++) {
      Row = (int_type) GRID64(i);
      EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumEntries, values, Indices)); // Set pointers
            ierr = InsertOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumSameIDs; i++) {
      Row = (int_type) GRID64(i);
      EPETRA_CHK_ERR(A.ExtractGlobalRowView(Row, NumEntries, values, Indices)); // Set pointers
            ierr = InsertGlobalValues(Row, NumEntries, values, Indices);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
      } 
    }
  }
  
  // Do local permutation next
  int_type FromRow, ToRow;
  if (NumPermuteIDs>0) {
    if (A.IndicesAreLocal()) {
      if (StaticGraph() || IndicesAreLocal()) {
        if(Indexor) {
          for (i=0; i<NumPermuteIDs; i++) {
      FromRow = (int_type) A.GRID64(PermuteFromLIDs[i]);
      ToRow = (int_type) GRID64(PermuteToLIDs[i]);
      EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices)); // Set pointers
      ierr = ReplaceOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumPermuteIDs; i++) {
      FromRow = (int_type) A.GRID64(PermuteFromLIDs[i]);
      ToRow = (int_type) GRID64(PermuteToLIDs[i]);
      EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices)); // Set pointers
      ierr = ReplaceGlobalValues(ToRow, NumEntries, values, Indices);
            if( ierr<0 ) EPETRA_CHK_ERR(ierr);
          }
        }
      }
      else {
        if(Indexor) {
          for (i=0; i<NumPermuteIDs; i++) {
      FromRow = (int_type) A.GRID64(PermuteFromLIDs[i]);
      ToRow = (int_type) GRID64(PermuteToLIDs[i]);
      EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices)); // Set pointers
      ierr = InsertOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]);
      if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumPermuteIDs; i++) {
      FromRow = (int_type) A.GRID64(PermuteFromLIDs[i]);
      ToRow = (int_type) GRID64(PermuteToLIDs[i]);
      EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, maxNumEntries, NumEntries, values, Indices)); // Set pointers
      ierr = InsertGlobalValues(ToRow, NumEntries, values, Indices);
      if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
      }
    }
    else { // A.IndicesAreGlobal()
      if (StaticGraph() || IndicesAreLocal()) {
        if(Indexor) {
          for (i=0; i<NumPermuteIDs; i++) {
      FromRow = (int_type) A.GRID64(PermuteFromLIDs[i]);
      ToRow = (int_type) GRID64(PermuteToLIDs[i]);
      EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumEntries, values, Indices)); // Set pointers
      ierr = ReplaceOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]);
      if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumPermuteIDs; i++) {
      FromRow = (int_type) A.GRID64(PermuteFromLIDs[i]);
      ToRow = (int_type) GRID64(PermuteToLIDs[i]);
      EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumEntries, values, Indices)); // Set pointers
      ierr = ReplaceGlobalValues(ToRow, NumEntries, values, Indices);
      if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
      }
      else {
        if(Indexor) {
          for (i=0; i<NumPermuteIDs; i++) {
      FromRow = (int_type) A.GRID64(PermuteFromLIDs[i]);
      ToRow = (int_type) GRID64(PermuteToLIDs[i]);
      EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumEntries, values, Indices)); // Set pointers
      ierr = InsertOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]);
      if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
        else {
          for (i=0; i<NumPermuteIDs; i++) {
      FromRow = (int_type) A.GRID64(PermuteFromLIDs[i]);
      ToRow = (int_type) GRID64(PermuteToLIDs[i]);
      EPETRA_CHK_ERR(A.ExtractGlobalRowView(FromRow, NumEntries, values, Indices)); // Set pointers
      ierr = InsertGlobalValues(ToRow, NumEntries, values, Indices);
      if (ierr<0) EPETRA_CHK_ERR(ierr);
          }
        }
      }
    }
  }

  if (maxNumEntries>0 && A.IndicesAreLocal() ) { // Delete Temp Space
    delete [] values;
    delete [] Indices;
  }
  
  return(0);
}

int Epetra_CrsMatrix::CopyAndPermuteCrsMatrix(const Epetra_CrsMatrix & A,
                                              int NumSameIDs, 
                int NumPermuteIDs,
                                              int * PermuteToLIDs,
                int *PermuteFromLIDs,
                                              const Epetra_OffsetIndex * Indexor)
{
  if(!A.RowMap().GlobalIndicesTypeMatch(RowMap()))
    throw ReportError("Epetra_CrsMatrix::CopyAndPermuteCrsMatrix: Incoming global index type does not match the one for *this",-1);

  if(A.RowMap().GlobalIndicesInt())
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    return TCopyAndPermuteCrsMatrix<int>(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs, Indexor);
#else
    throw ReportError("Epetra_CrsMatrix::CopyAndPermuteCrsMatrix: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif

  if(A.RowMap().GlobalIndicesLongLong())
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    return TCopyAndPermuteCrsMatrix<long long>(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs, Indexor);
#else
    throw ReportError("Epetra_CrsMatrix::CopyAndPermuteCrsMatrix: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif

  throw ReportError("Epetra_CrsMatrix::CopyAndPermuteCrsMatrix: Internal error, unable to determine global index type of maps", -1);
}

//=========================================================================
template<typename int_type>
int Epetra_CrsMatrix::TCopyAndPermuteRowMatrix(const Epetra_RowMatrix & A,
                                              int NumSameIDs, 
                int NumPermuteIDs,
                                              int * PermuteToLIDs,
                int *PermuteFromLIDs,
                                              const Epetra_OffsetIndex * Indexor ) {
  int i, j, ierr;
  
  int_type Row, ToRow;
  int NumEntries;
  int FromRow;
  int maxNumEntries = A.MaxNumEntries();
  int_type * gen_Indices = 0; // gen = general (int or long long)
  int * int_Indices = 0;
  double * values = 0;

  if (maxNumEntries>0) {
    if(StaticGraph() || IndicesAreLocal() || Indexor) {
      int_Indices = new int[maxNumEntries];
    }
    else {
      gen_Indices = new int_type[maxNumEntries];
      int_Indices = reinterpret_cast<int*>(gen_Indices);
    }

    values = new double[maxNumEntries]; // Must extract values even though we discard them
  }
  
  const Epetra_Map & rowMap = A.RowMatrixRowMap();
  const Epetra_Map & colMap = A.RowMatrixColMap();

  // Do copy first
  if (NumSameIDs>0) {
    if (StaticGraph() || IndicesAreLocal()) {
      if( Indexor ) {
        for (i=0; i<NumSameIDs; i++) {
          Row = (int_type) GRID64(i);
          int AlocalRow = rowMap.LID(Row);
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(AlocalRow, maxNumEntries, NumEntries, values, int_Indices));
    ierr = ReplaceOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]);
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
      else {
        for (i=0; i<NumSameIDs; i++) {
          Row = (int_type) GRID64(i);
          int AlocalRow = rowMap.LID(Row);
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(AlocalRow, maxNumEntries, NumEntries, values, int_Indices));
          for(j=0; j<NumEntries; ++j) {
            int_Indices[j] = LCID((int_type) colMap.GID64(int_Indices[j]));
          }
    ierr = ReplaceMyValues(i, NumEntries, values, int_Indices);
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
    }
    else {
      if( Indexor ) {
        for (i=0; i<NumSameIDs; i++) {
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(i, maxNumEntries, NumEntries, values, int_Indices));
          Row = (int_type) GRID64(i);
    ierr = InsertOffsetValues(Row, NumEntries, values, Indexor->SameOffsets()[i]); 
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
      else {
        for (i=0; i<NumSameIDs; i++) {
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(i, maxNumEntries, NumEntries, values, int_Indices));
          Row = (int_type) GRID64(i);

          // convert to GIDs, start from right.
          for(j = NumEntries; j > 0;) {
            --j;
           gen_Indices[j] = (int_type) colMap.GID64(int_Indices[j]);
          }
    ierr = InsertGlobalValues(Row, NumEntries, values, gen_Indices); 
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
    }
  }
  
  // Do local permutation next
  if (NumPermuteIDs>0) {
    if (StaticGraph() || IndicesAreLocal()) {
      if( Indexor ) {
        for (i=0; i<NumPermuteIDs; i++) {
          FromRow = PermuteFromLIDs[i];
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, maxNumEntries, NumEntries, values, int_Indices));
          ToRow = (int_type) GRID64(PermuteToLIDs[i]);
    ierr = ReplaceOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]);
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
      else {
        for (i=0; i<NumPermuteIDs; i++) {
          FromRow = PermuteFromLIDs[i];
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, maxNumEntries, NumEntries, values, int_Indices));
          ToRow = (int_type) GRID64(PermuteToLIDs[i]);
          for(j=0; j<NumEntries; ++j) {
            int_Indices[j] = LCID((int_type) colMap.GID64(int_Indices[j]));
          }
    ierr = ReplaceMyValues((int) ToRow, NumEntries, values, int_Indices);
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
    }
    else {
      if( Indexor ) {
        for (i=0; i<NumPermuteIDs; i++) {
          FromRow = PermuteFromLIDs[i];
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, maxNumEntries, NumEntries, values, int_Indices));
          ToRow = (int_type) GRID64(PermuteToLIDs[i]);
    ierr = InsertOffsetValues(ToRow, NumEntries, values, Indexor->PermuteOffsets()[i]); 
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
      else {
        for (i=0; i<NumPermuteIDs; i++) {
          FromRow = PermuteFromLIDs[i];
          EPETRA_CHK_ERR(A.ExtractMyRowCopy(FromRow, maxNumEntries, NumEntries, values, int_Indices));

          // convert to GIDs, start from right.
          for(j = NumEntries; j > 0;) {
            --j;
           gen_Indices[j] = (int_type) colMap.GID64(int_Indices[j]);
          }

          ToRow = (int_type) GRID64(PermuteToLIDs[i]);
    ierr = InsertGlobalValues(ToRow, NumEntries, values, gen_Indices); 
          if (ierr<0) EPETRA_CHK_ERR(ierr);
        }
      }
    }
  }  

  if (maxNumEntries>0) {
    delete [] values;
    if(StaticGraph() || IndicesAreLocal() || Indexor) {
      delete [] int_Indices;
    }
    else {
      delete [] gen_Indices;
    }
  }
  return(0);
}

int Epetra_CrsMatrix::CopyAndPermuteRowMatrix(const Epetra_RowMatrix & A,
                                              int NumSameIDs, 
                int NumPermuteIDs,
                                              int * PermuteToLIDs,
                int *PermuteFromLIDs,
                                              const Epetra_OffsetIndex * Indexor )
{
  if(!A.Map().GlobalIndicesTypeMatch(RowMap()))
    throw ReportError("Epetra_CrsMatrix::CopyAndPermuteRowMatrix: Incoming global index type does not match the one for *this",-1);

  if(A.RowMatrixRowMap().GlobalIndicesInt())
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    return TCopyAndPermuteRowMatrix<int>(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs, Indexor);
#else
    throw ReportError("Epetra_CrsMatrix::CopyAndPermuteRowMatrix: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif

  if(A.RowMatrixRowMap().GlobalIndicesLongLong())
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    return TCopyAndPermuteRowMatrix<long long>(A, NumSameIDs, NumPermuteIDs, PermuteToLIDs, PermuteFromLIDs, Indexor);
#else
    throw ReportError("Epetra_CrsMatrix::CopyAndPermuteRowMatrix: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif

  throw ReportError("Epetra_CrsMatrix::CopyAndPermuteRowMatrix: Internal error, unable to determine global index type of maps", -1);
}

//=========================================================================
int Epetra_CrsMatrix::PackAndPrepare(const Epetra_SrcDistObject & Source, 
             int NumExportIDs,
                                     int * ExportLIDs,
             int & LenExports,
                                     char *& Exports,
             int & SizeOfPacket,
                                     int * Sizes,
                                     bool & VarSizes,
                                     Epetra_Distributor & Distor)
{
  if(!Source.Map().GlobalIndicesTypeMatch(RowMap()))
    throw ReportError("Epetra_CrsMatrix::PackAndPrepare: Incoming global index type does not match the one for *this",-1);

  (void)Distor;  
  // Rest of work can be done using RowMatrix only  
  const Epetra_RowMatrix & A = dynamic_cast<const Epetra_RowMatrix &>(Source);

  VarSizes = true; //enable variable block size data comm

  int TotalSendLength = 0;
  int * IntSizes = 0; 
  if( NumExportIDs>0 ) IntSizes = new int[NumExportIDs];

  int SizeofIntType = -1;
  if(Source.Map().GlobalIndicesInt())
    SizeofIntType = (int)sizeof(int); 
  else if(Source.Map().GlobalIndicesLongLong())
    SizeofIntType = (int)sizeof(long long); 
  else
    throw ReportError("Epetra_CrsMatrix::PackAndPrepare: Unable to determine source global index type",-1);

  for( int i = 0; i < NumExportIDs; ++i )
  {    
    int NumEntries;
    A.NumMyRowEntries( ExportLIDs[i], NumEntries );
    // Will have NumEntries doubles, NumEntries +2 ints, pack them interleaved     Sizes[i] = NumEntries;
    Sizes[i] = NumEntries;
    IntSizes[i] = 1 + (((NumEntries+2)*SizeofIntType)/(int)sizeof(double));
    TotalSendLength += (Sizes[i]+IntSizes[i]);
  }    
         
  double * DoubleExports = 0; 
  SizeOfPacket = (int)sizeof(double);
       
  //setup buffer locally for memory management by this object
  if( TotalSendLength*SizeOfPacket > LenExports )
  {
    if( LenExports > 0 ) delete [] Exports;
    LenExports = TotalSendLength*SizeOfPacket;
    DoubleExports = new double[TotalSendLength];
    for( int i = 0; i < TotalSendLength; ++i ) DoubleExports[i] = 0.0;
    Exports = (char *) DoubleExports;
  } 
 
  int NumEntries;
  double * values;
  double * valptr, * dintptr; 
 
  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumEntries, Number of indices in row.
  // next NumEntries: The actual indices for the row.
 
  const Epetra_Map & rowMap = A.RowMatrixRowMap();
  const Epetra_Map & colMap = A.RowMatrixColMap();
 
  if( NumExportIDs > 0 )
  {
    if(Source.Map().GlobalIndicesInt()) {
      int * Indices;
      int FromRow;
      int * intptr;

      int maxNumEntries = A.MaxNumEntries();
      dintptr = (double *) Exports;
      valptr = dintptr + IntSizes[0];
      intptr = (int *) dintptr;
      for (int i=0; i<NumExportIDs; i++)
      {
        FromRow = (int) rowMap.GID64(ExportLIDs[i]);
        intptr[0] = FromRow;
        values = valptr;
        Indices = intptr + 2;
        EPETRA_CHK_ERR(A.ExtractMyRowCopy(ExportLIDs[i], maxNumEntries, NumEntries, values, Indices));
        for (int j=0; j<NumEntries; j++) Indices[j] = (int) colMap.GID64(Indices[j]); // convert to GIDs
        intptr[1] = NumEntries; // Load second slot of segment
        if( i < (NumExportIDs-1) )
        {
          dintptr += (IntSizes[i]+Sizes[i]);
          valptr = dintptr + IntSizes[i+1];
          intptr = (int *) dintptr;
        }
      }
    }
    else if(Source.Map().GlobalIndicesLongLong()) {
      long long * LL_Indices;
      long long FromRow;
      long long * LLptr;

      int maxNumEntries = A.MaxNumEntries();
      dintptr = (double *) Exports;
      valptr = dintptr + IntSizes[0];
      LLptr = (long long *) dintptr;
      for (int i=0; i<NumExportIDs; i++)
      {
        FromRow = rowMap.GID64(ExportLIDs[i]);
        LLptr[0] = FromRow;
        values = valptr;
        LL_Indices = LLptr + 2;
        int * int_indices = reinterpret_cast<int*>(LL_Indices);
        EPETRA_CHK_ERR(A.ExtractMyRowCopy(ExportLIDs[i], maxNumEntries, NumEntries, values, int_indices));

        // convert to GIDs, start from right.
        for(int j = NumEntries; j > 0;) {
           --j;
           LL_Indices[j] = colMap.GID64(int_indices[j]);
        }

        LLptr[1] = NumEntries; // Load second slot of segment
        if( i < (NumExportIDs-1) )
        {
          dintptr += (IntSizes[i]+Sizes[i]);
          valptr = dintptr + IntSizes[i+1];
          LLptr = (long long *) dintptr;
        }
      }
    }

    for( int i = 0; i < NumExportIDs; ++i )
      Sizes[i] += IntSizes[i];
  }
 
  if( IntSizes ) delete [] IntSizes;
 
  return(0);
}

//=========================================================================
template<typename int_type>
int Epetra_CrsMatrix::TUnpackAndCombine(const Epetra_SrcDistObject & Source, 
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
  if (NumImportIDs<=0) return(0);
  
  if (   CombineMode != Add
   && CombineMode != Insert
   && CombineMode != Zero )
    EPETRA_CHK_ERR(-1); //Unsupported CombineMode, defaults to Zero

  int NumEntries;
  int_type * Indices;
  double * values;
  int_type ToRow;
  int i, ierr;
  int IntSize;
  
  double * valptr, *dintptr;
  int_type * intptr;

  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumEntries, Number of indices in row.
  // next NumEntries: The actual indices for the row.

  dintptr = (double *) Imports;
  intptr = (int_type *) dintptr;
  NumEntries = (int) intptr[1];
  IntSize = 1 + (((NumEntries+2)*(int)sizeof(int_type))/(int)sizeof(double));
  valptr = dintptr + IntSize;
 
  for (i=0; i<NumImportIDs; i++)
  {
    ToRow = (int_type) GRID64(ImportLIDs[i]);
    assert((intptr[0])==ToRow); // Sanity check
    values = valptr;
    Indices = intptr + 2;
 
    if (CombineMode==Add) {
      if (StaticGraph() || IndicesAreLocal()) {
        if( Indexor )
          ierr = SumIntoOffsetValues(ToRow, NumEntries, values, Indexor->RemoteOffsets()[i]);
        else
          ierr = SumIntoGlobalValues(ToRow, NumEntries, values, Indices);
      }
      else {
        if( Indexor )
          ierr = InsertOffsetValues(ToRow, NumEntries, values, Indexor->RemoteOffsets()[i]);
        else
          ierr = InsertGlobalValues(ToRow, NumEntries, values, Indices);
      }
      if (ierr<0) EPETRA_CHK_ERR(ierr);
    }
    else if (CombineMode==Insert) {
      if (StaticGraph() || IndicesAreLocal()) {
        if( Indexor )
          ierr = ReplaceOffsetValues(ToRow, NumEntries, values, Indexor->RemoteOffsets()[i]);
        else
          ierr = ReplaceGlobalValues(ToRow, NumEntries, values, Indices);
      }
      else {
        if( Indexor )
          ierr = InsertOffsetValues(ToRow, NumEntries, values, Indexor->RemoteOffsets()[i]);
        else
          ierr = InsertGlobalValues(ToRow, NumEntries, values, Indices);
      }
      if (ierr<0) EPETRA_CHK_ERR(ierr);
    }
 
    if( i < (NumImportIDs-1) )
    {
      dintptr += IntSize + NumEntries;
      intptr = (int_type *) dintptr;
      NumEntries = (int) intptr[1];
      IntSize = 1 + (((NumEntries+2)*(int)sizeof(int_type))/(int)sizeof(double));
      valptr = dintptr + IntSize;
    }
  }

  return(0);
}

int Epetra_CrsMatrix::UnpackAndCombine(const Epetra_SrcDistObject & Source, 
               int NumImportIDs,
                                       int * ImportLIDs, 
                                       int LenImports,
               char * Imports,
                                       int & SizeOfPacket, 
               Epetra_Distributor & Distor, 
               Epetra_CombineMode CombineMode,
                                       const Epetra_OffsetIndex * Indexor )
{
  if(!Source.Map().GlobalIndicesTypeMatch(RowMap()))
    throw ReportError("Epetra_CrsMatrix::UnpackAndCombine: Incoming global index type does not match the one for *this",-1);

  if(Source.Map().GlobalIndicesInt())
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    return TUnpackAndCombine<int>(Source, NumImportIDs, ImportLIDs, LenImports,
               Imports, SizeOfPacket, Distor, CombineMode, Indexor);
#else
    throw ReportError("Epetra_CrsMatrix::UnpackAndCombine: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif

  if(Source.Map().GlobalIndicesLongLong())
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    return TUnpackAndCombine<long long>(Source, NumImportIDs, ImportLIDs, LenImports,
               Imports, SizeOfPacket, Distor, CombineMode, Indexor);
#else
    throw ReportError("Epetra_CrsMatrix::UnpackAndCombine: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif

  throw ReportError("Epetra_CrsMatrix::UnpackAndCombine: Internal error, unable to determine global index type of maps", -1);
}

//=========================================================================

void Epetra_CrsMatrix::Print(ostream& os) const {
  int MyPID = RowMap().Comm().MyPID();
  int NumProc = RowMap().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      /*      const Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
        const Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
        const int             oldp = os.precision(12); */
      if (MyPID==0) {
  os <<  "\nNumber of Global Rows        = "; os << NumGlobalRows64(); os << endl;
  os <<    "Number of Global Cols        = "; os << NumGlobalCols64(); os << endl;
  os <<    "Number of Global Diagonals   = "; os << NumGlobalDiagonals64(); os << endl;
  os <<    "Number of Global Nonzeros    = "; os << NumGlobalNonzeros64(); os << endl;
  os <<    "Global Maximum Num Entries   = "; os << GlobalMaxNumEntries(); os << endl;
  if (LowerTriangular()) os <<    " ** Matrix is Lower Triangular **"; os << endl;
  if (UpperTriangular()) os <<    " ** Matrix is Upper Triangular **"; os << endl;
  if (NoDiagonal())      os <<    " ** Matrix has no diagonal     **"; os << endl; os << endl;
      }
      
      os <<  "\nNumber of My Rows        = "; os << NumMyRows(); os << endl;
      os <<    "Number of My Cols        = "; os << NumMyCols(); os << endl;
      os <<    "Number of My Diagonals   = "; os << NumMyDiagonals(); os << endl;
      os <<    "Number of My Nonzeros    = "; os << NumMyNonzeros(); os << endl;
      os <<    "My Maximum Num Entries   = "; os << MaxNumEntries(); os << endl; os << endl;

      os << flush;
      
      // Reset os flags
      
      /*      os.setf(olda,ios::adjustfield);
        os.setf(oldf,ios::floatfield);
        os.precision(oldp); */
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }
  
  {for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumMyRows1 = NumMyRows();
      int MaxNumIndices = MaxNumEntries();

      int * Indices_int = 0;
      long long * Indices_LL = 0;
      if(RowMap().GlobalIndicesInt()) {
         Indices_int = new int[MaxNumIndices];
      }
      else if(RowMap().GlobalIndicesLongLong()) {
         Indices_LL = new long long[MaxNumIndices];
      }
      else
         throw ReportError("Epetra_CrsGraph::Print: Unable to determine source global index type",-1);

      double * values  = new double[MaxNumIndices];
#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
      int NumIndices;
      int j;
#endif
      int i;
      
      if (MyPID==0) {
  os.width(8);
  os <<  "   Processor ";
  os.width(10);
  os <<  "   Row Index ";
  os.width(10);
  os <<  "   Col Index ";
  os.width(20);
  os <<  "   Value     ";
  os << endl;
      }
      for (i=0; i<NumMyRows1; i++) {
        if(RowMap().GlobalIndicesInt()) {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
           int Row = (int) GRID64(i); // Get global row number
           ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, values, Indices_int);

           for (j = 0; j < NumIndices ; j++) {   
              os.width(8);
              os <<  MyPID ; os << "    ";  
              os.width(10);
              os <<  Row ; os << "    ";  
              os.width(10);
              os <<  Indices_int[j]; os << "    ";
              os.width(20);
              os <<  values[j]; os << "    ";
              os << endl;
           }
#else
    throw ReportError("Epetra_CrsMatrix::Print: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
        }
        else if(RowMap().GlobalIndicesLongLong()) {
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
           long long Row = GRID64(i); // Get global row number
           ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, values, Indices_LL);

           for (j = 0; j < NumIndices ; j++) {   
              os.width(8);
              os <<  MyPID ; os << "    ";  
              os.width(10);
              os <<  Row ; os << "    ";  
              os.width(10);
              os <<  Indices_LL[j]; os << "    ";
              os.width(20);
              os <<  values[j]; os << "    ";
              os << endl;
           }
#else
    throw ReportError("Epetra_CrsMatrix::Print: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
        }
      }

      if(RowMap().GlobalIndicesInt()) {
         delete [] Indices_int;
      }
      else if(RowMap().GlobalIndicesLongLong()) {
         delete [] Indices_LL;
      }
      delete [] values;
      
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    RowMap().Comm().Barrier();
    RowMap().Comm().Barrier();
    RowMap().Comm().Barrier();
  }}
  
  return;
}
//=============================================================================
int Epetra_CrsMatrix::Multiply(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Multiply(TransA,x,y)");
#endif

  //
  // This function forms the product y = A * x or y = A' * x
  //

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  double* xp = (double*) x.Values();
  double* yp = (double*) y.Values();

  Epetra_Vector * xcopy = 0;
  if (&x==&y && Importer()==0 && Exporter()==0) {
    xcopy = new Epetra_Vector(x);
    xp = (double *) xcopy->Values();
  }
  UpdateImportVector(1); // Refresh import and output vectors if needed
  UpdateExportVector(1);

  if(!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if(Importer() != 0) {
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*) ImportVector_->Values();
    }
    
    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if(Exporter() != 0)  yp = (double*) ExportVector_->Values();
    
    // Do actual computation
    GeneralMV(xp, yp);

    if(Exporter() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
  }
  
  else { // Transpose operation

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    if(Exporter() != 0) {
      EPETRA_CHK_ERR(ExportVector_->Import(x, *Exporter(), Insert));
      xp = (double*) ExportVector_->Values();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if(Importer() != 0) yp = (double*) ImportVector_->Values();

    // Do actual computation
    GeneralMTV(xp, yp);

    if(Importer() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
  }


  UpdateFlops(2 * NumGlobalNonzeros64());
  if (xcopy!=0) {
    delete xcopy;
    EPETRA_CHK_ERR(1); // Return positive code to alert the user about needing extra copy of x
    return(1);
  }
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Multiply(TransA,X,Y)");
#endif

#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  if(Epetra_CrsMatrixTraceDumpMultiply) {
    *out << std::boolalpha;
    *out << "\nEntering Epetra_CrsMatrix::Multipy("<<TransA<<",X,Y) ...\n";
    if(!TransA) {
      *out << "\nDomainMap =\n";
      this->DomainMap().Print(Teuchos::OSTab(out).o());
    }
    else {
      *out << "\nRangeMap =\n";
      this->RangeMap().Print(Teuchos::OSTab(out).o());
    }
    *out << "\nInitial input X with " << ( TransA ? "RangeMap" : "DomainMap" ) << " =\n\n";
    X.Print(Teuchos::OSTab(out).o());
  }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY

  //
  // This function forms the product Y = A * Y or Y = A' * X
  //
  if(!Filled()) {
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  }

  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) {
    EPETRA_CHK_ERR(-2); // Need same number of vectors in each MV
  }

  double** Xp = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();

  int LDX = X.ConstantStride() ? X.Stride() : 0;
  int LDY = Y.ConstantStride() ? Y.Stride() : 0;

  Epetra_MultiVector * Xcopy = 0;
  if (&X==&Y && Importer()==0 && Exporter()==0) {
    Xcopy = new Epetra_MultiVector(X);
    Xp = (double **) Xcopy->Pointers();
    LDX = Xcopy->ConstantStride() ? Xcopy->Stride() : 0;
  }
  UpdateImportVector(NumVectors); // Make sure Import and Export Vectors are compatible
  UpdateExportVector(NumVectors);

  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
      Xp = (double**)ImportVector_->Pointers();
      LDX = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
      if(Epetra_CrsMatrixTraceDumpMultiply) {
        *out << "\nColMap =\n";
        this->ColMap().Print(Teuchos::OSTab(out).o());
        *out << "\nX after import from DomainMap to ColMap =\n\n";
        ImportVector_->Print(Teuchos::OSTab(out).o());
      }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      Yp = (double**)ExportVector_->Pointers();
      LDY = ExportVector_->ConstantStride() ? ExportVector_->Stride() : 0;
    }

    // Do actual computation
    if (NumVectors==1)
      GeneralMV(*Xp, *Yp);
    else
      GeneralMM(Xp, LDX, Yp, LDY, NumVectors);
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    if(Epetra_CrsMatrixTraceDumpMultiply) {
      *out << "\nRowMap =\n";
      this->RowMap().Print(Teuchos::OSTab(out).o());
      *out << "\nY after local mat-vec where Y has RowMap =\n\n";
      if(Exporter()!=0)
        ExportVector_->Print(Teuchos::OSTab(out).o());
      else
        Y.Print(Teuchos::OSTab(out).o());
    }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    if (Exporter()!=0) {
      Y.PutScalar(0.0);  // Make sure target is zero
      Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
      if(Epetra_CrsMatrixTraceDumpMultiply) {
        *out << "\nRangeMap =\n";
        this->RangeMap().Print(Teuchos::OSTab(out).o());
        *out << "\nY after export from RowMap to RangeMap = \n\n";
        Y.Print(Teuchos::OSTab(out).o());
      }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());
  }
  else { // Transpose operation

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      EPETRA_CHK_ERR(ExportVector_->Import(X, *Exporter(), Insert));
      Xp = (double**)ExportVector_->Pointers();
      LDX = ExportVector_->ConstantStride() ? ExportVector_->Stride() : 0;
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
      if(Epetra_CrsMatrixTraceDumpMultiply) {
        *out << "\nRowMap =\n";
        this->RowMap().Print(Teuchos::OSTab(out).o());
        *out << "\nX after import from RangeMap to RowMap =\n\n";
        ExportVector_->Print(Teuchos::OSTab(out).o());
      }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      Yp = (double**)ImportVector_->Pointers();
      LDY = ImportVector_->ConstantStride() ? ImportVector_->Stride() : 0;
    }

    // Do actual computation
    if (NumVectors==1)
      GeneralMTV(*Xp, *Yp);
    else
      GeneralMTM(Xp, LDX, Yp, LDY, NumVectors);
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    if(Epetra_CrsMatrixTraceDumpMultiply) {
      *out << "\nColMap =\n";
      this->ColMap().Print(Teuchos::OSTab(out).o());
      *out << "\nY after local transpose mat-vec where Y has ColMap =\n\n";
      if(Importer()!=0)
        ImportVector_->Print(Teuchos::OSTab(out).o());
      else
        Y.Print(Teuchos::OSTab(out).o());
    }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    if (Importer()!=0) {
      Y.PutScalar(0.0);  // Make sure target is zero
      EPETRA_CHK_ERR(Y.Export(*ImportVector_, *Importer(), Add)); // Fill Y with Values from export vector
#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
      if(Epetra_CrsMatrixTraceDumpMultiply) {
        *out << "\nDomainMap =\n";
        this->DomainMap().Print(Teuchos::OSTab(out).o());
        *out << "\nY after export from ColMap to DomainMap =\n\n";
        Y.Print(Teuchos::OSTab(out).o());
      }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1)  EPETRA_CHK_ERR(Y.Reduce());
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros64());
  if (Xcopy!=0) {
    delete Xcopy;
    EPETRA_CHK_ERR(1); // Return positive code to alert the user about needing extra copy of X
    return(1);
  }

#ifdef EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY
  if(Epetra_CrsMatrixTraceDumpMultiply) {
    *out << "\nFinal output Y is the last Y printed above!\n";
    *out << "\nLeaving Epetra_CrsMatrix::Multipy("<<TransA<<",X,Y) ...\n";
  }
#endif // EPETRA_CRS_MATRIX_TRACE_DUMP_MULTIPLY

  return(0);
}
//=======================================================================================================
void Epetra_CrsMatrix::UpdateImportVector(int NumVectors) const {    
  if(Importer() != 0) {
    if(ImportVector_ != 0) {
      if(ImportVector_->NumVectors() != NumVectors) { 
  delete ImportVector_; 
  ImportVector_= 0;
      }
    }
    if(ImportVector_ == 0) 
      ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors); // Create import vector if needed
  }
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::UpdateExportVector(int NumVectors) const {    
  if(Exporter() != 0) {
    if(ExportVector_ != 0) {
      if(ExportVector_->NumVectors() != NumVectors) { 
  delete ExportVector_; 
  ExportVector_= 0;
      }
    }
    if(ExportVector_ == 0) 
      ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
  }
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralMV(double * x, double * y)  const {
  
if (StorageOptimized() && Graph().StorageOptimized()) {

  double * values = All_Values();
  int * Indices = Graph().All_Indices();
  int * IndexOffset = Graph().IndexOffset();
#ifdef EPETRA_HAVE_OMP
  const int numMyRows = NumMyRows_;
#pragma omp parallel for default(none) shared(IndexOffset,values,Indices,y,x)
     for (int row=0; row<numMyRows; ++row)
        {
     const int curOffset = IndexOffset[row];
          const double *val_ptr    = values+curOffset;
          const int    *colnum_ptr = Indices+curOffset;
     double s = 0.;
     const double *const val_end_of_row = &values[IndexOffset[row+1]];
     while (val_ptr != val_end_of_row)
       s += *val_ptr++ * x[*colnum_ptr++];
     y[row] = s;
  }
#else
#ifndef Epetra_ENABLE_CASK
       const double *val_ptr    = values;
       const int    *colnum_ptr = Indices;
       double       * dst_ptr = y;
       for (int row=0; row<NumMyRows_; ++row)
   {
     double s = 0.;
     const double *const val_end_of_row = &values[IndexOffset[row+1]];
     while (val_ptr != val_end_of_row)
       s += *val_ptr++ * x[*colnum_ptr++];
     *dst_ptr++ = s;
   }
#else
       cask_csr_dax_new(NumMyRows_, IndexOffset, Indices,
                        values, x, y, cask);
#endif // Epetra_ENABLE_CASK
#endif // EPETRA_HAVE_OMP

    return;
  }
  else if (!StorageOptimized() && !Graph().StorageOptimized()) {


    int* NumEntriesPerRow = Graph().NumIndicesPerRow();
    int** Indices = Graph().Indices();
    double** srcValues = Values();
        const int numMyRows = NumMyRows_;

    // Do actual computation
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(NumEntriesPerRow,Indices,srcValues,y,x)
#endif
    for(int i = 0; i < numMyRows; i++) {
      int     NumEntries = NumEntriesPerRow[i];
      int*    RowIndices = Indices[i];
      double* RowValues  = srcValues[i];
      double sum = 0.0;
      for(int j = 0; j < NumEntries; j++) 
  sum += *RowValues++ * x[*RowIndices++];
      
      y[i] = sum;
      
    }
  } 
  else { // Case where StorageOptimized is incompatible:  Use general accessors.

    const int numMyRows = NumMyRows_;

    // Do actual computation
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(x,y)
#endif
    for(int i = 0; i < numMyRows; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      double sum = 0.0;
      for(int j = 0; j < NumEntries; j++) 
  sum += *RowValues++ * x[*RowIndices++];
      
      y[i] = sum;
      
    }
  } 
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralMTV(double * x, double * y) const {

  int NumCols = NumMyCols();
#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
  if (StorageOptimized() && Graph().StorageOptimized()) {
    double * values = All_Values_;
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();
#ifndef Epetra_ENABLE_CASK
   int ione = 1;
   EPETRA_DCRSMV_F77(&ione, &NumMyRows_, &NumCols, values, Indices, IndexOffset, x, y);
#else
   cask_csr_datx( NumMyRows_, NumCols, IndexOffset,  Indices,  values,x ,y );
#endif

    return;
  }
#endif // FORTRAN_DISABLED
  for(int i = 0; i < NumCols; i++) 
    y[i] = 0.0; // Initialize y for transpose multiply

  if (StorageOptimized() && Graph().StorageOptimized()) {
    double * values = All_Values_;
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();
    for(int i = 0; i < NumMyRows_; ++i) {
      int prevOffset = *IndexOffset++;
      int NumEntries = *IndexOffset - prevOffset;
      double xi = x[i];
      for(int j = 0; j < NumEntries; j++) 
  y[*Indices++] += *values++ * xi;
    }
  }
  else if (!StorageOptimized() && !Graph().StorageOptimized()) {

    int* NumEntriesPerRow = Graph().NumIndicesPerRow();
    int** Indices = Graph().Indices();
    double** srcValues = Values();
    
    for(int i = 0; i < NumMyRows_; i++) {
      int     NumEntries = *NumEntriesPerRow++;
      int*    RowIndices = *Indices++;
      double* RowValues  = *srcValues++;
      double xi = x[i];
      for(int j = 0; j < NumEntries; j++) 
  y[*RowIndices++] += *RowValues++ * xi;
    }
  }
  else { // Case where StorageOptimized is incompatible:  Use general accessors.
  
    for(int i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      double xi = x[i];
      for(int j = 0; j < NumEntries; j++) 
  y[*RowIndices++] += *RowValues++ * xi;
    }
  }

  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralMM(double ** X, int LDX, double ** Y, int LDY, int NumVectors) const {

#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
  if (StorageOptimized() && Graph().StorageOptimized()) {
    double * values = All_Values_;
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();

    if (LDX!=0 && LDY!=0) {
#ifndef Epetra_ENABLE_CASK
    int izero = 0;
    EPETRA_DCRSMM_F77(&izero, &NumMyRows_, &NumMyRows_, values, Indices, IndexOffset, *X, &LDX, *Y, &LDY, &NumVectors);
#else
    cask_csr_dgesmm_new(0, 1.0, NumMyRows_, NumMyRows_,  NumVectors,
                    IndexOffset, Indices, values, *X, LDX, 0.0,  *Y, LDY,cask);
#endif
    return;
    }

    double ** const xp = X;
    double ** const yp = Y;
    const int numMyRows = NumMyRows_;
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(IndexOffset,Indices,values,NumVectors)
#endif
    for (int i=0; i < numMyRows; i++) {
      int prevOffset = IndexOffset[i];
      int NumEntries = IndexOffset[i+1] - prevOffset;
      int *    RowIndices = Indices+prevOffset;
      double * RowValues  = values+prevOffset;
      for (int k=0; k<NumVectors; k++) {
  double sum = 0.0;
  const double * const x = xp[k];
  double * const y = yp[k];
  for (int j=0; j < NumEntries; j++) sum += RowValues[j] * x[RowIndices[j]];
  y[i] = sum;
      }
    }
  }
  else if (!StorageOptimized() && !Graph().StorageOptimized()) {
#else
  if (!StorageOptimized() && !Graph().StorageOptimized()) {
#endif

    int* NumEntriesPerRow = Graph().NumIndicesPerRow();
    int** Indices = Graph().Indices();
    double** srcValues = Values();
    double ** const xp = X;
    double ** const yp = Y;
    const int numMyRows = NumMyRows_;

#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(NumEntriesPerRow,Indices,srcValues,NumVectors)
#endif
    for (int i=0; i < numMyRows; i++) {
      int      NumEntries = NumEntriesPerRow[i];
      int *    RowIndices = Indices[i];
      double * RowValues  = srcValues[i];
      for (int k=0; k<NumVectors; k++) {
  double sum = 0.0;
  const double * const x = xp[k];
  double * const y = yp[k];
  for (int j=0; j < NumEntries; j++) sum += RowValues[j] * x[RowIndices[j]];
  y[i] = sum;
      }
    }
  }
  else {

    double ** const xp = X;
    double ** const yp = Y;
    const int numMyRows = NumMyRows_;
#ifdef EPETRA_HAVE_OMP
#pragma omp parallel for default(none) shared(NumVectors)
#endif
    for (int i=0; i < numMyRows; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for (int k=0; k<NumVectors; k++) {
  double sum = 0.0;
  double * x = xp[k];
  for (int j=0; j < NumEntries; j++) sum += RowValues[j] * x[RowIndices[j]];
  yp[k][i] = sum;
      }
    }
  }
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralMTM(double ** X, int LDX, double ** Y, int LDY, int NumVectors)  const{

  int NumCols = NumMyCols();
#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
  if (StorageOptimized() && Graph().StorageOptimized()) {
    if (LDX!=0 && LDY!=0) {
      double * values = All_Values_;
      int * Indices = Graph().All_Indices();
      int * IndexOffset = Graph().IndexOffset();

#ifndef Epetra_ENABLE_CASK
      int ione = 1;
      EPETRA_DCRSMM_F77(&ione, &NumMyRows_, &NumCols, values, Indices, IndexOffset, *X, &LDX, *Y, &LDY, &NumVectors);
#else
      cask_csr_dgesmm_new(1, 1.0, NumMyRows_, NumCols,  NumVectors,
                          IndexOffset, Indices, values, *X, LDX, 0.0,
                          *Y, LDY, cask);
#endif
      return;
    }
  }
#endif
  for (int k=0; k<NumVectors; k++) 
    for (int i=0; i < NumCols; i++) 
      Y[k][i] = 0.0; // Initialize y for transpose multiply
  
  if (StorageOptimized() && Graph().StorageOptimized()) {
    double * values = All_Values_;
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();
    for (int i=0; i < NumMyRows_; i++) {
      int prevOffset = *IndexOffset++;
      int NumEntries = *IndexOffset - prevOffset;
      int *    RowIndices = Indices+prevOffset;
      double * RowValues  = values+prevOffset;
      
      for (int k=0; k<NumVectors; k++) {
  double * y = Y[k];
  double * x = X[k];
  for (int j=0; j < NumEntries; j++) 
    y[RowIndices[j]] += RowValues[j] * x[i];
      }
    }
  }
  else if (!StorageOptimized() && !Graph().StorageOptimized()) {
    
    int* NumEntriesPerRow = Graph().NumIndicesPerRow();
    int** Indices = Graph().Indices();
    double** srcValues = Values();

    for (int i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *srcValues++;
      for (int k=0; k<NumVectors; k++) {
  double * y = Y[k];
  double * x = X[k];
  for (int j=0; j < NumEntries; j++) 
    y[RowIndices[j]] += RowValues[j] * x[i];
      }
    }
  }
  else { // Case where StorageOptimized is incompatible:  Use general accessors.
    
    for (int i=0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for (int k=0; k<NumVectors; k++) {
  double * y = Y[k];
  double * x = X[k];
  for (int j=0; j < NumEntries; j++) 
    y[RowIndices[j]] += RowValues[j] * x[i];
      }
    }
  }
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralSV(bool Upper, bool Trans, bool UnitDiagonal, double * xp, double * yp)  const {


  int i, j, j0;

#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
  if (StorageOptimized() && Graph().StorageOptimized() && ((UnitDiagonal && NoDiagonal())|| (!UnitDiagonal && !NoDiagonal()))) {
    double * values = All_Values();
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();

    int iupper = Upper ? 1:0;
    int itrans = Trans ? 1:0;
    int udiag =  UnitDiagonal ? 1:0;
    int nodiag = NoDiagonal() ? 1:0;
    int xysame = (xp==yp) ? 1:0;

#ifndef Epetra_ENABLE_CASK
    EPETRA_DCRSSV_F77( &iupper, &itrans, &udiag, &nodiag, &NumMyRows_, &NumMyRows_, values, Indices, IndexOffset, xp, yp, &xysame);
#else
    cask_csr_dtrsv_new( iupper, itrans, udiag, nodiag, 0, xysame, NumMyRows_,
                    NumMyRows_, IndexOffset, Indices, values, xp, yp, cask);
#endif
    return;
  }
  //=================================================================
  else { // !StorageOptimized()
  //=================================================================
#endif    
    if (!Trans) {
      
      if (Upper) {
  
  j0 = 1;
  if (NoDiagonal()) 
    j0--; // Include first term if no diagonal
  for (i=NumMyRows_-1; i >=0; i--) {
    int      NumEntries = NumMyEntries(i);
    int *    RowIndices = Graph().Indices(i);
    double * RowValues  = Values(i);
    double sum = 0.0;
    for (j=j0; j < NumEntries; j++) 
      sum += RowValues[j] * yp[RowIndices[j]];
    
    if (UnitDiagonal) 
      yp[i] = xp[i] - sum;
    else 
      yp[i] = (xp[i] - sum)/RowValues[0];
    
  }
      }
      else {
  j0 = 1;
  if (NoDiagonal())
    j0--; // Include first term if no diagonal
  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = NumMyEntries(i) - j0;
    int *    RowIndices = Graph().Indices(i);
    double * RowValues  = Values(i);
    double sum = 0.0;
    for (j=0; j < NumEntries; j++) 
      sum += RowValues[j] * yp[RowIndices[j]];
    
    if (UnitDiagonal) 
      yp[i] = xp[i] - sum;
    else 
      yp[i] = (xp[i] - sum)/RowValues[NumEntries];
    
  }
      }
    }
    
    // ***********  Transpose case *******************************
    
    else {
      
      if (xp!=yp) 
  for (i=0; i < NumMyRows_; i++) 
    yp[i] = xp[i]; // Initialize y for transpose solve
      
      if (Upper) {
  
  j0 = 1;
  if (NoDiagonal()) 
    j0--; // Include first term if no diagonal
  
  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = NumMyEntries(i);
    int *    RowIndices = Graph().Indices(i);
    double * RowValues  = Values(i);
    if (!UnitDiagonal) 
      yp[i] = yp[i]/RowValues[0];
    double ytmp = yp[i];
    for (j=j0; j < NumEntries; j++) 
      yp[RowIndices[j]] -= RowValues[j] * ytmp;
  }
      }
      else {
  
  j0 = 1;
  if (NoDiagonal()) 
    j0--; // Include first term if no diagonal
  
  for (i=NumMyRows_-1; i >= 0; i--) {
    int      NumEntries = NumMyEntries(i) - j0;
    int *    RowIndices = Graph().Indices(i);
    double * RowValues  = Values(i);
    if (!UnitDiagonal) 
      yp[i] = yp[i]/RowValues[NumEntries];
    double ytmp = yp[i];
    for (j=0; j < NumEntries; j++) 
      yp[RowIndices[j]] -= RowValues[j] * ytmp;
  }
      }
      
    }
#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
  }
#endif
  return;
}
//=======================================================================================================
void Epetra_CrsMatrix::GeneralSM(bool Upper, bool Trans, bool UnitDiagonal, double ** Xp, int LDX, double ** Yp, int LDY, int NumVectors)  const{

  int i, j, j0, k;
  double diag = 0.0;

  if (StorageOptimized() && Graph().StorageOptimized()) {
    double * values = All_Values();
    int * Indices = Graph().All_Indices();
    int * IndexOffset = Graph().IndexOffset();
#if !defined(FORTRAN_DISABLED) || defined(Epetra_ENABLE_CASK)
    if (LDX!=0 && LDY!=0 && ((UnitDiagonal && NoDiagonal()) || (!UnitDiagonal && !NoDiagonal()))) {
      int iupper = Upper ? 1:0;
      int itrans = Trans ? 1:0;
      int udiag =  UnitDiagonal ? 1:0;
      int nodiag = NoDiagonal() ? 1:0;
      int xysame = (Xp==Yp) ? 1:0;

#ifndef Epetra_ENABLE_CASK
      EPETRA_DCRSSM_F77( &iupper, &itrans, &udiag, &nodiag, &NumMyRows_, &NumMyRows_, values, Indices, IndexOffset, 
       *Xp, &LDX, *Yp, &LDY, &xysame, &NumVectors);
#else
      cask_csr_dtrsm( iupper, itrans, udiag, nodiag, 0, xysame,  NumMyRows_, 
                      NumMyRows_, NumVectors, IndexOffset, Indices, values, 
                      *Xp, LDX, *Yp, LDY);
#endif
      return;
    }
#endif
    if(!Trans) {   
      if(Upper) {   
  j0 = 1;
  if(NoDiagonal()) 
    j0--; // Include first term if no diagonal
  for(i = NumMyRows_ - 1; i >= 0; i--) {
    int Offset = IndexOffset[i];
    int      NumEntries = IndexOffset[i+1]-Offset;
    int *    RowIndices = Indices+Offset;
    double * RowValues  = values+Offset;
    if(!UnitDiagonal) 
      diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
    for(k = 0; k < NumVectors; k++) {
      double sum = 0.0;
      for(j = j0; j < NumEntries; j++) 
        sum += RowValues[j] * Yp[k][RowIndices[j]];
          
      if(UnitDiagonal) 
        Yp[k][i] = Xp[k][i] - sum;
      else 
        Yp[k][i] = (Xp[k][i] - sum) * diag;
    }
  }
      }
      else {
  j0 = 1;
  if(NoDiagonal()) 
    j0--; // Include first term if no diagonal
  for(i = 0; i < NumMyRows_; i++) {
    int Offset = IndexOffset[i];
    int      NumEntries = IndexOffset[i+1]-Offset - j0;
    int *    RowIndices = Indices+Offset;
    double * RowValues  = values+Offset;
    if(!UnitDiagonal)
      diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
    for(k = 0; k < NumVectors; k++) {
      double sum = 0.0;
      for(j = 0; j < NumEntries; j++) 
        sum += RowValues[j] * Yp[k][RowIndices[j]];
          
      if(UnitDiagonal) 
        Yp[k][i] = Xp[k][i] - sum;
      else 
        Yp[k][i] = (Xp[k][i] - sum)*diag;
    }
  }
      }
    }
    // ***********  Transpose case *******************************

    else {
      for(k = 0; k < NumVectors; k++) 
  if(Yp[k] != Xp[k]) 
    for(i = 0; i < NumMyRows_; i++)
      Yp[k][i] = Xp[k][i]; // Initialize y for transpose multiply
    
      if(Upper) {
  j0 = 1;
  if(NoDiagonal()) 
    j0--; // Include first term if no diagonal
      
  for(i = 0; i < NumMyRows_; i++) {
    int Offset = IndexOffset[i];
    int      NumEntries = IndexOffset[i+1]-Offset;
    int *    RowIndices = Indices+Offset;
    double * RowValues  = values+Offset;
    if(!UnitDiagonal) 
      diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
    for(k = 0; k < NumVectors; k++) {
      if(!UnitDiagonal) 
        Yp[k][i] = Yp[k][i]*diag;
      double ytmp = Yp[k][i];
      for(j = j0; j < NumEntries; j++) 
        Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
    }
  }
      }
      else {
  j0 = 1;
  if(NoDiagonal()) 
    j0--; // Include first term if no diagonal  
  for(i = NumMyRows_ - 1; i >= 0; i--) {
    int Offset = IndexOffset[i];
    int      NumEntries = IndexOffset[i+1]-Offset - j0;
    int *    RowIndices = Indices+Offset;
    double * RowValues  = values+Offset;
    if(!UnitDiagonal) 
      diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
    for(k = 0; k < NumVectors; k++) {
      if(!UnitDiagonal)  
        Yp[k][i] = Yp[k][i]*diag;
      double ytmp = Yp[k][i];
      for(j = 0; j < NumEntries; j++)
        Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
    }
  }
      }
    }
  }
    // ========================================================
  else { // !StorageOptimized()
    // ========================================================

    if(!Trans) {   
      if(Upper) {   
  j0 = 1;
  if(NoDiagonal()) 
    j0--; // Include first term if no diagonal
  for(i = NumMyRows_ - 1; i >= 0; i--) {
    int     NumEntries = NumMyEntries(i);
    int*    RowIndices = Graph().Indices(i);
    double* RowValues  = Values(i);
    if(!UnitDiagonal) 
      diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
    for(k = 0; k < NumVectors; k++) {
      double sum = 0.0;
      for(j = j0; j < NumEntries; j++) 
        sum += RowValues[j] * Yp[k][RowIndices[j]];
          
      if(UnitDiagonal) 
        Yp[k][i] = Xp[k][i] - sum;
      else 
        Yp[k][i] = (Xp[k][i] - sum) * diag;
    }
  }
      }
      else {
  j0 = 1;
  if(NoDiagonal()) 
    j0--; // Include first term if no diagonal
  for(i = 0; i < NumMyRows_; i++) {
    int     NumEntries = NumMyEntries(i) - j0;
    int*    RowIndices = Graph().Indices(i);
    double* RowValues  = Values(i);
    if(!UnitDiagonal)
      diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
    for(k = 0; k < NumVectors; k++) {
      double sum = 0.0;
      for(j = 0; j < NumEntries; j++) 
        sum += RowValues[j] * Yp[k][RowIndices[j]];
          
      if(UnitDiagonal) 
        Yp[k][i] = Xp[k][i] - sum;
      else 
        Yp[k][i] = (Xp[k][i] - sum)*diag;
    }
  }
      }
    }
    // ***********  Transpose case *******************************

    else {
      for(k = 0; k < NumVectors; k++) 
  if(Yp[k] != Xp[k]) 
    for(i = 0; i < NumMyRows_; i++)
      Yp[k][i] = Xp[k][i]; // Initialize y for transpose multiply
    
      if(Upper) {
  j0 = 1;
  if(NoDiagonal()) 
    j0--; // Include first term if no diagonal
      
  for(i = 0; i < NumMyRows_; i++) {
    int     NumEntries = NumMyEntries(i);
    int*    RowIndices = Graph().Indices(i);
    double* RowValues  = Values(i);
    if(!UnitDiagonal) 
      diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
    for(k = 0; k < NumVectors; k++) {
      if(!UnitDiagonal) 
        Yp[k][i] = Yp[k][i]*diag;
      double ytmp = Yp[k][i];
      for(j = j0; j < NumEntries; j++) 
        Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
    }
  }
      }
      else {
  j0 = 1;
  if(NoDiagonal()) 
    j0--; // Include first term if no diagonal  
  for(i = NumMyRows_ - 1; i >= 0; i--) {
    int     NumEntries = NumMyEntries(i) - j0;
    int*    RowIndices = Graph().Indices(i);
    double* RowValues  = Values(i);
    if(!UnitDiagonal) 
      diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
    for(k = 0; k < NumVectors; k++) {
      if(!UnitDiagonal)  
        Yp[k][i] = Yp[k][i]*diag;
      double ytmp = Yp[k][i];
      for(j = 0; j < NumEntries; j++)
        Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
    }
  }
      }
    }
  }
  
  return;
}
//=============================================================================
int Epetra_CrsMatrix::Multiply1(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Multiply1(TransA,x,y)");
#endif

  //
  // This function forms the product y = A * x or y = A' * x
  //

  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  int i, j;
  double* xp = (double*) x.Values();
  double* yp = (double*) y.Values();
  int NumMyCols_ = NumMyCols();

  if(!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if(Importer() != 0) {
      if(ImportVector_ != 0) {
  if(ImportVector_->NumVectors() != 1) { 
    delete ImportVector_; 
    ImportVector_= 0;
  }
      }
      if(ImportVector_ == 0) 
  ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      EPETRA_CHK_ERR(ImportVector_->Import(x, *Importer(), Insert));
      xp = (double*) ImportVector_->Values();
    }
    
    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if(Exporter() != 0) {
      if(ExportVector_ != 0) {
  if(ExportVector_->NumVectors() != 1) { 
    delete ExportVector_; 
    ExportVector_= 0;
  }
      }
      if(ExportVector_ == 0) 
  ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      yp = (double*) ExportVector_->Values();
    }
    
    // Do actual computation
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      double sum = 0.0;
      for(j = 0; j < NumEntries; j++) 
  sum += RowValues[j] * xp[RowIndices[j]];
      
      yp[i] = sum;
      
    }
    if(Exporter() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ExportVector_, *Exporter(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
  }
  
  else { // Transpose operation

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    if(Exporter() != 0) {
      if(ExportVector_ != 0) {
  if(ExportVector_->NumVectors() != 1) { 
    delete ExportVector_; 
    ExportVector_= 0;
  }
      }
      if(ExportVector_ == 0) 
  ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      EPETRA_CHK_ERR(ExportVector_->Import(x, *Exporter(), Insert));
      xp = (double*) ExportVector_->Values();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if(Importer() != 0) {
      if(ImportVector_ != 0) {
  if(ImportVector_->NumVectors() != 1) { 
    delete ImportVector_; 
    ImportVector_= 0;
  }
      }
      if(ImportVector_ == 0) 
  ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      yp = (double*) ImportVector_->Values();
    }

    // Do actual computation
    for(i = 0; i < NumMyCols_; i++) 
      yp[i] = 0.0; // Initialize y for transpose multiply
        
    for(i = 0; i < NumMyRows_; i++) {
      int     NumEntries = NumMyEntries(i);
      int*    RowIndices = Graph().Indices(i);
      double* RowValues  = Values(i);
      for(j = 0; j < NumEntries; j++) 
  yp[RowIndices[j]] += RowValues[j] * xp[i];
    }
    if(Importer() != 0) {
      y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(y.Export(*ImportVector_, *Importer(), Add)); // Fill y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(y.Reduce());
  }

  UpdateFlops(2 * NumGlobalNonzeros64());
  return(0);
}
//=============================================================================
int Epetra_CrsMatrix::Multiply1(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Multiply1(TransA,X,Y)");
#endif

  //
  // This function forms the product Y = A * Y or Y = A' * X
  //
  if((X.NumVectors() == 1) && (Y.NumVectors() == 1)) {
    double* xp = (double*) X[0];
    double* yp = (double*) Y[0];
    Epetra_Vector x(View, X.Map(), xp);
    Epetra_Vector y(View, Y.Map(), yp);
    EPETRA_CHK_ERR(Multiply1(TransA, x, y));
    return(0);
  }
  if(!Filled()) {
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  }

  int i, j, k;

  double** Xp = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();

  int NumVectors = X.NumVectors();
  int NumMyCols_ = NumMyCols();


  // Need to better manage the Import and Export vectors:
  // - Need accessor functions
  // - Need to make the NumVector match (use a View to do this)
  // - Need to look at RightScale and ColSum routines too.

  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
  if (ImportVector_->NumVectors()!=NumVectors) { 
    delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) 
  ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors); // Create import vector if needed
      EPETRA_CHK_ERR(ImportVector_->Import(X, *Importer(), Insert));
      Xp = (double**)ImportVector_->Pointers();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
  if (ExportVector_->NumVectors()!=NumVectors) { 
    delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) 
  ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
      Yp = (double**)ExportVector_->Pointers();
    }

    // Do actual computation

    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = NumMyEntries(i);
      int *    RowIndices = Graph().Indices(i);
      double * RowValues  = Values(i);
      for (k=0; k<NumVectors; k++) {
  double sum = 0.0;
  for (j=0; j < NumEntries; j++) sum += RowValues[j] * Xp[k][RowIndices[j]];
  Yp[k][i] = sum;
      }
    }
    if (Exporter()!=0) {
      Y.PutScalar(0.0); // Make sure target is zero
      Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().RangeMap().DistributedGlobal() && Comm().NumProc()>1) EPETRA_CHK_ERR(Y.Reduce());
  }
  else { // Transpose operation
    

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      if (ExportVector_!=0) {
  if (ExportVector_->NumVectors()!=NumVectors) { 
    delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) 
  ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
      EPETRA_CHK_ERR(ExportVector_->Import(X, *Exporter(), Insert));
      Xp = (double**)ExportVector_->Pointers();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
  if (ImportVector_->NumVectors()!=NumVectors) { 
    delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) 
  ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors); // Create import vector if needed
      Yp = (double**)ImportVector_->Pointers();
    }

    // Do actual computation



    for (k=0; k<NumVectors; k++) 
      for (i=0; i < NumMyCols_; i++) 
  Yp[k][i] = 0.0; // Initialize y for transpose multiply
    
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = NumMyEntries(i);
      int *    RowIndices = Graph().Indices(i);
      double * RowValues  = Values(i);
      for (k=0; k<NumVectors; k++) {
  for (j=0; j < NumEntries; j++) 
    Yp[k][RowIndices[j]] += RowValues[j] * Xp[k][i];
      }
    }
    if (Importer()!=0) {
      Y.PutScalar(0.0); // Make sure target is zero
      EPETRA_CHK_ERR(Y.Export(*ImportVector_, *Importer(), Add)); // Fill Y with Values from export vector
    }
    // Handle case of rangemap being a local replicated map
    if (!Graph().DomainMap().DistributedGlobal() && Comm().NumProc()>1)  EPETRA_CHK_ERR(Y.Reduce());
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros64());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Solve1(bool Upper, bool Trans, bool UnitDiagonal,
          const Epetra_Vector& x, Epetra_Vector& y) const
{

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Solve1(Upper,Trans,UnitDiag,x,y)");
#endif

  //
  // This function finds y such that Ly = x or Uy = x or the transpose cases.
  //

  if (!Filled()) {
    EPETRA_CHK_ERR(-1); // Matrix must be filled.
  }

  if ((Upper) && (!UpperTriangular())) 
    EPETRA_CHK_ERR(-2);
  if ((!Upper) && (!LowerTriangular())) 
    EPETRA_CHK_ERR(-3);
  if ((!UnitDiagonal) && (NoDiagonal())) 
    EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) 
    EPETRA_CHK_ERR(-5); // Need each row to have a diagonal
      

  int i, j, j0;
  int * NumEntriesPerRow = Graph().NumIndicesPerRow();
  int ** Indices = Graph().Indices();
  double ** Vals = Values();
  int NumMyCols_ = NumMyCols();

  // If upper, point to last row
  if ((Upper && !Trans) || (!Upper && Trans)) {
    NumEntriesPerRow += NumMyRows_-1;
    Indices += NumMyRows_-1;
    Vals += NumMyRows_-1;
  }
    
  double *xp = (double*)x.Values();
  double *yp = (double*)y.Values();

  if (!Trans) {

    if (Upper) {

      j0 = 1;
      if (NoDiagonal()) 
  j0--; // Include first term if no diagonal
      for (i=NumMyRows_-1; i >=0; i--) {
  int      NumEntries = *NumEntriesPerRow--;
  int *    RowIndices = *Indices--;
  double * RowValues  = *Vals--;
  double sum = 0.0;
  for (j=j0; j < NumEntries; j++) 
    sum += RowValues[j] * yp[RowIndices[j]];
        
  if (UnitDiagonal) 
    yp[i] = xp[i] - sum;
  else 
    yp[i] = (xp[i] - sum)/RowValues[0];
        
      }
    }
    else {
      j0 = 1;
      if (NoDiagonal())
  j0--; // Include first term if no diagonal
      for (i=0; i < NumMyRows_; i++) {
  int      NumEntries = *NumEntriesPerRow++ - j0;
  int *    RowIndices = *Indices++;
  double * RowValues  = *Vals++;
  double sum = 0.0;
  for (j=0; j < NumEntries; j++) 
    sum += RowValues[j] * yp[RowIndices[j]];
        
  if (UnitDiagonal) 
    yp[i] = xp[i] - sum;
  else 
    yp[i] = (xp[i] - sum)/RowValues[NumEntries];
        
      }
    }
  }
  
  // ***********  Transpose case *******************************
  
  else {

    if (xp!=yp) 
      for (i=0; i < NumMyCols_; i++) 
  yp[i] = xp[i]; // Initialize y for transpose solve
    
    if (Upper) {

      j0 = 1;
      if (NoDiagonal()) 
  j0--; // Include first term if no diagonal
    
      for (i=0; i < NumMyRows_; i++) {
  int      NumEntries = *NumEntriesPerRow++;
  int *    RowIndices = *Indices++;
  double * RowValues  = *Vals++;
  if (!UnitDiagonal) 
    yp[i] = yp[i]/RowValues[0];
     double ytmp = yp[i];
  for (j=j0; j < NumEntries; j++) 
    yp[RowIndices[j]] -= RowValues[j] * ytmp;
      }
    }
    else {
      
      j0 = 1;
      if (NoDiagonal()) 
  j0--; // Include first term if no diagonal
    
      for (i=NumMyRows_-1; i >= 0; i--) {
  int      NumEntries = *NumEntriesPerRow-- - j0;
  int *    RowIndices = *Indices--;
  double * RowValues  = *Vals--;
  if (!UnitDiagonal) 
    yp[i] = yp[i]/RowValues[NumEntries];
     double ytmp = yp[i];
  for (j=0; j < NumEntries; j++) 
    yp[RowIndices[j]] -= RowValues[j] * ytmp;
      }
    }
    
  }
  UpdateFlops(2*NumGlobalNonzeros64());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Solve1(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

#ifdef EPETRA_CRSMATRIX_TEUCHOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Epetra_CrsMatrix::Solve(Upper,Trans,UnitDiag,X,Y)");
#endif

  //
  // This function find Y such that LY = X or UY = X or the transpose cases.
  //
  if((X.NumVectors() == 1) && (Y.NumVectors() == 1)) {
    double* xp = (double*) X[0];
    double* yp = (double*) Y[0];
    Epetra_Vector x(View, X.Map(), xp);
    Epetra_Vector y(View, Y.Map(), yp);
    EPETRA_CHK_ERR(Solve1(Upper, Trans, UnitDiagonal, x, y));
    return(0);
  }
  if(!Filled()) 
    EPETRA_CHK_ERR(-1); // Matrix must be filled.

  if((Upper) && (!UpperTriangular()))
    EPETRA_CHK_ERR(-2);
  if((!Upper) && (!LowerTriangular()))
    EPETRA_CHK_ERR(-3);
  if((!UnitDiagonal) && (NoDiagonal()))
    EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_))
    EPETRA_CHK_ERR(-5); // Need each row to have a diagonal

  int i, j, j0, k;
  int* NumEntriesPerRow = Graph().NumIndicesPerRow();
  int** Indices = Graph().Indices();
  double** Vals = Values();
  double diag = 0.0;

  // If upper, point to last row
  if((Upper && !Trans) || (!Upper && Trans)) {
    NumEntriesPerRow += NumMyRows_-1;
    Indices += NumMyRows_-1;
    Vals += NumMyRows_-1;
  }

  double** Xp = (double**) X.Pointers();
  double** Yp = (double**) Y.Pointers();

  int NumVectors = X.NumVectors();

  if(!Trans) {   
    if(Upper) {   
      j0 = 1;
      if(NoDiagonal()) 
  j0--; // Include first term if no diagonal
      for(i = NumMyRows_ - 1; i >= 0; i--) {
  int     NumEntries = *NumEntriesPerRow--;
  int*    RowIndices = *Indices--;
  double* RowValues  = *Vals--;
  if(!UnitDiagonal) 
    diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
  for(k = 0; k < NumVectors; k++) {
    double sum = 0.0;
    for(j = j0; j < NumEntries; j++) 
      sum += RowValues[j] * Yp[k][RowIndices[j]];
          
    if(UnitDiagonal) 
      Yp[k][i] = Xp[k][i] - sum;
    else 
      Yp[k][i] = (Xp[k][i] - sum) * diag;
  }
      }
    }
    else {
      j0 = 1;
      if(NoDiagonal()) 
  j0--; // Include first term if no diagonal
      for(i = 0; i < NumMyRows_; i++) {
  int     NumEntries = *NumEntriesPerRow++ - j0;
  int*    RowIndices = *Indices++;
  double* RowValues  = *Vals++;
  if(!UnitDiagonal)
    diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
  for(k = 0; k < NumVectors; k++) {
    double sum = 0.0;
    for(j = 0; j < NumEntries; j++) 
      sum += RowValues[j] * Yp[k][RowIndices[j]];
          
    if(UnitDiagonal) 
      Yp[k][i] = Xp[k][i] - sum;
    else 
      Yp[k][i] = (Xp[k][i] - sum)*diag;
  }
      }
    }
  }
  // ***********  Transpose case *******************************

  else {
    for(k = 0; k < NumVectors; k++) 
      if(Yp[k] != Xp[k]) 
  for(i = 0; i < NumMyRows_; i++)
    Yp[k][i] = Xp[k][i]; // Initialize y for transpose multiply
    
    if(Upper) {
      j0 = 1;
      if(NoDiagonal()) 
  j0--; // Include first term if no diagonal
      
      for(i = 0; i < NumMyRows_; i++) {
  int     NumEntries = *NumEntriesPerRow++;
  int*    RowIndices = *Indices++;
  double* RowValues  = *Vals++;
  if(!UnitDiagonal) 
    diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
  for(k = 0; k < NumVectors; k++) {
    if(!UnitDiagonal) 
      Yp[k][i] = Yp[k][i]*diag;
       double ytmp = Yp[k][i];
    for(j = j0; j < NumEntries; j++) 
      Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
  }
      }
    }
    else {
      j0 = 1;
      if(NoDiagonal()) 
  j0--; // Include first term if no diagonal  
      for(i = NumMyRows_ - 1; i >= 0; i--) {
  int     NumEntries = *NumEntriesPerRow-- - j0;
  int*    RowIndices = *Indices--;
  double* RowValues  = *Vals--;
     if (!UnitDiagonal)
       diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
  for(k = 0; k < NumVectors; k++) {
    if(!UnitDiagonal)  
      Yp[k][i] = Yp[k][i]*diag;
       double ytmp = Yp[k][i];
    for(j = 0; j < NumEntries; j++)
      Yp[k][RowIndices[j]] -= RowValues[j] * ytmp;
        }
      }
    }
  }
  
  UpdateFlops(2 * NumVectors * NumGlobalNonzeros64());
  return(0);
}



//=============================================================================    
Epetra_IntSerialDenseVector& Epetra_CrsMatrix::ExpertExtractIndexOffset(){ 
   return Graph_.CrsGraphData_->IndexOffset_;
 }
 
//=============================================================================    
Epetra_IntSerialDenseVector& Epetra_CrsMatrix::ExpertExtractIndices() {
  return Graph_.CrsGraphData_->data->All_Indices_;
 }

//=============================================================================    
int Epetra_CrsMatrix::ExpertMakeUniqueCrsGraphData(){
   if(Graph_.CrsGraphData_->ReferenceCount() > 1){
     Graph_.CrsGraphData_->DecrementReferenceCount();
     Graph_.CrsGraphData_=new Epetra_CrsGraphData(Copy, RowMap(),ColMap(),true); // true is for StaticProfile
   }
   return (0);
 }

//=============================================================================    
int Epetra_CrsMatrix::ExpertStaticFillComplete(const Epetra_Map & DomainMap,const Epetra_Map & RangeMap, const Epetra_Import * Importer,int NumMyDiagonals){

  Epetra_CrsGraphData& D=*Graph_.CrsGraphData_;
  int m=D.RowMap_.NumMyElements();

  bool UseLL=false;
  if(D.RowMap_.GlobalIndicesLongLong()) UseLL=true;
    
  // This only works for constant element size matrices w/ size=1
  if(!D.RowMap_.ConstantElementSize() || D.RowMap_.ElementSize()!=1 ||
     !D.ColMap_.ConstantElementSize() || D.ColMap_.ElementSize()!=1)
    EPETRA_CHK_ERR(-1);

  // Maps, import export
  D.DomainMap_ = DomainMap;
  D.RangeMap_  = RangeMap;

  if (!D.ColMap_.SameAs(D.DomainMap_)) {
    if (D.Importer_ != 0) {
      delete D.Importer_;
      D.Importer_ = 0;
    }
    if(Importer && Importer->SourceMap().SameAs(D.DomainMap_) && Importer->TargetMap().SameAs(D.ColMap_)){
      D.Importer_=Importer;
    }
    else {
      delete Importer;
      D.Importer_ = new Epetra_Import(D.ColMap_, D.DomainMap_,0,0);
    }
  }
  
  if (!D.RowMap_.SameAs(D.RangeMap_)) {
    if (D.Exporter_ != 0) {
      delete D.Exporter_;
      D.Exporter_ = 0;
    }
    D.Exporter_ = new Epetra_Export(D.RowMap_, D.RangeMap_); // Create Export object. 
  }

  // Matrix constants
  Allocated_                  = true;
  StaticGraph_                = true;
  UseTranspose_               = false;
  constructedWithFilledGraph_ = true;
  matrixFillCompleteCalled_   = true;
  StorageOptimized_           = true;
  squareFillCompleteCalled_   = false; // Always false for the full FillComplete
  // Note: Entries in D.Indices_ array point to memory we don't need to dealloc, but the D.Indices_
  // array itself needs deallocation.

  // Cleanup existing data
  for(int i=0;i<m;i++){
    if(Values_)            delete [] Values_[i];
    D.data->SortedEntries_[i].entries_.resize(0);
  }
  delete [] Values_;                 Values_=0;
  delete [] Values_alloc_lengths_;   Values_alloc_lengths_=0;
  delete [] D.data->Indices_;        D.data->Indices_=0; 
  D.NumAllocatedIndicesPerRow_.Resize(0);


  // GraphData constants
  D.Filled_                                = true;
  D.Allocated_                             = true;
  D.Sorted_                                = false;
  D.StorageOptimized_                      = true;
  D.NoRedundancies_                        = true;
  D.IndicesAreGlobal_                      = false;
  D.IndicesAreLocal_                       = true;
  D.IndicesAreContiguous_                  = true;
  D.GlobalConstantsComputed_               = true;
  D.StaticProfile_                         = true;
  D.SortGhostsAssociatedWithEachProcessor_ = true;
  D.HaveColMap_                            = true;

  // Don't change the index base

  // Local quantities (no computing)
  int nnz=D.IndexOffset_[m]-D.IndexOffset_[0];
  D.NumMyRows_       = D.NumMyBlockRows_ = m;
  D.NumMyCols_       = D.NumMyBlockCols_ = D.ColMap_.NumMyElements();
  D.NumMyNonzeros_   = D.NumMyEntries_   = nnz;
  D.MaxRowDim_       = 1;
  D.MaxColDim_       = 1;

  // A priori globals
  D.GlobalMaxRowDim_ = 1;
  D.GlobalMaxColDim_ = 1;

  // Compute max NZ per row, indices, diagonals, triangular
  D.MaxNumIndices_=0;
  for(int i=0; i<m; i++){
    int NumIndices=D.IndexOffset_[i+1]-D.IndexOffset_[i];
    D.MaxNumIndices_=EPETRA_MAX(D.MaxNumIndices_,NumIndices);

    // Code copied from Epetra_CrsGraph.Determine_Triangular()
    if(NumIndices > 0) {
      const int* col_indices = & D.data->All_Indices_[D.IndexOffset_[i]];

      int jl_0 = col_indices[0];
      int jl_n = col_indices[NumIndices-1];

      if(jl_n > i) D.LowerTriangular_ = false;
      if(jl_0 < i) D.UpperTriangular_ = false;


      if(NumMyDiagonals == -1) {
	// Binary search in GID space
	// NOTE: This turns out to be noticibly faster than doing a single LID lookup and then binary searching
	// on the LIDs directly.
	int insertPoint = -1;
	if(UseLL)  {
	  long long jg = D.RowMap_.GID64(i);
	  if (Epetra_Util_binary_search_aux(jg, col_indices, D.ColMap_.MyGlobalElements64(), NumIndices, insertPoint)>-1) {
	    D.NumMyBlockDiagonals_++;
	    D.NumMyDiagonals_ ++; 
	  }
	}
	else {
	  int jg = D.RowMap_.GID(i);
	  if (Epetra_Util_binary_search_aux(jg, col_indices, D.ColMap_.MyGlobalElements(), NumIndices, insertPoint)>-1) {
	    D.NumMyBlockDiagonals_++;
	    D.NumMyDiagonals_ ++; 
	  }
	}
      }
    }
  } // end DetermineTriangular

  if(NumMyDiagonals > -1) D.NumMyDiagonals_ = D.NumMyBlockDiagonals_ = NumMyDiagonals;

  D.MaxNumNonzeros_=D.MaxNumIndices_;

  // Compute global constants - Code copied from Epetra_CrsGraph.ComputeGlobalConstants()
  Epetra_IntSerialDenseVector tempvec(8); // Temp space
  tempvec[0] = D.NumMyEntries_;
  tempvec[1] = D.NumMyBlockDiagonals_;
  tempvec[2] = D.NumMyDiagonals_;
  tempvec[3] = D.NumMyNonzeros_;
  
  Comm().SumAll(&tempvec[0], &tempvec[4], 4);
    
  D.NumGlobalEntries_        = tempvec[4];
  D.NumGlobalBlockDiagonals_ = tempvec[5];
  D.NumGlobalDiagonals_      = tempvec[6];
  D.NumGlobalNonzeros_       = tempvec[7];

  tempvec[0] = D.MaxNumIndices_;
  tempvec[1] = D.MaxNumNonzeros_;

  Comm().MaxAll(&tempvec[0], &tempvec[2], 2);
  
  D.GlobalMaxNumIndices_  = tempvec[2];
  D.GlobalMaxNumNonzeros_ = tempvec[3];
  //end  Epetra_CrsGraph.ComputeGlobalConstants()


  // Global Info
  if(UseLL) {
    D.NumGlobalRows_ = D.RangeMap_.NumGlobalPoints64();
    D.NumGlobalCols_ = D.DomainMap_.NumGlobalPoints64();
  }
  else {
    D.NumGlobalRows_ = (int) D.RangeMap_.NumGlobalPoints64();
    D.NumGlobalCols_ = (int) D.DomainMap_.NumGlobalPoints64();
  }
  return 0;
}

//=========================================================================
int Epetra_CrsMatrix::PackAndPrepareWithOwningPIDs(const Epetra_SrcDistObject & Source, 
						   int NumExportIDs,
						   int * ExportLIDs,
						   int & LenExports,
						   char *& Exports,
						   int & SizeOfPacket,
						   int * Sizes,
						   bool & VarSizes,
						   Epetra_Distributor & Distor)
{
  (void)Distor;	

  // Rest of work can be done using CrsMatrix only  
  const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Source);

  // Grab the importer, if we have one
  const Epetra_Import *MyImporter= A.Importer();

  VarSizes = true; //enable variable block size data comm

  int TotalSendLength = 0;
  int * IntSizes = 0; 
  if( NumExportIDs>0 ) IntSizes = new int[NumExportIDs];

  int SizeofIntType = -1;
  if(Source.Map().GlobalIndicesInt())
    SizeofIntType = (int)sizeof(int); 
  else if(Source.Map().GlobalIndicesLongLong())
    SizeofIntType = (int)sizeof(long long); 
  else
    throw ReportError("Epetra_CrsMatrix::PackAndPrepareWithOwningPIDs: Unable to determine source global index type",-1);

  for( int i = 0; i < NumExportIDs; ++i )
  {    
    int NumEntries;
    A.NumMyRowEntries( ExportLIDs[i], NumEntries );
    // Will have NumEntries doubles, 2*NumEntries +2 ints pack them interleaved     Sizes[i] = NumEntries;
    // NTS: We add the owning PID as the SECOND int of the pair for each entry
    Sizes[i] = NumEntries;
    // NOTE: Mixing and matching Int Types would be more efficient, BUT what about variable alignment?
    IntSizes[i] = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
    TotalSendLength += (Sizes[i]+IntSizes[i]);
  }    
         
  double * DoubleExports = 0; 
  SizeOfPacket = (int)sizeof(double);
       
  //setup buffer locally for memory management by this object
  if( TotalSendLength*SizeOfPacket > LenExports )
  {
    if( LenExports > 0 ) delete [] Exports;
    LenExports = TotalSendLength*SizeOfPacket;
    DoubleExports = new double[TotalSendLength];
    for( int i = 0; i < TotalSendLength; ++i ) DoubleExports[i] = 0.0;
    Exports = (char *) DoubleExports;
  } 

  int NumEntries;
  double * values;
  double * valptr, * dintptr; 

  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumEntries, Number of indices in row
  // next 2*NumEntries: The actual indices and owning [1] PID each for the row in (GID,PID) pairs with the GID first.

  // [1] Owning is defined in the sense of "Who owns the GID in the DomainMap," aka, who sends the GID in the importer

  const Epetra_Map & rowMap = A.RowMatrixRowMap();
  const Epetra_Map & colMap = A.RowMatrixColMap();

  if( NumExportIDs > 0 )
  {
    // Grab the owning PIDs from the Importer (distributor, really)
    Epetra_Util util;
    std::vector< int > pids;
    if(MyImporter) util.GetPids(*MyImporter,pids,false);
    else {
      pids.resize(colMap.NumMyElements());
      pids.assign(colMap.NumMyElements(),A.Comm().MyPID());
    }

    if(Source.Map().GlobalIndicesInt()) { 
      int * Indices;
      int FromRow; 
      int * intptr;                         
        
      int maxNumEntries = A.MaxNumEntries();
      std::vector<int> MyIndices(maxNumEntries);
      dintptr = (double *) Exports;
      valptr = dintptr + IntSizes[0];
      intptr = (int *) dintptr;
      for (int i=0; i<NumExportIDs; i++)
	{
	  FromRow   = (int) rowMap.GID64(ExportLIDs[i]);
	  intptr[0] = FromRow;
	  values    = valptr;
	  Indices   = intptr + 2;
	  EPETRA_CHK_ERR(A.ExtractMyRowCopy(ExportLIDs[i], maxNumEntries, NumEntries, values, &MyIndices[0]));
	  for (int j=0; j<NumEntries; j++) {
	    Indices[2*j]   = (int)colMap.GID64(MyIndices[j]);   // convert to GIDs
	    Indices[2*j+1] = pids[MyIndices[j]];               // PID owning the entry.
	  }
	  intptr[1] = NumEntries; // Load second slot of segment
	  if( i < (NumExportIDs-1) )
	    {
	      dintptr += (IntSizes[i]+Sizes[i]);
	      valptr = dintptr + IntSizes[i+1];
	      intptr = (int *) dintptr;
	    }	
	}    
    }
    else if(Source.Map().GlobalIndicesLongLong()) {
      long long * LL_Indices;
      long long FromRow; 
      long long * LLptr;                         
  
      int maxNumEntries = A.MaxNumEntries();
      std::vector<int> MyIndices(maxNumEntries);
      
      dintptr = (double *) Exports;
      valptr = dintptr + IntSizes[0];
      LLptr = (long long *) dintptr;
      for (int i=0; i<NumExportIDs; i++)
	{
	  FromRow = rowMap.GID64(ExportLIDs[i]);
	  LLptr[0]   = FromRow;
	  values     = valptr;
	  LL_Indices = LLptr + 2;
	  EPETRA_CHK_ERR(A.ExtractMyRowCopy(ExportLIDs[i], maxNumEntries, NumEntries, values, &MyIndices[0]));
	  for (int j=0; j<NumEntries; j++) {
	    LL_Indices[2*j]   = colMap.GID64(MyIndices[j]);   // convert to GIDs
	    LL_Indices[2*j+1] = pids[MyIndices[j]];           // PID owning the entry.

	  }
	  LLptr[1] = NumEntries; // Load second slot of segment
	  if( i < (NumExportIDs-1) )
	    {
	      dintptr += (IntSizes[i]+Sizes[i]);
	      valptr = dintptr + IntSizes[i+1];
	      LLptr = (long long *) dintptr;
	    }	
	}    
    }
    
    for( int i = 0; i < NumExportIDs; ++i )
      Sizes[i] += IntSizes[i];
  }

  if( IntSizes ) delete [] IntSizes;

  return(0);
}


// private ===================================================================
 template<typename int_type>
 int Epetra_CrsMatrix::LowCommunicationMakeColMapAndReindex(const Epetra_Map& domainMap, const int * owningPIDs, std::vector<int>& RemotePIDs, const int_type *colind_LL)
   {
  int i,j;

  // Sanity checking - If int_type==int, colind (aka Graph_.CrsGraphData_->data->All_Indices_)  contains GIDs, if not, colind_LL does.
  bool UseLL=false;
  if(RowMap().GlobalIndicesLongLong()) UseLL=true;

  if(!UseLL && colind_LL)
    throw ReportError("LowCommunicationMakeColMapAndReindex can't determine whether it is in long-long or int mode due to argument conflict",-1);
  
  // Note: If we legitimately have no column indices, then colind_LL can be null even if we're in UseLL mode.

  // Scan all column indices and sort into two groups: 
  // Local:  those whose GID matches a GID of the domain map on this processor and
  // Remote: All others.
  int numDomainElements = domainMap.NumMyElements();
  bool * LocalGIDs  = 0;
  if (numDomainElements>0) LocalGIDs  = new bool[numDomainElements];
  for (i=0; i<numDomainElements; i++) LocalGIDs[i] = false; // Assume domain GIDs are not local

  bool DoSizes = !domainMap.ConstantElementSize(); // If not constant element size, then error
  if(DoSizes) EPETRA_CHK_ERR(-1);

  // Because these will always be useful
  int * rowptr = Graph_.CrsGraphData_->IndexOffset_.Values();
  int * colind = Graph_.CrsGraphData_->data->All_Indices_.Values(); 

  // Meta-pointer to get around GID/LID issue
  const int_type * colind_GID=0;
  if(UseLL) colind_GID = (int_type*) colind_LL;
  else colind_GID = (int_type*) colind;

  // In principle it is good to have RemoteGIDs and RemotGIDList be as long as the number of remote GIDs
  // on this processor, but this would require two passes through the column IDs, so we make it the max of 100
  // and the number of block rows.
  const int numMyBlockRows = NumMyRows();
  int  hashsize = numMyBlockRows; if (hashsize < 100) hashsize = 100;
  Epetra_HashTable<int_type> RemoteGIDs(hashsize); 
  std::vector<int_type> RemoteGIDList; RemoteGIDList.reserve(hashsize);
  std::vector<int> PIDList;            PIDList.reserve(hashsize);

  // Here we start using the *int* colind array.  If int_type==int this clobbers the GIDs, if
  // int_type==long long, then this is the first use of the colind array.
  // For *local* GID's set colind with with their LID in the domainMap.  For *remote* GIDs, 
  // we set colind with (numDomainElements+NumRemoteColGIDs) before the increment of
  // the remote count.  These numberings will be separate because no local LID is greater 
  // than numDomainElements. 

  int NumLocalColGIDs = 0;
  int NumRemoteColGIDs = 0;
  for(i = 0; i < numMyBlockRows; i++) {
    for(j = rowptr[i]; j < rowptr[i+1]; j++) {
      int_type GID = colind_GID[j];
      // Check if GID matches a row GID
      int LID = domainMap.LID(GID);
      if(LID != -1) {
	bool alreadyFound = LocalGIDs[LID];
	if (!alreadyFound) {
          LocalGIDs[LID] = true; // There is a column in the graph associated with this domain map GID
          NumLocalColGIDs++;
	}
	colind[j] = LID; 
      }
      else {
	int_type hash_value=RemoteGIDs.Get(GID);
	if(hash_value  == -1) { // This means its a new remote GID
	  int PID = owningPIDs[j];
	  if(PID==-1) throw ReportError("LowCommunicationMakeColMapAndReindex: Cannot figure out if PID is owned.",-1);
	  colind[j] = numDomainElements + NumRemoteColGIDs;
	  RemoteGIDs.Add(GID, NumRemoteColGIDs);
	  RemoteGIDList.push_back(GID);
	  PIDList.push_back(PID);
	  NumRemoteColGIDs++;
	}
	else
	  colind[j] = numDomainElements + hash_value;	  
      }
    }
  }

  // Possible short-circuit:  If all domain map GIDs are present as column indices, then set ColMap=domainMap and quit
  if (domainMap.Comm().NumProc()==1) { 
    
    if (NumRemoteColGIDs!=0) {
      throw ReportError("Some column IDs are not in domainMap.  If matrix is rectangular, you must pass in domainMap to FillComplete",-2); 
      // Sanity test: When one processor,there can be no remoteGIDs
    }
    if (NumLocalColGIDs==numDomainElements) {
      Graph_.CrsGraphData_->ColMap_ = domainMap;
      Graph_.CrsGraphData_->HaveColMap_ = true;
      if (LocalGIDs!=0) delete [] LocalGIDs; 
      // In this case, we just use the domainMap's indices, which is, not coincidently, what we clobbered colind with up above anyway. 
      // No further reindexing is needed.
      return(0); 
    }
  }
      
  // Now build integer array containing column GIDs
  // Build back end, containing remote GIDs, first
  int numMyBlockCols = NumLocalColGIDs + NumRemoteColGIDs;
  std::vector<int_type> ColIndices;
  int_type * RemoteColIndices=0;
  if(numMyBlockCols > 0) {
    ColIndices.resize(numMyBlockCols);
    RemoteColIndices = &ColIndices[NumLocalColGIDs]; // Points to back end of ColIndices
  }

  for(i = 0; i < NumRemoteColGIDs; i++) 
    RemoteColIndices[i] = RemoteGIDList[i]; 

  // Build permute array for *remote* reindexing.
  std::vector<int> RemotePermuteIDs(NumRemoteColGIDs);
  for(i=0; i<NumRemoteColGIDs; i++) RemotePermuteIDs[i]=i;

  // Sort External column indices so that all columns coming from a given remote processor are contiguous
  int NumListsInt=0;
  int NumListsLL =0;
  int* IntSortLists[2];
  long long * LLSortLists[2];
  if(!UseLL) {
    // int version
    IntSortLists[0] = (int*) RemoteColIndices;
    IntSortLists[1] = RemotePermuteIDs.data();
    NumListsInt=2;
  }
  else {
    //LL version
    LLSortLists[0]  = (long long*) RemoteColIndices;
    IntSortLists[0] = RemotePermuteIDs.data();
    NumListsInt = NumListsLL = 1;
  }
 
  Epetra_Util::Sort(true, NumRemoteColGIDs, PIDList.data(), 0, 0, NumListsInt, IntSortLists,NumListsLL,LLSortLists);

  // Stash the RemotePIDs  
  PIDList.resize(NumRemoteColGIDs);
  RemotePIDs = PIDList;

  if (Graph_.CrsGraphData_->SortGhostsAssociatedWithEachProcessor_) {
    // Sort external column indices so that columns from a given remote processor are not only contiguous
    // but also in ascending order. NOTE: I don't know if the number of externals associated
    // with a given remote processor is known at this point ... so I count them here.

    // NTS: Only sort the RemoteColIndices this time...
    int StartCurrent, StartNext;
    StartCurrent = 0; StartNext = 1;
    while ( StartNext < NumRemoteColGIDs ) {
      if (PIDList[StartNext]==PIDList[StartNext-1]) StartNext++;
      else {
	IntSortLists[0] =  &RemotePermuteIDs[StartCurrent];
	Epetra_Util::Sort(true,StartNext-StartCurrent, &(RemoteColIndices[StartCurrent]),0,0,1,IntSortLists,0,0);
        StartCurrent = StartNext; StartNext++;
      }
    }
    IntSortLists[0] =  &RemotePermuteIDs[StartCurrent];
    Epetra_Util::Sort(true, StartNext-StartCurrent, &(RemoteColIndices[StartCurrent]), 0, 0, 1,IntSortLists,0,0);
  }

  // Reverse the permutation to get the information we actually care about
  std::vector<int> ReverseRemotePermuteIDs(NumRemoteColGIDs);
  for(i=0; i<NumRemoteColGIDs; i++) ReverseRemotePermuteIDs[RemotePermuteIDs[i]]=i;

  // Build permute array for *local* reindexing.
  bool use_local_permute=false;
  std::vector<int> LocalPermuteIDs(numDomainElements);

  // Now fill front end. Two cases:
  // (1) If the number of Local column GIDs is the same as the number of Local domain GIDs, we
  //     can simply read the domain GIDs into the front part of ColIndices, otherwise 
  // (2) We step through the GIDs of the domainMap, checking to see if each domain GID is a column GID.
  //     we want to do this to maintain a consistent ordering of GIDs between the columns and the domain.

  if(NumLocalColGIDs == domainMap.NumMyElements()) {
    if(NumLocalColGIDs > 0) {
      domainMap.MyGlobalElements(&ColIndices[0]); // Load Global Indices into first numMyBlockCols elements column GID list
    }
  }
  else {
    int_type* MyGlobalElements;
    if (!UseLL) MyGlobalElements = (int_type*)domainMap.MyGlobalElements();
    else MyGlobalElements = (int_type*)domainMap.MyGlobalElements64();

    int* ElementSizeList = 0;
    if(DoSizes) 
      ElementSizeList = domainMap.ElementSizeList();
    int NumLocalAgain = 0;
    use_local_permute = true;    
    for(i = 0; i < numDomainElements; i++) {
      if(LocalGIDs[i]) {
	LocalPermuteIDs[i] = NumLocalAgain;
	ColIndices[NumLocalAgain++] = MyGlobalElements[i];
      }
    }
    assert(NumLocalAgain==NumLocalColGIDs); // Sanity test
  }

  // Done with this array
  if (LocalGIDs!=0) delete [] LocalGIDs; 

  // Make Column map with same element sizes as Domain map 
  Epetra_Map temp((int_type)(-1), numMyBlockCols, ColIndices.data(), (int)domainMap.IndexBase64(), domainMap.Comm());

  Graph_.CrsGraphData_->ColMap_ = temp;
  Graph_.CrsGraphData_->HaveColMap_ = true;

  // Low-cost reindex of the matrix
  for(i=0; i<numMyBlockRows; i++){
    for(j=rowptr[i]; j<rowptr[i+1]; j++){
      int ID=colind[j];
      if(ID < numDomainElements){
	if(use_local_permute) colind[j] = LocalPermuteIDs[colind[j]];
	// In the case where use_local_permute==false, we just copy the DomainMap's ordering, which it so happens
	// is what we put in colind to begin with.
      }
      else
	colind[j] =  NumLocalColGIDs + ReverseRemotePermuteIDs[colind[j]-numDomainElements];
    }
  }
  
  return(0);
}



// ===================================================================
template<typename int_type>
int Epetra_CrsMatrix::UnpackAndCombineIntoCrsArrays(const Epetra_SrcDistObject& Source, 
						    int NumSameIDs,
						    int NumRemoteIDs,
						    const int * RemoteLIDs,
						    int NumPermuteIDs,
						    const int *PermuteToLIDs,
						    const int *PermuteFromLIDs,
						    int LenImports,
						    char* Imports,
						    std::vector<int> & pids,
						    std::vector<int_type> &CSR_colind_LL)
{
  // What we really need to know is where in the CSR arrays each row should start (aka the rowptr).
  // We do that by (a) having each row record it's size in the rowptr (b) doing a cumulative sum to get the rowptr values correct and 
  // (c) Having each row copied into the right colind / values locations.

  // From Epetra_CrsMatrix UnpackAndCombine():
  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumEntries, Number of indices in row.
  // next NumEntries: The actual indices for the row.

  int i,j,rv;
  int N = NumMyRows();
  int MyPID = Comm().MyPID();

  // Rest of work can be done using CrsMatrix only  
  const Epetra_CrsMatrix & SourceMatrix = dynamic_cast<const Epetra_CrsMatrix &>(Source);

  // Pull and pointers & allocate memory (rowptr should be the right size already)
  Epetra_IntSerialDenseVector & CSR_rowptr = ExpertExtractIndexOffset();
  Epetra_IntSerialDenseVector & CSR_colind = ExpertExtractIndices();  
  double *&                     CSR_vals   = ExpertExtractValues();  
  // Presume the RowPtr is the right size

  bool UseLL=false;
  int SizeofIntType = -1;
  if(SourceMatrix.RowMap().GlobalIndicesInt()) {
    SizeofIntType = (int)sizeof(int); 
    UseLL=false;
  }
  else if(SourceMatrix.RowMap().GlobalIndicesLongLong()) {
    SizeofIntType = (int)sizeof(long long); 
    UseLL=true;
  }
  else
    throw ReportError("EpetraExt::npackAndCombineIntoStaticProfile: Unable to determine source global index type",-1);


  // SameIDs: Always first, always in the same place
  for(i=0; i<NumSameIDs; i++)
    CSR_rowptr[i]=SourceMatrix.NumMyEntries(i);
  
  // PermuteIDs: Still local, but reordered
  for(i=0; i<NumPermuteIDs; i++)
    CSR_rowptr[PermuteToLIDs[i]] = SourceMatrix.NumMyEntries(PermuteFromLIDs[i]);
  
  // RemoteIDs:  RemoteLIDs tells us the ID, we need to look up the length the hard way.  See UnpackAndCombine for where this code came from
  if(NumRemoteIDs > 0) {
    double * dintptr = (double *) Imports_;
    if(!UseLL) { 
      // Int version
      int    *  intptr = (int *) dintptr;
      int   NumEntries = intptr[1];
      int      IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
      for(i=0; i<NumRemoteIDs; i++) {
	CSR_rowptr[RemoteLIDs[i]] = NumEntries;
	
	if( i < (NumRemoteIDs-1) ) {
	  dintptr += IntSize + NumEntries;
	  intptr = (int *) dintptr;
	  NumEntries = intptr[1];
	  IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
	}
      }
    }
    else {
      // LongLong version
      long long     *  LLptr = (long long *) dintptr;
      int         NumEntries = (int) LLptr[1];
      int      IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
      for(i=0; i<NumRemoteIDs; i++) {
	CSR_rowptr[RemoteLIDs[i]] = NumEntries;
	
	if( i < (NumRemoteIDs-1) ) {
	  dintptr += IntSize + NumEntries;
	  LLptr = (long long *) dintptr;
	  NumEntries = (int) LLptr[1];
	  IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
	}
      }
    }
  }
  
  // Turn row length into a real CSR_rowptr
  int last_len = CSR_rowptr[0];
  CSR_rowptr[0] = 0;
  for(i=1; i<N+1; i++){
    int new_len = CSR_rowptr[i];
    CSR_rowptr[i] = last_len + CSR_rowptr[i-1];
    last_len=new_len;
  }

  // Allocate CSR_colind & CSR_values arrays
  int mynnz=CSR_rowptr[N];
  CSR_colind.Resize(mynnz);
  if(UseLL) CSR_colind_LL.resize(mynnz);
  delete [] CSR_vals; // should be a noop.
  CSR_vals = new double[mynnz];

  // Get the list of owning PIDs from SourceMatrix's importer
  std::vector<int> SourcePids(SourceMatrix.NumMyCols(),-1);
  if(SourceMatrix.Importer()) Epetra_Util::GetPids(*SourceMatrix.Importer(),SourcePids,true);

  // Get a list of GIDs ready for colmap, preseed with -1 for local
  pids.resize(mynnz); pids.assign(mynnz,-1);
 
  // Grab pointers for SourceMatrix
  int    * Source_rowptr, * Source_colind;
  double * Source_vals;
  rv=SourceMatrix.ExtractCrsDataPointers(Source_rowptr,Source_colind,Source_vals);
  if(rv) throw ReportError("Epetra_CrsMatrix: Fused copy constructor failed in ExtractCrsDataPointers",-4);

  // SameIDs: Copy the data over
  for(i=0; i<NumSameIDs; i++) {
    int FromRow = Source_rowptr[i];
    int ToRow   = CSR_rowptr[i];

    for(j=Source_rowptr[i]; j<Source_rowptr[i+1]; j++) {
      CSR_vals[ToRow + j - FromRow]                = Source_vals[j];      
      if(UseLL) CSR_colind_LL[ToRow + j - FromRow] = SourceMatrix.GCID64(Source_colind[j]);
      else      CSR_colind[ToRow + j - FromRow]    = SourceMatrix.GCID(Source_colind[j]);
      pids[ToRow + j - FromRow]                    = (SourcePids[Source_colind[j]] != MyPID) ? SourcePids[Source_colind[j]] : -1;
    }
  }

  // PermuteIDs: Copy the data over
  for(i=0; i<NumPermuteIDs; i++) {
    int FromLID  = PermuteFromLIDs[i];
    int FromRow = Source_rowptr[FromLID];
    int ToRow   = CSR_rowptr[PermuteToLIDs[i]];

    for(j=Source_rowptr[FromLID]; j<Source_rowptr[FromLID+1]; j++) {
      CSR_vals[ToRow + j - FromRow]                = Source_vals[j];      
      if(UseLL) CSR_colind_LL[ToRow + j - FromRow] = SourceMatrix.GCID64(Source_colind[j]);
      else      CSR_colind[ToRow + j - FromRow]    = SourceMatrix.GCID(Source_colind[j]);
      pids[ToRow + j - FromRow]                    = (SourcePids[Source_colind[j]] != MyPID) ? SourcePids[Source_colind[j]] : -1;
    }
  }

  // RemoteIDs: Loop structure following UnpackAndCombine  
  if(NumRemoteIDs > 0) {
    double * dintptr = (double *) Imports_;
    if(!UseLL) {
      // int version
      int *    intptr  = (int *) dintptr;
      int NumEntries   = intptr[1];
      int IntSize      = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
      double* valptr   = dintptr + IntSize;
      
      for (i=0; i<NumRemoteIDs; i++) {
	int ToLID    = RemoteLIDs[i];
	int StartRow = CSR_rowptr[ToLID];
	
	double * values  = valptr;
	int    * Indices = intptr + 2;
	for(j=0; j<NumEntries; j++){
	  CSR_colind[StartRow + j]  = Indices[2*j];
	  if(MyPID !=  Indices[2*j+1]) pids[StartRow+j]   = Indices[2*j+1];
	  CSR_vals[StartRow + j]   = values[j];	  
	}
	if( i < (NumRemoteIDs-1) ) {
	  dintptr += IntSize + NumEntries;
	  intptr = (int *) dintptr;
	  NumEntries = intptr[1];
	  IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
	  valptr = dintptr + IntSize;
	}
      }
    }
    else {
      // Long Long Version
      long long * LLptr  = (long long *) dintptr;
      int NumEntries     = (int) LLptr[1];
      int IntSize        = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
      double* valptr     = dintptr + IntSize;
      
      for (i=0; i<NumRemoteIDs; i++) {
	int ToLID    = RemoteLIDs[i];
	int StartRow = CSR_rowptr[ToLID];
	
	double * values        = valptr;
	long long    * Indices = LLptr + 2;
	for(j=0; j<NumEntries; j++){
	  CSR_colind_LL[StartRow + j]  = Indices[2*j];
	  if(MyPID !=  Indices[2*j+1]) pids[StartRow+j]   = (int) Indices[2*j+1];
	  CSR_vals[StartRow + j]   = values[j];	  
	}
	if( i < (NumRemoteIDs-1) ) {
	  dintptr += IntSize + NumEntries;
	  LLptr    = (long long *) dintptr;
	  NumEntries = (int) LLptr[1];
	  IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
	  valptr = dintptr + IntSize;
	}
      }    
    }
  }
  return 0;
}


// ===================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(const Epetra_CrsMatrix & SourceMatrix, const Epetra_Import & RowImporter,const Epetra_Map * RangeMap)
   : Epetra_DistObject(RowImporter.TargetMap(), "Epetra::CrsMatrix"),
   Epetra_CompObject(),
   Epetra_BLAS(),
   Graph_(Copy, RowImporter.TargetMap(), 0, false),
   Values_(0),
   Values_alloc_lengths_(0),
   All_Values_(0),
   NumMyRows_(RowImporter.TargetMap().NumMyPoints()),
   ImportVector_(0),
   ExportVector_(0),
   CV_(Copy)
{
   // Fused constructor, import & FillComplete
  int rv;
  int N = NumMyRows();
  Epetra_Util util;

  bool communication_needed = RowImporter.SourceMap().DistributedGlobal();

  bool UseLL=false;
  int SizeofIntType = -1;
  if(SourceMatrix.RowMap().GlobalIndicesInt()) {
    SizeofIntType = (int)sizeof(int); 
    UseLL=false;
  }
  else if(SourceMatrix.RowMap().GlobalIndicesLongLong()) {
    SizeofIntType = (int)sizeof(long long); 
    UseLL=true;
  }
  else
    throw ReportError("EpetraExt::LightweightCrsMatrix::PackAndPrepare: Unable to determine source global index type",-1);


  // The basic algorithm here is:
  // 1) Call Distor.Do to handle the import.
  // 2) Copy all the Imported and Copy/Permuted data into the raw Epetra_CrsMatrix / Epetra_CrsGraphData pointers, still using GIDs.
  // 3) Call an optimized version of MakeColMap that avoids the Directory lookups (since the importer knows who owns all the gids) AND
  //    reindexes to LIDs.
  // 4) Call ExpertStaticFillComplete(B.DomainMap(),A.RangeMap());

  // WARNING: This assumes both GIDs and LIDs are ints (we use the All_Indices array for both).  Will need to fix for Epetra64.

  // Sanity Check
  if (!SourceMatrix.RowMap().SameAs(RowImporter.SourceMap())) 
    throw ReportError("Epetra_CrsMatrix: Fused copy constructor requires Importer.SourceMap() to match SourceMatrix.RowMap()",-1);
  
  // Get information from the Importer
  int NumSameIDs             = RowImporter.NumSameIDs();
  int NumPermuteIDs          = RowImporter.NumPermuteIDs();
  int NumRemoteIDs           = RowImporter.NumRemoteIDs();
  int NumExportIDs           = RowImporter.NumExportIDs();
  int* ExportLIDs            = RowImporter.ExportLIDs();
  int* RemoteLIDs            = RowImporter.RemoteLIDs();
  int* PermuteToLIDs         = RowImporter.PermuteToLIDs();
  int* PermuteFromLIDs       = RowImporter.PermuteFromLIDs();
  Epetra_Distributor& Distor = RowImporter.Distributor();

  // Pull and pointers & allocate memory (rowptr should be the right size already)
  Epetra_IntSerialDenseVector & CSR_rowptr = ExpertExtractIndexOffset();
  Epetra_IntSerialDenseVector & CSR_colind = ExpertExtractIndices();  
  double *&                     CSR_vals   = ExpertExtractValues();  
  CSR_rowptr.Resize(N+1);

  // Unused if we're not in LL mode
  std::vector<long long> CSR_colind_LL;

  /***************************************************/
  /***** 1) From Epetra_DistObject::DoTransfer() *****/
  /***************************************************/
  rv=CheckSizes(SourceMatrix);
  if(rv) throw ReportError("Epetra_CrsMatrix: Fused copy constructor failed in CheckSizes()",-2);

  int SizeOfPacket; 
  bool VarSizes = false;
  if( NumExportIDs > 0) {
    delete [] Sizes_;
    Sizes_ = new int[NumExportIDs];
  }

  rv=PackAndPrepareWithOwningPIDs(SourceMatrix, NumExportIDs, ExportLIDs,LenExports_, Exports_, SizeOfPacket, Sizes_, VarSizes, Distor);
  if(rv) throw ReportError("Epetra_CrsMatrix: Fused copy constructor failed in PackAndPrepareWithOwningPIDs()",-3);

  if (communication_needed) {
    // Do the exchange of remote data
    if( VarSizes ) 
      rv=Distor.Do(Exports_, SizeOfPacket, Sizes_, LenImports_, Imports_);
    else
      rv=Distor.Do(Exports_, SizeOfPacket, LenImports_, Imports_);
  }
  if(rv) throw ReportError("Epetra_CrsMatrix: Fused copy constructor failed in Distor.Do",-3);


  /*********************************************************************/
  /**** 2) Copy all of the Same/Permute/Remote data into CSR_arrays ****/
  /*********************************************************************/

  std::vector<int> pids;
  if(UseLL) 
    UnpackAndCombineIntoCrsArrays<long long>(SourceMatrix,NumSameIDs,NumRemoteIDs,RemoteLIDs,NumPermuteIDs,PermuteToLIDs,PermuteFromLIDs,LenImports_,Imports_,pids,CSR_colind_LL);
  else {
    std::vector<int> dummy;
    UnpackAndCombineIntoCrsArrays<int>(SourceMatrix,NumSameIDs,NumRemoteIDs,RemoteLIDs,NumPermuteIDs,PermuteToLIDs,PermuteFromLIDs,LenImports_,Imports_,pids,dummy);
  }

  /**************************************************************/
  /**** 3) Call Optimized MakeColMap w/ no Directory Lookups ****/
  /**************************************************************/
  //Call an optimized version of MakeColMap that avoids the Directory lookups (since the importer knows who owns all the gids).
  std::vector<int> RemotePIDs;

  if(UseLL) {
    LowCommunicationMakeColMapAndReindex<long long>(SourceMatrix.DomainMap(),pids.data(),RemotePIDs,CSR_colind_LL.data());
  }
  else LowCommunicationMakeColMapAndReindex<int>(SourceMatrix.DomainMap(),pids.data(),RemotePIDs);  

  /********************************************/
  /**** 5) Call ExpertStaticFillComplete() ****/
  /********************************************/
  // Sort the entries
  Epetra_Util::SortCrsEntries(N, CSR_rowptr.Values(), CSR_colind.Values(), CSR_vals);

  // Pre-build the importer using the existing PIDs
  Epetra_Import * MyImport=0;
  int NumRemotePIDs = RemotePIDs.size();
  if(!SourceMatrix.DomainMap().SameAs(ColMap()))
    MyImport = new Epetra_Import(ColMap(),SourceMatrix.DomainMap(),NumRemotePIDs,RemotePIDs.data());

  if(RangeMap) ExpertStaticFillComplete(SourceMatrix.DomainMap(),*RangeMap,MyImport); 
  else {
    const Epetra_Map *MyRowMap = dynamic_cast<const Epetra_Map*>(&RowImporter.TargetMap());
    ExpertStaticFillComplete(SourceMatrix.DomainMap(),*MyRowMap,MyImport);
  }
  // Note: ExpertStaticFillComplete assumes ownership of the importer, if we made one...
  // We are not to deallocate that here.

}// end fused copy constructor

 


