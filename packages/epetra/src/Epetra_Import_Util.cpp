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
#include "Epetra_Import_Util.h"
#include "Epetra_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_HashTable.h"

#include <stdexcept>


namespace Epetra_Import_Util { 

// =========================================================================
// =========================================================================
template<typename int_type>
int TPackAndPrepareWithOwningPIDs(const Epetra_CrsMatrix & A, 
				 int NumExportIDs,
				 int * ExportLIDs,
				 int & LenExports,
				 char *& Exports,
				 int & SizeOfPacket,
				 int * Sizes,
				 bool & VarSizes,
				 std::vector<int>& pids)
{
  int i, j;

  VarSizes = true; //enable variable block size data comm

  int TotalSendLength = 0;
  int * IntSizes = 0; 
  if( NumExportIDs>0 ) IntSizes = new int[NumExportIDs];

  int SizeofIntType = sizeof(int_type);

  for(i = 0; i < NumExportIDs; ++i) {
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
  if( TotalSendLength*SizeOfPacket > LenExports ) {
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

  const Epetra_Map & rowMap = A.RowMap();
  const Epetra_Map & colMap = A.ColMap();
  
  if( NumExportIDs > 0 ) {
    int_type * Indices;
    int_type FromRow; 
    int_type * intptr;                         
    
    int maxNumEntries = A.MaxNumEntries();
    std::vector<int> MyIndices(maxNumEntries);
    dintptr = (double *) Exports;
    valptr = dintptr + IntSizes[0];
    intptr = (int_type *) dintptr;
    for (i=0; i<NumExportIDs; i++) {
      FromRow   = (int_type) rowMap.GID64(ExportLIDs[i]);
      intptr[0] = FromRow;
      values    = valptr;
      Indices   = intptr + 2;
      EPETRA_CHK_ERR(A.ExtractMyRowCopy(ExportLIDs[i], maxNumEntries, NumEntries, values, &MyIndices[0]));
      for (j=0; j<NumEntries; j++) {
	Indices[2*j]   = (int_type) colMap.GID64(MyIndices[j]);   // convert to GIDs
	Indices[2*j+1] = pids[MyIndices[j]];                      // PID owning the entry.
      }
      intptr[1] = NumEntries; // Load second slot of segment
      if( i < (NumExportIDs-1) ) {
	dintptr += (IntSizes[i]+Sizes[i]);
	valptr = dintptr + IntSizes[i+1];
	intptr = (int_type *) dintptr;
      }	
    }    

    for(i = 0; i < NumExportIDs; ++i )
      Sizes[i] += IntSizes[i];
  }
  
  if( IntSizes ) delete [] IntSizes;
  
  return(0);
}


// =========================================================================
// =========================================================================
int PackAndPrepareWithOwningPIDs(const Epetra_CrsMatrix & SourceMatrix, 
				 int NumExportIDs,
				 int * ExportLIDs,
				 int & LenExports,
				 char *& Exports,
				 int & SizeOfPacket,
				 int * Sizes,
				 bool & VarSizes,
				 std::vector<int>& SourcePids) 
{
  if(SourceMatrix.RowMap().GlobalIndicesInt()) {
    return TPackAndPrepareWithOwningPIDs<int>(SourceMatrix,NumExportIDs,ExportLIDs,LenExports,Exports,SizeOfPacket,Sizes,VarSizes,SourcePids);
  }
  else if(SourceMatrix.RowMap().GlobalIndicesLongLong()) {
    return TPackAndPrepareWithOwningPIDs<long long>(SourceMatrix,NumExportIDs,ExportLIDs,LenExports,Exports,SizeOfPacket,Sizes,VarSizes,SourcePids);
  }
  else
    throw std::runtime_error("UnpackWithOwningPIDsCount: Unable to determine source global index type");
}


// =========================================================================
// =========================================================================
template<typename int_type>
int TUnpackWithOwningPIDsCount(const Epetra_CrsMatrix& SourceMatrix, 
			       int NumSameIDs,
			       int NumRemoteIDs,
			       const int * RemoteLIDs,
			       int NumPermuteIDs,
			       const int *PermuteToLIDs,
			       const int *PermuteFromLIDs,
			       int LenImports,
			       char* Imports)
{  
  int i,nnz=0;
  int SizeofIntType = (int)sizeof(int_type);

  // SameIDs: Always first, always in the same place
  for(i=0; i<NumSameIDs; i++)
    nnz+=SourceMatrix.NumMyEntries(i);
  
  // PermuteIDs: Still local, but reordered
  for(i=0; i<NumPermuteIDs; i++)
    nnz += SourceMatrix.NumMyEntries(PermuteFromLIDs[i]);
  
  // RemoteIDs:  RemoteLIDs tells us the ID, we need to look up the length the hard way.  See UnpackAndCombine for where this code came from
  if(NumRemoteIDs > 0) {
    double * dintptr = (double *) Imports;
    // General version
    int_type *  intptr = (int_type *) dintptr;
    int   NumEntries = (int) intptr[1];
    int      IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
    for(i=0; i<NumRemoteIDs; i++) {
      nnz += NumEntries;
      
      if( i < (NumRemoteIDs-1) ) {
	dintptr   += IntSize + NumEntries;
	intptr     = (int_type *) dintptr;
	NumEntries = (int) intptr[1];
	IntSize    = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
      }
    }
  }

  return nnz;
}

  
// =========================================================================
// =========================================================================
int UnpackWithOwningPIDsCount(const Epetra_CrsMatrix& SourceMatrix, 
			      int NumSameIDs,
			      int NumRemoteIDs,
			      const int * RemoteLIDs,
			      int NumPermuteIDs,
			      const int *PermuteToLIDs,
			      const int *PermuteFromLIDs,
			      int LenImports,
			      char* Imports)	
{
  if(SourceMatrix.RowMap().GlobalIndicesInt()) {
    return TUnpackWithOwningPIDsCount<int>(SourceMatrix,NumSameIDs,NumRemoteIDs,RemoteLIDs,NumPermuteIDs,
					   PermuteToLIDs,PermuteFromLIDs,LenImports,Imports);
  }
  else if(SourceMatrix.RowMap().GlobalIndicesLongLong()) {
   return TUnpackWithOwningPIDsCount<long long>(SourceMatrix,NumSameIDs,NumRemoteIDs,RemoteLIDs,NumPermuteIDs,
						PermuteToLIDs,PermuteFromLIDs,LenImports,Imports);
  }
  else
    throw std::runtime_error("UnpackWithOwningPIDsCount: Unable to determine source global index type");

}

				
// =========================================================================
// =========================================================================
template<typename int_type>
int TUnpackAndCombineIntoCrsArrays(const Epetra_CrsMatrix& SourceMatrix, 
				  int NumSameIDs,
				  int NumRemoteIDs,
				  const int * RemoteLIDs,
				  int NumPermuteIDs,
				  const int *PermuteToLIDs,
				  const int *PermuteFromLIDs,
				  int LenImports,
				  char* Imports,
				  int TargetNumRows,
				  int TargetNumNonzeros,
				  int * CSR_rowptr,
				  int_type * CSR_colind,
				  double * CSR_vals,
				  const std::vector<int> &SourcePids,
				  std::vector<int> &TargetPids)
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
  int N = TargetNumRows;
  int mynnz = TargetNumNonzeros;
  int MyPID = SourceMatrix.Comm().MyPID();

  int SizeofIntType = sizeof(int_type);

  // Zero the rowptr
  for(i=0; i<N+1; i++) CSR_rowptr[i]=0;

  // SameIDs: Always first, always in the same place
  for(i=0; i<NumSameIDs; i++)
    CSR_rowptr[i]=SourceMatrix.NumMyEntries(i);
  
  // PermuteIDs: Still local, but reordered
  for(i=0; i<NumPermuteIDs; i++)
    CSR_rowptr[PermuteToLIDs[i]] = SourceMatrix.NumMyEntries(PermuteFromLIDs[i]);
  
  // RemoteIDs:  RemoteLIDs tells us the ID, we need to look up the length the hard way.  See UnpackAndCombine for where this code came from
  if(NumRemoteIDs > 0) {
    double * dintptr = (double *) Imports;    
    int_type    *  intptr = (int_type *) dintptr;
    int   NumEntries = (int) intptr[1];
    int      IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
    for(i=0; i<NumRemoteIDs; i++) {
      CSR_rowptr[RemoteLIDs[i]] += NumEntries;
      
      if( i < (NumRemoteIDs-1) ) {
	dintptr   += IntSize + NumEntries;
	intptr     = (int_type *) dintptr;
	NumEntries = (int) intptr[1];
	IntSize    = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
      }
    }
  }

  // If multiple procs contribute to a row;
  std::vector<int> NewStartRow(N+1);
 
  // Turn row length into a real CSR_rowptr
  int last_len = CSR_rowptr[0];
  CSR_rowptr[0] = 0;
  for(i=1; i<N+1; i++){
    int new_len    = CSR_rowptr[i];
    CSR_rowptr[i]  = last_len + CSR_rowptr[i-1];
    NewStartRow[i] = CSR_rowptr[i];
    last_len       = new_len;
  }

  // Preseed TargetPids with -1 for local
  if(TargetPids.size()!=(size_t)mynnz) TargetPids.resize(mynnz);
  TargetPids.assign(mynnz,-1);
 
  // Grab pointers for SourceMatrix
  int    * Source_rowptr, * Source_colind;
  double * Source_vals;
  rv=SourceMatrix.ExtractCrsDataPointers(Source_rowptr,Source_colind,Source_vals);
  if(rv) throw std::runtime_error("UnpackAndCombineIntoCrsArrays: failed in ExtractCrsDataPointers");

  // SameIDs: Copy the data over
  for(i=0; i<NumSameIDs; i++) {
    int FromRow = Source_rowptr[i];
    int ToRow   = CSR_rowptr[i];
    NewStartRow[i] += Source_rowptr[i+1]-Source_rowptr[i];

    for(j=Source_rowptr[i]; j<Source_rowptr[i+1]; j++) {
      CSR_vals[ToRow + j - FromRow]   = Source_vals[j];      
      CSR_colind[ToRow + j - FromRow] = (int_type) SourceMatrix.GCID64(Source_colind[j]);
      TargetPids[ToRow + j - FromRow] = (SourcePids[Source_colind[j]] != MyPID) ? SourcePids[Source_colind[j]] : -1;
    }
  }

  // PermuteIDs: Copy the data over
  for(i=0; i<NumPermuteIDs; i++) {
    int FromLID = PermuteFromLIDs[i];
    int FromRow = Source_rowptr[FromLID];
    int ToRow   = CSR_rowptr[PermuteToLIDs[i]];

    NewStartRow[PermuteToLIDs[i]] += Source_rowptr[FromLID+1]-Source_rowptr[FromLID];

    for(j=Source_rowptr[FromLID]; j<Source_rowptr[FromLID+1]; j++) {
      CSR_vals[ToRow + j - FromRow]   = Source_vals[j];      
      CSR_colind[ToRow + j - FromRow] = (int_type) SourceMatrix.GCID64(Source_colind[j]);
      TargetPids[ToRow + j - FromRow] = (SourcePids[Source_colind[j]] != MyPID) ? SourcePids[Source_colind[j]] : -1;
    }
  }

  // RemoteIDs: Loop structure following UnpackAndCombine  
  if(NumRemoteIDs > 0) {
    double * dintptr   = (double *) Imports;
    int_type * intptr  = (int_type *) dintptr;
    int NumEntries     = (int) intptr[1];
    int IntSize        = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
    double* valptr     = dintptr + IntSize;
      
    for (i=0; i<NumRemoteIDs; i++) {
      int ToLID    = RemoteLIDs[i];
      int StartRow = NewStartRow[ToLID];
      NewStartRow[ToLID]+=NumEntries;	
      
      double * values    = valptr;
      int_type * Indices = intptr + 2;
      for(j=0; j<NumEntries; j++){
	CSR_colind[StartRow + j]  = Indices[2*j];
	if(MyPID !=  Indices[2*j+1]) TargetPids[StartRow+j]   = Indices[2*j+1];
	CSR_vals[StartRow + j]   = values[j];	  
      }
      if( i < (NumRemoteIDs-1) ) {
	dintptr   += IntSize + NumEntries;
	intptr     = (int_type *) dintptr;
	NumEntries = (int) intptr[1];
	IntSize    = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
	valptr     = dintptr + IntSize;
      }
    }
  }

  return 0;
}



// =========================================================================
// =========================================================================
int UnpackAndCombineIntoCrsArrays(const Epetra_CrsMatrix& SourceMatrix, 
				  int NumSameIDs,
				  int NumRemoteIDs,
				  const int * RemoteLIDs,
				  int NumPermuteIDs,
				  const int *PermuteToLIDs,
				  const int *PermuteFromLIDs,
				  int LenImports,
				  char* Imports,
				  int TargetNumRows,
				  int TargetNumNonzeros,
				  int * CSR_rowptr,
				  int * CSR_colind,
				  double * CSR_values,
				  const std::vector<int> &SourcePids,
				  std::vector<int> &TargetPids)
{
  if(SourceMatrix.RowMap().GlobalIndicesInt()) {
    return TUnpackAndCombineIntoCrsArrays<int>(SourceMatrix,NumSameIDs,NumRemoteIDs,RemoteLIDs,NumPermuteIDs,
					       PermuteToLIDs,PermuteFromLIDs,LenImports,Imports,TargetNumRows,TargetNumNonzeros,
					       CSR_rowptr,CSR_colind,CSR_values,SourcePids,TargetPids);
  }
  else 
    throw std::runtime_error("UnpackAndCombineIntoCrsArrays: int type not matched.");
}

// =========================================================================
// =========================================================================
int UnpackAndCombineIntoCrsArrays(const Epetra_CrsMatrix& SourceMatrix, 
				  int NumSameIDs,
				  int NumRemoteIDs,
				  const int * RemoteLIDs,
				  int NumPermuteIDs,
				  const int *PermuteToLIDs,
				  const int *PermuteFromLIDs,
				  int LenImports,
				  char* Imports,
				  int TargetNumRows,
				  int TargetNumNonzeros,
				  int * CSR_rowptr,
				  long long * CSR_colind,
				  double * CSR_values,
				  const std::vector<int> &SourcePids,
				  std::vector<int> &TargetPids)
{
  if(SourceMatrix.RowMap().GlobalIndicesLongLong()) {
    return TUnpackAndCombineIntoCrsArrays<long long>(SourceMatrix,NumSameIDs,NumRemoteIDs,RemoteLIDs,NumPermuteIDs,
						     PermuteToLIDs,PermuteFromLIDs,LenImports,Imports,TargetNumRows,TargetNumNonzeros,
						     CSR_rowptr,CSR_colind,CSR_values,SourcePids,TargetPids);
  }
  else
    throw std::runtime_error("UnpackAndCombineIntoCrsArrays: int type not matched.");
}


// =========================================================================
// =========================================================================
  template<typename int_type, class MapType1, class MapType2>
 int TLowCommunicationMakeColMapAndReindex(int N, const int * rowptr, int * colind_LID, const int_type *colind_GID, const Epetra_Map& domainMap, const int * owningPIDs, bool SortGhostsAssociatedWithEachProcessor, std::vector<int>& RemotePIDs, MapType1 & NewColMap)
   {
  int i,j;


  // Sanity checks
  bool UseLL;
  if(domainMap.GlobalIndicesLongLong()) UseLL=true;
  else if(domainMap.GlobalIndicesInt()) UseLL=false;
  else throw std::runtime_error("LowCommunicationMakeColMapAndReindex: cannot detect int type.");

  // Scan all column indices and sort into two groups: 
  // Local:  those whose GID matches a GID of the domain map on this processor and
  // Remote: All others.
  int numDomainElements = domainMap.NumMyElements();
  bool * LocalGIDs  = 0;
  if (numDomainElements>0) LocalGIDs  = new bool[numDomainElements];
  for (i=0; i<numDomainElements; i++) LocalGIDs[i] = false; // Assume domain GIDs are not local

  bool DoSizes = !domainMap.ConstantElementSize(); // If not constant element size, then error
  if(DoSizes) throw std::runtime_error("LowCommunicationMakeColMapAndReindex: cannot handle non-constant sized domainMap.");


  // In principle it is good to have RemoteGIDs and RemotGIDList be as long as the number of remote GIDs
  // on this processor, but this would require two passes through the column IDs, so we make it the max of 100
  // and the number of block rows.
  const int numMyBlockRows = N;
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
	colind_LID[j] = LID; 
      }
      else {
	int_type hash_value=RemoteGIDs.Get(GID);
	if(hash_value  == -1) { // This means its a new remote GID
	  int PID = owningPIDs[j];
	  if(PID==-1) throw std::runtime_error("LowCommunicationMakeColMapAndReindex: Cannot figure out if PID is owned.");
	  colind_LID[j] = numDomainElements + NumRemoteColGIDs;
	  RemoteGIDs.Add(GID, NumRemoteColGIDs);
	  RemoteGIDList.push_back(GID);
	  PIDList.push_back(PID);
	  NumRemoteColGIDs++;
	}
	else
	  colind_LID[j] = numDomainElements + hash_value;	  
      }
    }
  }

  // Possible short-circuit:  If all domain map GIDs are present as column indices, then set ColMap=domainMap and quit
  if (domainMap.Comm().NumProc()==1) { 
    
    if (NumRemoteColGIDs!=0) {
      throw std::runtime_error("Some column IDs are not in domainMap.  If matrix is rectangular, you must pass in a domainMap"); 
      // Sanity test: When one processor,there can be no remoteGIDs
    }
    if (NumLocalColGIDs==numDomainElements) {
      if (LocalGIDs!=0) delete [] LocalGIDs; 
      // In this case, we just use the domainMap's indices, which is, not coincidently, what we clobbered colind with up above anyway. 
      // No further reindexing is needed.
      NewColMap = domainMap;
      return 0;
    }
  }
      
  // Now build integer array containing column GIDs
  // Build back end, containing remote GIDs, first
  int numMyBlockCols = NumLocalColGIDs + NumRemoteColGIDs;
  std::vector<int_type> ColIndices;
  int_type * RemoteColIndices=0;
  if(numMyBlockCols > 0) {
    ColIndices.resize(numMyBlockCols);
    if(NumLocalColGIDs!=numMyBlockCols) RemoteColIndices = &ColIndices[NumLocalColGIDs]; // Points to back end of ColIndices
    else RemoteColIndices=0;
  }

  for(i = 0; i < NumRemoteColGIDs; i++) 
    RemoteColIndices[i] = RemoteGIDList[i]; 

  // Build permute array for *remote* reindexing.
  std::vector<int> RemotePermuteIDs(NumRemoteColGIDs);
  for(i=0; i<NumRemoteColGIDs; i++) RemotePermuteIDs[i]=i;

  // Sort External column indices so that all columns coming from a given remote processor are contiguous
  int NumListsInt=0;
  int NumListsLL =0;
  int * IntSortLists[2];
  long long * LLSortLists[2];
  int * RemotePermuteIDs_ptr = RemotePermuteIDs.size() ? &RemotePermuteIDs[0] : 0;
  if(!UseLL) {
    // int version
    IntSortLists[0] = (int*) RemoteColIndices;
    IntSortLists[1] = RemotePermuteIDs_ptr;
    NumListsInt=2;
  }
  else {
    //LL version
    LLSortLists[0]  = (long long*) RemoteColIndices;
    IntSortLists[0] = RemotePermuteIDs_ptr;
    NumListsInt     = NumListsLL = 1;
  }

  int * PIDList_ptr = PIDList.size() ? &PIDList[0] : 0;
  Epetra_Util::Sort(true, NumRemoteColGIDs, PIDList_ptr, 0, 0, NumListsInt, IntSortLists,NumListsLL,LLSortLists);

  // Stash the RemotePIDs  
  PIDList.resize(NumRemoteColGIDs);
  RemotePIDs = PIDList;

  if (SortGhostsAssociatedWithEachProcessor) {
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
    int_type* MyGlobalElements = 0;
    domainMap.MyGlobalElementsPtr(MyGlobalElements);

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
  int_type * ColIndices_ptr  = ColIndices.size() ? &ColIndices[0] : 0;
  MapType2 temp((int_type)(-1), numMyBlockCols, ColIndices_ptr, (int_type)domainMap.IndexBase64(), domainMap.Comm());
  NewColMap = temp;

  // Low-cost reindex of the matrix
  for(i=0; i<numMyBlockRows; i++){
    for(j=rowptr[i]; j<rowptr[i+1]; j++){
      int ID=colind_LID[j];
      if(ID < numDomainElements){
	if(use_local_permute) colind_LID[j] = LocalPermuteIDs[colind_LID[j]];
	// In the case where use_local_permute==false, we just copy the DomainMap's ordering, which it so happens
	// is what we put in colind to begin with.
      }
      else
	colind_LID[j] =  NumLocalColGIDs + ReverseRemotePermuteIDs[colind_LID[j]-numDomainElements];
    }
  }
  
  return 0;
}


// =========================================================================
// =========================================================================
int LowCommunicationMakeColMapAndReindex(int N, const int * rowptr, int * colind, const Epetra_Map& domainMap, const int * owningPIDs, bool SortGhostsAssociatedWithEachProcessor, std::vector<int>& RemotePIDs, Epetra_BlockMap & NewColMap) {

  
  if(domainMap.GlobalIndicesInt()) 
    return TLowCommunicationMakeColMapAndReindex<int,Epetra_BlockMap,Epetra_Map>(N,rowptr,colind,colind,domainMap,owningPIDs,SortGhostsAssociatedWithEachProcessor,RemotePIDs,NewColMap); 
  else
    throw std::runtime_error("LowCommunicationMakeColMapAndReindex: GID type mismatch.");
  return -1;
}

// =========================================================================
// =========================================================================
int LowCommunicationMakeColMapAndReindex(int N, const int * rowptr, int * colind_LID, long long * colind_GID, const Epetra_Map& domainMap, const int * owningPIDs, bool SortGhostsAssociatedWithEachProcessor, std::vector<int>& RemotePIDs, Epetra_BlockMap & NewColMap) {

  
  if(domainMap.GlobalIndicesLongLong()) 
    return TLowCommunicationMakeColMapAndReindex<long long,Epetra_BlockMap,Epetra_Map>(N,rowptr,colind_LID,colind_GID,domainMap,owningPIDs,SortGhostsAssociatedWithEachProcessor,RemotePIDs,NewColMap); 
  else
    throw std::runtime_error("LowCommunicationMakeColMapAndReindex: GID type mismatch.");
  return -1;
}


}// end namespace Epetra_Import_Util
