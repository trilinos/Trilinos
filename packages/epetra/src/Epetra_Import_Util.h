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

#ifndef EPETRA_IMPORT_UTIL_H
#define EPETRA_IMPORT_UTIL_H

#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include <vector>
#include <stdexcept>
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;

//! Epetra_Import_Util:  The Epetra ImportUtil Wrapper Namespace.
/*! The Epetra_Import_Util namepsace is a collection of useful functions for data import.
   The goal of this is to rmevoe code duplication between Epetra and EpetraExt.

*/

namespace Epetra_Import_Util { 

//=========================================================================
//! PackAndPrepareWithOwningPIDs.  
 /*! Note: The pids vector should contain a list of owning PIDs for each column in the ColMap, as from Epetra_Util::GetPids,
   without the "-1 for local" option being used.

   \warning This method is intended for expert developer use only, and should never be called by user code.
   */
template<typename MatrixType>
int PackAndPrepareWithOwningPIDs(const MatrixType & A, 
				 int NumExportIDs,
				 int * ExportLIDs,
				 int & LenExports,
				 char *& Exports,
				 int & SizeOfPacket,
				 int * Sizes,
				 bool & VarSizes,
				 std::vector< int > pids)
{

  VarSizes = true; //enable variable block size data comm

  int TotalSendLength = 0;
  int * IntSizes = 0; 
  if( NumExportIDs>0 ) IntSizes = new int[NumExportIDs];

  int SizeofIntType = -1;
  if(A.RowMap().GlobalIndicesInt())
    SizeofIntType = (int)sizeof(int); 
  else if(A.RowMap().GlobalIndicesLongLong())
    SizeofIntType = (int)sizeof(long long); 
  else
    throw std::runtime_error("Epetra_Import_Util::PackAndPrepareWithOwningPIDs: Unable to determine source global index type");

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

  const Epetra_Map & rowMap = A.RowMap();
  const Epetra_Map & colMap = A.ColMap();

  if( NumExportIDs > 0 )
  {

    if(A.RowMap().GlobalIndicesInt()) { 
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
    else if(A.RowMap().GlobalIndicesLongLong()) {
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



}



#endif /* EPETRA_IMPORT_UTIL_H */
