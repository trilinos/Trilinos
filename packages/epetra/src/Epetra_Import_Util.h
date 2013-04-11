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
#include "Epetra_CrsMatrix.h"
#include <vector>
class Epetra_Map;
class Epetra_CrsMatrix;
class Epetra_Import;

//! Epetra_Import_Util:  The Epetra ImportUtil Wrapper Namespace.
/*! The Epetra_Import_Util namepsace is a collection of useful functions for data import.
   The goal of this is to rmevoe code duplication between Epetra and EpetraExt.

*/

namespace Epetra_Import_Util { 

//=========================================================================
//! PackAndPrepareWithOwningPIDs.  
 /*! Note: The SourcePids vector should contain a list of owning PIDs for each column in the ColMap, as from Epetra_Util::GetPids,
   without the "-1 for local" option being used.

   \warning This method is intended for expert developer use only, and should never be called by user code.
 */
int PackAndPrepareWithOwningPIDs(const Epetra_CrsMatrix & SourceMatrix, 
				 int NumExportIDs,
				 int * ExportLIDs,
				 int & LenExports,
				 char *& Exports,
				 int & SizeOfPacket,
				 int * Sizes,
				 bool & VarSizes,
				 std::vector<int>& SourcePids);

// ===================================================================
//! UnpackWithOwningPIDsCount
/*! Counts the number of non-zeros in the resulting matrix.  Call this before UnpackAndCombineIntoCrsArrays.

   \warning This method is intended for expert developer use only, and should never be called by user code.
*/
int UnpackWithOwningPIDsCount(const Epetra_CrsMatrix& SourceMatrix, 
			      int NumSameIDs,
			      int NumRemoteIDs,
			      const int * RemoteLIDs,
			      int NumPermuteIDs,
			      const int *PermuteToLIDs,
			      const int *PermuteFromLIDs,
			      int LenImports,
			      char* Imports);

// ===================================================================
//! UnpackAndCombineIntoCrsArrays
/*! You should call UnpackWithOwningPIDsCount first and allocate all arrays accordingly.
   
   Note: The SourcePids vector (on input) should contain of owning PIDs for each column in the ColMap, as from Epetra_Util::GetPids,
   with the "-1 for local" option being used.  

   Note: The TargetPids vector (on output) will contain of owning PIDs for each column in the ColMap, as from Epetra_Util::GetPids,
   with the "-1 for local" option being used.  

   \warning This method is intended for expert developer use only, and should never be called by user code.
   */
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
				  std::vector<int> &TargetPids);
// ===================================================================
//! UnpackAndCombineIntoCrsArrays
/*! You should call UnpackWithOwningPIDsCount first and allocate all arrays accordingly.
   
   Note: The SourcePids vector (on input) should contain of owning PIDs for each column in the ColMap, as from Epetra_Util::GetPids,
   with the "-1 for local" option being used.  

   Note: The TargetPids vector (on output) will contain of owning PIDs for each column in the ColMap, as from Epetra_Util::GetPids,
   with the "-1 for local" option being used.  

   \warning This method is intended for expert developer use only, and should never be called by user code.
   */
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
				  std::vector<int> &TargetPids);



// ===================================================================
//! LowCommunicationMakeColMapAndReindex
/*! If you know the owning PIDs already, you can make the colmap a lot less expensively.
   
   Note: The owningPids vector (on input) should contain of owning PIDs for each column in the ColMap, as from Epetra_Util::GetPids,
   with the "-1 for local" can be used, or not, here.

   Note: This method will return a std::vector of the RemotePIDs, used for construction of the importer.

   \warning This method is intended for expert developer use only, and should never be called by user code.
*/
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int LowCommunicationMakeColMapAndReindex(int N, const int *rowptr, int *colind, const Epetra_Map& domainMap, 
					 const int *owningPIDs, bool SortGhostsAssociatedWithEachProcessor, 
					 std::vector<int>& RemotePIDs, Epetra_BlockMap & NewColMap);
#endif 
// ===================================================================
//! LowCommunicationMakeColMapAndReindex
/*! If you know the owning PIDs already, you can make the colmap a lot less expensively.
   
   Note: The owningPids vector (on input) should contain of owning PIDs for each column in the ColMap, as from Epetra_Util::GetPids,
   with the "-1 for local" can be used, or not, here.

   Note: This method will return a std::vector of the RemotePIDs, used for construction of the importer.

   \warning This method is intended for expert developer use only, and should never be called by user code.
*/
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int LowCommunicationMakeColMapAndReindex(int N, const int *rowptr, int *colind_LID, long long *colind_GID, const Epetra_Map& domainMap, 
					 const int *owningPIDs, bool SortGhostsAssociatedWithEachProcessor, std::vector<int>& RemotePIDs, Epetra_BlockMap & NewColMap);
#endif

} /* Epetra_Import_Util namespace */

#endif /* EPETRA_IMPORT_UTIL_H */
