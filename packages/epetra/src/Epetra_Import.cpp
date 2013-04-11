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
#include "Epetra_Import.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Distributor.h"
#include "Epetra_Comm.h"
#include "Epetra_Util.h"
#include "Epetra_Export.h"

#include <algorithm>
#include <vector>


#ifdef HAVE_MPI
#include "Epetra_MpiDistributor.h"
#endif

#ifdef EPETRA_ENABLE_DEBUG
#include "Epetra_IntVector.h"
#endif

//==============================================================================
// Epetra_Import constructor function for a Epetra_BlockMap object
template<typename int_type>
void Epetra_Import::Construct_Expert( const Epetra_BlockMap &  targetMap, const Epetra_BlockMap & sourceMap, int NumRemotePIDs, const int * UserRemotePIDs,
				      const int & UserNumExportIDs, const int * UserExportLIDs,  const int * UserExportPIDs)
{
  int i,ierr;
  // Build three ID lists:
  // NumSameIDs - Number of IDs in TargetMap and SourceMap that are identical, up to the first
  //              nonidentical ID.
  // NumPermuteIDs - Number of IDs in SourceMap that must be indirectly loaded but are on this processor.
  // NumRemoteIDs - Number of IDs that are in SourceMap but not in TargetMap, and thus must be imported.
  
  int NumSourceIDs = sourceMap.NumMyElements();
  int NumTargetIDs = targetMap.NumMyElements();
  
  int_type *TargetGIDs = 0;
  if (NumTargetIDs>0) {
    TargetGIDs = new int_type[NumTargetIDs];
    targetMap.MyGlobalElements(TargetGIDs);
  }
  
  int_type * SourceGIDs = 0;
  if (NumSourceIDs>0) {
    SourceGIDs = new int_type[NumSourceIDs];
    sourceMap.MyGlobalElements(SourceGIDs);
  }
  
  int MinIDs = EPETRA_MIN(NumSourceIDs, NumTargetIDs);
    
  NumSameIDs_ = 0;
  for (i=0; i< MinIDs; i++) if (TargetGIDs[i]==SourceGIDs[i]) NumSameIDs_++; else break;
  
  // Find count of Target IDs that are truly remote and those that are local but permuted
  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) 
    if (sourceMap.MyGID(TargetGIDs[i])) NumPermuteIDs_++; // Check if Target GID is a local Source GID
    else NumRemoteIDs_++; // If not, then it is remote
     
  // Define remote and permutation lists  
  int_type * RemoteGIDs=0;
  RemoteLIDs_ = 0;
  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    RemoteGIDs = new int_type[NumRemoteIDs_];
  }
  if (NumPermuteIDs_>0)  {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
  }
  
  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) {
    if (sourceMap.MyGID(TargetGIDs[i])) {
      PermuteToLIDs_[NumPermuteIDs_] = i;
      PermuteFromLIDs_[NumPermuteIDs_++] = sourceMap.LID(TargetGIDs[i]);
    }
    else {
      //NumRecv_ +=TargetMap.ElementSize(i); // Count total number of entries to receive
      NumRecv_ +=targetMap.MaxElementSize(); // Count total number of entries to receive (currently need max)
      RemoteGIDs[NumRemoteIDs_] = TargetGIDs[i];
      RemoteLIDs_[NumRemoteIDs_++] = i;
    }
  }

  if( NumRemoteIDs_>0 && !sourceMap.DistributedGlobal() )
    ReportError("Warning in Epetra_Import: Serial Import has remote IDs. (Importing to Subset of Target Map)", 1);
  
  // Test for distributed cases
  int * RemotePIDs = 0;
  
  if (sourceMap.DistributedGlobal()) {
    if (NumRemoteIDs_>0)  RemotePIDs = new int[NumRemoteIDs_];
  
#ifdef EPETRA_ENABLE_DEBUG
    int myeq = (NumRemotePIDs==NumRemoteIDs_);
    int globaleq=0;
    sourceMap.Comm().MinAll(&myeq,&globaleq,1);
    if(globaleq!=1) { 
      printf("[%d] UserRemotePIDs count wrong %d != %d\n",sourceMap.Comm().MyPID(),NumRemotePIDs,NumRemoteIDs_);
      fflush(stdout);
      sourceMap.Comm().Barrier();
      sourceMap.Comm().Barrier();
      sourceMap.Comm().Barrier();
      throw ReportError("Epetra_Import: UserRemotePIDs count wrong");
    }
#endif

    if(NumRemotePIDs==NumRemoteIDs_){
      // Since I need to sort these, I'll copy them
      for(i=0; i<NumRemoteIDs_; i++)  RemotePIDs[i] = UserRemotePIDs[i];
    }

    //Get rid of IDs that don't exist in SourceMap
    if(NumRemoteIDs_>0) {
      int cnt = 0;
      for( i = 0; i < NumRemoteIDs_; ++i )
        if( RemotePIDs[i] == -1 ) ++cnt;
      if( cnt ) {
        if( NumRemoteIDs_-cnt ) {
          int_type * NewRemoteGIDs = new int_type[NumRemoteIDs_-cnt];
          int * NewRemotePIDs = new int[NumRemoteIDs_-cnt];
          int * NewRemoteLIDs = new int[NumRemoteIDs_-cnt];
          cnt = 0;
          for( i = 0; i < NumRemoteIDs_; ++i )
            if( RemotePIDs[i] != -1 ) {
              NewRemoteGIDs[cnt] = RemoteGIDs[i];
              NewRemotePIDs[cnt] = RemotePIDs[i];
              NewRemoteLIDs[cnt] = targetMap.LID(RemoteGIDs[i]);
              ++cnt;
            }
          NumRemoteIDs_ = cnt;
          delete [] RemoteGIDs;
          delete [] RemotePIDs;
          delete [] RemoteLIDs_;
          RemoteGIDs = NewRemoteGIDs;
          RemotePIDs = NewRemotePIDs;
          RemoteLIDs_ = NewRemoteLIDs;
          ReportError("Warning in Epetra_Import: Target IDs not found in Source Map (Do you want to import to subset of Target Map?)", 1);
        }
        else { //valid RemoteIDs empty
          NumRemoteIDs_ = 0;
          delete [] RemoteGIDs;
          RemoteGIDs = 0;
          delete [] RemotePIDs;
          RemotePIDs = 0;
        }
      }
    }

    //Sort Remote IDs by processor so DoReverses will work
    Epetra_Util util;
    
    if(targetMap.GlobalIndicesLongLong())
      {
	util.Sort(true,NumRemoteIDs_,RemotePIDs,0,0, 1,&RemoteLIDs_, 1,(long long**)&RemoteGIDs);
      }
    else if(targetMap.GlobalIndicesInt())
      {
	int* ptrs[2] = {RemoteLIDs_, (int*)RemoteGIDs};
	util.Sort(true,NumRemoteIDs_,RemotePIDs,0,0,2,&ptrs[0], 0, 0);
      }
    else
      {
	throw ReportError("Epetra_Import::Epetra_Import: GlobalIndices Internal Error", -1);
      }
    
    // Build distributor & Export lists
    Distor_ = sourceMap.Comm().CreateDistributor();    
    
    NumExportIDs_=UserNumExportIDs;
    ExportLIDs_ = new int[NumExportIDs_];
    ExportPIDs_ = new int[NumExportIDs_];
    for(i=0; i<NumExportIDs_; i++)  {
      ExportPIDs_[i] = UserExportPIDs[i];
      ExportLIDs_[i] = UserExportLIDs[i];
    }

#ifdef HAVE_MPI
    Epetra_MpiDistributor* MpiDistor = dynamic_cast< Epetra_MpiDistributor*>(Distor_);
    bool Deterministic = true;
    if(MpiDistor)
      ierr=MpiDistor->CreateFromSendsAndRecvs(NumExportIDs_,ExportPIDs_,					      
					      NumRemoteIDs_, RemoteGIDs, RemotePIDs,Deterministic);
    else ierr=-10;
#else
    ierr=-20;
#endif
    
    if (ierr!=0) throw ReportError("Error in Epetra_Distributor.CreateFromRecvs()", ierr);   
  }  

  if( NumRemoteIDs_>0 ) delete [] RemoteGIDs;
  if( NumRemoteIDs_>0 ) delete [] RemotePIDs;
  
  if (NumTargetIDs>0) delete [] TargetGIDs;
  if (NumSourceIDs>0) delete [] SourceGIDs;


#ifdef EPETRA_ENABLE_DEBUG
// Sanity check to make sure we got the import right
  Epetra_IntVector Source(sourceMap);
  Epetra_IntVector Target(targetMap);

  for(i=0; i<Source.MyLength(); i++)
    Source[i] = (int) (Source.Map().GID(i) % INT_MAX);
  Target.PutValue(-1);
 
  Target.Import(Source,*this,Insert);
  
  bool test_passed=true;
  for(i=0; i<Target.MyLength(); i++){
    if(Target[i] != Target.Map().GID(i) % INT_MAX) test_passed=false;
  }

  if(!test_passed) {  
    printf("[%d] PROCESSOR has a mismatch... prepearing to crash or hang!\n",sourceMap.Comm().MyPID());
    fflush(stdout);
    sourceMap.Comm().Barrier();
    sourceMap.Comm().Barrier();
    sourceMap.Comm().Barrier();
    throw ReportError("Epetra_Import: ERROR. User provided IDs do not match what an import generates.");
  }
#endif
  
  return;
}

//==============================================================================
// Epetra_Import constructor function for a Epetra_BlockMap object
template<typename int_type>
void Epetra_Import::Construct( const Epetra_BlockMap &  targetMap, const Epetra_BlockMap & sourceMap, int NumRemotePIDs, const int * UserRemotePIDs)
{
  int i,ierr;
  // Build three ID lists:
  // NumSameIDs - Number of IDs in TargetMap and SourceMap that are identical, up to the first
  //              nonidentical ID.
  // NumPermuteIDs - Number of IDs in SourceMap that must be indirectly loaded but are on this processor.
  // NumRemoteIDs - Number of IDs that are in SourceMap but not in TargetMap, and thus must be imported.
  
  int NumSourceIDs = sourceMap.NumMyElements();
  int NumTargetIDs = targetMap.NumMyElements();
  
  int_type *TargetGIDs = 0;
  if (NumTargetIDs>0) {
    TargetGIDs = new int_type[NumTargetIDs];
    targetMap.MyGlobalElements(TargetGIDs);
  }
  
  int_type * SourceGIDs = 0;
  if (NumSourceIDs>0) {
    SourceGIDs = new int_type[NumSourceIDs];
    sourceMap.MyGlobalElements(SourceGIDs);
  }
  
  int MinIDs = EPETRA_MIN(NumSourceIDs, NumTargetIDs);
  
  
  NumSameIDs_ = 0;
  for (i=0; i< MinIDs; i++) if (TargetGIDs[i]==SourceGIDs[i]) NumSameIDs_++; else break;
  
  // Find count of Target IDs that are truly remote and those that are local but permuted
  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) 
    if (sourceMap.MyGID(TargetGIDs[i])) NumPermuteIDs_++; // Check if Target GID is a local Source GID
    else NumRemoteIDs_++; // If not, then it is remote
  
  
  
  // Define remote and permutation lists  
  int_type * RemoteGIDs=0;
  RemoteLIDs_ = 0;
  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    RemoteGIDs = new int_type[NumRemoteIDs_];
  }
  if (NumPermuteIDs_>0)  {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
  }
  
  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) {
    if (sourceMap.MyGID(TargetGIDs[i])) {
      PermuteToLIDs_[NumPermuteIDs_] = i;
      PermuteFromLIDs_[NumPermuteIDs_++] = sourceMap.LID(TargetGIDs[i]);
    }
    else {
      //NumRecv_ +=TargetMap.ElementSize(i); // Count total number of entries to receive
      NumRecv_ +=targetMap.MaxElementSize(); // Count total number of entries to receive (currently need max)
      RemoteGIDs[NumRemoteIDs_] = TargetGIDs[i];
      RemoteLIDs_[NumRemoteIDs_++] = i;
    }
  }

  if( NumRemoteIDs_>0 && !sourceMap.DistributedGlobal() )
    ReportError("Warning in Epetra_Import: Serial Import has remote IDs. (Importing to Subset of Target Map)", 1);
  
  // Test for distributed cases
  int * RemotePIDs = 0;
  
  if (sourceMap.DistributedGlobal()) {
    if (NumRemoteIDs_>0)  RemotePIDs = new int[NumRemoteIDs_];

#ifdef EPETRA_ENABLE_DEBUG
    if(NumRemotePIDs!=-1){
      int myeq = (NumRemotePIDs==NumRemoteIDs_);
      int globaleq=0;
      sourceMap.Comm().MinAll(&myeq,&globaleq,1);
      if(globaleq!=1) { 
	printf("[%d] UserRemotePIDs count wrong %d != %d\n",sourceMap.Comm().MyPID(),NumRemotePIDs,NumRemoteIDs_);
	fflush(stdout);
	sourceMap.Comm().Barrier();
	sourceMap.Comm().Barrier();
	sourceMap.Comm().Barrier();
	throw ReportError("Epetra_Import: UserRemotePIDs count wrong",-1);
      }
    }
#endif

    if(NumRemotePIDs!=-1 && NumRemotePIDs==NumRemoteIDs_){
      // Since I need to sort these, I'll copy them
      for(i=0; i<NumRemoteIDs_; i++)  RemotePIDs[i] = UserRemotePIDs[i];
    }
    else{
      ierr=sourceMap.RemoteIDList(NumRemoteIDs_, RemoteGIDs, RemotePIDs, 0); // Get remote PIDs
      if (ierr) throw ReportError("Error in sourceMap.RemoteIDList call", ierr);
    }

    //Get rid of IDs that don't exist in SourceMap
    if(NumRemoteIDs_>0) {
      int cnt = 0;
      for( i = 0; i < NumRemoteIDs_; ++i )
        if( RemotePIDs[i] == -1 ) ++cnt;
      if( cnt ) {
        if( NumRemoteIDs_-cnt ) {
          int_type * NewRemoteGIDs = new int_type[NumRemoteIDs_-cnt];
          int * NewRemotePIDs = new int[NumRemoteIDs_-cnt];
          int * NewRemoteLIDs = new int[NumRemoteIDs_-cnt];
          cnt = 0;
          for( i = 0; i < NumRemoteIDs_; ++i )
            if( RemotePIDs[i] != -1 ) {
              NewRemoteGIDs[cnt] = RemoteGIDs[i];
              NewRemotePIDs[cnt] = RemotePIDs[i];
              NewRemoteLIDs[cnt] = targetMap.LID(RemoteGIDs[i]);
              ++cnt;
            }
          NumRemoteIDs_ = cnt;
          delete [] RemoteGIDs;
          delete [] RemotePIDs;
          delete [] RemoteLIDs_;
          RemoteGIDs = NewRemoteGIDs;
          RemotePIDs = NewRemotePIDs;
          RemoteLIDs_ = NewRemoteLIDs;
          ReportError("Warning in Epetra_Import: Target IDs not found in Source Map (Do you want to import to subset of Target Map?)", 1);
        }
        else { //valid RemoteIDs empty
          NumRemoteIDs_ = 0;
          delete [] RemoteGIDs;
          RemoteGIDs = 0;
          delete [] RemotePIDs;
          RemotePIDs = 0;
        }
      }
    }

    //Sort Remote IDs by processor so DoReverses will work
    if(targetMap.GlobalIndicesLongLong())
      {
	Epetra_Util::Sort(true,NumRemoteIDs_,RemotePIDs,0,0, 1,&RemoteLIDs_, 1,(long long**)&RemoteGIDs);
      }
    else if(targetMap.GlobalIndicesInt())
      {
	int* ptrs[2] = {RemoteLIDs_, (int*)RemoteGIDs};
	Epetra_Util::Sort(true,NumRemoteIDs_,RemotePIDs,0,0,2,&ptrs[0], 0, 0);
      }
    else
      {
	throw ReportError("Epetra_Import::Epetra_Import: GlobalIndices Internal Error", -1);
      }
    Distor_ = sourceMap.Comm().CreateDistributor();
    
    // Construct list of exports that calling processor needs to send as a result
    // of everyone asking for what it needs to receive.
    
    bool Deterministic = true;
    int_type* tmp_ExportLIDs; //Export IDs come in as GIDs
    ierr = Distor_->CreateFromRecvs( NumRemoteIDs_, RemoteGIDs, RemotePIDs,
				     Deterministic, NumExportIDs_, tmp_ExportLIDs, ExportPIDs_ );
    if (ierr!=0) throw ReportError("Error in Epetra_Distributor.CreateFromRecvs()", ierr);
    

    // Export IDs come in as GIDs, convert to LIDs
    if(targetMap.GlobalIndicesLongLong())
      {
	ExportLIDs_ = new int[NumExportIDs_];
	
	for (i=0; i< NumExportIDs_; i++) {
	  if (ExportPIDs_[i] < 0) throw ReportError("targetMap requested a GID that is not in the sourceMap.", -1);
	  ExportLIDs_[i] = sourceMap.LID(tmp_ExportLIDs[i]);
	  NumSend_ += sourceMap.MaxElementSize(); // Count total number of entries to send (currently need max)
	}
	
	delete[] tmp_ExportLIDs;
      }
    else if(targetMap.GlobalIndicesInt())
      {
	for (i=0; i< NumExportIDs_; i++) {
	  if (ExportPIDs_[i] < 0) throw ReportError("targetMap requested a GID that is not in the sourceMap.", -1);
	  tmp_ExportLIDs[i] = sourceMap.LID(tmp_ExportLIDs[i]);
	  NumSend_ += sourceMap.MaxElementSize(); // Count total number of entries to send (currently need max)
	}
	
	ExportLIDs_ = reinterpret_cast<int *>(tmp_ExportLIDs); // Can't reach here if tmp_ExportLIDs is long long.
      }
    else
      {
	throw ReportError("Epetra_Import::Epetra_Import: GlobalIndices Internal Error", -1);
      }
  }
  
  if( NumRemoteIDs_>0 ) delete [] RemoteGIDs;
  if( NumRemoteIDs_>0 ) delete [] RemotePIDs;

  if (NumTargetIDs>0) delete [] TargetGIDs;
  if (NumSourceIDs>0) delete [] SourceGIDs;

  return;
}

//==============================================================================
Epetra_Import::Epetra_Import( const Epetra_BlockMap &  targetMap, const Epetra_BlockMap & sourceMap,int NumRemotePIDs,const int * RemotePIDs)
  : Epetra_Object("Epetra::Import"),
    TargetMap_(targetMap),
    SourceMap_(sourceMap),
    NumSameIDs_(0),
    NumPermuteIDs_(0),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(0),
    RemoteLIDs_(0),
    NumExportIDs_(0),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(0),
    NumRecv_(0),
    Distor_(0)
{
  if(!targetMap.GlobalIndicesTypeMatch(sourceMap))
    throw ReportError("Epetra_Import::Epetra_Import: GlobalIndicesTypeMatch failed", -1);

  if(targetMap.GlobalIndicesInt())
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    Construct<int>(targetMap, sourceMap,NumRemotePIDs,RemotePIDs);
#else
    throw ReportError("Epetra_Import::Epetra_Import: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
  else if(targetMap.GlobalIndicesLongLong())
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    Construct<long long>(targetMap, sourceMap,NumRemotePIDs,RemotePIDs);
#else
    throw ReportError("Epetra_Import::Epetra_Import: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
  else
    throw ReportError("Epetra_Import::Epetra_Import: Bad global indices type", -1);
}



//==============================================================================
Epetra_Import::Epetra_Import( const Epetra_BlockMap &  targetMap, const Epetra_BlockMap & sourceMap,int NumRemotePIDs,const int * RemotePIDs,
			      const int & NumExportIDs, const int * ExportLIDs,  const int * ExportPIDs)
  : Epetra_Object("Epetra::Import"),
    TargetMap_(targetMap),
    SourceMap_(sourceMap),
    NumSameIDs_(0),
    NumPermuteIDs_(0),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(0),
    RemoteLIDs_(0),
    NumExportIDs_(0),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(0),
    NumRecv_(0),
    Distor_(0)
{
  if(!targetMap.GlobalIndicesTypeMatch(sourceMap))
    throw ReportError("Epetra_Import::Epetra_Import: GlobalIndicesTypeMatch failed", -1);

  if(targetMap.GlobalIndicesInt())
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    Construct_Expert<int>(targetMap, sourceMap,NumRemotePIDs,RemotePIDs,NumExportIDs,ExportLIDs,ExportPIDs);
#else
    throw ReportError("Epetra_Import::Epetra_Import: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
  else if(targetMap.GlobalIndicesLongLong())
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    Construct_Expert<long long>(targetMap, sourceMap,NumRemotePIDs,RemotePIDs,NumExportIDs,ExportLIDs,ExportPIDs);
#else
    throw ReportError("Epetra_Import::Epetra_Import: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
  else
    throw ReportError("Epetra_Import::Epetra_Import: Bad global indices type", -1);
}


//==============================================================================
Epetra_Import::Epetra_Import( const Epetra_BlockMap &  targetMap, const Epetra_BlockMap & sourceMap)
  : Epetra_Object("Epetra::Import"),
    TargetMap_(targetMap),
    SourceMap_(sourceMap),
    NumSameIDs_(0),
    NumPermuteIDs_(0),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(0),
    RemoteLIDs_(0),
    NumExportIDs_(0),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(0),
    NumRecv_(0),
    Distor_(0)
{
  if(!targetMap.GlobalIndicesTypeMatch(sourceMap))
    throw ReportError("Epetra_Import::Epetra_Import: GlobalIndicesTypeMatch failed", -1);

  if(targetMap.GlobalIndicesInt())
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    Construct<int>(targetMap, sourceMap);
#else
    throw ReportError("Epetra_Import::Epetra_Import: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
  else if(targetMap.GlobalIndicesLongLong())
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    Construct<long long>(targetMap, sourceMap);
#else
    throw ReportError("Epetra_Import::Epetra_Import: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
  else
    throw ReportError("Epetra_Import::Epetra_Import: Bad global indices type", -1);
}

//==============================================================================
// Epetra_Import copy constructor 
Epetra_Import::Epetra_Import(const Epetra_Import & Importer)
  : Epetra_Object(Importer),
    TargetMap_(Importer.TargetMap_),
    SourceMap_(Importer.SourceMap_),
    NumSameIDs_(Importer.NumSameIDs_),
    NumPermuteIDs_(Importer.NumPermuteIDs_),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(Importer.NumRemoteIDs_),
    RemoteLIDs_(0),
    NumExportIDs_(Importer.NumExportIDs_),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(Importer.NumSend_),
    NumRecv_(Importer.NumRecv_),
    Distor_(0)
{
  int i;
  if (NumPermuteIDs_>0) {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
    for (i=0; i< NumPermuteIDs_; i++) {
      PermuteToLIDs_[i] = Importer.PermuteToLIDs_[i];
      PermuteFromLIDs_[i] = Importer.PermuteFromLIDs_[i];
    }
  }

  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    for (i=0; i< NumRemoteIDs_; i++) RemoteLIDs_[i] = Importer.RemoteLIDs_[i];
  }

  if (NumExportIDs_>0) {
    ExportLIDs_ = new int[NumExportIDs_];
    ExportPIDs_ = new int[NumExportIDs_];
    for (i=0; i< NumExportIDs_; i++) {
      ExportLIDs_[i] = Importer.ExportLIDs_[i];
      ExportPIDs_[i] = Importer.ExportPIDs_[i];
    }
  }

  if (Importer.Distor_!=0) Distor_ = Importer.Distor_->Clone();

}

//==============================================================================
// Epetra_Import destructor 
Epetra_Import::~Epetra_Import()
{
  if( Distor_ != 0 ) delete Distor_;
  if (RemoteLIDs_ != 0) delete [] RemoteLIDs_;
  if (PermuteToLIDs_ != 0) delete [] PermuteToLIDs_;
  if (PermuteFromLIDs_ != 0) delete [] PermuteFromLIDs_;

  if( ExportPIDs_ != 0 ) delete [] ExportPIDs_; // These were created by Distor_
  if( ExportLIDs_ != 0 ) delete [] ExportLIDs_;
}

//==============================================================================
// Epetra_Import pseudo-copy constructor. 
Epetra_Import::Epetra_Import(const Epetra_Export& Exporter):
  TargetMap_(Exporter.SourceMap_), //reverse
  SourceMap_(Exporter.TargetMap_),//reverse
  NumSameIDs_(Exporter.NumSameIDs_),
  NumPermuteIDs_(Exporter.NumPermuteIDs_),
  PermuteToLIDs_(0),
  PermuteFromLIDs_(0),
  NumRemoteIDs_(Exporter.NumExportIDs_),//reverse
  RemoteLIDs_(0),
  NumExportIDs_(Exporter.NumRemoteIDs_),//reverse
  ExportLIDs_(0),
  ExportPIDs_(0),
  NumSend_(Exporter.NumRecv_),//reverse
  NumRecv_(Exporter.NumSend_),//revsese
  Distor_(0)
{
  int i;
  // Reverse the permutes
  if (NumPermuteIDs_>0) {
    PermuteToLIDs_   = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
    for (i=0; i< NumPermuteIDs_; i++) {
      PermuteFromLIDs_[i] = Exporter.PermuteToLIDs_[i];
      PermuteToLIDs_[i]   = Exporter.PermuteFromLIDs_[i];
    }
  }

  // Copy the exports to the remotes
  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    for (i=0; i< NumRemoteIDs_; i++) RemoteLIDs_[i] = Exporter.ExportLIDs_[i];
  }

  // Copy the remotes to the exports
  if (NumExportIDs_>0) {
    ExportLIDs_ = new int[NumExportIDs_];
    ExportPIDs_ = new int[NumExportIDs_];
    for (i=0; i< NumExportIDs_; i++) ExportLIDs_[i] = Exporter.RemoteLIDs_[i];

    
    // Extract the RemotePIDs from the Distributor
#ifdef HAVE_MPI
    Epetra_MpiDistributor *D=dynamic_cast<Epetra_MpiDistributor*>(&Exporter.Distributor());
    if(!D) throw ReportError("Epetra_Import: Can't have ExportPIDs w/o an Epetra::MpiDistributor.",-1);
    int i,j,k;
    
    // Get the distributor's data
    int NumReceives        = D->NumReceives();
    const int *ProcsFrom   = D->ProcsFrom();
    const int *LengthsFrom = D->LengthsFrom();
    
    // Now, for each remote ID, record who actually owns it.  This loop follows the operation order in the
    // MpiDistributor so it ought to duplicate that effect.
    for(i=0,j=0;i<NumReceives;i++){
      int pid=ProcsFrom[i];
      for(k=0;k<LengthsFrom[i];k++){
	ExportPIDs_[j]=pid;
	j++;
      }    
    }
#else
    throw ReportError("Epetra_Import: Can't have ExportPIDs w/o an Epetra::MpiDistributor.",-2);
#endif
  }//end NumExportIDs>0

  if (Exporter.Distor_!=0) Distor_ = Exporter.Distor_->ReverseClone();
}

//=============================================================================
void Epetra_Import::Print(ostream & os) const
{
  // mfh 14 Dec 2011: The implementation of Print() I found here
  // previously didn't print much at all, and it included a message
  // saying that it wasn't finished ("Epetra_Import::Print needs
  // attention!!!").  What you see below is a port of
  // Tpetra::Import::print, which does have a full implementation.
  // This should allow a side-by-side comparison of Epetra_Import with
  // Tpetra::Import.

  // If true, then copy the array data and sort it before printing.
  // Otherwise, leave the data in its original order.  
  //
  // NOTE: Do NOT sort the arrays in place!  Only sort in the copy.
  // Epetra depends on the order being preserved, and some arrays'
  // orders are coupled.
  const bool sortIDs = false;

  const Epetra_Comm& comm = SourceMap_.Comm();
  const int myRank = comm.MyPID();
  const int numProcs = comm.NumProc();
  
  if (myRank == 0) {
    os << "Import Data Members:" << endl;
  }
  // We don't need a barrier before this for loop, because Proc 0 is
  // the first one to do anything in the for loop anyway.
  for (int p = 0; p < numProcs; ++p) {
    if (myRank == p) {
      os << "Image ID       : " << myRank << endl;

      os << "permuteFromLIDs:";
      if (PermuteFromLIDs_ == NULL) {
  os << " NULL";
      } else {
  std::vector<int> permuteFromLIDs (NumPermuteIDs_);
  std::copy (PermuteFromLIDs_, PermuteFromLIDs_ + NumPermuteIDs_, 
       permuteFromLIDs.begin());
  if (sortIDs) {
    std::sort (permuteFromLIDs.begin(), permuteFromLIDs.end());
  }
  os << " {";
  for (int i = 0; i < NumPermuteIDs_; ++i) {
    os << permuteFromLIDs[i];
    if (i < NumPermuteIDs_ - 1) {
      os << ", ";
    }
  }
  os << "}";
      }
      os << endl;

      os << "permuteToLIDs  :";
      if (PermuteToLIDs_ == NULL) {
  os << " NULL";
      } else {
  std::vector<int> permuteToLIDs (NumPermuteIDs_);
  std::copy (PermuteToLIDs_, PermuteToLIDs_ + NumPermuteIDs_, 
       permuteToLIDs.begin());
  if (sortIDs) {
    std::sort (permuteToLIDs.begin(), permuteToLIDs.end());
  }
  os << " {";
  for (int i = 0; i < NumPermuteIDs_; ++i) {
    os << permuteToLIDs[i];
    if (i < NumPermuteIDs_ - 1) {
      os << ", ";
    }
  }
  os << "}";
      }
      os << endl;

      os << "remoteLIDs     :";
      if (RemoteLIDs_ == NULL) {
  os << " NULL";
      } else {
  std::vector<int> remoteLIDs (NumRemoteIDs_);
  std::copy (RemoteLIDs_, RemoteLIDs_ + NumRemoteIDs_, 
       remoteLIDs.begin());
  if (sortIDs) {
    std::sort (remoteLIDs.begin(), remoteLIDs.end());
  }
  os << " {";
  for (int i = 0; i < NumRemoteIDs_; ++i) {
    os << remoteLIDs[i];
    if (i < NumRemoteIDs_ - 1) {
      os << ", ";
    }
  }
  os << "}";
      }
      os << endl;

      // If sorting for output, the export LIDs and export PIDs have
      // to be sorted together.  We can use Epetra_Util::Sort, using
      // the PIDs as the keys to match Tpetra::Import.
      std::vector<int> exportLIDs (NumExportIDs_);
      std::vector<int> exportPIDs (NumExportIDs_);
      if (ExportLIDs_ != NULL) {
  std::copy (ExportLIDs_, ExportLIDs_ + NumExportIDs_, exportLIDs.begin());
  std::copy (ExportPIDs_, ExportPIDs_ + NumExportIDs_, exportPIDs.begin());

  if (sortIDs && NumExportIDs_ > 0) {
    int* intCompanions[1]; // Input for Epetra_Util::Sort().
    intCompanions[0] = &exportLIDs[0];
    Epetra_Util::Sort (true, NumExportIDs_, &exportPIDs[0], 
           0, (double**) NULL, 1, intCompanions, 0, 0);
  }
      }

      os << "exportLIDs     :";
      if (ExportLIDs_ == NULL) {
  os << " NULL";
      } else {
  os << " {";
  for (int i = 0; i < NumExportIDs_; ++i) {
    os << exportLIDs[i];
    if (i < NumExportIDs_ - 1) {
      os << ", ";
    }
  }
  os << "}";
      }
      os << endl;

      os << "exportImageIDs :";
      if (ExportPIDs_ == NULL) {
  os << " NULL";
      } else {
  os << " {";
  for (int i = 0; i < NumExportIDs_; ++i) {
    os << exportPIDs[i];
    if (i < NumExportIDs_ - 1) {
      os << ", ";
    }
  }
  os << "}";
      }
      os << endl;

      os << "numSameIDs     : " << NumSameIDs_ << endl;
      os << "numPermuteIDs  : " << NumPermuteIDs_ << endl;
      os << "numRemoteIDs   : " << NumRemoteIDs_ << endl;
      os << "numExportIDs   : " << NumExportIDs_ << endl;

      // Epetra keeps NumSend_ and NumRecv_, whereas in Tpetra, these
      // are stored in the Distributor object.  This is why we print
      // them here.
      os << "Number of sends: " << NumSend_ << endl;
      os << "Number of recvs: " << NumRecv_ << endl;
    } // if my rank is p

    // A few global barriers give I/O a chance to complete.
    comm.Barrier();
    comm.Barrier();
    comm.Barrier();
  } // for each rank p

  const bool printMaps = false;
  if (printMaps) {
    // The original implementation printed the Maps first.  We moved
    // printing the Maps to the end, for easy comparison with the
    // output of Tpetra::Import::print().
    if (myRank == 0) {
      os << endl << endl << "Source Map:" << endl << std::flush;
    }
    comm.Barrier();
    SourceMap_.Print(os);
    comm.Barrier();
  
    if (myRank == 0) {
      os << endl << endl << "Target Map:" << endl << std::flush;
    }
    comm.Barrier();
    TargetMap_.Print(os);
    comm.Barrier();
  }

  if (myRank == 0) {
    os << endl << endl << "Distributor:" << endl << std::flush;
  }
  comm.Barrier();
  if (Distor_ == NULL) {
    if (myRank == 0) {
      os << " is NULL." << endl;
    }
  } else {
    Distor_->Print(os); // Printing the Distributor is itself distributed.
  }
  comm.Barrier();
}

