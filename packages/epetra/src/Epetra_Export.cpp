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

#include "Epetra_Export.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Distributor.h"
#include "Epetra_Comm.h"
#include "Epetra_Util.h"

//==============================================================================
// Epetra_Export constructor for a Epetra_BlockMap object
Epetra_Export::Epetra_Export( const Epetra_BlockMap &  sourceMap, const Epetra_BlockMap & targetMap)
  : Epetra_Object("Epetra::Export"), 
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

  int i;

  // Build three ID lists:
  // NumSameIDs - Number of IDs in TargetMap and SourceMap that are identical, up to the first
  //              nonidentical ID.
  // NumPermuteIDs - Number of IDs in SourceMap that must be indirectly loaded but are on this processor.
  // NumExportIDs - Number of IDs that are in SourceMap but not in TargetMap, and thus must be exported.

  int NumSourceIDs = sourceMap.NumMyElements();
  int NumTargetIDs = targetMap.NumMyElements();

  int *TargetGIDs = 0;
  if (NumTargetIDs>0) {
    TargetGIDs = new int[NumTargetIDs];
    targetMap.MyGlobalElements(TargetGIDs);
  }

  int * SourceGIDs = 0;
  if (NumSourceIDs>0) {
    SourceGIDs = new int[NumSourceIDs];
    sourceMap.MyGlobalElements(SourceGIDs);
  }

  int MinIDs = EPETRA_MIN(NumSourceIDs, NumTargetIDs);

  NumSameIDs_ = 0;
  for (i=0; i< MinIDs; i++) if (TargetGIDs[i]==SourceGIDs[i]) NumSameIDs_++; else break;

  // Find count of Source IDs that are truly remote and those that are local but permuted

  NumPermuteIDs_ = 0;
  NumExportIDs_ = 0;
  for (i=NumSameIDs_; i< NumSourceIDs; i++) 
    if (targetMap.MyGID(SourceGIDs[i])) NumPermuteIDs_++; // Check if Source GID is a local Target GID
    else NumExportIDs_++; // If not, then it is remote

  // Define remote and permutation lists

  int * ExportGIDs = 0;
  if (NumExportIDs_>0) {
    ExportLIDs_ = new int[NumExportIDs_];
    ExportGIDs = new int[NumExportIDs_];
  }
  if (NumPermuteIDs_>0)  {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
  }

  NumPermuteIDs_ = 0;
  NumExportIDs_ = 0;
  for (i=NumSameIDs_; i< NumSourceIDs; i++) {
    if (targetMap.MyGID(SourceGIDs[i])) {
      PermuteFromLIDs_[NumPermuteIDs_] = i;
      PermuteToLIDs_[NumPermuteIDs_++] = targetMap.LID(SourceGIDs[i]);
    }
    else {
      //NumSend_ +=sourceMap.ElementSize(i); // Count total number of entries to send
      NumSend_ +=sourceMap.MaxElementSize(); // Count total number of entries to send (currently need max)
      ExportGIDs[NumExportIDs_] = SourceGIDs[i];
      ExportLIDs_[NumExportIDs_++] = i;
    }
  }
     
  if ( NumExportIDs_>0 && !sourceMap.DistributedGlobal()) 
    ReportError("Warning in Epetra_Export: Serial Export has remote IDs. (Exporting from Subset of Source Map)", 1);

  // Test for distributed cases
  int ierr = 0;

  if (sourceMap.DistributedGlobal()) {

    if (NumExportIDs_>0) ExportPIDs_ = new int[NumExportIDs_];
    ierr = targetMap.RemoteIDList(NumExportIDs_, ExportGIDs, ExportPIDs_, 0); // Get remote PIDs
    if( ierr ) throw ReportError("Error in Epetra_BlockMap::RemoteIDList", ierr);

    //Get rid of IDs not in Target Map
    if(NumExportIDs_>0) {
      int cnt = 0;
      for( i = 0; i < NumExportIDs_; ++i )
	if( ExportPIDs_[i] == -1 ) ++cnt;
      if( cnt ) {
	int * NewExportGIDs = 0;
	int * NewExportPIDs = 0;
	int * NewExportLIDs = 0;
	int cnt1 = NumExportIDs_-cnt;
	if (cnt1) {
	  NewExportGIDs = new int[cnt1];
	  NewExportPIDs = new int[cnt1];
	  NewExportLIDs = new int[cnt1];
	}
	cnt = 0;
	for( i = 0; i < NumExportIDs_; ++i )
	  if( ExportPIDs_[i] != -1 ) {
	    NewExportGIDs[cnt] = ExportGIDs[i];
	    NewExportPIDs[cnt] = ExportPIDs_[i];
	    NewExportLIDs[cnt] = ExportLIDs_[i];
	    ++cnt;
          }
	assert(cnt==cnt1); // Sanity test
	NumExportIDs_ = cnt;
	delete [] ExportGIDs;
	delete [] ExportPIDs_;
	delete [] ExportLIDs_;
	ExportGIDs = NewExportGIDs;
	ExportPIDs_ = NewExportPIDs;
	ExportLIDs_ = NewExportLIDs;
	ReportError("Warning in Epetra_Export: Source IDs not found in Target Map (Do you want to export from subset of Source Map?)", 1 );
      }
    }
    
    //Make sure Export IDs are ordered by processor
    Epetra_Util util;
    int * tmpPtr[2];
    tmpPtr[0] = ExportLIDs_, tmpPtr[1] = ExportGIDs;
    util.Sort(true,NumExportIDs_,ExportPIDs_,0,0,2,tmpPtr);

    Distor_ = sourceMap.Comm().CreateDistributor();
    
    // Construct list of exports that calling processor needs to send as a result
    // of everyone asking for what it needs to receive.
    
    ierr = Distor_->CreateFromSends( NumExportIDs_, ExportPIDs_, true, NumRemoteIDs_);
    if (ierr!=0) throw ReportError("Error in Epetra_Distributor.CreateFromSends()", ierr);
    
    // Use comm plan with ExportGIDs to find out who is sending to us and
    // get proper ordering of GIDs for remote entries 
    // (that we will convert to LIDs when done).
    
    if (NumRemoteIDs_>0) RemoteLIDs_ = new int[NumRemoteIDs_]; // Allocate space for LIDs in target that are
    // going to get something from off-processor.
    char * cRemoteGIDs = 0; //Do will alloc memory for this object
    int LenCRemoteGIDs = 0;
    ierr = Distor_->Do(reinterpret_cast<char *> (ExportGIDs), 
		sizeof( int ),
		LenCRemoteGIDs,
		cRemoteGIDs);
    if (ierr) throw ReportError("Error in Epetra_Distributor.Do()", ierr);
    int * RemoteGIDs = reinterpret_cast<int*>(cRemoteGIDs);

    // Remote IDs come in as GIDs, convert to LIDs
    for (i=0; i< NumRemoteIDs_; i++) {
      RemoteLIDs_[i] = targetMap.LID(RemoteGIDs[i]);
      //NumRecv_ += targetMap.ElementSize(RemoteLIDs_[i]); // Count total number of entries to receive
      NumRecv_ += targetMap.MaxElementSize(); // Count total number of entries to receive (currently need max)
    }

    if (LenCRemoteGIDs>0) delete [] cRemoteGIDs;
  }
  if (NumExportIDs_>0) delete [] ExportGIDs;
  if (NumTargetIDs>0) delete [] TargetGIDs;
  if (NumSourceIDs>0) delete [] SourceGIDs;
  
  return;
}

//==============================================================================
// Epetra_Export copy constructor 
Epetra_Export::Epetra_Export(const Epetra_Export & Exporter)
  : Epetra_Object(Exporter), 
     TargetMap_(Exporter.TargetMap_),
    SourceMap_(Exporter.SourceMap_),
    NumSameIDs_(Exporter.NumSameIDs_),
    NumPermuteIDs_(Exporter.NumPermuteIDs_),
    PermuteToLIDs_(0),
    PermuteFromLIDs_(0),
    NumRemoteIDs_(Exporter.NumRemoteIDs_),
    RemoteLIDs_(0),
    NumExportIDs_(Exporter.NumExportIDs_),
    ExportLIDs_(0),
    ExportPIDs_(0),
    NumSend_(Exporter.NumSend_),
    NumRecv_(Exporter.NumRecv_),
    Distor_(0)
{
  int i;
  if (NumPermuteIDs_>0) {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
    for (i=0; i< NumPermuteIDs_; i++) {
      PermuteToLIDs_[i] = Exporter.PermuteToLIDs_[i];
      PermuteFromLIDs_[i] = Exporter.PermuteFromLIDs_[i];
    }
  }

  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    for (i=0; i< NumRemoteIDs_; i++) RemoteLIDs_[i] = Exporter.RemoteLIDs_[i];
  }

  TargetMap().Comm().Barrier();
  if (NumExportIDs_>0) {
    ExportLIDs_ = new int[NumExportIDs_];
    ExportPIDs_ = new int[NumExportIDs_];
    for (i=0; i< NumExportIDs_; i++) {
      ExportLIDs_[i] = Exporter.ExportLIDs_[i];
      ExportPIDs_[i] = Exporter.ExportPIDs_[i];
    }
  }

  if (Exporter.Distor_!=0) Distor_ = Exporter.Distor_->Clone();

}

//==============================================================================
// Epetra_Export destructor 
Epetra_Export::~Epetra_Export()
{
  if( Distor_ != 0 ) delete Distor_;

  if (RemoteLIDs_ != 0) delete [] RemoteLIDs_;
  if (PermuteToLIDs_ != 0) delete [] PermuteToLIDs_;
  if (PermuteFromLIDs_ != 0) delete [] PermuteFromLIDs_;

  if( ExportPIDs_ != 0 ) delete [] ExportPIDs_; // These were created by GSPlan
  if( ExportLIDs_ != 0 ) delete [] ExportLIDs_;
}
//=============================================================================
void Epetra_Export::Print(ostream & os) const
{

  os << endl << endl << "Source Map:" << endl << endl;
  SourceMap_.Print(os);
  
  os << endl << endl << "Target Map:" << endl << endl;
  TargetMap_.Print(os);
  
  os << endl << endl << "Distributor:" << endl << endl;
  if (Distor_==0) os << "  Is empty...." << endl;
  else Distor_->Print(os);
  
  os << "Number of Same IDs = " << NumSameIDs_ << endl;
  os << "Number of Permute IDs = " << NumPermuteIDs_ << endl;
  os << "Number of Export IDs = " << NumExportIDs_ << endl;
  
  os << "Epetra_Export Print Needs attention!!!!" << endl;
  return;
}

//----------------------------------------------------------------------------
Epetra_Export& Epetra_Export::operator=(const Epetra_Export& src)
{
  (void)src;
  //not currently supported
  bool throw_err = true;
  if (throw_err) {
    throw ReportError("Epetra_Export::operator= not supported.",-1);
  }
  return(*this);
}
