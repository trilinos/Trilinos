
#include "Petra_Import.h"

#include "Petra_BlockMap.h"
#include "Petra_Map.h"

#ifdef PETRA_MPI
#include "GSComm_Plan.h"
#endif

//==============================================================================
// Petra_Import constructor for a Petra_BlockMap object
Petra_Import::Petra_Import( const Petra_BlockMap &  TargetMap, const Petra_BlockMap & SourceMap)
  : TargetMap_(TargetMap),
    SourceMap_(SourceMap),
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
#ifdef PETRA_MPI
    NumRecv_(0),
    GSPlan_(0)
#else
    NumRecv_(0)
#endif
{

  int i;
  
  // Build three ID lists:
  // NumSameIDs - Number of IDs in TargetMap and SourceMap that are identical, up to the first
  //              nonidentical ID.
  // NumPermuteIDs - Number of IDs in SourceMap that must be indirectly loaded but are on this processor.
  // NumRemoteIDs - Number of IDs that are in SourceMap but not in TargetMap, and thus must be imported.
  
  int NumSourceIDs = SourceMap.NumMyElements();
  int NumTargetIDs = TargetMap.NumMyElements();
  
  int *TargetGIDs = 0;
  if (NumTargetIDs>0) {
    TargetGIDs = new int[NumTargetIDs];
    TargetMap.MyGlobalElements(TargetGIDs);
  }
  
  int * SourceGIDs = 0;
  if (NumSourceIDs>0) {
    SourceGIDs = new int[NumSourceIDs];
    SourceMap.MyGlobalElements(SourceGIDs);
  }
  
  int MinIDs = minfn(NumSourceIDs, NumTargetIDs);
  
  
  NumSameIDs_ = 0;
  for (i=0; i< MinIDs; i++) if (TargetGIDs[i]==SourceGIDs[i]) NumSameIDs_++; else break;
  
  
  // Find count of Target IDs that are truly remote and those that are local but permuted

  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) 
    if (SourceMap.MyGID(TargetGIDs[i])) NumPermuteIDs_++; // Check if Target GID is a local Source GID
    else NumRemoteIDs_++; // If not, then it is remote
  
  
  
  // Define remote and permutation lists
  
  int * RemoteGIDs;
  if (NumRemoteIDs_>0) {
    RemoteLIDs_ = new int[NumRemoteIDs_];
    RemoteGIDs = new int[NumRemoteIDs_];
  }
  if (NumPermuteIDs_>0)  {
    PermuteToLIDs_ = new int[NumPermuteIDs_];
    PermuteFromLIDs_ = new int[NumPermuteIDs_];
  }
  
  NumPermuteIDs_ = 0;
  NumRemoteIDs_ = 0;
  for (i=NumSameIDs_; i< NumTargetIDs; i++) {
    if (SourceMap.MyGID(TargetGIDs[i])) {
      PermuteToLIDs_[NumPermuteIDs_] = i;
      PermuteFromLIDs_[NumPermuteIDs_++] = SourceMap.LID(TargetGIDs[i]);
    }
    else {
      //NumRecv_ +=TargetMap.ElementSize(i); // Count total number of entries to receive
      NumRecv_ +=TargetMap.MaxElementSize(); // Count total number of entries to receive (currently need max)
      RemoteGIDs[NumRemoteIDs_] = TargetGIDs[i];
      RemoteLIDs_[NumRemoteIDs_++] = i;
    }
  }
  
  if ( NumRemoteIDs_>0 && !SourceMap.DistributedGlobal()) {
    cout << "Error in Petra_Import: Serial Import has remote IDs." << endl; 
    abort();
  }
  
#ifdef PETRA_MPI
  // Test for distributed cases
  
  if (SourceMap.DistributedGlobal()) {
    
    int * RemotePIDs = 0;
    
    if (NumRemoteIDs_>0)  RemotePIDs = new int[NumRemoteIDs_];
    SourceMap.RemoteIDList(NumRemoteIDs_, RemoteGIDs, RemotePIDs, 0); // Get remote PIDs
    
    GSPlan_ = new GSComm_Plan();
    
    // Construct list of exports that calling processor needs to send as a result
    // of everyone asking for what it needs to receive.
    
    int msgtag = 32765; // Note: Should get message tags from Petra_Comm?
    bool Deterministic = true;
    
    bool comm_flag = GSPlan_->CreateFromRecvs( NumRemoteIDs_, RemoteGIDs,
					       RemotePIDs, SourceMap.Comm().Comm(), msgtag, Deterministic,
					       NumExportIDs_, ExportLIDs_, ExportPIDs_ );
    
    GSComm_Comm * GSComm = new GSComm_Comm;
    
    // Use comm plan with Export GIDs (stored in ExportLIDs_) to
    // get proper ordering of GIDs for remote entries 
    // (that we will convert to LIDs when done).
    
    GSComm->Do( *GSPlan_, msgtag, 
		reinterpret_cast<char *> (ExportLIDs_), 
		sizeof( int ),
		reinterpret_cast<char *> (RemoteGIDs) );
    
    delete GSComm;
    
    // Export IDs come in as GIDs, convert to LIDs
    for (i=0; i< NumExportIDs_; i++) {
      if (ExportPIDs_[i] < 0) { // Check for valid PID.  Abort if any bad PIDs.
	cout << "TargetMap requested a GID that is not in the SourceMap.  Must abort." << endl;
	abort();
      }
      
      ExportLIDs_[i] = SourceMap.LID(ExportLIDs_[i]);
      //NumSend_ += SourceMap.ElementSize(ExportLIDs_[i]); // Count total number of entries to send
      NumSend_ += SourceMap.MaxElementSize(); // Count total number of entries to send (currently need max)
    }
    
    // Remote IDs come in as GIDs, convert to LIDs in proper order

    // for (i=0; i< NumRemoteIDs_; i++) RemoteLIDs_[i] = TargetMap.LID(RemoteGIDs[i]); // Only works when target map has no repeated GIDs
    
    if (NumRemoteIDs_>0) {
      int * ReorderedRemoteLIDs = RemotePIDs; // Reuse some temp space
      for (i=0; i< NumRemoteIDs_; i++) {
	int CurrentGID = RemoteGIDs[i];
	bool Found = false;
	for (int j=0; j < NumRemoteIDs_; j++) {
	  if (RemoteLIDs_[j]!= -1) {
	    if (CurrentGID==TargetGIDs[RemoteLIDs_[j]]) {
	      ReorderedRemoteLIDs[i] = RemoteLIDs_[j];
	      RemoteLIDs_[j] = -1;
	      Found = true;
	      break;
	    }
	  }
	}
	if (!Found) {
	  cout << "Internal error.  Cannot map incoming GID to Target Map.  Must abort." << endl;
	  abort();
	}
      }
      
      // Clean up and leave....
      
      delete [] RemoteLIDs_;
      delete [] RemoteGIDs;
      RemoteLIDs_ = ReorderedRemoteLIDs;
    }
  }
#endif

  if (NumTargetIDs>0) delete [] TargetGIDs;
  if (NumSourceIDs>0) delete [] SourceGIDs;
  
  return;
}

//==============================================================================
// Petra_Import copy constructor 
Petra_Import::Petra_Import(const Petra_Import & Importer)
  : TargetMap_(Importer.TargetMap_),
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
#ifdef PETRA_MPI
    NumRecv_(Importer.NumRecv_),
    GSPlan_(0)
#else
    NumRecv_(Importer.NumRecv_)
#endif
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

#ifdef PETRA_MPI
  if (Importer.GSPlan_!=0) GSPlan_ = new GSComm_Plan(*Importer.GSPlan_);
#endif

}

//==============================================================================
// Petra_Import destructor 
Petra_Import::~Petra_Import()
{
#ifdef PETRA_MPI
  if( GSPlan_ != 0 ) delete GSPlan_;
#endif

  if (RemoteLIDs_ != 0) delete [] RemoteLIDs_;
  if (PermuteToLIDs_ != 0) delete [] PermuteToLIDs_;
  if (PermuteFromLIDs_ != 0) delete [] PermuteFromLIDs_;

  if( ExportPIDs_ != 0 ) delete [] ExportPIDs_; // These were created by GSPlan
  if( ExportLIDs_ != 0 ) delete [] ExportLIDs_;
}

//==============================================================================
// operator<<
/*
ostream & operator<<( ostream & os, const Petra_Import & pd )
{
 
  int MyPID;
  if( pd.DirectoryMap_ != 0 )
  {
    MyPID = pd.DirectoryMap_->Comm().MyPID();
    os << MyPID << " Petra_Import Object: "
      << pd.DirectoryMap_->NumMyElements() << endl;
    for( int i = 0; i < pd.DirectoryMap_->NumMyElements(); i++ )
      os << " " << i << " " << pd.ProcList_[i] << " "
        << pd.LocalIndexList_[i] << endl;
     os << endl;
//     os << "Directory " << *(pd.GSPlan_) << endl;
  }
  else
  {
    cout << "Petra_Import not setup<<<<<<" << endl;
  }

  return os;
}
*/
