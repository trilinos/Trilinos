
#include "Petra_Export.h"

#include "Petra_BlockMap.h"
#include "Petra_Map.h"

#ifdef PETRA_MPI
#include "GSComm_Plan.h"
#endif

//==============================================================================
// Petra_Export constructor for a Petra_BlockMap object
Petra_Export::Petra_Export( const Petra_BlockMap &  SourceMap, const Petra_BlockMap & TargetMap)
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
  // NumExportIDs - Number of IDs that are in SourceMap but not in TargetMap, and thus must be exported.

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


  // Find count of Source IDs that are truly remote and those that are local but permuted

  NumPermuteIDs_ = 0;
  NumExportIDs_ = 0;
  for (i=NumSameIDs_; i< NumSourceIDs; i++) 
    if (TargetMap.MyGID(SourceGIDs[i])) NumPermuteIDs_++; // Check if Source GID is a local Target GID
    else NumExportIDs_++; // If not, then it is remote



  // Define remote and permutation lists

  int * ExportGIDs;
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
    if (TargetMap.MyGID(SourceGIDs[i])) {
      PermuteFromLIDs_[NumPermuteIDs_] = i;
      PermuteToLIDs_[NumPermuteIDs_++] = TargetMap.LID(SourceGIDs[i]);
    }
    else {
      //NumSend_ +=SourceMap.ElementSize(i); // Count total number of entries to send
      NumSend_ +=SourceMap.MaxElementSize(); // Count total number of entries to send (currently need max)
      ExportGIDs[NumExportIDs_] = SourceGIDs[i];
      ExportLIDs_[NumExportIDs_++] = i;
    }
  }
     
  if ( NumExportIDs_>0 && !SourceMap.DistributedGlobal()) {
    cout << "Error in Petra_Export: Serial Export has remote IDs." << endl; 
    abort();
  }

#ifdef PETRA_MPI
  // Test for distributed cases

  if (SourceMap.DistributedGlobal()) {

    ExportPIDs_ = 0;
    
    if (NumExportIDs_>0) {
      ExportPIDs_ = new int[NumExportIDs_];
      TargetMap.RemoteIDList(NumExportIDs_, ExportGIDs, ExportPIDs_, 0); // Get remote PIDs
      for (i=0; i< NumExportIDs_; i++) { 
	if (ExportPIDs_[i] < 0) { // Check for valid PID.  Abort if any bad PIDs.
      cout << "SourceMap requested a GID that is not in the TargetMap.  Must abort." << endl;
      abort();
	}
      }
      
    }
    
    GSPlan_ = new GSComm_Plan();
    
    // Construct list of exports that calling processor needs to send as a result
    // of everyone asking for what it needs to receive.
    
    int msgtag = 32765; // Note: Should get message tags from Petra_Comm?
    bool Deterministic = true;
    
    bool comm_flag = GSPlan_->CreateFromSends( NumExportIDs_, ExportPIDs_, 
					       SourceMap.Comm().Comm(), msgtag, Deterministic,
					       NumRemoteIDs_);
    
    GSComm_Comm * GSComm = new GSComm_Comm;
    
    // Use comm plan with ExportGIDs to find out who is sending to us and
    // get proper ordering of GIDs for remote entries 
    // (that we will convert to LIDs when done).
    
    if (NumRemoteIDs_>0) RemoteLIDs_ = new int[NumRemoteIDs_]; // Allocate space for LIDs in target that are
    // going to get something from off-processor.
    GSComm->Do( *GSPlan_, msgtag, 
		reinterpret_cast<char *> (ExportGIDs), 
		sizeof( int ),
		reinterpret_cast<char *> (RemoteLIDs_) );
    
    delete GSComm;
    
  // Remote IDs come in as GIDs, convert to LIDs
  for (i=0; i< NumRemoteIDs_; i++) {
    RemoteLIDs_[i] = TargetMap.LID(RemoteLIDs_[i]);
    //NumRecv_ += TargetMap.ElementSize(RemoteLIDs_[i]); // Count total number of entries to receive
    NumRecv_ += TargetMap.MaxElementSize(); // Count total number of entries to receive (currently need max)
  }

    if (NumExportIDs_>0) delete [] ExportGIDs;
  }
#endif
  
  if (NumTargetIDs>0) delete [] TargetGIDs;
  if (NumSourceIDs>0) delete [] SourceGIDs;
  
  return;
}

//==============================================================================
// Petra_Export copy constructor 
Petra_Export::Petra_Export(const Petra_Export & Exporter)
  : TargetMap_(Exporter.TargetMap_),
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
#ifdef PETRA_MPI
    NumRecv_(Exporter.NumRecv_),
    GSPlan_(0)
#else
    NumRecv_(Exporter.NumRecv_)
#endif
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

#ifdef PETRA_MPI
  if (Exporter.GSPlan_!=0) GSPlan_ = new GSComm_Plan(*Exporter.GSPlan_);
#endif

}

//==============================================================================
// Petra_Export destructor 
Petra_Export::~Petra_Export()
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
ostream & operator<<( ostream & os, const Petra_Export & pd )
{
 
  int MyPID;
  if( pd.DirectoryMap_ != 0 )
  {
    MyPID = pd.DirectoryMap_->Comm().MyPID();
    os << MyPID << " Petra_Export Object: "
      << pd.DirectoryMap_->NumMyElements() << endl;
    for( int i = 0; i < pd.DirectoryMap_->NumMyElements(); i++ )
      os << " " << i << " " << pd.ProcList_[i] << " "
        << pd.LocalIndexList_[i] << endl;
     os << endl;
//     os << "Directory " << *(pd.GSPlan_) << endl;
  }
  else
  {
    cout << "Petra_Export not setup<<<<<<" << endl;
  }

  return os;
}
*/
