
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_BasicDirectory.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"


//==============================================================================
// Epetra_BasicDirectory constructor for a Epetra_BlockMap object
Epetra_BasicDirectory::Epetra_BasicDirectory(const Epetra_BlockMap & Map)
  : Map_(&Map),
    DirectoryMap_(0),
    ProcList_(0),
    LocalIndexList_(0),
    SizeList_(0),
    AllMinGIDs_(0)
{

  // Test for simple cases

  // Uniprocessor and local map cases (Nothing to set up)

  if (!(Map_->DistributedGlobal())) return;

  // Linear Map case

  else if (Map_->LinearMap()) {

    // Build a list of the Minimum global ids for all processors on each processor.
    // Since the map is linear, we know that all GIDs are contiguous on each processor
    // and can be found using the MinGIDs.

    int NumProc = Map_->Comm().NumProc();
    AllMinGIDs_ = new int[NumProc+1];
    int MinMyGID = Map_->MinMyGID();
    Map_->Comm().GatherAll(&MinMyGID, AllMinGIDs_, 1);
    AllMinGIDs_[NumProc] = 1 + Map_->MaxAllGID(); // Set max cap
  }

  // General case.  Need to build a directory via calls to communication functions
  else {
    
    int flag = Generate();
    assert(flag==0);
  }
}
//==============================================================================
// Epetra_BasicDirectory copy constructor
Epetra_BasicDirectory::Epetra_BasicDirectory(const Epetra_BasicDirectory & Directory)
  : Map_(Directory.Map_),
    DirectoryMap_(0),
    ProcList_(0),
    LocalIndexList_(0),
    SizeList_(0),
    AllMinGIDs_(0)
{
  int i;

  if (Directory.DirectoryMap_!=0) DirectoryMap_ = new Epetra_Map(Directory.DirectoryMap());

  int Dir_NumMyElements = DirectoryMap_->NumMyElements();

  if (Directory.ProcList_!=0) {
    ProcList_ = new int[Dir_NumMyElements];
    for (i=0; i<Dir_NumMyElements; i++) ProcList_[i] = Directory.ProcList_[i];
  }
  if (Directory.LocalIndexList_!=0) {
    LocalIndexList_ = new int[Dir_NumMyElements];
    for (int i=0; i<Dir_NumMyElements; i++) LocalIndexList_[i] = Directory.LocalIndexList_[i];
    }
  if (Directory.SizeList_!=0) {
    SizeList_ = new int[Dir_NumMyElements];
    for (int i=0; i<Dir_NumMyElements; i++) SizeList_[i] = Directory.SizeList_[i];
    }
  if (Directory.AllMinGIDs_!=0) {
    int NumProc = Map_->Comm().NumProc();
    AllMinGIDs_ = new int[NumProc+1];
    for (int i=0; i<NumProc+1; i++) AllMinGIDs_[i] = Directory.AllMinGIDs_[i];
    }

}
//==============================================================================
// Epetra_BasicDirectory destructor 
Epetra_BasicDirectory::~Epetra_BasicDirectory()
{
  if( DirectoryMap_ != 0 ) delete DirectoryMap_;
  if( ProcList_ != 0 ) delete [] ProcList_;
  if( LocalIndexList_ != 0 ) delete [] LocalIndexList_;
  if( SizeList_ != 0 ) delete [] SizeList_;
  if( AllMinGIDs_ != 0 ) delete [] AllMinGIDs_;

  DirectoryMap_ = 0;
  ProcList_ = 0 ;
  LocalIndexList_ = 0;
  SizeList_ = 0;
  AllMinGIDs_ = 0;
}

//==============================================================================
// Generate: Generates Directory Tables
int Epetra_BasicDirectory::Generate()
{
  int i;
  bool SizeIsConst = Map_->ConstantElementSize();
  int MinAllGID = Map_->MinAllGID();
  int MaxAllGID = Map_->MaxAllGID();
  // DirectoryMap will have a range of elements from the minimum to the maximum
  // GID of the user map, and an IndexBase of MinAllGID from the user map
  int Dir_NumGlobalElements = MaxAllGID - MinAllGID + 1;

  // Create a uniform linear map to contain the directory
  DirectoryMap_ = new Epetra_Map( Dir_NumGlobalElements, MinAllGID, Map_->Comm() );

  int Dir_NumMyElements = DirectoryMap_->NumMyElements(); // Get NumMyElements



  // Allocate Processor list and Local Index List.  Initialize to -1s.

  if (Dir_NumMyElements>0) {
    ProcList_ = new int[ Dir_NumMyElements ];
    LocalIndexList_ = new int[ Dir_NumMyElements ];
    if (!SizeIsConst) SizeList_ = new int[ Dir_NumMyElements ];
    // Initialize values to -1 in case the user global element list does
    // fill all IDs from MinAllGID to MaxAllGID (e.g., allows global indices to be 
    // all even integers.
    for (i=0; i<Dir_NumMyElements; i++) {
      ProcList_[i] = -1;
      LocalIndexList_[i] = -1;
      if (!SizeIsConst) SizeList_[i] = -1;
    }
  }

  
  // Get list of processors owning the directory entries for the Map GIDs

  int NumProc = Map_->Comm().NumProc();
  int MyPID = Map_->Comm().MyPID();

  int Map_NumMyElements = Map_->NumMyElements();
  int * send_procs = 0;
  if (Map_NumMyElements>0) send_procs = new int[Map_NumMyElements];
  int * Map_MyGlobalElements = Map_->MyGlobalElements();

  assert(DirectoryMap_->RemoteIDList(Map_NumMyElements, Map_MyGlobalElements, 
				     send_procs, 0)==0); 



  bool det_flag = true;

  int num_recvs=0;
    
  Epetra_Distributor * Distor = Map_->Comm().CreateDistributor();

  EPETRA_CHK_ERR(Distor->CreateFromSends( Map_NumMyElements, send_procs, det_flag, num_recvs ));

  if (Map_NumMyElements>0) delete [] send_procs;

  int * export_elements = 0;
  int * import_elements = 0;
  int * ElementSizeList = 0;
  int packetSize = 3; // Assume we will send GIDs, PIDs and LIDs (will increase to 4 if also sending sizes)
  if (!SizeIsConst) packetSize++; // Must send element size info also
 
  if (Map_NumMyElements>0) {
    if (!SizeIsConst) ElementSizeList = Map_->ElementSizeList();
    export_elements = new int[ packetSize * Map_NumMyElements ];
    int * ptr = export_elements;
    for( i = 0; i < Map_NumMyElements; i++ )
      {
	*ptr++ = Map_MyGlobalElements[i];
	*ptr++ = MyPID;
	*ptr++ = i;
	if (!SizeIsConst) *ptr++ = ElementSizeList[i];
      }
  }

  if (num_recvs>0) import_elements = new int[ packetSize * num_recvs ];
  //for (i=0; i< packetSize*num_recvs; i++) import_elements[i] = 0;

  EPETRA_CHK_ERR(Distor->Do(reinterpret_cast<char *> (export_elements), 
		  packetSize * sizeof( int ),
		  reinterpret_cast<char *> (import_elements) ));

  
  //bool MYPID = (Map_->Comm().MyPID()==0);
  int curr_LID;
  //if (MYPID) cout << "Processor " << Map_->Comm().MyPID()<< "  num_recvs = "<< num_recvs << endl << flush;
  int * ptr = import_elements;
  for( i = 0; i < num_recvs; i++ )
  {
    curr_LID = DirectoryMap_->LID(*ptr++); // Convert incoming GID to Directory LID
    //if (MYPID) cout << " Receive ID = " << i << "  GID = " << import_elements[3*i] << "  LID = " << curr_LID << endl << flush;
    assert(curr_LID !=-1); // Internal error
    ProcList_[ curr_LID ] = *ptr++;
    LocalIndexList_[ curr_LID ] = *ptr++;
    if (!SizeIsConst) SizeList_[ curr_LID ] = *ptr++;
  }

  if (import_elements!=0) delete [] import_elements;
  if (export_elements!=0) delete [] export_elements;
  
  delete Distor;
  return(0);
}
//==============================================================================
// GetDirectoryEntries: Get non-local GID references ( procID and localID )
// 			Space should already be allocated for Procs and
//     			LocalEntries.
int Epetra_BasicDirectory::GetDirectoryEntries( const int NumEntries,
					  const int * GlobalEntries, int * Procs,
					  int * LocalEntries, int * EntrySizes ) const {



  int ierr = 0;
  int j;
  int i;
  int MyPID = Map_->Comm().MyPID();
  int NumProc = Map_->Comm().NumProc();
  int n_over_p = Map_->NumGlobalElements() / NumProc;
  int remainder = Map_->NumGlobalElements() % NumProc;


  // Test for simple cases

  // Uniprocessor and local map cases

  if (!Map_->DistributedGlobal()) {
    for (i=0; i<NumEntries; i++) {
      if (LocalEntries!=0) LocalEntries[i] = Map_->LID(GlobalEntries[i]);

      // If GID is not valid, return -1 for both the proc and local entry info
      if (LocalEntries!=0)  if (LocalEntries[i]==-1) Procs[i] = -1; 
      else Procs[i] = MyPID;

    }
    if (EntrySizes!=0) {
      if (Map_->ConstantElementSize()) {
	int ElementSize = Map_->MaxElementSize();
	for (i=0; i<NumEntries; i++) EntrySizes[i] = ElementSize;
      }
      else {
	int * ElementSizeList = Map_->ElementSizeList();
	for (i=0; i<NumEntries; i++) if (LocalEntries[i]>-1) EntrySizes[i] = ElementSizeList[LocalEntries[i]];
      }
    }
    EPETRA_CHK_ERR(ierr);
    return(0);
  }

  // Linear Map case
  if (Map_->LinearMap()) {
    
    int MinAllGID = Map_->MinAllGID(); // Get Min of all GID
    int MaxAllGID = Map_->MaxAllGID(); // Get Max of all GID
    for (i=0; i<NumEntries; i++) {
      int LID = -1; // Assume not found
      int Proc = -1;
      int GID = GlobalEntries[i];
      if (GID<MinAllGID) ierr = 1;
      else if (GID>MaxAllGID) ierr = 1;
      else {
	// Guess uniform distribution and start a little above it
	int Proc1 = EPETRA_MIN(GID/EPETRA_MAX(n_over_p,1) + 2, NumProc-1);
	bool found = false;
	while (Proc1 >= 0 && Proc1< NumProc) {
	  if (AllMinGIDs_[Proc1]<=GID) {
	    if (GID <AllMinGIDs_[Proc1+1]) {
	    found = true;
	    break;
	    }
	    else Proc1++;
	  }
	  else Proc1--;
	}
	if (found) {
	  Proc = Proc1;
	  LID = GID - AllMinGIDs_[Proc];
	}
      }
      Procs[i] = Proc;
      if (LocalEntries!=0) LocalEntries[i] = LID;
    }
    if (EntrySizes!=0) {
      if (Map_->ConstantElementSize()) {
	int ElementSize = Map_->MaxElementSize();
	for (i=0; i<NumEntries; i++) EntrySizes[i] = ElementSize;
      }
      else {
	int * ElementSizeList = Map_->ElementSizeList(); // We know this exists
	
	
	Epetra_Distributor * Size_Distor = Map_->Comm().CreateDistributor();
	
	int Size_num_sends;
	int * Size_send_gids = 0;
	int * Size_send_procs = 0;

	
	EPETRA_CHK_ERR(Size_Distor->CreateFromRecvs( NumEntries, GlobalEntries, Procs, true,
						       Size_num_sends, Size_send_gids, Size_send_procs ));
	
	int * Size_exports = 0;
	int * Size_imports = 0;
	if (Size_num_sends>0) {
	  Size_exports = new int[ 2 * Size_num_sends ];
	  for( i = 0; i < Size_num_sends; i++ )
	    {
	      int Size_curr_GID = Size_send_gids[i];
	      int Size_curr_LID = Map_->LID(Size_curr_GID);
	      assert(Size_curr_LID!=-1); // Internal error 
	      Size_exports[2*i] = Size_curr_GID;
	      int Size_curr_size = ElementSizeList[Size_curr_LID];
	      Size_exports[2*i+1] = Size_curr_size;
	    }
	}
	
	if (NumEntries>0) Size_imports = new int[ 2 * NumEntries ];
	EPETRA_CHK_ERR(Size_Distor->Do( reinterpret_cast<char*> (Size_exports),
			 2 * sizeof( int ), reinterpret_cast<char*> (Size_imports)));
	
	for( i = 0; i < NumEntries; i++ )
	  {

	    // Need to change !!!!
	    //bool found = false;
	    int Size_curr_LID = Size_imports[2*i];
	    for( j = 0; j < NumEntries; j++ )
	      if( Size_curr_LID == GlobalEntries[j] )
		{
		  EntrySizes[j] = Size_imports[2*i+1];
		  // found = true;
		  break;
		}
	    //	if (!found) cout << "Internal error:  Epetra_BasicDirectory::GetDirectoryEntries: Global Index " << curr_LID
	    //	     << " not on processor " << MyPID << endl; abort();
	  }
	
	if( Size_send_gids != 0 ) delete [] Size_send_gids;
	if( Size_send_procs != 0 ) delete [] Size_send_procs;
	
	if( Size_imports != 0 ) delete [] Size_imports;
	if( Size_exports != 0 ) delete [] Size_exports;
	
	delete Size_Distor;
      }
    }
    EPETRA_CHK_ERR(ierr);
    return(0);
  }

  // General case (need to set up an actual directory structure)
  
  int * ElementSizeList = 0;
  int PacketSize = 2; // We will send at least the GID and PID.  Might also send LID and Size info
  bool DoSizes = false;
  if (EntrySizes!=0) {
    if (Map_->ConstantElementSize()) {
      int ElementSize = Map_->MaxElementSize();
	for (i=0; i<NumEntries; i++) EntrySizes[i] = ElementSize;
    }
    else {
      ElementSizeList = Map_->ElementSizeList(); // We know this exists
      DoSizes = true;
      PacketSize++; // Sending Size info
    }
  }

  bool DoLIDs = (LocalEntries!=0); // Do LIDs?
  if (DoLIDs) PacketSize++; // Sending LIDs also

  
  Epetra_Distributor * Distor = DirectoryMap_->Comm().CreateDistributor();
  
  
  int * dir_procs = 0;
  if (NumEntries>0) dir_procs = new int[ NumEntries ];
  
  // Get directory locations for the requested list of entries
  DirectoryMap_->RemoteIDList(NumEntries, GlobalEntries, dir_procs, 0);

  //Check for unfound GlobalEntries and set cooresponding Procs to -1
  int NumMissing = 0;
  for( int i = 0; i < NumEntries; ++i )
    if( dir_procs[i] == -1 )
    {
      Procs[i] = -1;
      if (DoLIDs) LocalEntries[i] = -1;
      ++NumMissing;
    }

  int num_sends;
  int * send_gids = 0;
  int * send_procs = 0;
  
  EPETRA_CHK_ERR(Distor->CreateFromRecvs( NumEntries, GlobalEntries, dir_procs, true,
					   num_sends, send_gids, send_procs));

  if (NumEntries>0) delete [] dir_procs;


  int curr_LID;
  int * exports = 0;
  int * imports = 0;
  if (num_sends>0) {
    exports = new int[ PacketSize * num_sends ];
    int * ptr = exports;
    for( i = 0; i < num_sends; i++ )
      {
	int curr_GID = send_gids[i];
	*ptr++ = curr_GID;
	curr_LID = DirectoryMap_->LID(curr_GID);
	assert(curr_LID!=-1); // Internal error 
	*ptr++ = ProcList_[ curr_LID ];
	if (DoLIDs) *ptr++ = LocalIndexList_[curr_LID];
	if (DoSizes) *ptr++ = SizeList_[curr_LID];
      }
  }

  int NumRecv = NumEntries - NumMissing;
  if (NumRecv>0) imports = new int[PacketSize*NumRecv];
  EPETRA_CHK_ERR(Distor->Do(reinterpret_cast<char*> (exports),
	     PacketSize * sizeof( int ), reinterpret_cast<char*> (imports)));

  int * ptr = imports;
  for( i = 0; i < NumRecv; i++ ) {
    curr_LID = *ptr++;
    for( j = 0; j < NumEntries; j++ )
      if( curr_LID == GlobalEntries[j] ) {
        Procs[j] = *ptr++;
        if (DoLIDs) LocalEntries[j] = *ptr++;
        if (DoSizes) EntrySizes[j] = *ptr++;
        break;
      }
  }
  
  if( send_gids ) delete [] send_gids;
  if( send_procs ) delete [] send_procs;
  
  if( imports ) delete [] imports;
  if( exports ) delete [] exports;
  
  delete Distor;
  return(0);
}

//==============================================================================
void Epetra_BasicDirectory::Print( ostream & os) const {
  
  int MyPID;
  if( DirectoryMap_ != 0 ) {
    bool SizeIsConst = Map_->ConstantElementSize();
    MyPID = DirectoryMap_->Comm().MyPID();
    os << MyPID << " Epetra_BasicDirectory Object: "
      << DirectoryMap_->NumMyElements() << endl;
    for( int i = 0; i < DirectoryMap_->NumMyElements(); i++ ) {
      os << " " << i << " " << ProcList_[i] << " "
	 << LocalIndexList_[i];
      if (!SizeIsConst)
	os  << " " <<  SizeList_[i];
      os << endl;
      os << endl;
    }
  }
  else
  {
    cout << "Epetra_BasicDirectory not setup<<<<<<" << endl;
  }

  return;
}
