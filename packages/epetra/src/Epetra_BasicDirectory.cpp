
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

#include "Epetra_ConfigDefs.h"
#include "Epetra_BasicDirectory.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include "Epetra_Util.h"

//==============================================================================
// Epetra_BasicDirectory constructor for a Epetra_BlockMap object
Epetra_BasicDirectory::Epetra_BasicDirectory(const Epetra_BlockMap & Map)
  : DirectoryMap_(0),
    ProcList_(0),
    ProcListLists_(0),
    ProcListLens_(0),
    numProcLists_(0),
    entryOnMultipleProcs_(false),
    LocalIndexList_(0),
    SizeList_(0),
    SizeIsConst_(true)
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    ,AllMinGIDs_int_(0)
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    ,AllMinGIDs_LL_(0)
#endif
{
  // Test for simple cases

  // Uniprocessor and local map cases (Nothing to set up)

  if (!(Map.DistributedGlobal())) return;

  // Linear Map case

  else if (Map.LinearMap()) {

    // Build a list of the Minimum global ids for all processors on each processor.
    // Since the map is linear, we know that all GIDs are contiguous on each processor
    // and can be found using the MinGIDs.

    int NumProc = Map.Comm().NumProc();

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
	if(Map.GlobalIndicesInt())
	{
       AllMinGIDs_int_ = new int[NumProc+1];
       int MinMyGID = (int) Map.MinMyGID64();
       Map.Comm().GatherAll(&MinMyGID, AllMinGIDs_int_, 1);
       AllMinGIDs_int_[NumProc] = (int) (1 + Map.MaxAllGID64()); // Set max cap
	}
	else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    if(Map.GlobalIndicesLongLong())
	{
       AllMinGIDs_LL_ = new long long[NumProc+1];
       long long MinMyGID = Map.MinMyGID64();
       Map.Comm().GatherAll(&MinMyGID, AllMinGIDs_LL_, 1);
       AllMinGIDs_LL_[NumProc] = 1 + Map.MaxAllGID64(); // Set max cap
	}
	else
#endif
		throw "Epetra_BasicDirectory::Epetra_BasicDirectory: Unknown map index type";
  }

  // General case.  Need to build a directory via calls to communication functions
  else {

	int flag = -1;
	if(Map.GlobalIndicesInt())
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      flag = Generate<int>(Map);
#else
      throw "Epetra_BasicDirectory::Epetra_BasicDirectory: ERROR, GlobalIndicesInt but no API for it.";
#endif
	else if(Map.GlobalIndicesLongLong())
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
      flag = Generate<long long>(Map);
#else
      throw "Epetra_BasicDirectory::Epetra_BasicDirectory: ERROR, GlobalIndicesLongLong but no API for it.";
#endif

	assert(flag==0);
  }
}

//==============================================================================
// Epetra_BasicDirectory copy constructor
Epetra_BasicDirectory::Epetra_BasicDirectory(const Epetra_BasicDirectory & Directory)
  : DirectoryMap_(0),
    ProcList_(0),
    ProcListLists_(0),
    ProcListLens_(0),
    numProcLists_(0),
    entryOnMultipleProcs_(false),
    LocalIndexList_(0),
    SizeList_(0),
    SizeIsConst_(Directory.SizeIsConst_)
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    ,AllMinGIDs_int_(0)
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    ,AllMinGIDs_LL_(0)
#endif
{
  if (Directory.DirectoryMap_!=0) DirectoryMap_ = new Epetra_Map(Directory.DirectoryMap());

  int Dir_NumMyElements = DirectoryMap_->NumMyElements();

  if (Directory.ProcList_!=0) {
    ProcList_ = new int[Dir_NumMyElements];
    for (int i=0; i<Dir_NumMyElements; i++) ProcList_[i] = Directory.ProcList_[i];
  }
  if (Directory.LocalIndexList_!=0) {
    LocalIndexList_ = new int[Dir_NumMyElements];
    for (int i=0; i<Dir_NumMyElements; i++) LocalIndexList_[i] = Directory.LocalIndexList_[i];
    }
  if (Directory.SizeList_!=0) {
    SizeList_ = new int[Dir_NumMyElements];
    for (int i=0; i<Dir_NumMyElements; i++) SizeList_[i] = Directory.SizeList_[i];
    }
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if (Directory.AllMinGIDs_int_!=0) {
       int NumProc = DirectoryMap_->Comm().NumProc();
       AllMinGIDs_int_ = new int[NumProc+1];
       for (int i=0; i<NumProc+1; i++) AllMinGIDs_int_[i] = Directory.AllMinGIDs_int_[i];
	}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if (Directory.AllMinGIDs_LL_!=0) {
       int NumProc = DirectoryMap_->Comm().NumProc();
       AllMinGIDs_LL_ = new long long[NumProc+1];
       for (int i=0; i<NumProc+1; i++) AllMinGIDs_LL_[i] = Directory.AllMinGIDs_LL_[i];
	}
#endif

  if (Directory.numProcLists_ > 0) {
    int num = Directory.numProcLists_;
    ProcListLens_ = new int[num];
    ProcListLists_ = new int*[num];
    numProcLists_ = num;

    for(int i=0; i<num; ++i) {
      int len = Directory.ProcListLens_[i];
      ProcListLens_[i] = len;

      if (len > 0) {
	ProcListLists_[i] = new int[len];
	const int* dir_list = Directory.ProcListLists_[i];
	for(int j=0; j<len; ++j) {
	  ProcListLists_[i][j] = dir_list[j];
	}
      }
      else ProcListLists_[i] = 0;
    }
  }

  entryOnMultipleProcs_ = Directory.entryOnMultipleProcs_;
}

//==============================================================================
// Epetra_BasicDirectory destructor 
Epetra_BasicDirectory::~Epetra_BasicDirectory()
{
  if (numProcLists_>0) {
    for(int i=0; i<numProcLists_; ++i) {
      if (ProcListLens_[i] > 0) delete [] ProcListLists_[i];
    }
    delete [] ProcListLists_; ProcListLists_ = 0;
    delete [] ProcListLens_;  ProcListLens_ = 0;
    numProcLists_ = 0;
  }

  if( DirectoryMap_ != 0 ) delete DirectoryMap_;
  if( ProcList_ != 0 ) delete [] ProcList_;
  if( LocalIndexList_ != 0 ) delete [] LocalIndexList_;
  if( SizeList_ != 0 ) delete [] SizeList_;
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if( AllMinGIDs_int_ != 0 ) delete [] AllMinGIDs_int_;
  AllMinGIDs_int_ = 0;
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if( AllMinGIDs_LL_ != 0 ) delete [] AllMinGIDs_LL_;
  AllMinGIDs_LL_ = 0;
#endif

  DirectoryMap_ = 0;
  ProcList_ = 0 ;
  LocalIndexList_ = 0;
  SizeList_ = 0;
}

//==============================================================================
void Epetra_BasicDirectory::create_ProcListArrays()
{
  numProcLists_ = DirectoryMap_->NumMyElements();
  ProcListLens_ = new int[numProcLists_];
  ProcListLists_ = new int*[numProcLists_];

  for(int i=0; i<numProcLists_; ++i) {
    ProcListLens_[i] = 0;
    ProcListLists_[i] = 0;
  }
}

//==============================================================================
void Epetra_BasicDirectory::addProcToList(int proc, int LID)
{
  int insertPoint = -1;
  int index = Epetra_Util_binary_search(proc, ProcListLists_[LID],
				    ProcListLens_[LID], insertPoint);
  if (index < 0) {
    int tmp = ProcListLens_[LID];
    Epetra_Util_insert(proc, insertPoint, ProcListLists_[LID],
		       ProcListLens_[LID], tmp, 1);
  }
}

//==============================================================================
// Generate: Generates Directory Tables
template<typename int_type>
int Epetra_BasicDirectory::Generate(const Epetra_BlockMap& Map)
{
  int i;
  SizeIsConst_ = Map.ConstantElementSize();
  int_type MinAllGID = (int_type) Map.MinAllGID64();
  int_type MaxAllGID = (int_type) Map.MaxAllGID64();
  // DirectoryMap will have a range of elements from the minimum to the maximum
  // GID of the user map, and an IndexBase of MinAllGID from the user map
  int_type Dir_NumGlobalElements = MaxAllGID - MinAllGID + 1;

  // Create a uniform linear map to contain the directory
  DirectoryMap_ = new Epetra_Map( Dir_NumGlobalElements, MinAllGID, Map.Comm() );

  int Dir_NumMyElements = DirectoryMap_->NumMyElements(); // Get NumMyElements



  // Allocate Processor list and Local Index List.  Initialize to -1s.

  if (Dir_NumMyElements>0) {
    ProcList_ = new int[ Dir_NumMyElements ];
    LocalIndexList_ = new int[ Dir_NumMyElements ];
    if (!SizeIsConst_) SizeList_ = new int[ Dir_NumMyElements ];
    // Initialize values to -1 in case the user global element list does
    // fill all IDs from MinAllGID to MaxAllGID (e.g., allows global indices to be 
    // all even integers.
    for (i=0; i<Dir_NumMyElements; i++) {
      ProcList_[i] = -1;
      LocalIndexList_[i] = -1;
      if (!SizeIsConst_) SizeList_[i] = -1;
    }
  }

  
  // Get list of processors owning the directory entries for the Map GIDs

  int MyPID = Map.Comm().MyPID();

  int Map_NumMyElements = Map.NumMyElements();
  int * send_procs = 0;
  if (Map_NumMyElements>0) send_procs = new int[Map_NumMyElements];
  int_type * Map_MyGlobalElements = 0;
  Map.MyGlobalElementsPtr(Map_MyGlobalElements);

  EPETRA_CHK_ERR(DirectoryMap_->RemoteIDList(Map_NumMyElements,
					     Map_MyGlobalElements, 
					     send_procs, 0));

  bool det_flag = true;

  int num_recvs=0;
    
  Epetra_Distributor * Distor = Map.Comm().CreateDistributor();

  EPETRA_CHK_ERR(Distor->CreateFromSends( Map_NumMyElements, send_procs, det_flag, num_recvs ));

  if (Map_NumMyElements>0) delete [] send_procs;

  int * export_elements = 0;
  char * c_import_elements = 0;
  int * import_elements = 0;
  int len_import_elements = 0;
  int * ElementSizeList = 0;

  int packetSize = (int) (sizeof(int_type) + 2*sizeof(int))/sizeof(int); // Assume we will send GIDs, PIDs and LIDs (will increase to 4 if also sending sizes)
  if (!SizeIsConst_) packetSize++; // Must send element size info also
 
  if (Map_NumMyElements>0) {
    if (!SizeIsConst_) ElementSizeList = Map.ElementSizeList();
    export_elements = new int[ packetSize * Map_NumMyElements ];
    int * ptr = export_elements;
    for( i = 0; i < Map_NumMyElements; i++ )
      {
	*(int_type*)ptr = Map_MyGlobalElements[i];
	ptr += sizeof(int_type)/sizeof(int);
	*ptr++ = MyPID;
	*ptr++ = i;
	if (!SizeIsConst_) *ptr++ = ElementSizeList[i];
      }
  }

  //if (num_recvs>0) import_elements = new int[ packetSize * num_recvs ];
  //for (i=0; i< packetSize*num_recvs; i++) import_elements[i] = 0;

  EPETRA_CHK_ERR(Distor->Do(reinterpret_cast<char *> (export_elements), 
			    packetSize * (int)sizeof( int ),
			    len_import_elements,
			    c_import_elements ));

  import_elements = reinterpret_cast<int *>(c_import_elements);
  
  //bool MYPID = (Map.Comm().MyPID()==0);
  int curr_LID;
  //if (MYPID) cout << "Processor " << Map.Comm().MyPID()<< "  num_recvs = "<< num_recvs << endl << flush;
  int * ptr = import_elements;
  for( i = 0; i < num_recvs; i++ )
  {
    curr_LID = DirectoryMap_->LID(*(int_type*)ptr); // Convert incoming GID to Directory LID
	ptr += sizeof(int_type)/sizeof(int);
    //if (MYPID) cout << " Receive ID = " << i << "  GID = " << import_elements[3*i] << "  LID = " << curr_LID << endl << flush;
    assert(curr_LID !=-1); // Internal error
    int proc = *ptr++;
    if (ProcList_[curr_LID] >= 0) {
      if (ProcList_[curr_LID] != proc) {
	if (numProcLists_ < 1) {
	  create_ProcListArrays();
	}

	addProcToList(ProcList_[curr_LID], curr_LID);
	addProcToList(proc, curr_LID);

	//leave the lowest-numbered proc in ProcList_[curr_LID].
	ProcList_[curr_LID] = ProcListLists_[curr_LID][0];
      }
    }
    else {
      ProcList_[curr_LID] = proc;
    }
    LocalIndexList_[ curr_LID ] = *ptr++;
    if (!SizeIsConst_) SizeList_[ curr_LID ] = *ptr++;
  }

  int localval, globalval;
  localval = numProcLists_;
  DirectoryMap_->Comm().MaxAll(&localval, &globalval, 1);
  entryOnMultipleProcs_ = globalval > 0 ? true : false;

  if (len_import_elements!=0) delete [] c_import_elements;
  if (export_elements!=0) delete [] export_elements;
  
  delete Distor;
  return(0);
}

//==============================================================================
bool Epetra_BasicDirectory::GIDsAllUniquelyOwned() const
{
  return( !entryOnMultipleProcs_ );
}

//==============================================================================
// GetDirectoryEntries: Get non-local GID references ( procID and localID )
// 			Space should already be allocated for Procs and
//     			LocalEntries.
template<typename int_type>
int Epetra_BasicDirectory::GetDirectoryEntries( const Epetra_BlockMap& Map,
						const int NumEntries,
						const int_type * GlobalEntries,
						int * Procs,
						int * LocalEntries,
						int * EntrySizes,
						bool high_rank_sharing_procs) const
{
  int ierr = 0;
  int j;
  int i;
  int MyPID = Map.Comm().MyPID();
  int NumProc = Map.Comm().NumProc();
  int_type n_over_p = (int_type) (Map.NumGlobalElements64() / NumProc);

  // Test for simple cases

  // Uniprocessor and local map cases

  if (!Map.DistributedGlobal()) {
    int ElementSize = 0;
    int * ElementSizeList = 0;
    bool ConstantElementSize = Map.ConstantElementSize();
    if (ConstantElementSize)
      ElementSize = Map.MaxElementSize();
    else
      ElementSizeList = Map.ElementSizeList();
    for (i=0; i<NumEntries; i++) {
      int LID = Map.LID(GlobalEntries[i]); // Get LID
      // Procs[i] will be MyPID, or -1 if the GID is not owned by this map
      if (LID==-1) {
	Procs[i] = -1; 
	ierr = 1; // Send warning error back that one of the GIDs is not part of this map
      }
      else Procs[i] = MyPID;

      // Put LID in return array if needed
      if (LocalEntries!=0) LocalEntries[i] = LID;
      
      // Fill EntrySizes if needed
      if (EntrySizes!=0) {
	if (ConstantElementSize)
	  EntrySizes[i] = ElementSize;
	else if (LID>-1) 
	  EntrySizes[i] = ElementSizeList[LID];
	else
	  EntrySizes[i] = 0;
      }
    }
    EPETRA_CHK_ERR(ierr);
    return(0);
  }

  // Linear Map case
  if (Map.LinearMap()) {
    
    int_type MinAllGID = (int_type) Map.MinAllGID64(); // Get Min of all GID
    int_type MaxAllGID = (int_type) Map.MaxAllGID64(); // Get Max of all GID
    for (i=0; i<NumEntries; i++) {
      int LID = -1; // Assume not found
      int Proc = -1;
      int_type GID = GlobalEntries[i];
      if (GID<MinAllGID) ierr = 1;
      else if (GID>MaxAllGID) ierr = 1;
      else {
	// Guess uniform distribution and start a little above it
	int Proc1 = (int) EPETRA_MIN(GID/EPETRA_MAX(n_over_p,(int_type)1) + 2, (int_type) NumProc-1);
	bool found = false;
	const int_type* AllMinGIDs_ptr = AllMinGIDs<int_type>();

	while (Proc1 >= 0 && Proc1< NumProc) {
	  if (AllMinGIDs_ptr[Proc1]<=GID) {
	    if (GID <AllMinGIDs_ptr[Proc1+1]) {
	    found = true;
	    break;
	    }
	    else Proc1++;
	  }
	  else Proc1--;
	}
	if (found) {
	  Proc = Proc1;
	  LID = (int) (GID - AllMinGIDs_ptr[Proc]);
	}
      }
      Procs[i] = Proc;
      if (LocalEntries!=0) LocalEntries[i] = LID;
    }
    if (EntrySizes!=0) {
      if (Map.ConstantElementSize()) {
	int ElementSize = Map.MaxElementSize();
	for (i=0; i<NumEntries; i++) EntrySizes[i] = ElementSize;
      }
      else {
	int * ElementSizeList = Map.ElementSizeList(); // We know this exists
	
	
	Epetra_Distributor * Size_Distor = Map.Comm().CreateDistributor();
	
	int Size_num_sends;
	int_type * Size_send_gids = 0;
	int * Size_send_procs = 0;

	
	EPETRA_CHK_ERR(Size_Distor->CreateFromRecvs( NumEntries, GlobalEntries, Procs, true,
						       Size_num_sends, Size_send_gids, Size_send_procs ));
	
	int * Size_exports = 0;
	char * c_Size_imports = 0;
	int * Size_imports = 0;
    int packetSize = (int) (sizeof(int_type) + sizeof(int))/sizeof(int);
	if (Size_num_sends>0) {
	  Size_exports = new int[ packetSize * Size_num_sends ];
	  for( i = 0; i < Size_num_sends; i++ )
	    {
	      int_type Size_curr_GID = Size_send_gids[i];
	      int Size_curr_LID = Map.LID(Size_curr_GID);
	      assert(Size_curr_LID!=-1); // Internal error 
	      *(int_type*)(Size_exports + packetSize*i) = Size_curr_GID;
	      int Size_curr_size = ElementSizeList[Size_curr_LID];
	      *(Size_exports + packetSize*i + (packetSize - 1)) = Size_curr_size;
	    }
	}
	
        int len_Size_imports = 0;
	EPETRA_CHK_ERR(Size_Distor->Do( reinterpret_cast<char*> (Size_exports),
                                        packetSize * (int)sizeof( int ),
                                        len_Size_imports,
                                        c_Size_imports));
	Size_imports = reinterpret_cast<int*>(c_Size_imports);
	
	for( i = 0; i < NumEntries; i++ )
	  {

	    // Need to change !!!!
	    //bool found = false;
	    int_type Size_curr_GID = *(int_type*)(Size_imports + packetSize*i);
	    for( j = 0; j < NumEntries; j++ )
	      if( Size_curr_GID == GlobalEntries[j] )
		{
		  EntrySizes[j] = *(Size_imports + packetSize*i + (packetSize - 1));
		  // found = true;
		  break;
		}
	    //	if (!found) cout << "Internal error:  Epetra_BasicDirectory::GetDirectoryEntries: Global Index " << curr_LID
	    //	     << " not on processor " << MyPID << endl; abort();
	  }
	
	if( Size_send_gids != 0 ) delete [] Size_send_gids;
	if( Size_send_procs != 0 ) delete [] Size_send_procs;
	
	if( len_Size_imports != 0 ) delete [] c_Size_imports;
	if( Size_exports != 0 ) delete [] Size_exports;
	
	delete Size_Distor;
      }
    }
    EPETRA_CHK_ERR(ierr);
    return(0);
  }

  // General case (need to set up an actual directory structure)
  
  int PacketSize = (int) (sizeof(int_type) + sizeof(int))/sizeof(int); // We will send at least the GID and PID.  Might also send LID and Size info
  bool DoSizes = false;
  if (EntrySizes!=0) {
    if (Map.ConstantElementSize()) {
      int ElementSize = Map.MaxElementSize();
	for (i=0; i<NumEntries; i++) EntrySizes[i] = ElementSize;
    }
    else {
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

  //Check for unfound GlobalEntries and set corresponding Procs to -1
  int NumMissing = 0;
  {for( i = 0; i < NumEntries; ++i )
    if( dir_procs[i] == -1 )
    {
      Procs[i] = -1;
      if (DoLIDs) LocalEntries[i] = -1;
      ++NumMissing;
  }}

  int num_sends;
  int_type * send_gids = 0;
  int * send_procs = 0;
  
  EPETRA_CHK_ERR(Distor->CreateFromRecvs( NumEntries, GlobalEntries, dir_procs, true,
					   num_sends, send_gids, send_procs));

  if (NumEntries>0) delete [] dir_procs;


  int curr_LID;
  int * exports = 0;
  char * c_imports = 0;
  int * imports = 0;
  int len_imports = 0;
  if (num_sends>0) {
    exports = new int[ PacketSize * num_sends ];
    int * ptr = exports;
    for( i = 0; i < num_sends; i++ )
      {
	int_type curr_GID = send_gids[i];
	*(int_type*)ptr = curr_GID;
	ptr += sizeof(int_type)/sizeof(int);
	curr_LID = DirectoryMap_->LID(curr_GID);
	assert(curr_LID!=-1); // Internal error 
	if (high_rank_sharing_procs==false) {
	  *ptr++ = ProcList_[ curr_LID ];
	}
	else {
	  //high_rank_sharing_procs==true means that if multiple procs share a
	  //GID, we want to use the proc with highest rank rather than the
	  //proc with lowest rank.
	  if (numProcLists_ > 0) {
	    int num = ProcListLens_[curr_LID];
	    if (num > 1) {
	      *ptr++ = ProcListLists_[curr_LID][num-1];
	    }
	    else {
	      *ptr++ = ProcList_[ curr_LID ];
	    }
	  }
	  else {
	    *ptr++ = ProcList_[ curr_LID ];
	  }
	}

	if (DoLIDs) *ptr++ = LocalIndexList_[curr_LID];
	if (DoSizes) *ptr++ = SizeList_[curr_LID];
      }
  }

  int NumRecv = NumEntries - NumMissing;
  EPETRA_CHK_ERR(Distor->Do(reinterpret_cast<char*> (exports),
                            PacketSize * (int)sizeof( int ),
                            len_imports,
                            c_imports));
  imports = reinterpret_cast<int*>(c_imports);

  //create a sorted copy of the GlobalEntries array, along with a companion
  //array that will allow us to put result arrays (Procs, LocalEntries &
  //EntrySizes) in the same order as the unsorted GlobalEntries array
  int* sortedGE_int = new int[NumEntries*(1 + sizeof(int_type)/sizeof(int))];
  int* offsets = sortedGE_int+NumEntries*sizeof(int_type)/sizeof(int);
  int_type* sortedGE = reinterpret_cast<int_type*>(sortedGE_int);

  for(i=0; i<NumEntries; ++i) {
    offsets[i] = i;
  }

  std::memcpy(sortedGE, GlobalEntries, NumEntries*sizeof(int_type));
  Epetra_Util Utils;
  Utils.Sort(true, NumEntries, sortedGE, 0, 0, 1, &offsets, 0, 0);

  int * ptr = imports;
  int insertPoint; //insertPoint won't be used, but is argument to binary_search

  for( i = 0; i < NumRecv; i++ ) {
    int_type curr_LID = *(int_type*)ptr;
	ptr += sizeof(int_type)/sizeof(int);
    j = Epetra_Util_binary_search(curr_LID, sortedGE, NumEntries, insertPoint);
    if (j > -1) {
      Procs[offsets[j]] = *ptr++;
      if (DoLIDs) LocalEntries[offsets[j]] = *ptr++;
      if (DoSizes) EntrySizes[offsets[j]] = *ptr++;
    }
  }

  delete [] sortedGE_int;

  if( send_gids ) delete [] send_gids;
  if( send_procs ) delete [] send_procs;
  
  if( len_imports ) delete [] c_imports;
  if( exports ) delete [] exports;

  delete Distor;
  return(0);
}

//==============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_BasicDirectory::GetDirectoryEntries( const Epetra_BlockMap& Map,
						const int NumEntries,
						const int * GlobalEntries,
						int * Procs,
						int * LocalEntries,
						int * EntrySizes,
						bool high_rank_sharing_procs) const
{
	if(!Map.GlobalIndicesInt())
		throw "Epetra_BasicDirectory::GetDirectoryEntries: int version can't be called for non int map";

	return GetDirectoryEntries<int>(Map, NumEntries, GlobalEntries, Procs,
		LocalEntries, EntrySizes, high_rank_sharing_procs);
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
// GetDirectoryEntries: Get non-local GID references ( procID and localID )
// 			Space should already be allocated for Procs and
//     			LocalEntries.
int Epetra_BasicDirectory::GetDirectoryEntries( const Epetra_BlockMap& Map,
						const int NumEntries,
						const long long * GlobalEntries,
						int * Procs,
						int * LocalEntries,
						int * EntrySizes,
						bool high_rank_sharing_procs) const
{
	if(!Map.GlobalIndicesLongLong())
		throw "Epetra_BasicDirectory::GetDirectoryEntries: long long version can't be called for non long long map";

	return GetDirectoryEntries<long long>(Map, NumEntries, GlobalEntries, Procs,
		LocalEntries, EntrySizes, high_rank_sharing_procs);
}
#endif

//==============================================================================
void Epetra_BasicDirectory::Print( ostream & os) const {
  
  int MyPID;
  if( DirectoryMap_ != 0 ) {;
    MyPID = DirectoryMap_->Comm().MyPID();
    os << MyPID << " Epetra_BasicDirectory Object: "
      << DirectoryMap_->NumMyElements() << endl;
    for( int i = 0; i < DirectoryMap_->NumMyElements(); i++ ) {
      os << " " << i << " " << ProcList_[i] << " "
	 << LocalIndexList_[i];
      if (!SizeIsConst_)
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

//--------------------------------------------------------------------------------
Epetra_BasicDirectory& Epetra_BasicDirectory::operator=(const Epetra_BasicDirectory& src)
{
  (void)src;
  //not currently supported
  bool throw_error = true;
  if (throw_error) {
    std::cerr << std::endl
	      << "Epetra_BasicDirectory::operator= not supported."
	      << std::endl;
    throw -1;
  }
  return( *this );
}
