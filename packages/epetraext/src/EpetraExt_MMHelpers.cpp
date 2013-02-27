//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include <EpetraExt_ConfigDefs.h>
#include <EpetraExt_MMHelpers.h>
#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>
#include <Epetra_Distributor.h>
#include <Epetra_HashTable.h>
#include <Epetra_Util.h>

#include <Teuchos_TimeMonitor.hpp>

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "Epetra_MpiDistributor.h"
#endif
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define MIN3(x,y,z) ((x)<(y)?(MIN(x,z)):(MIN(y,z)))

namespace EpetraExt {

CrsMatrixStruct::CrsMatrixStruct()
 : numRows(0), numEntriesPerRow(NULL), indices(NULL), values(NULL),
   remote(NULL), numRemote(0), importColMap(NULL), rowMap(NULL), colMap(NULL),
   domainMap(NULL), importMatrix(NULL), origMatrix(NULL)
{
}

CrsMatrixStruct::~CrsMatrixStruct()
{
  deleteContents();
}

void CrsMatrixStruct::deleteContents()
{
  numRows = 0;
  delete [] numEntriesPerRow; numEntriesPerRow = NULL;
  delete [] indices; indices = NULL;
  delete [] values; values = NULL;
  delete [] remote; remote = NULL;
  numRemote = 0;
  delete importMatrix; importMatrix=0;
  // origMatrix is not owned by me, so don't delete
  origMatrix=0;
  targetMapToOrigRow.resize(0);
  targetMapToImportRow.resize(0);
}

int dumpCrsMatrixStruct(const CrsMatrixStruct& M)
{
  cout << "proc " << M.rowMap->Comm().MyPID()<<endl;
  cout << "numRows: " << M.numRows<<endl;
  for(int i=0; i<M.numRows; ++i) {
    for(int j=0; j<M.numEntriesPerRow[i]; ++j) {
      if (M.remote[i]) {
        cout << "  *"<<M.rowMap->GID(i)<<"   "
             <<M.importColMap->GID(M.indices[i][j])<<"   "<<M.values[i][j]<<endl;
      }
      else {
        cout << "   "<<M.rowMap->GID(i)<<"   "
             <<M.colMap->GID(M.indices[i][j])<<"   "<<M.values[i][j]<<endl;
      }
    }
  }
  return(0);
}

CrsWrapper_Epetra_CrsMatrix::CrsWrapper_Epetra_CrsMatrix(Epetra_CrsMatrix& epetracrsmatrix)
 : ecrsmat_(epetracrsmatrix)
{
}

CrsWrapper_Epetra_CrsMatrix::~CrsWrapper_Epetra_CrsMatrix()
{
}

const Epetra_Map&
CrsWrapper_Epetra_CrsMatrix::RowMap() const
{
  return ecrsmat_.RowMap();
}

bool CrsWrapper_Epetra_CrsMatrix::Filled()
{
  return ecrsmat_.Filled();
}

int
CrsWrapper_Epetra_CrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
{
  return ecrsmat_.InsertGlobalValues(GlobalRow, NumEntries, Values, Indices);
}

int
CrsWrapper_Epetra_CrsMatrix::SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
{
  return ecrsmat_.SumIntoGlobalValues(GlobalRow, NumEntries, Values, Indices);
}


//------------------------------------

CrsWrapper_GraphBuilder::CrsWrapper_GraphBuilder(const Epetra_Map& emap)
 : graph_(),
   rowmap_(emap),
   max_row_length_(0)
{
  int num_rows = emap.NumMyElements();
  int* rows = emap.MyGlobalElements();

  for(int i=0; i<num_rows; ++i) {
    graph_[rows[i]] = new std::set<int>;
  }
}

CrsWrapper_GraphBuilder::~CrsWrapper_GraphBuilder()
{
  std::map<int,std::set<int>*>::iterator
    iter = graph_.begin(), iter_end = graph_.end();
  for(; iter!=iter_end; ++iter) {
    delete iter->second;
  }

  graph_.clear();
}

bool CrsWrapper_GraphBuilder::Filled()
{
  return false;
}

int
CrsWrapper_GraphBuilder::InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
{
  std::map<int,std::set<int>*>::iterator
    iter = graph_.find(GlobalRow);

  if (iter == graph_.end()) return(-1);

  std::set<int>& cols = *(iter->second);

  for(int i=0; i<NumEntries; ++i) {
    cols.insert(Indices[i]);
  }

  int row_length = cols.size();
  if (row_length > max_row_length_) max_row_length_ = row_length;

  return(0);
}

int
CrsWrapper_GraphBuilder::SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices)
{
  return InsertGlobalValues(GlobalRow, NumEntries, Values, Indices);
}

std::map<int,std::set<int>*>&
CrsWrapper_GraphBuilder::get_graph()
{
  return graph_;
}

void insert_matrix_locations(CrsWrapper_GraphBuilder& graphbuilder,
                              Epetra_CrsMatrix& C)
{
  int max_row_length = graphbuilder.get_max_row_length();
  if (max_row_length < 1) return;

  std::vector<int> indices(max_row_length);
  int* indices_ptr = &indices[0];
  std::vector<double> zeros(max_row_length, 0.0);
  double* zeros_ptr = &zeros[0];

  std::map<int,std::set<int>*>& graph = graphbuilder.get_graph();

  std::map<int,std::set<int>*>::iterator
    iter = graph.begin(), iter_end = graph.end();

  for(; iter!=iter_end; ++iter) {
    int row = iter->first;
    std::set<int>& cols = *(iter->second);
    int num_entries = cols.size();

    std::set<int>::iterator
      col_iter = cols.begin(), col_end = cols.end();
    for(int j=0; col_iter!=col_end; ++col_iter, ++j) {
      indices_ptr[j] = *col_iter;
    }

    C.InsertGlobalValues(row, num_entries, zeros_ptr, indices_ptr);
  }
}

void pack_outgoing_rows(const Epetra_CrsMatrix& mtx,
                        const std::vector<int>& proc_col_ranges,
                        std::vector<int>& send_rows,
                        std::vector<int>& rows_per_send_proc)
{
  const Epetra_Map& rowmap = mtx.RowMap();
  int numrows = mtx.NumMyRows();
  const Epetra_CrsGraph& graph = mtx.Graph();
  int rowlen = 0;
  int* col_indices = NULL;
  int num_col_ranges = proc_col_ranges.size()/2;
  rows_per_send_proc.resize(num_col_ranges);
  send_rows.clear();
  for(int nc=0; nc<num_col_ranges; ++nc) {
    int first_col = proc_col_ranges[nc*2];
    int last_col = proc_col_ranges[nc*2+1];
    int num_send_rows = 0;
    for(int i=0; i<numrows; ++i) {
      int grow = rowmap.GID(i);
      if (mtx.Filled()) {
        const Epetra_Map& colmap = mtx.ColMap();
        graph.ExtractMyRowView(i, rowlen, col_indices);
        for(int j=0; j<rowlen; ++j) {
          int col = colmap.GID(col_indices[j]);
          if (first_col <= col && last_col >= col) {
            ++num_send_rows;
            send_rows.push_back(grow);
            break;
          }
        }
      }
      else {
        graph.ExtractGlobalRowView(grow, rowlen, col_indices);
        for(int j=0; j<rowlen; ++j) {
          if (first_col <= col_indices[j] && last_col >= col_indices[j]) {
            ++num_send_rows;
            send_rows.push_back(grow);
            break;
          }
        }
      }
    }
    rows_per_send_proc[nc] = num_send_rows;
  }
}

std::pair<int,int> get_col_range(const Epetra_Map& emap)
{
  return std::make_pair(emap.MinMyGID(),emap.MaxMyGID());
}


std::pair<int,int> get_col_range(const Epetra_CrsMatrix& mtx)
{
  std::pair<int,int> col_range;
  if (mtx.Filled()) {
    col_range = get_col_range(mtx.ColMap());
  }
  else {
    const Epetra_Map& row_map = mtx.RowMap();
    col_range.first = row_map.MaxMyGID();
    col_range.second = row_map.MinMyGID();
    int rowlen = 0;
    int* col_indices = NULL;
    const Epetra_CrsGraph& graph = mtx.Graph();
    for(int i=0; i<row_map.NumMyElements(); ++i) {
      graph.ExtractGlobalRowView(row_map.GID(i), rowlen, col_indices);
      for(int j=0; j<rowlen; ++j) {
        if (col_indices[j] < col_range.first) col_range.first = col_indices[j];
        if (col_indices[j] > col_range.second) col_range.second = col_indices[j];
      }
    }
  }

  return col_range;
}

/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
#ifdef HAVE_MPI
template <typename MyType>
void boundary_exchange(const Epetra_MpiComm Comm, MPI_Datatype DataType, 
		       int NumSends, const int * SendProcs, const int * SendSizes, MyType* SendBuffer, 
		       int NumRecvs, const int * RecvProcs, const int * RecvSizes, MyType* RecvBuffer,int SizeOfPacket,int msg_tag)
{  

  MPI_Comm comm = Comm.Comm();
  std::vector<MPI_Request> requests(NumRecvs);
  std::vector<MPI_Status>  status(NumRecvs);

  int i,num_waits=0,MyPID=Comm.MyPID();
  int start, self_recv_len=-1,self_recv_start=-1, self_send_start=-1;
  
  // Default send/recv size if the Sizes arrays are NULL.
  int mysendsize=1, myrecvsize=1;

  // Post Recvs
  start=0;
  for(i=0; i<NumRecvs; i++){
    if(RecvSizes) myrecvsize=RecvSizes[i]*SizeOfPacket;
    if(RecvProcs[i] != MyPID) {
      MPI_Irecv(RecvBuffer + start, myrecvsize, DataType, RecvProcs[i], msg_tag, comm, &requests[num_waits]);
      num_waits++;
    }
    else {
      self_recv_len = myrecvsize;
      self_recv_start=start;
    }	
    start+=myrecvsize;
  }

  // Do sends
  start=0;
  for(i=0; i<NumSends; i++){
    if(SendSizes) mysendsize=SendSizes[i]*SizeOfPacket;
    if(SendProcs[i] != MyPID)
      MPI_Send(SendBuffer + start, mysendsize,DataType,SendProcs[i],msg_tag,comm);
    else 
      self_send_start=start;
    start+=mysendsize;
  }

  // Self-copy (if needed)
  if(self_recv_len != -1)
    memcpy(RecvBuffer+self_recv_start,SendBuffer+self_send_start,self_recv_len*sizeof(MyType)*SizeOfPacket);

  // Wait
  if(NumRecvs > 0)
    MPI_Waitall(num_waits, &requests[0],&status[0]);
}
#endif


#ifdef HAVE_MPI
template <typename MyType>
void boundary_exchange_varsize(const Epetra_MpiComm Comm, MPI_Datatype DataType, 
			       int NumSends, const int * SendProcs, const int * SendSizes, MyType* SendBuffer, 
			       int NumRecvs, const int * RecvProcs, int * RecvSizes, MyType*& RecvBuffer,int SizeOfPacket,int msg_tag)
{  

  int i,rbuffersize=0;

  // Do a first round of boundary exchange with the the SendBuffer sizes
  boundary_exchange<int>(Comm,MPI_INT,NumSends,SendProcs,(int*)0,const_cast<int*>(SendSizes),NumRecvs,RecvProcs,(int*)0,RecvSizes,1,msg_tag);

  // Allocate the RecvBuffer 
  for(i=0; i<NumRecvs; i++) rbuffersize+=RecvSizes[i]*SizeOfPacket;
  RecvBuffer = new MyType[rbuffersize];

  // Do a second round of boundary exchange to trade the actual values
  boundary_exchange<MyType>(Comm,DataType,NumSends,SendProcs,SendSizes,SendBuffer,NumRecvs,RecvProcs,RecvSizes,RecvBuffer,SizeOfPacket,msg_tag+100);
}
#endif



/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/
int sort_crs_entries(int NumRows, const int *CRS_rowptr, int *CRS_colind, double *CRS_vals){
  // For each row, sort column entries from smallest to largest.
  // Use shell sort. Stable sort so it is fast if indices are already sorted.
  // Code copied from  Epetra_CrsMatrix::SortEntries() 
  for(int i = 0; i < NumRows; i++){
    int start=CRS_rowptr[i];

    double* locValues = &CRS_vals[start];
    int NumEntries    = CRS_rowptr[i+1] - start;
    int* locIndices   = &CRS_colind[start];
		
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
  return(0);
}


//=========================================================================
//=========================================================================
//=========================================================================
LightweightMapData::LightweightMapData():
  Epetra_Data(),
  IndexBase_(0),
  LIDHash_(0),
  CopyMap_(0)
{
}
//=========================================================================
LightweightMapData::~LightweightMapData(){
  delete LIDHash_;
  delete CopyMap_;
}

//=========================================================================
LightweightMap::LightweightMap():Data_(0){;}

//=========================================================================
LightweightMap::LightweightMap(int NumGlobalElements,int NumMyElements, const int * MyGlobalElements, int IndexBase, bool GenerateHash)
{
  Data_=new LightweightMapData();
  Data_->MyGlobalElements_.resize(NumMyElements);

  // Build the hash table
  if(GenerateHash) Data_->LIDHash_ = new Epetra_HashTable<int>(NumMyElements + 1 );
    for(int i=0; i < NumMyElements; ++i ) {
      Data_->MyGlobalElements_[i]=MyGlobalElements[i];
      if(GenerateHash) Data_->LIDHash_->Add(MyGlobalElements[i], i);
    }  
}

//=========================================================================
LightweightMap::LightweightMap(const Epetra_Map & Map)
{
  Data_=new LightweightMapData();
  Data_->CopyMap_=new Epetra_Map(Map);
}

//=========================================================================
LightweightMap::LightweightMap(const LightweightMap& map)
  : Data_(map.Data_)
{
  Data_->IncrementReferenceCount();
}

//=========================================================================
LightweightMap::~LightweightMap(){
  CleanupData();
}

//=========================================================================
LightweightMap & LightweightMap::operator= (const LightweightMap & map)
{
  if((this != &map) && (Data_ != map.Data_)) {
    CleanupData();
    Data_ = map.Data_;
    Data_->IncrementReferenceCount();
  }
  return(*this);
}

//=========================================================================
void LightweightMap::CleanupData(){
  if(Data_){
    Data_->DecrementReferenceCount();
    if(Data_->ReferenceCount() == 0) {
      delete Data_;
    }
  }
}

//=========================================================================
int LightweightMap::NumMyElements() const {
  if(Data_->CopyMap_) return Data_->CopyMap_->NumMyElements();
  return Data_->MyGlobalElements_.size();
}

//=========================================================================
int LightweightMap::LID(int GID) const {
  if(Data_->CopyMap_) return Data_->CopyMap_->LID(GID);
  return Data_->LIDHash_->Get(GID);
}

//=========================================================================
int LightweightMap::GID(int LID) const {
  if(Data_->CopyMap_) return Data_->CopyMap_->GID(LID);
  if(LID < 0 || LID > (int)Data_->MyGlobalElements_.size()) return -1;
  return Data_->MyGlobalElements_[LID];
}

//=========================================================================
int* LightweightMap::MyGlobalElements() const {
  if(Data_->CopyMap_) return Data_->CopyMap_->MyGlobalElements();
  return const_cast<int*>(&Data_->MyGlobalElements_[0]);
}

//=========================================================================
int LightweightMap::MinLID() const {
  if(Data_->CopyMap_) return Data_->CopyMap_->MinLID();
  else return 0;
}

//=========================================================================
int LightweightMap::MaxLID() const {
  if(Data_->CopyMap_) return Data_->CopyMap_->MaxLID();
  else return (int)(Data_->MyGlobalElements_.size()-1);
}


//=========================================================================
//=========================================================================
//=========================================================================
RemoteOnlyImport::RemoteOnlyImport(const Epetra_Import & Importer, LightweightMap & RemoteOnlyTargetMap)
{
  int i;

  // Build an "Importer" that only takes the remote parts of the Importer.
  SourceMap_=&Importer.SourceMap();
  TargetMap_=&RemoteOnlyTargetMap;

  // Pull data from the Importer
  NumSend_            = Importer.NumSend();
  NumRemoteIDs_       = Importer.NumRemoteIDs();
  NumExportIDs_       = Importer.NumExportIDs();
  Distor_             = &Importer.Distributor();
  int * OldRemoteLIDs = Importer.RemoteLIDs();
  int * OldExportLIDs = Importer.ExportLIDs();
  int * OldExportPIDs = Importer.ExportPIDs();
  
  // Sanity Check
  if(NumRemoteIDs_ != RemoteOnlyTargetMap.NumMyElements())
    throw std::runtime_error("RemoteOnlyImport: Importer doesn't match RemoteOnlyTargetMap for number of remotes.");

  // Copy the ExportIDs_, since they don't change
  ExportLIDs_ = new int[NumExportIDs_];
  ExportPIDs_ = new int[NumExportIDs_];
  for(i=0; i<NumExportIDs_; i++) {
    ExportLIDs_[i] = OldExportLIDs[i];
    ExportPIDs_[i] = OldExportPIDs[i];
  }

  // The RemoteIDs, on the other hand, do change.  So let's do this right.
  // Note: We might be able to bypass the LID call by just indexing off the Same and Permute GIDs, but at the moment this
  // is fast enough not to worry about it.
  RemoteLIDs_ = new int[NumRemoteIDs_];
  for(i=0; i<NumRemoteIDs_; i++) 
    RemoteLIDs_[i] = TargetMap_->LID(Importer.TargetMap().GID(OldRemoteLIDs[i]));

  // Nowe we make sure these guys are in sorted order.  AztecOO, ML and all that jazz.
  for(i=0; i<NumRemoteIDs_-1; i++) 
    if(RemoteLIDs_[i] > RemoteLIDs_[i+1])
      throw std::runtime_error("RemoteOnlyImport: Importer and RemoteOnlyTargetMap order don't match.");
}

//=========================================================================
RemoteOnlyImport::~RemoteOnlyImport()
{
  delete [] ExportLIDs_;
  delete [] ExportPIDs_;
  delete [] RemoteLIDs_;
  // Don't delete the Distributor, SourceMap_ or TargetMap_ - those were shallow copies
}

//=========================================================================
//=========================================================================
//=========================================================================
template <class GO>
int LightweightCrsMatrix::MakeColMapAndReindex(std::vector<int> owningPIDs, std::vector<GO> Gcolind)
{
  int i,j;

#ifdef ENABLE_MMM_TIMINGS
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor MM(myTime);
  Teuchos::RCP<Teuchos::Time> mtime;
  mtime=MM.getNewTimer("LWCRS C-3.1");
  mtime->start();
#endif

  // Scan all column indices and sort into two groups: 
  // Local:  those whose GID matches a GID of the domain map on this processor and
  // Remote: All others.
  int numDomainElements = DomainMap_.NumMyElements();
  std::vector<bool> LocalGIDs(numDomainElements,false);

  bool DoSizes = !DomainMap_.ConstantElementSize(); // If not constant element size, then error
  if(DoSizes) EPETRA_CHK_ERR(-1);

  // In principle it is good to have RemoteGIDs and RemoteGIDList be as long as the number of remote GIDs
  // on this processor, but this would require two passes through the column IDs, so we make it the max of 100
  // and the number of block rows.
  int numMyBlockRows;
  if(use_lw) numMyBlockRows = RowMapLW_->NumMyElements();
  else numMyBlockRows = RowMapEP_->NumMyElements();

  int  hashsize = numMyBlockRows; if (hashsize < 100) hashsize = 100;
  Epetra_HashTable<GO> RemoteGIDs(hashsize); 
  std::vector<GO>  RemoteGIDList;     RemoteGIDList.reserve(hashsize);
  std::vector<int> RemoteOwningPIDs;  RemoteOwningPIDs.reserve(hashsize);
    
  // In order to do the map reindexing inexpensively, we clobber the GIDs during this pass.  For *local* GID's we clobber them
  // with their LID in the domainMap.  For *remote* GIDs, we clobber them with (numDomainElements+NumRemoteColGIDs) before the increment of
  // the remote count.  These numberings will be separate because no local LID is greater than numDomainElements. 
  int NumLocalColGIDs = 0;
  int NumRemoteColGIDs = 0;
  for(i = 0; i < numMyBlockRows; i++) {
    for(j = rowptr_[i]; j < rowptr_[i+1]; j++) {
      GO GID = Gcolind[j];
      // Check if GID matches a row GID
      int LID = DomainMap_.LID(GID);
      if(LID != -1) {
	bool alreadyFound = LocalGIDs[LID];
	if (!alreadyFound) {
          LocalGIDs[LID] = true; // There is a column in the graph associated with this domain map GID
          NumLocalColGIDs++;
	}
	colind_[j] = LID; 
      }
      else {
	int hash_value=RemoteGIDs.Get(GID);
	if(hash_value  == -1) { // This means its a new remote GID
	  int PID = owningPIDs[j];
	  if(PID==-1) printf("[%d] ERROR: Remote PID should not be -1\n",DomainMap_.Comm().MyPID());
	  colind_[j] = numDomainElements + NumRemoteColGIDs;
	  RemoteGIDs.Add(GID, NumRemoteColGIDs);
	  RemoteGIDList.push_back(GID);
	  RemoteOwningPIDs.push_back(PID);
	  NumRemoteColGIDs++;
	}
	else
	  colind_[j] = numDomainElements + hash_value;	  
      }
    }
  }

  // Possible short-circuit:  If all domain map GIDs are present as column indices, then set ColMap=domainMap and quit
  if (DomainMap_.Comm().NumProc()==1) {     
    if (NumRemoteColGIDs!=0) {
      throw "Some column IDs are not in domainMap.  If matrix is rectangular, you must pass in domainMap to FillComplete";
      // Sanity test: When one processor,there can be no remoteGIDs
    }
    if (NumLocalColGIDs==numDomainElements) {
      ColMap_ = DomainMap_;

      // In this case, we just use the domainMap's indices, which is, not coincidently, what we clobbered colind_ with up above anyway. 
      // No further reindexing is needed.
      return(0); 
    }
  }
      
  // Now build integer array containing column GIDs
  // Build back end, containing remote GIDs, first
  int numMyBlockCols = NumLocalColGIDs + NumRemoteColGIDs;
  Epetra_IntSerialDenseVector Colindices;
  if(numMyBlockCols > 0) 
    Colindices.Size(numMyBlockCols);
  int* RemoteColindices = Colindices.Values() + NumLocalColGIDs; // Points to back end of Colindices

  for(i = 0; i < NumRemoteColGIDs; i++) 
    RemoteColindices[i] = RemoteGIDList[i]; 

  int NLists = 2;

  // Build permute array for *remote* reindexing.
  std::vector<int> RemotePermuteIDs(NumRemoteColGIDs);
  for(i=0; i<NumRemoteColGIDs; i++) RemotePermuteIDs[i]=i;

  // Sort External column indices so that all columns coming from a given remote processor are contiguous
  int* SortLists[3]; // this array is allocated on the stack, and so we won't need to delete it.
  if(NumRemoteColGIDs > 0) {
    SortLists[0] = RemoteColindices;    
    SortLists[1] = &RemotePermuteIDs[0];
    Epetra_Util::Sort(true, NumRemoteColGIDs, &RemoteOwningPIDs[0], 0, 0, NLists, SortLists);
  }


  bool SortGhostsAssociatedWithEachProcessor_ = false;
  if (SortGhostsAssociatedWithEachProcessor_) {
    // Sort external column indices so that columns from a given remote processor are not only contiguous
    // but also in ascending order. NOTE: I don't know if the number of externals associated
    // with a given remote processor is known at this point ... so I count them here.

    NLists=1;
    int StartCurrent, StartNext;
    StartCurrent = 0; StartNext = 1;
    while ( StartNext < NumRemoteColGIDs ) {
      if (RemoteOwningPIDs[StartNext]==RemoteOwningPIDs[StartNext-1]) StartNext++;
      else {
	SortLists[0] =  &RemotePermuteIDs[StartCurrent];
	Epetra_Util::Sort(true,StartNext-StartCurrent, &(RemoteColindices[StartCurrent]),0,0,NLists,SortLists);
        StartCurrent = StartNext; StartNext++;
      }
    }
    SortLists[0] =  &RemotePermuteIDs[StartCurrent];
    Epetra_Util::Sort(true, StartNext-StartCurrent, &(RemoteColindices[StartCurrent]), 0, 0, NLists, SortLists);
  }

  // Reverse the permutation to get the information we actually care about
  std::vector<int> ReverseRemotePermuteIDs(NumRemoteColGIDs);
  for(i=0; i<NumRemoteColGIDs; i++) ReverseRemotePermuteIDs[RemotePermuteIDs[i]]=i;

  // Build permute array for *local* reindexing.
  bool use_local_permute=false;
  std::vector<int> LocalPermuteIDs(numDomainElements);

  // Now fill front end. Two cases:
  // (1) If the number of Local column GIDs is the same as the number of Local domain GIDs, we
  //     can simply read the domain GIDs into the front part of Colindices, otherwise 
  // (2) We step through the GIDs of the domainMap, checking to see if each domain GID is a column GID.
  //     we want to do this to maintain a consistent ordering of GIDs between the columns and the domain.

  if(NumLocalColGIDs == DomainMap_.NumMyElements()) {
    DomainMap_.MyGlobalElements(Colindices.Values()); // Load Global Indices into first numMyBlockCols elements column GID list
  }
  else {
    int* MyGlobalElements = DomainMap_.MyGlobalElements();
    int NumLocalAgain = 0;
    use_local_permute = true;    
    for(i = 0; i < numDomainElements; i++) {
      if(LocalGIDs[i]) {
	LocalPermuteIDs[i] = NumLocalAgain;
	Colindices[NumLocalAgain++] = MyGlobalElements[i];
      }
    }
    assert(NumLocalAgain==NumLocalColGIDs); // Sanity test
  }


  // Copy the remote PID list correctly
  ColMapOwningPIDs_.resize(numMyBlockCols);
  ColMapOwningPIDs_.assign(numMyBlockCols,DomainMap_.Comm().MyPID());
  for(int i=0;i<NumRemoteColGIDs;i++)
    ColMapOwningPIDs_[NumLocalColGIDs+i] = RemoteOwningPIDs[i];

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("LWCRS C-3.2");
  mtime->start();
#endif

  // Make Column map with same element sizes as Domain map 
  LightweightMap temp(-1, numMyBlockCols, Colindices.Values(), DomainMap_.IndexBase());
  ColMap_ = temp;


#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("LWCRS C-3.3");
  mtime->start();
#endif

  // Low-cost reindex of the matrix
  for(i=0; i<numMyBlockRows; i++){
    for(j=rowptr_[i]; j<rowptr_[i+1]; j++){
      int ID=colind_[j];
      if(ID < numDomainElements){
	if(use_local_permute) colind_[j] = LocalPermuteIDs[colind_[j]];
	// In the case where use_local_permute==false, we just copy the DomainMap's ordering, which it so happens
	// is what we put in colind_ to begin with.
      }
      else
	colind_[j] =  NumLocalColGIDs + ReverseRemotePermuteIDs[colind_[j]-numDomainElements];
    }
  }

  assert((size_t)ColMap_.NumMyElements() == ColMapOwningPIDs_.size());

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
#endif

  return(0);
}


//=========================================================================
int LightweightCrsMatrix::PackAndPrepareWithOwningPIDs(const Epetra_DistObject & Source, 
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

  // Rest of work can be done using RowMatrix only  
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
    throw std::runtime_error("EpetraExt::LightweightCrsMatrix::PackAndPrepare: Unable to determine source global index type");

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


//=========================================================================
// Template params are <PID,GID>
static inline bool lessthan12(std::pair<int,int> i, std::pair<int,int> j){
  return ((i.first<j.first) || (i.first==j.first && i.second <j.second));
}


template<typename ImportType>
int LightweightCrsMatrix::PackAndPrepareReverseComm(const Epetra_CrsMatrix & SourceMatrix, ImportType & RowImporter,
						    std::vector<int> &ReverseSendSizes, std::vector<int> &ReverseSendBuffer) {
#ifdef HAVE_MPI
  // Buffer pairs are in (PID,GID) order
  int i,j,k;
  const Epetra_Import *MyImporter= SourceMatrix.Importer();
  if(MyImporter == 0) return -1;
  const Epetra_MpiComm * MpiComm        = dynamic_cast<const Epetra_MpiComm*>(&SourceMatrix.Comm());
  int MyPID                             = MpiComm->MyPID();

  // Things related to messages I and sending in forward mode (RowImporter)
  int NumExportIDs                      = RowImporter.NumExportIDs();
  int* ExportLIDs                       = RowImporter.ExportLIDs();
  int* ExportPIDs                       = RowImporter.ExportPIDs();

  // Things related to messages I am sending in reverse mode (MyImporter)
  Epetra_Distributor& Distor            = MyImporter->Distributor();
  const Epetra_MpiDistributor * MDistor = dynamic_cast<Epetra_MpiDistributor*>(&Distor);
  int NumRecvs                          = MDistor->NumReceives();
  int* RemoteLIDs                       = MyImporter->RemoteLIDs();
  const int * ProcsFrom                 = MDistor->ProcsFrom();
  const int * LengthsFrom               = MDistor->LengthsFrom();

  
  // Get the owning pids in a special way, s.t. ProcsFrom[RemotePIDs[i]] is the guy who owns 
  // RemoteLIDs[j]....
  std::vector<int> RemotePIDOrder(SourceMatrix.NumMyCols(),-1);

  // Now, for each remote ID, record which processor (in ProcsFrom ordering) owns it.
  for(i=0,j=0;i<NumRecvs;i++){
    for(k=0;k<LengthsFrom[i];k++){
      int pid=ProcsFrom[i];
      if(pid!=MyPID) RemotePIDOrder[RemoteLIDs[j]]=i;
      j++;
    }    
  }

  // Step One: Start tacking the (GID,PID) pairs on the std sets
  std::vector<std::set<std::pair<int,int> > > ReversePGIDs(NumRecvs);
  int    *rowptr, *colind;
  double *vals;
  EPETRA_CHK_ERR(SourceMatrix.ExtractCrsDataPointers(rowptr,colind,vals));
  

  // Loop over each exported row and add to the temp list
  for(i=0; i < NumExportIDs; i++) {
    int lid = ExportLIDs[i];
    int exp_pid = ExportPIDs[i];    
    for(j=rowptr[lid]; j<rowptr[lid+1]; j++){
      int pid_order = RemotePIDOrder[colind[j]];
      if(pid_order!=-1) {
	int gid = SourceMatrix.GCID(colind[j]);
	// This GID is getting shipped off somewhere
	ReversePGIDs[pid_order].insert(std::pair<int,int>(exp_pid,gid));
      }
    }
  }

  // Step 2: Count sizes (one too large to avoid std::vector errors)
  ReverseSendSizes.resize(NumRecvs+1);
  int totalsize=0;
  for(i=0; i<NumRecvs; i++) {
    ReverseSendSizes[i] = 2*ReversePGIDs[i].size();
    totalsize += ReverseSendSizes[i];
  }

  // Step 3: Alloc and fill the send buffer (one too large to avoid std::vector errors)
  ReverseSendBuffer.resize(totalsize+1);
  for(i=0, j=0; i<NumRecvs; i++) {
    for(std::set<std::pair<int,int> >::iterator it=ReversePGIDs[i].begin(); it!=ReversePGIDs[i].end(); it++) {
      ReverseSendBuffer[j]   =  it->first;
      ReverseSendBuffer[j+1] =  it->second;
      j+=2;
    }
  }
#endif

  return 0;
}

//=========================================================================
int build_type3_exports(int MyPID,int Nrecv, Epetra_BlockMap &DomainMap, std::vector<int> &ReverseRecvSizes, const int *ReverseRecvBuffer,  std::vector<int> &ExportLID3, std::vector<int> &ExportPID3){
  int i,j;

  // Estimate total length of procs_to for Type 3
  int total_length3=0;
  for(i=0; i<Nrecv; i++)
    total_length3+=ReverseRecvSizes[i]/2;
  if(total_length3==0) return 0;

  std::vector<int> ExportGID3(total_length3);
  ExportLID3.resize(total_length3);
  ExportPID3.resize(total_length3);

  // Build a sorted colmap-style list for Type3 (removing any self-sends)    
  for(i=0,j=0; i<2*total_length3; i+=2) {
    if(ReverseRecvBuffer[i] != MyPID){
      ExportPID3[j]=ReverseRecvBuffer[i];
      ExportGID3[j]=ReverseRecvBuffer[i+1];
      j++;
    }    
  }
  total_length3=j;

  if(total_length3==0) return 0;

  // Sort (ala Epetra_CrsGraph)
  int * companion = &ExportGID3[0];
  Epetra_Util::Sort(true,total_length3,&ExportPID3[0],0,0,1,&companion,0,0);
  int StartCurrent, StartNext;
  StartCurrent = 0; StartNext = 1;
  while ( StartNext < total_length3 ) {
    if(ExportPID3[StartNext] == ExportPID3[StartNext-1]) StartNext++;
    else {
      Epetra_Util::Sort(true,StartNext-StartCurrent, &(ExportGID3[StartCurrent]),0,0,0,0, 0, 0);
      StartCurrent = StartNext; StartNext++;
    }
  }
  Epetra_Util::Sort(true,StartNext-StartCurrent, &(ExportGID3[StartCurrent]),0,0,0,0, 0, 0);
  

  /*  printf("[%d] Type 3 Sorted= ",MyComm.MyPID());
  for(i=0; i<total_length3; i++)
    printf("(--,%2d,%2d) ",ExportGID3[i],ExportPID3[i]);
  printf("\n");
  fflush(stdout);*/


  // Uniq & resize
  for(i=1,j=1; i<total_length3; i++){
    if(ExportPID3[i]!=ExportPID3[i-1] || ExportGID3[i]!=ExportGID3[i-1]){
      ExportPID3[j] = ExportPID3[i];
      ExportGID3[j] = ExportGID3[i];
      j++;
    }
  }
  ExportPID3.resize(j);
  ExportLID3.resize(j);
  total_length3=j;

  /*  printf("[%d] Type 3 UNIQ  = ",MyComm.MyPID());
  for(i=0; i<total_length3; i++)
    printf("(--,%2d,%2d) ",ExportGID3[i],ExportPID3[i]);
  printf("\n");
  fflush(stdout);*/


  
  // Now index down to LIDs   
  for(i=0; i<total_length3; i++) {
    ExportLID3[i]=DomainMap.LID(ExportGID3[i]);
    if(ExportLID3[i] < 0) throw std::runtime_error("LightweightCrsMatrix:MakeExportLists invalid LID");    
  }
    
  /*  printf("[%d] Type 3 FINAL = ",MyComm.MyPID());
  for(i=0; i<total_length3; i++)
    printf("(%2d,%2d,%2d) ",ExportLID3[i],ExportGID3[i],ExportPID3[i]);
  printf("\n");
  fflush(stdout);*/
 
  return total_length3;
}

//=========================================================================
template<typename ImportType>
int build_type2_exports(const Epetra_CrsMatrix & SourceMatrix, ImportType & MyImporter, std::vector<int> &ExportLID2, std::vector<int> &ExportPID2){
  int total_length2=0;
#ifdef HAVE_MPI
  int i,j;
  const Epetra_Import *SourceImporter= SourceMatrix.Importer();

  int    *rowptr, *colind;
  double *vals;
  EPETRA_CHK_ERR(SourceMatrix.ExtractCrsDataPointers(rowptr,colind,vals));

  // Things related to messages I am sending in forward mode (MyImporter)
  int NumExportIDs                      = MyImporter.NumExportIDs();
  const int* ExportLIDs                 = MyImporter.ExportLIDs();
  const int* ExportPIDs                 = MyImporter.ExportPIDs();
  if(NumExportIDs==0) return 0;

  // Assume I own all the cols, then flag any cols I don't own
  // This allows us to avoid LID calls later...
  std::vector<bool> IsOwned(SourceMatrix.NumMyCols(),true);
  if(SourceImporter) {
    const int * RemoteLIDs = SourceImporter->RemoteLIDs();
    // Now flag the cols I don't own
    for(i=0; i<SourceImporter->NumRemoteIDs(); i++) 
      IsOwned[RemoteLIDs[i]]=false;
  }

  // Status vector 
  std::vector<int> SentTo(SourceMatrix.NumMyCols(),-1);

  // Initial allocation: Too large (assume send involves a max size row)
  total_length2 = MyImporter.NumSend() * SourceMatrix.MaxNumEntries() + 1;
  std::vector<int> ExportGID2(total_length2);

  ExportLID2.resize(total_length2);
  ExportPID2.resize(total_length2);

  int current=0, last_start=0, last_pid=ExportPIDs[0];
  for(i=0; i<NumExportIDs; i++){
    // For each row I have to send via MyImporter...
    int row=ExportLIDs[i];
    int pid=ExportPIDs[i];

    if(i!=0 && pid>last_pid) {
      // We have a new PID, so lets finish up the current one      
      if(current!=last_start){
       	int *lids = &ExportLID2[last_start];
	Epetra_Util::Sort(true,current-last_start,&ExportGID2[last_start],0,0,1,&lids,0,0);
	// Note: we don't need to sort the ExportPIDs since they're the same since last_start
      }
      // Reset the list
      last_pid=pid;
      last_start=current;		      
    }
    else if(pid < last_pid) {
      throw std::runtime_error("build_type2_exports: ExportPIDs are not sorted!");
    }

    for(j=rowptr[row]; j<rowptr[row+1]; j++) {
      // For each column in that row...      
      int col=colind[j];
      if(IsOwned[col] && SentTo[col]!=pid){
	// We haven't added this guy to the list yet.
	SentTo[col]         = pid;
	ExportGID2[current] = SourceMatrix.GCID(col);
	ExportLID2[current] = SourceMatrix.DomainMap().LID(ExportGID2[current]);
	ExportPID2[current] = pid;
	current++;
      }
    }
  }//end main loop

  // Final Sort
  int *lids = &ExportLID2[last_start];
  Epetra_Util::Sort(true,current-last_start,&ExportGID2[last_start],0,0,1,&lids,0,0);
  // Note: we don't need to sort the ExportPIDs since they're the same since last_start

  total_length2=current;
  ExportLID2.resize(total_length2);
  ExportPID2.resize(total_length2);
#endif
  return total_length2;
}

//=========================================================================
int build_type1_exports(const Epetra_Import * Importer1, std::vector<int> &ExportLID1, std::vector<int> &ExportPID1){
  int i, total_length1=0;
  if(!Importer1) return 0;
  total_length1 = Importer1->NumSend();
  if(total_length1==0) return 0;

  std::vector<int> ExportGID1(total_length1);
  ExportLID1.resize(total_length1);
  ExportPID1.resize(total_length1);
  const int * ExportLID1Base = Importer1->ExportLIDs();
  const int * ExportPID1Base = Importer1->ExportPIDs();

  for(i=0; i<total_length1; i++){
    ExportLID1[i] = ExportLID1Base[i];
    ExportPID1[i] = ExportPID1Base[i];
    ExportGID1[i] = Importer1->SourceMap().GID(ExportLID1Base[i]);
  }

  // Sort (ala Epetra_CrsGraph)
  int * companion[2] = {&ExportLID1[0],&ExportGID1[0]};
  Epetra_Util::Sort(true,total_length1,&ExportPID1[0],0,0,2,&companion[0],0,0);

  int StartCurrent, StartNext;
  StartCurrent = 0; StartNext = 1;
  while ( StartNext < total_length1 ) {
    if(ExportPID1[StartNext] == ExportPID1[StartNext-1]) StartNext++;
    else {
      int *new_companion = {&ExportLID1[StartCurrent]};
      Epetra_Util::Sort(true,StartNext-StartCurrent, &(ExportGID1[StartCurrent]),0,0,1,&new_companion, 0, 0);
      StartCurrent = StartNext; StartNext++;
    }
  }
  int *new_companion = {&ExportLID1[StartCurrent]};
  Epetra_Util::Sort(true,StartNext-StartCurrent, &(ExportGID1[StartCurrent]),0,0,1,&new_companion, 0, 0);
  return total_length1;
}

//=========================================================================
template<typename ImportType>
int LightweightCrsMatrix::MakeExportLists(const Epetra_CrsMatrix & SourceMatrix, ImportType & Importer2,
					  std::vector<int> &ReverseRecvSizes, const int *ReverseRecvBuffer,
					  std::vector<int> & ExportPIDs, std::vector<int> & ExportLIDs) {
#ifdef HAVE_MPI
  int MyPID = SourceMatrix.Comm().MyPID();

  // This could (legitimately) be zero, in which case we don't have any ReverseRecvs either.
  const Epetra_Import *Importer1= SourceMatrix.Importer();

  // So for each entry in the DomainMap, I have to answer the question: Do I need to send this to anybody?  And if so, to whom?
  //
  // This information comes from three sources:
  // 1) IDs in my DomainMap that are in someone else's ColMap (e.g. SourceMatrix.Importer()).
  // 2) IDs that I own that I sent to another proc in the Forward communication round (e.g. RowImporter).
  // 3) IDs that someone else sent on during the Forward round (and they told me about duing the reverse round).
  //
  // Any of these could legitimately be null.
  Epetra_MpiDistributor * Distor1 = (Importer1)?(dynamic_cast<Epetra_MpiDistributor*>(&Importer1->Distributor())):0;

  int Nsend1 = (Distor1)?(Distor1->NumSends()):0; // Also the number of messages we'll need to parse through in build_type3_exports

  std::vector<int> ExportPID3;
  std::vector<int> ExportLID3;

  std::vector<int> ExportPID2;
  std::vector<int> ExportLID2;

  std::vector<int> ExportPID1;
  std::vector<int> ExportLID1;

  // Build (sorted) ExportLID / ExportGID list for Type
  int Len1=build_type1_exports(Importer1, ExportLID1, ExportPID1);
  int Len2=build_type2_exports<ImportType>(SourceMatrix, Importer2, ExportLID2, ExportPID2);
  int Len3=build_type3_exports(MyPID,Nsend1,DomainMap_,ReverseRecvSizes, ReverseRecvBuffer, ExportLID3, ExportPID3);

  // Since everything should now be sorted correctly, we can do a streaming merge of the three Export lists...
#ifdef HAVE_EPETRAEXT_DEBUG
  {
    int i;
    // Now we sanity check the 1 & 2 lists
    bool test_passed=true;
    for(i=1; i<Len1; i++) {
      if(ExportPID1[i] < ExportPID1[i-1] || (ExportPID1[i] == ExportPID1[i-1] && DomainMap_.GID(ExportLID1[i]) < DomainMap_.GID(ExportLID1[i-1])))
	test_passed=false;
    }
    SourceMatrix.Comm().Barrier();
    if(!test_passed) {
      printf("[%d] Type1 ERRORS  = ",SourceMatrix.Comm().MyPID());
      for(int i=0; i<Len1; i++)
	printf("(%2d,%2d,%2d) ",ExportLID1[i],DomainMap_.GID(ExportLID1[i]),ExportPID1[i]);
      printf("\n");
      fflush(stdout);      
      throw std::runtime_error("Importer1 fails the sanity test");      
    }

    for(i=1; i<Len2; i++) {
      if(ExportPID2[i] < ExportPID2[i-1]  || (ExportPID2[i] == ExportPID2[i-1] && DomainMap_.GID(ExportLID2[i]) < DomainMap_.GID(ExportLID2[i-1])))
	test_passed=false;
    }

    SourceMatrix.Comm().Barrier();
    if(!test_passed) {
      printf("[%d] Type2 ERRORS  = ",SourceMatrix.Comm().MyPID());
      for(int i=0; i<Len2; i++)
	printf("(%2d,%2d,%2d) ",ExportLID2[i],DomainMap_.GID(ExportLID2[i]),ExportPID2[i]);
      printf("\n");
      fflush(stdout);      
      throw std::runtime_error("Importer2 fails the sanity test");      
    }

    for(i=1; i<Len3; i++) {
      if(ExportPID3[i] < ExportPID3[i-1]  || (ExportPID3[i] == ExportPID3[i-1] && DomainMap_.GID(ExportLID3[i]) < DomainMap_.GID(ExportLID3[i-1])))
	test_passed=false;
    }

    SourceMatrix.Comm().Barrier();
    if(!test_passed) {
      printf("[%d] Type3 ERRORS  = ",SourceMatrix.Comm().MyPID());
      for(int i=0; i<Len3; i++)
	printf("(%2d,%2d,%2d) ",ExportLID3[i],DomainMap_.GID(ExportLID3[i]),ExportPID3[i]);
      printf("\n");
      fflush(stdout);      
      throw std::runtime_error("Importer3 fails the sanity test");      
    }


  }
#endif

  if(Importer1 && !Importer1->SourceMap().SameAs(DomainMap_))
    throw std::runtime_error("ERROR: Map Mismatch Importer1");

  if(!Importer2.SourceMap().SameAs(SourceMatrix.RowMap()))
    throw std::runtime_error("ERROR: Map Mismatch Importer2");

  int InfGID = INT_MAX;
  int InfPID = INT_MAX;

  int i1=0, i2=0, i3=0, current=0;

  int MyLen=Len1+Len2+Len3;  
  ExportLIDs.resize(MyLen);
  ExportPIDs.resize(MyLen);

  while(i1 < Len1 || i2 < Len2 || i3 < Len3){
    int PID1 = (i1<Len1)?(ExportPID1[i1]):InfPID;
    int PID2 = (i2<Len2)?(ExportPID2[i2]):InfPID;
    int PID3 = (i3<Len3)?(ExportPID3[i3]):InfPID;

    int GID1 = (i1<Len1)?(DomainMap_.GID(ExportLID1[i1])):InfGID;
    int GID2 = (i2<Len2)?(DomainMap_.GID(ExportLID2[i2])):InfGID;
    int GID3 = (i3<Len3)?(DomainMap_.GID(ExportLID3[i3])):InfGID;

    int MIN_PID = MIN3(PID1,PID2,PID3);
    int MIN_GID = MIN3( ((PID1==MIN_PID)?GID1:InfGID), ((PID2==MIN_PID)?GID2:InfGID), ((PID3==MIN_PID)?GID3:InfGID));
    bool added_entry=false;

    // Case 1: Add off list 1   
    if(PID1 == MIN_PID && GID1 == MIN_GID){
      ExportLIDs[current] = ExportLID1[i1];
      ExportPIDs[current] = ExportPID1[i1];
      current++;
      i1++;
      added_entry=true;
    }

    // Case 2: Add off list 2
    if(PID2 == MIN_PID && GID2 == MIN_GID){
      if(!added_entry) {
	ExportLIDs[current] = ExportLID2[i2];
	ExportPIDs[current] = ExportPID2[i2];
	current++;
	added_entry=true;
      }
      i2++; 
    }

    // Case 3: Add off list 3
    if(PID3 == MIN_PID && GID3 == MIN_GID){
      if(!added_entry) {
	ExportLIDs[current] = ExportLID3[i3];
	ExportPIDs[current] = ExportPID3[i3];
	current++;
      }
      i3++; 
    }
  }// end while
  if(current!=MyLen) {
    ExportLIDs.resize(current);
    ExportPIDs.resize(current);
  }
#endif
  return 0;
}

//=========================================================================
template<typename ImportType>
void LightweightCrsMatrix::Construct(const Epetra_CrsMatrix & SourceMatrix, ImportType & RowImporter)
{  
  // Do we need to use long long for GCIDs?
  bool UseLL=false;
  int SizeofIntType=-1;

  // Get the size of int typpe
  if(SourceMatrix.Map().GlobalIndicesInt()) {
    SizeofIntType = (int)sizeof(int); 
    UseLL=false;
  }
  else if(SourceMatrix.Map().GlobalIndicesLongLong()) {
    SizeofIntType = (int)sizeof(long long); 
    UseLL=true;
  }
  else
    throw std::runtime_error("EpetraExt::LightweightCrsMatrix::PackAndPrepare: Unable to determine source global index type");

#ifdef ENABLE_MMM_TIMINGS
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor MM(myTime);
  Teuchos::RCP<Teuchos::Time> mtime;
  mtime=MM.getNewTimer("LWCRS C-1");
  mtime->start();
#endif


  // Fused constructor, import & FillComplete
  int i,j,rv=0;
  int N;
  if(use_lw) N = RowMapLW_->NumMyElements();
  else N = RowMapEP_->NumMyElements();

  Epetra_Util util;
  int MyPID = SourceMatrix.Comm().MyPID();

#ifdef HAVE_MPI
  std::vector<int> ReverseSendSizes, ReverseSendBuffer;
  std::vector<int> ReverseRecvSizes;
  int * ReverseRecvBuffer=0;
#endif

  bool communication_needed = RowImporter.SourceMap().DistributedGlobal();

  // The basic algorithm here is:
  // 1) Call Distor.Do to handle the import.
  // 2) Copy all the Imported and Copy/Permuted data into the raw Epetra_CrsMatrix / Epetra_CrsGraphData pointers, still using GIDs.
  // 3) Call an optimized version of MakeColMap that avoids the Directory lookups (since the importer knows who owns all the gids) AND
  //    reindexes to LIDs.

  // Sanity Check
  if (!SourceMatrix.RowMap().SameAs(RowImporter.SourceMap())) 
    throw "LightweightCrsMatrix: Fused copy constructor requires Importer.SourceMap() to match SourceMatrix.RowMap()";

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

#ifdef HAVE_MPI
  const Epetra_MpiComm * MpiComm        = dynamic_cast<const Epetra_MpiComm*>(&SourceMatrix.Comm());
  const Epetra_MpiDistributor * MDistor = dynamic_cast<Epetra_MpiDistributor*>(&Distor);
#endif

  // Allocate memory
  rowptr_.resize(N+1);

  /***************************************************/
  /***** 1) From Epetra_DistObject::DoTransfer() *****/
  /***************************************************/
  //  rv=SourceMatrix.CheckSizes(SourceMatrix);

  // NTS: Add CheckSizes stuff here.
  if(rv) throw "LightweightCrsMatrix: Fused copy constructor failed in CheckSizes()";
  
  // Buffers & Other Relevant Info
  char* Exports_  = 0;
  char* Imports_  = 0;
  int LenExports_ = 0;
  int LenImports_ = 0;
  int *Sizes_     = 0;

  int SizeOfPacket; 
  bool VarSizes = false;
  if( NumExportIDs > 0) {
    Sizes_ = new int[NumExportIDs];
  }
  
  rv=PackAndPrepareWithOwningPIDs(SourceMatrix, NumExportIDs, ExportLIDs,LenExports_, Exports_, SizeOfPacket, Sizes_, VarSizes, Distor);
  if(rv) throw "LightweightCrsMatrix: Fused copy constructor failed in PackAndPrepare()";

  if (communication_needed) {
#ifdef HAVE_MPI 
    // Do the exchange of remote data
    int curr_pid;
    const int * ExportPIDs = RowImporter.ExportPIDs();

    // Use the fact that the export procs are sorted to avoid building a hash table.
    // NOTE: The +1's on the message size lists are to avoid std::vector problems if a proc has no sends or recvs. 
    std::vector<int> SendSizes(MDistor->NumSends()+1,0);
    for(i=0, curr_pid=0; i<NumExportIDs; i++) { 
      if(i>0 &&  ExportPIDs[i] > ExportPIDs[i-1]) curr_pid++;
      SendSizes[curr_pid] +=Sizes_[i];

      // sanity check
      if(i>0 &&  ExportPIDs[i] < ExportPIDs[i-1]) throw "ExportPIDs not sorted";
    }

    LenImports_=0; // remove compiler warning.

    std::vector<int> RecvSizes(MDistor->NumReceives()+1);
    int msg_tag=MpiComm->GetMpiTag();
    boundary_exchange_varsize<char>(*MpiComm,MPI_CHAR,MDistor->NumSends(),MDistor->ProcsTo(),&SendSizes[0],Exports_,
				    MDistor->NumReceives(),MDistor->ProcsFrom(),&RecvSizes[0],Imports_,SizeOfPacket,msg_tag);

    // If the  source matrix doesn't have an importer, then nobody sent data belonging to me in the forward round.
    if(SourceMatrix.Importer()) {
      Epetra_Import* SourceImporter=const_cast<Epetra_Import*>(SourceMatrix.Importer());
      const Epetra_MpiDistributor * MyDistor = dynamic_cast<Epetra_MpiDistributor*>(&SourceImporter->Distributor());

      // Setup the reverse communication
      // Note: Buffer pairs are in (PID,GID) order
      PackAndPrepareReverseComm<ImportType>(SourceMatrix,RowImporter,ReverseSendSizes,ReverseSendBuffer);
      
      // Do the reverse communication
      // NOTE: Make the vector one too large to avoid std::vector errors
      ReverseRecvSizes.resize(MyDistor->NumSends()+1);
      int msg_tag=MpiComm->GetMpiTag();
      boundary_exchange_varsize<int>(*MpiComm,MPI_INT,MyDistor->NumReceives(),MyDistor->ProcsFrom(),&ReverseSendSizes[0],&ReverseSendBuffer[0],
				      MyDistor->NumSends(),MyDistor->ProcsTo(),&ReverseRecvSizes[0],ReverseRecvBuffer,1,msg_tag);      
    }
#endif
  }

  if(rv) throw "LightweightCrsMatrix: Fused copy constructor failed in Distor.Do";

  /*********************************************************************/
  /**** 2) Copy all of the Same/Permute/Remote data into CSR_arrays ****/
  /*********************************************************************/
#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("LWCRS C-2");
  mtime->start();
#endif
  // What we really need to know is where in the CSR arrays each row should start (aka the rowptr_).
  // We do that by (a) having each row record it's size in the rowptr_ (b) doing a cumulative sum to get the rowptr_ values correct and 
  // (c) Having each row copied into the right colind_ / values locations.

  // From Epetra_CrsMatrix UnpackAndCombine():
  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumEntries, Number of indices in row.
  // next NumEntries: The actual indices for the row.


  // SameIDs: Always first, always in the same place
  for(i=0; i<NumSameIDs; i++)
    rowptr_[i]=SourceMatrix.NumMyEntries(i);
  
  // PermuteIDs: Still local, but reordered
  for(i=0; i<NumPermuteIDs; i++)
    rowptr_[PermuteToLIDs[i]] = SourceMatrix.NumMyEntries(PermuteFromLIDs[i]);
  
  // RemoteIDs:  RemoteLIDs tells us the ID, we need to look up the length the hard way.  See UnpackAndCombine for where this code came from
  if(NumRemoteIDs > 0) {
    double * dintptr = (double *) Imports_;

    if(UseLL) {
      long long * LLptr = (long long *) dintptr;
      int   NumEntries  = LLptr[1];
      int      IntSize  = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
      for(i=0; i<NumRemoteIDs; i++) {
	rowptr_[RemoteLIDs[i]] = NumEntries;
	
	if( i < (NumRemoteIDs-1) ) {
	  dintptr += IntSize + NumEntries;
	  LLptr = (long long *) dintptr;
	  NumEntries = LLptr[1];
	  IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
	}
      }
    }
    else {
      int    *  intptr = (int *) dintptr;
      int   NumEntries = intptr[1];
      int      IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
      for(i=0; i<NumRemoteIDs; i++) {
	rowptr_[RemoteLIDs[i]] = NumEntries;
	
	if( i < (NumRemoteIDs-1) ) {
	  dintptr += IntSize + NumEntries;
	  intptr = (int *) dintptr;
	  NumEntries = intptr[1];
	  IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
	}
      }
    }
  }

  // Turn row length into a real rowptr_
  int last_len = rowptr_[0];
  rowptr_[0] = 0;
  for(i=1; i<N+1; i++){
    int new_len = rowptr_[i];
    rowptr_[i] = last_len + rowptr_[i-1];
    last_len=new_len;
  }

  // Allocate CSR_colind_ & CSR_values arrays
  int mynnz=rowptr_[N];
  colind_.resize(mynnz);
  vals_.resize(mynnz);

  // Get a list of PIDs ready for colmap, preseed with -1 for local
  std::vector<int> pids(mynnz,-1);

  // Find the PIDs from the Source Matrix (whatever SourceMatrix.Importer() knows)
  std::vector<int> SourcePIDs(SourceMatrix.NumMyCols(),-1);
  if(SourceMatrix.Importer()) util.GetPids(*SourceMatrix.Importer(),SourcePIDs,true);

  // Grab pointers for SourceMatrix
  int    * Source_rowptr, * Source_colind;
  double * Source_vals;
  rv=SourceMatrix.ExtractCrsDataPointers(Source_rowptr,Source_colind,Source_vals);
  if(rv) throw "LightweightCrsMatrix: Fused copy constructor failed in ExtractCrsDataPointers";

  // SameIDs: Copy the data over
  for(i=0; i<NumSameIDs; i++) {
    int FromRow = Source_rowptr[i];
    int ToRow   = rowptr_[i];

    for(j=Source_rowptr[i]; j<Source_rowptr[i+1]; j++) {
      if (UseLL)  colind_LL_[ToRow + j - FromRow] = SourceMatrix.GCID64(Source_colind[j]);
      else        colind_[ToRow + j - FromRow]    = (int) SourceMatrix.GCID64(Source_colind[j]);
      vals_[ToRow + j - FromRow]   = Source_vals[j];      
      pids[ToRow + j - FromRow]    = SourcePIDs[Source_colind[j]];
    }
  }

  // PermuteIDs: Copy the data over
  for(i=0; i<NumPermuteIDs; i++) {
    int FromLID  = PermuteFromLIDs[i];
    int FromRow = Source_rowptr[FromLID];
    int ToRow   = rowptr_[PermuteToLIDs[i]];

    for(j=Source_rowptr[FromLID]; j<Source_rowptr[FromLID+1]; j++) {
      if (UseLL)  colind_LL_[ToRow + j - FromRow] = SourceMatrix.GCID64(Source_colind[j]);
      else        colind_[ToRow + j - FromRow]    = (int) SourceMatrix.GCID64(Source_colind[j]);
      vals_[ToRow + j - FromRow]   = Source_vals[j];      
      pids[ToRow + j - FromRow]    = SourcePIDs[Source_colind[j]];
    }
  }

  // RemoteIDs: Loop structure following UnpackAndCombine  
  if(NumRemoteIDs > 0) {
    double * dintptr   = (double *) Imports_;    
    
    if(UseLL) {
      long long * LLptr      =  (long long *) dintptr;
      int NumEntries         = LLptr[1];
      int IntSize            = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
      double* valptr         = dintptr + IntSize;
      for (i=0; i<NumRemoteIDs; i++) {	  	  
	int ToLID            = RemoteLIDs[i];
	int StartRow         = rowptr_[ToLID];      
	double * values      = valptr;
	long long * Indices  = LLptr + 2;
	for(j=0; j<NumEntries; j++){
	  colind_LL_[StartRow+j] = Indices[2*j];
	  if(MyPID !=  Indices[2*j+1]) pids[StartRow+j]   = (int) Indices[2*j+1];
	  vals_[StartRow + j]    = values[j];	  
	}
	if( i < (NumRemoteIDs-1) ) {
	  dintptr += IntSize + NumEntries;
	  LLptr = (long long *) dintptr;
	  NumEntries = LLptr[1];
	  IntSize = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
	  valptr = dintptr + IntSize;
	}	
      }
    }
    else {
      int * intptr       =  (int *) dintptr;
      int NumEntries     = intptr[1];
      int IntSize        = 1 + (((2*NumEntries+2)*SizeofIntType)/(int)sizeof(double));
      double* valptr     = dintptr + IntSize;

      for (i=0; i<NumRemoteIDs; i++) {	  	
	int ToLID        = RemoteLIDs[i];
	int StartRow     = rowptr_[ToLID];       
	double * values  = valptr;
	int * Indices    = intptr + 2;
	for(j=0; j<NumEntries; j++){
	  colind_[StartRow + j]  = Indices[2*j];
	  if(MyPID !=  Indices[2*j+1]) pids[StartRow+j]   = Indices[2*j+1];
	  vals_[StartRow + j]    = values[j];
	  
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
  }

  /**************************************************************/
  /**** 3) Call Optimized MakeColMap w/ no Directory Lookups ****/
  /**************************************************************/
#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("LWCRS C-3");
  mtime->start();
#endif

  //Call an optimized version of MakeColMap that avoids the Directory lookups (since the importer knows who owns all the gids).
  if(UseLL) MakeColMapAndReindex<long long>(pids,colind_LL_);
  else MakeColMapAndReindex<int>(pids,colind_);

  /********************************************/
  /**** 4) Make Export Lists for Import    ****/
  /********************************************/
#ifdef HAVE_MPI
  MakeExportLists<ImportType>(SourceMatrix,RowImporter,ReverseRecvSizes,ReverseRecvBuffer,ExportPIDs_,ExportLIDs_);
#endif 

  /********************************************/
  /**** 5) Call sort the entries           ****/
  /********************************************/
  // NOTE: If we have no entries the &blah[0] will cause the STL to die in debug mode
#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("LWCRS C-4");
  mtime->start();
#endif
  if(N>0) sort_crs_entries(N, &rowptr_[0], &colind_[0], &vals_[0]);

  /********************************************/
  /**** 6) Cleanup                         ****/
  /********************************************/
#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("LWCRS C-5");
  mtime->start();
#endif

#ifdef HAVE_MPI
  delete [] ReverseRecvBuffer;
#endif

  delete [] Exports_;
  delete [] Imports_;
  delete [] Sizes_;

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
#endif
 }// end fused copy constructor




//=========================================================================
LightweightCrsMatrix::LightweightCrsMatrix(const Epetra_CrsMatrix & SourceMatrix, RemoteOnlyImport & RowImporter):
  use_lw(true),
  RowMapLW_(0),
  RowMapEP_(0),
  DomainMap_(SourceMatrix.DomainMap())
{ 
#ifdef ENABLE_MMM_TIMINGS
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor MM(myTime);
  Teuchos::RCP<Teuchos::Time> mtime;
  mtime=MM.getNewTimer("LWCRS Total");
  mtime->start();
#endif
  RowMapLW_= new LightweightMap(RowImporter.TargetMap());
   
  Construct<RemoteOnlyImport>(SourceMatrix,RowImporter);
#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
#endif

}


//=========================================================================
LightweightCrsMatrix::LightweightCrsMatrix(const Epetra_CrsMatrix & SourceMatrix, Epetra_Import & RowImporter):
  use_lw(false),
  RowMapLW_(0),
  RowMapEP_(0),
  DomainMap_(SourceMatrix.DomainMap())
{
  RowMapEP_= new Epetra_BlockMap(RowImporter.TargetMap());
  Construct<Epetra_Import>(SourceMatrix,RowImporter);
}


LightweightCrsMatrix::~LightweightCrsMatrix(){
  delete RowMapLW_;
  delete RowMapEP_;
}


}//namespace EpetraExt

