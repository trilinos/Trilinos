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

#include "Epetra_MpiDistributor.h"
#include "Epetra_MpiComm.h"

#include <stdexcept>
#include <vector>

//==============================================================================
Epetra_MpiDistributor::Epetra_MpiDistributor(const Epetra_MpiComm & Comm): 
  Epetra_Object("Epetra::MpiDistributor"),
  lengths_to_(0),
  procs_to_(0),
  indices_to_(0),
  size_indices_to_(0),
  lengths_from_(0),
  procs_from_(0),
  indices_from_(0),
  size_indices_from_(0),
  resized_(false),
  sizes_(0),
  sizes_to_(0),
  starts_to_(0),
  starts_to_ptr_(0),
  indices_to_ptr_(0),
  sizes_from_(0),
  starts_from_(0),
  starts_from_ptr_(0),
  indices_from_ptr_(0),
  nrecvs_(0),
  nsends_(0),
  nexports_(0),
  self_msg_(0),
  max_send_length_(0),
  total_recv_length_(0),
  tag_(Comm.GetMpiTag()),
  epComm_(&Comm),
  comm_(Comm.GetMpiComm()),
  request_(0),
  status_(0),
  no_delete_(false),
  send_array_(0),
  send_array_size_(0),
  comm_plan_reverse_(0)
{
}

//==============================================================================
Epetra_MpiDistributor::Epetra_MpiDistributor(const Epetra_MpiDistributor & Distributor):
  Epetra_Object("Epetra::MpiDistributor"),
  lengths_to_(0),
  procs_to_(0),
  indices_to_(0),
  size_indices_to_(Distributor.size_indices_to_),
  lengths_from_(0),
  procs_from_(0),
  indices_from_(0),
  size_indices_from_(Distributor.size_indices_from_),
  resized_(false),
  sizes_(0),
  sizes_to_(0),
  starts_to_(0),
  starts_to_ptr_(0),
  indices_to_ptr_(0),
  sizes_from_(0),
  starts_from_(0),
  starts_from_ptr_(0),
  indices_from_ptr_(0),
  nrecvs_(Distributor.nrecvs_),
  nsends_(Distributor.nsends_),
  nexports_(Distributor.nexports_),
  self_msg_(Distributor.self_msg_),
  max_send_length_(Distributor.max_send_length_),
  total_recv_length_(Distributor.total_recv_length_),
  tag_(Distributor.tag_),
  epComm_(Distributor.epComm_),
  comm_(Distributor.comm_),
  request_(0),
  status_(0),
  no_delete_(Distributor.no_delete_),
  send_array_(0),
  send_array_size_(0),
  comm_plan_reverse_(0)
{
  int i;
  if (nsends_>0) {
    lengths_to_ = new int[nsends_];
    procs_to_ = new int[nsends_];
    for (i=0; i<nsends_; i++) {
      lengths_to_[i] = Distributor.lengths_to_[i];
      procs_to_[i] = Distributor.procs_to_[i];
    }
  }
  if (size_indices_to_>0) {
    indices_to_ = new int[size_indices_to_];
    for (i=0; i<size_indices_to_; i++) {
      indices_to_[i] = Distributor.indices_to_[i];
    }
  }

  if (nrecvs_>0) {
    lengths_from_ = new int[nrecvs_];
    procs_from_ = new int[nrecvs_];
    request_ = new MPI_Request[ nrecvs_ ];
    status_ = new MPI_Status[ nrecvs_ ];
    for (i=0; i<nrecvs_; i++) {
      lengths_from_[i] = Distributor.lengths_from_[i];
      procs_from_[i] = Distributor.procs_from_[i];
    }
  }
  if (size_indices_from_>0) {
    indices_from_ = new int[size_indices_from_];
    for (i=0; i<size_indices_from_; i++) {
      indices_from_[i] = Distributor.indices_from_[i];
    }
  }
}

//==============================================================================
Epetra_MpiDistributor::~Epetra_MpiDistributor()
{
  if( !no_delete_ )
  {
    if( lengths_to_ != 0 ) delete [] lengths_to_;
    if( procs_to_ != 0 ) delete [] procs_to_;
    if( indices_to_ != 0 ) delete [] indices_to_;
    if( lengths_from_ != 0 ) delete [] lengths_from_;
    if( procs_from_ != 0 ) delete [] procs_from_;
    if( indices_from_ != 0 ) delete [] indices_from_;
    if( starts_to_ != 0 ) delete [] starts_to_;
    if( starts_from_ != 0 ) delete [] starts_from_;
  }

  if( sizes_ != 0 ) delete [] sizes_;
  if( sizes_to_ != 0 ) delete [] sizes_to_;
  if( sizes_from_ != 0 ) delete [] sizes_from_;
  if( starts_to_ptr_ != 0 ) delete [] starts_to_ptr_;
  if( starts_from_ptr_ != 0 ) delete [] starts_from_ptr_;
  if( indices_to_ptr_ != 0 ) delete [] indices_to_ptr_;
  if( indices_from_ptr_ != 0 ) delete [] indices_from_ptr_;

  if( request_ != 0 ) delete [] request_;
  if( status_ != 0 ) delete [] status_;

  if( send_array_size_ != 0 ) { delete [] send_array_; send_array_size_ = 0; }

  if( comm_plan_reverse_ != 0 ) delete comm_plan_reverse_;
}


//==============================================================================
int Epetra_MpiDistributor::CreateFromSends( const int & NumExportIDs,
                                            const int * ExportPIDs,
                                            bool Deterministic,
                                            int & NumRemoteIDs )
{
 (void)Deterministic; // Prevent compiler warnings for unused argument.
  int my_proc;
  MPI_Comm_rank( comm_, &my_proc );

  int nprocs;
  MPI_Comm_size( comm_, &nprocs );

  // Do the forward map component
  CreateSendStructures_(my_proc,nprocs,NumExportIDs,ExportPIDs);

  //Invert map to see what msgs are received and what length
  EPETRA_CHK_ERR( ComputeRecvs_( my_proc, nprocs ) );

  if (nrecvs_>0) {
    if( !request_ ) {
      request_ = new MPI_Request[ nrecvs_ ];
      status_ = new MPI_Status[ nrecvs_ ];
    }
  }

  NumRemoteIDs = total_recv_length_;

  return 0;
}

//==============================================================================
int Epetra_MpiDistributor::CreateFromRecvs( const int & NumRemoteIDs,
				   const int * RemoteGIDs,
			           const int * RemotePIDs,
				   bool Deterministic,
			           int & NumExportIDs,
				   int *& ExportGIDs,
				   int *& ExportPIDs )
{
  int my_proc;
  MPI_Comm_rank( comm_, &my_proc );

  int nprocs;
  MPI_Comm_size( comm_, &nprocs );

  EPETRA_CHK_ERR( ComputeSends_( NumRemoteIDs, RemoteGIDs, RemotePIDs, NumExportIDs,
				 ExportGIDs, ExportPIDs, my_proc) );

  int testNumRemoteIDs;
  EPETRA_CHK_ERR( CreateFromSends( NumExportIDs, ExportPIDs,
				   Deterministic, testNumRemoteIDs ) );

  return(0);
}

//==============================================================================
//---------------------------------------------------------------------------
//CreateFromRecvs Method
// - create communication plan given a known list of procs to recv from
//---------------------------------------------------------------------------
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_MpiDistributor::CreateFromRecvs( const int & NumRemoteIDs,
				   const long long * RemoteGIDs,
			           const int * RemotePIDs,
				   bool Deterministic,
			           int & NumExportIDs,
				   long long *& ExportGIDs,
				   int *& ExportPIDs )
{
  int my_proc;
  MPI_Comm_rank( comm_, &my_proc );

  int nprocs;
  MPI_Comm_size( comm_, &nprocs );

  EPETRA_CHK_ERR( ComputeSends_( NumRemoteIDs, RemoteGIDs, RemotePIDs, NumExportIDs,
				 ExportGIDs, ExportPIDs, my_proc) );

  int testNumRemoteIDs;
  EPETRA_CHK_ERR( CreateFromSends( NumExportIDs, ExportPIDs,
				   Deterministic, testNumRemoteIDs ) );

  return(0);
}
#endif



//==============================================================================
//---------------------------------------------------------------------------
//CreateFromSendsAndRecvs Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::CreateFromSendsAndRecvs( const int & NumExportIDs,
						    const int * ExportPIDs,
						    const int & NumRemoteIDs,
						    const int * RemoteGIDs,
						    const int * RemotePIDs,
						    bool Deterministic)
{
  (void)RemoteGIDs;
  (void)Deterministic; // Prevent compiler warnings for unused argument.
  nexports_ = NumExportIDs;

  int my_proc;
  MPI_Comm_rank( comm_, &my_proc );
  int nprocs;
  MPI_Comm_size( comm_, &nprocs );

  // Do the forward map component
  CreateSendStructures_(my_proc,nprocs,NumExportIDs,ExportPIDs);

  // Do the reverse map component
  CreateRecvStructures_(NumRemoteIDs,RemotePIDs);

  return 0;
}
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int  Epetra_MpiDistributor::CreateFromSendsAndRecvs( const int & NumExportIDs,
						     const int * ExportPIDs,
						     const int & NumRemoteIDs,
						     const long long * RemoteGIDs,
						     const int * RemotePIDs,
						     bool Deterministic)
{
  (void)RemoteGIDs;
  (void)Deterministic; // Prevent compiler warnings for unused argument.
  nexports_ = NumExportIDs;

  int my_proc;
  MPI_Comm_rank( comm_, &my_proc );
  int nprocs;
  MPI_Comm_size( comm_, &nprocs );

  // Do the forward map component
  CreateSendStructures_(my_proc,nprocs,NumExportIDs,ExportPIDs);

  // Do the reverse map component
  CreateRecvStructures_(NumRemoteIDs,RemotePIDs);

  return 0;


}
#endif



//==============================================================================
int Epetra_MpiDistributor::CreateSendStructures_(int my_proc,
						 int nprocs,
						 const int & NumExportIDs,
						 const int * ExportPIDs)
{
  nexports_ = NumExportIDs;

  int i;

  // Check to see if items are grouped by processor w/o gaps
  // If so, indices_to -> 0

  // Setup data structures for quick traversal of arrays
  int * starts = new int[ nprocs + 1 ];
  for( i = 0; i < nprocs; i++ )
    starts[i] = 0;

  int nactive = 0;
  bool no_send_buff = true;
  int numDeadIndices = 0; // In some cases the GIDs will not be owned by any processors and the PID will be -1

  for( i = 0; i < NumExportIDs; i++ )
  {
    if( no_send_buff && i && (ExportPIDs[i] < ExportPIDs[i-1]) )
      no_send_buff = false;
    if( ExportPIDs[i] >= 0 )
    {
      ++starts[ ExportPIDs[i] ];
      ++nactive;
    }
    else numDeadIndices++; // Increase the number of dead indices.  Used below to leave these out of the analysis
  }

  self_msg_ = ( starts[my_proc] != 0 ) ? 1 : 0;

  nsends_ = 0;

  if( no_send_buff ) //grouped by processor, no send buffer or indices_to_ needed
  {
    for( i = 0; i < nprocs; ++i )
      if( starts[i] ) ++nsends_;

    if( nsends_ )
    {
      procs_to_ = new int[nsends_];
      starts_to_ = new int[nsends_];
      lengths_to_ = new int[nsends_];
    }

    int index = numDeadIndices;  // Leave off the dead indices (PID = -1)
    int proc;
    for( i = 0; i < nsends_; ++i )
    {
      starts_to_[i] = index;
      proc = ExportPIDs[index];
      procs_to_[i] = proc;
      index += starts[proc];
    }

    if( nsends_ )
      Sort_ints_( procs_to_, starts_to_, nsends_ );

    max_send_length_ = 0;

    for( i = 0; i < nsends_; ++i )
    {
      proc = procs_to_[i];
      lengths_to_[i] = starts[proc];
      if( (proc != my_proc) && (lengths_to_[i] > max_send_length_) )
        max_send_length_ = lengths_to_[i];
    }
  }
  else //not grouped by processor, need send buffer and indices_to_
  {
    if( starts[0] != 0 ) nsends_ = 1;

    for( i = 1; i < nprocs; i++ )
    {
      if( starts[i] != 0 ) ++nsends_;
      starts[i] += starts[i-1];
    }

    for( i = nprocs-1; i != 0; i-- )
      starts[i] = starts[i-1];

    starts[0] = 0;

    if (nactive>0) {
      indices_to_ = new int[ nactive ];
      size_indices_to_ = nactive;
    }

    for( i = 0; i < NumExportIDs; i++ )
    if( ExportPIDs[i] >= 0 )
    {
      indices_to_[ starts[ ExportPIDs[i] ] ] = i;
      ++starts[ ExportPIDs[i] ];
    }

    //Reconstuct starts array to index into indices_to.

    for( i = nprocs-1; i != 0; i-- )
      starts[i] = starts[i-1];
    starts[0] = 0;
    starts[nprocs] = nactive;

    if (nsends_>0) {
      lengths_to_ = new int[ nsends_ ];
      procs_to_ = new int[ nsends_ ];
      starts_to_ = new int[ nsends_ ];
    }

    int j = 0;
    max_send_length_ = 0;

    for( i = 0; i < nprocs; i++ )
      if( starts[i+1] != starts[i] )
      {
        lengths_to_[j] = starts[i+1] - starts[i];
        starts_to_[j] = starts[i];
        if( ( i != my_proc ) && ( lengths_to_[j] > max_send_length_ ) )
          max_send_length_ = lengths_to_[j];
        procs_to_[j] = i;
        j++;
      }
  }
    
  delete [] starts;

  nsends_ -= self_msg_;

  return 0;
}

//==============================================================================
int Epetra_MpiDistributor::CreateRecvStructures_(const int & NumRemoteIDs,
						 const int * RemotePIDs)
{
  int i, j;

  // Since the RemotePIDs should be sorted, counting the total number of recvs should be easy...
  // use nsends as an initial guess for space.
  std::vector<int> recv_list;
  recv_list.reserve(nsends_);
  
  int last_pid=-2;
  for(i=0; i<NumRemoteIDs; i++) { 
    if(RemotePIDs[i]>last_pid) {
      recv_list.push_back(RemotePIDs[i]);
      last_pid = RemotePIDs[i];
    }
    else if (RemotePIDs[i]<last_pid)
      throw std::runtime_error("Epetra_MpiDistributor::CreateRecvStructures_ expected RemotePIDs to be in sorted order");    
  }
  nrecvs_=recv_list.size();

  if (nrecvs_>0) {
    starts_from_  = new int[nrecvs_];
    procs_from_   = new int[nrecvs_];
    lengths_from_ = new int[nrecvs_];
    request_      = new MPI_Request[ nrecvs_ ];
    status_       = new MPI_Status[ nrecvs_ ];
  }

  for(i=0,j=0; i<nrecvs_; ++i) {
    int jlast=j;
    procs_from_[i]  = recv_list[i];
    starts_from_[i] = j;
    for( ; RemotePIDs[jlast]==RemotePIDs[j] && j<NumRemoteIDs ; j++){;}
    lengths_from_[i]=j-jlast;
  }
  total_recv_length_=NumRemoteIDs;

  nrecvs_ -= self_msg_;

  return 0;
}


//==============================================================================
//---------------------------------------------------------------------------
//ComputeRecvs Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::ComputeRecvs_( int my_proc, 
			        int nprocs )
{
  int * msg_count = new int[ nprocs ];
  int * counts = new int[ nprocs ];

  int i;
  MPI_Status status;

  for( i = 0; i < nprocs; i++ )
  {
    msg_count[i] = 0;
    counts[i] = 1;
  }

  for( i = 0; i < nsends_+self_msg_; i++ )
    msg_count[ procs_to_[i] ] = 1;

#if defined(REDUCE_SCATTER_BUG)
// the bug is found in mpich on linux platforms
  MPI_Reduce(msg_count, counts, nprocs, MPI_INT, MPI_SUM, 0, comm_);
  MPI_Scatter(counts, 1, MPI_INT, &nrecvs_, 1, MPI_INT, 0, comm_);
#else
  MPI_Reduce_scatter( msg_count, &nrecvs_, counts, MPI_INT, MPI_SUM, comm_ );
#endif

  delete [] msg_count;
  delete [] counts;

  if (nrecvs_>0) {
    lengths_from_ = new int[nrecvs_];
    procs_from_ = new int[nrecvs_];
    for(i=0; i<nrecvs_; ++i) {
      lengths_from_[i] = 0;
      procs_from_[i] = 0;
    }
  }

#ifndef NEW_COMM_PATTERN
  for( i = 0; i < (nsends_+self_msg_); i++ )
    if( procs_to_[i] != my_proc ) {
      MPI_Send( &(lengths_to_[i]), 1, MPI_INT, procs_to_[i], tag_, comm_ );
    }
    else
    {
      //set self_msg_ to end block of recv arrays
      lengths_from_[nrecvs_-1] = lengths_to_[i];
      procs_from_[nrecvs_-1] = my_proc;
    }

  for( i = 0; i < (nrecvs_-self_msg_); i++ )
  {
    MPI_Recv( &(lengths_from_[i]), 1, MPI_INT, MPI_ANY_SOURCE, tag_, comm_, &status );
    procs_from_[i] = status.MPI_SOURCE;
  }

  MPI_Barrier( comm_ );
#else
  if (nrecvs_>0) {
    if( !request_ ) {
      request_ = new MPI_Request[nrecvs_-self_msg_];
      status_ = new MPI_Status[nrecvs_-self_msg_];
    }
  }

  for( i = 0; i < (nrecvs_-self_msg_); i++ )
    MPI_Irecv( &(lengths_from_[i]), 1, MPI_INT, MPI_ANY_SOURCE, tag_, comm_, &(request_[i]) );

  MPI_Barrier( comm_ );

  for( i = 0; i < (nsends_+self_msg_); i++ )
    if( procs_to_[i] != my_proc ) {
      MPI_Rsend( &(lengths_to_[i]), 1, MPI_INT, procs_to_[i], tag_, comm_ );
    }
    else
    {
      //set self_msg_ to end block of recv arrays
      lengths_from_[nrecvs_-1] = lengths_to_[i];
      procs_from_[nrecvs_-1] = my_proc;
    }

  if( (nrecvs_-self_msg_) > 0 ) MPI_Waitall( (nrecvs_-self_msg_), request_, status_ );

  for( i = 0; i < (nrecvs_-self_msg_); i++ )
    procs_from_[i] = status_[i].MPI_SOURCE;
#endif

  Sort_ints_( procs_from_, lengths_from_, nrecvs_ );

  // Compute indices_from_
  // Seems to break some rvs communication
/* Not necessary since rvs communication is always blocked
  size_indices_from_ = 0;
  if( nrecvs_ > 0 )
  {
    for( i = 0; i < nrecvs_; i++ )  size_indices_from_ += lengths_from_[i];
    indices_from_ = new int[ size_indices_from_ ];

    for (i=0; i<size_indices_from_; i++) indices_from_[i] = i;
  }
*/

  if (nrecvs_>0) starts_from_ = new int[nrecvs_];
  int j = 0;
  for( i=0; i<nrecvs_; ++i )
  {
    starts_from_[i] = j;
    j += lengths_from_[i];
  }

  total_recv_length_ = 0;
  for( i = 0; i < nrecvs_; i++ )
    total_recv_length_ += lengths_from_[i];

  nrecvs_ -= self_msg_;

  MPI_Barrier( comm_ );
  
  return false;
}

//==============================================================================
//---------------------------------------------------------------------------
//ComputeSends Method
//---------------------------------------------------------------------------
template<typename id_type>
int Epetra_MpiDistributor::ComputeSends_( int num_imports,
				const id_type *& import_ids,
				const int *& import_procs,
				int & num_exports,
				id_type *& export_ids,
				int *& export_procs,
				int my_proc ) {
 
  Epetra_MpiDistributor tmp_plan(*epComm_);
  int i;

  int * proc_list = 0;
  int * import_objs = 0;
  char * c_export_objs = 0;
  const int pack_size = (1 + sizeof(id_type)/sizeof(int));

  if( num_imports > 0 )
  {
    proc_list = new int[ num_imports ];
    import_objs = new int[ num_imports * pack_size];

    for( i = 0; i < num_imports; i++ )
    {
      proc_list[i] = import_procs[i];

      *(id_type*)(import_objs + pack_size*i) = import_ids[i];
      *(import_objs + pack_size*i + (pack_size-1)) = my_proc;
    }
  }

  EPETRA_CHK_ERR(tmp_plan.CreateFromSends( num_imports, proc_list,
					   true, num_exports) );
  if( num_exports > 0 )
  {
    //export_objs = new int[ 2 * num_exports ];
    export_ids = new id_type[ num_exports ];
    export_procs = new int[ num_exports ];
  }
  else
  {
    export_ids = 0;
    export_procs = 0;
  }

  int len_c_export_objs = 0;
  EPETRA_CHK_ERR( tmp_plan.Do(reinterpret_cast<char *> (import_objs),
			      pack_size * (int)sizeof( int ), 
			      len_c_export_objs,
			      c_export_objs) );
  int * export_objs = reinterpret_cast<int *>(c_export_objs);

  for( i = 0; i < num_exports; i++ ) {
    export_ids[i] = *(id_type*)(export_objs + pack_size*i);
    export_procs[i] = *(export_objs + pack_size*i + (pack_size-1));
  }

  if( proc_list != 0 ) delete [] proc_list;
  if( import_objs != 0 ) delete [] import_objs;
  if( len_c_export_objs != 0 ) delete [] c_export_objs;

  return(0);

}

//==============================================================================
int Epetra_MpiDistributor::Do( char * export_objs,
                               int obj_size,
                               int & len_import_objs,
                               char *& import_objs )
{
  EPETRA_CHK_ERR( DoPosts(export_objs, obj_size, len_import_objs, import_objs) );
  EPETRA_CHK_ERR( DoWaits() );
  return(0);
}

//==============================================================================
int Epetra_MpiDistributor::DoReverse( char * export_objs,
                                      int obj_size,
                                      int & len_import_objs,
                                      char *& import_objs )
{
  EPETRA_CHK_ERR( DoReversePosts(export_objs, obj_size,
				 len_import_objs, import_objs) );
  EPETRA_CHK_ERR( DoReverseWaits() );
  return(0);
}

//==============================================================================
int Epetra_MpiDistributor::DoPosts( char * export_objs,
				    int obj_size,
                                    int & len_import_objs,
				    char *& import_objs )
{
  int i, j, k;

  int my_proc = 0;
  int self_recv_address = 0;

  MPI_Comm_rank( comm_, &my_proc );

  if( len_import_objs < (total_recv_length_*obj_size) )
  {
    if( import_objs!=0 ) {delete [] import_objs; import_objs = 0;}
    len_import_objs = total_recv_length_*obj_size;
    if (len_import_objs>0) import_objs = new char[len_import_objs];
    for( i=0; i<len_import_objs; ++i ) import_objs[i]=0;
  }

  k = 0;

  j = 0;
  for( i = 0; i < (nrecvs_+self_msg_); i++ )
  {
    if( procs_from_[i] != my_proc )
    {
      MPI_Irecv( &(import_objs[j]),
                 lengths_from_[i] * obj_size,
                 MPI_CHAR, procs_from_[i],
                 tag_, comm_,
                 &(request_[k]) );
      k++;
    }
    else
      self_recv_address = j;

    j += lengths_from_[i] * obj_size;
  }

#ifndef EPETRA_NO_READY_SEND_IN_DO_POSTS
  // NOTE (mfh 19 Mar 2012):
  //
  // The ready-sends below require that each ready-send's matching
  // receive (see above) has already been posted.  We ensure this with
  // a barrier.  (Otherwise, some process that doesn't need to post
  // receives might post its ready-send before the receiving process
  // gets to post its receive.)  If you want to remove the barrier,
  // you'll have to replace the ready-sends below with standard sends
  // or Isends.
  MPI_Barrier( comm_ );
#endif // EPETRA_NO_READY_SEND_IN_DO_POSTS

  //setup scan through procs_to list starting w/ higher numbered procs 
  //Should help balance msg traffic
  int nblocks = nsends_ + self_msg_; 
  int proc_index = 0; 
  while( proc_index < nblocks && procs_to_[proc_index] < my_proc )
    ++proc_index;                    
  if( proc_index == nblocks ) proc_index = 0;
   
  int self_num = 0, self_index = 0;
  int p;

  if( !indices_to_ ) //data already blocked by processor
  {
    for( i = 0; i < nblocks; ++i )
    {
      p = i + proc_index;
      if( p > (nblocks-1) ) p -= nblocks;

      if( procs_to_[p] != my_proc ) {

#ifndef EPETRA_NO_READY_SEND_IN_DO_POSTS
        MPI_Rsend( &export_objs[starts_to_[p]*obj_size],
                   lengths_to_[p]*obj_size,
                   MPI_CHAR,
                   procs_to_[p],
                   tag_,
                   comm_ );
#else
        MPI_Send( &export_objs[starts_to_[p]*obj_size],
                   lengths_to_[p]*obj_size,
                   MPI_CHAR,
                   procs_to_[p],
                   tag_,
                   comm_ );
#endif // EPETRA_NO_READY_SEND_IN_DO_POSTS
      }
      else {
       self_num = p;
      }
    }

    if( self_msg_ )
      memcpy( &import_objs[self_recv_address],
              &export_objs[starts_to_[self_num]*obj_size],
              lengths_to_[self_num]*obj_size );
  }
  else //data not blocked by proc, use send buffer
  {
    if( send_array_size_ < (max_send_length_*obj_size) )
    {
      if( send_array_!=0 ) {delete [] send_array_; send_array_ = 0;}
      send_array_size_ = max_send_length_*obj_size;
      if (send_array_size_>0) send_array_ = new char[send_array_size_];
    }

    j = 0;
    for( i = 0; i < nblocks; i++ )
    {
      p = i + proc_index;
      if( p > (nblocks-1) ) p -= nblocks;
      if( procs_to_[p] != my_proc )
      {
        int offset = 0;
        j = starts_to_[p];
        for( k = 0; k < lengths_to_[p]; k++ )
        {
          memcpy( &(send_array_[offset]), 
                  &(export_objs[indices_to_[j]*obj_size]),
                  obj_size );
          ++j;
          offset += obj_size;
        }
#ifndef EPETRA_NO_READY_SEND_IN_DO_POSTS
        MPI_Rsend( send_array_,
                   lengths_to_[p] * obj_size,
                   MPI_CHAR,
                   procs_to_[p],
                   tag_, comm_ );
#else
        MPI_Send( send_array_,
		  lengths_to_[p] * obj_size,
		  MPI_CHAR,
		  procs_to_[p],
		  tag_, comm_ );
#endif // EPETRA_NO_READY_SEND_IN_DO_POSTS
      }
      else
      {
        self_num = p;
        self_index = starts_to_[p];
      }
    }

    if( self_msg_ )
      for( k = 0; k < lengths_to_[self_num]; k++ )
      {
        memcpy( &(import_objs[self_recv_address]),
                &(export_objs[indices_to_[self_index]*obj_size]),
                obj_size );
        self_index++;
        self_recv_address += obj_size;
      }
  }
  return(0);
}

//==============================================================================
int Epetra_MpiDistributor::DoWaits()
{
  if( nrecvs_ > 0 ) MPI_Waitall( nrecvs_, request_, status_ );

  return(0);
}

//==============================================================================
int Epetra_MpiDistributor::DoReversePosts( char * export_objs,
                                           int obj_size,
                                           int & len_import_objs,
                                           char *& import_objs )
{
  assert(indices_to_==0); //Can only do reverse comm when original data
                          // is blocked by processor

  int i;
  int my_proc = 0;

  MPI_Comm_rank( comm_, &my_proc );

  if( comm_plan_reverse_ == 0 )
  {
    int total_send_length = 0;
    for( i = 0; i < nsends_+self_msg_; i++ )
      total_send_length += lengths_to_[i];

    int max_recv_length = 0;
    for( i = 0; i < nrecvs_; i++ )
      if( procs_from_[i] != my_proc )
        if( lengths_from_[i] > max_recv_length )
          max_recv_length = lengths_from_[i];

    comm_plan_reverse_ = new Epetra_MpiDistributor(*epComm_);

    comm_plan_reverse_->lengths_to_ = lengths_from_;
    comm_plan_reverse_->procs_to_ = procs_from_;
    comm_plan_reverse_->indices_to_ = indices_from_;
    comm_plan_reverse_->starts_to_ = starts_from_;

    comm_plan_reverse_->lengths_from_ = lengths_to_;
    comm_plan_reverse_->procs_from_ = procs_to_;
    comm_plan_reverse_->indices_from_ = indices_to_;
    comm_plan_reverse_->starts_from_ = starts_to_;

    comm_plan_reverse_->nsends_ = nrecvs_;
    comm_plan_reverse_->nrecvs_ = nsends_;
    comm_plan_reverse_->self_msg_ = self_msg_;

    comm_plan_reverse_->max_send_length_ = max_recv_length;
    comm_plan_reverse_->total_recv_length_ = total_send_length;

    comm_plan_reverse_->request_ = new MPI_Request[ comm_plan_reverse_->nrecvs_ ];
    comm_plan_reverse_->status_= new MPI_Status[ comm_plan_reverse_->nrecvs_ ];

    comm_plan_reverse_->no_delete_ = true;
  }

  int comm_flag = comm_plan_reverse_->DoPosts(export_objs, obj_size, len_import_objs, import_objs);

  return(comm_flag);
}

//==============================================================================
int Epetra_MpiDistributor::DoReverseWaits()
{
  if( comm_plan_reverse_ == 0 ) return (-1);

  int comm_flag = comm_plan_reverse_->DoWaits();

  return(comm_flag);
}

//==============================================================================
//---------------------------------------------------------------------------
//Resize Method                      (Heaphy) 
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::Resize_( int * sizes )
{ 
  int  i, j, k;         // loop counters
  int sum; 
 
  //if (sizes == 0) return 0; 
     
  int my_proc; 
  MPI_Comm_rank (comm_, &my_proc);  
  int nprocs;
  MPI_Comm_size( comm_, &nprocs );
     
  if( resized_ )
  {  
    //test and see if we are already setup for these sizes
    bool match = true; 
    for( i = 0; i < nexports_; ++i )
      match = match && (sizes_[i]==sizes[i]);
    int matched = match?1:0;
    int match_count = 0;
    MPI_Allreduce( &matched, &match_count, 1, MPI_INT, MPI_SUM, comm_ );
    if( match_count == nprocs )
      return 0;
    else //reset existing sizing arrays
      max_send_length_ = 0;
  } 
 
  if( !sizes_ && nexports_ ) sizes_ = new int[nexports_]; 
  for (i = 0; i < nexports_; i++)
    sizes_[i] = sizes[i]; 
 
  if( !sizes_to_ && (nsends_+self_msg_) ) sizes_to_ = new int[nsends_+self_msg_];
  for (i = 0; i < (nsends_+self_msg_); ++i) 
    sizes_to_[i] = 0;                       
 
  if( !starts_to_ptr_ && (nsends_+self_msg_) ) starts_to_ptr_ = new int[nsends_+self_msg_];

  if( !indices_to_ ) //blocked sends
  {
    
    int * index = 0;
    int * sort_val = 0;
    if (nsends_+self_msg_>0) {
      index = new int[nsends_+self_msg_];
      sort_val = new int[nsends_+self_msg_];
    }
    for( i = 0; i < (nsends_+self_msg_); ++i )
    {
      j = starts_to_[i];
      for( k = 0; k < lengths_to_[i]; ++k )
        sizes_to_[i] += sizes[j++];
      if( (sizes_to_[i] > max_send_length_) && (procs_to_[i] != my_proc) )
        max_send_length_ = sizes_to_[i];
    }

    for( i = 0; i < (nsends_+self_msg_); ++i )
    {
      sort_val[i] = starts_to_[i];
      index[i] = i;
    }

    if( nsends_+self_msg_ )
      Sort_ints_( sort_val, index, (nsends_+self_msg_) );

    sum = 0;
    for( i = 0; i < (nsends_+self_msg_); ++i )
    {
      starts_to_ptr_[ index[i] ] = sum;
      sum += sizes_to_[ index[i] ];
    }

    if (index!=0) {delete [] index; index = 0;}
    if (sort_val!=0) {delete [] sort_val; sort_val = 0;}
  }
  else //Sends not blocked, so have to do more work
  {
    if( !indices_to_ptr_ && nexports_ ) indices_to_ptr_ = new int[nexports_]; 
    int * offset = 0;
    if( nexports_ ) offset = new int[nexports_];
   
    //Compute address for every item in send array
    sum = 0; 
    for( i = 0; i < nexports_; ++i )
    {  
      offset[i] = sum; 
      sum += sizes_[i];
    }
   
    sum = 0;
    max_send_length_ = 0;
    for( i = 0; i < (nsends_+self_msg_); ++i )
    {
      starts_to_ptr_[i] = sum;
      for( j = starts_to_[i]; j < (starts_to_[i]+lengths_to_[i]); ++j )
      {
        indices_to_ptr_[j] = offset[ indices_to_[j] ];
        sizes_to_[i] += sizes_[ indices_to_[j] ];
      }
      if( sizes_to_[i] > max_send_length_ && procs_to_[i] != my_proc )
        max_send_length_ = sizes_to_[i];
      sum += sizes_to_[i];
    }
 
    if (offset!=0) {delete [] offset; offset = 0;}
  }
 
  //  Exchange sizes routine inserted here:
  int self_index_to = -1;
  total_recv_length_ = 0;
  if( !sizes_from_ && (nrecvs_+self_msg_) ) sizes_from_ = new int [nrecvs_+self_msg_];

#ifndef EPETRA_NEW_COMM_PATTERN
  for (i = 0; i < (nsends_+self_msg_); i++)
  {
    if(procs_to_[i] != my_proc)
      MPI_Send ((void *) &(sizes_to_[i]), 1, MPI_INT, procs_to_[i], tag_, comm_);
    else
      self_index_to = i;
  }

  MPI_Status status;
  for (i = 0; i < (nrecvs_+self_msg_); ++i)
  {
    sizes_from_[i] = 0;
    if (procs_from_[i] != my_proc)
      MPI_Recv((void *) &(sizes_from_[i]), 1, MPI_INT, procs_from_[i], tag_, comm_, &status);
    else
      sizes_from_[i] = sizes_to_[self_index_to];
    total_recv_length_ += sizes_from_[i];
  }
#else
  if (nrecvs_>0 && !request_) {
    request_ = new MPI_Request[ nrecvs_-self_msg_ ];
    status_ = new MPI_Status[ nrecvs_-self_msg_ ];
  }

  for (i = 0; i < (nsends_+self_msg_); i++)
  {
    if(procs_to_[i] == my_proc)
      self_index_to = i;
  }

  for (i = 0; i < (nrecvs_+self_msg_); ++i)
  {
    sizes_from_[i] = 0;
    if (procs_from_[i] != my_proc)
      MPI_Irecv((void *) &(sizes_from_[i]), 1, MPI_INT, procs_from_[i], tag_, comm_, &(request_[i]));
    else
    {
      sizes_from_[i] = sizes_to_[self_index_to];
      total_recv_length_ += sizes_from_[i];
    }
  }

  MPI_Barrier( comm_ );

  for (i = 0; i < (nsends_+self_msg_); i++)
  {
    if(procs_to_[i] != my_proc)
      MPI_Rsend ((void *) &(sizes_to_[i]), 1, MPI_INT, procs_to_[i], tag_, comm_);
  }

  if( nrecvs_ > 0 ) MPI_Waitall( nrecvs_, request_, status_ );

  for (i = 0; i < (nrecvs_+self_msg_); ++i)
  {
    if (procs_from_[i] != my_proc)
      total_recv_length_ += sizes_from_[i];
  }
#endif
  // end of exchanges sizes insert
 
  sum = 0;
  if( !starts_from_ptr_ ) starts_from_ptr_  = new int[nrecvs_+self_msg_]; 
  for (i = 0; i < (nrecvs_+self_msg_); ++i)
  {
     starts_from_ptr_[i] = sum;
     sum += sizes_from_[i];
  }

  resized_ = true;

  return 0;
}

//==============================================================================
int Epetra_MpiDistributor::Do( char * export_objs,
                               int obj_size,
                               int *& sizes,
                               int & len_import_objs,
                               char *& import_objs )
{
  EPETRA_CHK_ERR( DoPosts(export_objs, obj_size, sizes,
			  len_import_objs, import_objs) );
  EPETRA_CHK_ERR( DoWaits() );

  return(0);
}

//==============================================================================
int Epetra_MpiDistributor::DoReverse( char * export_objs,
                                      int obj_size,
                                      int *& sizes,
                                      int & len_import_objs,
                                      char *& import_objs )
{
  EPETRA_CHK_ERR( DoReversePosts(export_objs, obj_size, sizes,
				 len_import_objs, import_objs) );
  EPETRA_CHK_ERR( DoReverseWaits() );

  return(0);
}

//==============================================================================
int Epetra_MpiDistributor::DoPosts( char * export_objs,
                                    int obj_size,
                                    int *& sizes,
                                    int & len_import_objs,
                                    char *& import_objs )
{
  int ierr = Resize_(sizes);
  if (ierr != 0) {
    return(ierr);
  }

  MPI_Barrier( comm_ );

  int i, j, k;

  int my_proc = 0;
  int self_recv_address = 0;

  MPI_Comm_rank( comm_, &my_proc );

  if( len_import_objs < (total_recv_length_*obj_size) )
  {
    if( import_objs!=0 ) {delete [] import_objs; import_objs = 0;}
    len_import_objs = total_recv_length_*obj_size;
    if (len_import_objs>0) import_objs = new char[len_import_objs];
  }

  k = 0;

  for( i = 0; i < (nrecvs_+self_msg_); ++i )
  {
    if( procs_from_[i] != my_proc )
    {
      MPI_Irecv( &(import_objs[starts_from_ptr_[i] * obj_size]),
                 sizes_from_[i] * obj_size,
                 MPI_CHAR, procs_from_[i], tag_, comm_, &(request_[k]) );
      k++;
    }
    else
      self_recv_address = starts_from_ptr_[i] * obj_size;
  }

  MPI_Barrier( comm_ );

  //setup scan through procs_to list starting w/ higher numbered procs 
  //Should help balance msg traffic
  int nblocks = nsends_ + self_msg_; 
  int proc_index = 0; 
  while( proc_index < nblocks && procs_to_[proc_index] < my_proc )
    ++proc_index;                    
  if( proc_index == nblocks ) proc_index = 0;
   
  int self_num = 0;
  int p;

  if( !indices_to_ ) //data already blocked by processor
  {
    for( i = 0; i < nblocks; ++i )
    {
      p = i + proc_index;
      if( p > (nblocks-1) ) p -= nblocks;

      if( procs_to_[p] != my_proc )
        MPI_Rsend( &export_objs[starts_to_ptr_[p]*obj_size],
                   sizes_to_[p]*obj_size,
                   MPI_CHAR,
                   procs_to_[p],
                   tag_,
                   comm_ );
      else
        self_num = p;
    }

    if( self_msg_ )
      memcpy( &import_objs[self_recv_address],
              &export_objs[starts_to_ptr_[self_num]*obj_size],
              sizes_to_[self_num]*obj_size );
  }
  else //data not blocked by proc, need to copy to buffer
  {
    if( send_array_size_ && send_array_size_ < (max_send_length_*obj_size) )
    {
      if (send_array_!=0) {delete [] send_array_; send_array_ = 0;}
      send_array_ = 0;
      send_array_size_ = 0;
    }
    if( !send_array_size_ )
    {
      send_array_size_ = max_send_length_*obj_size;
      if (send_array_size_>0) send_array_ = new char[send_array_size_];
    }

    for( i=0; i<nblocks; ++i )
    {
      p = i + proc_index;
      if( p > (nblocks-1) ) p -= nblocks;

      if( procs_to_[p] != my_proc )
      {
        int offset = 0;
        j = starts_to_[p];
        for( k=0; k<lengths_to_[p]; ++k )
        {
          memcpy( &send_array_[offset],
                  &export_objs[indices_to_ptr_[j]*obj_size],
                  sizes_[indices_to_[j]]*obj_size );
          offset += sizes_[indices_to_[j]]*obj_size;
          ++j;
        }
        MPI_Rsend( send_array_, sizes_to_[p]*obj_size,
                   MPI_CHAR, procs_to_[p], tag_, comm_ );
      }
      else
        self_num = p;
    }

    if( self_msg_ )
    {
      int jj;
      j = starts_to_[self_num];
      for( k=0; k<lengths_to_[self_num]; ++k )
      {
        jj = indices_to_ptr_[j];
        memcpy( &import_objs[self_recv_address],
                &export_objs[jj*obj_size],
                sizes_[indices_to_[j]*obj_size] );
        self_recv_address += (obj_size*sizes_[indices_to_[j]]);
        ++jj;
      }
    }
  }

  return(0);
}

//==============================================================================
int Epetra_MpiDistributor::DoReversePosts( char * export_objs,
				           int obj_size,
                                           int *& sizes,
                                           int & len_import_objs,
				           char *& import_objs )
{
  assert(indices_to_==0); //Can only do reverse comm when original data
                          // is blocked by processor

  int i;
  int my_proc = 0;

  MPI_Comm_rank( comm_, &my_proc );

  if( comm_plan_reverse_ == 0 )
  {
    int total_send_length = 0;
    for( i = 0; i < nsends_+self_msg_; i++ )
      total_send_length += lengths_to_[i];

    int max_recv_length = 0;
    for( i = 0; i < nrecvs_; i++ )
      if( procs_from_[i] != my_proc )
        if( lengths_from_[i] > max_recv_length )
          max_recv_length = lengths_from_[i];

    comm_plan_reverse_ = new Epetra_MpiDistributor(*epComm_);

    comm_plan_reverse_->lengths_to_ = lengths_from_;
    comm_plan_reverse_->procs_to_ = procs_from_;
    comm_plan_reverse_->indices_to_ = indices_from_;
    comm_plan_reverse_->starts_to_ = starts_from_;

    comm_plan_reverse_->lengths_from_ = lengths_to_;
    comm_plan_reverse_->procs_from_ = procs_to_;
    comm_plan_reverse_->indices_from_ = indices_to_;
    comm_plan_reverse_->starts_from_ = starts_to_;

    comm_plan_reverse_->nsends_ = nrecvs_;
    comm_plan_reverse_->nrecvs_ = nsends_;
    comm_plan_reverse_->self_msg_ = self_msg_;

    comm_plan_reverse_->max_send_length_ = max_recv_length;
    comm_plan_reverse_->total_recv_length_ = total_send_length;

    comm_plan_reverse_->request_ = new MPI_Request[ comm_plan_reverse_->nrecvs_ ];
    comm_plan_reverse_->status_= new MPI_Status[ comm_plan_reverse_->nrecvs_ ];

    comm_plan_reverse_->no_delete_ = true;
  }

  int comm_flag = comm_plan_reverse_->DoPosts(export_objs, obj_size, sizes, len_import_objs, import_objs);

  return(comm_flag);
}

//==============================================================================
void Epetra_MpiDistributor::Print( ostream & os) const
{
  using std::endl;

  int myRank = 0, numProcs = 1;
  MPI_Comm_rank (comm_, &myRank);
  MPI_Comm_size (comm_, &numProcs);

  if (myRank == 0) {
    os << "Epetra_MpiDistributor (implements Epetra_Distributor)" << endl;
  }
  // Let each MPI process print its data.  We assume that all
  // processes can print to the given output stream, and execute
  // barriers to make it more likely that the output will be in the
  // right order.
  for (int p = 0; p < numProcs; ++p) {
    if (myRank == p) {
      os << "[Node " << p << " of " << numProcs << "]" << endl;
      os << " selfMessage: " << self_msg_ << endl;
      os << " numSends: " << nsends_ << endl;

      os << " imagesTo: [";
      for (int i = 0; i < nsends_; ++i) {
	os << procs_to_[i];
	if (i < nsends_ - 1) {
	  os << " ";
	}
      }
      os << "]" << endl;

      os << " lengthsTo: [";
      for (int i = 0; i < nsends_; ++i) {
	os << lengths_to_[i];
	if (i < nsends_ - 1) {
	  os << " ";
	}
      }
      os << "]" << endl;

      os << " maxSendLength: " << max_send_length_ << endl;

      os << " startsTo: ";
      if (starts_to_ == NULL) {
	os << "(NULL)" << endl;
      } else {
	os << "[";
	for (int i = 0; i < nsends_; ++i) {
	  os << starts_to_[i];
	  if (i < nsends_ - 1) {
	    os << " ";
	  }
	}
	os << "]" << endl;
      }

      os << " indicesTo: ";
      if (indices_to_ == NULL) {
	os << "(NULL)" << endl;
      } else {
	os << "[";
	int k = 0;
	for (int i = 0; i < nsends_; ++i) {
	  for (int j = 0; j < lengths_to_[i]; ++j) {
	    os << " " << indices_to_[j+k];
	  }
	  k += lengths_to_[i];
	}
	os << "]" << endl;
      }

      os << " numReceives: " << nrecvs_ << endl;
      os << " totalReceiveLength: " << total_recv_length_ << endl;

      os << " lengthsFrom: [";
      for (int i = 0; i < nrecvs_; ++i) {
	os << lengths_from_[i];
	if (i < nrecvs_ - 1) {
	  os << " ";
	}
      }
      os << "]" << endl;

      os << " startsFrom: [";
      for (int i = 0; i < nrecvs_; ++i) {
	os << starts_from_[i];
	if (i < nrecvs_ - 1) {
	  os << " ";
	}
      }
      os << "]" << endl;

      os << " imagesFrom: [";
      for (int i = 0; i < nrecvs_; ++i) {
	os << procs_from_[i];
	if (i < nrecvs_ - 1) {
	  os << " ";
	}
      }
      os << "]" << endl;

      // mfh 16 Dec 2011: I found this commented out here; not sure if
      // we want to print this, so I'm leaving it commented out.
      /*
	os << "indices_from: ";
	k = 0;
	for( i = 0; i < nrecvs_; i++ )
	{
	for( j = 0; j < lengths_from_[i]; j++ )
	os << " " << indices_from_[j+k];
	k += lengths_from_[i];
	}
      */

      // Last output is a flush; it leaves a space and also 
      // helps synchronize output.
      os << std::flush;
    } // if it's my process' turn to print

    // Execute barriers to give output time to synchronize.
    // One barrier generally isn't enough.
    MPI_Barrier (comm_);
    MPI_Barrier (comm_);
    MPI_Barrier (comm_);
  }
}

//---------------------------------------------------------------------------
int Epetra_MpiDistributor::Sort_ints_(
 int *vals_sort,     //  values to be sorted  
 int *vals_other,    // other array to be reordered with sort
 int  nvals)         // length of these two arrays
{
// It is primarily used to sort messages to improve communication flow. 
// This routine will also insure that the ordering produced by the invert_map
// routines is deterministic.  This should make bugs more reproducible.  This
// is accomplished by sorting the message lists by processor ID.
// This is a distribution count sort algorithm (see Knuth)
//  This version assumes non negative integers. 
   
    if (nvals <= 1) return 0;
        
    int i;                        // loop counter        
     
    // find largest int, n, to size sorting array, then allocate and clear it
    int n = 0;  
    for (i = 0; i < nvals; i++)  
       if (n < vals_sort[i]) n = vals_sort[i]; 
    int *pos = new int [n+2];  
    for (i = 0; i < n+2; i++) pos[i] = 0;
 
    // copy input arrays into temporary copies to allow sorting original arrays
    int *copy_sort  = new int [nvals]; 
    int *copy_other = new int [nvals]; 
    for (i = 0; i < nvals; i++)
    { 
      copy_sort[i]  = vals_sort[i]; 
      copy_other[i] = vals_other[i];  
    }                           
 
    // count the occurances of integers ("distribution count")
    int *p = pos+1;
    for (i = 0; i < nvals; i++) p[copy_sort[i]]++;
 
    // create the partial sum of distribution counts 
    for (i = 1; i < n; i++) p[i] += p[i-1]; 
 
    // the shifted partitial sum is the index to store the data  in sort order
    p = pos; 
    for (i = 0; i < nvals; i++)         
    {                                   
      vals_sort  [p[copy_sort [i]]]   = copy_sort[i];
      vals_other [p[copy_sort [i]]++] = copy_other[i]; 
    } 
 
    delete [] copy_sort;
    delete [] copy_other; 
    delete [] pos; 
 
    return 0;
}

//-------------------------------------------------------------------------
Epetra_MpiDistributor& Epetra_MpiDistributor::operator=(const Epetra_MpiDistributor& src)
{
  (void)src;
  //not currently supported
  bool throw_error = true;
  if (throw_error) {
    throw ReportError("Epetra_MpiDistributor::operator= not supported.",-1);
  }
  return( *this );
}
