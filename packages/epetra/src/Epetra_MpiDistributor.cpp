//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_MpiDistributor.h"
#include "Epetra_MpiComm.h"


//==============================================================================
// Epetra_MpiDistributor constructor
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
// Epetra_MpiDistributor destructor
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
//---------------------------------------------------------------------------
//CreateFromSends Method
// - create communication plan given a known list of procs to send to
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::CreateFromSends( const int & NumExportIDs,
                                            const int * ExportPIDs,
                                            bool Deterministic,
                                            int & NumRemoteIDs )
{
  nexports_ = NumExportIDs;

  int i;

  int my_proc;
  MPI_Comm_rank( comm_, &my_proc );

  int nprocs;
  MPI_Comm_size( comm_, &nprocs );

  // Check to see if items are grouped by processor w/o gaps
  // If so, indices_to -> 0

  // Setup data structures for quick traversal of arrays
  int * starts = new int[ nprocs + 1 ];
  for( i = 0; i < nprocs; i++ )
    starts[i] = 0;

  int nactive = 0;
  bool no_send_buff = true;

  for( i = 0; i < NumExportIDs; i++ )
  {
    if( no_send_buff && i && (ExportPIDs[i] < ExportPIDs[i-1]) )
      no_send_buff = false;
    if( ExportPIDs[i] >= 0 )
    {
      ++starts[ ExportPIDs[i] ];
      ++nactive;
    }
  }

  self_msg_ = ( starts[my_proc] != 0 );

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

    int index = 0;
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

  //Invert map to see what msgs are received and what length
  EPETRA_CHK_ERR( ComputeRecvs_( my_proc, nprocs ) );

  if (nrecvs_>0) {
    request_ = new MPI_Request[ nrecvs_ ];
    status_ = new MPI_Status[ nrecvs_ ];
  }

  NumRemoteIDs = total_recv_length_;

  return 0;
}

//==============================================================================
//---------------------------------------------------------------------------
//CreateFromRecvs Method
// - create communication plan given a known list of procs to recv from
//---------------------------------------------------------------------------
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

  MPI_Reduce_scatter( msg_count, &nrecvs_, counts, MPI_INT, MPI_SUM, comm_ );

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

  Sort_ints_( procs_from_, lengths_from_, nrecvs_ );

  // Compute indices_from_
  size_indices_from_ = 0;
  if( nrecvs_ > 0 )
  {
    for( i = 0; i < nrecvs_; i++ )  size_indices_from_ += lengths_from_[i];
    indices_from_ = new int[ size_indices_from_ ];

    for (i=0; i<size_indices_from_; i++) indices_from_[i] = i;
  }

  starts_from_ = new int[nrecvs_];
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
int Epetra_MpiDistributor::ComputeSends_( int num_imports,
				const int *& import_ids,
				const int *& import_procs,
				int & num_exports,
				int *& export_ids,
				int *& export_procs,
				int my_proc ) {
 
  Epetra_MpiDistributor tmp_plan(*epComm_);
  int i;

  int * proc_list = 0;
  int * import_objs = 0;
  char * c_export_objs = 0;

  if( num_imports > 0 )
  {
    proc_list = new int[ num_imports ];
    import_objs = new int[ 2 * num_imports ];

    for( i = 0; i < num_imports; i++ )
    {
      proc_list[i] = import_procs[i];

      import_objs[2*i] = import_ids[i];
      import_objs[2*i+1] = my_proc;
    }
  }

  EPETRA_CHK_ERR(tmp_plan.CreateFromSends( num_imports, proc_list,
					   true, num_exports) );
  if( num_exports > 0 )
  {
    //export_objs = new int[ 2 * num_exports ];
    export_ids = new int[ num_exports ];
    export_procs = new int[ num_exports ];
  }
  else
  {
    export_ids = 0;
    export_procs = 0;
  }

  int len_c_export_objs = 0;
  EPETRA_CHK_ERR( tmp_plan.Do(reinterpret_cast<char *> (import_objs),
			      2 * sizeof( int ), 
			      len_c_export_objs,
			      c_export_objs) );
  int * export_objs = reinterpret_cast<int *>(c_export_objs);

  for( i = 0; i < num_exports; i++ ) {
    export_ids[i] = export_objs[2*i];
    export_procs[i] = export_objs[2*i+1];
  }

  if( proc_list != 0 ) delete [] proc_list;
  if( import_objs != 0 ) delete [] import_objs;
  if( len_c_export_objs != 0 ) delete [] c_export_objs;

  return(0);

}

//==============================================================================
// Do method
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
// DoReverse method
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
//---------------------------------------------------------------------------
//Do_Posts Method
//---------------------------------------------------------------------------
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
    if( len_import_objs ) delete [] import_objs;
    len_import_objs = total_recv_length_*obj_size;
    import_objs = new char[len_import_objs];
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

  MPI_Barrier( comm_ );

  //setup scan through procs_to list starting w/ higher numbered procs 
  //Should help balance msg traffic
  int nblocks = nsends_ + self_msg_; 
  int proc_index = 0; 
  while( proc_index < nblocks && procs_to_[proc_index] < my_proc )
    ++proc_index;                    
  if( proc_index == nblocks ) proc_index = 0;
   
  int self_num, self_index;
  int p;

  if( !indices_to_ ) //data already blocked by processor
  {
    for( i = 0; i < nblocks; ++i )
    {
      p = i + proc_index;
      if( p > (nblocks-1) ) p -= nblocks;

      if( procs_to_[p] != my_proc )
        MPI_Rsend( &export_objs[starts_to_[p]*obj_size],
                   lengths_to_[p]*obj_size,
                   MPI_CHAR,
                   procs_to_[p],
                   tag_,
                   comm_ );
      else
       self_num = p;
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
      if( send_array_size_ ) delete [] send_array_;
      send_array_size_ = max_send_length_*obj_size;
      send_array_ = new char[send_array_size_];
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
        MPI_Rsend( send_array_,
                   lengths_to_[p] * obj_size,
                   MPI_CHAR,
                   procs_to_[p],
                   tag_, comm_ );
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
//---------------------------------------------------------------------------
//Do_Waits Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::DoWaits()
{
  if( nrecvs_ > 0 ) MPI_Waitall( nrecvs_, request_, status_ );

  return(0);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Posts Method
//---------------------------------------------------------------------------
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
//---------------------------------------------------------------------------
//DoReverse_Waits Method
//---------------------------------------------------------------------------
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
    int * index = new int[nsends_+self_msg_];
    int * sort_val = new int[nsends_+self_msg_];

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

    int sum = 0;
    for( i = 0; i < (nsends_+self_msg_); ++i )
    {
      starts_to_ptr_[ index[i] ] = sum;
      sum += sizes_to_[ index[i] ];
    }

    delete [] index;
    delete [] sort_val;
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
 
    delete [] offset;
  }
 
  //  Exchange sizes routine inserted here:
  int self_index_to = -1;
  for (i = 0; i < (nsends_+self_msg_); i++)
  {
    if(procs_to_[i] != my_proc)
      MPI_Send ((void *) &(sizes_to_[i]), 1, MPI_INT, procs_to_[i], tag_, comm_);
    else
      self_index_to = i;
  }

  total_recv_length_ = 0;
  MPI_Status status;
  if( !sizes_from_ && (nrecvs_+self_msg_) ) sizes_from_ = new int [nrecvs_+self_msg_];
  for (i = 0; i < (nrecvs_+self_msg_); ++i)
  {
    sizes_from_[i] = 0;
    if (procs_from_[i] != my_proc)
      MPI_Recv((void *) &(sizes_from_[i]), 1, MPI_INT, procs_from_[i], tag_, comm_, &status);
    else
      sizes_from_[i] = sizes_to_[self_index_to];
    total_recv_length_ += sizes_from_[i];
  }
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
// GSComm_Comm Do method
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
// GSComm_Comm DoReverse method
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
//---------------------------------------------------------------------------
//Do_Posts Method (Variable Block Size)
//---------------------------------------------------------------------------
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
    if( len_import_objs ) delete [] import_objs;
    len_import_objs = total_recv_length_*obj_size;
    import_objs = new char[len_import_objs];
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
   
  int self_num;
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
      delete [] send_array_;
      send_array_size_ = 0;
    }
    if( !send_array_size_ )
    {
      send_array_size_ = max_send_length_*obj_size;
      send_array_ = new char[send_array_size_];
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
//---------------------------------------------------------------------------
//DoReverse_Posts Method
//---------------------------------------------------------------------------
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
  int i, j;
  os << "nsends: " << nsends_ << endl;
  os << "procs_to: ";
  for( i = 0; i < nsends_; i++ )
    os << " " << procs_to_[i];
  os << endl;
  os<< "lengths_to: ";
  for( i = 0; i < nsends_; i++ )
    os << " " << lengths_to_[i];
  os << endl;
  os << "indices_to: ";
  int k = 0;
  if( indices_to_ )
  {
    for( i = 0; i < nsends_; i++ )
    {
      for( j = 0; j < lengths_to_[i]; j++ )
        os << " " << indices_to_[j+k];
      k += lengths_to_[i];
    }
  }
  os << endl;
  os << "nrecvs: " << nrecvs_ << endl;
  os << "procs_from: ";
  for( i = 0; i < nrecvs_; i++ )
    os << " " << procs_from_[i];
  os << endl;
  os << "lengths_from: ";
  for( i = 0; i < nrecvs_; i++ )
    os << " " << lengths_from_[i];
  os << endl;
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
  os << "self_msg: " << self_msg_ << endl;
  os << "max_send_length: " << max_send_length_ << endl;
  os << "total_recv_length: " << total_recv_length_ << endl;
  os << endl;

  return;
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
