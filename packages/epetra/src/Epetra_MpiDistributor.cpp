
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
  nrecvs_(0),
  nsends_(0),
  self_msg_(0),
  max_send_length_(0),
  total_recv_length_(0),
  tag_(Comm.GetMpiTag()),
  epComm_(&Comm),
  comm_(Comm.GetMpiComm()),
  request_(0),
  status_(0),
  no_delete_(false),
  recv_array_(0),
  send_array_(0),
  comm_plan_reverse_(0){
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
  nrecvs_(Distributor.nrecvs_),
  nsends_(Distributor.nsends_),
  self_msg_(Distributor.self_msg_),
  max_send_length_(Distributor.max_send_length_),
  total_recv_length_(Distributor.total_recv_length_),
  tag_(Distributor.tag_),
  epComm_(Distributor.epComm_),
  comm_(Distributor.comm_),
  request_(0),
  status_(0),
  no_delete_(Distributor.no_delete_),
  recv_array_(0),
  send_array_(0),
  comm_plan_reverse_(0){
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
Epetra_MpiDistributor::~Epetra_MpiDistributor() {
  if( !no_delete_ )
  {
    if( lengths_to_ != 0 ) delete [] lengths_to_;
    if( procs_to_ != 0 ) delete [] procs_to_;
    if( indices_to_ != 0 ) delete [] indices_to_;
    if( lengths_from_ != 0 ) delete [] lengths_from_;
    if( procs_from_ != 0 ) delete [] procs_from_;
    if( indices_from_ != 0 ) delete [] indices_from_;
  }

  if( request_ != 0 ) delete [] request_;
  if( status_ != 0 ) delete [] status_;

  if( send_array_ != 0 ) delete [] send_array_;

  if( comm_plan_reverse_ != 0 ) delete comm_plan_reverse_;
}


//==============================================================================
//---------------------------------------------------------------------------
//CreateFromSends Method
// - create communication plan given a known list of procs to send to
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::CreateFromSends( const int & NumExportIDs,
			           const int * ExportPIDs,
			           const bool & Deterministic,
			           int & NumRemoteIDs ) {
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

  for( i = 0; i < NumExportIDs; i++ )
  {
    if( ExportPIDs[i] >= 0 )
    {
      ++starts[ ExportPIDs[i] ];
      ++nactive;
    }
  }

  self_msg_ = ( starts[my_proc] != 0 );

  nsends_ = 0;
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
  }

  int j = 0;
  max_send_length_ = 0;

  for( i = 0; i < nprocs; i++ )
    if( starts[i+1] != starts[i] )
    {
      lengths_to_[j] = starts[i+1] - starts[i];
      if( ( i != my_proc ) && ( lengths_to_[j] > max_send_length_ ) )
        max_send_length_ = lengths_to_[j];
      procs_to_[j] = i;
      j++;
    }
    
  delete [] starts;

  //Invert map to see what msgs are received and what length
  EPETRA_CHK_ERR(ComputeRecvs( my_proc, nprocs, Deterministic ));

  total_recv_length_ = 0;
  for( i = 0; i < nrecvs_; i++ )
    total_recv_length_ += lengths_from_[i];

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
			           const bool & Deterministic,
			           int & NumExportIDs,
				   int *& ExportGIDs,
				   int *& ExportPIDs )
{
  int i;

  int my_proc;
  MPI_Comm_rank( comm_, &my_proc );

  int nprocs;
  MPI_Comm_size( comm_, &nprocs );

  EPETRA_CHK_ERR(ComputeSends( NumRemoteIDs, RemoteGIDs, RemotePIDs, NumExportIDs,
			    ExportGIDs, ExportPIDs, my_proc));

  // Setup data structures for quick traversal of arrays
  int * starts = new int[ nprocs + 1 ];
  for( i = 0; i < nprocs; i++ )
    starts[i] = 0;

  int nactive = 0;

  for( i = 0; i < NumExportIDs; i++ )
  {
    if( ExportPIDs[i] >= 0 )
    {
      ++starts[ ExportPIDs[i] ];
      ++nactive;
    }
  }

  self_msg_ = ( starts[my_proc] != 0 );

  nsends_ = 0;
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
  }

  int j = 0;
  max_send_length_ = 0;

  for( i = 0; i < nprocs; i++ )
    if( starts[i+1] != starts[i] )
    {
      lengths_to_[j] = starts[i+1] - starts[i];
      if( ( i != my_proc ) && ( lengths_to_[j] > max_send_length_ ) )
        max_send_length_ = lengths_to_[j];
      procs_to_[j] = i;
      j++;
    }
    
  delete [] starts;

  //Invert map to see what msgs are received and what length
  EPETRA_CHK_ERR(ComputeRecvs( my_proc, nprocs, Deterministic));

  total_recv_length_ = 0;
  for( i = 0; i < ( nrecvs_ + self_msg_ ); i++ )
    total_recv_length_ += lengths_from_[i];

  if (nrecvs_>0) {
    request_ = new MPI_Request[ nrecvs_ ];
    status_ = new MPI_Status[ nrecvs_ ];
  }

  return(0);
}

//==============================================================================
//---------------------------------------------------------------------------
//ComputeRecvs Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::ComputeRecvs( const int & my_proc, 
			        const int & nprocs, 
			        const bool & Deterministic )
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

  for( i = 0; i < nsends_; i++ )
    msg_count[ procs_to_[i] ] = 1;

  MPI_Reduce_scatter( msg_count, &nrecvs_, counts, MPI_INT, MPI_SUM, comm_ );

  delete [] msg_count;
  delete [] counts;

  if (nrecvs_>0) {
    lengths_from_ = new int[ nrecvs_ ];
    procs_from_ = new int[ nrecvs_ ];
  }
  for( i = 0; i < nsends_; i++ )
    if( procs_to_[i] != my_proc ) {
      MPI_Send( &(lengths_to_[i]), 1, MPI_INT, procs_to_[i], tag_, comm_ );
    }
    else
    {
      assert(nrecvs_>0);
      lengths_from_[nrecvs_-1] = lengths_to_[i];
      procs_from_[nrecvs_-1] = my_proc;
    }

  for( i = 0; i < ( nrecvs_ - self_msg_ ); i++ )
  {
    MPI_Recv( &(lengths_from_[i]), 1, MPI_INT, MPI_ANY_SOURCE, tag_, comm_, &status );
    procs_from_[i] = status.MPI_SOURCE;
  }

  MPI_Barrier( comm_ );

  if( Deterministic )
  {
    int j;
    int temp;

    for( i = 1; i < ( nrecvs_ - self_msg_ ); i++ )
    {
      j = i;

      while( ( j > 0 ) && ( procs_from_[j] < procs_from_[j-1] ) )
      {
        temp = procs_from_[j];
        procs_from_[j] = procs_from_[j-1];
        procs_from_[j-1] = temp;

        temp = lengths_from_[j];
        lengths_from_[j] = lengths_from_[j-1];
        lengths_from_[j-1] = temp;

        j--;
      }
    }
  }

  // Compute indices_from_

  size_indices_from_ = 0;
  for( i = 0; i < nrecvs_; i++ )  size_indices_from_ += lengths_from_[i];
  indices_from_ = new int[ size_indices_from_ ];

  for (i=0; i<size_indices_from_; i++) indices_from_[i] = i;

  MPI_Barrier( comm_ );

  
  return false;

}

//==============================================================================
//---------------------------------------------------------------------------
//ComputeSends Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::ComputeSends( const int & num_imports,
				const int * import_ids,
				const int * import_procs,
				int & num_exports,
				int *& export_ids,
				int *& export_procs,
				const int & my_proc ) {
 
 Epetra_MpiDistributor tmp_plan(*epComm_);
  int i;

  int * proc_list = 0;
  int * import_objs = 0;
  int * export_objs = 0;

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

  EPETRA_CHK_ERR(tmp_plan.CreateFromSends( num_imports, proc_list, true, num_exports)); 
  if( num_exports > 0 )
  {
    export_objs = new int[ 2 * num_exports ];
    export_ids = new int[ num_exports ]; // Note: export_ids and export_procs must be deleted by the calling routine
    export_procs = new int[ num_exports ];
  }
  else
  {
    export_ids = 0;
    export_procs = 0;
  }

  EPETRA_CHK_ERR(tmp_plan.Do(reinterpret_cast<char *> (import_objs), 
			   2 * sizeof( int ), 
			   reinterpret_cast<char *> (export_objs)));


  for( i = 0; i < num_exports; i++ ) {
    export_ids[i] = export_objs[2*i];
    export_procs[i] = export_objs[2*i+1];
  }

  if( proc_list != 0 ) delete [] proc_list;
  if( import_objs != 0 ) delete [] import_objs;
  if( export_objs != 0 ) delete [] export_objs;

  return(0);

}

//==============================================================================
// Do method
int Epetra_MpiDistributor::Do(char * export_objs, const int & obj_size, char * import_objs )
{
  EPETRA_CHK_ERR(DoPosts(export_objs, obj_size, import_objs));
  EPETRA_CHK_ERR(DoWaits(export_objs, obj_size, import_objs));
 return(0);
}

//==============================================================================
// DoReverse method
int Epetra_MpiDistributor::DoReverse(char * export_objs,const int & obj_size, char * import_objs )
{

  EPETRA_CHK_ERR(DoReversePosts(export_objs, obj_size, import_objs));
  EPETRA_CHK_ERR(DoReverseWaits(export_objs, obj_size, import_objs));
  return(0);
}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Posts Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::DoPosts(char * export_objs,
				const int & obj_size,
				char * import_objs ) {
  int i, j, k;

  int my_proc = 0;
  int self_recv_address = 0;

  MPI_Comm_rank( comm_, &my_proc );

  recv_array_ = import_objs;

  j = 0;
  k = 0;

  for( i = 0; i < nrecvs_; i++ )
  {
    if( procs_from_[i] != my_proc )
    {
      //cout << "Processor " << my_proc << "lengths_from_["<<i<<"] = " << lengths_from_[i] << " obj_size = " <<obj_size<<endl;
      MPI_Irecv( &(recv_array_[j]), lengths_from_[i] * obj_size,
	MPI_CHAR, procs_from_[i], tag_, comm_,
	&(request_[k]) );
      k++;
    }
    else
      self_recv_address = j;

    j += lengths_from_[i] * obj_size;
  }

  MPI_Barrier( comm_ );

  int self_num, self_index;

  if (max_send_length_* obj_size>0) 
    send_array_ = new char[ max_send_length_ * obj_size ];

  j = 0;
  for( i = 0; i < nsends_; i++ )
  {
    if( procs_to_[i] != my_proc )
    {
      int offset = 0;
      for( k = 0; k < lengths_to_[i]; k++ )
      {
        memcpy( &(send_array_[offset]), 
 	  &(export_objs[indices_to_[j]*obj_size]), obj_size );
        j++;
        offset += obj_size;
      }
      //   cout << "my_proc = " << my_proc << " length = " << lengths_to_[i] * obj_size 
      //   << " send to = " << procs_to_[i] << " tag = " << tag << endl;
      //cout << "Processor " << my_proc << "lengths_to_["<<i<<"] = " << lengths_to_[i] << " obj_size = " <<obj_size<<endl;
      MPI_Rsend( send_array_, lengths_to_[i] * obj_size,
	MPI_CHAR, procs_to_[i], tag_, comm_ );
    }
    else
    {
      self_num = i;
      self_index = j;
      j += lengths_to_[i];
    }
  }

  if( self_msg_ )
    for( k = 0; k < lengths_to_[self_num]; k++ )
    {
      memcpy( &(recv_array_[self_recv_address]),
	 &(export_objs[indices_to_[self_index]*obj_size]),
	  obj_size );
      self_index++;
      self_recv_address += obj_size;
    }

  if (send_array_!=0) {
    delete [] send_array_;
    send_array_ = 0;
  }
  return(0);
}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Waits Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::DoWaits(char * export_objs,
			       const int & obj_size,
			       char * import_objs )
{
  if( nrecvs_ - self_msg_ > 0 )
    MPI_Waitall( nrecvs_ - self_msg_, 
		 request_, status_ );

  return(0);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Posts Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::DoReversePosts(char * export_objs,
				       const int & obj_size,
				       char * import_objs )
{
  int i;
  int my_proc = 0;

  MPI_Comm_rank( comm_, &my_proc );

  int total_send_length = 0;
  for( i = 0; i < nsends_; i++ )
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
  comm_plan_reverse_->lengths_from_ = lengths_to_;
  comm_plan_reverse_->procs_from_ = procs_to_;
  comm_plan_reverse_->indices_from_ = indices_to_;
  comm_plan_reverse_->nrecvs_ = nsends_;
  comm_plan_reverse_->nsends_ = nrecvs_;
  comm_plan_reverse_->self_msg_ = self_msg_;
  comm_plan_reverse_->max_send_length_ = max_recv_length;
  comm_plan_reverse_->total_recv_length_ = total_send_length;

  comm_plan_reverse_->request_ = new MPI_Request[ comm_plan_reverse_->nrecvs_ ];
  comm_plan_reverse_->status_= new MPI_Status[ comm_plan_reverse_->nrecvs_ ];

  comm_plan_reverse_->no_delete_ = true;

  int comm_flag = comm_plan_reverse_->DoPosts(export_objs, obj_size, import_objs);

  return(comm_flag);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Waits Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::DoReverseWaits(char * export_objs,
		 	           const int & obj_size,
			           char * import_objs )
{
  if( comm_plan_reverse_ == 0 ) return (-1);

  int comm_flag = comm_plan_reverse_->DoWaits(export_objs, obj_size, import_objs );

  if (comm_plan_reverse_!=0) {
    delete comm_plan_reverse_;
    comm_plan_reverse_ = 0;
  }

  return(comm_flag);
}

//==============================================================================
// GSComm_Comm Do method
int Epetra_MpiDistributor::Do(char * export_objs, const int * & obj_size, char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
// GSComm_Comm DoReverse method
int Epetra_MpiDistributor::DoReverse(char * export_objs,const int * & obj_size, char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Posts Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::DoPosts(char * export_objs,
				const int * & obj_size,
				char * import_objs ) {
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);

}
//==============================================================================
//---------------------------------------------------------------------------
//Do_Waits Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::DoWaits(char * export_objs,
			       const int * & obj_size,
			       char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Posts Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::DoReversePosts(char * export_objs,
				       const int * & obj_size,
				       char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
}

//==============================================================================
//---------------------------------------------------------------------------
//DoReverse_Waits Method
//---------------------------------------------------------------------------
int Epetra_MpiDistributor::DoReverseWaits(char * export_objs,
		 	           const int * & obj_size,
			           char * import_objs )
{
  EPETRA_CHK_ERR(-1); // This method should never be called 
  return(-1);
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
  for( i = 0; i < nsends_; i++ )
  {
    for( j = 0; j < lengths_to_[i]; j++ )
      os << " " << indices_to_[j+k];
    k += lengths_to_[i];
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

