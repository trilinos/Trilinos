#ifdef PETRA_MPI
#include "GSComm_Plan.h"

#include "GSComm_Comm.h"

GSComm_Plan::GSComm_Plan(const GSComm_Plan & Plan) 
: lengths_to_(0),
  procs_to_(0),
  indices_to_(0),
  size_indices_to_(Plan.size_indices_to_),
  lengths_from_(0),
  procs_from_(0),
  indices_from_(0),
  size_indices_from_(Plan.size_indices_from_),
  nrecvs_(Plan.nrecvs_),
  nsends_(Plan.nsends_),
  self_msg_(Plan.self_msg_),
  max_send_length_(Plan.max_send_length_),
  total_recv_length_(Plan.total_recv_length_),
  comm_(Plan.comm_),
  request_(0),
  status_(0),
  no_delete_(Plan.no_delete_)
{
  int i;
  if (nsends_>0) {
    lengths_to_ = new int[nsends_];
    procs_to_ = new int[nsends_];
    for (i=0; i<nsends_; i++) {
      lengths_to_[i] = Plan.lengths_to_[i];
      procs_to_[i] = Plan.procs_to_[i];
    }
  }
  if (size_indices_to_>0) {
    indices_to_ = new int[size_indices_to_];
    for (i=0; i<size_indices_to_; i++) {
      indices_to_[i] = Plan.indices_to_[i];
    }
  }

  if (nrecvs_>0) {
    lengths_from_ = new int[nrecvs_];
    procs_from_ = new int[nrecvs_];
    request_ = new MPI_Request[ nrecvs_ ];
    status_ = new MPI_Status[ nrecvs_ ];
    for (i=0; i<nrecvs_; i++) {
      lengths_from_[i] = Plan.lengths_from_[i];
      procs_from_[i] = Plan.procs_from_[i];
    }
  }
  if (size_indices_from_>0) {
    indices_from_ = new int[size_indices_from_];
    for (i=0; i<size_indices_from_; i++) {
      indices_from_[i] = Plan.indices_from_[i];
    }
  }
}



//---------------------------------------------------------------------------
//CreateFromSends Method
// - create communication plan given a known list of procs to send to
//---------------------------------------------------------------------------
bool GSComm_Plan::CreateFromSends( const int & nvals,
			           const int * assign,
			           MPI_Comm comm,
			           const int & tag,
			           const bool & deterministic,
			           int & pnrecv )
{
  int i;

  comm_ = comm;

  int my_proc;
  MPI_Comm_rank( comm, &my_proc );

  int nprocs;
  MPI_Comm_size( comm, &nprocs );

  // Check to see if items are grouped by processor w/o gaps
  // If so, indices_to -> 0

  // Setup data structures for quick traversal of arrays
  int * starts = new int[ nprocs + 1 ];
  for( i = 0; i < nprocs; i++ )
    starts[i] = 0;

  int nactive = 0;

  for( i = 0; i < nvals; i++ )
  {
    if( assign[i] >= 0 )
    {
      ++starts[ assign[i] ];
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

  for( i = 0; i < nvals; i++ )
    if( assign[i] >= 0 )
    {
      indices_to_[ starts[ assign[i] ] ] = i;
      ++starts[ assign[i] ];
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
  bool invert_flag = ComputeRecvs( my_proc, nprocs, tag, deterministic );

  if( invert_flag ) return invert_flag;

  total_recv_length_ = 0;
  for( i = 0; i < nrecvs_; i++ )
    total_recv_length_ += lengths_from_[i];

  if (nrecvs_>0) {
    request_ = new MPI_Request[ nrecvs_ ];
    status_ = new MPI_Status[ nrecvs_ ];
  }
  pnrecv = total_recv_length_;

  return false;
}

//---------------------------------------------------------------------------
//CreateFromRecvs Method
// - create communication plan given a known list of procs to recv from
//---------------------------------------------------------------------------
bool GSComm_Plan::CreateFromRecvs( const int & nvals,
				   const int * recv_gids,
			           const int * assign,
			           MPI_Comm comm,
			           const int & tag,
			           const bool & deterministic,
			           int & pnsends,
				   int *& send_gids,
				   int *& send_procs )
{
  bool comm_flag;

  int i;

  comm_ = comm;

  int my_proc;
  MPI_Comm_rank( comm, &my_proc );

  int nprocs;
  MPI_Comm_size( comm, &nprocs );

  comm_flag = ComputeSends( nvals, recv_gids, assign, pnsends,
	send_gids, send_procs, my_proc );

  // Setup data structures for quick traversal of arrays
  int * starts = new int[ nprocs + 1 ];
  for( i = 0; i < nprocs; i++ )
    starts[i] = 0;

  int nactive = 0;

  for( i = 0; i < pnsends; i++ )
  {
    if( send_procs[i] >= 0 )
    {
      ++starts[ send_procs[i] ];
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

  for( i = 0; i < pnsends; i++ )
    if( send_procs[i] >= 0 )
    {
      indices_to_[ starts[ send_procs[i] ] ] = i;
      ++starts[ send_procs[i] ];
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
  bool invert_flag = ComputeRecvs( my_proc, nprocs, tag, deterministic );

  if( invert_flag ) return invert_flag;

  total_recv_length_ = 0;
  for( i = 0; i < ( nrecvs_ + self_msg_ ); i++ )
    total_recv_length_ += lengths_from_[i];

  if (nrecvs_>0) {
    request_ = new MPI_Request[ nrecvs_ ];
    status_ = new MPI_Status[ nrecvs_ ];
  }

  return comm_flag;
}

//---------------------------------------------------------------------------
//ComputeRecvs Method
//---------------------------------------------------------------------------
bool GSComm_Plan::ComputeRecvs( const int & my_proc, 
			        const int & nprocs,
			        const int & tag, 
			        const bool & deterministic )
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
      MPI_Send( &(lengths_to_[i]), 1, MPI_INT, procs_to_[i], tag, comm_ );
    }
    else
    {
      assert(nrecvs_>0);
      lengths_from_[nrecvs_-1] = lengths_to_[i];
      procs_from_[nrecvs_-1] = my_proc;
    }

  for( i = 0; i < ( nrecvs_ - self_msg_ ); i++ )
  {
    MPI_Recv( &(lengths_from_[i]), 1, MPI_INT, MPI_ANY_SOURCE, tag, comm_, &status );
    procs_from_[i] = status.MPI_SOURCE;
  }

  MPI_Barrier( comm_ );

  if( deterministic )
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

//---------------------------------------------------------------------------
//ComputeSends Method
//---------------------------------------------------------------------------
bool GSComm_Plan::ComputeSends( const int & num_imports,
				const int * import_ids,
				const int * import_procs,
				int & num_exports,
				int *& export_ids,
				int *& export_procs,
				const int & my_proc )
{
  GSComm_Plan comm_plan;
  GSComm_Comm comm_doer;

  int i;

  int msgtag = 32767;
  int msgtag2 = 32766;

  int * proc_list = 0;
  int * import_objs = 0;
  int * export_objs = 0;

  bool comm_flag;

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

  comm_flag = comm_plan.CreateFromSends( num_imports, proc_list,
		comm_, msgtag, true, num_exports ); 

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

  comm_doer.Do( comm_plan, msgtag2, 
	reinterpret_cast<char *> (import_objs), 
	2 * sizeof( int ), 
	reinterpret_cast<char *> (export_objs) );


  for( i = 0; i < num_exports; i++ )
  {
    export_ids[i] = export_objs[2*i];
    export_procs[i] = export_objs[2*i+1];
  }

  if( proc_list != 0 ) delete [] proc_list;
  if( import_objs != 0 ) delete [] import_objs;
  if( export_objs != 0 ) delete [] export_objs;

  return comm_flag;

}

ostream & operator<< ( ostream & os, const GSComm_Plan & plan )
{
  int i, j;
  os << "GSComm_Plan Object" << endl;
  os << "------------------" << endl;
  os << "nsends: " << plan.nsends_ << endl;
  os << "procs_to: ";
  for( i = 0; i < plan.nsends_; i++ )
    os << " " << plan.procs_to_[i];
  os << endl;
  os<< "lengths_to: ";
  for( i = 0; i < plan.nsends_; i++ )
    os << " " << plan.lengths_to_[i];
  os << endl;
  os << "indices_to: ";
  int k = 0;
  for( i = 0; i < plan.nsends_; i++ )
  {
    for( j = 0; j < plan.lengths_to_[i]; j++ )
      os << " " << plan.indices_to_[j+k];
    k += plan.lengths_to_[i];
  }
  os << endl;
  os << "nrecvs: " << plan.nrecvs_ << endl;
  os << "procs_from: ";
  for( i = 0; i < plan.nrecvs_; i++ )
    os << " " << plan.procs_from_[i];
  os << endl;
  os << "lengths_from: ";
  for( i = 0; i < plan.nrecvs_; i++ )
    os << " " << plan.lengths_from_[i];
  os << endl;
/*
  os << "indices_from: ";
  k = 0;
  for( i = 0; i < plan.nrecvs_; i++ )
  {
    for( j = 0; j < plan.lengths_from_[i]; j++ )
      os << " " << plan.indices_from_[j+k];
    k += plan.lengths_from_[i];
  }
*/
  os << "self_msg: " << plan.self_msg_ << endl;
  os << "max_send_length: " << plan.max_send_length_ << endl;
  os << "total_recv_length: " << plan.total_recv_length_ << endl;
  os << endl;

  return os;
}
#endif /* PETRA_MPI */
