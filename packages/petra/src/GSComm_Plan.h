#ifdef PETRA_MPI
#ifndef _GSCOMM_PLAN_H_
#define _GSCOMM_PLAN_H_

//
// This plan is constructed for gather/scatter operations
// OO version derived from Zoltan's LB_comm_create
//

#include "Petra_Petra.h"

class GSComm_Plan {
    
  public:

    GSComm_Plan();

    GSComm_Plan(const GSComm_Plan & Plan);

    ~GSComm_Plan();

    bool CreateFromSends( const int & nvals,
	    	          const int * assign,
	    	          MPI_Comm comm,
	                  const int & tag,
                          const bool & deterministic,
	                  int & pnrecv );

    bool CreateFromRecvs( const int & nvals,
			  const int * recv_gids,
	    	          const int * assign,
	    	          MPI_Comm comm,
	                  const int & tag,
                          const bool & deterministic,
	                  int & pnsend,
			  int *& send_gids,
			  int *& send_procs );

    void Clear();

  private:

    bool ComputeRecvs( const int & my_proc,
	               const int & nprocs,
	               const int & tag,
	               const bool & deterministic );

    bool ComputeSends( const int & num_imports,
		       const int * import_ids,
		       const int * import_procs,
		       int & num_exports,
		       int *& export_ids,
		       int *& export_procs,
		       const int & my_proc );

  private:

    int * lengths_to_;
    int * procs_to_;
    int * indices_to_;
    int  size_indices_to_;
    int * lengths_from_;
    int * procs_from_;
    int * indices_from_;
    int  size_indices_from_;
    int   nrecvs_;
    int   nsends_;
    int   self_msg_;
    int   max_send_length_;
    int   total_recv_length_;

    MPI_Comm      comm_;
    MPI_Request * request_;
    MPI_Status  * status_;

    bool no_delete_;

    friend class GSComm_Comm;

    friend ostream & operator<<( ostream & os, const GSComm_Plan & plan );

};

// GSComm_Plan constructor
inline GSComm_Plan::GSComm_Plan()
	: lengths_to_(0),
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
	  comm_(MPI_COMM_WORLD),
	  request_(0),
          status_(0),
	  no_delete_(false)
{
}

// GSComm_Plan destructor
inline GSComm_Plan::~GSComm_Plan()
{
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
}

// Clear method to clear attributes.
inline void GSComm_Plan::Clear()
{
  if( !no_delete_ )
  {
    if( lengths_to_ != 0 )   { delete [] lengths_to_;   lengths_to_ = 0; }
    if( procs_to_ != 0 )     { delete [] procs_to_;     procs_to_  = 0; }
    if( indices_to_ != 0 )   { delete [] indices_to_;   indices_to_  = 0; }
    if( lengths_from_ != 0 ) { delete [] lengths_from_; lengths_from_ = 0; }
    if( procs_from_ != 0 )   { delete [] procs_from_;   procs_from_  = 0; }
    if( indices_from_ != 0 ) { delete [] indices_from_; indices_from_ = 0; }
  }

  if( request_ != 0 ) { delete [] request_; request_ = 0; }
  if( status_ != 0 )  { delete [] status_;  status_  = 0; }

  nrecvs_ = 0;
  nsends_ = 0;
  self_msg_ = 0;
  max_send_length_ = 0;
  total_recv_length_ = 0;
}

#endif /* _GSCOMM_PLAN_H_ */
#endif /* PETRA_MPI */
