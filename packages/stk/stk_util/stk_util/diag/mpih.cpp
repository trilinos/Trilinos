/**   ------------------------------------------------------------
 *    Copyright 2002 - 2009 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

/** \file
    Purpose:

    These functions provide a C++ framework interface to the
    MPIH collective communication functions.

    Written by James R. Overfelt
    jroverf@sandia.gov
    Oct 2002
*/

#include <stk_util/diag/mpih.hpp>
#include <MPIH_Include.h>               // for MPIH_Has_Exception_Flag, etc
#include <map>                          // for map<>::mapped_type
#include <ostream>                      // for endl, operator<<, etc
#include <stk_util/diag/Env.hpp>        // for parallel_comm, output
#include <stk_util/diag/ExceptionReport.hpp>  // for RuntimeWarning
#include <stk_util/diag/Trace.hpp>      // for Traceback
#include <stk_util/environment/ProductRegistry.hpp>
#include <string>                       // for basic_string
#include "mpi.h"                        // for MPI_SUCCESS, MPI_Datatype, etc
#include "stk_util/diag/Exception.hpp"  // for ExTemp1, RuntimeError, etc
#include "stk_util/environment/ReportHandler.hpp"  // for StackTrace
#include "stk_util/environment/RuntimeWarning.hpp"

using sierra::Diag::Traceback;

extern "C" {

  /**  Handler compete function needed by MPIH.
   *
   *  This particular Compete function does nothing.   The compete
   *  function is used in two places.  One is when different
   *  processors throw different exceptions.  Then the compete
   *  function is called to resolve the conflict about which
   *  exception to throw on all processors in parallel.
   *  The second use for the compete function is when
   *  Set_Exception is called multiple times on a single
   *  processor, Compete will be called to resolve which.
   *  exception is kept.
   *
   *  The Execute handler is called on each exception during
   *  the parallel exception handling.  After all of the exceptions
   *  have been distributed to all of the processors, the exception
   *  list is iterated and Execute called on each non-NULL
   *  exception.  The Execute function is called in parallel
   *  so collective communication calls are valid in an Execute
   *  handle.
   *  Besides the Execute handle, there is a call back function
   *  registered in the standard Framework exception class called
   *  Parallel_Handler that will also be called in parallel
   *  on each exception, allowing class specific handling.
   */
  static void Handler_Compete (void **h_exists, void *h_in) throw () {
    /* %TRACE[TRACEBACK]% */ Traceback trace__("Handler_Compete(void **h_exists, void *h_in)"); /* %TRACE% */}


#if 0
  /** Execution compete function needed by MPIH.
   *
   * The execution of the exception handler only depends on the
   * handler having the virtual function what() to return an error
   * string.  This is just an example of what could be done.
   * For Framework, each class defines it's own callback function
   * and Handler_Execute is not needed.
   */

  static void Handler_Execute (void *X,
			       int  processor_num) {
    /* %TRACE[TRACEBACK]% */ Traceback trace__("Handler_Execute(void *X, int processor_num)"); /* %TRACE% */

    /* Determining if a void pointer points to a certain class is tricky.
     * It is not possible to call typeid or dynamic_cast on a void pointer
     * directly.  Even after casting to a class pointer, typeid will
     * cause a segmentation fault if the pointer is not actually
     * derived from that class.  So I cast to std::exception then
     * dynamic cast to all known derived classes.  This avoids
     * segmentation errors, but limits the handling to just the
     * known exception classes.
     */

    std::exception *std_X = static_cast<std::exception *>(X);
    std_X = dynamic_cast<std::exception *>(std_X);

    if (std_X) {
      std::cerr
	<< "Internal Exception handling error." << endl
	<< "Function Handler_Execute recieved an exception" << endl
	<< "not derived from standard library base class." << endl
	<< "Therefore no specific exception error message is available."
	<< endl;
      Env::output()
	<< "Internal Exception handling error." << endl
	<< "Function Handler_Execute recieved an exception" << endl
	<< "not derived from standard library base class." << endl
	<< "Therefore no specific exception error message is available."
	<< endl;
    } else {
      int me;
      try {
	me = Env::parallel_rank();
      }
      catch (...) {
	std::cerr << "Exception thrown during processing of exception handle" << endl
		  << "Error messages will be printed for processor 0 only." << endl;
	me =0;
      }
      if (processor_num==me) {
	Env::output()
	  << "*******************************************" << endl
	  << "Exception thrown on processor " << me << ":" << endl
	  << std_X->what() << endl
	  << "*******************************************" << endl;
      }
    }
    return;
  }

#endif

}

namespace sierra {
namespace mpih {

void
Sub_Communicator(
  MPI_Comm		old_comm,
  MPI_Comm		new_comm)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Sub_Communicator(MPI_Comm old_comm, MPI_Comm new_comm)"); /* %TRACE% */
  int       result=MPI_SUCCESS;

  result = MPIH_Copy_Key(old_comm, new_comm);
  if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Sub_Communicator :" << std::endl
      << "Error return from MPIH_Copy_Key:" << result << std::endl << StackTrace;
  }
}


void
Enable()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Enable()"); /* %TRACE% */

  /* A note about these tag numbers.  They are just random
   * numbers.
   * They are sequential to make it easier to keep
   * track of them.  Every MPIH application has to choose
   * tag numbers to use, but there is no way to register
   * or check which numbers are being used by other
   * functions.  Since there are lots of numbers and only
   * a few used as tags, the chances of accidentally using
   * the same tags is small, but non-zero.  Tag collisions
   * would probably cause erratic behavior that would
   * be hard to diagnose.
   * Within Sierra only Framework uses MPIH, so there
   * can be no problem with duplicate tag numbers.
   */
  const int MPIH_TAG_SPARSE    = 12131;
  const int MPIH_TAG_NOMINAL   = 12132;
  const int MPIH_TAG_EXCEPTION = 12133;

  int       result=MPI_SUCCESS;

  result = MPIH_Comm_Enable(Env::parallel_comm(),
			    MPIH_TAG_SPARSE,
			    MPIH_TAG_NOMINAL,
			    MPIH_TAG_EXCEPTION,
			    Handler_Compete,
			    NULL);
  if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Enable :" << std::endl
      << "Error return from MPIH_Comm_Enable:" << result << std::endl << StackTrace;
  }
}

void Keyval_delete(MPI_Comm comm)
{
  MPIH_Keyval_delete(comm);
}

void
Map_Query(
  MPIH_Map          map         /* in:  Map for the sparse operation */,
  int              *symmetric   /* out:                         */,
  size_t           *nsend       /* out: Number of sends         */,
  std::vector<int> *sendlist    /* out: List of destinations    */,
  std::vector<int> *sendlength  /* out:                         */,
  std::vector<int> *sendbuflen  /* out:                         */,
  size_t           *nrecv       /* out: Number of receives      */,
  std::vector<int> *recvlist    /* out: List of destinations    */,
  std::vector<int> *recvlength  /* out:                         */,
  std::vector<int> *recvbuflen  /* out:                      */ )
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Map_Query(MPIH_Map map , int *symmetric , size_t *nsend , std::vector<Int> *sendlist , std::vector<Int> *sendlength , std::vector<Int> *sendbuflen , size_t *nrecv , std::vector<Int> *recvlist , std::vector<Int> *recvlength , std::vector<Int> *recvbuflen  )"); /* %TRACE% */
  int
    c_symmetric = 0,
    c_nsend     = 0,
    c_nrecv     = 0;

  if (!map) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Map_Query." << std::endl
      << "MPIH_Map to query is a NULL pointer." << std::endl << StackTrace;
  }
  MPIH_Map_query(   map     ,
		    &c_symmetric ,
		    &c_nsend   ,
		    NULL      ,
		    NULL      ,
		    NULL      ,
		    &c_nrecv   ,
		    NULL      ,
		    NULL      ,
		    NULL      );
  if ( symmetric  ) *symmetric = c_symmetric;
  if ( nsend      ) *nsend     = c_nsend;
  if ( nrecv      ) *nrecv     = c_nrecv;
  if ( sendlist   ) sendlist  ->resize(c_nsend);
  if ( sendlength ) sendlength->resize(c_nsend);
  if ( sendbuflen ) sendbuflen->resize(c_nsend);
  if ( recvlist   ) recvlist  ->resize(c_nrecv);
  if ( recvlength ) recvlength->resize(c_nrecv);
  if ( recvbuflen ) recvbuflen->resize(c_nrecv);
  const bool non_null =  sendlist || sendlength || sendbuflen  ||
    recvlist || recvlength || recvbuflen;
  if (non_null) {
    int
      *c_sendlist   = ( sendlist   && c_nsend ) ? &(*sendlist  )[0] : NULL,
      *c_sendlength = ( sendlength && c_nsend ) ? &(*sendlength)[0] : NULL,
      *c_sendbuflen = ( sendbuflen && c_nsend ) ? &(*sendbuflen)[0] : NULL,
      *c_recvlist   = ( recvlist   && c_nrecv ) ? &(*recvlist  )[0] : NULL,
      *c_recvlength = ( recvlength && c_nrecv ) ? &(*recvlength)[0] : NULL,
      *c_recvbuflen = ( recvbuflen && c_nrecv ) ? &(*recvbuflen)[0] : NULL;

    MPIH_Map_query(   map         ,
		      NULL        ,
		      NULL        ,
		      c_sendlist  ,
		      c_sendlength,
		      c_sendbuflen,
		      NULL        ,
		      c_recvlist  ,
		      c_recvlength,
		      c_recvbuflen);
  }
}


void
Sparse_Map(
  const std::vector<int> &lengths  /* in: Byte length of each message  */ ,
  const std::vector<int> &buflens  /* in: Byte length of each message buffer*/,
  const std::vector<int> &sendlist /* in: Destination processors       */ ,
  MPIH_Map * map                  /* out: Map for the sparse operation */ )
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Sparse_Map(const std::vector<Int> lengths  , const std::vector<Int> buflens , const std::vector<Int> sendlist  , MPIH_Map * map  )"); /* %TRACE% */
  int result=MPI_SUCCESS;
  const unsigned int nsend = lengths.size();

  if (nsend != buflens.size()) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Sparse_Map." << std::endl
      << "Length of size array:" << nsend << std::endl
      << "Not equal to length of buffer size array:" << buflens.size() << std::endl << StackTrace;
  }
  if (nsend != sendlist.size()) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Sparse_Map." << std::endl
      << "Length of size array:" << nsend << std::endl
      << "Not equal to length of sendlist array:" << sendlist.size() << std::endl << StackTrace;
  }

  const int* const c_lengths  = ( nsend ) ? &lengths [0] : NULL;
  const int* const c_buflens  = ( nsend ) ? &buflens [0] : NULL;
  const int* const c_sendlist = ( nsend ) ? &sendlist[0] : NULL;
  result = MPIH_Sparse_Map( nsend     ,
			    c_lengths ,
			    c_buflens ,
			    c_sendlist,
			    Env::parallel_comm(),
			    map       );

  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Sparse_Map." << std::endl
      << "MPIH_Sparse_Map returned error:" << result << std::endl << StackTrace;
  }
}


void
Sparse_Symmetric_Map(
  const std::vector<int> &	lengths		/* in:  Length of each message       */ ,
  const std::vector<int> &	buflens		/* in:  Length of each message       */ ,
  const std::vector<int> &	sendlist	/* in:  Destination processors       */ ,
  MPIH_Map *			map)		/* out: Map for the sparse operation */
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Sparse_Symmetric_Map(const std::vector<Int> lengths  , const std::vector<Int> buflens  , const std::vector<Int> sendlist  , MPIH_Map * map  )"); /* %TRACE% */
  int result=MPI_SUCCESS;
  const unsigned int nsend = lengths.size();

  if (nsend != buflens.size()) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Sparse_Map." << std::endl
      << "Length of size array:" << nsend << std::endl
      << "Not equal to length of buffer size array:" << buflens.size() << std::endl << StackTrace;
  }
  if (nsend != sendlist.size()) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Sparse_Map." << std::endl
      << "Length of size array:" << nsend << std::endl
      << "Not equal to length of sendlist array:" << sendlist.size() << std::endl << StackTrace;
  }

  const int* const c_lengths  = ( nsend ) ? &lengths [0] : NULL;
  const int* const c_buflens  = ( nsend ) ? &buflens [0] : NULL;
  const int* const c_sendlist = ( nsend ) ? &sendlist[0] : NULL;
  result = MPIH_Sparse_Symmetric_Map( nsend     ,
				      c_lengths   ,
				      c_buflens   ,
				      c_sendlist  ,
				      Env::parallel_comm(),
				      map       );
  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Sparse_Symmetric_Map." << std::endl
      << "MPIH_Sparse_Symmetric_Map returned error:" << result << std::endl << StackTrace;
  }
}


void
Add_Handle(
  const ExParallel &		X)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Add_Handle(const ExParallel &X)"); /* %TRACE% */
//  ExParallel &Ex = X.singleton();
  int result = MPIH_Comm_add_handle(Env::parallel_comm(),
				    static_cast<void *>(const_cast<ExParallel *>(&X)));
  if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Add_Handle:" << std::endl
      << "Error returned from MPIH_Comm_add_handle:" << result << std::endl;
  }
}


void
Delete_Handles()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Delete_Handles()"); /* %TRACE% */
  int result = MPIH_Comm_delete_handles(Env::parallel_comm());
  if (MPI_SUCCESS != result) {
    RuntimeWarning()
      << "Error in function sierra::mpih::Delete_Handles:" << std::endl
      << "Error returned from MPIH_Comm_delete_handles:" << result << std::endl
      << " Non-Fatal Error." << std::endl;
  }
}


void
Activate_Handles()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Activate_Handles()"); /* %TRACE% */
  int result=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_activate_handles(Env::parallel_comm()))) {
    Env::output()
      << "Error in function sierra::mpih::Activate_Handles." << std::endl
      << "MPIH_Comm_activate_handles returned error:" << result << std::endl
      << "Parallel exception handling failed." << std::endl;
  }
}


void
Deactivate_Handles()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Deactivate_Handles()"); /* %TRACE% */
  int result=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_deactivate_handles(Env::parallel_comm()))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Deactivate_Handles." << std::endl
      << "MPIH_Comm_deactivate_handles returned error:" << result << std::endl << StackTrace;
  }
}


ExParallel *
Get_Local_Handle()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Get_Local_Handle()"); /* %TRACE% */

  void * p;

  int result = MPIH_Comm_get_local_handle(Env::parallel_comm(), &p);
  if (result != MPI_SUCCESS)
    throw RuntimeError()
      << "Error in function sierra::mpih::Get_Local_Handle." << std::endl
      << "MPIH_Comm_get_local_handle returned error:" << result << std::endl << StackTrace;

  ExParallel *handle = reinterpret_cast<ExParallel*>(p);

  return handle;
}


void
Set_Local_Handle(
  ExParallel &handle)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Set_Local_Handle(ExParallel &handle)"); /* %TRACE% */
  int result=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_set_local_handle(Env::parallel_comm(),
					   static_cast<void *>(&handle)))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Set_Local_Handle." << std::endl
      << "MPIH_Comm_set_local_handle returned error:" << result << std::endl << StackTrace;
  }
}


void
Reset_Local_Handle()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Reset_Local_Handle()"); /* %TRACE% */
  int result=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_reset_status(Env::parallel_comm()))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Reset_Local_Handle." << std::endl
      << "MPIH_Comm_reset_status returned error:" << result << std::endl << StackTrace;
  }
}


int
Get_Global_Status()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Get_Global_Status()"); /* %TRACE% */
  int result=-1, status=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_get_global_status(Env::parallel_comm(), &status))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Get_Global_Status." << std::endl
      << "MPIH_Comm_get_global_status returned error:" << result << std::endl << StackTrace;
  }
  return status;
}


void
Set_Status_Check()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Set_Status_Check()"); /* %TRACE% */
  int result=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_set_status_check(Env::parallel_comm()))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Set_Status_Check." << std::endl
      << "MPIH_Comm_set_status_check returned error:" << result << std::endl << StackTrace;
  }
}


void
Reset_Status_Check()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Reset_Status_Check()"); /* %TRACE% */
  int result=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_reset_status_check(Env::parallel_comm()))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Reset_Status_Check." << std::endl
      << "MPIH_Comm_reset_status_check returned error:" << result << std::endl << StackTrace;
  }
}


int
Get_Status_Check()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Get_Status_Check()"); /* %TRACE% */
  int result=-1, status=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_get_status_check(Env::parallel_comm(),
					   &status))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Get_Status_Check." << std::endl
      << "MPIH_Comm_get_status_check returned error:" << result << std::endl << StackTrace;
  }
  return status;
}


int
Get_Nhandles()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Get_Nhandles()"); /* %TRACE% */
  int result=-1, nhandles=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_get_nhandles(Env::parallel_comm(),
				       &nhandles))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Get_Nhandles." << std::endl
      << "MPIH_Comm_get_nhandles returned error:" << result << std::endl << StackTrace;
  }
  return nhandles;
}


void
Get_Handles(
  ExParallel **handles)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Get_Handles(ExParallel **handles)"); /* %TRACE% */
  int result=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_get_handles(Env::parallel_comm(),
				      reinterpret_cast<void **>(handles)))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Get_Handles." << std::endl
      << "MPIH_Comm_get_handles returned error:" << result << std::endl << StackTrace;
  }
}


void
Get_Tags(
  int *			active,
  int *			tag_sparse,
  int *			tag_nominal,
  int *			tag_message)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Get_Tags(int *active, int *tag_sparse, int *tag_nominal, int *tag_message)"); /* %TRACE% */
  int result=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_get_tags(Env::parallel_comm(),
				   active,
				   tag_sparse,
				   tag_nominal,
				   tag_message))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Get_Tags." << std::endl
      << "MPIH_Comm_get_tags returned error:" << result << std::endl << StackTrace;
  }
}


void
Get_Functions(
  MPIH_Handler_compete *	handler_compete_fn ,
  MPIH_Handler_execute *	handler_execute_fn)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Get_Functions(MPIH_Handler_compete *handler_compete_fn , MPIH_Handler_execute *handler_execute_fn)"); /* %TRACE% */
  int result=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_get_functions(Env::parallel_comm(),
					handler_compete_fn ,
					handler_execute_fn))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Get_Functions." << std::endl
      << "MPIH_Comm_get_functions returned error:" << result << std::endl << StackTrace;
  }
}


int
Get_Control_Message()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Get_Control_Message()"); /* %TRACE% */
  return MPIH_Comm_get_control_message(Env::parallel_comm());
}


void
Get_Global_Handles(
  ExParallel **		handles)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Get_Global_Handles(ExParallel **handles)"); /* %TRACE% */
  int result=-1;
  if (MPI_SUCCESS !=
      (result = MPIH_Comm_get_global_handles(Env::parallel_comm(),
					     reinterpret_cast<void **>(handles)))) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Get_Global_Handles." << std::endl
      << "MPIH_Comm_get_global_handles returned error:" << result << std::endl << StackTrace;
  }
}


void
Bcast(
  void *		buffer,
  int			count,
  MPI_Datatype		datatype,
  int			root)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Bcast(void *buffer, int count, MPI_Datatype datatype, int root)"); /* %TRACE% */
  int result=MPI_SUCCESS;
  result = MPIH_Bcast(buffer,
		      count,
		      datatype,
		      root,
		      Env::parallel_comm());

  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Bcast." << std::endl
      << "MPIH_Bcast returned error:" << result << std::endl << StackTrace;
  }
}


void
Allreduce(
  void *		in_buffer,
  void *		out_buffer,
  int			count,
  MPI_Datatype		datatype,
  MPI_Op		op )
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Allreduce( void *in_buffer, void *out_buffer, int count, MPI_Datatype datatype, MPI_Op op )"); /* %TRACE% */
  int result=MPI_SUCCESS;
  result = MPIH_Allreduce(in_buffer,
			  out_buffer,
			  count,
			  datatype,
			  op,
			  Env::parallel_comm());

  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Allreduce." << std::endl
      << "MPIH_Allreduce returned error:" << result << std::endl << StackTrace;
  }
}


void
Gather(
  void *		send_buf,
  int			send_size,
  MPI_Datatype		send_datatype,
  void *		recv_buf,
  int			recv_size,
  MPI_Datatype		recv_datatype,
  int			root)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Gather( void *send_buf, int send_size, MPI_Datatype send_datatype, void *recv_buf, int recv_size, MPI_Datatype recv_datatype, int root)"); /* %TRACE% */
  int result=MPI_SUCCESS;
  result = MPIH_Gather(send_buf,
		       send_size,
		       send_datatype,
		       recv_buf,
		       recv_size,
		       recv_datatype,
		       root,
		       Env::parallel_comm());

  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Gather." << std::endl
      << "MPIH_Gather returned error:" << result << std::endl << StackTrace;
  }
}


void
Reduce(
  void *		in_buffer,
  void *		out_buffer,
  int			count,
  MPI_Datatype		datatype,
  MPI_Op		op,
  int			root)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Reduce( void *in_buffer, void *out_buffer, int count, MPI_Datatype datatype, MPI_Op op, int root)"); /* %TRACE% */
  int result=MPI_SUCCESS;
  result = MPIH_Reduce( in_buffer,
			out_buffer,
			count,
			datatype,
			op,
			root,
			Env::parallel_comm());

  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Reduce." << std::endl
      << "MPIH_Reduce returned error:" << result << std::endl << StackTrace;
  }
}


void
Reduce_Scatter(
  void *		in_buffer,
  void *		out_buffer,
  int			recv_count,
  MPI_Datatype		datatype,
  MPI_Op		op)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Reduce_Scatter( void *in_buffer, void *out_buffer, int recv_count, MPI_Datatype datatype, MPI_Op op)"); /* %TRACE% */
  int result=MPI_SUCCESS;
  result = MPIH_Reduce_scatter(in_buffer,
			       out_buffer,
			       recv_count,
			       datatype,
			       op,
			       Env::parallel_comm());

  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Reduce_Scatter." << std::endl
      << "MPIH_Reduce_scatter returned error:" << result << std::endl << StackTrace;
  }
}


void
Scatter(
  void *		send_buf,
  int			send_size,
  MPI_Datatype		send_datatype,
  void *		recv_buf,
  int			recv_size,
  MPI_Datatype		recv_datatype,
  int			root)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Scatter( void *send_buf, int send_size, MPI_Datatype send_datatype, void *recv_buf, int recv_size, MPI_Datatype recv_datatype, int root)"); /* %TRACE% */
  int result=MPI_SUCCESS;
  result = MPIH_Scatter(send_buf,
			send_size,
			send_datatype,
			recv_buf,
			recv_size,
			recv_datatype,
			root,
			Env::parallel_comm());

  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Scatter." << std::endl
      << "MPIH_Scatter returned error:" << result << std::endl << StackTrace;
  }
}


void
Map_Free(
  MPIH_Map *		map)
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Map_Free(MPIH_Map * map)"); /* %TRACE% */
  MPIH_Map_free(map);
}


void
Sparse(
  void *		sendbuf,			///< in:  address of send buffer
  MPI_Datatype		sendtype,			///< in:  datatype of the send messages
  void *		recvbuf,			///< in:  address of receive buffer
  MPI_Datatype		recvtype,			///< in:  datatype of the recv messages
  int			transpose,			///< in:  whether to reverse communic
  MPIH_Map		map)				///< in:  communication map
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Sparse(void * sendbuf  , MPI_Datatype sendtype  , void * recvbuf  , MPI_Datatype recvtype  , int transpose  , MPIH_Map map  )"); /* %TRACE% */
  int result=MPI_SUCCESS;
  result = MPIH_Sparse( sendbuf   ,
			sendtype  ,
			recvbuf   ,
			recvtype  ,
			transpose ,
			map       );
  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Sparse." << std::endl
      << "MPIH_Sparse returned error:" << result << std::endl << StackTrace;
  }
}

void
Initialize_Sparse(
  void  *		sendbuf,			///< in:  address of send buffer
  MPI_Datatype		sendtype,			///< in:  datatype of the send messages
  void *		recvbuf,			///< in:  address of receive buffer
  MPI_Datatype		recvtype,			///< in:  datatype of the recv messages
  int			transpose,			///< in:  whether to reverse communication
  MPIH_Map		map)				///< in:  communication map
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Initialize_Sparse(void * sendbuf, MPI_Datatype sendtype  , void * recvbuf  , MPI_Datatype recvtype  , int transpose  , MPIH_Map map )"); /* %TRACE% */
  int result=MPI_SUCCESS;
  result = MPIH_Isparse( sendbuf   ,
			 sendtype  ,
			 recvbuf   ,
			 recvtype  ,
			 transpose ,
			 map       );
  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Initialize_Sparse." << std::endl
      << "MPIH_Isparse returned error:" << result << std::endl << StackTrace;
  }
}

void
Wait_Sparse(
  MPIH_Map		map)				///< in:  communication map
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::Wait_Sparse(MPIH_Map map )"); /* %TRACE% */
  int result=MPI_SUCCESS;
  result = MPIH_Wait(map);
  if (MPIH_Has_Exception_Flag() == result) {
    parallel_throw(Env::parallel_comm());
  }
  else if (MPI_SUCCESS != result) {
    throw RuntimeError()
      << "Error in function sierra::mpih::Wait_Sparse." << std::endl
      << "MPIH_Wait returned error:" << result << std::endl << StackTrace;
  }
}


const char *
get_product_name()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::get_product_name()"); /* %TRACE% */
  return "MPIH";
}


void
register_product()
{
  /* %TRACE[TRACEBACK]% */ Traceback trace__("sierra::mpih::register_product()"); /* %TRACE% */

  // Register mpih
  ProductRegistry::AttributeMap &attr_map = ProductRegistry::instance().addProduct(get_product_name());
  attr_map[ProductRegistry::VERSION] = ProductRegistry::version();
  attr_map[ProductRegistry::TITLE] = "MPI parallel exception handling";
  attr_map[ProductRegistry::CONTACT] = "framework-developers@sourceforge.sandia.gov";

  // Register TPL's and other things which may not be properly registered but used directly.

}

} // namespace mpih
} // namespace sierra
