/*--------------------------------------------------------------------*/
/*    Copyright 2002, 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
/** \file
    Framework interface to MPIH system.

    sierra::mpih provides a C++ interface to the C based MPIH system.
    It handles all of the explicit exception handling and run time
    type identification that is not possible in C.  Sierra applications
    should not have to call MPIH directly, but should call sierra::mpih.

*/

#ifndef STK_UTIL_PARALLEL_mpih_h
#define STK_UTIL_PARALLEL_mpih_h

#include <vector>
#include <mpi.h>

#include <stk_util/parallel/Exception.hpp>

#include <MPIH_Include.h>

namespace sierra {
namespace mpih {

/**
 * @brief Member function <b>Enable</b> initializes exception handling.
 *
 * <b>Enable()</b> wraps the MPIH_enable fuction call, defaulting all of the parameters.
 * <b>Enable()</b> should be called after MPI_Init and before Register_Handles.
 *
 * @see MPIH_Comm_Enable
 */
void Enable();

void Keyval_delete(MPI_Comm comm);

/**
 * @brief Member function <b>Sub_Communicator</b> initialize a new communicator from
 * an old one.  This is an unusual call in that the communicators have to be explicitly
 * passed.  Normally the default communicator is taken as sierra::Fmwk::Env::comm().
 * Sub_Communicator is used when Env is reset with a partition of communicators.  Each of
 * the new sub communicators need to be initialized from the old one.  This avoids having
 * to register all of the exception handlers again.  The new communicator will inherit
 * everything from the old and even use the same memory to store flags and such.
 *
 * @param old_comm	a <b>MPI_Comm</b> value of the old MPI communicator.
 *
 * @param new_comm	a <b>MPI_Comm</b> value of the new MPI communicator.
 *
 */
void Sub_Communicator (MPI_Comm old_comm,
		       MPI_Comm new_comm);

/**
 * @brief Member function <b>Register_Handles</b> is called during creation of this
 * singleton.  All of the exception subclasses that can be thrown by Fmwk are registered
 * as handles with the MPIH subsystem.  \see Register_Handles
 *
 */
void Register_Handles ();

/**
 * @brief Member function <b>Add_Handle</b> additional exceptions can be added with
 * Add_Handle.  This would be neccessary if an application defined its own exception
 * classes.
 *
 * @see MPIH_Comm_add_handle
 *
 * @param X		an <b>ExParallel</b> reference to the exception singleton to
 *			add to the registered exceptions list.
 */
void Add_Handle (const ExParallel &X);

/**
 * @brief Member function <b>Delete_Handles</b> deletes all exception handles.
 * Should only be called in destructor.
 *
 */
void Delete_Handles ();

/** Collective Communication routine.
 * Bcast, Allreduce, ... are the standard collective
 * communication calls.  These are simple wrappers
 * that default the communicator to sierra::Env::mpi::comm()
 * and check for error return codes.  All of these
 * functions sierra::Env::abort() on error, so there is
 * no return code to check.  If an exception has
 * been thrown, that is automatically handled also.
 * \see MPIH_Allreduce MPIH_Bcast
 * MPIH_Gather MPIH_Reduce MPIH_Reduce_scatter
 * MPIH_Scatter
 */
void Bcast(void        * buffer,
	   int           count,
	   MPI_Datatype  datatype,
	   int           root);

void Allreduce( void        * in_buffer,
		void        * out_buffer,
		int           count,
		MPI_Datatype  datatype,
		MPI_Op        op );

void Gather( void        * send_buf,
	     int           send_size,
	     MPI_Datatype  send_datatype,
	     void        * recv_buf,
	     int           recv_size,
	     MPI_Datatype  recv_datatype,
	     int           root);

void Reduce( void        * in_buffer,
	     void        * out_buffer,
	     int           count,
	     MPI_Datatype  datatype,
	     MPI_Op        op,
	     int           root);

void Reduce_Scatter( void        * in_buffer,
		     void        * out_buffer,
		     int           recv_count,
		     MPI_Datatype  datatype,
		     MPI_Op        op);

void Scatter( void        * send_buf,
	      int           send_size,
	      MPI_Datatype  send_datatype,
	      void        * recv_buf,
	      int           recv_size,
	      MPI_Datatype  recv_datatype,
	      int           root);

void Map_Free(MPIH_Map * map);

void Map_Query
(MPIH_Map          map         /* in:  Map for the sparse operation */,
 int              *symmetric   /* out:                         */,
 size_t           *nsend       /* out: Number of sends         */,
 std::vector<int> *sendlist    /* out: List of destinations    */,
 std::vector<int> *sendlength  /* out:                         */,
 std::vector<int> *sendbuflen  /* out:                         */,
 size_t           *nrecv       /* out: Number of receives      */,
 std::vector<int> *recvlist    /* out: List of destinations    */,
 std::vector<int> *recvlength  /* out:                      */,
 std::vector<int> *recvbuflen  /* out:                      */ );

void Sparse
(void       * sendbuf   /* in:  address of send buffer        */ ,
 MPI_Datatype sendtype  /* in:  datatype of the send messages */ ,
 void       * recvbuf   /* in:  address of receive buffer     */ ,
 MPI_Datatype recvtype  /* in:  datatype of the recv messages */ ,
 int          transpose /* in:  whether to reverse communic.  */ ,
 MPIH_Map     map       /* in:  communication map             */ );

/* The Sparse is a combination of Initialize_Sparse and Wait_Sparse.
 * This allows computation to happen while the communication is
 * begin processed.  After the call to Initialize_Sparse the
 * buffers should NOT be reused until the completion of the
 * corresponding Wait_Sparse call.  The map passed to Wait_Sparse
 * must be the same map that was used in Initialize_Sparse.
 * Also, it is not possible to call Initialize_Sparse twice followed
 * by two calls to Wait_Sparse.  Because MPI does not guarentee
 * the order that messages are received, the first Wait_Sparse
 * could get the messages from the second Initialize_Sparse.
 */
void Initialize_Sparse
(void       * sendbuf   /* in:  address of send buffer        */ ,
 MPI_Datatype sendtype  /* in:  datatype of the send messages */ ,
 void       * recvbuf   /* in:  address of receive buffer     */ ,
 MPI_Datatype recvtype  /* in:  datatype of the recv messages */ ,
 int          transpose /* in:  whether to reverse communic.  */ ,
 MPIH_Map     map       /* in:  communication map             */ );

void Wait_Sparse
(MPIH_Map     map       /* in:  communication map             */ );

void Sparse_Map
(const std::vector<int> &lengths   /* in:  Byte length of each message  */ ,
 const std::vector<int> &buflens   /* in:  Byte length of each message buffer */ ,
 const std::vector<int> &sendlist  /* in:  Destination processors       */ ,
 MPIH_Map * map              /* out: Map for the sparse operation */ );

void Sparse_Symmetric_Map
(const std::vector<int> &lengths   /* in:  Length of each message       */ ,
 const std::vector<int> &buflens   /* in:  Length of each message buffe */ ,
 const std::vector<int> &sendlist  /* in:  Destination processors       */ ,
 MPIH_Map * map       /* out: Map for the sparse operation */ );

inline void ParallelExceptionCheck()
{
  int dummy = 0;
  Bcast(&dummy, 0, MPI_INT, 0);
}


/*
 * The rest of the functions are unlikely to be used
 * outside of mpih.C.
 */

/* Local handle, Global handles, and just plain handles:
 *
 * Just plain handles are the handles that are registered
 *    on this processor.  The number of them is unknown
 *    in advance.  These are required to be consistant
 *    across all processors.
 *
 * local_handle is singular, it refers to the exception
 *    handle that has been set on the local processor.
 *    This is the handle that will be propagated on
 *    the next collective communication.
 *
 * Global handles is the collection of local handles
 *    from all of the processors after the collective
 *    communication.  There will be one for each
 *    processor, but some will be NULL if there was not
 *    a local handle set for that processor.
 */
void Activate_Handles   ();
void Deactivate_Handles ();
int  Get_Nhandles       ();
void Get_Handles        (ExParallel **handles);

ExParallel *Get_Local_Handle ();
void Set_Local_Handle   (ExParallel &handle);
void Reset_Local_Handle ();

void Get_Global_Handles (ExParallel ** handles);
int  Get_Global_Status  ();

void Set_Status_Check   ();
void Reset_Status_Check ();
int  Get_Status_Check   ();

void Get_Tags (int *active,
	       int *tag_sparse,
	       int *tag_normal,
	       int *tag_message);

void Get_Functions (MPIH_Handler_compete *handler_compete_fn ,
		    MPIH_Handler_execute *handler_execute_fn );

int  Get_Control_Message();

/* interface for the product versioning. */
const char *get_product_name();
const char *get_product_version();
const char *get_product_qualifier();
void register_product();

} // namespace mpih
} // namespace sierra

#endif //  STK_UTIL_PARALLEL_mpih_h
