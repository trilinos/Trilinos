/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "msg_const.h"


/*
 * interface to message-passing mechanism
 */
static void msg_abort(MPI_Comm, int, int errcode);

/*****************************************************************************/
/*
 * Zoltan_Oct_msg_int_scan
 *
 * perform "exclusive" scan
 *
 * in    v0   v1   v2     v3       ...
 * out    0   v0   v0+v1  v0+v1+v2 ...
 *
 */
int Zoltan_Oct_msg_int_scan(MPI_Comm communicator, int proc, int value)
{
  int recvbuf;
  int ret;

  ret=MPI_Scan(&value,&recvbuf,1,MPI_INT,MPI_SUM,communicator);

  if (ret!= MPI_SUCCESS) {
    fprintf(stderr,"%d(%d) Zoltan_Oct_msg_int_scan: Scan ret=%d\n",proc, proc,ret);
    msg_abort(communicator, proc, ret);
  }
  
  return(recvbuf-value);
}

/*****************************************************************************/
/*
 * Zoltan_Oct_msg_float_scan
 *
 * perform "exclusive" scan
 *
 * in    v0   v1   v2     v3       ...
 * out    0   v0   v0+v1  v0+v1+v2 ...
 *
 */
float Zoltan_Oct_msg_float_scan(MPI_Comm communicator, int proc, float value)
{
  float recvbuf;
  int ret;

  ret=MPI_Scan(&value,&recvbuf,1,MPI_FLOAT,MPI_SUM,communicator);

  if (ret!= MPI_SUCCESS) {
    fprintf(stderr,"%d(%d) Zoltan_Oct_msg_float_scan: Scan ret=%d\n",proc,proc,ret);
    msg_abort(communicator, proc, ret);
  }

  return(recvbuf-value);
}

/*****************************************************************************/
/*
 * msg_abort(int errcode)
 *
 * Try to abort all tasks
 *
 */
static void msg_abort(MPI_Comm communicator, int proc, int errcode)
{
  char errmsg[MPI_MAX_ERROR_STRING];
  int errsize;

  MPI_Error_string(errcode,errmsg,&errsize);
  fprintf(stderr,"%d  error string: %s\n",proc,errmsg);
  MPI_Abort(communicator,errcode);
  abort();
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
