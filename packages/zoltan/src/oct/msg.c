/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_msgc_id = "$Id$";
#endif

#include "lb_const.h"
#include "msg_const.h"


/*
 * interface to message-passing mechanism
 */
static void msg_abort(int errcode);

/*****************************************************************************/
/*
 * msg_int_scan
 *
 * perform "exclusive" scan
 *
 * in    v0   v1   v2     v3       ...
 * out    0   v0   v0+v1  v0+v1+v2 ...
 *
 */
int msg_int_scan(int value)
{
  int recvbuf;
  int ret;

  ret=MPI_Scan(&value,&recvbuf,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  if (ret!= MPI_SUCCESS) {
    fprintf(stderr,"%d msg_int_scan: Scan ret=%d\n",LB_Proc,ret);
    msg_abort(ret);
  }
  
  return(recvbuf-value);
}

/*****************************************************************************/
/*
 * msg_float_scan
 *
 * perform "exclusive" scan
 *
 * in    v0   v1   v2     v3       ...
 * out    0   v0   v0+v1  v0+v1+v2 ...
 *
 */
float msg_float_scan(float value)
{
  float recvbuf;
  int ret;

  ret=MPI_Scan(&value,&recvbuf,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

  if (ret!= MPI_SUCCESS) {
    fprintf(stderr,"%d msg_float_scan: Scan ret=%d\n",LB_Proc,ret);
    msg_abort(ret);
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
static void msg_abort(int errcode)
{
  char errmsg[MPI_MAX_ERROR_STRING];
  int errsize;

  MPI_Error_string(errcode,errmsg,&errsize);
  fprintf(stderr,"%d  error string: %s\n",LB_Proc,errmsg);
  MPI_Abort(MPI_COMM_WORLD,errcode);
  abort();
}
