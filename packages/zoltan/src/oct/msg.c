#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "msg_const.h"


/*
 * interface to message-passing mechanism
 *
 */
int msg_mypid;                 /* GLOBALS !!! */
int msg_temp;
int msg_nprocs;

static int *receives;
static int bytes_sent;

/*
 * msg_init()
 *
 * call this for initialization before calling any other msg_xxx() routines
 * Among other things, will set globals msg_nprocs and msg_mypid
 *
 */
void msg_init()
{
  int ret;
  void *bufptr;

  /*
  ret=MPI_Init(argc,argv); 
  
  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"msg_inti: MPI_Init() failed\n");
    msg_abort(ret);
  }
  */

  ret=MPI_Comm_rank(MPI_COMM_WORLD,&msg_mypid);
  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"msg_init: MPI_Comm_rank() failed\n");
    msg_abort(ret);
  }

  ret=MPI_Comm_size(MPI_COMM_WORLD,&msg_nprocs);
  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"msg_init: MPI_Comm_size() failed\n");
    msg_abort(ret);
  }

  bufptr=(int *)malloc(MSG_BUFSIZE);

  ret=MPI_Buffer_attach(bufptr,MSG_BUFSIZE);
  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"msg_init: MPI_Buffer_attach() failed\n");
    msg_abort(ret);
  }

  receives=(int *)calloc(msg_nprocs, sizeof(int));
  if (!receives) {
    fprintf(stderr,"msg_init: could not malloc receive array\n");
    abort();
  }
  
}


void msg_end(void)
{ MPI_Finalize(); }

/*
 * msg_send_init()
 *
 * call this to initialize the data structure that keeps track
 * of sends.  After calling msg_bsend() as much as desired, this
 * info is used by msg_nreceives() to figure out how many messages
 * each processor will receive.
 *
 */
void msg_send_init()
{
  int i;

  for (i=0; i<msg_nprocs; i++)
    receives[i]=0;

  bytes_sent=0;
}

/*
 * msg_barrier();
 *
 */
void msg_barrier()
{ MPI_Barrier(MPI_COMM_WORLD); }


/*
 * msg_bsend(buf,size,dest,type)
 *
 * Do a send, blocking until buf is free.
 * Checks error codes.
 *
 */
void msg_bsend(void *buf, int size, int dest, int type) {
  int ret;

  if (dest==msg_mypid) {
    fprintf(stderr,"%d msg_bsend: attempt to send to self\n",msg_mypid);
    abort();
  }

  if (receives)
    receives[dest]++;
  bytes_sent+=size;
  
  ret=MPI_Bsend(buf,size,MPI_BYTE,dest,type,MPI_COMM_WORLD);

  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"%d msg_bsend: size=%d dest=%d type=%d errno=%d\n",
	    msg_mypid,size,dest,type,ret);
    msg_abort(ret);
  }
}


/*
 * msg_breceive(buf,size,from,type)
 *
 * receive a message of the specified size and type from any other
 * processor.  Size and type must match.
 * Checks for errors.
 *
 */
void msg_breceive(void *buf, int size, int *from, int type) {
  int ret;
  int nbytes;
  MPI_Status status;

  ret=MPI_Recv(buf,size,MPI_BYTE,MPI_ANY_SOURCE,type,MPI_COMM_WORLD,&status);

  *from=status.MPI_SOURCE;
  /* MPI_ANY_TAG, status.MPI_TAG */
  
  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"%d msg_breceive: receive error ret=%d\n",msg_mypid,ret);
    msg_abort(ret);
  }

  ret=MPI_Get_count(&status,MPI_BYTE,&nbytes);

  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"%d msg_breceive: get count failed\n",msg_mypid);
    msg_abort(ret);
  }

  if (nbytes!=size) {
    fprintf(stderr,"%d msg_breceive: expected %d bytes but got %d bytes\n",
	    msg_mypid,size,nbytes);
    abort();
  }
}


/*
 * msg_breceive_any(buf,size,from,type)
 *
 * receive a variable-sized message of the specified 
 * type from any other processor.
 *
 * Checks for errors.
 * Caller must free the buffer when done.
 *
 */
void msg_breceive_anysize(void **buffer_ret,
                          int *size_ret,
                          int *from_ret,
                          int type)
{
  void *buffer;
  int ret;
  int probesize, recvsize;
  MPI_Status probestat, recvstat;

  ret=MPI_Probe(MPI_ANY_SOURCE,type,MPI_COMM_WORLD,&probestat);

  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"msg_breceive_anysize: probe error\n");
    msg_abort(ret);
  }

  ret=MPI_Get_count(&probestat,MPI_BYTE,&probesize);

  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"msg_breceive_anysize: get probesize failed\n");
    msg_abort(ret);
  }

  buffer= (void *) malloc(probesize);

  if (!buffer) {
    fprintf(stderr,"msg_breceive_anysize: malloc failure %d bytes\n",
	    probesize);
    abort();
  }

  ret=MPI_Recv(buffer,probesize,MPI_BYTE,
	       probestat.MPI_SOURCE,probestat.MPI_TAG,
	       MPI_COMM_WORLD,&recvstat);

  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"msg_breceive_anysize: receive error\n");
    msg_abort(ret);
  }

  ret=MPI_Get_count(&recvstat,MPI_BYTE,&recvsize);

  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"%d msg_breceive: get recvsize failed\n",msg_mypid);
    msg_abort(ret);
  }

  if (probesize!=recvsize || recvstat.MPI_TAG != type) {
    fprintf(stderr,"msg_breceive: recv doesn't match probe\n");
    abort();
  }

  *buffer_ret=buffer;
  *size_ret=recvsize;
  *from_ret=recvstat.MPI_SOURCE;
}

/*
 * msg_nreceives()
 *
 * all processors should call this to see what they are
 * sending each other.
 *
 * Usage: call msg_send_init(), one or more msg_bsend(),
 * then this to determine how many each processor will receive.
 *
 */
int msg_nreceives()
{
  int rectot,inmsg;
  int i;
  int ret;

  rectot= -1;


  for (i=0; i<msg_nprocs; i++) {
    inmsg=0;
    
    ret=MPI_Reduce((void *)&(receives[i]),
		   (void *)&inmsg,
		   1,                     /* one MPI_INT */
		   MPI_INT,
		   MPI_SUM,
		   i,                     /* root */
		   MPI_COMM_WORLD);
    
    if (ret!= MPI_SUCCESS) {
      fprintf(stderr,"%d msg_nreceives: reduce ret=%d\n",msg_mypid,ret);
      msg_abort(ret);
    }

    if (msg_mypid==i)
      rectot=inmsg;
  }


  /* MLS: re-initialize processor array once done reducing it 5/29/97 */
  msg_send_init();
  
  /* printf("msg_nreceives: %7d bytes sent\n",bytes_sent); */
  return(rectot);

}


/* Written by: MLS
 *
 * msg_pending(int *my_receives)
 *
 * All processors should call this to see how many messages they
 * are receiving. The return value is TRUE if there are still messages
 * in the system (even if my_receives == 0), and FALSE otherwise.
 * Use when the number of communication cycles is unknown.
 *
 */ 

int msg_pending(int *my_receives) {

  int total;         /* total number of messages in the system */

  *my_receives = msg_nreceives();

  MPI_Allreduce(my_receives, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  return(total);
}


/*
 * msg_int_sum
 *
 * Sum a value from each processor, and return sum to all
 * processors.
 *
 */
int msg_int_sum(int value)
{
  int recvbuf;
  int ret;

  ret=MPI_Allreduce(&value,&recvbuf,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  if (ret!= MPI_SUCCESS) {
    fprintf(stderr,"%d msg_int_reduce: Allreduce ret=%d\n",msg_mypid,ret);
    msg_abort(ret);
  }

  return(recvbuf);
}




/*
 * msg_double_sum
 *
 * Sum a value from each processor, and return sum to all
 * processors.
 *
 */

double msg_double_sum(double value)
{
  double recvbuf;
  int ret;

  ret=MPI_Allreduce(&value,&recvbuf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if (ret!= MPI_SUCCESS) {
    fprintf(stderr,"%d msg_double_reduce: Allreduce ret=%d\n",msg_mypid,ret);
    msg_abort(ret);
  }

  return(recvbuf);
}

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
    fprintf(stderr,"%d msg_int_scan: Scan ret=%d\n",msg_mypid,ret);
    msg_abort(ret);
  }
  
  return(recvbuf-value);
}

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
    fprintf(stderr,"%d msg_float_scan: Scan ret=%d\n",msg_mypid,ret);
    msg_abort(ret);
  }

  return(recvbuf-value);
}

/*
 * msg_int_bcast
 *
 * broadcast from one proc to all
 *
 */

int msg_int_bcast(int value, int sender)
{
  int ret;

  ret=MPI_Bcast(&value,1,MPI_INT,sender,MPI_COMM_WORLD);

  if (ret!= MPI_SUCCESS) {
    fprintf(stderr,"%d msg_int_bcast: Bcast ret=%d\n",msg_mypid,ret);
    msg_abort(ret);
  }

  return(value);
}

/*
 * msg_double_min
 *
 * Find min value across procs
 *
 */
double msg_double_min(double value)
{
  double recvbuf;
  int ret;

  ret=MPI_Allreduce(&value,&recvbuf,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

  if (ret!= MPI_SUCCESS) {
    fprintf(stderr,"%d msg_double_min: Allreduce ret=%d\n",msg_mypid,ret);
    msg_abort(ret);
  }

  return(recvbuf);
}

/*
 * msg_double_max
 *
 * Find max value across procs
 *
 */
double msg_double_max(double value)
{
  double recvbuf;
  int ret;

  ret=MPI_Allreduce(&value,&recvbuf,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (ret!= MPI_SUCCESS) {
    fprintf(stderr,"%d msg_double_max: Allreduce ret=%d\n",msg_mypid,ret);
    msg_abort(ret);
  }

  return(recvbuf);
}

/*
 * msg_sync()
 *
 * synchronize all processors
 *
 */
void msg_sync()
{
  int ret;

  ret=MPI_Barrier(MPI_COMM_WORLD);

  if (ret!=MPI_SUCCESS) {
    fprintf(stderr,"%d msg_sync: failed, ret=%d\n",msg_mypid,ret);
    msg_abort(ret);
  }
}

/*
 * msg_inorder_begin()
 *
 * start a block of code that the processors should execute
 * one-at-a-time, in order of processor number.
 * End the block of code with msg_inorder_end()
 * 
 */
void msg_inorder_begin()
{
  int from;

  msg_sync();

  if (msg_mypid>0)
    msg_breceive(NULL,0,&from,MTYPE_INORDER);
}


/* 
 * msg_inorder_end()
 *
 * end a block of code that the processors should execute
 * one-at-a-time, in order of processor number.
 * Start the block of code with msg_inorder_begin()
 *
 */
void msg_inorder_end()
{

  if (msg_mypid+1<msg_nprocs)
    msg_bsend(NULL,0,msg_mypid+1,MTYPE_INORDER);

  msg_sync();

}

/*
 * msg_abort(int errcode)
 *
 * Try to abort all tasks
 *
 */
void msg_abort(int errcode)
{
  char errmsg[MPI_MAX_ERROR_STRING];
  int errsize;

  MPI_Error_string(errcode,errmsg,&errsize);
  fprintf(stderr,"%d  error string: %s\n",msg_mypid,errmsg);
  MPI_Abort(MPI_COMM_WORLD,errcode);
  abort();
}




#if MSGMAIN

main(argc,argv)
int argc;
char *argv[];
{
  int i,nr,from,data;

  msg_init(&argc,&argv);
  msg_send_init();

  for (i=0; i<msg_mypid; i++)
    {
      printf("sending to %d\n",i);
      data=100*msg_mypid+i;
      msg_bsend(&data,sizeof(int),i,0);
    }

  /* sleep(1); */
  msg_sync();

  nr=msg_nreceives();
  printf("receives: %d\n",nr);


  sleep(1);
  msg_sync();

  for (i=0; i<nr; i++)
    {
      msg_breceive(&data,sizeof(int),&from,0);
      printf("received from %d : %d\n",from,data);
    }

  msg_sync();

}
#endif

