#define MSG_BUFSIZE 2024000
/*#define MSG_BUFSIZE 65024000*/
/*#define MSG_BUFSIZE 90024000*/

extern int msg_nprocs;
extern int msg_mypid;
extern int msg_temp;

extern void   msg_init();
extern void   msg_endBuffer();
extern void   msg_send_init(void);
extern void   msg_barrier(void);
extern void   msg_bsend(void *buf, int size, int dest, int type);
extern void   msg_breceive(void *buf, int size, int *from, int type);
extern void   msg_breceive_anysize(void **buf, int *size, int *from, int type);
extern int    msg_nreceives(void);
extern int    msg_pending(int *my_receives);
extern int    msg_int_sum(int value);
extern double msg_double_sum(double value);
extern int    msg_int_scan(int value);
extern float  msg_float_scan(float value);
extern int    msg_int_bcast(int value, int sender);
extern double msg_double_min(double value);
extern double msg_double_max(double value);
extern void   msg_sync(void);
extern void   msg_inorder_begin(void);
extern void   msg_inorder_end(void);
extern void   msg_end(void);
extern void   msg_abort(int errcode);

#define PRINT_IN_ORDER() \
  msg_sync(); \
  if (msg_mypid==0) \
    printf("************************************************************\n"); \
  for (msg_temp=0 /*, mpc_flush(1)*/; \
       msg_temp<msg_nprocs; \
       msg_temp++ /*, mpc_flush(1)*/ ) \
    if (msg_mypid==msg_temp)

/*
 * Message types
 *
 */
#define MTYPE_INORDER    0
#define MTYPE_OREQ       1
#define MTYPE_OREP       2
#define MTYPE_OUPDATE    3
#define MTYPE_OMIGRATE   4
#define MTYPE_DMIGRATE   5
