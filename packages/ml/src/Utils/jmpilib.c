#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#ifdef DEBUG
#	define	ASSERT(expr)						\
		if (!(expr))						\
			jerror("assertion failed; line %d, file %s - "	\
				#expr, __LINE__, __FILE__);
#	define	MONITOR(expr)	expr
#else
#	define	ASSERT(expr)		;
#	define	MONITOR(expr)		;
#endif

#define	CALLOC(_ptr_,_type_,_l_) \
		_ptr_ = (_type_*) jalloc(_l_, sizeof(_type_), 0, \
			1, #_type_, __LINE__, __FILE__);
#define	FREE(_ptr_,_type_,_l_) \
		jfree((void **) &_ptr_, _l_, sizeof(_type_), 0, \
			#_type_, __LINE__, __FILE__);

#define	MAXIMUM(_a_,_b_)	((_a_)>(_b_)?(_a_):(_b_))
#define	MINIMUM(_a_,_b_)	((_a_)<(_b_)?(_a_):(_b_))

enum	global_tags	{ SUM, MAX, MAXBYTE, MAXLOC, MIN, MINBYTE, MINLOC, MULTI, DSUM, DMAX, DMIN };

enum	tags	{ ANY_TAG = 8665 };

/*
#define	MESSAGE_LOG
*/

extern	void	jerror(char*,...);
extern	void	jwarning(int,char*,...);
extern	void	*jalloc(int,int,int,int,char*,int,char*);
extern	void	jfree(void**,int,int,int,char*,int,char*);

extern	int	int2dbl_init(void);
extern	void	dbl2int(double,int*,int*);
extern	double	int2dbl(int*,int*);

extern	int	jpid, j_nprocs, j_nparts;

#ifdef DEBUG
static	int	success;
#endif

static	MPI_Request	*request;
static	MPI_Op	multi_op_handle;
static	MPI_Group	non_idle_group;
static	MPI_Group	global_group = NULL;
static	MPI_Comm	original_global_comm, global_comm;

static	int	gs_initialised = 0;
static	int	*request_no;
#ifdef MESSAGE_LOG
static	FILE	*logfile;

static	void	print_message(int n, int *buffer)
{	/* print_message: O(? - MESSAGE_LOG)
	 * write a message to the logfile */
	int	i;

	for (i = 0; i < n; ++i) {
		/* (void) fprintf(logfile," %.12d",buffer[i]); */
		(void) fprintf(logfile," %d",buffer[i]);
		if ((i+1)%6 == 0) (void) fprintf(logfile,"\n");
	}
	if (i%6 != 0) (void) fprintf(logfile,"\n");
	(void) fflush(logfile);
}
#endif

#ifdef DEVELOPMENT
static	int	max_byte(int i, int j)
{	/* max_byte: O(?)
	 * never used or tested */
	int	k = 0;
	int	two_n = 1;

	while (i && j) {
		if (i%2 || j%2)
			k += two_n;
		two_n *= 2;
		i /= 2;
		j /= 2;
	}
	return(k);
}

static	int	min_byte(int i, int j)
{	/* min_byte: O(?)
	 * never used or tested */
	int	k = 0;
	int	two_n = 1;

	while (i && j) {
		if (i%2 && j%2)
			k += two_n;
		two_n *= 2;
		i /= 2;
		j /= 2;
	}
	return(k);
}
#endif

static	void	multi_op(void *ivec, void *iovec, int *n, MPI_Datatype *dptr)
{	/* multi_op: O(G*I)
	 * routine to handle multiple operations for a global reduction */
	int	j;
	int	*invec = (int *) ivec;
	int	*inoutvec = (int *) iovec;
	int	ctr;
	double	value, value1, value2;
	int	k;
	int	nval;

	ASSERT(*dptr == MPI_INT);

	for (j = 0; j < *n;) {
		inoutvec[j] = invec[j];
		nval = invec[j++];
		inoutvec[j] = invec[j];
		switch (invec[j++]) {
		case SUM:
			for (k = 0; k < nval; ++k) {
				inoutvec[j] += invec[j];
				j += 1;
			}
			break;
		case MAX:
			for (k = 0; k < nval; ++k) {
				inoutvec[j] = MAXIMUM(inoutvec[j],invec[j]);
				j += 1;
			}
			break;
#ifdef DEVELOPMENT
		case MAXBYTE:
			for (k = 0; k < nval; ++k) {
				inoutvec[j] = max_byte(inoutvec[j],invec[j]);
				j += 1;
			}
			break;
#endif
		case MAXLOC:
			for (k = 0; k < nval; ++k) {
				if (invec[j] > inoutvec[j]) {
					inoutvec[j] = invec[j];
					inoutvec[j+1] = invec[j+1];
				} else if (invec[j] == inoutvec[j]
				 && invec[j+1] < inoutvec[j+1]) {
					inoutvec[j+1] = invec[j+1];
				}
				j += 2;
			}
			break;
		case MIN:
			for (k = 0; k < nval; ++k) {
				inoutvec[j] = MINIMUM(inoutvec[j],invec[j]);
				j += 1;
			}
			break;
#ifdef DEVELOPMENT
		case MINBYTE:
			for (k = 0; k < nval; ++k) {
				inoutvec[j] = min_byte(inoutvec[j],invec[j]);
				j += 1;
			}
			break;
#endif
		case MINLOC:
			for (k = 0; k < nval; ++k) {
				if (invec[j] < inoutvec[j]) {
					inoutvec[j] = invec[j];
					inoutvec[j+1] = invec[j+1];
				} else if (invec[j] == inoutvec[j]
				 && invec[j+1] < inoutvec[j+1]) {
					inoutvec[j+1] = invec[j+1];
				}
				j += 2;
			}
			break;
		case DSUM:
			for (k = 0; k < nval; ++k) {
				ctr = 0;
				value = int2dbl(&invec[j],&ctr);
				ctr = 0;
				value += int2dbl(&inoutvec[j],&ctr);
				dbl2int(value,inoutvec,&j);
			}
			break;
		case DMAX:
			for (k = 0; k < nval; ++k) {
				ctr = 0;
				value1 = int2dbl(&invec[j],&ctr);
				ctr = 0;
				value2 = int2dbl(&inoutvec[j],&ctr);
				value = MAXIMUM(value1,value2);
				dbl2int(value,inoutvec,&j);
			}
			break;
		case DMIN:
			for (k = 0; k < nval; ++k) {
				ctr = 0;
				value1 = int2dbl(&invec[j],&ctr);
				ctr = 0;
				value2 = int2dbl(&inoutvec[j],&ctr);
				value = MINIMUM(value1,value2);
				dbl2int(value,inoutvec,&j);
			}
			break;
		}
	}
}

void	pjostle_init(int *P, int *p)
{	/* pjostle_init: O(1)
	 * initialise the communications for pjostle */
#ifdef MESSAGE_LOG
	char	fname[99];
#endif

	if (gs_initialised) {
		jwarning(0,"pjostle_init already called");
		return;
	} else
		gs_initialised = 1;

	if (*p < 0 || *p >= *P)
		jerror("pid (= %d) < 0 || pid >= nprocs (= %d)",*p,*P);

	jpid = *p;
	j_nprocs = *P;

	original_global_comm = global_comm = MPI_COMM_WORLD;

	if (global_group == NULL) {
		MONITOR(success =) MPI_Comm_group(MPI_COMM_WORLD, &global_group);
		ASSERT(success >= 0);
	}

	MONITOR(success =) MPI_Op_create(multi_op,1,&multi_op_handle);
	ASSERT(success >= 0);

	CALLOC(request,MPI_Request,j_nprocs);

	CALLOC(request_no,int,j_nprocs);

	if (int2dbl_init()) jerror("int2dbl will not work");

#ifdef MESSAGE_LOG
	(void) sprintf(fname,".messages/%.3d.log",jpid);
	if ((logfile = fopen(fname,"w")) == NULL)
		jerror("can't open %s",fname);
	(void) fprintf(logfile,"INITIALISE\n");
	(void) fflush(logfile);
#endif
}

void	pjostle_comm(MPI_Comm *comm)
{

	ASSERT(gs_initialised);

	if (*comm != original_global_comm) {
		MONITOR(success =) MPI_Comm_group(MPI_COMM_WORLD, &global_group);
		ASSERT(success >= 0);
	}
	original_global_comm = global_comm = *comm;
}

int	jinit(int *arg, int *argc, char ***argv)
{	/* jinit: O(1)
	 * initialise the communications for standalone pjostle */
	/* pp 196 */
	int	nparts;
	int	nP;
	int	pid;
MPI_Comm	comm;

	MONITOR(success =) MPI_Init(argc, argv);
	ASSERT(success >= 0);

	/* pp 147-8 */
	MONITOR(success =) MPI_Comm_group(MPI_COMM_WORLD, &global_group);
	ASSERT(success >= 0);

	MONITOR(success =) MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	ASSERT(success >= 0);

	MONITOR(success =) MPI_Group_size(global_group,&nP);
	ASSERT(success >= 0);

	nparts = atoi((*argv)[(*arg)++]);

	pjostle_init(&nP, &pid);

comm = MPI_COMM_WORLD;
pjostle_comm(&comm);

	return(nparts);
}

void	jexcl(int n_idle)
{
	int	range[1][3];

	ASSERT(gs_initialised);

	if (n_idle == 0) {
		global_comm = original_global_comm;
	} else {
		range[0][0] = j_nprocs - n_idle;
		range[0][1] = j_nprocs - 1;
		range[0][2] = 1;

		MONITOR(success =) MPI_Group_range_excl(global_group, 1, range, &non_idle_group);
		ASSERT(success >= 0);

		j_nprocs -= n_idle;

		/* if (jpid < j_nprocs) { */
			MONITOR(success =) MPI_Comm_create(original_global_comm, non_idle_group, &global_comm);
			ASSERT(success >= 0);
		/* } */
	}
}

void	jincl(int n_idle)
{

	ASSERT(gs_initialised);

	if (n_idle) {

		if (jpid < j_nprocs) {
			/* pp 148 */
			MONITOR(success =) MPI_Comm_free(&global_comm);
			ASSERT(success >= 0);
		}

		MONITOR(success =) MPI_Group_free(&non_idle_group);
		ASSERT(success >= 0);

		j_nprocs += n_idle;

	}

	global_comm = original_global_comm;
}

void	jstop(void)
{	/* jstop: O(1)
	 * stops the parallel execution by hanging all processors */
	int	dummy = 0;

	ASSERT(gs_initialised);

	MONITOR(success =) MPI_Ssend(&dummy, 1, MPI_INT, (jpid+1)%j_nprocs, 4000, global_comm);
	ASSERT(success >= 0);
}

void	jfinalise(void)
{	/* jfinalise: O(1)
	 * finalise the communications for pjostle */
	/* pp 197 */

	ASSERT(gs_initialised);

#ifdef MESSAGE_LOG
	(void) fprintf(logfile,"FINALISE\n");
	(void) fclose(logfile);
#endif

	FREE(request_no,int,j_nprocs);

	FREE(request,MPI_Request,j_nprocs);

	/* pp 120 */
	MONITOR(success =) MPI_Op_free(&multi_op_handle);
	ASSERT(success >= 0);

	ASSERT(global_comm == MPI_COMM_WORLD);

	MONITOR(success =) MPI_Group_free(&global_group);
	ASSERT(success >= 0);

	MONITOR(success =) MPI_Finalize();
	ASSERT(success >= 0);
}

void	jsend(int to, int n, int *buffer, int tag)
{	/* jsend: O(?)
	 * blocking send */

	ASSERT(gs_initialised);

#ifdef MESSAGE_LOG
	(void) fprintf(logfile,"send (tag = %d, to = %d, n = %d)\n",tag,to,n);
	print_message(n,buffer);
	(void) fflush(logfile);
#endif

	ASSERT(n == 1 || n == buffer[0]);

	MONITOR(success =) MPI_Send(buffer, n, MPI_INT, to, tag, global_comm);
	ASSERT(success >= 0);
}

void	jrecv_direct(int from, int n, int *buffer, int tag)
{	/* jrecv_direct: O(?)
	 * receive a direct message */
	MPI_Status	status;

	ASSERT(gs_initialised);

	MONITOR(success =) MPI_Recv(buffer, n, MPI_INT, from, tag, global_comm, &status);
	ASSERT(success >= 0);

	ASSERT(n == 1 || n == buffer[0]);

#ifdef MESSAGE_LOG
	(void) fprintf(logfile,"recv (tag = %d, from = 0, n = %d)\n",tag,n);
	print_message(n,buffer);
	(void) fflush(logfile);
#endif
}

void	jasend(int to, int n, int *buffer, int tag, int *rqust)
{	/* jasend: O(?)
	 * non-blocking send */
	/* section 3.7 */

	ASSERT(gs_initialised);

#ifdef MESSAGE_LOG
	(void) fprintf(logfile,"send (tag = %d, to = %d, n = %d)\n",tag,to,n);
	print_message(n,buffer);
	(void) fflush(logfile);
#endif

	ASSERT(n == 1 || n == buffer[0]);

	*rqust = ++request_no[to];
	MONITOR(success =) MPI_Isend(buffer, n, MPI_INT, to, tag, global_comm, &request[to]);
	ASSERT(success >= 0);
}

void	jwait(int to, int rqust)
{	/* jwait: O(?)
	 * wait for a non blocking send to complete */
	MPI_Status	status;

	ASSERT(gs_initialised);

#ifdef MESSAGE_LOG
	(void) fprintf(logfile,"wait (to = %d)\n",to);
	(void) fflush(logfile);
#endif

	ASSERT(rqust == request_no[to]);

	MONITOR(success =) MPI_Wait(&request[to], &status);
	ASSERT(success >= 0);
}

void	jrecv(int from, int n, int *buffer, int tag)
{	/* jrecv: O(?)
	 * receive a message */
	MPI_Status	status;

	ASSERT(gs_initialised);

	MONITOR(success =) MPI_Recv(buffer, n, MPI_INT, from, tag, global_comm, &status);
	ASSERT(success >= 0);

	ASSERT(n == 1 || n == buffer[0]);

#ifdef MESSAGE_LOG
	(void) fprintf(logfile,"recv (tag = %d, from = %d, n = %d)\n",tag,from,n);
	print_message(n,buffer);
	(void) fflush(logfile);
#endif
}

void	jsetup_recvs(int first, int limit, void* processor, int scheduled)
/*ARGSUSED0*/ {}

void	jsetup_recvs_input(int first, int limit, int *nghbr_proc, int scheduled)
/*ARGSUSED0*/ {}

void	jprobe(int *from, int *n, int *type, int tag)
{	/* jprobe: O(?)
	 * look for an arriving message */
	/* pp 48 */
	MPI_Status	status;

	ASSERT(gs_initialised);

#ifdef MESSAGE_LOG
	(void) fprintf(logfile,"probing (tag = %d, from = %d)\n",tag,*from);
	(void) fflush(logfile);
#endif

	if (*from < 0) { /* could be from any source */
		MONITOR(success =) MPI_Probe(MPI_ANY_SOURCE, tag, global_comm, &status);
		ASSERT(success >= 0);
		*from = status.MPI_SOURCE;
	} else {
		if (tag == ANY_TAG) {
			MONITOR(success =) MPI_Probe(*from, MPI_ANY_TAG, global_comm, &status);
			ASSERT(success >= 0);
			*type = status.MPI_TAG;
		} else {
			MONITOR(success =) MPI_Probe(*from, tag, global_comm, &status);
			ASSERT(success >= 0);
		}
	}
	MONITOR(success =) MPI_Get_count(&status, MPI_INT, n);
	ASSERT(success >= 0);
#ifdef MESSAGE_LOG
	(void) fprintf(logfile,"probe (tag = %d, from = %d, n = %d)\n",tag,*from,*n);
	(void) fflush(logfile);
#endif
}

int	*jglobal(int n, int *buffer, int tag)
{	/* jglobal: O(G*I)
	 * global reduction operation */
	/* pp 90-93, 111-123 */
	int	recv_value;
	int	*recvbuf;

	ASSERT(gs_initialised);

#ifdef MESSAGE_LOG
	(void) fprintf(logfile,"\nglobal (tag = %d)\n",tag);
	(void) fflush(logfile);
#endif

	ASSERT(tag != MULTI || n == buffer[0]);

	if (tag == MULTI) {
		recvbuf = &buffer[n];
		recvbuf[0] = n;
		recvbuf += 1;
		buffer += 1;
		n -= 1;
	} else if (n > 1) {
		recvbuf = &buffer[n];
	} else {
		recvbuf = &recv_value;
	}

	switch (tag) {
	case MAX:
		MONITOR(success =) MPI_Allreduce(buffer, recvbuf, n, MPI_INT, MPI_MAX, global_comm);
		ASSERT(success >= 0);
		break;
	case MIN:
		MONITOR(success =) MPI_Allreduce(buffer, recvbuf, n, MPI_INT, MPI_MIN, global_comm);
		ASSERT(success >= 0);
		break;
	case SUM:
		MONITOR(success =) MPI_Allreduce(buffer, recvbuf, n, MPI_INT, MPI_SUM, global_comm);
		ASSERT(success >= 0);
		break;
	case MULTI:
		MONITOR(success =) MPI_Allreduce(buffer, recvbuf, n, MPI_INT, multi_op_handle, global_comm);
/*
		MONITOR(success =) MPI_Reduce(buffer, recvbuf, n, MPI_INT, multi_op_handle, 0, global_comm);
		MONITOR(success =) MPI_Bcast(recvbuf, n, MPI_INT, 0, global_comm);
*/
		ASSERT(success >= 0);
		break;
	}

	if (tag == MULTI) {
		recvbuf -= 1;
		n += 1;
	} else if (n == 1) {
		*buffer = recv_value;
		recvbuf = buffer;
	}

#ifdef MESSAGE_LOG
	print_message(n,recvbuf);
	(void) fflush(logfile);
#endif

	return(recvbuf);
}

void	jabort(void)
{	/* jabort: O(1)
	 * abort a parallel run */
	/* pp 197 */

	ASSERT(gs_initialised);

#ifdef MESSAGE_LOG
	(void) fprintf(logfile,"ABORT\n");
	(void) fflush(logfile);
	(void) fclose(logfile);
#endif

	MONITOR(success =) MPI_Abort(MPI_COMM_WORLD, 1);
	ASSERT(success >= 0);
#ifndef PACKAGE
	(void) system("/home/wc06/bin/nkill pjostle");
#endif
}

void	jbarrier(void)
{	/* jbarrier: O(1)
	 * barrier until all process reach this point */

	if (j_nprocs == 1) return;

	ASSERT(gs_initialised);

	MONITOR(success =) MPI_Barrier(global_comm);
	ASSERT(success >= 0);
}

double	jtime(void)
{	/* jtime: O(?)
	 * parallel timer returning wall clock time */
	/* pp 195 */

	ASSERT(gs_initialised);

#ifdef T3D
	return((double) 6.6e-9 * rtclock());
#else
	return(MPI_Wtime());
#endif
}

