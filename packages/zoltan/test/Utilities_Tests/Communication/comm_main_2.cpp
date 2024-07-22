// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "zoltan_types.h"
#include "zoltan_util.h"
#include "zoltan_mem.h"
#include "zoltan_comm_cpp.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

static int out_level = 1;	/* how much output to generate? */

struct Params {
    int       nvals;		/* total number of values on all procs */
    int       same_size;	/* sizes of values (if not variable) */
    int       variable_sizes;	/* are values of variable sizes? */
    int       blocked;		/* should data be blocked by dest proc? */
    int       seed;		/* random number seed */
    float     drop_freq;	/* probability a value is skipped */
};

struct Data {
    int       nvals;		/* total number of values on all procs */
    int       same_size;	/* sizes of values (if not variable) */
    int      *proc_owner;	/* current owner of each value */
    int      *proc_dest;	/* eventual owner of each value */
    int      *sizes;		/* sizes of values (if variable) */
    float    *vals;		/* actual values */
};

struct Answer {
    int       nvals;		/* total number of values on all procs */
    int      *proc_sender;	/* current owner of each value */
    int       same_size;	/* sizes of values (if not variable) */
    int      *sizes;		/* sizes of values (if variable) */
    float    *vals;		/* actual values */
};

void check_comm_info( Zoltan_Comm *plan, struct Data *my_send_data, 
  int nvals_recv, int my_proc);
void print_plan( char *s, Zoltan_Comm *planObject, int my_proc);
void print_data( char *s, struct Data *data);
void free_comm_data( struct Data *data, struct Data *my_data, 
                     struct Answer *true_answer);
void check_comm_answer_reverse( struct Data *my_data, float *reverse_data, int my_proc);
void check_comm_answer( struct Answer *answer, float *recv_data, int my_proc);
void gen_comm_data( struct Params *params, struct Data *data, int nprocs);
void extract_comm_answer( struct Data *data, struct Answer *answer, int my_proc,
  int nprocs);
void set_up_comm_from_send( struct Data *data, struct Data *my_data,
  int    blocked, int    my_proc, int    nprocs);
int read_comm_problem( FILE *in_file, struct Params *params,
                       int *out_level, int my_proc);
double drandom();

int main(int argc, char *argv[])
{
    FILE     *in_file = NULL;	/* file with data for problems */
    Zoltan_Comm *plan;	        /* communication object*/
    struct Params params;	/* parameters describing a problem */
    struct Data data;		/* data describing a problem instance */
    struct Data my_send_data;	/* data I initially own */
    struct Answer true_answer;	/* expected outcome of exchange */
    int       nvals_recv;	/* number of vals I own after exchange */
    float    *recv_data;	/* values I own after exchange */
    float    *reverse_data;	/* values I own after reversing communication */
    char      file_name[100];	/* name of input file */
    int       nbytes;		/* size of objects */
    int       my_proc;		/* my processor ID */
    int       nprocs;		/* total number of processors */
    int       flag;		/* return code from comm operations */
    int       more_problems;	/* are there more problems to do? */
    int       i, j;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    Zoltan_Memory_Debug(2);

    if (argc > 1) strcpy(file_name, argv[1]);
    else strcpy(file_name,"comm_input.dat");
    
    if (my_proc == 0) {
	in_file = fopen(file_name, "r");
	if (in_file == NULL) {
	    printf("No input file `%s' found.\n", file_name);
	}
    }

    /* Read some problem descriptors from a file */
    more_problems = read_comm_problem(in_file, &params, &out_level, my_proc);

    while (more_problems) {

	/* Generate full data at random on each proc */
	gen_comm_data(&params, &data, nprocs);

if (out_level > 2) if (my_proc == 0) print_data("DATA", &data);

	/* Figure out from the data what what to expect to receive */
	extract_comm_answer(&data, &true_answer, my_proc, nprocs);

	/* Extract what to send */
	set_up_comm_from_send(&data, &my_send_data, params.blocked, my_proc, nprocs);

if (out_level > 2) print_data("MY_DATA", &my_send_data);

if (out_level > 1) printf("%d: About to call comm_create\n", my_proc);
	/* Call comm routines */
        plan = new Zoltan_Comm(my_send_data.nvals, my_send_data.proc_dest,
            MPI_COMM_WORLD, 1, &nvals_recv);

if (out_level > 1) printf("%d: About to call comm_info\n", my_proc);
        check_comm_info(plan, &my_send_data, nvals_recv, my_proc);

if (out_level > 2) print_plan("BEFORE RESIZE", plan, my_proc);
	/* "4" reflects the max_sizes value in gen_comm_data */
	recv_data = (float *) ZOLTAN_MALLOC(nvals_recv * 4 * sizeof(float));
	reverse_data = (float *) ZOLTAN_MALLOC(my_send_data.nvals * 4 * sizeof(float));

	if (my_send_data.sizes != NULL) {
if (out_level > 1) printf("%d: About to call comm_resize\n", my_proc);
	    plan->Resize(my_send_data.sizes, 43, NULL);
	    plan->Resize(NULL, 43, NULL);
	    plan->Resize(my_send_data.sizes, 43, NULL);
if (out_level > 2) print_plan("AFTER RESIZE", plan, my_proc);
	}

	if (my_send_data.sizes == NULL) nbytes = my_send_data.same_size * sizeof(float);
	else nbytes = sizeof(float);

if (out_level > 1) printf("%d: About to call comm_do\n", my_proc);
	flag = plan->Do(2, (char *) my_send_data.vals, nbytes,
	    (char *) recv_data);

	if (flag == ZOLTAN_OK) {
if (out_level > 1) printf("%d: About to call check_answer\n", my_proc);
	    /* Check answers */
	    check_comm_answer(&true_answer, recv_data, my_proc);
	}
	else {
	    printf("%d: Comm_Do returned error code %d\n", my_proc, flag);
	}

if (out_level > 1) printf("%d: About to call comm_do_reverse\n", my_proc);

        ZOLTAN_COMM_OBJ *planStruct = plan->Get_ZOLTAN_COMM_OBJ();

	i = (true_answer.sizes != NULL && planStruct->indices_to != NULL);
	MPI_Allreduce(&i, &j, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (j == 0) flag = plan->Do_Reverse(2, (char *) recv_data,
	    nbytes, true_answer.sizes, (char *) reverse_data);
	else {
	    if (my_proc == 0)
		printf(">> Non-blocked, variable-sized recvs not supported\n");
	    flag = ZOLTAN_FATAL;
	}

	if (flag == ZOLTAN_OK) {
if (out_level > 1) printf("%d: About to call check_answer_reverse\n", my_proc);
	    /* Check answers */
	    check_comm_answer_reverse(&my_send_data, reverse_data, my_proc);
	}
	else {
if (out_level > 1) printf("%d: Comm_Do_Reverse returned error code %d\n", my_proc, flag);
	}


	/* Free up data structures */
	ZOLTAN_FREE((void *) &reverse_data);
	ZOLTAN_FREE((void *) &recv_data);

	free_comm_data(&data, &my_send_data, &true_answer);

        delete plan;

	/* Read some problem descriptors from a file */
	more_problems = read_comm_problem(in_file, &params, &out_level, my_proc);

    }

    Zoltan_Memory_Stats();
    MPI_Finalize();

    return(0);
}


int read_comm_problem(
FILE *in_file,
struct Params *params,
int *level,
int my_proc)
{
    char  line[121];
    char *cptr;
    int   in_level;
    int   i;
    int   flag;

    /* Skip comment lines, then read set of ints */

    if (my_proc == 0) {
	flag = 1;
	line[0] = '%';
	while (line[0] == '%' || line[0] == '#') {
	    cptr = fgets(line, 120, in_file);
	    if (cptr == NULL) {
		flag = 0;
		break;
	    }
	}
	if (flag) {
	    i = sscanf(line, "%d%d%d%d%d%f%d", &params->nvals, &params->same_size,
	         &params->variable_sizes, &params->blocked, &params->seed,
		 &params->drop_freq, &in_level);
	    if (i != 6 && i != 7) {
	        printf("Input error, text is '%s'\n", line);
	        flag = 0;
	    }
	    if (i == 7) {
		*level = in_level;
	    }
	}
        if (flag == 0) params->nvals = -1;
    }

    MPI_Bcast((void *) params, sizeof(struct Params), MPI_BYTE, 0,
	MPI_COMM_WORLD);
    MPI_Bcast((void *) level, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (params->nvals == -1) return(0);
    else return(1);
}


void gen_comm_data(
struct Params *params,
struct Data *data,
int nprocs)
{
    int       sizes_max = 4;	/* maximum size of variable sized data */
    int       index;		/* pointer into vals array */
    int       i, j;		/* loop counters */

    data->proc_owner = (int *) ZOLTAN_MALLOC(params->nvals * sizeof(int));
    data->proc_dest = (int *) ZOLTAN_MALLOC(params->nvals * sizeof(int));
    if (params->variable_sizes) {
	data->sizes = (int *) ZOLTAN_MALLOC((params->nvals + 1) * sizeof(int));
	data->vals = (float *)
	   ZOLTAN_MALLOC(params->nvals * sizes_max * sizeof(float));
    }
    else {
	data->sizes = NULL;
	data->same_size = params->same_size;
	data->vals = (float *)
	   ZOLTAN_MALLOC(params->nvals * data->same_size * sizeof(float));
    }

    data->nvals = params->nvals;

    srand(params->seed);
    index = 0;
    for (i = 0; i < params->nvals; i++) {
	data->proc_owner[i] = (int)(drandom() * nprocs);
	data->proc_dest[i] = (int)(drandom() * nprocs);
	if (drandom() < params->drop_freq)
	    data->proc_dest[i] = -1;
	if (params->variable_sizes) {
	    data->sizes[i] = (int)(1 + drandom() * (sizes_max - 1));
	    for (j = 0; j < data->sizes[i]; j++) {
		data->vals[index++] = (int)drandom();
	    }
	}
	else {
	    for (j = 0; j < params->same_size; j++) {
		data->vals[index++] = (int)drandom();
	    }
	}
    }
}


void extract_comm_answer(
struct Data *data,
struct Answer *answer,
int my_proc,
int nprocs)
{
    int       recv_count;	/* number of items I'll recv */
    int       recv_size;	/* total size of data I'll recv */
    int       nfrom, sfrom;	/* count/size of items I'll recv */
    int       i, j, k;		/* loop counters */
    int       ii, jj, kk;	/* loop counters */

    recv_size = 0;
    recv_count = 0;
    for (i = 0; i < data->nvals; i++) {
	if (data->proc_dest[i] == my_proc) {
	    recv_count++;
	    if (data->sizes != NULL)
		recv_size += data->sizes[i];
	}
    }
    if (data->sizes == NULL)
	recv_size = recv_count * data->same_size;

    answer->nvals = recv_count;
    answer->proc_sender = (int *) ZOLTAN_MALLOC(recv_count * sizeof(int));
    if (data->sizes == NULL) {
	answer->same_size = data->same_size;
	answer->sizes = NULL;
    }
    else {
	answer->sizes = (int *) ZOLTAN_MALLOC((recv_count + 1) * sizeof(int));
    }

    answer->vals = (float *) ZOLTAN_MALLOC(recv_size * sizeof(float));

    /* Now copy correct answer into vals.  Use dumb algorithm since */
    /* it's just for testing. */
    kk = 0;
    jj = 0;
    for (i = 0; i < nprocs; i++) {
        ii = 0;
	nfrom = 0;
	sfrom = 0;
	for (j = 0; j < data->nvals; j++) {
	    if (data->proc_dest[j] == my_proc && data->proc_owner[j] == i) {
		/* found next val I'll recv */
		++nfrom;
		if (data->sizes != NULL) {
		    sfrom += data->sizes[j];
		    answer->sizes[kk] = data->sizes[j];
		    for (k = 0; k < data->sizes[j]; k++) {
			answer->vals[jj++] = data->vals[ii + k];
		    }
		}
		else {
		    sfrom += data->same_size;
		    for (k = 0; k < data->same_size; k++) {
			answer->vals[jj++] = data->vals[ii + k];
		    }
		}

		answer->proc_sender[kk++] = i;
	    }
	    if (data->sizes != NULL)
		ii += data->sizes[j];
	    else
		ii += data->same_size;
	}
/*printf("%d will receive %d (%d) from %d\n", my_proc, nfrom, sfrom, i);*/
    }
}


void set_up_comm_from_send(
struct Data *data,
struct Data *my_data,
int    blocked,
int    my_proc,
int    nprocs)
{
    int my_size;	/* total size of data I own */
    int proc;		/* processor number */
    int jstart;		/* saved index pointer */
    int ii, jj, kk;	/* loop counters */
    int i, j, k;	/* loop counters */


    my_data->nvals = 0;
    my_size = 0;
    for (i = 0; i < data->nvals; i++) {
	if (data->proc_owner[i] == my_proc) {
	    ++my_data->nvals;
	    if (data->sizes != NULL) my_size += data->sizes[i];
	}
    }
    if (data->sizes == NULL) {
	my_size = my_data->nvals * data->same_size;
	my_data->sizes = NULL;
	my_data->same_size = data->same_size;
    }
    else {
	my_data->sizes = (int *) ZOLTAN_MALLOC((my_size + 1) * sizeof(int));
    }

    my_data->proc_owner = NULL;
    my_data->proc_dest = (int *) ZOLTAN_MALLOC(my_data->nvals * sizeof(int));
    my_data->vals = (float *) ZOLTAN_MALLOC(my_size * sizeof(float));

    /* Now I can copy data into my_data */
    /* Two cases - data blocked by destination processor or not */

    if (blocked) {	/* Make all data to a proc adjacent */
	jstart = 0;
	ii = 0;
	jj = 0;
	for (i = 0; i < nprocs; i++) {
	    proc = -nprocs - 2;
	    for (j = jstart; j < data->nvals && proc < 0; j++) {
		if (data->proc_owner[j] == my_proc && data->proc_dest[j] >= 0) {
		    proc = data->proc_dest[j];
		    jstart = j;
		}
	    }

	    kk = 0;
	    for (j = 0; j < data->nvals; j++) {
		if (data->proc_owner[j] == my_proc && data->proc_dest[j] == proc) {
		    if (data->sizes != NULL) {
		        my_data->sizes[ii] = data->sizes[j];
			for (k = 0; k < data->sizes[j]; k++ ) {
			    my_data->vals[jj++] = data->vals[kk++];
			}
		    }
		    else {
			for (k = 0; k < data->same_size; k++ ) {
			    my_data->vals[jj++] = data->vals[kk++];
			}
		    }
		    my_data->proc_dest[ii++] = proc;
		    data->proc_dest[j] = -proc - 2;
		}
		else {
		    if (data->sizes != NULL) kk += data->sizes[j];
		    else kk += data->same_size;
		}
	    }
	}
	for (j = 0; j < data->nvals; j++) {	/* Restore dest values */
	    if (data->proc_dest[j] < -1) 
		data->proc_dest[j] = -data->proc_dest[j] - 2;
	}
	/* Don't count negative values since data is being blocked */
	my_data->nvals = ii;
    }

    else {		/* Generate data scattered by processor */
      ii = 0;
      jj = 0;
      kk = 0;
      for (i = 0; i < data->nvals; i++) {
	if (data->proc_owner[i] == my_proc) {
	    if (data->sizes != NULL) {
	        my_data->sizes[ii] = data->sizes[i];
		for (j = 0; j < data->sizes[i]; j++ ) {
		    my_data->vals[jj++] = data->vals[kk++];
		}
	    }
	    else {
		for (j = 0; j < data->same_size; j++ ) {
		    my_data->vals[jj++] = data->vals[kk++];
		}
	    }
	    my_data->proc_dest[ii++] = data->proc_dest[i];
	}
	else {
	    if (data->sizes != NULL) kk += data->sizes[i];
	    else kk += data->same_size;
	}
      }
    }
}


void check_comm_answer(
struct Answer *answer,
float *recv_data,
int my_proc)
{
    int size;			/* size of an object */
    int err;			/* anyone have a problem? */
    int i, ii, j;		/* loop counter */

    err = 0;
    ii = 0;

    for (i = 0; i < answer->nvals; i++) {
	if (answer->sizes == NULL) size = answer->same_size;
	else size = answer->sizes[i];
	for (j = 0; j < size; j++) {
	    if (recv_data[ii] != answer->vals[ii]) {
		if (out_level > 0) printf("%d: *** val %d, item %d, from %d -- Answer = %f, got %f\n",
		    my_proc, i, j, answer->proc_sender[i],
		    answer->vals[ii], recv_data[ii]);
		err = 1;
		goto End;
	    }
	    ii++;
	}
    }	

End:
    MPI_Reduce(&err, &i, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (my_proc == 0) {
	if (i == 0) printf(">> Right answer for forward communication!\n");
        else printf(">> Wrong answer for forward communication on %d procs\n", i);
    }
}


void check_comm_answer_reverse(
struct Data *my_data,
float *reverse_data,
int my_proc)
{
    int size;			/* size of an object */
    int err;			/* any errors? */
    int i, ii, j, jj;		/* loop counter */

    ii = 0;
    jj = 0;
    err = 0;
    for (i = 0; i < my_data->nvals; i++) {
	if (my_data->sizes == NULL) size = my_data->same_size;
	else size = my_data->sizes[i];
	if (my_data->proc_dest[i] >= 0) {
	    for (j = 0; j < size; j++) {
	        if (reverse_data[ii] != my_data->vals[ii]) {
		    if (out_level > 0) printf("%d (R): *** val %d, item %d, via %d -- Answer = %f, got %f\n",
		        my_proc, i, j, my_data->proc_dest[i],
		        my_data->vals[ii], reverse_data[jj]);
		    err = 1;
		    goto End;
	        }
	        ii++;
		jj++;
	    }
	}
	else {
	    ii += size;
	}
    }	

End:
    MPI_Reduce(&err, &i, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (my_proc == 0) {
	if (i == 0) printf(">> Right answer for reverse communication!\n");
        else printf(">> Wrong answer for reverse communication on %d procs\n", i);
    }
}


void free_comm_data(
struct Data *data,
struct Data *my_data,
struct Answer *true_answer)
{
    ZOLTAN_FREE((void *) &data->vals);
    ZOLTAN_FREE((void *) &data->sizes);
    ZOLTAN_FREE((void *) &data->proc_dest);
    ZOLTAN_FREE((void *) &data->proc_owner);

    ZOLTAN_FREE((void *) &my_data->vals);
    ZOLTAN_FREE((void *) &my_data->sizes);
    ZOLTAN_FREE((void *) &my_data->proc_dest);
    ZOLTAN_FREE((void *) &my_data->proc_owner);

    ZOLTAN_FREE((void *) &true_answer->vals);
    ZOLTAN_FREE((void *) &true_answer->sizes);
    ZOLTAN_FREE((void *) &true_answer->proc_sender);
}


double    drandom()
{
    int       val;

    val = rand();

    return (((double) val) / (.01 + RAND_MAX));
}



void print_data(
char *s,
struct Data *data)
{
    int my_proc;
    int i, j, k;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);

    k = 0;
    for (i = 0; i < data->nvals; i++) {
	if (data->proc_owner != NULL) 
	    printf("%d: %s item %d, owner %d, dest = %d, vals =", my_proc, s,
		i, data->proc_owner[i], data->proc_dest[i]);
	else
	    printf("%d: %s item %d, owner %d, dest = %d, vals =", my_proc, s,
		i, my_proc, data->proc_dest[i]);
	if (data->sizes == NULL) {
	    for (j = 0; j < data->same_size; j++) {
		printf(" %f", data->vals[k++]);
	    }
  	
	}
  	else {
	    for (j = 0; j < data->sizes[i]; j++) {
		printf(" %f", data->vals[k++]);
	    }
	}
	printf("\n");
    }
    printf("\n");
}

void print_plan(
char *s,
Zoltan_Comm *planObject,	/* communication object */
int my_proc)
{
    int i;

    ZOLTAN_COMM_OBJ *plan = planObject->Get_ZOLTAN_COMM_OBJ();

    if (!plan) {
      return;
      }

    if (plan->sizes == NULL) {	/* Constant sizes */
        printf("%s %d: Sending (constant size) ", s, my_proc);
	if (plan->indices_to != NULL)
	    printf("(scattered)\n");
	else 
            printf("(blocked)\n");

	for (i = 0; i < plan->nsends + plan->self_msg; i++) {
	    printf("    %d vals to %d, starts_to = %d\n", plan->lengths_to[i],
	        plan->procs_to[i], plan->starts_to[i]);
	    
        }
        printf("\n");

        printf("%s %d: Expecting (constant size) ", s, my_proc);
        if (plan->indices_from != NULL)
	    printf("(scattered)\n");
	else
            printf("(blocked)\n");

        for (i = 0; i < plan->nrecvs + plan->self_msg; i++) {
	    printf("    %d vals from %d, starts_from = %d\n", plan->lengths_from[i],
	        plan->procs_from[i], plan->starts_from[i]);
        }
        printf("\n");
    }
    else {		/* Variable sizes */
        printf("%s %d: Sending (variable sizes) ", s, my_proc);
	if (plan->indices_to != NULL)
	    printf("(scattered)\n");
	else 
            printf("(blocked)\n");

	for (i = 0; i < plan->nsends + plan->self_msg; i++) {
	    printf("    %d vals to %d, starts_to_ptr = %d\n", plan->sizes_to[i],
	        plan->procs_to[i], plan->starts_to_ptr[i]);
	    
        }
        printf("\n");

        printf("%s %d: Expecting (variable sizes) ", s, my_proc);
        if (plan->indices_from != NULL)
	    printf("(scattered)\n");
	else
            printf("(blocked)\n");

        for (i = 0; i < plan->nrecvs + plan->self_msg; i++) {
	    printf("    %d vals from %d, starts_from_ptr = %d\n",
		plan->sizes_from[i], plan->procs_from[i], plan->starts_from_ptr[i]);
        }
        printf("\n");
    }
}

void check_comm_info(
  Zoltan_Comm *plan, 
  struct Data *my_send_data, 
  int nvals_recv,
  int my_proc
)
{
int *info_tmp = NULL,  /* Temp bufs for verifying Zoltan_Comm_Info */
    *info_tmp_send = NULL,
    *info_tmp_recv = NULL;
int *info_send_list = NULL,
    *info_send_procs = NULL,
    *info_send_lengths = NULL;
int  info_nsend = 0,
     info_send_nvals = 0,
     info_send_size = 0,
     info_max_send_size = 0;
int *info_recv_list = NULL,
    *info_recv_procs = NULL,
    *info_recv_lengths = NULL;
int  info_nrecv = 0,
     info_recv_nvals = 0,
     info_recv_size = 0,
     info_total_recv_size = 0;
int  info_self_msg = 0;
int  i;


    plan->Info(&info_nsend, NULL, NULL,
               &info_send_nvals, NULL, NULL, &info_nrecv, NULL,
               NULL, &info_recv_nvals, NULL, NULL, &info_self_msg);

    if (info_send_nvals != my_send_data->nvals)
        printf("%d Error in Zoltan_Comm_Info info_send_nvals %d != %d\n",
             my_proc, info_send_nvals, my_send_data->nvals);

    if (info_recv_nvals != nvals_recv)
        printf("%d Error in Zoltan_Comm_Info info_recv_nvals %d != %d\n",
             my_proc, info_recv_nvals, nvals_recv);

    info_send_size = 2 * (info_nsend + info_self_msg) + info_send_nvals;
    info_send_procs = (int *) ZOLTAN_MALLOC(info_send_size * sizeof(int));
    info_send_lengths = info_send_procs + (info_nsend + info_self_msg);
    info_send_list = info_send_lengths + (info_nsend + info_self_msg);

    info_recv_size = 2 * (info_nrecv + info_self_msg) + info_recv_nvals;
    info_recv_procs = (int *) ZOLTAN_MALLOC(info_recv_size * sizeof(int));
    info_recv_lengths = info_recv_procs + (info_nrecv + info_self_msg);
    info_recv_list = info_recv_lengths + (info_nrecv + info_self_msg);

    plan->Info(&info_nsend, info_send_procs, info_send_lengths,
                &info_send_nvals, &info_max_send_size, info_send_list, 
                &info_nrecv, info_recv_procs, info_recv_lengths,
                &info_recv_nvals, &info_total_recv_size,
                info_recv_list, &info_self_msg);

    for (i = 0; i < info_send_nvals; i++) 
        if (info_send_list[i] != my_send_data->proc_dest[i])
            printf("%d Error in Zoltan_Comm_Info send_list[%d]: %d != %d\n",
                 my_proc, i, info_send_list[i], my_send_data->proc_dest[i]);

    info_tmp = (int *) ZOLTAN_MALLOC((info_send_nvals + info_recv_nvals)
                                      * sizeof(int));
    info_tmp_send = info_tmp;
    info_tmp_recv = info_tmp_send + info_send_nvals;
    for (i = 0; i < info_send_nvals; i++)
        info_tmp_send[i] = my_proc;

    plan->Do(12, (char *) info_tmp_send, sizeof(int), 
              (char *) info_tmp_recv);

    for (i = 0; i < info_recv_nvals; i++)
        if (info_recv_list[i] != info_tmp_recv[i])
            printf("%d Error in Zoltan_Comm_Info recv_list[%d]: %d != %d\n",
                 my_proc, i, info_recv_list[i], info_tmp_recv[i]);

    ZOLTAN_FREE(&info_tmp);
    ZOLTAN_FREE(&info_recv_procs);
    ZOLTAN_FREE(&info_send_procs);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
