/*
// @HEADER
// ***********************************************************************
// 
//                Pliris: Parallel Dense Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

Title:
  LU factorization under the MPI message passing protocol
 
General Description:
  This code contains all the routines to do LU factorization with
  pivoting double precision real matrices.  The matrix L is not
  available explicitly, but rather as a product of Gauss matrices and
  permutation matrices.  Routines for triangular solves are included.
  Results are reported for the linpack benchmark

Files Required:
  BLAS_prototypes.h
  README
  factor.c           factorization of matrix
  factor.h
  init.c             initialization of matrix and right hand side
  init.h
  macros.h           macros for processor assignment and decomposition
  pcomm.c            paragon communication routines
  pcomm.h
  mat.c              main program
  solve.c            triangular solve
  solve.h
  time.c             timing routines

Modification History:
  09/30/93 - initial port to Intel Paragon completed
  10/28/93 - linpack residual checks incorporated
  10/29/93 - out-of-core facilities removed from code
  04/04/97 - message passing changed to MPI, ported to TFLOPS
Authors: 
  1. David E. Womble
     Sandia National Laboratories
     P.O. Box 5800
     Albuquerque, NM 87185-5800
 
     dewombl@cs.sandia.gov
     (505) 845-7471

  2. David S. Greenberg
 
  3. Joe Kotulski (primary point of contact)
     Sandia National Laboratories
     P. O. Box 5800
     Albuquerque, NM 87185-1166

     jdkotul@sandia.gov
     (505)-845-7955
     (Ported the solver to the Tflop machine, MPI changed to pass bytes
       for the matrix entries)
 
  4.Brian Driessen
Authors' conditions of release:
  1. The authors receive acknowledgement of authorship in any publication 
     or in any distributions of code in which this code or any part of this
     code is used;
  2. The authors is notified of any problems encountered with the
     code;
*/

#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "defines.h"
#include "init.h"
/* double  seconds(double start);
double  timing(double secs, int type);  */
#include "factor.h"
#include "solve.h"
#include "permute_rhs.h"
#include "permute_mat.h"
#include "forward.h"
#include "macros.h"
#include "pcomm.h"
#include "block.h"

#define TIMERTYPE1 2048
#define TIMERTYPE2 2049
#define TIMERTYPE3 2050
#define TIMERTYPE4 2051
#define TIMERTYPE5 2052
#define TIMERTYPE6 2053
#define TIMERTYPE7 2054
#define TIMERTYPE8 2055
#define TIMERTYPE9 2056

#define BCASTTYPE 2060

#define RESIDTYPE0 2070
#define RESIDTYPE1 2071
#define RESIDTYPE2 2072

#define NUMARGS 8
#define MINHEAP 14000000

/* the number of integer arguments passed by node 0 to rest of nodes 
   currently are matrix size, row cube size, segment size, flags, 
   BLAS blocksize, number of disk controllers */

int   dsg;
int   me, me_proc, host;	/* processor id information */

int   nprocs_cube;		/* num of procs in the allocated cube */
int   nprocs_row;		/* num of procs to which a row is assigned */
int   nprocs_col;		/* num of procs to which a col is assigned */
int   max_procs;                /* max num of procs in any dimension */

int  myrow,mycol;
MPI_Comm row_comm,col_comm;


int   nrows_matrix;		/* number of rows in the matrix */
int   ncols_matrix;		/* number of cols in the matrix */
int   matrix_size;		/* order of matrix=nrows_matrix=ncols_matrix */
int   DISKS;		        /* number of disk controllers */

int   my_first_row;		/* proc position in a row */
int   my_first_col;		/* proc position in a col */

int   my_rows;			/* num of rows I own */
int   my_cols;			/* num of cols I own */

int   nrhs;                     /* number of right hand sides in the matrix */
int   my_rhs;                   /* number of right hand sides that I own */

int   nsegs_row;		/* num of segs to which each row is assigned */
int   ncols_seg;		/* num of cols in each segment */
int   ncols_last;		/* num of cols in the rightmost segment */

int i_one = 1;

int   my_cols_seg;		/* num of cols I own in a seg */
int   my_cols_last;		/* num of cols I own in the rightmost seg */

int   mat_stride;               /* stride to second dimension of mat */
int   col1_stride;              /* stride to second dimension of col1 */
int   row1_stride;              /* stride to second dimension of row1 */
int   size_of_entry;            /* number of bytes in a double or complex */

int   *pivot_vec;               /* stores pivot information */

DATA_TYPE *mat_tmp;
int kstart;
DATA_TYPE *mat;			/* incore storage for col being factored */
DATA_TYPE *rhs;                 /* storage for right hand sides */
DATA_TYPE *rhs_copy;            /* store the original rhs for comparison */
DATA_TYPE *col1,*col1a;	        /* ptrs to col used for updating a col */
DATA_TYPE *col2;		/* ptr to col received in message buf */
DATA_TYPE *row1,*row1a;		/* ptr to diagonal row */
DATA_TYPE *row2;		/* ptr to pivot row */
DATA_TYPE *row3;                /* ptr to row used for pivots */
double    resid, tresid;        /* residuals */

int zpflag;

int   blksz;			/* block size for BLAS 3 operations */
int   rhs_blksz;                /* agglomeration block size for backsolve */
int   colcnt;			/* number of columns stored for BLAS 3 ops */
int   totmem = 0;		/* Total memory (heap) used */
int k,kk,kkk,checkrow,checkcol;

#ifdef COMPLEX
  volatile int   MSPLIT=100;      /* ZGEMM splitting parameter */
#else
  volatile int   MSPLIT=20 ;      /* DGEMM splitting parameter */
#endif

volatile int goflag=0;
volatile int goflag1=0;

extern double seconds( double);
extern double timing(double , int);

void 
main(int argc, char *argv[])
{
    extern char *optarg;
    int   error;
    char  ch;
    static int buf[NUMARGS];
    int   nnodes;		/* number of nodes */
    int   mlen;
    int   i, j;

    double normA1;              /* 1-norm of the matrix */
    double normAinf;            /* inf-norm of the matrix */
    double macheps;             /* machine epsilon */
    double max_error, gmax_error;
    double max_resid, gmax_resid;
    double max_x, gmax_x;

    int   temp;

    double timer0, timer1, timer2, timer3, timer4, timer5, timer6;
    double timer10, timer11, timer12;
    double avg_run, avg_lu, avg_comp, avg_write, avg_read, ops;	
                                /* Calculate MFLOPS */
    DATA_TYPE *p,*ptr1,*ptr2,*ptr3;/* ptr to a segment for printing */
    int *permutations;
  
/*  Initialize into the parallel process */

    /* initialize MPI and find out number of nodes and rank */

    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &me) ;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_cube) ;

    for(i=0;i<NUMARGS;i++) buf[i]=-1;
    error = 0;

    if (me == 0) {
	/* check command line args */

	while ((ch = getopt(argc, argv, "n:r:s:b:d:v:m:CPRW")) != EOF) {
	    switch (ch) {
		case 'n':
		    sscanf(optarg, "%d", &buf[0]);
		    break;
		case 'r':
		    sscanf(optarg, "%d", &buf[1]);
		    break;
		case 'b':
		    sscanf(optarg, "%d", &buf[4]);
		    break;
		case 'v':
		    sscanf(optarg, "%d", &buf[6]);
		    break;
		case 'm':
		    sscanf(optarg, "%d", &buf[7]);
		    break;
		default:
		    error = 1;
	    }
	}

	if (error) {
	    fprintf(stderr, 
                "Usage: %s [-n ms] [-r pr] [-v rhs] [-b blk] [-m split]\n",
		argv[0]);
	    fprintf(stderr, "       ms    matrix size\n");
	    fprintf(stderr, "       pr    processors per row\n");
	    fprintf(stderr, "       rhs   number of right hand sides\n");
	    fprintf(stderr, "       blk   BLAS3 block size\n");
	    fprintf(stderr, "       split cache blocking size\n");
	    buf[0] = -1;
	    exit(-1);
	}
        fprintf(stderr, "  ----  MPI version ----\n");
	if (buf[0] < 0) {
	    fprintf(stderr, "Enter size of matrix: ");
	    scanf("%d", &buf[0]);
	}
	if (buf[1] < 0) {
	    fprintf(stderr, 
              "Enter number of processors to which each row is assigned ");
	    scanf("%d", &buf[1]);
	}
        if ( (nprocs_cube/buf[1]) * buf[1] != nprocs_cube)
        {
     		if (me == 0)
     		{
       			printf("nprocs_row must go into numprocs perfectly!\n");
       			printf("Try a different value of nprocs_row.\n");
     		}
     		MPI_Barrier(MPI_COMM_WORLD);
     		exit(0);
    	}

        buf[2] = buf[0];
	if (buf[4] < 0) {
            buf[4] = DEFBLKSZ;
	}
	if (buf[6] < 0) {
            buf[6] = 1;
	}
	if (buf[6] > 1) {
            buf[6] = 1;
	}
        if (buf[7] < 0) {
          buf[7] = MSPLIT;
        }
    }
    mlen = NUMARGS * sizeof(int);

    /* Send the initilization data to each processor    */

    MPI_Bcast(buf,mlen,MPI_CHAR,0,MPI_COMM_WORLD);
    if (buf[0] <= 0) {
	/* Node 0 received garbage input and terminated */

	exit(-1);
    }

    blksz = buf[4];

    matrix_size = buf[0];
    nrows_matrix = matrix_size;
    ncols_matrix = matrix_size;

    nprocs_row = buf[1];
    nprocs_col = nprocs_cube/nprocs_row;
    max_procs = (nprocs_row < nprocs_col) ? nprocs_col : nprocs_row;

    /* set up communicators for rows and columns */

    myrow = mesh_row(me);
    mycol = mesh_col(me);
    MPI_Comm_split(MPI_COMM_WORLD,myrow,mycol,&row_comm);
    MPI_Comm_split(MPI_COMM_WORLD,mycol,myrow,&col_comm);
    MPI_Comm_rank(row_comm,&checkcol);
    MPI_Comm_rank(col_comm,&checkrow);
    if (myrow != checkrow ||  mycol != checkcol)
    {
	printf("myrow or mycol error!\n");fflush(NULL);
	exit(0);
    } 

    MPI_Barrier(MPI_COMM_WORLD);

    {int checkcol,checkrow;
     MPI_Comm_rank(col_comm, &checkrow) ;
     MPI_Comm_rank(row_comm, &checkcol) ;
     if (myrow != checkrow)
       printf("Node %d: my row = %d but rank in col = %d\n",me,myrow,checkrow);
     if (mycol != checkcol)
       printf("Node %d: my col = %d but rank in row = %d\n",me,mycol,checkcol);
   }

    /* Set up the matrix distribution */

    my_first_col = mesh_col(me);
    my_first_row = mesh_row(me);

    my_rows = nrows_matrix / nprocs_col;
    if (my_first_row < nrows_matrix % nprocs_col)
	++my_rows;
    my_cols = ncols_matrix / nprocs_row;
    if (my_first_col < ncols_matrix % nprocs_row)
	++my_cols;

    nrhs = buf[6];

    blksz = buf[4];
    MSPLIT = buf[7];

    ncols_seg =  ncols_matrix;
    nsegs_row = 1;
    my_cols_seg = my_cols;
    ncols_last = ncols_matrix;
    my_cols_last = my_cols;

    nrhs = (nrhs < ncols_seg) ? nrhs : ncols_seg;

    /* Distribute the RHS per processor */

    my_rhs = nrhs / nprocs_row;
    if (my_first_col < nrhs % nprocs_row) ++my_rhs;
    
    /* Test case for 1 RHS */

    my_rhs = 1;
    dsg = my_rhs;
    
    /* Compute operation counts */
    
#ifdef COMPLEX
    ops = 8.0 * matrix_size * matrix_size * matrix_size / 3.0;
    ops += nrhs * 8.0 * matrix_size * matrix_size;
#else
    ops = 2.0 * matrix_size * matrix_size * matrix_size / 3.0;
    ops += nrhs * 2.0 * matrix_size * matrix_size;
#endif

    if (me == 0) {
      fprintf(stderr, "\nMatrix size: %d x %d\n", nrows_matrix, ncols_matrix);
#ifdef ZCPLX
      fprintf(stderr, "Data type is double precision, complex\n");
#endif
#ifdef SCPLX
      fprintf(stderr, "Data type is single precision, complex\n");
#endif

#ifdef DREAL 
      fprintf(stderr, "Data type is double precision, real\n");
#endif
#ifdef SREAL
      fprintf(stderr, "Data type is single precision, real\n");
#endif
      fprintf(stderr, "Number of right hand sides is %d\n",nrhs);
      fprintf(stderr, "Block size is %d\n", blksz);
      fprintf(stderr, "Number of processors: %d\n", nprocs_cube);
      fprintf(stderr, 
              "Number of processors used: %d\n", nprocs_row*nprocs_col);
      fprintf(stderr, 
              "Processor decomposition is %d x %d\n", nprocs_col, nprocs_row);
      fprintf(stderr, "Using libkmath gemms\n");
    }

    /* Initialize the memory that needs that needs to be used       */

    permutations = (int *)malloc(nrows_matrix * sizeof(int));
    totmem += (nrows_matrix)*sizeof(int);
    if (permutations == NULL) {
       fprintf(stderr, "Node %d: Out of memory in setting the Permutations \n",me);
    }

    if (nrhs > 0) {
      temp = (my_rhs > 1) ? my_rhs : 1;
      rhs_copy = (DATA_TYPE *) malloc((my_rows+1)*temp*sizeof(DATA_TYPE));
      totmem += (my_rows+1)*temp*sizeof(DATA_TYPE);
      if (rhs_copy == NULL) {
	fprintf(stderr, "Node %d: Out of memory in setting the RHS\n", me);
	exit(-1);
      }
    }

    pivot_vec = (int *) malloc(my_cols_last * sizeof(int));
    totmem += my_cols_last * sizeof(int);
    if (pivot_vec == NULL) {
	fprintf(stderr, "Node %d: Out of memory for the pivot vector\n", me);
	exit(-1);
    }

    row3 = (DATA_TYPE *) malloc((my_cols_last + (2*blksz) + dsg) 
                                    * sizeof(DATA_TYPE));
    totmem += (my_cols_seg + (2*blksz) + dsg) * sizeof(DATA_TYPE);
    if (row3 == NULL) {
	fprintf(stderr, "Node %d: Out of memory for row3 vector\n", me);
	exit(-1);
    }

    row2 = (DATA_TYPE *) malloc((my_cols_last + (2*blksz) + dsg) 
                                    * sizeof(DATA_TYPE));
    totmem += (my_cols_seg + (2*blksz) + dsg) * sizeof(DATA_TYPE);
    if (row2 == NULL) {
	fprintf(stderr, "Node %d: Out of memory for row2 vector\n", me);
	exit(-1);
    }

    col2 = (DATA_TYPE *) malloc((my_rows + 1) * sizeof(DATA_TYPE));
    totmem += (my_rows + 1) * sizeof(DATA_TYPE);
    if (col2 == NULL) {
	fprintf(stderr, "Node %d: Out of memory for col2 vector\n", me);
	exit(-1);
    }

    /* this is a two dimensional array, its leading dimension is 
       coded into routines change with care */

    col1 = (DATA_TYPE *) malloc(my_rows*blksz* sizeof(DATA_TYPE));
    col1_stride = my_rows;
    totmem += my_rows*blksz* sizeof(DATA_TYPE);
    if (col1 == NULL) {
      fprintf(stderr, "Node %d: Out of memory for col1 vector\n", me);
      exit(-1);
    }

    col1a = (DATA_TYPE *) malloc(my_rows*blksz* sizeof(DATA_TYPE));
    col1_stride = my_rows;
    totmem += my_rows*blksz* sizeof(DATA_TYPE);
    if (col1a == NULL) {
      fprintf(stderr, "Node %d: Out of memory for cola vector\n", me);
      exit(-1);
    }

    /* this is a two dimensional array, its leading dimension is 
       coded into routines change with care */

    row1 = (DATA_TYPE *) malloc(blksz*(my_cols_last+blksz+dsg)
                                   * sizeof(DATA_TYPE));
    row1_stride = my_cols_last+blksz+dsg;
    totmem += blksz * (my_cols_last + blksz + dsg) * sizeof(DATA_TYPE);
    if (row1 == NULL) {
	fprintf(stderr, "Node %d: Out of memory for row1 vector\n", me);
	exit(-1);
    }

    row1a = (DATA_TYPE *) malloc(blksz*(my_cols_last+blksz+dsg)
                                   * sizeof(DATA_TYPE));
    row1_stride = my_cols_last+blksz+dsg;
    totmem += blksz * (my_cols_last + blksz + dsg) * sizeof(DATA_TYPE);
    if (row1a == NULL) {
	fprintf(stderr, "Node %d: Out of memory for row1a\n", me);
	exit(-1);
    }

    /* this is a two dimensional array, its leading dimension is 
       coded into routines change with care */

    mat = (DATA_TYPE *) malloc(my_rows * (my_cols_last+dsg) 
                                  * sizeof(DATA_TYPE));
    mat_stride = my_rows;
    totmem += my_rows * (my_cols_last+dsg) * sizeof(DATA_TYPE);
    if (mat == NULL) {
      fprintf(stderr, "Node %d: Out of memory for mat vector\n", me);
      exit(-1);
    }


    if (me == 0) {
	fprintf(stderr, "Memory allocated from heap: %d\n", totmem);
    }

    /*
     * the matrix will be stored on a per processor basis.  The processor is
     * responsible for its portion of each segment. The size of its portion
     * of each segment is in bytes[] so its total file size is the sum of the
     * entries of bytes[].
     */

    /* Initialize the matrix */

    timer0 = (double) 0.0;
    timer1 = (double) 0.0;

    MPI_Barrier(MPI_COMM_WORLD);

    timer0 = seconds(timer0);
    init_seg(mat, 0); 

    /* Initialize the RHS and pack into the matrix */
     
    if ((nrhs > 0)) {
      rhs = mat + my_rows*my_cols_last;
      init_rhs(rhs, mat, 0);
      XCOPY(my_rows, rhs, i_one, rhs_copy, i_one);
    }
    timer0 = seconds(timer0);

    /* initialize communication routines */

    initcomm();
    
    if (me == 0)
	fprintf(stderr, "Initialization time = ");
    timing(timer0, TIMERTYPE1);

    if (me == 0)
	fprintf(stderr, "\n");

    timer0 = (double) 0.0;
    timer2 = (double) 0.0;
    timer4 = (double) 0.0;
    timer6 = (double) 0.0;

    MPI_Barrier(MPI_COMM_WORLD);

    if ((me == 0)) {
      fprintf(stderr, "Beginning Factor\n");
    }
    timer0 = seconds(timer0);

    timer4 = seconds(timer4);
    /* Factor and forward solve the matrix */

    factor(mat);
    timer4 = seconds(timer4);
    /* Insert the actual right hand side that we want to solve
        into "rhs" (in this case it is rhs_copy): */
    /*    XCOPY(my_rows, rhs_copy, i_one, rhs, i_one);  */
    /* Forward solve the matrix for a right hand side of rhs_copy.  
	Note:  The solution will be stored in rhs spread across
	the processors that share column 0 */
    /* forwardall(mat,permutations,rhs_copy,rhs);  */
    /* solve system if there are right hand sides */

    if ((nrhs > 0)) {
      if ((me == 0)) {
        fprintf(stderr, "Beginning solve\n");
      }
      timer6 = seconds(timer6);
      back_solve6(mat, rhs); 
      timer6 = seconds(timer6);
    }

    timer0 = seconds(timer0);

    if (me == 0)
	fprintf(stderr, "\nTotal run time = ");
    avg_run = timing(timer0, TIMERTYPE3);
    if (me == 0)
	fprintf(stderr, "Time in factor() = ");
    avg_lu = timing(timer4, TIMERTYPE7);
    if (me == 0)
	fprintf(stderr, "Time in solve() = ");
    timing(timer6, TIMERTYPE8);
    if (me == 0)
	fprintf(stderr, "\nMFLOPS per node = %f, total = %f\n",
		((ops / avg_run) / 1000000.0) / nprocs_cube,
		(ops / avg_run) / 1000000.0);
    if (me == 0)
	fprintf(stderr, "MFLOPS per node in factorization = %f, total = %f\n",
		((ops / avg_lu) / 1000000.0) / nprocs_cube,
		(ops / avg_lu) / 1000000.0);


    /* Compute quantities for the linpack check */

     MPI_Bcast(rhs,my_rows,MPI_DOUBLE,0,row_comm);

     /* nrhs = 0; */   /* Circumvent the checking  */
    if ((nrhs > 0) ) {

    /* Initialize the matrix   */
      
      init_seg(mat, 0);

    /* Perform the matrix vector multiplication Ax    */
 
      mat_vec(mat, 0, rhs);

      ptr3 = row1;
      max_error = (double) 0.0;
      if (myrow == 0) {
        for (i=0; i<my_cols; i++) {
#ifdef COMPLEX
          tresid = fabs((*ptr3).r - (double) 1.0);
          max_error = (tresid > max_error) ? tresid : max_error;
          tresid = fabs((*ptr3).i - (double) 0.0);
          max_error = (tresid > max_error) ? tresid : max_error;
#else
          tresid = fabs(*ptr3 - (double) 1.0);
          max_error = (tresid > max_error) ? tresid : max_error;
#endif
          ptr3++;
        }
      }
  
      /* Find maximum error across the processors  */
      gmax_error=max_error;
      MPI_Allreduce(&max_error,&gmax_error,1,MPI_DATA_TYPE,MPI_MAX,MPI_COMM_WORLD);
      ptr1 = rhs;
      ptr2 = rhs_copy;
      max_resid = (double) 0.0;
      if (me == col_owner(0)) {
        for (i=0; i<my_rows; i++) {
#ifdef COMPLEX
          tresid = fabs((*ptr1).r - (*ptr2).r);
          max_resid = (tresid > max_resid) ? tresid : max_resid;
          tresid = fabs((*ptr1).i - (*ptr2).i);
          max_resid = (tresid > max_resid) ? tresid : max_resid;
#else
          tresid = fabs(*ptr1 - *ptr2);
          max_resid = (tresid > max_resid) ? tresid : max_resid;
#endif
          ptr1++;
          ptr2++;
        }
      }

      /* Find max residual */
 
      MPI_Allreduce(&max_resid,&gmax_resid,1,MPI_DATA_TYPE,MPI_MAX,
		MPI_COMM_WORLD);

      max_x = (double) 0.0;    /* norm(x,inf) */
      ptr3 = row1;
      if (me == row_owner(0)) {
        for (i=0; i<my_cols; i++){
#ifdef COMPLEX
          tresid = fabs((*ptr3).r);
          max_x = (tresid > max_x) ? tresid : max_x;
          tresid = fabs((*ptr3).i);
          max_x = (tresid > max_x) ? tresid : max_x;
#else
          tresid = fabs(*ptr3);
          max_x = (tresid > max_x) ? tresid : max_x;
#endif
          ptr3++;
        }
      }

      /* Find maximum norm(x,inf) */

      MPI_Allreduce(&max_x,&gmax_x,1,MPI_DATA_TYPE,MPI_MAX,
	MPI_COMM_WORLD);
      normA1 = one_norm(mat, 0);
      normAinf = inf_norm(mat, 0);
      macheps = init_eps();

      if (me == 0){
        fprintf(stderr, "\nCheck error in solution\n");
        /* fprintf(stderr, "  Norm(error,inf) = %12.6e\n",gmax_error); */
        fprintf(stderr, "  Norm(resid,inf) = %12.6e\n",gmax_resid); 
        /* fprintf(stderr, "  Norm(x,inf) = %12.6e\n",gmax_x); */
        fprintf(stderr, "  Norm(A,1) = %12.6e\n",normA1);
        fprintf(stderr, "  Norm(A,inf) = %12.6e\n",normAinf);
        fprintf(stderr, "  Machine epsilon = %12.6e\n",macheps);
        fprintf(stderr, "  Size of matrix = %d\n",ncols_matrix);
        fprintf(stderr, "  Norm(resid,inf)/(norm(A,inf)*norm(x,inf)) = %12.6e\n",gmax_resid/normAinf/gmax_x);
        /* fprintf(stderr, "  Norm(error,inf)/norm(x,inf) = %12.6e\n", 
		gmax_error/gmax_x); */
        fprintf(stderr, "  n * norm(A,1) * macheps = %12.6e\n",ncols_matrix*macheps*normA1);

        if ((gmax_resid/normAinf/gmax_x) > (ncols_matrix*macheps*normA1))
          fprintf(stderr,"SOLUTION FAILS RESIDUAL CHECK\n");
        else
          fprintf(stderr,"Solution passes residual check\n");
      }
    }
    free(permutations);

    MPI_Finalize();
}
