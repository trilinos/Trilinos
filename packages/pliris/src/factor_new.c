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
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "defines.h"
#include "BLAS_prototypes.h"

#include "factor.h"
#include "macros.h"
#include "pcomm.h"
#include "mytime.h"


extern int myrow;
extern int mycol; 
extern int me;	               /* processor id information */
extern int nprocs_row;         /* num of procs to which a row is assigned */
extern int nprocs_col;         /* num of procs to which a col is assigned */
extern int nrows_matrix;       /* number of rows in the matrix */
extern int ncols_matrix;       /* number of cols in the matrix */
extern int my_rows;            /* num of rows I own */
extern int my_cols;            /* num of cols I own */
extern int my_rhs;             /* num of right hand side I own */

extern DATA_TYPE *col1;           /* ptrs to col used for updating a col */
extern DATA_TYPE *row1;           /* ptr to diagonal row */
extern DATA_TYPE *row2;           /* ptr to pivot row */
extern DATA_TYPE *row3;           /* ptr to temporary vector for rows */
extern int *pivot_vec;         /* ptr to vector storing list of pivot rows */
extern int mat_stride,col1_stride,row1_stride;  /* strides for 2nd dim of 2d mats */

extern int blksz;              /* block size for BLAS 3 operations */
extern int ringnext,ringprev,ringnex2,ringpre2,ringnex3,ringpre3,ringnex4,ringpre4;
#define LUSTATUSINT 64

extern MPI_Comm col_comm; 

#define LUPIVOTTYPE (1<<19)
#define LUCOLTYPE (1<<20)
#define LUROWTYPE (1<<21)
#define LUPIVROWTYPE ((1<<21) + (1<<20))
#define LUSENDTYPE (1<<22)

#define rowplus(I) proc_num(mesh_row(me),(mesh_col(me)+(I) < nprocs_row) ? mesh_col(me)+(I) : mesh_col(me)+(I)-nprocs_row)
#define rowminus(I) rowplus(nprocs_row - (I))
#define MAXDIST 1

void
factor(DATA_TYPE *seg)
{ int j,k;                        /* loop counter */

  struct pivot_type {           /* a stucture for storing pivot info */
    DATA_TYPE entry;            /*   pivot entry */
    DATA_TYPE current;          /*   current row entry */
    int row;                    /*   pivot row number */
  } pivot;              /* pivot info and some temporary pivot infor */
  
  DATA_TYPE d_zero = CONST_ZERO;/* constant for initializations */
  int one = 1;                  /* constant for the daxpy routines */
  DATA_TYPE d_one = CONST_ONE;  /* constant for the daxpy routines */

  DATA_TYPE d_min_one = CONST_MINUS_ONE; /* constant for the daxpy routines */
  char transA = 'N';            /* all dgemm's don't transpose matrix A  */
  char transB = 'N';            /* nor transpose matrix B */
  int colcnt,cols_in_blk_owned; /* number of columns stored for BLAS 3 ops */
  int c_owner,r_owner,pivot_owner;
  int col_len,row_len,length,row_size,rows_used,cols_used;
  int rel_lpivot_row,lpivot_row,*sav_pivot_ptr;
 
  double pivot_mag;
  DATA_TYPE invpiv;
  DATA_TYPE *cur_col_ptr,*cur_row_ptr,*update_ptr,*piv_row_ptr;
  DATA_TYPE *sav_col_ptr,*sav_row_ptr,*sav_piv_row_ptr;
  DATA_TYPE *cur_col1_row_ptr,*piv_col1_row_ptr;
  DATA_TYPE *temp_row_ptr;
  DATA_TYPE *act_col_ptr,*act_row_ptr,*act_cur_row_ptr,*act_piv_row_ptr;
  int gpivot_row; /* make sure that this is well aligned */

  int dest,ringdist,rdist;
  long type,bytes;


#ifdef DREAL 
  struct {
    double  val;
    int proc;
  } pivot_in,pivot_out;
#endif
#ifdef ZCPLX
  struct {
    double  val;
    int proc;
  } pivot_in,pivot_out;
#endif
#ifdef SREAL
  struct {
    float  val;
    int proc;
  } pivot_in,pivot_out;
#endif
#ifdef SCPLX
  struct {
    float  val;
    int proc;
  } pivot_in,pivot_out;
#endif
  DATA_TYPE entry,current;
  int row;
 
  int numprocs;
 
#ifdef TIMING0
  double updatetime,colupdtime,rowupdtime,scaltime;
  double xpivmsgtime,bcastpivstime,bcastpivrtime,bcastcolstime,bcastcolrtime,bcastcolwtime,bcastrowtime,sendrowtime,recvrowtime;
  double copycoltime,copyrowtime,copyrow1time,copypivrowtime,copypivrow1time;
  double t1,t2;
  double msgtime,copytime,dgemmtime,totaltime;
#endif

  MPI_Request msgrequest;
  MPI_Status msgstatus;
  /* Distribution for the matrix on me */

  MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
  if ( (numprocs/nprocs_row) * nprocs_row != numprocs )
  {
     if (me == 0)
     {
       printf("nprocs_row must go into numprocs perfectly!\n");
       printf("Try a different value of nprocs_row.\n");
     }
     MPI_Barrier(MPI_COMM_WORLD);
     exit(0);
  }
#ifdef TIMING0
  t2 = MPI_Wtime();
#endif
  colcnt = 0;           /* number of column's currently saved for update */
  col_len = my_rows;    /* length of column in remaining local matrix */
  cur_col_ptr = seg;    /* location of first column in remaining local matrix */
  sav_col_ptr = col1;   /* location to store next active column */
  act_col_ptr = col1;   /* location of matrix of columns being saved for dgemm update */

  row_len = my_cols + my_rhs;  /* length of row in local matrix including 
			          rhs's*/

  rows_used = 0;      /* haven't used any local rows yet */
  cols_used = 0;
  cur_row_ptr = seg;  /* location of first row in local matrix */
  cur_col1_row_ptr = col1;  /* location of first row in col1 matrix */
  act_row_ptr = row1; /* location of matrix of rows being saved for dgemm 
			 update */

  sav_piv_row_ptr = row1; /* location for next row being saved for dgemm 
	                 update */

  temp_row_ptr = row3; /* location for pivot row while being sent and 
		         before transposing */
 
  sav_row_ptr = row2;  /* location to save current row and send to 
			 owner of pivot row */

  sav_pivot_ptr = pivot_vec; /* location to store name of pivot row */

  update_ptr = seg;     /* location of remaining local matrix */

#ifdef TIMING0
  xpivmsgtime=bcastpivstime=bcastpivrtime=bcastcolstime=bcastcolrtime=bcastcolwtime=bcastrowtime=sendrowtime=recvrowtime=0.0;
  copycoltime=copyrowtime=copyrow1time=copypivrowtime=copypivrow1time=0.0;
  updatetime=colupdtime=rowupdtime=scaltime=0.0;
#endif

#ifdef PRINT_STATUS
  if (me == 0) {
    printf("Pivots sent separate from column \n");
    printf("Send to %4d neighbors \n",MAXDIST);
    printf("Attempt to do row work while waiting for column to arrive \n");
#ifdef OVERLAP
    printf("Column updates within block done one column at a time \n");
#endif
  }
#endif  

  for (j=0; j<ncols_matrix; j++) {
    c_owner = col_owner(j); r_owner = row_owner(j);
    ringdist = mesh_col(me) - mesh_col(c_owner);
    if (ringdist < 0) ringdist += nprocs_row;
    if (me == c_owner) {
      if (col_len > 0) {

#ifndef OVERLAP	
	if (colcnt != 0){ /* update current column with saved columns */
#ifdef TIMING0
	  t1 = MPI_Wtime();
#endif
          XGEMMS_(&transA, &transB,
                 &col_len, &one, &colcnt, &d_min_one,
                 act_col_ptr, &col1_stride,
                 act_row_ptr, &blksz, &d_one,
                 cur_col_ptr, &mat_stride);

#ifdef TIMING0
	  colupdtime += (MPI_Wtime()-t1);
#endif
	}
#endif
	/* find maximum local pivot */

        pivot.entry = *(cur_col_ptr);
#ifdef CBLAS 
	rel_lpivot_row = IXAMAX(col_len, cur_col_ptr, one);
#else
    /*   Shift for the fortran to C array definitions   */

        rel_lpivot_row = IXAMAX(col_len, cur_col_ptr, one);
	rel_lpivot_row = rel_lpivot_row -1;
#endif
	pivot.entry = *(cur_col_ptr + rel_lpivot_row);
	pivot.row = lrow_to_grow(rows_used+rel_lpivot_row);
	if (me == r_owner) 
          pivot.current = *(cur_col_ptr);  
        else 
          pivot.current = d_zero;
      } else {
	pivot.row = 0;  pivot.entry = d_zero;  pivot.current = d_zero;
      }
#ifdef TIMING0
      t1 = MPI_Wtime();
#endif

       /* Pivot for column j on me */

       entry = pivot.entry;
       row = pivot.row;
       current = pivot.current;
       pivot_in.val = ABS_VAL(pivot.entry);
       pivot_in.proc = me;
       bytes = sizeof(DATA_TYPE);
#ifdef DEBUG
       printf("Node %d: pivot val %g, pivot row %d \n",me,pivot_in.val,pivot.row);
#endif  

       /* Exchange to find global pivot value */
       MPI_Allreduce(&pivot_in,&pivot_out,1,MPI_DATA_TYPE2,MPI_MAXLOC,col_comm);
       /* Changed the mesh_row argument  */
       MPI_Bcast(&current,bytes,MPI_BYTE,mesh_row(r_owner),col_comm);
       MPI_Barrier(col_comm);
       pivot.current = current;
       MPI_Bcast(&row,1,MPI_INT,mesh_row(pivot_out.proc),col_comm);
       MPI_Barrier(col_comm);
       pivot.row = row;
       MPI_Bcast(&entry,bytes,MPI_BYTE,mesh_row(pivot_out.proc),col_comm);
       MPI_Barrier(col_comm);
       pivot.entry = entry;  
  
      
#ifdef TIMING0
      xpivmsgtime += (MPI_Wtime()-t1);
#endif
      *sav_pivot_ptr = pivot.row;
      gpivot_row = pivot.row;
      pivot_mag = ABS_VAL(pivot.entry);
      if (pivot_mag == 0.0) {
        printf("Node %d error -- zero pivot found in column %d -- exiting\n",me,j);
	return; }
      /* divide everything including the diagonal by the pivot entry. */

#ifdef TIMING0
      t1 = MPI_Wtime();
#endif
      INVERSE(pivot.entry,pivot_mag,invpiv);
      XSCAL(col_len, invpiv, cur_col_ptr, one);
#ifdef TIMING0
      scaltime += (MPI_Wtime()-t1);
#endif
      /* restore the pivot entry */

      /* swap pivot and current row in current column */

      if (me == row_owner(gpivot_row)){
        MULTIPLY(pivot.current, invpiv, *(cur_col_ptr+rel_lpivot_row));
      }
      if (me == r_owner)  {
        *(cur_col_ptr) = pivot.entry;
      }
#ifdef TIMING0
      t1 = MPI_Wtime();
#endif
      XCOPY(col_len, cur_col_ptr, one, sav_col_ptr, one);
#ifdef TIMING0
      copycoltime += (MPI_Wtime()-t1);
#endif
#ifdef TIMING0
      t1 = MPI_Wtime();
#endif
      /* send column and pivot down one's row for column j */

      for (rdist = 1;rdist <= MAXDIST;rdist++){
	if (rowplus(rdist) == c_owner) break;
        bytes = sizeof(gpivot_row);
        MPI_Send(&gpivot_row,bytes,MPI_BYTE,rowplus(rdist),
                 LUPIVROWTYPE+j,MPI_COMM_WORLD);

      }
#ifdef TIMING0
      bcastpivstime += (MPI_Wtime()-t1);
#endif
#ifdef TIMING0
      t1 = MPI_Wtime();
#endif
      for (rdist = 1;rdist <= MAXDIST;rdist++){
	if (rowplus(rdist) == c_owner) break;
        bytes=sizeof(DATA_TYPE)*col_len;
        MPI_Send(sav_col_ptr,bytes,MPI_BYTE,
                 rowplus(rdist),LUROWTYPE+j,MPI_COMM_WORLD);
      }
#ifdef TIMING0
      bcastcolstime += (MPI_Wtime()-t1);
#endif
      /* if own active column don't include it in row work anymore */

      row_len--;
      update_ptr += mat_stride; cur_col_ptr += mat_stride;
      cur_row_ptr += mat_stride; 
      sav_pivot_ptr++;
      act_row_ptr += blksz;  sav_piv_row_ptr += blksz;
      cols_used++;
    }else{
      /* recv column and pivot */

     bytes=col_len*sizeof(DATA_TYPE);
       MPI_Irecv(sav_col_ptr,bytes,MPI_BYTE,
                MPI_ANY_SOURCE,LUROWTYPE+j,MPI_COMM_WORLD,&msgrequest);
  

#ifdef TIMING0
      t1 = MPI_Wtime();
#endif
      bytes = 0; dest = -1; type = LUPIVROWTYPE+j;
      bytes=4;
      bytes = sizeof(gpivot_row);
      MPI_Recv(&gpivot_row,bytes,MPI_BYTE,MPI_ANY_SOURCE,type,MPI_COMM_WORLD,&msgstatus);
#ifdef TIMING0
      bcastpivrtime += (MPI_Wtime()-t1);
#endif
 /*      bytes=col_len*sizeof(DATA_TYPE);
          MPI_Recv(sav_col_ptr,bytes,MPI_BYTE,
         MPI_ANY_SOURCE,LUROWTYPE+j,MPI_COMM_WORLD,&msgstatus);  */
      /* if necessary forward column and pivot */

      if ((ringdist % MAXDIST) == 0) {
#ifdef TIMING0
	t1 = MPI_Wtime();
#endif
	for (rdist = 1;rdist <= MAXDIST;rdist++){
	  if (rowplus(rdist) == c_owner) break;
          bytes = sizeof(gpivot_row);
          MPI_Send(&gpivot_row,bytes,MPI_BYTE,rowplus(rdist),
                  LUPIVROWTYPE+j,MPI_COMM_WORLD);
	}
#ifdef TIMING0
	bcastpivstime += (MPI_Wtime()-t1);
#endif
#ifdef TIMING0
	t1 = MPI_Wtime();
#endif
        MPI_Wait(&msgrequest,&msgstatus);    
#ifdef TIMING0
	bcastcolrtime += (MPI_Wtime()-t1);
#endif
#ifdef TIMING0
	t1 = MPI_Wtime();
#endif
	for (rdist = 1;rdist <= MAXDIST;rdist++){
	  if (rowplus(rdist) == c_owner) break;
          bytes=col_len*sizeof(DATA_TYPE);
          MPI_Send(sav_col_ptr,bytes,MPI_BYTE,
                   rowplus(rdist),LUROWTYPE+j,MPI_COMM_WORLD);
	}
#ifdef TIMING0
	bcastcolstime += (MPI_Wtime()-t1);
#endif
      }
    }
    pivot_owner = row_owner(gpivot_row); lpivot_row = grow_to_lrow(gpivot_row);
    act_cur_row_ptr = col1 + rows_used;  act_piv_row_ptr = col1 + lpivot_row;
    piv_row_ptr = cur_col_ptr + (lpivot_row - rows_used);
    piv_col1_row_ptr = act_col_ptr + (lpivot_row - rows_used);
    row_size = (row_len + colcnt)*sizeof(DATA_TYPE);

    /* send current row to owner of pivot row, skip this if pivot row is current row */

    if (gpivot_row != j){
      if (me == r_owner){ /* own current row so pack it up*/

	/* first take current row then portion in stored active columns */
	/* stored portion is at act_cur_row_ptr with col1's stride */

#ifdef TIMING0
	t1 = MPI_Wtime();
#endif
	/* copy from matrix and col1 */

        XCOPY(row_len, cur_row_ptr, mat_stride, sav_row_ptr, one);
        XCOPY(colcnt, cur_col1_row_ptr, col1_stride, sav_row_ptr+row_len, one);
#ifdef TIMING0
 copyrowtime += (MPI_Wtime()-t1);
#endif

      }
    }
    /* update pivot row and save in active row matrix */

    if (me == pivot_owner){
      if (colcnt > 0) {
#ifdef TIMING0
	t1 = MPI_Wtime();
#endif
#ifdef OVERLAP
	/* don't need to update columns which were already updated */

	cols_in_blk_owned = 0;
	for (k = 1; k < (blksz - colcnt); k++)
	  if (j+k < ncols_matrix)
	    if (me == col_owner(j+k)) cols_in_blk_owned++;
	length = row_len - cols_in_blk_owned;
	if (length > 0)
          XGEMMS_(&transA, &transB, &one, &length, &colcnt, &d_min_one,
                 act_piv_row_ptr, &col1_stride,
                 act_row_ptr+(cols_in_blk_owned*blksz), &blksz, &d_one,
                 piv_row_ptr+(cols_in_blk_owned*mat_stride), &mat_stride);
#else
          XGEMMS_(&transA, &transB, &one, &length, &colcnt, &d_min_one,
                 act_piv_row_ptr, &col1_stride,
                 act_row_ptr+(cols_in_blk_owned*blksz), &blksz, &d_one,
                 piv_row_ptr+(cols_in_blk_owned*mat_stride), &mat_stride);
#endif	
#ifdef TIMING0
	rowupdtime += (MPI_Wtime()-t1);
#endif
      }
      /* copy pivot row to temp holder */

#ifdef TIMING0
      t1 = MPI_Wtime();
#endif

	/* copy from matrix and col1 */

      XCOPY(row_len, piv_row_ptr, mat_stride, temp_row_ptr, one);
      XCOPY(colcnt, piv_col1_row_ptr, col1_stride, temp_row_ptr+row_len, one);
      
#ifdef TIMING0
      copypivrowtime += (MPI_Wtime()-t1);
#endif
    }
    /* broadcast pivot row */

#ifdef TIMING0
    t1 = MPI_Wtime();
#endif

    bytes=sizeof(DATA_TYPE)*row_size ;
    MPI_Bcast((char *) temp_row_ptr, row_size, MPI_CHAR, mesh_row(pivot_owner), col_comm);   
    MPI_Barrier(col_comm);  
#ifdef TIMING0
    bcastrowtime += (MPI_Wtime()-t1);
#endif

#ifdef TIMING0
    t1 = MPI_Wtime();
#endif

    XCOPY(row_len, temp_row_ptr, one, sav_piv_row_ptr, blksz);

#ifdef TIMING0
    copypivrowtime += (MPI_Wtime()-t1);
#endif
    if (gpivot_row != j){
        if (me != pivot_owner && me == r_owner) 
	  {
#ifdef TIMING0
	    t1 = MPI_Wtime();
#endif
            bytes=(row_len+colcnt)*sizeof(DATA_TYPE);
	    MPI_Send(sav_row_ptr,bytes,MPI_BYTE,pivot_owner,
                     LUSENDTYPE+j,MPI_COMM_WORLD);
#ifdef TIMING0
	    sendrowtime += (MPI_Wtime()-t1);
#endif
	  }
      if (me == pivot_owner){
	/* receive top row and copy into pivot row */

#ifdef TIMING0
	t1 = MPI_Wtime();
#endif
	if (me != r_owner) {
          bytes=(row_len+colcnt)*sizeof(DATA_TYPE);
          MPI_Recv(sav_row_ptr,bytes,MPI_BYTE,r_owner,
                   LUSENDTYPE+j,MPI_COMM_WORLD,&msgstatus);
        }
#ifdef TIMING0
	recvrowtime += (MPI_Wtime()-t1);
#endif
#ifdef TIMING0
	t1 = MPI_Wtime();
#endif
	/* copy from matrix and col1 */

	XCOPY(row_len, sav_row_ptr, one, piv_row_ptr, mat_stride);
	XCOPY(colcnt, sav_row_ptr+row_len, one, piv_col1_row_ptr, col1_stride);
#ifdef TIMING0
	copyrow1time += (MPI_Wtime()-t1);
#endif
      }
      if (me == r_owner) { /* copy pivot row into current row */
#ifdef TIMING0
	t1 = MPI_Wtime();
#endif

	/* copy from matrix and col1 */

        XCOPY(row_len, temp_row_ptr, one, cur_row_ptr, mat_stride);
        XCOPY(colcnt, temp_row_ptr+row_len, one, cur_col1_row_ptr, col1_stride);

#ifdef TIMING0
	copypivrow1time += (MPI_Wtime()-t1);
#endif
      }
    }
    if ((me != c_owner) && ((ringdist % MAXDIST) != 0)) {
#ifdef TIMING0
      t1 = MPI_Wtime();
#endif
      MPI_Wait(&msgrequest,&msgstatus);
#ifdef TIMING0
      bcastcolrtime += (MPI_Wtime()-t1);
#endif
    }
    /* saved this active row and column so get ready for next ones */

    if (me == r_owner) { /* finished with this row so update all 
	                    column pointers */

      col_len--; rows_used++; update_ptr++; cur_row_ptr++; cur_col1_row_ptr++;
      cur_col_ptr++; sav_col_ptr++;  act_col_ptr++;
    }
    colcnt++;

#ifdef OVERLAP
    cols_in_blk_owned = 0;
    for (k = 1; k <= (blksz - colcnt); k++)
      if (j+k < ncols_matrix)
	if (me == col_owner(j+k)) cols_in_blk_owned++;
    if (cols_in_blk_owned > 0){ /* update current column with latest column */
#ifdef TIMING0
      t1 = MPI_Wtime();
#endif
      XGEMMS_(&transA, &transB,
             &col_len, &cols_in_blk_owned, &one, &d_min_one,
             sav_col_ptr, &col1_stride,
             sav_piv_row_ptr, &blksz, &d_one,
             cur_col_ptr, &mat_stride);
      
#ifdef TIMING0
      colupdtime += (MPI_Wtime()-t1);
#endif
    }
#endif

    sav_col_ptr += col1_stride;
    sav_piv_row_ptr++;

    /* if we have saved up enough columns, we do the outer product update. */

    if (colcnt == blksz)
      if (j != ncols_matrix-1){
#ifdef TIMING0
	t1 = MPI_Wtime();
#endif
        XGEMM_(&transA, &transB, &col_len, &row_len, &colcnt, &d_min_one,
               act_col_ptr, &col1_stride, act_row_ptr, &blksz, &d_one,
               update_ptr, &mat_stride);
#ifdef TIMING0
	updatetime += (MPI_Wtime()-t1);
#endif
	/* reset active matrix pointers */

	colcnt = 0;
	act_col_ptr = sav_col_ptr = col1 + rows_used;	
	act_row_ptr = sav_piv_row_ptr = row1;
      }
#ifdef PRINT_STATUS
    if (((j%1000) == 0) && (me == 0)) {
      fprintf(stderr," Column %d completed\n",j);
    }
#endif
  }
#ifdef TIMING0
  totaltime = MPI_Wtime() - t2;
#endif
#ifdef TIMING0
  msgtime = xpivmsgtime+bcastpivstime+bcastpivrtime+bcastcolstime+bcastcolrtime+bcastrowtime+sendrowtime+recvrowtime;
  copytime = copycoltime+copyrowtime+copyrow1time+copypivrowtime+copypivrow1time;
  dgemmtime = updatetime+colupdtime+rowupdtime+scaltime;
  showtime("Time to xchgpivot",&xpivmsgtime);
  showtime("Time to do send in bcast pivot",&bcastpivstime);
  showtime("Time to do recv in bcast pivot",&bcastpivrtime);
  tmp = bcastpivrtime+bcastpivstime; 
  showtime("Time to do bcast pivot",&tmp);
  showtime("Time to do send in bcast cur col",&bcastcolstime);
  showtime("Time to do recv bcast cur col",&bcastcolrtime);
  tmp = bcastcolrtime+bcastcolstime;
  showtime("Time to do bcast cur col",&tmp);
  tmp = bcastcolrtime+bcastcolstime+bcastpivrtime+bcastpivstime;
  showtime("Time to do bcast cur col and pivot",&tmp);
  showtime("Time to bcast piv row",&bcastrowtime);
  showtime("Time to send cur row",&sendrowtime);
  showtime("Time to recv cur row",&recvrowtime);
  showtime("Total msg passing time",&msgtime);
  tmp = 100*msgtime/totaltime;
  showtime("Percent msg passing time",&tmp);
  showtime("Time to copy cur col",&copycoltime);
  showtime("Time to copy cur row to sav row",&copyrowtime);
  showtime("Time to copy piv row to sav piv",&copypivrowtime);
  showtime("Time to copy sav row to cur row",&copyrow1time);
  showtime("Time to copy sav piv  to piv row",&copypivrow1time);
  showtime("Total copying time",&copytime);
  tmp = 100*copytime/totaltime;
  showtime("Percent copying time",&tmp);
  showtime("Time to scale cur col",&scaltime);
  showtime("Time to update cur col",&colupdtime);
  showtime("Time to update piv row",&rowupdtime);
  showtime("Time to update matrix",&updatetime);
  showtime("Total update time",&dgemmtime);
  tmp = 100*dgemmtime/totaltime;
  showtime("Percent update time",&tmp);
  showtime("Total time in factor",&totaltime);
#endif   
}
