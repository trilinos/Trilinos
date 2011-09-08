/*
//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#include <stdlib.h>
#include <stdio.h>
#include "az_aztec.h"

void create_vbr(char *partition_file, int *proc_config,
		int *N_global, int *N_blk_global, 
		int *n_nonzeros, int *n_blk_nonzeros,
		int *N_update, int **update,
		int *bindx_msr, double *val_msr,
		double **val, int **indx, int **rpntr, int **cpntr,
		int **bpntr, int **bindx)
#undef DEBUG 
     /*  read ASCII data file:
	 line 1: N_global, number of entries (%d,%d)
	 line 2-...: i,j,real (%d, %d, %f)
     */

{
  FILE *data ;


  int i, n_entries, N_columns;
  int ii, jj ;
  int kk = 0;
  int max_ii = 0, max_jj = 0;
  int ione = 1;
  double value;
  double *cnt;
  int *pntr, *indx1, *pntr1;
  double *val1;
  int blocksize, N_blk_equations, N_block_entries, n_vbr_nonzeros, n_msr_nonzeros;
  int total_msr_storage, total_vbr_storage;
  int variable_block_size, blk_type;
  int cur_blk_ptr=0, prev_blk_ptr;

  if(proc_config[AZ_node] == 0) 
    { 

      /* Do case where command line argument is an integer.
	 Interpret integer as the constant block size */
      printf("***************************************************************\n");
      if (partition_file[0] >='0' && partition_file[0] <='9')
	{
	  blocksize = atoi(partition_file);
	  printf("Using block size of %d to convert from MSR to VBR\n",blocksize);
	  N_blk_equations = *N_global/blocksize;

	  /* Allocate memory for cpntr */
	  *cpntr = (int   *) calloc(N_blk_equations+2,sizeof(int)) ;

	  /* Define block sizes for all but last equation */
	  for (i=0; i<N_blk_equations; i++) (*cpntr)[i] = blocksize;

	  /* Check if number of equations is multiple of blocksize */
	  variable_block_size = *N_global%blocksize;
	  blk_type = blocksize;

	  if (variable_block_size)
	    {
	      N_blk_equations ++;
	      (*cpntr)[N_blk_equations-1] = variable_block_size;
	      blk_type = -blocksize;
	    }
	}
      else
	{
	  /* Otherwise command line arg is a file name containing partition 
	     information.  
	     The first line of the file must be the integer value zero.
	     The last line of the file must equal the number of global equations,
	     i.e., N_global.
	     Lines in between are incremented by the number of equations per
             block row.
	  */
	  /* This should be a short file, so read once to get number of block
	     equations, then read again to fill values */
	  printf("Using partition from %s to convert from MSR to VBR\n",
		 partition_file);
	  data = fopen(partition_file,"r") ;
	  N_blk_equations = 0;
	  while(cur_blk_ptr !=*N_global)
	    {
	      fscanf(data, "%ld", cur_blk_ptr);
	      N_blk_equations++;
	    }
	  close(data);

	  /* Allocate memory for cpntr */
	  *cpntr = (int   *) calloc(N_blk_equations+1,sizeof(int)) ;

	  N_blk_equations = 0;
	  data = fopen(partition_file,"r") ;
	  fscanf(data, "%ld", prev_blk_ptr);
	  while(cur_blk_ptr !=*N_global)
	    {
	      fscanf(data, "%ld", cur_blk_ptr);
	      (*cpntr)[N_blk_equations] = cur_blk_ptr - prev_blk_ptr;
	      prev_blk_ptr = cur_blk_ptr;
	      N_blk_equations++;
	    }
	  close(data);
	  blk_type = -1; /* assume variable block for now */
    
	}

      /* Estimate storage needed for VBR and allocate space */

      N_block_entries = *n_nonzeros;
      n_vbr_nonzeros = min(abs(*n_nonzeros * blocksize * blocksize),
			   420000000/8);
      *N_blk_global = N_blk_equations;

      printf("\nEstimated Storage parameters for VBR:\n");
      printf("   Number of block  equations = %d\n",N_blk_equations);
      printf("   Number of block  entries   = %d\n",N_block_entries);
      printf("   Number of scalar entries   = %d\n",n_vbr_nonzeros);
	
      
      *bpntr = (int   *) calloc(N_blk_equations+1,sizeof(int)) ;
      *rpntr = (int   *) calloc(N_blk_equations+1,sizeof(int)) ;
      *bindx = (int   *) calloc(N_block_entries+1,sizeof(int)) ;
      *indx  = (int   *) calloc(N_block_entries+1,sizeof(int)) ;
      *val = (double *) calloc(n_vbr_nonzeros+1,  sizeof(double)) ;
  
      
      while (n_vbr_nonzeros >= *n_nonzeros && (*val) == NULL)
	{
	  printf("Error: Unable to allocate %ld bytes to create VBR matrix.\n",
		 n_vbr_nonzeros*sizeof(double));
	  printf("       Trying to allocate %ld bytes.\n",
		 n_vbr_nonzeros*sizeof(double)/2);
	  n_vbr_nonzeros /= 2;
	  *val = (double *) calloc(n_vbr_nonzeros+1,  sizeof(double)) ;
	}

      AZ_msr2vbr(*val, *indx, *rpntr, *cpntr, *bpntr, *bindx, bindx_msr,val_msr, 
		 N_blk_equations, N_blk_equations, N_block_entries,
		 n_vbr_nonzeros, blk_type);

      n_msr_nonzeros = *n_nonzeros;

      *n_nonzeros = (*indx)[(*bpntr)[*N_blk_global]];
      *n_blk_nonzeros = (*bpntr)[*N_blk_global];
      
      *bindx = (int   *) realloc((void *) (*bindx),
				 (*n_blk_nonzeros+1)*sizeof(int)) ;
      *indx  = (int   *) realloc((void *) (*indx),
				 (*n_blk_nonzeros+1)*sizeof(int)) ;
      *val = (double *) realloc((void *) (*val),
				(*n_nonzeros+1)*sizeof(double)) ;
      printf("\nActual Storage parameters for VBR:\n");
      printf("   Number of block  equations = %d\n",N_blk_equations);
      printf("   Number of block  entries   = %d\n",*n_blk_nonzeros);
      printf("   Number of scalar entries   = %d\n",*n_nonzeros);
      
      total_msr_storage = 4*  (n_msr_nonzeros+1)  +   8*(n_msr_nonzeros+1);
      total_vbr_storage = 4*3*(N_blk_equations+1) + 4*2*(*n_blk_nonzeros+1) + 
	                  8*(*n_nonzeros);
      printf("\nTotal MSR storage (bytes)   = %d\n",total_msr_storage);
      printf(  "Total VBR storage (bytes)   = %d\n",total_vbr_storage);
      printf(  "Ratio of VBR to MSR storage = %5.2f\n",
	     (float)total_vbr_storage/(float)total_msr_storage);

  
      printf("***************************************************************\n");
    }
  /* end create_vbr */
}
