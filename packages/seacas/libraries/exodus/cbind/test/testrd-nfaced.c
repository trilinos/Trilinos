/*
 * Copyright (c) 2010 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.  
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/*****************************************************************************
*
* testrd - read exodus file test-nsided.exo created by testwt-nsided
*
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "netcdf.h"
#include <assert.h>
#include "exodusII.h"

int main (int argc, char **argv)
{
  int exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int *num_elem_in_block, *num_face_in_block, *num_nodes_per_elem, *num_edges_per_elem, *num_faces_per_elem, *num_attr;
  int error, nnodes;
  int i, j, k;
  int *connect, *fconnect;
  int *ids, *nnpe, *nnpf; 
  int num_qa_rec, num_info;
  int CPU_word_size,IO_word_size;
  int idum;

  float *x, *y, *z;
  float version, fdum;

  char *coord_names[3], *qa_record[2][4], *info[3];
  char *block_names[10];
  char *elem_type[10];
  char name[MAX_STR_LENGTH+1];
  char *cdum = 0;

  CPU_word_size = 0;                    /* sizeof(float) */
  IO_word_size = 0;                     /* use what is stored in file */

  ex_opts (EX_VERBOSE | EX_ABORT );

  /* open EXODUS II files */
  exoid = ex_open ("test-nfaced.exo",  /* filename path */
                   EX_READ,             /* access mode = READ */
                   &CPU_word_size,      /* CPU word size */
                   &IO_word_size,       /* IO word size */
                   &version);           /* ExodusII library version */

  printf ("\nafter ex_open\n");
  if (exoid < 0) exit(1);

  printf ("test.exo is an EXODUSII file; version %4.2f\n",
          version);
  printf ("         I/O word size %1d\n",IO_word_size);

  ex_inquire(exoid,EX_INQ_LIB_VERS, &idum, &version, cdum);
  printf ("EXODUSII Library API; version %4.2f (%d)\n", version, idum);

  /* read database parameters */
  {
    ex_init_params par;
    error = ex_get_init_ext (exoid, &par);

    printf ("after ex_get_init, error = %3d\n", error);

    printf ("database parameters:\n");
    printf ("title =  '%s'\n",par.title);
    printf ("num_dim = %3d\n",par.num_dim);
    printf ("num_nodes = %3d\n",par.num_nodes);
    printf ("num_edge = %3d\n",par.num_edge);
    printf ("num_face = %3d\n",par.num_face);
    printf ("num_elem = %3d\n",par.num_elem);
    printf ("num_elem_blk = %3d\n",par.num_elem_blk);
    printf ("num_node_sets = %3d\n",par.num_node_sets);
    printf ("num_side_sets = %3d\n",par.num_side_sets);

    num_dim = par.num_dim;
    num_nodes = par.num_nodes;
    num_elem = par.num_elem;
    num_elem_blk = par.num_elem_blk;
  }
    
  assert(num_dim == 3);

  /* read nodal coordinates values and names from database */

  x = (float *) calloc(num_nodes, sizeof(float));
  y = (float *) calloc(num_nodes, sizeof(float));
  z = (float *) calloc(num_nodes, sizeof(float));

  error = ex_get_coord (exoid, x, y, z);
  printf ("\nafter ex_get_coord, error = %3d\n", error);

  printf ("x, y, z coords = \n");
  for (i=0; i<num_nodes; i++)
    {
      printf ("%5.1f\t%5.1f\t%5.1f\n", x[i], y[i], z[i]);
    }

  free (x);
  free (y);
  free (z);

  for (i=0; i<num_dim; i++) {
    coord_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
  }
  
  error = ex_get_coord_names (exoid, coord_names);
  printf ("\nafter ex_get_coord_names, error = %3d\n", error);

  printf ("x coord name = '%s'\n", coord_names[0]);
  printf ("y coord name = '%s'\n", coord_names[1]);
  printf ("z coord name = '%s'\n", coord_names[2]);

  for (i=0; i<num_dim; i++)
    free(coord_names[i]);

  /* read element block parameters */
  if (num_elem_blk > 0) {
    ids = (int *) calloc(num_elem_blk, sizeof(int));
    num_elem_in_block = (int *) calloc(num_elem_blk, sizeof(int));
    num_face_in_block = (int *) calloc(num_elem_blk, sizeof(int));
    num_nodes_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
    num_edges_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
    num_faces_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
    num_attr = (int *) calloc(num_elem_blk, sizeof(int));
     
    for (i=0; i<num_elem_blk; i++) {
      elem_type[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
      block_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
    }

    error = ex_get_elem_blk_ids (exoid, ids);
    printf ("\nafter ex_get_elem_blk_ids, error = %3d\n", error);
     
    error = ex_get_names(exoid, EX_ELEM_BLOCK, block_names);
    printf ("\nafter ex_get_names, error = %3d\n", error);
    
    for (i=0; i<num_elem_blk; i++) {
      ex_get_name(exoid, EX_ELEM_BLOCK, ids[i], name);
      if (strcmp(name, block_names[i]) != 0) {
	printf ("error in ex_get_name for block id %d\n", ids[i]);
      }
      error = ex_get_block (exoid, EX_ELEM_BLOCK, ids[i], elem_type[i],
			    &(num_elem_in_block[i]),
			    &(num_nodes_per_elem[i]),
			    &(num_edges_per_elem[i]),
			    &(num_faces_per_elem[i]),
			    &(num_attr[i]));
      printf ("\nafter ex_get_elem_block, error = %d\n", error);
         
      printf ("element block id = %2d\n",ids[i]);
      printf ("element block type = '%s'\n", elem_type[i]);
      printf ("num_elem_in_block = %2d\n",num_elem_in_block[i]);
      printf ("num_total_nodes_per_block = %2d\n",num_nodes_per_elem[i]);
      printf ("num_total_edges_per_block = %2d\n",num_edges_per_elem[i]);
      printf ("num_total_faces_per_block = %2d\n",num_faces_per_elem[i]);
      printf ("num_attr = %2d\n",num_attr[i]);
      printf ("name = '%s'\n",block_names[i]);
    }
     
  }
   
  /* read connectivity */
  for (i=0; i<num_elem_blk; i++) {
    if (num_elem_in_block[i] > 0) {
      if (strcmp(elem_type[i], "NFACED") == 0 || strcmp(elem_type[i], "nfaced") == 0) {
	int nfaces = 0;
	connect = (int *) calloc((num_faces_per_elem[i]), sizeof(int));

	nnpe = (int *) calloc(num_elem_in_block[i], sizeof(int));
	error = ex_get_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK, ids[i], nnpe);
	printf ("\nafter ex_get_entity_count_per_polyhedra, error = %d\n", error);
	
	for (j=0; j < num_elem_in_block[i]; j++) {
	  nfaces += nnpe[j];
	}
	assert(nfaces == num_faces_per_elem[i]);
	
	error = ex_get_conn (exoid, EX_ELEM_BLOCK, ids[i], NULL, NULL, connect);
	printf ("\nafter ex_get_conn, error = %d\n", error);
	
	printf ("face connectivity array for elem block %2d\n", ids[i]);
	nfaces = 0;
	for (j=0; j < num_elem_in_block[i]; j++) {
	  printf("Element %d, %d faces:\t", j+1, nnpe[j]);
	  for (k=0; k < nnpe[j]; k++) {
	    printf("%3d ", connect[nfaces+k]);
	  }
	  printf("\n");
	  nfaces += nnpe[j];
	}

	/* Now get the faces and their connectivity... */
	/*
	 * Convention is that the faces for an nfaced block are in a
	 * face block which has the same id as the element block...
	 * (Or, at least let's try that for awhile and see if it works...)
	 */

	/* NOTE: We are overwriting the element block data here... */
	error = ex_get_block (exoid, EX_FACE_BLOCK, ids[i], elem_type[i],
			      &(num_face_in_block[i]),
			      &(num_nodes_per_elem[i]),
			      NULL, NULL,
			      &(num_attr[i]));

	printf ("\nafter ex_get_block (EX_FACE_BLOCK), error = %d\n", error);
	
	error = ex_get_names(exoid, EX_FACE_BLOCK, block_names);
	printf ("\nafter ex_get_names, error = %3d\n", error);

	printf ("\tface block id = %2d\n",ids[i]);
	printf ("\tface block type = '%s'\n", elem_type[i]);
	printf ("\tnum_face_in_block = %2d\n",num_face_in_block[i]);
	printf ("\tnum_total_nodes_per_block = %2d\n",num_nodes_per_elem[i]);
	printf ("\tnum_attr = %2d\n",num_attr[i]);
	printf ("\tname = '%s'\n",block_names[i]);

	fconnect = (int *) calloc((num_nodes_per_elem[i]), sizeof(int));
	nnpf = (int *) calloc(num_face_in_block[i], sizeof(int));
	error = ex_get_entity_count_per_polyhedra(exoid, EX_FACE_BLOCK, ids[i], nnpf);
	printf ("\nafter ex_get_entity_count_per_polyhedra, error = %d\n", error);
	
	nnodes = 0;
	for (j=0; j < num_face_in_block[i]; j++) {
	  nnodes += nnpf[j];
	}
	assert(nnodes == num_nodes_per_elem[i]);
	
	error = ex_get_conn (exoid, EX_FACE_BLOCK, ids[i], fconnect, NULL, NULL);
	printf ("\nafter ex_get_conn, error = %d\n", error);
	
	printf ("node connectivity array for face block %2d\n", ids[i]);
	nnodes = 0;
	for (j=0; j < num_face_in_block[i]; j++) {
	  printf("Face %d, %d nodes:\t", j+1, nnpf[j]);
	  for (k=0; k < nnpf[j]; k++) {
	    printf("%3d ", fconnect[nnodes+k]);
	  }
	  printf("\n");
	  nnodes += nnpf[j];
	}
	free(fconnect);
	free(nnpe);
	free(nnpf);
      } else {
	connect = (int *) calloc((num_nodes_per_elem[i] * num_elem_in_block[i]), 
				 sizeof(int));
	error = ex_get_elem_conn (exoid, ids[i], connect);
	printf ("\nafter ex_get_elem_conn, error = %d\n", error);
	
	printf ("connect array for elem block %2d\n", ids[i]);
	
	for (j=0; j<num_nodes_per_elem[i]; j++) {
	  printf ("%3d\n", connect[j]);
	}
      }
      free (connect);
    }
  }

  for (i=0; i<num_elem_blk; i++) {
    free(elem_type[i]);
    free(block_names[i]);
  }
  if (num_elem_blk > 0) {
    free (ids);
    free (num_nodes_per_elem);
    free (num_edges_per_elem);
    free (num_faces_per_elem);
    free (num_attr);
  }

  /* read QA records */
  ex_inquire (exoid, EX_INQ_QA, &num_qa_rec, &fdum, cdum);

  for (i=0; i<num_qa_rec; i++)
    {
      for (j=0; j<4; j++)
	{
	  qa_record[i][j] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
	}
    }

  error = ex_get_qa (exoid, qa_record); 
  printf ("\nafter ex_get_qa, error = %3d\n", error);

  printf ("QA records = \n");
  for (i=0; i<num_qa_rec; i++) 
    {
      for (j=0; j<4; j++)
	{
	  printf (" '%s'\n", qa_record[i][j]);
	  free(qa_record[i][j]);
	}
    }

  /* read information records */

  error = ex_inquire (exoid, EX_INQ_INFO, &num_info, &fdum, cdum);
  printf ("\nafter ex_inquire, error = %3d\n", error);

  for (i=0; i<num_info; i++)
    {
      info[i] = (char *) calloc ((MAX_LINE_LENGTH+1), sizeof(char));
    }

  error = ex_get_info (exoid, info); 
  printf ("\nafter ex_get_info, error = %3d\n", error);

  printf ("info records = \n");
  for (i=0; i<num_info; i++)
    {
      printf (" '%s'\n", info[i]);
      free(info[i]);
    }

  error = ex_close (exoid);
  printf ("\nafter ex_close, error = %3d\n", error);
  return 0;
}
