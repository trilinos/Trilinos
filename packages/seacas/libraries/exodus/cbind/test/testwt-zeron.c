/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
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
* testwt - test write an ExodusII database file
*
* author - Sandia National Laboratories
*          Larry A. Schoof - Original
*          Vic Yarberry    - Added headers and error logging
*               7/7/93          Modified for use with Exodus 2.00
*
*          
* environment - UNIX
*
* entry conditions - 
*
* exit conditions - 
*
* revision history - 
*
*  This is a test program for the C binding of the EXODUS II 
*  database write routines.
*
*
*****************************************************************************/


#include <stdlib.h>
#include <stdio.h>
/* #include "netcdf.h" */
#include "exodusII.h"

int main (int argc, char **argv)
{
   int exoid, num_dim, num_nodes, num_elem, num_elem_blk;
   int num_elem_in_block[10], num_nodes_per_elem[10];
   int num_node_sets, num_side_sets, error;
   int i, j, *elem_map;
   int ebids[10];
   int  num_qa_rec, num_info;
   int num_glo_vars;
   int whole_time_step, num_time_steps;
   int CPU_word_size,IO_word_size;

   float *glob_var_vals, *nodal_var_vals, *elem_var_vals;
   float time_value;
   float x[100], y[100], z[100];
   char *coord_names[3], *qa_record[2][4], *info[3], *var_names[3];

   ex_opts (EX_VERBOSE | EX_ABORT );

/* Specify compute and i/o word size */

   CPU_word_size = 0;                   /* sizeof(float) */
   IO_word_size = 4;                    /* (4 bytes) */

/* create EXODUS II file */

   exoid = ex_create ("test.exo",       /* filename path */
                       EX_CLOBBER,      /* create mode */
                       &CPU_word_size,  /* CPU float word size in bytes */
                       &IO_word_size);  /* I/O float word size in bytes */
   printf ("after ex_create for test.exo, exoid = %d\n", exoid);
   printf (" cpu word size: %d io word size: %d\n",CPU_word_size,IO_word_size);

   /* ncopts = NC_VERBOSE; */

/* initialize file with parameters */

   num_dim = 3;
   num_nodes = 0;
   num_elem = 0;
   num_elem_blk = 7;
   num_node_sets = 0;
   num_side_sets = 0;

   error = ex_put_init (exoid, "This is a test", num_dim, num_nodes, num_elem,
                        num_elem_blk, num_node_sets, num_side_sets);

   printf ("after ex_put_init, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

/* write nodal varcoordinates values and names to database */

/* Quad #1 */
   x[0] = 0.0; y[0] = 0.0; z[0] = 0.0;
   x[1] = 1.0; y[1] = 0.0; z[1] = 0.0;
   x[2] = 1.0; y[2] = 1.0; z[2] = 0.0;
   x[3] = 0.0; y[3] = 1.0; z[3] = 0.0;

/* Quad #2 */
   x[4]  =  1.0; y[4]  =  0.0; z[4]  =  0.0;
   x[5]  =  2.0; y[5]  =  0.0; z[5]  =  0.0;
   x[6]  =  2.0; y[6]  =  1.0; z[6]  =  0.0;
   x[7]  =  1.0; y[7]  =  1.0; z[7]  =  0.0;

/* Hex #1 */
   x[8]  =  0.0; y[8]  =  0.0; z[8]  =  0.0;
   x[9]  = 10.0; y[9]  =  0.0; z[9]  =  0.0;
   x[10] = 10.0; y[10] =  0.0; z[10] =-10.0;
   x[11] =  1.0; y[11] =  0.0; z[11] =-10.0;
   x[12] =  1.0; y[12] = 10.0; z[12] =  0.0;
   x[13] = 10.0; y[13] = 10.0; z[13] =  0.0;
   x[14] = 10.0; y[14] = 10.0; z[14] =-10.0;
   x[15] =  1.0; y[15] = 10.0; z[15] =-10.0;

/* Tetra #1 */
   x[16] =  0.0; y[16] =  0.0; z[16] =  0.0;
   x[17] =  1.0; y[17] =  0.0; z[17] =  5.0;
   x[18] = 10.0; y[18] =  0.0; z[18] =  2.0;
   x[19] =  7.0; y[19] =  5.0; z[19] =  3.0;

/* Wedge #1 */
   x[20] =  3.0; y[20] =  0.0; z[20] =  6.0;
   x[21] =  6.0; y[21] =  0.0; z[21] =  0.0;
   x[22] =  0.0; y[22] =  0.0; z[22] =  0.0;
   x[23] =  3.0; y[23] =  2.0; z[23] =  6.0;
   x[24] =  6.0; y[24] =  2.0; z[24] =  2.0;
   x[25] =  0.0; y[25] =  2.0; z[25] =  0.0;

/* Tetra #2 */
   x[26] =  2.7; y[26] =  1.7; z[26] =  2.7;
   x[27] =  6.0; y[27] =  1.7; z[27] =  3.3;
   x[28] =  5.7; y[28] =  1.7; z[28] =  1.7;
   x[29] =  3.7; y[29] =  0.0; z[29] =  2.3;

/* 3d Tri */
   x[30] =  0.0; y[30] =  0.0; z[30] =  0.0;
   x[31] = 10.0; y[31] =  0.0; z[31] =  0.0;
   x[32] = 10.0; y[32] = 10.0; z[32] = 10.0;

   /*   error = ex_put_coord (exoid, x, y, z); */
   printf ("after ex_put_coord, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   coord_names[0] = "xcoor";
   coord_names[1] = "ycoor";
   coord_names[2] = "zcoor";

   error = ex_put_coord_names (exoid, coord_names);
   printf ("after ex_put_coord_names, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

/* write element order map */

   elem_map = (int *) calloc(num_elem, sizeof(int));

   for (i=1; i<=num_elem; i++)
   {
      elem_map[i-1] = i;
   }

   error = ex_put_map (exoid, elem_map);
   printf ("after ex_put_map, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   free (elem_map);


/* write element block parameters */

   num_elem_in_block[0] = 0;
   num_elem_in_block[1] = 0;
   num_elem_in_block[2] = 0;
   num_elem_in_block[3] = 0;
   num_elem_in_block[4] = 0;
   num_elem_in_block[5] = 0;
   num_elem_in_block[6] = 0;

   num_nodes_per_elem[0] = 4; /* elements in block #1 are 4-node quads  */
   num_nodes_per_elem[1] = 4; /* elements in block #2 are 4-node quads  */
   num_nodes_per_elem[2] = 8; /* elements in block #3 are 8-node hexes  */
   num_nodes_per_elem[3] = 4; /* elements in block #4 are 4-node tetras */
   num_nodes_per_elem[4] = 6; /* elements in block #5 are 6-node wedges */
   num_nodes_per_elem[5] = 8; /* elements in block #6 are 8-node tetras */
   num_nodes_per_elem[6] = 3; /* elements in block #7 are 3-node tris   */

   ebids[0] = 10;
   ebids[1] = 11;
   ebids[2] = 12;
   ebids[3] = 13;
   ebids[4] = 14;
   ebids[5] = 15;
   ebids[6] = 16;

   error = ex_put_elem_block (exoid, ebids[0], "quad", num_elem_in_block[0],
                              num_nodes_per_elem[0], 1);
   printf ("after ex_put_elem_block, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   error = ex_put_elem_block (exoid, ebids[1], "quad", num_elem_in_block[1],
                               num_nodes_per_elem[1], 1);
   printf ("after ex_put_elem_block, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   error = ex_put_elem_block (exoid, ebids[2], "hex", num_elem_in_block[2],
                               num_nodes_per_elem[2], 1);
   printf ("after ex_put_elem_block, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   error = ex_put_elem_block (exoid, ebids[3], "tetra", num_elem_in_block[3],
                               num_nodes_per_elem[3], 1);
   printf ("after ex_put_elem_block, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   error = ex_put_elem_block (exoid, ebids[4], "wedge", num_elem_in_block[4],
                               num_nodes_per_elem[4], 1);
   printf ("after ex_put_elem_block, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   error = ex_put_elem_block (exoid, ebids[5], "tetra", num_elem_in_block[5],
                               num_nodes_per_elem[5], 1);
   printf ("after ex_put_elem_block, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   error = ex_put_elem_block (exoid, ebids[6], "tri", num_elem_in_block[6],
                               num_nodes_per_elem[6], 1);
   printf ("after ex_put_elem_block, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }

/* write QA records; test empty and just blank-filled records */

   num_qa_rec = 2;


   qa_record[0][0] = "TESTWT";
   qa_record[0][1] = "testwt";
   qa_record[0][2] = "07/07/93";
   qa_record[0][3] = "15:41:33";
   qa_record[1][0] = "";
   qa_record[1][1] = "                            ";
   qa_record[1][2] = "";
   qa_record[1][3] = "                        ";

   error = ex_put_qa (exoid, num_qa_rec, qa_record);
   printf ("after ex_put_qa, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }


/* write information records; test empty and just blank-filled records */

   num_info = 3;


   info[0] = "This is the first information record.";
   info[1] = "";
   info[2] = "                                     ";

   error = ex_put_info (exoid, num_info, info);
   printf ("after ex_put_info, error = %d\n", error);

   if (error) {
     ex_close (exoid);
     exit(-1);
   }



/* write results variables parameters and names */

   num_glo_vars = 1;

   var_names[0] = "glo_vars";

   error = ex_put_var_param (exoid, "g", num_glo_vars);
   printf ("after ex_put_var_param, error = %d\n", error);
   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   error = ex_put_var_names (exoid, "g", num_glo_vars, var_names);
   printf ("after ex_put_var_names, error = %d\n", error);
   if (error) {
     ex_close (exoid);
     exit(-1);
   }



/* for each time step, write the analysis results;
 * the code below fills the arrays glob_var_vals, 
 * nodal_var_vals, and elem_var_vals with values for debugging purposes;
 * obviously the analysis code will populate these arrays
 */

   whole_time_step = 1;
   num_time_steps = 10;

   glob_var_vals = (float *) calloc (num_glo_vars, CPU_word_size);
   nodal_var_vals = (float *) calloc (num_nodes, CPU_word_size);
   elem_var_vals = (float *) calloc (4, CPU_word_size);

   for (i=0; i<num_time_steps; i++)
   {
     time_value = (float)(i+1)/100.;

/* write time value */

     error = ex_put_time (exoid, whole_time_step, &time_value);
     printf ("after ex_put_time, error = %d\n", error);

     if (error) {
       ex_close (exoid);
       exit(-1);
     }

/* write global variables */

     for (j=0; j<num_glo_vars; j++)
     {
       glob_var_vals[j] = (float)(j+2) * time_value;
     }

     error = ex_put_glob_vars (exoid, whole_time_step, num_glo_vars, 
                               glob_var_vals);
     printf ("after ex_put_glob_vars, error = %d\n", error);

     if (error) {
       ex_close (exoid);
       exit(-1);
     }

     whole_time_step++;

/* update the data file; this should be done at the end of every time step
 * to ensure that no data is lost if the analysis dies
 */
     error = ex_update (exoid);
     printf ("after ex_update, error = %d\n", error);
     if (error) {
       ex_close (exoid);
       exit(-1);
     }
   }
   free(glob_var_vals);
   free(nodal_var_vals);
   free(elem_var_vals);


/* close the EXODUS files
 */
   error = ex_close (exoid);
   printf ("after ex_close, error = %d\n", error);
   if (error) {
     ex_close (exoid);
     exit(-1);
   }
   return 0;
}
