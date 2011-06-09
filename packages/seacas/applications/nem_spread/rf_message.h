/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
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
/*
*
*
*  MESSAGE TYPE FILE -- ALL COMMUNICATION TYPES SHOULD BE SET HERE
*
*                    -- DEFAULTS FOR SET-UP OF COMMUNICATIONS BUFFERS
*
*     In the future, the range of message types that each message operation
* takes should also be included here.  This file is meant to be a
* "clearinghouse" for message types within the program.
*     Also suggest that all message types have a common style.  They should
* being with "MT_", and should only be 19 characters long.
*
*     At the present, there is a possible problem with the overlap of
* message types with the krylib library.  This should probably be resolved
* by restricting the krylib library to have a limited and specified range of
* message types in the high part of the allowed range (e.g., 20000 to 24000).
*
*/


/*****************************************************************************/
/*      	SETUP OF COMMUNICATIONS BUFFERS	AND MESSAGE SIZES            */
/*****************************************************************************/

#ifndef COMM_BUFF_SIZE
#define COMM_BUFF_SIZE 130000
#endif

/* This is the amount of communications buffer size that the
 * broadcast routine will be allowed to use.
 */
#define BRCST_COMM_BUFF_SIZE  (COMM_BUFF_SIZE - 3000)
/*
 * This is the max size of all brdcst writes
 */
#define MAX_MESG_SIZE    8192
/*
 * This is the max size of all fannin messages
 */
#define MSIZE_FANNIN    8192
/*
 * This is number of broadcast messages that are assured to fit in
 * the buffer at one time.
 */
#define MAX_MESG_IN_BUF        (BRCST_COMM_BUFF_SIZE / MAX_MESG_SIZE)


/*****************************************************************************/
/*             			MESSAGE TYPES			             */
/*****************************************************************************/

/* NOTE: messages for the sl_ routines are typed in sl_krylib.h */

/*
 *    Message type for use in the nwrite_big and nread_big routines.
 */

#define MT_BIG_WRITE_ACK     1  /* There is no range to this message type    */

#define SYNC                10	/* sync routine in rf_comm.c		     */
				/* sync messages will have types :
					10, 11,..., 10+dim_hypercube,
				        20, 21 ..., 20+dim_hypercube,
					       ...
				        90, 91 ..., 90+dim_hypercube	     */

#define GSUM_DOUBLE 250  /* total sum of a scalar in rf_comm.c           */
#define GMAX_DOUBLE 2000 /* total search for a max scalar in rf_comm.c   */
#define GSUM_INT    3000 /* total sum of a scalar in rf_comm.c           */

#define COMM_SETUP_1       100	/* local communication setup routine, used   */
#define COMM_SETUP_2       200	/* in exchange_local_info (rf_comm.c)        */


#define ACKNOWLEDGE       3000  /* used in print_global_vec */
#define INT_PRINT_GLOBAL_VEC    3100    /* see rf_comm.c            */
#define DOUBLE_PRINT_GLOBAL_VEC 3200

#define BRDCST            4000  /* logrithmic fan out in comm.c (brdcst) */

#define PRINT_SYNC        5000  /* print_sync_start() and print_sync_end()
				   in rf_comm.c  types:
					5001, 5002, 5099		     */

#define MT_FANNIN         9000  /* Used by fannin routines
                                   Range: 9000 - 9999                        */

#define EXOII_MESH_INT   10000  /* mesh loading message (el_exoII_io.c)      */
#define EXOII_IV_LENGTH  10001  /* mesh loading message (el_exoII_io.c)      */
#define EXOII_MESH_FLOAT 20000  /* mesh loading message (el_exoII_io.c)      */
#define EXOII_FV_LENGTH  20001  /* mesh loading message (el_exoII_io.c)      */




/*****************************************************************************/
/*         PROTOTYPES FOR COMMON MESSAGE PASSING ROUTINES - rf_comm.c        */
/*****************************************************************************/

extern void brdcst        (int node, int Num_Proc, char *buf, int len,
                           int sender);
extern void brdcst_maxlen (const int node, const int Num_Proc,
                           char *buf, int len, int sender);
extern void psync         (int node, int Num_Proc);
extern void print_sync_start (int node, int Dim, int do_print_line);
extern void print_sync_end   (int node, int Dim, int do_print_line);

extern void gather_global_vec (double sol_vec [], int gnodes [], int var_no,
			       int k, int proc, int num_procs, int gindex[],
			       double gvec[], int *N);

extern int  nwrite_big (char *buffer, int nbytes, int dest, int type,
			int *flag);
extern int  nread_big (char *buffer, int nbytes, int *source, int *type,
			int *flag);

/*------------- rf_fanin.c ---------------*/

extern int  fanin_void_slab (int word_size, int gindex[], void *gvec,
                             int gsize, int start_pos, int end_pos,
                             void *gvec_ord);
extern int  fanin_int_slab (int gindex[], int *gvec, int gsize,
                            int start_pos, int end_pos, int *gvec_ord);
extern int fanin_int_void_slab (int word_size, int gindex[], int *gvec,
                                void *fvec, int gsize, int start_pos,
                                int end_pos, int *gvec_ord, void *fvec_ord);
extern int fanin_void_2d_slab (int word_size, int gindex[], void *gvec,
                               int gsize, int block_size, int start_pos,
                               int end_pos, void *gvec_ord);
extern int fanin_int_2d_slab (int gindex[], int **gvec, int gsize,
                              int block_size, int start_pos, int end_pos,
                              int *gvec_ord);
extern int fanin_iblk_slab (int gindex[], int *gvec, int gsize, int gblk_size,
                           int start_pos, int end_pos, int *gvec_ord);
extern int fanin_iblk_ptr_slab (int word_size, int gindex[], int *gvec,
                                void *fvec, int gsize, int gblk_size,
                                int fblk_size, int start_pos, int end_pos,
                                int *gvec_ord, void *fvec_ord);

extern void fanin_float (float *vector, int len, int sum, int proc, int Dim);
extern void fanin_double (double *vector, int len, int sum, int proc, int Dim);
extern void fanin_int (int *vector, int len, int sum, int proc, int Dim);
extern void fanin_str (char *vector, int len, int proc, int Dim);

/*****************************************************************************/
