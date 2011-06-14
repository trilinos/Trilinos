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
*  rf_salsa.h:
*
*	This include file contains common preprocessor definitions, typedef
* declarations,  and function PROTO declarations for the salsa project.
* This include file is included in all source code files starting with the
* "rf_" and "el_" prefix.
*
*/

#define UTIL_NAME "nem_spread"
#define VER_STR   "5.14 (2009/10/22)"

/*****************************************************************************/
/*		PROTOTYPES FOR COMMON MP FUNCTIONS		  	     */
/*****************************************************************************/

extern int md_write(char *buffer, int nbytes, int dest,    int  type, int *flag);
extern int md_read( char *buffer, int nbytes, int *source, int *type, int *flag);

extern void whoami (int *Proc, int *proc, int *host, int *Dim);

extern void check_exodus_error (int, char *);

/*****************************************************************************/
/*     EXTERN STATEMENTS FOR GLOBAL FUNCTIONS CALLED DIRECTLY BY rf_salsa.c  */
/*****************************************************************************/

extern int    check_inp            (void);
extern void   pre_calc             (void);
extern void   exch_init_info       (void),
              init_info            (void),
              read_input_file      (char *cmd_file_name),
              read_mesh_param      (int *io_ws),
  	      read_restart_params  (int  io_ws),
              load_lb_info         (void),
              load_mesh 	   (int  io_ws),
  	      read_restart_data    (int  io_ws),
              pre_process          (void),
              solve_problem        (int *, int *, int *, int *, double *,
				    double *),
              get_parallel_info    (int *, int *, int *),
              handle_ieee          (void),
              setup_actual_problem (void),
              set_local_map_comm   (void),
              cleanup              (int *indx, int *bindx, int *rpntr,
                                    int *bpntr, double *a, double *x);
extern double second               (void);

/*****************************************************************************/
