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

#ifndef _ELB_INP_CONST_H_
#define _ELB_INP_CONST_H_

#include "elb_const.h"

/* Prototype for command-line parsing function */
extern
int cmd_line_arg_parse(
  int argc, 		/* The command line argument count */
  char *argv[],		/* The command line arguments array */
  char ex_inp_file[],	/* The ExodusII input FEM file name */
  char as_inp_file[],	/* The ASCII input file name */
  char ne_out_file[],	/* The output NemesisI file name */
  MACHINE_PTR machine,	/* Pointer to structure in which to place machine
                         * information */
  LB_INFO_PTR lb,	/* Pointer to structure in which to place load
                         * balance parameters */
  PROB_INFO_PTR prob,	/* Pointer to structure in which to place general
                         * information about the run */
  SOLVE_INFO_PTR sol,	/* Pointer to structure in which to place parameters
                         * for the eigensolver */
  WEIGHT_INFO_PTR wgh	/* Pointer to structure in which to place parameters
                         * for the graph weighting scheme */
  );

/* Prototype for function which reads in the ASCII input file */
extern
int read_cmd_file(
  char as_inp_file[],	/* The ASCII input file name */
  char ex_inp_file[],	/* The ExodusII input FEM file name */
  char ne_out_file[],	/* The output NemesisI file name */
  MACHINE_PTR machine,	/* Pointer to structure in which to place machine
                         * information */
  LB_INFO_PTR lb,	/* Pointer to structure in which to place load
                         * balance parameters */
  PROB_INFO_PTR prob,	/* Pointer to structure in which to place general
                         * information about the run */
  SOLVE_INFO_PTR sol,	/* Pointer to structure in which to place parameters
                         * for the eigensolver */
  WEIGHT_INFO_PTR wght	/* Pointer to structure in which to place parameters
                         * for the eigensolver */
);

/* Prototype for function which checks the user specified input */
extern
int check_inp_specs(
  char ex_inp_file[],	/* The ExodusII input FEM file name */
  char ne_out_file[],	/* The output NemesisI file name */
  MACHINE_PTR machine,	/* Pointer to structure in which to place machine
                         * information */
  LB_INFO_PTR lb,	/* Pointer to structure in which to place load
                         * balance parameters */
  PROB_INFO_PTR prob,	/* Pointer to structure in which to place general
                         * information about the run */
  SOLVE_INFO_PTR sol,	/* Pointer to structure in which to place parameters
                         * for the eigensolver */
  WEIGHT_INFO_PTR wght	/* Pointer to structure in which to place parameters
                         * for the weighting scheme */
);

/* Various defines used by the input routines */
#define NONE -1

#endif /* _ELB_INP_CONST_H_ */
