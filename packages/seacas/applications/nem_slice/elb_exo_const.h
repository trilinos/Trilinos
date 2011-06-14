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

#ifndef _ELB_EXO_CONST_H_
#define _ELB_EXO_CONST_H_

#include "elb_const.h"

/* Function prototypes */
extern
int read_exo_weights(
  PROB_INFO_PTR prob,		/* Pointer to problem info structure */
  WEIGHT_INFO_PTR weight	/* Pointer to weight info structure */
);

extern
int read_mesh_params(
  char exo_file[],		/* Name of ExodusII geometry file */
  PROB_INFO_PTR prob,		/* Pointer to problem info structure */
  MESH_INFO_PTR mesh,		/* Mesh information structure */
  SPHERE_INFO_PTR sphere	/* Sphere element info structure */
);

extern
int read_mesh(
  char exo_file[],		/* Name of ExodusII geometry file */
  PROB_INFO_PTR prob,		/* Problem information */
  MESH_INFO_PTR mesh,		/* Mesh information structure */
  WEIGHT_INFO_PTR weight	/* Weight specification structure */
);

extern
int init_weight_struct(PROB_INFO_PTR problem,	/* Problem information */
                       MESH_INFO_PTR mesh,      /* Mesh information structure */
                       WEIGHT_INFO_PTR weight); /* Weight specification structure */

#endif /* _ELB_EXO_CONST_H_ */
