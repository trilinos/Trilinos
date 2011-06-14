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
*
*	Include file for I/O global varibales used in FEM
*	problem specification
*
*/

char ExoFile[MAX_FNL+1];    /* Exodus II File containing problem definition. */
                            /* This name is the root name.                   */
char Exo_LB_File[MAX_FNL+1];/* Exodus II file containing the mesh
                             * load-balanceinformation                       */
char Exo_Res_File[MAX_FNL+1]; /* Exodus II file containing the mesh results  */
                              /* information                                   */

int  Debug_Flag = 1;	    /* Flag to specify debug info is to be printed out.
			       The value of this flag determines the level of
			       diagnostic info which is printed to stdout
			       Debug_Flag == 0 	No debug output
			                     .
					     .
					     9	maximum output               */
int Gen_Flag = 1;           /* Flag used by nem_join to determine if the user
                               wants to use an existing genesis file rather
                               than generating one from the parallel files */

int Num_Nod_Var  = 0;		/* The number of nodal variables to reserve */
				/* space in the output file for. */
int Num_Elem_Var = 0;		/* The number of elemental variables to */
				/* reserve space in the output file for. */
int Num_Glob_Var = 0;		/* The number of global variables to reserve */
				/* space in the output file for. */
int Num_Nset_Var = 0;		/* The number of nodeset variables to reserve */
				/* space in the output file for. */
int Num_Sset_Var = 0;		/* The number of sideset variables to reserve */
				/* space in the output file for. */
RESTART Restart_Info;		/* Restart information structure */
