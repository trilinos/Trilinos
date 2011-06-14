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
#include "exodusII.h"

/*-----------------------------------------------------------------------------
*
*  Include file containing variable definitions used for reading and setting up
*  the mesh information on the processors.  These variables are only used
*  in the program in the file, el_exoII_io.c.  Use of these variables outside
*  this file constitutes an error.  (There are far too many global .h files
*  anyway to have global definitions in this file.)
*
*
*----------------------------------------------------------------------------*/

/*
 *         Information from the Mesh File
 */

static int *Node_Set_Ids = NULL;       /* Vector of node set ids             *
			                * (ptr to vector of length, 	     *
			                *  Num_Node_Set)    		     */
static int *Side_Set_Ids = NULL;       /* Side set ids  		     *
			                * (ptr to vector of length, 	     *
				        *  Num_Side_Set)      		     */
static int *Num_Elem_In_Blk = NULL;    /* Number of elements in each element *
					* block (ptr to vector of length,    *
				        *        Num_Elem_Blk)     	     */
static int *Num_Nodes_Per_Elem = NULL; /* Number of nodes per element in each*
					* elem block (ptr to vector of       *
					* length, Num_Elem_Blk)              */
static int *Num_Attr_Per_Elem = NULL;  /* Number of attributes per element in*
					* each block (ptr to vector of       *
					* length, Num_Elem_Blk)              */
static int *Elem_Blk_Ids = NULL;       /* Element block id's of each element *
					* block (ptr to vector of length,    *
					*        Num_Elem_Blk)               */
static char **Elem_Blk_Types = NULL;   /* Element block types for each       *
					* element block (ptr to vector of    *
					* length, Num_Elem_Blk, of vectors of*
					* char of length MAX_STR_LENGTH+1)   */
static char **Elem_Blk_Names = NULL;   /* Element block names for each       *
					* element block (ptr to vector of    *
					* length, Num_Elem_Blk, of vectors of*
					* char of length MAX_STR_LENGTH+1)   */
static char **Node_Set_Names = NULL;   /* Nodeset names for each             *
					* nodeset (ptr to vector of          *
					* length, Num_Node_Set, of vectors of*
					* char of length MAX_STR_LENGTH+1)   */
static char **Side_Set_Names = NULL;   /* Sideset names for each             *
					* sideset (ptr to vector of          *
					* length, Num_Side_Set, of vectors of*
					* char of length MAX_STR_LENGTH+1)   */
static char ***Elem_Blk_Attr_Names = NULL;  /* Element block attribute names for each 
					     * attribute in each element block
					     * ragged array (size varies for each
					     * element block          	     */

int *GM_Elem_Types = NULL;             /* This is a list of the element      *
                                        * in the entire global mesh. It is   *
                                        * stored on Processor 0 to           *
                                        * facilitate the reading of the side *
                                        * set distribution factors.          */


/*****************************************************************************/
/*          VARIABLES FOR TIMING THE MESH LOADING ROUTINES                   */
/*****************************************************************************/

#ifdef DEBUG_TIME

extern void    sync (int, int);
static double  Start_Time = 0.0;
static double  End_Time   = 0.0;
static double  Time_Exo_Read = 0.0;       /* Time to do large exoII reads    */
static double  Time_File_Ops = 0.0;       /* Time to do file operations      */
static double  Time_Big_Brdcst = 0.0;     /* Time to do large broadcasts     */
static double  Time_Find_Int = 0.0;       /* Time to find intersections      */
static double  Time_Elem_Blk = 0.0;       /* Time to process element block   *
					   * level information               */
static double  Time_Dist_Ids = 0.0;       /* Time to read/distribute/check   *
					   * Ids for Element Blocks, Node    *
					   * Sets, and side sets             */
static double  Time_Set_Up = 0.0;         /* Time to do initial setup
       				           * - involves reading small amounts
                                           *   of exoII info and broadcasting
                                           it */
#endif
