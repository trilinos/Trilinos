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
/*-----------------------------------------------------------------------------
*
*  Include file containing element geometry variable declarations
*----------------------------------------------------------------------------*/

#ifndef NULL
#include <stdio.h>
#endif

/*---------------------------------------------------------------------------*/
/*	GLOBAL QUANTITITES THAT ARE THE SAME ON ALL PROCESSORS		     */
/*---------------------------------------------------------------------------*/

int   Num_Dim = 0;            /* Number of physical dimensions in the problem*/
int   Num_Node = 0;           /* Total number of nodes in the entire mesh    *
		               * - this is a global quantity                 */
int   Num_Elem = 0;           /* Total number of elements in the entire mesh *
		               * - this is a global quantity                 */
int   Num_Elem_Blk = 0;       /* Total number of element blocks              */
int   Num_Node_Set = 0;       /* Total number of node sets defined in the    *
			       * mesh exoII file			     */
int   Num_Side_Set = 0;       /* Total number of side sets defined in the    *
			       * mesh exoII file			     */
int   Max_NP_Elem = 0;        /* Maximum number of nodes in any element      *
		               *  - this is a global quantity                */
int   Num_QA_Recs = 0;	      /* Number of QA records in original file       */
int   Num_Info_Recs = 0;      /* Number of Info records in original file     */

int   Num_Coordinate_Frames = 0;

/*---------------------------------------------------------------------------*/
/*    VARIABLES THAT DEAL WITH SPECIFICATION OF LOAD BALANCE PROPERTIES      */
/*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
/*---------------------------------------------------------------------------*/

int  *Num_Internal_Nodes = NULL; /* Number internal nodes on the local proc  */
int  *Num_Border_Nodes   = NULL;  /* Number border nodes on the local proc   */
int  *Num_External_Nodes = NULL;  /* Number external nodes on the local proc */
int  *Num_Internal_Elems = NULL;  /* Number Elements on the local proc       */
int  *Num_Border_Elems   = NULL;  /* Number Elements on the local proc       *
                                   * but shared by other processors          */

int  *Num_N_Comm_Maps = NULL;   /* Number of nodal communication maps */

int  *Num_E_Comm_Maps = NULL;   /* Number of elemental communication maps */

ELEM_COMM_MAP **E_Comm_Map=NULL;    /* Elemental communication map structure */
NODE_COMM_MAP **N_Comm_Map=NULL;    /* Nodal communication map structure */

int **GNodes = NULL;          /* Data structure which contains the internal, *
			       * border, and external nodes on each processor*
			       * They are structured in that order, and      *
			       * monotonically within each category          *
			       *  Type: int vector of length                 *
 			       *       (Num_Internal_Nodes + Num_Border_Nodes*
			                Num_External_Nodes) 		     */

int **GElems = NULL;          /* Data structure which contains the internal  *
			       * elements on each processor.  It is a map    *
			       * from the local element number to the global *
			       * element number.  			     *
			       *  Type: int vector of length                 *
 			       *        Num_Internal_Elems                   */

int **Proc_Global_Elem_Id_Map = NULL; /* Data structure which contains the internal  *
			       * elements on each processor.  It is a map    *
			       * from the local element number to the global *
			       * element id (as found in elem_num_map)	     *
			       *  Type: int vector of length                 *
 			       *        Num_Internal_Elems                   */

int **Proc_Global_Node_Id_Map = NULL; /* Data structure which contains the internal  *
			       * nodes on each processor.  It is a map    *
			       * from the local node number to the global *
			       * node id (as found in node_num_map)	     *
			       *  Type: int vector of length                 *
 			       *        Num_Internal_Nodes                   */

int **Elem_Map = NULL;	      /* Map for Nemesis output */

int   **GElem_Blks = NULL;    /* Data structure which contains the mapping   *
                               * from the local element block number to the  *
	                       * global element block number                 *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Elem_Blk                    */

int   **GNode_Sets = NULL;    /* Data structure which contains the mapping   *
                               * from the local node set number to the       *
	                       * global node set number                      *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Node_Sets                   */

int   **GSide_Sets = NULL;    /* Data structure which contains the mapping   *
                               * from the local side set number to the       *
	                       * global side set number                      *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Side_Sets                   */

/*---------------------------------------------------------------------------*/
/*    VARIABLES THAT DEAL WITH SPECIFICATION OF NODAL PROPERTIES   	     */
/*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
/*---------------------------------------------------------------------------*/

float  ***Coor_sp = NULL;
double ***Coor_dp = NULL;     /* 2d dynamically allocated array containing   *
                               * the physical space coordinates for each     *
			       * node, defined on the local processor.       *
			       *  Type: double vector of length              *
			       *        Num_Dim   by			     *
 			       *       (Num_Internal_Nodes + Num_Border_Nodes*
			       *        Num_External_Nodes) 		     */

/*---------------------------------------------------------------------------*/
/*    VARIABLES THAT DEAL WITH SPECIFICATION OF ELEMENT PROPERTIES   	     */
/*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
/*---------------------------------------------------------------------------*/

int  **Proc_Connect_Ptr = NULL;/* Vector of pointers to the start of each    *
                                * element in Proc_Elem_Connect               *
                                *  Type: int vector of length                *
                                *        Num_Internal_Elems                  */

int  **Elem_Type = NULL;      /* Vector which contains the element code for  *
 			       * each defined on the current processor       *
			       *  Type: int vector of length                 *
 			       *        Num_Internal_Elems                   */

int **Proc_Elem_Connect = NULL;/* Connectivity lists for the elements        *
                                * required by the current processor          *
                                *  Type: int vector of variable length       */

float  **Proc_Elem_Attr_sp = NULL;
double **Proc_Elem_Attr_dp = NULL; /* Attribute list for the elements        *
                                    * required by the current processor      *
                                    *  Type: float vector of variable length */

/*---------------------------------------------------------------------------*/
/*    VARIABLES THAT DEAL WITH SPECIFICATION OF ELEMENT BLOCK PROPERTIES     */
/*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
/*---------------------------------------------------------------------------*/

int *Proc_Num_Elem_Blk = NULL;  /* Number of element blocks on this proc     */

int **Proc_Num_Elem_In_Blk = NULL;
                              /* Number of elements in the processor's       *
			       * element blocks                              *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Elem_Blk                    */

int **Proc_Elem_Blk_Ids = NULL;
			      /* Element block id's for the processor's      *
			       * element blocks                              *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Elem_Blk                    */

int **Proc_Elem_Blk_Types = NULL;
			      /* Element block types for the processor's     *
			       * element blocks                              *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Elem_Blk                    */

int **Proc_Nodes_Per_Elem = NULL;
			      /* Number of nodes per element for each        *
			       * block on the current processor              *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Elem_Blk                    */

int **Proc_Num_Attr = NULL;   /* Number of attributes for each block on the  *
                               * current processor                           *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Elem_Blk                    */

int  *Elem_Blk_2_Matls = NULL;
                              /* Mapping of element block NUMBERS to material*
			       * IDs  	                	             *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Elem_Blk                    */

/*---------------------------------------------------------------------------*/
/*	VARIABLES THAT DEAL WITH SPECIFICATION OF NODE SETS		     */
/*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
/*---------------------------------------------------------------------------*/

int *Proc_Num_Node_Sets = NULL;  /* Number of node sets on the current proc  */

int *Proc_NS_List_Length = NULL;  /* Total length of all the node set lists  */

int **Proc_NS_Ids = NULL;     /* Node sets ids for the node sets in a given  *
			       * processor                                   *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Node_Sets                   */

int **Proc_NS_Count = NULL;   /* Number of nodes in each node set defined on *
			       * the current processor.                      *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Node_Sets                   */

int **Proc_NS_DF_Count = NULL;
                              /* Number of distribution factors associated   *
                               * with a node set.                            *
                               *  Type: int vector of length                 *
                               *        Proc_Num_Node_Sets                   */

int **Proc_NS_Pointers = NULL; /* Vector of pointers into Proc_NS_List for   *
			       * each node set defined on the current proc   *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Node_Sets                   */

int **Proc_NS_List = NULL;    /* Node sets list record for the nodes sets in *
			       * a given processor                           *
			       *  Type: int vector of length                 *
 			       *        Proc_NS_List_Length                  */

float  **Proc_NS_Dist_Fact_sp = NULL;
double **Proc_NS_Dist_Fact_dp = NULL;
			      /* Node sets distribution factors for the node *
			       * sets on a given processor                   *
			       *  Type: float vector of length               *
 			       *        Proc_NS_List_Length                  */

int **Proc_NS_GNMap_List = NULL;
                              /* Map local node position to global nodeset node list position *
			       * in a given processor                        *
			       *  Type: int vector of length                 *
 			       *        Proc_NS_List_Length             */

/*---------------------------------------------------------------------------*/
/*	VARIABLES THAT DEAL WITH SPECIFICATION OF SIDE SETS		     */
/*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
/*---------------------------------------------------------------------------*/


int *Proc_Num_Side_Sets = NULL;  /* Number of side sets on the current proc  */

int *Proc_SS_Elem_List_Length = NULL;
		              /* Total number of sides in all side sets      */

int **Proc_SS_Ids = NULL;     /* Side sets ids for side sets in a given      *
			       * processor                                   *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Side_Sets                   */
int **Proc_SS_Elem_Count = NULL;
                              /* A count of the number of sides in each side *
                               * set for a given processor.                  *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Side_Sets                   */

int **Proc_SS_DF_Count = NULL;
                              /* Count of the number of distribution         *
                               * factors in this side set.                   *
                               *  Type: int vector of length                 *
                               *        Proc_Num_Side_Sets                   */

int **Proc_SS_Elem_Pointers = NULL;
			      /* Side sets pointer record for elements (see  *
			       * the EXODUS manual) for the side sets in a   *
			       * given processor                             *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Side_Sets                   */

int **Proc_SS_DF_Pointers = NULL;
                              /* Pointer into the distribution factors for   *
                               * the side sets.                              *
                               *  Type: int array of dimensions              *
                               *    Num_Proc * Proc_Num_Side_Sets            */

float  **Proc_SS_Dist_Fact_sp = NULL;
double **Proc_SS_Dist_Fact_dp = NULL;
                              /* Pointer for storage of the distribution     *
                               * factors.                                    */

int **Proc_SS_Elem_List = NULL;
			      /* Side sets element list record for the side  *
			       * sets in a given processor                   *
			       *  Type: int vector of length                 *
 			       *        Proc_Num_Side_Sets                   */

int **Proc_SS_Side_List = NULL;
			      /* Side ID list for the side sets on a given   *
                               * processor                                   *
			       *  Type: int vector of length                 *
 			       *        Proc_SS_Node_List_Length             */

int **Proc_SS_GEMap_List = NULL;
                              /* Map local element position to global sideset element list position *
			       * in a given processor                        *
			       *  Type: int vector of length                 *
 			       *        Proc_SS_Elem_List_Length             */

/*---------------------------------------------------------------------------*/
/*		VARIABLES THAT DEAL WITH GENERAL INFORMATION THAT IS         */
/*			THE SAME ON EVERY PROCESSOR                          */
/*---------------------------------------------------------------------------*/

char **QA_Record = NULL;	/* The QA Records from the original file     */

char **Info_Record = NULL;	/* The Information Records from the original *
                                 * file                                      */

int    *Coordinate_Frame_Ids            = NULL;
float  *Coordinate_Frame_Coordinates_sp = NULL;
double *Coordinate_Frame_Coordinates_dp = NULL;
char   *Coordinate_Frame_Tags           = NULL;

/*---------------------------------------------------------------------------*/
/*			END OF el_geom.h				     */
/*---------------------------------------------------------------------------*/
