/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


#ifndef _DR_ELM_CONST_H
#define _DR_ELM_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* Define element types */
typedef enum {E_TYPE_ERROR=-1, SPHERE, BAR1, BAR2, QUAD1, S_QUAD2, QUAD2,
              SHELL1, SHELL2, TRI1, TRI2, TSHELL1, TSHELL2, HEX1,
              S_HEX2, HEX2, TET1, TET2, WEDGE1, WEDGE2,
              HEXSHELL, NULL_EL} E_Type;

extern
E_Type get_elem_type(
  char       *elem_name,	/* ExodusII element name */
  const int   num_nodes,	/* Number of nodes in the element */
  const int   num_dim		/* Number of dimensions of the mesh */
);

extern 
const char *get_elem_name(
  int         itype             /* ExodusII element type */
);

extern
int get_elem_info(
  const int info,		/* The requested information */
  const E_Type elem_type,	/* The element type */
  const ZOLTAN_ID_TYPE sid			/* side id (to get number of side nodes) */
);

extern
int get_side_id(
  E_Type     etype,		/* The element type */
  const ZOLTAN_ID_TYPE *conn,		/* The element connectivity */
  const int  nsnodes,	/* The number of side nodes */
  int side_nodes[]	/* The list of side node IDs */
);

extern
int ss_to_node_list(
  const E_Type  etype,		/* The element type */
  const ZOLTAN_ID_TYPE *connect,	/* The element connectivity */
  int  side_num,		 /* The element side number */
  int ss_node_list[]	/* The list of side node IDs */
);

extern
ZOLTAN_ID_TYPE get_ss_mirror(
  const E_Type etype,		/* The element type */
  const ZOLTAN_ID_TYPE *ss_node_list,	/* The list of side node IDs */
  ZOLTAN_ID_TYPE side_num,			/* The element side number */
  ZOLTAN_ID_TYPE mirror_node_list[]	/* The list of the mirror side node IDs */
);


/* Define element info requests */
#define NNODES		0
#define NQUAD		1
#define NDIM		2
#define NQUAD_SURF	3
#define NSNODES		4
#define NSIDES		5

/* Define for the maximum number of nodes on an element side/face */
#define MAX_SIDE_NODES	9
/*
 * Define for the maximum number of sides (and hence communications
 * entries) that an element can have
 */
#define MAX_ELEM_SIDES	6


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
#endif /* _DR_ELM_CONST_H */
