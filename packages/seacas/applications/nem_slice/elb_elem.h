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

#ifndef _ELB_ELM_CONST_H
#define _ELB_ELM_CONST_H

/* Define element types */
typedef enum {SPHERE, BAR2, BAR3, QUAD4, QUAD8, QUAD9,
              SHELL4, SHELL8, TRI3, TRI6, TSHELL3, TSHELL6, HEX8,
              HEX20, HEX27, HEXSHELL, TET4, TET10, TET8, WEDGE6, WEDGE15,
              WEDGE16, PYRAMID5, PYRAMID13, SHELL2, SHELL3, NULL_EL} E_Type;

extern const char* elem_name_from_enum(const E_Type elem_type);

extern
E_Type get_elem_type(
  const char *elem_name,	/* ExodusII element name */
  const int   num_nodes,	/* Number of nodes in the element */
  const int   num_dim		/* Number of dimensions of the mesh */
);

extern
int get_elem_info(
  const int info,		/* The requested information */
  const E_Type elem_type	/* The element type */
);

template <typename INT>
int get_side_id(const E_Type etype, const INT *connect, const int nsnodes,
                INT side_nodes[], const int skip_check, const int partial_adj);

template <typename INT>
int get_side_id_hex_tet(
  const E_Type  etype,          /* The element type */
  const INT *conn,              /* The element connectivity */
  const int  nsnodes,           /* The number of side nodes */
  const INT  side_nodes[]       /* The list of side node IDs */
);

template <typename INT>
int ss_to_node_list(
  const E_Type  etype,		/* The element type */
  const INT *connect,		/* The element connectivity */
  int  side_num,		/* The element side number */
  INT  ss_node_list[]		/* The list of side node IDs */
);

template <typename INT>
int get_ss_mirror(
  const E_Type etype,		/* The element type */
  const INT *ss_node_list,	/* The list of side node IDs */
  int side_num,			/* The element side number */
  INT mirror_node_list[]	/* The list of the mirror side node IDs */
);


/* Define element info requests */
#define NNODES		0
#define NQUAD		1
#define NDIM		2
#define NQUAD_SURF	3
#define NSIDE_NODES	4
#define NSIDES		5

/* Define for the maximum number of nodes on an element side/face */
#define MAX_SIDE_NODES	9
/*
 * Define for the maximum number of sides (and hence communications
 * entries) that an element can have
 */
#define MAX_ELEM_SIDES	6

#endif /* _ELB_ELM_CONST_H */
