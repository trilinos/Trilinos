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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "elb_const.h"
#include "elb_elem_const.h"
#include "elb_err_const.h"
#include "elb_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function get_elem_type() begins:
 *----------------------------------------------------------------------------
 * This function returns the type of element based on the ExodusII element
 * string and number of nodes.
 *
 * Need the number of dimensions in order to distinguish between
 * TRI elements in a 2d mesh from TRI elements in a 3d mesh.
 *****************************************************************************/
E_Type get_elem_type(const char *elem_name, const int num_nodes,
                     const int num_dim)
{

  E_Type answer = NULL_EL;
  switch (elem_name[0]) {
  case 'h':
  case 'H':
  if(strncasecmp(elem_name, "HEX", 3) == 0)
  {
    switch(num_nodes)
    {
    case 8:
      answer = HEX8;
      break;
    case 12:
      answer = HEXSHELL;
      break;
    case 20:
      answer = HEX20;
      break;
    case 27:
      answer = HEX27;
      break;
    default:
      Gen_Error(0, "fatal: unsupported HEX element");
      error_report();
      exit(1);
    }
  }
  break;

  case 's':
  case 'S':
  if(strncasecmp(elem_name, "SPHERE", 6) == 0)
    answer = SPHERE;
  else if(strncasecmp(elem_name, "SHELL", 5) == 0)
  {
    switch(num_nodes)
    {
    case 2:
      if(num_dim == 2) answer = SHELL2;
      else {
        Gen_Error(0, "fatal: unsupported SHELL element");
        error_report();
        exit(1);
      }
      break;
    case 3:
      if(num_dim == 2) answer = SHELL3;
      else {
        Gen_Error(0, "fatal: unsupported SHELL element");
        error_report();
        exit(1);
      }
      break;
    case 4:
      answer = SHELL4;
      break;
    case 8:
      answer = SHELL8;
      break;
    default:
      Gen_Error(0, "fatal: unsupported SHELL element");
      error_report();
      exit(1);
    }
  }
  break;
  
  case 'b':
  case 'B':
  case 't':
  case 'T':
  case 'r':
  case 'R':
  if(strncasecmp(elem_name, "BEAM", 4) == 0 ||
     strncasecmp(elem_name, "TRUSS", 5) == 0 ||
     strncasecmp(elem_name, "ROD", 3) == 0 ||
     strncasecmp(elem_name, "BAR", 3) == 0)
  {
    switch(num_nodes)
    {
    case 2:
      answer=BAR2;
      break;
    case 3:
      answer=BAR3;
      break;
    default:
      Gen_Error(0, "fatal: unsupported BAR/BEAM/TRUSS element");
      error_report();
      exit(1);
    }
  }
  else if(strncasecmp(elem_name, "TRI", 3) == 0)
  {
    switch(num_nodes)
    {
    case 3:
      if (num_dim == 2) answer = TRI3;
      else              answer = TSHELL3;
      break;
    case 6:
      if (num_dim == 2) answer = TRI6;
      else              answer = TSHELL6;
      break;
    default:
      Gen_Error(0, "fatal: unsupported TRI element");
      error_report();
      exit(1);
    }
  }
  else if(strncasecmp(elem_name, "TET", 3) == 0)
  {
    switch(num_nodes)
    {
    case 4:
      answer = TET4;
      break;
    case 8:
      answer = TET8;
      break;
    case 10:
      answer = TET10;
      break;
    default:
      Gen_Error(0, "fatal: unsupported TET element");
      error_report();
      exit(1);
    }
  }
  break;
  
  case 'q':
  case 'Q':
  if(strncasecmp(elem_name, "QUAD", 4) == 0)
  {
    switch(num_nodes)
    {
    case 4:
      if(num_dim == 2) answer = QUAD4;
      else             answer = SHELL4;
      break;
    case 8:
      if(num_dim == 2) answer = QUAD8;
      else             answer = SHELL8;
      break;
    case 9:
      if(num_dim == 2) answer = QUAD9;
      else {
        Gen_Error(0, "fatal: unsupported SHELL element");
        error_report();
        exit(1);
      }
      break;
    default:
      Gen_Error(0, "fatal: unsupported QUAD element");
      error_report();
      exit(1);
    }
  }
  break;

  case 'w':
  case 'W':
  if(strncasecmp(elem_name, "WEDGE", 5) == 0)
  {
    switch(num_nodes)
    {
    case 6:
      answer = WEDGE6;
      break;
    case 15:
      answer = WEDGE16;
      break;
    case 16:
      answer = WEDGE15;
      break;
    default:
      Gen_Error(0, "fatal: unsupported WEDGE element");
      error_report();
      exit(1);
    }
  }
  break;
  
  case 'p':
  case 'P':
  if(strncasecmp(elem_name, "PYR", 3) == 0)
  {
    switch(num_nodes)
    {
    case 5:
      answer = PYRAMID5;
      break;
    case 13:
      answer = PYRAMID13;
      break;
    default:
      Gen_Error(0, "fatal: unsupported PYRAMID element");
      error_report();
      exit(1);
    }
  }
  break;
  
  default:
    break;
  }
  
  if (answer == NULL_EL) {
    Gen_Error(0, "fatal: unknown element type read");
    error_report();
    exit(1);
  }
  
  return answer;

} /*---------------------------End get_elem_type()---------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function get_elem_info() begins:
 *----------------------------------------------------------------------------
 * This function returns various information about the input element type.
 *****************************************************************************/
int get_elem_info(const int req, const E_Type etype)
{

  int answer=0;

  switch(etype)		/* Switch over the element type */
  {
  case BAR2:
    switch(req)
    {
    case NNODES:
      answer = 2;
      break;
    case NSIDE_NODES:
      answer = 2;
      break;
    case NSIDES:
      answer = 1;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 1;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case SHELL2:
    switch(req)
    {
    case NNODES:
      answer = 2;
      break;
    case NSIDE_NODES:
      answer = 2;
      break;
    case NSIDES:
      answer = 1;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 1;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case SHELL3:
    switch(req)
    {
    case NNODES:
      answer = 3;
      break;
    case NSIDE_NODES:
      answer = 2;
      break;
    case NSIDES:
      answer = 1;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 1;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case BAR3:
    switch(req)
    {
    case NNODES:
      answer = 3;
      break;
    case NSIDE_NODES:
      answer = 2;
      break;
    case NSIDES:
      answer = 1;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 1;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case SPHERE:
    switch(req)
    {
    case NNODES:
      answer = 1;
      break;
    case NSIDE_NODES:
      answer = 0;
      break;
    case NSIDES:
      answer = 0;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    }
    break;

  case QUAD4:		/* First order quad */
    switch(req)		/* select type of information required*/
    {
    case NNODES:	/* number of nodes */
      answer = 4;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NSIDE_NODES:
      answer = 2;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal:unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case QUAD8:		/* 2nd order serendipity quad */
    switch(req)		/* select type of information required */
    {
    case NNODES:		/* number of nodes */
      answer = 8;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NSIDE_NODES:
      answer = 3;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case QUAD9:	/* biquadratic quadrilateral */
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 9;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NSIDE_NODES:
      answer = 3;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  /* NOTE: cannot determine NSIDE_NODES for SHELL element */
  case SHELL4:
    switch(req)
    {
    case NNODES:
      answer = 4;
      break;
    case NSIDES:
      answer = 6;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case SHELL8:
    switch(req)
    {
    case NNODES:
      answer = 8;
      break;
    case NSIDES:
      answer = 6;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TRI3:
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 3;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NSIDE_NODES:
      answer = 2;
      break;
    case NSIDES:
      answer = 3;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TRI6:
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 6;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NSIDE_NODES:
      answer = 3;
      break;
    case NSIDES:
      answer = 3;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  /* NOTE: cannot determine NSIDE_NODES for TSHELL element */
  case TSHELL3:
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 3;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NSIDES:
      answer = 5;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TSHELL6:
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 6;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NSIDES:
      answer = 5;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case HEX8:		/* trilinear hexahedron */
    switch(req)	/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 8;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSIDE_NODES:
      answer = 4;
      break;
    case NSIDES:
      answer = 6;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case HEX20:		/* serendipity triquadratic hexahedron */
    switch(req)	/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 20;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSIDE_NODES:
      answer = 8;
      break;
    case NSIDES:
      answer = 6;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case HEX27:		/* triquadratic hexahedron */
    switch(req)		/* select type of information required*/
    {
    case NNODES:	/* number of nodes */
      answer = 27;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSIDE_NODES:
      answer = 9;
      break;
    case NSIDES:
      answer = 6;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  /* NOTE: cannot determine NSIDE_NODES for HEXSHELL element */
  case HEXSHELL:
    switch(req)
    {
    case NNODES:
      answer = 12;
      break;
    case NSIDES:
      answer = 6;
      break;
    case NDIM:          /* number of physical dimensions */
      answer = 3;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TET4:		/* trilinear tetrahedron */
    switch(req)		/* select type of information required*/
    {
    case NNODES:	/* number of nodes */
      answer = 4;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSIDE_NODES:
      answer = 3;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TET10:		/* triquadradic tetrahedron */
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 10;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSIDE_NODES:
      answer = 6;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TET8:		/* 8-node (midface nodes) tetrahedron */
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 8;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSIDE_NODES:
      answer = 4;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  /* NOTE: cannot determine NSIDE_NODES for WEDGE elements */
  case WEDGE6:
    switch(req)
    {
    case NNODES:
      answer = 6;
      break;
    case NSIDES:
      answer = 5;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case WEDGE15:
    switch(req)
    {
    case NNODES:
      answer = 16;
      break;
    case NSIDES:
      answer = 5;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case WEDGE16:
    switch(req)
    {
    case NNODES:
      answer = 15;
      break;
    case NSIDES:
      answer = 5;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  /* NOTE: cannot determine NSIDE_NODES for PYRAMID element */
  case PYRAMID5:
    switch(req)
    {
    case NNODES:
      answer = 5;
      break;
    case NSIDES:
      answer = 5;
      break;
    case NDIM:          /* number of physical dimensions */
      answer = 3;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case PYRAMID13:
    switch(req)
    {
    case NNODES:
      answer = 13;
      break;
    case NSIDES:
      answer = 5;
      break;
    case NDIM:          /* number of physical dimensions */
      answer = 3;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  default:
    Gen_Error(0, "fatal: unknown or unimplemented element type");
    error_report();
    exit(1);
  }

  return answer;

} /*---------------------------End get_elem_info()---------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function get_side_id() begins:
 *----------------------------------------------------------------------------
 * This function returns the Side ID (as used in ExodusII) given a list of
 * nodes on that side.
 *
 * Changed so that it is now order dependent, but independent of starting
 * node for 3-D sides. On 2-D sides (lines), the starting node is important.
 *
 * Now supoports degenrate faces in HEX elements.
 *****************************************************************************/
int get_side_id(const E_Type etype, const int *connect, const int nsnodes,
                int side_nodes[], const int skip_check, const int partial_adj)
{
  char *func_name="get_side_id";
  char  err_buff[300];

  int nnodes, i, j, num;
  int dup, location[9];
  int count;
  /*  min_match for hex elements means that min_match+1 nodes 
      on a face of a hex must match to return the side of the
      hex on which the nodes exist, i.e., if 3/4 nodes of a hex
      match, then it might be considered connected.  If it is 
      connected, then min_match states if 3/4 nodes is considered a face of 
      a hex or not */

  /*  Default for min_match is 3, 2 is trial and error stuff */

  const int min_match = 3; 
  /* const int min_match = 2; */

  /* check if this is a degenerate face */
  dup = 0;
  for (i = 0; i < (nsnodes - 1); i++) {
    for (j = (i + 1); j < nsnodes; j++) {
      if (side_nodes[i] == side_nodes[j]) {
        location[dup++] = i; /* location of duplicated node */
      }
    }
  }

  nnodes = get_elem_info(NNODES, etype);

  /* Find all of the side nodes in the connect table */
  num = 0;
  for(i=0; i < nnodes; i++)
  {
    for(j=0; j < nsnodes; j++)
    {
      if(connect[i] == side_nodes[j])
      {
        num++;
        break;
      }
    }
    if(num == nsnodes)
      break;
  }

  /* printf("%s::%d num = %d and nsnodes = %d\n",__FILE__,__LINE__,num,nsnodes); */

  /* I commented out the conditional statement causing the 
     error if 2 hexes only share 3 out of 4 nodes.  I replaced
     this with what is seen below.  It works, but only for
     this particular case */


  /* the following ifdef is used to determine face adjacency 
     old way:  numnodes on face must match on both elements
     new way:  only 3/4 of hex nodes have to match to be face adjacent */

  if (((partial_adj == 1) && (num < nsnodes - 1) &&  (num >= 2)) ||
      ((partial_adj != 1) && (num != nsnodes )))
  {
    if (skip_check)
      Gen_Error(0, "warning: not all side nodes in connect table for element");
    else
    {
      Gen_Error(0, "fatal: not all side nodes in connect table for element");
      return -1;
    }
  }

  if ((partial_adj == 1) && ( num != nsnodes ))
    {
      return 0;
    }

  /* Find the side ID */
  switch(etype)
  {
  case BAR2: 
  case BAR3:
  case SHELL2:
  case SHELL3:
    /* SIDE 1 */
    if (side_nodes[0] == connect[0] &&
        side_nodes[1] == connect[1]) return 1;
    break;
  case QUAD4:
  case QUAD8:
  case QUAD9:
    /* SIDE 1 */
    if (side_nodes[0] == connect[0] &&
        side_nodes[1] == connect[1]) return 1;

    /* SIDE 2 */
    if (side_nodes[0] == connect[1] &&
        side_nodes[1] == connect[2]) return 2;

    /* SIDE 3 */
    if (side_nodes[0] == connect[2] &&
        side_nodes[1] == connect[3]) return 3;

    /* SIDE 4 */
    if (side_nodes[0] == connect[3] &&
        side_nodes[1] == connect[0]) return 4;

    break;

  case TRI3:
  case TRI6:
    /* SIDE 1 */
    if (side_nodes[0] == connect[0] &&
        side_nodes[1] == connect[1]) return 1;

    /* SIDE 2 */
    if (side_nodes[0] == connect[1] &&
        side_nodes[1] == connect[2]) return 2;

    /* SIDE 3 */
    if (side_nodes[0] == connect[2] &&
        side_nodes[1] == connect[0]) return 3;

    break;

  case TET4:
  case TET10:
  case TET8:
    /* check the # of side nodes */
    if (nsnodes < 3) return 0;

    /* SIDE 1 */
    if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 3] == connect[1] &&
          side_nodes[(2 + num) % 3] == connect[3]) return 1;
    }

    /* SIDE 2 */
    if((num = in_list(connect[1], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 3] == connect[2] &&
          side_nodes[(2 + num) % 3] == connect[3]) return 2;
    }

    /* SIDE 3 */
    if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 3] == connect[3] &&
          side_nodes[(2 + num) % 3] == connect[2]) return 3;
    }

    /* SIDE 4 */
    if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 3] == connect[2] &&
          side_nodes[(2 + num) % 3] == connect[1]) return 4;
    }

    break;

  case HEX8:
  case HEX20:
  case HEX27:
  case HEXSHELL:  /* this should be the same as a HEX element */
    /* check the # of side nodes */
    if (nsnodes < 4) return 0;

    /* SIDE 1 */
    if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
      count = 0;
      count += numbermatch(side_nodes,1,num,4,connect[1]);
      count += numbermatch(side_nodes,2,num,4,connect[5]);
      count += numbermatch(side_nodes,3,num,4,connect[4]);
      if ( count >= min_match ) return 1;  

      /* if this is the duplicated node, then find the next occurence */
      if (dup) {
	for (i=0; i < dup; i++) {
	  if (connect[0] == side_nodes[location[i]]) {
	    num = in_list(connect[0], (nsnodes-num), &(side_nodes[num+1])) + location[i] + 1;
	    count = 0;
	    count += numbermatch(side_nodes,1,num,4,connect[1]);
	    count += numbermatch(side_nodes,2,num,4,connect[5]);
	    count += numbermatch(side_nodes,3,num,4,connect[4]); 
	    if (count >=min_match ) return 1;       
	  }
	}
      }
    }

    /* SIDE 2 */
    if((num = in_list(connect[1], nsnodes, side_nodes)) >= 0) {
      count = 0;
      count += numbermatch(side_nodes,1,num,4,connect[2]);
      count += numbermatch(side_nodes,2,num,4,connect[6]);
      count += numbermatch(side_nodes,3,num,4,connect[5]);
      if (count >= min_match) return 2;
      
      /* if this is the duplicated node, then find the next occurence */
      if (dup) {
	for (i=0; i < dup; i++) {
	  if (connect[1] == side_nodes[location[i]]) {
	    num = in_list(connect[1], (nsnodes-num), &(side_nodes[num+1])) + location[i] + 1;
	    count = 0;
	    count += numbermatch(side_nodes,1,num,4,connect[2]);
	    count += numbermatch(side_nodes,2,num,4,connect[6]);
	    count += numbermatch(side_nodes,3,num,4,connect[5]);
	    if (count >= min_match) return 2;
	  }
	}
      }
    }

    /* SIDE 3 */
    if((num = in_list(connect[2], nsnodes, side_nodes)) >= 0) {
      count = 0;
      count += numbermatch(side_nodes,1,num,4,connect[3]);
      count += numbermatch(side_nodes,2,num,4,connect[7]);
      count += numbermatch(side_nodes,3,num,4,connect[6]);
      if (count >= min_match) return 3;
   
      /* if this is the duplicated node, then find the next occurence */
      if (dup) {
	for (i=0; i < dup; i++) {
	  if (connect[2] == side_nodes[location[i]]) {
	    num = in_list(connect[2], (nsnodes-num), &(side_nodes[num+1])) + location[i] + 1;
	    count = 0;
	    count += numbermatch(side_nodes,1,num,4,connect[3]);
	    count += numbermatch(side_nodes,2,num,4,connect[7]);
	    count += numbermatch(side_nodes,3,num,4,connect[6]);
	    if (count >= min_match) return 3;
	  }
	}
      }
    }

    /* SIDE 4 */
    if((num = in_list(connect[3], nsnodes, side_nodes)) >= 0) {
      count = 0;
      count += numbermatch(side_nodes,1,num,4,connect[0]);
      count += numbermatch(side_nodes,2,num,4,connect[4]);
      count += numbermatch(side_nodes,3,num,4,connect[7]);
      if (count >= min_match) return 4;

      /* if this is the duplicated node, then find the next occurence */
      if (dup) {
	for (i=0; i < dup; i++) {
	  if (connect[3] == side_nodes[location[i]]) {
	    num = in_list(connect[3], (nsnodes-num), &(side_nodes[num+1])) + location[i] + 1;
	    count = 0;
	    count += numbermatch(side_nodes,1,num,4,connect[0]);
	    count += numbermatch(side_nodes,2,num,4,connect[4]);
	    count += numbermatch(side_nodes,3,num,4,connect[7]);
	    if (count >= min_match) return 4;
	  }
	}
      }
    }

    /* SIDE 5 */
    if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
      count = 0;
      count += numbermatch(side_nodes,1,num,4,connect[3]);
      count += numbermatch(side_nodes,2,num,4,connect[2]);
      count += numbermatch(side_nodes,3,num,4,connect[1]);
      if (count >= min_match) return 5;

      /* if this is the duplicated node, then find the next occurence */
      if (dup) {
	for (i=0; i < dup; i++) {
	  if (connect[0] == side_nodes[location[i]]) {
	    num = in_list(connect[0], (nsnodes-num), &(side_nodes[num+1])) + location[i] + 1;
	    count = 0;
	    count += numbermatch(side_nodes,1,num,4,connect[3]);
	    count += numbermatch(side_nodes,2,num,4,connect[2]);
	    count += numbermatch(side_nodes,3,num,4,connect[1]);
	    if (count >= min_match) return 5;
	  }
	}
      } 
    }

    /* SIDE 6 */
    if((num = in_list(connect[4], nsnodes, side_nodes)) >= 0) {
      count = 0;
      count += numbermatch(side_nodes,1,num,4,connect[5]);
      count += numbermatch(side_nodes,2,num,4,connect[6]);
      count += numbermatch(side_nodes,3,num,4,connect[7]);
      if (count >= min_match) return 6;

      /* if this is the duplicated node, then find the next occurence */
      if (dup) {
	for (i=0; i < dup; i++) {
	  if (connect[4] == side_nodes[location[i]]) {
	    num = in_list(connect[4], (nsnodes-num), &(side_nodes[num+1])) + location[i] + 1;
	    count = 0;
	    count += numbermatch(side_nodes,1,num,4,connect[5]);
	    count += numbermatch(side_nodes,2,num,4,connect[6]);
	    count += numbermatch(side_nodes,3,num,4,connect[7]);
	    if (count >= min_match) return 6;
	  }
	}
      }
    }

    break;

  case SHELL4:
  case SHELL8:

    /* 2D sides */
    if(nsnodes == 2 || nsnodes == 3) {
      /* SIDE 3 */
      if (side_nodes[0] == connect[0] &&
          side_nodes[1] == connect[1]) return 3;

      /* SIDE 4 */
      if (side_nodes[0] == connect[1] &&
          side_nodes[1] == connect[2]) return 4;

      /* SIDE 5 */
      if (side_nodes[0] == connect[2] &&
          side_nodes[1] == connect[3]) return 5;

      /* SIDE 6 */
      if (side_nodes[0] == connect[3] &&
          side_nodes[1] == connect[0]) return 6;
    }

    /* 3D faces */
    else if (nsnodes == 4 || nsnodes == 8) {

      /* SIDE 1 */
      if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == connect[1] &&
            side_nodes[(2 + num) % 4] == connect[2] &&
            side_nodes[(3 + num) % 4] == connect[3]) return 1;
      }

      /* SIDE 2 */
      if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == connect[3] &&
            side_nodes[(2 + num) % 4] == connect[2] &&
            side_nodes[(3 + num) % 4] == connect[1]) return 2;
      }
    }

    break;

  case WEDGE6:
  case WEDGE15:
  case WEDGE16:

    /* quad sides */
    if (nsnodes == 4 || nsnodes == 8) {
      /* SIDE 1 */
      if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == connect[1] &&
            side_nodes[(2 + num) % 4] == connect[4] &&
            side_nodes[(3 + num) % 4] == connect[3]) return 1;
      }

      /* SIDE 2 */
      if((num = in_list(connect[1], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == connect[2] &&
            side_nodes[(2 + num) % 4] == connect[5] &&
            side_nodes[(3 + num) % 4] == connect[4]) return 2;
      }

      /* SIDE 3 */
      if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == connect[3] &&
            side_nodes[(2 + num) % 4] == connect[5] &&
            side_nodes[(3 + num) % 4] == connect[2]) return 3;
      }
    }

    /* triangle sides */
    else if (nsnodes == 3 || nsnodes == 6) {
      /* SIDE 4 */
      if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == connect[2] &&
            side_nodes[(2 + num) % 3] == connect[1]) return 4;
      }

      /* SIDE 5 */
      if((num = in_list(connect[3], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == connect[4] &&
            side_nodes[(2 + num) % 3] == connect[5]) return 5;
      }
    }

    break;

  case TSHELL3:
  case TSHELL6:

    /* 2D sides */
    if(nsnodes == 2 || (etype == TSHELL6 && nsnodes == 3)) {
      /* SIDE 3 */
      if (side_nodes[0] == connect[0] &&
          side_nodes[1] == connect[1]) return 3;

      /* SIDE 4 */
      if (side_nodes[0] == connect[1] &&
          side_nodes[1] == connect[2]) return 4;

      /* SIDE 5 */
      if (side_nodes[0] == connect[2] &&
          side_nodes[1] == connect[0]) return 5;

    }

    /* 3D faces */
    else if (nsnodes == 3 || nsnodes == 6) {

      /* SIDE 1 */
      if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == connect[1] &&
            side_nodes[(2 + num) % 3] == connect[2]) return 1;
      }

      /* SIDE 2 */
      if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == connect[2] &&
            side_nodes[(2 + num) % 3] == connect[1]) return 2;
      }
    }

    break;

  case PYRAMID5:
  case PYRAMID13:
    /* triangular sides */
    if(nsnodes == 3 || nsnodes == 6) {
      /* SIDE 1 */
      if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == connect[1] &&
            side_nodes[(2 + num) % 3] == connect[4]) return 1;
      }

      /* SIDE 2 */
      if((num = in_list(connect[1], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == connect[2] &&
            side_nodes[(2 + num) % 3] == connect[4]) return 2;
      }

      /* SIDE 3 */
      if((num = in_list(connect[2], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == connect[3] &&
            side_nodes[(2 + num) % 3] == connect[4]) return 3;
      }

      /* SIDE 4 */
      if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == connect[4] &&
            side_nodes[(2 + num) % 3] == connect[3]) return 4;
      }
    }

    else if (nsnodes == 4 || nsnodes == 8) {
      /* SIDE 5 */
      if((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == connect[3] &&
            side_nodes[(2 + num) % 4] == connect[2] &&
            side_nodes[(3 + num) % 4] == connect[1]) return 5;
      }
    }

    break;

  case SPHERE:
    break;

  default:
    sprintf(err_buff, "fatal: unknown element type %d in function %s",
            etype, func_name);
    Gen_Error(0, err_buff);
    error_report();
    exit(1);

  } /* End "switch(etype)" */

  return 0;
} /*---------------------------End get_side_id()-----------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function get_side_id_hex_tet() begins:
 *----------------------------------------------------------------------------
 * This function returns the Side ID (as used in ExodusII) of a HEX element
 * given a list of at least 3 nodes on that side. This function works on
 * the fact that any three nodes define a side of a hex. When a tet is
 * connected to a side of a hex, there are only three nodes connecting
 * the two. In this case a side id can be found.
 *****************************************************************************/
int get_side_id_hex_tet(const E_Type etype, const int *connect,
                        const int nsnodes, const int side_nodes[])
{
  char *func_name="get_side_id_hex";
  char  err_buff[300];

  int nnodes, lcnt, i1, i2;
  int loc_node_ids[MAX_SIDE_NODES];

  nnodes = get_elem_info(NNODES, etype);

  /* Find the local node numbers for nodes forming the side */
  lcnt = 0;
  for(i1=0; i1 < nnodes; i1++)
  {
    for(i2=0; i2 < nsnodes; i2++)
    {
      if(connect[i1] == side_nodes[i2])
      {
        loc_node_ids[lcnt++] = i1+1;
        break;
      }
    }
    if(lcnt == nsnodes)
      break;
  }

  switch (etype) {
  case TET4:
  case TET10:
  case TET8:
    /* SIDE 1 */
    if(in_list(1, lcnt, loc_node_ids) >= 0 &&
       in_list(2, lcnt, loc_node_ids) >= 0 &&
       in_list(4, lcnt, loc_node_ids) >= 0) return 1;

    /* SIDE 2 */
    if(in_list(2, lcnt, loc_node_ids) >= 0 &&
       in_list(3, lcnt, loc_node_ids) >= 0 &&
       in_list(4, lcnt, loc_node_ids) >= 0) return 2;

    /* SIDE 3 */
    if(in_list(1, lcnt, loc_node_ids) >= 0 &&
       in_list(3, lcnt, loc_node_ids) >= 0 &&
       in_list(4, lcnt, loc_node_ids) >= 0) return 3;

    /* SIDE 4 */
    if(in_list(1, lcnt, loc_node_ids) >= 0 &&
       in_list(2, lcnt, loc_node_ids) >= 0 &&
       in_list(3, lcnt, loc_node_ids) >= 0) return 4;

    break;

  case HEX8:
  case HEX20:
  case HEX27:
    /* SIDE 1 */
    nnodes = 0;
    if(in_list(1, lcnt, loc_node_ids) >= 0) nnodes++; 
    if(in_list(2, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(5, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(6, lcnt, loc_node_ids) >= 0) nnodes++;
    if (nnodes > 2) return 1;

    /* SIDE 2 */
    nnodes = 0;
    if(in_list(2, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(3, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(6, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(7, lcnt, loc_node_ids) >= 0) nnodes++;
    if (nnodes > 2) return 2;

    /* SIDE 3 */
    nnodes = 0;
    if(in_list(3, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(4, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(7, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(8, lcnt, loc_node_ids) >= 0) nnodes++;
    if (nnodes > 2) return 3;

    /* SIDE 4 */
    nnodes = 0;
    if(in_list(1, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(4, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(5, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(8, lcnt, loc_node_ids) >= 0) nnodes++;
    if (nnodes > 2) return 4;

    /* SIDE 5 */
    nnodes = 0;
    if(in_list(1, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(2, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(3, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(4, lcnt, loc_node_ids) >= 0) nnodes++;
    if (nnodes > 2) return 5;

    /* SIDE 6 */
    nnodes = 0;
    if(in_list(5, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(6, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(7, lcnt, loc_node_ids) >= 0) nnodes++;
    if(in_list(8, lcnt, loc_node_ids) >= 0) nnodes++;
    if (nnodes > 2) return 6;

    break;

  default:
    sprintf(err_buff, "fatal: unknown element type %d in function %s",
            etype, func_name);
    Gen_Error(0, err_buff);
    error_report();
    exit(1);

  } /* End "switch(etype)" */

  return 0;
} /*-------------------------End get_side_id_hex()---------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function ss_to_node_list() begins:
 *----------------------------------------------------------------------------
 * This function returns the list of nodes in a side of an element given
 * the element type, and the side id. It also returns the number of nodes
 * in that side.
 *****************************************************************************/
int ss_to_node_list(const E_Type etype, const int *connect, int side_num,
                    int ss_node_list[])

{
  int i=0;

  /*
   * This function returns a list of global node numbers forming a
   * side set.
   */

  /* triangle */
  static int tri_table[3][3] = {
  /*   1        2        3                                            side   */
    {1,2,4}, {2,3,5}, {3,1,6}                                      /* nodes  */
  };

  /* tshell */
  static int tshell_table[2][6] = {
  /*        1                  2                                      side   */
    {1,2,3,4,5,6,}, {1,3,2,6,5,4}                                  /* nodes  */
  };

  /* quad */
  static int quad_table[4][3] = {
  /*   1        2        3        4                                   side   */
    {1,2,5}, {2,3,6}, {3,4,7}, {4,1,8}                             /* nodes  */
  };

  /* shell */
  static int shell_table[2][8] = {
  /*        1                  2                                      side   */
    {1,2,3,4,5,6,7,8}, {1,4,3,2,8,7,6,5}                           /* nodes  */
  };

  /* tetra */
  static int tetra_table[4][6] = {
  /*      1              2               3               4            side   */
    {1,2,4,5,9,8}, {2,3,4,6,10,9}, {1,4,3,8,10,7}, {1,3,2,7,6,5}   /* nodes  */
  };

  /* wedge */
  static int wedge_table[5][8] = {
  /*        1                     2                     3             side   */
    {1,2,5,4,7,11,13,10}, {2,3,6,5,8,12,14,11}, {1,4,6,3,10,15,12,9},
  /*        4                  5                                      side   */
    {1,3,2,9,8,7,0,0}, {4,5,6,13,14,15,0,0}                        /* nodes  */
  };

  /* hex */
  static int hex_table[6][9] = {
  /*        1                     2                                   side   */
    {1,2,6,5,9,14,17,13,26},  {2,3,7,6,10,15,18,14,25},
  /*        3                     4                                   side   */
    {3,4,8,7,11,16,19,15,27}, {4,1,5,8,13,20,16,12,24},
  /*        5                     6                                   side   */
    {1,4,3,2,12,11,10,9,22},  {5,6,7,8,17,18,19,20,23}                /*nodes*/
  };

  /* hexshell */
  static int hexshell_table[6][6] = {
  /*      1               2                3                4         side   */
    {1,2,6,5,10,9}, {2,3,7,6,11,10}, {3,4,8,7,12,11}, {4,1,5,8,9,12},
  /*      5               6                                           side   */
    {1,4,3,2,0,0},  {5,6,7,8,0,0}                                     /*nodes*/
  };

  /* pyramid */
  static int pyramid_table[5][8] = {
  /*          1                   2                    3              side   */
    {1,2,5,6,11,10,0,0}, {2,3,5,7,12,11,0,0}, {3,4,5,8,13,12,0,0},
  /*          4                  5                                    side   */
    {1,5,4,10,13,9,0,0}, {1,4,3,2,9,8,7,6}                         /* nodes  */
  };

  static int bar_table[1][3] = { {1, 2, 3} };

/* {2, 0} , {1,2} */
/***************************** execution begins ******************************/

  /* Locally decrement side_num */

  side_num--;

  /* Switch over the element type. */
  switch (etype) {
  case BAR2:
  case SHELL2:
    /* Bar1 has 1 side */
        for (i=0;i<2;i++)
           ss_node_list[i] = connect[(bar_table[side_num][i]-1)];
    break;

  case BAR3:
  case SHELL3:
    /* Bar has 1 side */
        for (i=0;i<3;i++)
	  ss_node_list[i] = connect[(bar_table[side_num][i]-1)];
    break;
    
  case QUAD4:
    for (i = 0; i < 2; i++)
      ss_node_list[i] = connect[(quad_table[side_num][i] - 1)];
    break;

  case QUAD8:
  case QUAD9:
    for (i = 0; i < 3; i++)
      ss_node_list[i] = connect[(quad_table[side_num][i] - 1)];
    break;

  case SHELL4:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 4; i++)
        ss_node_list[i] = connect[(shell_table[side_num][i] - 1)];
      break;

    default:
      /*
       * sides 3, 4, 5, & 6 correspond to sides 1, 2, 3, & 4
       * of the quad element.
       */
      for (i = 0; i < 2; i++)
        ss_node_list[i] = connect[(quad_table[(side_num-2)][i] - 1)];
      break;
    }
    break;

  case SHELL8:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 8; i++)
        ss_node_list[i] = connect[(shell_table[side_num][i] - 1)];
      break;

    default:
      /*
       * sides 3, 4, 5, & 6 correspond to sides 1, 2, 3, & 4
       * of the quad element.
       */
      for (i = 0; i < 3; i++)
        ss_node_list[i] = connect[(quad_table[(side_num-2)][i] - 1)];
      break;
    }
    break;

  case TRI3:
    for (i = 0; i < 2; i++)
      ss_node_list[i] = connect[(tri_table[side_num][i] - 1)];
    break;

  case TRI6:
    for (i = 0; i < 3; i++)
      ss_node_list[i] = connect[(tri_table[side_num][i] - 1)];
    break;

  case TSHELL3:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 3; i++)
        ss_node_list[i] = connect[(tshell_table[side_num][i] - 1)];
      break;

    default:
      /*
       * sides 3, 4 & 5 correspond to sides 1, 2 & 3
       * of the tri element.
       */
      for (i = 0; i < 2; i++)
        ss_node_list[i] = connect[(tri_table[(side_num-2)][i] - 1)];
      break;
    }
    break;

  case TSHELL6:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 6; i++)
        ss_node_list[i] = connect[(tshell_table[side_num][i] - 1)];
      break;

    default:
      /*
       * sides 3, 4 & 5 correspond to sides 1, 2 & 3
       * of the tri element.
       */
      for (i = 0; i < 3; i++)
        ss_node_list[i] = connect[(tri_table[(side_num-2)][i] - 1)];
      break;
    }
    break;

  case HEX8:
    for (i = 0; i < 4; i++)
      ss_node_list[i] = connect[(hex_table[side_num][i] - 1)];
    break;

  case HEX20:
    for (i = 0; i < 8; i++)
      ss_node_list[i] = connect[(hex_table[side_num][i] - 1)];
    break;

  case HEX27:
    for (i = 0; i < 9; i++)
      ss_node_list[i] = connect[(hex_table[side_num][i] - 1)];
    break;

  case TET4:
    for (i = 0; i < 3; i++)
      ss_node_list[i] = connect[(tetra_table[side_num][i] - 1)];
    break;

  case TET10:
    for (i = 0; i < 6; i++)
      ss_node_list[i] = connect[(tetra_table[side_num][i] - 1)];
    break;

  case TET8:
    for (i = 0; i < 4; i++)
      ss_node_list[i] = connect[(tetra_table[side_num][i] - 1)];
    break;

  case WEDGE6:
    switch (side_num) {
    case 3:
    case 4:
      for (i = 0; i < 3; i++)
        ss_node_list[i] = connect[(wedge_table[side_num][i] - 1)];
      break;

    default:
      for (i = 0; i < 4; i++)
        ss_node_list[i] = connect[(wedge_table[side_num][i] - 1)];
      break;
    }
    break;

  case WEDGE15:
    switch (side_num) {
    case 3:
    case 4:
      for (i = 0; i < 6; i++)
        ss_node_list[i] = connect[(wedge_table[side_num][i] - 1)];
      break;

    default:
      for (i = 0; i < 8; i++)
        ss_node_list[i] = connect[(wedge_table[side_num][i] - 1)];
      break;
    }
    break;

  case WEDGE16:
    switch (side_num) {
    case 3:
    case 4:
      for (i = 0; i < 6; i++)
        ss_node_list[i] = connect[(wedge_table[side_num][i] - 1)];
      break;

    default:
      for (i = 0; i < 8; i++)
        ss_node_list[i] = connect[(wedge_table[side_num][i] - 1)];
      break;
    }
    break;

  case HEXSHELL:
    switch (side_num) {
    case 5:
    case 6:
      for (i = 0; i < 4; i++)
        ss_node_list[i] = connect[(hexshell_table[side_num][i] - 1)];
      break;

    default:
      for (i = 0; i < 6; i++)
        ss_node_list[i] = connect[(hexshell_table[side_num][i] - 1)];
      break;
    }
    break;

  case PYRAMID5:
    switch (side_num) {
    case 4:
      for (i = 0; i < 4; i++)
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      break;

    default:
      for (i = 0; i < 3; i++)
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      break;
    }
    break;

  case PYRAMID13:
    switch (side_num) {
    case 4:
      for (i = 0; i < 8; i++)
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      break;

    default:
      for (i = 0; i < 6; i++)
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      break;
    }
    break;

  case SPHERE: /* SHPERE's have no side sets */
  case NULL_EL:
    i = 0;
    break;

  } /* End "switch (etype)" */

  /* the variable "i" should be the number of positions that I filled */
  return (i);

} /*-------------------------End ss_to_node_list()---------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function get_ss_mirror() begins:
 *----------------------------------------------------------------------------
 * This function returns the node list for the mirror of the list
 * given. This will be the node list of a face that is connected
 * to this element on this face.
 *****************************************************************************/
int get_ss_mirror(const E_Type etype, const int *ss_node_list, int side_num,
                  int mirror_node_list[])

{
  int i=0;

  /*
   * the following arrays are the conversion from the side to
   * an opposing face
   */

  /* line (1-d) */
  static int line_table[3] = {1,0,2};

  /* square (2-d) */
  static int sqr_table[9] = {0,3,2,1,7,6,5,4,8};

  /* square hexshell (2-d) */
  static int hs_table[6] = {0,3,2,1,5,4};

  /* triangle (2-d) */
  static int tri_table[6] = {1,0,2,3,5,4};

/***************************** execution begins ******************************/

  /* Switch over the element type. */
  switch (etype) {
  case BAR2:
  case SHELL2:
    for (i = 0; i < 2; i++)
      mirror_node_list[i] = ss_node_list[line_table[i]];
    break;
  case BAR3:
  case SHELL3:
    for (i = 0; i < 3; i++)
      mirror_node_list[i] = ss_node_list[line_table[i]];
    break;
  case QUAD4:
    for (i = 0; i < 2; i++)
      mirror_node_list[i] = ss_node_list[line_table[i]];
    break;

  case QUAD8:
  case QUAD9:
    for (i = 0; i < 3; i++)
      mirror_node_list[i] = ss_node_list[line_table[i]];
    break;

  case SHELL4:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 4; i++)
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      break;

    default:
      for (i = 0; i < 2; i++)
        mirror_node_list[i] = ss_node_list[line_table[i]];
      break;
    }
    break;

  case SHELL8:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 8; i++)
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      break;

    default:
      for (i = 0; i < 3; i++)
        mirror_node_list[i] = ss_node_list[line_table[i]];
      break;
    }
    break;

  case TRI3:
    for (i = 0; i < 2; i++)
      mirror_node_list[i] = ss_node_list[line_table[i]];
    break;

  case TRI6:
    for (i = 0; i < 3; i++)
      mirror_node_list[i] = ss_node_list[line_table[i]];
    break;

  case TSHELL3:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 3; i++)
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      break;

    default:
      for (i = 0; i < 2; i++)
        mirror_node_list[i] = ss_node_list[line_table[i]];
      break;
    }
    break;

  case TSHELL6:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 6; i++)
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      break;

    default:
      for (i = 0; i < 3; i++)
        mirror_node_list[i] = ss_node_list[line_table[i]];
      break;
    }
    break;

  case HEX8:
    for (i = 0; i < 4; i++)
      mirror_node_list[i] = ss_node_list[sqr_table[i]];
    break;

  case HEX27:
    for (i = 0; i < 9; i++)
      mirror_node_list[i] = ss_node_list[sqr_table[i]];
    break;

  case HEX20:
    for (i = 0; i < 8; i++)
      mirror_node_list[i] = ss_node_list[sqr_table[i]];
    break;

  case TET4:
    for (i = 0; i < 3; i++)
      mirror_node_list[i] = ss_node_list[tri_table[i]];
    break;

  case TET8:
    for (i = 0; i < 4; i++)
      mirror_node_list[i] = ss_node_list[tri_table[i]];
    break;

  case TET10:
    for (i = 0; i < 6; i++)
      mirror_node_list[i] = ss_node_list[tri_table[i]];
    break;

  case WEDGE6:
    switch (side_num) {
    case 4:
    case 5:
      for (i = 0; i < 3; i++)
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      break;

    default:
      for (i = 0; i < 4; i++)
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      break;
    }
    break;

  case WEDGE15:
    switch (side_num) {
    case 4:
    case 5:
      for (i = 0; i < 6; i++)
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      break;

    default:
      for (i = 0; i < 8; i++)
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      break;
    }
    break;

  case WEDGE16:
    switch (side_num) {
    case 4:
    case 5:
      for (i = 0; i < 6; i++)
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      break;

    default:
      for (i = 0; i < 8; i++)
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      break;
    }
    break;

  case HEXSHELL:
    switch (side_num) {
    case 5:
    case 6:
      for (i = 0; i < 4; i++)
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      break;

    default:
      for (i = 0; i < 6; i++)
        mirror_node_list[i] = ss_node_list[hs_table[i]];
      break;
    }
    break;

  case PYRAMID5:
    switch (side_num) {
    case 5:
      for (i = 0; i < 4; i++)
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      break;

    default:
      for (i = 0; i < 3; i++)
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      break;
    }
    break;

  case PYRAMID13:
    switch (side_num) {
    case 5:
      for (i = 0; i < 8; i++)
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      break;

    default:
      for (i = 0; i < 6; i++)
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      break;
    }
    break;

  case SPHERE: /* SHPERE's have no side sets */
  case NULL_EL:
    i = 0;
    break;

  } /* End "switch (etype)" */

  /* the variable "i" should be the number of positions that I filled */
  return (i);

} /*-------------------------Ed get_ss_mirror()---------------------------*/

int numbermatch(int* sidenodes, int i, int j, int k, int value )
{
  if ( sidenodes[(i+j)%k] == value ) return 1;
  return 0;
}

/* #define numbermatch( sidenodes, i, j, k, value ) \
   ( sidenode[(i+j)%k] == value ) ? 1 : 0 */
