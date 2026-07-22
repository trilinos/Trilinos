// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER



#include "dr_const.h"
#include "dr_elem_const.h"
#include "dr_util_const.h"
#include "dr_err_const.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Author(s): Gary L. Hennigan (SNL 9221)
 *            Scott A. Hutchinson (SNL 9221)
 *            Matthew M. St. John (SNL 9226)
 *----------------------------------------------------------------------------
 * Functions contained in this file:
 *    get_elem_type()
 *    get_elem_info()
 *    get_side_id()
 *    ss_to_node_list()
 *    get_ss_mirror()
 *    get_elem_name()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

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
E_Type get_elem_type(char *elem_name, const int num_nodes,
                     const int num_dim)
{

  E_Type answer;
  int i;

  /* Convert elem_name to upper case. */
  for (i = strlen(elem_name); i>= 0; i--) {
    elem_name[i] = toupper(elem_name[i]);
  }

  if(strncmp(elem_name, "SPHERE", 6) == 0)
    answer = SPHERE;
  else if(strncmp(elem_name, "BEAM", 4) == 0 ||
          strncmp(elem_name, "TRUSS", 5) == 0 ||
          strncmp(elem_name, "BAR", 3) == 0)
  {
    switch(num_nodes)
    {
    case 2:
      answer=BAR1;
      break;
    case 3:
      answer=BAR2;
      break;
    default:
      Gen_Error(0, "fatal: unsupported BAR/BEAM/TRUSS element");
      answer = E_TYPE_ERROR;
    }
  }
  else if(strncmp(elem_name, "QUAD", 4) == 0)
  {
    switch(num_nodes)
    {
    case 4:
      answer = QUAD1;
      break;
    case 8:
      answer = S_QUAD2;
      break;
    case 9:
      answer = QUAD2;
      break;
    default:
      Gen_Error(0, "fatal: unsupported QUAD element");
      answer = E_TYPE_ERROR;
    }
  }
  else if(strncmp(elem_name, "HEX", 3) == 0)
  {
    switch(num_nodes)
    {
    case 8:
      answer = HEX1;
      break;
    case 12:
      answer = HEXSHELL;
      break;
    case 20:
      answer = S_HEX2;
      break;
    case 27:
      answer = HEX2;
      break;
    default:
      Gen_Error(0, "fatal: unsupported HEX element");
      answer = E_TYPE_ERROR;
    }
  }
  else if(strncmp(elem_name, "TRI", 3) == 0)
  {
    switch(num_nodes)
    {
    case 3:
      if (num_dim == 2) answer = TRI1;
      else              answer = TSHELL1;
      break;
    case 6:
      if (num_dim == 2) answer = TRI2;
      else              answer = TSHELL2;
      break;
    default:
      Gen_Error(0, "fatal: unsupported TRI element");
      answer = E_TYPE_ERROR;
    }
  }
  else if(strncmp(elem_name, "TET", 3) == 0)
  {
    switch(num_nodes)
    {
    case 4:
      answer = TET1;
      break;
    case 10:
      answer = TET2;
      break;
    default:
      Gen_Error(0, "fatal: unsupported TET element");
      answer = E_TYPE_ERROR;
    }
  }
  else if(strncmp(elem_name, "SHELL", 5) == 0)
  {
    switch(num_nodes)
    {
    case 4:
      answer = SHELL1;
      break;
    case 8:
      answer = SHELL2;
      break;
    default:
      Gen_Error(0, "fatal: unsupported SHELL element");
      answer = E_TYPE_ERROR;
    }
  }
  else if(strncmp(elem_name, "WEDGE", 5) == 0)
  {
    switch(num_nodes)
    {
    case 6:
      answer = WEDGE1;
      break;
    case 16:
      answer = WEDGE2;
      break;
    default:
      Gen_Error(0, "fatal: unsupported WEDGE element");
      answer = E_TYPE_ERROR;
    }
  }
  else
  {
    Gen_Error(0, "fatal: unknown element type read");
    answer = E_TYPE_ERROR;
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
int get_elem_info(const int req, const E_Type etype, const ZOLTAN_ID_TYPE sid)
{

  int answer;

  switch(etype)		/* Switch over the element type */
  {
  case BAR1:
    switch(req)
    {
    case NNODES:
      answer = 2;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 1;
      break;
    case NSIDES:
      answer = 0;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case BAR2:
    switch(req)
    {
    case NNODES:
      answer = 3;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 1;
      break;
    case NSIDES:
      answer = 0;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case SPHERE:
    switch(req)
    {
    case NNODES:
      answer = 1;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 0;
      break;
    case NSIDES:
      answer = 0;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD:
    case NQUAD_SURF:
    default:
      Gen_Error(0, "fatal:unknown quantity");
      answer = -1;
    }
    break;

  case QUAD1:		/* First order quad */
    switch(req)		/* select type of information required*/
    {
    case NNODES:	/* number of nodes */
      answer = 4;
      break;
    case NQUAD:		/* number of quadrature points */
      answer = 4;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:	/* number of surface quad points */
      answer = 2;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 2;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal:unknown quantity");
      answer = -1;
    }
    break;

  case S_QUAD2:		/* 2nd order serendipity quad */
    switch(req)		/* select type of information required */
    {
    case NNODES:		/* number of nodes */
      answer = 8;
      break;
    case NQUAD:		/* number of quadrature points */
      answer = 9;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:	/* number of surface quad points */
      answer = 3;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 3;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case QUAD2:	/* biquadratic quadrilateral */
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 9;
      break;
    case NQUAD:		/* number of quadrature points */
      answer = 9;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:	/* number of surface quad points */
      answer = 3;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 3;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case SHELL1:
    switch(req)
    {
    case NNODES:
      answer = 4;
      break;
    case NSIDES:
      answer = 6;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSNODES:	/* number of side nodes */
      switch (sid)
      {
      case 1:
      case 2:
        answer = 4;
        break;
      case 3:
      case 4:
      case 5:
      case 6:
        answer = 2;
        break;
      default:
        Gen_Error(0, "fatal: unknown side id");
        answer = -1;
      }
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case SHELL2:
    switch(req)
    {
    case NNODES:
      answer = 8;
      break;
    case NSIDES:
      answer = 6;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSNODES:	/* number of side nodes */
      switch (sid)
      {
      case 1:
      case 2:
        answer = 8;
        break;
      case 3:
      case 4:
      case 5:
      case 6:
        answer = 3;
        break;
      default:
        Gen_Error(0, "fatal: unknown side id");
        answer = -1;
      }
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case TRI1:
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 3;
      break;
    case NQUAD:		/* number of quadrature points */
      answer = 4;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:	/* number of surface quad points */
      answer = 2;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 2;
      break;
    case NSIDES:
      answer = 3;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case TRI2:
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 6;
      break;
    case NQUAD:		/* number of quadrature points */
      answer = 7;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 2;
      break;
    case NQUAD_SURF:	/* number of surface quad points */
      answer = 3;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 3;
      break;
    case NSIDES:
      answer = 3;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case TSHELL1:
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 3;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSIDES:
      answer = 5;
      break;
    case NSNODES:	/* number of side nodes */
      switch (sid)
      {
      case 1:
      case 2:
        answer = 3;
        break;
      case 3:
      case 4:
      case 5:
        answer = 2;
        break;
      default:
        Gen_Error(0, "fatal: unknown side id");
        answer = -1;
      }
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case TSHELL2:
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 6;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSIDES:
      answer = 5;
      break;
    case NSNODES:	/* number of side nodes */
      switch (sid)
      {
      case 1:
      case 2:
        answer = 6;
        break;
      case 3:
      case 4:
      case 5:
        answer = 3;
        break;
      default:
        Gen_Error(0, "fatal: unknown side id");
        answer = -1;
      }
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case HEX1:		/* trilinear hexahedron */
    switch(req)	/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 8;
      break;
    case NQUAD:		/* number of quadrature points */
      answer = 8;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:	/* number of surface quad points */
      answer = 4;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 4;
      break;
    case NSIDES:
      answer = 6;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case S_HEX2:		/* serendipity triquadratic hexahedron */
    switch(req)	/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 20;
      break;
    case NQUAD:		/* number of quadrature points */
      answer = 27;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:	/* number of surface quad points */
      answer = 9;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 8;
      break;
    case NSIDES:
      answer = 6;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case HEX2:		/* triquadratic hexahedron */
    switch(req)		/* select type of information required*/
    {
    case NNODES:	/* number of nodes */
      answer = 27;
      break;
    case NQUAD:		/* number of quadrature points */
      answer = 27;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:	/* number of surface quad points */
      answer = 9;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 9;
      break;
    case NSIDES:
      answer = 6;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case TET1:		/* trilinear tetrahedron */
    switch(req)		/* select type of information required*/
    {
    case NNODES:	/* number of nodes */
      answer = 4;
      break;
    case NQUAD:		/* number of quadrature points */
      answer = 5;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:	/* number of surface quad points */
      answer = 4;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 3;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case TET2:		/* triquadradic tetrahedron */
    switch(req)		/* select type of information required */
    {
    case NNODES:	/* number of nodes */
      answer = 10;
      break;
    case NQUAD:		/* number of quadrature points */
      answer = 15;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NQUAD_SURF:	/* number of surface quad points */
      answer = 7;
      break;
    case NSNODES:	/* number of side nodes */
      answer = 6;
      break;
    case NSIDES:
      answer = 4;
      break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case WEDGE1:
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
    case NSNODES:	/* number of side nodes */
      switch (sid)
      {
      case 1:
      case 2:
      case 3:
        answer = 4;
        break;
      case 4:
      case 5:
        answer = 3;
        break;
      default:
        Gen_Error(0, "fatal: unknown side id");
        answer = -1;
      }
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case WEDGE2:
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
    case NSNODES:	/* number of side nodes */
      switch (sid)
      {
      case 1:
      case 2:
      case 3:
        answer = 8;
        break;
      case 4:
      case 5:
        answer = 6;
        break;
      default:
        Gen_Error(0, "fatal: unknown side id");
        answer = -1;
      }
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;

  case HEXSHELL:
    switch(req)
    {
    case NNODES:
      answer = 12;
      break;
    case NSIDES:
      answer = 6;
      break;
    case NDIM:		/* number of physical dimensions */
      answer = 3;
      break;
    case NSNODES:	/* number of side nodes */
      switch (sid)
      {
      case 1:
      case 2:
      case 3:
      case 4:
        answer = 6;
        break;
      case 5:
      case 6:
        answer = 4;
        break;
      default:
        Gen_Error(0, "fatal: unknown side id");
        answer = -1;
      }
    default:
      Gen_Error(0, "fatal: unknown quantity");
      answer = -1;
    }
    break;


  default:
    Gen_Error(0, "fatal: unknown or unimplemented element type");
    answer = -1;
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
int get_side_id(E_Type etype, const ZOLTAN_ID_TYPE *connect, const int nsnodes,
                int side_nodes[])
{
  const char *func_name="get_side_id";
  char  err_buff[300];

  int nnodes, i, j, num;
  int dup, location = 0;

  /* check if this is a degenerate face */
  dup = 0;
  for (i = 0; i < (nsnodes - 1); i++) {
    for (j = (i + 1); j < nsnodes; j++)
      if (side_nodes[i] == side_nodes[j]) {
        dup = 1;
        location = i; /* location of duplicated node */
        break;
      }
    if (dup) break;
  }

  nnodes = get_elem_info(NNODES, etype, 0);

  /* Find all of the side nodes in the connect table */
  num = 0;
  for(i=0; i < nnodes; i++)
  {
    for(j=0; j < nsnodes; j++)
    {
      if((int)connect[i] == side_nodes[j])
      {
        num++;
        break;
      }
    }
    if(num == nsnodes)
      break;
  }

  if(num != nsnodes)
  {
    Gen_Error(0, "fatal: not all side nodes in connect table for element");
    return -1;
  }

  /* Find the side ID */
  switch(etype)
  {
  case QUAD1:
  case S_QUAD2:
  case QUAD2:
    /* SIDE 1 */
    if (side_nodes[0] == (int)connect[0] &&
        side_nodes[1] == (int)connect[1]) return 1;

    /* SIDE 2 */
    if (side_nodes[0] == (int)connect[1] &&
        side_nodes[1] == (int)connect[2]) return 2;

    /* SIDE 3 */
    if (side_nodes[0] == (int)connect[2] &&
        side_nodes[1] == (int)connect[3]) return 3;

    /* SIDE 4 */
    if (side_nodes[0] == (int)connect[3] &&
        side_nodes[1] == (int)connect[0]) return 4;

    break;

  case TRI1:
  case TRI2:
    /* SIDE 1 */
    if (side_nodes[0] == (int)connect[0] &&
        side_nodes[1] == (int)connect[1]) return 1;

    /* SIDE 2 */
    if (side_nodes[0] == (int)connect[1] &&
        side_nodes[1] == (int)connect[2]) return 2;

    /* SIDE 3 */
    if (side_nodes[0] == (int)connect[2] &&
        side_nodes[1] == (int)connect[0]) return 3;

    break;

  case TET1:
  case TET2:
    /* SIDE 1 */
    if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 3] == (int)connect[1] &&
          side_nodes[(2 + num) % 3] == (int)connect[3]) return 1;
    }

    /* SIDE 2 */
    if((num = in_list2((int)connect[1], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 3] == (int)connect[2] &&
          side_nodes[(2 + num) % 3] == (int)connect[3]) return 2;
    }

    /* SIDE 3 */
    if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 3] == (int)connect[3] &&
          side_nodes[(2 + num) % 3] == (int)connect[2]) return 3;
    }

    /* SIDE 4 */
    if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 3] == (int)connect[2] &&
          side_nodes[(2 + num) % 3] == (int)connect[1]) return 4;
    }

    break;

  case HEX1:
  case S_HEX2:
  case HEX2:
  case HEXSHELL:  /* this should be the same as a HEX element */
    /* SIDE 1 */
    if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 4] == (int)connect[1] &&
          side_nodes[(2 + num) % 4] == (int)connect[5] &&
          side_nodes[(3 + num) % 4] == (int)connect[4]) return 1;

      /* if this is the duplicated node, then find the next occurence */
      if (dup && (int)connect[0] == side_nodes[location]) {
        num = in_list2((int)connect[0], (nsnodes-num), &(side_nodes[num+1]));
        if (side_nodes[(1 + num) % 4] == (int)connect[1] &&
            side_nodes[(2 + num) % 4] == (int)connect[5] &&
            side_nodes[(3 + num) % 4] == (int)connect[4]) return 1;
      }
    }

    /* SIDE 2 */
    if((num = in_list2((int)connect[1], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 4] == (int)connect[2] &&
          side_nodes[(2 + num) % 4] == (int)connect[6] &&
          side_nodes[(3 + num) % 4] == (int)connect[5]) return 2;

      /* if this is the duplicated node, then find the next occurence */
      if (dup && (int)connect[1] == side_nodes[location]) {
        num = in_list2((int)connect[1], (nsnodes-num), &(side_nodes[num+1]));
        if (side_nodes[(1 + num) % 4] == (int)connect[2] &&
            side_nodes[(2 + num) % 4] == (int)connect[6] &&
            side_nodes[(3 + num) % 4] == (int)connect[5]) return 2;
      }
    }

    /* SIDE 3 */
    if((num = in_list2((int)connect[2], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 4] == (int)connect[3] &&
          side_nodes[(2 + num) % 4] == (int)connect[7] &&
          side_nodes[(3 + num) % 4] == (int)connect[6]) return 3;

      /* if this is the duplicated node, then find the next occurence */
      if (dup && (int)connect[2] == side_nodes[location]) {
        num = in_list2((int)connect[2], (nsnodes-num), &(side_nodes[num+1]));
        if (side_nodes[(1 + num) % 4] == (int)connect[3] &&
            side_nodes[(2 + num) % 4] == (int)connect[7] &&
            side_nodes[(3 + num) % 4] == (int)connect[6]) return 3;
      }
    }

    /* SIDE 4 */
    if((num = in_list2((int)connect[3], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 4] == (int)connect[0] &&
          side_nodes[(2 + num) % 4] == (int)connect[4] &&
          side_nodes[(3 + num) % 4] == (int)connect[7]) return 4;

      /* if this is the duplicated node, then find the next occurence */
      if (dup && (int)connect[3] == side_nodes[location]) {
        num = in_list2((int)connect[3], (nsnodes-num), &(side_nodes[num+1]));
        if (side_nodes[(1 + num) % 4] == (int)connect[0] &&
            side_nodes[(2 + num) % 4] == (int)connect[4] &&
            side_nodes[(3 + num) % 4] == (int)connect[7]) return 4;
      }
    }

    /* SIDE 5 */
    if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 4] == (int)connect[3] &&
          side_nodes[(2 + num) % 4] == (int)connect[2] &&
          side_nodes[(3 + num) % 4] == (int)connect[1]) return 5;

      /* if this is the duplicated node, then find the next occurence */
      if (dup && (int)connect[0] == side_nodes[location]) {
        num = in_list2((int)connect[0], (nsnodes-num), &(side_nodes[num+1]));
        if (side_nodes[(1 + num) % 4] == (int)connect[3] &&
            side_nodes[(2 + num) % 4] == (int)connect[2] &&
            side_nodes[(3 + num) % 4] == (int)connect[1]) return 5;
      }
    }

    /* SIDE 6 */
    if((num = in_list2((int)connect[4], nsnodes, side_nodes)) >= 0) {
      if (side_nodes[(1 + num) % 4] == (int)connect[5] &&
          side_nodes[(2 + num) % 4] == (int)connect[6] &&
          side_nodes[(3 + num) % 4] == (int)connect[7]) return 6;

      /* if this is the duplicated node, then find the next occurence */
      if (dup && (int)connect[4] == side_nodes[location]) {
        num = in_list2((int)connect[4], (nsnodes-num), &(side_nodes[num+1]));
        if (side_nodes[(1 + num) % 4] == (int)connect[5] &&
            side_nodes[(2 + num) % 4] == (int)connect[6] &&
            side_nodes[(3 + num) % 4] == (int)connect[7]) return 6;
      }
    }

    break;

  case SHELL1:
  case SHELL2:

    /* 2D sides */
    if(nsnodes == 2 || nsnodes == 3) {
      /* SIDE 3 */
      if (side_nodes[0] == (int)connect[0] &&
          side_nodes[1] == (int)connect[1]) return 3;

      /* SIDE 4 */
      if (side_nodes[0] == (int)connect[1] &&
          side_nodes[1] == (int)connect[2]) return 4;

      /* SIDE 5 */
      if (side_nodes[0] == (int)connect[2] &&
          side_nodes[1] == (int)connect[3]) return 5;

      /* SIDE 6 */
      if (side_nodes[0] == (int)connect[3] &&
          side_nodes[1] == (int)connect[0]) return 6;
    }

    /* 3D faces */
    else if (nsnodes == 4 || nsnodes == 8) {

      /* SIDE 1 */
      if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == (int)connect[1] &&
            side_nodes[(2 + num) % 4] == (int)connect[2] &&
            side_nodes[(3 + num) % 4] == (int)connect[3]) return 1;
      }

      /* SIDE 2 */
      if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == (int)connect[3] &&
            side_nodes[(2 + num) % 4] == (int)connect[2] &&
            side_nodes[(3 + num) % 4] == (int)connect[1]) return 2;
      }
    }

    break;

  case WEDGE1:
  case WEDGE2:

    /* quad sides */
    if (nsnodes == 4 || nsnodes == 8) {
      /* SIDE 1 */
      if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == (int)connect[1] &&
            side_nodes[(2 + num) % 4] == (int)connect[4] &&
            side_nodes[(3 + num) % 4] == (int)connect[3]) return 1;
      }

      /* SIDE 2 */
      if((num = in_list2((int)connect[1], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == (int)connect[2] &&
            side_nodes[(2 + num) % 4] == (int)connect[5] &&
            side_nodes[(3 + num) % 4] == (int)connect[4]) return 2;
      }

      /* SIDE 3 */
      if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 4] == (int)connect[3] &&
            side_nodes[(2 + num) % 4] == (int)connect[5] &&
            side_nodes[(3 + num) % 4] == (int)connect[2]) return 3;
      }
    }

    /* triangle sides */
    else if (nsnodes == 3 || nsnodes == 6) {
      /* SIDE 4 */
      if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == (int)connect[2] &&
            side_nodes[(3 + num) % 3] == (int)connect[1]) return 4;
      }

      /* SIDE 5 */
      if((num = in_list2((int)connect[3], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == (int)connect[4] &&
            side_nodes[(3 + num) % 3] == (int)connect[5]) return 5;
      }
    }

    break;

  case TSHELL1:
  case TSHELL2:

    /* 2D sides */
    if(nsnodes == 2 || ((etype == TSHELL2) && (nsnodes == 3))) {
      /* SIDE 3 */
      if (side_nodes[0] == (int)connect[0] &&
          side_nodes[1] == (int)connect[1]) return 3;

      /* SIDE 4 */
      if (side_nodes[0] == (int)connect[1] &&
          side_nodes[1] == (int)connect[2]) return 4;

      /* SIDE 5 */
      if (side_nodes[0] == (int)connect[2] &&
          side_nodes[1] == (int)connect[0]) return 5;

    }

    /* 3D faces */
    else if (nsnodes == 3 || nsnodes == 6) {

      /* SIDE 1 */
      if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == (int)connect[1] &&
            side_nodes[(2 + num) % 3] == (int)connect[2]) return 1;
      }

      /* SIDE 2 */
      if((num = in_list2((int)connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == (int)connect[2] &&
            side_nodes[(2 + num) % 3] == (int)connect[1]) return 2;
      }
    }

    break;


  default:
    sprintf(err_buff, "fatal: unknown element type %d in function %s",
            etype, func_name);
    Gen_Error(0, err_buff);
    return -1;

  } /* End "switch(etype)" */

  return 0;
} /*---------------------------End get_side_id()-----------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function ss_to_node_list() begins:
 *----------------------------------------------------------------------------
 * This function returns the list of nodes in a side of an element given
 * the element type, and the side id. It also returns the number of nodes
 * in that side.
 *****************************************************************************/
int ss_to_node_list(const E_Type etype, const ZOLTAN_ID_TYPE *connect, int side_num,
                    int ss_node_list[])

{
  ZOLTAN_ID_TYPE i;
  char msg[80];

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
  {1,4,3,2,12,11,10,9,22},  {5,6,7,8,17,18,19,20,23}                  /*nodes*/
 };

  /* hexshell */
  static int hexshell_table[6][6] = {
  /*      1               2                3                4         side   */
    {1,2,6,5,10,9}, {2,3,7,6,11,10}, {3,4,8,7,12,11}, {4,1,5,8,9,12},
  /*      5               6                                           side   */
    {1,4,3,2,0,0},  {5,6,7,8,0,0}                                     /*nodes*/
 };

/***************************** execution begins ******************************/

  /* Locally decrement side_num */

  side_num--;

  /* Switch over the element type. */
  switch (etype) {
  case QUAD1:
    for (i = 0; i < 2; i++)
      ss_node_list[i] = (int)connect[(quad_table[side_num][i] - 1)];
    break;

  case S_QUAD2:
  case QUAD2:
    for (i = 0; i < 3; i++)
      ss_node_list[i] = (int)connect[(quad_table[side_num][i] - 1)];
    break;

  case SHELL1:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 4; i++)
        ss_node_list[i] = (int)connect[(shell_table[side_num][i] - 1)];
      break;

    default:
      /*
       * sides 3, 4, 5, & 6 correspond to sides 1, 2, 3, & 4
       * of the quad element.
       */
      for (i = 0; i < 2; i++)
        ss_node_list[i] = (int)connect[(quad_table[(side_num-2)][i] - 1)];
      break;
    }
    break;

  case SHELL2:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 8; i++)
        ss_node_list[i] = (int)connect[(shell_table[side_num][i] - 1)];
      break;

    default:
      /*
       * sides 3, 4, 5, & 6 correspond to sides 1, 2, 3, & 4
       * of the quad element.
       */
      for (i = 0; i < 3; i++)
        ss_node_list[i] = (int)connect[(quad_table[(side_num-2)][i] - 1)];
      break;
    }
    break;

  case TRI1:
    for (i = 0; i < 2; i++)
      ss_node_list[i] = (int)connect[(tri_table[side_num][i] - 1)];
    break;

  case TRI2:
    for (i = 0; i < 3; i++)
      ss_node_list[i] = (int)connect[(tri_table[side_num][i] - 1)];
    break;

  case TSHELL1:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 3; i++)
        ss_node_list[i] = (int)connect[(tshell_table[side_num][i] - 1)];
      break;

    default:
      /*
       * sides 3, 4 & 5 correspond to sides 1, 2 & 3
       * of the tri element.
       */
      for (i = 0; i < 2; i++)
        ss_node_list[i] = (int)connect[(tri_table[(side_num-2)][i] - 1)];
      break;
    }
    break;

  case TSHELL2:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 6; i++)
        ss_node_list[i] = (int)connect[(tshell_table[side_num][i] - 1)];
      break;

    default:
      /*
       * sides 3, 4 & 5 correspond to sides 1, 2 & 3
       * of the tri element.
       */
      for (i = 0; i < 3; i++)
        ss_node_list[i] = (int)connect[(tri_table[(side_num-2)][i] - 1)];
      break;
    }
    break;

  case HEX1:
    for (i = 0; i < 4; i++)
      ss_node_list[i] = (int)connect[(hex_table[side_num][i] - 1)];
    break;

  case S_HEX2:
    for (i = 0; i < 8; i++)
      ss_node_list[i] = (int)connect[(hex_table[side_num][i] - 1)];
    break;

  case HEX2:
    for (i = 0; i < 9; i++)
      ss_node_list[i] = (int)connect[(hex_table[side_num][i] - 1)];
    break;

  case TET1:
    for (i = 0; i < 3; i++)
      ss_node_list[i] = (int)connect[(tetra_table[side_num][i] - 1)];
    break;

  case TET2:
    for (i = 0; i < 6; i++)
      ss_node_list[i] = (int)connect[(tetra_table[side_num][i] - 1)];
    break;

  case WEDGE1:
    switch (side_num) {
    case 3:
    case 4:
      for (i = 0; i < 3; i++)
        ss_node_list[i] = (int)connect[(wedge_table[side_num][i] - 1)];
      break;

    default:
      for (i = 0; i < 4; i++)
        ss_node_list[i] = (int)connect[(wedge_table[side_num][i] - 1)];
      break;
    }
    break;

  case WEDGE2:
    switch (side_num) {
    case 3:
    case 4:
      for (i = 0; i < 6; i++)
        ss_node_list[i] = (int)connect[(wedge_table[side_num][i] - 1)];
      break;

    default:
      for (i = 0; i < 8; i++)
        ss_node_list[i] = (int)connect[(wedge_table[side_num][i] - 1)];
      break;
    }
    break;

  case HEXSHELL:
    switch (side_num) {
    case 5:
    case 6:
      for (i = 0; i < 4; i++)
        ss_node_list[i] = (int)connect[(hexshell_table[side_num][i] - 1)];
      break;

    default:
      for (i = 0; i < 6; i++)
        ss_node_list[i] = (int)connect[(hexshell_table[side_num][i] - 1)];
      break;
    }
    break;

  case SPHERE: /* SHPERE's have no side sets */
    i = 0;
    break;

  default:  /* Yowsa!  More element types? */
    sprintf(msg, "Warning: ss_to_node_list not supported for element type %d",
	    etype);
    Gen_Error(0, msg);
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
ZOLTAN_ID_TYPE get_ss_mirror(const E_Type etype, const ZOLTAN_ID_TYPE *ss_node_list, ZOLTAN_ID_TYPE side_num,
                  ZOLTAN_ID_TYPE mirror_node_list[])

{
  ZOLTAN_ID_TYPE   i;
  char msg[80];

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
  case QUAD1:
    for (i = 0; i < 2; i++)
      mirror_node_list[i] = ss_node_list[line_table[i]];
    break;

  case S_QUAD2:
  case QUAD2:
    for (i = 0; i < 3; i++)
      mirror_node_list[i] = ss_node_list[line_table[i]];
    break;

  case SHELL1:
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

  case SHELL2:
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

  case TRI1:
    for (i = 0; i < 2; i++)
      mirror_node_list[i] = ss_node_list[line_table[i]];
    break;

  case TRI2:
    for (i = 0; i < 3; i++)
      mirror_node_list[i] = ss_node_list[line_table[i]];
    break;

  case TSHELL1:
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

  case TSHELL2:
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

  case HEX1:
    for (i = 0; i < 4; i++)
      mirror_node_list[i] = ss_node_list[sqr_table[i]];
    break;

  case HEX2:
    for (i = 0; i < 9; i++)
      mirror_node_list[i] = ss_node_list[sqr_table[i]];
    break;

  case S_HEX2:
    for (i = 0; i < 8; i++)
      mirror_node_list[i] = ss_node_list[sqr_table[i]];
    break;

  case TET1:
    for (i = 0; i < 3; i++)
      mirror_node_list[i] = ss_node_list[tri_table[i]];
    break;

  case TET2:
    for (i = 0; i < 6; i++)
      mirror_node_list[i] = ss_node_list[tri_table[i]];
    break;

  case WEDGE1:
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

  case WEDGE2:
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

  case SPHERE: /* SHPERE's have no side sets */
    i = 0;
    break;

  default:  /* Yowsa!  More element types? */
    sprintf(msg, "Warning: get_ss_mirror not supported for element type %d",
	    etype);
    Gen_Error(0, msg);
    i = 0;
    break;
  } /* End "switch (etype)" */

  /* the variable "i" should be the number of positions that I filled */
  return (i);

} /*-------------------------Ed get_ss_mirror()---------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

const char *get_elem_name(int itype) {
/* Function to return the name of an element given its type.
 * Inverse of get_elem_type().
 */
E_Type etype = (E_Type) itype;

  switch (etype) {
  case SPHERE:
    return("SPHERE");
  case BAR1:
  case BAR2:
    return("BAR");
  case QUAD1:
  case S_QUAD2:
  case QUAD2:
    return("QUAD");
  case HEX1:
  case HEXSHELL:
  case S_HEX2:
  case HEX2:
    return("HEX");
  case TRI1:
  case TSHELL1:
  case TRI2:
  case TSHELL2:
    return("TRI");
  case TET1:
  case TET2:
    return("TET");
  case SHELL1:
  case SHELL2:
    return("SHELL");
  case WEDGE1:
  case WEDGE2:
    return("WEDGE");
  default:
    Gen_Error(0, "fatal: unknown element type read");
    return("NULL");
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
