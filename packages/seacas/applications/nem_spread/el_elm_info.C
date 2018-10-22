/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
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
 *     * Neither the name of NTESS nor the names of its
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

#include "el_elm.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define LINEAR 1
#define QUADRATIC 2

/*************** R O U T I N E S   I N   T H I S   F I L E ********************
 *
 *  NAME				TYPE		CALL_BY
 * ---------------		-------		------------------------
 *  elem_info ()			int		"rf_fill.c" matrix_fill
 *  get_type  ()			int
 *  calc_elem_vol()              double          multiple routines
 *
 ******************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int elem_info(int info, int ielem_type, int supp)

/*
 * Function which returns the various parameters * for the elements, e.g.,
 * polynomial order, number of * Gaussian-quadrature points, etc, based upon
 * what a code * passed from the calling routine.
 *
 * Author:          Scott Hutchinson (1421)
 * Date:            15 May 1992
 *
 * The routine currently handles the following requests for information
 * about element type, ielem_type:
 *
 *	NNODES     = Number of nodes in the element
 *	NDIM	   = Dimension of the element
 *      NN_SIDE    = Number of nodes on a side/face of the element.
 */

{

  int         answer = 0;
  const char *yo     = "elem_info: ";

  /* return desired element information */

  switch (ielem_type) { /* select type of element */

  case SPHERE:
    switch (info) {
    case NNODES: answer = 1; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case SHELL4:
    switch (info) {
    case NNODES: answer = 4; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 1:
      case 2: answer = 4; break;
      default: answer = 2; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case SHELL8:
    switch (info) {
    case NNODES: answer = 8; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 1:
      case 2: answer = 8; break;
      default: answer = 3; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case SHELL9:
    switch (info) {
    case NNODES: answer = 9; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 1:
      case 2: answer = 9; break;
      default: answer = 3; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case TSHELL3:
    switch (info) {
    case NNODES: answer = 3; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 1:
      case 2: answer = 3; break;
      default: answer = 2; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case TSHELL4:
    switch (info) {
    case NNODES: answer = 4; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 1:
      case 2: answer = 4; break;
      default: answer = 2; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case TSHELL6:
    switch (info) {
    case NNODES: answer = 6; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 1:
      case 2: answer = 6; break;
      default: answer = 3; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case TSHELL7:
    switch (info) {
    case NNODES: answer = 7; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 1:
      case 2: answer = 7; break;
      default: answer = 3; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case BAR2:
  case SHELL2:
    switch (info) { /* select type of information required */
    case NNODES: answer = 2; break;
    case NDIM: answer = 1; break;
    case NN_SIDE: answer = 1; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case BAR3:
  case SHELL3:
    switch (info) { /* select type of information required */
    case NNODES: answer = 3; break;
    case NDIM: answer = 1; break;
    case NN_SIDE: answer = 1; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case QUAD4:       /* bilinear quadrilateral */
    switch (info) { /* select type of information required */
    case NNODES: answer = 4; break;
    case NDIM: answer = 2; break;
    case NN_SIDE: answer = 2; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case QUAD8:       /* biquadratic serendipity quadrilateral */
    switch (info) { /* select type of information required */
    case NNODES: answer = 8; break;
    case NDIM: answer = 2; break;
    case NN_SIDE: answer = 3; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case QUAD9:       /* biquadratic quadrilateral */
    switch (info) { /* select type of information required */
    case NNODES: answer = 9; break;
    case NDIM: answer = 2; break;
    case NN_SIDE: answer = 3; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case TRI3:
    switch (info) { /* select type of information required */
    case NNODES: answer = 3; break;
    case NDIM: answer = 2; break;
    case NN_SIDE: answer = 2; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case TRI4:
    switch (info) { /* select type of information required */
    case NNODES: answer = 4; break;
    case NDIM: answer = 2; break;
    case NN_SIDE: answer = 2; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case TRI6:
    switch (info) { /* select type of information required */
    case NNODES: answer = 6; break;
    case NDIM: answer = 2; break;
    case NN_SIDE: answer = 3; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case TRI7:
    switch (info) { /* select type of information required */
    case NNODES: answer = 7; break;
    case NDIM: answer = 2; break;
    case NN_SIDE: answer = 3; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case HEX8:        /* trilinear hexahedron */
    switch (info) { /* select type of information required */
    case NNODES: answer = 8; break;
    case NDIM: answer = 3; break;
    case NN_SIDE: answer = 4; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case HEX16: /* localization element */
    switch (info) {
    case NNODES: answer = 16; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 5:
      case 6: answer = 8; break;
      default: answer = 6; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case HEX20:       /* serendipity triquadratic hexahedron */
    switch (info) { /* select type of information required */
    case NNODES: answer = 20; break;
    case NDIM: answer = 3; break;
    case NN_SIDE: answer = 8; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case HEX27:       /* triquadratic hexahedron */
    switch (info) { /* select type of information required */
    case NNODES: answer = 27; break;
    case NDIM: answer = 3; break;
    case NN_SIDE: answer = 9; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity.\n", yo); exit(1);
    }
    break;

  case TET4:        /* trilinear tetrahedron */
    switch (info) { /* select type of information required */
    case NNODES: answer = 4; break;
    case NDIM: answer = 3; break;
    case NN_SIDE: answer = 3; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case TET10:       /* triquadradic tetrahedron */
    switch (info) { /* select type of information required */
    case NNODES: answer = 10; break;
    case NDIM: answer = 3; break;
    case NN_SIDE: answer = 6; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity.\n", yo); exit(1);
    }
    break;

  case TET14:       /* triquadradic tetrahedron */
    switch (info) { /* select type of information required */
    case NNODES: answer = 14; break;
    case NDIM: answer = 3; break;
    case NN_SIDE: answer = 7; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity.\n", yo); exit(1);
    }
    break;

  case TET15:       /* triquadradic tetrahedron */
    switch (info) { /* select type of information required */
    case NNODES: answer = 15; break;
    case NDIM: answer = 3; break;
    case NN_SIDE: answer = 7; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity.\n", yo); exit(1);
    }
    break;

  case TET8:        /* triquadradic tetrahedron */
    switch (info) { /* select type of information required */
    case NNODES: answer = 8; break;
    case NDIM: answer = 3; break;
    case NN_SIDE: answer = 4; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity.\n", yo); exit(1);
    }
    break;

  case WEDGE6:
    switch (info) {
    case NNODES: answer = 6; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 4:
      case 5: answer = 3; break;
      default: answer = 4; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case WEDGE12:
    switch (info) {
    case NNODES: answer = 12; break;
    case NDIM: answer = 3; break;
    case NN_SIDE: answer = 6; break;
    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case WEDGE16:
    switch (info) {
    case NNODES: answer = 16; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 4:
      case 5: answer = 6; break;
      default: answer = 8; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case WEDGE15:
    switch (info) {
    case NNODES: answer = 15; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 4:
      case 5: answer = 6; break;

      default: answer = 8; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case WEDGE20:
    switch (info) {
    case NNODES: answer = 20; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 4:
      case 5: answer = 7; break;

      default: answer = 9; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case WEDGE21:
    switch (info) {
    case NNODES: answer = 21; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 4:
      case 5: answer = 7; break;

      default: answer = 9; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case HEXSHELL:
    switch (info) {
    case NNODES: answer = 12; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 5:
      case 6: answer = 4; break;
      default: answer = 6; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case PYRAMID5:
    switch (info) {
    case NNODES: answer = 5; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 5: answer = 4; break;
      default: answer = 3; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case PYRAMID13:
    switch (info) {
    case NNODES: answer = 13; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 5: answer = 8; break;
      default: answer = 6; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case PYRAMID14:
    switch (info) {
    case NNODES: answer = 14; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 5: answer = 9; break;
      default: answer = 6; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case PYRAMID18:
    switch (info) {
    case NNODES: answer = 18; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 5: answer = 9; break;
      default: answer = 7; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  case PYRAMID19:
    switch (info) {
    case NNODES: answer = 19; break;
    case NDIM: answer = 3; break;
    case NN_SIDE:
      switch (supp) {
      case 5: answer = 9; break;
      default: answer = 7; break;
      }
      break;

    default: fprintf(stderr, "%sERROR: Unknown quantity\n", yo); exit(1);
    }
    break;

  default: fprintf(stderr, "%sERROR: Unimplemented element type.\n", yo); exit(1);
  }

  return answer;

} /* elem_info */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int get_type(char string[], int nodes, int num_dim)

/*
 * Function which returns the element type according to this analysis code
 * based on the EXODUS element type string and the number of nodes in the
 * element.
 *
 * Author:          Scott Hutchinson (1421)
 * Date:            29 June 1992
 */

{

  int         answer = 0;
  const char *yo     = "get_type: ";

  /* Precondition: string is lower case */
  switch (string[0]) {
  case 'q':

    if (strncmp(string, "quad", 4) == 0) { /* select element shape */
      switch (nodes) {                     /* select number of nodes in this element */
      case 4: /* bilinear quadralateral */ answer = QUAD4; break;
      case 8: /* serendipity biquadratic quadralateral */ answer = QUAD8; break;
      case 9: /* biquadratic quadrilateral */ answer = QUAD9; break;
      default:
        fprintf(stderr,
                "%sERROR: Quadralateral element with %d nodes "
                "not valid.\n",
                yo, nodes);
        exit(1);
      }
    }
    break;

  case 's':
    if (strncmp(string, "sphere", 6) == 0) {
      answer = SPHERE;
    }
    else if (strncmp(string, "shell", 5) == 0) {
      switch (nodes) {
      case 2:
        if (num_dim == 2) {
          answer = SHELL2;
        }
        else {
          fprintf(stderr,
                  "%sERROR: Shell element with %d nodes "
                  "only valid in 2D.\n",
                  yo, nodes);
          exit(1);
        }
        break;
      case 3:
        if (num_dim == 2) {
          answer = SHELL3;
        }
        else {
          fprintf(stderr,
                  "%sERROR: Shell element with %d nodes "
                  "only valid in 2D.\n",
                  yo, nodes);
          exit(1);
        }
        break;
      case 4: answer = SHELL4; break;
      case 8: answer = SHELL8; break;
      case 9: answer = SHELL9; break;
      default:
        fprintf(stderr, "%sERROR: Shell element with %d nodes unknown.\n", yo, nodes);
        exit(1);
      }
    }
    break;

  case 'c':
    if (strncmp(string, "circle", 6) == 0) {
      answer = SPHERE;
    }
    break;

  case 'b':
  case 't':
  case 'r':
    if (strncmp(string, "bar", 3) == 0 || strncmp(string, "beam", 3) == 0 ||
        strncmp(string, "rod", 3) == 0 || strncmp(string, "truss", 3) == 0) {
      switch (nodes) {
      case 2: answer = BAR2; break;
      case 3: answer = BAR3; break;
      default:
        fprintf(stderr, "%sERROR: Bar/beam/truss elements with %d nodes unknown.\n", yo, nodes);
        exit(1);
      }
    }
    else if (strncmp(string, "tri", 3) == 0) { /* select element shape */
      switch (nodes) {                         /* select number of nodes in this element */
      case 3:                                  /* bilinear triangle */
        if (num_dim == 2) {
          answer = TRI3;
        }
        else {
          answer = TSHELL3;
        }
        break;
      case 4: /* bilinear triangle */
        if (num_dim == 2) {
          answer = TRI4;
        }
        else {
          answer = TSHELL4;
        }
        break;
      case 6: /* biquadratic triangle */
        if (num_dim == 2) {
          answer = TRI6;
        }
        else {
          answer = TSHELL6;
        }
        break;
      case 7: /* biquadratic triangle */
        if (num_dim == 2) {
          answer = TRI7;
        }
        else {
          answer = TSHELL7;
        }
        break;
      default:
        if (num_dim == 2) {
          fprintf(stderr, "%sERROR: triangle element with %d nodes not valid.\n", yo, nodes);
        }
        else {
          fprintf(stderr, "%sERROR: triangle shell element with %d nodes not valid.\n", yo, nodes);
        }
        exit(1);
      }
    }

    else if (strncmp(string, "tet", 3) == 0) { /* select element shape */
      switch (nodes) {                         /* select number of nodes in this element */
      case 4: /* trilinear tetrahedron */ answer = TET4; break;
      case 8: /* 8-node (mid-face) tetrahedron */ answer = TET8; break;
      case 10: /* triquadratic tetrahedron */ answer = TET10; break;
      default:
        fprintf(stderr, "%sERROR: tetrahedral element with %d nodes not valid.\n", yo, nodes);
        exit(1);
      }
    }

    break;

  case 'h':
    /* must check for this before checking for HEX */
    if (strncmp(string, "hexshell", 8) == 0) { /* select element shape */
      switch (nodes) {                         /* select number of nodes in this element */
      case 12: /* only one hexshell */ answer = HEXSHELL; break;
      default:
        fprintf(stderr, "%sERROR: hexshell element with %d nodes not valid.\n", yo, nodes);
        exit(1);
      }
    }

    else if (strncmp(string, "hex", 3) == 0) { /* select element shape */
      switch (nodes) {                         /* select number of nodes in this element */
      case 8: /* trilinear hexahedron */ answer = HEX8; break;
      case 16: /* localization element */ answer = HEX16; break;
      case 20: /* serendipity triquadratic hexahedron */ answer = HEX20; break;
      case 27: /* triquadratic hexahedron */ answer = HEX27; break;
      default:
        fprintf(stderr, "%sERROR: Hexahedron element with %d nodes not valid.\n", yo, nodes);
        exit(1);
      }
    }
    break;

  case 'p':
    if (strncmp(string, "pyra", 4) == 0) { /* select element shape */
      switch (nodes) {                     /* select number of nodes in this element */
      case 5: answer = PYRAMID5; break;
      case 13: answer = PYRAMID13; break;
      case 14: answer = PYRAMID14; break;
      case 18: answer = PYRAMID18; break;
      case 19: answer = PYRAMID19; break;
      default:
        fprintf(stderr, "%sERROR: pyramid element with %d nodes not valid.\n", yo, nodes);
        exit(1);
      }
    }
    break;

  case 'w':
    if (strncmp(string, "wedge", 5) == 0) { /* select element shape */
      switch (nodes) {                      /* select number of nodes in this element */
      case 6: /* trilinear wedge */ answer = WEDGE6; break;
      case 12: /* localization element */ answer = WEDGE12; break;
      case 15: /* triquadratic wedge */ answer = WEDGE15; break;
      case 16: /* triquadratic wedge */ answer = WEDGE16; break;
      default:
        fprintf(stderr, "%sERROR: wedge element with %d nodes not valid.\n", yo, nodes);
        exit(1);
      }
    }
    break;

  default: fprintf(stderr, "%sERROR: Element type %s not supported!\n", yo, string); exit(1);
  }
  /* return desired element information */
  return answer;

} /* get_type */

/*****************************************************************************/
/*			END of el_elm_info.c				     */
/*****************************************************************************/
