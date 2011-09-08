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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pe_common.h"
#include "el_geom_const.h"
#include "el_elm.h"
#include "rf_mp_const.h"

/************* R O U T I N E S   I N   T H I S   F I L E **********************

    Name_of_Routine		type		     Called by
---------------------------------------------------------------

check_exodus_error ()	      void    	        (many routines)
find_message_info ()	      void		several_routines
ss_to_node_list ()             int		el_exoII_io:read_side_sets()
find_message_info2 ()         void              (not currently used)
******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void check_exodus_error(int error, char function_name[])

/*
 *       Function which checks the error code returned from and EXODUS II API
 *   call and exits if it is -1.
 *
 *       Author:          Scott Hutchinson (1421)
 *       Date:            30 June 1992
 */

{

  if (error == -1) {
    fprintf(stderr, "ERROR returned from %s on Processor %d!\n",
            function_name, Proc);
    exit(1);
  }

} /* check_exodus_error */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void find_message_info(int iunit_size, int max_units,
		       int *num_units_per_message, int *num_messages,
		       int *num_units_left_over)

/*
  Function which determines information regarding the size and number of
  messages to be used by the functions which broadcast the mesh information.

  Author:          Scott Hutchinson (1421)
  Date:            13 March 1993

  */

{

  /**************************** execution begins ******************************/

  if (iunit_size > 0)
    *num_units_per_message = MAX_CHUNK_SIZE/(2*iunit_size);
  else {
    (void) fprintf(stderr,
		   "ERROR:  find_message_info called with unit_size = 0.\n");
    exit(1);
  }

  if (*num_units_per_message > max_units)
    *num_units_per_message = max_units;

  *num_messages        = max_units/(*num_units_per_message);
  *num_units_left_over = max_units - *num_messages*(*num_units_per_message);

  if (max_units%(*num_units_per_message) != 0)
    (*num_messages)++;

} /* find_message_info */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ss_to_node_list(int elem_type, int elem_id, int side_num,
		    int ss_node_list[], int iproc)

{
  int i, itmp1;

  /*
   * This function returns a list of global node numbers forming a
   * side set.
   */

  /* triangle */
  static int tri_table[3][3] = {
  /*   1        2        3                                            side   */
    {1,2,4}, {2,3,5}, {3,1,6}                                      /* nodes  */
  };

  /* quad */
  static int quad_table[4][3] = {
  /*   1        2        3        4                                   side   */
    {1,2,5}, {2,3,6}, {3,4,7}, {4,1,8}                             /* nodes  */
  };

 /* tshell */
  static int tshell_table[2][6] = {
  /*        1                  2                                      side   */
    {1,2,3,4,5,6,}, {1,3,2,6,5,4}                                  /* nodes  */
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

/***************************** execution begins ******************************/

  /* Locally decrement side_num */

  side_num--;

  /* Switch over the element type. */

  switch (elem_type) {

  case SHELL4:
    switch(side_num) {

    case 0:
    case 1:
      for (i = 0; i < 4; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (shell_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;

    default:
      /*
       * sides 2,3,4,5 correspond to
       * sides 0,1,2,3 of the quad element.
       */
      for (i = 0; i < 2; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (quad_table[(side_num-2)][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;
    }
    break;

  case SHELL8:
    switch(side_num) {

    case 0:
    case 1:
      for(i=0; i < 8; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (shell_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;

    default:
      /*
       * sides 2, 3, 4, 5 correspond to
       * sides 0, 1, 2, 3 of the quad element.
       */
      for(i=0; i < 3; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (quad_table[(side_num-2)][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;
    }
    break;

  case TSHELL3:
    switch(side_num) {

    case 0:
    case 1:
      for (i = 0; i < 3; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (tshell_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;

    default:
      /*
       * sides 2, 3, 4 correspond to
       * sides 0, 1, 2 of the tri element.
       */
      for (i = 0; i < 2; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (tri_table[(side_num-2)][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;
    }
    break;

  case TSHELL6:
    switch(side_num) {

    case 0:
    case 1:
      for(i=0; i < 6; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (tshell_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;

    default:
      /*
       * sides 3, 4, & 5 correspond to sides 1, 2, & 3
       * of the tri element.
       */
      for(i=0; i < 3; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (tri_table[(side_num-2)][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;
    }
    break;

  case QUAD4:
    for (i = 0; i < 2; i++) {
      itmp1 = Proc_Connect_Ptr[iproc][elem_id] + (quad_table[side_num][i] - 1);
      ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
    }
    break;

  case QUAD8:
  case QUAD9:
    for (i = 0; i < 3; i++) {
      itmp1 = Proc_Connect_Ptr[iproc][elem_id] + (quad_table[side_num][i] - 1);
      ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
    }
    break;

  case TRI3:
    for (i = 0; i < 2; i++) {
      itmp1 = Proc_Connect_Ptr[iproc][elem_id] + (tri_table[side_num][i] - 1);
      ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
    }
    break;

  case TRI6:
    for (i = 0; i < 3; i++) {
      itmp1 = Proc_Connect_Ptr[iproc][elem_id] + (tri_table[side_num][i] - 1);
      ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
    }
    break;

  case HEX8:
    for (i = 0; i < 4; i++) {
      itmp1 = Proc_Connect_Ptr[iproc][elem_id] + (hex_table[side_num][i] - 1);
      ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
    }
    break;

  case HEX20:
    for (i = 0; i < 8; i++) {
      itmp1 = Proc_Connect_Ptr[iproc][elem_id] + (hex_table[side_num][i] - 1);
      ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
    }
    break;

  case HEX27:
    for (i = 0; i < 9; i++) {
      itmp1 = Proc_Connect_Ptr[iproc][elem_id] + (hex_table[side_num][i] - 1);
      ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
    }
    break;

  case TET4:
    for (i = 0; i < 3; i++) {
      itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
	(tetra_table[side_num][i] - 1);
      ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
    }
    break;

  case TET10:
    for (i = 0; i < 6; i++) {
      itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
	(tetra_table[side_num][i] - 1);
      ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
    }
    break;

  case TET8:
    for (i = 0; i < 4; i++) {
      itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
	(tetra_table[side_num][i] - 1);
      ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
    }
    break;

  case WEDGE6:
    /* NOTE: side_num is decremented above, runs from 0..4 */
    switch(side_num) {

    case 3:
    case 4:
      for(i=0; i < 3; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (wedge_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;

    default:
      for(i=0; i < 4; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (wedge_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;
    }
    break;

  case WEDGE16:
    switch(side_num) {

    case 3:
    case 4:
      for(i=0; i < 6; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (wedge_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;

    default:
      for(i=0; i < 8; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (wedge_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;
    }
    break;

  case WEDGE15:
    switch(side_num) {

    case 3:
    case 4:
      for(i=0; i < 6; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (wedge_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;

    default:
      for(i=0; i < 8; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (wedge_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;
    }
    break;

  case PYRAMID5:
    /* NOTE: side_num is decremented above, runs from 0..4 */
    switch(side_num) {

    case 4:
      for(i=0; i < 4; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (pyramid_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;

    default:
      for(i=0; i < 3; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (pyramid_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;
    }
    break;

  case PYRAMID13:
    switch(side_num) {

    case 4:
      for(i=0; i < 8; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (pyramid_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;

    default:
      for(i=0; i < 6; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (pyramid_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;
    }
    break;

  case HEXSHELL:
    switch (side_num) {
    case 4:
    case 5:
      for (i = 0; i < 4; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (hexshell_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;

    default:
      for (i = 0; i < 6; i++) {
        itmp1 = Proc_Connect_Ptr[iproc][elem_id] +
          (hexshell_table[side_num][i] - 1);
        ss_node_list[i] = Proc_Elem_Connect[iproc][itmp1];
      }
      break;
    }
    break;

  default:
    (void) fprintf(stderr,
		   "ss_to_node_list: ERROR: Unsupported element type!\n");
    exit(1);
  }

  return elem_info(NN_SIDE, elem_type, side_num+1);

} /* End of function "ss_to_node_list" */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void find_message_info2(int iunit_size, int max_units,
                        int *num_units_per_message, int *num_messages,
                        int *num_units_left_over, int max_size)

/*
  Function which determines information regarding the size and number of
  messages to be used by the functions which broadcast the mesh information.

  Author:          Scott Hutchinson (1421)
  Date:            13 March 1993

  */

{

 /**************************** execution begins ******************************/

  if (iunit_size > 0)
    *num_units_per_message = max_size/iunit_size;
  else {
    (void) fprintf(stderr,
		   "ERROR:  find_message_info called with unit_size = 0.\n");
    exit(1);
  }

  if (*num_units_per_message > max_units)
    *num_units_per_message = max_units;

  *num_messages        = max_units/(*num_units_per_message);
  *num_units_left_over = max_units - *num_messages*(*num_units_per_message);

  if (max_units%(*num_units_per_message) != 0)
    (*num_messages)++;

} /* find_message_info2 */

/*****************************************************************************/
/*			END of el_exoII_io_util.c			     */
/*****************************************************************************/
