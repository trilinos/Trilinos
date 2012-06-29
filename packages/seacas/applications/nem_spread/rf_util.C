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

#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for printf, fprintf, NULL, etc
#include <stdlib.h>                     // for exit
#include "rf_allo.h"                    // for array_alloc

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
void check_exodus_error(int error, const char function_name[])
{

  if (error == -1) {
    fprintf(stderr, "ERROR returned from %s!\n",
            function_name);
    exit(1);
  }

} /* check_exodus_error */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void print_line (const char *charstr, int ntimes)

{
  int i;
  for (i = 0; i < ntimes; i++)
    printf("%c", *charstr);
  printf("\n");
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int break_message_up (size_t unit_size, size_t num_units, size_t max_bytes,
		      int *start_pos[])

/*
 *   break_message_up:
 *
 *         This routine will break a long message up into chunks.  It returns
 * the number of chunks and the starting and ending positions, in terms of
 * unit numbers, of each chunk.
 * It assumes that the first unit of the message has a positional value of 0.
 *
 * NOTE:
 *   This routine does malloc memory.  This memory must be freed in the
 * calling routine, when it is no longer needed.  Use the following statement
 * to free the memory:
 *
 *       safe_free ((void **) &start_pos);
 *
 * Input
 * -------
 *
 * unit_size      = Size in bytes of each "unit" in the message.
 * num_units      = Number of units in the message.
 * max_bytes      = Max_bytes allowed in each message.
 *
 * Output
 * -------
 *
 * return_value   = Number of messages
 * start_pos []   = Starting positions of each message
 *                      (note: start_pos[0] = 0, by definition).
 *                    It is defined to have a length of num_mesg+1
 *                    start_pos[num_mesg] = num_units
 *
 *
 * Usage
 *----------
 *       To get the length of message i:
 *
 *                   mesg_length = start_pos[i+1] - start_pos[i]
 *
 */

{
  size_t  i, num_units_per_message, remainder, num_mesg;

 /*------------------------- begin execution --------------------------------*/

  if (num_units <= 0) {
    num_mesg = 0;
    *start_pos = NULL;
    return (num_mesg);
  }

  num_units_per_message = max_bytes / unit_size ;
  if (num_units < num_units_per_message) num_units_per_message = num_units;
  num_mesg   = num_units / num_units_per_message;
  remainder  = num_units % num_units_per_message;
  if (remainder != 0) num_mesg++;

  *start_pos = (int *) array_alloc (__FILE__, __LINE__, 1, (num_mesg + 1),
                                    sizeof (int));

  for (i = 0; i < num_mesg; i++) {
    (*start_pos)[i] =  i * num_units_per_message;
  }
  (*start_pos) [num_mesg] = num_units;

  return num_mesg;
}

/*****************************************************************************/
/*                END OF FILE rf_util.c					     */
/*****************************************************************************/
