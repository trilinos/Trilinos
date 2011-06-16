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

#include "elb_err_const.h"

static int error_cnt = 0;
int error_lev = 1;

static ERROR_MSG_PTR error_info;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function adds the specified error message to the array of error
 * structures.
 *
 * A level 0 error indicates a fatal error, otherwise it's a warning.
 *****************************************************************************/
void error_add(int level, char *message, char *filename, int line_no)
{

  if(error_cnt == 0)
  {
    error_info = (ERROR_MSG_PTR)malloc(sizeof(ERROR_MSG));
    if(!error_info)
    {
      fprintf(stderr, "no memory for error message\n");
      return;
    }
  }
  else
  {
    error_info = (ERROR_MSG_PTR)realloc(error_info,
                                        (error_cnt+1)*sizeof(ERROR_MSG));
    if(!error_info)
    {
      fprintf(stderr, "no memory for error message\n");
      return;
    }
  }

  /* Store the requested error message */
  (error_info+error_cnt)->level    = level;
  (error_info+error_cnt)->err_mesg = (char *)malloc((strlen(message)+1)*
                                                    sizeof(char));
  if(!((error_info+error_cnt)->err_mesg))
  {
    error_info = (ERROR_MSG_PTR)realloc(error_info,
                                        error_cnt*sizeof(ERROR_MSG));
    fprintf(stderr, "no memory for error message\n");
    return;
  }

  strcpy((error_info+error_cnt)->err_mesg, message);

  /* Store the line number info */
  (error_info+error_cnt)->line_no = line_no;

  /* Store the name of the file in which the error occured */
  (error_info+error_cnt)->filename = malloc((strlen(filename)+1)*sizeof(char));
  if(!((error_info+error_cnt)->filename))
  {
    fprintf(stderr, "insufficient memory for entire error message\n");
    error_cnt++;
    return;
  }
  strcpy((error_info+error_cnt)->filename, filename);

  error_cnt++;
  return;

} /*----------------End error_msg()------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function outputs the accumulated error messages to stderr
 *****************************************************************************/
void error_report(void)
{
  int i, iflag=0;


  if(error_lev >= 1)
  {
    for(i=0; i < error_cnt; i++)
    {
      if((error_info+i)->level == 0 || error_lev > 1)
      {
        if(iflag == 0)
        {
          fprintf(stderr, "================================");
          fprintf(stderr, "messages");
          fprintf(stderr, "================================\n");
          iflag = 1;
        }

        fprintf(stderr, "\t%s\n", (error_info+i)->err_mesg);
        if(error_lev >= 2)
          fprintf(stderr, "\t\tin file %s\n", (error_info+i)->filename);

        if(error_lev >= 3)
          fprintf(stderr, "\t\t\tat line %d\n", (error_info+i)->line_no);
      }

    }
  }

  return;
}
