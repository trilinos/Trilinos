/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dr_err_const.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Author(s): Gary L. Hennigan (SNL 9221)
 *----------------------------------------------------------------------------
 * Functions contained in this file:
 *    error_add()
 *    error_report()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static int error_cnt = 0;
static int error_lev = 3;

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
  (error_info+error_cnt)->filename =
                          (char *) malloc((strlen(filename)+1)*sizeof(char));
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
void error_report(int Proc)
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
          fprintf(stderr, "========================="
                          "messages from Proc %d"
                          "=========================\n", Proc);
          iflag = 1;
        }

        fprintf(stderr, "Proc %d: \t%s\n", Proc, (error_info+i)->err_mesg);
        if(error_lev >= 2)
          fprintf(stderr, "Proc %d: \t\tin file %s\n", 
                           Proc, (error_info+i)->filename);

        if(error_lev >= 3)
          fprintf(stderr, "Proc %d: \t\t\tat line %d\n", 
                           Proc, (error_info+i)->line_no);
      }

    }
  }

  return;
}
