/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/**************************************************************************/
/* FILE   **************          mpirun.c         ************************/
/**************************************************************************/
/* Author: Patrick Miller July 19 2002					  */
/* Copyright (C) 2002 University of California Regents			  */
/**************************************************************************/
/*                                                                        */
/**************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

int main(int argc, char** argv) {
  char* arguments[1000];
  int i;

  /* ----------------------------------------------- */
  /* Copy over original arguments, replacing mpirun  */
  /* with the executable                             */
  /* ----------------------------------------------- */
  arguments[0] = argv[3];
  for(i=1;i<argc;++i) arguments[i] = argv[i+3];

  execvp(arguments[0],arguments);
  perror(arguments[0]);
  return 1;
}
