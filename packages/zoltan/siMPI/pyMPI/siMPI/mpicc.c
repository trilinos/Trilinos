/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/**************************************************************************/
/* FILE   **************          mpicc.c          ************************/
/**************************************************************************/
/* Author: Patrick Miller July 15 2002					  */
/* Copyright (C) 2002 University of California Regents			  */
/**************************************************************************/
/*  */
/**************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "mpi_path.h"

#ifdef COMPILER_PATH
char* cc = COMPILER_PATH;
#else
char* cc = CCPATH;
#endif

char* include_directory = IFLAG;
char* lib_directory = LFLAG;
char* lib_library = "-lsimpi";

int main(int argc, char** argv) {
  char* arguments[1000];
  int i;
  int verbose = 0;
  int dash_c = 0;

  /* ----------------------------------------------- */
  /* Copy over original arguments, replacing mpicc   */
  /* with real C compiler                            */
  /* ----------------------------------------------- */
  arguments[0] = cc;
  for(i=1;i<argc;++i) arguments[i] = argv[i];


  /* ----------------------------------------------- */
  /* Check for special flags                         */
  /* ----------------------------------------------- */
  for(i=1;i<argc;++i) {
    if ( strcmp(argv[i],"-c") == 0 ) dash_c = 1;
    if ( strcmp(argv[i],"-v") == 0 ) verbose = 1;
    if ( strcmp(argv[i],"--verbose") == 0 ) verbose = 1;
  }


  /* ----------------------------------------------- */
  /* Add in include directories....                  */
  /* ----------------------------------------------- */
  arguments[argc++] = "-I.";
  arguments[argc++] = include_directory;


  /* ----------------------------------------------- */
  /* If the -c switch is not there, add in link flag */
  /* ----------------------------------------------- */
  if ( !dash_c ) {
    arguments[argc++] = "-L.";
    arguments[argc++] = lib_directory;
    arguments[argc++] = lib_library;
  }
  arguments[argc] = 0;
  

  /* ----------------------------------------------- */
  /* In verbose mode, echo the command               */
  /* ----------------------------------------------- */
  if ( verbose ) {
    for(i=0;i<argc;++i) {
      fprintf(stderr,"%s ",arguments[i]);
    }
    fprintf(stderr,"\n");
  }

  execvp(arguments[0],arguments);
  perror(arguments[0]);
  return 1;
}
