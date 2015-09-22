/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/**************************************************************************/
/* FILE   **************   fortran_parameterize.c  ************************/
/**************************************************************************/
/* Author: Patrick Miller July 19 2002					  */
/* Copyright (C) 2002 University of California Regents			  */
/**************************************************************************/
/*  */
/**************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <ctype.h>

#include "mpi_config.h"

#ifndef BUFSIZE
#define BUFSIZE 4096
#endif


/**************************************************************************/
/* GLOBAL **************            main           ************************/
/************************************************************************ **/
/*  */
/**************************************************************************/
int main(int argc, char** argv) {
  char* filename = 0;
  FILE* file = 0;
  char line[BUFSIZE];
  char ident[BUFSIZE];
  char args[BUFSIZE];
  char cast[BUFSIZE];

  if ( argc <= 1 ) {
    fprintf(stderr,"Usage %s filename\n",argv[0]);
    exit(1);
  }
  filename = argv[1];

  file = fopen(filename,"r");
  if ( !file ) {
    perror(argv[0]);
    return 1;
  }

  /* ----------------------------------------------- */
  /* look through the file for #define MPI_xxx(....) */
  /* ----------------------------------------------- */
  while( fgets(line,BUFSIZE,file) && line[0] ) {
    if ( isspace(line[0]) ) continue;

    /* ----------------------------------------------- */
    /* Look for straight integers                      */
    /* ----------------------------------------------- */
    ident[0] = 0;
    args[0] = 0;
    sscanf(line,"#define %[a-zA-Z0-9_] ( %[-+0-9] )",ident,args);
    if ( ident[0] && args[0]) {
      printf("      INTEGER %s\n",ident);
      printf("      PARAMETER ( %s= (%s) )\n",ident,args);
      continue;
    }

    /* ----------------------------------------------- */
    /* Look for cast integers                          */
    /* ----------------------------------------------- */
    ident[0] = 0;
    cast[0] = 0;
    args[0] = 0;
    sscanf(line,"#define %[a-zA-Z0-9_] ( ( %[a-zA-Z0-9_] ) %[-+0-9] )",ident,cast,args);
    if ( ident[0] && cast[0] && args[0]) {
      printf("      INTEGER %s\n",ident);
      printf("      PARAMETER ( %s= (%s) )\n",ident,args);
      continue;
    }

    /* ----------------------------------------------- */
    /* Look for FORTRAN special lines                  */
    /* ----------------------------------------------- */
    args[0] = 0;
    sscanf(line,"FORTRAN: %[ ,a-zA-Z0-9_]",args);
    if ( args[0] ) {
      printf("      %s\n",args);
    }
  }
  

  return 0;  
}
