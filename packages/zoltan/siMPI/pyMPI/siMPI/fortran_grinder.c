/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/**************************************************************************/
/* FILE   **************     fortran_grinder.c     ************************/
/**************************************************************************/
/* Author: Patrick Miller July 19 2002					  */
/* Copyright (C) 2002 University of California Regents			  */
/**************************************************************************/
/* Convert C code into FORTRAN77 interface code                           */
/**************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include "mpi_config.h"

#ifndef BUFSIZE
#define BUFSIZE 4096
#endif

typedef struct {
  char* ctype;
  char* ftype;
  char* protocol;
} interface_description_type;

/**************************************************************************/
/* LOCAL  **************   interface_description   ************************/
/************************************************************************ **/
/*  */
/**************************************************************************/
static interface_description_type interface_description[] = {
  {"int","MPI_FORTRAN_INTEGER*","*(arg%d)"},
  {"MPI_Comm","MPI_FORTRAN_INTEGER*","*(arg%d)"},
  {"MPI_Comm*","MPI_FORTRAN_INTEGER*","arg%d"},
  {"MPI_Group","MPI_FORTRAN_INTEGER*","*(arg%d)"},
  {"MPI_Group*","MPI_FORTRAN_INTEGER*","arg%d"},
  {"MPI_Request","MPI_FORTRAN_INTEGER*","*(arg%d)"},
  {"MPI_Request*","MPI_FORTRAN_INTEGER*","arg%d"},
  {"MPI_Datatype","MPI_FORTRAN_INTEGER*","*(arg%d)"},
  {"MPI_Datatype*","MPI_FORTRAN_INTEGER*","arg%d"},
  {"MPI_Op","MPI_FORTRAN_INTEGER*","*(arg%d)"},
  {"void*","void*","arg%d"},
  {0,0,0}
};

/**************************************************************************/
/* GLOBAL **************           lookup          ************************/
/************************************************************************ **/
/*  */
/**************************************************************************/
interface_description_type* lookup(char* name) {
  interface_description_type* record = 0;
  for( record = interface_description; record->ctype; ++record) {
    if ( strcmp(record->ctype,name) == 0 ) return record;
  }
  fprintf(stderr,"Error on lookup(\"%s\")\n",name);
  exit(1);
}

/**************************************************************************/
/* LOCAL  **************           mangle          ************************/
/************************************************************************ **/
/*  */
/**************************************************************************/
static void mangle(char* from, char* to) {
  char pattern[BUFSIZE];
  char cased[BUFSIZE];
  char* format = 0;
  char* p = 0;
  char* bar = 0;


  /* Figure out if we have an _ in our name */
  bar = index(MPI_FORTRAN_MANGLING_PATTERNS,'|');
  if ( ! bar ) {
    fprintf(stderr,"Configuration problem\n");
    exit(1);
  }
  if ( index(from,'_') ) {
    strcpy(pattern,bar+1);
  } else {
    strcpy(pattern,MPI_FORTRAN_MANGLING_PATTERNS);
    *bar = 0;
  }

  format = pattern+5;
  strcpy(cased,from);
  if (strncmp(pattern,"lower",5) == 0) {
    for(p=cased;*p;++p) *p = tolower(*p);
  } else {
    for(p=cased;*p;++p) *p = toupper(*p);
  }

  sprintf(to,format,cased);
}

/**************************************************************************/
/* LOCAL  **************         interface         ************************/
/************************************************************************ **/
/*  */
/**************************************************************************/
static void interface(char* rtype, char* ident, char* args) {
  char* parameter[20];
  char* p;
  char* last;
  int i,n;
  interface_description_type* record;
  char mangled[BUFSIZE];

  /* ----------------------------------------------- */
  /* Break out parameters			     */
  /* ----------------------------------------------- */
  for(n=0,last=args,p=args; *p; ++p) {
    if ( *p == ',' ) {
      parameter[n++] = last;
      *p = 0;
      last = p+1;
      while(last && isspace(*last)) ++last;
    }
    if ( isspace(*p) ) *p = 0;
    if ( p > args && *p == '*' && *(p-1) == 0 ) {
      fprintf(stderr,"Bad * placement on argument %d\n",n+1);
      exit(1);
    }
  }
  parameter[n++] = last;

  mangle(ident,mangled);

  printf("void %s(\n",mangled);
  for(i=0;i<n;++i) {
    record = lookup(parameter[i]);
    printf("  %s arg%d,\n",record->ftype,i);
  }
  puts("  MPI_FORTRAN_INTEGER* status)");
  puts("{");
  puts("  int istatus;");
  printf("  istatus = %s(\n",ident);
  for(i=0;i<n;++i) {
    record = lookup(parameter[i]);
    printf("    ");
    printf(record->protocol,i);
    if ( i < n-1) putchar(',');
    putchar('\n');
  }
  puts("  );");
  puts("  *status = istatus;");
  puts("}");

}


/**************************************************************************/
/* GLOBAL **************            main           ************************/
/************************************************************************ **/
/*  */
/**************************************************************************/
int main(int argc, char** argv) {
  char* filename = 0;
  char* body = 0;
  FILE* file = 0;
  char line[BUFSIZE];
  char rtype[BUFSIZE];
  char ident[BUFSIZE];
  char args[BUFSIZE];
  char mangled[BUFSIZE];

  if ( argc <= 1 ) {
    fprintf(stderr,"Usage %s filename {optional body}\n",argv[0]);
    exit(1);
  }
  filename = argv[1];
  if ( argc > 2 ) {
    body = argv[2];
  }


  file = fopen(filename,"r");
  if ( !file ) {
    perror(argv[0]);
    return 1;
  }

  /* ----------------------------------------------- */
  /* look through the file for xxx MPI_xxx(....)     */
  /* ----------------------------------------------- */
  while( fgets(line,BUFSIZE,file) && line[0] ) {
    if ( isspace(line[0]) ) continue;
    rtype[0] = 0;
    ident[0] = 0;
    sscanf(line,"%[a-zA-Z0-9_] %[a-zA-Z0-9_] ( %[^)] )",rtype,ident,args);
    if ( rtype[0] && ident[0] ) {
      puts("#include \"mpi.h\"");
      puts("#include \"mpi_config.h\"");
      puts("#include \"mpi_implementation.h\"");
      
      if ( body ) {
	mangle(ident,mangled);
	printf("void %s%s\n",mangled,body);
      } else {
	interface(rtype,ident,args);
      }
      return 0;
    }
  }
  

  return 0;  
}
