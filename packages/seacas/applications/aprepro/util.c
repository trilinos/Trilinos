/* 
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
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

/***
   NAME
     util
   PURPOSE
     Utility routines for Aprepro
   NOTES
     do_pathopen - from GAWK source code.  GNU license applies.
        Copyright (C) 1986, 1988, 1989, 1991 the Free Software Foundation, Inc.

   HISTORY
     gdsjaar - Nov 14, 1991: Created.
***/

#ifndef ENVSEP
#if defined(MSDOS) || defined(VMS) || defined(AMIGA)
#define ENVSEP  ';'
#else
#define ENVSEP	':'
#endif
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>

#include <my_aprepro.h>
extern aprepro_options ap_options;

char *pathopen (char *file);
void  conv_string(char *string);
FILE *open_file(char *file, char *mode);
FILE *check_open_file(char *file, char *mode);

static char trypath[BUFSIZ];

char *pathopen (char *file)
{
  static char *savepath;
  static int first = 1;
  char *matspath, *accesspath, *cp, *mats = "/mats";
  FILE *fd;

  if (first)
    {
      first = 0;
      /* get environment variable */
      if ((accesspath = (char *) getenv("ACCESS")) == NULL )
        {
          printf ("Your ACCESS environment variable is not set\n");
          fprintf (stderr, "%s\n", "Please set your ACCESS environment variable and try again");
          exit(1);
        }
      matspath = (char *)malloc (sizeof(char)*(strlen(accesspath)+strlen(mats)+1));
      strcpy(matspath, accesspath);
      strcat(matspath, mats);
      savepath = matspath;	/* used for restarting */
    }
  matspath = savepath;

  do
    {
      trypath[0] = '\0';
      /* this should take into account limits on size of trypath */
      for (cp = trypath; *matspath && *matspath != ENVSEP;)
	*cp++ = *matspath++;

      if (cp != trypath)
	{			/* nun-null element in path */
	  /* add directory punctuation only if needed */
	  if (*(cp - 1) != '/')
	    *cp++ = '/';
	  /* append filename */
	  strcpy (cp, file);
	}
      else
	strcpy (trypath, file);
#ifdef DEBUG
      fprintf (stderr, "Path Tried: %s\n", trypath);
#endif
      if ((fd = fopen (trypath, "r")))
	{
	  fclose (fd);
	  return (trypath);
	}

      /* no luck, keep going */
      if (*matspath == ENVSEP && matspath[1] != '\0')
	matspath++;		/* skip colon */
    }
  while (*matspath);
  /*
	 * You might have one of the mats
	 * paths defined, WITHOUT the current working directory in it.
	 * Therefore try to open the file in the current directory.
	 */
#ifdef DEBUG
  fprintf (stderr, "Path Tried: %s\n", file);
#endif

  if ((fd = fopen (file, "r")))
    {
      fclose (fd);
      return (file);
    }
  else
    return ('\0');
}

/* Convert string to all lower-case and replace all spaces with '_' */
void conv_string (char *string)
{
  char *p = string;
  while (*p != '\0')
    {
      if (*p == ' ')
	*p = '_';
      else if (isupper ((int)*p))
	*p = tolower ((int)*p);
      p++;
    }
}

/* Two methods for opening files. In OPEN_FILE, the file must exist
   or else the code will exit. If CHECK_OPEN_FILE, it is OK if the
   file does not exist. A Null 'pointer' will be returned.
 */
FILE* open_file(char *file, char *mode)
{
  FILE *pointer;
  errno = 0;
  pointer = NULL;
  /* See if file exists in current directory (or as specified) */
  pointer = fopen(file, mode);
  if (pointer == NULL && ap_options.include_path != NULL) {
    /* If there is an include path specified, try opening file there */
    char *file_path;
    int   length;
    length = strlen(ap_options.include_path) + strlen(file) + 2;
    file_path = (char*) malloc(length * sizeof(char));
    strcpy(file_path, ap_options.include_path);
    strcat(file_path, "/");
    strcat(file_path, file);
    
    pointer = fopen(file_path, mode);
    
    free(file_path);
  }

  /* If pointer still null, print error message */
  if (pointer == NULL) {
    char tmpstr[128];
    sprintf(tmpstr, "Aprepro: ERR:  Can't open '%s'",file); 
    perror(tmpstr);
    exit(EXIT_FAILURE);
  }
  return pointer;
}

FILE *check_open_file(char *file, char *mode)
{
  FILE *pointer = NULL;
  pointer = fopen(file, mode);

  if (pointer == NULL && ap_options.include_path != NULL) {
    /* If there is an include path specified, try opening file there */
    char *file_path;
    int   length;
    length = strlen(ap_options.include_path) + strlen(file) + 2;
    file_path = (char*) malloc(length * sizeof(char));
    strcpy(file_path, ap_options.include_path);
    strcat(file_path, "/");
    strcat(file_path, file);
    
    pointer = fopen(file_path, mode);
    
    free(file_path);
  }
  return pointer;
}

/* This function returns a pointer to a static character array.
 * If called multiple times, the name will change on each call,
 * so don't save the pointer; copy the data it points to if you
 * need to keep it for awhile.
 */
char *get_temp_filename()
{
  static char template[] = "./aprepro_temp_XXXXXX";
  int fd;

  strcpy(template, "./aprepro_temp_XXXXXX");
#if defined(__CYGWIN__) && defined(__NO_CYGWIN_OPTION__) 
  fd = mkstemps(template, 0);
#else
  fd = mkstemp(template);
#endif
  close(fd);
  return template;
}  
