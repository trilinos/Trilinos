/* $Id: util.c,v 1.10 2009/03/25 12:46:02 gdsjaar Exp $ */

/* GNUPLOT - util.c */
/*
 * Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.
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
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>

#include "help.h"

#undef EXTERN
#define EXTERN 
#include "util.h"
#undef EXTERN

BOOLEAN screen_ok;
/* TRUE if command just typed; becomes FALSE whenever we
   send some other output to screen.  If FALSE, the command line
   will be echoed to the screen before the ^ error message. */

extern char input_line[];
extern struct lexical_unit token[];
extern BOOLEAN interactive;	/* from plot.c */

/*
 * equals() compares string value of token number t_num with str[], and
 *   returns TRUE if they are identical.
 */
int equals(int t_num, char *str)
{
  register int i;

  if (!token[t_num].is_token)
    return(FALSE);	/* must be a value--can't be equal */
  for (i = 0; i < token[t_num].length; i++) {
    if (input_line[token[t_num].start_index+i] != str[i])
      return(FALSE);
  }
  /* now return TRUE if at end of str[], FALSE if not */
  return(str[i] == '\0');
}


/*
 *	capture() copies into str[] the part of input_line[] which lies between
 *	the begining of token[start] and end of token[end].
 */
void capture(char *str, int start, int end)
{
  register int i,e;

  e = token[end].start_index + token[end].length;
  for (i = token[start].start_index; i < e && input_line[i] != '\0'; i++)
    *str++ = input_line[i];
  *str = '\0';
}

/* Lower-case the given string (DFK) */
/* Done in place. */
void lower_case(char *s)
{
  register char *p = s;

  while (*p != '\0') {
    if (isupper((int)*p))
      *p = tolower((int)*p);
    p++;
  }
}

/* Squash spaces in the given string (DFK) */
/* That is, reduce all multiple white-space chars to single spaces */
/* Done in place. */
void squash_spaces(char *s)
{
  register char *r = s;		/* reading point */
  register char *w = s;		/* writing point */
  BOOLEAN space = FALSE;	/* TRUE if we've already copied a space */

  for (w = r = s; *r != '\0'; r++) {
    if (isspace((int)*r)) {
      /* white space; only copy if we haven't just copied a space */
      if (!space) {
	space = TRUE;
	*w++ = ' ';
      }				/* else ignore multiple spaces */
    } else {
      /* non-space character; copy it and clear flag */
      *w++ = *r;
      space = FALSE;
    }
  }
  *w = '\0';				/* null terminate string */
}

void int_error(char *str, int t_num)
{
  register int i;

  /* reprint line if screen has been written to */

  if (t_num != NO_CARET) {		/* put caret under error */
    if (!screen_ok)
      fprintf(stderr,"\n%s%s\n", PROMPT, input_line);

    for (i = 0; i < sizeof(PROMPT) - 1; i++)
      (void) putc(' ',stderr);
    for (i = 0; i < token[t_num].start_index; i++) {
      (void) putc((input_line[i] == '\t') ? '\t' : ' ',stderr);
    }
    (void) putc('^',stderr);
    (void) putc('\n',stderr);
  }

  for (i = 0; i < sizeof(PROMPT) - 1; i++)
    (void) putc(' ',stderr);
  fprintf(stderr,"%s\n\n", str);
}

void os_error(char *str, int t_num)
{
  register int i;
  int save_error;
  save_error = errno; 

  /* reprint line if screen has been written to */

  if (t_num != NO_CARET) {		/* put caret under error */
    if (!screen_ok)
      fprintf(stderr,"\n%s%s\n", PROMPT, input_line);

    for (i = 0; i < sizeof(PROMPT) - 1; i++)
      (void) putc(' ',stderr);
    for (i = 0; i < token[t_num].start_index; i++) {
      (void) putc((input_line[i] == '\t') ? '\t' : ' ',stderr);
    }
    (void) putc('^',stderr);
    (void) putc('\n',stderr);
  }

  for (i = 0; i < sizeof(PROMPT) - 1; i++)
    (void) putc(' ',stderr);
  fprintf(stderr,"%s\n",str);

  for (i = 0; i < sizeof(PROMPT) - 1; i++)
    (void) putc(' ',stderr);

  errno = save_error;
  perror(""); 
  exit(1);
}

void read_line(char *prompt)
{
  int start = 0;
  BOOLEAN more;
  int last = 0;

  if (interactive)
    fputs(prompt,stderr);
  do {
    /* grab some input */
    if ( fgets(&(input_line[start]), MAX_LINE_LEN - start, stdin) 
	 == (char *)NULL ) {
      /* end-of-file */
      if (interactive)
	(void) putc('\n',stderr);
      input_line[start] = '\0';
      if (start > 0)	/* don't quit yet - process what we have */
	more = FALSE;
      else
	exit(IO_SUCCESS); /* no return */
    } else {
      /* normal line input */
      last = strlen(input_line) - 1;
      if (input_line[last] == '\n') { /* remove any newline */
	input_line[last] = '\0';
	/* Watch out that we don't backup beyond 0 (1-1-1) */
	if (last > 0) --last;
      } else if (last+1 >= MAX_LINE_LEN)
	int_error("Input line too long",NO_CARET);
				 
      if (input_line[last] == '\\') { /* line continuation */
	start = last;
	more = TRUE;
      } else
	more = FALSE;
    }
    if (more && interactive)
      fputs("> ", stderr);
  } while(more);
}

/* find char c in string str; return p such that str[p]==c;
 * if c not in str then p=strlen(str)
 */
int instring(char *str, char c)
{
  int pos = 0;

  while (str != NULL && *str != '\0' && c != *str) {
    str++; 
    pos++;
  }
  return (pos);
}
