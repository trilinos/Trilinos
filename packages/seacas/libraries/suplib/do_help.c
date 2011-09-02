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
/* $Id: do_help.c,v 1.4 2009/03/25 12:46:01 gdsjaar Exp $ */

#include <stdlib.h>
#include <string.h>
#include "help.h"
#include "util.h"

typedef int boolean;
extern int help_c(char *keyword, char *path, boolean *subtopics);
extern int scanner(char *expression);

/* do_help: (not VMS, although it would work)
 * Give help to the user. 
 * It parses the command line into helpbuf and supplies help for that 
 * string. Then, if there are subtopics available for that key,
 * it prompts the user with this string. If more input is
 * given, do_help is called recursively, with the argument the index of 
 * null character in the string. Thus a more specific help can be 
 * supplied. This can be done repeatedly. 
 * If null input is given, the function returns, effecting a
 * backward climb up the tree.
 * David Kotz (dfk@cs.duke.edu) 10/89
 */
void do_help()
{
  static char helpbuf[MAX_LINE_LEN] = "";
  static char prompt[MAX_LINE_LEN] = "";
  int base;				/* index of first char AFTER help string */
  int len;				/* length of current help string */
  BOOLEAN more_help;
  BOOLEAN only;			/* TRUE if only printing subtopics */
  int subtopics;			/* 0 if no subtopics for this topic */
  int start;				/* starting token of help string */
  char *help_ptr;			/* name of help file */

  if ( (help_ptr = (char *)getenv("XHELP")) == (char *)NULL )
    /* if can't find environment variable then just use HELPFILE */
    help_ptr = HELPFILE;

  len = base = strlen(helpbuf);

  /* find the end of the help command */
  for (start = c_token; !(END_OF_COMMAND); c_token++)
    ;
  /* copy new help input into helpbuf */
  if (len > 0)
    helpbuf[len++] = ' ';	/* add a space */
  capture(helpbuf+len, start, c_token-1);
  squash_spaces(helpbuf+base); /* only bother with new stuff */
  lower_case(helpbuf+base); /* only bother with new stuff */
  len = strlen(helpbuf);

  /* now, a lone ? will print subtopics only */
  if (strcmp(helpbuf + (base ? base+1 : 0), "?") == 0) {
    /* subtopics only */
    subtopics = 1;
    only = TRUE;
    helpbuf[base] = '\0';	/* cut off question mark */
  } else {
    /* normal help request */
    subtopics = 0;
    only = FALSE;
  }

  switch (help_c(helpbuf, help_ptr, &subtopics)) {
  case H_FOUND: {
    /* already printed the help info */
    /* subtopics now is true if there were any subtopics */
    screen_ok = FALSE;
    
    do {
      if (subtopics && !only) {
				/* prompt for subtopic with current help string */
	if (len > 0)
	  (void) sprintf(prompt, "Subtopic of `%s': ", helpbuf);
	else
	  (void) strcpy(prompt, "Help topic: ");
	read_line(prompt);
	num_tokens = scanner(input_line);
	c_token = 0;
	more_help = !(END_OF_COMMAND);
	if (more_help)
	  /* base for next level is all of current helpbuf */
	  do_help();
      } else 
	more_help = FALSE;
    } while(more_help);
    
    break;
  }
  case H_NOTFOUND: {
    printf("Sorry, no help for '%s'\n", helpbuf);
    break;
  }
  case H_ERROR: {
    perror(help_ptr);
    break;
  }
  default: {		/* defensive programming */
    int_error("Impossible case in switch\n", NO_CARET);
    /* NOTREACHED */
  }
  }
    
  helpbuf[base] = '\0';	/* cut it off where we started */
}
