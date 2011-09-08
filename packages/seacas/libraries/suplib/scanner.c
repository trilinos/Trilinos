/* $Id: scanner.c,v 1.3 2009/03/25 12:46:02 gdsjaar Exp $ */

/* GNUPLOT - scanner.c */
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

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "help.h"
#include "util.h"

void substitute(char *str, int max); 	/* substitute output from ` ` */
int get_num(char *str);

#define isident(c) (isalnum(c) || (c) == '_')

#ifndef STDOUT
#define STDOUT 1
#endif

#define LBRACE '{'
#define RBRACE '}'

#define APPEND_TOKEN {token[t_num].length++; current++;}

#define SCAN_IDENTIFIER while (isident((int)expression[current + 1]))\
				APPEND_TOKEN

extern struct lexical_unit token[MAX_TOKENS];

static int t_num;	/* number of token I'm working on */

/*
 * scanner() breaks expression[] into lexical units, storing them in token[].
 *   The total number of tokens found is returned as the function value.
 *   Scanning will stop when '\0' is found in expression[], or when token[]
 *     is full.
 *
 *	 Scanning is performed by following rules:
 *
 *	Current char	token should contain
 *     -------------    -----------------------
 *	1.  alpha	all following alpha-numerics
 *	2.  digit	0 or more following digits, 0 or 1 decimal point,
 *			0 or more digits, 0 or 1 'e' or 'E',
 *			0 or more digits.
 *	3.  ^,+,-,/	only current char
 *	    %,~,(,)
 *	    [,],;,:,
 *	    ?,comma
 *	4.  &,|,=,*	current char; also next if next is same
 *	5.  !,<,>	current char; also next if next is =
 *	6.  ", '	all chars up until matching quote
 *	7.  #          this token cuts off scanning of the line (DFK).
 *
 *	white space between tokens is ignored
 */
int scanner(char *expression)
{
  register int current;	/* index of current char in expression[] */
  register int quote;
  char brace;

  for (current = t_num = 0;
       t_num < MAX_TOKENS && expression[current] != '\0';
       current++) {
       again:
    if (isspace((int)expression[current]))
      continue;				/* skip the whitespace */
    token[t_num].start_index = current;
    token[t_num].length = 1;
    token[t_num].is_token = TRUE;	/* to start with...*/

    if (expression[current] == '`') {
      substitute(&expression[current],MAX_LINE_LEN - current);
      goto again;
    }
    if (isalpha((int)expression[current])) {
      SCAN_IDENTIFIER;
    } else if (isdigit((int)expression[current]) || expression[current] == '.'){
      token[t_num].is_token = FALSE;
      token[t_num].length = get_num(&expression[current]);
      current += (token[t_num].length - 1);
    } else if (expression[current] == LBRACE) {
      token[t_num].is_token = FALSE;
      token[t_num].l_val.type = CMPLX;
      if ((sscanf(&expression[++current],"%lf , %lf %c",
		  &token[t_num].l_val.v.cmplx_val.real,
		  &token[t_num].l_val.v.cmplx_val.imag,
		  &brace) != 3) || (brace != RBRACE))
	int_error("invalid complex constant",t_num);
      token[t_num].length += 2;
      while (expression[++current] != RBRACE) {
	token[t_num].length++;
	if (expression[current] == '\0')			/* { for vi % */
	  int_error("no matching '}'", t_num);
      }
    } else if (expression[current] == '\'' || expression[current] == '\"'){
      token[t_num].length++;
      quote = expression[current];
      while (expression[++current] != quote) {
	if (!expression[current]) {
	  expression[current] = quote;
	  expression[current+1] = '\0';
	  break;
	} else
	  token[t_num].length++;
      }
    } else switch (expression[current]) {
    case '#':		/* DFK: add comments to gnuplot */
      goto endline; /* ignore the rest of the line */
    case '^':
    case '+':
    case '-':
    case '/':
    case '%':
    case '~':
    case '(':
    case ')':
    case '[':
    case ']':
    case ';':
    case ':':
    case '?':
    case ',':
      break;
    case '&':
    case '|':
    case '=':
    case '*':
      if (expression[current] == expression[current + 1])
	APPEND_TOKEN;
      break;
    case '!':
    case '<':
    case '>':
      if (expression[current + 1] == '=')
	APPEND_TOKEN;
      break;
    default:
      int_error("invalid character",t_num);
    }
    ++t_num;	/* next token if not white space */
  }

 endline:					/* comments jump here to ignore line */

  /* Now kludge an extra token which points to '\0' at end of expression[].
     This is useful so printerror() looks nice even if we've fallen off the
     line. */

  token[t_num].start_index = current;
  token[t_num].length = 0;
  return(t_num);
}


int get_num(char *str)
{
  register int count = 0;
  register long lval;

  token[t_num].is_token = FALSE;
  token[t_num].l_val.type = INT;		/* assume unless . or E found */
  while (isdigit((int)str[count]))
    count++;
  if (str[count] == '.') {
    token[t_num].l_val.type = CMPLX;
    while (isdigit((int)str[++count]))
      /* swallow up digits until non-digit */
      ;
    /* now str[count] is other than a digit */
  }
  if (str[count] == 'e' || str[count] == 'E') {
    token[t_num].l_val.type = CMPLX;
    /* modified if statement to allow + sign in exponent
       rjl 26 July 1988 */
    count++;
    if (str[count] == '-' || str[count] == '+')
      count++;
    if (!isdigit((int)str[count])) {
      token[t_num].start_index += count;
      int_error("expecting exponent",t_num);
    }
    while (isdigit((int)str[++count]))
      ;
  }
  if (token[t_num].l_val.type == INT) {
    lval = atol(str);
    if ((token[t_num].l_val.v.int_val = lval) != lval)
      int_error("integer overflow; change to floating point",t_num);
  } else {
    token[t_num].l_val.v.cmplx_val.imag = 0.0;
    token[t_num].l_val.v.cmplx_val.real = atof(str);
  }
  return(count);
}

void substitute(char *str, int max) 	/* substitute output from ` ` */
{
  register char *last;
  register int i,c;
  register FILE *f;
  static char pgm[MAX_LINE_LEN+1],output[MAX_LINE_LEN+1];

  i = 0;
  last = str;
  while (*(++last) != '`') {
    if (*last == '\0')
      int_error("unmatched `",t_num);
    pgm[i++] = *last;
  }
  pgm[i] = '\0';		/* end with null */
  max -= strlen(last);	/* max is now the max length of output sub. */
  
  if ((f = popen(pgm,"r")) == NULL)
    os_error("popen failed",NO_CARET);

  i = 0;
  while ((c = getc(f)) != EOF) {
    output[i++] = ((c == '\n') ? ' ' : c);	/* newlines become blanks*/
    if (i == max) {
      (void) pclose(f);
      int_error("substitution overflow", t_num);
    }
  }
  (void) pclose(f);
  if (i + strlen(last) > max)
    int_error("substitution overflowed rest of line", t_num);
  (void) strncpy(output+i,last+1,MAX_LINE_LEN-i);
  /* tack on rest of line to output */
  (void) strcpy(str,output);		/* now replace ` ` with output */
  screen_ok = FALSE;
}
