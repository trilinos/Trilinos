/* 
 * Copyright (c) 2006 Sandia Corporation. Under the terms of Contract
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

#ifndef APREPRO_H
#define APREPRO_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#ifndef EXIT_FAILURE
# define EXIT_FAILURE 1
# define EXIT_SUCCESS 0
#endif

#define True  (1)
#define False (0)

#ifndef bool
#define bool	char
#endif

extern FILE *yyin;	  /* Input file  */
extern FILE *yyout;       /* Output file */

/* Global options */
struct aprepro_options
{
  char comment;
  char *include_path;
  
  int end_on_exit;
  int warning_msg;
  int info_msg;
  int copyright;
  int quiet;
  int debugging;
  int statistics;
  int interactive;
};

typedef struct aprepro_options aprepro_options;

/* Data type for links in the chain of symbols. */

struct symrec
{
  char *name;
  char *syntax;
  char *info;
  int   type;
  int   isInternal;  /* Only need a bit here; combine with type? */
  union {
    double var;
    double (*fnctptr)();
    char *svar;
    char *(*strfnct)();
  } value;
  struct symrec *next;
};

typedef struct symrec symrec;

/* Structure for holding file names and line counters */
struct file_rec
{
  char 	*name;
  int	lineno;
  bool	tmp_file;
  int	loop_count;
};

typedef struct file_rec file_rec;
extern file_rec ap_file_list[];

symrec	*putsym(char *sym_name, int sym_type, int isInternal);
symrec	*getsym(char*);
void     version(char*);

extern void push_file(char *filename, int is_tmp, int loop_count);
extern FILE* open_file(char *filename, char *mode);

/* Code for saving quoted strings into allocated memory 
   From: "Lex and Yacc," O'Reilly and Associates, p. 95
*/
void yyerror(char *);
#define ALLOC(x,s,t)	do { x = (t)calloc(1, (s)); \
			       if (x == (t)0)  { \
				 yyerror("memory allocation failed");\
					   exit(EXIT_FAILURE); }} while(0)

/* NOTE: the while(0) loop causes lint to return a warning: constant in 
         conditional context.  However, it is valid C and is used to 
	 create a block so we don't get into bad nesting problems.  It
	 also lets the user terminate the macro with a ;
*/
#define NEWSTR(from,to)	do { int len=strlen(from);\
			       ALLOC(to, len+1, char *);\
			       memcpy(to, from, len+1);	\
			       } while(0)

#define SET_FILE_LIST(file_num, line, temp_file, loop_cnt) do { \
				ap_file_list[file_num].lineno = line; \
				ap_file_list[file_num].tmp_file = temp_file; \
				if (loop_cnt >= 0) \
				ap_file_list[file_num].loop_count = loop_cnt; \
			     } while(0)

#define MATH_ERROR(function) do { \
				  if (errno != 0) \
				    yyerror(function); \
				  if (errno == EDOM) \
				    perror("	DOMAIN error"); \
				  else if (errno == ERANGE) \
				    perror("	RANGE error"); \
				  else if (errno != 0) \
				    perror("	Unknown error"); \
				    } while (0)

#endif /* APREPRO_H */



