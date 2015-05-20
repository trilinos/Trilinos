/* -*- Mode: c++ -*- */
/* 
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
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

%{

#include <stdio.h>
#include <string.h>
#include <fcntl.h> 
#include "ap_getline_int.h"

#define show(x)   printf("<%s>", x);
#undef yywrap
int yywrap(void);
#define MAX_INCLUDE_DEPTH 128
YY_BUFFER_STATE file_pointer[MAX_INCLUDE_DEPTH];
file_rec ap_file_list[MAX_INCLUDE_DEPTH];
extern int nfile;

int ifdef;
int file_must_exist = 0; /* Global used by include/conditional include */
#define LINENUM_INC	ap_file_list[nfile].lineno++

static int get_input(char *buf, int size);
#define YY_INPUT(buf,result,max_size) \
  if ( YY_CURRENT_BUFFER_LVALUE->yy_is_interactive ) {	\
    result = get_input(buf, max_size); \
  } else { 			       \
    errno=0;								\
    while ( (result = fread(buf, 1, max_size, yyin))==0 && ferror(yyin)) \
      {									\
	if( errno != EINTR)						\
	  {								\
	    YY_FATAL_ERROR( "input in flex scanner failed" );		\
	    break;							\
	  }								\
	errno=0;							\
	clearerr(yyin);							\
      }									\
  }
 
extern aprepro_options ap_options;

/* Global variables used by the looping mechanism */
int loop_lvl = 0;
FILE *tmp_file;
char  *temp_f, *tmpstr;

/* Routines used by material database routines */
extern char *get_temp_filename(void);
extern FILE *open_file(char *file, char *mode);
extern FILE *check_open_file(char *file, char *mode);
extern int   check_valid_var(char *var);
void bad_variable_syntax(char *string);
 
#define MAX_IF_NESTING 64

int if_state[MAX_IF_NESTING] = {0}; /* INITIAL */
int if_case_run[MAX_IF_NESTING] = {False}; /* Has any if or elseif condition executed */
int if_lvl = 0;
int if_skip_level = 0;
int suppress_nl = 0; /* Used by some constructs to temporarily suppress a newline */

int switch_active = 0; /* Are we in a switch */
int switch_case_run = 0; /* has there been a case which matched condition run? */
int switch_skip_to_endcase = 0;
double switch_condition = 0.0; /* Value specified in "switch(condition)" */

%}

%option nounput

qstring	\"[^\"\n]*[\"\n]
mlstring \'[^\']*[\']
D     [0-9]
E     [Ee][+-]?{D}+
L     [A-Za-z_]
id    {L}(:|{L}|{D})*
WS    [ \t\f]*
NL    "\n"
number {D}*\.({D}+)?({E})?
integer {D}+({E})?

%START PARSING GET_FILENAME IF_SKIP GET_VAR VERBATIM IF_WHILE_SKIP GET_LOOP_VAR LOOP LOOP_SKIP END_CASE_SKIP

%%
<VERBATIM>{
  "{VERBATIM(OFF)}"   { BEGIN(INITIAL);   }
  [A-Za-z0-9_: ]* | 
  .                   { if (echo) ECHO; }
  "\n"                { if (echo) ECHO; LINENUM_INC; }
}

<INITIAL>{
  "{VERBATIM(ON)}"   { BEGIN(VERBATIM);  }
  {WS}"{ECHO}"   | 
  {WS}"{ECHO(ON)}"	    { echo = True;	}
  {WS}"{NOECHO}" | 
  {WS}"{ECHO(OFF)}"	    { echo = False;	}

  {WS}"{IMMUTABLE(ON)}"	    { state_immutable = True;	}
  {WS}"{IMMUTABLE(OFF)}"    { state_immutable = ap_options.immutable; }

  {WS}"{"[Ll]"oop"{WS}"(" { BEGIN(GET_LOOP_VAR); }
}

<GET_LOOP_VAR>{
  {number}")".*"\n" | 
	    {integer}")}".*"\n" {
    /* Loop control defined by integer */
    char *pt = strchr(yytext, ')');
    *pt = '\0';
    sscanf (yytext, "%lf", &yylval.val);
    ap_file_list[nfile+1].loop_count = (int)yylval.val;
    if (ap_file_list[nfile+1].loop_count <= 0)
      {
	BEGIN(LOOP_SKIP);
      }
    else /* Value defined and != 0. */
      {
	temp_f = get_temp_filename();
	tmp_file = open_file(temp_f, "w");
	loop_lvl++;
	BEGIN(LOOP);
      }
    LINENUM_INC;
  }

  .+")}".*"\n"  {
    /* Loop control defined by variable */
    symrec *s;
    char *pt = strchr(yytext, ')');
    *pt = '\0';
    if (!check_valid_var(yytext)) {
      bad_variable_syntax(yytext);
      BEGIN(LOOP_SKIP);
    } else {
      s = getsym(yytext);
    
      if (s == 0 || (s->type != SVAR && s->type != IMMSVAR && s->value.var == 0.))
	{
	  BEGIN(LOOP_SKIP);
	}
      else /* Value defined and != 0. */
	{
	  ap_file_list[nfile+1].loop_count = s->value.var;
	  temp_f = get_temp_filename();
	  tmp_file = open_file(temp_f, "w");
	  loop_lvl++;
	  BEGIN(LOOP);
	}
    }
    LINENUM_INC;
  }
} /* End of <GET_LOOP_VAR> start condition block */

<LOOP>{
  {WS}"{"[Ee]"nd"[Ll]"oop".*"\n" {
    LINENUM_INC;
    if (--loop_lvl == 0)
      {
	BEGIN(INITIAL);
	fclose(tmp_file);
	file_pointer[nfile++] = YY_CURRENT_BUFFER;
	NEWSTR(temp_f, ap_file_list[nfile].name);
	SET_FILE_LIST(nfile, 0, True, -1);
	yyin = open_file(temp_f, "r");
	yy_switch_to_buffer (
			     yy_create_buffer( yyin, YY_BUF_SIZE));
      }
    else
      fprintf(tmp_file, "%s", yytext);
  }

  {WS}"{"[Ll]"oop"{WS}"(".*"\n"  { 
    loop_lvl++; /* Nested Loop */
    fprintf(tmp_file, "%s", yytext);
    LINENUM_INC;
  }

  .*"\n"  { 
    fprintf(tmp_file, "%s", yytext);
    LINENUM_INC;
  }
} /* End of <LOOP> start condition block */


<LOOP_SKIP>{
  /* IF in a loop that has 0 or negative iteration count, then
   * just need to skip until the end of that loop.  If 
   * there is a nested loop in the skipped loop body,
   * detect that so the correct "endloop" is recognized
   * to turn off the loop skipping
   */
  {WS}"{"[Ee]"nd"[Ll]"oop".*"\n" {
    LINENUM_INC;
    if (--loop_lvl == 0)
      BEGIN(INITIAL);
  }

  {WS}"{"[Ll]"oop"{WS}"(".*"\n" {
    loop_lvl++; /* Nested Loop */
    LINENUM_INC;
  }

  .*"\n" {
    LINENUM_INC; 
  }
} /* End of <LOOP_SKIP> start condition block */


<END_CASE_SKIP>"{case".*"\n"  {
  yyless(0);
  BEGIN(INITIAL);
  switch_skip_to_endcase = 0;
}

<INITIAL,END_CASE_SKIP>"{default}".*"\n"     {
  LINENUM_INC;
  if (!switch_active) {
    yyerror("default statement found outside switch statement.");
  }

  if (!switch_case_run) {
    switch_case_run = 1;
    BEGIN(INITIAL);
    switch_skip_to_endcase = 0;
    if (ap_options.debugging) 
      fprintf (stderr, "DEBUG SWITCH: 'default' code executing at line %d\n",
	       ap_file_list[nfile].lineno);
  } 
  else {
    if (ap_options.debugging) 
      fprintf (stderr, "DEBUG SWITCH: 'default' not executing since a previous case already ran at line %d\n",
	       ap_file_list[nfile].lineno);
    
    /* Need to skip all code until end of case */
    BEGIN(END_CASE_SKIP);
  }
}

<END_CASE_SKIP>"{endswitch}".*"\n"  {
  LINENUM_INC;
  BEGIN(INITIAL);
  switch_active = 0;
  switch_skip_to_endcase = 0;
}

<END_CASE_SKIP>.*"\n" { LINENUM_INC; }

<INITIAL>{WS}"{endswitch}".*"\n"        {
  LINENUM_INC;
  if (!switch_active) {
    yyerror("endswitch statement found without matching switch.");
  }
  switch_active = 0;
}

<IF_WHILE_SKIP>{
  /* If an if was found while skipping, then eat
   * that entire if block until endif
   * found since there is no way that
   * any of the code in that if block could be executed.
   * Make sure to handle multiple levels of skipped ifs...
   *
   * NOTE: if_lvl was not incremented, so don't need to decrement when
   *       endif found.
   */
  {WS}"{"[Ee]"ndif}".*"\n"     { 
    LINENUM_INC;  
    if (--if_skip_level == 0)
      BEGIN(IF_SKIP);
  }

  {WS}"{"[Ii]"fdef"{WS}"(".*"\n"  { 
    LINENUM_INC;  
    if_skip_level++;
  }

  {WS}"{if"{WS}"(".*"\n"  { 
    LINENUM_INC;  
    if_skip_level++;
  }

  {WS}"{"[Ii]"fndef"{WS}"(".*"\n" {
    LINENUM_INC;  
    if_skip_level++;
  }

  .*"\n" {
    LINENUM_INC;
  }
}

<IF_SKIP>{
  /* IF an if, ifdef, or ifndef found while skipping, then
   * skip the entire block up and including the endif.
   * The (IF_WHILE_SKIP) start condition handles this skipping.
   */
  {WS}"{"[Ii]"fdef"{WS}"("  { 
    if (ap_options.debugging) 
      fprintf (stderr, "DEBUG IF: 'ifdef'  found while skipping at line %d\n",
	       ap_file_list[nfile].lineno);
    if_skip_level = 1;
    BEGIN(IF_WHILE_SKIP);
  }

  {WS}"{if"{WS}"("  { 
    if (ap_options.debugging) 
      fprintf (stderr, "DEBUG IF: 'ifdef'  found while skipping at line %d\n",
	       ap_file_list[nfile].lineno);
    if_skip_level = 1;
    BEGIN(IF_WHILE_SKIP);
  }

  {WS}"{"[Ii]"fndef"{WS}"(" {
    if (ap_options.debugging) 
      fprintf (stderr, "DEBUG IF: 'ifndef'  found while skipping at line %d\n",
	       ap_file_list[nfile].lineno);
    if_skip_level = 1;
    BEGIN(IF_WHILE_SKIP);
  }
}

{WS}"{"[Ee]"lse}".*"\n"  { 
  LINENUM_INC; 
  if (ap_options.debugging) 
    fprintf (stderr, "DEBUG IF: 'else'   at level = %d at line %d\n",
	     if_lvl, ap_file_list[nfile].lineno);
  if (if_state[if_lvl] == IF_SKIP) {
    if (!if_case_run[if_lvl]) {
      BEGIN(INITIAL);
      if_state[if_lvl] = INITIAL;
      if_case_run[if_lvl] = True;
    } else {
      BEGIN(IF_SKIP);
      if_state[if_lvl] = IF_SKIP;
    }
  }
  else if (if_state[if_lvl] == INITIAL) {
    BEGIN(IF_SKIP);
    if_state[if_lvl] = IF_SKIP;
  }
  
  /* If neither is true, this is a nested 
     if that should be skipped */
}

<IF_SKIP>{WS}"{"{WS}"elseif".*"\n"  { 
  /* If any previous 'block' of this if has executed, then
   * just skip this block; otherwise see if condition is
   * true and execute this block
   */
  if (ap_options.debugging) 
    fprintf (stderr, "DEBUG IF: 'elseif'   at level = %d at line %d\n",
	     if_lvl, ap_file_list[nfile].lineno);

  if (if_case_run[if_lvl]) { /* A previous else/elseif has run */
    LINENUM_INC; 
    /* Already in IF_SKIP, so don't need to change state */
  } else {
    /* Need to check the elseif condition; push back and parse */
    yyless(0);
    BEGIN(INITIAL);
    if_state[if_lvl] = INITIAL;
  }
}

{WS}"{"[Ee]"ndif}".*"\n"     { 
  if (if_state[if_lvl] == IF_SKIP || if_state[if_lvl] == INITIAL)
    BEGIN(INITIAL);
  /* If neither is true, this is a nested if that should be skipped */
  if (ap_options.debugging) 
    printf ("DEBUG IF: 'endif'  at level = %d at line %d\n",
	    if_lvl, ap_file_list[nfile].lineno);
  if (--if_lvl < 0) {
    if_lvl = 0;
    yyerror("Improperly Nested ifdef/ifndef statements");
  }
  LINENUM_INC;  
  /* Ignore endif if not skipping */
}

<IF_SKIP>[A-Za-z0-9_ ]* |
<IF_SKIP>\\\{           |
<IF_SKIP>\\\}           |
<IF_SKIP>.                 { ; }
<IF_SKIP>"\n"              { LINENUM_INC; }

<INITIAL>{WS}"{"[Ii]"nclude"{WS}"("           { BEGIN(GET_FILENAME); 
                             file_must_exist = True; }
<INITIAL>{WS}"{"[Cc]"include"{WS}"("          { BEGIN(GET_FILENAME);
                             file_must_exist = !True; }
<GET_FILENAME>.+")"{WS}"}"{NL}* {
  BEGIN(INITIAL); 
  {
    symrec *s;
    int quoted = False;
    FILE *yytmp;
    char *pt = strchr(yytext, ')');
    *pt = '\0';
    /* Check to see if surrounded by double quote */ 
    if ((pt = strchr(yytext, '"')) != NULL) {
      yytext++;
      quoted = True;
    }
    if ((pt = strrchr(yytext, '"')) != NULL) {
      *pt = '\0';
      quoted = True;
    }
    
    if (quoted == False) {
      /* See if this is an aprepro variable referring to a name */
      s = getsym(yytext);
      if (s == 0 || (s->type != SVAR && s->type != IMMSVAR)) {
	pt = yytext;
      } else {
	pt = s->value.svar;
      }
    } else {
      pt = yytext;
    }
    
    if (file_must_exist)
      yytmp = open_file(pt, "r");
    else
      yytmp = check_open_file(pt, "r");
    if (yytmp != NULL) {
      yyin = yytmp;
      if (ap_options.info_msg == True) {
	fprintf(stderr,
		"Aprepro: INFO: Included File: '%s' (%s, line %d)\n", pt,
		ap_file_list[nfile].name,
		ap_file_list[nfile].lineno);
      }
      if (nfile+1 >= MAX_INCLUDE_DEPTH) {
	yyerror("Too Many Includes");
	exit(EXIT_FAILURE);
      }
      file_pointer[nfile++] = YY_CURRENT_BUFFER;
      NEWSTR(pt, ap_file_list[nfile].name);
      SET_FILE_LIST(nfile, 0, False, 0);
      yy_switch_to_buffer (
			   yy_create_buffer( yyin, YY_BUF_SIZE));
    } else {
      if (ap_options.warning_msg == True) {
	fprintf(stderr,
		"Aprepro: WARN: Can't open '%s'\n",
		yytext);
      }
    }
    LINENUM_INC;
  }
}

<PARSING>{
  {integer} |        
  {number} {
    sscanf (yytext, "%lf", &yylval.val);
    return(NUM); 
  }

  {WS}              ; /* Empty Rule */

  {id} {
    symrec *s;
    s = getsym(yytext);
    if (s == 0)
      s = putsym (yytext, UNDVAR, 0);
    yylval.tptr = s;
    return(s->type);
  }

  "="               return('=');
  "+="              return(EQ_PLUS);
  "-="              return(EQ_MINUS);
  "*="              return(EQ_TIME);
  "/="              return(EQ_DIV);
  "^="              return(EQ_POW);
  "**="             return(EQ_POW);
  "++"              return(INC);
  "--"              return(DEC);
  "+"               return('+');
  "-"               return('-');
  "*"               return('*');
  "~"               return('*');         /* ~ is same as multiply */
  "//"              return(CONCAT);      /* String concatenation */
  "/"               return('/');
  "%"               return('%');
  "^"               return(POW);
  "**"              return(POW);
  "\n"              LINENUM_INC;
  "("               return('(');
  ")"               return(')');
  ","               return(',');
  ";"               return(';');
  ":"               return(':');
  "?"               return('?');
  "<"               return('<');
  ">"               return('>');
  "<="              return(LE);
  ">="              return(GE);
  "=="              return(EQ);
  "!="              return(NE);
  "&&"              return(LAND);
  "||"              return(LOR);
  "!"               return(NOT);
  "["               return('[');
  "]"               return(']');
  {qstring} {
    char *pt = strrchr(yytext, '"');
    *pt = '\0';
    NEWSTR(yytext+1, yylval.string);
    return QSTRING; 
  }
  
  {mlstring} {
    char *pt = strrchr(yytext, '\'');
    *pt = '\0';
    NEWSTR(yytext+1, yylval.string);
    return QSTRING; 
  }

  "}" {
    if (switch_skip_to_endcase)
      BEGIN(END_CASE_SKIP);
    else
      BEGIN(if_state[if_lvl]);
    return('}');
  }
} /* End of <PARSING> start condition block */

\\\{                      { if (echo) fwrite( "{", 1, 1, yyout ); }

\\\}                      { if (echo) fwrite( "}", 1, 1, yyout ); }

"{"                        { BEGIN(PARSING); return('{'); }

[Ee][Xx][Ii][Tt] |
[Qq][Uu][Ii][Tt]           { if (ap_options.end_on_exit)
			       {
				 if (echo) ECHO;
				 return(-1);  
			       }
                              else 
                               if (echo) ECHO;
			   }


\$			   { if (echo) ECHO; }


{id} |
.                          { if (echo) ECHO; }

"\n"                       { if (echo && !suppress_nl) ECHO; suppress_nl = 0; LINENUM_INC; }

%%
int yywrap(void)
{
  if (nfile == 0)		/* End of main file, not in nested include */
    {
      fclose (yyin);
      return (1);
    }
  else
    /* We are in an included or looping file */
    {
      yy_delete_buffer (YY_CURRENT_BUFFER);
      if (ap_file_list[nfile].tmp_file)
	{
          if (ap_options.debugging)
	      fprintf (stderr, "DEBUG LOOP: Loop count = %d\n",
		       ap_file_list[nfile].loop_count);
	  if (--ap_file_list[nfile].loop_count <= 0)
	    {
	      if (strcmp("_string_", ap_file_list[nfile].name) != 0) {
		
		fclose (yyin);
		if (!ap_options.debugging)
		  remove(ap_file_list[nfile].name);	/* Delete file if temporary */
	      } else {
		/* Free the string buffer */
		/* free(file_pointer[nfile]); */
	      }
	      free (ap_file_list[nfile].name);
	      yy_switch_to_buffer (file_pointer[--nfile]);
	    }
	  else
	    {
	      rewind (yyin);
	      yy_switch_to_buffer (yy_create_buffer (yyin, YY_BUF_SIZE));
	      ap_file_list[nfile].lineno = 0;
	    }
	}
      else
	{
	  fclose (yyin);
	  free (ap_file_list[nfile].name);
	  yy_switch_to_buffer (file_pointer[--nfile]);
	  /* Turn echoing back on at end of included files. */
	  echo = True;
	  /* Set immutability back to global immutable state */
	  state_immutable = ap_options.immutable;	  
	}
      return (0);
    }
}

void set_type(symrec *var, int type)
{
  if (var->name[0] == '_' || !state_immutable) {
    var->type = type;
  } else {
    if (type == VAR)
      var->type = IMMVAR;
    else if (type == SVAR)
      var->type = IMMSVAR;
    else
    var->type = type;
  }
}

/* Print error message to standard error and return.  Note: internally
 *   'lineno' starts at zero.  To avoid confusion, we add 1 to value
 *   when it is output.
 */

void yyerror (char *s)
{
  fprintf (stderr,
	   "Aprepro: ERROR:  %s (%s, line %d)\n",
	   s, ap_file_list[nfile].name, ap_file_list[nfile].lineno + 1);
}

void undefined_warning (char *var)
{
  if (ap_options.warning_msg == True)
    fprintf (stderr,
	     "Aprepro: WARN: Undefined variable '%s' (%s, line %d)\n",
	     var, ap_file_list[nfile].name, ap_file_list[nfile].lineno + 1);
}

void immutable_modify(symrec *var)
{
  fprintf (stderr,
	   "Aprepro: (IMMUTABLE) Variable '%s' is immutable and cannot be modified (%s, line %d)\n",
	   var->name, ap_file_list[nfile].name, ap_file_list[nfile].lineno + 1);
}

void redefined_warning (symrec *var)
{
  if (*var->name != '_' && ap_options.warning_msg == True)
    fprintf (stderr,
             "Aprepro: WARN: %s Variable '%s' redefined (%s, line %d)\n",
	     var->isInternal ? "Pre-defined" : "User-defined",
             var->name, ap_file_list[nfile].name, ap_file_list[nfile].lineno + 1);
}

void bad_variable_syntax(char *string)
{
  if (ap_options.warning_msg == True)
    fprintf (stderr,
	     "Aprepro: WARN: Invalid variable name syntax '%s' (%s, line %d)\n",
	     string, ap_file_list[nfile].name, ap_file_list[nfile].lineno + 1);
}

void warning (char *string)
{
  if (ap_options.warning_msg == True)
    fprintf (stderr,
	     "Aprepro: WARN: '%s' (%s, line %d)\n",
	     string, ap_file_list[nfile].name, ap_file_list[nfile].lineno + 1);
}


char *do_execute (char string[])
{
  /* Push the contents of 'string' onto the stack to be reread.
   * 'string' will be surrounded by {} so it must be a valid expression.
   */
  
  int i;
  char *new_string = NULL;

  /*
   * NOTE: The closing } has not yet been scanned in the call to execute();
   *       therefore, we read it ourselves using input(), then we push:
   *       '}{' + our string + '}' 
   */
  while ((i = input()) != '}' && i != EOF)
    ;				/* eat up values */
  if (i == EOF) {
    yyerror ("End-of-file in execute() command");
  }

  /* Allocate space for string + '}' + '{' + end_of_string */
  new_string = malloc(strlen(string)+4); 
  strcpy(new_string, "}{");
  strcat(new_string, string);
  strcat(new_string, "}");

  file_pointer[nfile++] = YY_CURRENT_BUFFER;
  NEWSTR("_string_", ap_file_list[nfile].name);
  SET_FILE_LIST(nfile, 0, True, -1);
  
  {
    YY_BUFFER_STATE bp = yy_scan_string(new_string);
    file_pointer[nfile] = bp;
    yy_switch_to_buffer(bp);
    free(new_string);
  }

  return (NULL);
}

/* Push the contents of 'string' onto the stack to be reread.
 * 'string' will not be surrounded by {}.
 */

char *do_rescan (char *string)
{
  int i;
/*
 * NOTE: The closing } has not yet been scanned in the call to rescan();
 *       therefore, we read it ourselves using input(), then we push our
 *       string and then put the closing } back on the stack last
 *       (to be read first),
 */
  while ((i = input ()) != '}' && i != EOF)
    ;				/* eat up values */
  if (i == EOF) {
    yyerror ("End-of-file in rescan() command");
  }
  
  {
    YY_BUFFER_STATE bp;
    char *new_string = malloc(strlen(string)+2);
    file_pointer[nfile++] = YY_CURRENT_BUFFER;
    NEWSTR("_string_", ap_file_list[nfile].name);
    SET_FILE_LIST(nfile, 0, True, -1);
    new_string[0] = '}';
    new_string[1] = '\0';
    strcat(new_string, string);
    bp = yy_scan_string(new_string);
    file_pointer[nfile] = bp;
    yy_switch_to_buffer(bp);
    free(new_string);
  }
  return (NULL);
}

char *do_switch(double x)
{
  /* save that we are in a switch statement */
  /* save the value of 'x' for use in deciding which case to execute */
  if (switch_active) {
    yyerror("switch statement found while switch already active. Nested switch not supported.");
  }

  switch_active = 1;
  switch_case_run = 0;
  switch_condition = x;
  switch_skip_to_endcase = 1; /* Skip everything until first case */
  suppress_nl = 1;

  if (ap_options.debugging) 
    fprintf (stderr, "DEBUG SWITCH: 'switch' with condition = %g at line %d\n",
	     switch_condition, ap_file_list[nfile].lineno);

  return(NULL);
}

char *do_case(double x)
{
  /* make sure we are in a switch statement */
  /* if 'x' matches the value saved in the switch statement 
   * and no other case has been executed, then
   * execute the code in the case and set a flag indicating
   * the switch has run;
   * if 'x' does not match the value saved, then skip to endcase
   */
  suppress_nl = 1;

  if (!switch_active) {
    yyerror("case statement found outside switch statement.");
  }

  if (!switch_case_run && x == switch_condition) {
    switch_case_run = 1;
    if (ap_options.debugging) 
      fprintf (stderr, "DEBUG SWITCH: 'case' condition = %g matches switch condition = %g at line %d\n",
	       x, switch_condition, ap_file_list[nfile].lineno);
  } 
  else {
    if (ap_options.debugging) 
      fprintf (stderr, "DEBUG SWITCH: 'case' condition = %g does not match switch condition = %g (or case already matched) at line %d\n",
	       x, switch_condition, ap_file_list[nfile].lineno);

    /* Need to skip all code until end of case */
    switch_skip_to_endcase = 1;
  }
  return(NULL);
}

char *do_if(double x)
{
  if_lvl++;
  if (if_lvl >= MAX_IF_NESTING) {
    yyerror("Too many nested if statements");
  }
  else {
    if (x == 0) {
      if_state[if_lvl] = IF_SKIP;
      if_case_run[if_lvl] = False;
    } else {
      suppress_nl = 1;
      if_state[if_lvl] = INITIAL;
      if_case_run[if_lvl] = True;
    }
    if (ap_options.debugging) 
      fprintf(stderr, "DEBUG IF: If level %d %d\n", if_lvl, if_state[if_lvl]);
  }
  return(NULL);
}

char *do_elseif(double x)
{
  /* Check if in an if */
  /* See if any of the 'cases' in the if have executed */

  if (x == 0 || if_case_run[if_lvl]) {
    if_state[if_lvl] = IF_SKIP;
  } else {
    suppress_nl = 1;
    if_state[if_lvl] = INITIAL;
    if_case_run[if_lvl] = True;
  }
  if (ap_options.debugging) 
    fprintf(stderr, "DEBUG IF: elseif at level %d %d\n", if_lvl, if_state[if_lvl]); 
  return(NULL);
}

char *do_notif(double x)
{
  return do_if(x == 0);
}

char *do_error (char *error_string)
{
  /* Print error message (to stderr) and exit */
  yyerror(error_string);
  exit(EXIT_FAILURE);
  /* NOTREACHED */
  return(NULL);
}

void push_file(char *tmp_file, int is_tmp, int loop_count)
{
    file_pointer[nfile++] = YY_CURRENT_BUFFER;
    NEWSTR(tmp_file, ap_file_list[nfile].name);
    SET_FILE_LIST(nfile, 0, is_tmp, loop_count);
    yyin = open_file(tmp_file, "r");
    yy_switch_to_buffer(yy_create_buffer( yyin, YY_BUF_SIZE));        
}


void add_input_file(char *filename)
{
  yyin = open_file(filename, "r");
  file_pointer[nfile] = yy_create_buffer(yyin, YY_BUF_SIZE);        
  NEWSTR(filename, ap_file_list[nfile].name);
  SET_FILE_LIST(nfile, 0, False, 0);
}

static int get_input(char *buf, int size)
{
  char *line;
  if (feof(yyin))
    return YY_NULL;

  line = getline_int(NULL);

  if (strlen(line) == 0)
    return YY_NULL;

  gl_histadd(line);
  
  if (strlen(line) > size - 2) {
    yyerror("input line too long");
    return YY_NULL;
  }
  strcpy(buf, line);
  strcat(buf, "\n");
  
  return strlen(buf);
}

