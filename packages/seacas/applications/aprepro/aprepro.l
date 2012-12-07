/* -*- Mode: c++ -*- */
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
extern void  conv_string(char *string);
extern FILE *open_file(char *file, char *mode);
extern FILE *check_open_file(char *file, char *mode);
 
#define MAX_IF_NESTING 64

int if_state[MAX_IF_NESTING];
int if_lvl = 0;

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

%START PARSING GET_FILENAME IF_SKIP GET_VAR VERBATIM IF_WHILE_SKIP GET_LOOP_VAR LOOP LOOP_SKIP

%%
<INITIAL>"{VERBATIM(ON)}"   { BEGIN(VERBATIM);  }
<VERBATIM>"{VERBATIM(OFF)}" { BEGIN(INITIAL);   }
<VERBATIM>[A-Za-z0-9_: ]* |
<VERBATIM>.                 { if (echo) ECHO; }
<VERBATIM>"\n"              { if (echo) ECHO; LINENUM_INC;   }

<INITIAL>{WS}"{ECHO}" |
{WS}"{ECHO(ON)}"	    { echo = True;	}
<INITIAL>{WS}"{NOECHO}" |
{WS}"{ECHO(OFF)}"	    { echo = False;	}

<INITIAL>{WS}"{IMMUTABLE(ON)}"	    { state_immutable = True;	}
<INITIAL>{WS}"{IMMUTABLE(OFF)}"	    { state_immutable = ap_options.immutable; }

<INITIAL>{WS}"{"[Ll]"oop(" { BEGIN(GET_LOOP_VAR); }

<GET_LOOP_VAR>{number}")".*"\n" |
<GET_LOOP_VAR>{integer}")}".*"\n" {/* Loop control defined by integer */
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
<GET_LOOP_VAR>.+")}".*"\n"  { /* Loop control defined by variable */
                              symrec *s;
			      char *pt = strchr(yytext, ')');
			      *pt = '\0';
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
				 LINENUM_INC;
			       }
<LOOP>{WS}"{"[Ee]"nd"[Ll]"oop".*"\n" { LINENUM_INC;
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
<LOOP>{WS}"{"[Ll]"oop(".*"\n"  { loop_lvl++; /* Nested Loop */
			     fprintf(tmp_file, "%s", yytext);
			     LINENUM_INC;
			   }
<LOOP>.*"\n"		   { fprintf(tmp_file, "%s", yytext);
			     LINENUM_INC;
			   }


<LOOP_SKIP>{WS}"{"[Ee]"nd"[Ll]"oop".*"\n" { LINENUM_INC;
					if (--loop_lvl == 0)
					  BEGIN(INITIAL);
				      }
<LOOP_SKIP>{WS}"{"[Ll]"oop(".*"\n"        { loop_lvl++; /* Nested Loop */
					LINENUM_INC;
				      }
<LOOP_SKIP>.*"\n"		      { LINENUM_INC; }

<IF_SKIP>{WS}"{"[Ii]"fdef("  { if_lvl++; 
    if (ap_options.debugging) 
	fprintf (stderr, "DEBUG IF: 'ifdef'  at level = %d at line %d\n",
		 if_lvl, ap_file_list[nfile].lineno);
			   if (if_lvl >= MAX_IF_NESTING)
			     yyerror("Too many nested if statements");
			   if_state[if_lvl] = IF_WHILE_SKIP; }
<IF_SKIP>{WS}"{"[Ii]"fndef(" { if_lvl++; 
    if (ap_options.debugging) 
	fprintf (stderr, "DEBUG IF: 'ifndef' at level = %d at line %d\n",
		 if_lvl, ap_file_list[nfile].lineno);
			   if (if_lvl >= MAX_IF_NESTING)
			     yyerror("Too many nested if statements");
			   if_state[if_lvl] = IF_WHILE_SKIP; }
<INITIAL>{WS}"{"[Ii]"fdef("  { if_lvl++; 
    if (ap_options.debugging) 
	fprintf (stderr, "DEBUG IF: 'ifdef'  at level = %d at line %d\n",
		 if_lvl, ap_file_list[nfile].lineno);
			   if (if_lvl >= MAX_IF_NESTING)
			     yyerror("Too many nested if statements");
			   ifdef = 1; BEGIN(GET_VAR); }
<INITIAL>{WS}"{"[Ii]"fndef(" { if_lvl++; 
    if (ap_options.debugging)
	fprintf (stderr, "DEBUG IF: 'ifndef' at level = %d at line %d\n",
		 if_lvl, ap_file_list[nfile].lineno);
			   if (if_lvl >= MAX_IF_NESTING)
			     yyerror("Too many nested if statements");
			   ifdef = 0; BEGIN(GET_VAR); }

<GET_VAR>.+")}".*"\n"     { symrec *s;
			      char *pt = strchr(yytext, ')');
			      *pt = '\0';
			      s = getsym(yytext);
			      if (s == 0 || (s->type != SVAR && s->type != IMMSVAR && s->value.var == 0.))
				{
				  if (ifdef == 1) {
				    BEGIN(IF_SKIP);
				    if_state[if_lvl] = IF_SKIP;
				  }
				  else {
				    BEGIN(INITIAL);
				    if_state[if_lvl] = INITIAL;
				  }
				}
			      else /* Value defined and != 0. */
				{
				  if (ifdef == 1) {
				    BEGIN(INITIAL);
				    if_state[if_lvl] = INITIAL;
				  }
				  else {
				    BEGIN(IF_SKIP);
				    if_state[if_lvl] = IF_SKIP;
				  }
				}
			      LINENUM_INC;
			    }

{WS}"{"[Ee]"lse}".*"\n"     { LINENUM_INC; 
    if (ap_options.debugging) 
	fprintf (stderr, "DEBUG IF: 'else'   at level = %d at line %d\n",
		 if_lvl, ap_file_list[nfile].lineno);
			    if (if_state[if_lvl] == IF_SKIP) 
			      BEGIN(INITIAL);
			    if (if_state[if_lvl] == INITIAL)
			      BEGIN(IF_SKIP);
			    /* If neither is true, this is a nested 
			       if that should be skipped */
			  }
{WS}"{"[Ee]"ndif}".*"\n"     { if (if_state[if_lvl] == IF_SKIP ||
			       if_state[if_lvl] == INITIAL)
			     BEGIN(INITIAL);
			   /* If neither is true, this is a nested 
			      if that should be skipped */
    if (ap_options.debugging) 
	printf ("DEBUG IF: 'endif'  at level = %d at line %d\n",
		if_lvl, ap_file_list[nfile].lineno);
			   if (--if_lvl < 0) {
			     if_lvl = 0;
			     yyerror("Improperly Nested ifdef/ifndef statements");
			   }
			   LINENUM_INC;  
			   /* Ignore endif if not skipping */ }
<IF_SKIP>[A-Za-z0-9_ ]* |
<IF_SKIP>\\\{           |
<IF_SKIP>\\\}           |
<IF_SKIP>.                 { ; }
<IF_SKIP>"\n"              { LINENUM_INC; }

<INITIAL>{WS}"{"[Ii]"nclude("           { BEGIN(GET_FILENAME); 
                             file_must_exist = True; }
<INITIAL>{WS}"{"[Cc]"include("          { BEGIN(GET_FILENAME);
                             file_must_exist = !True; }
<GET_FILENAME>.+")"{WS}"}"{NL}* { BEGIN(INITIAL); 
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
				 if (nfile+1 >= MAX_INCLUDE_DEPTH)
				   {
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


<PARSING>{integer}  |        
<PARSING>{number}	   { sscanf (yytext, "%lf", &yylval.val);
                             return(NUM); }

<PARSING>{WS}              ; /* Empty Rule */

<PARSING>{id}              { symrec *s;
			     s = getsym(yytext);
			     if (s == 0)
			       s = putsym (yytext, UNDVAR, 0);
			     yylval.tptr = s;
			     return(s->type);
			   }
<PARSING>"="               return('=');
<PARSING>"+="              return(EQ_PLUS);
<PARSING>"-="              return(EQ_MINUS);
<PARSING>"*="              return(EQ_TIME);
<PARSING>"/="              return(EQ_DIV);
<PARSING>"^="              return(EQ_POW);
<PARSING>"**="             return(EQ_POW);
<PARSING>"++"              return(INC);
<PARSING>"--"              return(DEC);
<PARSING>"+"               return('+');
<PARSING>"-"               return('-');
<PARSING>"*"               return('*');
<PARSING>"~"		   return('*');		/* ~ is same as multiply */
<PARSING>"//"		   return(CONCAT);	/* String concatenation */
<PARSING>"/"               return('/');
<PARSING>"%"               return('%');
<PARSING>"^"               return(POW);
<PARSING>"**"              return(POW);
<PARSING>"\n"              LINENUM_INC;
<PARSING>"("               return('(');
<PARSING>")"               return(')');
<PARSING>","               return(',');
<PARSING>";"               return(';');
<PARSING>":"               return(':');
<PARSING>"?"               return('?');
<PARSING>"<"               return('<');
<PARSING>">"               return('>');
<PARSING>"<="              return(LE);
<PARSING>">="              return(GE);
<PARSING>"=="              return(EQ);
<PARSING>"!="              return(NE);
<PARSING>"&&"              return(LAND);
<PARSING>"||"              return(LOR);
<PARSING>"!"               return(NOT);
<PARSING>"["               return('[');
<PARSING>"]"               return(']');
<PARSING>{qstring}	   { char *pt = strrchr(yytext, '"');
			     *pt = '\0';
                             NEWSTR(yytext+1, yylval.string);
			     return QSTRING; }

<PARSING>{mlstring}	   { char *pt = strrchr(yytext, '\'');
			     *pt = '\0';
                             NEWSTR(yytext+1, yylval.string);
			     return QSTRING; }

<PARSING>"}"               { BEGIN(INITIAL); return('}'); }

\\\{                      { if (echo) fwrite( "{", 1, 1, yyout ); }

\\\}                      { if (echo) fwrite( "}", 1, 1, yyout ); }

"{"                        { BEGIN(PARSING); return('{');  }

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

"\n"                       { if (echo) ECHO; LINENUM_INC; }

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

void redefined_warning (char *var)
{
  if (*var != '_' && ap_options.warning_msg == True)
    fprintf (stderr,
             "Aprepro: WARN: Variable '%s' redefined (%s, line %d)\n",
             var, ap_file_list[nfile].name, ap_file_list[nfile].lineno + 1);
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
  while ((i = input ()) != '}' && i != EOF)
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

