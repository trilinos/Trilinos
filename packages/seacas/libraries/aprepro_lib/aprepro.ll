/* -*- Mode: c++ -*- */
%{

#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>
#include <stdio.h>
#include <string.h>
#include <fcntl.h> 

#include "apr_scanner.h"
#include "aprepro.h"
#include "apr_util.h"

/* import the parser's token type into a local typedef */
typedef SEAMS::Parser::token token;
typedef SEAMS::Parser::token_type token_type;

/* By default yylex returns int, we use token_type. Unfortunately yyterminate
 * by default returns 0, which is not of token_type. */
#define yyterminate() return token::END

#define show(x)   printf("<%s>", x);
 namespace SEAMS {
   extern int echo;
   extern char *get_temp_filename(void);
   extern char *pathopen(const char *file);
   extern void  conv_string(const char *string);
   void yyerror(const char *s);
 }
 
int ifdef;
int file_must_exist = 0; /* Global used by include/conditional include */

/* Global variables used by the looping mechanism */
int loop_lvl = 0;
std::fstream *tmp_file;
char  *temp_f;

#define MAX_IF_NESTING 64

int if_state[MAX_IF_NESTING];
int if_lvl = 0;

%}
/*** Flex Declarations and Options ***/

/* enable c++ scanner class generation */
%option c++

/* change the name of the scanner class. results in "SEAMSFlexLexer" */
%option prefix="SEAMS"

/* enable scanner to generate debug output. disable this for release
 * versions. */
%option debug

/* enables the use of start condition stacks */
%option stack

qstring	\"[^\"\n]*[\"\n]
mlstring \'[^\']*[\']
D     [0-9]
E     [Ee][+-]?{D}+
L     [A-Za-z_]
id    {L}({L}|{D}|:)*
WS    [ \t\f]*
NL    "\n"
number {D}*\.({D}+)?({E})?
integer {D}+({E})?

%START PARSING GET_FILENAME IF_SKIP GET_VAR VERBATIM IF_WHILE_SKIP GET_LOOP_VAR LOOP LOOP_SKIP

%%
<INITIAL>"{VERBATIM(ON)}"   { BEGIN(VERBATIM);  }
<VERBATIM>"{VERBATIM(OFF)}" { BEGIN(INITIAL);   }
<VERBATIM>[A-Za-z0-9_ ]* |
<VERBATIM>.                 { if (echo) ECHO; }
<VERBATIM>"\n"              { if (echo) ECHO; aprepro.ap_file_list.top().lineno++;   }

<INITIAL>{WS}"{ECHO}" |
{WS}"{ECHO(ON)}"	    { echo = true;	}
<INITIAL>{WS}"{NOECHO}" |
{WS}"{ECHO(OFF)}"	    { echo = false;	}

<INITIAL>{WS}"{"[Ll]"oop(" { BEGIN(GET_LOOP_VAR);
			      if (aprepro.ap_options.debugging) 
				std::cerr << "DEBUG LOOP - Found loop begin test " << yytext << " in file "
					  << aprepro.ap_file_list.top().name << "\n";

                           }

<GET_LOOP_VAR>{number}")".*"\n" |
<GET_LOOP_VAR>{integer}")}".*"\n" {/* Loop control defined by integer */
                              char *pt = strchr(yytext, ')');
			      *pt = '\0';
			      sscanf (yytext, "%lf", &yylval->val);

			      if (yylval->val <= 0) {
				BEGIN(LOOP_SKIP);
			      }
			      else {/* Value defined and != 0. */
				temp_f = get_temp_filename();
				SEAMS::file_rec new_file(temp_f, 0, true, (int)yylval->val);
				aprepro.ap_file_list.push(new_file);

				if (aprepro.ap_options.debugging) 
				  std::cerr << "DEBUG LOOP VAR = " << aprepro.ap_file_list.top().loop_count
					    << " in file " << aprepro.ap_file_list.top().name
					    << " at line " << aprepro.ap_file_list.top().lineno << "\n";

				tmp_file = new std::fstream(temp_f, std::ios::out);
				loop_lvl++;
				BEGIN(LOOP);
			      }
			      aprepro.ap_file_list.top().lineno++;
                            }
<GET_LOOP_VAR>.+")}".*"\n"  { /* Loop control defined by variable */
                              symrec *s;
			      char *pt = strchr(yytext, ')');
			      *pt = '\0';
			      s = aprepro.getsym(yytext);

			      if (s == 0 || (s->type != token::SVAR && s->value.var == 0.)) {
				BEGIN(LOOP_SKIP);
			      }
			      else { /* Value defined and != 0. */
				temp_f = get_temp_filename();
				SEAMS::file_rec new_file(temp_f, 0, true, (int)s->value.var);
				aprepro.ap_file_list.push(new_file);
				
				if (aprepro.ap_options.debugging) 
				  std::cerr << "DEBUG LOOP VAR = " << aprepro.ap_file_list.top().loop_count
					    << " in file " << aprepro.ap_file_list.top().name
					    << " at line " << aprepro.ap_file_list.top().lineno << "\n";

				tmp_file = new std::fstream(temp_f, std::ios::out);
				loop_lvl++;
				BEGIN(LOOP);
			      }
			      aprepro.ap_file_list.top().lineno++;
                             }
<LOOP>{WS}"{"[Ee]"nd"[Ll]"oop".*"\n" { aprepro.ap_file_list.top().lineno++;
				   if (--loop_lvl == 0) {
				     BEGIN(INITIAL);
				     tmp_file->close();
				     delete tmp_file;
				     
				     yyin = aprepro.open_file(aprepro.ap_file_list.top().name, "r");
				     yyFlexLexer::yypush_buffer_state (yyFlexLexer::yy_create_buffer( yyin, YY_BUF_SIZE));
				   }
				   else {
				     (*tmp_file) << yytext;
				   }
				 }
<LOOP>{WS}"{"[Ll]"oop(".*"\n"  { loop_lvl++; /* Nested Loop */
	                         (*tmp_file) << yytext;
			         aprepro.ap_file_list.top().lineno++;
			        }
<LOOP>.*"\n"		        { (*tmp_file) << yytext;
 			          aprepro.ap_file_list.top().lineno++;
			        }


<LOOP_SKIP>{WS}"{"[Ee]"nd"[Ll]"oop".*"\n" { aprepro.ap_file_list.top().lineno++;
					if (--loop_lvl == 0)
					  BEGIN(INITIAL);
				      }
<LOOP_SKIP>{WS}"{"[Ll]"oop(".*"\n"        { loop_lvl++; /* Nested Loop */
					aprepro.ap_file_list.top().lineno++;
				      }
<LOOP_SKIP>.*"\n"		      { aprepro.ap_file_list.top().lineno++; }

<IF_SKIP>{WS}"{"[Ii]"fdef("  { if_lvl++; 
    if (aprepro.ap_options.debugging) 
	fprintf (stderr, "DEBUG IF: 'ifdef'  at level = %d at line %d\n",
		 if_lvl, aprepro.ap_file_list.top().lineno);
			   if (if_lvl >= MAX_IF_NESTING)
			     yyerror("Too many nested if statements");
			   if_state[if_lvl] = IF_WHILE_SKIP; }
<IF_SKIP>{WS}"{"[Ii]"fndef(" { if_lvl++; 
    if (aprepro.ap_options.debugging) 
	fprintf (stderr, "DEBUG IF: 'ifndef' at level = %d at line %d\n",
		 if_lvl, aprepro.ap_file_list.top().lineno);
			   if (if_lvl >= MAX_IF_NESTING)
			     yyerror("Too many nested if statements");
			   if_state[if_lvl] = IF_WHILE_SKIP; }
<INITIAL>{WS}"{"[Ii]"fdef("  { if_lvl++; 
    if (aprepro.ap_options.debugging) 
	fprintf (stderr, "DEBUG IF: 'ifdef'  at level = %d at line %d\n",
		 if_lvl, aprepro.ap_file_list.top().lineno);
			   if (if_lvl >= MAX_IF_NESTING)
			     yyerror("Too many nested if statements");
			   ifdef = 1; BEGIN(GET_VAR); }
<INITIAL>{WS}"{"[Ii]"fndef(" { if_lvl++; 
    if (aprepro.ap_options.debugging)
	fprintf (stderr, "DEBUG IF: 'ifndef' at level = %d at line %d\n",
		 if_lvl, aprepro.ap_file_list.top().lineno);
			   if (if_lvl >= MAX_IF_NESTING)
			     yyerror("Too many nested if statements");
			   ifdef = 0; BEGIN(GET_VAR); }

<GET_VAR>.+")}".*"\n"     { symrec *s;
			      char *pt = strchr(yytext, ')');
			      *pt = '\0';
			      s = aprepro.getsym(yytext);
			      if (s == 0 || (s->type != token::SVAR && s->value.var == 0.))
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
			      aprepro.ap_file_list.top().lineno++;
			    }

"{"[Ee]"lse}".*"\n"     { aprepro.ap_file_list.top().lineno++; 
    if (aprepro.ap_options.debugging) 
	fprintf (stderr, "DEBUG IF: 'else'   at level = %d at line %d\n",
		 if_lvl, aprepro.ap_file_list.top().lineno);
			    if (if_state[if_lvl] == IF_SKIP) 
			      BEGIN(INITIAL);
			    if (if_state[if_lvl] == INITIAL)
			      BEGIN(IF_SKIP);
			    /* If neither is true, this is a nested 
			       if that should be skipped */
			  }
"{"[Ee]"ndif}".*"\n"     { if (if_state[if_lvl] == IF_SKIP ||
			       if_state[if_lvl] == INITIAL)
			     BEGIN(INITIAL);
			   /* If neither is true, this is a nested 
			      if that should be skipped */
    if (aprepro.ap_options.debugging) 
	printf ("DEBUG IF: 'endif'  at level = %d at line %d\n",
		if_lvl, aprepro.ap_file_list.top().lineno);
			   if (--if_lvl < 0) {
			     if_lvl = 0;
			     yyerror("Improperly Nested ifdef/ifndef statements");
			   }
			   aprepro.ap_file_list.top().lineno++;  
			   /* Ignore endif if not skipping */ }
<IF_SKIP>[A-Za-z0-9_ ]* |
<IF_SKIP>\\\{           |
<IF_SKIP>\\\}           |
<IF_SKIP>.                 { ; }
<IF_SKIP>"\n"              { aprepro.ap_file_list.top().lineno++; }

<INITIAL>{WS}"{"[Ii]"nclude("           { BEGIN(GET_FILENAME); 
                             file_must_exist = true; }
<INITIAL>{WS}"{"[Cc]"include("          { BEGIN(GET_FILENAME);
                             file_must_exist = !true; }
<GET_FILENAME>.+")"{WS}"}"{NL}* { BEGIN(INITIAL); 
			     {
			       symrec *s;
			       int quoted = false;
			       std::fstream *yytmp;
			       char *pt = strchr(yytext, ')');
			       *pt = '\0';
			       /* Check to see if surrounded by double quote */ 
			       if ((pt = strchr(yytext, '"')) != NULL) {
				 yytext++;
				 quoted = true;
			       }
			       if ((pt = strrchr(yytext, '"')) != NULL) {
				 *pt = '\0';
				 quoted = true;
			       }

			       if (quoted == false) {
				 /* See if this is an aprepro variable referring to a name */
				 s = aprepro.getsym(yytext);
				 if (s == 0 || s->type != token::SVAR) {
				   pt = yytext;
				 } else {
				   pt = (char*)s->value.svar;
				 }
			       } else {
				 pt = yytext;
			       }
			       
			       if (file_must_exist)
				 yytmp = aprepro.open_file(pt, "r");
			       else
				 yytmp = aprepro.check_open_file(pt, "r");
			       if (yytmp != NULL) {
				 yyin = yytmp;
				 if (aprepro.ap_options.info_msg == true) {
				   std::cerr << "Aprepro: INFO: Included File: '"
					     << pt << "' (" << aprepro.ap_file_list.top().name
					     << ", line " << aprepro.ap_file_list.top().lineno
					     << ")\n";
				 }
				 SEAMS::file_rec new_file(pt, 0, false, 0);
				 aprepro.ap_file_list.push(new_file);

				 yyFlexLexer::yypush_buffer_state (
				    yyFlexLexer::yy_create_buffer( yyin, YY_BUF_SIZE));
			       } else {
				 if (aprepro.ap_options.warning_msg == true) {
				   std::cerr << "Aprepro: WARN: Can't open '"
					     << yytext << "'\n";
				 }
			       }
			       aprepro.ap_file_list.top().lineno++;
			     }
			   }


<PARSING>{integer}  |        
<PARSING>{number}	   { sscanf (yytext, "%lf", &yylval->val);
                             return(token::NUM); }

<PARSING>{WS}              ; /* Empty Rule */

<PARSING>{id}              { symrec *s;
			     s = aprepro.getsym(yytext);
			     if (s == 0)
			       s = aprepro.putsym (yytext, SEAMS::Aprepro::UNDEFINED_VARIABLE, 0);
			     yylval->tptr = s;
			     return((token::yytokentype)s->type);
			   }
<PARSING>"="               return(token::EQUAL);
<PARSING>"+="              return(token::EQ_PLUS);
<PARSING>"-="              return(token::EQ_MINUS);
<PARSING>"*="              return(token::EQ_TIME);
<PARSING>"/="              return(token::EQ_DIV);
<PARSING>"^="              return(token::EQ_POW);
<PARSING>"**="             return(token::EQ_POW);
<PARSING>"++"              return(token::INC);
<PARSING>"--"              return(token::DEC);
<PARSING>"+"               return(token::PLU);
<PARSING>"-"               return(token::SUB);
<PARSING>"*"               return(token::TIM);
<PARSING>"~"		   return(token::TIM);		/* ~ is same as multiply */
<PARSING>"//"		   return(token::CONCAT);	/* String concatenation */
<PARSING>"/"               return(token::DIV);
<PARSING>"%"               return(token::MOD);
<PARSING>"^"               return(token::POW);
<PARSING>"**"              return(token::POW);
<PARSING>"\n"              aprepro.ap_file_list.top().lineno++;
<PARSING>"("               return(token::LPAR);
<PARSING>")"               return(token::RPAR);
<PARSING>","               return(token::COMMA);
<PARSING>";"               return(token::SEMI);
<PARSING>":"               return(token::COLON);
<PARSING>"?"               return(token::QUEST);
<PARSING>"<"               return(token::LT);
<PARSING>">"               return(token::GT);
<PARSING>"<="              return(token::LE);
<PARSING>">="              return(token::GE);
<PARSING>"=="              return(token::EQ);
<PARSING>"!="              return(token::NE);
<PARSING>"&&"              return(token::LAND);
<PARSING>"||"              return(token::LOR);
<PARSING>"!"               return(token::NOT);
<PARSING>"["               return(token::LBRACK);
<PARSING>"]"               return(token::RBRACK);
<PARSING>{qstring}	   { char *pt = strrchr(yytext, '"');
			     *pt = '\0';
                             new_string(yytext+1, &yylval->string);
			     return token::QSTRING; }

<PARSING>{mlstring}	   { char *pt = strrchr(yytext, '\'');
			     *pt = '\0';
                             new_string(yytext+1, &yylval->string);
			     return token::QSTRING; }

<PARSING>"}"               { BEGIN(INITIAL); return(token::RBRACE); }

\\\{                      { if (echo) LexerOutput("{", 1); }

\\\}                      { if (echo) LexerOutput("}", 1); }

"{"                        { BEGIN(PARSING); return(token::LBRACE);  }

[Ee][Xx][Ii][Tt] |
[Qq][Uu][Ii][Tt]           { if (aprepro.ap_options.end_on_exit)
			       {
				 if (echo) ECHO;
				 return((token::yytokentype)-1);  
			       }
                              else 
                               if (echo) ECHO;
			   }


\$			   { if (echo) ECHO; }


{id} |
.                          { if (echo) ECHO; }

"\n"                       { if (echo) ECHO; aprepro.ap_file_list.top().lineno++; }

%%

/* When the scanner receives an end-of-file indication from YY_INPUT, it then
 * checks the yywrap() function. If yywrap() returns false (zero), then it is
 * assumed that the function has gone ahead and set up `yyin' to point to
 * another input file, and scanning continues. If it returns true (non-zero),
 * then the scanner terminates, returning 0 to its caller. */

namespace SEAMS {

  Scanner::Scanner(Aprepro& aprepro_yyarg,
		   std::istream* in,
		   std::ostream* out)
    : SEAMSFlexLexer(in, out), aprepro(aprepro_yyarg)
  {
    aprepro.outputStream.push(out);
  }

  Scanner::~Scanner()
  { }

  void Scanner::LexerOutput(const char* buf, int size )
  {
    aprepro.outputStream.top()->write( buf, size );
    if (aprepro.ap_options.interactive && aprepro.outputStream.size() == 1) {
      // In interactive mode, output to stdout in addition to the
      // output stream, unless user has redirected output...
      std::cout << buf;
    }
  }

  int Scanner::yywrap()
  {
    if (aprepro.ap_file_list.size() <= 1) {		/* End of main file, not in nested include */
      return (1);
    }
    else {
      /* We are in an included or looping file */
      if (aprepro.ap_file_list.top().tmp_file) {
	if (aprepro.ap_options.debugging)
	  std::cerr << "DEBUG LOOP: Loop count = " << aprepro.ap_file_list.top().loop_count << "\n";
	if (--aprepro.ap_file_list.top().loop_count <= 0)  {
	  if (strcmp("_string_", aprepro.ap_file_list.top().name.c_str()) != 0) {
	    if (!aprepro.ap_options.debugging)
	      remove(aprepro.ap_file_list.top().name.c_str());	/* Delete file if temporary */
	  }
	  delete yyin;
	  aprepro.ap_file_list.pop(); 
	  yyFlexLexer::yypop_buffer_state();
	}
	else {
	  // Do not pop ap_file_list; we are rereading that file...
	  delete yyin;
	  yyFlexLexer::yypop_buffer_state();
	  yyin = aprepro.open_file(aprepro.ap_file_list.top().name, "r");
	  yyFlexLexer::yypush_buffer_state (yyFlexLexer::yy_create_buffer(yyin, YY_BUF_SIZE));
	  aprepro.ap_file_list.top().lineno = 0;
	}
      }
      else {
	delete yyin;
	yyFlexLexer::yypop_buffer_state();
	aprepro.ap_file_list.pop();
	/* Turn echoing back on at end of included files. */
	echo = true;
      }
      return (0);
    }
  }

  /* Print error message to standard error and return.  Note: internally
   *   'lineno' starts at zero.  To avoid confusion, we add 1 to value
   *   when it is output.
   */

  void Scanner::yyerror (const char *s)
  {
    std::cerr << "Aprepro: ERROR:  " << s << " ("
	      << aprepro.ap_file_list.top().name<< ", line "
	      << aprepro.ap_file_list.top().lineno + 1 << ")\n";
  }

  char *Scanner::execute (char string[])
  {
    /* Push the contents of 'string' onto the stack to be reread.
     * 'string' will be surrounded by {} so it must be a valid expression.
     */

    /*
     * NOTE: The closing } has not yet been scanned in the call to execute();
     *       therefore, we read it ourselves using input(), then we push:
     *       '}{' + our string + '}' 
     */
    int i;
    while ((i = yyFlexLexer::yyinput ()) != '}' && i != EOF)
      ;				/* eat up values */

    /* Allocate space for string + '}' + '{' + end_of_string */
    std::string new_string;
    new_string += "}{";
    new_string += string;
    new_string += "}";

    aprepro.ap_file_list.push(SEAMS::file_rec("_string_", 0, true, -1));
  
    std::istringstream *ins = new std::istringstream(new_string); // Declare an input string stream.
    yyFlexLexer::yypush_buffer_state(yyFlexLexer::yy_create_buffer(ins, new_string.size()));
    return (NULL);
  }

  /* Push the contents of 'string' onto the stack to be reread.
   * 'string' will not be surrounded by {}.
   */

  char *Scanner::rescan (char *string)
  {
    int i;
    /*
     * NOTE: The closing } has not yet been scanned in the call to rescan();
     *       therefore, we read it ourselves using input(), then we push our
     *       string and then put the closing } back on the stack last
     *       (to be read first),
     */
    while ((i = yyFlexLexer::yyinput ()) != '}' && i != EOF)
      ;				/* eat up values */
    {
      aprepro.ap_file_list.push(SEAMS::file_rec("_string_", 0, true, -1));
      std::string new_string("}");
      new_string += string;

      std::istringstream *ins = new std::istringstream(new_string); // Declare an input string stream.
      yyFlexLexer::yypush_buffer_state(yyFlexLexer::yy_create_buffer(ins, new_string.size()));
    }
    return (NULL);
  }
}

/* This implementation of ExampleFlexLexer::yylex() is required to fill the
 * vtable of the class ExampleFlexLexer. We define the scanner's main yylex
 * function via YY_DECL to reside in the Scanner class instead. */

#ifdef yylex
#undef yylex
#endif
int SEAMSFlexLexer::yylex()
{
    std::cerr << "in ExampleFlexLexer::yylex() !" << std::endl;
    return 0;
}

/* When the scanner receives an end-of-file indication from YY_INPUT, it then
 * checks the yywrap() function. If yywrap() returns false (zero), then it is
 * assumed that the function has gone ahead and set up `yyin' to point to
 * another input file, and scanning continues. If it returns true (non-zero),
 * then the scanner terminates, returning 0 to its caller. */

int SEAMSFlexLexer::yywrap()
{
    return 1;
}



