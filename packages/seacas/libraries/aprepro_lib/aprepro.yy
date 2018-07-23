// Copyright (c) 2014-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.  
// 
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
%{
#include "aprepro.h"
#include "apr_util.h"
#include "apr_array.h"

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <cstdio>
#include <cfenv>

namespace {
  void reset_error()
  {
    if (math_errhandling & MATH_ERREXCEPT) {
      std::feclearexcept(FE_ALL_EXCEPT);
    }
    if (math_errhandling & MATH_ERRNO) {
      errno = 0;
    }
  }
}

namespace SEAMS {
   extern bool echo;
 }

 %}

/* Require bison 2.4 or later */
%require "3.0"

/* add debug output code to generated parser. disable this for release
 * versions. */
%debug

/* write out a header file containing the token defines */
%defines

/* use newer C++ skeleton file */
%skeleton "lalr1.cc"

/* namespace to enclose parser in */
%name-prefix "SEAMS"

/* set the parser's class identifier */
%define "parser_class_name" {Parser}

%error-verbose

/* aprepro is passed by reference to the parser and to the scanner. This
 * provides a simple but effective pure interface, not relying on global
 * variables. */
%parse-param { class Aprepro& aprepro }

%union {
  double  val;		/* For returning numbers.		*/
  struct symrec *tptr;	/* For returning symbol-table pointers	*/
  char   *string;	/* For returning quoted strings		*/
  struct array  *arrval;       /* For returning arrays                 */
}

%token	<val>	NUM		/* Simple double precision number	*/
%token	<string> QSTRING	/* Quoted string			*/
%token	<tptr>	UNDVAR 	/* Variable and function		*/
%token  <tptr>  VAR
%token	<tptr>	SVAR 	/* String Variable */
%token  <tptr>  IMMVAR  /* Immutable Variable */
%token	<tptr>	IMMSVAR /* Immutable String Variable */
%token  <tptr>  AVAR    /* array data [i,j] */
%token  <tptr>  FNCT
%token  <tptr>  SFNCT
%token  <tptr>  AFNCT
%type	<val>	exp 
%type   <arrval> aexp 
%type   <val>   bool
%type	<string>	sexp

%token END 0 "end of file" 
%token COMMA LPAR RPAR LBRACK RBRACK LBRACE RBRACE SEMI
/* Precedence (Lowest to Highest) and associativity */
%right	EQUAL
%right  EQ_PLUS EQ_MINUS
%right  EQ_TIME EQ_DIV
%right  EQ_POW
%right  QUEST COLON
%left   LOR                /* Logical OR     */
%left   LAND               /* Logical AND    */
%left   LT GT LE GE EQ NE  /* <=, >=, ==, != */
%left	PLU SUB
%left	DIV TIM MOD
%left	UNARY NOT 	/* Negation--unary minus/plus 		*/
%right	POW	  	/* Exponentiation	      		*/
%left	INC DEC   	/* increment (++), decrement (--) 	*/
%left	CONCAT	  	/* Concatenate Strings			*/

%{

#include "aprepro.h"
#include "apr_scanner.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex

%}

/* Grammar Rules: */
%%
input:	/* empty rule */
	| input line
;

line:	  '\n'			{ if (echo) aprepro.lexer->LexerOutput("\n", 1); }
	| LBRACE exp RBRACE 	{ if (echo) {
	                             static char tmpstr[512];
				     SEAMS::symrec *format = aprepro.getsym("_FORMAT");
				     int len = sprintf(tmpstr, format->value.svar.c_str(), $2);
				     aprepro.lexer->LexerOutput(tmpstr, len);
				   }
                                }
	| LBRACE sexp RBRACE	{ if (echo && $2 != NULL) {
				    aprepro.lexer->LexerOutput($2, strlen($2));
                                  }
                                }
        | LBRACE aexp RBRACE    {                                       }
	| LBRACE RBRACE	        {                                       }
	| error RBRACE	        { yyerrok;				}
;

bool:     exp LT exp            { $$ = $1 < $3;                         }
        | exp GT exp            { $$ = $1 > $3;                         }
        | NOT exp               { $$ = !($2);                           }
        | exp LE  exp           { $$ = $1 <= $3;                        }
        | exp GE  exp           { $$ = $1 >= $3;                        }
        | exp EQ  exp           { $$ = $1 == $3;                        }
        | exp NE  exp           { $$ = $1 != $3;                        }
        | exp LOR exp           { $$ = $1 || $3;                        }
        | exp LAND exp          { $$ = $1 && $3;                        }
        | bool LOR bool         { $$ = $1 || $3;                        }
        | bool LAND bool        { $$ = $1 && $3;                        }
        | LPAR bool RPAR        { $$ = $2;                              }
;

bool:     sexp LT sexp          { $$ = (strcmp($1,$3) <  0 ? 1 : 0);	}
        | sexp GT sexp          { $$ = (strcmp($1,$3) >  0 ? 1 : 0);	}
        | sexp LE  sexp         { $$ = (strcmp($1,$3) <= 0 ? 1 : 0);	}
        | sexp GE  sexp         { $$ = (strcmp($1,$3) >= 0 ? 1 : 0);	}
        | sexp EQ  sexp         { $$ = (strcmp($1,$3) == 0 ? 1 : 0);	}
        | sexp NE  sexp         { $$ = (strcmp($1,$3) != 0 ? 1 : 0);	}

aexp:   AVAR                    { $$ = $1->value.avar;}
        | AFNCT LPAR sexp RPAR  {
	  if (arg_check($1, $1->value.arrfnct_c == NULL))
	    $$ = (*($1->value.arrfnct_c))($3);
	  else
	    yyerrok;
	}
        | AFNCT LPAR sexp COMMA exp RPAR  {
	  if (arg_check($1, $1->value.arrfnct_cd == NULL))
	    $$ = (*($1->value.arrfnct_cd))($3,$5);
	  else
	    yyerrok;
	}
        | AFNCT LPAR sexp COMMA sexp RPAR  {
	  if (arg_check($1, $1->value.arrfnct_cc == NULL))
	    $$ = (*($1->value.arrfnct_cc))($3,$5);
	  else
	    yyerrok;
	}
        | AFNCT LPAR exp COMMA exp RPAR {
	  if (arg_check($1, $1->value.arrfnct_dd == NULL))
	    $$ = (*($1->value.arrfnct_dd))($3,$5);
	  else
	    yyerrok;
	}
        | AFNCT LPAR exp RPAR {
	  if (arg_check($1, $1->value.arrfnct_d == NULL))
	    $$ = (*($1->value.arrfnct_d))($3);
	  else
	    yyerrok;
	}
        | AFNCT LPAR aexp RPAR  {
	  if (arg_check($1, $1->value.arrfnct_a == NULL))
	    $$ = (*($1->value.arrfnct_a))($3);
	  else
	    yyerrok;
	}
        | AVAR EQUAL aexp       { $$ = $3; delete $1->value.avar; $1->value.avar = $3; 
                                  redefined_warning(aprepro, $1);
                                  set_type(aprepro, $1, token::AVAR); }
        | UNDVAR EQUAL aexp     { $$ = $3; $1->value.avar = $3; 
                                  set_type(aprepro, $1, token::AVAR); }
        | aexp PLU aexp         { if ($1->cols == $3->cols && $1->rows == $3->rows ) {
                                     $$ = array_add($1, $3); 
                                  }
                                  else {
                                    yyerror(aprepro, "Arrays do not have same row and column count"); 
                                    yyerrok;
                                  }
                                }
        | SUB aexp %prec UNARY  { $$ = array_scale($2, -1.0);           }

        | aexp SUB aexp         { if ($1->cols == $3->cols && $1->rows == $3->rows ) {
                                     $$ = array_sub($1, $3); 
                                  }
                                  else {
                                    yyerror(aprepro, "Arrays do not have same row and column count"); 
                                    yyerrok;
                                  }
                                }
        | aexp TIM exp          { $$ = array_scale($1, $3);             }
        | aexp DIV exp          { $$ = array_scale($1, 1.0/$3);         }
        | exp  TIM aexp         { $$ = array_scale($3, $1);             }
        | aexp TIM aexp         { if ($1->cols == $3->rows) {
                                    $$ = array_mult($1, $3);
                                  }
                                  else {
                                    yyerror(aprepro, "Column count of first array does not match row count of second array"); 
                                    yyerrok;
                                  }
				}

sexp:     QSTRING		{ $$ = $1;				}
        | SVAR			{ $$ = (char*)$1->value.svar.c_str();			}
        | IMMSVAR		{ $$ = (char*)$1->value.svar.c_str();			}
    	| UNDVAR EQUAL sexp	{ $$ = $3; $1->value.svar = $3;
		                  set_type(aprepro, $1, Parser::token::SVAR);	}
        | SVAR EQUAL sexp	{ $$ = $3; 
				  $1->value.svar = $3;
				  redefined_warning(aprepro, $1);          }
        | VAR EQUAL sexp	{ $$ = $3; 
				  $1->value.svar= $3;
				  redefined_warning(aprepro, $1);          
		                  set_type(aprepro, $1, token::SVAR);		}
	| IMMSVAR EQUAL sexp	{ $$ = (char*)$1->value.svar.c_str(); immutable_modify(aprepro, $1); }
        | IMMVAR EQUAL sexp	{ immutable_modify(aprepro, $1); YYERROR; }
        | SFNCT LPAR sexp RPAR	{
	  if (arg_check($1, $1->value.strfnct_c == NULL))
	    $$ = (char*)(*($1->value.strfnct_c))($3);
	  else
	    $$ = (char*)"";
	}
	| SFNCT LPAR RPAR	{
	  if (arg_check($1, $1->value.strfnct == NULL))
	    $$ = (char*)(*($1->value.strfnct))();
	  else
	    $$ = (char*)"";
	}
        | SFNCT LPAR exp  RPAR	{
	  if (arg_check($1, $1->value.strfnct_d == NULL))
	    $$ = (char*)(*($1->value.strfnct_d))($3);
	  else
	    $$ = (char*)"";
	}
        | SFNCT LPAR aexp RPAR	{
	  if (arg_check($1, $1->value.strfnct_a == NULL))
	    $$ = (char*)(*($1->value.strfnct_a))($3);
	  else
	    $$ = (char*)"";
	}
        | sexp CONCAT sexp	{ concat_string($1, $3, &$$); }
        | SFNCT LPAR exp COMMA exp RPAR {
	  if (arg_check($1, $1->value.strfnct_dd == NULL))
	    $$ = (char*)(*($1->value.strfnct_dd))($3, $5);
	  else
	    $$ = (char*)"";
	}
        | SFNCT LPAR exp COMMA sexp COMMA sexp COMMA sexp COMMA sexp RPAR {
	  if (arg_check($1, $1->value.strfnct_dcccc == NULL))
	    $$ = (char*)(*($1->value.strfnct_dcccc))($3, $5, $7, $9, $11);
	  else
	    $$ = (char*)"";
	}
        | SFNCT LPAR exp COMMA sexp COMMA sexp  RPAR {
	  if (arg_check($1, $1->value.strfnct_dcc == NULL))
	    $$ = (char*)(*($1->value.strfnct_dcc))($3, $5, $7);
	  else
	    $$ = (char*)"";
	}
        | SFNCT LPAR sexp COMMA sexp COMMA sexp  RPAR {
	  if (arg_check($1, $1->value.strfnct_ccc == NULL))
	    $$ = (char*)(*($1->value.strfnct_ccc))($3, $5, $7);
	  else
	    $$ = (char*)"";
	}
        | SFNCT LPAR sexp COMMA sexp RPAR {
	  if (arg_check($1, $1->value.strfnct_cc == NULL))
	    $$ = (char*)(*($1->value.strfnct_cc))($3, $5);
	  else
	    $$ = (char*)"";
	}
        | bool QUEST sexp COLON sexp  { $$ = ($1) ? ($3) : ($5);              }

exp:	  NUM			{ $$ = $1; 				}
        | INC NUM		{ $$ = $2 + 1;				}
        | DEC NUM		{ $$ = $2 - 1;				}
        | VAR			{ $$ = $1->value.var;			}
        | IMMVAR		{ $$ = $1->value.var;			}
	| INC VAR		{ $$ = ++($2->value.var);		}
	| DEC VAR		{ $$ = --($2->value.var);		}
	| VAR INC		{ $$ = ($1->value.var)++;		}
	| VAR DEC		{ $$ = ($1->value.var)--;		}
	| VAR EQUAL exp		{ $$ = $3; $1->value.var = $3;
				  redefined_warning(aprepro, $1);          }
	| SVAR EQUAL exp		{ $$ = $3; $1->value.var = $3;
				  redefined_warning(aprepro, $1);          
		                  set_type(aprepro, $1, token::VAR);			}
	| VAR EQ_PLUS exp	{ $1->value.var += $3; $$ = $1->value.var; }
	| VAR EQ_MINUS exp	{ $1->value.var -= $3; $$ = $1->value.var; }
	| VAR EQ_TIME exp	{ $1->value.var *= $3; $$ = $1->value.var; }
	| VAR EQ_DIV exp	{ $1->value.var /= $3; $$ = $1->value.var; }
	| VAR EQ_POW exp	{ reset_error();
				  $1->value.var = std::pow($1->value.var,$3); 
				  $$ = $1->value.var; 
				  SEAMS::math_error(aprepro, "Power");
				}
        | INC IMMVAR		{ $$ = $2->value.var; immutable_modify(aprepro, $2);  }
	| DEC IMMVAR		{ $$ = $2->value.var; immutable_modify(aprepro, $2);  }
	| IMMVAR INC		{ $$ = $1->value.var; immutable_modify(aprepro, $1);  }
	| IMMVAR DEC		{ $$ = $1->value.var; immutable_modify(aprepro, $1);  }
        | IMMVAR EQUAL exp	{ $$ = $1->value.var; immutable_modify(aprepro, $1);  }
	| IMMSVAR EQUAL exp	{ immutable_modify(aprepro, $1); YYERROR; }
	| IMMVAR EQ_PLUS exp	{ $$ = $1->value.var; immutable_modify(aprepro, $1);  }
	| IMMVAR EQ_MINUS exp	{ $$ = $1->value.var; immutable_modify(aprepro, $1);  }
	| IMMVAR EQ_TIME exp	{ $$ = $1->value.var; immutable_modify(aprepro, $1);  }
	| IMMVAR EQ_DIV exp	{ $$ = $1->value.var; immutable_modify(aprepro, $1);  }
	| IMMVAR EQ_POW exp	{ $$ = $1->value.var; immutable_modify(aprepro, $1);  }

	| UNDVAR		{ $$ = $1->value.var;
				  undefined_error(aprepro, $1->name);          }
	| INC UNDVAR		{ $$ = ++($2->value.var);		
		                  set_type(aprepro, $2, token::VAR);
				  undefined_error(aprepro, $2->name);          }
	| DEC UNDVAR		{ $$ = --($2->value.var);		
		                  set_type(aprepro, $2, token::VAR);
				  undefined_error(aprepro, $2->name);          }
	| UNDVAR INC		{ $$ = ($1->value.var)++;		
		                  set_type(aprepro, $1, token::VAR);
				  undefined_error(aprepro, $1->name);          }
	| UNDVAR DEC		{ $$ = ($1->value.var)--;		
		                  set_type(aprepro, $1, token::VAR);
				  undefined_error(aprepro, $1->name);          }
	| UNDVAR EQUAL exp	{ $$ = $3; $1->value.var = $3;
		                  set_type(aprepro, $1, token::VAR);                      }
	| UNDVAR EQ_PLUS exp	{ $1->value.var += $3; $$ = $1->value.var; 
		                  set_type(aprepro, $1, token::VAR);
				  undefined_error(aprepro, $1->name);          }
	| UNDVAR EQ_MINUS exp	{ $1->value.var -= $3; $$ = $1->value.var; 
		                  set_type(aprepro, $1, token::VAR);
				  undefined_error(aprepro, $1->name);          }
	| UNDVAR EQ_TIME exp	{ $1->value.var *= $3; $$ = $1->value.var; 
		                  set_type(aprepro, $1, token::VAR);
				  undefined_error(aprepro, $1->name);          }
	| UNDVAR EQ_DIV exp	{ $1->value.var /= $3; $$ = $1->value.var; 
		                  set_type(aprepro, $1, token::VAR);
				  undefined_error(aprepro, $1->name);          }
	| UNDVAR EQ_POW exp	{ reset_error();
				  $1->value.var = std::pow($1->value.var,$3); 
				  $$ = $1->value.var; 
		                  set_type(aprepro, $1, token::VAR);
				  SEAMS::math_error(aprepro, "Power");
				  undefined_error(aprepro, $1->name);          }

        | FNCT LPAR RPAR	{
	  if (arg_check($1, $1->value.fnctptr == NULL))
	    $$ = (*($1->value.fnctptr))();
	  else 
	    $$ = 0.0;
	  }

	| FNCT LPAR exp RPAR	{
	  if (arg_check($1, $1->value.fnctptr_d == NULL))
	    $$ = (*($1->value.fnctptr_d))($3);
	  else
	    $$ = 0.0;
	  }

	| FNCT LPAR sexp RPAR	{
	  if (arg_check($1, $1->value.fnctptr_c == NULL))
	    $$ = (*($1->value.fnctptr_c))($3);
	  else
	    $$ = 0.0;
	  }

	| FNCT LPAR aexp RPAR	{
	  if (arg_check($1, $1->value.fnctptr_a == NULL))
	    $$ = (*($1->value.fnctptr_a))($3);
	  else
	    $$ = 0.0;
	  }

	| FNCT LPAR sexp COMMA exp RPAR {
	    if (arg_check($1, $1->value.fnctptr_cd == NULL))
	      $$ = (*($1->value.fnctptr_cd))($3, $5);
	    else
	      $$ = 0.0;
	  }	    

        | FNCT LPAR exp COMMA sexp RPAR {
	  if (arg_check($1, $1->value.fnctptr_dc == NULL))
	    $$ = (*($1->value.fnctptr_dc))($3, $5);
	  else
	    $$ = 0.0;
	  }

	| FNCT LPAR sexp COMMA sexp RPAR {
	    if (arg_check($1, $1->value.fnctptr_cc == NULL))
	      $$ = (*($1->value.fnctptr_cc))($3, $5);
	    else
	      $$ = 0.0;
	  }

        | FNCT LPAR sexp COMMA sexp COMMA sexp RPAR  {
	  if (arg_check($1, $1->value.fnctptr_ccc == NULL))
	    $$ = (*($1->value.fnctptr_ccc))($3,$5,$7);
	  else
	    yyerrok;
	}

        | FNCT LPAR exp COMMA exp RPAR {
	    if (arg_check($1, $1->value.fnctptr_dd == NULL))
	      $$ = (*($1->value.fnctptr_dd))($3, $5);
	    else
	      $$ = 0.0;
	  }
	| FNCT LPAR exp COMMA exp COMMA exp RPAR {
	    if (arg_check($1, $1->value.fnctptr_ddd == NULL))
	      $$ = (*($1->value.fnctptr_ddd))($3, $5, $7);
	    else
	      $$ = 0.0;
	  }
	| FNCT LPAR sexp COMMA sexp COMMA exp RPAR {
	    if (arg_check($1, $1->value.fnctptr_ccd == NULL))
	      $$ = (*($1->value.fnctptr_ccd))($3, $5, $7);
	    else
	      $$ = 0.0;
	  }
	| FNCT LPAR exp COMMA exp SEMI exp COMMA exp RPAR {
	    if (arg_check($1, $1->value.fnctptr_dddd == NULL))
	      $$ = (*($1->value.fnctptr_dddd))($3, $5, $7, $9);
	    else
	      $$ = 0.0;
	  }
	| FNCT LPAR exp COMMA exp COMMA exp COMMA exp RPAR {
	    if (arg_check($1, $1->value.fnctptr_dddd == NULL))
	      $$ = (*($1->value.fnctptr_dddd))($3, $5, $7, $9);
	    else
	      $$ = 0.0;
	  }
        | FNCT LPAR exp COMMA exp COMMA exp COMMA exp COMMA sexp RPAR {
	    if (arg_check($1, $1->value.fnctptr_ddddc == NULL))
	      $$ = (*($1->value.fnctptr_ddddc))($3, $5, $7, $9, $11);
	    else
	      $$ = 0.0;
	  }
	| FNCT LPAR exp COMMA exp COMMA exp COMMA exp COMMA exp COMMA exp RPAR {
	    if (arg_check($1, $1->value.fnctptr_dddddd == NULL))
	      $$ = (*($1->value.fnctptr_dddddd))($3, $5, $7, $9, $11, $13);
	    else
	      $$ = 0.0;
	  }
	| exp PLU exp		{ $$ = $1 + $3; 			}
	| exp SUB exp		{ $$ = $1 - $3; 			}
	| exp TIM exp		{ $$ = $1 * $3; 			}
	| exp DIV exp		{ if ($3 == 0.)
				    {
				      yyerror(aprepro, "Zero divisor"); 
				      yyerrok;
				    }
				  else
				    $$ = $1 / $3; 			}
	| exp MOD exp		{ if ($3 == 0.)
				    {
				      yyerror(aprepro, "Zero divisor");
				      yyerrok;
				    }
				  else
				    $$ = (int)$1 % (int)$3;		}  
	| SUB exp %prec UNARY	{ $$ = -$2;				}
	| PLU exp %prec UNARY	{ $$ =  $2;				}
	| exp POW exp 		{ reset_error();
				  $$ = std::pow($1, $3); 
				  SEAMS::math_error(aprepro, "Power");			}
	| LPAR exp RPAR		{ $$ = $2;				}
	| LBRACK exp RBRACK     { reset_error();
				  $$ = (double)($2 < 0 ? -floor(-($2)): floor($2) );
				  SEAMS::math_error(aprepro, "floor (int)");		}
        | bool                   { $$ = ($1) ? 1 : 0; }
        | bool QUEST exp COLON exp   { $$ = ($1) ? ($3) : ($5);              }
        | AVAR LBRACK exp RBRACK { $$ = array_value($1->value.avar, $3, 0); }
        | AVAR LBRACK exp COMMA exp RBRACK { $$ = array_value($1->value.avar, $3, $5); }
        | AVAR LBRACK exp RBRACK EQUAL exp 
                                  { $$ = $6;
                                    array *arr = $1->value.avar;
                                    int cols = arr->cols;
				    if (cols > 1) {
                                      yyerror(aprepro, "Cannot use [index] array access with multi-column array"); 
                                      yyerrok;
				    }
                                    int rows = arr->rows;
				    int row = $3;
				    if (aprepro.ap_options.one_based_index) {
				      row--;
				    }
				    if (row < rows) {
                                      int offset = row*cols;
                                      $1->value.avar->data[offset] = $6;
                                    }
                                    else {
                                      yyerror(aprepro, "Row or Column index out of range"); 
                                      yyerrok;
                                    }
                                  }
        | AVAR LBRACK exp COMMA exp RBRACK EQUAL exp 
                                  { $$ = $8;
                                    array *arr = $1->value.avar;
                                    int cols = arr->cols;
                                    int rows = arr->rows;
				    int row = $3;
				    int col = $5;
				    if (aprepro.ap_options.one_based_index) {
				      row--;
				      col--;
				    }
				    if (row < rows && col < cols) {
                                      int offset = row*cols+col;
                                      $1->value.avar->data[offset] = $8;
                                    }
                                    else {
                                      yyerror(aprepro, "Row or Column index out of range"); 
                                      yyerrok;
                                    }
                                  }


/* End of grammar */
%%

void SEAMS::Parser::error(const std::string& m)
{
    aprepro.error(m);
}

