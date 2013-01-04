%{
#include "aprepro.h"
#include "apr_util.h"

#include <stdlib.h>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <cstdio>

namespace SEAMS {
  int   echo = true;
}

%}

/* Require bison 2.4 or later */
%require "2.4"

/* add debug output code to generated parser. disable this for release
 * versions. */
%debug

/* write out a header file containing the token defines */
%defines

/* use newer C++ skeleton file */
%skeleton "lalr1.cc"

/* namespace to enclose parser in */
%name-prefix="SEAMS"

/* set the parser's class identifier */
%define "parser_class_name" "Parser"

%error-verbose

/* aprepro is passed by reference to the parser and to the scanner. This
 * provides a simple but effective pure interface, not relying on global
 * variables. */
%parse-param { class Aprepro& aprepro }

%union {
  double  val;		/* For returning numbers.		*/
  struct symrec *tptr;	/* For returning symbol-table pointers	*/
  char   *string;	/* For returning quoted strings		*/
}

%token	<val>	NUM		/* Simple double precision number	*/
%token	<string> QSTRING	/* Quoted string			*/
%token	<tptr>	UNDVAR 	/* Variable and function		*/
%token  <tptr>  VAR
%token	<tptr>	SVAR 	/* String Variable */
%token  <tptr>  IMMVAR  /* Immutable Variable */
%token	<tptr>	IMMSVAR /* Immutable String Variable */
%token  <tptr>  FNCT
%token  <tptr>  SFNCT
%type	<val>	exp 
%type   <val>   bool
%type	<string>	sexp

%token END 0 "end of file" 
%token COMMA RT LPAR RPAR LBRACK RBRACK LBRACE RBRACE SEMI
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
				     int len = sprintf(tmpstr, format->value.svar, $2);
				     aprepro.lexer->LexerOutput(tmpstr, len);
				   }
                                }
	| LBRACE sexp RBRACE	{ if (echo && $2 != NULL) {
				    aprepro.lexer->LexerOutput($2, strlen($2));
                                  }
                                }
	| LBRACE RBRACE	        {                                       }
	| error RBRACE	        { yyerrok;				}
;

bool:     exp LT exp            { $$ = $1 < $3;                         }
        | exp RT exp            { $$ = $1 > $3;                         }
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

sexp:     QSTRING		{ $$ = $1;				}
        | SVAR			{ $$ = (char*)$1->value.svar;			}
        | IMMSVAR		{ $$ = (char*)$1->value.svar;			}
    	| UNDVAR EQUAL sexp	{ $$ = $3; $1->value.svar = $3;
		                  set_type(aprepro, $1, Parser::token::SVAR);			}
        | SVAR EQUAL sexp	{ $$ = $3; 
				  $1->value.svar = $3;
				  redefined_warning(aprepro, $1->name);          }
        | VAR EQUAL sexp	{ $$ = $3; 
				  $1->value.svar= $3;
				  redefined_warning(aprepro, $1->name);          
		                  set_type(aprepro, $1, token::SVAR);			}
	| IMMSVAR EQUAL sexp	{ immutable_modify(aprepro, $1); YYERROR; }
        | IMMVAR EQUAL sexp	{ immutable_modify(aprepro, $1); YYERROR; }
        | SFNCT LPAR sexp RPAR	{ $$ = (char*)(*($1->value.strfnct_c))($3);	}
	| SFNCT LPAR RPAR	{ $$ = (char*)(*($1->value.strfnct))();	}
        | SFNCT LPAR exp  RPAR	{ $$ = (char*)(*($1->value.strfnct_d))($3);	}
        | sexp CONCAT sexp	{ int len1 = strlen($1);
				  int len3 = strlen($3);
				  $$ = (char*)calloc(1, (len1+len3+1));
				  memcpy($$, $1, len1+1);
				  strcat($$, $3); }
        | SFNCT LPAR exp COMMA sexp COMMA sexp COMMA sexp COMMA sexp RPAR
				{ $$ = (char*)(*($1->value.strfnct_dcccc))($3, $5, $7, $9, $11); }
        | SFNCT LPAR exp COMMA sexp COMMA sexp  RPAR
				{ $$ = (char*)(*($1->value.strfnct_dcc))($3, $5, $7); }
        | SFNCT LPAR sexp COMMA sexp COMMA sexp  RPAR
				{ $$ = (char*)(*($1->value.strfnct_ccc))($3, $5, $7); }
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
				  redefined_warning(aprepro, $1->name);          }
	| SVAR EQUAL exp		{ $$ = $3; $1->value.var = $3;
				  redefined_warning(aprepro, $1->name);          
		                  set_type(aprepro, $1, token::VAR);			}
	| VAR EQ_PLUS exp	{ $1->value.var += $3; $$ = $1->value.var; }
	| VAR EQ_MINUS exp	{ $1->value.var -= $3; $$ = $1->value.var; }
	| VAR EQ_TIME exp	{ $1->value.var *= $3; $$ = $1->value.var; }
	| VAR EQ_DIV exp	{ $1->value.var /= $3; $$ = $1->value.var; }
	| VAR EQ_POW exp	{ errno = 0;
				  $1->value.var = std::pow($1->value.var,$3); 
				  $$ = $1->value.var; 
				  SEAMS::math_error(aprepro, "Power");
				}
        | INC IMMVAR		{ immutable_modify(aprepro, $2); YYERROR; }
	| DEC IMMVAR		{ immutable_modify(aprepro, $2); YYERROR; }
	| IMMVAR INC		{ immutable_modify(aprepro, $1); YYERROR; }
	| IMMVAR DEC		{ immutable_modify(aprepro, $1); YYERROR; }
        | IMMVAR EQUAL exp	{ immutable_modify(aprepro, $1); YYERROR; }
	| IMMSVAR EQUAL exp	{ immutable_modify(aprepro, $1); YYERROR; }
	| IMMVAR EQ_PLUS exp	{ immutable_modify(aprepro, $1); YYERROR; }
	| IMMVAR EQ_MINUS exp	{ immutable_modify(aprepro, $1); YYERROR; }
	| IMMVAR EQ_TIME exp	{ immutable_modify(aprepro, $1); YYERROR; }
	| IMMVAR EQ_DIV exp	{ immutable_modify(aprepro, $1); YYERROR; }
	| IMMVAR EQ_POW exp	{ immutable_modify(aprepro, $1); YYERROR; }

	| UNDVAR		{ $$ = $1->value.var;
				  undefined_warning(aprepro, $1->name);          }
	| INC UNDVAR		{ $$ = ++($2->value.var);		
		                  set_type(aprepro, $2, token::VAR);
				  undefined_warning(aprepro, $2->name);          }
	| DEC UNDVAR		{ $$ = --($2->value.var);		
		                  set_type(aprepro, $2, token::VAR);
				  undefined_warning(aprepro, $2->name);          }
	| UNDVAR INC		{ $$ = ($1->value.var)++;		
		                  set_type(aprepro, $1, token::VAR);
				  undefined_warning(aprepro, $1->name);          }
	| UNDVAR DEC		{ $$ = ($1->value.var)--;		
		                  set_type(aprepro, $1, token::VAR);
				  undefined_warning(aprepro, $1->name);          }
	| UNDVAR EQUAL exp	{ $$ = $3; $1->value.var = $3;
		                  set_type(aprepro, $1, token::VAR);                      }
	| UNDVAR EQ_PLUS exp	{ $1->value.var += $3; $$ = $1->value.var; 
		                  set_type(aprepro, $1, token::VAR);
				  undefined_warning(aprepro, $1->name);          }
	| UNDVAR EQ_MINUS exp	{ $1->value.var -= $3; $$ = $1->value.var; 
		                  set_type(aprepro, $1, token::VAR);
				  undefined_warning(aprepro, $1->name);          }
	| UNDVAR EQ_TIME exp	{ $1->value.var *= $3; $$ = $1->value.var; 
		                  set_type(aprepro, $1, token::VAR);
				  undefined_warning(aprepro, $1->name);          }
	| UNDVAR EQ_DIV exp	{ $1->value.var /= $3; $$ = $1->value.var; 
		                  set_type(aprepro, $1, token::VAR);
				  undefined_warning(aprepro, $1->name);          }
	| UNDVAR EQ_POW exp	{ errno = 0;
				  $1->value.var = std::pow($1->value.var,$3); 
				  $$ = $1->value.var; 
		                  set_type(aprepro, $1, token::VAR);
				  SEAMS::math_error(aprepro, "Power");
				  undefined_warning(aprepro, $1->name);          }
	| FNCT LPAR RPAR		{ $$ = (*($1->value.fnctptr))();	}
	| FNCT LPAR exp RPAR	{ $$ = (*($1->value.fnctptr_d))($3); 	}
	| FNCT LPAR sexp RPAR	{ $$ = (*($1->value.fnctptr_c))($3); 	}
	| FNCT LPAR sexp COMMA exp RPAR
                                { $$ = (*($1->value.fnctptr_cd))($3, $5); 	}
	| FNCT LPAR sexp COMMA sexp RPAR
                                { $$ = (*($1->value.fnctptr_cc))($3, $5); 	}
   	| FNCT LPAR exp COMMA exp RPAR
				{ $$ = (*($1->value.fnctptr_dd))($3, $5);	}
	| FNCT LPAR exp COMMA exp COMMA exp RPAR
				{ $$ = (*($1->value.fnctptr_ddd))($3, $5, $7); }
	| FNCT LPAR exp COMMA exp SEMI exp COMMA exp RPAR
				{ $$ = (*($1->value.fnctptr_dddd))($3, $5, $7, $9); }
	| FNCT LPAR exp COMMA exp COMMA exp COMMA exp RPAR
				{ $$ = (*($1->value.fnctptr_dddd))($3, $5, $7, $9); }
	| FNCT LPAR exp COMMA exp COMMA exp COMMA exp COMMA exp COMMA exp RPAR
  		     { $$ = (*($1->value.fnctptr_dddddd))($3, $5, $7, $9, $11, $13); }
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
	| exp POW exp 		{ errno = 0;
				  $$ = std::pow($1, $3); 
				  SEAMS::math_error(aprepro, "Power");			}
	| LPAR exp RPAR		{ $$ = $2;				}
	| LBRACK exp RBRACK     { errno = 0;
				  $$ = (double)($2 < 0 ? -floor(-($2)): floor($2) );
				  SEAMS::math_error(aprepro, "floor (int)");		}
        | bool                   { $$ = ($1) ? 1 : 0; }
        | bool '?' exp ':' exp   { $$ = ($1) ? ($3) : ($5);              }


/* End of grammar */
%%

void SEAMS::Parser::error(const Parser::location_type&, const std::string& m)
{
    aprepro.error(m);
}

