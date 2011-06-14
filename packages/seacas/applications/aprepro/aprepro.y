%{
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
#include "my_aprepro.h"
#include <stdlib.h>

void undefined_warning(char* var);
void redefined_warning(char* var);
void yyerror(char* var);
void warning(char *string);
int  yylex(void);

int   echo = True;

symrec *format;
%}

%union {
  double  val;		/* For returning numbers.		*/
  symrec *tptr;		/* For returning symbol-table pointers	*/
  char   *string;	/* For returning quoted strings		*/
}

%token	<val>	NUM		/* Simple double precision number	*/
%token	<string> QSTRING	/* Quoted string			*/
%token	<tptr>	UNDVAR 	/* Variable and function		*/
%token  <tptr>  VAR
%token	<tptr>	SVAR 	/* String Variable and function		*/
%token  <tptr>  FNCT
%token  <tptr>  SFNCT
%type	<val>	exp 
%type   <val>   bool
%type	<string>	sexp

/* Precedence (Lowest to Highest) and associativity */
%right	'=' 
%right  EQ_PLUS EQ_MINUS
%right  EQ_TIME EQ_DIV
%right  EQ_POW
%right '?' ':'
%left   LOR                /* Logical OR     */
%left   LAND               /* Logical AND    */
%left '<' '>' LE GE EQ NE  /* <=, >=, ==, != */
%left	'-' '+'
%left	'/' '*' '%'
%left	UNARY NOT 	/* Negation--unary minus/plus 		*/
%right	POW	  	/* Exponentiation	      		*/
%left	INC DEC   	/* increment (++), decrement (--) 	*/
%left	CONCAT	  	/* Concatenate Strings			*/
/* Grammar Rules: */

%%
input:	/* empty rule */
	| input line
;

line:	  '\n'			{ if (echo) fprintf(yyout,"\n");	}
	| '{' exp  '}' 		{ if (echo) {
	                             format = getsym("_FORMAT");
	                             fprintf(yyout, format->value.svar, $2);
				   }                                    }
	| '{' sexp '}' 		{ if (echo && $2 != NULL)
				    fprintf(yyout, "%s", $2);	}
	| error '}'		{ yyerrok;				}
;

bool:     exp '<' exp           { $$ = $1 < $3;                         }
        | exp '>' exp           { $$ = $1 > $3;                         }
        | NOT exp               { $$ = !($2);                           }
        | exp LE  exp           { $$ = $1 <= $3;                        }
        | exp GE  exp           { $$ = $1 >= $3;                        }
        | exp EQ  exp           { $$ = $1 == $3;                        }
        | exp NE  exp           { $$ = $1 != $3;                        }
        | exp LOR exp           { $$ = $1 || $3;                        }
        | exp LAND exp          { $$ = $1 && $3;                        }
        | bool LOR bool         { $$ = $1 || $3;                        }
        | bool LAND bool        { $$ = $1 && $3;                        }
        | '(' bool ')'          { $$ = $2;                              }
;

bool:     sexp '<' sexp         { $$ = (strcmp($1,$3) <  0 ? 1 : 0);	}
        | sexp '>' sexp         { $$ = (strcmp($1,$3) >  0 ? 1 : 0);	}
        | sexp LE  sexp         { $$ = (strcmp($1,$3) <= 0 ? 1 : 0);	}
        | sexp GE  sexp         { $$ = (strcmp($1,$3) >= 0 ? 1 : 0);	}
        | sexp EQ  sexp         { $$ = (strcmp($1,$3) == 0 ? 1 : 0);	}
        | sexp NE  sexp         { $$ = (strcmp($1,$3) != 0 ? 1 : 0);	}

sexp:     QSTRING		{ $$ = $1;				}
        | SVAR			{ $$ = $1->value.svar;			}
    	| UNDVAR '=' sexp	{ $$ = $3; $1->value.svar = $3;
				  $1->type = SVAR;			}
        | SVAR '=' sexp		{ $$ = $3; 
				  $1->value.svar = $3;
				  redefined_warning($1->name);          }
        | VAR '=' sexp		{ $$ = $3; 
				  $1->value.svar= $3;
				  redefined_warning($1->name);          
				  $1->type = SVAR; 		}
        | SFNCT '(' sexp ')'	{ $$ = (*($1->value.strfnct))($3);	}
	| SFNCT '(' ')'		{ $$ = (*($1->value.strfnct))();	}
        | SFNCT '(' exp  ')'	{ $$ = (*($1->value.strfnct))($3);	}
        | sexp CONCAT sexp	{ int len1 = strlen($1);
				  int len3 = strlen($3);
				  ALLOC($$, len1+len3+1, char *);
				  memcpy($$, $1, len1+1);
				  strcat($$, $3); }
        | SFNCT '(' exp ',' sexp ',' sexp ',' sexp ',' sexp ')'
				{ $$ = (*($1->value.strfnct))($3, $5, $7, $9, $11); }
        | SFNCT '(' exp ',' sexp ',' sexp  ')'
				{ $$ = (*($1->value.strfnct))($3, $5, $7); }
        | SFNCT '(' sexp ',' sexp ',' sexp  ')'
				{ $$ = (*($1->value.strfnct))($3, $5, $7); }
        | bool '?' sexp ':' sexp  { $$ = ($1) ? ($3) : ($5);              }

exp:	  NUM			{ $$ = $1; 				}
        | INC NUM		{ $$ = $2 + 1;				}
        | DEC NUM		{ $$ = $2 - 1;				}
        | VAR			{ $$ = $1->value.var;			}
	| INC VAR		{ $$ = ++($2->value.var);		}
	| DEC VAR		{ $$ = --($2->value.var);		}
	| VAR INC		{ $$ = ($1->value.var)++;		}
	| VAR DEC		{ $$ = ($1->value.var)--;		}
	| VAR '=' exp		{ $$ = $3; $1->value.var = $3;
				  redefined_warning($1->name);          }
	| SVAR '=' exp		{ $$ = $3; $1->value.var = $3;
				  redefined_warning($1->name);          
				  $1->type = VAR;			}
	| VAR EQ_PLUS exp	{ $1->value.var += $3; $$ = $1->value.var; }
	| VAR EQ_MINUS exp	{ $1->value.var -= $3; $$ = $1->value.var; }
	| VAR EQ_TIME exp	{ $1->value.var *= $3; $$ = $1->value.var; }
	| VAR EQ_DIV exp	{ $1->value.var /= $3; $$ = $1->value.var; }
	| VAR EQ_POW exp	{ errno = 0;
				  $1->value.var = pow($1->value.var,$3); 
				  $$ = $1->value.var; 
				  MATH_ERROR("Power");
				}
	| UNDVAR		{ $$ = $1->value.var;
				  undefined_warning($1->name);          }
	| INC UNDVAR		{ $$ = ++($2->value.var);		
				  $2->type = VAR;                       
				  undefined_warning($2->name);          }
	| DEC UNDVAR		{ $$ = --($2->value.var);		
				  $2->type = VAR;                       
				  undefined_warning($2->name);          }
	| UNDVAR INC		{ $$ = ($1->value.var)++;		
				  $1->type = VAR;                       
				  undefined_warning($1->name);          }
	| UNDVAR DEC		{ $$ = ($1->value.var)--;		
				  $1->type = VAR;                       
				  undefined_warning($1->name);          }
	| UNDVAR '=' exp	{ $$ = $3; $1->value.var = $3;
				  $1->type = VAR;                       }
	| UNDVAR EQ_PLUS exp	{ $1->value.var += $3; $$ = $1->value.var; 
				  $1->type = VAR;                       
				  undefined_warning($1->name);          }
	| UNDVAR EQ_MINUS exp	{ $1->value.var -= $3; $$ = $1->value.var; 
				  $1->type = VAR;                       
				  undefined_warning($1->name);          }
	| UNDVAR EQ_TIME exp	{ $1->value.var *= $3; $$ = $1->value.var; 
				  $1->type = VAR;                       
				  undefined_warning($1->name);          }
	| UNDVAR EQ_DIV exp	{ $1->value.var /= $3; $$ = $1->value.var; 
				  $1->type = VAR;                       
				  undefined_warning($1->name);          }
	| UNDVAR EQ_POW exp	{ errno = 0;
				  $1->value.var = pow($1->value.var,$3); 
				  $$ = $1->value.var; 
				  $1->type = VAR;                       
				  MATH_ERROR("Power");
				  undefined_warning($1->name);          }
	| FNCT '(' ')'		{ $$ = (*($1->value.fnctptr))();	}
	| FNCT '(' exp ')'	{ $$ = (*($1->value.fnctptr))($3); 	}
	| FNCT '(' sexp ')'	{ $$ = (*($1->value.fnctptr))($3); 	}
	| FNCT '(' sexp ',' exp ')'
                                { $$ = (*($1->value.fnctptr))($3, $5); 	}
	| FNCT '(' sexp ',' sexp ')'
                                { $$ = (*($1->value.fnctptr))($3, $5); 	}
   	| FNCT '(' exp ',' exp ')'
				{ $$ = (*($1->value.fnctptr))($3, $5);	}
	| FNCT '(' exp ',' exp ',' exp')'
				{ $$ = (*($1->value.fnctptr))($3, $5, $7); }
	| FNCT '(' exp ',' exp ';' exp ',' exp ')'
				{ $$ = (*($1->value.fnctptr))($3, $5, $7, $9); }
	| FNCT '(' exp ',' exp ',' exp ',' exp ')'
				{ $$ = (*($1->value.fnctptr))($3, $5, $7, $9); }
	| FNCT '(' exp ',' exp ',' exp ',' exp ',' exp ',' exp ')'
  		     { $$ = (*($1->value.fnctptr))($3, $5, $7, $9, $11, $13); }
	| exp '+' exp		{ $$ = $1 + $3; 			}
	| exp '-' exp		{ $$ = $1 - $3; 			}
	| exp '*' exp		{ $$ = $1 * $3; 			}
	| exp '/' exp		{ if ($3 == 0.)
				    {
				      yyerror("Zero divisor"); 
				      yyerrok;
				    }
				  else
				    $$ = $1 / $3; 			}
	| exp '%' exp		{ if ($3 == 0.)
				    {
				      yyerror("Zero divisor");
				      yyerrok;
				    }
				  else
				    $$ = (int)$1 % (int)$3;		}  
	| '-' exp %prec UNARY	{ $$ = -$2;				}
	| '+' exp %prec UNARY	{ $$ =  $2;				}
	| exp POW exp 		{ errno = 0;
				  $$ = pow($1, $3); 
				  MATH_ERROR("Power");			}
	| '(' exp ')'		{ $$ = $2;				}
    | '[' exp ']'           { errno = 0;
				  $$ = (double)($2 < 0 ? 
					-floor(-($2)): floor($2) );
				  MATH_ERROR("floor (int)");		}

    | bool { $$ = ($1) ? 1 : 0; }
    | bool '?' exp ':' exp  { $$ = ($1) ? ($3) : ($5);              }


/* End of grammar */
%%
# include "lex.yy.c"


