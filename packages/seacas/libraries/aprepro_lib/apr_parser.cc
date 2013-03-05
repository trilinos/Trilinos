/* A Bison parser, made by GNU Bison 2.4.2.  */

/* Skeleton implementation for Bison LALR(1) parsers in C++
   
      Copyright (C) 2002-2010 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

// Take the name prefix into account.
#define yylex   SEAMSlex

/* First part of user declarations.  */

/* Line 310 of lalr1.cc  */
#line 1 "aprepro.yy"

#include "aprepro.h"
#include "apr_util.h"

#include <stdlib.h>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <cstdio>

namespace SEAMS {
   extern int echo;
 }
 

/* Line 310 of lalr1.cc  */
#line 57 "apr_parser.cc"


#include "aprepro_parser.h"

/* User implementation prologue.  */

/* Line 316 of lalr1.cc  */
#line 79 "aprepro.yy"


#include "aprepro.h"
#include "apr_scanner.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex



/* Line 316 of lalr1.cc  */
#line 80 "apr_parser.cc"

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* FIXME: INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#define YYUSE(e) ((void) (e))

/* Enable debugging if requested.  */
#if YYDEBUG

/* A pseudo ostream that takes yydebug_ into account.  */
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)	\
do {							\
  if (yydebug_)						\
    {							\
      *yycdebug_ << Title << ' ';			\
      yy_symbol_print_ ((Type), (Value), (Location));	\
      *yycdebug_ << std::endl;				\
    }							\
} while (false)

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug_)				\
    yy_reduce_print_ (Rule);		\
} while (false)

# define YY_STACK_PRINT()		\
do {					\
  if (yydebug_)				\
    yystack_print_ ();			\
} while (false)

#else /* !YYDEBUG */

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_REDUCE_PRINT(Rule)
# define YY_STACK_PRINT()

#endif /* !YYDEBUG */

#define yyerrok		(yyerrstatus_ = 0)
#define yyclearin	(yychar = yyempty_)

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)


namespace SEAMS {

/* Line 379 of lalr1.cc  */
#line 146 "apr_parser.cc"
#if YYERROR_VERBOSE

  /* Return YYSTR after stripping away unnecessary quotes and
     backslashes, so that it's suitable for yyerror.  The heuristic is
     that double-quoting is unnecessary unless the string contains an
     apostrophe, a comma, or backslash (other than backslash-backslash).
     YYSTR is taken from yytname.  */
  std::string
  Parser::yytnamerr_ (const char *yystr)
  {
    if (*yystr == '"')
      {
        std::string yyr = "";
        char const *yyp = yystr;

        for (;;)
          switch (*++yyp)
            {
            case '\'':
            case ',':
              goto do_not_strip_quotes;

            case '\\':
              if (*++yyp != '\\')
                goto do_not_strip_quotes;
              /* Fall through.  */
            default:
              yyr += *yyp;
              break;

            case '"':
              return yyr;
            }
      do_not_strip_quotes: ;
      }

    return yystr;
  }

#endif

  /// Build a parser object.
  Parser::Parser (class Aprepro& aprepro_yyarg)
    :
#if YYDEBUG
      yydebug_ (false),
      yycdebug_ (&std::cerr),
#endif
      aprepro (aprepro_yyarg)
  {
  }

  Parser::~Parser ()
  {
  }

#if YYDEBUG
  /*--------------------------------.
  | Print this symbol on YYOUTPUT.  |
  `--------------------------------*/

  inline void
  Parser::yy_symbol_value_print_ (int yytype,
			   const semantic_type* yyvaluep, const location_type* yylocationp)
  {
    YYUSE (yylocationp);
    YYUSE (yyvaluep);
    switch (yytype)
      {
         default:
	  break;
      }
  }


  void
  Parser::yy_symbol_print_ (int yytype,
			   const semantic_type* yyvaluep, const location_type* yylocationp)
  {
    *yycdebug_ << (yytype < yyntokens_ ? "token" : "nterm")
	       << ' ' << yytname_[yytype] << " ("
	       << *yylocationp << ": ";
    yy_symbol_value_print_ (yytype, yyvaluep, yylocationp);
    *yycdebug_ << ')';
  }
#endif

  void
  Parser::yydestruct_ (const char* yymsg,
			   int yytype, semantic_type* yyvaluep, location_type* yylocationp)
  {
    YYUSE (yylocationp);
    YYUSE (yymsg);
    YYUSE (yyvaluep);

    YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

    switch (yytype)
      {
  
	default:
	  break;
      }
  }

  void
  Parser::yypop_ (unsigned int n)
  {
    yystate_stack_.pop (n);
    yysemantic_stack_.pop (n);
    yylocation_stack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
  Parser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  Parser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  Parser::debug_level_type
  Parser::debug_level () const
  {
    return yydebug_;
  }

  void
  Parser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif

  int
  Parser::parse ()
  {
    /// Lookahead and lookahead in internal form.
    int yychar = yyempty_;
    int yytoken = 0;

    /* State.  */
    int yyn;
    int yylen = 0;
    int yystate = 0;

    /* Error handling.  */
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// Semantic value of the lookahead.
    semantic_type yylval;
    /// Location of the lookahead.
    location_type yylloc;
    /// The locations where the error started and ended.
    location_type yyerror_range[2];

    /// $$.
    semantic_type yyval;
    /// @$.
    location_type yyloc;

    int yyresult;

    YYCDEBUG << "Starting parse" << std::endl;


    /* Initialize the stacks.  The initial state will be pushed in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystate_stack_ = state_stack_type (0);
    yysemantic_stack_ = semantic_stack_type (0);
    yylocation_stack_ = location_stack_type (0);
    yysemantic_stack_.push (yylval);
    yylocation_stack_.push (yylloc);

    /* New state.  */
  yynewstate:
    yystate_stack_.push (yystate);
    YYCDEBUG << "Entering state " << yystate << std::endl;

    /* Accept?  */
    if (yystate == yyfinal_)
      goto yyacceptlab;

    goto yybackup;

    /* Backup.  */
  yybackup:

    /* Try to take a decision without lookahead.  */
    yyn = yypact_[yystate];
    if (yyn == yypact_ninf_)
      goto yydefault;

    /* Read a lookahead token.  */
    if (yychar == yyempty_)
      {
	YYCDEBUG << "Reading a token: ";
	yychar = yylex (&yylval);
      }


    /* Convert token to internal form.  */
    if (yychar <= yyeof_)
      {
	yychar = yytoken = yyeof_;
	YYCDEBUG << "Now at end of input." << std::endl;
      }
    else
      {
	yytoken = yytranslate_ (yychar);
	YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
      }

    /* If the proper action on seeing token YYTOKEN is to reduce or to
       detect an error, take that action.  */
    yyn += yytoken;
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yytoken)
      goto yydefault;

    /* Reduce or error.  */
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
	if (yyn == 0 || yyn == yytable_ninf_)
	goto yyerrlab;
	yyn = -yyn;
	goto yyreduce;
      }

    /* Shift the lookahead token.  */
    YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

    /* Discard the token being shifted.  */
    yychar = yyempty_;

    yysemantic_stack_.push (yylval);
    yylocation_stack_.push (yylloc);

    /* Count tokens shifted since error; after three, turn off error
       status.  */
    if (yyerrstatus_)
      --yyerrstatus_;

    yystate = yyn;
    goto yynewstate;

  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[yystate];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;

  /*-----------------------------.
  | yyreduce -- Do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    /* If YYLEN is nonzero, implement the default value of the action:
       `$$ = $1'.  Otherwise, use the top of the stack.

       Otherwise, the following line sets YYVAL to garbage.
       This behavior is undocumented and Bison
       users should not rely upon it.  */
    if (yylen)
      yyval = yysemantic_stack_[yylen - 1];
    else
      yyval = yysemantic_stack_[0];

    {
      slice<location_type, location_stack_type> slice (yylocation_stack_, yylen);
      YYLLOC_DEFAULT (yyloc, slice, yylen);
    }
    YY_REDUCE_PRINT (yyn);
    switch (yyn)
      {
	  case 4:

/* Line 677 of lalr1.cc  */
#line 98 "aprepro.yy"
    { if (echo) aprepro.lexer->LexerOutput("\n", 1); }
    break;

  case 5:

/* Line 677 of lalr1.cc  */
#line 99 "aprepro.yy"
    { if (echo) {
	                             static char tmpstr[512];
				     SEAMS::symrec *format = aprepro.getsym("_FORMAT");
				     int len = sprintf(tmpstr, format->value.svar, (yysemantic_stack_[(3) - (2)].val));
				     aprepro.lexer->LexerOutput(tmpstr, len);
				   }
                                }
    break;

  case 6:

/* Line 677 of lalr1.cc  */
#line 106 "aprepro.yy"
    { if (echo && (yysemantic_stack_[(3) - (2)].string) != NULL) {
				    aprepro.lexer->LexerOutput((yysemantic_stack_[(3) - (2)].string), strlen((yysemantic_stack_[(3) - (2)].string)));
                                  }
                                }
    break;

  case 7:

/* Line 677 of lalr1.cc  */
#line 110 "aprepro.yy"
    {                                       }
    break;

  case 8:

/* Line 677 of lalr1.cc  */
#line 111 "aprepro.yy"
    { yyerrok;				}
    break;

  case 9:

/* Line 677 of lalr1.cc  */
#line 114 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) < (yysemantic_stack_[(3) - (3)].val);                         }
    break;

  case 10:

/* Line 677 of lalr1.cc  */
#line 115 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) > (yysemantic_stack_[(3) - (3)].val);                         }
    break;

  case 11:

/* Line 677 of lalr1.cc  */
#line 116 "aprepro.yy"
    { (yyval.val) = !((yysemantic_stack_[(2) - (2)].val));                           }
    break;

  case 12:

/* Line 677 of lalr1.cc  */
#line 117 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) <= (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 13:

/* Line 677 of lalr1.cc  */
#line 118 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) >= (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 14:

/* Line 677 of lalr1.cc  */
#line 119 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) == (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 15:

/* Line 677 of lalr1.cc  */
#line 120 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) != (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 16:

/* Line 677 of lalr1.cc  */
#line 121 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) || (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 17:

/* Line 677 of lalr1.cc  */
#line 122 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) && (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 18:

/* Line 677 of lalr1.cc  */
#line 123 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) || (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 19:

/* Line 677 of lalr1.cc  */
#line 124 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) && (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 20:

/* Line 677 of lalr1.cc  */
#line 125 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (2)].val);                              }
    break;

  case 21:

/* Line 677 of lalr1.cc  */
#line 128 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) <  0 ? 1 : 0);	}
    break;

  case 22:

/* Line 677 of lalr1.cc  */
#line 129 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) >  0 ? 1 : 0);	}
    break;

  case 23:

/* Line 677 of lalr1.cc  */
#line 130 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) <= 0 ? 1 : 0);	}
    break;

  case 24:

/* Line 677 of lalr1.cc  */
#line 131 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) >= 0 ? 1 : 0);	}
    break;

  case 25:

/* Line 677 of lalr1.cc  */
#line 132 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) == 0 ? 1 : 0);	}
    break;

  case 26:

/* Line 677 of lalr1.cc  */
#line 133 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) != 0 ? 1 : 0);	}
    break;

  case 27:

/* Line 677 of lalr1.cc  */
#line 135 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(1) - (1)].string);				}
    break;

  case 28:

/* Line 677 of lalr1.cc  */
#line 136 "aprepro.yy"
    { (yyval.string) = (char*)(yysemantic_stack_[(1) - (1)].tptr)->value.svar;			}
    break;

  case 29:

/* Line 677 of lalr1.cc  */
#line 137 "aprepro.yy"
    { (yyval.string) = (char*)(yysemantic_stack_[(1) - (1)].tptr)->value.svar;			}
    break;

  case 30:

/* Line 677 of lalr1.cc  */
#line 138 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); (yysemantic_stack_[(3) - (1)].tptr)->value.svar = (yysemantic_stack_[(3) - (3)].string);
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), Parser::token::SVAR);	}
    break;

  case 31:

/* Line 677 of lalr1.cc  */
#line 140 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); 
				  (yysemantic_stack_[(3) - (1)].tptr)->value.svar = (yysemantic_stack_[(3) - (3)].string);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 32:

/* Line 677 of lalr1.cc  */
#line 143 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); 
				  (yysemantic_stack_[(3) - (1)].tptr)->value.svar= (yysemantic_stack_[(3) - (3)].string);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::SVAR);		}
    break;

  case 33:

/* Line 677 of lalr1.cc  */
#line 147 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 34:

/* Line 677 of lalr1.cc  */
#line 148 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 35:

/* Line 677 of lalr1.cc  */
#line 149 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_c))((yysemantic_stack_[(4) - (3)].string));	}
    break;

  case 36:

/* Line 677 of lalr1.cc  */
#line 150 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(3) - (1)].tptr)->value.strfnct))();	}
    break;

  case 37:

/* Line 677 of lalr1.cc  */
#line 151 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_d))((yysemantic_stack_[(4) - (3)].val));	}
    break;

  case 38:

/* Line 677 of lalr1.cc  */
#line 152 "aprepro.yy"
    { concat_string((yysemantic_stack_[(3) - (1)].string), (yysemantic_stack_[(3) - (3)].string), &(yyval.string)); }
    break;

  case 39:

/* Line 677 of lalr1.cc  */
#line 154 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(12) - (1)].tptr)->value.strfnct_dcccc))((yysemantic_stack_[(12) - (3)].val), (yysemantic_stack_[(12) - (5)].string), (yysemantic_stack_[(12) - (7)].string), (yysemantic_stack_[(12) - (9)].string), (yysemantic_stack_[(12) - (11)].string)); }
    break;

  case 40:

/* Line 677 of lalr1.cc  */
#line 156 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(8) - (1)].tptr)->value.strfnct_dcc))((yysemantic_stack_[(8) - (3)].val), (yysemantic_stack_[(8) - (5)].string), (yysemantic_stack_[(8) - (7)].string)); }
    break;

  case 41:

/* Line 677 of lalr1.cc  */
#line 158 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(8) - (1)].tptr)->value.strfnct_ccc))((yysemantic_stack_[(8) - (3)].string), (yysemantic_stack_[(8) - (5)].string), (yysemantic_stack_[(8) - (7)].string)); }
    break;

  case 42:

/* Line 677 of lalr1.cc  */
#line 159 "aprepro.yy"
    { (yyval.string) = ((yysemantic_stack_[(5) - (1)].val)) ? ((yysemantic_stack_[(5) - (3)].string)) : ((yysemantic_stack_[(5) - (5)].string));              }
    break;

  case 43:

/* Line 677 of lalr1.cc  */
#line 161 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].val); 				}
    break;

  case 44:

/* Line 677 of lalr1.cc  */
#line 162 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(2) - (2)].val) + 1;				}
    break;

  case 45:

/* Line 677 of lalr1.cc  */
#line 163 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(2) - (2)].val) - 1;				}
    break;

  case 46:

/* Line 677 of lalr1.cc  */
#line 164 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;			}
    break;

  case 47:

/* Line 677 of lalr1.cc  */
#line 165 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;			}
    break;

  case 48:

/* Line 677 of lalr1.cc  */
#line 166 "aprepro.yy"
    { (yyval.val) = ++((yysemantic_stack_[(2) - (2)].tptr)->value.var);		}
    break;

  case 49:

/* Line 677 of lalr1.cc  */
#line 167 "aprepro.yy"
    { (yyval.val) = --((yysemantic_stack_[(2) - (2)].tptr)->value.var);		}
    break;

  case 50:

/* Line 677 of lalr1.cc  */
#line 168 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)++;		}
    break;

  case 51:

/* Line 677 of lalr1.cc  */
#line 169 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)--;		}
    break;

  case 52:

/* Line 677 of lalr1.cc  */
#line 170 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 53:

/* Line 677 of lalr1.cc  */
#line 172 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);			}
    break;

  case 54:

/* Line 677 of lalr1.cc  */
#line 175 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var += (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 55:

/* Line 677 of lalr1.cc  */
#line 176 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var -= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 56:

/* Line 677 of lalr1.cc  */
#line 177 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var *= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 57:

/* Line 677 of lalr1.cc  */
#line 178 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var /= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 58:

/* Line 677 of lalr1.cc  */
#line 179 "aprepro.yy"
    { errno = 0;
				  (yysemantic_stack_[(3) - (1)].tptr)->value.var = std::pow((yysemantic_stack_[(3) - (1)].tptr)->value.var,(yysemantic_stack_[(3) - (3)].val)); 
				  (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
				  SEAMS::math_error(aprepro, "Power");
				}
    break;

  case 59:

/* Line 677 of lalr1.cc  */
#line 184 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (2)].tptr)); YYERROR; }
    break;

  case 60:

/* Line 677 of lalr1.cc  */
#line 185 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (2)].tptr)); YYERROR; }
    break;

  case 61:

/* Line 677 of lalr1.cc  */
#line 186 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (1)].tptr)); YYERROR; }
    break;

  case 62:

/* Line 677 of lalr1.cc  */
#line 187 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (1)].tptr)); YYERROR; }
    break;

  case 63:

/* Line 677 of lalr1.cc  */
#line 188 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 64:

/* Line 677 of lalr1.cc  */
#line 189 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 65:

/* Line 677 of lalr1.cc  */
#line 190 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 66:

/* Line 677 of lalr1.cc  */
#line 191 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 67:

/* Line 677 of lalr1.cc  */
#line 192 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 68:

/* Line 677 of lalr1.cc  */
#line 193 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 69:

/* Line 677 of lalr1.cc  */
#line 194 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 70:

/* Line 677 of lalr1.cc  */
#line 196 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;
				  undefined_warning(aprepro, (yysemantic_stack_[(1) - (1)].tptr)->name);          }
    break;

  case 71:

/* Line 677 of lalr1.cc  */
#line 198 "aprepro.yy"
    { (yyval.val) = ++((yysemantic_stack_[(2) - (2)].tptr)->value.var);		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (2)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (2)].tptr)->name);          }
    break;

  case 72:

/* Line 677 of lalr1.cc  */
#line 201 "aprepro.yy"
    { (yyval.val) = --((yysemantic_stack_[(2) - (2)].tptr)->value.var);		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (2)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (2)].tptr)->name);          }
    break;

  case 73:

/* Line 677 of lalr1.cc  */
#line 204 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)++;		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (1)].tptr)->name);          }
    break;

  case 74:

/* Line 677 of lalr1.cc  */
#line 207 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)--;		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (1)].tptr)->name);          }
    break;

  case 75:

/* Line 677 of lalr1.cc  */
#line 210 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);                      }
    break;

  case 76:

/* Line 677 of lalr1.cc  */
#line 212 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var += (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 77:

/* Line 677 of lalr1.cc  */
#line 215 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var -= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 78:

/* Line 677 of lalr1.cc  */
#line 218 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var *= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 79:

/* Line 677 of lalr1.cc  */
#line 221 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var /= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 80:

/* Line 677 of lalr1.cc  */
#line 224 "aprepro.yy"
    { errno = 0;
				  (yysemantic_stack_[(3) - (1)].tptr)->value.var = std::pow((yysemantic_stack_[(3) - (1)].tptr)->value.var,(yysemantic_stack_[(3) - (3)].val)); 
				  (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  SEAMS::math_error(aprepro, "Power");
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 81:

/* Line 677 of lalr1.cc  */
#line 230 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(3) - (1)].tptr)->value.fnctptr))();	}
    break;

  case 82:

/* Line 677 of lalr1.cc  */
#line 231 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_d))((yysemantic_stack_[(4) - (3)].val)); 	}
    break;

  case 83:

/* Line 677 of lalr1.cc  */
#line 232 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_c))((yysemantic_stack_[(4) - (3)].string)); 	}
    break;

  case 84:

/* Line 677 of lalr1.cc  */
#line 234 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_cd))((yysemantic_stack_[(6) - (3)].string), (yysemantic_stack_[(6) - (5)].val)); 	}
    break;

  case 85:

/* Line 677 of lalr1.cc  */
#line 236 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_cc))((yysemantic_stack_[(6) - (3)].string), (yysemantic_stack_[(6) - (5)].string)); 	}
    break;

  case 86:

/* Line 677 of lalr1.cc  */
#line 238 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_dd))((yysemantic_stack_[(6) - (3)].val), (yysemantic_stack_[(6) - (5)].val));	}
    break;

  case 87:

/* Line 677 of lalr1.cc  */
#line 240 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(8) - (1)].tptr)->value.fnctptr_ddd))((yysemantic_stack_[(8) - (3)].val), (yysemantic_stack_[(8) - (5)].val), (yysemantic_stack_[(8) - (7)].val)); }
    break;

  case 88:

/* Line 677 of lalr1.cc  */
#line 242 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(10) - (1)].tptr)->value.fnctptr_dddd))((yysemantic_stack_[(10) - (3)].val), (yysemantic_stack_[(10) - (5)].val), (yysemantic_stack_[(10) - (7)].val), (yysemantic_stack_[(10) - (9)].val)); }
    break;

  case 89:

/* Line 677 of lalr1.cc  */
#line 244 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(10) - (1)].tptr)->value.fnctptr_dddd))((yysemantic_stack_[(10) - (3)].val), (yysemantic_stack_[(10) - (5)].val), (yysemantic_stack_[(10) - (7)].val), (yysemantic_stack_[(10) - (9)].val)); }
    break;

  case 90:

/* Line 677 of lalr1.cc  */
#line 246 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(14) - (1)].tptr)->value.fnctptr_dddddd))((yysemantic_stack_[(14) - (3)].val), (yysemantic_stack_[(14) - (5)].val), (yysemantic_stack_[(14) - (7)].val), (yysemantic_stack_[(14) - (9)].val), (yysemantic_stack_[(14) - (11)].val), (yysemantic_stack_[(14) - (13)].val)); }
    break;

  case 91:

/* Line 677 of lalr1.cc  */
#line 247 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) + (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 92:

/* Line 677 of lalr1.cc  */
#line 248 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) - (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 93:

/* Line 677 of lalr1.cc  */
#line 249 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) * (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 94:

/* Line 677 of lalr1.cc  */
#line 250 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (3)].val) == 0.)
				    {
				      yyerror(aprepro, "Zero divisor"); 
				      yyerrok;
				    }
				  else
				    (yyval.val) = (yysemantic_stack_[(3) - (1)].val) / (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 95:

/* Line 677 of lalr1.cc  */
#line 257 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (3)].val) == 0.)
				    {
				      yyerror(aprepro, "Zero divisor");
				      yyerrok;
				    }
				  else
				    (yyval.val) = (int)(yysemantic_stack_[(3) - (1)].val) % (int)(yysemantic_stack_[(3) - (3)].val);		}
    break;

  case 96:

/* Line 677 of lalr1.cc  */
#line 264 "aprepro.yy"
    { (yyval.val) = -(yysemantic_stack_[(2) - (2)].val);				}
    break;

  case 97:

/* Line 677 of lalr1.cc  */
#line 265 "aprepro.yy"
    { (yyval.val) =  (yysemantic_stack_[(2) - (2)].val);				}
    break;

  case 98:

/* Line 677 of lalr1.cc  */
#line 266 "aprepro.yy"
    { errno = 0;
				  (yyval.val) = std::pow((yysemantic_stack_[(3) - (1)].val), (yysemantic_stack_[(3) - (3)].val)); 
				  SEAMS::math_error(aprepro, "Power");			}
    break;

  case 99:

/* Line 677 of lalr1.cc  */
#line 269 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (2)].val);				}
    break;

  case 100:

/* Line 677 of lalr1.cc  */
#line 270 "aprepro.yy"
    { errno = 0;
				  (yyval.val) = (double)((yysemantic_stack_[(3) - (2)].val) < 0 ? -floor(-((yysemantic_stack_[(3) - (2)].val))): floor((yysemantic_stack_[(3) - (2)].val)) );
				  SEAMS::math_error(aprepro, "floor (int)");		}
    break;

  case 101:

/* Line 677 of lalr1.cc  */
#line 273 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(1) - (1)].val)) ? 1 : 0; }
    break;

  case 102:

/* Line 677 of lalr1.cc  */
#line 274 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(5) - (1)].val)) ? ((yysemantic_stack_[(5) - (3)].val)) : ((yysemantic_stack_[(5) - (5)].val));              }
    break;



/* Line 677 of lalr1.cc  */
#line 1191 "apr_parser.cc"
	default:
          break;
      }
    YY_SYMBOL_PRINT ("-> $$ =", yyr1_[yyn], &yyval, &yyloc);

    yypop_ (yylen);
    yylen = 0;
    YY_STACK_PRINT ();

    yysemantic_stack_.push (yyval);
    yylocation_stack_.push (yyloc);

    /* Shift the result of the reduction.  */
    yyn = yyr1_[yyn];
    yystate = yypgoto_[yyn - yyntokens_] + yystate_stack_[0];
    if (0 <= yystate && yystate <= yylast_
	&& yycheck_[yystate] == yystate_stack_[0])
      yystate = yytable_[yystate];
    else
      yystate = yydefgoto_[yyn - yyntokens_];
    goto yynewstate;

  /*------------------------------------.
  | yyerrlab -- here on detecting error |
  `------------------------------------*/
  yyerrlab:
    /* If not already recovering from an error, report this error.  */
    if (!yyerrstatus_)
      {
	++yynerrs_;
	error (yylloc, yysyntax_error_ (yystate, yytoken));
      }

    yyerror_range[0] = yylloc;
    if (yyerrstatus_ == 3)
      {
	/* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

	if (yychar <= yyeof_)
	  {
	  /* Return failure if at end of input.  */
	  if (yychar == yyeof_)
	    YYABORT;
	  }
	else
	  {
	    yydestruct_ ("Error: discarding", yytoken, &yylval, &yylloc);
	    yychar = yyempty_;
	  }
      }

    /* Else will try to reuse lookahead token after shifting the error
       token.  */
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:

    /* Pacify compilers like GCC when the user code never invokes
       YYERROR and the label yyerrorlab therefore never appears in user
       code.  */
    if (false)
      goto yyerrorlab;

    yyerror_range[0] = yylocation_stack_[yylen - 1];
    /* Do not reclaim the symbols of the rule which action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    yystate = yystate_stack_[0];
    goto yyerrlab1;

  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;	/* Each real token shifted decrements this.  */

    for (;;)
      {
	yyn = yypact_[yystate];
	if (yyn != yypact_ninf_)
	{
	  yyn += yyterror_;
	  if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == yyterror_)
	    {
	      yyn = yytable_[yyn];
	      if (0 < yyn)
		break;
	    }
	}

	/* Pop the current state because it cannot handle the error token.  */
	if (yystate_stack_.height () == 1)
	YYABORT;

	yyerror_range[0] = yylocation_stack_[0];
	yydestruct_ ("Error: popping",
		     yystos_[yystate],
		     &yysemantic_stack_[0], &yylocation_stack_[0]);
	yypop_ ();
	yystate = yystate_stack_[0];
	YY_STACK_PRINT ();
      }

    yyerror_range[1] = yylloc;
    // Using YYLLOC is tempting, but would change the location of
    // the lookahead.  YYLOC is available though.
    YYLLOC_DEFAULT (yyloc, (yyerror_range - 1), 2);
    yysemantic_stack_.push (yylval);
    yylocation_stack_.push (yyloc);

    /* Shift the error token.  */
    YY_SYMBOL_PRINT ("Shifting", yystos_[yyn],
		     &yysemantic_stack_[0], &yylocation_stack_[0]);

    yystate = yyn;
    goto yynewstate;

    /* Accept.  */
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;

    /* Abort.  */
  yyabortlab:
    yyresult = 1;
    goto yyreturn;

  yyreturn:
    if (yychar != yyempty_)
      yydestruct_ ("Cleanup: discarding lookahead", yytoken, &yylval, &yylloc);

    /* Do not reclaim the symbols of the rule which action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (yystate_stack_.height () != 1)
      {
	yydestruct_ ("Cleanup: popping",
		   yystos_[yystate_stack_[0]],
		   &yysemantic_stack_[0],
		   &yylocation_stack_[0]);
	yypop_ ();
      }

    return yyresult;
  }

  // Generate an error message.
  std::string
  Parser::yysyntax_error_ (int yystate, int tok)
  {
    std::string res;
    YYUSE (yystate);
#if YYERROR_VERBOSE
    int yyn = yypact_[yystate];
    if (yypact_ninf_ < yyn && yyn <= yylast_)
      {
	/* Start YYX at -YYN if negative to avoid negative indexes in
	   YYCHECK.  */
	int yyxbegin = yyn < 0 ? -yyn : 0;

	/* Stay within bounds of both yycheck and yytname.  */
	int yychecklim = yylast_ - yyn + 1;
	int yyxend = yychecklim < yyntokens_ ? yychecklim : yyntokens_;
	int count = 0;
	for (int x = yyxbegin; x < yyxend; ++x)
	  if (yycheck_[x + yyn] == x && x != yyterror_)
	    ++count;

	// FIXME: This method of building the message is not compatible
	// with internationalization.  It should work like yacc.c does it.
	// That is, first build a string that looks like this:
	// "syntax error, unexpected %s or %s or %s"
	// Then, invoke YY_ on this string.
	// Finally, use the string as a format to output
	// yytname_[tok], etc.
	// Until this gets fixed, this message appears in English only.
	res = "syntax error, unexpected ";
	res += yytnamerr_ (yytname_[tok]);
	if (count < 5)
	  {
	    count = 0;
	    for (int x = yyxbegin; x < yyxend; ++x)
	      if (yycheck_[x + yyn] == x && x != yyterror_)
		{
		  res += (!count++) ? ", expecting " : " or ";
		  res += yytnamerr_ (yytname_[x]);
		}
	  }
      }
    else
#endif
      res = YY_("syntax error");
    return res;
  }


  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
  const signed char Parser::yypact_ninf_ = -42;
  const short int
  Parser::yypact_[] =
  {
       -42,     4,   -42,   -13,   120,   -42,   -42,   -42,   -42,   -42,
     164,   330,    -9,   360,    -2,     6,     9,   277,   277,   -42,
     277,   277,   277,   101,   111,   -19,   784,   582,   277,   277,
     277,   277,   277,   277,   -42,   -42,   277,   277,   277,   277,
     277,   277,   -42,   -42,   277,   277,   277,   277,   277,   277,
     277,   -42,   -42,   277,   197,   263,   122,   376,   612,   641,
      -6,    -6,    -6,   -42,   -42,   -42,   -42,   -42,   -42,   -42,
     -42,   277,   277,   277,   277,   -42,   277,   277,   277,   277,
     277,   277,   277,   277,   -42,   277,   277,   277,   277,   277,
     277,   277,   277,   277,   277,   277,   277,   277,   376,   753,
     753,   753,   753,   753,   753,   376,   753,   753,   753,   753,
     753,   753,   376,   753,   376,   753,   753,   753,   753,   753,
     753,   376,   753,   -42,   183,   388,   -42,   327,   418,   -42,
     -42,   -42,   794,   753,   -28,   -41,   297,   -32,   -32,   -32,
     -32,   -32,   -32,   -42,   753,   768,   492,   214,   214,   214,
     214,   214,   128,   128,    -6,    -6,    -6,    -6,   277,   -42,
     277,   -42,   277,   -42,   277,   -42,   277,   277,    33,   671,
     358,    27,    66,   376,   753,   -42,   -42,   277,   -42,   277,
     277,   277,   448,   508,   145,   333,   277,   -42,   277,   -42,
     277,   -42,   478,   701,   228,   277,   -42,   -42,   277,   538,
     558,   277,   -42,   731,   -42
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned char
  Parser::yydefact_[] =
  {
         2,     0,     1,     0,     0,     4,     3,     8,    43,    27,
      70,    46,    28,    47,    29,     0,     0,     0,     0,     7,
       0,     0,     0,     0,     0,   101,     0,     0,     0,     0,
       0,     0,     0,     0,    74,    73,     0,     0,     0,     0,
       0,     0,    51,    50,     0,     0,     0,     0,     0,     0,
       0,    62,    61,     0,     0,     0,   101,     0,     0,     0,
      96,    97,    11,    45,    72,    49,    60,    44,    71,    48,
      59,     0,     0,     0,     0,     6,     0,     0,     0,     0,
       0,     0,     0,     0,     5,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    30,    75,
      77,    76,    79,    78,    80,    32,    52,    55,    54,    57,
      56,    58,    31,    53,    34,    63,    66,    65,    68,    67,
      69,    33,    64,    81,     0,     0,    36,     0,     0,    20,
      99,   100,     0,     0,    18,    19,     0,    26,    25,    24,
      23,    22,    21,    38,    10,    16,    17,    15,    14,    13,
      12,     9,    92,    91,    95,    93,    94,    98,     0,    83,
       0,    82,     0,    35,     0,    37,     0,     0,     0,     0,
       0,     0,     0,    42,   102,    85,    84,     0,    86,     0,
       0,     0,     0,     0,     0,     0,     0,    87,     0,    41,
       0,    40,     0,     0,     0,     0,    89,    88,     0,     0,
       0,     0,    39,     0,    90
  };

  /* YYPGOTO[NTERM-NUM].  */
  const signed char
  Parser::yypgoto_[] =
  {
       -42,   -42,   -42,   -16,    67,    -4
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  Parser::yydefgoto_[] =
  {
        -1,     1,     6,    25,    57,   133
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char Parser::yytable_ninf_ = -1;
  const unsigned char
  Parser::yytable_[] =
  {
        27,    56,    73,    80,     2,     3,     7,    83,    74,    71,
      72,    73,    44,    58,    59,    82,    60,    61,    62,    53,
      54,    74,     4,    55,    99,   100,   101,   102,   103,   104,
      74,     0,   106,   107,   108,   109,   110,   111,    97,   180,
     113,   115,   116,   117,   118,   119,   120,     0,   175,   122,
     125,   128,     5,     0,     0,     0,   134,   135,    76,    77,
      78,    79,    80,    81,    76,    77,    78,    79,    80,    81,
     136,    26,     0,     0,    82,     0,     0,     0,   181,   144,
      82,   145,   146,   147,   148,   149,   150,   151,   152,   153,
     154,   155,   156,   157,     0,    98,     0,    76,    77,    78,
      79,    80,    81,   105,    63,     0,    64,    65,     0,    66,
       0,   112,   114,    82,    67,     0,    68,    69,     0,    70,
     121,   124,   127,     8,     9,    10,    11,    12,    13,    14,
      15,    16,     0,     0,    17,     0,    18,   129,   132,    19,
       0,    83,     0,   137,   138,   139,   140,   141,   142,   143,
      71,    72,    73,     0,   169,     0,   170,    20,    21,     0,
     189,     0,    22,   174,     0,    23,    24,    94,    95,    96,
       0,    74,    97,   182,     0,   183,    76,    77,    78,    79,
      80,    81,   192,     0,   193,    28,    29,    30,    31,    32,
      33,   199,    82,     0,     0,   158,     0,   203,   159,     0,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    34,
      35,    17,   123,    18,    76,    77,    78,    79,    80,    81,
       0,     0,     0,     0,     0,   168,     0,    83,     0,   171,
      82,   172,     0,   173,    20,    21,     0,     0,     0,    22,
     198,     0,    23,    24,     0,     0,     0,   184,   185,     0,
       0,    92,    93,    94,    95,    96,     0,   194,    97,    76,
      77,    78,    79,    80,    81,   200,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    82,     0,    17,   126,    18,
       8,     9,    10,    11,    12,    13,    14,    15,    16,     0,
       0,    17,     0,    18,     0,     0,     0,     0,     0,     0,
      20,    21,     0,     0,     0,    22,     0,     0,    23,    24,
      83,     0,     0,     0,    20,    21,     0,     0,     0,    22,
       0,     0,    23,    24,     0,     0,    85,    86,    87,    88,
      89,    90,     0,    91,    92,    93,    94,    95,    96,   162,
       0,    97,   163,     0,     0,   190,     0,   167,   191,     0,
       0,    36,    37,    38,    39,    40,    41,     0,    76,    77,
      78,    79,    80,    81,    76,    77,    78,    79,    80,    81,
     177,    83,     0,   178,    82,    42,    43,     0,   179,     0,
      82,    45,    46,    47,    48,    49,    50,    85,    86,    87,
      88,    89,    90,     0,    91,    92,    93,    94,    95,    96,
     160,    83,    97,   161,     0,    51,    52,    76,    77,    78,
      79,    80,    81,     0,     0,     0,     0,    85,    86,    87,
      88,    89,    90,    82,    91,    92,    93,    94,    95,    96,
     164,    83,    97,   165,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    85,    86,    87,
      88,    89,    90,     0,    91,    92,    93,    94,    95,    96,
     186,    83,    97,   187,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    85,    86,    87,
      88,    89,    90,     0,    91,    92,    93,    94,    95,    96,
     195,    83,    97,   196,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    83,     0,    85,    86,    87,
      88,    89,    90,     0,    91,    92,    93,    94,    95,    96,
     188,    83,    97,    87,    88,    89,    90,     0,    91,    92,
      93,    94,    95,    96,     0,     0,    97,    85,    86,    87,
      88,    89,    90,     0,    91,    92,    93,    94,    95,    96,
     201,    83,    97,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    85,    86,    87,
      88,    89,    90,   202,    91,    92,    93,    94,    95,    96,
       0,     0,    97,     0,     0,     0,     0,     0,     0,    76,
      77,    78,    79,    80,    81,    83,     0,     0,     0,     0,
       0,    84,     0,     0,     0,    82,     0,     0,     0,     0,
       0,    85,    86,    87,    88,    89,    90,     0,    91,    92,
      93,    94,    95,    96,     0,    83,    97,   130,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    85,    86,    87,    88,    89,    90,     0,    91,    92,
      93,    94,    95,    96,    83,     0,    97,     0,   131,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      85,    86,    87,    88,    89,    90,     0,    91,    92,    93,
      94,    95,    96,     0,    83,    97,   176,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      85,    86,    87,    88,    89,    90,     0,    91,    92,    93,
      94,    95,    96,     0,    83,    97,   197,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      85,    86,    87,    88,    89,    90,     0,    91,    92,    93,
      94,    95,    96,     0,    83,    97,   204,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      85,    86,    87,    88,    89,    90,    83,    91,    92,    93,
      94,    95,    96,     0,     0,    97,     0,     0,     0,     0,
       0,    83,    85,    86,    87,    88,    89,    90,     0,    91,
      92,    93,    94,    95,    96,     0,     0,    97,    86,    87,
      88,    89,    90,    75,    91,    92,    93,    94,    95,    96,
       0,     0,    97,     0,     0,    76,    77,    78,    79,    80,
      81,   166,     0,     0,     0,    76,    77,    78,    79,    80,
      81,    82,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    82
  };

  /* YYCHECK.  */
  const short int
  Parser::yycheck_[] =
  {
         4,    17,    30,    35,     0,     1,    19,    13,    49,    28,
      29,    30,    21,    17,    18,    47,    20,    21,    22,    21,
      14,    49,    18,    14,    28,    29,    30,    31,    32,    33,
      49,    -1,    36,    37,    38,    39,    40,    41,    44,    12,
      44,    45,    46,    47,    48,    49,    50,    -1,    15,    53,
      54,    55,    48,    -1,    -1,    -1,    72,    73,    31,    32,
      33,    34,    35,    36,    31,    32,    33,    34,    35,    36,
      74,     4,    -1,    -1,    47,    -1,    -1,    -1,    12,    83,
      47,    85,    86,    87,    88,    89,    90,    91,    92,    93,
      94,    95,    96,    97,    -1,    28,    -1,    31,    32,    33,
      34,    35,    36,    36,     3,    -1,     5,     6,    -1,     8,
      -1,    44,    45,    47,     3,    -1,     5,     6,    -1,     8,
      53,    54,    55,     3,     4,     5,     6,     7,     8,     9,
      10,    11,    -1,    -1,    14,    -1,    16,    15,    71,    19,
      -1,    13,    -1,    76,    77,    78,    79,    80,    81,    82,
      28,    29,    30,    -1,   158,    -1,   160,    37,    38,    -1,
      15,    -1,    42,   167,    -1,    45,    46,    39,    40,    41,
      -1,    49,    44,   177,    -1,   179,    31,    32,    33,    34,
      35,    36,   186,    -1,   188,    21,    22,    23,    24,    25,
      26,   195,    47,    -1,    -1,    12,    -1,   201,    15,    -1,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    45,
      46,    14,    15,    16,    31,    32,    33,    34,    35,    36,
      -1,    -1,    -1,    -1,    -1,   158,    -1,    13,    -1,   162,
      47,   164,    -1,   166,    37,    38,    -1,    -1,    -1,    42,
      12,    -1,    45,    46,    -1,    -1,    -1,   180,   181,    -1,
      -1,    37,    38,    39,    40,    41,    -1,   190,    44,    31,
      32,    33,    34,    35,    36,   198,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    47,    -1,    14,    15,    16,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    -1,
      -1,    14,    -1,    16,    -1,    -1,    -1,    -1,    -1,    -1,
      37,    38,    -1,    -1,    -1,    42,    -1,    -1,    45,    46,
      13,    -1,    -1,    -1,    37,    38,    -1,    -1,    -1,    42,
      -1,    -1,    45,    46,    -1,    -1,    29,    30,    31,    32,
      33,    34,    -1,    36,    37,    38,    39,    40,    41,    12,
      -1,    44,    15,    -1,    -1,    12,    -1,    50,    15,    -1,
      -1,    21,    22,    23,    24,    25,    26,    -1,    31,    32,
      33,    34,    35,    36,    31,    32,    33,    34,    35,    36,
      12,    13,    -1,    15,    47,    45,    46,    -1,    20,    -1,
      47,    21,    22,    23,    24,    25,    26,    29,    30,    31,
      32,    33,    34,    -1,    36,    37,    38,    39,    40,    41,
      12,    13,    44,    15,    -1,    45,    46,    31,    32,    33,
      34,    35,    36,    -1,    -1,    -1,    -1,    29,    30,    31,
      32,    33,    34,    47,    36,    37,    38,    39,    40,    41,
      12,    13,    44,    15,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    29,    30,    31,
      32,    33,    34,    -1,    36,    37,    38,    39,    40,    41,
      12,    13,    44,    15,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    29,    30,    31,
      32,    33,    34,    -1,    36,    37,    38,    39,    40,    41,
      12,    13,    44,    15,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    13,    -1,    29,    30,    31,
      32,    33,    34,    -1,    36,    37,    38,    39,    40,    41,
      12,    13,    44,    31,    32,    33,    34,    -1,    36,    37,
      38,    39,    40,    41,    -1,    -1,    44,    29,    30,    31,
      32,    33,    34,    -1,    36,    37,    38,    39,    40,    41,
      12,    13,    44,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    29,    30,    31,
      32,    33,    34,    15,    36,    37,    38,    39,    40,    41,
      -1,    -1,    44,    -1,    -1,    -1,    -1,    -1,    -1,    31,
      32,    33,    34,    35,    36,    13,    -1,    -1,    -1,    -1,
      -1,    19,    -1,    -1,    -1,    47,    -1,    -1,    -1,    -1,
      -1,    29,    30,    31,    32,    33,    34,    -1,    36,    37,
      38,    39,    40,    41,    -1,    13,    44,    15,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    29,    30,    31,    32,    33,    34,    -1,    36,    37,
      38,    39,    40,    41,    13,    -1,    44,    -1,    17,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      29,    30,    31,    32,    33,    34,    -1,    36,    37,    38,
      39,    40,    41,    -1,    13,    44,    15,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      29,    30,    31,    32,    33,    34,    -1,    36,    37,    38,
      39,    40,    41,    -1,    13,    44,    15,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      29,    30,    31,    32,    33,    34,    -1,    36,    37,    38,
      39,    40,    41,    -1,    13,    44,    15,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      29,    30,    31,    32,    33,    34,    13,    36,    37,    38,
      39,    40,    41,    -1,    -1,    44,    -1,    -1,    -1,    -1,
      -1,    13,    29,    30,    31,    32,    33,    34,    -1,    36,
      37,    38,    39,    40,    41,    -1,    -1,    44,    30,    31,
      32,    33,    34,    19,    36,    37,    38,    39,    40,    41,
      -1,    -1,    44,    -1,    -1,    31,    32,    33,    34,    35,
      36,    27,    -1,    -1,    -1,    31,    32,    33,    34,    35,
      36,    47,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    47
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned char
  Parser::yystos_[] =
  {
         0,    52,     0,     1,    18,    48,    53,    19,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    14,    16,    19,
      37,    38,    42,    45,    46,    54,    55,    56,    21,    22,
      23,    24,    25,    26,    45,    46,    21,    22,    23,    24,
      25,    26,    45,    46,    21,    21,    22,    23,    24,    25,
      26,    45,    46,    21,    14,    14,    54,    55,    56,    56,
      56,    56,    56,     3,     5,     6,     8,     3,     5,     6,
       8,    28,    29,    30,    49,    19,    31,    32,    33,    34,
      35,    36,    47,    13,    19,    29,    30,    31,    32,    33,
      34,    36,    37,    38,    39,    40,    41,    44,    55,    56,
      56,    56,    56,    56,    56,    55,    56,    56,    56,    56,
      56,    56,    55,    56,    55,    56,    56,    56,    56,    56,
      56,    55,    56,    15,    55,    56,    15,    55,    56,    15,
      15,    17,    55,    56,    54,    54,    56,    55,    55,    55,
      55,    55,    55,    55,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    12,    15,
      12,    15,    12,    15,    12,    15,    27,    50,    55,    56,
      56,    55,    55,    55,    56,    15,    15,    12,    15,    20,
      12,    12,    56,    56,    55,    55,    12,    15,    12,    15,
      12,    15,    56,    56,    55,    12,    15,    15,    12,    56,
      55,    12,    15,    56,    15
  };

#if YYDEBUG
  /* TOKEN_NUMBER_[YYLEX-NUM] -- Internal symbol number corresponding
     to YYLEX-NUM.  */
  const unsigned short int
  Parser::yytoken_number_[] =
  {
         0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,    10,    63,
      58
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned char
  Parser::yyr1_[] =
  {
         0,    51,    52,    52,    53,    53,    53,    53,    53,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    55,    55,    55,
      55,    55,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    55,    55,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  Parser::yyr2_[] =
  {
         0,     2,     0,     2,     1,     3,     3,     2,     2,     3,
       3,     2,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     1,     1,
       3,     3,     3,     3,     3,     4,     3,     4,     3,    12,
       8,     8,     5,     1,     2,     2,     1,     1,     2,     2,
       2,     2,     3,     3,     3,     3,     3,     3,     3,     2,
       2,     2,     2,     3,     3,     3,     3,     3,     3,     3,
       1,     2,     2,     2,     2,     3,     3,     3,     3,     3,
       3,     3,     4,     4,     6,     6,     6,     8,    10,    10,
      14,     3,     3,     3,     3,     3,     2,     2,     3,     3,
       3,     1,     5
  };

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
  /* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
     First, the terminals, then, starting at \a yyntokens_, nonterminals.  */
  const char*
  const Parser::yytname_[] =
  {
    "\"end of file\"", "error", "$undefined", "NUM", "QSTRING", "UNDVAR",
  "VAR", "SVAR", "IMMVAR", "IMMSVAR", "FNCT", "SFNCT", "COMMA", "RT",
  "LPAR", "RPAR", "LBRACK", "RBRACK", "LBRACE", "RBRACE", "SEMI", "EQUAL",
  "EQ_MINUS", "EQ_PLUS", "EQ_DIV", "EQ_TIME", "EQ_POW", "COLON", "QUEST",
  "LOR", "LAND", "NE", "EQ", "GE", "LE", "GT", "LT", "SUB", "PLU", "MOD",
  "TIM", "DIV", "NOT", "UNARY", "POW", "DEC", "INC", "CONCAT", "'\\n'",
  "'?'", "':'", "$accept", "input", "line", "bool", "sexp", "exp", 0
  };
#endif

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const Parser::rhs_number_type
  Parser::yyrhs_[] =
  {
        52,     0,    -1,    -1,    52,    53,    -1,    48,    -1,    18,
      56,    19,    -1,    18,    55,    19,    -1,    18,    19,    -1,
       1,    19,    -1,    56,    36,    56,    -1,    56,    13,    56,
      -1,    42,    56,    -1,    56,    34,    56,    -1,    56,    33,
      56,    -1,    56,    32,    56,    -1,    56,    31,    56,    -1,
      56,    29,    56,    -1,    56,    30,    56,    -1,    54,    29,
      54,    -1,    54,    30,    54,    -1,    14,    54,    15,    -1,
      55,    36,    55,    -1,    55,    35,    55,    -1,    55,    34,
      55,    -1,    55,    33,    55,    -1,    55,    32,    55,    -1,
      55,    31,    55,    -1,     4,    -1,     7,    -1,     9,    -1,
       5,    21,    55,    -1,     7,    21,    55,    -1,     6,    21,
      55,    -1,     9,    21,    55,    -1,     8,    21,    55,    -1,
      11,    14,    55,    15,    -1,    11,    14,    15,    -1,    11,
      14,    56,    15,    -1,    55,    47,    55,    -1,    11,    14,
      56,    12,    55,    12,    55,    12,    55,    12,    55,    15,
      -1,    11,    14,    56,    12,    55,    12,    55,    15,    -1,
      11,    14,    55,    12,    55,    12,    55,    15,    -1,    54,
      28,    55,    27,    55,    -1,     3,    -1,    46,     3,    -1,
      45,     3,    -1,     6,    -1,     8,    -1,    46,     6,    -1,
      45,     6,    -1,     6,    46,    -1,     6,    45,    -1,     6,
      21,    56,    -1,     7,    21,    56,    -1,     6,    23,    56,
      -1,     6,    22,    56,    -1,     6,    25,    56,    -1,     6,
      24,    56,    -1,     6,    26,    56,    -1,    46,     8,    -1,
      45,     8,    -1,     8,    46,    -1,     8,    45,    -1,     8,
      21,    56,    -1,     9,    21,    56,    -1,     8,    23,    56,
      -1,     8,    22,    56,    -1,     8,    25,    56,    -1,     8,
      24,    56,    -1,     8,    26,    56,    -1,     5,    -1,    46,
       5,    -1,    45,     5,    -1,     5,    46,    -1,     5,    45,
      -1,     5,    21,    56,    -1,     5,    23,    56,    -1,     5,
      22,    56,    -1,     5,    25,    56,    -1,     5,    24,    56,
      -1,     5,    26,    56,    -1,    10,    14,    15,    -1,    10,
      14,    56,    15,    -1,    10,    14,    55,    15,    -1,    10,
      14,    55,    12,    56,    15,    -1,    10,    14,    55,    12,
      55,    15,    -1,    10,    14,    56,    12,    56,    15,    -1,
      10,    14,    56,    12,    56,    12,    56,    15,    -1,    10,
      14,    56,    12,    56,    20,    56,    12,    56,    15,    -1,
      10,    14,    56,    12,    56,    12,    56,    12,    56,    15,
      -1,    10,    14,    56,    12,    56,    12,    56,    12,    56,
      12,    56,    12,    56,    15,    -1,    56,    38,    56,    -1,
      56,    37,    56,    -1,    56,    40,    56,    -1,    56,    41,
      56,    -1,    56,    39,    56,    -1,    37,    56,    -1,    38,
      56,    -1,    56,    44,    56,    -1,    14,    56,    15,    -1,
      16,    56,    17,    -1,    54,    -1,    54,    49,    56,    50,
      56,    -1
  };

  /* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
     YYRHS.  */
  const unsigned short int
  Parser::yyprhs_[] =
  {
         0,     0,     3,     4,     7,     9,    13,    17,    20,    23,
      27,    31,    34,    38,    42,    46,    50,    54,    58,    62,
      66,    70,    74,    78,    82,    86,    90,    94,    96,    98,
     100,   104,   108,   112,   116,   120,   125,   129,   134,   138,
     151,   160,   169,   175,   177,   180,   183,   185,   187,   190,
     193,   196,   199,   203,   207,   211,   215,   219,   223,   227,
     230,   233,   236,   239,   243,   247,   251,   255,   259,   263,
     267,   269,   272,   275,   278,   281,   285,   289,   293,   297,
     301,   305,   309,   314,   319,   326,   333,   340,   349,   360,
     371,   386,   390,   394,   398,   402,   406,   409,   412,   416,
     420,   424,   426
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  Parser::yyrline_[] =
  {
         0,    94,    94,    95,    98,    99,   106,   110,   111,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   128,   129,   130,   131,   132,   133,   135,   136,   137,
     138,   140,   143,   147,   148,   149,   150,   151,   152,   153,
     155,   157,   159,   161,   162,   163,   164,   165,   166,   167,
     168,   169,   170,   172,   175,   176,   177,   178,   179,   184,
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     196,   198,   201,   204,   207,   210,   212,   215,   218,   221,
     224,   230,   231,   232,   233,   235,   237,   239,   241,   243,
     245,   247,   248,   249,   250,   257,   264,   265,   266,   269,
     270,   273,   274
  };

  // Print the state stack on the debug stream.
  void
  Parser::yystack_print_ ()
  {
    *yycdebug_ << "Stack now";
    for (state_stack_type::const_iterator i = yystate_stack_.begin ();
	 i != yystate_stack_.end (); ++i)
      *yycdebug_ << ' ' << *i;
    *yycdebug_ << std::endl;
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void
  Parser::yy_reduce_print_ (int yyrule)
  {
    unsigned int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    /* Print the symbols being reduced, and their result.  */
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
	       << " (line " << yylno << "):" << std::endl;
    /* The symbols being reduced.  */
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
		       yyrhs_[yyprhs_[yyrule] + yyi],
		       &(yysemantic_stack_[(yynrhs) - (yyi + 1)]),
		       &(yylocation_stack_[(yynrhs) - (yyi + 1)]));
  }
#endif // YYDEBUG

  /* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
  Parser::token_number_type
  Parser::yytranslate_ (int t)
  {
    static
    const token_number_type
    translate_table[] =
    {
           0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      48,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    50,     2,
       2,     2,     2,    49,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int Parser::yyeof_ = 0;
  const int Parser::yylast_ = 841;
  const int Parser::yynnts_ = 6;
  const int Parser::yyempty_ = -2;
  const int Parser::yyfinal_ = 2;
  const int Parser::yyterror_ = 1;
  const int Parser::yyerrcode_ = 256;
  const int Parser::yyntokens_ = 51;

  const unsigned int Parser::yyuser_token_number_max_ = 302;
  const Parser::token_number_type Parser::yyundef_token_ = 2;


} // SEAMS

/* Line 1053 of lalr1.cc  */
#line 1921 "apr_parser.cc"


/* Line 1055 of lalr1.cc  */
#line 278 "aprepro.yy"


void SEAMS::Parser::error(const Parser::location_type&, const std::string& m)
{
    aprepro.error(m);
}


