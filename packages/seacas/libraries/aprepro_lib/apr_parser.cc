/* A Bison parser, made by GNU Bison 2.7.12-4996.  */

/* Skeleton implementation for Bison LALR(1) parsers in C++
   
      Copyright (C) 2002-2013 Free Software Foundation, Inc.
   
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
/* Line 283 of lalr1.cc  */
#line 1 "aprepro.yy"

#include "aprepro.h"
#include "apr_util.h"
#include "apr_array.h"

#include <stdlib.h>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <cstdio>

namespace SEAMS {
   extern int echo;
 }
 
/* Line 283 of lalr1.cc  */
#line 56 "apr_parser.cc"


#include "aprepro_parser.h"

/* User implementation prologue.  */
/* Line 289 of lalr1.cc  */
#line 84 "aprepro.yy"


#include "aprepro.h"
#include "apr_scanner.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex


/* Line 289 of lalr1.cc  */
#line 77 "apr_parser.cc"


# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

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

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (/*CONSTCOND*/ false)
# endif


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
# define YY_SYMBOL_PRINT(Title, Type, Value, Location) YYUSE(Type)
# define YY_REDUCE_PRINT(Rule)        static_cast<void>(0)
# define YY_STACK_PRINT()             static_cast<void>(0)

#endif /* !YYDEBUG */

#define yyerrok		(yyerrstatus_ = 0)
#define yyclearin	(yychar = yyempty_)

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)


namespace SEAMS {
/* Line 357 of lalr1.cc  */
#line 172 "apr_parser.cc"

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
    std::ostream& yyo = debug_stream ();
    std::ostream& yyoutput = yyo;
    YYUSE (yyoutput);
    YYUSE (yytype);
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

    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

    YYUSE (yytype);
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

  inline bool
  Parser::yy_pact_value_is_default_ (int yyvalue)
  {
    return yyvalue == yypact_ninf_;
  }

  inline bool
  Parser::yy_table_value_is_error_ (int yyvalue)
  {
    return yyvalue == yytable_ninf_;
  }

  int
  Parser::parse ()
  {
    /// Lookahead and lookahead in internal form.
    int yychar = yyempty_;
    int yytoken = 0;

    // State.
    int yyn;
    int yylen = 0;
    int yystate = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// Semantic value of the lookahead.
    static semantic_type yyval_default;
    semantic_type yylval = yyval_default;
    /// Location of the lookahead.
    location_type yylloc;
    /// The locations where the error started and ended.
    location_type yyerror_range[3];

    /// $$.
    semantic_type yyval;
    /// @$.
    location_type yyloc;

    int yyresult;

    // FIXME: This shoud be completely indented.  It is not yet to
    // avoid gratuitous conflicts when merging into the master branch.
    try
      {
    YYCDEBUG << "Starting parse" << std::endl;


    /* Initialize the stacks.  The initial state will be pushed in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystate_stack_.clear ();
    yysemantic_stack_.clear ();
    yylocation_stack_.clear ();
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
    if (yy_pact_value_is_default_ (yyn))
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
	if (yy_table_value_is_error_ (yyn))
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

    // Compute the default @$.
    {
      slice<location_type, location_stack_type> slice (yylocation_stack_, yylen);
      YYLLOC_DEFAULT (yyloc, slice, yylen);
    }

    // Perform the reduction.
    YY_REDUCE_PRINT (yyn);
    switch (yyn)
      {
          case 4:
/* Line 664 of lalr1.cc  */
#line 103 "aprepro.yy"
    { if (echo) aprepro.lexer->LexerOutput("\n", 1); }
    break;

  case 5:
/* Line 664 of lalr1.cc  */
#line 104 "aprepro.yy"
    { if (echo) {
	                             static char tmpstr[512];
				     SEAMS::symrec *format = aprepro.getsym("_FORMAT");
				     int len = sprintf(tmpstr, format->value.svar, (yysemantic_stack_[(3) - (2)].val));
				     aprepro.lexer->LexerOutput(tmpstr, len);
				   }
                                }
    break;

  case 6:
/* Line 664 of lalr1.cc  */
#line 111 "aprepro.yy"
    { if (echo && (yysemantic_stack_[(3) - (2)].string) != NULL) {
				    aprepro.lexer->LexerOutput((yysemantic_stack_[(3) - (2)].string), strlen((yysemantic_stack_[(3) - (2)].string)));
                                  }
                                }
    break;

  case 7:
/* Line 664 of lalr1.cc  */
#line 115 "aprepro.yy"
    {                                       }
    break;

  case 8:
/* Line 664 of lalr1.cc  */
#line 116 "aprepro.yy"
    {                                       }
    break;

  case 9:
/* Line 664 of lalr1.cc  */
#line 117 "aprepro.yy"
    { yyerrok;				}
    break;

  case 10:
/* Line 664 of lalr1.cc  */
#line 120 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) < (yysemantic_stack_[(3) - (3)].val);                         }
    break;

  case 11:
/* Line 664 of lalr1.cc  */
#line 121 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) > (yysemantic_stack_[(3) - (3)].val);                         }
    break;

  case 12:
/* Line 664 of lalr1.cc  */
#line 122 "aprepro.yy"
    { (yyval.val) = !((yysemantic_stack_[(2) - (2)].val));                           }
    break;

  case 13:
/* Line 664 of lalr1.cc  */
#line 123 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) <= (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 14:
/* Line 664 of lalr1.cc  */
#line 124 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) >= (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 15:
/* Line 664 of lalr1.cc  */
#line 125 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) == (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 16:
/* Line 664 of lalr1.cc  */
#line 126 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) != (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 17:
/* Line 664 of lalr1.cc  */
#line 127 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) || (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 18:
/* Line 664 of lalr1.cc  */
#line 128 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) && (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 19:
/* Line 664 of lalr1.cc  */
#line 129 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) || (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 20:
/* Line 664 of lalr1.cc  */
#line 130 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) && (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 21:
/* Line 664 of lalr1.cc  */
#line 131 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (2)].val);                              }
    break;

  case 22:
/* Line 664 of lalr1.cc  */
#line 134 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) <  0 ? 1 : 0);	}
    break;

  case 23:
/* Line 664 of lalr1.cc  */
#line 135 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) >  0 ? 1 : 0);	}
    break;

  case 24:
/* Line 664 of lalr1.cc  */
#line 136 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) <= 0 ? 1 : 0);	}
    break;

  case 25:
/* Line 664 of lalr1.cc  */
#line 137 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) >= 0 ? 1 : 0);	}
    break;

  case 26:
/* Line 664 of lalr1.cc  */
#line 138 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) == 0 ? 1 : 0);	}
    break;

  case 27:
/* Line 664 of lalr1.cc  */
#line 139 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) != 0 ? 1 : 0);	}
    break;

  case 28:
/* Line 664 of lalr1.cc  */
#line 141 "aprepro.yy"
    { (yyval.arrval) = (yysemantic_stack_[(1) - (1)].tptr)->value.avar;}
    break;

  case 29:
/* Line 664 of lalr1.cc  */
#line 142 "aprepro.yy"
    { (yyval.arrval) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.arrfnct_c))((yysemantic_stack_[(4) - (3)].string));      }
    break;

  case 30:
/* Line 664 of lalr1.cc  */
#line 144 "aprepro.yy"
    { (yyval.arrval) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.arrfnct_dd))((yysemantic_stack_[(6) - (3)].val),(yysemantic_stack_[(6) - (5)].val));   }
    break;

  case 31:
/* Line 664 of lalr1.cc  */
#line 146 "aprepro.yy"
    { (yyval.arrval) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.arrfnct_d))((yysemantic_stack_[(4) - (3)].val));      }
    break;

  case 32:
/* Line 664 of lalr1.cc  */
#line 148 "aprepro.yy"
    { (yyval.arrval) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.arrfnct_a))((yysemantic_stack_[(4) - (3)].arrval));      }
    break;

  case 33:
/* Line 664 of lalr1.cc  */
#line 149 "aprepro.yy"
    { (yyval.arrval) = (yysemantic_stack_[(3) - (3)].arrval); (yysemantic_stack_[(3) - (1)].tptr)->value.avar = (yysemantic_stack_[(3) - (3)].arrval); 
                                  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr));
                                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::AVAR); }
    break;

  case 34:
/* Line 664 of lalr1.cc  */
#line 152 "aprepro.yy"
    { (yyval.arrval) = (yysemantic_stack_[(3) - (3)].arrval); (yysemantic_stack_[(3) - (1)].tptr)->value.avar = (yysemantic_stack_[(3) - (3)].arrval); 
                                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::AVAR); }
    break;

  case 35:
/* Line 664 of lalr1.cc  */
#line 154 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (1)].arrval)->cols == (yysemantic_stack_[(3) - (3)].arrval)->cols && (yysemantic_stack_[(3) - (1)].arrval)->rows == (yysemantic_stack_[(3) - (3)].arrval)->rows ) {
                                     (yyval.arrval) = array_add((yysemantic_stack_[(3) - (1)].arrval), (yysemantic_stack_[(3) - (3)].arrval)); 
                                  }
                                  else {
                                    yyerror(aprepro, "Arrays do not have same row and column count"); 
                                    yyerrok;
                                  }
                                }
    break;

  case 36:
/* Line 664 of lalr1.cc  */
#line 162 "aprepro.yy"
    { (yyval.arrval) = array_scale((yysemantic_stack_[(2) - (2)].arrval), -1.0);           }
    break;

  case 37:
/* Line 664 of lalr1.cc  */
#line 164 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (1)].arrval)->cols == (yysemantic_stack_[(3) - (3)].arrval)->cols && (yysemantic_stack_[(3) - (1)].arrval)->rows == (yysemantic_stack_[(3) - (3)].arrval)->rows ) {
                                     (yyval.arrval) = array_sub((yysemantic_stack_[(3) - (1)].arrval), (yysemantic_stack_[(3) - (3)].arrval)); 
                                  }
                                  else {
                                    yyerror(aprepro, "Arrays do not have same row and column count"); 
                                    yyerrok;
                                  }
                                }
    break;

  case 38:
/* Line 664 of lalr1.cc  */
#line 172 "aprepro.yy"
    { (yyval.arrval) = array_scale((yysemantic_stack_[(3) - (1)].arrval), (yysemantic_stack_[(3) - (3)].val));             }
    break;

  case 39:
/* Line 664 of lalr1.cc  */
#line 173 "aprepro.yy"
    { (yyval.arrval) = array_scale((yysemantic_stack_[(3) - (1)].arrval), 1.0/(yysemantic_stack_[(3) - (3)].val));         }
    break;

  case 40:
/* Line 664 of lalr1.cc  */
#line 174 "aprepro.yy"
    { (yyval.arrval) = array_scale((yysemantic_stack_[(3) - (3)].arrval), (yysemantic_stack_[(3) - (1)].val));             }
    break;

  case 41:
/* Line 664 of lalr1.cc  */
#line 175 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (1)].arrval)->cols == (yysemantic_stack_[(3) - (3)].arrval)->rows) {
                                    (yyval.arrval) = array_mult((yysemantic_stack_[(3) - (1)].arrval), (yysemantic_stack_[(3) - (3)].arrval));
                                  }
                                  else {
                                    yyerror(aprepro, "Column count of first array does not match row count of second array"); 
                                    yyerrok;
                                  }
				}
    break;

  case 42:
/* Line 664 of lalr1.cc  */
#line 184 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(1) - (1)].string);				}
    break;

  case 43:
/* Line 664 of lalr1.cc  */
#line 185 "aprepro.yy"
    { (yyval.string) = (char*)(yysemantic_stack_[(1) - (1)].tptr)->value.svar;			}
    break;

  case 44:
/* Line 664 of lalr1.cc  */
#line 186 "aprepro.yy"
    { (yyval.string) = (char*)(yysemantic_stack_[(1) - (1)].tptr)->value.svar;			}
    break;

  case 45:
/* Line 664 of lalr1.cc  */
#line 187 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); (yysemantic_stack_[(3) - (1)].tptr)->value.svar = (yysemantic_stack_[(3) - (3)].string);
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), Parser::token::SVAR);	}
    break;

  case 46:
/* Line 664 of lalr1.cc  */
#line 189 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); 
				  (yysemantic_stack_[(3) - (1)].tptr)->value.svar = (yysemantic_stack_[(3) - (3)].string);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr));          }
    break;

  case 47:
/* Line 664 of lalr1.cc  */
#line 192 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); 
				  (yysemantic_stack_[(3) - (1)].tptr)->value.svar= (yysemantic_stack_[(3) - (3)].string);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr));          
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::SVAR);		}
    break;

  case 48:
/* Line 664 of lalr1.cc  */
#line 196 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 49:
/* Line 664 of lalr1.cc  */
#line 197 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 50:
/* Line 664 of lalr1.cc  */
#line 198 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_c))((yysemantic_stack_[(4) - (3)].string));	}
    break;

  case 51:
/* Line 664 of lalr1.cc  */
#line 199 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(3) - (1)].tptr)->value.strfnct))();	}
    break;

  case 52:
/* Line 664 of lalr1.cc  */
#line 200 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_d))((yysemantic_stack_[(4) - (3)].val));	}
    break;

  case 53:
/* Line 664 of lalr1.cc  */
#line 201 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_a))((yysemantic_stack_[(4) - (3)].arrval));	}
    break;

  case 54:
/* Line 664 of lalr1.cc  */
#line 202 "aprepro.yy"
    { concat_string((yysemantic_stack_[(3) - (1)].string), (yysemantic_stack_[(3) - (3)].string), &(yyval.string)); }
    break;

  case 55:
/* Line 664 of lalr1.cc  */
#line 204 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(12) - (1)].tptr)->value.strfnct_dcccc))((yysemantic_stack_[(12) - (3)].val), (yysemantic_stack_[(12) - (5)].string), (yysemantic_stack_[(12) - (7)].string), (yysemantic_stack_[(12) - (9)].string), (yysemantic_stack_[(12) - (11)].string)); }
    break;

  case 56:
/* Line 664 of lalr1.cc  */
#line 206 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(8) - (1)].tptr)->value.strfnct_dcc))((yysemantic_stack_[(8) - (3)].val), (yysemantic_stack_[(8) - (5)].string), (yysemantic_stack_[(8) - (7)].string)); }
    break;

  case 57:
/* Line 664 of lalr1.cc  */
#line 208 "aprepro.yy"
    { (yyval.string) = (char*)(*((yysemantic_stack_[(8) - (1)].tptr)->value.strfnct_ccc))((yysemantic_stack_[(8) - (3)].string), (yysemantic_stack_[(8) - (5)].string), (yysemantic_stack_[(8) - (7)].string)); }
    break;

  case 58:
/* Line 664 of lalr1.cc  */
#line 209 "aprepro.yy"
    { (yyval.string) = ((yysemantic_stack_[(5) - (1)].val)) ? ((yysemantic_stack_[(5) - (3)].string)) : ((yysemantic_stack_[(5) - (5)].string));              }
    break;

  case 59:
/* Line 664 of lalr1.cc  */
#line 211 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].val); 				}
    break;

  case 60:
/* Line 664 of lalr1.cc  */
#line 212 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(2) - (2)].val) + 1;				}
    break;

  case 61:
/* Line 664 of lalr1.cc  */
#line 213 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(2) - (2)].val) - 1;				}
    break;

  case 62:
/* Line 664 of lalr1.cc  */
#line 214 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;			}
    break;

  case 63:
/* Line 664 of lalr1.cc  */
#line 215 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;			}
    break;

  case 64:
/* Line 664 of lalr1.cc  */
#line 216 "aprepro.yy"
    { (yyval.val) = ++((yysemantic_stack_[(2) - (2)].tptr)->value.var);		}
    break;

  case 65:
/* Line 664 of lalr1.cc  */
#line 217 "aprepro.yy"
    { (yyval.val) = --((yysemantic_stack_[(2) - (2)].tptr)->value.var);		}
    break;

  case 66:
/* Line 664 of lalr1.cc  */
#line 218 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)++;		}
    break;

  case 67:
/* Line 664 of lalr1.cc  */
#line 219 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)--;		}
    break;

  case 68:
/* Line 664 of lalr1.cc  */
#line 220 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr));          }
    break;

  case 69:
/* Line 664 of lalr1.cc  */
#line 222 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr));          
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);			}
    break;

  case 70:
/* Line 664 of lalr1.cc  */
#line 225 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var += (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 71:
/* Line 664 of lalr1.cc  */
#line 226 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var -= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 72:
/* Line 664 of lalr1.cc  */
#line 227 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var *= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 73:
/* Line 664 of lalr1.cc  */
#line 228 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var /= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 74:
/* Line 664 of lalr1.cc  */
#line 229 "aprepro.yy"
    { errno = 0;
				  (yysemantic_stack_[(3) - (1)].tptr)->value.var = std::pow((yysemantic_stack_[(3) - (1)].tptr)->value.var,(yysemantic_stack_[(3) - (3)].val)); 
				  (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
				  SEAMS::math_error(aprepro, "Power");
				}
    break;

  case 75:
/* Line 664 of lalr1.cc  */
#line 234 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (2)].tptr)); YYERROR; }
    break;

  case 76:
/* Line 664 of lalr1.cc  */
#line 235 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (2)].tptr)); YYERROR; }
    break;

  case 77:
/* Line 664 of lalr1.cc  */
#line 236 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (1)].tptr)); YYERROR; }
    break;

  case 78:
/* Line 664 of lalr1.cc  */
#line 237 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (1)].tptr)); YYERROR; }
    break;

  case 79:
/* Line 664 of lalr1.cc  */
#line 238 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 80:
/* Line 664 of lalr1.cc  */
#line 239 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 81:
/* Line 664 of lalr1.cc  */
#line 240 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 82:
/* Line 664 of lalr1.cc  */
#line 241 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 83:
/* Line 664 of lalr1.cc  */
#line 242 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 84:
/* Line 664 of lalr1.cc  */
#line 243 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 85:
/* Line 664 of lalr1.cc  */
#line 244 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 86:
/* Line 664 of lalr1.cc  */
#line 246 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;
				  undefined_warning(aprepro, (yysemantic_stack_[(1) - (1)].tptr)->name);          }
    break;

  case 87:
/* Line 664 of lalr1.cc  */
#line 248 "aprepro.yy"
    { (yyval.val) = ++((yysemantic_stack_[(2) - (2)].tptr)->value.var);		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (2)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (2)].tptr)->name);          }
    break;

  case 88:
/* Line 664 of lalr1.cc  */
#line 251 "aprepro.yy"
    { (yyval.val) = --((yysemantic_stack_[(2) - (2)].tptr)->value.var);		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (2)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (2)].tptr)->name);          }
    break;

  case 89:
/* Line 664 of lalr1.cc  */
#line 254 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)++;		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (1)].tptr)->name);          }
    break;

  case 90:
/* Line 664 of lalr1.cc  */
#line 257 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)--;		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (1)].tptr)->name);          }
    break;

  case 91:
/* Line 664 of lalr1.cc  */
#line 260 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);                      }
    break;

  case 92:
/* Line 664 of lalr1.cc  */
#line 262 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var += (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 93:
/* Line 664 of lalr1.cc  */
#line 265 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var -= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 94:
/* Line 664 of lalr1.cc  */
#line 268 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var *= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 95:
/* Line 664 of lalr1.cc  */
#line 271 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var /= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 96:
/* Line 664 of lalr1.cc  */
#line 274 "aprepro.yy"
    { errno = 0;
				  (yysemantic_stack_[(3) - (1)].tptr)->value.var = std::pow((yysemantic_stack_[(3) - (1)].tptr)->value.var,(yysemantic_stack_[(3) - (3)].val)); 
				  (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  SEAMS::math_error(aprepro, "Power");
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 97:
/* Line 664 of lalr1.cc  */
#line 280 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(3) - (1)].tptr)->value.fnctptr))();	}
    break;

  case 98:
/* Line 664 of lalr1.cc  */
#line 281 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_d))((yysemantic_stack_[(4) - (3)].val)); 	}
    break;

  case 99:
/* Line 664 of lalr1.cc  */
#line 282 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_c))((yysemantic_stack_[(4) - (3)].string)); 	}
    break;

  case 100:
/* Line 664 of lalr1.cc  */
#line 283 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_a))((yysemantic_stack_[(4) - (3)].arrval)); 	}
    break;

  case 101:
/* Line 664 of lalr1.cc  */
#line 285 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_cd))((yysemantic_stack_[(6) - (3)].string), (yysemantic_stack_[(6) - (5)].val)); 	}
    break;

  case 102:
/* Line 664 of lalr1.cc  */
#line 287 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_cc))((yysemantic_stack_[(6) - (3)].string), (yysemantic_stack_[(6) - (5)].string)); 	}
    break;

  case 103:
/* Line 664 of lalr1.cc  */
#line 289 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_dd))((yysemantic_stack_[(6) - (3)].val), (yysemantic_stack_[(6) - (5)].val));	}
    break;

  case 104:
/* Line 664 of lalr1.cc  */
#line 291 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(8) - (1)].tptr)->value.fnctptr_ddd))((yysemantic_stack_[(8) - (3)].val), (yysemantic_stack_[(8) - (5)].val), (yysemantic_stack_[(8) - (7)].val)); }
    break;

  case 105:
/* Line 664 of lalr1.cc  */
#line 293 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(10) - (1)].tptr)->value.fnctptr_dddd))((yysemantic_stack_[(10) - (3)].val), (yysemantic_stack_[(10) - (5)].val), (yysemantic_stack_[(10) - (7)].val), (yysemantic_stack_[(10) - (9)].val)); }
    break;

  case 106:
/* Line 664 of lalr1.cc  */
#line 295 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(10) - (1)].tptr)->value.fnctptr_dddd))((yysemantic_stack_[(10) - (3)].val), (yysemantic_stack_[(10) - (5)].val), (yysemantic_stack_[(10) - (7)].val), (yysemantic_stack_[(10) - (9)].val)); }
    break;

  case 107:
/* Line 664 of lalr1.cc  */
#line 297 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(14) - (1)].tptr)->value.fnctptr_dddddd))((yysemantic_stack_[(14) - (3)].val), (yysemantic_stack_[(14) - (5)].val), (yysemantic_stack_[(14) - (7)].val), (yysemantic_stack_[(14) - (9)].val), (yysemantic_stack_[(14) - (11)].val), (yysemantic_stack_[(14) - (13)].val)); }
    break;

  case 108:
/* Line 664 of lalr1.cc  */
#line 298 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) + (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 109:
/* Line 664 of lalr1.cc  */
#line 299 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) - (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 110:
/* Line 664 of lalr1.cc  */
#line 300 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) * (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 111:
/* Line 664 of lalr1.cc  */
#line 301 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (3)].val) == 0.)
				    {
				      yyerror(aprepro, "Zero divisor"); 
				      yyerrok;
				    }
				  else
				    (yyval.val) = (yysemantic_stack_[(3) - (1)].val) / (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 112:
/* Line 664 of lalr1.cc  */
#line 308 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (3)].val) == 0.)
				    {
				      yyerror(aprepro, "Zero divisor");
				      yyerrok;
				    }
				  else
				    (yyval.val) = (int)(yysemantic_stack_[(3) - (1)].val) % (int)(yysemantic_stack_[(3) - (3)].val);		}
    break;

  case 113:
/* Line 664 of lalr1.cc  */
#line 315 "aprepro.yy"
    { (yyval.val) = -(yysemantic_stack_[(2) - (2)].val);				}
    break;

  case 114:
/* Line 664 of lalr1.cc  */
#line 316 "aprepro.yy"
    { (yyval.val) =  (yysemantic_stack_[(2) - (2)].val);				}
    break;

  case 115:
/* Line 664 of lalr1.cc  */
#line 317 "aprepro.yy"
    { errno = 0;
				  (yyval.val) = std::pow((yysemantic_stack_[(3) - (1)].val), (yysemantic_stack_[(3) - (3)].val)); 
				  SEAMS::math_error(aprepro, "Power");			}
    break;

  case 116:
/* Line 664 of lalr1.cc  */
#line 320 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (2)].val);				}
    break;

  case 117:
/* Line 664 of lalr1.cc  */
#line 321 "aprepro.yy"
    { errno = 0;
				  (yyval.val) = (double)((yysemantic_stack_[(3) - (2)].val) < 0 ? -floor(-((yysemantic_stack_[(3) - (2)].val))): floor((yysemantic_stack_[(3) - (2)].val)) );
				  SEAMS::math_error(aprepro, "floor (int)");		}
    break;

  case 118:
/* Line 664 of lalr1.cc  */
#line 324 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(1) - (1)].val)) ? 1 : 0; }
    break;

  case 119:
/* Line 664 of lalr1.cc  */
#line 325 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(5) - (1)].val)) ? ((yysemantic_stack_[(5) - (3)].val)) : ((yysemantic_stack_[(5) - (5)].val));              }
    break;

  case 120:
/* Line 664 of lalr1.cc  */
#line 326 "aprepro.yy"
    { array *arr = (yysemantic_stack_[(6) - (1)].tptr)->value.avar;
                                      int cols = arr->cols;
                                      int rows = arr->rows;
                                      if ((yysemantic_stack_[(6) - (3)].val) < rows && (yysemantic_stack_[(6) - (5)].val) < cols) {
                                        int offset = (yysemantic_stack_[(6) - (3)].val)*cols+(yysemantic_stack_[(6) - (5)].val);
                                        (yyval.val) = (yysemantic_stack_[(6) - (1)].tptr)->value.avar->data[offset];
                                      }
                                      else {
                                        yyerror(aprepro, "Row or Column index out of range"); 
                                        yyerrok;
                                      }
                                    }
    break;

  case 121:
/* Line 664 of lalr1.cc  */
#line 339 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(8) - (8)].val);
                                    array *arr = (yysemantic_stack_[(8) - (1)].tptr)->value.avar;
                                    int cols = arr->cols;
                                    int rows = arr->rows;
                                    if ((yysemantic_stack_[(8) - (3)].val) < rows && (yysemantic_stack_[(8) - (5)].val) < cols) {
                                      int offset = (yysemantic_stack_[(8) - (3)].val)*cols+(yysemantic_stack_[(8) - (5)].val);
                                      (yysemantic_stack_[(8) - (1)].tptr)->value.avar->data[offset] = (yysemantic_stack_[(8) - (8)].val);
                                    }
                                    else {
                                      yyerror(aprepro, "Row or Column index out of range"); 
                                      yyerrok;
                                    }
                                  }
    break;


/* Line 664 of lalr1.cc  */
#line 1290 "apr_parser.cc"
      default:
        break;
      }

    /* User semantic actions sometimes alter yychar, and that requires
       that yytoken be updated with the new translation.  We take the
       approach of translating immediately before every use of yytoken.
       One alternative is translating here after every semantic action,
       but that translation would be missed if the semantic action
       invokes YYABORT, YYACCEPT, or YYERROR immediately after altering
       yychar.  In the case of YYABORT or YYACCEPT, an incorrect
       destructor might then be invoked immediately.  In the case of
       YYERROR, subsequent parser actions might lead to an incorrect
       destructor call or verbose syntax error message before the
       lookahead is translated.  */
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
    /* Make sure we have latest lookahead translation.  See comments at
       user semantic actions for why this is necessary.  */
    yytoken = yytranslate_ (yychar);

    /* If not already recovering from an error, report this error.  */
    if (!yyerrstatus_)
      {
	++yynerrs_;
	if (yychar == yyempty_)
	  yytoken = yyempty_;
	error (yylloc, yysyntax_error_ (yystate, yytoken));
      }

    yyerror_range[1] = yylloc;
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

    yyerror_range[1] = yylocation_stack_[yylen - 1];
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
	if (!yy_pact_value_is_default_ (yyn))
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

	yyerror_range[1] = yylocation_stack_[0];
	yydestruct_ ("Error: popping",
		     yystos_[yystate],
		     &yysemantic_stack_[0], &yylocation_stack_[0]);
	yypop_ ();
	yystate = yystate_stack_[0];
	YY_STACK_PRINT ();
      }

    yyerror_range[2] = yylloc;
    // Using YYLLOC is tempting, but would change the location of
    // the lookahead.  YYLOC is available though.
    YYLLOC_DEFAULT (yyloc, yyerror_range, 2);
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
      {
        /* Make sure we have latest lookahead translation.  See comments
           at user semantic actions for why this is necessary.  */
        yytoken = yytranslate_ (yychar);
        yydestruct_ ("Cleanup: discarding lookahead", yytoken, &yylval,
                     &yylloc);
      }

    /* Do not reclaim the symbols of the rule which action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (1 < yystate_stack_.height ())
      {
        yydestruct_ ("Cleanup: popping",
                     yystos_[yystate_stack_[0]],
                     &yysemantic_stack_[0],
                     &yylocation_stack_[0]);
        yypop_ ();
      }

    return yyresult;
    }
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack"
                 << std::endl;
        // Do not try to display the values of the reclaimed symbols,
        // as their printer might throw an exception.
        if (yychar != yyempty_)
          {
            /* Make sure we have latest lookahead translation.  See
               comments at user semantic actions for why this is
               necessary.  */
            yytoken = yytranslate_ (yychar);
            yydestruct_ (YY_NULL, yytoken, &yylval, &yylloc);
          }

        while (1 < yystate_stack_.height ())
          {
            yydestruct_ (YY_NULL,
                         yystos_[yystate_stack_[0]],
                         &yysemantic_stack_[0],
                         &yylocation_stack_[0]);
            yypop_ ();
          }
        throw;
      }
  }

  // Generate an error message.
  std::string
  Parser::yysyntax_error_ (int yystate, int yytoken)
  {
    std::string yyres;
    // Number of reported tokens (one for the "unexpected", one per
    // "expected").
    size_t yycount = 0;
    // Its maximum.
    enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
    // Arguments of yyformat.
    char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];

    /* There are many possibilities here to consider:
       - If this state is a consistent state with a default action, then
         the only way this function was invoked is if the default action
         is an error action.  In that case, don't check for expected
         tokens because there are none.
       - The only way there can be no lookahead present (in yytoken) is
         if this state is a consistent state with a default action.
         Thus, detecting the absence of a lookahead is sufficient to
         determine that there is no unexpected or expected token to
         report.  In that case, just report a simple "syntax error".
       - Don't assume there isn't a lookahead just because this state is
         a consistent state with a default action.  There might have
         been a previous inconsistent state, consistent state with a
         non-default action, or user semantic action that manipulated
         yychar.
       - Of course, the expected token list depends on states to have
         correct lookahead information, and it depends on the parser not
         to perform extra reductions after fetching a lookahead from the
         scanner and before detecting a syntax error.  Thus, state
         merging (from LALR or IELR) and default reductions corrupt the
         expected token list.  However, the list is correct for
         canonical LR with one exception: it will still contain any
         token that will not be accepted due to an error action in a
         later state.
    */
    if (yytoken != yyempty_)
      {
        yyarg[yycount++] = yytname_[yytoken];
        int yyn = yypact_[yystate];
        if (!yy_pact_value_is_default_ (yyn))
          {
            /* Start YYX at -YYN if negative to avoid negative indexes in
               YYCHECK.  In other words, skip the first -YYN actions for
               this state because they are default actions.  */
            int yyxbegin = yyn < 0 ? -yyn : 0;
            /* Stay within bounds of both yycheck and yytname.  */
            int yychecklim = yylast_ - yyn + 1;
            int yyxend = yychecklim < yyntokens_ ? yychecklim : yyntokens_;
            for (int yyx = yyxbegin; yyx < yyxend; ++yyx)
              if (yycheck_[yyx + yyn] == yyx && yyx != yyterror_
                  && !yy_table_value_is_error_ (yytable_[yyx + yyn]))
                {
                  if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                    {
                      yycount = 1;
                      break;
                    }
                  else
                    yyarg[yycount++] = yytname_[yyx];
                }
          }
      }

    char const* yyformat = YY_NULL;
    switch (yycount)
      {
#define YYCASE_(N, S)                         \
        case N:                               \
          yyformat = S;                       \
        break
        YYCASE_(0, YY_("syntax error"));
        YYCASE_(1, YY_("syntax error, unexpected %s"));
        YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
        YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
        YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
        YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
#undef YYCASE_
      }

    // Argument number.
    size_t yyi = 0;
    for (char const* yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount)
        {
          yyres += yytnamerr_ (yyarg[yyi++]);
          ++yyp;
        }
      else
        yyres += *yyp;
    return yyres;
  }


  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
  const signed char Parser::yypact_ninf_ = -36;
  const short int
  Parser::yypact_[] =
  {
       -36,     2,   -36,   -12,   245,   -36,   -36,   -36,   -36,   -36,
     -13,   278,    -5,   498,    19,    53,    27,    35,    44,   351,
     351,   -36,   371,   351,   351,    71,   162,    37,   -16,   862,
     930,   371,   351,   351,   351,   351,   351,   -36,   -36,   351,
     351,   351,   351,   351,   351,   -36,   -36,   351,   351,   351,
     351,   351,   351,   351,   -36,   -36,   351,   351,   371,   306,
     326,   371,   528,    45,   351,    42,   966,   731,   882,   -36,
      16,    16,    16,   -36,   -36,   -36,   -36,   -36,   -36,   -36,
     -36,   351,   351,   351,   351,   -36,   371,   371,   371,   351,
     -36,   351,   351,   351,   351,   351,   351,   351,   -36,   351,
     351,   351,   351,   351,   351,   351,   351,   351,   351,   351,
     371,   351,   351,   160,   966,   985,  1001,  1001,  1001,  1001,
    1001,   966,  1001,  1001,  1001,  1001,  1001,  1001,   966,  1001,
     966,  1001,  1001,  1001,  1001,  1001,  1001,   966,  1001,   646,
     160,   985,   -36,   205,   100,   497,   -36,   225,   148,   527,
     255,   123,   557,   351,    16,   -36,   -36,   351,   -36,   855,
    1001,   -26,    13,   946,   -35,   -35,   -36,  1017,  1017,    21,
      21,    21,    21,    21,    21,   -36,  1032,  1046,    86,    86,
      86,    86,    86,    86,    70,    70,    16,   -36,    16,    16,
      16,   351,   -36,   351,   -36,   351,   -36,   -36,   351,   -36,
     351,   -36,   -36,   -36,   351,   -36,  1001,    16,   351,   351,
     907,   415,   758,   467,   180,   387,   785,   966,  1001,    56,
     -36,   -36,   351,   -36,   351,   351,   351,   -36,   351,   587,
     675,   421,   174,  1001,   351,   -36,   351,   -36,   351,   -36,
     617,   812,   393,   351,   -36,   -36,   351,   704,   443,   351,
     -36,   839,   -36
  };

  /* YYDEFACT[S] -- default reduction number in state S.  Performed when
     YYTABLE doesn't specify something else to do.  Zero means the
     default is an error.  */
  const unsigned char
  Parser::yydefact_[] =
  {
         2,     0,     1,     0,     0,     4,     3,     9,    59,    42,
      86,    62,    43,    63,    44,    28,     0,     0,     0,     0,
       0,     8,     0,     0,     0,     0,     0,   118,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    90,    89,     0,
       0,     0,     0,     0,     0,    67,    66,     0,     0,     0,
       0,     0,     0,     0,    78,    77,     0,     0,     0,     0,
       0,     0,    86,     0,     0,   118,     0,     0,     0,    36,
     113,   114,    12,    61,    88,    65,    76,    60,    87,    64,
      75,     0,     0,     0,     0,     7,     0,     0,     0,     0,
       6,     0,     0,     0,     0,     0,     0,     0,     5,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    34,    45,    91,    93,    92,    95,    94,
      96,    47,    68,    71,    70,    73,    72,    74,    46,    69,
      49,    79,    82,    81,    84,    83,    85,    48,    80,     0,
      33,     0,    97,     0,     0,     0,    51,     0,     0,     0,
       0,     0,     0,     0,   113,    21,   116,     0,   117,     0,
       0,    19,    20,     0,    37,    35,    41,    38,    39,    27,
      26,    25,    24,    23,    22,    54,    17,    18,    16,    15,
      14,    13,    11,    10,   109,   108,   112,    40,   110,   111,
     115,     0,   100,     0,    99,     0,    98,    53,     0,    50,
       0,    52,    32,    29,     0,    31,    91,   110,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    58,   119,   120,
     102,   101,     0,   103,     0,     0,     0,    30,     0,     0,
       0,     0,     0,   121,     0,   104,     0,    57,     0,    56,
       0,     0,     0,     0,   106,   105,     0,     0,     0,     0,
      55,     0,   107
  };

  /* YYPGOTO[NTERM-NUM].  */
  const signed char
  Parser::yypgoto_[] =
  {
       -36,   -36,   -36,   -18,    59,    82,    -4
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  Parser::yydefgoto_[] =
  {
        -1,     1,     6,    27,    28,    66,   160
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If YYTABLE_NINF_, syntax error.  */
  const signed char Parser::yytable_ninf_ = -1;
  const unsigned char
  Parser::yytable_[] =
  {
        30,    65,     2,     3,    85,    83,    88,    89,     7,    31,
      32,    33,    34,    35,    36,    67,    68,    47,    70,    71,
      72,     4,    86,    87,    84,    88,    89,   115,   116,   117,
     118,   119,   120,    37,    38,   122,   123,   124,   125,   126,
     127,    56,    59,   129,   131,   132,   133,   134,   135,   136,
      60,     5,   138,   139,   141,   145,   149,   152,   155,    61,
     154,   112,    57,    84,   161,   162,    81,    82,    83,    97,
      57,    81,    82,    83,    73,    58,    74,    75,   228,    76,
     163,    69,   141,   141,   167,   168,    29,    84,     0,     0,
     113,     0,    84,     0,     0,   176,   177,   178,   179,   180,
     181,   182,   183,   184,   185,   186,   188,   189,   190,     0,
     109,   157,   111,   114,   193,   112,   194,   140,   143,   147,
     150,   121,     0,     0,   107,   108,   109,   157,   111,   128,
     130,   112,    91,    92,    93,    94,    95,    96,   137,   203,
       0,   144,   148,   151,     0,   164,   165,   166,    97,   206,
       0,     0,     0,   207,     0,    91,    92,    93,    94,    95,
      96,     0,   198,   159,   199,    77,     0,    78,    79,   187,
      80,    97,     0,   169,   170,   171,   172,   173,   174,   175,
      91,    92,    93,    94,    95,    96,     0,   210,   238,   212,
     239,   213,     0,     0,   225,     0,    97,     0,    86,    87,
     216,    88,    89,     0,     0,   218,    91,    92,    93,    94,
      95,    96,    91,    92,    93,    94,    95,    96,   229,     0,
     230,   192,    97,     0,   233,     0,     0,     0,    97,     0,
     240,     0,   241,     0,     0,   114,     0,     0,     0,   247,
       0,   197,     0,    86,    87,   251,    88,    89,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,     0,
      19,     0,    20,    86,    87,    21,    88,    89,     0,     0,
       0,   202,     0,     0,     0,   211,     0,     0,     0,     0,
     214,     0,   215,    22,    23,     0,     0,     0,    24,     0,
     217,    25,    26,    86,    87,     0,    88,    89,     0,     0,
      39,    40,    41,    42,    43,    44,     0,   231,   232,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
     242,    19,   142,    20,    45,    46,     0,     0,   248,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
       0,    19,   146,    20,    22,    23,     0,     0,     0,    24,
       0,     0,    25,    26,     8,     9,    62,    11,    12,    13,
      14,    63,    16,    17,    22,    23,    19,     0,    20,    24,
       0,     0,    25,    26,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,     0,    19,     0,    20,    64,
      23,     0,     0,     0,    24,     0,     0,    25,    26,     0,
       0,   226,     0,     0,     0,     0,     0,   246,     0,    22,
      23,     0,     0,     0,    24,     0,     0,    25,    26,    91,
      92,    93,    94,    95,    96,    91,    92,    93,    94,    95,
      96,   220,     0,     0,     0,    97,     0,   237,     0,     0,
       0,    97,     0,     0,     0,     0,     0,    91,    92,    93,
      94,    95,    96,    91,    92,    93,    94,    95,    96,   250,
       0,     0,     0,    97,     0,     0,     0,     0,     0,    97,
       0,     0,     0,     0,     0,    91,    92,    93,    94,    95,
      96,   222,     0,   223,     0,     0,     0,     0,   224,     0,
       0,    97,     0,     0,     0,     0,     0,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   109,   157,   111,
       0,   195,   112,   196,     0,     0,     0,     0,     0,     0,
      48,    49,    50,    51,    52,    53,     0,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   109,   110,   111,
       0,   200,   112,   201,    54,    55,     0,     0,     0,     0,
     153,    32,    33,    34,    35,    36,     0,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   109,   110,   111,
       0,   204,   112,   205,    37,    38,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   109,   110,   111,
       0,   234,   112,   235,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   109,   157,   111,
       0,   243,   112,   244,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   109,   157,   111,
     191,     0,   112,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    99,   100,   101,   102,
     103,   104,   105,   106,   107,   108,   109,   157,   111,   236,
       0,   112,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   157,   111,   249,     0,
     112,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   157,   111,   156,     0,   112,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    99,   100,   101,   102,   103,   104,   105,   106,   107,
     108,   109,   157,   111,   221,     0,   112,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    99,   100,
     101,   102,   103,   104,   105,   106,   107,   108,   109,   157,
     111,   227,     0,   112,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   157,   111,   245,     0,
     112,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    99,   100,   101,   102,   103,   104,   105,   106,
     107,   108,   109,   157,   111,   252,     0,   112,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,   109,
     157,   111,    90,   208,   112,     0,     0,    91,    92,    93,
      94,    95,    96,     0,    91,    92,    93,    94,    95,    96,
     158,     0,     0,    97,     0,     0,     0,     0,     0,     0,
      97,     0,    99,   100,   101,   102,   103,   104,   105,   106,
     107,   108,   109,   157,   111,   219,     0,   112,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   109,   157,   111,
      98,     0,   112,     0,     0,     0,     0,     0,     0,     0,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     109,   110,   111,     0,     0,   112,    99,   100,   101,   102,
     103,   104,   105,   106,   107,   108,   109,   157,   111,     0,
       0,   112,     0,     0,     0,     0,     0,   209,    91,    92,
      93,    94,    95,    96,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    97,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,   111,     0,     0,
     112,    99,   100,   101,   102,   103,   104,   105,   106,   107,
     108,   109,   157,   111,     0,     0,   112,    99,   100,   101,
     102,   103,   104,   105,   106,     0,     0,   109,     0,     0,
       0,     0,   112,   100,   101,   102,   103,   104,   105,   106,
     107,   108,   109,   157,   111,     0,     0,   112,   101,   102,
     103,   104,   105,   106,   107,   108,   109,   157,   111,     0,
       0,   112
  };

  /* YYCHECK.  */
  const short int
  Parser::yycheck_[] =
  {
         4,    19,     0,     1,    20,    31,    41,    42,    20,    22,
      23,    24,    25,    26,    27,    19,    20,    22,    22,    23,
      24,    19,    38,    39,    50,    41,    42,    31,    32,    33,
      34,    35,    36,    46,    47,    39,    40,    41,    42,    43,
      44,    22,    15,    47,    48,    49,    50,    51,    52,    53,
      15,    49,    56,    57,    58,    59,    60,    61,    16,    15,
      64,    45,    17,    50,    82,    83,    29,    30,    31,    48,
      17,    29,    30,    31,     3,    22,     5,     6,    22,     8,
      84,    22,    86,    87,    88,    89,     4,    50,    -1,    -1,
      31,    -1,    50,    -1,    -1,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,   111,   112,    -1,
      40,    41,    42,    31,    14,    45,    16,    58,    59,    60,
      61,    39,    -1,    -1,    38,    39,    40,    41,    42,    47,
      48,    45,    32,    33,    34,    35,    36,    37,    56,    16,
      -1,    59,    60,    61,    -1,    86,    87,    88,    48,   153,
      -1,    -1,    -1,   157,    -1,    32,    33,    34,    35,    36,
      37,    -1,    14,    81,    16,     3,    -1,     5,     6,   110,
       8,    48,    -1,    91,    92,    93,    94,    95,    96,    97,
      32,    33,    34,    35,    36,    37,    -1,   191,    14,   193,
      16,   195,    -1,    -1,    14,    -1,    48,    -1,    38,    39,
     204,    41,    42,    -1,    -1,   209,    32,    33,    34,    35,
      36,    37,    32,    33,    34,    35,    36,    37,   222,    -1,
     224,    16,    48,    -1,   228,    -1,    -1,    -1,    48,    -1,
     234,    -1,   236,    -1,    -1,   153,    -1,    -1,    -1,   243,
      -1,    16,    -1,    38,    39,   249,    41,    42,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    -1,
      15,    -1,    17,    38,    39,    20,    41,    42,    -1,    -1,
      -1,    16,    -1,    -1,    -1,   193,    -1,    -1,    -1,    -1,
     198,    -1,   200,    38,    39,    -1,    -1,    -1,    43,    -1,
     208,    46,    47,    38,    39,    -1,    41,    42,    -1,    -1,
      22,    23,    24,    25,    26,    27,    -1,   225,   226,     3,
       4,     5,     6,     7,     8,     9,    10,    11,    12,    13,
     238,    15,    16,    17,    46,    47,    -1,    -1,   246,     3,
       4,     5,     6,     7,     8,     9,    10,    11,    12,    13,
      -1,    15,    16,    17,    38,    39,    -1,    -1,    -1,    43,
      -1,    -1,    46,    47,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    38,    39,    15,    -1,    17,    43,
      -1,    -1,    46,    47,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    -1,    15,    -1,    17,    38,
      39,    -1,    -1,    -1,    43,    -1,    -1,    46,    47,    -1,
      -1,    14,    -1,    -1,    -1,    -1,    -1,    14,    -1,    38,
      39,    -1,    -1,    -1,    43,    -1,    -1,    46,    47,    32,
      33,    34,    35,    36,    37,    32,    33,    34,    35,    36,
      37,    16,    -1,    -1,    -1,    48,    -1,    16,    -1,    -1,
      -1,    48,    -1,    -1,    -1,    -1,    -1,    32,    33,    34,
      35,    36,    37,    32,    33,    34,    35,    36,    37,    16,
      -1,    -1,    -1,    48,    -1,    -1,    -1,    -1,    -1,    48,
      -1,    -1,    -1,    -1,    -1,    32,    33,    34,    35,    36,
      37,    14,    -1,    16,    -1,    -1,    -1,    -1,    21,    -1,
      -1,    48,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      -1,    14,    45,    16,    -1,    -1,    -1,    -1,    -1,    -1,
      22,    23,    24,    25,    26,    27,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      -1,    14,    45,    16,    46,    47,    -1,    -1,    -1,    -1,
      22,    23,    24,    25,    26,    27,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      -1,    14,    45,    16,    46,    47,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      -1,    14,    45,    16,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      -1,    14,    45,    16,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      14,    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    14,
      -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    14,    -1,
      45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    16,    -1,    45,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    16,    -1,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    16,    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    16,    -1,
      45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    16,    -1,    45,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    20,    28,    45,    -1,    -1,    32,    33,    34,
      35,    36,    37,    -1,    32,    33,    34,    35,    36,    37,
      18,    -1,    -1,    48,    -1,    -1,    -1,    -1,    -1,    -1,
      48,    -1,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    18,    -1,    45,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      20,    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    -1,    -1,    45,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    -1,
      -1,    45,    -1,    -1,    -1,    -1,    -1,    51,    32,    33,
      34,    35,    36,    37,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    48,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    -1,    -1,
      45,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    -1,    -1,    45,    30,    31,    32,
      33,    34,    35,    36,    37,    -1,    -1,    40,    -1,    -1,
      -1,    -1,    45,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    -1,    -1,    45,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    -1,
      -1,    45
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned char
  Parser::yystos_[] =
  {
         0,    53,     0,     1,    19,    49,    54,    20,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    15,
      17,    20,    38,    39,    43,    46,    47,    55,    56,    57,
      58,    22,    23,    24,    25,    26,    27,    46,    47,    22,
      23,    24,    25,    26,    27,    46,    47,    22,    22,    23,
      24,    25,    26,    27,    46,    47,    22,    17,    22,    15,
      15,    15,     5,    10,    38,    55,    57,    58,    58,    56,
      58,    58,    58,     3,     5,     6,     8,     3,     5,     6,
       8,    29,    30,    31,    50,    20,    38,    39,    41,    42,
      20,    32,    33,    34,    35,    36,    37,    48,    20,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    45,    56,    57,    58,    58,    58,    58,    58,
      58,    57,    58,    58,    58,    58,    58,    58,    57,    58,
      57,    58,    58,    58,    58,    58,    58,    57,    58,    58,
      56,    58,    16,    56,    57,    58,    16,    56,    57,    58,
      56,    57,    58,    22,    58,    16,    16,    41,    18,    57,
      58,    55,    55,    58,    56,    56,    56,    58,    58,    57,
      57,    57,    57,    57,    57,    57,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    56,    58,    58,
      58,    14,    16,    14,    16,    14,    16,    16,    14,    16,
      14,    16,    16,    16,    14,    16,    58,    58,    28,    51,
      58,    57,    58,    58,    57,    57,    58,    57,    58,    18,
      16,    16,    14,    16,    21,    14,    14,    16,    22,    58,
      58,    57,    57,    58,    14,    16,    14,    16,    14,    16,
      58,    58,    57,    14,    16,    16,    14,    58,    57,    14,
      16,    58,    16
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
     295,   296,   297,   298,   299,   300,   301,   302,   303,    10,
      63,    58
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned char
  Parser::yyr1_[] =
  {
         0,    52,    53,    53,    54,    54,    54,    54,    54,    54,
      55,    55,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    55,    55,    55,    55,    55,    55,    55,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    57,    57,    57,    57,    57,    57,    57,    57,
      57,    57,    57,    57,    57,    57,    57,    57,    57,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  Parser::yyr2_[] =
  {
         0,     2,     0,     2,     1,     3,     3,     3,     2,     2,
       3,     3,     2,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     1,     4,
       6,     4,     4,     3,     3,     3,     2,     3,     3,     3,
       3,     3,     1,     1,     1,     3,     3,     3,     3,     3,
       4,     3,     4,     4,     3,    12,     8,     8,     5,     1,
       2,     2,     1,     1,     2,     2,     2,     2,     3,     3,
       3,     3,     3,     3,     3,     2,     2,     2,     2,     3,
       3,     3,     3,     3,     3,     3,     1,     2,     2,     2,
       2,     3,     3,     3,     3,     3,     3,     3,     4,     4,
       4,     6,     6,     6,     8,    10,    10,    14,     3,     3,
       3,     3,     3,     2,     2,     3,     3,     3,     1,     5,
       6,     8
  };


  /* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
     First, the terminals, then, starting at \a yyntokens_, nonterminals.  */
  const char*
  const Parser::yytname_[] =
  {
    "\"end of file\"", "error", "$undefined", "NUM", "QSTRING", "UNDVAR",
  "VAR", "SVAR", "IMMVAR", "IMMSVAR", "AVAR", "FNCT", "SFNCT", "AFNCT",
  "COMMA", "LPAR", "RPAR", "LBRACK", "RBRACK", "LBRACE", "RBRACE", "SEMI",
  "EQUAL", "EQ_MINUS", "EQ_PLUS", "EQ_DIV", "EQ_TIME", "EQ_POW", "COLON",
  "QUEST", "LOR", "LAND", "NE", "EQ", "GE", "LE", "GT", "LT", "SUB", "PLU",
  "MOD", "TIM", "DIV", "NOT", "UNARY", "POW", "DEC", "INC", "CONCAT",
  "'\\n'", "'?'", "':'", "$accept", "input", "line", "bool", "aexp",
  "sexp", "exp", YY_NULL
  };

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const Parser::rhs_number_type
  Parser::yyrhs_[] =
  {
        53,     0,    -1,    -1,    53,    54,    -1,    49,    -1,    19,
      58,    20,    -1,    19,    57,    20,    -1,    19,    56,    20,
      -1,    19,    20,    -1,     1,    20,    -1,    58,    37,    58,
      -1,    58,    36,    58,    -1,    43,    58,    -1,    58,    35,
      58,    -1,    58,    34,    58,    -1,    58,    33,    58,    -1,
      58,    32,    58,    -1,    58,    30,    58,    -1,    58,    31,
      58,    -1,    55,    30,    55,    -1,    55,    31,    55,    -1,
      15,    55,    16,    -1,    57,    37,    57,    -1,    57,    36,
      57,    -1,    57,    35,    57,    -1,    57,    34,    57,    -1,
      57,    33,    57,    -1,    57,    32,    57,    -1,    10,    -1,
      13,    15,    57,    16,    -1,    13,    15,    58,    14,    58,
      16,    -1,    13,    15,    58,    16,    -1,    13,    15,    56,
      16,    -1,    10,    22,    56,    -1,     5,    22,    56,    -1,
      56,    39,    56,    -1,    38,    56,    -1,    56,    38,    56,
      -1,    56,    41,    58,    -1,    56,    42,    58,    -1,    58,
      41,    56,    -1,    56,    41,    56,    -1,     4,    -1,     7,
      -1,     9,    -1,     5,    22,    57,    -1,     7,    22,    57,
      -1,     6,    22,    57,    -1,     9,    22,    57,    -1,     8,
      22,    57,    -1,    12,    15,    57,    16,    -1,    12,    15,
      16,    -1,    12,    15,    58,    16,    -1,    12,    15,    56,
      16,    -1,    57,    48,    57,    -1,    12,    15,    58,    14,
      57,    14,    57,    14,    57,    14,    57,    16,    -1,    12,
      15,    58,    14,    57,    14,    57,    16,    -1,    12,    15,
      57,    14,    57,    14,    57,    16,    -1,    55,    29,    57,
      28,    57,    -1,     3,    -1,    47,     3,    -1,    46,     3,
      -1,     6,    -1,     8,    -1,    47,     6,    -1,    46,     6,
      -1,     6,    47,    -1,     6,    46,    -1,     6,    22,    58,
      -1,     7,    22,    58,    -1,     6,    24,    58,    -1,     6,
      23,    58,    -1,     6,    26,    58,    -1,     6,    25,    58,
      -1,     6,    27,    58,    -1,    47,     8,    -1,    46,     8,
      -1,     8,    47,    -1,     8,    46,    -1,     8,    22,    58,
      -1,     9,    22,    58,    -1,     8,    24,    58,    -1,     8,
      23,    58,    -1,     8,    26,    58,    -1,     8,    25,    58,
      -1,     8,    27,    58,    -1,     5,    -1,    47,     5,    -1,
      46,     5,    -1,     5,    47,    -1,     5,    46,    -1,     5,
      22,    58,    -1,     5,    24,    58,    -1,     5,    23,    58,
      -1,     5,    26,    58,    -1,     5,    25,    58,    -1,     5,
      27,    58,    -1,    11,    15,    16,    -1,    11,    15,    58,
      16,    -1,    11,    15,    57,    16,    -1,    11,    15,    56,
      16,    -1,    11,    15,    57,    14,    58,    16,    -1,    11,
      15,    57,    14,    57,    16,    -1,    11,    15,    58,    14,
      58,    16,    -1,    11,    15,    58,    14,    58,    14,    58,
      16,    -1,    11,    15,    58,    14,    58,    21,    58,    14,
      58,    16,    -1,    11,    15,    58,    14,    58,    14,    58,
      14,    58,    16,    -1,    11,    15,    58,    14,    58,    14,
      58,    14,    58,    14,    58,    14,    58,    16,    -1,    58,
      39,    58,    -1,    58,    38,    58,    -1,    58,    41,    58,
      -1,    58,    42,    58,    -1,    58,    40,    58,    -1,    38,
      58,    -1,    39,    58,    -1,    58,    45,    58,    -1,    15,
      58,    16,    -1,    17,    58,    18,    -1,    55,    -1,    55,
      50,    58,    51,    58,    -1,    10,    17,    58,    14,    58,
      18,    -1,    10,    17,    58,    14,    58,    18,    22,    58,
      -1
  };

  /* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
     YYRHS.  */
  const unsigned short int
  Parser::yyprhs_[] =
  {
         0,     0,     3,     4,     7,     9,    13,    17,    21,    24,
      27,    31,    35,    38,    42,    46,    50,    54,    58,    62,
      66,    70,    74,    78,    82,    86,    90,    94,    98,   100,
     105,   112,   117,   122,   126,   130,   134,   137,   141,   145,
     149,   153,   157,   159,   161,   163,   167,   171,   175,   179,
     183,   188,   192,   197,   202,   206,   219,   228,   237,   243,
     245,   248,   251,   253,   255,   258,   261,   264,   267,   271,
     275,   279,   283,   287,   291,   295,   298,   301,   304,   307,
     311,   315,   319,   323,   327,   331,   335,   337,   340,   343,
     346,   349,   353,   357,   361,   365,   369,   373,   377,   382,
     387,   392,   399,   406,   413,   422,   433,   444,   459,   463,
     467,   471,   475,   479,   482,   485,   489,   493,   497,   499,
     505,   512
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  Parser::yyrline_[] =
  {
         0,    99,    99,   100,   103,   104,   111,   115,   116,   117,
     120,   121,   122,   123,   124,   125,   126,   127,   128,   129,
     130,   131,   134,   135,   136,   137,   138,   139,   141,   142,
     143,   145,   147,   149,   152,   154,   162,   164,   172,   173,
     174,   175,   184,   185,   186,   187,   189,   192,   196,   197,
     198,   199,   200,   201,   202,   203,   205,   207,   209,   211,
     212,   213,   214,   215,   216,   217,   218,   219,   220,   222,
     225,   226,   227,   228,   229,   234,   235,   236,   237,   238,
     239,   240,   241,   242,   243,   244,   246,   248,   251,   254,
     257,   260,   262,   265,   268,   271,   274,   280,   281,   282,
     283,   284,   286,   288,   290,   292,   294,   296,   298,   299,
     300,   301,   308,   315,   316,   317,   320,   321,   324,   325,
     326,   338
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
      49,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    51,     2,
       2,     2,     2,    50,     2,     2,     2,     2,     2,     2,
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
      45,    46,    47,    48
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int Parser::yyeof_ = 0;
  const int Parser::yylast_ = 1091;
  const int Parser::yynnts_ = 7;
  const int Parser::yyempty_ = -2;
  const int Parser::yyfinal_ = 2;
  const int Parser::yyterror_ = 1;
  const int Parser::yyerrcode_ = 256;
  const int Parser::yyntokens_ = 52;

  const unsigned int Parser::yyuser_token_number_max_ = 303;
  const Parser::token_number_type Parser::yyundef_token_ = 2;


} // SEAMS
/* Line 1135 of lalr1.cc  */
#line 2196 "apr_parser.cc"
/* Line 1136 of lalr1.cc  */
#line 355 "aprepro.yy"


void SEAMS::Parser::error(const Parser::location_type&, const std::string& m)
{
    aprepro.error(m);
}

