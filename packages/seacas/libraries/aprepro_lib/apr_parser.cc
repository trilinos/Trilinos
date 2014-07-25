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

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <cstdio>

namespace SEAMS {
   extern int echo;
 }

 
/* Line 283 of lalr1.cc  */
#line 58 "apr_parser.cc"


#include "aprepro_parser.h"

/* User implementation prologue.  */
/* Line 289 of lalr1.cc  */
#line 86 "aprepro.yy"


#include "aprepro.h"
#include "apr_scanner.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex


/* Line 289 of lalr1.cc  */
#line 79 "apr_parser.cc"


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
#line 174 "apr_parser.cc"

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
#line 105 "aprepro.yy"
    { if (echo) aprepro.lexer->LexerOutput("\n", 1); }
    break;

  case 5:
/* Line 664 of lalr1.cc  */
#line 106 "aprepro.yy"
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
#line 113 "aprepro.yy"
    { if (echo && (yysemantic_stack_[(3) - (2)].string) != NULL) {
				    aprepro.lexer->LexerOutput((yysemantic_stack_[(3) - (2)].string), strlen((yysemantic_stack_[(3) - (2)].string)));
                                  }
                                }
    break;

  case 7:
/* Line 664 of lalr1.cc  */
#line 117 "aprepro.yy"
    {                                       }
    break;

  case 8:
/* Line 664 of lalr1.cc  */
#line 118 "aprepro.yy"
    {                                       }
    break;

  case 9:
/* Line 664 of lalr1.cc  */
#line 119 "aprepro.yy"
    { yyerrok;				}
    break;

  case 10:
/* Line 664 of lalr1.cc  */
#line 122 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) < (yysemantic_stack_[(3) - (3)].val);                         }
    break;

  case 11:
/* Line 664 of lalr1.cc  */
#line 123 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) > (yysemantic_stack_[(3) - (3)].val);                         }
    break;

  case 12:
/* Line 664 of lalr1.cc  */
#line 124 "aprepro.yy"
    { (yyval.val) = !((yysemantic_stack_[(2) - (2)].val));                           }
    break;

  case 13:
/* Line 664 of lalr1.cc  */
#line 125 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) <= (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 14:
/* Line 664 of lalr1.cc  */
#line 126 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) >= (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 15:
/* Line 664 of lalr1.cc  */
#line 127 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) == (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 16:
/* Line 664 of lalr1.cc  */
#line 128 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) != (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 17:
/* Line 664 of lalr1.cc  */
#line 129 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) || (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 18:
/* Line 664 of lalr1.cc  */
#line 130 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) && (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 19:
/* Line 664 of lalr1.cc  */
#line 131 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) || (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 20:
/* Line 664 of lalr1.cc  */
#line 132 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) && (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 21:
/* Line 664 of lalr1.cc  */
#line 133 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (2)].val);                              }
    break;

  case 22:
/* Line 664 of lalr1.cc  */
#line 136 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) <  0 ? 1 : 0);	}
    break;

  case 23:
/* Line 664 of lalr1.cc  */
#line 137 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) >  0 ? 1 : 0);	}
    break;

  case 24:
/* Line 664 of lalr1.cc  */
#line 138 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) <= 0 ? 1 : 0);	}
    break;

  case 25:
/* Line 664 of lalr1.cc  */
#line 139 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) >= 0 ? 1 : 0);	}
    break;

  case 26:
/* Line 664 of lalr1.cc  */
#line 140 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) == 0 ? 1 : 0);	}
    break;

  case 27:
/* Line 664 of lalr1.cc  */
#line 141 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) != 0 ? 1 : 0);	}
    break;

  case 28:
/* Line 664 of lalr1.cc  */
#line 143 "aprepro.yy"
    { (yyval.arrval) = (yysemantic_stack_[(1) - (1)].tptr)->value.avar;}
    break;

  case 29:
/* Line 664 of lalr1.cc  */
#line 144 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(4) - (1)].tptr), (yysemantic_stack_[(4) - (1)].tptr)->value.arrfnct_c == NULL))
	    (yyval.arrval) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.arrfnct_c))((yysemantic_stack_[(4) - (3)].string));
	  else
	    yyerrok;
	}
    break;

  case 30:
/* Line 664 of lalr1.cc  */
#line 150 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(6) - (1)].tptr), (yysemantic_stack_[(6) - (1)].tptr)->value.arrfnct_cd == NULL))
	    (yyval.arrval) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.arrfnct_cd))((yysemantic_stack_[(6) - (3)].string),(yysemantic_stack_[(6) - (5)].val));
	  else
	    yyerrok;
	}
    break;

  case 31:
/* Line 664 of lalr1.cc  */
#line 156 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(6) - (1)].tptr), (yysemantic_stack_[(6) - (1)].tptr)->value.arrfnct_cc == NULL))
	    (yyval.arrval) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.arrfnct_cc))((yysemantic_stack_[(6) - (3)].string),(yysemantic_stack_[(6) - (5)].string));
	  else
	    yyerrok;
	}
    break;

  case 32:
/* Line 664 of lalr1.cc  */
#line 162 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(6) - (1)].tptr), (yysemantic_stack_[(6) - (1)].tptr)->value.arrfnct_dd == NULL))
	    (yyval.arrval) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.arrfnct_dd))((yysemantic_stack_[(6) - (3)].val),(yysemantic_stack_[(6) - (5)].val));
	  else
	    yyerrok;
	}
    break;

  case 33:
/* Line 664 of lalr1.cc  */
#line 168 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(4) - (1)].tptr), (yysemantic_stack_[(4) - (1)].tptr)->value.arrfnct_d == NULL))
	    (yyval.arrval) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.arrfnct_d))((yysemantic_stack_[(4) - (3)].val));
	  else
	    yyerrok;
	}
    break;

  case 34:
/* Line 664 of lalr1.cc  */
#line 174 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(4) - (1)].tptr), (yysemantic_stack_[(4) - (1)].tptr)->value.arrfnct_a == NULL))
	    (yyval.arrval) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.arrfnct_a))((yysemantic_stack_[(4) - (3)].arrval));
	  else
	    yyerrok;
	}
    break;

  case 35:
/* Line 664 of lalr1.cc  */
#line 180 "aprepro.yy"
    { (yyval.arrval) = (yysemantic_stack_[(3) - (3)].arrval); delete (yysemantic_stack_[(3) - (1)].tptr)->value.avar; (yysemantic_stack_[(3) - (1)].tptr)->value.avar = (yysemantic_stack_[(3) - (3)].arrval); 
                                  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr));
                                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::AVAR); }
    break;

  case 36:
/* Line 664 of lalr1.cc  */
#line 183 "aprepro.yy"
    { (yyval.arrval) = (yysemantic_stack_[(3) - (3)].arrval); (yysemantic_stack_[(3) - (1)].tptr)->value.avar = (yysemantic_stack_[(3) - (3)].arrval); 
                                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::AVAR); }
    break;

  case 37:
/* Line 664 of lalr1.cc  */
#line 185 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (1)].arrval)->cols == (yysemantic_stack_[(3) - (3)].arrval)->cols && (yysemantic_stack_[(3) - (1)].arrval)->rows == (yysemantic_stack_[(3) - (3)].arrval)->rows ) {
                                     (yyval.arrval) = array_add((yysemantic_stack_[(3) - (1)].arrval), (yysemantic_stack_[(3) - (3)].arrval)); 
                                  }
                                  else {
                                    yyerror(aprepro, "Arrays do not have same row and column count"); 
                                    yyerrok;
                                  }
                                }
    break;

  case 38:
/* Line 664 of lalr1.cc  */
#line 193 "aprepro.yy"
    { (yyval.arrval) = array_scale((yysemantic_stack_[(2) - (2)].arrval), -1.0);           }
    break;

  case 39:
/* Line 664 of lalr1.cc  */
#line 195 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (1)].arrval)->cols == (yysemantic_stack_[(3) - (3)].arrval)->cols && (yysemantic_stack_[(3) - (1)].arrval)->rows == (yysemantic_stack_[(3) - (3)].arrval)->rows ) {
                                     (yyval.arrval) = array_sub((yysemantic_stack_[(3) - (1)].arrval), (yysemantic_stack_[(3) - (3)].arrval)); 
                                  }
                                  else {
                                    yyerror(aprepro, "Arrays do not have same row and column count"); 
                                    yyerrok;
                                  }
                                }
    break;

  case 40:
/* Line 664 of lalr1.cc  */
#line 203 "aprepro.yy"
    { (yyval.arrval) = array_scale((yysemantic_stack_[(3) - (1)].arrval), (yysemantic_stack_[(3) - (3)].val));             }
    break;

  case 41:
/* Line 664 of lalr1.cc  */
#line 204 "aprepro.yy"
    { (yyval.arrval) = array_scale((yysemantic_stack_[(3) - (1)].arrval), 1.0/(yysemantic_stack_[(3) - (3)].val));         }
    break;

  case 42:
/* Line 664 of lalr1.cc  */
#line 205 "aprepro.yy"
    { (yyval.arrval) = array_scale((yysemantic_stack_[(3) - (3)].arrval), (yysemantic_stack_[(3) - (1)].val));             }
    break;

  case 43:
/* Line 664 of lalr1.cc  */
#line 206 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (1)].arrval)->cols == (yysemantic_stack_[(3) - (3)].arrval)->rows) {
                                    (yyval.arrval) = array_mult((yysemantic_stack_[(3) - (1)].arrval), (yysemantic_stack_[(3) - (3)].arrval));
                                  }
                                  else {
                                    yyerror(aprepro, "Column count of first array does not match row count of second array"); 
                                    yyerrok;
                                  }
				}
    break;

  case 44:
/* Line 664 of lalr1.cc  */
#line 215 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(1) - (1)].string);				}
    break;

  case 45:
/* Line 664 of lalr1.cc  */
#line 216 "aprepro.yy"
    { (yyval.string) = (char*)(yysemantic_stack_[(1) - (1)].tptr)->value.svar;			}
    break;

  case 46:
/* Line 664 of lalr1.cc  */
#line 217 "aprepro.yy"
    { (yyval.string) = (char*)(yysemantic_stack_[(1) - (1)].tptr)->value.svar;			}
    break;

  case 47:
/* Line 664 of lalr1.cc  */
#line 218 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); (yysemantic_stack_[(3) - (1)].tptr)->value.svar = (yysemantic_stack_[(3) - (3)].string);
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), Parser::token::SVAR);	}
    break;

  case 48:
/* Line 664 of lalr1.cc  */
#line 220 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); 
				  (yysemantic_stack_[(3) - (1)].tptr)->value.svar = (yysemantic_stack_[(3) - (3)].string);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr));          }
    break;

  case 49:
/* Line 664 of lalr1.cc  */
#line 223 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); 
				  (yysemantic_stack_[(3) - (1)].tptr)->value.svar= (yysemantic_stack_[(3) - (3)].string);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr));          
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::SVAR);		}
    break;

  case 50:
/* Line 664 of lalr1.cc  */
#line 227 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 51:
/* Line 664 of lalr1.cc  */
#line 228 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 52:
/* Line 664 of lalr1.cc  */
#line 229 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(4) - (1)].tptr), (yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_c == NULL))
	    (yyval.string) = (char*)(*((yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_c))((yysemantic_stack_[(4) - (3)].string));
	  else
	    (yyval.string) = (char*)"";
	}
    break;

  case 53:
/* Line 664 of lalr1.cc  */
#line 235 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(3) - (1)].tptr), (yysemantic_stack_[(3) - (1)].tptr)->value.strfnct == NULL))
	    (yyval.string) = (char*)(*((yysemantic_stack_[(3) - (1)].tptr)->value.strfnct))();
	  else
	    (yyval.string) = (char*)"";
	}
    break;

  case 54:
/* Line 664 of lalr1.cc  */
#line 241 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(4) - (1)].tptr), (yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_d == NULL))
	    (yyval.string) = (char*)(*((yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_d))((yysemantic_stack_[(4) - (3)].val));
	  else
	    (yyval.string) = (char*)"";
	}
    break;

  case 55:
/* Line 664 of lalr1.cc  */
#line 247 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(4) - (1)].tptr), (yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_a == NULL))
	    (yyval.string) = (char*)(*((yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_a))((yysemantic_stack_[(4) - (3)].arrval));
	  else
	    (yyval.string) = (char*)"";
	}
    break;

  case 56:
/* Line 664 of lalr1.cc  */
#line 253 "aprepro.yy"
    { concat_string((yysemantic_stack_[(3) - (1)].string), (yysemantic_stack_[(3) - (3)].string), &(yyval.string)); }
    break;

  case 57:
/* Line 664 of lalr1.cc  */
#line 254 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(6) - (1)].tptr), (yysemantic_stack_[(6) - (1)].tptr)->value.strfnct_dd == NULL))
	    (yyval.string) = (char*)(*((yysemantic_stack_[(6) - (1)].tptr)->value.strfnct_dd))((yysemantic_stack_[(6) - (3)].val), (yysemantic_stack_[(6) - (5)].val));
	  else
	    (yyval.string) = (char*)"";
	}
    break;

  case 58:
/* Line 664 of lalr1.cc  */
#line 260 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(12) - (1)].tptr), (yysemantic_stack_[(12) - (1)].tptr)->value.strfnct_dcccc == NULL))
	    (yyval.string) = (char*)(*((yysemantic_stack_[(12) - (1)].tptr)->value.strfnct_dcccc))((yysemantic_stack_[(12) - (3)].val), (yysemantic_stack_[(12) - (5)].string), (yysemantic_stack_[(12) - (7)].string), (yysemantic_stack_[(12) - (9)].string), (yysemantic_stack_[(12) - (11)].string));
	  else
	    (yyval.string) = (char*)"";
	}
    break;

  case 59:
/* Line 664 of lalr1.cc  */
#line 266 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(8) - (1)].tptr), (yysemantic_stack_[(8) - (1)].tptr)->value.strfnct_dcc == NULL))
	    (yyval.string) = (char*)(*((yysemantic_stack_[(8) - (1)].tptr)->value.strfnct_dcc))((yysemantic_stack_[(8) - (3)].val), (yysemantic_stack_[(8) - (5)].string), (yysemantic_stack_[(8) - (7)].string));
	  else
	    (yyval.string) = (char*)"";
	}
    break;

  case 60:
/* Line 664 of lalr1.cc  */
#line 272 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(8) - (1)].tptr), (yysemantic_stack_[(8) - (1)].tptr)->value.strfnct_ccc == NULL))
	    (yyval.string) = (char*)(*((yysemantic_stack_[(8) - (1)].tptr)->value.strfnct_ccc))((yysemantic_stack_[(8) - (3)].string), (yysemantic_stack_[(8) - (5)].string), (yysemantic_stack_[(8) - (7)].string));
	  else
	    (yyval.string) = (char*)"";
	}
    break;

  case 61:
/* Line 664 of lalr1.cc  */
#line 278 "aprepro.yy"
    { (yyval.string) = ((yysemantic_stack_[(5) - (1)].val)) ? ((yysemantic_stack_[(5) - (3)].string)) : ((yysemantic_stack_[(5) - (5)].string));              }
    break;

  case 62:
/* Line 664 of lalr1.cc  */
#line 280 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].val); 				}
    break;

  case 63:
/* Line 664 of lalr1.cc  */
#line 281 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(2) - (2)].val) + 1;				}
    break;

  case 64:
/* Line 664 of lalr1.cc  */
#line 282 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(2) - (2)].val) - 1;				}
    break;

  case 65:
/* Line 664 of lalr1.cc  */
#line 283 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;			}
    break;

  case 66:
/* Line 664 of lalr1.cc  */
#line 284 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;			}
    break;

  case 67:
/* Line 664 of lalr1.cc  */
#line 285 "aprepro.yy"
    { (yyval.val) = ++((yysemantic_stack_[(2) - (2)].tptr)->value.var);		}
    break;

  case 68:
/* Line 664 of lalr1.cc  */
#line 286 "aprepro.yy"
    { (yyval.val) = --((yysemantic_stack_[(2) - (2)].tptr)->value.var);		}
    break;

  case 69:
/* Line 664 of lalr1.cc  */
#line 287 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)++;		}
    break;

  case 70:
/* Line 664 of lalr1.cc  */
#line 288 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)--;		}
    break;

  case 71:
/* Line 664 of lalr1.cc  */
#line 289 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr));          }
    break;

  case 72:
/* Line 664 of lalr1.cc  */
#line 291 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr));          
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);			}
    break;

  case 73:
/* Line 664 of lalr1.cc  */
#line 294 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var += (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 74:
/* Line 664 of lalr1.cc  */
#line 295 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var -= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 75:
/* Line 664 of lalr1.cc  */
#line 296 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var *= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 76:
/* Line 664 of lalr1.cc  */
#line 297 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var /= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 77:
/* Line 664 of lalr1.cc  */
#line 298 "aprepro.yy"
    { errno = 0;
				  (yysemantic_stack_[(3) - (1)].tptr)->value.var = std::pow((yysemantic_stack_[(3) - (1)].tptr)->value.var,(yysemantic_stack_[(3) - (3)].val)); 
				  (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
				  SEAMS::math_error(aprepro, "Power");
				}
    break;

  case 78:
/* Line 664 of lalr1.cc  */
#line 303 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (2)].tptr)); YYERROR; }
    break;

  case 79:
/* Line 664 of lalr1.cc  */
#line 304 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (2)].tptr)); YYERROR; }
    break;

  case 80:
/* Line 664 of lalr1.cc  */
#line 305 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (1)].tptr)); YYERROR; }
    break;

  case 81:
/* Line 664 of lalr1.cc  */
#line 306 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(2) - (1)].tptr)); YYERROR; }
    break;

  case 82:
/* Line 664 of lalr1.cc  */
#line 307 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 83:
/* Line 664 of lalr1.cc  */
#line 308 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 84:
/* Line 664 of lalr1.cc  */
#line 309 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 85:
/* Line 664 of lalr1.cc  */
#line 310 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 86:
/* Line 664 of lalr1.cc  */
#line 311 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 87:
/* Line 664 of lalr1.cc  */
#line 312 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 88:
/* Line 664 of lalr1.cc  */
#line 313 "aprepro.yy"
    { immutable_modify(aprepro, (yysemantic_stack_[(3) - (1)].tptr)); YYERROR; }
    break;

  case 89:
/* Line 664 of lalr1.cc  */
#line 315 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;
				  undefined_warning(aprepro, (yysemantic_stack_[(1) - (1)].tptr)->name);          }
    break;

  case 90:
/* Line 664 of lalr1.cc  */
#line 317 "aprepro.yy"
    { (yyval.val) = ++((yysemantic_stack_[(2) - (2)].tptr)->value.var);		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (2)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (2)].tptr)->name);          }
    break;

  case 91:
/* Line 664 of lalr1.cc  */
#line 320 "aprepro.yy"
    { (yyval.val) = --((yysemantic_stack_[(2) - (2)].tptr)->value.var);		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (2)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (2)].tptr)->name);          }
    break;

  case 92:
/* Line 664 of lalr1.cc  */
#line 323 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)++;		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (1)].tptr)->name);          }
    break;

  case 93:
/* Line 664 of lalr1.cc  */
#line 326 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)--;		
		                  set_type(aprepro, (yysemantic_stack_[(2) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (1)].tptr)->name);          }
    break;

  case 94:
/* Line 664 of lalr1.cc  */
#line 329 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);                      }
    break;

  case 95:
/* Line 664 of lalr1.cc  */
#line 331 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var += (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 96:
/* Line 664 of lalr1.cc  */
#line 334 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var -= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 97:
/* Line 664 of lalr1.cc  */
#line 337 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var *= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 98:
/* Line 664 of lalr1.cc  */
#line 340 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var /= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 99:
/* Line 664 of lalr1.cc  */
#line 343 "aprepro.yy"
    { errno = 0;
				  (yysemantic_stack_[(3) - (1)].tptr)->value.var = std::pow((yysemantic_stack_[(3) - (1)].tptr)->value.var,(yysemantic_stack_[(3) - (3)].val)); 
				  (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
		                  set_type(aprepro, (yysemantic_stack_[(3) - (1)].tptr), token::VAR);
				  SEAMS::math_error(aprepro, "Power");
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 100:
/* Line 664 of lalr1.cc  */
#line 350 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(3) - (1)].tptr), (yysemantic_stack_[(3) - (1)].tptr)->value.fnctptr == NULL))
	    (yyval.val) = (*((yysemantic_stack_[(3) - (1)].tptr)->value.fnctptr))();
	  else 
	    (yyval.val) = 0.0;
	  }
    break;

  case 101:
/* Line 664 of lalr1.cc  */
#line 357 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(4) - (1)].tptr), (yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_d == NULL))
	    (yyval.val) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_d))((yysemantic_stack_[(4) - (3)].val));
	  else
	    (yyval.val) = 0.0;
	  }
    break;

  case 102:
/* Line 664 of lalr1.cc  */
#line 364 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(4) - (1)].tptr), (yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_c == NULL))
	    (yyval.val) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_c))((yysemantic_stack_[(4) - (3)].string));
	  else
	    (yyval.val) = 0.0;
	  }
    break;

  case 103:
/* Line 664 of lalr1.cc  */
#line 371 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(4) - (1)].tptr), (yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_a == NULL))
	    (yyval.val) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_a))((yysemantic_stack_[(4) - (3)].arrval));
	  else
	    (yyval.val) = 0.0;
	  }
    break;

  case 104:
/* Line 664 of lalr1.cc  */
#line 378 "aprepro.yy"
    {
	    if (arg_check((yysemantic_stack_[(6) - (1)].tptr), (yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_cd == NULL))
	      (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_cd))((yysemantic_stack_[(6) - (3)].string), (yysemantic_stack_[(6) - (5)].val));
	    else
	      (yyval.val) = 0.0;
	  }
    break;

  case 105:
/* Line 664 of lalr1.cc  */
#line 385 "aprepro.yy"
    {
	  if (arg_check((yysemantic_stack_[(6) - (1)].tptr), (yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_dc == NULL))
	    (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_dc))((yysemantic_stack_[(6) - (3)].val), (yysemantic_stack_[(6) - (5)].string));
	  else
	    (yyval.val) = 0.0;
	  }
    break;

  case 106:
/* Line 664 of lalr1.cc  */
#line 392 "aprepro.yy"
    {
	    if (arg_check((yysemantic_stack_[(6) - (1)].tptr), (yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_cc == NULL))
	      (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_cc))((yysemantic_stack_[(6) - (3)].string), (yysemantic_stack_[(6) - (5)].string));
	    else
	      (yyval.val) = 0.0;
	  }
    break;

  case 107:
/* Line 664 of lalr1.cc  */
#line 399 "aprepro.yy"
    {
	    if (arg_check((yysemantic_stack_[(6) - (1)].tptr), (yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_dd == NULL))
	      (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_dd))((yysemantic_stack_[(6) - (3)].val), (yysemantic_stack_[(6) - (5)].val));
	    else
	      (yyval.val) = 0.0;
	  }
    break;

  case 108:
/* Line 664 of lalr1.cc  */
#line 405 "aprepro.yy"
    {
	    if (arg_check((yysemantic_stack_[(8) - (1)].tptr), (yysemantic_stack_[(8) - (1)].tptr)->value.fnctptr_ddd == NULL))
	      (yyval.val) = (*((yysemantic_stack_[(8) - (1)].tptr)->value.fnctptr_ddd))((yysemantic_stack_[(8) - (3)].val), (yysemantic_stack_[(8) - (5)].val), (yysemantic_stack_[(8) - (7)].val));
	    else
	      (yyval.val) = 0.0;
	  }
    break;

  case 109:
/* Line 664 of lalr1.cc  */
#line 411 "aprepro.yy"
    {
	    if (arg_check((yysemantic_stack_[(8) - (1)].tptr), (yysemantic_stack_[(8) - (1)].tptr)->value.fnctptr_ccd == NULL))
	      (yyval.val) = (*((yysemantic_stack_[(8) - (1)].tptr)->value.fnctptr_ccd))((yysemantic_stack_[(8) - (3)].string), (yysemantic_stack_[(8) - (5)].string), (yysemantic_stack_[(8) - (7)].val));
	    else
	      (yyval.val) = 0.0;
	  }
    break;

  case 110:
/* Line 664 of lalr1.cc  */
#line 417 "aprepro.yy"
    {
	    if (arg_check((yysemantic_stack_[(10) - (1)].tptr), (yysemantic_stack_[(10) - (1)].tptr)->value.fnctptr_dddd == NULL))
	      (yyval.val) = (*((yysemantic_stack_[(10) - (1)].tptr)->value.fnctptr_dddd))((yysemantic_stack_[(10) - (3)].val), (yysemantic_stack_[(10) - (5)].val), (yysemantic_stack_[(10) - (7)].val), (yysemantic_stack_[(10) - (9)].val));
	    else
	      (yyval.val) = 0.0;
	  }
    break;

  case 111:
/* Line 664 of lalr1.cc  */
#line 423 "aprepro.yy"
    {
	    if (arg_check((yysemantic_stack_[(10) - (1)].tptr), (yysemantic_stack_[(10) - (1)].tptr)->value.fnctptr_dddd == NULL))
	      (yyval.val) = (*((yysemantic_stack_[(10) - (1)].tptr)->value.fnctptr_dddd))((yysemantic_stack_[(10) - (3)].val), (yysemantic_stack_[(10) - (5)].val), (yysemantic_stack_[(10) - (7)].val), (yysemantic_stack_[(10) - (9)].val));
	    else
	      (yyval.val) = 0.0;
	  }
    break;

  case 112:
/* Line 664 of lalr1.cc  */
#line 429 "aprepro.yy"
    {
	    if (arg_check((yysemantic_stack_[(12) - (1)].tptr), (yysemantic_stack_[(12) - (1)].tptr)->value.fnctptr_ddddc == NULL))
	      (yyval.val) = (*((yysemantic_stack_[(12) - (1)].tptr)->value.fnctptr_ddddc))((yysemantic_stack_[(12) - (3)].val), (yysemantic_stack_[(12) - (5)].val), (yysemantic_stack_[(12) - (7)].val), (yysemantic_stack_[(12) - (9)].val), (yysemantic_stack_[(12) - (11)].string));
	    else
	      (yyval.val) = 0.0;
	  }
    break;

  case 113:
/* Line 664 of lalr1.cc  */
#line 435 "aprepro.yy"
    {
	    if (arg_check((yysemantic_stack_[(14) - (1)].tptr), (yysemantic_stack_[(14) - (1)].tptr)->value.fnctptr_dddddd == NULL))
	      (yyval.val) = (*((yysemantic_stack_[(14) - (1)].tptr)->value.fnctptr_dddddd))((yysemantic_stack_[(14) - (3)].val), (yysemantic_stack_[(14) - (5)].val), (yysemantic_stack_[(14) - (7)].val), (yysemantic_stack_[(14) - (9)].val), (yysemantic_stack_[(14) - (11)].val), (yysemantic_stack_[(14) - (13)].val));
	    else
	      (yyval.val) = 0.0;
	  }
    break;

  case 114:
/* Line 664 of lalr1.cc  */
#line 441 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) + (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 115:
/* Line 664 of lalr1.cc  */
#line 442 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) - (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 116:
/* Line 664 of lalr1.cc  */
#line 443 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) * (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 117:
/* Line 664 of lalr1.cc  */
#line 444 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (3)].val) == 0.)
				    {
				      yyerror(aprepro, "Zero divisor"); 
				      yyerrok;
				    }
				  else
				    (yyval.val) = (yysemantic_stack_[(3) - (1)].val) / (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 118:
/* Line 664 of lalr1.cc  */
#line 451 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (3)].val) == 0.)
				    {
				      yyerror(aprepro, "Zero divisor");
				      yyerrok;
				    }
				  else
				    (yyval.val) = (int)(yysemantic_stack_[(3) - (1)].val) % (int)(yysemantic_stack_[(3) - (3)].val);		}
    break;

  case 119:
/* Line 664 of lalr1.cc  */
#line 458 "aprepro.yy"
    { (yyval.val) = -(yysemantic_stack_[(2) - (2)].val);				}
    break;

  case 120:
/* Line 664 of lalr1.cc  */
#line 459 "aprepro.yy"
    { (yyval.val) =  (yysemantic_stack_[(2) - (2)].val);				}
    break;

  case 121:
/* Line 664 of lalr1.cc  */
#line 460 "aprepro.yy"
    { errno = 0;
				  (yyval.val) = std::pow((yysemantic_stack_[(3) - (1)].val), (yysemantic_stack_[(3) - (3)].val)); 
				  SEAMS::math_error(aprepro, "Power");			}
    break;

  case 122:
/* Line 664 of lalr1.cc  */
#line 463 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (2)].val);				}
    break;

  case 123:
/* Line 664 of lalr1.cc  */
#line 464 "aprepro.yy"
    { errno = 0;
				  (yyval.val) = (double)((yysemantic_stack_[(3) - (2)].val) < 0 ? -floor(-((yysemantic_stack_[(3) - (2)].val))): floor((yysemantic_stack_[(3) - (2)].val)) );
				  SEAMS::math_error(aprepro, "floor (int)");		}
    break;

  case 124:
/* Line 664 of lalr1.cc  */
#line 467 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(1) - (1)].val)) ? 1 : 0; }
    break;

  case 125:
/* Line 664 of lalr1.cc  */
#line 468 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(5) - (1)].val)) ? ((yysemantic_stack_[(5) - (3)].val)) : ((yysemantic_stack_[(5) - (5)].val));              }
    break;

  case 126:
/* Line 664 of lalr1.cc  */
#line 469 "aprepro.yy"
    { (yyval.val) = array_value((yysemantic_stack_[(6) - (1)].tptr)->value.avar, (yysemantic_stack_[(6) - (3)].val), (yysemantic_stack_[(6) - (5)].val)); }
    break;

  case 127:
/* Line 664 of lalr1.cc  */
#line 471 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(8) - (8)].val);
                                    array *arr = (yysemantic_stack_[(8) - (1)].tptr)->value.avar;
                                    int cols = arr->cols;
                                    int rows = arr->rows;
				    int row = (yysemantic_stack_[(8) - (3)].val);
				    int col = (yysemantic_stack_[(8) - (5)].val);
				    if (aprepro.ap_options.one_based_index) {
				      row--;
				      col--;
				    }
				    if (row < rows && col < cols) {
                                      int offset = row*cols+col;
                                      (yysemantic_stack_[(8) - (1)].tptr)->value.avar->data[offset] = (yysemantic_stack_[(8) - (8)].val);
                                    }
                                    else {
                                      yyerror(aprepro, "Row or Column index out of range"); 
                                      yyerrok;
                                    }
                                  }
    break;


/* Line 664 of lalr1.cc  */
#line 1463 "apr_parser.cc"
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
  const signed char Parser::yypact_ninf_ = -39;
  const short int
  Parser::yypact_[] =
  {
       -39,    22,   -39,    -9,   258,   -39,   -39,   -39,   -39,   -39,
     213,   403,     2,   565,     4,    52,    18,    27,    36,   337,
     337,   -39,   357,   337,   337,   141,   186,    49,   268,    91,
    1060,   357,   337,   337,   337,   337,   337,   -39,   -39,   337,
     337,   337,   337,   337,   337,   -39,   -39,   337,   337,   337,
     337,   337,   337,   337,   -39,   -39,   337,   337,   357,   206,
     312,   357,   595,    41,   337,   103,  1090,   798,  1012,   -39,
      14,    14,    14,   -39,   -39,   -39,   -39,   -39,   -39,   -39,
     -39,   337,   337,   337,   -39,   357,   357,   357,   337,   -39,
     337,   337,   337,   337,   337,   337,   337,   -39,   337,   337,
     337,   337,   337,   337,   337,   337,   337,   337,   337,   357,
     337,   337,   165,  1090,  1109,  1125,  1125,  1125,  1125,  1125,
    1090,  1125,  1125,  1125,  1125,  1125,  1125,  1090,  1125,  1090,
    1125,  1125,  1125,  1125,  1125,  1125,  1090,  1125,   713,   165,
    1109,   -39,    34,   121,   564,   -39,    76,   148,   594,   129,
     373,   624,   337,    14,   -39,   -39,   337,   -39,   -27,  1076,
      37,  1125,   -39,   -38,   -38,   -39,  1141,  1141,    19,    19,
      19,    19,    19,    19,   -39,  1156,  1170,   253,   253,   253,
     253,   253,   253,   -28,   -28,    14,   -39,    14,    14,    14,
     337,   -39,   337,   -39,   337,   -39,   -39,   337,   -39,   337,
     -39,   -39,   337,   -39,   337,   -39,  1125,    14,   337,   337,
    1037,   383,   825,   454,   534,   400,   425,   852,   460,   879,
     906,  1090,  1125,    48,   337,   -39,   -39,   -39,   337,   -39,
     337,   337,   337,   -39,   -39,   -39,   -39,   337,   933,   654,
     742,   482,   408,  1125,   -39,   337,   -39,   337,   -39,   337,
     -39,   684,   960,   432,   337,   -39,   -39,   337,   488,   771,
     510,   -39,   337,   -39,   987,   -39
  };

  /* YYDEFACT[S] -- default reduction number in state S.  Performed when
     YYTABLE doesn't specify something else to do.  Zero means the
     default is an error.  */
  const unsigned char
  Parser::yydefact_[] =
  {
         2,     0,     1,     0,     0,     4,     3,     9,    62,    44,
      89,    65,    45,    66,    46,    28,     0,     0,     0,     0,
       0,     8,     0,     0,     0,     0,     0,   124,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    93,    92,     0,
       0,     0,     0,     0,     0,    70,    69,     0,     0,     0,
       0,     0,     0,     0,    81,    80,     0,     0,     0,     0,
       0,     0,    89,     0,     0,   124,     0,     0,     0,    38,
     119,   120,    12,    64,    91,    68,    79,    63,    90,    67,
      78,     0,     0,     0,     7,     0,     0,     0,     0,     6,
       0,     0,     0,     0,     0,     0,     0,     5,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    36,    47,    94,    96,    95,    98,    97,    99,
      49,    71,    74,    73,    76,    75,    77,    48,    72,    51,
      82,    85,    84,    87,    86,    88,    50,    83,     0,    35,
       0,   100,     0,     0,     0,    53,     0,     0,     0,     0,
       0,     0,     0,   119,    21,   122,     0,   123,     0,     0,
      19,     0,    20,    39,    37,    43,    40,    41,    27,    26,
      25,    24,    23,    22,    56,    17,    18,    16,    15,    14,
      13,    11,    10,   115,   114,   118,    42,   116,   117,   121,
       0,   103,     0,   102,     0,   101,    55,     0,    52,     0,
      54,    34,     0,    29,     0,    33,    94,   116,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    61,   125,   126,     0,   106,   104,   105,     0,   107,
       0,     0,     0,    57,    31,    30,    32,     0,     0,     0,
       0,     0,     0,   127,   109,     0,   108,     0,    60,     0,
      59,     0,     0,     0,     0,   111,   110,     0,     0,     0,
       0,   112,     0,    58,     0,   113
  };

  /* YYPGOTO[NTERM-NUM].  */
  const signed char
  Parser::yypgoto_[] =
  {
       -39,   -39,   -39,   -17,     3,    82,    -4
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  Parser::yydefgoto_[] =
  {
        -1,     1,     6,    27,    28,    66,   161
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If YYTABLE_NINF_, syntax error.  */
  const signed char Parser::yytable_ninf_ = -1;
  const unsigned short int
  Parser::yytable_[] =
  {
        30,   208,    65,    87,    88,    90,    91,    92,    93,    94,
      95,     7,   108,   156,   110,    67,    68,   111,    70,    71,
      72,    96,     2,     3,    47,    69,    56,   114,   115,   116,
     117,   118,   119,    59,   112,   121,   122,   123,   124,   125,
     126,     4,    60,   128,   130,   131,   132,   133,   134,   135,
     191,    61,   137,   138,   140,   144,   148,   151,    57,   111,
     153,   139,   142,   146,   149,   160,   162,    96,    83,    57,
     237,     5,    85,    86,    58,    87,    88,   159,    81,    82,
      83,   140,   140,   166,   167,     0,    29,     0,   163,   164,
     165,     0,   196,     0,   175,   176,   177,   178,   179,   180,
     181,   182,   183,   184,   185,   187,   188,   189,     0,     0,
       0,    89,   186,   113,    85,    86,     0,    87,    88,   154,
       0,   120,     0,    90,    91,    92,    93,    94,    95,   127,
     129,     0,    81,    82,    83,   192,     0,   193,   136,    96,
       0,   143,   147,   150,    73,   201,    74,    75,   206,    76,
       0,     0,   207,    90,    91,    92,    93,    94,    95,     0,
       0,     0,   197,   158,   198,     0,     0,    85,    86,    96,
      87,    88,   168,   169,   170,   171,   172,   173,   174,     0,
      90,    91,    92,    93,    94,    95,   210,     0,   212,    77,
     214,    78,    79,     0,    80,   217,    96,     0,   219,     0,
     220,     0,     0,    85,    86,   222,    87,    88,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
     238,    19,   141,    20,   239,     0,   240,     0,     0,     0,
       0,     0,     0,   243,   113,    31,    32,    33,    34,    35,
      36,   251,     0,   252,    22,    23,     0,     0,     0,    24,
     259,     0,    25,    26,     0,     0,     0,     0,   264,    37,
      38,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,     0,    19,   211,    20,   213,     0,    21,   215,
       0,   216,     0,     0,   218,     0,     0,     0,    84,     0,
     221,   106,   107,   108,   156,   110,    22,    23,   111,     0,
       0,    24,     0,     0,    25,    26,    85,    86,     0,    87,
      88,     0,     0,   241,   242,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,     0,    19,   145,    20,
       0,   253,     0,     0,     0,     0,   258,     0,     0,   260,
       8,     9,    62,    11,    12,    13,    14,    63,    16,    17,
      22,    23,    19,     0,    20,    24,     0,     0,    25,    26,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,     0,    19,     0,    20,    64,    23,     0,     0,     0,
      24,     0,     0,    25,    26,     0,     0,   202,     0,   203,
       0,     0,     0,     0,     0,    22,    23,   224,     0,   225,
      24,     0,     0,    25,    26,    90,    91,    92,    93,    94,
      95,     0,     0,     0,   231,    90,    91,    92,    93,    94,
      95,    96,   249,     0,   250,    39,    40,    41,    42,    43,
      44,    96,    90,    91,    92,    93,    94,    95,     0,   232,
      90,    91,    92,    93,    94,    95,   257,     0,    96,    45,
      46,     0,     0,     0,     0,     0,    96,    90,    91,    92,
      93,    94,    95,     0,    90,    91,    92,    93,    94,    95,
     227,     0,     0,    96,     0,     0,   234,     0,     0,     0,
      96,     0,     0,     0,     0,     0,    90,    91,    92,    93,
      94,    95,    90,    91,    92,    93,    94,    95,   248,     0,
       0,     0,    96,     0,   261,     0,     0,     0,    96,     0,
       0,     0,     0,     0,    90,    91,    92,    93,    94,    95,
      90,    91,    92,    93,    94,    95,   263,     0,     0,     0,
      96,     0,     0,     0,     0,     0,    96,     0,     0,     0,
       0,     0,    90,    91,    92,    93,    94,    95,   228,     0,
     229,     0,     0,     0,     0,   230,     0,     0,    96,     0,
       0,     0,     0,     0,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   156,   110,     0,   194,   111,
     195,     0,     0,     0,     0,     0,     0,    48,    49,    50,
      51,    52,    53,     0,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,     0,   199,   111,
     200,    54,    55,     0,     0,     0,     0,   152,    32,    33,
      34,    35,    36,     0,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,     0,   204,   111,
     205,    37,    38,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,     0,   245,   111,
     246,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   156,   110,     0,   254,   111,
     255,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   156,   110,   190,     0,   111,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   156,   110,   247,     0,   111,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   156,   110,   262,     0,   111,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     107,   108,   156,   110,   155,     0,   111,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,   156,
     110,   226,     0,   111,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    98,    99,   100,   101,   102,
     103,   104,   105,   106,   107,   108,   156,   110,   233,     0,
     111,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   156,   110,   235,     0,   111,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     156,   110,   236,     0,   111,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   156,   110,   244,
       0,   111,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   156,   110,   256,     0,   111,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      98,    99,   100,   101,   102,   103,   104,   105,   106,   107,
     108,   156,   110,   265,     0,   111,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    98,    99,   100,
     101,   102,   103,   104,   105,   106,   107,   108,   156,   110,
     157,     0,   111,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   156,   110,   223,     0,   111,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    98,    99,   100,
     101,   102,   103,   104,   105,   106,   107,   108,   156,   110,
      97,     0,   111,     0,     0,     0,     0,     0,     0,     0,
      98,    99,   100,   101,   102,   103,   104,   105,   106,   107,
     108,   109,   110,     0,   209,   111,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   156,   110,     0,
       0,   111,    90,    91,    92,    93,    94,    95,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    96,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     109,   110,     0,     0,   111,    98,    99,   100,   101,   102,
     103,   104,   105,   106,   107,   108,   156,   110,     0,     0,
     111,    98,    99,   100,   101,   102,   103,   104,   105,     0,
       0,   108,     0,     0,     0,     0,   111,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   156,   110,     0,
       0,   111,   100,   101,   102,   103,   104,   105,   106,   107,
     108,   156,   110,     0,     0,   111
  };

  /* YYCHECK.  */
  const short int
  Parser::yycheck_[] =
  {
         4,    28,    19,    41,    42,    32,    33,    34,    35,    36,
      37,    20,    40,    41,    42,    19,    20,    45,    22,    23,
      24,    48,     0,     1,    22,    22,    22,    31,    32,    33,
      34,    35,    36,    15,    31,    39,    40,    41,    42,    43,
      44,    19,    15,    47,    48,    49,    50,    51,    52,    53,
      16,    15,    56,    57,    58,    59,    60,    61,    17,    45,
      64,    58,    59,    60,    61,    82,    83,    48,    31,    17,
      22,    49,    38,    39,    22,    41,    42,    81,    29,    30,
      31,    85,    86,    87,    88,    -1,     4,    -1,    85,    86,
      87,    -1,    16,    -1,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,   111,    -1,    -1,
      -1,    20,   109,    31,    38,    39,    -1,    41,    42,    16,
      -1,    39,    -1,    32,    33,    34,    35,    36,    37,    47,
      48,    -1,    29,    30,    31,    14,    -1,    16,    56,    48,
      -1,    59,    60,    61,     3,    16,     5,     6,   152,     8,
      -1,    -1,   156,    32,    33,    34,    35,    36,    37,    -1,
      -1,    -1,    14,    81,    16,    -1,    -1,    38,    39,    48,
      41,    42,    90,    91,    92,    93,    94,    95,    96,    -1,
      32,    33,    34,    35,    36,    37,   190,    -1,   192,     3,
     194,     5,     6,    -1,     8,   199,    48,    -1,   202,    -1,
     204,    -1,    -1,    38,    39,   209,    41,    42,    -1,     3,
       4,     5,     6,     7,     8,     9,    10,    11,    12,    13,
     224,    15,    16,    17,   228,    -1,   230,    -1,    -1,    -1,
      -1,    -1,    -1,   237,   152,    22,    23,    24,    25,    26,
      27,   245,    -1,   247,    38,    39,    -1,    -1,    -1,    43,
     254,    -1,    46,    47,    -1,    -1,    -1,    -1,   262,    46,
      47,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    -1,    15,   192,    17,   194,    -1,    20,   197,
      -1,   199,    -1,    -1,   202,    -1,    -1,    -1,    20,    -1,
     208,    38,    39,    40,    41,    42,    38,    39,    45,    -1,
      -1,    43,    -1,    -1,    46,    47,    38,    39,    -1,    41,
      42,    -1,    -1,   231,   232,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    13,    -1,    15,    16,    17,
      -1,   249,    -1,    -1,    -1,    -1,   254,    -1,    -1,   257,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      38,    39,    15,    -1,    17,    43,    -1,    -1,    46,    47,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    -1,    15,    -1,    17,    38,    39,    -1,    -1,    -1,
      43,    -1,    -1,    46,    47,    -1,    -1,    14,    -1,    16,
      -1,    -1,    -1,    -1,    -1,    38,    39,    14,    -1,    16,
      43,    -1,    -1,    46,    47,    32,    33,    34,    35,    36,
      37,    -1,    -1,    -1,    14,    32,    33,    34,    35,    36,
      37,    48,    14,    -1,    16,    22,    23,    24,    25,    26,
      27,    48,    32,    33,    34,    35,    36,    37,    -1,    14,
      32,    33,    34,    35,    36,    37,    14,    -1,    48,    46,
      47,    -1,    -1,    -1,    -1,    -1,    48,    32,    33,    34,
      35,    36,    37,    -1,    32,    33,    34,    35,    36,    37,
      16,    -1,    -1,    48,    -1,    -1,    16,    -1,    -1,    -1,
      48,    -1,    -1,    -1,    -1,    -1,    32,    33,    34,    35,
      36,    37,    32,    33,    34,    35,    36,    37,    16,    -1,
      -1,    -1,    48,    -1,    16,    -1,    -1,    -1,    48,    -1,
      -1,    -1,    -1,    -1,    32,    33,    34,    35,    36,    37,
      32,    33,    34,    35,    36,    37,    16,    -1,    -1,    -1,
      48,    -1,    -1,    -1,    -1,    -1,    48,    -1,    -1,    -1,
      -1,    -1,    32,    33,    34,    35,    36,    37,    14,    -1,
      16,    -1,    -1,    -1,    -1,    21,    -1,    -1,    48,    -1,
      -1,    -1,    -1,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    -1,    14,    45,
      16,    -1,    -1,    -1,    -1,    -1,    -1,    22,    23,    24,
      25,    26,    27,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    -1,    14,    45,
      16,    46,    47,    -1,    -1,    -1,    -1,    22,    23,    24,
      25,    26,    27,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    -1,    14,    45,
      16,    46,    47,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    -1,    14,    45,
      16,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    -1,    14,    45,
      16,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    14,    -1,    45,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    14,    -1,    45,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    14,    -1,    45,    -1,    -1,
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
      41,    42,    16,    -1,    45,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    16,
      -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    16,    -1,    45,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    16,    -1,    45,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      18,    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    18,    -1,    45,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      20,    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    -1,    28,    45,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    -1,
      -1,    45,    32,    33,    34,    35,    36,    37,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    48,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    -1,    -1,    45,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    -1,    -1,
      45,    30,    31,    32,    33,    34,    35,    36,    37,    -1,
      -1,    40,    -1,    -1,    -1,    -1,    45,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    -1,
      -1,    45,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    -1,    -1,    45
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned char
  Parser::yystos_[] =
  {
         0,    51,     0,     1,    19,    49,    52,    20,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    15,
      17,    20,    38,    39,    43,    46,    47,    53,    54,    55,
      56,    22,    23,    24,    25,    26,    27,    46,    47,    22,
      23,    24,    25,    26,    27,    46,    47,    22,    22,    23,
      24,    25,    26,    27,    46,    47,    22,    17,    22,    15,
      15,    15,     5,    10,    38,    53,    55,    56,    56,    54,
      56,    56,    56,     3,     5,     6,     8,     3,     5,     6,
       8,    29,    30,    31,    20,    38,    39,    41,    42,    20,
      32,    33,    34,    35,    36,    37,    48,    20,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    45,    54,    55,    56,    56,    56,    56,    56,    56,
      55,    56,    56,    56,    56,    56,    56,    55,    56,    55,
      56,    56,    56,    56,    56,    56,    55,    56,    56,    54,
      56,    16,    54,    55,    56,    16,    54,    55,    56,    54,
      55,    56,    22,    56,    16,    16,    41,    18,    55,    56,
      53,    56,    53,    54,    54,    54,    56,    56,    55,    55,
      55,    55,    55,    55,    55,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    54,    56,    56,    56,
      14,    16,    14,    16,    14,    16,    16,    14,    16,    14,
      16,    16,    14,    16,    14,    16,    56,    56,    28,    28,
      56,    55,    56,    55,    56,    55,    55,    56,    55,    56,
      56,    55,    56,    18,    14,    16,    16,    16,    14,    16,
      21,    14,    14,    16,    16,    16,    16,    22,    56,    56,
      56,    55,    55,    56,    16,    14,    16,    14,    16,    14,
      16,    56,    56,    55,    14,    16,    16,    14,    55,    56,
      55,    16,    14,    16,    56,    16
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
     295,   296,   297,   298,   299,   300,   301,   302,   303,    10
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned char
  Parser::yyr1_[] =
  {
         0,    50,    51,    51,    52,    52,    52,    52,    52,    52,
      53,    53,    53,    53,    53,    53,    53,    53,    53,    53,
      53,    53,    53,    53,    53,    53,    53,    53,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    55,    55,    55,    55,    55,    55,
      55,    55,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    55,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  Parser::yyr2_[] =
  {
         0,     2,     0,     2,     1,     3,     3,     3,     2,     2,
       3,     3,     2,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     1,     4,
       6,     6,     6,     4,     4,     3,     3,     3,     2,     3,
       3,     3,     3,     3,     1,     1,     1,     3,     3,     3,
       3,     3,     4,     3,     4,     4,     3,     6,    12,     8,
       8,     5,     1,     2,     2,     1,     1,     2,     2,     2,
       2,     3,     3,     3,     3,     3,     3,     3,     2,     2,
       2,     2,     3,     3,     3,     3,     3,     3,     3,     1,
       2,     2,     2,     2,     3,     3,     3,     3,     3,     3,
       3,     4,     4,     4,     6,     6,     6,     6,     8,     8,
      10,    10,    12,    14,     3,     3,     3,     3,     3,     2,
       2,     3,     3,     3,     1,     5,     6,     8
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
  "'\\n'", "$accept", "input", "line", "bool", "aexp", "sexp", "exp", YY_NULL
  };

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const Parser::rhs_number_type
  Parser::yyrhs_[] =
  {
        51,     0,    -1,    -1,    51,    52,    -1,    49,    -1,    19,
      56,    20,    -1,    19,    55,    20,    -1,    19,    54,    20,
      -1,    19,    20,    -1,     1,    20,    -1,    56,    37,    56,
      -1,    56,    36,    56,    -1,    43,    56,    -1,    56,    35,
      56,    -1,    56,    34,    56,    -1,    56,    33,    56,    -1,
      56,    32,    56,    -1,    56,    30,    56,    -1,    56,    31,
      56,    -1,    53,    30,    53,    -1,    53,    31,    53,    -1,
      15,    53,    16,    -1,    55,    37,    55,    -1,    55,    36,
      55,    -1,    55,    35,    55,    -1,    55,    34,    55,    -1,
      55,    33,    55,    -1,    55,    32,    55,    -1,    10,    -1,
      13,    15,    55,    16,    -1,    13,    15,    55,    14,    56,
      16,    -1,    13,    15,    55,    14,    55,    16,    -1,    13,
      15,    56,    14,    56,    16,    -1,    13,    15,    56,    16,
      -1,    13,    15,    54,    16,    -1,    10,    22,    54,    -1,
       5,    22,    54,    -1,    54,    39,    54,    -1,    38,    54,
      -1,    54,    38,    54,    -1,    54,    41,    56,    -1,    54,
      42,    56,    -1,    56,    41,    54,    -1,    54,    41,    54,
      -1,     4,    -1,     7,    -1,     9,    -1,     5,    22,    55,
      -1,     7,    22,    55,    -1,     6,    22,    55,    -1,     9,
      22,    55,    -1,     8,    22,    55,    -1,    12,    15,    55,
      16,    -1,    12,    15,    16,    -1,    12,    15,    56,    16,
      -1,    12,    15,    54,    16,    -1,    55,    48,    55,    -1,
      12,    15,    56,    14,    56,    16,    -1,    12,    15,    56,
      14,    55,    14,    55,    14,    55,    14,    55,    16,    -1,
      12,    15,    56,    14,    55,    14,    55,    16,    -1,    12,
      15,    55,    14,    55,    14,    55,    16,    -1,    53,    29,
      55,    28,    55,    -1,     3,    -1,    47,     3,    -1,    46,
       3,    -1,     6,    -1,     8,    -1,    47,     6,    -1,    46,
       6,    -1,     6,    47,    -1,     6,    46,    -1,     6,    22,
      56,    -1,     7,    22,    56,    -1,     6,    24,    56,    -1,
       6,    23,    56,    -1,     6,    26,    56,    -1,     6,    25,
      56,    -1,     6,    27,    56,    -1,    47,     8,    -1,    46,
       8,    -1,     8,    47,    -1,     8,    46,    -1,     8,    22,
      56,    -1,     9,    22,    56,    -1,     8,    24,    56,    -1,
       8,    23,    56,    -1,     8,    26,    56,    -1,     8,    25,
      56,    -1,     8,    27,    56,    -1,     5,    -1,    47,     5,
      -1,    46,     5,    -1,     5,    47,    -1,     5,    46,    -1,
       5,    22,    56,    -1,     5,    24,    56,    -1,     5,    23,
      56,    -1,     5,    26,    56,    -1,     5,    25,    56,    -1,
       5,    27,    56,    -1,    11,    15,    16,    -1,    11,    15,
      56,    16,    -1,    11,    15,    55,    16,    -1,    11,    15,
      54,    16,    -1,    11,    15,    55,    14,    56,    16,    -1,
      11,    15,    56,    14,    55,    16,    -1,    11,    15,    55,
      14,    55,    16,    -1,    11,    15,    56,    14,    56,    16,
      -1,    11,    15,    56,    14,    56,    14,    56,    16,    -1,
      11,    15,    55,    14,    55,    14,    56,    16,    -1,    11,
      15,    56,    14,    56,    21,    56,    14,    56,    16,    -1,
      11,    15,    56,    14,    56,    14,    56,    14,    56,    16,
      -1,    11,    15,    56,    14,    56,    14,    56,    14,    56,
      14,    55,    16,    -1,    11,    15,    56,    14,    56,    14,
      56,    14,    56,    14,    56,    14,    56,    16,    -1,    56,
      39,    56,    -1,    56,    38,    56,    -1,    56,    41,    56,
      -1,    56,    42,    56,    -1,    56,    40,    56,    -1,    38,
      56,    -1,    39,    56,    -1,    56,    45,    56,    -1,    15,
      56,    16,    -1,    17,    56,    18,    -1,    53,    -1,    53,
      29,    56,    28,    56,    -1,    10,    17,    56,    14,    56,
      18,    -1,    10,    17,    56,    14,    56,    18,    22,    56,
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
     105,   112,   119,   126,   131,   136,   140,   144,   148,   151,
     155,   159,   163,   167,   171,   173,   175,   177,   181,   185,
     189,   193,   197,   202,   206,   211,   216,   220,   227,   240,
     249,   258,   264,   266,   269,   272,   274,   276,   279,   282,
     285,   288,   292,   296,   300,   304,   308,   312,   316,   319,
     322,   325,   328,   332,   336,   340,   344,   348,   352,   356,
     358,   361,   364,   367,   370,   374,   378,   382,   386,   390,
     394,   398,   403,   408,   413,   420,   427,   434,   441,   450,
     459,   470,   481,   494,   509,   513,   517,   521,   525,   529,
     532,   535,   539,   543,   547,   549,   555,   562
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  Parser::yyrline_[] =
  {
         0,   101,   101,   102,   105,   106,   113,   117,   118,   119,
     122,   123,   124,   125,   126,   127,   128,   129,   130,   131,
     132,   133,   136,   137,   138,   139,   140,   141,   143,   144,
     150,   156,   162,   168,   174,   180,   183,   185,   193,   195,
     203,   204,   205,   206,   215,   216,   217,   218,   220,   223,
     227,   228,   229,   235,   241,   247,   253,   254,   260,   266,
     272,   278,   280,   281,   282,   283,   284,   285,   286,   287,
     288,   289,   291,   294,   295,   296,   297,   298,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   315,
     317,   320,   323,   326,   329,   331,   334,   337,   340,   343,
     350,   357,   364,   371,   378,   385,   392,   399,   405,   411,
     417,   423,   429,   435,   441,   442,   443,   444,   451,   458,
     459,   460,   463,   464,   467,   468,   469,   470
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
  const int Parser::yylast_ = 1215;
  const int Parser::yynnts_ = 7;
  const int Parser::yyempty_ = -2;
  const int Parser::yyfinal_ = 2;
  const int Parser::yyterror_ = 1;
  const int Parser::yyerrcode_ = 256;
  const int Parser::yyntokens_ = 50;

  const unsigned int Parser::yyuser_token_number_max_ = 303;
  const Parser::token_number_type Parser::yyundef_token_ = 2;


} // SEAMS
/* Line 1135 of lalr1.cc  */
#line 2399 "apr_parser.cc"
/* Line 1136 of lalr1.cc  */
#line 493 "aprepro.yy"


void SEAMS::Parser::error(const Parser::location_type&, const std::string& m)
{
    aprepro.error(m);
}

