/* A Bison parser, made by GNU Bison 2.4.3.  */

/* Skeleton implementation for Bison LALR(1) parsers in C++
   
      Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 Free
   Software Foundation, Inc.
   
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

/* Line 311 of lalr1.cc  */
#line 1 "aprepro.yy"

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



/* Line 311 of lalr1.cc  */
#line 59 "apr_parser.cc"


#include "aprepro_parser.h"

/* User implementation prologue.  */

/* Line 317 of lalr1.cc  */
#line 78 "aprepro.yy"


#include "aprepro.h"
#include "apr_scanner.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex



/* Line 317 of lalr1.cc  */
#line 82 "apr_parser.cc"

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

/* Line 380 of lalr1.cc  */
#line 148 "apr_parser.cc"
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
    location_type yyerror_range[3];

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

/* Line 678 of lalr1.cc  */
#line 97 "aprepro.yy"
    { if (echo) aprepro.lexer->LexerOutput("\n", 1); }
    break;

  case 5:

/* Line 678 of lalr1.cc  */
#line 98 "aprepro.yy"
    { if (echo) {
	                             static char tmpstr[512];
				     SEAMS::symrec *format = aprepro.getsym("_FORMAT");
				     int len = sprintf(tmpstr, format->value.svar, (yysemantic_stack_[(3) - (2)].val));
				     aprepro.lexer->LexerOutput(tmpstr, len);
				   }
                                }
    break;

  case 6:

/* Line 678 of lalr1.cc  */
#line 105 "aprepro.yy"
    { if (echo && (yysemantic_stack_[(3) - (2)].string) != NULL) {
				    aprepro.lexer->LexerOutput((yysemantic_stack_[(3) - (2)].string), strlen((yysemantic_stack_[(3) - (2)].string)));
                                  }
                                }
    break;

  case 7:

/* Line 678 of lalr1.cc  */
#line 109 "aprepro.yy"
    {                                       }
    break;

  case 8:

/* Line 678 of lalr1.cc  */
#line 110 "aprepro.yy"
    { yyerrok;				}
    break;

  case 9:

/* Line 678 of lalr1.cc  */
#line 113 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) < (yysemantic_stack_[(3) - (3)].val);                         }
    break;

  case 10:

/* Line 678 of lalr1.cc  */
#line 114 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) > (yysemantic_stack_[(3) - (3)].val);                         }
    break;

  case 11:

/* Line 678 of lalr1.cc  */
#line 115 "aprepro.yy"
    { (yyval.val) = !((yysemantic_stack_[(2) - (2)].val));                           }
    break;

  case 12:

/* Line 678 of lalr1.cc  */
#line 116 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) <= (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 13:

/* Line 678 of lalr1.cc  */
#line 117 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) >= (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 14:

/* Line 678 of lalr1.cc  */
#line 118 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) == (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 15:

/* Line 678 of lalr1.cc  */
#line 119 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) != (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 16:

/* Line 678 of lalr1.cc  */
#line 120 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) || (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 17:

/* Line 678 of lalr1.cc  */
#line 121 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) && (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 18:

/* Line 678 of lalr1.cc  */
#line 122 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) || (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 19:

/* Line 678 of lalr1.cc  */
#line 123 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) && (yysemantic_stack_[(3) - (3)].val);                        }
    break;

  case 20:

/* Line 678 of lalr1.cc  */
#line 124 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (2)].val);                              }
    break;

  case 21:

/* Line 678 of lalr1.cc  */
#line 127 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) <  0 ? 1 : 0);	}
    break;

  case 22:

/* Line 678 of lalr1.cc  */
#line 128 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) >  0 ? 1 : 0);	}
    break;

  case 23:

/* Line 678 of lalr1.cc  */
#line 129 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) <= 0 ? 1 : 0);	}
    break;

  case 24:

/* Line 678 of lalr1.cc  */
#line 130 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) >= 0 ? 1 : 0);	}
    break;

  case 25:

/* Line 678 of lalr1.cc  */
#line 131 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) == 0 ? 1 : 0);	}
    break;

  case 26:

/* Line 678 of lalr1.cc  */
#line 132 "aprepro.yy"
    { (yyval.val) = (strcmp((yysemantic_stack_[(3) - (1)].string),(yysemantic_stack_[(3) - (3)].string)) != 0 ? 1 : 0);	}
    break;

  case 27:

/* Line 678 of lalr1.cc  */
#line 134 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(1) - (1)].string);				}
    break;

  case 28:

/* Line 678 of lalr1.cc  */
#line 135 "aprepro.yy"
    { (yyval.string) = (char*)(yysemantic_stack_[(1) - (1)].tptr)->value.svar;			}
    break;

  case 29:

/* Line 678 of lalr1.cc  */
#line 136 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); (yysemantic_stack_[(3) - (1)].tptr)->value.svar = (yysemantic_stack_[(3) - (3)].string);
				  (yysemantic_stack_[(3) - (1)].tptr)->type = token::SVAR;			}
    break;

  case 30:

/* Line 678 of lalr1.cc  */
#line 138 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); 
				  (yysemantic_stack_[(3) - (1)].tptr)->value.svar = (yysemantic_stack_[(3) - (3)].string);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 31:

/* Line 678 of lalr1.cc  */
#line 141 "aprepro.yy"
    { (yyval.string) = (yysemantic_stack_[(3) - (3)].string); 
				  (yysemantic_stack_[(3) - (1)].tptr)->value.svar= (yysemantic_stack_[(3) - (3)].string);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          
				  (yysemantic_stack_[(3) - (1)].tptr)->type = token::SVAR; 		}
    break;

  case 32:

/* Line 678 of lalr1.cc  */
#line 145 "aprepro.yy"
    { (yyval.string) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_c))((yysemantic_stack_[(4) - (3)].string));	}
    break;

  case 33:

/* Line 678 of lalr1.cc  */
#line 146 "aprepro.yy"
    { (yyval.string) = (*((yysemantic_stack_[(3) - (1)].tptr)->value.strfnct))();	}
    break;

  case 34:

/* Line 678 of lalr1.cc  */
#line 147 "aprepro.yy"
    { (yyval.string) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.strfnct_d))((yysemantic_stack_[(4) - (3)].val));	}
    break;

  case 35:

/* Line 678 of lalr1.cc  */
#line 148 "aprepro.yy"
    { int len1 = strlen((yysemantic_stack_[(3) - (1)].string));
				  int len3 = strlen((yysemantic_stack_[(3) - (3)].string));
				  (yyval.string) = (char*)calloc(1, (len1+len3+1));
				  memcpy((yyval.string), (yysemantic_stack_[(3) - (1)].string), len1+1);
				  (void *)strcat((yyval.string), (yysemantic_stack_[(3) - (3)].string)); }
    break;

  case 36:

/* Line 678 of lalr1.cc  */
#line 154 "aprepro.yy"
    { (yyval.string) = (*((yysemantic_stack_[(12) - (1)].tptr)->value.strfnct_dcccc))((yysemantic_stack_[(12) - (3)].val), (yysemantic_stack_[(12) - (5)].string), (yysemantic_stack_[(12) - (7)].string), (yysemantic_stack_[(12) - (9)].string), (yysemantic_stack_[(12) - (11)].string)); }
    break;

  case 37:

/* Line 678 of lalr1.cc  */
#line 156 "aprepro.yy"
    { (yyval.string) = (*((yysemantic_stack_[(8) - (1)].tptr)->value.strfnct_dcc))((yysemantic_stack_[(8) - (3)].val), (yysemantic_stack_[(8) - (5)].string), (yysemantic_stack_[(8) - (7)].string)); }
    break;

  case 38:

/* Line 678 of lalr1.cc  */
#line 158 "aprepro.yy"
    { (yyval.string) = (*((yysemantic_stack_[(8) - (1)].tptr)->value.strfnct_ccc))((yysemantic_stack_[(8) - (3)].string), (yysemantic_stack_[(8) - (5)].string), (yysemantic_stack_[(8) - (7)].string)); }
    break;

  case 39:

/* Line 678 of lalr1.cc  */
#line 159 "aprepro.yy"
    { (yyval.string) = ((yysemantic_stack_[(5) - (1)].val)) ? ((yysemantic_stack_[(5) - (3)].string)) : ((yysemantic_stack_[(5) - (5)].string));              }
    break;

  case 40:

/* Line 678 of lalr1.cc  */
#line 161 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].val); 				}
    break;

  case 41:

/* Line 678 of lalr1.cc  */
#line 162 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(2) - (2)].val) + 1;				}
    break;

  case 42:

/* Line 678 of lalr1.cc  */
#line 163 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(2) - (2)].val) - 1;				}
    break;

  case 43:

/* Line 678 of lalr1.cc  */
#line 164 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;			}
    break;

  case 44:

/* Line 678 of lalr1.cc  */
#line 165 "aprepro.yy"
    { (yyval.val) = ++((yysemantic_stack_[(2) - (2)].tptr)->value.var);		}
    break;

  case 45:

/* Line 678 of lalr1.cc  */
#line 166 "aprepro.yy"
    { (yyval.val) = --((yysemantic_stack_[(2) - (2)].tptr)->value.var);		}
    break;

  case 46:

/* Line 678 of lalr1.cc  */
#line 167 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)++;		}
    break;

  case 47:

/* Line 678 of lalr1.cc  */
#line 168 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)--;		}
    break;

  case 48:

/* Line 678 of lalr1.cc  */
#line 169 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 49:

/* Line 678 of lalr1.cc  */
#line 171 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
				  redefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          
				  (yysemantic_stack_[(3) - (1)].tptr)->type = token::VAR;			}
    break;

  case 50:

/* Line 678 of lalr1.cc  */
#line 174 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var += (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 51:

/* Line 678 of lalr1.cc  */
#line 175 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var -= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 52:

/* Line 678 of lalr1.cc  */
#line 176 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var *= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 53:

/* Line 678 of lalr1.cc  */
#line 177 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var /= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; }
    break;

  case 54:

/* Line 678 of lalr1.cc  */
#line 178 "aprepro.yy"
    { errno = 0;
				  (yysemantic_stack_[(3) - (1)].tptr)->value.var = std::pow((yysemantic_stack_[(3) - (1)].tptr)->value.var,(yysemantic_stack_[(3) - (3)].val)); 
				  (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
				  SEAMS::math_error(aprepro, "Power");
				}
    break;

  case 55:

/* Line 678 of lalr1.cc  */
#line 183 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(1) - (1)].tptr)->value.var;
				  undefined_warning(aprepro, (yysemantic_stack_[(1) - (1)].tptr)->name);          }
    break;

  case 56:

/* Line 678 of lalr1.cc  */
#line 185 "aprepro.yy"
    { (yyval.val) = ++((yysemantic_stack_[(2) - (2)].tptr)->value.var);		
				  (yysemantic_stack_[(2) - (2)].tptr)->type = token::VAR;                       
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (2)].tptr)->name);          }
    break;

  case 57:

/* Line 678 of lalr1.cc  */
#line 188 "aprepro.yy"
    { (yyval.val) = --((yysemantic_stack_[(2) - (2)].tptr)->value.var);		
				  (yysemantic_stack_[(2) - (2)].tptr)->type = token::VAR;                       
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (2)].tptr)->name);          }
    break;

  case 58:

/* Line 678 of lalr1.cc  */
#line 191 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)++;		
				  (yysemantic_stack_[(2) - (1)].tptr)->type = token::VAR;                       
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (1)].tptr)->name);          }
    break;

  case 59:

/* Line 678 of lalr1.cc  */
#line 194 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(2) - (1)].tptr)->value.var)--;		
				  (yysemantic_stack_[(2) - (1)].tptr)->type = token::VAR;                       
				  undefined_warning(aprepro, (yysemantic_stack_[(2) - (1)].tptr)->name);          }
    break;

  case 60:

/* Line 678 of lalr1.cc  */
#line 197 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (3)].val); (yysemantic_stack_[(3) - (1)].tptr)->value.var = (yysemantic_stack_[(3) - (3)].val);
				  (yysemantic_stack_[(3) - (1)].tptr)->type = token::VAR;                       }
    break;

  case 61:

/* Line 678 of lalr1.cc  */
#line 199 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var += (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
				  (yysemantic_stack_[(3) - (1)].tptr)->type = token::VAR;                       
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 62:

/* Line 678 of lalr1.cc  */
#line 202 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var -= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
				  (yysemantic_stack_[(3) - (1)].tptr)->type = token::VAR;                       
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 63:

/* Line 678 of lalr1.cc  */
#line 205 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var *= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
				  (yysemantic_stack_[(3) - (1)].tptr)->type = token::VAR;                       
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 64:

/* Line 678 of lalr1.cc  */
#line 208 "aprepro.yy"
    { (yysemantic_stack_[(3) - (1)].tptr)->value.var /= (yysemantic_stack_[(3) - (3)].val); (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
				  (yysemantic_stack_[(3) - (1)].tptr)->type = token::VAR;                       
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 65:

/* Line 678 of lalr1.cc  */
#line 211 "aprepro.yy"
    { errno = 0;
				  (yysemantic_stack_[(3) - (1)].tptr)->value.var = std::pow((yysemantic_stack_[(3) - (1)].tptr)->value.var,(yysemantic_stack_[(3) - (3)].val)); 
				  (yyval.val) = (yysemantic_stack_[(3) - (1)].tptr)->value.var; 
				  (yysemantic_stack_[(3) - (1)].tptr)->type = token::VAR;                       
				  SEAMS::math_error(aprepro, "Power");
				  undefined_warning(aprepro, (yysemantic_stack_[(3) - (1)].tptr)->name);          }
    break;

  case 66:

/* Line 678 of lalr1.cc  */
#line 217 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(3) - (1)].tptr)->value.fnctptr))();	}
    break;

  case 67:

/* Line 678 of lalr1.cc  */
#line 218 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_d))((yysemantic_stack_[(4) - (3)].val)); 	}
    break;

  case 68:

/* Line 678 of lalr1.cc  */
#line 219 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(4) - (1)].tptr)->value.fnctptr_c))((yysemantic_stack_[(4) - (3)].string)); 	}
    break;

  case 69:

/* Line 678 of lalr1.cc  */
#line 221 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_cd))((yysemantic_stack_[(6) - (3)].string), (yysemantic_stack_[(6) - (5)].val)); 	}
    break;

  case 70:

/* Line 678 of lalr1.cc  */
#line 223 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_cc))((yysemantic_stack_[(6) - (3)].string), (yysemantic_stack_[(6) - (5)].string)); 	}
    break;

  case 71:

/* Line 678 of lalr1.cc  */
#line 225 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(6) - (1)].tptr)->value.fnctptr_dd))((yysemantic_stack_[(6) - (3)].val), (yysemantic_stack_[(6) - (5)].val));	}
    break;

  case 72:

/* Line 678 of lalr1.cc  */
#line 227 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(8) - (1)].tptr)->value.fnctptr_ddd))((yysemantic_stack_[(8) - (3)].val), (yysemantic_stack_[(8) - (5)].val), (yysemantic_stack_[(8) - (7)].val)); }
    break;

  case 73:

/* Line 678 of lalr1.cc  */
#line 229 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(10) - (1)].tptr)->value.fnctptr_dddd))((yysemantic_stack_[(10) - (3)].val), (yysemantic_stack_[(10) - (5)].val), (yysemantic_stack_[(10) - (7)].val), (yysemantic_stack_[(10) - (9)].val)); }
    break;

  case 74:

/* Line 678 of lalr1.cc  */
#line 231 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(10) - (1)].tptr)->value.fnctptr_dddd))((yysemantic_stack_[(10) - (3)].val), (yysemantic_stack_[(10) - (5)].val), (yysemantic_stack_[(10) - (7)].val), (yysemantic_stack_[(10) - (9)].val)); }
    break;

  case 75:

/* Line 678 of lalr1.cc  */
#line 233 "aprepro.yy"
    { (yyval.val) = (*((yysemantic_stack_[(14) - (1)].tptr)->value.fnctptr_dddddd))((yysemantic_stack_[(14) - (3)].val), (yysemantic_stack_[(14) - (5)].val), (yysemantic_stack_[(14) - (7)].val), (yysemantic_stack_[(14) - (9)].val), (yysemantic_stack_[(14) - (11)].val), (yysemantic_stack_[(14) - (13)].val)); }
    break;

  case 76:

/* Line 678 of lalr1.cc  */
#line 234 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) + (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 77:

/* Line 678 of lalr1.cc  */
#line 235 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) - (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 78:

/* Line 678 of lalr1.cc  */
#line 236 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (1)].val) * (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 79:

/* Line 678 of lalr1.cc  */
#line 237 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (3)].val) == 0.)
				    {
				      yyerror(aprepro, "Zero divisor"); 
				      yyerrok;
				    }
				  else
				    (yyval.val) = (yysemantic_stack_[(3) - (1)].val) / (yysemantic_stack_[(3) - (3)].val); 			}
    break;

  case 80:

/* Line 678 of lalr1.cc  */
#line 244 "aprepro.yy"
    { if ((yysemantic_stack_[(3) - (3)].val) == 0.)
				    {
				      yyerror(aprepro, "Zero divisor");
				      yyerrok;
				    }
				  else
				    (yyval.val) = (int)(yysemantic_stack_[(3) - (1)].val) % (int)(yysemantic_stack_[(3) - (3)].val);		}
    break;

  case 81:

/* Line 678 of lalr1.cc  */
#line 251 "aprepro.yy"
    { (yyval.val) = -(yysemantic_stack_[(2) - (2)].val);				}
    break;

  case 82:

/* Line 678 of lalr1.cc  */
#line 252 "aprepro.yy"
    { (yyval.val) =  (yysemantic_stack_[(2) - (2)].val);				}
    break;

  case 83:

/* Line 678 of lalr1.cc  */
#line 253 "aprepro.yy"
    { errno = 0;
				  (yyval.val) = std::pow((yysemantic_stack_[(3) - (1)].val), (yysemantic_stack_[(3) - (3)].val)); 
				  SEAMS::math_error(aprepro, "Power");			}
    break;

  case 84:

/* Line 678 of lalr1.cc  */
#line 256 "aprepro.yy"
    { (yyval.val) = (yysemantic_stack_[(3) - (2)].val);				}
    break;

  case 85:

/* Line 678 of lalr1.cc  */
#line 257 "aprepro.yy"
    { errno = 0;
				  (yyval.val) = (double)((yysemantic_stack_[(3) - (2)].val) < 0 ? -floor(-((yysemantic_stack_[(3) - (2)].val))): floor((yysemantic_stack_[(3) - (2)].val)) );
				  SEAMS::math_error(aprepro, "floor (int)");		}
    break;

  case 86:

/* Line 678 of lalr1.cc  */
#line 260 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(1) - (1)].val)) ? 1 : 0; }
    break;

  case 87:

/* Line 678 of lalr1.cc  */
#line 261 "aprepro.yy"
    { (yyval.val) = ((yysemantic_stack_[(5) - (1)].val)) ? ((yysemantic_stack_[(5) - (3)].val)) : ((yysemantic_stack_[(5) - (5)].val));              }
    break;



/* Line 678 of lalr1.cc  */
#line 1092 "apr_parser.cc"
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
  const signed char Parser::yypact_ninf_ = -35;
  const short int
  Parser::yypact_[] =
  {
       -35,     1,   -35,   -14,   100,   -35,   -35,   -35,   -35,   -35,
     231,   315,   -13,     8,     9,   233,   233,   -35,   233,   233,
     233,    58,    87,   -18,   337,   592,   233,   233,   233,   233,
     233,   233,   -35,   -35,   233,   233,   233,   233,   233,   233,
     -35,   -35,   233,   175,   221,    15,    22,   622,   651,    -6,
      -6,    -6,   -35,   -35,   -35,   -35,   -35,   -35,   233,   233,
     233,   233,   -35,   233,   233,   233,   233,   233,   233,   233,
     233,   -35,   233,   233,   233,   233,   233,   233,   233,   233,
     233,   233,   233,   233,   233,    22,   763,   763,   763,   763,
     763,   763,    22,   763,   763,   763,   763,   763,   763,    22,
     763,   -35,   116,   398,   -35,   281,   428,   -35,   -35,   -35,
     793,   763,   -10,   -34,   251,   -26,   -26,   -26,   -26,   -26,
     -26,   -35,   763,   778,   502,    49,    49,    49,    49,    49,
     126,   126,    -6,    -6,    -6,    -6,   233,   -35,   233,   -35,
     233,   -35,   233,   -35,   233,   233,   125,   681,   368,   312,
     318,    22,   763,   -35,   -35,   233,   -35,   233,   233,   233,
     458,   518,   172,   287,   233,   -35,   233,   -35,   233,   -35,
     488,   711,   343,   233,   -35,   -35,   233,   548,   568,   233,
     -35,   741,   -35
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned char
  Parser::yydefact_[] =
  {
         2,     0,     1,     0,     0,     4,     3,     8,    40,    27,
      55,    43,    28,     0,     0,     0,     0,     7,     0,     0,
       0,     0,     0,    86,     0,     0,     0,     0,     0,     0,
       0,     0,    59,    58,     0,     0,     0,     0,     0,     0,
      47,    46,     0,     0,     0,    86,     0,     0,     0,    81,
      82,    11,    42,    57,    45,    41,    56,    44,     0,     0,
       0,     0,     6,     0,     0,     0,     0,     0,     0,     0,
       0,     5,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    29,    60,    62,    61,    64,
      63,    65,    31,    48,    51,    50,    53,    52,    54,    30,
      49,    66,     0,     0,    33,     0,     0,    20,    84,    85,
       0,     0,    18,    19,     0,    26,    25,    24,    23,    22,
      21,    35,    10,    16,    17,    15,    14,    13,    12,     9,
      77,    76,    80,    78,    79,    83,     0,    68,     0,    67,
       0,    32,     0,    34,     0,     0,     0,     0,     0,     0,
       0,    39,    87,    70,    69,     0,    71,     0,     0,     0,
       0,     0,     0,     0,     0,    72,     0,    38,     0,    37,
       0,     0,     0,     0,    74,    73,     0,     0,     0,     0,
      36,     0,    75
  };

  /* YYPGOTO[NTERM-NUM].  */
  const signed char
  Parser::yypgoto_[] =
  {
       -35,   -35,   -35,   -11,    55,    -4
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const signed char
  Parser::yydefgoto_[] =
  {
        -1,     1,     6,    23,    46,   111
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char Parser::yytable_ninf_ = -1;
  const unsigned char
  Parser::yytable_[] =
  {
        25,     2,     3,     7,    45,    70,    42,    67,    58,    59,
      60,    47,    48,    61,    49,    50,    51,     4,    60,    69,
      43,    44,    86,    87,    88,    89,    90,    91,   107,    61,
      93,    94,    95,    96,    97,    98,    84,    61,   100,   103,
     106,    58,    59,    60,     0,     0,     0,     5,   112,   113,
       0,    63,    64,    65,    66,    67,    68,   114,     0,    24,
      70,    52,    61,    53,    54,     0,   122,    69,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,    85,     0,     0,    79,    80,    81,    82,    83,    92,
      55,    84,    56,    57,     0,     0,     0,    99,   102,   105,
       0,     0,     0,     8,     9,    10,    11,    12,    13,    14,
       0,     0,    15,   110,    16,     0,     0,    17,   115,   116,
     117,   118,   119,   120,   121,     0,   136,     0,     0,   137,
       0,     0,   147,     0,   148,    18,    19,    70,   153,     0,
      20,   152,     0,    21,    22,    63,    64,    65,    66,    67,
      68,   160,     0,   161,    63,    64,    65,    66,    67,    68,
     170,    69,   171,    81,    82,    83,     0,     0,    84,   177,
      69,     0,     0,     0,     0,   181,     0,     0,     8,     9,
      10,    11,    12,    13,    14,   167,     0,    15,   101,    16,
       0,   146,     0,     0,     0,   149,     0,   150,     0,   151,
       0,    63,    64,    65,    66,    67,    68,     0,     0,     0,
      18,    19,     0,   162,   163,    20,     0,    69,    21,    22,
       0,     0,     0,   172,     8,     9,    10,    11,    12,    13,
      14,   178,     0,    15,   104,    16,     8,     9,    10,    11,
      12,    13,    14,     0,     0,    15,     0,    16,     0,     0,
      26,    27,    28,    29,    30,    31,    18,    19,     0,     0,
       0,    20,    70,     0,    21,    22,     0,     0,    18,    19,
       0,     0,     0,    20,    32,    33,    21,    22,    72,    73,
      74,    75,    76,    77,     0,    78,    79,    80,    81,    82,
      83,   140,     0,    84,   141,     0,     0,   168,     0,   145,
     169,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      63,    64,    65,    66,    67,    68,    63,    64,    65,    66,
      67,    68,   158,     0,     0,     0,    69,     0,   159,     0,
       0,     0,    69,     0,    34,    35,    36,    37,    38,    39,
       0,    63,    64,    65,    66,    67,    68,    63,    64,    65,
      66,    67,    68,   176,    62,     0,     0,    69,    40,    41,
       0,     0,     0,    69,     0,     0,    63,    64,    65,    66,
      67,    68,    63,    64,    65,    66,    67,    68,   155,    70,
       0,   156,    69,     0,     0,     0,   157,     0,    69,     0,
       0,     0,     0,     0,     0,    72,    73,    74,    75,    76,
      77,     0,    78,    79,    80,    81,    82,    83,   138,    70,
      84,   139,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    72,    73,    74,    75,    76,
      77,     0,    78,    79,    80,    81,    82,    83,   142,    70,
      84,   143,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    72,    73,    74,    75,    76,
      77,     0,    78,    79,    80,    81,    82,    83,   164,    70,
      84,   165,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    72,    73,    74,    75,    76,
      77,     0,    78,    79,    80,    81,    82,    83,   173,    70,
      84,   174,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    70,     0,    72,    73,    74,    75,    76,
      77,     0,    78,    79,    80,    81,    82,    83,   166,    70,
      84,    74,    75,    76,    77,     0,    78,    79,    80,    81,
      82,    83,     0,     0,    84,    72,    73,    74,    75,    76,
      77,     0,    78,    79,    80,    81,    82,    83,   179,    70,
      84,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    72,    73,    74,    75,    76,
      77,   180,    78,    79,    80,    81,    82,    83,     0,     0,
      84,     0,     0,     0,     0,     0,     0,    63,    64,    65,
      66,    67,    68,    70,     0,     0,     0,     0,     0,    71,
       0,     0,     0,    69,     0,     0,     0,     0,     0,    72,
      73,    74,    75,    76,    77,     0,    78,    79,    80,    81,
      82,    83,     0,    70,    84,   108,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    72,
      73,    74,    75,    76,    77,     0,    78,    79,    80,    81,
      82,    83,    70,     0,    84,     0,   109,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    72,    73,
      74,    75,    76,    77,     0,    78,    79,    80,    81,    82,
      83,     0,    70,    84,   154,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    72,    73,
      74,    75,    76,    77,     0,    78,    79,    80,    81,    82,
      83,     0,    70,    84,   175,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    72,    73,
      74,    75,    76,    77,     0,    78,    79,    80,    81,    82,
      83,     0,    70,    84,   182,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    72,    73,
      74,    75,    76,    77,    70,    78,    79,    80,    81,    82,
      83,     0,     0,    84,     0,     0,     0,     0,     0,    70,
      72,    73,    74,    75,    76,    77,     0,    78,    79,    80,
      81,    82,    83,     0,     0,    84,    73,    74,    75,    76,
      77,     0,    78,    79,    80,    81,    82,    83,   144,     0,
      84,     0,    63,    64,    65,    66,    67,    68,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    69
  };

  /* YYCHECK.  */
  const short int
  Parser::yycheck_[] =
  {
         4,     0,     1,    17,    15,    11,    19,    33,    26,    27,
      28,    15,    16,    47,    18,    19,    20,    16,    28,    45,
      12,    12,    26,    27,    28,    29,    30,    31,    13,    47,
      34,    35,    36,    37,    38,    39,    42,    47,    42,    43,
      44,    26,    27,    28,    -1,    -1,    -1,    46,    59,    60,
      -1,    29,    30,    31,    32,    33,    34,    61,    -1,     4,
      11,     3,    47,     5,     6,    -1,    70,    45,    72,    73,
      74,    75,    76,    77,    78,    79,    80,    81,    82,    83,
      84,    26,    -1,    -1,    35,    36,    37,    38,    39,    34,
       3,    42,     5,     6,    -1,    -1,    -1,    42,    43,    44,
      -1,    -1,    -1,     3,     4,     5,     6,     7,     8,     9,
      -1,    -1,    12,    58,    14,    -1,    -1,    17,    63,    64,
      65,    66,    67,    68,    69,    -1,    10,    -1,    -1,    13,
      -1,    -1,   136,    -1,   138,    35,    36,    11,    13,    -1,
      40,   145,    -1,    43,    44,    29,    30,    31,    32,    33,
      34,   155,    -1,   157,    29,    30,    31,    32,    33,    34,
     164,    45,   166,    37,    38,    39,    -1,    -1,    42,   173,
      45,    -1,    -1,    -1,    -1,   179,    -1,    -1,     3,     4,
       5,     6,     7,     8,     9,    13,    -1,    12,    13,    14,
      -1,   136,    -1,    -1,    -1,   140,    -1,   142,    -1,   144,
      -1,    29,    30,    31,    32,    33,    34,    -1,    -1,    -1,
      35,    36,    -1,   158,   159,    40,    -1,    45,    43,    44,
      -1,    -1,    -1,   168,     3,     4,     5,     6,     7,     8,
       9,   176,    -1,    12,    13,    14,     3,     4,     5,     6,
       7,     8,     9,    -1,    -1,    12,    -1,    14,    -1,    -1,
      19,    20,    21,    22,    23,    24,    35,    36,    -1,    -1,
      -1,    40,    11,    -1,    43,    44,    -1,    -1,    35,    36,
      -1,    -1,    -1,    40,    43,    44,    43,    44,    27,    28,
      29,    30,    31,    32,    -1,    34,    35,    36,    37,    38,
      39,    10,    -1,    42,    13,    -1,    -1,    10,    -1,    48,
      13,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      29,    30,    31,    32,    33,    34,    29,    30,    31,    32,
      33,    34,    10,    -1,    -1,    -1,    45,    -1,    10,    -1,
      -1,    -1,    45,    -1,    19,    20,    21,    22,    23,    24,
      -1,    29,    30,    31,    32,    33,    34,    29,    30,    31,
      32,    33,    34,    10,    17,    -1,    -1,    45,    43,    44,
      -1,    -1,    -1,    45,    -1,    -1,    29,    30,    31,    32,
      33,    34,    29,    30,    31,    32,    33,    34,    10,    11,
      -1,    13,    45,    -1,    -1,    -1,    18,    -1,    45,    -1,
      -1,    -1,    -1,    -1,    -1,    27,    28,    29,    30,    31,
      32,    -1,    34,    35,    36,    37,    38,    39,    10,    11,
      42,    13,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    27,    28,    29,    30,    31,
      32,    -1,    34,    35,    36,    37,    38,    39,    10,    11,
      42,    13,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    27,    28,    29,    30,    31,
      32,    -1,    34,    35,    36,    37,    38,    39,    10,    11,
      42,    13,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    27,    28,    29,    30,    31,
      32,    -1,    34,    35,    36,    37,    38,    39,    10,    11,
      42,    13,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    11,    -1,    27,    28,    29,    30,    31,
      32,    -1,    34,    35,    36,    37,    38,    39,    10,    11,
      42,    29,    30,    31,    32,    -1,    34,    35,    36,    37,
      38,    39,    -1,    -1,    42,    27,    28,    29,    30,    31,
      32,    -1,    34,    35,    36,    37,    38,    39,    10,    11,
      42,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    27,    28,    29,    30,    31,
      32,    13,    34,    35,    36,    37,    38,    39,    -1,    -1,
      42,    -1,    -1,    -1,    -1,    -1,    -1,    29,    30,    31,
      32,    33,    34,    11,    -1,    -1,    -1,    -1,    -1,    17,
      -1,    -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,    27,
      28,    29,    30,    31,    32,    -1,    34,    35,    36,    37,
      38,    39,    -1,    11,    42,    13,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    27,
      28,    29,    30,    31,    32,    -1,    34,    35,    36,    37,
      38,    39,    11,    -1,    42,    -1,    15,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    27,    28,
      29,    30,    31,    32,    -1,    34,    35,    36,    37,    38,
      39,    -1,    11,    42,    13,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    27,    28,
      29,    30,    31,    32,    -1,    34,    35,    36,    37,    38,
      39,    -1,    11,    42,    13,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    27,    28,
      29,    30,    31,    32,    -1,    34,    35,    36,    37,    38,
      39,    -1,    11,    42,    13,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    27,    28,
      29,    30,    31,    32,    11,    34,    35,    36,    37,    38,
      39,    -1,    -1,    42,    -1,    -1,    -1,    -1,    -1,    11,
      27,    28,    29,    30,    31,    32,    -1,    34,    35,    36,
      37,    38,    39,    -1,    -1,    42,    28,    29,    30,    31,
      32,    -1,    34,    35,    36,    37,    38,    39,    25,    -1,
      42,    -1,    29,    30,    31,    32,    33,    34,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    45
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned char
  Parser::yystos_[] =
  {
         0,    50,     0,     1,    16,    46,    51,    17,     3,     4,
       5,     6,     7,     8,     9,    12,    14,    17,    35,    36,
      40,    43,    44,    52,    53,    54,    19,    20,    21,    22,
      23,    24,    43,    44,    19,    20,    21,    22,    23,    24,
      43,    44,    19,    12,    12,    52,    53,    54,    54,    54,
      54,    54,     3,     5,     6,     3,     5,     6,    26,    27,
      28,    47,    17,    29,    30,    31,    32,    33,    34,    45,
      11,    17,    27,    28,    29,    30,    31,    32,    34,    35,
      36,    37,    38,    39,    42,    53,    54,    54,    54,    54,
      54,    54,    53,    54,    54,    54,    54,    54,    54,    53,
      54,    13,    53,    54,    13,    53,    54,    13,    13,    15,
      53,    54,    52,    52,    54,    53,    53,    53,    53,    53,
      53,    53,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    10,    13,    10,    13,
      10,    13,    10,    13,    25,    48,    53,    54,    54,    53,
      53,    53,    54,    13,    13,    10,    13,    18,    10,    10,
      54,    54,    53,    53,    10,    13,    10,    13,    10,    13,
      54,    54,    53,    10,    13,    13,    10,    54,    53,    10,
      13,    54,    13
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
     295,   296,   297,   298,   299,   300,    10,    63,    58
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned char
  Parser::yyr1_[] =
  {
         0,    49,    50,    50,    51,    51,    51,    51,    51,    52,
      52,    52,    52,    52,    52,    52,    52,    52,    52,    52,
      52,    52,    52,    52,    52,    52,    52,    53,    53,    53,
      53,    53,    53,    53,    53,    53,    53,    53,    53,    53,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  Parser::yyr2_[] =
  {
         0,     2,     0,     2,     1,     3,     3,     2,     2,     3,
       3,     2,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     1,     3,
       3,     3,     4,     3,     4,     3,    12,     8,     8,     5,
       1,     2,     2,     1,     2,     2,     2,     2,     3,     3,
       3,     3,     3,     3,     3,     1,     2,     2,     2,     2,
       3,     3,     3,     3,     3,     3,     3,     4,     4,     6,
       6,     6,     8,    10,    10,    14,     3,     3,     3,     3,
       3,     2,     2,     3,     3,     3,     1,     5
  };

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
  /* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
     First, the terminals, then, starting at \a yyntokens_, nonterminals.  */
  const char*
  const Parser::yytname_[] =
  {
    "\"end of file\"", "error", "$undefined", "NUM", "QSTRING", "UNDVAR",
  "VAR", "SVAR", "FNCT", "SFNCT", "COMMA", "RT", "LPAR", "RPAR", "LBRACK",
  "RBRACK", "LBRACE", "RBRACE", "SEMI", "EQUAL", "EQ_MINUS", "EQ_PLUS",
  "EQ_DIV", "EQ_TIME", "EQ_POW", "COLON", "QUEST", "LOR", "LAND", "NE",
  "EQ", "GE", "LE", "GT", "LT", "SUB", "PLU", "MOD", "TIM", "DIV", "NOT",
  "UNARY", "POW", "DEC", "INC", "CONCAT", "'\\n'", "'?'", "':'", "$accept",
  "input", "line", "bool", "sexp", "exp", 0
  };
#endif

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const Parser::rhs_number_type
  Parser::yyrhs_[] =
  {
        50,     0,    -1,    -1,    50,    51,    -1,    46,    -1,    16,
      54,    17,    -1,    16,    53,    17,    -1,    16,    17,    -1,
       1,    17,    -1,    54,    34,    54,    -1,    54,    11,    54,
      -1,    40,    54,    -1,    54,    32,    54,    -1,    54,    31,
      54,    -1,    54,    30,    54,    -1,    54,    29,    54,    -1,
      54,    27,    54,    -1,    54,    28,    54,    -1,    52,    27,
      52,    -1,    52,    28,    52,    -1,    12,    52,    13,    -1,
      53,    34,    53,    -1,    53,    33,    53,    -1,    53,    32,
      53,    -1,    53,    31,    53,    -1,    53,    30,    53,    -1,
      53,    29,    53,    -1,     4,    -1,     7,    -1,     5,    19,
      53,    -1,     7,    19,    53,    -1,     6,    19,    53,    -1,
       9,    12,    53,    13,    -1,     9,    12,    13,    -1,     9,
      12,    54,    13,    -1,    53,    45,    53,    -1,     9,    12,
      54,    10,    53,    10,    53,    10,    53,    10,    53,    13,
      -1,     9,    12,    54,    10,    53,    10,    53,    13,    -1,
       9,    12,    53,    10,    53,    10,    53,    13,    -1,    52,
      26,    53,    25,    53,    -1,     3,    -1,    44,     3,    -1,
      43,     3,    -1,     6,    -1,    44,     6,    -1,    43,     6,
      -1,     6,    44,    -1,     6,    43,    -1,     6,    19,    54,
      -1,     7,    19,    54,    -1,     6,    21,    54,    -1,     6,
      20,    54,    -1,     6,    23,    54,    -1,     6,    22,    54,
      -1,     6,    24,    54,    -1,     5,    -1,    44,     5,    -1,
      43,     5,    -1,     5,    44,    -1,     5,    43,    -1,     5,
      19,    54,    -1,     5,    21,    54,    -1,     5,    20,    54,
      -1,     5,    23,    54,    -1,     5,    22,    54,    -1,     5,
      24,    54,    -1,     8,    12,    13,    -1,     8,    12,    54,
      13,    -1,     8,    12,    53,    13,    -1,     8,    12,    53,
      10,    54,    13,    -1,     8,    12,    53,    10,    53,    13,
      -1,     8,    12,    54,    10,    54,    13,    -1,     8,    12,
      54,    10,    54,    10,    54,    13,    -1,     8,    12,    54,
      10,    54,    18,    54,    10,    54,    13,    -1,     8,    12,
      54,    10,    54,    10,    54,    10,    54,    13,    -1,     8,
      12,    54,    10,    54,    10,    54,    10,    54,    10,    54,
      10,    54,    13,    -1,    54,    36,    54,    -1,    54,    35,
      54,    -1,    54,    38,    54,    -1,    54,    39,    54,    -1,
      54,    37,    54,    -1,    35,    54,    -1,    36,    54,    -1,
      54,    42,    54,    -1,    12,    54,    13,    -1,    14,    54,
      15,    -1,    52,    -1,    52,    47,    54,    48,    54,    -1
  };

  /* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
     YYRHS.  */
  const unsigned short int
  Parser::yyprhs_[] =
  {
         0,     0,     3,     4,     7,     9,    13,    17,    20,    23,
      27,    31,    34,    38,    42,    46,    50,    54,    58,    62,
      66,    70,    74,    78,    82,    86,    90,    94,    96,    98,
     102,   106,   110,   115,   119,   124,   128,   141,   150,   159,
     165,   167,   170,   173,   175,   178,   181,   184,   187,   191,
     195,   199,   203,   207,   211,   215,   217,   220,   223,   226,
     229,   233,   237,   241,   245,   249,   253,   257,   262,   267,
     274,   281,   288,   297,   308,   319,   334,   338,   342,   346,
     350,   354,   357,   360,   364,   368,   372,   374
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  Parser::yyrline_[] =
  {
         0,    93,    93,    94,    97,    98,   105,   109,   110,   113,
     114,   115,   116,   117,   118,   119,   120,   121,   122,   123,
     124,   127,   128,   129,   130,   131,   132,   134,   135,   136,
     138,   141,   145,   146,   147,   148,   153,   155,   157,   159,
     161,   162,   163,   164,   165,   166,   167,   168,   169,   171,
     174,   175,   176,   177,   178,   183,   185,   188,   191,   194,
     197,   199,   202,   205,   208,   211,   217,   218,   219,   220,
     222,   224,   226,   228,   230,   232,   234,   235,   236,   237,
     244,   251,   252,   253,   256,   257,   260,   261
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
      46,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    48,     2,
       2,     2,     2,    47,     2,     2,     2,     2,     2,     2,
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
      45
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int Parser::yyeof_ = 0;
  const int Parser::yylast_ = 838;
  const int Parser::yynnts_ = 6;
  const int Parser::yyempty_ = -2;
  const int Parser::yyfinal_ = 2;
  const int Parser::yyterror_ = 1;
  const int Parser::yyerrcode_ = 256;
  const int Parser::yyntokens_ = 49;

  const unsigned int Parser::yyuser_token_number_max_ = 300;
  const Parser::token_number_type Parser::yyundef_token_ = 2;


} // SEAMS

/* Line 1054 of lalr1.cc  */
#line 1799 "apr_parser.cc"


/* Line 1056 of lalr1.cc  */
#line 265 "aprepro.yy"


void SEAMS::Parser::error(const Parser::location_type&, const std::string& m)
{
    aprepro.error(m);
}


