// A Bison parser, made by GNU Bison 3.0.4.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.

// Take the name prefix into account.
#define yylex   SEAMSlex

// First part of user declarations.
#line 33 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:404

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

 

#line 57 "apr_parser.cc" // lalr1.cc:404

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

#include "aprepro_parser.h"

// User implementation prologue.
#line 118 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:412


#include "aprepro.h"
#include "apr_scanner.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex


#line 83 "apr_parser.cc" // lalr1.cc:412


#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif



// Suppress unused-variable warnings by "using" E.
#define YYUSE(E) ((void) (E))

// Enable debugging if requested.
#if YYDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << '\n';                  \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yystack_print_ ();                \
  } while (false)

#else // !YYDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YYUSE(Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void>(0)
# define YY_STACK_PRINT()                static_cast<void>(0)

#endif // !YYDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)


namespace SEAMS {
#line 150 "apr_parser.cc" // lalr1.cc:479

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
              // Fall through.
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
  {}

  Parser::~Parser ()
  {}


  /*---------------.
  | Symbol types.  |
  `---------------*/

  inline
  Parser::syntax_error::syntax_error (const std::string& m)
    : std::runtime_error (m)
  {}

  // basic_symbol.
  template <typename Base>
  inline
  Parser::basic_symbol<Base>::basic_symbol ()
    : value ()
  {}

  template <typename Base>
  inline
  Parser::basic_symbol<Base>::basic_symbol (const basic_symbol& other)
    : Base (other)
    , value ()
  {
    value = other.value;
  }


  template <typename Base>
  inline
  Parser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, semantic_type  v)
    : Base (t)
    , value (std::move(v))
  {}


  /// Constructor for valueless symbols.
  template <typename Base>
  inline
  Parser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t)
    : Base (t)
    , value ()
  {}

  template <typename Base>
  inline
  Parser::basic_symbol<Base>::~basic_symbol ()
  {
    clear ();
  }

  template <typename Base>
  inline
  void
  Parser::basic_symbol<Base>::clear ()
  {
    Base::clear ();
  }

  template <typename Base>
  inline
  bool
  Parser::basic_symbol<Base>::empty () const
  {
    return Base::type_get () == empty_symbol;
  }

  template <typename Base>
  inline
  void
  Parser::basic_symbol<Base>::move (basic_symbol& s)
  {
    super_type::move(s);
    value = s.value;
  }

  // by_type.
  inline
  Parser::by_type::by_type ()
    : type (empty_symbol)
  {}

  inline
  Parser::by_type::by_type (const by_type& other)
    : type (other.type)
  {}

  inline
  Parser::by_type::by_type (token_type t)
    : type (yytranslate_ (t))
  {}

  inline
  void
  Parser::by_type::clear ()
  {
    type = empty_symbol;
  }

  inline
  void
  Parser::by_type::move (by_type& that)
  {
    type = that.type;
    that.clear ();
  }

  inline
  int
  Parser::by_type::type_get () const
  {
    return type;
  }


  // by_state.
  inline
  Parser::by_state::by_state ()
    : state (empty_state)
  {}

  inline
  Parser::by_state::by_state (const by_state& other)
    : state (other.state)
  {}

  inline
  void
  Parser::by_state::clear ()
  {
    state = empty_state;
  }

  inline
  void
  Parser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  inline
  Parser::by_state::by_state (state_type s)
    : state (s)
  {}

  inline
  Parser::symbol_number_type
  Parser::by_state::type_get () const
  {
    if (state == empty_state)
      return empty_symbol;
    else
      return yystos_[state];
  }

  inline
  Parser::stack_symbol_type::stack_symbol_type ()
  {}


  inline
  Parser::stack_symbol_type::stack_symbol_type (state_type s, symbol_type& that)
    : super_type (s)
  {
    value = that.value;
    // that is emptied.
    that.type = empty_symbol;
  }

  inline
  Parser::stack_symbol_type&
  Parser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
    value = that.value;
    return *this;
  }


  template <typename Base>
  inline
  void
  Parser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);

    // User destructor.
    YYUSE (yysym.type_get ());
  }

#if YYDEBUG
  template <typename Base>
  void
  Parser::yy_print_ (std::ostream& yyo,
                                     const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YYUSE (yyoutput);
    symbol_number_type yytype = yysym.type_get ();
    // Avoid a (spurious) G++ 4.8 warning about "array subscript is
    // below array bounds".
    if (yysym.empty ())
      std::abort ();
    yyo << (yytype < yyntokens_ ? "token" : "nterm")
        << ' ' << yytname_[yytype] << " (";
    YYUSE (yytype);
    yyo << ')';
  }
#endif

  inline
  void
  Parser::yypush_ (const char* m, state_type s, symbol_type& sym)
  {
    stack_symbol_type t (s, sym);
    yypush_ (m, t);
  }

  inline
  void
  Parser::yypush_ (const char* m, stack_symbol_type& s)
  {
    if (m)
      YY_SYMBOL_PRINT (m, s);
    yystack_.push (s);
  }

  inline
  void
  Parser::yypop_ (unsigned int n)
  {
    yystack_.pop (n);
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
#endif // YYDEBUG

  inline Parser::state_type
  Parser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - yyntokens_] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - yyntokens_];
  }

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
    // State.
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The return value of parse ().
    int yyresult;

    // FIXME: This shoud be completely indented.  It is not yet to
    // avoid gratuitous conflicts when merging into the master branch.
    try
      {
    YYCDEBUG << "Starting parse" << '\n';


    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, yyla);

    // A new symbol was pushed on the stack.
  yynewstate:
    YYCDEBUG << "Entering state " << yystack_[0].state << '\n';

    // Accept?
    if (yystack_[0].state == yyfinal_)
      goto yyacceptlab;

    goto yybackup;

    // Backup.
  yybackup:

    // Try to take a decision without lookahead.
    yyn = yypact_[yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token: ";
        try
          {
            yyla.type = yytranslate_ (yylex (&yyla.value));
          }
        catch (const syntax_error& yyexc)
          {
            error (yyexc);
            goto yyerrlab1;
          }
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.type_get ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.type_get ())
      goto yydefault;

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", yyn, yyla);
    goto yynewstate;

  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;

  /*-----------------------------.
  | yyreduce -- Do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_(yystack_[yylen].state, yyr1_[yyn]);
      /* If YYLEN is nonzero, implement the default value of the
         action: '$$ = $1'.  Otherwise, use the top of the stack.

         Otherwise, the following line sets YYLHS.VALUE to garbage.
         This behavior is undocumented and Bison users should not rely
         upon it.  */
      if (yylen)
        yylhs.value = yystack_[yylen - 1].value;
      else
        yylhs.value = yystack_[0].value;


      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
      try
        {
          switch (yyn)
            {
  case 4:
#line 137 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { if (echo) aprepro.lexer->LexerOutput("\n", 1); }
#line 615 "apr_parser.cc" // lalr1.cc:859
    break;

  case 5:
#line 138 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { if (echo) {
	                             static char tmpstr[512];
				     SEAMS::symrec *format = aprepro.getsym("_FORMAT");
				     int len = sprintf(tmpstr, format->value.svar, (yystack_[1].value.val));
				     aprepro.lexer->LexerOutput(tmpstr, len);
				   }
                                }
#line 627 "apr_parser.cc" // lalr1.cc:859
    break;

  case 6:
#line 145 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { if (echo && (yystack_[1].value.string) != nullptr) {
				    aprepro.lexer->LexerOutput((yystack_[1].value.string), strlen((yystack_[1].value.string)));
                                  }
                                }
#line 636 "apr_parser.cc" // lalr1.cc:859
    break;

  case 7:
#line 149 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {                                       }
#line 642 "apr_parser.cc" // lalr1.cc:859
    break;

  case 8:
#line 150 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {                                       }
#line 648 "apr_parser.cc" // lalr1.cc:859
    break;

  case 9:
#line 151 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { yyerrok;				}
#line 654 "apr_parser.cc" // lalr1.cc:859
    break;

  case 10:
#line 154 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) < (yystack_[0].value.val);                         }
#line 660 "apr_parser.cc" // lalr1.cc:859
    break;

  case 11:
#line 155 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) > (yystack_[0].value.val);                         }
#line 666 "apr_parser.cc" // lalr1.cc:859
    break;

  case 12:
#line 156 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = !((yystack_[0].value.val));                           }
#line 672 "apr_parser.cc" // lalr1.cc:859
    break;

  case 13:
#line 157 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) <= (yystack_[0].value.val);                        }
#line 678 "apr_parser.cc" // lalr1.cc:859
    break;

  case 14:
#line 158 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) >= (yystack_[0].value.val);                        }
#line 684 "apr_parser.cc" // lalr1.cc:859
    break;

  case 15:
#line 159 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) == (yystack_[0].value.val);                        }
#line 690 "apr_parser.cc" // lalr1.cc:859
    break;

  case 16:
#line 160 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) != (yystack_[0].value.val);                        }
#line 696 "apr_parser.cc" // lalr1.cc:859
    break;

  case 17:
#line 161 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);                        }
#line 702 "apr_parser.cc" // lalr1.cc:859
    break;

  case 18:
#line 162 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);                        }
#line 708 "apr_parser.cc" // lalr1.cc:859
    break;

  case 19:
#line 163 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);                        }
#line 714 "apr_parser.cc" // lalr1.cc:859
    break;

  case 20:
#line 164 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);                        }
#line 720 "apr_parser.cc" // lalr1.cc:859
    break;

  case 21:
#line 165 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[1].value.val);                              }
#line 726 "apr_parser.cc" // lalr1.cc:859
    break;

  case 22:
#line 168 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (strcmp((yystack_[2].value.string),(yystack_[0].value.string)) <  0 ? 1 : 0);	}
#line 732 "apr_parser.cc" // lalr1.cc:859
    break;

  case 23:
#line 169 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (strcmp((yystack_[2].value.string),(yystack_[0].value.string)) >  0 ? 1 : 0);	}
#line 738 "apr_parser.cc" // lalr1.cc:859
    break;

  case 24:
#line 170 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (strcmp((yystack_[2].value.string),(yystack_[0].value.string)) <= 0 ? 1 : 0);	}
#line 744 "apr_parser.cc" // lalr1.cc:859
    break;

  case 25:
#line 171 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (strcmp((yystack_[2].value.string),(yystack_[0].value.string)) >= 0 ? 1 : 0);	}
#line 750 "apr_parser.cc" // lalr1.cc:859
    break;

  case 26:
#line 172 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (strcmp((yystack_[2].value.string),(yystack_[0].value.string)) == 0 ? 1 : 0);	}
#line 756 "apr_parser.cc" // lalr1.cc:859
    break;

  case 27:
#line 173 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (strcmp((yystack_[2].value.string),(yystack_[0].value.string)) != 0 ? 1 : 0);	}
#line 762 "apr_parser.cc" // lalr1.cc:859
    break;

  case 28:
#line 175 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.arrval) = (yystack_[0].value.tptr)->value.avar;}
#line 768 "apr_parser.cc" // lalr1.cc:859
    break;

  case 29:
#line 176 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[3].value.tptr), (yystack_[3].value.tptr)->value.arrfnct_c == nullptr))
	    (yylhs.value.arrval) = (*((yystack_[3].value.tptr)->value.arrfnct_c))((yystack_[1].value.string));
	  else
	    yyerrok;
	}
#line 779 "apr_parser.cc" // lalr1.cc:859
    break;

  case 30:
#line 182 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[5].value.tptr), (yystack_[5].value.tptr)->value.arrfnct_cd == nullptr))
	    (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cd))((yystack_[3].value.string),(yystack_[1].value.val));
	  else
	    yyerrok;
	}
#line 790 "apr_parser.cc" // lalr1.cc:859
    break;

  case 31:
#line 188 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[5].value.tptr), (yystack_[5].value.tptr)->value.arrfnct_cc == nullptr))
	    (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cc))((yystack_[3].value.string),(yystack_[1].value.string));
	  else
	    yyerrok;
	}
#line 801 "apr_parser.cc" // lalr1.cc:859
    break;

  case 32:
#line 194 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[5].value.tptr), (yystack_[5].value.tptr)->value.arrfnct_dd == nullptr))
	    (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_dd))((yystack_[3].value.val),(yystack_[1].value.val));
	  else
	    yyerrok;
	}
#line 812 "apr_parser.cc" // lalr1.cc:859
    break;

  case 33:
#line 200 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[3].value.tptr), (yystack_[3].value.tptr)->value.arrfnct_d == nullptr))
	    (yylhs.value.arrval) = (*((yystack_[3].value.tptr)->value.arrfnct_d))((yystack_[1].value.val));
	  else
	    yyerrok;
	}
#line 823 "apr_parser.cc" // lalr1.cc:859
    break;

  case 34:
#line 206 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[3].value.tptr), (yystack_[3].value.tptr)->value.arrfnct_a == nullptr))
	    (yylhs.value.arrval) = (*((yystack_[3].value.tptr)->value.arrfnct_a))((yystack_[1].value.arrval));
	  else
	    yyerrok;
	}
#line 834 "apr_parser.cc" // lalr1.cc:859
    break;

  case 35:
#line 212 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.arrval) = (yystack_[0].value.arrval); delete (yystack_[2].value.tptr)->value.avar; (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval); 
                                  redefined_warning(aprepro, (yystack_[2].value.tptr));
                                  set_type(aprepro, (yystack_[2].value.tptr), token::AVAR); }
#line 842 "apr_parser.cc" // lalr1.cc:859
    break;

  case 36:
#line 215 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.arrval) = (yystack_[0].value.arrval); (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval); 
                                  set_type(aprepro, (yystack_[2].value.tptr), token::AVAR); }
#line 849 "apr_parser.cc" // lalr1.cc:859
    break;

  case 37:
#line 217 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { if ((yystack_[2].value.arrval)->cols == (yystack_[0].value.arrval)->cols && (yystack_[2].value.arrval)->rows == (yystack_[0].value.arrval)->rows ) {
                                     (yylhs.value.arrval) = array_add((yystack_[2].value.arrval), (yystack_[0].value.arrval)); 
                                  }
                                  else {
                                    yyerror(aprepro, "Arrays do not have same row and column count"); 
                                    yyerrok;
                                  }
                                }
#line 862 "apr_parser.cc" // lalr1.cc:859
    break;

  case 38:
#line 225 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), -1.0);           }
#line 868 "apr_parser.cc" // lalr1.cc:859
    break;

  case 39:
#line 227 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { if ((yystack_[2].value.arrval)->cols == (yystack_[0].value.arrval)->cols && (yystack_[2].value.arrval)->rows == (yystack_[0].value.arrval)->rows ) {
                                     (yylhs.value.arrval) = array_sub((yystack_[2].value.arrval), (yystack_[0].value.arrval)); 
                                  }
                                  else {
                                    yyerror(aprepro, "Arrays do not have same row and column count"); 
                                    yyerrok;
                                  }
                                }
#line 881 "apr_parser.cc" // lalr1.cc:859
    break;

  case 40:
#line 235 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.arrval) = array_scale((yystack_[2].value.arrval), (yystack_[0].value.val));             }
#line 887 "apr_parser.cc" // lalr1.cc:859
    break;

  case 41:
#line 236 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.arrval) = array_scale((yystack_[2].value.arrval), 1.0/(yystack_[0].value.val));         }
#line 893 "apr_parser.cc" // lalr1.cc:859
    break;

  case 42:
#line 237 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), (yystack_[2].value.val));             }
#line 899 "apr_parser.cc" // lalr1.cc:859
    break;

  case 43:
#line 238 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { if ((yystack_[2].value.arrval)->cols == (yystack_[0].value.arrval)->rows) {
                                    (yylhs.value.arrval) = array_mult((yystack_[2].value.arrval), (yystack_[0].value.arrval));
                                  }
                                  else {
                                    yyerror(aprepro, "Column count of first array does not match row count of second array"); 
                                    yyerrok;
                                  }
				}
#line 912 "apr_parser.cc" // lalr1.cc:859
    break;

  case 44:
#line 247 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.string) = (yystack_[0].value.string);				}
#line 918 "apr_parser.cc" // lalr1.cc:859
    break;

  case 45:
#line 248 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.string) = (char*)(yystack_[0].value.tptr)->value.svar;			}
#line 924 "apr_parser.cc" // lalr1.cc:859
    break;

  case 46:
#line 249 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.string) = (char*)(yystack_[0].value.tptr)->value.svar;			}
#line 930 "apr_parser.cc" // lalr1.cc:859
    break;

  case 47:
#line 250 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.string) = (yystack_[0].value.string); (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
		                  set_type(aprepro, (yystack_[2].value.tptr), Parser::token::SVAR);	}
#line 937 "apr_parser.cc" // lalr1.cc:859
    break;

  case 48:
#line 252 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.string) = (yystack_[0].value.string); 
				  (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
				  redefined_warning(aprepro, (yystack_[2].value.tptr));          }
#line 945 "apr_parser.cc" // lalr1.cc:859
    break;

  case 49:
#line 255 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.string) = (yystack_[0].value.string); 
				  (yystack_[2].value.tptr)->value.svar= (yystack_[0].value.string);
				  redefined_warning(aprepro, (yystack_[2].value.tptr));          
		                  set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);		}
#line 954 "apr_parser.cc" // lalr1.cc:859
    break;

  case 50:
#line 259 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[2].value.tptr)); YYERROR; }
#line 960 "apr_parser.cc" // lalr1.cc:859
    break;

  case 51:
#line 260 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[2].value.tptr)); YYERROR; }
#line 966 "apr_parser.cc" // lalr1.cc:859
    break;

  case 52:
#line 261 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[3].value.tptr), (yystack_[3].value.tptr)->value.strfnct_c == nullptr))
	    (yylhs.value.string) = (char*)(*((yystack_[3].value.tptr)->value.strfnct_c))((yystack_[1].value.string));
	  else
	    (yylhs.value.string) = (char*)"";
	}
#line 977 "apr_parser.cc" // lalr1.cc:859
    break;

  case 53:
#line 267 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[2].value.tptr), (yystack_[2].value.tptr)->value.strfnct == nullptr))
	    (yylhs.value.string) = (char*)(*((yystack_[2].value.tptr)->value.strfnct))();
	  else
	    (yylhs.value.string) = (char*)"";
	}
#line 988 "apr_parser.cc" // lalr1.cc:859
    break;

  case 54:
#line 273 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[3].value.tptr), (yystack_[3].value.tptr)->value.strfnct_d == nullptr))
	    (yylhs.value.string) = (char*)(*((yystack_[3].value.tptr)->value.strfnct_d))((yystack_[1].value.val));
	  else
	    (yylhs.value.string) = (char*)"";
	}
#line 999 "apr_parser.cc" // lalr1.cc:859
    break;

  case 55:
#line 279 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[3].value.tptr), (yystack_[3].value.tptr)->value.strfnct_a == nullptr))
	    (yylhs.value.string) = (char*)(*((yystack_[3].value.tptr)->value.strfnct_a))((yystack_[1].value.arrval));
	  else
	    (yylhs.value.string) = (char*)"";
	}
#line 1010 "apr_parser.cc" // lalr1.cc:859
    break;

  case 56:
#line 285 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { concat_string((yystack_[2].value.string), (yystack_[0].value.string), &(yylhs.value.string)); }
#line 1016 "apr_parser.cc" // lalr1.cc:859
    break;

  case 57:
#line 286 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[5].value.tptr), (yystack_[5].value.tptr)->value.strfnct_dd == nullptr))
	    (yylhs.value.string) = (char*)(*((yystack_[5].value.tptr)->value.strfnct_dd))((yystack_[3].value.val), (yystack_[1].value.val));
	  else
	    (yylhs.value.string) = (char*)"";
	}
#line 1027 "apr_parser.cc" // lalr1.cc:859
    break;

  case 58:
#line 292 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[11].value.tptr), (yystack_[11].value.tptr)->value.strfnct_dcccc == nullptr))
	    (yylhs.value.string) = (char*)(*((yystack_[11].value.tptr)->value.strfnct_dcccc))((yystack_[9].value.val), (yystack_[7].value.string), (yystack_[5].value.string), (yystack_[3].value.string), (yystack_[1].value.string));
	  else
	    (yylhs.value.string) = (char*)"";
	}
#line 1038 "apr_parser.cc" // lalr1.cc:859
    break;

  case 59:
#line 298 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[7].value.tptr), (yystack_[7].value.tptr)->value.strfnct_dcc == nullptr))
	    (yylhs.value.string) = (char*)(*((yystack_[7].value.tptr)->value.strfnct_dcc))((yystack_[5].value.val), (yystack_[3].value.string), (yystack_[1].value.string));
	  else
	    (yylhs.value.string) = (char*)"";
	}
#line 1049 "apr_parser.cc" // lalr1.cc:859
    break;

  case 60:
#line 304 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[7].value.tptr), (yystack_[7].value.tptr)->value.strfnct_ccc == nullptr))
	    (yylhs.value.string) = (char*)(*((yystack_[7].value.tptr)->value.strfnct_ccc))((yystack_[5].value.string), (yystack_[3].value.string), (yystack_[1].value.string));
	  else
	    (yylhs.value.string) = (char*)"";
	}
#line 1060 "apr_parser.cc" // lalr1.cc:859
    break;

  case 61:
#line 310 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.string) = ((yystack_[4].value.val)) ? ((yystack_[2].value.string)) : ((yystack_[0].value.string));              }
#line 1066 "apr_parser.cc" // lalr1.cc:859
    break;

  case 62:
#line 312 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[0].value.val); 				}
#line 1072 "apr_parser.cc" // lalr1.cc:859
    break;

  case 63:
#line 313 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[0].value.val) + 1;				}
#line 1078 "apr_parser.cc" // lalr1.cc:859
    break;

  case 64:
#line 314 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[0].value.val) - 1;				}
#line 1084 "apr_parser.cc" // lalr1.cc:859
    break;

  case 65:
#line 315 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;			}
#line 1090 "apr_parser.cc" // lalr1.cc:859
    break;

  case 66:
#line 316 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;			}
#line 1096 "apr_parser.cc" // lalr1.cc:859
    break;

  case 67:
#line 317 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);		}
#line 1102 "apr_parser.cc" // lalr1.cc:859
    break;

  case 68:
#line 318 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);		}
#line 1108 "apr_parser.cc" // lalr1.cc:859
    break;

  case 69:
#line 319 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;		}
#line 1114 "apr_parser.cc" // lalr1.cc:859
    break;

  case 70:
#line 320 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;		}
#line 1120 "apr_parser.cc" // lalr1.cc:859
    break;

  case 71:
#line 321 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[0].value.val); (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
				  redefined_warning(aprepro, (yystack_[2].value.tptr));          }
#line 1127 "apr_parser.cc" // lalr1.cc:859
    break;

  case 72:
#line 323 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[0].value.val); (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
				  redefined_warning(aprepro, (yystack_[2].value.tptr));          
		                  set_type(aprepro, (yystack_[2].value.tptr), token::VAR);			}
#line 1135 "apr_parser.cc" // lalr1.cc:859
    break;

  case 73:
#line 326 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val); (yylhs.value.val) = (yystack_[2].value.tptr)->value.var; }
#line 1141 "apr_parser.cc" // lalr1.cc:859
    break;

  case 74:
#line 327 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val); (yylhs.value.val) = (yystack_[2].value.tptr)->value.var; }
#line 1147 "apr_parser.cc" // lalr1.cc:859
    break;

  case 75:
#line 328 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val); (yylhs.value.val) = (yystack_[2].value.tptr)->value.var; }
#line 1153 "apr_parser.cc" // lalr1.cc:859
    break;

  case 76:
#line 329 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val); (yylhs.value.val) = (yystack_[2].value.tptr)->value.var; }
#line 1159 "apr_parser.cc" // lalr1.cc:859
    break;

  case 77:
#line 330 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { errno = 0;
				  (yystack_[2].value.tptr)->value.var = std::pow((yystack_[2].value.tptr)->value.var,(yystack_[0].value.val)); 
				  (yylhs.value.val) = (yystack_[2].value.tptr)->value.var; 
				  SEAMS::math_error(aprepro, "Power");
				}
#line 1169 "apr_parser.cc" // lalr1.cc:859
    break;

  case 78:
#line 335 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[0].value.tptr)); YYERROR; }
#line 1175 "apr_parser.cc" // lalr1.cc:859
    break;

  case 79:
#line 336 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[0].value.tptr)); YYERROR; }
#line 1181 "apr_parser.cc" // lalr1.cc:859
    break;

  case 80:
#line 337 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[1].value.tptr)); YYERROR; }
#line 1187 "apr_parser.cc" // lalr1.cc:859
    break;

  case 81:
#line 338 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[1].value.tptr)); YYERROR; }
#line 1193 "apr_parser.cc" // lalr1.cc:859
    break;

  case 82:
#line 339 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[2].value.tptr)); YYERROR; }
#line 1199 "apr_parser.cc" // lalr1.cc:859
    break;

  case 83:
#line 340 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[2].value.tptr)); YYERROR; }
#line 1205 "apr_parser.cc" // lalr1.cc:859
    break;

  case 84:
#line 341 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[2].value.tptr)); YYERROR; }
#line 1211 "apr_parser.cc" // lalr1.cc:859
    break;

  case 85:
#line 342 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[2].value.tptr)); YYERROR; }
#line 1217 "apr_parser.cc" // lalr1.cc:859
    break;

  case 86:
#line 343 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[2].value.tptr)); YYERROR; }
#line 1223 "apr_parser.cc" // lalr1.cc:859
    break;

  case 87:
#line 344 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[2].value.tptr)); YYERROR; }
#line 1229 "apr_parser.cc" // lalr1.cc:859
    break;

  case 88:
#line 345 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { immutable_modify(aprepro, (yystack_[2].value.tptr)); YYERROR; }
#line 1235 "apr_parser.cc" // lalr1.cc:859
    break;

  case 89:
#line 347 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
				  undefined_error(aprepro, (yystack_[0].value.tptr)->name);          }
#line 1242 "apr_parser.cc" // lalr1.cc:859
    break;

  case 90:
#line 349 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);		
		                  set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
				  undefined_error(aprepro, (yystack_[0].value.tptr)->name);          }
#line 1250 "apr_parser.cc" // lalr1.cc:859
    break;

  case 91:
#line 352 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);		
		                  set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
				  undefined_error(aprepro, (yystack_[0].value.tptr)->name);          }
#line 1258 "apr_parser.cc" // lalr1.cc:859
    break;

  case 92:
#line 355 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;		
		                  set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
				  undefined_error(aprepro, (yystack_[1].value.tptr)->name);          }
#line 1266 "apr_parser.cc" // lalr1.cc:859
    break;

  case 93:
#line 358 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;		
		                  set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
				  undefined_error(aprepro, (yystack_[1].value.tptr)->name);          }
#line 1274 "apr_parser.cc" // lalr1.cc:859
    break;

  case 94:
#line 361 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[0].value.val); (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
		                  set_type(aprepro, (yystack_[2].value.tptr), token::VAR);                      }
#line 1281 "apr_parser.cc" // lalr1.cc:859
    break;

  case 95:
#line 363 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val); (yylhs.value.val) = (yystack_[2].value.tptr)->value.var; 
		                  set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
				  undefined_error(aprepro, (yystack_[2].value.tptr)->name);          }
#line 1289 "apr_parser.cc" // lalr1.cc:859
    break;

  case 96:
#line 366 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val); (yylhs.value.val) = (yystack_[2].value.tptr)->value.var; 
		                  set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
				  undefined_error(aprepro, (yystack_[2].value.tptr)->name);          }
#line 1297 "apr_parser.cc" // lalr1.cc:859
    break;

  case 97:
#line 369 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val); (yylhs.value.val) = (yystack_[2].value.tptr)->value.var; 
		                  set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
				  undefined_error(aprepro, (yystack_[2].value.tptr)->name);          }
#line 1305 "apr_parser.cc" // lalr1.cc:859
    break;

  case 98:
#line 372 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val); (yylhs.value.val) = (yystack_[2].value.tptr)->value.var; 
		                  set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
				  undefined_error(aprepro, (yystack_[2].value.tptr)->name);          }
#line 1313 "apr_parser.cc" // lalr1.cc:859
    break;

  case 99:
#line 375 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { errno = 0;
				  (yystack_[2].value.tptr)->value.var = std::pow((yystack_[2].value.tptr)->value.var,(yystack_[0].value.val)); 
				  (yylhs.value.val) = (yystack_[2].value.tptr)->value.var; 
		                  set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
				  SEAMS::math_error(aprepro, "Power");
				  undefined_error(aprepro, (yystack_[2].value.tptr)->name);          }
#line 1324 "apr_parser.cc" // lalr1.cc:859
    break;

  case 100:
#line 382 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[2].value.tptr), (yystack_[2].value.tptr)->value.fnctptr == nullptr))
	    (yylhs.value.val) = (*((yystack_[2].value.tptr)->value.fnctptr))();
	  else 
	    (yylhs.value.val) = 0.0;
	  }
#line 1335 "apr_parser.cc" // lalr1.cc:859
    break;

  case 101:
#line 389 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[3].value.tptr), (yystack_[3].value.tptr)->value.fnctptr_d == nullptr))
	    (yylhs.value.val) = (*((yystack_[3].value.tptr)->value.fnctptr_d))((yystack_[1].value.val));
	  else
	    (yylhs.value.val) = 0.0;
	  }
#line 1346 "apr_parser.cc" // lalr1.cc:859
    break;

  case 102:
#line 396 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[3].value.tptr), (yystack_[3].value.tptr)->value.fnctptr_c == nullptr))
	    (yylhs.value.val) = (*((yystack_[3].value.tptr)->value.fnctptr_c))((yystack_[1].value.string));
	  else
	    (yylhs.value.val) = 0.0;
	  }
#line 1357 "apr_parser.cc" // lalr1.cc:859
    break;

  case 103:
#line 403 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[3].value.tptr), (yystack_[3].value.tptr)->value.fnctptr_a == nullptr))
	    (yylhs.value.val) = (*((yystack_[3].value.tptr)->value.fnctptr_a))((yystack_[1].value.arrval));
	  else
	    (yylhs.value.val) = 0.0;
	  }
#line 1368 "apr_parser.cc" // lalr1.cc:859
    break;

  case 104:
#line 410 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	    if (arg_check((yystack_[5].value.tptr), (yystack_[5].value.tptr)->value.fnctptr_cd == nullptr))
	      (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cd))((yystack_[3].value.string), (yystack_[1].value.val));
	    else
	      (yylhs.value.val) = 0.0;
	  }
#line 1379 "apr_parser.cc" // lalr1.cc:859
    break;

  case 105:
#line 417 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	  if (arg_check((yystack_[5].value.tptr), (yystack_[5].value.tptr)->value.fnctptr_dc == nullptr))
	    (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dc))((yystack_[3].value.val), (yystack_[1].value.string));
	  else
	    (yylhs.value.val) = 0.0;
	  }
#line 1390 "apr_parser.cc" // lalr1.cc:859
    break;

  case 106:
#line 424 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	    if (arg_check((yystack_[5].value.tptr), (yystack_[5].value.tptr)->value.fnctptr_cc == nullptr))
	      (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cc))((yystack_[3].value.string), (yystack_[1].value.string));
	    else
	      (yylhs.value.val) = 0.0;
	  }
#line 1401 "apr_parser.cc" // lalr1.cc:859
    break;

  case 107:
#line 431 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	    if (arg_check((yystack_[5].value.tptr), (yystack_[5].value.tptr)->value.fnctptr_dd == nullptr))
	      (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dd))((yystack_[3].value.val), (yystack_[1].value.val));
	    else
	      (yylhs.value.val) = 0.0;
	  }
#line 1412 "apr_parser.cc" // lalr1.cc:859
    break;

  case 108:
#line 437 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	    if (arg_check((yystack_[7].value.tptr), (yystack_[7].value.tptr)->value.fnctptr_ddd == nullptr))
	      (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ddd))((yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
	    else
	      (yylhs.value.val) = 0.0;
	  }
#line 1423 "apr_parser.cc" // lalr1.cc:859
    break;

  case 109:
#line 443 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	    if (arg_check((yystack_[7].value.tptr), (yystack_[7].value.tptr)->value.fnctptr_ccd == nullptr))
	      (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccd))((yystack_[5].value.string), (yystack_[3].value.string), (yystack_[1].value.val));
	    else
	      (yylhs.value.val) = 0.0;
	  }
#line 1434 "apr_parser.cc" // lalr1.cc:859
    break;

  case 110:
#line 449 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	    if (arg_check((yystack_[9].value.tptr), (yystack_[9].value.tptr)->value.fnctptr_dddd == nullptr))
	      (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))((yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
	    else
	      (yylhs.value.val) = 0.0;
	  }
#line 1445 "apr_parser.cc" // lalr1.cc:859
    break;

  case 111:
#line 455 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	    if (arg_check((yystack_[9].value.tptr), (yystack_[9].value.tptr)->value.fnctptr_dddd == nullptr))
	      (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))((yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
	    else
	      (yylhs.value.val) = 0.0;
	  }
#line 1456 "apr_parser.cc" // lalr1.cc:859
    break;

  case 112:
#line 461 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	    if (arg_check((yystack_[11].value.tptr), (yystack_[11].value.tptr)->value.fnctptr_ddddc == nullptr))
	      (yylhs.value.val) = (*((yystack_[11].value.tptr)->value.fnctptr_ddddc))((yystack_[9].value.val), (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.string));
	    else
	      (yylhs.value.val) = 0.0;
	  }
#line 1467 "apr_parser.cc" // lalr1.cc:859
    break;

  case 113:
#line 467 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    {
	    if (arg_check((yystack_[13].value.tptr), (yystack_[13].value.tptr)->value.fnctptr_dddddd == nullptr))
	      (yylhs.value.val) = (*((yystack_[13].value.tptr)->value.fnctptr_dddddd))((yystack_[11].value.val), (yystack_[9].value.val), (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
	    else
	      (yylhs.value.val) = 0.0;
	  }
#line 1478 "apr_parser.cc" // lalr1.cc:859
    break;

  case 114:
#line 473 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) + (yystack_[0].value.val); 			}
#line 1484 "apr_parser.cc" // lalr1.cc:859
    break;

  case 115:
#line 474 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) - (yystack_[0].value.val); 			}
#line 1490 "apr_parser.cc" // lalr1.cc:859
    break;

  case 116:
#line 475 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[2].value.val) * (yystack_[0].value.val); 			}
#line 1496 "apr_parser.cc" // lalr1.cc:859
    break;

  case 117:
#line 476 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { if ((yystack_[0].value.val) == 0.)
				    {
				      yyerror(aprepro, "Zero divisor"); 
				      yyerrok;
				    }
				  else
				    (yylhs.value.val) = (yystack_[2].value.val) / (yystack_[0].value.val); 			}
#line 1508 "apr_parser.cc" // lalr1.cc:859
    break;

  case 118:
#line 483 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { if ((yystack_[0].value.val) == 0.)
				    {
				      yyerror(aprepro, "Zero divisor");
				      yyerrok;
				    }
				  else
				    (yylhs.value.val) = (int)(yystack_[2].value.val) % (int)(yystack_[0].value.val);		}
#line 1520 "apr_parser.cc" // lalr1.cc:859
    break;

  case 119:
#line 490 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = -(yystack_[0].value.val);				}
#line 1526 "apr_parser.cc" // lalr1.cc:859
    break;

  case 120:
#line 491 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) =  (yystack_[0].value.val);				}
#line 1532 "apr_parser.cc" // lalr1.cc:859
    break;

  case 121:
#line 492 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { errno = 0;
				  (yylhs.value.val) = std::pow((yystack_[2].value.val), (yystack_[0].value.val)); 
				  SEAMS::math_error(aprepro, "Power");			}
#line 1540 "apr_parser.cc" // lalr1.cc:859
    break;

  case 122:
#line 495 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[1].value.val);				}
#line 1546 "apr_parser.cc" // lalr1.cc:859
    break;

  case 123:
#line 496 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { errno = 0;
				  (yylhs.value.val) = (double)((yystack_[1].value.val) < 0 ? -floor(-((yystack_[1].value.val))): floor((yystack_[1].value.val)) );
				  SEAMS::math_error(aprepro, "floor (int)");		}
#line 1554 "apr_parser.cc" // lalr1.cc:859
    break;

  case 124:
#line 499 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = ((yystack_[0].value.val)) ? 1 : 0; }
#line 1560 "apr_parser.cc" // lalr1.cc:859
    break;

  case 125:
#line 500 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = ((yystack_[4].value.val)) ? ((yystack_[2].value.val)) : ((yystack_[0].value.val));              }
#line 1566 "apr_parser.cc" // lalr1.cc:859
    break;

  case 126:
#line 501 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = array_value((yystack_[5].value.tptr)->value.avar, (yystack_[3].value.val), (yystack_[1].value.val)); }
#line 1572 "apr_parser.cc" // lalr1.cc:859
    break;

  case 127:
#line 503 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:859
    { (yylhs.value.val) = (yystack_[0].value.val);
                                    array *arr = (yystack_[7].value.tptr)->value.avar;
                                    int cols = arr->cols;
                                    int rows = arr->rows;
				    int row = (yystack_[5].value.val);
				    int col = (yystack_[3].value.val);
				    if (aprepro.ap_options.one_based_index) {
				      row--;
				      col--;
				    }
				    if (row < rows && col < cols) {
                                      int offset = row*cols+col;
                                      (yystack_[7].value.tptr)->value.avar->data[offset] = (yystack_[0].value.val);
                                    }
                                    else {
                                      yyerror(aprepro, "Row or Column index out of range"); 
                                      yyerrok;
                                    }
                                  }
#line 1596 "apr_parser.cc" // lalr1.cc:859
    break;


#line 1600 "apr_parser.cc" // lalr1.cc:859
            default:
              break;
            }
        }
      catch (const syntax_error& yyexc)
        {
          error (yyexc);
          YYERROR;
        }
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;
      YY_STACK_PRINT ();

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, yylhs);
    }
    goto yynewstate;

  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        error (yysyntax_error_ (yystack_[0].state, yyla));
      }


    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.type_get () == yyeof_)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
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
    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    goto yyerrlab1;

  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    {
      stack_symbol_type error_token;
      for (;;)
        {
          yyn = yypact_[yystack_[0].state];
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

          // Pop the current state because it cannot handle the error token.
          if (yystack_.size () == 1)
            YYABORT;

          yy_destroy_ ("Error: popping", yystack_[0]);
          yypop_ ();
          YY_STACK_PRINT ();
        }


      // Shift the error token.
      error_token.state = yyn;
      yypush_ ("Shifting", error_token);
    }
    goto yynewstate;

    // Accept.
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;

    // Abort.
  yyabortlab:
    yyresult = 1;
    goto yyreturn;

  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack"
                 << '\n';
        // Do not try to display the values of the reclaimed symbols,
        // as their printer might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
  }

  void
  Parser::error (const syntax_error& yyexc)
  {
    error (yyexc.what());
  }

  // Generate an error message.
  std::string
  Parser::yysyntax_error_ (state_type yystate, const symbol_type& yyla) const
  {
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
       - The only way there can be no lookahead present (in yyla) is
         if this state is a consistent state with a default action.
         Thus, detecting the absence of a lookahead is sufficient to
         determine that there is no unexpected or expected token to
         report.  In that case, just report a simple "syntax error".
       - Don't assume there isn't a lookahead just because this state is
         a consistent state with a default action.  There might have
         been a previous inconsistent state, consistent state with a
         non-default action, or user semantic action that manipulated
         yyla.  (However, yyla is currently not documented for users.)
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
    if (!yyla.empty ())
      {
        int yytoken = yyla.type_get ();
        yyarg[yycount++] = yytname_[yytoken];
        int yyn = yypact_[yystate];
        if (!yy_pact_value_is_default_ (yyn))
          {
            /* Start YYX at -YYN if negative to avoid negative indexes in
               YYCHECK.  In other words, skip the first -YYN actions for
               this state because they are default actions.  */
            int yyxbegin = yyn < 0 ? -yyn : 0;
            // Stay within bounds of both yycheck and yytname.
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

    char const* yyformat = YY_NULLPTR;
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

    std::string yyres;
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


  const signed char Parser::yypact_ninf_ = -25;

  const signed char Parser::yytable_ninf_ = -1;

  const short int
  Parser::yypact_[] =
  {
     -25,    22,   -25,   -13,   258,   -25,   -25,   -25,   -25,   -25,
     574,   604,    -8,   634,    -5,    63,     6,    10,    27,   417,
     417,   -25,   417,   372,   417,    -2,   210,    48,   127,    91,
    1069,   372,   417,   417,   417,   417,   417,   -25,   -25,   417,
     417,   417,   417,   417,   417,   -25,   -25,   417,   417,   417,
     417,   417,   417,   417,   -25,   -25,   417,   417,   372,   312,
     357,   372,   664,    16,   417,   171,   309,   807,  1021,    42,
     -25,    42,    42,   -25,   -25,   -25,   -25,   -25,   -25,   -25,
     -25,   417,   417,   417,   -25,   372,   372,   417,   372,   -25,
     417,   417,   417,   417,   417,   417,   417,   -25,   417,   417,
     417,   417,   417,   417,   417,   417,   417,   417,   417,   372,
     417,   417,    29,   309,  1102,  1118,  1118,  1118,  1118,  1118,
     309,  1118,  1118,  1118,  1118,  1118,  1118,   309,  1118,   309,
    1118,  1118,  1118,  1118,  1118,  1118,   309,  1118,   722,    29,
    1102,   -25,    34,   121,   573,   -25,    76,   148,   603,    93,
     175,   633,   417,    42,   -25,   -25,   417,   -25,   -24,  1086,
      20,  1118,   -25,    18,    18,  1134,   -25,  1134,    45,    45,
      45,    45,    45,    45,   -25,  1149,  1163,   214,   214,   214,
     214,   214,   214,   104,   104,    42,   -25,    42,    42,    42,
     417,   -25,   417,   -25,   417,   -25,   -25,   417,   -25,   417,
     -25,   -25,   417,   -25,   417,   -25,  1118,    42,   417,   417,
    1046,   203,   834,   212,   543,   403,   439,   861,   469,   888,
     915,   309,  1118,    68,   417,   -25,   -25,   -25,   417,   -25,
     417,   417,   417,   -25,   -25,   -25,   -25,   417,   942,   663,
     751,   491,   433,  1118,   -25,   417,   -25,   417,   -25,   417,
     -25,   693,   969,   463,   417,   -25,   -25,   417,   497,   780,
     519,   -25,   417,   -25,   996,   -25
  };

  const unsigned char
  Parser::yydefact_[] =
  {
       2,     0,     1,     0,     0,     4,     3,     9,    62,    44,
      89,    65,    45,    66,    46,    28,     0,     0,     0,     0,
       0,     8,     0,     0,     0,     0,     0,   124,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    92,    93,     0,
       0,     0,     0,     0,     0,    69,    70,     0,     0,     0,
       0,     0,     0,     0,    80,    81,     0,     0,     0,     0,
       0,     0,    89,     0,     0,   124,     0,     0,     0,   120,
      38,   119,    12,    63,    90,    67,    78,    64,    91,    68,
      79,     0,     0,     0,     7,     0,     0,     0,     0,     6,
       0,     0,     0,     0,     0,     0,     0,     5,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    36,    47,    94,    95,    96,    97,    98,    99,
      49,    71,    73,    74,    75,    76,    77,    48,    72,    51,
      82,    84,    85,    86,    87,    88,    50,    83,     0,    35,
       0,   100,     0,     0,     0,    53,     0,     0,     0,     0,
       0,     0,     0,   119,    21,   122,     0,   123,     0,     0,
      19,     0,    20,    37,    39,    41,    43,    40,    22,    23,
      24,    25,    26,    27,    56,    17,    18,    10,    11,    13,
      14,    15,    16,   114,   115,   117,    42,   116,   118,   121,
       0,   103,     0,   102,     0,   101,    55,     0,    52,     0,
      54,    34,     0,    29,     0,    33,    94,   116,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    61,   125,   126,     0,   106,   104,   105,     0,   107,
       0,     0,     0,    57,    31,    30,    32,     0,     0,     0,
       0,     0,     0,   127,   109,     0,   108,     0,    60,     0,
      59,     0,     0,     0,     0,   111,   110,     0,     0,     0,
       0,   112,     0,    58,     0,   113
  };

  const signed char
  Parser::yypgoto_[] =
  {
     -25,   -25,   -25,   -17,     3,    82,    -4
  };

  const short int
  Parser::yydefgoto_[] =
  {
      -1,     1,     6,    27,    28,    66,   161
  };

  const unsigned short int
  Parser::yytable_[] =
  {
      30,    73,    65,    74,    75,   208,    76,     7,    90,    91,
      92,    93,    94,    95,    47,    67,    68,    56,    69,    71,
      72,    59,     2,     3,    96,    60,    70,   114,   115,   116,
     117,   118,   119,    57,   112,   121,   122,   123,   124,   125,
     126,     4,    61,   128,   130,   131,   132,   133,   134,   135,
     191,    83,   137,   138,   140,   144,   148,   151,    87,    88,
     153,   139,   142,   146,   149,   160,   162,    85,    86,    87,
      88,     5,    85,    86,    87,    88,    81,   159,    82,    83,
      57,   140,   140,   165,   167,    58,    29,   111,   163,   164,
     237,   166,   196,    96,   175,   176,   177,   178,   179,   180,
     181,   182,   183,   184,   185,   187,   188,   189,     0,   201,
       0,    89,   186,   113,    85,    86,    87,    88,     0,     0,
       0,   120,     0,    90,    91,    92,    93,    94,    95,   127,
     129,    85,    86,    87,    88,   192,     0,   193,   136,    96,
       0,   143,   147,   150,   108,   156,   110,    84,   206,   111,
       0,     0,   207,    90,    91,    92,    93,    94,    95,     0,
       0,     0,   197,   158,   198,    85,    86,    87,    88,    96,
       0,     0,   168,   169,   170,   171,   172,   173,   174,     0,
      90,    91,    92,    93,    94,    95,   210,   154,   212,   202,
     214,   203,     0,     0,     0,   217,    96,     0,   219,    81,
     220,    82,    83,     0,     0,   222,     0,    90,    91,    92,
      93,    94,    95,    77,     0,    78,    79,   224,    80,   225,
     238,     0,     0,    96,   239,     0,   240,     0,   227,     0,
       0,     0,     0,   243,   113,    90,    91,    92,    93,    94,
      95,   251,     0,   252,    90,    91,    92,    93,    94,    95,
     259,    96,   106,   107,   108,   156,   110,     0,   264,   111,
      96,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,     0,    19,   211,    20,   213,     0,    21,   215,
       0,   216,     0,     0,   218,     0,     0,     0,     0,     0,
     221,     0,     0,     0,     0,     0,    22,    23,     0,     0,
       0,     0,    24,     0,    25,    26,     0,     0,     0,     0,
       0,     0,     0,   241,   242,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,     0,    19,   141,    20,
       0,   253,     0,     0,     0,     0,   258,     0,     0,   260,
       0,    90,    91,    92,    93,    94,    95,     0,     0,     0,
      22,    23,     0,     0,     0,     0,    24,    96,    25,    26,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,     0,    19,   145,    20,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,     0,    19,     0,    20,
       0,     0,     0,     0,     0,    22,    23,     0,     0,     0,
       0,    24,     0,    25,    26,     0,     0,     0,     0,     0,
      22,    23,     0,     0,     0,     0,    24,   231,    25,    26,
       8,     9,    62,    11,    12,    13,    14,    63,    16,    17,
       0,     0,    19,     0,    20,    90,    91,    92,    93,    94,
      95,     0,     0,     0,     0,     0,     0,   249,     0,   250,
       0,    96,     0,   232,     0,    22,    64,     0,     0,     0,
       0,    24,     0,    25,    26,    90,    91,    92,    93,    94,
      95,    90,    91,    92,    93,    94,    95,   257,     0,     0,
       0,    96,     0,     0,     0,   234,     0,    96,     0,     0,
       0,     0,     0,     0,     0,    90,    91,    92,    93,    94,
      95,    90,    91,    92,    93,    94,    95,   248,     0,     0,
       0,    96,     0,   261,     0,     0,     0,    96,     0,     0,
       0,     0,     0,    90,    91,    92,    93,    94,    95,    90,
      91,    92,    93,    94,    95,   263,     0,     0,     0,    96,
       0,     0,     0,     0,     0,    96,     0,     0,     0,     0,
       0,    90,    91,    92,    93,    94,    95,   228,     0,   229,
       0,     0,     0,     0,   230,     0,     0,    96,     0,     0,
       0,     0,     0,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   156,   110,     0,   194,   111,   195,
       0,     0,     0,     0,     0,     0,    31,    32,    33,    34,
      35,    36,     0,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,     0,   199,   111,   200,
      37,    38,     0,     0,     0,     0,    39,    40,    41,    42,
      43,    44,     0,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,     0,   204,   111,   205,
      45,    46,     0,     0,     0,     0,    48,    49,    50,    51,
      52,    53,     0,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,     0,   245,   111,   246,
      54,    55,     0,     0,     0,     0,   152,    32,    33,    34,
      35,    36,     0,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   156,   110,     0,   254,   111,   255,
      37,    38,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   156,   110,   190,     0,   111,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   156,   110,   247,     0,   111,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     107,   108,   156,   110,   262,     0,   111,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      98,    99,   100,   101,   102,   103,   104,   105,   106,   107,
     108,   156,   110,   155,     0,   111,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    98,    99,   100,
     101,   102,   103,   104,   105,   106,   107,   108,   156,   110,
     226,     0,   111,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   156,   110,   233,     0,   111,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     107,   108,   156,   110,   235,     0,   111,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,   156,
     110,   236,     0,   111,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    98,    99,   100,   101,   102,
     103,   104,   105,   106,   107,   108,   156,   110,   244,     0,
     111,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   156,   110,   256,     0,   111,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     156,   110,   265,     0,   111,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   156,   110,   157,
       0,   111,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     107,   108,   156,   110,   223,     0,   111,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   156,   110,    97,
       0,   111,     0,     0,     0,     0,     0,     0,     0,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     109,   110,     0,     0,   111,   209,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   156,   110,     0,
       0,   111,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   109,   110,     0,     0,   111,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,   156,
     110,     0,     0,   111,    98,    99,   100,   101,   102,   103,
     104,   105,     0,     0,     0,     0,   110,     0,     0,   111,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     156,   110,     0,     0,   111,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   156,   110,     0,     0,   111
  };

  const short int
  Parser::yycheck_[] =
  {
       4,     3,    19,     5,     6,    29,     8,    20,    32,    33,
      34,    35,    36,    37,    22,    19,    20,    22,    22,    23,
      24,    15,     0,     1,    48,    15,    23,    31,    32,    33,
      34,    35,    36,    17,    31,    39,    40,    41,    42,    43,
      44,    19,    15,    47,    48,    49,    50,    51,    52,    53,
      16,    31,    56,    57,    58,    59,    60,    61,    40,    41,
      64,    58,    59,    60,    61,    82,    83,    38,    39,    40,
      41,    49,    38,    39,    40,    41,    28,    81,    30,    31,
      17,    85,    86,    87,    88,    22,     4,    45,    85,    86,
      22,    88,    16,    48,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,   111,    -1,    16,
      -1,    20,   109,    31,    38,    39,    40,    41,    -1,    -1,
      -1,    39,    -1,    32,    33,    34,    35,    36,    37,    47,
      48,    38,    39,    40,    41,    14,    -1,    16,    56,    48,
      -1,    59,    60,    61,    40,    41,    42,    20,   152,    45,
      -1,    -1,   156,    32,    33,    34,    35,    36,    37,    -1,
      -1,    -1,    14,    81,    16,    38,    39,    40,    41,    48,
      -1,    -1,    90,    91,    92,    93,    94,    95,    96,    -1,
      32,    33,    34,    35,    36,    37,   190,    16,   192,    14,
     194,    16,    -1,    -1,    -1,   199,    48,    -1,   202,    28,
     204,    30,    31,    -1,    -1,   209,    -1,    32,    33,    34,
      35,    36,    37,     3,    -1,     5,     6,    14,     8,    16,
     224,    -1,    -1,    48,   228,    -1,   230,    -1,    16,    -1,
      -1,    -1,    -1,   237,   152,    32,    33,    34,    35,    36,
      37,   245,    -1,   247,    32,    33,    34,    35,    36,    37,
     254,    48,    38,    39,    40,    41,    42,    -1,   262,    45,
      48,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    -1,    15,   192,    17,   194,    -1,    20,   197,
      -1,   199,    -1,    -1,   202,    -1,    -1,    -1,    -1,    -1,
     208,    -1,    -1,    -1,    -1,    -1,    38,    39,    -1,    -1,
      -1,    -1,    44,    -1,    46,    47,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   231,   232,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    13,    -1,    15,    16,    17,
      -1,   249,    -1,    -1,    -1,    -1,   254,    -1,    -1,   257,
      -1,    32,    33,    34,    35,    36,    37,    -1,    -1,    -1,
      38,    39,    -1,    -1,    -1,    -1,    44,    48,    46,    47,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    -1,    15,    16,    17,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    13,    -1,    15,    -1,    17,
      -1,    -1,    -1,    -1,    -1,    38,    39,    -1,    -1,    -1,
      -1,    44,    -1,    46,    47,    -1,    -1,    -1,    -1,    -1,
      38,    39,    -1,    -1,    -1,    -1,    44,    14,    46,    47,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      -1,    -1,    15,    -1,    17,    32,    33,    34,    35,    36,
      37,    -1,    -1,    -1,    -1,    -1,    -1,    14,    -1,    16,
      -1,    48,    -1,    14,    -1,    38,    39,    -1,    -1,    -1,
      -1,    44,    -1,    46,    47,    32,    33,    34,    35,    36,
      37,    32,    33,    34,    35,    36,    37,    14,    -1,    -1,
      -1,    48,    -1,    -1,    -1,    16,    -1,    48,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    32,    33,    34,    35,    36,
      37,    32,    33,    34,    35,    36,    37,    16,    -1,    -1,
      -1,    48,    -1,    16,    -1,    -1,    -1,    48,    -1,    -1,
      -1,    -1,    -1,    32,    33,    34,    35,    36,    37,    32,
      33,    34,    35,    36,    37,    16,    -1,    -1,    -1,    48,
      -1,    -1,    -1,    -1,    -1,    48,    -1,    -1,    -1,    -1,
      -1,    32,    33,    34,    35,    36,    37,    14,    -1,    16,
      -1,    -1,    -1,    -1,    21,    -1,    -1,    48,    -1,    -1,
      -1,    -1,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    -1,    14,    45,    16,
      -1,    -1,    -1,    -1,    -1,    -1,    22,    23,    24,    25,
      26,    27,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    -1,    14,    45,    16,
      46,    47,    -1,    -1,    -1,    -1,    22,    23,    24,    25,
      26,    27,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    -1,    14,    45,    16,
      46,    47,    -1,    -1,    -1,    -1,    22,    23,    24,    25,
      26,    27,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    -1,    14,    45,    16,
      46,    47,    -1,    -1,    -1,    -1,    22,    23,    24,    25,
      26,    27,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    -1,    14,    45,    16,
      46,    47,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    14,    -1,    45,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    14,    -1,    45,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    14,    -1,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    16,    -1,    45,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      16,    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
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
      41,    42,    16,    -1,    45,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    18,
      -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    18,    -1,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    20,
      -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    -1,    -1,    45,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    -1,
      -1,    45,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    -1,    -1,    45,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    -1,    -1,    45,    30,    31,    32,    33,    34,    35,
      36,    37,    -1,    -1,    -1,    -1,    42,    -1,    -1,    45,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    -1,    -1,    45,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    -1,    -1,    45
  };

  const unsigned char
  Parser::yystos_[] =
  {
       0,    51,     0,     1,    19,    49,    52,    20,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    15,
      17,    20,    38,    39,    44,    46,    47,    53,    54,    55,
      56,    22,    23,    24,    25,    26,    27,    46,    47,    22,
      23,    24,    25,    26,    27,    46,    47,    22,    22,    23,
      24,    25,    26,    27,    46,    47,    22,    17,    22,    15,
      15,    15,     5,    10,    39,    53,    55,    56,    56,    56,
      54,    56,    56,     3,     5,     6,     8,     3,     5,     6,
       8,    28,    30,    31,    20,    38,    39,    40,    41,    20,
      32,    33,    34,    35,    36,    37,    48,    20,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    45,    54,    55,    56,    56,    56,    56,    56,    56,
      55,    56,    56,    56,    56,    56,    56,    55,    56,    55,
      56,    56,    56,    56,    56,    56,    55,    56,    56,    54,
      56,    16,    54,    55,    56,    16,    54,    55,    56,    54,
      55,    56,    22,    56,    16,    16,    41,    18,    55,    56,
      53,    56,    53,    54,    54,    56,    54,    56,    55,    55,
      55,    55,    55,    55,    55,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    54,    56,    56,    56,
      14,    16,    14,    16,    14,    16,    16,    14,    16,    14,
      16,    16,    14,    16,    14,    16,    56,    56,    29,    29,
      56,    55,    56,    55,    56,    55,    55,    56,    55,    56,
      56,    55,    56,    18,    14,    16,    16,    16,    14,    16,
      21,    14,    14,    16,    16,    16,    16,    22,    56,    56,
      56,    55,    55,    56,    16,    14,    16,    14,    16,    14,
      16,    56,    56,    55,    14,    16,    16,    14,    55,    56,
      55,    16,    14,    16,    56,    16
  };

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



  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a yyntokens_, nonterminals.
  const char*
  const Parser::yytname_[] =
  {
  "\"end of file\"", "error", "$undefined", "NUM", "QSTRING", "UNDVAR",
  "VAR", "SVAR", "IMMVAR", "IMMSVAR", "AVAR", "FNCT", "SFNCT", "AFNCT",
  "COMMA", "LPAR", "RPAR", "LBRACK", "RBRACK", "LBRACE", "RBRACE", "SEMI",
  "EQUAL", "EQ_PLUS", "EQ_MINUS", "EQ_TIME", "EQ_DIV", "EQ_POW", "QUEST",
  "COLON", "LOR", "LAND", "LT", "GT", "LE", "GE", "EQ", "NE", "PLU", "SUB",
  "DIV", "TIM", "MOD", "UNARY", "NOT", "POW", "INC", "DEC", "CONCAT",
  "'\\n'", "$accept", "input", "line", "bool", "aexp", "sexp", "exp", YY_NULLPTR
  };

#if YYDEBUG
  const unsigned short int
  Parser::yyrline_[] =
  {
       0,   133,   133,   134,   137,   138,   145,   149,   150,   151,
     154,   155,   156,   157,   158,   159,   160,   161,   162,   163,
     164,   165,   168,   169,   170,   171,   172,   173,   175,   176,
     182,   188,   194,   200,   206,   212,   215,   217,   225,   227,
     235,   236,   237,   238,   247,   248,   249,   250,   252,   255,
     259,   260,   261,   267,   273,   279,   285,   286,   292,   298,
     304,   310,   312,   313,   314,   315,   316,   317,   318,   319,
     320,   321,   323,   326,   327,   328,   329,   330,   335,   336,
     337,   338,   339,   340,   341,   342,   343,   344,   345,   347,
     349,   352,   355,   358,   361,   363,   366,   369,   372,   375,
     382,   389,   396,   403,   410,   417,   424,   431,   437,   443,
     449,   455,   461,   467,   473,   474,   475,   476,   483,   490,
     491,   492,   495,   496,   499,   500,   501,   502
  };

  // Print the state stack on the debug stream.
  void
  Parser::yystack_print_ ()
  {
    *yycdebug_ << "Stack now";
    for (const auto & elem : yystack_)
      *yycdebug_ << ' ' << elem.state;
    *yycdebug_ << '\n';
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void
  Parser::yy_reduce_print_ (int yyrule)
  {
    unsigned int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):" << '\n';
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // YYDEBUG

  // Symbol number corresponding to token number t.
  inline
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
    const unsigned int user_token_number_max_ = 303;
    const token_number_type undef_token_ = 2;

    if (static_cast<int>(t) <= yyeof_)
      return yyeof_;
    else if (static_cast<unsigned int> (t) <= user_token_number_max_)
      return translate_table[t];
    else
      return undef_token_;
  }


} // SEAMS
#line 2368 "apr_parser.cc" // lalr1.cc:1167
#line 525 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy" // lalr1.cc:1168


void SEAMS::Parser::error(const std::string& m)
{
    aprepro.error(m);
}

