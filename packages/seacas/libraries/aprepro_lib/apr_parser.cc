// A Bison parser, made by GNU Bison 3.3.2.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2019 Free Software Foundation, Inc.

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

// Undocumented macros, especially those whose name start with YY_,
// are private implementation details.  Do not rely on them.

// Take the name prefix into account.
#define yylex SEAMSlex

// First part of user prologue.
#line 33 "aprepro.yy" // lalr1.cc:429

#include "apr_array.h"
#include "apr_util.h"
#include "aprepro.h"

#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdlib.h>

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
} // namespace

namespace SEAMS {
  extern bool echo;
}

#line 74 "apr_parser.cc" // lalr1.cc:429

#include "aprepro_parser.h"

// Second part of user prologue.
#line 130 "aprepro.yy" // lalr1.cc:434

#include "apr_scanner.h"
#include "aprepro.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex

#line 92 "apr_parser.cc" // lalr1.cc:434

#ifndef YY_
#if defined YYENABLE_NLS && YYENABLE_NLS
#if ENABLE_NLS
#include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#define YY_(msgid) dgettext("bison-runtime", msgid)
#endif
#endif
#ifndef YY_
#define YY_(msgid) msgid
#endif
#endif

// Whether we are compiled with exception support.
#ifndef YY_EXCEPTIONS
#if defined __GNUC__ && !defined __EXCEPTIONS
#define YY_EXCEPTIONS 0
#else
#define YY_EXCEPTIONS 1
#endif
#endif

// Suppress unused-variable warnings by "using" E.
#define YYUSE(E) ((void)(E))

// Enable debugging if requested.
#if SEAMSDEBUG

// A pseudo ostream that takes yydebug_ into account.
#define YYCDEBUG                                                                                   \
  if (yydebug_)                                                                                    \
  (*yycdebug_)

#define YY_SYMBOL_PRINT(Title, Symbol)                                                             \
  do {                                                                                             \
    if (yydebug_) {                                                                                \
      *yycdebug_ << Title << ' ';                                                                  \
      yy_print_(*yycdebug_, Symbol);                                                               \
      *yycdebug_ << '\n';                                                                          \
    }                                                                                              \
  } while (false)

#define YY_REDUCE_PRINT(Rule)                                                                      \
  do {                                                                                             \
    if (yydebug_)                                                                                  \
      yy_reduce_print_(Rule);                                                                      \
  } while (false)

#define YY_STACK_PRINT()                                                                           \
  do {                                                                                             \
    if (yydebug_)                                                                                  \
      yystack_print_();                                                                            \
  } while (false)

#else // !SEAMSDEBUG

#define YYCDEBUG                                                                                   \
  if (false)                                                                                       \
  std::cerr
#define YY_SYMBOL_PRINT(Title, Symbol) YYUSE(Symbol)
#define YY_REDUCE_PRINT(Rule) static_cast<void>(0)
#define YY_STACK_PRINT() static_cast<void>(0)

#endif // !SEAMSDEBUG

#define yyerrok (yyerrstatus_ = 0)
#define yyclearin (yyla.clear())

#define YYACCEPT goto yyacceptlab
#define YYABORT goto yyabortlab
#define YYERROR goto yyerrorlab
#define YYRECOVERING() (!!yyerrstatus_)

namespace SEAMS {
#line 168 "apr_parser.cc" // lalr1.cc:510

  /* Return YYSTR after stripping away unnecessary quotes and
     backslashes, so that it's suitable for yyerror.  The heuristic is
     that double-quoting is unnecessary unless the string contains an
     apostrophe, a comma, or backslash (other than backslash-backslash).
     YYSTR is taken from yytname.  */
  std::string Parser::yytnamerr_(const char *yystr)
  {
    if (*yystr == '"') {
      std::string yyr;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp) {
        case '\'':
        case ',': goto do_not_strip_quotes;

        case '\\':
          if (*++yyp != '\\')
            goto do_not_strip_quotes;
          else
            goto append;

        append:
        default: yyr += *yyp; break;

        case '"': return yyr;
        }
    do_not_strip_quotes:;
    }

    return yystr;
  }

  /// Build a parser object.
  Parser::Parser(class Aprepro &aprepro_yyarg)
      :
#if SEAMSDEBUG
        yydebug_(false), yycdebug_(&std::cerr),
#endif
        aprepro(aprepro_yyarg)
  {
  }

  Parser::~Parser() {}

  Parser::syntax_error::~syntax_error() YY_NOEXCEPT YY_NOTHROW {}

  /*---------------.
  | Symbol types.  |
  `---------------*/

  // basic_symbol.
#if 201103L <= YY_CPLUSPLUS
  template <typename Base>
  Parser::basic_symbol<Base>::basic_symbol(basic_symbol &&that)
      : Base(std::move(that)), value(std::move(that.value))
  {
  }
#endif

  template <typename Base>
  Parser::basic_symbol<Base>::basic_symbol(const basic_symbol &that) : Base(that), value(that.value)
  {
  }

  /// Constructor for valueless symbols.
  template <typename Base>
  Parser::basic_symbol<Base>::basic_symbol(typename Base::kind_type t) : Base(t), value()
  {
  }

  template <typename Base>
  Parser::basic_symbol<Base>::basic_symbol(typename Base::kind_type t, YY_RVREF(semantic_type) v)
      : Base(t), value(YY_MOVE(v))
  {
  }

  template <typename Base> bool Parser::basic_symbol<Base>::empty() const YY_NOEXCEPT
  {
    return Base::type_get() == empty_symbol;
  }

  template <typename Base> void Parser::basic_symbol<Base>::move(basic_symbol &s)
  {
    super_type::move(s);
    value = YY_MOVE(s.value);
  }

  // by_type.
  Parser::by_type::by_type() : type(empty_symbol) {}

#if 201103L <= YY_CPLUSPLUS
  Parser::by_type::by_type(by_type &&that) : type(that.type) { that.clear(); }
#endif

  Parser::by_type::by_type(const by_type &that) : type(that.type) {}

  Parser::by_type::by_type(token_type t) : type(yytranslate_(t)) {}

  void Parser::by_type::clear() { type = empty_symbol; }

  void Parser::by_type::move(by_type &that)
  {
    type = that.type;
    that.clear();
  }

  int Parser::by_type::type_get() const YY_NOEXCEPT { return type; }

  // by_state.
  Parser::by_state::by_state() YY_NOEXCEPT : state(empty_state) {}

  Parser::by_state::by_state(const by_state &that) YY_NOEXCEPT : state(that.state) {}

  void Parser::by_state::clear() YY_NOEXCEPT { state = empty_state; }

  void Parser::by_state::move(by_state &that)
  {
    state = that.state;
    that.clear();
  }

  Parser::by_state::by_state(state_type s) YY_NOEXCEPT : state(s) {}

  Parser::symbol_number_type Parser::by_state::type_get() const YY_NOEXCEPT
  {
    if (state == empty_state)
      return empty_symbol;
    else
      return yystos_[state];
  }

  Parser::stack_symbol_type::stack_symbol_type() {}

  Parser::stack_symbol_type::stack_symbol_type(YY_RVREF(stack_symbol_type) that)
      : super_type(YY_MOVE(that.state), YY_MOVE(that.value))
  {
#if 201103L <= YY_CPLUSPLUS
    // that is emptied.
    that.state = empty_state;
#endif
  }

  Parser::stack_symbol_type::stack_symbol_type(state_type s, YY_MOVE_REF(symbol_type) that)
      : super_type(s, YY_MOVE(that.value))
  {
    // that is emptied.
    that.type = empty_symbol;
  }

#if YY_CPLUSPLUS < 201103L
  Parser::stack_symbol_type &Parser::stack_symbol_type::operator=(stack_symbol_type &that)
  {
    state = that.state;
    value = that.value;
    // that is emptied.
    that.state = empty_state;
    return *this;
  }
#endif

  template <typename Base>
  void Parser::yy_destroy_(const char *yymsg, basic_symbol<Base> &yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT(yymsg, yysym);

    // User destructor.
    YYUSE(yysym.type_get());
  }

#if SEAMSDEBUG
  template <typename Base>
  void Parser::yy_print_(std::ostream &yyo, const basic_symbol<Base> &yysym) const
  {
    std::ostream &yyoutput = yyo;
    YYUSE(yyoutput);
    symbol_number_type yytype = yysym.type_get();
#if defined __GNUC__ && !defined __clang__ && !defined __ICC &&                                    \
    __GNUC__ * 100 + __GNUC_MINOR__ <= 408
    // Avoid a (spurious) G++ 4.8 warning about "array subscript is
    // below array bounds".
    if (yysym.empty())
      std::abort();
#endif
    yyo << (yytype < yyntokens_ ? "token" : "nterm") << ' ' << yytname_[yytype] << " (";
    YYUSE(yytype);
    yyo << ')';
  }
#endif

  void Parser::yypush_(const char *m, YY_MOVE_REF(stack_symbol_type) sym)
  {
    if (m)
      YY_SYMBOL_PRINT(m, sym);
    yystack_.push(YY_MOVE(sym));
  }

  void Parser::yypush_(const char *m, state_type s, YY_MOVE_REF(symbol_type) sym)
  {
#if 201103L <= YY_CPLUSPLUS
    yypush_(m, stack_symbol_type(s, std::move(sym)));
#else
    stack_symbol_type ss(s, sym);
    yypush_(m, ss);
#endif
  }

  void Parser::yypop_(int n) { yystack_.pop(n); }

#if SEAMSDEBUG
  std::ostream &Parser::debug_stream() const { return *yycdebug_; }

  void Parser::set_debug_stream(std::ostream &o) { yycdebug_ = &o; }

  Parser::debug_level_type Parser::debug_level() const { return yydebug_; }

  void Parser::set_debug_level(debug_level_type l) { yydebug_ = l; }
#endif // SEAMSDEBUG

  Parser::state_type Parser::yy_lr_goto_state_(state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - yyntokens_] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - yyntokens_];
  }

  bool Parser::yy_pact_value_is_default_(int yyvalue) { return yyvalue == yypact_ninf_; }

  bool Parser::yy_table_value_is_error_(int yyvalue) { return yyvalue == yytable_ninf_; }

  int Parser::operator()() { return parse(); }

  int Parser::parse()
  {
    // State.
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_     = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The return value of parse ().
    int yyresult;

#if YY_EXCEPTIONS
    try
#endif // YY_EXCEPTIONS
    {
      YYCDEBUG << "Starting parse\n";

      /* Initialize the stack.  The initial state will be set in
         yynewstate, since the latter expects the semantical and the
         location values to have been already stored, initialize these
         stacks with a primary value.  */
      yystack_.clear();
      yypush_(YY_NULLPTR, 0, YY_MOVE(yyla));

    /*-----------------------------------------------.
    | yynewstate -- push a new symbol on the stack.  |
    `-----------------------------------------------*/
    yynewstate:
      YYCDEBUG << "Entering state " << yystack_[0].state << '\n';

      // Accept?
      if (yystack_[0].state == yyfinal_)
        YYACCEPT;

      goto yybackup;

    /*-----------.
    | yybackup.  |
    `-----------*/
    yybackup:
      // Try to take a decision without lookahead.
      yyn = yypact_[yystack_[0].state];
      if (yy_pact_value_is_default_(yyn))
        goto yydefault;

      // Read a lookahead token.
      if (yyla.empty()) {
        YYCDEBUG << "Reading a token: ";
#if YY_EXCEPTIONS
        try
#endif // YY_EXCEPTIONS
        {
          yyla.type = yytranslate_(yylex(&yyla.value));
        }
#if YY_EXCEPTIONS
        catch (const syntax_error &yyexc) {
          YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
          error(yyexc);
          goto yyerrlab1;
        }
#endif // YY_EXCEPTIONS
      }
      YY_SYMBOL_PRINT("Next token is", yyla);

      /* If the proper action on seeing token YYLA.TYPE is to reduce or
         to detect an error, take that action.  */
      yyn += yyla.type_get();
      if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.type_get())
        goto yydefault;

      // Reduce or error.
      yyn = yytable_[yyn];
      if (yyn <= 0) {
        if (yy_table_value_is_error_(yyn))
          goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
      }

      // Count tokens shifted since error; after three, turn off error status.
      if (yyerrstatus_)
        --yyerrstatus_;

      // Shift the lookahead token.
      yypush_("Shifting", yyn, YY_MOVE(yyla));
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
    | yyreduce -- do a reduction.  |
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
        YY_REDUCE_PRINT(yyn);
#if YY_EXCEPTIONS
        try
#endif // YY_EXCEPTIONS
        {
          switch (yyn) {
          case 4:
#line 149 "aprepro.yy" // lalr1.cc:919
          {
            if (echo)
              aprepro.lexer->LexerOutput("\n", 1);
          }
#line 636 "apr_parser.cc" // lalr1.cc:919
          break;

          case 5:
#line 150 "aprepro.yy" // lalr1.cc:919
          {
            if (echo) {
              static char    tmpstr[512];
              SEAMS::symrec *format = aprepro.getsym("_FORMAT");
              int len = sprintf(tmpstr, format->value.svar.c_str(), (yystack_[1].value.val));
              aprepro.lexer->LexerOutput(tmpstr, len);
            }
          }
#line 648 "apr_parser.cc" // lalr1.cc:919
          break;

          case 6:
#line 157 "aprepro.yy" // lalr1.cc:919
          {
            if (echo && (yystack_[1].value.string) != NULL) {
              aprepro.lexer->LexerOutput((yystack_[1].value.string),
                                         strlen((yystack_[1].value.string)));
            }
          }
#line 657 "apr_parser.cc" // lalr1.cc:919
          break;

          case 7:
#line 161 "aprepro.yy" // lalr1.cc:919
          {
          }
#line 663 "apr_parser.cc" // lalr1.cc:919
          break;

          case 8:
#line 162 "aprepro.yy" // lalr1.cc:919
          {
          }
#line 669 "apr_parser.cc" // lalr1.cc:919
          break;

          case 9:
#line 163 "aprepro.yy" // lalr1.cc:919
          {
            yyerrok;
          }
#line 675 "apr_parser.cc" // lalr1.cc:919
          break;

          case 10:
#line 166 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) < (yystack_[0].value.val);
          }
#line 681 "apr_parser.cc" // lalr1.cc:919
          break;

          case 11:
#line 167 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) > (yystack_[0].value.val);
          }
#line 687 "apr_parser.cc" // lalr1.cc:919
          break;

          case 12:
#line 168 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = !((yystack_[0].value.val));
          }
#line 693 "apr_parser.cc" // lalr1.cc:919
          break;

          case 13:
#line 169 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) <= (yystack_[0].value.val);
          }
#line 699 "apr_parser.cc" // lalr1.cc:919
          break;

          case 14:
#line 170 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) >= (yystack_[0].value.val);
          }
#line 705 "apr_parser.cc" // lalr1.cc:919
          break;

          case 15:
#line 171 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) == (yystack_[0].value.val);
          }
#line 711 "apr_parser.cc" // lalr1.cc:919
          break;

          case 16:
#line 172 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) != (yystack_[0].value.val);
          }
#line 717 "apr_parser.cc" // lalr1.cc:919
          break;

          case 17:
#line 173 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 723 "apr_parser.cc" // lalr1.cc:919
          break;

          case 18:
#line 174 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 729 "apr_parser.cc" // lalr1.cc:919
          break;

          case 19:
#line 175 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 735 "apr_parser.cc" // lalr1.cc:919
          break;

          case 20:
#line 176 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 741 "apr_parser.cc" // lalr1.cc:919
          break;

          case 21:
#line 177 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          }
#line 747 "apr_parser.cc" // lalr1.cc:919
          break;

          case 22:
#line 180 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) < 0 ? 1 : 0);
          }
#line 753 "apr_parser.cc" // lalr1.cc:919
          break;

          case 23:
#line 181 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) > 0 ? 1 : 0);
          }
#line 759 "apr_parser.cc" // lalr1.cc:919
          break;

          case 24:
#line 182 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) <= 0 ? 1 : 0);
          }
#line 765 "apr_parser.cc" // lalr1.cc:919
          break;

          case 25:
#line 183 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) >= 0 ? 1 : 0);
          }
#line 771 "apr_parser.cc" // lalr1.cc:919
          break;

          case 26:
#line 184 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) == 0 ? 1 : 0);
          }
#line 777 "apr_parser.cc" // lalr1.cc:919
          break;

          case 27:
#line 185 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) != 0 ? 1 : 0);
          }
#line 783 "apr_parser.cc" // lalr1.cc:919
          break;

          case 28:
#line 187 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) = (yystack_[0].value.tptr)->value.avar;
          }
#line 789 "apr_parser.cc" // lalr1.cc:919
          break;

          case 29:
#line 188 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_c == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_c))((yystack_[1].value.string));
            else
              yyerrok;
          }
#line 800 "apr_parser.cc" // lalr1.cc:919
          break;

          case 30:
#line 194 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_cd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 811 "apr_parser.cc" // lalr1.cc:919
          break;

          case 31:
#line 200 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_cc == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              yyerrok;
          }
#line 822 "apr_parser.cc" // lalr1.cc:919
          break;

          case 32:
#line 206 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.arrfnct_ddd == NULL))
              (yylhs.value.arrval) = (*((yystack_[7].value.tptr)->value.arrfnct_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 833 "apr_parser.cc" // lalr1.cc:919
          break;

          case 33:
#line 212 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_dd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 844 "apr_parser.cc" // lalr1.cc:919
          break;

          case 34:
#line 218 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_d == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_d))((yystack_[1].value.val));
            else
              yyerrok;
          }
#line 855 "apr_parser.cc" // lalr1.cc:919
          break;

          case 35:
#line 224 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_a == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_a))((yystack_[1].value.arrval));
            else
              yyerrok;
          }
#line 866 "apr_parser.cc" // lalr1.cc:919
          break;

          case 36:
#line 230 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) = (yystack_[0].value.arrval);
            delete (yystack_[2].value.tptr)->value.avar;
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 874 "apr_parser.cc" // lalr1.cc:919
          break;

          case 37:
#line 233 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 881 "apr_parser.cc" // lalr1.cc:919
          break;

          case 38:
#line 235 "aprepro.yy" // lalr1.cc:919
          {
            if ((yystack_[2].value.arrval)->cols == (yystack_[0].value.arrval)->cols &&
                (yystack_[2].value.arrval)->rows == (yystack_[0].value.arrval)->rows) {
              (yylhs.value.arrval) =
                  array_add((yystack_[2].value.arrval), (yystack_[0].value.arrval));
            }
            else {
              yyerror(aprepro, "Arrays do not have same row and column count");
              yyerrok;
            }
          }
#line 894 "apr_parser.cc" // lalr1.cc:919
          break;

          case 39:
#line 243 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), -1.0);
          }
#line 900 "apr_parser.cc" // lalr1.cc:919
          break;

          case 40:
#line 245 "aprepro.yy" // lalr1.cc:919
          {
            if ((yystack_[2].value.arrval)->cols == (yystack_[0].value.arrval)->cols &&
                (yystack_[2].value.arrval)->rows == (yystack_[0].value.arrval)->rows) {
              (yylhs.value.arrval) =
                  array_sub((yystack_[2].value.arrval), (yystack_[0].value.arrval));
            }
            else {
              yyerror(aprepro, "Arrays do not have same row and column count");
              yyerrok;
            }
          }
#line 913 "apr_parser.cc" // lalr1.cc:919
          break;

          case 41:
#line 253 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) = array_scale((yystack_[2].value.arrval), (yystack_[0].value.val));
          }
#line 919 "apr_parser.cc" // lalr1.cc:919
          break;

          case 42:
#line 254 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) =
                array_scale((yystack_[2].value.arrval), 1.0 / (yystack_[0].value.val));
          }
#line 925 "apr_parser.cc" // lalr1.cc:919
          break;

          case 43:
#line 255 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), (yystack_[2].value.val));
          }
#line 931 "apr_parser.cc" // lalr1.cc:919
          break;

          case 44:
#line 256 "aprepro.yy" // lalr1.cc:919
          {
            if ((yystack_[2].value.arrval)->cols == (yystack_[0].value.arrval)->rows) {
              (yylhs.value.arrval) =
                  array_mult((yystack_[2].value.arrval), (yystack_[0].value.arrval));
            }
            else {
              yyerror(aprepro,
                      "Column count of first array does not match row count of second array");
              yyerrok;
            }
          }
#line 944 "apr_parser.cc" // lalr1.cc:919
          break;

          case 45:
#line 265 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string) = (yystack_[0].value.string);
          }
#line 950 "apr_parser.cc" // lalr1.cc:919
          break;

          case 46:
#line 266 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          }
#line 956 "apr_parser.cc" // lalr1.cc:919
          break;

          case 47:
#line 267 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          }
#line 962 "apr_parser.cc" // lalr1.cc:919
          break;

          case 48:
#line 268 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            set_type(aprepro, (yystack_[2].value.tptr), Parser::token::SVAR);
          }
#line 969 "apr_parser.cc" // lalr1.cc:919
          break;

          case 49:
#line 270 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          }
#line 977 "apr_parser.cc" // lalr1.cc:919
          break;

          case 50:
#line 273 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);
          }
#line 986 "apr_parser.cc" // lalr1.cc:919
          break;

          case 51:
#line 277 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string) = (char *)(yystack_[2].value.tptr)->value.svar.c_str();
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 992 "apr_parser.cc" // lalr1.cc:919
          break;

          case 52:
#line 278 "aprepro.yy" // lalr1.cc:919
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          }
#line 998 "apr_parser.cc" // lalr1.cc:919
          break;

          case 53:
#line 279 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_c == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_c))(
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1009 "apr_parser.cc" // lalr1.cc:919
          break;

          case 54:
#line 285 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.strfnct == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[2].value.tptr)->value.strfnct))();
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1020 "apr_parser.cc" // lalr1.cc:919
          break;

          case 55:
#line 291 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_d == NULL))
              (yylhs.value.string) =
                  (char *)(*((yystack_[3].value.tptr)->value.strfnct_d))((yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1031 "apr_parser.cc" // lalr1.cc:919
          break;

          case 56:
#line 297 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_a == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_a))(
                  (yystack_[1].value.arrval));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1042 "apr_parser.cc" // lalr1.cc:919
          break;

          case 57:
#line 303 "aprepro.yy" // lalr1.cc:919
          {
            concat_string((yystack_[2].value.string), (yystack_[0].value.string),
                          &(yylhs.value.string));
          }
#line 1048 "apr_parser.cc" // lalr1.cc:919
          break;

          case 58:
#line 304 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_dd == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1059 "apr_parser.cc" // lalr1.cc:919
          break;

          case 59:
#line 310 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.strfnct_dcccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[11].value.tptr)->value.strfnct_dcccc))(
                  (yystack_[9].value.val), (yystack_[7].value.string), (yystack_[5].value.string),
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1070 "apr_parser.cc" // lalr1.cc:919
          break;

          case 60:
#line 316 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_dcc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_dcc))(
                  (yystack_[5].value.val), (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1081 "apr_parser.cc" // lalr1.cc:919
          break;

          case 61:
#line 322 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_ccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1092 "apr_parser.cc" // lalr1.cc:919
          break;

          case 62:
#line 328 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_cc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1103 "apr_parser.cc" // lalr1.cc:919
          break;

          case 63:
#line 334 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string) = ((yystack_[4].value.val)) ? ((yystack_[2].value.string))
                                                             : ((yystack_[0].value.string));
          }
#line 1109 "apr_parser.cc" // lalr1.cc:919
          break;

          case 64:
#line 336 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          }
#line 1115 "apr_parser.cc" // lalr1.cc:919
          break;

          case 65:
#line 337 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.val) + 1;
          }
#line 1121 "apr_parser.cc" // lalr1.cc:919
          break;

          case 66:
#line 338 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.val) - 1;
          }
#line 1127 "apr_parser.cc" // lalr1.cc:919
          break;

          case 67:
#line 339 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          }
#line 1133 "apr_parser.cc" // lalr1.cc:919
          break;

          case 68:
#line 340 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          }
#line 1139 "apr_parser.cc" // lalr1.cc:919
          break;

          case 69:
#line 341 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
          }
#line 1145 "apr_parser.cc" // lalr1.cc:919
          break;

          case 70:
#line 342 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
          }
#line 1151 "apr_parser.cc" // lalr1.cc:919
          break;

          case 71:
#line 343 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
          }
#line 1157 "apr_parser.cc" // lalr1.cc:919
          break;

          case 72:
#line 344 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
          }
#line 1163 "apr_parser.cc" // lalr1.cc:919
          break;

          case 73:
#line 345 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          }
#line 1170 "apr_parser.cc" // lalr1.cc:919
          break;

          case 74:
#line 347 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1178 "apr_parser.cc" // lalr1.cc:919
          break;

          case 75:
#line 350 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1184 "apr_parser.cc" // lalr1.cc:919
          break;

          case 76:
#line 351 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1190 "apr_parser.cc" // lalr1.cc:919
          break;

          case 77:
#line 352 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1196 "apr_parser.cc" // lalr1.cc:919
          break;

          case 78:
#line 353 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1202 "apr_parser.cc" // lalr1.cc:919
          break;

          case 79:
#line 354 "aprepro.yy" // lalr1.cc:919
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            SEAMS::math_error(aprepro, "Power");
          }
#line 1212 "apr_parser.cc" // lalr1.cc:919
          break;

          case 80:
#line 359 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          }
#line 1218 "apr_parser.cc" // lalr1.cc:919
          break;

          case 81:
#line 360 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          }
#line 1224 "apr_parser.cc" // lalr1.cc:919
          break;

          case 82:
#line 361 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          }
#line 1230 "apr_parser.cc" // lalr1.cc:919
          break;

          case 83:
#line 362 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          }
#line 1236 "apr_parser.cc" // lalr1.cc:919
          break;

          case 84:
#line 363 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1242 "apr_parser.cc" // lalr1.cc:919
          break;

          case 85:
#line 364 "aprepro.yy" // lalr1.cc:919
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          }
#line 1248 "apr_parser.cc" // lalr1.cc:919
          break;

          case 86:
#line 365 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1254 "apr_parser.cc" // lalr1.cc:919
          break;

          case 87:
#line 366 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1260 "apr_parser.cc" // lalr1.cc:919
          break;

          case 88:
#line 367 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1266 "apr_parser.cc" // lalr1.cc:919
          break;

          case 89:
#line 368 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1272 "apr_parser.cc" // lalr1.cc:919
          break;

          case 90:
#line 369 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1278 "apr_parser.cc" // lalr1.cc:919
          break;

          case 91:
#line 371 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1285 "apr_parser.cc" // lalr1.cc:919
          break;

          case 92:
#line 373 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1293 "apr_parser.cc" // lalr1.cc:919
          break;

          case 93:
#line 376 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1301 "apr_parser.cc" // lalr1.cc:919
          break;

          case 94:
#line 379 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          }
#line 1309 "apr_parser.cc" // lalr1.cc:919
          break;

          case 95:
#line 382 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          }
#line 1317 "apr_parser.cc" // lalr1.cc:919
          break;

          case 96:
#line 385 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1324 "apr_parser.cc" // lalr1.cc:919
          break;

          case 97:
#line 387 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1332 "apr_parser.cc" // lalr1.cc:919
          break;

          case 98:
#line 390 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1340 "apr_parser.cc" // lalr1.cc:919
          break;

          case 99:
#line 393 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1348 "apr_parser.cc" // lalr1.cc:919
          break;

          case 100:
#line 396 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1356 "apr_parser.cc" // lalr1.cc:919
          break;

          case 101:
#line 399 "aprepro.yy" // lalr1.cc:919
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            SEAMS::math_error(aprepro, "Power");
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1367 "apr_parser.cc" // lalr1.cc:919
          break;

          case 102:
#line 406 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.fnctptr == NULL))
              (yylhs.value.val) = (*((yystack_[2].value.tptr)->value.fnctptr))();
            else
              (yylhs.value.val) = 0.0;
          }
#line 1378 "apr_parser.cc" // lalr1.cc:919
          break;

          case 103:
#line 413 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_d == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_d))((yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1389 "apr_parser.cc" // lalr1.cc:919
          break;

          case 104:
#line 420 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_c == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_c))((yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1400 "apr_parser.cc" // lalr1.cc:919
          break;

          case 105:
#line 427 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_a == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_a))((yystack_[1].value.arrval));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1411 "apr_parser.cc" // lalr1.cc:919
          break;

          case 106:
#line 434 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1422 "apr_parser.cc" // lalr1.cc:919
          break;

          case 107:
#line 441 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dc))(
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1433 "apr_parser.cc" // lalr1.cc:919
          break;

          case 108:
#line 448 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1444 "apr_parser.cc" // lalr1.cc:919
          break;

          case 109:
#line 455 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccc == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              yyerrok;
          }
#line 1455 "apr_parser.cc" // lalr1.cc:919
          break;

          case 110:
#line 462 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1466 "apr_parser.cc" // lalr1.cc:919
          break;

          case 111:
#line 468 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ddd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1477 "apr_parser.cc" // lalr1.cc:919
          break;

          case 112:
#line 474 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccd))(
                  (yystack_[5].value.string), (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1488 "apr_parser.cc" // lalr1.cc:919
          break;

          case 113:
#line 480 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1499 "apr_parser.cc" // lalr1.cc:919
          break;

          case 114:
#line 486 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1510 "apr_parser.cc" // lalr1.cc:919
          break;

          case 115:
#line 492 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.fnctptr_ddddc == NULL))
              (yylhs.value.val) = (*((yystack_[11].value.tptr)->value.fnctptr_ddddc))(
                  (yystack_[9].value.val), (yystack_[7].value.val), (yystack_[5].value.val),
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1521 "apr_parser.cc" // lalr1.cc:919
          break;

          case 116:
#line 498 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[13].value.tptr),
                          (yystack_[13].value.tptr)->value.fnctptr_dddddd == NULL))
              (yylhs.value.val) = (*((yystack_[13].value.tptr)->value.fnctptr_dddddd))(
                  (yystack_[11].value.val), (yystack_[9].value.val), (yystack_[7].value.val),
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1532 "apr_parser.cc" // lalr1.cc:919
          break;

          case 117:
#line 504 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) + (yystack_[0].value.val);
          }
#line 1538 "apr_parser.cc" // lalr1.cc:919
          break;

          case 118:
#line 505 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) - (yystack_[0].value.val);
          }
#line 1544 "apr_parser.cc" // lalr1.cc:919
          break;

          case 119:
#line 506 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) * (yystack_[0].value.val);
          }
#line 1550 "apr_parser.cc" // lalr1.cc:919
          break;

          case 120:
#line 507 "aprepro.yy" // lalr1.cc:919
          {
            if ((yystack_[0].value.val) == 0.) {
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (yystack_[2].value.val) / (yystack_[0].value.val);
          }
#line 1562 "apr_parser.cc" // lalr1.cc:919
          break;

          case 121:
#line 514 "aprepro.yy" // lalr1.cc:919
          {
            if ((yystack_[0].value.val) == 0.) {
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (int)(yystack_[2].value.val) % (int)(yystack_[0].value.val);
          }
#line 1574 "apr_parser.cc" // lalr1.cc:919
          break;

          case 122:
#line 521 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = -(yystack_[0].value.val);
          }
#line 1580 "apr_parser.cc" // lalr1.cc:919
          break;

          case 123:
#line 522 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          }
#line 1586 "apr_parser.cc" // lalr1.cc:919
          break;

          case 124:
#line 523 "aprepro.yy" // lalr1.cc:919
          {
            reset_error();
            (yylhs.value.val) = std::pow((yystack_[2].value.val), (yystack_[0].value.val));
            SEAMS::math_error(aprepro, "Power");
          }
#line 1594 "apr_parser.cc" // lalr1.cc:919
          break;

          case 125:
#line 526 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          }
#line 1600 "apr_parser.cc" // lalr1.cc:919
          break;

          case 126:
#line 527 "aprepro.yy" // lalr1.cc:919
          {
            reset_error();
            (yylhs.value.val) =
                (double)((yystack_[1].value.val) < 0 ? -floor(-((yystack_[1].value.val)))
                                                     : floor((yystack_[1].value.val)));
            SEAMS::math_error(aprepro, "floor (int)");
          }
#line 1608 "apr_parser.cc" // lalr1.cc:919
          break;

          case 127:
#line 530 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ((yystack_[0].value.val)) ? 1 : 0;
          }
#line 1614 "apr_parser.cc" // lalr1.cc:919
          break;

          case 128:
#line 531 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                ((yystack_[4].value.val)) ? ((yystack_[2].value.val)) : ((yystack_[0].value.val));
          }
#line 1620 "apr_parser.cc" // lalr1.cc:919
          break;

          case 129:
#line 532 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                array_value((yystack_[3].value.tptr)->value.avar, (yystack_[1].value.val), 0);
          }
#line 1626 "apr_parser.cc" // lalr1.cc:919
          break;

          case 130:
#line 533 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = array_value((yystack_[5].value.tptr)->value.avar,
                                            (yystack_[3].value.val), (yystack_[1].value.val));
          }
#line 1632 "apr_parser.cc" // lalr1.cc:919
          break;

          case 131:
#line 535 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.val);
            array *arr        = (yystack_[5].value.tptr)->value.avar;
            int    cols       = arr->cols;
            if (cols > 1) {
              yyerror(aprepro, "Cannot use [index] array access with multi-column array");
              yyerrok;
            }
            int rows = arr->rows;
            int row  = (yystack_[3].value.val);
            if (aprepro.ap_options.one_based_index) {
              row--;
            }
            if (row < rows) {
              int offset                                         = row * cols;
              (yystack_[5].value.tptr)->value.avar->data[offset] = (yystack_[0].value.val);
            }
            else {
              yyerror(aprepro, "Row or Column index out of range");
              yyerrok;
            }
          }
#line 1658 "apr_parser.cc" // lalr1.cc:919
          break;

          case 132:
#line 557 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.val);
            array *arr        = (yystack_[7].value.tptr)->value.avar;
            int    cols       = arr->cols;
            int    rows       = arr->rows;
            int    row        = (yystack_[5].value.val);
            int    col        = (yystack_[3].value.val);
            if (aprepro.ap_options.one_based_index) {
              row--;
              col--;
            }
            if (row < rows && col < cols) {
              int offset                                         = row * cols + col;
              (yystack_[7].value.tptr)->value.avar->data[offset] = (yystack_[0].value.val);
            }
            else {
              yyerror(aprepro, "Row or Column index out of range");
              yyerrok;
            }
          }
#line 1682 "apr_parser.cc" // lalr1.cc:919
          break;

#line 1686 "apr_parser.cc" // lalr1.cc:919
          default: break;
          }
        }
#if YY_EXCEPTIONS
        catch (const syntax_error &yyexc) {
          YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
          error(yyexc);
          YYERROR;
        }
#endif // YY_EXCEPTIONS
        YY_SYMBOL_PRINT("-> $$ =", yylhs);
        yypop_(yylen);
        yylen = 0;
        YY_STACK_PRINT();

        // Shift the result of the reduction.
        yypush_(YY_NULLPTR, YY_MOVE(yylhs));
      }
      goto yynewstate;

    /*--------------------------------------.
    | yyerrlab -- here on detecting error.  |
    `--------------------------------------*/
    yyerrlab:
      // If not already recovering from an error, report this error.
      if (!yyerrstatus_) {
        ++yynerrs_;
        error(yysyntax_error_(yystack_[0].state, yyla));
      }

      if (yyerrstatus_ == 3) {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.type_get() == yyeof_)
          YYABORT;
        else if (!yyla.empty()) {
          yy_destroy_("Error: discarding", yyla);
          yyla.clear();
        }
      }

      // Else will try to reuse lookahead token after shifting the error token.
      goto yyerrlab1;

    /*---------------------------------------------------.
    | yyerrorlab -- error raised explicitly by YYERROR.  |
    `---------------------------------------------------*/
    yyerrorlab:
      /* Pacify compilers when the user code never invokes YYERROR and
         the label yyerrorlab therefore never appears in user code.  */
      if (false)
        YYERROR;

      /* Do not reclaim the symbols of the rule whose action triggered
         this YYERROR.  */
      yypop_(yylen);
      yylen = 0;
      goto yyerrlab1;

    /*-------------------------------------------------------------.
    | yyerrlab1 -- common code for both syntax error and YYERROR.  |
    `-------------------------------------------------------------*/
    yyerrlab1:
      yyerrstatus_ = 3; // Each real token shifted decrements this.
      {
        stack_symbol_type error_token;
        for (;;) {
          yyn = yypact_[yystack_[0].state];
          if (!yy_pact_value_is_default_(yyn)) {
            yyn += yyterror_;
            if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == yyterror_) {
              yyn = yytable_[yyn];
              if (0 < yyn)
                break;
            }
          }

          // Pop the current state because it cannot handle the error token.
          if (yystack_.size() == 1)
            YYABORT;

          yy_destroy_("Error: popping", yystack_[0]);
          yypop_();
          YY_STACK_PRINT();
        }

        // Shift the error token.
        error_token.state = yyn;
        yypush_("Shifting", YY_MOVE(error_token));
      }
      goto yynewstate;

    /*-------------------------------------.
    | yyacceptlab -- YYACCEPT comes here.  |
    `-------------------------------------*/
    yyacceptlab:
      yyresult = 0;
      goto yyreturn;

    /*-----------------------------------.
    | yyabortlab -- YYABORT comes here.  |
    `-----------------------------------*/
    yyabortlab:
      yyresult = 1;
      goto yyreturn;

    /*-----------------------------------------------------.
    | yyreturn -- parsing is finished, return the result.  |
    `-----------------------------------------------------*/
    yyreturn:
      if (!yyla.empty())
        yy_destroy_("Cleanup: discarding lookahead", yyla);

      /* Do not reclaim the symbols of the rule whose action triggered
         this YYABORT or YYACCEPT.  */
      yypop_(yylen);
      while (1 < yystack_.size()) {
        yy_destroy_("Cleanup: popping", yystack_[0]);
        yypop_();
      }

      return yyresult;
    }
#if YY_EXCEPTIONS
    catch (...) {
      YYCDEBUG << "Exception caught: cleaning lookahead and stack\n";
      // Do not try to display the values of the reclaimed symbols,
      // as their printers might throw an exception.
      if (!yyla.empty())
        yy_destroy_(YY_NULLPTR, yyla);

      while (1 < yystack_.size()) {
        yy_destroy_(YY_NULLPTR, yystack_[0]);
        yypop_();
      }
      throw;
    }
#endif // YY_EXCEPTIONS
  }

  void Parser::error(const syntax_error &yyexc) { error(yyexc.what()); }

  // Generate an error message.
  std::string Parser::yysyntax_error_(state_type yystate, const symbol_type &yyla) const
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
    if (!yyla.empty()) {
      int yytoken      = yyla.type_get();
      yyarg[yycount++] = yytname_[yytoken];
      int yyn          = yypact_[yystate];
      if (!yy_pact_value_is_default_(yyn)) {
        /* Start YYX at -YYN if negative to avoid negative indexes in
           YYCHECK.  In other words, skip the first -YYN actions for
           this state because they are default actions.  */
        int yyxbegin = yyn < 0 ? -yyn : 0;
        // Stay within bounds of both yycheck and yytname.
        int yychecklim = yylast_ - yyn + 1;
        int yyxend     = yychecklim < yyntokens_ ? yychecklim : yyntokens_;
        for (int yyx = yyxbegin; yyx < yyxend; ++yyx)
          if (yycheck_[yyx + yyn] == yyx && yyx != yyterror_ &&
              !yy_table_value_is_error_(yytable_[yyx + yyn])) {
            if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM) {
              yycount = 1;
              break;
            }
            else
              yyarg[yycount++] = yytname_[yyx];
          }
      }
    }

    char const *yyformat = YY_NULLPTR;
    switch (yycount) {
#define YYCASE_(N, S)                                                                              \
  case N: yyformat = S; break
    default: // Avoid compiler warnings.
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
    for (char const *yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount) {
        yyres += yytnamerr_(yyarg[yyi++]);
        ++yyp;
      }
      else
        yyres += *yyp;
    return yyres;
  }

  const signed char Parser::yypact_ninf_ = -34;

  const signed char Parser::yytable_ninf_ = -1;

  const short Parser::yypact_[] = {
      -34,  2,    -34,  -3,   289,  -34,  -34,  -34,  -34,  -34,  -13,  296,  4,    600,  28,
      45,   44,   46,   53,   405,  405,  -34,  405,  360,  405,  111,  179,  48,   -16,  38,
      1126, 360,  405,  405,  405,  405,  405,  -34,  -34,  405,  405,  405,  405,  405,  405,
      -34,  -34,  405,  405,  405,  405,  405,  405,  405,  -34,  -34,  405,  405,  360,  209,
      345,  360,  630,  52,   405,  93,   1174, 864,  1078, 18,   -34,  18,   18,   -34,  -34,
      -34,  -34,  -34,  -34,  -34,  -34,  405,  405,  405,  -34,  360,  360,  405,  360,  -34,
      405,  405,  405,  405,  405,  405,  405,  -34,  405,  405,  405,  405,  405,  405,  405,
      405,  405,  405,  405,  360,  405,  405,  -33,  1174, 1193, 1209, 1209, 1209, 1209, 1209,
      1174, 1209, 1209, 1209, 1209, 1209, 1209, 1174, 1209, 1174, 1209, 1209, 1209, 1209, 1209,
      1209, 1174, 1209, 599,  -33,  1193, -34,  50,   99,   629,  -34,  154,  229,  659,  230,
      421,  689,  405,  18,   -34,  -34,  405,  -34,  1140, 1160, 49,   1209, -34,  1,    1,
      1225, -34,  1225, 39,   39,   39,   39,   39,   39,   -34,  1240, 1254, 340,  340,  340,
      340,  340,  340,  190,  190,  18,   -34,  18,   18,   18,   405,  70,   -34,  405,  -34,
      405,  -34,  -34,  405,  -34,  405,  -34,  -34,  405,  -34,  405,  -34,  1209, 18,   405,
      405,  1103, 405,  431,  891,  484,  570,  456,  131,  918,  490,  945,  719,  1174, 1209,
      71,   1209, 405,  -34,  -34,  -34,  405,  -34,  405,  405,  -34,  405,  -34,  -34,  -34,
      405,  -34,  405,  512,  972,  749,  808,  518,  462,  999,  1209, -34,  -34,  405,  -34,
      405,  -34,  405,  -34,  -34,  779,  1026, 391,  405,  -34,  -34,  405,  540,  837,  546,
      -34,  405,  -34,  1053, -34};

  const unsigned char Parser::yydefact_[] = {
      2,   0,   1,   0,   0,   4,   3,   9,   64,  45, 91, 67, 46,  68,  47,  28,  0,   0,   0,
      0,   0,   8,   0,   0,   0,   0,   0,   127, 0,  0,  0,  0,   0,   0,   0,   0,   0,   94,
      95,  0,   0,   0,   0,   0,   0,   71,  72,  0,  0,  0,  0,   0,   0,   0,   82,  83,  0,
      0,   0,   0,   0,   0,   91,  0,   0,   127, 0,  0,  0,  123, 39,  122, 12,  65,  92,  69,
      80,  66,  93,  70,  81,  0,   0,   0,   7,   0,  0,  0,  0,   6,   0,   0,   0,   0,   0,
      0,   0,   5,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,   0,   0,   0,   0,   37,  48,
      96,  97,  98,  99,  100, 101, 50,  73,  75,  76, 77, 78, 79,  49,  74,  52,  84,  86,  87,
      88,  89,  90,  51,  85,  0,   36,  0,   102, 0,  0,  0,  54,  0,   0,   0,   0,   0,   0,
      0,   122, 21,  125, 0,   126, 0,   0,   19,  0,  20, 38, 40,  42,  44,  41,  22,  23,  24,
      25,  26,  27,  57,  17,  18,  10,  11,  13,  14, 15, 16, 117, 118, 120, 43,  119, 121, 124,
      0,   129, 105, 0,   104, 0,   103, 56,  0,   53, 0,  55, 35,  0,   29,  0,   34,  96,  119,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,   0,   63,  128, 130, 131, 0,
      108, 106, 107, 0,   110, 0,   0,   62,  0,   58, 31, 30, 0,   33,  0,   0,   0,   0,   0,
      0,   0,   0,   132, 109, 112, 0,   111, 0,   61, 0,  60, 32,  0,   0,   0,   0,   114, 113,
      0,   0,   0,   0,   115, 0,   59,  0,   116};

  const signed char Parser::yypgoto_[] = {-34, -34, -34, -18, 95, 81, -4};

  const short Parser::yydefgoto_[] = {-1, 1, 6, 27, 28, 66, 161};

  const unsigned short Parser::yytable_[] = {
      30,  65,  2,   3,   84,  85,  86,  87,  88,  31,  32,  33,  34,  35,  36,  67,  68,  7,   69,
      71,  72,  4,   85,  86,  87,  88,  47,  114, 115, 116, 117, 118, 119, 37,  38,  121, 122, 123,
      124, 125, 126, 87,  88,  128, 130, 131, 132, 133, 134, 135, 56,  5,   137, 138, 140, 144, 148,
      151, 89,  59,  153, 60,  57,  111, 160, 162, 192, 58,  61,  57,  90,  91,  92,  93,  94,  95,
      81,  159, 82,  83,  83,  140, 140, 165, 167, 29,  96,  96,  85,  86,  87,  88,  212, 242, 175,
      176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 187, 188, 189, 0,   154, 0,   0,   113, 193,
      73,  194, 74,  75,  70,  76,  120, 81,  0,   82,  83,  0,   112, 0,   127, 129, 0,   90,  91,
      92,  93,  94,  95,  136, 0,   0,   143, 147, 150, 0,   0,   236, 0,   96,  207, 0,   0,   0,
      208, 139, 142, 146, 149, 0,   0,   0,   0,   0,   158, 90,  91,  92,  93,  94,  95,  0,   197,
      168, 169, 170, 171, 172, 173, 174, 0,   96,  163, 164, 77,  166, 78,  79,  211, 80,  0,   214,
      0,   216, 85,  86,  87,  88,  219, 0,   0,   221, 0,   222, 0,   0,   186, 0,   224, 0,   226,
      0,   0,   0,   8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  244, 19,  141, 20,  245,
      0,   246, 108, 156, 110, 113, 0,   111, 249, 0,   250, 0,   0,   0,   0,   198, 0,   199, 202,
      22,  23,  260, 0,   261, 0,   24,  0,   25,  26,  0,   0,   268, 0,   90,  91,  92,  93,  94,
      95,  273, 85,  86,  87,  88,  0,   0,   213, 0,   215, 96,  0,   217, 0,   218, 0,   0,   220,
      0,   0,   0,   0,   0,   223, 0,   8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,
      19,  0,   20,  0,   243, 21,  0,   0,   0,   0,   0,   247, 0,   248, 39,  40,  41,  42,  43,
      44,  0,   0,   0,   22,  23,  0,   0,   0,   0,   24,  0,   25,  26,  0,   262, 0,   0,   0,
      45,  46,  267, 0,   0,   269, 8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,
      145, 20,  8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  0,   20,  106, 107,
      108, 156, 110, 22,  23,  111, 0,   0,   0,   24,  0,   25,  26,  0,   0,   0,   0,   0,   22,
      23,  0,   0,   0,   0,   24,  266, 25,  26,  8,   9,   62,  11,  12,  13,  14,  63,  16,  17,
      0,   0,   19,  0,   20,  90,  91,  92,  93,  94,  95,  0,   0,   0,   0,   0,   0,   203, 0,
      204, 0,   96,  0,   0,   0,   22,  64,  227, 0,   228, 0,   24,  0,   25,  26,  90,  91,  92,
      93,  94,  95,  0,   0,   0,   0,   90,  91,  92,  93,  94,  95,  96,  234, 0,   235, 0,   0,
      0,   257, 0,   258, 96,  0,   0,   0,   0,   0,   0,   0,   0,   90,  91,  92,  93,  94,  95,
      90,  91,  92,  93,  94,  95,  230, 0,   0,   0,   96,  0,   238, 0,   0,   0,   96,  0,   0,
      0,   0,   0,   90,  91,  92,  93,  94,  95,  90,  91,  92,  93,  94,  95,  251, 0,   0,   0,
      96,  0,   256, 0,   0,   0,   96,  0,   0,   0,   0,   0,   90,  91,  92,  93,  94,  95,  90,
      91,  92,  93,  94,  95,  270, 0,   0,   0,   96,  0,   272, 0,   0,   0,   96,  0,   0,   0,
      0,   0,   90,  91,  92,  93,  94,  95,  90,  91,  92,  93,  94,  95,  231, 0,   232, 0,   96,
      0,   0,   233, 0,   0,   96,  0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105,
      106, 107, 108, 156, 110, 190, 0,   111, 0,   191, 0,   0,   0,   0,   48,  49,  50,  51,  52,
      53,  0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 0,   195, 111, 196,
      54,  55,  0,   0,   0,   0,   152, 32,  33,  34,  35,  36,  0,   98,  99,  100, 101, 102, 103,
      104, 105, 106, 107, 108, 109, 110, 0,   200, 111, 201, 37,  38,  0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 0,
      205, 111, 206, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100,
      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 0,   240, 111, 241, 0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108,
      156, 110, 0,   253, 111, 254, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 0,   263, 111, 264, 0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105,
      106, 107, 108, 156, 110, 255, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 271, 0,   111, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104,
      105, 106, 107, 108, 156, 110, 155, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 229, 0,   111, 0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107,
      108, 156, 110, 237, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,
      100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 239, 0,   111, 0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110,
      252, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102,
      103, 104, 105, 106, 107, 108, 156, 110, 259, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 265, 0,   111,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105,
      106, 107, 108, 156, 110, 274, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 157, 0,   111, 0,   0,   0,
      0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110,
      225, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104,
      105, 106, 107, 108, 156, 110, 97,  0,   111, 0,   0,   0,   0,   0,   0,   0,   98,  99,  100,
      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 209, 0,   111, 90,  91,  92,  93,  94,  95,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   96,  210, 98,  99,  100, 101, 102, 103, 104,
      105, 106, 107, 108, 156, 110, 0,   0,   111, 90,  91,  92,  93,  94,  95,  0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   96,  98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
      110, 0,   0,   111, 98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 0,   0,
      111, 98,  99,  100, 101, 102, 103, 104, 105, 0,   0,   0,   0,   110, 0,   0,   111, 99,  100,
      101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 0,   0,   111, 100, 101, 102, 103, 104, 105,
      106, 107, 108, 156, 110, 0,   0,   111};

  const short Parser::yycheck_[] = {
      4,   19,  0,   1,   20,  38,  39,  40,  41,  22,  23,  24,  25,  26,  27,  19,  20,  20, 22,
      23,  24,  19,  38,  39,  40,  41,  22,  31,  32,  33,  34,  35,  36,  46,  47,  39,  40, 41,
      42,  43,  44,  40,  41,  47,  48,  49,  50,  51,  52,  53,  22,  49,  56,  57,  58,  59, 60,
      61,  20,  15,  64,  15,  17,  45,  82,  83,  16,  22,  15,  17,  32,  33,  34,  35,  36, 37,
      28,  81,  30,  31,  31,  85,  86,  87,  88,  4,   48,  48,  38,  39,  40,  41,  22,  22, 98,
      99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, -1,  16,  -1,  -1,  31, 14,
      3,   16,  5,   6,   23,  8,   39,  28,  -1,  30,  31,  -1,  31,  -1,  47,  48,  -1,  32, 33,
      34,  35,  36,  37,  56,  -1,  -1,  59,  60,  61,  -1,  -1,  14,  -1,  48,  152, -1,  -1, -1,
      156, 58,  59,  60,  61,  -1,  -1,  -1,  -1,  -1,  81,  32,  33,  34,  35,  36,  37,  -1, 16,
      90,  91,  92,  93,  94,  95,  96,  -1,  48,  85,  86,  3,   88,  5,   6,   190, 8,   -1, 193,
      -1,  195, 38,  39,  40,  41,  200, -1,  -1,  203, -1,  205, -1,  -1,  109, -1,  210, -1, 212,
      -1,  -1,  -1,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  227, 15,  16,  17, 231,
      -1,  233, 40,  41,  42,  152, -1,  45,  240, -1,  242, -1,  -1,  -1,  -1,  14,  -1,  16, 16,
      38,  39,  253, -1,  255, -1,  44,  -1,  46,  47,  -1,  -1,  263, -1,  32,  33,  34,  35, 36,
      37,  271, 38,  39,  40,  41,  -1,  -1,  193, -1,  195, 48,  -1,  198, -1,  200, -1,  -1, 203,
      -1,  -1,  -1,  -1,  -1,  209, -1,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13, -1,
      15,  -1,  17,  -1,  227, 20,  -1,  -1,  -1,  -1,  -1,  234, -1,  236, 22,  23,  24,  25, 26,
      27,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  257, -1,  -1, -1,
      46,  47,  263, -1,  -1,  266, 3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1, 15,
      16,  17,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,  -1,  17,  38, 39,
      40,  41,  42,  38,  39,  45,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,  -1,  -1, 38,
      39,  -1,  -1,  -1,  -1,  44,  14,  46,  47,  3,   4,   5,   6,   7,   8,   9,   10,  11, 12,
      -1,  -1,  15,  -1,  17,  32,  33,  34,  35,  36,  37,  -1,  -1,  -1,  -1,  -1,  -1,  14, -1,
      16,  -1,  48,  -1,  -1,  -1,  38,  39,  14,  -1,  16,  -1,  44,  -1,  46,  47,  32,  33, 34,
      35,  36,  37,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  48,  14,  -1,  16,  -1, -1,
      -1,  14,  -1,  16,  48,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36, 37,
      32,  33,  34,  35,  36,  37,  16,  -1,  -1,  -1,  48,  -1,  16,  -1,  -1,  -1,  48,  -1, -1,
      -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  32,  33,  34,  35,  36,  37,  16,  -1,  -1, -1,
      48,  -1,  16,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37, 32,
      33,  34,  35,  36,  37,  16,  -1,  -1,  -1,  48,  -1,  16,  -1,  -1,  -1,  48,  -1,  -1, -1,
      -1,  -1,  32,  33,  34,  35,  36,  37,  32,  33,  34,  35,  36,  37,  14,  -1,  16,  -1, 48,
      -1,  -1,  21,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36, 37,
      38,  39,  40,  41,  42,  14,  -1,  45,  -1,  18,  -1,  -1,  -1,  -1,  22,  23,  24,  25, 26,
      27,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45, 16,
      46,  47,  -1,  -1,  -1,  -1,  22,  23,  24,  25,  26,  27,  -1,  30,  31,  32,  33,  34, 35,
      36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  46,  47,  -1,  -1,  -1,  -1,  -1, -1,
      -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42, -1,
      14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31, 32,
      33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1, -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39, 40,
      41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1,
      30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1, -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36, 37,
      38,  39,  40,  41,  42,  14,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1,
      -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  14,  -1,  45, -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35, 36,
      37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1,
      -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1, -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38, 39,
      40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30, 31,
      32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1, -1,
      -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41, 42,
      16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33, 34,
      35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1,
      -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1, 45,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36, 37,
      38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1,
      30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  18,  -1,  45,  -1,  -1, -1,
      -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41, 42,
      18,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35, 36,
      37,  38,  39,  40,  41,  42,  20,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31, 32,
      33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  29,  -1,  45,  32,  33,  34,  35,  36, 37,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  48,  29,  30,  31,  32,  33,  34,  35, 36,
      37,  38,  39,  40,  41,  42,  -1,  -1,  45,  32,  33,  34,  35,  36,  37,  -1,  -1,  -1, -1,
      -1,  -1,  -1,  -1,  -1,  -1,  48,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40, 41,
      42,  -1,  -1,  45,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1, -1,
      45,  30,  31,  32,  33,  34,  35,  36,  37,  -1,  -1,  -1,  -1,  42,  -1,  -1,  45,  31, 32,
      33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  32,  33,  34,  35,  36, 37,
      38,  39,  40,  41,  42,  -1,  -1,  45};

  const unsigned char Parser::yystos_[] = {
      0,  51, 0,  1,  19, 49, 52, 20, 3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 15, 17, 20, 38,
      39, 44, 46, 47, 53, 54, 55, 56, 22, 23, 24, 25, 26, 27, 46, 47, 22, 23, 24, 25, 26, 27, 46,
      47, 22, 22, 23, 24, 25, 26, 27, 46, 47, 22, 17, 22, 15, 15, 15, 5,  10, 39, 53, 55, 56, 56,
      56, 54, 56, 56, 3,  5,  6,  8,  3,  5,  6,  8,  28, 30, 31, 20, 38, 39, 40, 41, 20, 32, 33,
      34, 35, 36, 37, 48, 20, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 45, 54, 55, 56,
      56, 56, 56, 56, 56, 55, 56, 56, 56, 56, 56, 56, 55, 56, 55, 56, 56, 56, 56, 56, 56, 55, 56,
      56, 54, 56, 16, 54, 55, 56, 16, 54, 55, 56, 54, 55, 56, 22, 56, 16, 16, 41, 18, 55, 56, 53,
      56, 53, 54, 54, 56, 54, 56, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 54, 56, 56, 56, 14, 18, 16, 14, 16, 14, 16, 16, 14, 16, 14, 16, 16, 14, 16, 14, 16,
      56, 56, 29, 29, 56, 22, 55, 56, 55, 56, 55, 55, 56, 55, 56, 56, 55, 56, 18, 56, 14, 16, 16,
      16, 14, 16, 21, 14, 16, 14, 16, 16, 16, 14, 16, 22, 55, 56, 56, 56, 55, 55, 56, 56, 16, 16,
      14, 16, 14, 16, 14, 16, 16, 56, 56, 55, 14, 16, 16, 14, 55, 56, 55, 16, 14, 16, 56, 16};

  const unsigned char Parser::yyr1_[] = {
      0,  50, 51, 51, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
      53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 55,
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56};

  const unsigned char Parser::yyr2_[] = {
      0, 2, 0, 2, 1, 3,  3,  3,  2,  2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
      3, 1, 4, 6, 6, 8,  6,  4,  4,  3, 3, 3, 2, 3, 3, 3, 3, 3, 1, 1, 1, 3, 3, 3, 3, 3, 4,
      3, 4, 4, 3, 6, 12, 8,  8,  6,  5, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2,
      2, 2, 2, 3, 3, 3,  3,  3,  3,  3, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 6, 6,
      6, 8, 6, 8, 8, 10, 10, 12, 14, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 1, 5, 4, 6, 6, 8};

  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a yyntokens_, nonterminals.
  const char *const Parser::yytname_[] = {"\"end of file\"",
                                          "error",
                                          "$undefined",
                                          "NUM",
                                          "QSTRING",
                                          "UNDVAR",
                                          "VAR",
                                          "SVAR",
                                          "IMMVAR",
                                          "IMMSVAR",
                                          "AVAR",
                                          "FNCT",
                                          "SFNCT",
                                          "AFNCT",
                                          "COMMA",
                                          "LPAR",
                                          "RPAR",
                                          "LBRACK",
                                          "RBRACK",
                                          "LBRACE",
                                          "RBRACE",
                                          "SEMI",
                                          "EQUAL",
                                          "EQ_PLUS",
                                          "EQ_MINUS",
                                          "EQ_TIME",
                                          "EQ_DIV",
                                          "EQ_POW",
                                          "QUEST",
                                          "COLON",
                                          "LOR",
                                          "LAND",
                                          "LT",
                                          "GT",
                                          "LE",
                                          "GE",
                                          "EQ",
                                          "NE",
                                          "PLU",
                                          "SUB",
                                          "DIV",
                                          "TIM",
                                          "MOD",
                                          "UNARY",
                                          "NOT",
                                          "POW",
                                          "INC",
                                          "DEC",
                                          "CONCAT",
                                          "'\\n'",
                                          "$accept",
                                          "input",
                                          "line",
                                          "bool",
                                          "aexp",
                                          "sexp",
                                          "exp",
                                          YY_NULLPTR};

#if SEAMSDEBUG
  const unsigned short Parser::yyrline_[] = {
      0,   145, 145, 146, 149, 150, 157, 161, 162, 163, 166, 167, 168, 169, 170, 171, 172,
      173, 174, 175, 176, 177, 180, 181, 182, 183, 184, 185, 187, 188, 194, 200, 206, 212,
      218, 224, 230, 233, 235, 243, 245, 253, 254, 255, 256, 265, 266, 267, 268, 270, 273,
      277, 278, 279, 285, 291, 297, 303, 304, 310, 316, 322, 328, 334, 336, 337, 338, 339,
      340, 341, 342, 343, 344, 345, 347, 350, 351, 352, 353, 354, 359, 360, 361, 362, 363,
      364, 365, 366, 367, 368, 369, 371, 373, 376, 379, 382, 385, 387, 390, 393, 396, 399,
      406, 413, 420, 427, 434, 441, 448, 455, 462, 468, 474, 480, 486, 492, 498, 504, 505,
      506, 507, 514, 521, 522, 523, 526, 527, 530, 531, 532, 533, 534, 556};

  // Print the state stack on the debug stream.
  void Parser::yystack_print_()
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator i = yystack_.begin(), i_end = yystack_.end(); i != i_end; ++i)
      *yycdebug_ << ' ' << i->state;
    *yycdebug_ << '\n';
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void Parser::yy_reduce_print_(int yyrule)
  {
    unsigned yylno  = yyrline_[yyrule];
    int      yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1 << " (line " << yylno << "):\n";
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT("   $" << yyi + 1 << " =", yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // SEAMSDEBUG

  Parser::token_number_type Parser::yytranslate_(int t)
  {
    // YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to
    // TOKEN-NUM as returned by yylex.
    static const token_number_type translate_table[] = {
        0,  2,  2,  2,  2,  2,  2,  2,  2,  2,  49, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  1,  2,  3,  4,  5,  6,  7,  8,
        9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
        31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48};
    const unsigned          user_token_number_max_ = 303;
    const token_number_type undef_token_           = 2;

    if (static_cast<int>(t) <= yyeof_)
      return yyeof_;
    else if (static_cast<unsigned>(t) <= user_token_number_max_)
      return translate_table[t];
    else
      return undef_token_;
  }

} // namespace SEAMS
#line 2494 "apr_parser.cc" // lalr1.cc:1242
#line 579 "aprepro.yy"     // lalr1.cc:1243

void SEAMS::Parser::error(const std::string &m) { aprepro.error(m); }
