// A Bison parser, made by GNU Bison 3.8.2.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2021, 2023 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

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

// DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
// especially those whose name start with YY_ or yy_.  They are
// private implementation details that can be changed or removed.

// Take the name prefix into account.
#define yylex SEAMSlex

// First part of user prologue.
#line 6 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"

#include "apr_array.h"
#include "apr_util.h"
#include "aprepro.h"

#if defined FMT_SUPPORT
#include <fmt/format.h>
#endif
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace {
  void reset_error()
  {
#if !defined(WIN32) && !defined(__WIN32__) && !defined(_WIN32) && !defined(_MSC_VER) &&            \
    !defined(__MINGW32__) && !defined(_WIN64) && !defined(__MINGW64__)
#ifndef math_errhandling
#define math_errhandling MATH_ERRNO
#endif

    if (math_errhandling & MATH_ERREXCEPT) {
      std::feclearexcept(FE_ALL_EXCEPT);
    }
    if (math_errhandling & MATH_ERRNO) {
      errno = 0;
    }
#endif
  }
} // namespace

namespace SEAMS {
  extern bool echo;
}

#line 82 "apr_parser.cc"

#include "aprepro_parser.h"

// Second part of user prologue.
#line 110 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"

#include "apr_scanner.h"
#include "aprepro.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex

#line 101 "apr_parser.cc"

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
      yy_stack_print_();                                                                           \
  } while (false)

#else // !SEAMSDEBUG

#define YYCDEBUG                                                                                   \
  if (false)                                                                                       \
  std::cerr
#define YY_SYMBOL_PRINT(Title, Symbol) YY_USE(Symbol)
#define YY_REDUCE_PRINT(Rule)          static_cast<void>(0)
#define YY_STACK_PRINT()               static_cast<void>(0)

#endif // !SEAMSDEBUG

#define yyerrok   (yyerrstatus_ = 0)
#define yyclearin (yyla.clear())

#define YYACCEPT       goto yyacceptlab
#define YYABORT        goto yyabortlab
#define YYERROR        goto yyerrorlab
#define YYRECOVERING() (!!yyerrstatus_)

namespace SEAMS {
#line 175 "apr_parser.cc"

  /// Build a parser object.
  Parser::Parser(class Aprepro &aprepro_yyarg)
#if SEAMSDEBUG
      : yydebug_(false), yycdebug_(&std::cerr),
#else
      :
#endif
        aprepro(aprepro_yyarg)
  {
  }

  Parser::~Parser() = default;

  Parser::syntax_error::~syntax_error() YY_NOEXCEPT YY_NOTHROW {}

  /*---------.
  | symbol.  |
  `---------*/

  // basic_symbol.
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
  Parser::basic_symbol<Base>::basic_symbol(typename Base::kind_type t, YY_RVREF(value_type) v)
      : Base(t), value(YY_MOVE(v))
  {
  }

  template <typename Base>
  Parser::symbol_kind_type Parser::basic_symbol<Base>::type_get() const YY_NOEXCEPT
  {
    return this->kind();
  }

  template <typename Base> bool Parser::basic_symbol<Base>::empty() const YY_NOEXCEPT
  {
    return this->kind() == symbol_kind::S_YYEMPTY;
  }

  template <typename Base> void Parser::basic_symbol<Base>::move(basic_symbol &s)
  {
    super_type::move(s);
    value = YY_MOVE(s.value);
  }

  // by_kind.
  Parser::by_kind::by_kind() YY_NOEXCEPT : kind_(symbol_kind::S_YYEMPTY) {}

#if 201103L <= YY_CPLUSPLUS
  Parser::by_kind::by_kind(by_kind &&that) YY_NOEXCEPT : kind_(that.kind_) { that.clear(); }
#endif

  Parser::by_kind::by_kind(const by_kind &that) YY_NOEXCEPT : kind_(that.kind_) {}

  Parser::by_kind::by_kind(token_kind_type t) YY_NOEXCEPT : kind_(yytranslate_(t)) {}

  void Parser::by_kind::clear() YY_NOEXCEPT { kind_ = symbol_kind::S_YYEMPTY; }

  void Parser::by_kind::move(by_kind &that)
  {
    kind_ = that.kind_;
    that.clear();
  }

  Parser::symbol_kind_type Parser::by_kind::kind() const YY_NOEXCEPT { return kind_; }

  Parser::symbol_kind_type Parser::by_kind::type_get() const YY_NOEXCEPT { return this->kind(); }

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

  Parser::symbol_kind_type Parser::by_state::kind() const YY_NOEXCEPT
  {
    if (state == empty_state)
      return symbol_kind::S_YYEMPTY;
    else
      return YY_CAST(symbol_kind_type, yystos_[+state]);
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
    that.kind_ = symbol_kind::S_YYEMPTY;
  }

#if YY_CPLUSPLUS < 201103L
  Parser::stack_symbol_type &Parser::stack_symbol_type::operator=(const stack_symbol_type &that)
  {
    state = that.state;
    value = that.value;
    return *this;
  }

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
    YY_USE(yysym.kind());
  }

#if SEAMSDEBUG
  template <typename Base>
  void Parser::yy_print_(std::ostream &yyo, const basic_symbol<Base> &yysym) const
  {
    std::ostream &yyoutput = yyo;
    YY_USE(yyoutput);
    if (yysym.empty())
      yyo << "empty symbol";
    else {
      symbol_kind_type yykind = yysym.kind();
      yyo << (yykind < YYNTOKENS ? "token" : "nterm") << ' ' << yysym.name() << " (";
      YY_USE(yykind);
      yyo << ')';
    }
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

  void Parser::yypop_(int n) YY_NOEXCEPT { yystack_.pop(n); }

#if SEAMSDEBUG
  std::ostream &Parser::debug_stream() const { return *yycdebug_; }

  void Parser::set_debug_stream(std::ostream &o) { yycdebug_ = &o; }

  Parser::debug_level_type Parser::debug_level() const { return yydebug_; }

  void Parser::set_debug_level(debug_level_type l) { yydebug_ = l; }
#endif // SEAMSDEBUG

  Parser::state_type Parser::yy_lr_goto_state_(state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - YYNTOKENS] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - YYNTOKENS];
  }

  bool Parser::yy_pact_value_is_default_(int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yypact_ninf_;
  }

  bool Parser::yy_table_value_is_error_(int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yytable_ninf_;
  }

  int Parser::operator()() { return parse(); }

  int Parser::parse()
  {
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
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
      YYCDEBUG << "Entering state " << int(yystack_[0].state) << '\n';
      YY_STACK_PRINT();

      // Accept?
      if (yystack_[0].state == yyfinal_)
        YYACCEPT;

      goto yybackup;

    /*-----------.
    | yybackup.  |
    `-----------*/
    yybackup:
      // Try to take a decision without lookahead.
      yyn = yypact_[+yystack_[0].state];
      if (yy_pact_value_is_default_(yyn))
        goto yydefault;

      // Read a lookahead token.
      if (yyla.empty()) {
        YYCDEBUG << "Reading a token\n";
#if YY_EXCEPTIONS
        try
#endif // YY_EXCEPTIONS
        {
          yyla.kind_ = yytranslate_(yylex(&yyla.value));
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

      if (yyla.kind() == symbol_kind::S_YYerror) {
        // The scanner already issued an error message, process directly
        // to error recovery.  But do not keep the error token as
        // lookahead, it is too special and may lead us to an endless
        // loop in error recovery. */
        yyla.kind_ = symbol_kind::S_YYUNDEF;
        goto yyerrlab1;
      }

      /* If the proper action on seeing token YYLA.TYPE is to reduce or
         to detect an error, take that action.  */
      yyn += yyla.kind();
      if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.kind()) {
        goto yydefault;
      }

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
      yypush_("Shifting", state_type(yyn), YY_MOVE(yyla));
      goto yynewstate;

    /*-----------------------------------------------------------.
    | yydefault -- do the default action for the current state.  |
    `-----------------------------------------------------------*/
    yydefault:
      yyn = yydefact_[+yystack_[0].state];
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
          case 4: // line: '\n'
#line 129 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (echo)
              aprepro.lexer->LexerOutput("\n", 1);
          }
#line 632 "apr_parser.cc"
          break;

          case 5: // line: LBRACE exp RBRACE
#line 130 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (echo) {
              SEAMS::symrec *format = aprepro.getsym("_FORMAT");
              if (format->value.svar.empty()) {
#if defined FMT_SUPPORT
                auto tmpstr = fmt::format("{}", (yystack_[1].value.val));
                aprepro.lexer->LexerOutput(tmpstr.c_str(), tmpstr.size());
#else
                yyerror(aprepro, "Empty _FORMAT string -- no output will be printed. Optional "
                                 "Lib::FMT dependency is not enabled.");
#endif
              }
              else {
                static char tmpstr[512];
                int         len =
                    snprintf(tmpstr, 512, format->value.svar.c_str(), (yystack_[1].value.val));
                aprepro.lexer->LexerOutput(tmpstr, len);
              }
            }
          }
#line 644 "apr_parser.cc"
          break;

          case 6: // line: LBRACE sexp RBRACE
#line 137 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (echo && (yystack_[1].value.string) != NULL) {
              aprepro.lexer->LexerOutput((yystack_[1].value.string),
                                         strlen((yystack_[1].value.string)));
            }
          }
#line 653 "apr_parser.cc"
          break;

          case 7: // line: LBRACE aexp RBRACE
#line 141 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
          }
#line 659 "apr_parser.cc"
          break;

          case 8: // line: LBRACE RBRACE
#line 142 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
          }
#line 665 "apr_parser.cc"
          break;

          case 9: // line: error RBRACE
#line 143 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            yyerrok;
          }
#line 671 "apr_parser.cc"
          break;

          case 10: // bool: exp LT exp
#line 146 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) < (yystack_[0].value.val);
          }
#line 677 "apr_parser.cc"
          break;

          case 11: // bool: exp GT exp
#line 147 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) > (yystack_[0].value.val);
          }
#line 683 "apr_parser.cc"
          break;

          case 12: // bool: NOT exp
#line 148 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = !((yystack_[0].value.val));
          }
#line 689 "apr_parser.cc"
          break;

          case 13: // bool: exp LE exp
#line 149 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) <= (yystack_[0].value.val);
          }
#line 695 "apr_parser.cc"
          break;

          case 14: // bool: exp GE exp
#line 150 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) >= (yystack_[0].value.val);
          }
#line 701 "apr_parser.cc"
          break;

          case 15: // bool: exp EQ exp
#line 151 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) == (yystack_[0].value.val);
          }
#line 707 "apr_parser.cc"
          break;

          case 16: // bool: exp NE exp
#line 152 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) != (yystack_[0].value.val);
          }
#line 713 "apr_parser.cc"
          break;

          case 17: // bool: exp LOR exp
#line 153 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 719 "apr_parser.cc"
          break;

          case 18: // bool: exp LAND exp
#line 154 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 725 "apr_parser.cc"
          break;

          case 19: // bool: bool LOR bool
#line 155 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 731 "apr_parser.cc"
          break;

          case 20: // bool: bool LAND bool
#line 156 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 737 "apr_parser.cc"
          break;

          case 21: // bool: bool LOR exp
#line 157 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 743 "apr_parser.cc"
          break;

          case 22: // bool: bool LAND exp
#line 158 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 749 "apr_parser.cc"
          break;

          case 23: // bool: exp LOR bool
#line 159 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 755 "apr_parser.cc"
          break;

          case 24: // bool: exp LAND bool
#line 160 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 761 "apr_parser.cc"
          break;

          case 25: // bool: LPAR bool RPAR
#line 161 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          }
#line 767 "apr_parser.cc"
          break;

          case 26: // bool: sexp LT sexp
#line 164 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) < 0 ? 1 : 0);
          }
#line 773 "apr_parser.cc"
          break;

          case 27: // bool: sexp GT sexp
#line 165 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) > 0 ? 1 : 0);
          }
#line 779 "apr_parser.cc"
          break;

          case 28: // bool: sexp LE sexp
#line 166 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) <= 0 ? 1 : 0);
          }
#line 785 "apr_parser.cc"
          break;

          case 29: // bool: sexp GE sexp
#line 167 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) >= 0 ? 1 : 0);
          }
#line 791 "apr_parser.cc"
          break;

          case 30: // bool: sexp EQ sexp
#line 168 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) == 0 ? 1 : 0);
          }
#line 797 "apr_parser.cc"
          break;

          case 31: // bool: sexp NE sexp
#line 169 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) != 0 ? 1 : 0);
          }
#line 803 "apr_parser.cc"
          break;

          case 32: // aexp: AVAR
#line 171 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) = aprepro.make_array(*((yystack_[0].value.tptr)->value.avar));
          }
#line 809 "apr_parser.cc"
          break;

          case 33: // aexp: AFNCT LPAR sexp RPAR
#line 172 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_c == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_c))((yystack_[1].value.string));
            else
              yyerrok;
          }
#line 820 "apr_parser.cc"
          break;

          case 34: // aexp: AFNCT LPAR sexp COMMA exp RPAR
#line 178 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_cd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 831 "apr_parser.cc"
          break;

          case 35: // aexp: AFNCT LPAR sexp COMMA sexp RPAR
#line 184 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_cc == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              yyerrok;
          }
#line 842 "apr_parser.cc"
          break;

          case 36: // aexp: AFNCT LPAR exp COMMA exp COMMA exp RPAR
#line 190 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.arrfnct_ddd == NULL))
              (yylhs.value.arrval) = (*((yystack_[7].value.tptr)->value.arrfnct_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 853 "apr_parser.cc"
          break;

          case 37: // aexp: AFNCT LPAR exp COMMA exp RPAR
#line 196 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_dd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 864 "apr_parser.cc"
          break;

          case 38: // aexp: AFNCT LPAR exp RPAR
#line 202 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_d == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_d))((yystack_[1].value.val));
            else
              yyerrok;
          }
#line 875 "apr_parser.cc"
          break;

          case 39: // aexp: AFNCT LPAR aexp RPAR
#line 208 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_a == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_a))((yystack_[1].value.arrval));
            else
              yyerrok;
          }
#line 886 "apr_parser.cc"
          break;

          case 40: // aexp: SVAR EQUAL aexp
#line 214 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 895 "apr_parser.cc"
          break;

          case 41: // aexp: VAR EQUAL aexp
#line 218 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 904 "apr_parser.cc"
          break;

          case 42: // aexp: AVAR EQUAL aexp
#line 222 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) = (yystack_[0].value.arrval);
            aprepro.redefine_array((yystack_[2].value.tptr)->value.avar);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 912 "apr_parser.cc"
          break;

          case 43: // aexp: UNDVAR EQUAL aexp
#line 225 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 919 "apr_parser.cc"
          break;

          case 44: // aexp: aexp PLU aexp
#line 227 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
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
#line 932 "apr_parser.cc"
          break;

          case 45: // aexp: SUB aexp
#line 235 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), -1.0);
          }
#line 938 "apr_parser.cc"
          break;

          case 46: // aexp: aexp SUB aexp
#line 237 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
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
#line 951 "apr_parser.cc"
          break;

          case 47: // aexp: aexp TIM exp
#line 245 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) = array_scale((yystack_[2].value.arrval), (yystack_[0].value.val));
          }
#line 957 "apr_parser.cc"
          break;

          case 48: // aexp: aexp DIV exp
#line 246 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) =
                array_scale((yystack_[2].value.arrval), 1.0 / (yystack_[0].value.val));
          }
#line 963 "apr_parser.cc"
          break;

          case 49: // aexp: exp TIM aexp
#line 247 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), (yystack_[2].value.val));
          }
#line 969 "apr_parser.cc"
          break;

          case 50: // aexp: aexp TIM aexp
#line 248 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
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
#line 982 "apr_parser.cc"
          break;

          case 51: // sexp: QSTRING
#line 257 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = (yystack_[0].value.string);
          }
#line 988 "apr_parser.cc"
          break;

          case 52: // sexp: SVAR
#line 258 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          }
#line 994 "apr_parser.cc"
          break;

          case 53: // sexp: IMMSVAR
#line 259 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          }
#line 1000 "apr_parser.cc"
          break;

          case 54: // sexp: UNDVAR EQUAL sexp
#line 260 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            set_type(aprepro, (yystack_[2].value.tptr), Parser::token::SVAR);
          }
#line 1007 "apr_parser.cc"
          break;

          case 55: // sexp: SVAR EQUAL sexp
#line 262 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          }
#line 1015 "apr_parser.cc"
          break;

          case 56: // sexp: VAR EQUAL sexp
#line 265 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);
          }
#line 1024 "apr_parser.cc"
          break;

          case 57: // sexp: AVAR EQUAL sexp
#line 269 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = (yystack_[0].value.string);
            aprepro.redefine_array((yystack_[2].value.tptr)->value.avar);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);
          }
#line 1034 "apr_parser.cc"
          break;

          case 58: // sexp: IMMSVAR EQUAL sexp
#line 274 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = (char *)(yystack_[2].value.tptr)->value.svar.c_str();
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1040 "apr_parser.cc"
          break;

          case 59: // sexp: IMMVAR EQUAL sexp
#line 275 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          }
#line 1046 "apr_parser.cc"
          break;

          case 60: // sexp: SFNCT LPAR sexp RPAR
#line 276 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_c == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_c))(
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1057 "apr_parser.cc"
          break;

          case 61: // sexp: SFNCT LPAR RPAR
#line 282 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.strfnct == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[2].value.tptr)->value.strfnct))();
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1068 "apr_parser.cc"
          break;

          case 62: // sexp: SFNCT LPAR exp RPAR
#line 288 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_d == NULL))
              (yylhs.value.string) =
                  (char *)(*((yystack_[3].value.tptr)->value.strfnct_d))((yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1079 "apr_parser.cc"
          break;

          case 63: // sexp: SFNCT LPAR aexp RPAR
#line 294 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_a == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_a))(
                  (yystack_[1].value.arrval));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1090 "apr_parser.cc"
          break;

          case 64: // sexp: sexp CONCAT sexp
#line 300 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            concat_string((yystack_[2].value.string), (yystack_[0].value.string),
                          &(yylhs.value.string));
          }
#line 1096 "apr_parser.cc"
          break;

          case 65: // sexp: SFNCT LPAR exp COMMA exp RPAR
#line 301 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_dd == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1107 "apr_parser.cc"
          break;

          case 66: // sexp: SFNCT LPAR exp COMMA sexp COMMA sexp COMMA sexp COMMA sexp RPAR
#line 307 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.strfnct_dcccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[11].value.tptr)->value.strfnct_dcccc))(
                  (yystack_[9].value.val), (yystack_[7].value.string), (yystack_[5].value.string),
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1118 "apr_parser.cc"
          break;

          case 67: // sexp: SFNCT LPAR exp COMMA sexp COMMA sexp RPAR
#line 313 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_dcc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_dcc))(
                  (yystack_[5].value.val), (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1129 "apr_parser.cc"
          break;

          case 68: // sexp: SFNCT LPAR sexp COMMA sexp COMMA sexp RPAR
#line 319 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_ccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1140 "apr_parser.cc"
          break;

          case 69: // sexp: SFNCT LPAR sexp COMMA sexp RPAR
#line 325 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_cc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1151 "apr_parser.cc"
          break;

          case 70: // sexp: bool QUEST sexp COLON sexp
#line 331 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = ((yystack_[4].value.val)) ? ((yystack_[2].value.string))
                                                             : ((yystack_[0].value.string));
          }
#line 1157 "apr_parser.cc"
          break;

          case 71: // exp: NUM
#line 333 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          }
#line 1163 "apr_parser.cc"
          break;

          case 72: // exp: INC NUM
#line 334 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val) + 1;
          }
#line 1169 "apr_parser.cc"
          break;

          case 73: // exp: DEC NUM
#line 335 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val) - 1;
          }
#line 1175 "apr_parser.cc"
          break;

          case 74: // exp: VAR
#line 336 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          }
#line 1181 "apr_parser.cc"
          break;

          case 75: // exp: IMMVAR
#line 337 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          }
#line 1187 "apr_parser.cc"
          break;

          case 76: // exp: INC VAR
#line 338 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
          }
#line 1193 "apr_parser.cc"
          break;

          case 77: // exp: DEC VAR
#line 339 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
          }
#line 1199 "apr_parser.cc"
          break;

          case 78: // exp: VAR INC
#line 340 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
          }
#line 1205 "apr_parser.cc"
          break;

          case 79: // exp: VAR DEC
#line 341 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
          }
#line 1211 "apr_parser.cc"
          break;

          case 80: // exp: VAR EQUAL exp
#line 342 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          }
#line 1218 "apr_parser.cc"
          break;

          case 81: // exp: SVAR EQUAL exp
#line 344 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1226 "apr_parser.cc"
          break;

          case 82: // exp: AVAR EQUAL exp
#line 347 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val);
            aprepro.redefine_array((yystack_[2].value.tptr)->value.avar);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1236 "apr_parser.cc"
          break;

          case 83: // exp: VAR EQ_PLUS exp
#line 352 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1242 "apr_parser.cc"
          break;

          case 84: // exp: VAR EQ_MINUS exp
#line 353 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1248 "apr_parser.cc"
          break;

          case 85: // exp: VAR EQ_TIME exp
#line 354 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1254 "apr_parser.cc"
          break;

          case 86: // exp: VAR EQ_DIV exp
#line 355 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1260 "apr_parser.cc"
          break;

          case 87: // exp: VAR EQ_POW exp
#line 356 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            SEAMS::math_error(aprepro, "Power");
          }
#line 1270 "apr_parser.cc"
          break;

          case 88: // exp: INC IMMVAR
#line 361 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          }
#line 1276 "apr_parser.cc"
          break;

          case 89: // exp: DEC IMMVAR
#line 362 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          }
#line 1282 "apr_parser.cc"
          break;

          case 90: // exp: IMMVAR INC
#line 363 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          }
#line 1288 "apr_parser.cc"
          break;

          case 91: // exp: IMMVAR DEC
#line 364 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          }
#line 1294 "apr_parser.cc"
          break;

          case 92: // exp: IMMVAR EQUAL exp
#line 365 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1300 "apr_parser.cc"
          break;

          case 93: // exp: IMMSVAR EQUAL exp
#line 366 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          }
#line 1306 "apr_parser.cc"
          break;

          case 94: // exp: IMMVAR EQ_PLUS exp
#line 367 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1312 "apr_parser.cc"
          break;

          case 95: // exp: IMMVAR EQ_MINUS exp
#line 368 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1318 "apr_parser.cc"
          break;

          case 96: // exp: IMMVAR EQ_TIME exp
#line 369 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1324 "apr_parser.cc"
          break;

          case 97: // exp: IMMVAR EQ_DIV exp
#line 370 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1330 "apr_parser.cc"
          break;

          case 98: // exp: IMMVAR EQ_POW exp
#line 371 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1336 "apr_parser.cc"
          break;

          case 99: // exp: UNDVAR
#line 373 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1343 "apr_parser.cc"
          break;

          case 100: // exp: INC UNDVAR
#line 375 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1351 "apr_parser.cc"
          break;

          case 101: // exp: DEC UNDVAR
#line 378 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1359 "apr_parser.cc"
          break;

          case 102: // exp: UNDVAR INC
#line 381 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          }
#line 1367 "apr_parser.cc"
          break;

          case 103: // exp: UNDVAR DEC
#line 384 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          }
#line 1375 "apr_parser.cc"
          break;

          case 104: // exp: UNDVAR EQUAL exp
#line 387 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1382 "apr_parser.cc"
          break;

          case 105: // exp: UNDVAR EQ_PLUS exp
#line 389 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1390 "apr_parser.cc"
          break;

          case 106: // exp: UNDVAR EQ_MINUS exp
#line 392 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1398 "apr_parser.cc"
          break;

          case 107: // exp: UNDVAR EQ_TIME exp
#line 395 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1406 "apr_parser.cc"
          break;

          case 108: // exp: UNDVAR EQ_DIV exp
#line 398 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1414 "apr_parser.cc"
          break;

          case 109: // exp: UNDVAR EQ_POW exp
#line 401 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            SEAMS::math_error(aprepro, "Power");
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1425 "apr_parser.cc"
          break;

          case 110: // exp: FNCT LPAR RPAR
#line 408 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.fnctptr == NULL))
              (yylhs.value.val) = (*((yystack_[2].value.tptr)->value.fnctptr))();
            else
              (yylhs.value.val) = 0.0;
          }
#line 1436 "apr_parser.cc"
          break;

          case 111: // exp: FNCT LPAR exp RPAR
#line 415 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_d == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_d))((yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1447 "apr_parser.cc"
          break;

          case 112: // exp: FNCT LPAR sexp RPAR
#line 422 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_c == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_c))((yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1458 "apr_parser.cc"
          break;

          case 113: // exp: FNCT LPAR aexp RPAR
#line 429 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_a == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_a))((yystack_[1].value.arrval));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1469 "apr_parser.cc"
          break;

          case 114: // exp: FNCT LPAR sexp COMMA exp RPAR
#line 436 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1480 "apr_parser.cc"
          break;

          case 115: // exp: FNCT LPAR exp COMMA sexp RPAR
#line 443 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dc))(
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1491 "apr_parser.cc"
          break;

          case 116: // exp: FNCT LPAR sexp COMMA sexp RPAR
#line 450 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1502 "apr_parser.cc"
          break;

          case 117: // exp: FNCT LPAR sexp COMMA sexp COMMA sexp RPAR
#line 457 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccc == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              yyerrok;
          }
#line 1513 "apr_parser.cc"
          break;

          case 118: // exp: FNCT LPAR exp COMMA exp RPAR
#line 464 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1524 "apr_parser.cc"
          break;

          case 119: // exp: FNCT LPAR exp COMMA exp COMMA exp RPAR
#line 470 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ddd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1535 "apr_parser.cc"
          break;

          case 120: // exp: FNCT LPAR sexp COMMA sexp COMMA exp RPAR
#line 476 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccd))(
                  (yystack_[5].value.string), (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1546 "apr_parser.cc"
          break;

          case 121: // exp: FNCT LPAR exp COMMA exp SEMI exp COMMA exp RPAR
#line 482 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1557 "apr_parser.cc"
          break;

          case 122: // exp: FNCT LPAR exp COMMA exp COMMA exp COMMA exp RPAR
#line 488 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1568 "apr_parser.cc"
          break;

          case 123: // exp: FNCT LPAR exp COMMA exp COMMA exp COMMA exp COMMA sexp RPAR
#line 494 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.fnctptr_ddddc == NULL))
              (yylhs.value.val) = (*((yystack_[11].value.tptr)->value.fnctptr_ddddc))(
                  (yystack_[9].value.val), (yystack_[7].value.val), (yystack_[5].value.val),
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1579 "apr_parser.cc"
          break;

          case 124: // exp: FNCT LPAR exp COMMA exp COMMA exp COMMA exp COMMA exp COMMA exp RPAR
#line 500 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[13].value.tptr),
                          (yystack_[13].value.tptr)->value.fnctptr_dddddd == NULL))
              (yylhs.value.val) = (*((yystack_[13].value.tptr)->value.fnctptr_dddddd))(
                  (yystack_[11].value.val), (yystack_[9].value.val), (yystack_[7].value.val),
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1590 "apr_parser.cc"
          break;

          case 125: // exp: exp PLU exp
#line 506 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) + (yystack_[0].value.val);
          }
#line 1596 "apr_parser.cc"
          break;

          case 126: // exp: exp SUB exp
#line 507 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) - (yystack_[0].value.val);
          }
#line 1602 "apr_parser.cc"
          break;

          case 127: // exp: exp TIM exp
#line 508 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) * (yystack_[0].value.val);
          }
#line 1608 "apr_parser.cc"
          break;

          case 128: // exp: exp DIV exp
#line 509 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if ((yystack_[0].value.val) == 0.) {
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (yystack_[2].value.val) / (yystack_[0].value.val);
          }
#line 1620 "apr_parser.cc"
          break;

          case 129: // exp: exp MOD exp
#line 516 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if ((yystack_[0].value.val) == 0.) {
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (int)(yystack_[2].value.val) % (int)(yystack_[0].value.val);
          }
#line 1632 "apr_parser.cc"
          break;

          case 130: // exp: SUB exp
#line 523 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = -(yystack_[0].value.val);
          }
#line 1638 "apr_parser.cc"
          break;

          case 131: // exp: PLU exp
#line 524 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          }
#line 1644 "apr_parser.cc"
          break;

          case 132: // exp: exp POW exp
#line 525 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            reset_error();
            (yylhs.value.val) = std::pow((yystack_[2].value.val), (yystack_[0].value.val));
            SEAMS::math_error(aprepro, "Power");
          }
#line 1652 "apr_parser.cc"
          break;

          case 133: // exp: LPAR exp RPAR
#line 528 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          }
#line 1658 "apr_parser.cc"
          break;

          case 134: // exp: LBRACK exp RBRACK
#line 529 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            reset_error();
            (yylhs.value.val) =
                (double)((yystack_[1].value.val) < 0 ? -floor(-((yystack_[1].value.val)))
                                                     : floor((yystack_[1].value.val)));
            SEAMS::math_error(aprepro, "floor (int)");
          }
#line 1666 "apr_parser.cc"
          break;

          case 135: // exp: bool
#line 532 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[0].value.val)) ? 1 : 0;
          }
#line 1672 "apr_parser.cc"
          break;

          case 136: // exp: bool QUEST exp COLON exp
#line 533 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                ((yystack_[4].value.val)) ? ((yystack_[2].value.val)) : ((yystack_[0].value.val));
          }
#line 1678 "apr_parser.cc"
          break;

          case 137: // exp: AVAR LBRACK exp RBRACK
#line 534 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                array_value((yystack_[3].value.tptr)->value.avar, (yystack_[1].value.val), 0);
          }
#line 1684 "apr_parser.cc"
          break;

          case 138: // exp: AVAR LBRACK exp COMMA exp RBRACK
#line 535 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = array_value((yystack_[5].value.tptr)->value.avar,
                                            (yystack_[3].value.val), (yystack_[1].value.val));
          }
#line 1690 "apr_parser.cc"
          break;

          case 139: // exp: AVAR LBRACK exp RBRACK EQUAL exp
#line 537 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
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
#line 1716 "apr_parser.cc"
          break;

          case 140: // exp: AVAR LBRACK exp COMMA exp RBRACK EQUAL exp
#line 559 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
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
#line 1740 "apr_parser.cc"
          break;

#line 1744 "apr_parser.cc"

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
        context     yyctx(*this, yyla);
        std::string msg = yysyntax_error_(yyctx);
        error(YY_MOVE(msg));
      }

      if (yyerrstatus_ == 3) {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.kind() == symbol_kind::S_YYEOF)
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
      YY_STACK_PRINT();
      goto yyerrlab1;

    /*-------------------------------------------------------------.
    | yyerrlab1 -- common code for both syntax error and YYERROR.  |
    `-------------------------------------------------------------*/
    yyerrlab1:
      yyerrstatus_ = 3; // Each real token shifted decrements this.
      // Pop stack until we find a state that shifts the error token.
      for (;;) {
        yyn = yypact_[+yystack_[0].state];
        if (!yy_pact_value_is_default_(yyn)) {
          yyn += symbol_kind::S_YYerror;
          if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == symbol_kind::S_YYerror) {
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
      {
        stack_symbol_type error_token;

        // Shift the error token.
        error_token.state = state_type(yyn);
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
      YY_STACK_PRINT();
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

  std::string Parser::symbol_name(symbol_kind_type yysymbol)
  {
    return yytnamerr_(yytname_[yysymbol]);
  }

  // Parser::context.
  Parser::context::context(const Parser &yyparser, const symbol_type &yyla)
      : yyparser_(yyparser), yyla_(yyla)
  {
  }

  int Parser::context::expected_tokens(symbol_kind_type yyarg[], int yyargn) const
  {
    // Actual number of expected tokens
    int yycount = 0;

    const int yyn = yypact_[+yyparser_.yystack_[0].state];
    if (!yy_pact_value_is_default_(yyn)) {
      /* Start YYX at -YYN if negative to avoid negative indexes in
         YYCHECK.  In other words, skip the first -YYN actions for
         this state because they are default actions.  */
      const int yyxbegin = yyn < 0 ? -yyn : 0;
      // Stay within bounds of both yycheck and yytname.
      const int yychecklim = yylast_ - yyn + 1;
      const int yyxend     = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      for (int yyx = yyxbegin; yyx < yyxend; ++yyx)
        if (yycheck_[yyx + yyn] == yyx && yyx != symbol_kind::S_YYerror &&
            !yy_table_value_is_error_(yytable_[yyx + yyn])) {
          if (!yyarg)
            ++yycount;
          else if (yycount == yyargn)
            return 0;
          else
            yyarg[yycount++] = YY_CAST(symbol_kind_type, yyx);
        }
    }

    if (yyarg && yycount == 0 && 0 < yyargn)
      yyarg[0] = symbol_kind::S_YYEMPTY;
    return yycount;
  }

  int Parser::yy_syntax_error_arguments_(const context &yyctx, symbol_kind_type yyarg[],
                                         int yyargn) const
  {
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
         scanner and before detecting a syntax error.  Thus, state merging
         (from LALR or IELR) and default reductions corrupt the expected
         token list.  However, the list is correct for canonical LR with
         one exception: it will still contain any token that will not be
         accepted due to an error action in a later state.
    */

    if (!yyctx.lookahead().empty()) {
      if (yyarg)
        yyarg[0] = yyctx.token();
      int yyn = yyctx.expected_tokens(yyarg ? yyarg + 1 : yyarg, yyargn - 1);
      return yyn + 1;
    }
    return 0;
  }

  // Generate an error message.
  std::string Parser::yysyntax_error_(const context &yyctx) const
  {
    // Its maximum.
    enum { YYARGS_MAX = 5 };
    // Arguments of yyformat.
    symbol_kind_type yyarg[YYARGS_MAX];
    int              yycount = yy_syntax_error_arguments_(yyctx, yyarg, YYARGS_MAX);

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
    std::ptrdiff_t yyi = 0;
    for (char const *yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount) {
        yyres += symbol_name(yyarg[yyi++]);
        ++yyp;
      }
      else
        yyres += *yyp;
    return yyres;
  }

  const signed char Parser::yypact_ninf_ = -37;

  const signed char Parser::yytable_ninf_ = -1;

  const short Parser::yypact_[] = {
      -37,  2,    -37,  -7,   281,  -37,  -37,  -37,  -37,  -37,  233,  284,  0,    311,  11,
      -5,   19,   27,   35,   399,  399,  -37,  399,  384,  399,  116,  163,  175,  178,  41,
      1151, 384,  399,  399,  399,  399,  399,  -37,  -37,  384,  399,  399,  399,  399,  399,
      -37,  -37,  384,  399,  399,  399,  399,  399,  399,  -37,  -37,  399,  399,  384,  225,
      339,  384,  482,  625,  36,   73,   399,  101,  1199, 889,  1103, 14,   -37,  14,   14,
      -37,  -37,  -37,  -37,  -37,  -37,  -37,  -37,  399,  399,  399,  -37,  384,  384,  399,
      384,  -37,  399,  399,  399,  399,  399,  399,  399,  -37,  399,  399,  399,  399,  399,
      399,  399,  399,  399,  399,  399,  384,  399,  399,  29,   1199, 1218, 1234, 1234, 1234,
      1234, 1234, 29,   1199, 1218, 1234, 1234, 1234, 1234, 1234, 29,   1199, 1218, 1199, 1234,
      1234, 1234, 1234, 1234, 1234, 1199, 1234, 624,  29,   1199, 1218, -37,  -15,  78,   654,
      -37,  25,   418,  684,  100,  425,  714,  399,  399,  399,  399,  14,   -37,  -37,  399,
      -37,  1165, 1185, 51,   1265, -37,  1279, -36,  1218, -36,  1250, -37,  1250, 12,   1234,
      12,   12,   12,   12,   12,   -37,  51,   1265, -37,  1279, -31,  -31,  -31,  -31,  -31,
      -31,  155,  155,  14,   -37,  14,   14,   14,   399,  69,   -37,  399,  -37,  399,  -37,
      -37,  399,  -37,  399,  -37,  -37,  399,  -37,  399,  -37,  1234, 1234, 1234, 1234, 14,
      399,  399,  1128, 399,  449,  916,  141,  595,  455,  485,  943,  515,  970,  744,  1199,
      1234, 96,   1234, 399,  -37,  -37,  -37,  399,  -37,  399,  399,  -37,  399,  -37,  -37,
      -37,  399,  -37,  399,  537,  997,  774,  833,  543,  479,  1024, 1234, -37,  -37,  399,
      -37,  399,  -37,  399,  -37,  -37,  804,  1051, 509,  399,  -37,  -37,  399,  565,  862,
      571,  -37,  399,  -37,  1078, -37};

  const unsigned char Parser::yydefact_[] = {
      2,   0,   1,   0,   0,   4,   3,   9,   71,  51,  99,  74,  52,  75,  53,  32,  0,   0,   0,
      0,   0,   8,   0,   0,   0,   0,   0,   135, 0,   0,   0,   0,   0,   0,   0,   0,   0,   102,
      103, 0,   0,   0,   0,   0,   0,   78,  79,  0,   0,   0,   0,   0,   0,   0,   90,  91,  0,
      0,   0,   0,   0,   0,   99,  74,  52,  0,   0,   135, 0,   0,   0,   131, 45,  130, 12,  72,
      100, 76,  88,  73,  101, 77,  89,  0,   0,   0,   7,   0,   0,   0,   0,   6,   0,   0,   0,
      0,   0,   0,   0,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      43,  54,  104, 105, 106, 107, 108, 109, 41,  56,  80,  83,  84,  85,  86,  87,  40,  55,  81,
      59,  92,  94,  95,  96,  97,  98,  58,  93,  0,   42,  57,  82,  110, 0,   0,   0,   61,  0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   130, 25,  133, 0,   134, 0,   0,   19,  21,  20,
      22,  44,  0,   46,  48,  50,  47,  26,  0,   27,  28,  29,  30,  31,  64,  23,  17,  24,  18,
      10,  11,  13,  14,  15,  16,  125, 126, 128, 49,  127, 129, 132, 0,   137, 113, 0,   112, 0,
      111, 63,  0,   60,  0,   62,  39,  0,   33,  0,   38,  104, 80,  81,  82,  127, 0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   70,  136, 138, 139, 0,   116, 114, 115,
      0,   118, 0,   0,   69,  0,   65,  35,  34,  0,   37,  0,   0,   0,   0,   0,   0,   0,   0,
      140, 117, 120, 0,   119, 0,   68,  0,   67,  36,  0,   0,   0,   0,   122, 121, 0,   0,   0,
      0,   123, 0,   66,  0,   124};

  const signed char Parser::yypgoto_[] = {-37, -37, -37, -13, 104, 89, -4};

  const unsigned char Parser::yydefgoto_[] = {0, 1, 6, 27, 28, 68, 179};

  const short Parser::yytable_[] = {
      30,  205, 2,   3,   89,  90,  67,  108, 109, 110, 164, 112, 57,  7,   113, 69,  70,  58,  71,
      73,  74,  4,   47,  87,  88,  89,  90,  116, 117, 118, 119, 120, 121, 56,  59,  124, 125, 126,
      127, 128, 129, 210, 60,  132, 134, 135, 136, 137, 138, 139, 61,  5,   141, 142, 145, 149, 153,
      156, 159, 113, 98,  91,  161, 87,  88,  89,  90,  87,  88,  89,  90,  168, 170, 92,  93,  94,
      95,  96,  97,  167, 169, 171, 85,  173, 173, 175, 177, 186, 188, 98,  57,  228, 206, 29,  207,
      160, 187, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 200, 201, 202, 92,  93,  94,  95,
      96,  97,  215, 162, 258, 75,  115, 76,  77,  0,   78,  0,   98,  72,  123, 83,  0,   84,  85,
      0,   0,   114, 131, 133, 87,  88,  89,  90,  0,   122, 0,   140, 0,   144, 148, 152, 155, 130,
      0,   220, 221, 222, 223, 246, 0,   0,   224, 0,   143, 147, 151, 154, 79,  0,   80,  81,  0,
      82,  166, 92,  93,  94,  95,  96,  97,  0,   0,   178, 180, 181, 182, 183, 184, 185, 0,   98,
      0,   172, 174, 0,   176, 110, 164, 112, 86,  227, 113, 0,   230, 83,  232, 84,  85,  0,   0,
      235, 0,   0,   237, 0,   238, 199, 87,  88,  89,  90,  0,   0,   240, 0,   242, 0,   0,   0,
      8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  260, 19,  146, 20,  261, 0,   262, 115,
      123, 131, 144, 0,   0,   265, 0,   266, 31,  32,  33,  34,  35,  36,  0,   0,   22,  23,  276,
      0,   277, 0,   24,  0,   25,  26,  0,   0,   284, 0,   0,   0,   37,  38,  0,   0,   289, 8,
      9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  229, 19,  231, 20,  0,   233, 21,  234, 0,
      0,   236, 39,  40,  41,  42,  43,  44,  0,   0,   239, 0,   0,   0,   0,   22,  23,  0,   0,
      0,   0,   24,  0,   25,  26,  0,   45,  46,  259, 48,  49,  50,  51,  52,  53,  263, 0,   264,
      8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  150, 20,  54,  55,  0,   0,
      0,   278, 0,   0,   0,   0,   0,   283, 0,   0,   285, 0,   0,   0,   0,   0,   22,  23,  0,
      0,   0,   0,   24,  0,   25,  26,  8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,
      19,  0,   20,  8,   9,   62,  63,  64,  13,  14,  65,  16,  17,  0,   0,   19,  0,   20,  0,
      0,   0,   0,   0,   22,  23,  0,   0,   0,   0,   24,  0,   25,  26,  211, 0,   212, 0,   0,
      22,  66,  216, 0,   217, 0,   24,  0,   25,  26,  0,   0,   0,   92,  93,  94,  95,  96,  97,
      0,   92,  93,  94,  95,  96,  97,  243, 0,   244, 98,  0,   0,   250, 0,   251, 0,   98,  0,
      0,   0,   0,   0,   0,   0,   92,  93,  94,  95,  96,  97,  92,  93,  94,  95,  96,  97,  273,
      0,   274, 0,   98,  0,   252, 0,   0,   0,   98,  157, 32,  33,  34,  35,  36,  0,   92,  93,
      94,  95,  96,  97,  92,  93,  94,  95,  96,  97,  282, 0,   0,   0,   98,  37,  38,  0,   254,
      0,   98,  0,   0,   0,   0,   0,   0,   0,   92,  93,  94,  95,  96,  97,  92,  93,  94,  95,
      96,  97,  267, 0,   0,   0,   98,  0,   272, 0,   0,   0,   98,  0,   0,   0,   0,   0,   92,
      93,  94,  95,  96,  97,  92,  93,  94,  95,  96,  97,  286, 0,   0,   0,   98,  0,   288, 0,
      0,   0,   98,  0,   0,   0,   0,   0,   92,  93,  94,  95,  96,  97,  92,  93,  94,  95,  96,
      97,  247, 0,   248, 0,   98,  0,   0,   249, 0,   0,   98,  0,   0,   0,   0,   0,   100, 101,
      102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 203, 0,   113, 0,   204, 0,   0,   0,
      0,   158, 40,  41,  42,  43,  44,  0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
      164, 112, 0,   208, 113, 209, 45,  46,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 0,   213, 113, 214, 0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107,
      108, 109, 110, 111, 112, 0,   218, 113, 219, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 0,   256, 113,
      257, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104,
      105, 106, 107, 108, 109, 110, 164, 112, 0,   269, 113, 270, 0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112,
      0,   279, 113, 280, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101,
      102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 271, 0,   113, 0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
      164, 112, 287, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100,
      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 163, 0,   113, 0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164,
      112, 245, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103,
      104, 105, 106, 107, 108, 109, 110, 164, 112, 253, 0,   113, 0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 255, 0,
      113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106,
      107, 108, 109, 110, 164, 112, 268, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 275, 0,   113, 0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
      110, 164, 112, 281, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101,
      102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 290, 0,   113, 0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112,
      165, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106,
      107, 108, 109, 110, 164, 112, 241, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   100,
      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 99,  0,   113, 0,   0,   0,   0,
      0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 225, 0,   113,
      92,  93,  94,  95,  96,  97,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  226, 100,
      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 0,   0,   113, 92,  93,  94,  95,
      96,  97,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  100, 101, 102, 103, 104, 105,
      106, 107, 108, 109, 110, 111, 112, 0,   0,   113, 100, 101, 102, 103, 104, 105, 106, 107, 108,
      109, 110, 164, 112, 0,   0,   113, 100, 101, 102, 103, 104, 105, 106, 107, 0,   0,   0,   0,
      112, 0,   0,   113, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 0,   0,   113,
      102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 0,   0,   113};

  const short Parser::yycheck_[] = {
      4,   16,  0,   1,   40,  41,  19,  38,  39,  40,  41,  42,  17,  20,  45,  19,  20,  22,  22,
      23,  24,  19,  22,  38,  39,  40,  41,  31,  32,  33,  34,  35,  36,  22,  15,  39,  40,  41,
      42,  43,  44,  16,  15,  47,  48,  49,  50,  51,  52,  53,  15,  49,  56,  57,  58,  59,  60,
      61,  22,  45,  48,  20,  66,  38,  39,  40,  41,  38,  39,  40,  41,  84,  85,  32,  33,  34,
      35,  36,  37,  83,  84,  85,  31,  87,  88,  89,  90,  100, 101, 48,  17,  22,  14,  4,   16,
      22,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 32,  33,  34,  35,
      36,  37,  16,  16,  22,  3,   31,  5,   6,   -1,  8,   -1,  48,  23,  39,  28,  -1,  30,  31,
      -1,  -1,  31,  47,  48,  38,  39,  40,  41,  -1,  39,  -1,  56,  -1,  58,  59,  60,  61,  47,
      -1,  157, 158, 159, 160, 16,  -1,  -1,  164, -1,  58,  59,  60,  61,  3,   -1,  5,   6,   -1,
      8,   83,  32,  33,  34,  35,  36,  37,  -1,  -1,  92,  93,  94,  95,  96,  97,  98,  -1,  48,
      -1,  87,  88,  -1,  90,  40,  41,  42,  20,  203, 45,  -1,  206, 28,  208, 30,  31,  -1,  -1,
      213, -1,  -1,  216, -1,  218, 111, 38,  39,  40,  41,  -1,  -1,  226, -1,  228, -1,  -1,  -1,
      3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  243, 15,  16,  17,  247, -1,  249, 157,
      158, 159, 160, -1,  -1,  256, -1,  258, 22,  23,  24,  25,  26,  27,  -1,  -1,  38,  39,  269,
      -1,  271, -1,  44,  -1,  46,  47,  -1,  -1,  279, -1,  -1,  -1,  46,  47,  -1,  -1,  287, 3,
      4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  206, 15,  208, 17,  -1,  211, 20,  213, -1,
      -1,  216, 22,  23,  24,  25,  26,  27,  -1,  -1,  225, -1,  -1,  -1,  -1,  38,  39,  -1,  -1,
      -1,  -1,  44,  -1,  46,  47,  -1,  46,  47,  243, 22,  23,  24,  25,  26,  27,  250, -1,  252,
      3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,  16,  17,  46,  47,  -1,  -1,
      -1,  273, -1,  -1,  -1,  -1,  -1,  279, -1,  -1,  282, -1,  -1,  -1,  -1,  -1,  38,  39,  -1,
      -1,  -1,  -1,  44,  -1,  46,  47,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,
      15,  -1,  17,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  -1,
      -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  14,  -1,  16,  -1,  -1,
      38,  39,  14,  -1,  16,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,
      -1,  32,  33,  34,  35,  36,  37,  14,  -1,  16,  48,  -1,  -1,  14,  -1,  16,  -1,  48,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  32,  33,  34,  35,  36,  37,  14,
      -1,  16,  -1,  48,  -1,  14,  -1,  -1,  -1,  48,  22,  23,  24,  25,  26,  27,  -1,  32,  33,
      34,  35,  36,  37,  32,  33,  34,  35,  36,  37,  14,  -1,  -1,  -1,  48,  46,  47,  -1,  16,
      -1,  48,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  32,  33,  34,  35,
      36,  37,  16,  -1,  -1,  -1,  48,  -1,  16,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  32,
      33,  34,  35,  36,  37,  32,  33,  34,  35,  36,  37,  16,  -1,  -1,  -1,  48,  -1,  16,  -1,
      -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  32,  33,  34,  35,  36,
      37,  14,  -1,  16,  -1,  48,  -1,  -1,  21,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  30,  31,
      32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  14,  -1,  45,  -1,  18,  -1,  -1,  -1,
      -1,  22,  23,  24,  25,  26,  27,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
      41,  42,  -1,  14,  45,  16,  46,  47,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,
      38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,
      16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,
      35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
      -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,
      32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  14,  -1,  45,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
      41,  42,  14,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
      42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,
      34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,
      45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,
      37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
      40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,
      32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
      18,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,
      37,  38,  39,  40,  41,  42,  18,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  20,  -1,  45,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  29,  -1,  45,
      32,  33,  34,  35,  36,  37,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  48,  29,  30,
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  32,  33,  34,  35,
      36,  37,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  48,  30,  31,  32,  33,  34,  35,
      36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  30,  31,  32,  33,  34,  35,  36,  37,  38,
      39,  40,  41,  42,  -1,  -1,  45,  30,  31,  32,  33,  34,  35,  36,  37,  -1,  -1,  -1,  -1,
      42,  -1,  -1,  45,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,
      32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45};

  const signed char Parser::yystos_[] = {
      0,  51, 0,  1,  19, 49, 52, 20, 3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 15, 17, 20, 38,
      39, 44, 46, 47, 53, 54, 55, 56, 22, 23, 24, 25, 26, 27, 46, 47, 22, 23, 24, 25, 26, 27, 46,
      47, 22, 22, 23, 24, 25, 26, 27, 46, 47, 22, 17, 22, 15, 15, 15, 5,  6,  7,  10, 39, 53, 55,
      56, 56, 56, 54, 56, 56, 3,  5,  6,  8,  3,  5,  6,  8,  28, 30, 31, 20, 38, 39, 40, 41, 20,
      32, 33, 34, 35, 36, 37, 48, 20, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 45, 54,
      55, 56, 56, 56, 56, 56, 56, 54, 55, 56, 56, 56, 56, 56, 56, 54, 55, 56, 55, 56, 56, 56, 56,
      56, 56, 55, 56, 56, 54, 55, 56, 16, 54, 55, 56, 16, 54, 55, 56, 54, 55, 56, 22, 22, 22, 22,
      56, 16, 16, 41, 18, 55, 56, 53, 56, 53, 56, 54, 56, 54, 56, 54, 56, 55, 56, 55, 55, 55, 55,
      55, 55, 53, 56, 53, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 54, 56, 56, 56, 14, 18, 16, 14,
      16, 14, 16, 16, 14, 16, 14, 16, 16, 14, 16, 14, 16, 56, 56, 56, 56, 56, 29, 29, 56, 22, 55,
      56, 55, 56, 55, 55, 56, 55, 56, 56, 55, 56, 18, 56, 14, 16, 16, 16, 14, 16, 21, 14, 16, 14,
      16, 16, 16, 14, 16, 22, 55, 56, 56, 56, 55, 55, 56, 56, 16, 16, 14, 16, 14, 16, 14, 16, 16,
      56, 56, 55, 14, 16, 16, 14, 55, 56, 55, 16, 14, 16, 56, 16};

  const signed char Parser::yyr1_[] = {
      0,  50, 51, 51, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
      54, 54, 54, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55,
      55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56};

  const signed char Parser::yyr2_[] = {
      0, 2, 0, 2, 1, 3,  3,  3,  2,  2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 1, 4, 6,  6,  8,  6,  4, 4, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 1, 1, 1, 3, 3, 3, 3,
      3, 3, 4, 3, 4, 4,  3,  6,  12, 8, 8, 6, 5, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
      3, 2, 2, 2, 2, 3,  3,  3,  3,  3, 3, 3, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 6, 6,
      6, 8, 6, 8, 8, 10, 10, 12, 14, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 1, 5, 4, 6, 6, 8};

#if SEAMSDEBUG || 1
  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a YYNTOKENS, nonterminals.
  const char *const Parser::yytname_[] = {"\"end of file\"",
                                          "error",
                                          "\"invalid token\"",
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
#endif

#if SEAMSDEBUG
  const short Parser::yyrline_[] = {
      0,   125, 125, 126, 129, 130, 137, 141, 142, 143, 146, 147, 148, 149, 150, 151, 152, 153,
      154, 155, 156, 157, 158, 159, 160, 161, 164, 165, 166, 167, 168, 169, 171, 172, 178, 184,
      190, 196, 202, 208, 214, 218, 222, 225, 227, 235, 237, 245, 246, 247, 248, 257, 258, 259,
      260, 262, 265, 269, 274, 275, 276, 282, 288, 294, 300, 301, 307, 313, 319, 325, 331, 333,
      334, 335, 336, 337, 338, 339, 340, 341, 342, 344, 347, 352, 353, 354, 355, 356, 361, 362,
      363, 364, 365, 366, 367, 368, 369, 370, 371, 373, 375, 378, 381, 384, 387, 389, 392, 395,
      398, 401, 408, 415, 422, 429, 436, 443, 450, 457, 464, 470, 476, 482, 488, 494, 500, 506,
      507, 508, 509, 516, 523, 524, 525, 528, 529, 532, 533, 534, 535, 536, 558};

  void Parser::yy_stack_print_() const
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator i = yystack_.begin(), i_end = yystack_.end(); i != i_end; ++i)
      *yycdebug_ << ' ' << int(i->state);
    *yycdebug_ << '\n';
  }

  void Parser::yy_reduce_print_(int yyrule) const
  {
    int yylno  = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1 << " (line " << yylno << "):\n";
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT("   $" << yyi + 1 << " =", yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // SEAMSDEBUG

  Parser::symbol_kind_type Parser::yytranslate_(int t) YY_NOEXCEPT
  {
    // YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to
    // TOKEN-NUM as returned by yylex.
    static const signed char translate_table[] = {
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
    // Last valid token kind.
    const int code_max = 303;

    if (t <= 0)
      return symbol_kind::S_YYEOF;
    else if (t <= code_max)
      return static_cast<symbol_kind_type>(translate_table[t]);
    else
      return symbol_kind::S_YYUNDEF;
  }

} // namespace SEAMS
#line 2649 "apr_parser.cc"

#line 581 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"

void SEAMS::Parser::error(const std::string &m) { aprepro.error(m); }
