// A Bison parser, made by GNU Bison 3.8.2.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2021 Free Software Foundation, Inc.

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

#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fmt/format.h>
#include <fmt/printf.h>
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

#line 84 "apr_parser.cc"

#include "aprepro_parser.h"

// Second part of user prologue.
#line 112 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"

#include "apr_scanner.h"
#include "aprepro.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex

#line 103 "apr_parser.cc"

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
#line 177 "apr_parser.cc"

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

  Parser::~Parser() {}

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
#line 131 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (echo)
              aprepro.lexer->LexerOutput("\n", 1);
          }
#line 634 "apr_parser.cc"
          break;

          case 5: // line: LBRACE exp RBRACE
#line 132 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (echo) {
              SEAMS::symrec *format = aprepro.getsym("_FORMAT");
              if (format->value.svar.empty()) {
                auto tmpstr = fmt::format("{}", (yystack_[1].value.val));
                aprepro.lexer->LexerOutput(tmpstr.c_str(), tmpstr.size());
              }
              else {
                static char tmpstr[512];
                int         len =
                    snprintf(tmpstr, 512, format->value.svar.c_str(), (yystack_[1].value.val));
                aprepro.lexer->LexerOutput(tmpstr, len);
              }
            }
          }
#line 652 "apr_parser.cc"
          break;

          case 6: // line: LBRACE sexp RBRACE
#line 145 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (echo && (yystack_[1].value.string) != NULL) {
              aprepro.lexer->LexerOutput((yystack_[1].value.string),
                                         strlen((yystack_[1].value.string)));
            }
          }
#line 661 "apr_parser.cc"
          break;

          case 7: // line: LBRACE aexp RBRACE
#line 149 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
          }
#line 667 "apr_parser.cc"
          break;

          case 8: // line: LBRACE RBRACE
#line 150 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
          }
#line 673 "apr_parser.cc"
          break;

          case 9: // line: error RBRACE
#line 151 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            yyerrok;
          }
#line 679 "apr_parser.cc"
          break;

          case 10: // bool: exp LT exp
#line 154 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) < (yystack_[0].value.val);
          }
#line 685 "apr_parser.cc"
          break;

          case 11: // bool: exp GT exp
#line 155 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) > (yystack_[0].value.val);
          }
#line 691 "apr_parser.cc"
          break;

          case 12: // bool: NOT exp
#line 156 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = !((yystack_[0].value.val));
          }
#line 697 "apr_parser.cc"
          break;

          case 13: // bool: exp LE exp
#line 157 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) <= (yystack_[0].value.val);
          }
#line 703 "apr_parser.cc"
          break;

          case 14: // bool: exp GE exp
#line 158 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) >= (yystack_[0].value.val);
          }
#line 709 "apr_parser.cc"
          break;

          case 15: // bool: exp EQ exp
#line 159 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) == (yystack_[0].value.val);
          }
#line 715 "apr_parser.cc"
          break;

          case 16: // bool: exp NE exp
#line 160 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) != (yystack_[0].value.val);
          }
#line 721 "apr_parser.cc"
          break;

          case 17: // bool: exp LOR exp
#line 161 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 727 "apr_parser.cc"
          break;

          case 18: // bool: exp LAND exp
#line 162 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 733 "apr_parser.cc"
          break;

          case 19: // bool: bool LOR bool
#line 163 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 739 "apr_parser.cc"
          break;

          case 20: // bool: bool LAND bool
#line 164 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 745 "apr_parser.cc"
          break;

          case 21: // bool: bool LOR exp
#line 165 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 751 "apr_parser.cc"
          break;

          case 22: // bool: bool LAND exp
#line 166 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 757 "apr_parser.cc"
          break;

          case 23: // bool: exp LOR bool
#line 167 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 763 "apr_parser.cc"
          break;

          case 24: // bool: exp LAND bool
#line 168 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 769 "apr_parser.cc"
          break;

          case 25: // bool: LPAR bool RPAR
#line 169 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          }
#line 775 "apr_parser.cc"
          break;

          case 26: // bool: UNDVAR LOR exp
#line 171 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = 0 || (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 781 "apr_parser.cc"
          break;

          case 27: // bool: UNDVAR LAND exp
#line 172 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = 0 && (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 787 "apr_parser.cc"
          break;

          case 28: // bool: exp LOR UNDVAR
#line 173 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 793 "apr_parser.cc"
          break;

          case 29: // bool: exp LAND UNDVAR
#line 174 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 799 "apr_parser.cc"
          break;

          case 30: // bool: bool LOR UNDVAR
#line 175 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 805 "apr_parser.cc"
          break;

          case 31: // bool: bool LAND UNDVAR
#line 176 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 811 "apr_parser.cc"
          break;

          case 32: // bool: UNDVAR LOR bool
#line 177 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = 0 || (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 817 "apr_parser.cc"
          break;

          case 33: // bool: UNDVAR LAND bool
#line 178 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = 0 && (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 823 "apr_parser.cc"
          break;

          case 34: // bool: sexp LT sexp
#line 181 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) < 0 ? 1 : 0);
          }
#line 829 "apr_parser.cc"
          break;

          case 35: // bool: sexp GT sexp
#line 182 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) > 0 ? 1 : 0);
          }
#line 835 "apr_parser.cc"
          break;

          case 36: // bool: sexp LE sexp
#line 183 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) <= 0 ? 1 : 0);
          }
#line 841 "apr_parser.cc"
          break;

          case 37: // bool: sexp GE sexp
#line 184 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) >= 0 ? 1 : 0);
          }
#line 847 "apr_parser.cc"
          break;

          case 38: // bool: sexp EQ sexp
#line 185 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) == 0 ? 1 : 0);
          }
#line 853 "apr_parser.cc"
          break;

          case 39: // bool: sexp NE sexp
#line 186 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) != 0 ? 1 : 0);
          }
#line 859 "apr_parser.cc"
          break;

          case 40: // bool: UNDVAR LT sexp
#line 188 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) < 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 865 "apr_parser.cc"
          break;

          case 41: // bool: UNDVAR GT sexp
#line 189 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) > 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 871 "apr_parser.cc"
          break;

          case 42: // bool: UNDVAR LE sexp
#line 190 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) <= 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 877 "apr_parser.cc"
          break;

          case 43: // bool: UNDVAR GE sexp
#line 191 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) >= 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 883 "apr_parser.cc"
          break;

          case 44: // bool: UNDVAR EQ sexp
#line 192 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) == 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 889 "apr_parser.cc"
          break;

          case 45: // bool: UNDVAR NE sexp
#line 193 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) != 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 895 "apr_parser.cc"
          break;

          case 46: // bool: sexp LT UNDVAR
#line 195 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") < 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 901 "apr_parser.cc"
          break;

          case 47: // bool: sexp GT UNDVAR
#line 196 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") > 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 907 "apr_parser.cc"
          break;

          case 48: // bool: sexp LE UNDVAR
#line 197 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") <= 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 913 "apr_parser.cc"
          break;

          case 49: // bool: sexp GE UNDVAR
#line 198 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") >= 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 919 "apr_parser.cc"
          break;

          case 50: // bool: sexp EQ UNDVAR
#line 199 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") == 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 925 "apr_parser.cc"
          break;

          case 51: // bool: sexp NE UNDVAR
#line 200 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") != 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 931 "apr_parser.cc"
          break;

          case 52: // bool: UNDVAR LT exp
#line 202 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = 0 < (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 937 "apr_parser.cc"
          break;

          case 53: // bool: UNDVAR GT exp
#line 203 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = 0 > (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 943 "apr_parser.cc"
          break;

          case 54: // bool: UNDVAR LE exp
#line 204 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = 0 <= (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 949 "apr_parser.cc"
          break;

          case 55: // bool: UNDVAR GE exp
#line 205 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = 0 >= (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 955 "apr_parser.cc"
          break;

          case 56: // bool: UNDVAR EQ exp
#line 206 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = 0 == (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 961 "apr_parser.cc"
          break;

          case 57: // bool: UNDVAR NE exp
#line 207 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = 0 != (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 967 "apr_parser.cc"
          break;

          case 58: // bool: exp LT UNDVAR
#line 209 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) < 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 973 "apr_parser.cc"
          break;

          case 59: // bool: exp GT UNDVAR
#line 210 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) > 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 979 "apr_parser.cc"
          break;

          case 60: // bool: exp LE UNDVAR
#line 211 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) <= 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 985 "apr_parser.cc"
          break;

          case 61: // bool: exp GE UNDVAR
#line 212 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) >= 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 991 "apr_parser.cc"
          break;

          case 62: // bool: exp EQ UNDVAR
#line 213 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) == 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 997 "apr_parser.cc"
          break;

          case 63: // bool: exp NE UNDVAR
#line 214 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) != 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1003 "apr_parser.cc"
          break;

          case 64: // aexp: AVAR
#line 216 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) = aprepro.make_array(*((yystack_[0].value.tptr)->value.avar));
          }
#line 1009 "apr_parser.cc"
          break;

          case 65: // aexp: AFNCT LPAR sexp RPAR
#line 217 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_c == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_c))((yystack_[1].value.string));
            else
              yyerrok;
          }
#line 1020 "apr_parser.cc"
          break;

          case 66: // aexp: AFNCT LPAR sexp COMMA exp RPAR
#line 223 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_cd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 1031 "apr_parser.cc"
          break;

          case 67: // aexp: AFNCT LPAR sexp COMMA sexp RPAR
#line 229 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_cc == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              yyerrok;
          }
#line 1042 "apr_parser.cc"
          break;

          case 68: // aexp: AFNCT LPAR exp COMMA exp COMMA exp RPAR
#line 235 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.arrfnct_ddd == NULL))
              (yylhs.value.arrval) = (*((yystack_[7].value.tptr)->value.arrfnct_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 1053 "apr_parser.cc"
          break;

          case 69: // aexp: AFNCT LPAR exp COMMA exp RPAR
#line 241 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_dd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 1064 "apr_parser.cc"
          break;

          case 70: // aexp: AFNCT LPAR exp RPAR
#line 247 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_d == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_d))((yystack_[1].value.val));
            else
              yyerrok;
          }
#line 1075 "apr_parser.cc"
          break;

          case 71: // aexp: AFNCT LPAR aexp RPAR
#line 253 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_a == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_a))((yystack_[1].value.arrval));
            else
              yyerrok;
          }
#line 1086 "apr_parser.cc"
          break;

          case 72: // aexp: SVAR EQUAL aexp
#line 259 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 1095 "apr_parser.cc"
          break;

          case 73: // aexp: VAR EQUAL aexp
#line 263 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 1104 "apr_parser.cc"
          break;

          case 74: // aexp: AVAR EQUAL aexp
#line 267 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) = (yystack_[0].value.arrval);
            aprepro.redefine_array((yystack_[2].value.tptr)->value.avar);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 1112 "apr_parser.cc"
          break;

          case 75: // aexp: UNDVAR EQUAL aexp
#line 270 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 1119 "apr_parser.cc"
          break;

          case 76: // aexp: aexp PLU aexp
#line 272 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
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
#line 1132 "apr_parser.cc"
          break;

          case 77: // aexp: SUB aexp
#line 280 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), -1.0);
          }
#line 1138 "apr_parser.cc"
          break;

          case 78: // aexp: aexp SUB aexp
#line 282 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
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
#line 1151 "apr_parser.cc"
          break;

          case 79: // aexp: aexp TIM exp
#line 290 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) = array_scale((yystack_[2].value.arrval), (yystack_[0].value.val));
          }
#line 1157 "apr_parser.cc"
          break;

          case 80: // aexp: aexp DIV exp
#line 291 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) =
                array_scale((yystack_[2].value.arrval), 1.0 / (yystack_[0].value.val));
          }
#line 1163 "apr_parser.cc"
          break;

          case 81: // aexp: exp TIM aexp
#line 292 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), (yystack_[2].value.val));
          }
#line 1169 "apr_parser.cc"
          break;

          case 82: // aexp: aexp TIM aexp
#line 293 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
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
#line 1182 "apr_parser.cc"
          break;

          case 83: // sexp: QSTRING
#line 302 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = (yystack_[0].value.string);
          }
#line 1188 "apr_parser.cc"
          break;

          case 84: // sexp: SVAR
#line 303 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          }
#line 1194 "apr_parser.cc"
          break;

          case 85: // sexp: IMMSVAR
#line 304 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          }
#line 1200 "apr_parser.cc"
          break;

          case 86: // sexp: UNDVAR EQUAL sexp
#line 305 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            set_type(aprepro, (yystack_[2].value.tptr), Parser::token::SVAR);
          }
#line 1207 "apr_parser.cc"
          break;

          case 87: // sexp: SVAR EQUAL sexp
#line 307 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          }
#line 1215 "apr_parser.cc"
          break;

          case 88: // sexp: VAR EQUAL sexp
#line 310 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);
          }
#line 1224 "apr_parser.cc"
          break;

          case 89: // sexp: AVAR EQUAL sexp
#line 314 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = (yystack_[0].value.string);
            aprepro.redefine_array((yystack_[2].value.tptr)->value.avar);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);
          }
#line 1234 "apr_parser.cc"
          break;

          case 90: // sexp: IMMSVAR EQUAL sexp
#line 319 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = (char *)(yystack_[2].value.tptr)->value.svar.c_str();
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1240 "apr_parser.cc"
          break;

          case 91: // sexp: IMMVAR EQUAL sexp
#line 320 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          }
#line 1246 "apr_parser.cc"
          break;

          case 92: // sexp: SFNCT LPAR sexp RPAR
#line 321 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_c == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_c))(
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1257 "apr_parser.cc"
          break;

          case 93: // sexp: SFNCT LPAR RPAR
#line 327 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.strfnct == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[2].value.tptr)->value.strfnct))();
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1268 "apr_parser.cc"
          break;

          case 94: // sexp: SFNCT LPAR exp RPAR
#line 333 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_d == NULL))
              (yylhs.value.string) =
                  (char *)(*((yystack_[3].value.tptr)->value.strfnct_d))((yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1279 "apr_parser.cc"
          break;

          case 95: // sexp: SFNCT LPAR aexp RPAR
#line 339 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_a == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_a))(
                  (yystack_[1].value.arrval));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1290 "apr_parser.cc"
          break;

          case 96: // sexp: sexp CONCAT sexp
#line 345 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            concat_string((yystack_[2].value.string), (yystack_[0].value.string),
                          &(yylhs.value.string));
          }
#line 1296 "apr_parser.cc"
          break;

          case 97: // sexp: SFNCT LPAR exp COMMA exp RPAR
#line 346 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_dd == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1307 "apr_parser.cc"
          break;

          case 98: // sexp: SFNCT LPAR exp COMMA sexp COMMA sexp COMMA sexp COMMA sexp RPAR
#line 352 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.strfnct_dcccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[11].value.tptr)->value.strfnct_dcccc))(
                  (yystack_[9].value.val), (yystack_[7].value.string), (yystack_[5].value.string),
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1318 "apr_parser.cc"
          break;

          case 99: // sexp: SFNCT LPAR exp COMMA sexp COMMA sexp RPAR
#line 358 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_dcc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_dcc))(
                  (yystack_[5].value.val), (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1329 "apr_parser.cc"
          break;

          case 100: // sexp: SFNCT LPAR sexp COMMA sexp COMMA sexp RPAR
#line 364 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_ccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1340 "apr_parser.cc"
          break;

          case 101: // sexp: SFNCT LPAR sexp COMMA sexp RPAR
#line 370 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_cc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1351 "apr_parser.cc"
          break;

          case 102: // sexp: bool QUEST sexp COLON sexp
#line 376 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.string) = ((yystack_[4].value.val)) ? ((yystack_[2].value.string))
                                                             : ((yystack_[0].value.string));
          }
#line 1357 "apr_parser.cc"
          break;

          case 103: // exp: NUM
#line 378 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          }
#line 1363 "apr_parser.cc"
          break;

          case 104: // exp: INC NUM
#line 379 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val) + 1;
          }
#line 1369 "apr_parser.cc"
          break;

          case 105: // exp: DEC NUM
#line 380 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val) - 1;
          }
#line 1375 "apr_parser.cc"
          break;

          case 106: // exp: VAR
#line 381 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          }
#line 1381 "apr_parser.cc"
          break;

          case 107: // exp: IMMVAR
#line 382 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          }
#line 1387 "apr_parser.cc"
          break;

          case 108: // exp: INC VAR
#line 383 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
          }
#line 1393 "apr_parser.cc"
          break;

          case 109: // exp: DEC VAR
#line 384 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
          }
#line 1399 "apr_parser.cc"
          break;

          case 110: // exp: VAR INC
#line 385 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
          }
#line 1405 "apr_parser.cc"
          break;

          case 111: // exp: VAR DEC
#line 386 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
          }
#line 1411 "apr_parser.cc"
          break;

          case 112: // exp: VAR EQUAL exp
#line 387 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          }
#line 1418 "apr_parser.cc"
          break;

          case 113: // exp: SVAR EQUAL exp
#line 389 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1426 "apr_parser.cc"
          break;

          case 114: // exp: AVAR EQUAL exp
#line 392 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val);
            aprepro.redefine_array((yystack_[2].value.tptr)->value.avar);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1436 "apr_parser.cc"
          break;

          case 115: // exp: VAR EQ_PLUS exp
#line 397 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1442 "apr_parser.cc"
          break;

          case 116: // exp: VAR EQ_MINUS exp
#line 398 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1448 "apr_parser.cc"
          break;

          case 117: // exp: VAR EQ_TIME exp
#line 399 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1454 "apr_parser.cc"
          break;

          case 118: // exp: VAR EQ_DIV exp
#line 400 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1460 "apr_parser.cc"
          break;

          case 119: // exp: VAR EQ_POW exp
#line 401 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            SEAMS::math_error(aprepro, "Power");
          }
#line 1470 "apr_parser.cc"
          break;

          case 120: // exp: INC IMMVAR
#line 406 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          }
#line 1476 "apr_parser.cc"
          break;

          case 121: // exp: DEC IMMVAR
#line 407 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          }
#line 1482 "apr_parser.cc"
          break;

          case 122: // exp: IMMVAR INC
#line 408 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          }
#line 1488 "apr_parser.cc"
          break;

          case 123: // exp: IMMVAR DEC
#line 409 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          }
#line 1494 "apr_parser.cc"
          break;

          case 124: // exp: IMMVAR EQUAL exp
#line 410 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1500 "apr_parser.cc"
          break;

          case 125: // exp: IMMSVAR EQUAL exp
#line 411 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          }
#line 1506 "apr_parser.cc"
          break;

          case 126: // exp: IMMVAR EQ_PLUS exp
#line 412 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1512 "apr_parser.cc"
          break;

          case 127: // exp: IMMVAR EQ_MINUS exp
#line 413 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1518 "apr_parser.cc"
          break;

          case 128: // exp: IMMVAR EQ_TIME exp
#line 414 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1524 "apr_parser.cc"
          break;

          case 129: // exp: IMMVAR EQ_DIV exp
#line 415 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1530 "apr_parser.cc"
          break;

          case 130: // exp: IMMVAR EQ_POW exp
#line 416 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1536 "apr_parser.cc"
          break;

          case 131: // exp: UNDVAR
#line 418 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1543 "apr_parser.cc"
          break;

          case 132: // exp: INC UNDVAR
#line 420 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1551 "apr_parser.cc"
          break;

          case 133: // exp: DEC UNDVAR
#line 423 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1559 "apr_parser.cc"
          break;

          case 134: // exp: UNDVAR INC
#line 426 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          }
#line 1567 "apr_parser.cc"
          break;

          case 135: // exp: UNDVAR DEC
#line 429 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          }
#line 1575 "apr_parser.cc"
          break;

          case 136: // exp: UNDVAR EQUAL exp
#line 432 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1582 "apr_parser.cc"
          break;

          case 137: // exp: UNDVAR EQ_PLUS exp
#line 434 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1590 "apr_parser.cc"
          break;

          case 138: // exp: UNDVAR EQ_MINUS exp
#line 437 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1598 "apr_parser.cc"
          break;

          case 139: // exp: UNDVAR EQ_TIME exp
#line 440 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1606 "apr_parser.cc"
          break;

          case 140: // exp: UNDVAR EQ_DIV exp
#line 443 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1614 "apr_parser.cc"
          break;

          case 141: // exp: UNDVAR EQ_POW exp
#line 446 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            SEAMS::math_error(aprepro, "Power");
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1625 "apr_parser.cc"
          break;

          case 142: // exp: FNCT LPAR RPAR
#line 453 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.fnctptr == NULL))
              (yylhs.value.val) = (*((yystack_[2].value.tptr)->value.fnctptr))();
            else
              (yylhs.value.val) = 0.0;
          }
#line 1636 "apr_parser.cc"
          break;

          case 143: // exp: FNCT LPAR exp RPAR
#line 460 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_d == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_d))((yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1647 "apr_parser.cc"
          break;

          case 144: // exp: FNCT LPAR sexp RPAR
#line 467 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_c == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_c))((yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1658 "apr_parser.cc"
          break;

          case 145: // exp: FNCT LPAR aexp RPAR
#line 474 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_a == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_a))((yystack_[1].value.arrval));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1669 "apr_parser.cc"
          break;

          case 146: // exp: FNCT LPAR sexp COMMA exp RPAR
#line 481 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1680 "apr_parser.cc"
          break;

          case 147: // exp: FNCT LPAR exp COMMA sexp RPAR
#line 488 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dc))(
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1691 "apr_parser.cc"
          break;

          case 148: // exp: FNCT LPAR sexp COMMA sexp RPAR
#line 495 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1702 "apr_parser.cc"
          break;

          case 149: // exp: FNCT LPAR sexp COMMA sexp COMMA sexp RPAR
#line 502 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccc == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              yyerrok;
          }
#line 1713 "apr_parser.cc"
          break;

          case 150: // exp: FNCT LPAR exp COMMA exp RPAR
#line 509 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1724 "apr_parser.cc"
          break;

          case 151: // exp: FNCT LPAR exp COMMA exp COMMA exp RPAR
#line 515 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ddd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1735 "apr_parser.cc"
          break;

          case 152: // exp: FNCT LPAR sexp COMMA sexp COMMA exp RPAR
#line 521 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccd))(
                  (yystack_[5].value.string), (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1746 "apr_parser.cc"
          break;

          case 153: // exp: FNCT LPAR exp COMMA exp SEMI exp COMMA exp RPAR
#line 527 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1757 "apr_parser.cc"
          break;

          case 154: // exp: FNCT LPAR exp COMMA exp COMMA exp COMMA exp RPAR
#line 533 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1768 "apr_parser.cc"
          break;

          case 155: // exp: FNCT LPAR exp COMMA exp COMMA exp COMMA exp COMMA sexp RPAR
#line 539 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.fnctptr_ddddc == NULL))
              (yylhs.value.val) = (*((yystack_[11].value.tptr)->value.fnctptr_ddddc))(
                  (yystack_[9].value.val), (yystack_[7].value.val), (yystack_[5].value.val),
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1779 "apr_parser.cc"
          break;

          case 156: // exp: FNCT LPAR exp COMMA exp COMMA exp COMMA exp COMMA exp COMMA exp RPAR
#line 545 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if (arg_check((yystack_[13].value.tptr),
                          (yystack_[13].value.tptr)->value.fnctptr_dddddd == NULL))
              (yylhs.value.val) = (*((yystack_[13].value.tptr)->value.fnctptr_dddddd))(
                  (yystack_[11].value.val), (yystack_[9].value.val), (yystack_[7].value.val),
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1790 "apr_parser.cc"
          break;

          case 157: // exp: exp PLU exp
#line 551 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) + (yystack_[0].value.val);
          }
#line 1796 "apr_parser.cc"
          break;

          case 158: // exp: exp SUB exp
#line 552 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) - (yystack_[0].value.val);
          }
#line 1802 "apr_parser.cc"
          break;

          case 159: // exp: exp TIM exp
#line 553 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) * (yystack_[0].value.val);
          }
#line 1808 "apr_parser.cc"
          break;

          case 160: // exp: exp DIV exp
#line 554 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if ((yystack_[0].value.val) == 0.) {
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (yystack_[2].value.val) / (yystack_[0].value.val);
          }
#line 1820 "apr_parser.cc"
          break;

          case 161: // exp: exp MOD exp
#line 561 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            if ((yystack_[0].value.val) == 0.) {
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (int)(yystack_[2].value.val) % (int)(yystack_[0].value.val);
          }
#line 1832 "apr_parser.cc"
          break;

          case 162: // exp: SUB exp
#line 568 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = -(yystack_[0].value.val);
          }
#line 1838 "apr_parser.cc"
          break;

          case 163: // exp: PLU exp
#line 569 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          }
#line 1844 "apr_parser.cc"
          break;

          case 164: // exp: exp POW exp
#line 570 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            reset_error();
            (yylhs.value.val) = std::pow((yystack_[2].value.val), (yystack_[0].value.val));
            SEAMS::math_error(aprepro, "Power");
          }
#line 1852 "apr_parser.cc"
          break;

          case 165: // exp: LPAR exp RPAR
#line 573 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          }
#line 1858 "apr_parser.cc"
          break;

          case 166: // exp: LBRACK exp RBRACK
#line 574 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            reset_error();
            (yylhs.value.val) =
                (double)((yystack_[1].value.val) < 0 ? -floor(-((yystack_[1].value.val)))
                                                     : floor((yystack_[1].value.val)));
            SEAMS::math_error(aprepro, "floor (int)");
          }
#line 1866 "apr_parser.cc"
          break;

          case 167: // exp: bool
#line 577 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[0].value.val)) ? 1 : 0;
          }
#line 1872 "apr_parser.cc"
          break;

          case 168: // exp: bool QUEST exp COLON exp
#line 578 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                ((yystack_[4].value.val)) ? ((yystack_[2].value.val)) : ((yystack_[0].value.val));
          }
#line 1878 "apr_parser.cc"
          break;

          case 169: // exp: AVAR LBRACK exp RBRACK
#line 579 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) =
                array_value((yystack_[3].value.tptr)->value.avar, (yystack_[1].value.val), 0);
          }
#line 1884 "apr_parser.cc"
          break;

          case 170: // exp: AVAR LBRACK exp COMMA exp RBRACK
#line 580 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
          {
            (yylhs.value.val) = array_value((yystack_[5].value.tptr)->value.avar,
                                            (yystack_[3].value.val), (yystack_[1].value.val));
          }
#line 1890 "apr_parser.cc"
          break;

          case 171: // exp: AVAR LBRACK exp RBRACK EQUAL exp
#line 582 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
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
#line 1916 "apr_parser.cc"
          break;

          case 172: // exp: AVAR LBRACK exp COMMA exp RBRACK EQUAL exp
#line 604 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"
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
#line 1940 "apr_parser.cc"
          break;

#line 1944 "apr_parser.cc"

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

  const signed char Parser::yypact_ninf_ = -25;

  const signed char Parser::yytable_ninf_ = -1;

  const short Parser::yypact_[] = {
      -25,  22,   -25,  -18,  209,  -25,  -25,  -25,  -25,  -25,  1629, 163,  -8,   211,  -5,
      154,  6,    51,   58,   471,  471,  -25,  471,  426,  471,  -2,   379,  71,   80,   1560,
      1665, 426,  471,  471,  471,  471,  471,  471,  471,  471,  471,  471,  471,  471,  471,
      -25,  -25,  426,  471,  471,  471,  471,  471,  -25,  -25,  426,  471,  471,  471,  471,
      471,  471,  -25,  -25,  471,  471,  426,  354,  411,  426,  1647, 326,  20,   181,  471,
      120,  343,  1348, 1580, 14,   -25,  14,   14,   -25,  -25,  -25,  -25,  -25,  -25,  -25,
      -25,  471,  486,  531,  -25,  426,  426,  471,  426,  -25,  546,  591,  606,  651,  666,
      711,  471,  -25,  726,  771,  786,  831,  846,  891,  906,  951,  471,  471,  471,  426,
      471,  471,  44,   343,  1698, 1714, 1714, 1714, 1714, 1714, 55,   1745, -25,  1759, 50,
      367,  50,   367,  50,   367,  50,   367,  50,   367,  50,   367,  44,   343,  1698, 1714,
      1714, 1714, 1714, 1714, 44,   343,  1698, 343,  1714, 1714, 1714, 1714, 1714, 1714, 343,
      1714, 1083, 44,   343,  1698, -25,  87,   254,  1113, -25,  167,  264,  1143, 279,  292,
      1173, 471,  471,  471,  471,  14,   -25,  -25,  471,  -25,  -24,  1682, 1647, 55,   1745,
      1647, -25,  1759, 9,    1698, 9,    1730, -25,  1730, 1647, 50,   1714, 1647, 50,   1647,
      50,   1647, 50,   1647, 50,   1647, 50,   -25,  1647, 55,   1745, 1647, -25,  1759, 1647,
      367,  1647, 367,  1647, 367,  1647, 367,  1647, 367,  1647, 367,  27,   27,   14,   -25,
      14,   14,   14,   471,  101,  -25,  471,  -25,  471,  -25,  -25,  471,  -25,  471,  -25,
      -25,  471,  -25,  471,  -25,  1714, 1714, 1714, 1714, 14,   471,  471,  1605, 471,  967,
      1375, 42,   1054, 977,  937,  1402, 108,  1429, 1203, 343,  1714, 117,  1714, 471,  -25,
      -25,  -25,  471,  -25,  471,  471,  -25,  471,  -25,  -25,  -25,  471,  -25,  471,  133,
      1456, 1233, 1292, 298,  1002, 1483, 1714, -25,  -25,  471,  -25,  471,  -25,  471,  -25,
      -25,  1263, 1510, 994,  471,  -25,  -25,  471,  1024, 1321, 1030, -25,  471,  -25,  1537,
      -25};

  const unsigned char Parser::yydefact_[] = {
      2,   0,   1,   0,   0,   4,   3,   9,   103, 83,  131, 106, 84,  107, 85,  64,  0,   0,   0,
      0,   0,   8,   0,   0,   0,   0,   0,   167, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   134, 135, 0,   0,   0,   0,   0,   0,   110, 111, 0,   0,
      0,   0,   0,   0,   0,   122, 123, 0,   0,   0,   0,   0,   0,   131, 106, 84,  0,   0,   167,
      0,   0,   0,   163, 77,  162, 12,  104, 132, 108, 120, 105, 133, 109, 121, 0,   0,   0,   7,
      0,   0,   0,   0,   6,   0,   0,   0,   0,   0,   0,   0,   5,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   75,  86,  136, 137, 138, 139, 140, 141, 32,  26,  33,
      27,  40,  52,  41,  53,  42,  54,  43,  55,  44,  56,  45,  57,  73,  88,  112, 115, 116, 117,
      118, 119, 72,  87,  113, 91,  124, 126, 127, 128, 129, 130, 90,  125, 0,   74,  89,  114, 142,
      0,   0,   0,   93,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   162, 25,  165, 0,   166,
      0,   0,   30,  19,  21,  31,  20,  22,  76,  0,   78,  80,  82,  79,  46,  34,  0,   47,  35,
      48,  36,  49,  37,  50,  38,  51,  39,  96,  28,  23,  17,  29,  24,  18,  58,  10,  59,  11,
      60,  13,  61,  14,  62,  15,  63,  16,  157, 158, 160, 81,  159, 161, 164, 0,   169, 145, 0,
      144, 0,   143, 95,  0,   92,  0,   94,  71,  0,   65,  0,   70,  136, 112, 113, 114, 159, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   102, 168, 170, 171, 0,   148,
      146, 147, 0,   150, 0,   0,   101, 0,   97,  67,  66,  0,   69,  0,   0,   0,   0,   0,   0,
      0,   0,   172, 149, 152, 0,   151, 0,   100, 0,   99,  68,  0,   0,   0,   0,   154, 153, 0,
      0,   0,   0,   155, 0,   98,  0,   156};

  const signed char Parser::yypgoto_[] = {-25, -25, -25, -12, 106, 91, -4};

  const unsigned char Parser::yydefgoto_[] = {0, 1, 6, 27, 28, 76, 206};

  const short Parser::yytable_[] = {
      30,  83,  7,   84,  85,  265, 86,  75,  100, 101, 102, 103, 104, 105, 55,  77,  78,  64,
      79,  81,  82,  67,  2,   3,   106, 130, 132, 124, 125, 126, 127, 128, 129, 131, 133, 135,
      137, 139, 141, 143, 145, 4,   183, 148, 149, 150, 151, 152, 153, 97,  98,  156, 158, 159,
      160, 161, 162, 163, 286, 121, 165, 166, 169, 173, 177, 180, 68,  118, 188, 120, 185, 5,
      121, 69,  100, 101, 102, 103, 104, 105, 193, 196, 95,  96,  97,  98,  93,  191, 194, 197,
      106, 199, 199, 201, 203, 29,  219, 222, 106, 91,  94,  92,  93,  245, 220, 223, 225, 227,
      229, 231, 233, 235, 236, 237, 238, 240, 241, 242, 95,  96,  97,  98,  123, 268, 294, 95,
      96,  97,  98,  80,  134, 136, 138, 140, 142, 144, 186, 122, 147, 298, 100, 101, 102, 103,
      104, 105, 155, 157, 91,  307, 92,  93,  0,   146, 0,   164, 106, 168, 172, 176, 179, 154,
      0,   0,   0,   100, 101, 102, 103, 104, 105, 65,  167, 171, 175, 178, 66,  260, 261, 262,
      263, 106, 190, 250, 264, 47,  48,  49,  50,  51,  52,  205, 208, 210, 212, 214, 216, 217,
      65,  0,   0,   198, 200, 184, 202, 95,  96,  97,  98,  53,  54,  0,   8,   9,   10,  11,
      12,  13,  14,  15,  16,  17,  18,  0,   19,  239, 20,  0,   0,   21,  0,   0,   0,   56,
      57,  58,  59,  60,  61,  267, 0,   0,   270, 0,   272, 0,   0,   22,  23,  275, 0,   0,
      277, 24,  278, 25,  26,  62,  63,  0,   0,   0,   280, 0,   282, 0,   0,   0,   246, 0,
      247, 0,   123, 147, 155, 168, 0,   0,   251, 300, 252, 0,   0,   301, 0,   302, 100, 101,
      102, 103, 104, 105, 305, 0,   306, 255, 100, 101, 102, 103, 104, 105, 106, 0,   0,   316,
      256, 317, 257, 0,   0,   0,   106, 0,   312, 324, 0,   95,  96,  97,  98,  0,   0,   329,
      100, 101, 102, 103, 104, 105, 100, 101, 102, 103, 104, 105, 0,   269, 0,   271, 106, 0,
      273, 0,   274, 0,   106, 276, 182, 48,  49,  50,  51,  52,  0,   0,   279, 8,   9,   10,
      11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  170, 20,  53,  54,  299, 100, 101, 102,
      103, 104, 105, 303, 87,  304, 88,  89,  0,   90,  0,   0,   0,   106, 22,  23,  0,   0,
      0,   0,   24,  0,   25,  26,  0,   0,   318, 116, 117, 118, 188, 120, 323, 0,   121, 325,
      8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  174, 20,  8,   9,   10,
      11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  0,   20,  0,   0,   0,   0,   0,   22,
      23,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,   0,   0,   22,  23,  0,   0,
      0,   0,   24,  0,   25,  26,  8,   9,   70,  71,  72,  13,  14,  73,  16,  17,  0,   0,
      19,  0,   20,  8,   9,   192, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,
      0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,
      0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  8,   9,   195, 71,  72,  13,
      14,  73,  16,  17,  0,   0,   19,  0,   20,  8,   9,   204, 71,  72,  13,  14,  73,  16,
      17,  0,   0,   19,  0,   20,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,
      0,   25,  26,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,
      8,   9,   207, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  8,   9,   209,
      71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  0,   0,   0,   0,   0,   22,
      74,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,   0,   0,   22,  74,  0,   0,
      0,   0,   24,  0,   25,  26,  8,   9,   211, 71,  72,  13,  14,  73,  16,  17,  0,   0,
      19,  0,   20,  8,   9,   213, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,
      0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,
      0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  8,   9,   215, 71,  72,  13,
      14,  73,  16,  17,  0,   0,   19,  0,   20,  8,   9,   218, 71,  72,  13,  14,  73,  16,
      17,  0,   0,   19,  0,   20,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,
      0,   25,  26,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,
      8,   9,   221, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  8,   9,   224,
      71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  0,   0,   0,   0,   0,   22,
      74,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,   0,   0,   22,  74,  0,   0,
      0,   0,   24,  0,   25,  26,  8,   9,   226, 71,  72,  13,  14,  73,  16,  17,  0,   0,
      19,  0,   20,  8,   9,   228, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,
      0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,
      0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  8,   9,   230, 71,  72,  13,
      14,  73,  16,  17,  0,   0,   19,  0,   20,  8,   9,   232, 71,  72,  13,  14,  73,  16,
      17,  0,   0,   19,  0,   20,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,
      0,   25,  26,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  292, 25,  26,
      8,   9,   234, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  100, 101, 102,
      103, 104, 105, 0,   0,   0,   0,   0,   0,   283, 0,   284, 0,   106, 0,   0,   0,   22,
      74,  290, 0,   291, 0,   24,  0,   25,  26,  100, 101, 102, 103, 104, 105, 0,   0,   0,
      322, 100, 101, 102, 103, 104, 105, 106, 313, 0,   314, 0,   0,   0,   0,   0,   0,   106,
      100, 101, 102, 103, 104, 105, 0,   0,   100, 101, 102, 103, 104, 105, 326, 0,   106, 0,
      0,   0,   328, 0,   0,   0,   106, 0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105,
      100, 101, 102, 103, 104, 105, 287, 0,   288, 0,   106, 0,   0,   289, 0,   0,   106, 0,
      0,   0,   0,   0,   108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 243,
      0,   121, 0,   244, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110,
      111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 0,   248, 121, 249, 0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110, 111, 112, 113, 114, 115, 116,
      117, 118, 119, 120, 0,   253, 121, 254, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 0,   258,
      121, 259, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110,
      111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 0,   296, 121, 297, 0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110, 111, 112, 113, 114, 115, 116,
      117, 118, 188, 120, 0,   309, 121, 310, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 0,   319,
      121, 320, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110,
      111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 311, 0,   121, 0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110, 111, 112, 113, 114, 115, 116, 117,
      118, 188, 120, 327, 0,   121, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 187, 0,   121, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110, 111, 112, 113, 114, 115,
      116, 117, 118, 188, 120, 285, 0,   121, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 293, 0,   121, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110, 111, 112, 113, 114, 115,
      116, 117, 118, 188, 120, 295, 0,   121, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 308, 0,   121, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110, 111, 112, 113, 114, 115,
      116, 117, 118, 188, 120, 315, 0,   121, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 321, 0,   121, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110, 111, 112, 113, 114, 115,
      116, 117, 118, 188, 120, 330, 0,   121, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 99,  0,   121, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 189, 0,   0,   0,
      0,   0,   0,   0,   0,   0,   106, 0,   108, 109, 110, 111, 112, 113, 114, 115, 116, 117,
      118, 188, 120, 281, 0,   121, 0,   0,   0,   0,   0,   0,   0,   0,   0,   108, 109, 110,
      111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 0,   0,   121, 31,  32,  33,  34,  35,
      36,  0,   0,   37,  38,  39,  40,  41,  42,  43,  44,  0,   0,   181, 32,  33,  34,  35,
      36,  45,  46,  37,  38,  39,  40,  41,  42,  43,  44,  107, 0,   0,   0,   0,   0,   0,
      0,   45,  46,  108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 0,   0,
      121, 266, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 0,   0,   121,
      108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 0,   0,   121, 108, 109,
      110, 111, 112, 113, 114, 115, 116, 117, 118, 188, 120, 0,   0,   121, 108, 109, 110, 111,
      112, 113, 114, 115, 0,   0,   0,   0,   120, 0,   0,   121, 109, 110, 111, 112, 113, 114,
      115, 116, 117, 118, 188, 120, 0,   0,   121, 110, 111, 112, 113, 114, 115, 116, 117, 118,
      188, 120, 0,   0,   121};

  const short Parser::yycheck_[] = {
      4,   3,   20,  5,   6,   29,  8,   19,  32,  33,  34,  35,  36,  37,  22,  19,  20,  22,  22,
      23,  24,  15,  0,   1,   48,  37,  38,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
      42,  43,  44,  19,  22,  47,  48,  49,  50,  51,  52,  40,  41,  55,  56,  57,  58,  59,  60,
      61,  16,  45,  64,  65,  66,  67,  68,  69,  15,  40,  41,  42,  74,  49,  45,  15,  32,  33,
      34,  35,  36,  37,  92,  93,  38,  39,  40,  41,  31,  91,  92,  93,  48,  95,  96,  97,  98,
      4,   108, 109, 48,  28,  20,  30,  31,  16,  108, 109, 110, 111, 112, 113, 114, 115, 116, 117,
      118, 119, 120, 121, 38,  39,  40,  41,  31,  22,  16,  38,  39,  40,  41,  23,  39,  40,  41,
      42,  43,  44,  16,  31,  47,  22,  32,  33,  34,  35,  36,  37,  55,  56,  28,  16,  30,  31,
      -1,  47,  -1,  64,  48,  66,  67,  68,  69,  55,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,
      17,  66,  67,  68,  69,  22,  181, 182, 183, 184, 48,  91,  16,  188, 22,  23,  24,  25,  26,
      27,  100, 101, 102, 103, 104, 105, 106, 17,  -1,  -1,  95,  96,  22,  98,  38,  39,  40,  41,
      46,  47,  -1,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,  119, 17,  -1,
      -1,  20,  -1,  -1,  -1,  22,  23,  24,  25,  26,  27,  243, -1,  -1,  246, -1,  248, -1,  -1,
      38,  39,  253, -1,  -1,  256, 44,  258, 46,  47,  46,  47,  -1,  -1,  -1,  266, -1,  268, -1,
      -1,  -1,  14,  -1,  16,  -1,  181, 182, 183, 184, -1,  -1,  14,  283, 16,  -1,  -1,  287, -1,
      289, 32,  33,  34,  35,  36,  37,  296, -1,  298, 16,  32,  33,  34,  35,  36,  37,  48,  -1,
      -1,  309, 14,  311, 16,  -1,  -1,  -1,  48,  -1,  16,  319, -1,  38,  39,  40,  41,  -1,  -1,
      327, 32,  33,  34,  35,  36,  37,  32,  33,  34,  35,  36,  37,  -1,  246, -1,  248, 48,  -1,
      251, -1,  253, -1,  48,  256, 22,  23,  24,  25,  26,  27,  -1,  -1,  265, 3,   4,   5,   6,
      7,   8,   9,   10,  11,  12,  13,  -1,  15,  16,  17,  46,  47,  283, 32,  33,  34,  35,  36,
      37,  290, 3,   292, 5,   6,   -1,  8,   -1,  -1,  -1,  48,  38,  39,  -1,  -1,  -1,  -1,  44,
      -1,  46,  47,  -1,  -1,  313, 38,  39,  40,  41,  42,  319, -1,  45,  322, 3,   4,   5,   6,
      7,   8,   9,   10,  11,  12,  13,  -1,  15,  16,  17,  3,   4,   5,   6,   7,   8,   9,   10,
      11,  12,  13,  -1,  15,  -1,  17,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,
      -1,  46,  47,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  3,
      4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  3,   4,   5,   6,   7,
      8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,
      -1,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,
      46,  47,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  3,   4,
      5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  -1,  -1,  -1,  -1,  -1,  38,
      39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,
      -1,  44,  -1,  46,  47,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,
      17,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  -1,  -1,  -1,
      -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,  -1,  -1,  38,  39,
      -1,  -1,  -1,  -1,  44,  -1,  46,  47,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,
      -1,  15,  -1,  17,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,
      -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,  -1,
      -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  3,   4,   5,   6,   7,   8,   9,   10,
      11,  12,  -1,  -1,  15,  -1,  17,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,
      15,  -1,  17,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,
      -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  3,   4,   5,   6,   7,
      8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  3,   4,   5,   6,   7,   8,   9,   10,  11,
      12,  -1,  -1,  15,  -1,  17,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,
      46,  47,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  3,   4,
      5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  3,   4,   5,   6,   7,   8,
      9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,
      -1,  44,  -1,  46,  47,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,
      47,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  3,   4,   5,
      6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  -1,  -1,  -1,  -1,  -1,  38,  39,
      -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,
      44,  14,  46,  47,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,
      32,  33,  34,  35,  36,  37,  -1,  -1,  -1,  -1,  -1,  -1,  14,  -1,  16,  -1,  48,  -1,  -1,
      -1,  38,  39,  14,  -1,  16,  -1,  44,  -1,  46,  47,  32,  33,  34,  35,  36,  37,  -1,  -1,
      -1,  14,  32,  33,  34,  35,  36,  37,  48,  14,  -1,  16,  -1,  -1,  -1,  -1,  -1,  -1,  48,
      32,  33,  34,  35,  36,  37,  -1,  -1,  32,  33,  34,  35,  36,  37,  16,  -1,  48,  -1,  -1,
      -1,  16,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  32,  33,
      34,  35,  36,  37,  14,  -1,  16,  -1,  48,  -1,  -1,  21,  -1,  -1,  48,  -1,  -1,  -1,  -1,
      -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  14,  -1,  45,  -1,  18,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,
      38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,
      16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,
      35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
      -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,
      32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
      40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  14,  -1,  45,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,
      38,  39,  40,  41,  42,  14,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,
      39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
      42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,
      34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,
      45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,
      37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
      40,  41,  42,  20,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,
      36,  37,  18,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  48,  -1,  30,  31,  32,  33,  34,
      35,  36,  37,  38,  39,  40,  41,  42,  18,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  22,  23,
      24,  25,  26,  27,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  -1,  -1,  22,  23,  24,
      25,  26,  27,  46,  47,  30,  31,  32,  33,  34,  35,  36,  37,  20,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  46,  47,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,
      45,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  30,
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  30,  31,  32,  33,
      34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  30,  31,  32,  33,  34,  35,  36,
      37,  -1,  -1,  -1,  -1,  42,  -1,  -1,  45,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
      41,  42,  -1,  -1,  45,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45};

  const signed char Parser::yystos_[] = {
      0,  51, 0,  1,  19, 49, 52, 20, 3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 15, 17, 20, 38,
      39, 44, 46, 47, 53, 54, 55, 56, 22, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 35, 36, 37, 46,
      47, 22, 23, 24, 25, 26, 27, 46, 47, 22, 22, 23, 24, 25, 26, 27, 46, 47, 22, 17, 22, 15, 15,
      15, 5,  6,  7,  10, 39, 53, 55, 56, 56, 56, 54, 56, 56, 3,  5,  6,  8,  3,  5,  6,  8,  28,
      30, 31, 20, 38, 39, 40, 41, 20, 32, 33, 34, 35, 36, 37, 48, 20, 30, 31, 32, 33, 34, 35, 36,
      37, 38, 39, 40, 41, 42, 45, 54, 55, 56, 56, 56, 56, 56, 56, 53, 56, 53, 56, 55, 56, 55, 56,
      55, 56, 55, 56, 55, 56, 55, 56, 54, 55, 56, 56, 56, 56, 56, 56, 54, 55, 56, 55, 56, 56, 56,
      56, 56, 56, 55, 56, 56, 54, 55, 56, 16, 54, 55, 56, 16, 54, 55, 56, 54, 55, 56, 22, 22, 22,
      22, 56, 16, 16, 41, 18, 55, 56, 5,  53, 56, 5,  53, 56, 54, 56, 54, 56, 54, 56, 5,  55, 56,
      5,  55, 5,  55, 5,  55, 5,  55, 5,  55, 55, 5,  53, 56, 5,  53, 56, 5,  56, 5,  56, 5,  56,
      5,  56, 5,  56, 5,  56, 56, 56, 56, 54, 56, 56, 56, 14, 18, 16, 14, 16, 14, 16, 16, 14, 16,
      14, 16, 16, 14, 16, 14, 16, 56, 56, 56, 56, 56, 29, 29, 56, 22, 55, 56, 55, 56, 55, 55, 56,
      55, 56, 56, 55, 56, 18, 56, 14, 16, 16, 16, 14, 16, 21, 14, 16, 14, 16, 16, 16, 14, 16, 22,
      55, 56, 56, 56, 55, 55, 56, 56, 16, 16, 14, 16, 14, 16, 14, 16, 16, 56, 56, 55, 14, 16, 16,
      14, 55, 56, 55, 16, 14, 16, 56, 16};

  const signed char Parser::yyr1_[] = {
      0,  50, 51, 51, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 54, 54,
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55,
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56};

  const signed char Parser::yyr2_[] = {
      0, 2, 0, 2, 1, 3, 3, 3, 2,  2,  3,  3,  2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3, 3,  3,  3,  3,  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 1, 4, 6,  6,  8,  6,  4, 4, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 1, 1, 1, 3,
      3, 3, 3, 3, 3, 4, 3, 4, 4,  3,  6,  12, 8, 8, 6, 5, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
      3, 3, 3, 3, 2, 2, 2, 2, 3,  3,  3,  3,  3, 3, 3, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4,
      4, 6, 6, 6, 8, 6, 8, 8, 10, 10, 12, 14, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 1, 5, 4, 6, 6, 8};

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
      0,   127, 127, 128, 131, 132, 145, 149, 150, 151, 154, 155, 156, 157, 158, 159, 160, 161,
      162, 163, 164, 165, 166, 167, 168, 169, 171, 172, 173, 174, 175, 176, 177, 178, 181, 182,
      183, 184, 185, 186, 188, 189, 190, 191, 192, 193, 195, 196, 197, 198, 199, 200, 202, 203,
      204, 205, 206, 207, 209, 210, 211, 212, 213, 214, 216, 217, 223, 229, 235, 241, 247, 253,
      259, 263, 267, 270, 272, 280, 282, 290, 291, 292, 293, 302, 303, 304, 305, 307, 310, 314,
      319, 320, 321, 327, 333, 339, 345, 346, 352, 358, 364, 370, 376, 378, 379, 380, 381, 382,
      383, 384, 385, 386, 387, 389, 392, 397, 398, 399, 400, 401, 406, 407, 408, 409, 410, 411,
      412, 413, 414, 415, 416, 418, 420, 423, 426, 429, 432, 434, 437, 440, 443, 446, 453, 460,
      467, 474, 481, 488, 495, 502, 509, 515, 521, 527, 533, 539, 545, 551, 552, 553, 554, 561,
      568, 569, 570, 573, 574, 577, 578, 579, 580, 581, 603};

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
#line 2966 "apr_parser.cc"

#line 626 "/Users/gdsjaar/src/seacas/packages/seacas/libraries/aprepro_lib/aprepro.yy"

void SEAMS::Parser::error(const std::string &m) { aprepro.error(m); }
