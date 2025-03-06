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

#include "aprepro_parser.h"

// Second part of user prologue.

#include "apr_scanner.h"
#include "aprepro.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex

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
          {
            if (echo)
              aprepro.lexer->LexerOutput("\n", 1);
          } break;

          case 5: // line: LBRACE exp RBRACE
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
          } break;

          case 6: // line: LBRACE sexp RBRACE
          {
            if (echo && (yystack_[1].value.string) != NULL) {
              aprepro.lexer->LexerOutput((yystack_[1].value.string),
                                         strlen((yystack_[1].value.string)));
            }
          } break;

          case 7: // line: LBRACE aexp RBRACE
          {
          } break;

          case 8: // line: LBRACE RBRACE
          {
          } break;

          case 9: // line: error RBRACE
          {
            yyerrok;
          } break;

          case 10: // bool: exp LT exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) < (yystack_[0].value.val);
          } break;

          case 11: // bool: exp GT exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) > (yystack_[0].value.val);
          } break;

          case 12: // bool: NOT exp
          {
            (yylhs.value.val) = !((yystack_[0].value.val));
          } break;

          case 13: // bool: exp LE exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) <= (yystack_[0].value.val);
          } break;

          case 14: // bool: exp GE exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) >= (yystack_[0].value.val);
          } break;

          case 15: // bool: exp EQ exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) == (yystack_[0].value.val);
          } break;

          case 16: // bool: exp NE exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) != (yystack_[0].value.val);
          } break;

          case 17: // bool: exp LOR exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          } break;

          case 18: // bool: exp LAND exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          } break;

          case 19: // bool: bool LOR bool
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          } break;

          case 20: // bool: bool LAND bool
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          } break;

          case 21: // bool: bool LOR exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          } break;

          case 22: // bool: bool LAND exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          } break;

          case 23: // bool: exp LOR bool
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          } break;

          case 24: // bool: exp LAND bool
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          } break;

          case 25: // bool: LPAR bool RPAR
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          } break;

          case 26: // bool: UNDVAR LOR exp
          {
            (yylhs.value.val) = 0 || (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 27: // bool: UNDVAR LAND exp
          {
            (yylhs.value.val) = 0 && (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 28: // bool: exp LOR UNDVAR
          {
            (yylhs.value.val) = (yystack_[2].value.val) || 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 29: // bool: exp LAND UNDVAR
          {
            (yylhs.value.val) = (yystack_[2].value.val) && 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 30: // bool: bool LOR UNDVAR
          {
            (yylhs.value.val) = (yystack_[2].value.val) || 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 31: // bool: bool LAND UNDVAR
          {
            (yylhs.value.val) = (yystack_[2].value.val) && 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 32: // bool: UNDVAR LOR bool
          {
            (yylhs.value.val) = 0 || (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 33: // bool: UNDVAR LAND bool
          {
            (yylhs.value.val) = 0 && (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 34: // bool: sexp LT sexp
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) < 0 ? 1 : 0);
          } break;

          case 35: // bool: sexp GT sexp
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) > 0 ? 1 : 0);
          } break;

          case 36: // bool: sexp LE sexp
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) <= 0 ? 1 : 0);
          } break;

          case 37: // bool: sexp GE sexp
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) >= 0 ? 1 : 0);
          } break;

          case 38: // bool: sexp EQ sexp
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) == 0 ? 1 : 0);
          } break;

          case 39: // bool: sexp NE sexp
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) != 0 ? 1 : 0);
          } break;

          case 40: // bool: exp LT sexp
          {
            (yylhs.value.val) = false;
            yyerror(aprepro, "Comparison of arithmetic with string not defined");
            yyerrok;
          } break;

          case 41: // bool: exp GT sexp
          {
            (yylhs.value.val) = false;
            yyerror(aprepro, "Comparison of arithmetic with string not defined");
            yyerrok;
          } break;

          case 42: // bool: exp LE sexp
          {
            (yylhs.value.val) = false;
            yyerror(aprepro, "Comparison of arithmetic with string not defined");
            yyerrok;
          } break;

          case 43: // bool: exp GE sexp
          {
            (yylhs.value.val) = false;
            yyerror(aprepro, "Comparison of arithmetic with string not defined");
            yyerrok;
          } break;

          case 44: // bool: exp EQ sexp
          {
            (yylhs.value.val) = false;
          } break;

          case 45: // bool: exp NE sexp
          {
            (yylhs.value.val) = true;
          } break;

          case 46: // bool: sexp LT exp
          {
            (yylhs.value.val) = false;
            yyerror(aprepro, "Comparison of string with arithmetic not defined");
            yyerrok;
          } break;

          case 47: // bool: sexp GT exp
          {
            (yylhs.value.val) = false;
            yyerror(aprepro, "Comparison of string with arithmetic not defined");
            yyerrok;
          } break;

          case 48: // bool: sexp LE exp
          {
            (yylhs.value.val) = false;
            yyerror(aprepro, "Comparison of string with arithmetic not defined");
            yyerrok;
          } break;

          case 49: // bool: sexp GE exp
          {
            (yylhs.value.val) = false;
            yyerror(aprepro, "Comparison of string with arithmetic not defined");
            yyerrok;
          } break;

          case 50: // bool: sexp EQ exp
          {
            (yylhs.value.val) = false;
          } break;

          case 51: // bool: sexp NE exp
          {
            (yylhs.value.val) = true;
          } break;

          case 52: // bool: UNDVAR LT sexp
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) < 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 53: // bool: UNDVAR GT sexp
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) > 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 54: // bool: UNDVAR LE sexp
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) <= 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 55: // bool: UNDVAR GE sexp
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) >= 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 56: // bool: UNDVAR EQ sexp
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) == 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 57: // bool: UNDVAR NE sexp
          {
            (yylhs.value.val) = (strcmp("", (yystack_[0].value.string)) != 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 58: // bool: sexp LT UNDVAR
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") < 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 59: // bool: sexp GT UNDVAR
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") > 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 60: // bool: sexp LE UNDVAR
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") <= 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 61: // bool: sexp GE UNDVAR
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") >= 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 62: // bool: sexp EQ UNDVAR
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") == 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 63: // bool: sexp NE UNDVAR
          {
            (yylhs.value.val) = (strcmp((yystack_[2].value.string), "") != 0 ? 1 : 0);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 64: // bool: UNDVAR LT exp
          {
            (yylhs.value.val) = 0 < (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 65: // bool: UNDVAR GT exp
          {
            (yylhs.value.val) = 0 > (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 66: // bool: UNDVAR LE exp
          {
            (yylhs.value.val) = 0 <= (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 67: // bool: UNDVAR GE exp
          {
            (yylhs.value.val) = 0 >= (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 68: // bool: UNDVAR EQ exp
          {
            (yylhs.value.val) = 0 == (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 69: // bool: UNDVAR NE exp
          {
            (yylhs.value.val) = 0 != (yystack_[0].value.val);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 70: // bool: exp LT UNDVAR
          {
            (yylhs.value.val) = (yystack_[2].value.val) < 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 71: // bool: exp GT UNDVAR
          {
            (yylhs.value.val) = (yystack_[2].value.val) > 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 72: // bool: exp LE UNDVAR
          {
            (yylhs.value.val) = (yystack_[2].value.val) <= 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 73: // bool: exp GE UNDVAR
          {
            (yylhs.value.val) = (yystack_[2].value.val) >= 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 74: // bool: exp EQ UNDVAR
          {
            (yylhs.value.val) = (yystack_[2].value.val) == 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 75: // bool: exp NE UNDVAR
          {
            (yylhs.value.val) = (yystack_[2].value.val) != 0;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 76: // aexp: AVAR
          {
            (yylhs.value.arrval) = aprepro.make_array(*((yystack_[0].value.tptr)->value.avar));
          } break;

          case 77: // aexp: AFNCT LPAR sexp RPAR
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_c == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_c))((yystack_[1].value.string));
            else
              yyerrok;
          } break;

          case 78: // aexp: AFNCT LPAR sexp COMMA exp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_cd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              yyerrok;
          } break;

          case 79: // aexp: AFNCT LPAR sexp COMMA sexp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_cc == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              yyerrok;
          } break;

          case 80: // aexp: AFNCT LPAR exp COMMA exp COMMA exp RPAR
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.arrfnct_ddd == NULL))
              (yylhs.value.arrval) = (*((yystack_[7].value.tptr)->value.arrfnct_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          } break;

          case 81: // aexp: AFNCT LPAR exp COMMA exp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_dd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          } break;

          case 82: // aexp: AFNCT LPAR exp RPAR
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_d == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_d))((yystack_[1].value.val));
            else
              yyerrok;
          } break;

          case 83: // aexp: AFNCT LPAR aexp RPAR
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_a == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_a))((yystack_[1].value.arrval));
            else
              yyerrok;
          } break;

          case 84: // aexp: SVAR EQUAL aexp
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          } break;

          case 85: // aexp: VAR EQUAL aexp
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          } break;

          case 86: // aexp: AVAR EQUAL aexp
          {
            (yylhs.value.arrval) = (yystack_[0].value.arrval);
            aprepro.redefine_array((yystack_[2].value.tptr)->value.avar);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          } break;

          case 87: // aexp: UNDVAR EQUAL aexp
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          } break;

          case 88: // aexp: aexp PLU aexp
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
          } break;

          case 89: // aexp: SUB aexp
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), -1.0);
          } break;

          case 90: // aexp: aexp SUB aexp
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
          } break;

          case 91: // aexp: aexp TIM exp
          {
            (yylhs.value.arrval) = array_scale((yystack_[2].value.arrval), (yystack_[0].value.val));
          } break;

          case 92: // aexp: aexp DIV exp
          {
            (yylhs.value.arrval) =
                array_scale((yystack_[2].value.arrval), 1.0 / (yystack_[0].value.val));
          } break;

          case 93: // aexp: exp TIM aexp
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), (yystack_[2].value.val));
          } break;

          case 94: // aexp: aexp TIM aexp
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
          } break;

          case 95: // sexp: QSTRING
          {
            (yylhs.value.string) = (yystack_[0].value.string);
          } break;

          case 96: // sexp: SVAR
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          } break;

          case 97: // sexp: IMMSVAR
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          } break;

          case 98: // sexp: UNDVAR EQUAL sexp
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            set_type(aprepro, (yystack_[2].value.tptr), Parser::token::SVAR);
          } break;

          case 99: // sexp: SVAR EQUAL sexp
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          } break;

          case 100: // sexp: VAR EQUAL sexp
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);
          } break;

          case 101: // sexp: AVAR EQUAL sexp
          {
            (yylhs.value.string) = (yystack_[0].value.string);
            aprepro.redefine_array((yystack_[2].value.tptr)->value.avar);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);
          } break;

          case 102: // sexp: IMMSVAR EQUAL sexp
          {
            (yylhs.value.string) = (char *)(yystack_[2].value.tptr)->value.svar.c_str();
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          } break;

          case 103: // sexp: IMMVAR EQUAL sexp
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          } break;

          case 104: // sexp: SFNCT LPAR sexp RPAR
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_c == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_c))(
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 105: // sexp: SFNCT LPAR RPAR
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.strfnct == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[2].value.tptr)->value.strfnct))();
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 106: // sexp: SFNCT LPAR exp RPAR
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_d == NULL))
              (yylhs.value.string) =
                  (char *)(*((yystack_[3].value.tptr)->value.strfnct_d))((yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 107: // sexp: SFNCT LPAR aexp RPAR
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_a == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_a))(
                  (yystack_[1].value.arrval));
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 108: // sexp: sexp CONCAT sexp
          {
            concat_string((yystack_[2].value.string), (yystack_[0].value.string),
                          &(yylhs.value.string));
          } break;

          case 109: // sexp: sexp PLU sexp
          {
            concat_string((yystack_[2].value.string), (yystack_[0].value.string),
                          &(yylhs.value.string));
          } break;

          case 110: // sexp: SFNCT LPAR exp COMMA exp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_dd == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 111: // sexp: SFNCT LPAR exp COMMA sexp COMMA sexp COMMA sexp COMMA sexp RPAR
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.strfnct_dcccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[11].value.tptr)->value.strfnct_dcccc))(
                  (yystack_[9].value.val), (yystack_[7].value.string), (yystack_[5].value.string),
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 112: // sexp: SFNCT LPAR exp COMMA sexp COMMA sexp RPAR
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_dcc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_dcc))(
                  (yystack_[5].value.val), (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 113: // sexp: SFNCT LPAR exp COMMA sexp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_dc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_dc))(
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 114: // sexp: SFNCT LPAR sexp COMMA exp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_cd == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 115: // sexp: SFNCT LPAR sexp COMMA sexp COMMA sexp RPAR
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_ccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 116: // sexp: SFNCT LPAR sexp COMMA sexp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_cc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          } break;

          case 117: // sexp: bool QUEST sexp COLON sexp
          {
            (yylhs.value.string) = ((yystack_[4].value.val)) ? ((yystack_[2].value.string))
                                                             : ((yystack_[0].value.string));
          } break;

          case 118: // sexp: exp TIM sexp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Multiplying an arithmetic with a string is not defined");
            yyerrok;
          } break;

          case 119: // sexp: sexp TIM exp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Multiplying a string with an arithmetic is not defined");
            yyerrok;
          } break;

          case 120: // sexp: sexp TIM sexp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Multiplying a string with a string is not defined");
            yyerrok;
          } break;

          case 121: // sexp: sexp DIV exp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Dividing a string by an arithmetic is not defined");
            yyerrok;
          } break;

          case 122: // sexp: exp DIV sexp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Dividing an arithmetic by a string is not defined");
            yyerrok;
          } break;

          case 123: // sexp: sexp DIV sexp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Dividing a string by a string is not defined");
            yyerrok;
          } break;

          case 124: // sexp: exp PLU sexp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Adding an arithmetic and a string is not defined");
            yyerrok;
          } break;

          case 125: // sexp: sexp PLU exp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Adding a string and an arithmetic is not defined");
            yyerrok;
          } break;

          case 126: // sexp: exp SUB sexp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Subtracting an arithmetic and a string is not defined");
            yyerrok;
          } break;

          case 127: // sexp: sexp SUB exp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Subtracting a string and an arithmetic is not defined");
            yyerrok;
          } break;

          case 128: // sexp: sexp SUB sexp
          {
            (yylhs.value.string) = (char *)"";
            yyerror(aprepro, "Subtracting a string from a string is not defined");
            yyerrok;
          } break;

          case 129: // exp: NUM
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          } break;

          case 130: // exp: INC NUM
          {
            (yylhs.value.val) = (yystack_[0].value.val) + 1;
          } break;

          case 131: // exp: DEC NUM
          {
            (yylhs.value.val) = (yystack_[0].value.val) - 1;
          } break;

          case 132: // exp: VAR
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          } break;

          case 133: // exp: IMMVAR
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          } break;

          case 134: // exp: INC VAR
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
          } break;

          case 135: // exp: DEC VAR
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
          } break;

          case 136: // exp: VAR INC
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
          } break;

          case 137: // exp: VAR DEC
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
          } break;

          case 138: // exp: VAR EQUAL exp
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          } break;

          case 139: // exp: SVAR EQUAL exp
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          } break;

          case 140: // exp: AVAR EQUAL exp
          {
            (yylhs.value.val) = (yystack_[0].value.val);
            aprepro.redefine_array((yystack_[2].value.tptr)->value.avar);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          } break;

          case 141: // exp: VAR EQ_PLUS exp
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          } break;

          case 142: // exp: VAR EQ_MINUS exp
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          } break;

          case 143: // exp: VAR EQ_TIME exp
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          } break;

          case 144: // exp: VAR EQ_DIV exp
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          } break;

          case 145: // exp: VAR EQ_POW exp
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            SEAMS::math_error(aprepro, "Power");
          } break;

          case 146: // exp: INC IMMVAR
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          } break;

          case 147: // exp: DEC IMMVAR
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          } break;

          case 148: // exp: IMMVAR INC
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          } break;

          case 149: // exp: IMMVAR DEC
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          } break;

          case 150: // exp: IMMVAR EQUAL exp
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          } break;

          case 151: // exp: IMMSVAR EQUAL exp
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          } break;

          case 152: // exp: IMMVAR EQ_PLUS exp
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          } break;

          case 153: // exp: IMMVAR EQ_MINUS exp
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          } break;

          case 154: // exp: IMMVAR EQ_TIME exp
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          } break;

          case 155: // exp: IMMVAR EQ_DIV exp
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          } break;

          case 156: // exp: IMMVAR EQ_POW exp
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          } break;

          case 157: // exp: UNDVAR
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 158: // exp: INC UNDVAR
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 159: // exp: DEC UNDVAR
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          } break;

          case 160: // exp: UNDVAR INC
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          } break;

          case 161: // exp: UNDVAR DEC
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          } break;

          case 162: // exp: UNDVAR EQUAL exp
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          } break;

          case 163: // exp: UNDVAR EQ_PLUS exp
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 164: // exp: UNDVAR EQ_MINUS exp
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 165: // exp: UNDVAR EQ_TIME exp
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 166: // exp: UNDVAR EQ_DIV exp
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 167: // exp: UNDVAR EQ_POW exp
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            SEAMS::math_error(aprepro, "Power");
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          } break;

          case 168: // exp: FNCT LPAR RPAR
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.fnctptr == NULL))
              (yylhs.value.val) = (*((yystack_[2].value.tptr)->value.fnctptr))();
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 169: // exp: FNCT LPAR exp RPAR
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_d == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_d))((yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 170: // exp: FNCT LPAR sexp RPAR
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_c == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_c))((yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 171: // exp: FNCT LPAR aexp RPAR
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_a == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_a))((yystack_[1].value.arrval));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 172: // exp: FNCT LPAR sexp COMMA exp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 173: // exp: FNCT LPAR exp COMMA sexp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dc))(
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 174: // exp: FNCT LPAR sexp COMMA sexp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 175: // exp: FNCT LPAR sexp COMMA sexp COMMA sexp RPAR
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccc == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              yyerrok;
          } break;

          case 176: // exp: FNCT LPAR exp COMMA exp RPAR
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 177: // exp: FNCT LPAR exp COMMA exp COMMA exp RPAR
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ddd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 178: // exp: FNCT LPAR sexp COMMA sexp COMMA exp RPAR
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccd))(
                  (yystack_[5].value.string), (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 179: // exp: FNCT LPAR exp COMMA exp SEMI exp COMMA exp RPAR
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 180: // exp: FNCT LPAR exp COMMA exp COMMA exp COMMA exp RPAR
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 181: // exp: FNCT LPAR exp COMMA exp COMMA exp COMMA exp COMMA sexp RPAR
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.fnctptr_ddddc == NULL))
              (yylhs.value.val) = (*((yystack_[11].value.tptr)->value.fnctptr_ddddc))(
                  (yystack_[9].value.val), (yystack_[7].value.val), (yystack_[5].value.val),
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 182: // exp: FNCT LPAR exp COMMA exp COMMA exp COMMA exp COMMA exp COMMA exp RPAR
          {
            if (arg_check((yystack_[13].value.tptr),
                          (yystack_[13].value.tptr)->value.fnctptr_dddddd == NULL))
              (yylhs.value.val) = (*((yystack_[13].value.tptr)->value.fnctptr_dddddd))(
                  (yystack_[11].value.val), (yystack_[9].value.val), (yystack_[7].value.val),
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          } break;

          case 183: // exp: exp PLU exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) + (yystack_[0].value.val);
          } break;

          case 184: // exp: exp SUB exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) - (yystack_[0].value.val);
          } break;

          case 185: // exp: exp TIM exp
          {
            (yylhs.value.val) = (yystack_[2].value.val) * (yystack_[0].value.val);
          } break;

          case 186: // exp: exp DIV exp
          {
            if ((yystack_[0].value.val) == 0.) {
              (yylhs.value.val) = std::numeric_limits<double>::infinity();
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (yystack_[2].value.val) / (yystack_[0].value.val);
          } break;

          case 187: // exp: exp MOD exp
          {
            if ((int)(yystack_[0].value.val) == 0.) {
              (yylhs.value.val) = (int)(yystack_[2].value.val);
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (int)(yystack_[2].value.val) % (int)(yystack_[0].value.val);
          } break;

          case 188: // exp: SUB exp
          {
            (yylhs.value.val) = -(yystack_[0].value.val);
          } break;

          case 189: // exp: PLU exp
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          } break;

          case 190: // exp: exp POW exp
          {
            reset_error();
            (yylhs.value.val) = std::pow((yystack_[2].value.val), (yystack_[0].value.val));
            SEAMS::math_error(aprepro, "Power");
          } break;

          case 191: // exp: LPAR exp RPAR
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          } break;

          case 192: // exp: LBRACK exp RBRACK
          {
            reset_error();
            (yylhs.value.val) =
                (double)((yystack_[1].value.val) < 0 ? -floor(-((yystack_[1].value.val)))
                                                     : floor((yystack_[1].value.val)));
            SEAMS::math_error(aprepro, "floor (int)");
          } break;

          case 193: // exp: bool
          {
            (yylhs.value.val) = ((yystack_[0].value.val)) ? 1 : 0;
          } break;

          case 194: // exp: bool QUEST exp COLON exp
          {
            (yylhs.value.val) =
                ((yystack_[4].value.val)) ? ((yystack_[2].value.val)) : ((yystack_[0].value.val));
          } break;

          case 195: // exp: AVAR LBRACK exp RBRACK
          {
            (yylhs.value.val) =
                array_value((yystack_[3].value.tptr)->value.avar, (yystack_[1].value.val), 0);
          } break;

          case 196: // exp: AVAR LBRACK exp COMMA exp RBRACK
          {
            (yylhs.value.val) = array_value((yystack_[5].value.tptr)->value.avar,
                                            (yystack_[3].value.val), (yystack_[1].value.val));
          } break;

          case 197: // exp: AVAR LBRACK exp RBRACK EQUAL exp
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
          } break;

          case 198: // exp: AVAR LBRACK exp COMMA exp RBRACK EQUAL exp
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
          } break;

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

  const signed char Parser::yypact_ninf_ = -32;

  const signed char Parser::yytable_ninf_ = -1;

  const short Parser::yypact_[] = {
      -32,  22,   -32,  4,    353,  -32,  -32,  -32,  -32,  -32,  1894, -21,  20,   178,  28,
      140,  44,   51,   54,   121,  121,  -32,  121,  496,  121,  251,  334,  162,  293,  115,
      1930, 496,  121,  121,  121,  121,  121,  121,  121,  121,  121,  121,  121,  121,  121,
      -32,  -32,  496,  121,  121,  121,  121,  121,  -32,  -32,  496,  121,  121,  121,  121,
      121,  121,  -32,  -32,  121,  121,  496,  421,  481,  496,  1912, 325,  50,   158,  121,
      148,  1961, 1604, 1845, 39,   -32,  39,   39,   -32,  -32,  -32,  -32,  -32,  -32,  -32,
      -32,  121,  541,  556,  -32,  496,  496,  121,  496,  -32,  601,  616,  661,  676,  721,
      736,  121,  121,  121,  121,  121,  -32,  781,  796,  841,  856,  901,  916,  961,  976,
      121,  121,  121,  496,  121,  121,  188,  1961, 1980, 1996, 1996, 1996, 1996, 1996, 55,
      2027, -32,  2041, -31,  264,  -31,  264,  -31,  264,  -31,  264,  -31,  264,  -31,  264,
      188,  1961, 1980, 1996, 1996, 1996, 1996, 1996, 188,  1961, 1980, 1961, 1996, 1996, 1996,
      1996, 1996, 1996, 1961, 1996, 1339, 188,  1961, 1980, -32,  179,  175,  1369, -32,  245,
      995,  1399, 255,  1023, 1429, 121,  121,  121,  121,  39,   -32,  -32,  121,  -32,  345,
      1947, 1912, 55,   2027, 1912, -32,  2041, -29,  1980, -29,  2012, -32,  2012, 1912, -31,
      264,  1912, -31,  264,  1912, -31,  264,  1912, -31,  264,  1912, -31,  264,  1912, -31,
      264,  -27,  99,   -27,  99,   25,   39,   25,   39,   -32,  1996, 1912, 55,   2027, 1912,
      -32,  2041, 1912, -31,  264,  1912, -31,  264,  1912, -31,  264,  1912, -31,  264,  1912,
      -31,  264,  1912, -31,  264,  -27,  99,   -27,  99,   25,   39,   -32,  25,   39,   39,
      39,   121,  63,   -32,  121,  -32,  121,  -32,  -32,  121,  -32,  121,  -32,  -32,  121,
      -32,  121,  -32,  1996, 1996, 1996, 1996, 39,   121,  121,  1870, 121,  1051, 1631, 42,
      1310, 1079, 1658, 1107, 1685, 1178, 1712, 1459, 1961, 1996, 73,   1996, 121,  -32,  -32,
      -32,  121,  -32,  121,  121,  -32,  -32,  121,  -32,  -32,  -32,  -32,  121,  -32,  121,
      1204, 1739, 1489, 1548, 1230, 1135, 1766, 1996, -32,  -32,  121,  -32,  121,  -32,  121,
      -32,  -32,  1519, 1793, 1152, 121,  -32,  -32,  121,  1256, 1577, 1282, -32,  121,  -32,
      1820, -32};

  const unsigned char Parser::yydefact_[] = {
      2,   0,   1,   0,   0,   4,   3,   9,   129, 95,  157, 132, 96,  133, 97,  76,  0,   0,   0,
      0,   0,   8,   0,   0,   0,   0,   0,   193, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   160, 161, 0,   0,   0,   0,   0,   0,   136, 137, 0,   0,
      0,   0,   0,   0,   0,   148, 149, 0,   0,   0,   0,   0,   0,   157, 132, 96,  0,   0,   193,
      0,   0,   0,   189, 89,  188, 12,  130, 158, 134, 146, 131, 159, 135, 147, 0,   0,   0,   7,
      0,   0,   0,   0,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   5,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   87,  98,  162, 163, 164, 165, 166,
      167, 32,  26,  33,  27,  52,  64,  53,  65,  54,  66,  55,  67,  56,  68,  57,  69,  85,  100,
      138, 141, 142, 143, 144, 145, 84,  99,  139, 103, 150, 152, 153, 154, 155, 156, 102, 151, 0,
      86,  101, 140, 168, 0,   0,   0,   105, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   188,
      25,  191, 0,   192, 0,   0,   30,  19,  21,  31,  20,  22,  88,  0,   90,  92,  94,  91,  58,
      34,  46,  59,  35,  47,  60,  36,  48,  61,  37,  49,  62,  38,  50,  63,  39,  51,  109, 125,
      128, 127, 123, 121, 120, 119, 108, 0,   28,  23,  17,  29,  24,  18,  70,  40,  10,  71,  41,
      11,  72,  42,  13,  73,  43,  14,  74,  44,  15,  75,  45,  16,  124, 183, 126, 184, 122, 186,
      93,  118, 185, 187, 190, 0,   195, 171, 0,   170, 0,   169, 107, 0,   104, 0,   106, 83,  0,
      77,  0,   82,  162, 138, 139, 140, 185, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   117, 194, 196, 197, 0,   174, 172, 173, 0,   176, 0,   0,   116, 114, 0,
      113, 110, 79,  78,  0,   81,  0,   0,   0,   0,   0,   0,   0,   0,   198, 175, 178, 0,   177,
      0,   115, 0,   112, 80,  0,   0,   0,   0,   180, 179, 0,   0,   0,   0,   181, 0,   111, 0,
      182};

  const short Parser::yypgoto_[] = {-32, -32, -32, 30, 232, 130, -4};

  const unsigned char Parser::yydefgoto_[] = {0, 1, 6, 27, 28, 76, 235};

  const short Parser::yytable_[] = {
      30,  47,  48,  49,  50,  51,  52,  106, 107, 108, 109, 97,  98,  108, 109, 77,  78,  110, 79,
      81,  82,  110, 2,   3,   7,   53,  54,  128, 129, 130, 131, 132, 133, 135, 137, 139, 141, 143,
      145, 147, 149, 4,   55,  152, 153, 154, 155, 156, 157, 75,  64,  160, 162, 163, 164, 165, 166,
      167, 315, 67,  169, 170, 173, 177, 181, 184, 68,  134, 136, 69,  189, 5,   187, 110, 100, 101,
      102, 103, 104, 105, 106, 107, 108, 109, 125, 296, 93,  195, 198, 201, 110, 203, 203, 205, 207,
      329, 210, 213, 216, 219, 222, 225, 227, 229, 231, 233, 0,   0,   238, 241, 244, 247, 250, 253,
      256, 259, 261, 263, 265, 268, 269, 270, 197, 200, 8,   9,   70,  71,  72,  13,  14,  73,  16,
      17,  29,  99,  19,  0,   20,  122, 192, 124, 237, 240, 125, 0,   0,   100, 101, 102, 103, 104,
      105, 106, 107, 108, 109, 65,  0,   22,  74,  127, 66,  110, 190, 24,  0,   25,  26,  138, 140,
      142, 144, 146, 148, 65,  91,  151, 92,  93,  188, 288, 289, 290, 291, 159, 161, 0,   292, 274,
      91,  275, 92,  93,  168, 273, 172, 176, 180, 183, 56,  57,  58,  59,  60,  61,  0,   100, 101,
      102, 103, 104, 105, 106, 107, 108, 109, 95,  96,  97,  98,  194, 0,   110, 62,  63,  95,  96,
      97,  98,  209, 212, 215, 218, 221, 224, 226, 228, 230, 232, 234, 0,   0,   0,   243, 246, 249,
      252, 255, 258, 260, 262, 264, 267, 83,  80,  84,  85,  0,   86,  0,   278, 0,   126, 0,   0,
      0,   295, 0,   0,   298, 283, 300, 0,   0,   302, 0,   304, 0,   150, 306, 0,   307, 95,  96,
      97,  98,  158, 0,   0,   309, 0,   311, 95,  96,  97,  98,  0,   171, 175, 179, 182, 120, 121,
      122, 192, 124, 0,   331, 125, 0,   0,   332, 94,  333, 127, 151, 159, 172, 0,   0,   0,   267,
      336, 0,   337, 0,   202, 204, 0,   206, 95,  96,  97,  98,  0,   347, 87,  348, 88,  89,  0,
      90,  0,   0,   0,   355, 186, 48,  49,  50,  51,  52,  0,   360, 266, 8,   9,   10,  11,  12,
      13,  14,  15,  16,  17,  18,  0,   19,  0,   20,  53,  54,  21,  293, 0,   0,   100, 101, 102,
      103, 104, 105, 106, 107, 108, 109, 0,   0,   0,   0,   22,  23,  110, 0,   0,   0,   24,  0,
      25,  26,  0,   0,   0,   297, 0,   299, 0,   0,   301, 0,   303, 0,   0,   305, 0,   0,   0,
      0,   0,   0,   0,   0,   308, 8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,
      174, 20,  0,   0,   0,   330, 0,   0,   0,   0,   0,   0,   334, 0,   0,   335, 0,   0,   0,
      0,   0,   0,   22,  23,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,   0,   0,   349,
      0,   0,   0,   0,   0,   354, 0,   0,   356, 8,   9,   10,  11,  12,  13,  14,  15,  16,  17,
      18,  0,   19,  178, 20,  8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  0,
      20,  0,   0,   0,   0,   0,   22,  23,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,
      0,   0,   22,  23,  0,   0,   0,   0,   24,  0,   25,  26,  8,   9,   196, 71,  72,  13,  14,
      73,  16,  17,  0,   0,   19,  0,   20,  8,   9,   199, 71,  72,  13,  14,  73,  16,  17,  0,
      0,   19,  0,   20,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,
      0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  8,   9,   208, 71,
      72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  8,   9,   211, 71,  72,  13,  14,  73,
      16,  17,  0,   0,   19,  0,   20,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,
      0,   25,  26,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  8,
      9,   214, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  8,   9,   217, 71,  72,
      13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  0,   0,   0,   0,   0,   22,  74,  0,   0,
      0,   0,   24,  0,   25,  26,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,
      25,  26,  8,   9,   220, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  8,   9,
      223, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  0,   0,   0,   0,   0,   22,
      74,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,
      0,   24,  0,   25,  26,  8,   9,   236, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,
      20,  8,   9,   239, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  0,   0,   0,
      0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,   0,   0,   22,  74,
      0,   0,   0,   0,   24,  0,   25,  26,  8,   9,   242, 71,  72,  13,  14,  73,  16,  17,  0,
      0,   19,  0,   20,  8,   9,   245, 71,  72,  13,  14,  73,  16,  17,  0,   0,   19,  0,   20,
      0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,   0,
      0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  8,   9,   248, 71,  72,  13,  14,  73,
      16,  17,  0,   0,   19,  0,   20,  8,   9,   251, 71,  72,  13,  14,  73,  16,  17,  0,   0,
      19,  0,   20,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  0,
      0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  8,   9,   254, 71,  72,
      13,  14,  73,  16,  17,  0,   0,   19,  0,   20,  8,   9,   257, 71,  72,  13,  14,  73,  16,
      17,  0,   0,   19,  0,   20,  0,   0,   0,   0,   0,   22,  74,  0,   0,   0,   0,   24,  0,
      25,  26,  279, 0,   280, 0,   0,   22,  74,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,
      0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 284, 0,   285, 0,   0,   0,   110, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108,
      109, 312, 0,   313, 0,   0,   0,   110, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 319, 0,   320, 0,   0,   0,   110, 0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
      322, 0,   323, 0,   0,   0,   110, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100,
      101, 102, 103, 104, 105, 106, 107, 108, 109, 344, 0,   345, 0,   0,   0,   110, 0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   353, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 0,
      0,   0,   0,   0,   0,   110, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 325, 0,   0,
      0,   0,   0,   110, 0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105,
      106, 107, 108, 109, 338, 0,   0,   0,   0,   0,   110, 0,   0,   0,   0,   0,   0,   0,   0,
      0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 343, 0,   0,   0,   0,   0,   110, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 357,
      0,   0,   0,   0,   0,   110, 0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103,
      104, 105, 106, 107, 108, 109, 359, 0,   0,   0,   0,   0,   110, 0,   0,   0,   0,   0,   0,
      0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 316, 0,   317, 0,   0,   0,
      110, 318, 0,   0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117, 118, 119, 120,
      121, 122, 192, 124, 271, 0,   125, 0,   272, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 0,   276, 125, 277, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117, 118,
      119, 120, 121, 122, 123, 124, 0,   281, 125, 282, 0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 0,   286,
      125, 287, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115,
      116, 117, 118, 119, 120, 121, 122, 123, 124, 0,   327, 125, 328, 0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192,
      124, 0,   340, 125, 341, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   112,
      113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 0,   350, 125, 351, 0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117, 118, 119, 120,
      121, 122, 192, 124, 342, 0,   125, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 358, 0,   125, 0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117, 118, 119,
      120, 121, 122, 192, 124, 191, 0,   125, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 314, 0,   125, 0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
      192, 124, 321, 0,   125, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   112, 113, 114,
      115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 324, 0,   125, 0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 326,
      0,   125, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117,
      118, 119, 120, 121, 122, 192, 124, 339, 0,   125, 0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 346, 0,   125, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117, 118, 119, 120,
      121, 122, 192, 124, 352, 0,   125, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   112,
      113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 361, 0,   125, 0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192,
      124, 193, 0,   125, 0,   0,   0,   0,   0,   0,   0,   0,   0,   112, 113, 114, 115, 116, 117,
      118, 119, 120, 121, 122, 192, 124, 310, 0,   125, 0,   0,   0,   0,   0,   0,   0,   0,   0,
      112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 0,   0,   125, 31,  32,  33,
      34,  35,  36,  0,   0,   37,  38,  39,  40,  41,  42,  43,  44,  0,   0,   185, 32,  33,  34,
      35,  36,  45,  46,  37,  38,  39,  40,  41,  42,  43,  44,  111, 0,   0,   0,   0,   0,   0,
      0,   45,  46,  112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 0,   0,   125,
      294, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 0,   0,   125, 100, 101,
      102, 103, 104, 105, 106, 107, 108, 109, 0,   0,   0,   0,   0,   0,   110, 112, 113, 114, 115,
      116, 117, 118, 119, 120, 121, 122, 123, 124, 0,   0,   125, 112, 113, 114, 115, 116, 117, 118,
      119, 120, 121, 122, 192, 124, 0,   0,   125, 112, 113, 114, 115, 116, 117, 118, 119, 0,   0,
      0,   0,   124, 0,   0,   125, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 0,
      0,   125, 114, 115, 116, 117, 118, 119, 120, 121, 122, 192, 124, 0,   0,   125};

  const short Parser::yycheck_[] = {
      4,   22,  23,  24,  25,  26,  27,  38,  39,  40,  41,  40,  41,  40,  41,  19,  20,  48,  22,
      23,  24,  48,  0,   1,   20,  46,  47,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
      42,  43,  44,  19,  22,  47,  48,  49,  50,  51,  52,  19,  22,  55,  56,  57,  58,  59,  60,
      61,  16,  15,  64,  65,  66,  67,  68,  69,  15,  37,  38,  15,  74,  49,  22,  48,  32,  33,
      34,  35,  36,  37,  38,  39,  40,  41,  45,  22,  31,  91,  92,  93,  48,  95,  96,  97,  98,
      22,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, -1,  -1,  112, 113, 114, 115, 116, 117,
      118, 119, 120, 121, 122, 123, 124, 125, 92,  93,  3,   4,   5,   6,   7,   8,   9,   10,  11,
      12,  4,   20,  15,  -1,  17,  40,  41,  42,  112, 113, 45,  -1,  -1,  32,  33,  34,  35,  36,
      37,  38,  39,  40,  41,  17,  -1,  38,  39,  31,  22,  48,  16,  44,  -1,  46,  47,  39,  40,
      41,  42,  43,  44,  17,  28,  47,  30,  31,  22,  185, 186, 187, 188, 55,  56,  -1,  192, 14,
      28,  16,  30,  31,  64,  16,  66,  67,  68,  69,  22,  23,  24,  25,  26,  27,  -1,  32,  33,
      34,  35,  36,  37,  38,  39,  40,  41,  38,  39,  40,  41,  91,  -1,  48,  46,  47,  38,  39,
      40,  41,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, -1,  -1,  -1,  114, 115, 116,
      117, 118, 119, 120, 121, 122, 123, 3,   23,  5,   6,   -1,  8,   -1,  16,  -1,  31,  -1,  -1,
      -1,  271, -1,  -1,  274, 16,  276, -1,  -1,  279, -1,  281, -1,  47,  284, -1,  286, 38,  39,
      40,  41,  55,  -1,  -1,  294, -1,  296, 38,  39,  40,  41,  -1,  66,  67,  68,  69,  38,  39,
      40,  41,  42,  -1,  312, 45,  -1,  -1,  316, 20,  318, 185, 186, 187, 188, -1,  -1,  -1,  192,
      327, -1,  329, -1,  95,  96,  -1,  98,  38,  39,  40,  41,  -1,  340, 3,   342, 5,   6,   -1,
      8,   -1,  -1,  -1,  350, 22,  23,  24,  25,  26,  27,  -1,  358, 123, 3,   4,   5,   6,   7,
      8,   9,   10,  11,  12,  13,  -1,  15,  -1,  17,  46,  47,  20,  29,  -1,  -1,  32,  33,  34,
      35,  36,  37,  38,  39,  40,  41,  -1,  -1,  -1,  -1,  38,  39,  48,  -1,  -1,  -1,  44,  -1,
      46,  47,  -1,  -1,  -1,  274, -1,  276, -1,  -1,  279, -1,  281, -1,  -1,  284, -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  293, 3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,
      16,  17,  -1,  -1,  -1,  312, -1,  -1,  -1,  -1,  -1,  -1,  319, -1,  -1,  322, -1,  -1,  -1,
      -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,  -1,  -1,  344,
      -1,  -1,  -1,  -1,  -1,  350, -1,  -1,  353, 3,   4,   5,   6,   7,   8,   9,   10,  11,  12,
      13,  -1,  15,  16,  17,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,  -1,
      17,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,
      -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  3,   4,   5,   6,   7,   8,   9,
      10,  11,  12,  -1,  -1,  15,  -1,  17,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,
      -1,  15,  -1,  17,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,
      -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  3,   4,   5,   6,
      7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,  17,  3,   4,   5,   6,   7,   8,   9,   10,
      11,  12,  -1,  -1,  15,  -1,  17,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,
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
      46,  47,  14,  -1,  16,  -1,  -1,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  -1,
      -1,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  14,  -1,  16,  -1,  -1,  -1,  48,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  38,  39,  40,
      41,  14,  -1,  16,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  14,  -1,  16,  -1,  -1,  -1,  48,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
      14,  -1,  16,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  32,
      33,  34,  35,  36,  37,  38,  39,  40,  41,  14,  -1,  16,  -1,  -1,  -1,  48,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  14,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  -1,
      -1,  -1,  -1,  -1,  -1,  48,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  16,  -1,  -1,
      -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,
      38,  39,  40,  41,  16,  -1,  -1,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  16,  -1,  -1,  -1,  -1,  -1,  48,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  16,
      -1,  -1,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,
      36,  37,  38,  39,  40,  41,  16,  -1,  -1,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  14,  -1,  16,  -1,  -1,  -1,
      48,  21,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,
      39,  40,  41,  42,  14,  -1,  45,  -1,  18,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,
      37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,
      45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,
      34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
      42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,
      39,  40,  41,  42,  14,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  14,  -1,  45,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,
      38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
      41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,
      33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,
      -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,
      36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,
      39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
      42,  18,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,
      36,  37,  38,  39,  40,  41,  42,  18,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  22,  23,  24,
      25,  26,  27,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  -1,  -1,  22,  23,  24,  25,
      26,  27,  46,  47,  30,  31,  32,  33,  34,  35,  36,  37,  20,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  46,  47,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,
      29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  32,  33,
      34,  35,  36,  37,  38,  39,  40,  41,  -1,  -1,  -1,  -1,  -1,  -1,  48,  30,  31,  32,  33,
      34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  30,  31,  32,  33,  34,  35,  36,
      37,  38,  39,  40,  41,  42,  -1,  -1,  45,  30,  31,  32,  33,  34,  35,  36,  37,  -1,  -1,
      -1,  -1,  42,  -1,  -1,  45,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,
      -1,  45,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45};

  const signed char Parser::yystos_[] = {
      0,  51, 0,  1,  19, 49, 52, 20, 3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 15, 17, 20, 38,
      39, 44, 46, 47, 53, 54, 55, 56, 22, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 35, 36, 37, 46,
      47, 22, 23, 24, 25, 26, 27, 46, 47, 22, 22, 23, 24, 25, 26, 27, 46, 47, 22, 17, 22, 15, 15,
      15, 5,  6,  7,  10, 39, 53, 55, 56, 56, 56, 54, 56, 56, 3,  5,  6,  8,  3,  5,  6,  8,  28,
      30, 31, 20, 38, 39, 40, 41, 20, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 48, 20, 30, 31, 32,
      33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 45, 54, 55, 56, 56, 56, 56, 56, 56, 53, 56, 53, 56,
      55, 56, 55, 56, 55, 56, 55, 56, 55, 56, 55, 56, 54, 55, 56, 56, 56, 56, 56, 56, 54, 55, 56,
      55, 56, 56, 56, 56, 56, 56, 55, 56, 56, 54, 55, 56, 16, 54, 55, 56, 16, 54, 55, 56, 54, 55,
      56, 22, 22, 22, 22, 56, 16, 16, 41, 18, 55, 56, 5,  53, 56, 5,  53, 56, 54, 56, 54, 56, 54,
      56, 5,  55, 56, 5,  55, 56, 5,  55, 56, 5,  55, 56, 5,  55, 56, 5,  55, 56, 55, 56, 55, 56,
      55, 56, 55, 56, 55, 56, 5,  53, 56, 5,  53, 56, 5,  55, 56, 5,  55, 56, 5,  55, 56, 5,  55,
      56, 5,  55, 56, 5,  55, 56, 55, 56, 55, 56, 55, 56, 54, 55, 56, 56, 56, 14, 18, 16, 14, 16,
      14, 16, 16, 14, 16, 14, 16, 16, 14, 16, 14, 16, 56, 56, 56, 56, 56, 29, 29, 56, 22, 55, 56,
      55, 56, 55, 56, 55, 56, 55, 56, 56, 55, 56, 18, 56, 14, 16, 16, 16, 14, 16, 21, 14, 16, 16,
      14, 16, 16, 16, 16, 14, 16, 22, 55, 56, 56, 56, 55, 55, 56, 56, 16, 16, 14, 16, 14, 16, 14,
      16, 16, 56, 56, 55, 14, 16, 16, 14, 55, 56, 55, 16, 14, 16, 56, 16};

  const signed char Parser::yyr1_[] = {
      0,  50, 51, 51, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
      53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
      54, 54, 54, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55,
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56};

  const signed char Parser::yyr2_[] = {
      0, 2, 0, 2, 1, 3,  3,  3,  2,  2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  3, 3, 3, 3,
      3, 3, 3, 3, 3, 3,  3,  3,  3,  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  3, 3, 3, 3,
      3, 3, 3, 3, 3, 3,  3,  3,  3,  3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 4, 6, 6, 8, 6, 4,  4, 3, 3, 3,
      3, 3, 2, 3, 3, 3,  3,  3,  1,  1, 1, 3, 3, 3, 3, 3, 3, 4, 3, 4, 4, 3, 3, 6, 12, 8, 6, 6, 8,
      6, 5, 3, 3, 3, 3,  3,  3,  3,  3, 3, 3, 3, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 3,  3, 3, 3, 3,
      3, 2, 2, 2, 2, 3,  3,  3,  3,  3, 3, 3, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4,  4, 4, 6, 6,
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
      0,   127, 127, 128, 131, 132, 145, 149, 150, 151, 154, 155, 156, 157, 158, 159, 160, 161, 162,
      163, 164, 165, 166, 167, 168, 169, 171, 172, 173, 174, 175, 176, 177, 178, 181, 182, 183, 184,
      185, 186, 188, 189, 190, 191, 192, 193, 195, 196, 197, 198, 199, 200, 202, 203, 204, 205, 206,
      207, 209, 210, 211, 212, 213, 214, 216, 217, 218, 219, 220, 221, 223, 224, 225, 226, 227, 228,
      230, 231, 237, 243, 249, 255, 261, 267, 273, 277, 281, 284, 286, 294, 296, 304, 305, 306, 307,
      316, 317, 318, 319, 321, 324, 328, 333, 334, 335, 341, 347, 353, 359, 360, 361, 367, 373, 379,
      385, 391, 397, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 416, 417, 418, 419,
      420, 421, 422, 423, 424, 425, 427, 430, 435, 436, 437, 438, 439, 444, 445, 446, 447, 448, 449,
      450, 451, 452, 453, 454, 456, 458, 461, 464, 467, 470, 472, 475, 478, 481, 484, 491, 498, 505,
      512, 519, 526, 533, 540, 547, 553, 559, 565, 571, 577, 583, 589, 590, 591, 592, 600, 608, 609,
      610, 613, 614, 617, 618, 619, 620, 621, 643};

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

void SEAMS::Parser::error(const std::string &m) { aprepro.error(m); }
