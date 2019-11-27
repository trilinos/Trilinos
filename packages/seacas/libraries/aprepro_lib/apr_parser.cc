// A Bison parser, made by GNU Bison 3.4.

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
#line 33 "aprepro.yy"

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

#line 74 "apr_parser.cc"

#include "aprepro_parser.h"

// Second part of user prologue.
#line 130 "aprepro.yy"

#include "apr_scanner.h"
#include "aprepro.h"

/* this "connects" the bison parser in aprepro to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the aprepro context. */
#undef yylex
#define yylex aprepro.lexer->lex

#line 93 "apr_parser.cc"

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
#line 169 "apr_parser.cc"

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
#line 149 "aprepro.yy"
          {
            if (echo)
              aprepro.lexer->LexerOutput("\n", 1);
          }
#line 638 "apr_parser.cc"
          break;

          case 5:
#line 150 "aprepro.yy"
          {
            if (echo) {
              static char    tmpstr[512];
              SEAMS::symrec *format = aprepro.getsym("_FORMAT");
              int len = sprintf(tmpstr, format->value.svar.c_str(), (yystack_[1].value.val));
              aprepro.lexer->LexerOutput(tmpstr, len);
            }
          }
#line 650 "apr_parser.cc"
          break;

          case 6:
#line 157 "aprepro.yy"
          {
            if (echo && (yystack_[1].value.string) != NULL) {
              aprepro.lexer->LexerOutput((yystack_[1].value.string),
                                         strlen((yystack_[1].value.string)));
            }
          }
#line 659 "apr_parser.cc"
          break;

          case 7:
#line 161 "aprepro.yy"
          {
          }
#line 665 "apr_parser.cc"
          break;

          case 8:
#line 162 "aprepro.yy"
          {
          }
#line 671 "apr_parser.cc"
          break;

          case 9:
#line 163 "aprepro.yy"
          {
            yyerrok;
          }
#line 677 "apr_parser.cc"
          break;

          case 10:
#line 166 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) < (yystack_[0].value.val);
          }
#line 683 "apr_parser.cc"
          break;

          case 11:
#line 167 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) > (yystack_[0].value.val);
          }
#line 689 "apr_parser.cc"
          break;

          case 12:
#line 168 "aprepro.yy"
          {
            (yylhs.value.val) = !((yystack_[0].value.val));
          }
#line 695 "apr_parser.cc"
          break;

          case 13:
#line 169 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) <= (yystack_[0].value.val);
          }
#line 701 "apr_parser.cc"
          break;

          case 14:
#line 170 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) >= (yystack_[0].value.val);
          }
#line 707 "apr_parser.cc"
          break;

          case 15:
#line 171 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) == (yystack_[0].value.val);
          }
#line 713 "apr_parser.cc"
          break;

          case 16:
#line 172 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) != (yystack_[0].value.val);
          }
#line 719 "apr_parser.cc"
          break;

          case 17:
#line 173 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 725 "apr_parser.cc"
          break;

          case 18:
#line 174 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 731 "apr_parser.cc"
          break;

          case 19:
#line 175 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 737 "apr_parser.cc"
          break;

          case 20:
#line 176 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 743 "apr_parser.cc"
          break;

          case 21:
#line 177 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          }
#line 749 "apr_parser.cc"
          break;

          case 22:
#line 180 "aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) < 0 ? 1 : 0);
          }
#line 755 "apr_parser.cc"
          break;

          case 23:
#line 181 "aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) > 0 ? 1 : 0);
          }
#line 761 "apr_parser.cc"
          break;

          case 24:
#line 182 "aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) <= 0 ? 1 : 0);
          }
#line 767 "apr_parser.cc"
          break;

          case 25:
#line 183 "aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) >= 0 ? 1 : 0);
          }
#line 773 "apr_parser.cc"
          break;

          case 26:
#line 184 "aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) == 0 ? 1 : 0);
          }
#line 779 "apr_parser.cc"
          break;

          case 27:
#line 185 "aprepro.yy"
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) != 0 ? 1 : 0);
          }
#line 785 "apr_parser.cc"
          break;

          case 28:
#line 187 "aprepro.yy"
          {
            (yylhs.value.arrval) = (yystack_[0].value.tptr)->value.avar;
          }
#line 791 "apr_parser.cc"
          break;

          case 29:
#line 188 "aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_c == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_c))((yystack_[1].value.string));
            else
              yyerrok;
          }
#line 802 "apr_parser.cc"
          break;

          case 30:
#line 194 "aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_cd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 813 "apr_parser.cc"
          break;

          case 31:
#line 200 "aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_cc == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              yyerrok;
          }
#line 824 "apr_parser.cc"
          break;

          case 32:
#line 206 "aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.arrfnct_ddd == NULL))
              (yylhs.value.arrval) = (*((yystack_[7].value.tptr)->value.arrfnct_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 835 "apr_parser.cc"
          break;

          case 33:
#line 212 "aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_dd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 846 "apr_parser.cc"
          break;

          case 34:
#line 218 "aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_d == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_d))((yystack_[1].value.val));
            else
              yyerrok;
          }
#line 857 "apr_parser.cc"
          break;

          case 35:
#line 224 "aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_a == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_a))((yystack_[1].value.arrval));
            else
              yyerrok;
          }
#line 868 "apr_parser.cc"
          break;

          case 36:
#line 230 "aprepro.yy"
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 877 "apr_parser.cc"
          break;

          case 37:
#line 234 "aprepro.yy"
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 886 "apr_parser.cc"
          break;

          case 38:
#line 238 "aprepro.yy"
          {
            (yylhs.value.arrval) = (yystack_[0].value.arrval);
            delete (yystack_[2].value.tptr)->value.avar;
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 894 "apr_parser.cc"
          break;

          case 39:
#line 241 "aprepro.yy"
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 901 "apr_parser.cc"
          break;

          case 40:
#line 243 "aprepro.yy"
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
#line 914 "apr_parser.cc"
          break;

          case 41:
#line 251 "aprepro.yy"
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), -1.0);
          }
#line 920 "apr_parser.cc"
          break;

          case 42:
#line 253 "aprepro.yy"
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
#line 933 "apr_parser.cc"
          break;

          case 43:
#line 261 "aprepro.yy"
          {
            (yylhs.value.arrval) = array_scale((yystack_[2].value.arrval), (yystack_[0].value.val));
          }
#line 939 "apr_parser.cc"
          break;

          case 44:
#line 262 "aprepro.yy"
          {
            (yylhs.value.arrval) =
                array_scale((yystack_[2].value.arrval), 1.0 / (yystack_[0].value.val));
          }
#line 945 "apr_parser.cc"
          break;

          case 45:
#line 263 "aprepro.yy"
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), (yystack_[2].value.val));
          }
#line 951 "apr_parser.cc"
          break;

          case 46:
#line 264 "aprepro.yy"
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
#line 964 "apr_parser.cc"
          break;

          case 47:
#line 273 "aprepro.yy"
          {
            (yylhs.value.string) = (yystack_[0].value.string);
          }
#line 970 "apr_parser.cc"
          break;

          case 48:
#line 274 "aprepro.yy"
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          }
#line 976 "apr_parser.cc"
          break;

          case 49:
#line 275 "aprepro.yy"
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          }
#line 982 "apr_parser.cc"
          break;

          case 50:
#line 276 "aprepro.yy"
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            set_type(aprepro, (yystack_[2].value.tptr), Parser::token::SVAR);
          }
#line 989 "apr_parser.cc"
          break;

          case 51:
#line 278 "aprepro.yy"
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          }
#line 997 "apr_parser.cc"
          break;

          case 52:
#line 281 "aprepro.yy"
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);
          }
#line 1006 "apr_parser.cc"
          break;

          case 53:
#line 285 "aprepro.yy"
          {
            (yylhs.value.string) = (yystack_[0].value.string);
            delete (yystack_[2].value.tptr)->value.avar;
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);
          }
#line 1016 "apr_parser.cc"
          break;

          case 54:
#line 290 "aprepro.yy"
          {
            (yylhs.value.string) = (char *)(yystack_[2].value.tptr)->value.svar.c_str();
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1022 "apr_parser.cc"
          break;

          case 55:
#line 291 "aprepro.yy"
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          }
#line 1028 "apr_parser.cc"
          break;

          case 56:
#line 292 "aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_c == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_c))(
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1039 "apr_parser.cc"
          break;

          case 57:
#line 298 "aprepro.yy"
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.strfnct == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[2].value.tptr)->value.strfnct))();
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1050 "apr_parser.cc"
          break;

          case 58:
#line 304 "aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_d == NULL))
              (yylhs.value.string) =
                  (char *)(*((yystack_[3].value.tptr)->value.strfnct_d))((yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1061 "apr_parser.cc"
          break;

          case 59:
#line 310 "aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_a == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_a))(
                  (yystack_[1].value.arrval));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1072 "apr_parser.cc"
          break;

          case 60:
#line 316 "aprepro.yy"
          {
            concat_string((yystack_[2].value.string), (yystack_[0].value.string),
                          &(yylhs.value.string));
          }
#line 1078 "apr_parser.cc"
          break;

          case 61:
#line 317 "aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_dd == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1089 "apr_parser.cc"
          break;

          case 62:
#line 323 "aprepro.yy"
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.strfnct_dcccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[11].value.tptr)->value.strfnct_dcccc))(
                  (yystack_[9].value.val), (yystack_[7].value.string), (yystack_[5].value.string),
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1100 "apr_parser.cc"
          break;

          case 63:
#line 329 "aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_dcc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_dcc))(
                  (yystack_[5].value.val), (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1111 "apr_parser.cc"
          break;

          case 64:
#line 335 "aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_ccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1122 "apr_parser.cc"
          break;

          case 65:
#line 341 "aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_cc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1133 "apr_parser.cc"
          break;

          case 66:
#line 347 "aprepro.yy"
          {
            (yylhs.value.string) = ((yystack_[4].value.val)) ? ((yystack_[2].value.string))
                                                             : ((yystack_[0].value.string));
          }
#line 1139 "apr_parser.cc"
          break;

          case 67:
#line 349 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          }
#line 1145 "apr_parser.cc"
          break;

          case 68:
#line 350 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val) + 1;
          }
#line 1151 "apr_parser.cc"
          break;

          case 69:
#line 351 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val) - 1;
          }
#line 1157 "apr_parser.cc"
          break;

          case 70:
#line 352 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          }
#line 1163 "apr_parser.cc"
          break;

          case 71:
#line 353 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          }
#line 1169 "apr_parser.cc"
          break;

          case 72:
#line 354 "aprepro.yy"
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
          }
#line 1175 "apr_parser.cc"
          break;

          case 73:
#line 355 "aprepro.yy"
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
          }
#line 1181 "apr_parser.cc"
          break;

          case 74:
#line 356 "aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
          }
#line 1187 "apr_parser.cc"
          break;

          case 75:
#line 357 "aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
          }
#line 1193 "apr_parser.cc"
          break;

          case 76:
#line 358 "aprepro.yy"
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          }
#line 1200 "apr_parser.cc"
          break;

          case 77:
#line 360 "aprepro.yy"
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1208 "apr_parser.cc"
          break;

          case 78:
#line 363 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val);
            delete (yystack_[2].value.tptr)->value.avar;
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1218 "apr_parser.cc"
          break;

          case 79:
#line 368 "aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1224 "apr_parser.cc"
          break;

          case 80:
#line 369 "aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1230 "apr_parser.cc"
          break;

          case 81:
#line 370 "aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1236 "apr_parser.cc"
          break;

          case 82:
#line 371 "aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1242 "apr_parser.cc"
          break;

          case 83:
#line 372 "aprepro.yy"
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            SEAMS::math_error(aprepro, "Power");
          }
#line 1252 "apr_parser.cc"
          break;

          case 84:
#line 377 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          }
#line 1258 "apr_parser.cc"
          break;

          case 85:
#line 378 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          }
#line 1264 "apr_parser.cc"
          break;

          case 86:
#line 379 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          }
#line 1270 "apr_parser.cc"
          break;

          case 87:
#line 380 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          }
#line 1276 "apr_parser.cc"
          break;

          case 88:
#line 381 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1282 "apr_parser.cc"
          break;

          case 89:
#line 382 "aprepro.yy"
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          }
#line 1288 "apr_parser.cc"
          break;

          case 90:
#line 383 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1294 "apr_parser.cc"
          break;

          case 91:
#line 384 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1300 "apr_parser.cc"
          break;

          case 92:
#line 385 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1306 "apr_parser.cc"
          break;

          case 93:
#line 386 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1312 "apr_parser.cc"
          break;

          case 94:
#line 387 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1318 "apr_parser.cc"
          break;

          case 95:
#line 389 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1325 "apr_parser.cc"
          break;

          case 96:
#line 391 "aprepro.yy"
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1333 "apr_parser.cc"
          break;

          case 97:
#line 394 "aprepro.yy"
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1341 "apr_parser.cc"
          break;

          case 98:
#line 397 "aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          }
#line 1349 "apr_parser.cc"
          break;

          case 99:
#line 400 "aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          }
#line 1357 "apr_parser.cc"
          break;

          case 100:
#line 403 "aprepro.yy"
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1364 "apr_parser.cc"
          break;

          case 101:
#line 405 "aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1372 "apr_parser.cc"
          break;

          case 102:
#line 408 "aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1380 "apr_parser.cc"
          break;

          case 103:
#line 411 "aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1388 "apr_parser.cc"
          break;

          case 104:
#line 414 "aprepro.yy"
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1396 "apr_parser.cc"
          break;

          case 105:
#line 417 "aprepro.yy"
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            SEAMS::math_error(aprepro, "Power");
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1407 "apr_parser.cc"
          break;

          case 106:
#line 424 "aprepro.yy"
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.fnctptr == NULL))
              (yylhs.value.val) = (*((yystack_[2].value.tptr)->value.fnctptr))();
            else
              (yylhs.value.val) = 0.0;
          }
#line 1418 "apr_parser.cc"
          break;

          case 107:
#line 431 "aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_d == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_d))((yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1429 "apr_parser.cc"
          break;

          case 108:
#line 438 "aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_c == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_c))((yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1440 "apr_parser.cc"
          break;

          case 109:
#line 445 "aprepro.yy"
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_a == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_a))((yystack_[1].value.arrval));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1451 "apr_parser.cc"
          break;

          case 110:
#line 452 "aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1462 "apr_parser.cc"
          break;

          case 111:
#line 459 "aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dc))(
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1473 "apr_parser.cc"
          break;

          case 112:
#line 466 "aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1484 "apr_parser.cc"
          break;

          case 113:
#line 473 "aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccc == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              yyerrok;
          }
#line 1495 "apr_parser.cc"
          break;

          case 114:
#line 480 "aprepro.yy"
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1506 "apr_parser.cc"
          break;

          case 115:
#line 486 "aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ddd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1517 "apr_parser.cc"
          break;

          case 116:
#line 492 "aprepro.yy"
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccd))(
                  (yystack_[5].value.string), (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1528 "apr_parser.cc"
          break;

          case 117:
#line 498 "aprepro.yy"
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1539 "apr_parser.cc"
          break;

          case 118:
#line 504 "aprepro.yy"
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1550 "apr_parser.cc"
          break;

          case 119:
#line 510 "aprepro.yy"
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.fnctptr_ddddc == NULL))
              (yylhs.value.val) = (*((yystack_[11].value.tptr)->value.fnctptr_ddddc))(
                  (yystack_[9].value.val), (yystack_[7].value.val), (yystack_[5].value.val),
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1561 "apr_parser.cc"
          break;

          case 120:
#line 516 "aprepro.yy"
          {
            if (arg_check((yystack_[13].value.tptr),
                          (yystack_[13].value.tptr)->value.fnctptr_dddddd == NULL))
              (yylhs.value.val) = (*((yystack_[13].value.tptr)->value.fnctptr_dddddd))(
                  (yystack_[11].value.val), (yystack_[9].value.val), (yystack_[7].value.val),
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1572 "apr_parser.cc"
          break;

          case 121:
#line 522 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) + (yystack_[0].value.val);
          }
#line 1578 "apr_parser.cc"
          break;

          case 122:
#line 523 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) - (yystack_[0].value.val);
          }
#line 1584 "apr_parser.cc"
          break;

          case 123:
#line 524 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[2].value.val) * (yystack_[0].value.val);
          }
#line 1590 "apr_parser.cc"
          break;

          case 124:
#line 525 "aprepro.yy"
          {
            if ((yystack_[0].value.val) == 0.) {
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (yystack_[2].value.val) / (yystack_[0].value.val);
          }
#line 1602 "apr_parser.cc"
          break;

          case 125:
#line 532 "aprepro.yy"
          {
            if ((yystack_[0].value.val) == 0.) {
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (int)(yystack_[2].value.val) % (int)(yystack_[0].value.val);
          }
#line 1614 "apr_parser.cc"
          break;

          case 126:
#line 539 "aprepro.yy"
          {
            (yylhs.value.val) = -(yystack_[0].value.val);
          }
#line 1620 "apr_parser.cc"
          break;

          case 127:
#line 540 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          }
#line 1626 "apr_parser.cc"
          break;

          case 128:
#line 541 "aprepro.yy"
          {
            reset_error();
            (yylhs.value.val) = std::pow((yystack_[2].value.val), (yystack_[0].value.val));
            SEAMS::math_error(aprepro, "Power");
          }
#line 1634 "apr_parser.cc"
          break;

          case 129:
#line 544 "aprepro.yy"
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          }
#line 1640 "apr_parser.cc"
          break;

          case 130:
#line 545 "aprepro.yy"
          {
            reset_error();
            (yylhs.value.val) =
                (double)((yystack_[1].value.val) < 0 ? -floor(-((yystack_[1].value.val)))
                                                     : floor((yystack_[1].value.val)));
            SEAMS::math_error(aprepro, "floor (int)");
          }
#line 1648 "apr_parser.cc"
          break;

          case 131:
#line 548 "aprepro.yy"
          {
            (yylhs.value.val) = ((yystack_[0].value.val)) ? 1 : 0;
          }
#line 1654 "apr_parser.cc"
          break;

          case 132:
#line 549 "aprepro.yy"
          {
            (yylhs.value.val) =
                ((yystack_[4].value.val)) ? ((yystack_[2].value.val)) : ((yystack_[0].value.val));
          }
#line 1660 "apr_parser.cc"
          break;

          case 133:
#line 550 "aprepro.yy"
          {
            (yylhs.value.val) =
                array_value((yystack_[3].value.tptr)->value.avar, (yystack_[1].value.val), 0);
          }
#line 1666 "apr_parser.cc"
          break;

          case 134:
#line 551 "aprepro.yy"
          {
            (yylhs.value.val) = array_value((yystack_[5].value.tptr)->value.avar,
                                            (yystack_[3].value.val), (yystack_[1].value.val));
          }
#line 1672 "apr_parser.cc"
          break;

          case 135:
#line 553 "aprepro.yy"
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
#line 1698 "apr_parser.cc"
          break;

          case 136:
#line 575 "aprepro.yy"
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
#line 1722 "apr_parser.cc"
          break;

#line 1726 "apr_parser.cc"

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
      -34,  2,    -34, -3,   305,  -34,  -34,  -34,  -34,  -34,  -13,  68,   4,    312,  28,   47,
      53,   62,   67,  421,  421,  -34,  421,  376,  421,  55,   181,  50,   -16,  39,   1142, 376,
      421,  421,  421, 421,  421,  -34,  -34,  376,  421,  421,  421,  421,  421,  -34,  -34,  376,
      421,  421,  421, 421,  421,  421,  -34,  -34,  421,  421,  376,  221,  361,  376,  616,  646,
      95,   48,   421, 107,  1190, 880,  1094, 82,   -34,  82,   82,   -34,  -34,  -34,  -34,  -34,
      -34,  -34,  -34, 421,  421,  421,  -34,  376,  376,  421,  376,  -34,  421,  421,  421,  421,
      421,  421,  421, -34,  421,  421,  421,  421,  421,  421,  421,  421,  421,  421,  421,  376,
      421,  421,  -33, 1190, 1209, 1225, 1225, 1225, 1225, 1225, -33,  1190, 1209, 1225, 1225, 1225,
      1225, 1225, -33, 1190, 1209, 1190, 1225, 1225, 1225, 1225, 1225, 1225, 1190, 1225, 615,  -33,
      1190, 1209, -34, 72,   240,  645,  -34,  262,  248,  675,  288,  437,  705,  421,  421,  421,
      421,  82,   -34, -34,  421,  -34,  1156, 1176, 97,   1225, -34,  1,    1209, 1,    1241, -34,
      1241, 81,   81,  81,   81,   81,   81,   -34,  1256, 1270, 80,   80,   80,   80,   80,   80,
      161,  161,  82,  -34,  82,   82,   82,   421,  108,  -34,  421,  -34,  421,  -34,  -34,  421,
      -34,  421,  -34, -34,  421,  -34,  421,  -34,  1225, 1225, 1225, 1225, 82,   421,  421,  1119,
      421,  447,  907, 500,  586,  472,  137,  934,  506,  961,  735,  1190, 1225, 109,  1225, 421,
      -34,  -34,  -34, 421,  -34,  421,  421,  -34,  421,  -34,  -34,  -34,  421,  -34,  421,  528,
      988,  765,  824, 534,  478,  1015, 1225, -34,  -34,  421,  -34,  421,  -34,  421,  -34,  -34,
      795,  1042, 407, 421,  -34,  -34,  421,  556,  853,  562,  -34,  421,  -34,  1069, -34};

  const unsigned char Parser::yydefact_[] = {
      2,   0,   1,   0,   0,   4,   3,   9,   67,  47,  95,  70,  48,  71,  49,  28, 0,   0,
      0,   0,   0,   8,   0,   0,   0,   0,   0,   131, 0,   0,   0,   0,   0,   0,  0,   0,
      0,   98,  99,  0,   0,   0,   0,   0,   0,   74,  75,  0,   0,   0,   0,   0,  0,   0,
      86,  87,  0,   0,   0,   0,   0,   0,   95,  70,  48,  0,   0,   131, 0,   0,  0,   127,
      41,  126, 12,  68,  96,  72,  84,  69,  97,  73,  85,  0,   0,   0,   7,   0,  0,   0,
      0,   6,   0,   0,   0,   0,   0,   0,   0,   5,   0,   0,   0,   0,   0,   0,  0,   0,
      0,   0,   0,   0,   0,   0,   39,  50,  100, 101, 102, 103, 104, 105, 37,  52, 76,  79,
      80,  81,  82,  83,  36,  51,  77,  55,  88,  90,  91,  92,  93,  94,  54,  89, 0,   38,
      53,  78,  106, 0,   0,   0,   57,  0,   0,   0,   0,   0,   0,   0,   0,   0,  0,   126,
      21,  129, 0,   130, 0,   0,   19,  0,   20,  40,  0,   42,  44,  46,  43,  22, 23,  24,
      25,  26,  27,  60,  17,  18,  10,  11,  13,  14,  15,  16,  121, 122, 124, 45, 123, 125,
      128, 0,   133, 109, 0,   108, 0,   107, 59,  0,   56,  0,   58,  35,  0,   29, 0,   34,
      100, 76,  77,  78,  123, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,   0,
      0,   66,  132, 134, 135, 0,   112, 110, 111, 0,   114, 0,   0,   65,  0,   61, 31,  30,
      0,   33,  0,   0,   0,   0,   0,   0,   0,   0,   136, 113, 116, 0,   115, 0,  64,  0,
      63,  32,  0,   0,   0,   0,   118, 117, 0,   0,   0,   0,   119, 0,   62,  0,  120};

  const signed char Parser::yypgoto_[] = {-34, -34, -34, -18, 103, 85, -4};

  const short Parser::yydefgoto_[] = {-1, 1, 6, 27, 28, 68, 169};

  const unsigned short Parser::yytable_[] = {
      30,  67,  2,   3,   86,  87,  88,  89,  90,  31,  32,  33,  34,  35,  36,  69,  70,  7,   71,
      73,  74,  4,   87,  88,  89,  90,  47,  116, 117, 118, 119, 120, 121, 37,  38,  124, 125, 126,
      127, 128, 129, 89,  90,  132, 134, 135, 136, 137, 138, 139, 56,  5,   141, 142, 145, 149, 153,
      156, 75,  91,  76,  77,  161, 78,  57,  57,  168, 170, 59,  58,  160, 92,  93,  94,  95,  96,
      97,  60,  83,  167, 84,  85,  61,  172, 172, 174, 176, 98,  201, 29,  39,  40,  41,  42,  43,
      44,  184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 196, 197, 198, 87,  88,  89,  90,
      45,  46,  115, 159, 108, 109, 110, 164, 112, 162, 123, 113, 72,  113, 85,  98,  224, 254, 131,
      133, 114, 83,  0,   84,  85,  0,   0,   140, 122, 144, 148, 152, 155, 0,   0,   0,   130, 248,
      0,   216, 217, 218, 219, 0,   0,   0,   220, 143, 147, 151, 154, 0,   0,   0,   166, 92,  93,
      94,  95,  96,  97,  0,   0,   177, 178, 179, 180, 181, 182, 183, 79,  98,  80,  81,  0,   82,
      171, 173, 0,   175, 0,   223, 0,   0,   226, 0,   228, 110, 164, 112, 0,   231, 113, 0,   233,
      0,   234, 0,   0,   0,   195, 0,   0,   0,   236, 0,   238, 0,   0,   0,   8,   9,   10,  11,
      12,  13,  14,  15,  16,  17,  18,  256, 19,  146, 20,  257, 0,   258, 115, 123, 131, 144, 0,
      0,   261, 0,   262, 0,   0,   0,   202, 0,   203, 0,   0,   22,  23,  272, 207, 273, 208, 24,
      0,   25,  26,  0,   0,   280, 92,  93,  94,  95,  96,  97,  206, 285, 92,  93,  94,  95,  96,
      97,  0,   225, 98,  227, 0,   0,   229, 0,   230, 0,   98,  232, 0,   0,   87,  88,  89,  90,
      211, 0,   235, 0,   8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  0,   20,
      0,   255, 21,  87,  88,  89,  90,  0,   259, 0,   260, 48,  49,  50,  51,  52,  53,  0,   0,
      0,   22,  23,  0,   0,   0,   0,   24,  0,   25,  26,  0,   274, 0,   0,   0,   54,  55,  279,
      0,   0,   281, 8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  150, 20,  8,
      9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  0,   20,  0,   0,   0,   0,   0,
      22,  23,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,   0,   0,   22,  23,  0,   0,
      0,   0,   24,  278, 25,  26,  8,   9,   62,  63,  64,  13,  14,  65,  16,  17,  0,   0,   19,
      0,   20,  92,  93,  94,  95,  96,  97,  0,   0,   0,   0,   0,   0,   212, 0,   213, 0,   98,
      0,   0,   0,   22,  66,  239, 0,   240, 0,   24,  0,   25,  26,  92,  93,  94,  95,  96,  97,
      0,   0,   0,   0,   92,  93,  94,  95,  96,  97,  98,  246, 0,   247, 0,   0,   0,   269, 0,
      270, 98,  0,   0,   0,   0,   0,   0,   0,   0,   92,  93,  94,  95,  96,  97,  92,  93,  94,
      95,  96,  97,  242, 0,   0,   0,   98,  0,   250, 0,   0,   0,   98,  0,   0,   0,   0,   0,
      92,  93,  94,  95,  96,  97,  92,  93,  94,  95,  96,  97,  263, 0,   0,   0,   98,  0,   268,
      0,   0,   0,   98,  0,   0,   0,   0,   0,   92,  93,  94,  95,  96,  97,  92,  93,  94,  95,
      96,  97,  282, 0,   0,   0,   98,  0,   284, 0,   0,   0,   98,  0,   0,   0,   0,   0,   92,
      93,  94,  95,  96,  97,  92,  93,  94,  95,  96,  97,  243, 0,   244, 0,   98,  0,   0,   245,
      0,   0,   98,  0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
      164, 112, 199, 0,   113, 0,   200, 0,   0,   0,   0,   157, 32,  33,  34,  35,  36,  0,   100,
      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 0,   204, 113, 205, 37,  38,  0,
      0,   0,   0,   158, 40,  41,  42,  43,  44,  0,   100, 101, 102, 103, 104, 105, 106, 107, 108,
      109, 110, 111, 112, 0,   209, 113, 210, 45,  46,  0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 0,   214, 113, 215,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105,
      106, 107, 108, 109, 110, 111, 112, 0,   252, 113, 253, 0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 0,
      265, 113, 266, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102,
      103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 0,   275, 113, 276, 0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
      164, 112, 267, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100,
      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 283, 0,   113, 0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
      110, 164, 112, 163, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101,
      102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 241, 0,   113, 0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112,
      249, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104,
      105, 106, 107, 108, 109, 110, 164, 112, 251, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 264, 0,   113,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107,
      108, 109, 110, 164, 112, 271, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 277, 0,   113, 0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
      164, 112, 286, 0,   113, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102,
      103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 165, 0,   113, 0,   0,   0,   0,   0,   0,
      0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 237, 0,   113,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
      110, 164, 112, 99,  0,   113, 0,   0,   0,   0,   0,   0,   0,   100, 101, 102, 103, 104, 105,
      106, 107, 108, 109, 110, 111, 112, 221, 0,   113, 92,  93,  94,  95,  96,  97,  0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   98,  222, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
      110, 164, 112, 0,   0,   113, 92,  93,  94,  95,  96,  97,  0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   98,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 0,   0,
      113, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 164, 112, 0,   0,   113, 100, 101,
      102, 103, 104, 105, 106, 107, 0,   0,   0,   0,   112, 0,   0,   113, 101, 102, 103, 104, 105,
      106, 107, 108, 109, 110, 164, 112, 0,   0,   113, 102, 103, 104, 105, 106, 107, 108, 109, 110,
      164, 112, 0,   0,   113};

  const short Parser::yycheck_[] = {
      4,  19,  0,   1,   20,  38,  39,  40,  41,  22,  23,  24,  25,  26,  27,  19,  20,  20,  22,
      23, 24,  19,  38,  39,  40,  41,  22,  31,  32,  33,  34,  35,  36,  46,  47,  39,  40,  41,
      42, 43,  44,  40,  41,  47,  48,  49,  50,  51,  52,  53,  22,  49,  56,  57,  58,  59,  60,
      61, 3,   20,  5,   6,   66,  8,   17,  17,  84,  85,  15,  22,  22,  32,  33,  34,  35,  36,
      37, 15,  28,  83,  30,  31,  15,  87,  88,  89,  90,  48,  16,  4,   22,  23,  24,  25,  26,
      27, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 38,  39,  40,  41,
      46, 47,  31,  22,  38,  39,  40,  41,  42,  16,  39,  45,  23,  45,  31,  48,  22,  22,  47,
      48, 31,  28,  -1,  30,  31,  -1,  -1,  56,  39,  58,  59,  60,  61,  -1,  -1,  -1,  47,  14,
      -1, 157, 158, 159, 160, -1,  -1,  -1,  164, 58,  59,  60,  61,  -1,  -1,  -1,  83,  32,  33,
      34, 35,  36,  37,  -1,  -1,  92,  93,  94,  95,  96,  97,  98,  3,   48,  5,   6,   -1,  8,
      87, 88,  -1,  90,  -1,  199, -1,  -1,  202, -1,  204, 40,  41,  42,  -1,  209, 45,  -1,  212,
      -1, 214, -1,  -1,  -1,  111, -1,  -1,  -1,  222, -1,  224, -1,  -1,  -1,  3,   4,   5,   6,
      7,  8,   9,   10,  11,  12,  13,  239, 15,  16,  17,  243, -1,  245, 157, 158, 159, 160, -1,
      -1, 252, -1,  254, -1,  -1,  -1,  14,  -1,  16,  -1,  -1,  38,  39,  265, 14,  267, 16,  44,
      -1, 46,  47,  -1,  -1,  275, 32,  33,  34,  35,  36,  37,  16,  283, 32,  33,  34,  35,  36,
      37, -1,  202, 48,  204, -1,  -1,  207, -1,  209, -1,  48,  212, -1,  -1,  38,  39,  40,  41,
      16, -1,  221, -1,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,  -1,  17,
      -1, 239, 20,  38,  39,  40,  41,  -1,  246, -1,  248, 22,  23,  24,  25,  26,  27,  -1,  -1,
      -1, 38,  39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  269, -1,  -1,  -1,  46,  47,  275,
      -1, -1,  278, 3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,  16,  17,  3,
      4,  5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,  -1,  17,  -1,  -1,  -1,  -1,  -1,
      38, 39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,
      -1, -1,  44,  14,  46,  47,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,
      -1, 17,  32,  33,  34,  35,  36,  37,  -1,  -1,  -1,  -1,  -1,  -1,  14,  -1,  16,  -1,  48,
      -1, -1,  -1,  38,  39,  14,  -1,  16,  -1,  44,  -1,  46,  47,  32,  33,  34,  35,  36,  37,
      -1, -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  48,  14,  -1,  16,  -1,  -1,  -1,  14,  -1,
      16, 48,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  32,  33,  34,
      35, 36,  37,  16,  -1,  -1,  -1,  48,  -1,  16,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,
      32, 33,  34,  35,  36,  37,  32,  33,  34,  35,  36,  37,  16,  -1,  -1,  -1,  48,  -1,  16,
      -1, -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  32,  33,  34,  35,
      36, 37,  16,  -1,  -1,  -1,  48,  -1,  16,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  32,
      33, 34,  35,  36,  37,  32,  33,  34,  35,  36,  37,  14,  -1,  16,  -1,  48,  -1,  -1,  21,
      -1, -1,  48,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
      41, 42,  14,  -1,  45,  -1,  18,  -1,  -1,  -1,  -1,  22,  23,  24,  25,  26,  27,  -1,  30,
      31, 32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  46,  47,  -1,
      -1, -1,  -1,  22,  23,  24,  25,  26,  27,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,
      39, 40,  41,  42,  -1,  14,  45,  16,  46,  47,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1, -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,
      -1, -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,
      36, 37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1, -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,
      14, 45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,
      33, 34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,
      -1, -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
      41, 42,  14,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,
      31, 32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  14,  -1,  45,  -1,  -1,  -1,  -1,
      -1, -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
      40, 41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,
      32, 33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,
      -1, -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
      16, -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,
      35, 36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1, -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,
      -1, -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,
      38, 39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      30, 31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,
      -1, -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
      41, 42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,
      33, 34,  35,  36,  37,  38,  39,  40,  41,  42,  18,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,
      -1, -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  18,  -1,  45,
      -1, -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
      40, 41,  42,  20,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,
      36, 37,  38,  39,  40,  41,  42,  29,  -1,  45,  32,  33,  34,  35,  36,  37,  -1,  -1,  -1,
      -1, -1,  -1,  -1,  -1,  -1,  -1,  48,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
      40, 41,  42,  -1,  -1,  45,  32,  33,  34,  35,  36,  37,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1, -1,  -1,  48,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,
      45, 30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  30,  31,
      32, 33,  34,  35,  36,  37,  -1,  -1,  -1,  -1,  42,  -1,  -1,  45,  31,  32,  33,  34,  35,
      36, 37,  38,  39,  40,  41,  42,  -1,  -1,  45,  32,  33,  34,  35,  36,  37,  38,  39,  40,
      41, 42,  -1,  -1,  45};

  const unsigned char Parser::yystos_[] = {
      0,  51, 0,  1,  19, 49, 52, 20, 3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 15, 17, 20, 38,
      39, 44, 46, 47, 53, 54, 55, 56, 22, 23, 24, 25, 26, 27, 46, 47, 22, 23, 24, 25, 26, 27, 46,
      47, 22, 22, 23, 24, 25, 26, 27, 46, 47, 22, 17, 22, 15, 15, 15, 5,  6,  7,  10, 39, 53, 55,
      56, 56, 56, 54, 56, 56, 3,  5,  6,  8,  3,  5,  6,  8,  28, 30, 31, 20, 38, 39, 40, 41, 20,
      32, 33, 34, 35, 36, 37, 48, 20, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 45, 54,
      55, 56, 56, 56, 56, 56, 56, 54, 55, 56, 56, 56, 56, 56, 56, 54, 55, 56, 55, 56, 56, 56, 56,
      56, 56, 55, 56, 56, 54, 55, 56, 16, 54, 55, 56, 16, 54, 55, 56, 54, 55, 56, 22, 22, 22, 22,
      56, 16, 16, 41, 18, 55, 56, 53, 56, 53, 54, 56, 54, 56, 54, 56, 55, 55, 55, 55, 55, 55, 55,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 54, 56, 56, 56, 14, 18, 16, 14, 16, 14, 16, 16,
      14, 16, 14, 16, 16, 14, 16, 14, 16, 56, 56, 56, 56, 56, 29, 29, 56, 22, 55, 56, 55, 56, 55,
      55, 56, 55, 56, 56, 55, 56, 18, 56, 14, 16, 16, 16, 14, 16, 21, 14, 16, 14, 16, 16, 16, 14,
      16, 22, 55, 56, 56, 56, 55, 55, 56, 56, 16, 16, 14, 16, 14, 16, 14, 16, 16, 56, 56, 55, 14,
      16, 16, 14, 55, 56, 55, 16, 14, 16, 56, 16};

  const unsigned char Parser::yyr1_[] = {
      0,  50, 51, 51, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
      53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
      54, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56};

  const unsigned char Parser::yyr2_[] = {
      0, 2, 0, 2, 1, 3,  3,  3,  2,  2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
      1, 4, 6, 6, 8, 6,  4,  4,  3,  3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 1, 1, 1, 3, 3, 3, 3, 3, 3,
      4, 3, 4, 4, 3, 6,  12, 8,  8,  6, 5, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
      2, 2, 2, 2, 3, 3,  3,  3,  3,  3, 3, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 6, 6,
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
      0,   145, 145, 146, 149, 150, 157, 161, 162, 163, 166, 167, 168, 169, 170, 171, 172, 173,
      174, 175, 176, 177, 180, 181, 182, 183, 184, 185, 187, 188, 194, 200, 206, 212, 218, 224,
      230, 234, 238, 241, 243, 251, 253, 261, 262, 263, 264, 273, 274, 275, 276, 278, 281, 285,
      290, 291, 292, 298, 304, 310, 316, 317, 323, 329, 335, 341, 347, 349, 350, 351, 352, 353,
      354, 355, 356, 357, 358, 360, 363, 368, 369, 370, 371, 372, 377, 378, 379, 380, 381, 382,
      383, 384, 385, 386, 387, 389, 391, 394, 397, 400, 403, 405, 408, 411, 414, 417, 424, 431,
      438, 445, 452, 459, 466, 473, 480, 486, 492, 498, 504, 510, 516, 522, 523, 524, 525, 532,
      539, 540, 541, 544, 545, 548, 549, 550, 551, 552, 574};

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
#line 2541 "apr_parser.cc"

#line 597 "aprepro.yy"

void SEAMS::Parser::error(const std::string &m) { aprepro.error(m); }
