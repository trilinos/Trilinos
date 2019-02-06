// A Bison parser, made by GNU Bison 3.3.1.

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
#line 131 "aprepro.yy" // lalr1.cc:434

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
#if YYDEBUG

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

#else // !YYDEBUG

#define YYCDEBUG                                                                                   \
  if (false)                                                                                       \
  std::cerr
#define YY_SYMBOL_PRINT(Title, Symbol) YYUSE(Symbol)
#define YY_REDUCE_PRINT(Rule) static_cast<void>(0)
#define YY_STACK_PRINT() static_cast<void>(0)

#endif // !YYDEBUG

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
#if YYDEBUG
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

#if YYDEBUG
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

#if YYDEBUG
  std::ostream &Parser::debug_stream() const { return *yycdebug_; }

  void Parser::set_debug_stream(std::ostream &o) { yycdebug_ = &o; }

  Parser::debug_level_type Parser::debug_level() const { return yydebug_; }

  void Parser::set_debug_level(debug_level_type l) { yydebug_ = l; }
#endif // YYDEBUG

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
#line 150 "aprepro.yy" // lalr1.cc:919
          {
            if (echo)
              aprepro.lexer->LexerOutput("\n", 1);
          }
#line 636 "apr_parser.cc" // lalr1.cc:919
          break;

          case 5:
#line 151 "aprepro.yy" // lalr1.cc:919
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
#line 158 "aprepro.yy" // lalr1.cc:919
          {
            if (echo && (yystack_[1].value.string) != NULL) {
              aprepro.lexer->LexerOutput((yystack_[1].value.string),
                                         strlen((yystack_[1].value.string)));
            }
          }
#line 657 "apr_parser.cc" // lalr1.cc:919
          break;

          case 7:
#line 162 "aprepro.yy" // lalr1.cc:919
          {
          }
#line 663 "apr_parser.cc" // lalr1.cc:919
          break;

          case 8:
#line 163 "aprepro.yy" // lalr1.cc:919
          {
          }
#line 669 "apr_parser.cc" // lalr1.cc:919
          break;

          case 9:
#line 164 "aprepro.yy" // lalr1.cc:919
          {
            yyerrok;
          }
#line 675 "apr_parser.cc" // lalr1.cc:919
          break;

          case 10:
#line 167 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) < (yystack_[0].value.val);
          }
#line 681 "apr_parser.cc" // lalr1.cc:919
          break;

          case 11:
#line 168 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) > (yystack_[0].value.val);
          }
#line 687 "apr_parser.cc" // lalr1.cc:919
          break;

          case 12:
#line 169 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = !((yystack_[0].value.val));
          }
#line 693 "apr_parser.cc" // lalr1.cc:919
          break;

          case 13:
#line 170 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) <= (yystack_[0].value.val);
          }
#line 699 "apr_parser.cc" // lalr1.cc:919
          break;

          case 14:
#line 171 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) >= (yystack_[0].value.val);
          }
#line 705 "apr_parser.cc" // lalr1.cc:919
          break;

          case 15:
#line 172 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) == (yystack_[0].value.val);
          }
#line 711 "apr_parser.cc" // lalr1.cc:919
          break;

          case 16:
#line 173 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) != (yystack_[0].value.val);
          }
#line 717 "apr_parser.cc" // lalr1.cc:919
          break;

          case 17:
#line 174 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 723 "apr_parser.cc" // lalr1.cc:919
          break;

          case 18:
#line 175 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 729 "apr_parser.cc" // lalr1.cc:919
          break;

          case 19:
#line 176 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) || (yystack_[0].value.val);
          }
#line 735 "apr_parser.cc" // lalr1.cc:919
          break;

          case 20:
#line 177 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) && (yystack_[0].value.val);
          }
#line 741 "apr_parser.cc" // lalr1.cc:919
          break;

          case 21:
#line 178 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          }
#line 747 "apr_parser.cc" // lalr1.cc:919
          break;

          case 22:
#line 181 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) < 0 ? 1 : 0);
          }
#line 753 "apr_parser.cc" // lalr1.cc:919
          break;

          case 23:
#line 182 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) > 0 ? 1 : 0);
          }
#line 759 "apr_parser.cc" // lalr1.cc:919
          break;

          case 24:
#line 183 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) <= 0 ? 1 : 0);
          }
#line 765 "apr_parser.cc" // lalr1.cc:919
          break;

          case 25:
#line 184 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) >= 0 ? 1 : 0);
          }
#line 771 "apr_parser.cc" // lalr1.cc:919
          break;

          case 26:
#line 185 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) == 0 ? 1 : 0);
          }
#line 777 "apr_parser.cc" // lalr1.cc:919
          break;

          case 27:
#line 186 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                (strcmp((yystack_[2].value.string), (yystack_[0].value.string)) != 0 ? 1 : 0);
          }
#line 783 "apr_parser.cc" // lalr1.cc:919
          break;

          case 28:
#line 188 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) = (yystack_[0].value.tptr)->value.avar;
          }
#line 789 "apr_parser.cc" // lalr1.cc:919
          break;

          case 29:
#line 189 "aprepro.yy" // lalr1.cc:919
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
#line 195 "aprepro.yy" // lalr1.cc:919
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
#line 201 "aprepro.yy" // lalr1.cc:919
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
#line 207 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.arrfnct_dd == NULL))
              (yylhs.value.arrval) = (*((yystack_[5].value.tptr)->value.arrfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              yyerrok;
          }
#line 833 "apr_parser.cc" // lalr1.cc:919
          break;

          case 33:
#line 213 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_d == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_d))((yystack_[1].value.val));
            else
              yyerrok;
          }
#line 844 "apr_parser.cc" // lalr1.cc:919
          break;

          case 34:
#line 219 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.arrfnct_a == NULL))
              (yylhs.value.arrval) =
                  (*((yystack_[3].value.tptr)->value.arrfnct_a))((yystack_[1].value.arrval));
            else
              yyerrok;
          }
#line 855 "apr_parser.cc" // lalr1.cc:919
          break;

          case 35:
#line 225 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) = (yystack_[0].value.arrval);
            delete (yystack_[2].value.tptr)->value.avar;
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 863 "apr_parser.cc" // lalr1.cc:919
          break;

          case 36:
#line 228 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval)                 = (yystack_[0].value.arrval);
            (yystack_[2].value.tptr)->value.avar = (yystack_[0].value.arrval);
            set_type(aprepro, (yystack_[2].value.tptr), token::AVAR);
          }
#line 870 "apr_parser.cc" // lalr1.cc:919
          break;

          case 37:
#line 230 "aprepro.yy" // lalr1.cc:919
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
#line 883 "apr_parser.cc" // lalr1.cc:919
          break;

          case 38:
#line 238 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), -1.0);
          }
#line 889 "apr_parser.cc" // lalr1.cc:919
          break;

          case 39:
#line 240 "aprepro.yy" // lalr1.cc:919
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
#line 902 "apr_parser.cc" // lalr1.cc:919
          break;

          case 40:
#line 248 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) = array_scale((yystack_[2].value.arrval), (yystack_[0].value.val));
          }
#line 908 "apr_parser.cc" // lalr1.cc:919
          break;

          case 41:
#line 249 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) =
                array_scale((yystack_[2].value.arrval), 1.0 / (yystack_[0].value.val));
          }
#line 914 "apr_parser.cc" // lalr1.cc:919
          break;

          case 42:
#line 250 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.arrval) = array_scale((yystack_[0].value.arrval), (yystack_[2].value.val));
          }
#line 920 "apr_parser.cc" // lalr1.cc:919
          break;

          case 43:
#line 251 "aprepro.yy" // lalr1.cc:919
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
#line 933 "apr_parser.cc" // lalr1.cc:919
          break;

          case 44:
#line 260 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string) = (yystack_[0].value.string);
          }
#line 939 "apr_parser.cc" // lalr1.cc:919
          break;

          case 45:
#line 261 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          }
#line 945 "apr_parser.cc" // lalr1.cc:919
          break;

          case 46:
#line 262 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string) = (char *)(yystack_[0].value.tptr)->value.svar.c_str();
          }
#line 951 "apr_parser.cc" // lalr1.cc:919
          break;

          case 47:
#line 263 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            set_type(aprepro, (yystack_[2].value.tptr), Parser::token::SVAR);
          }
#line 958 "apr_parser.cc" // lalr1.cc:919
          break;

          case 48:
#line 265 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          }
#line 966 "apr_parser.cc" // lalr1.cc:919
          break;

          case 49:
#line 268 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string)                 = (yystack_[0].value.string);
            (yystack_[2].value.tptr)->value.svar = (yystack_[0].value.string);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::SVAR);
          }
#line 975 "apr_parser.cc" // lalr1.cc:919
          break;

          case 50:
#line 272 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string) = (char *)(yystack_[2].value.tptr)->value.svar.c_str();
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 981 "apr_parser.cc" // lalr1.cc:919
          break;

          case 51:
#line 273 "aprepro.yy" // lalr1.cc:919
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          }
#line 987 "apr_parser.cc" // lalr1.cc:919
          break;

          case 52:
#line 274 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_c == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_c))(
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 998 "apr_parser.cc" // lalr1.cc:919
          break;

          case 53:
#line 280 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.strfnct == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[2].value.tptr)->value.strfnct))();
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1009 "apr_parser.cc" // lalr1.cc:919
          break;

          case 54:
#line 286 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_d == NULL))
              (yylhs.value.string) =
                  (char *)(*((yystack_[3].value.tptr)->value.strfnct_d))((yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1020 "apr_parser.cc" // lalr1.cc:919
          break;

          case 55:
#line 292 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.strfnct_a == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[3].value.tptr)->value.strfnct_a))(
                  (yystack_[1].value.arrval));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1031 "apr_parser.cc" // lalr1.cc:919
          break;

          case 56:
#line 298 "aprepro.yy" // lalr1.cc:919
          {
            concat_string((yystack_[2].value.string), (yystack_[0].value.string),
                          &(yylhs.value.string));
          }
#line 1037 "apr_parser.cc" // lalr1.cc:919
          break;

          case 57:
#line 299 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_dd == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1048 "apr_parser.cc" // lalr1.cc:919
          break;

          case 58:
#line 305 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.strfnct_dcccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[11].value.tptr)->value.strfnct_dcccc))(
                  (yystack_[9].value.val), (yystack_[7].value.string), (yystack_[5].value.string),
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1059 "apr_parser.cc" // lalr1.cc:919
          break;

          case 59:
#line 311 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_dcc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_dcc))(
                  (yystack_[5].value.val), (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1070 "apr_parser.cc" // lalr1.cc:919
          break;

          case 60:
#line 317 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.strfnct_ccc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[7].value.tptr)->value.strfnct_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1081 "apr_parser.cc" // lalr1.cc:919
          break;

          case 61:
#line 323 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.strfnct_cc == NULL))
              (yylhs.value.string) = (char *)(*((yystack_[5].value.tptr)->value.strfnct_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.string) = (char *)"";
          }
#line 1092 "apr_parser.cc" // lalr1.cc:919
          break;

          case 62:
#line 329 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.string) = ((yystack_[4].value.val)) ? ((yystack_[2].value.string))
                                                             : ((yystack_[0].value.string));
          }
#line 1098 "apr_parser.cc" // lalr1.cc:919
          break;

          case 63:
#line 331 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          }
#line 1104 "apr_parser.cc" // lalr1.cc:919
          break;

          case 64:
#line 332 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.val) + 1;
          }
#line 1110 "apr_parser.cc" // lalr1.cc:919
          break;

          case 65:
#line 333 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.val) - 1;
          }
#line 1116 "apr_parser.cc" // lalr1.cc:919
          break;

          case 66:
#line 334 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          }
#line 1122 "apr_parser.cc" // lalr1.cc:919
          break;

          case 67:
#line 335 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
          }
#line 1128 "apr_parser.cc" // lalr1.cc:919
          break;

          case 68:
#line 336 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
          }
#line 1134 "apr_parser.cc" // lalr1.cc:919
          break;

          case 69:
#line 337 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
          }
#line 1140 "apr_parser.cc" // lalr1.cc:919
          break;

          case 70:
#line 338 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
          }
#line 1146 "apr_parser.cc" // lalr1.cc:919
          break;

          case 71:
#line 339 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
          }
#line 1152 "apr_parser.cc" // lalr1.cc:919
          break;

          case 72:
#line 340 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
          }
#line 1159 "apr_parser.cc" // lalr1.cc:919
          break;

          case 73:
#line 342 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            redefined_warning(aprepro, (yystack_[2].value.tptr));
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1167 "apr_parser.cc" // lalr1.cc:919
          break;

          case 74:
#line 345 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1173 "apr_parser.cc" // lalr1.cc:919
          break;

          case 75:
#line 346 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1179 "apr_parser.cc" // lalr1.cc:919
          break;

          case 76:
#line 347 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1185 "apr_parser.cc" // lalr1.cc:919
          break;

          case 77:
#line 348 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
          }
#line 1191 "apr_parser.cc" // lalr1.cc:919
          break;

          case 78:
#line 349 "aprepro.yy" // lalr1.cc:919
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            SEAMS::math_error(aprepro, "Power");
          }
#line 1201 "apr_parser.cc" // lalr1.cc:919
          break;

          case 79:
#line 354 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          }
#line 1207 "apr_parser.cc" // lalr1.cc:919
          break;

          case 80:
#line 355 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[0].value.tptr));
          }
#line 1213 "apr_parser.cc" // lalr1.cc:919
          break;

          case 81:
#line 356 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          }
#line 1219 "apr_parser.cc" // lalr1.cc:919
          break;

          case 82:
#line 357 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[1].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[1].value.tptr));
          }
#line 1225 "apr_parser.cc" // lalr1.cc:919
          break;

          case 83:
#line 358 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1231 "apr_parser.cc" // lalr1.cc:919
          break;

          case 84:
#line 359 "aprepro.yy" // lalr1.cc:919
          {
            immutable_modify(aprepro, (yystack_[2].value.tptr));
            YYERROR;
          }
#line 1237 "apr_parser.cc" // lalr1.cc:919
          break;

          case 85:
#line 360 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1243 "apr_parser.cc" // lalr1.cc:919
          break;

          case 86:
#line 361 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1249 "apr_parser.cc" // lalr1.cc:919
          break;

          case 87:
#line 362 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1255 "apr_parser.cc" // lalr1.cc:919
          break;

          case 88:
#line 363 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1261 "apr_parser.cc" // lalr1.cc:919
          break;

          case 89:
#line 364 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            immutable_modify(aprepro, (yystack_[2].value.tptr));
          }
#line 1267 "apr_parser.cc" // lalr1.cc:919
          break;

          case 90:
#line 366 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.tptr)->value.var;
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1274 "apr_parser.cc" // lalr1.cc:919
          break;

          case 91:
#line 368 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ++((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1282 "apr_parser.cc" // lalr1.cc:919
          break;

          case 92:
#line 371 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = --((yystack_[0].value.tptr)->value.var);
            set_type(aprepro, (yystack_[0].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[0].value.tptr)->name);
          }
#line 1290 "apr_parser.cc" // lalr1.cc:919
          break;

          case 93:
#line 374 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)++;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          }
#line 1298 "apr_parser.cc" // lalr1.cc:919
          break;

          case 94:
#line 377 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ((yystack_[1].value.tptr)->value.var)--;
            set_type(aprepro, (yystack_[1].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[1].value.tptr)->name);
          }
#line 1306 "apr_parser.cc" // lalr1.cc:919
          break;

          case 95:
#line 380 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val)                   = (yystack_[0].value.val);
            (yystack_[2].value.tptr)->value.var = (yystack_[0].value.val);
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
          }
#line 1313 "apr_parser.cc" // lalr1.cc:919
          break;

          case 96:
#line 382 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var += (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1321 "apr_parser.cc" // lalr1.cc:919
          break;

          case 97:
#line 385 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var -= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1329 "apr_parser.cc" // lalr1.cc:919
          break;

          case 98:
#line 388 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var *= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1337 "apr_parser.cc" // lalr1.cc:919
          break;

          case 99:
#line 391 "aprepro.yy" // lalr1.cc:919
          {
            (yystack_[2].value.tptr)->value.var /= (yystack_[0].value.val);
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1345 "apr_parser.cc" // lalr1.cc:919
          break;

          case 100:
#line 394 "aprepro.yy" // lalr1.cc:919
          {
            reset_error();
            (yystack_[2].value.tptr)->value.var =
                std::pow((yystack_[2].value.tptr)->value.var, (yystack_[0].value.val));
            (yylhs.value.val) = (yystack_[2].value.tptr)->value.var;
            set_type(aprepro, (yystack_[2].value.tptr), token::VAR);
            SEAMS::math_error(aprepro, "Power");
            undefined_error(aprepro, (yystack_[2].value.tptr)->name);
          }
#line 1356 "apr_parser.cc" // lalr1.cc:919
          break;

          case 101:
#line 401 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[2].value.tptr),
                          (yystack_[2].value.tptr)->value.fnctptr == NULL))
              (yylhs.value.val) = (*((yystack_[2].value.tptr)->value.fnctptr))();
            else
              (yylhs.value.val) = 0.0;
          }
#line 1367 "apr_parser.cc" // lalr1.cc:919
          break;

          case 102:
#line 408 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_d == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_d))((yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1378 "apr_parser.cc" // lalr1.cc:919
          break;

          case 103:
#line 415 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_c == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_c))((yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1389 "apr_parser.cc" // lalr1.cc:919
          break;

          case 104:
#line 422 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[3].value.tptr),
                          (yystack_[3].value.tptr)->value.fnctptr_a == NULL))
              (yylhs.value.val) =
                  (*((yystack_[3].value.tptr)->value.fnctptr_a))((yystack_[1].value.arrval));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1400 "apr_parser.cc" // lalr1.cc:919
          break;

          case 105:
#line 429 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cd))(
                  (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1411 "apr_parser.cc" // lalr1.cc:919
          break;

          case 106:
#line 436 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dc))(
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1422 "apr_parser.cc" // lalr1.cc:919
          break;

          case 107:
#line 443 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_cc == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_cc))(
                  (yystack_[3].value.string), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1433 "apr_parser.cc" // lalr1.cc:919
          break;

          case 108:
#line 450 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccc == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccc))(
                  (yystack_[5].value.string), (yystack_[3].value.string),
                  (yystack_[1].value.string));
            else
              yyerrok;
          }
#line 1444 "apr_parser.cc" // lalr1.cc:919
          break;

          case 109:
#line 457 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[5].value.tptr),
                          (yystack_[5].value.tptr)->value.fnctptr_dd == NULL))
              (yylhs.value.val) = (*((yystack_[5].value.tptr)->value.fnctptr_dd))(
                  (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1455 "apr_parser.cc" // lalr1.cc:919
          break;

          case 110:
#line 463 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ddd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ddd))(
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1466 "apr_parser.cc" // lalr1.cc:919
          break;

          case 111:
#line 469 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[7].value.tptr),
                          (yystack_[7].value.tptr)->value.fnctptr_ccd == NULL))
              (yylhs.value.val) = (*((yystack_[7].value.tptr)->value.fnctptr_ccd))(
                  (yystack_[5].value.string), (yystack_[3].value.string), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1477 "apr_parser.cc" // lalr1.cc:919
          break;

          case 112:
#line 475 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[9].value.tptr),
                          (yystack_[9].value.tptr)->value.fnctptr_dddd == NULL))
              (yylhs.value.val) = (*((yystack_[9].value.tptr)->value.fnctptr_dddd))(
                  (yystack_[7].value.val), (yystack_[5].value.val), (yystack_[3].value.val),
                  (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1488 "apr_parser.cc" // lalr1.cc:919
          break;

          case 113:
#line 481 "aprepro.yy" // lalr1.cc:919
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
#line 487 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[11].value.tptr),
                          (yystack_[11].value.tptr)->value.fnctptr_ddddc == NULL))
              (yylhs.value.val) = (*((yystack_[11].value.tptr)->value.fnctptr_ddddc))(
                  (yystack_[9].value.val), (yystack_[7].value.val), (yystack_[5].value.val),
                  (yystack_[3].value.val), (yystack_[1].value.string));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1510 "apr_parser.cc" // lalr1.cc:919
          break;

          case 115:
#line 493 "aprepro.yy" // lalr1.cc:919
          {
            if (arg_check((yystack_[13].value.tptr),
                          (yystack_[13].value.tptr)->value.fnctptr_dddddd == NULL))
              (yylhs.value.val) = (*((yystack_[13].value.tptr)->value.fnctptr_dddddd))(
                  (yystack_[11].value.val), (yystack_[9].value.val), (yystack_[7].value.val),
                  (yystack_[5].value.val), (yystack_[3].value.val), (yystack_[1].value.val));
            else
              (yylhs.value.val) = 0.0;
          }
#line 1521 "apr_parser.cc" // lalr1.cc:919
          break;

          case 116:
#line 499 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) + (yystack_[0].value.val);
          }
#line 1527 "apr_parser.cc" // lalr1.cc:919
          break;

          case 117:
#line 500 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) - (yystack_[0].value.val);
          }
#line 1533 "apr_parser.cc" // lalr1.cc:919
          break;

          case 118:
#line 501 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[2].value.val) * (yystack_[0].value.val);
          }
#line 1539 "apr_parser.cc" // lalr1.cc:919
          break;

          case 119:
#line 502 "aprepro.yy" // lalr1.cc:919
          {
            if ((yystack_[0].value.val) == 0.) {
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (yystack_[2].value.val) / (yystack_[0].value.val);
          }
#line 1551 "apr_parser.cc" // lalr1.cc:919
          break;

          case 120:
#line 509 "aprepro.yy" // lalr1.cc:919
          {
            if ((yystack_[0].value.val) == 0.) {
              yyerror(aprepro, "Zero divisor");
              yyerrok;
            }
            else
              (yylhs.value.val) = (int)(yystack_[2].value.val) % (int)(yystack_[0].value.val);
          }
#line 1563 "apr_parser.cc" // lalr1.cc:919
          break;

          case 121:
#line 516 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = -(yystack_[0].value.val);
          }
#line 1569 "apr_parser.cc" // lalr1.cc:919
          break;

          case 122:
#line 517 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[0].value.val);
          }
#line 1575 "apr_parser.cc" // lalr1.cc:919
          break;

          case 123:
#line 518 "aprepro.yy" // lalr1.cc:919
          {
            reset_error();
            (yylhs.value.val) = std::pow((yystack_[2].value.val), (yystack_[0].value.val));
            SEAMS::math_error(aprepro, "Power");
          }
#line 1583 "apr_parser.cc" // lalr1.cc:919
          break;

          case 124:
#line 521 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = (yystack_[1].value.val);
          }
#line 1589 "apr_parser.cc" // lalr1.cc:919
          break;

          case 125:
#line 522 "aprepro.yy" // lalr1.cc:919
          {
            reset_error();
            (yylhs.value.val) =
                (double)((yystack_[1].value.val) < 0 ? -floor(-((yystack_[1].value.val)))
                                                     : floor((yystack_[1].value.val)));
            SEAMS::math_error(aprepro, "floor (int)");
          }
#line 1597 "apr_parser.cc" // lalr1.cc:919
          break;

          case 126:
#line 525 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = ((yystack_[0].value.val)) ? 1 : 0;
          }
#line 1603 "apr_parser.cc" // lalr1.cc:919
          break;

          case 127:
#line 526 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                ((yystack_[4].value.val)) ? ((yystack_[2].value.val)) : ((yystack_[0].value.val));
          }
#line 1609 "apr_parser.cc" // lalr1.cc:919
          break;

          case 128:
#line 527 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) =
                array_value((yystack_[3].value.tptr)->value.avar, (yystack_[1].value.val), 0);
          }
#line 1615 "apr_parser.cc" // lalr1.cc:919
          break;

          case 129:
#line 528 "aprepro.yy" // lalr1.cc:919
          {
            (yylhs.value.val) = array_value((yystack_[5].value.tptr)->value.avar,
                                            (yystack_[3].value.val), (yystack_[1].value.val));
          }
#line 1621 "apr_parser.cc" // lalr1.cc:919
          break;

          case 130:
#line 530 "aprepro.yy" // lalr1.cc:919
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
#line 1647 "apr_parser.cc" // lalr1.cc:919
          break;

          case 131:
#line 552 "aprepro.yy" // lalr1.cc:919
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
#line 1671 "apr_parser.cc" // lalr1.cc:919
          break;

#line 1675 "apr_parser.cc" // lalr1.cc:919
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
      -34,  2,    -34,  -3,   315,  -34,  -34,  -34, -34,  -34,  -13,  585,  4,    615,  28,
      45,   44,   46,   53,   420,  420,  -34,  420, 207,  420,  111,  179,  48,   -16,  38,
      1081, 207,  420,  420,  420,  420,  420,  -34, -34,  420,  420,  420,  420,  420,  420,
      -34,  -34,  420,  420,  420,  420,  420,  420, 420,  -34,  -34,  420,  420,  207,  360,
      375,  207,  645,  52,   420,  93,   268,  819, 1033, 18,   -34,  18,   18,   -34,  -34,
      -34,  -34,  -34,  -34,  -34,  -34,  420,  420, 420,  -34,  207,  207,  420,  207,  -34,
      420,  420,  420,  420,  420,  420,  420,  -34, 420,  420,  420,  420,  420,  420,  420,
      420,  420,  420,  420,  207,  420,  420,  -33, 268,  1131, 1147, 1147, 1147, 1147, 1147,
      268,  1147, 1147, 1147, 1147, 1147, 1147, 268, 1147, 268,  1147, 1147, 1147, 1147, 1147,
      1147, 268,  1147, 584,  -33,  1131, -34,  50,  99,   614,  -34,  154,  225,  644,  247,
      234,  674,  420,  18,   -34,  -34,  420,  -34, 1095, 1115, 49,   1147, -34,  1,    1,
      1163, -34,  1163, 39,   39,   39,   39,   39,  39,   -34,  1178, 1192, 307,  307,  307,
      307,  307,  307,  190,  190,  18,   -34,  18,  18,   18,   420,  70,   -34,  420,  -34,
      420,  -34,  -34,  420,  -34,  420,  -34,  -34, 420,  -34,  420,  -34,  1147, 18,   420,
      420,  1058, 420,  261,  846,  469,  555,  436, 131,  873,  475,  900,  927,  268,  1147,
      71,   1147, 420,  -34,  -34,  -34,  420,  -34, 420,  420,  -34,  420,  -34,  -34,  -34,
      -34,  420,  497,  954,  704,  763,  503,  446, 1147, -34,  -34,  420,  -34,  420,  -34,
      420,  -34,  734,  981,  406,  420,  -34,  -34, 420,  525,  792,  531,  -34,  420,  -34,
      1008, -34};

  const unsigned char Parser::yydefact_[] = {
      2,   0,   1,   0,   0,   4,   3,   9,  63,  44, 90, 66, 45,  67,  46,  28,  0,   0,   0,
      0,   0,   8,   0,   0,   0,   0,   0,  126, 0,  0,  0,  0,   0,   0,   0,   0,   0,   93,
      94,  0,   0,   0,   0,   0,   0,   70, 71,  0,  0,  0,  0,   0,   0,   0,   81,  82,  0,
      0,   0,   0,   0,   0,   90,  0,   0,  126, 0,  0,  0,  122, 38,  121, 12,  64,  91,  68,
      79,  65,  92,  69,  80,  0,   0,   0,  7,   0,  0,  0,  0,   6,   0,   0,   0,   0,   0,
      0,   0,   5,   0,   0,   0,   0,   0,  0,   0,  0,  0,  0,   0,   0,   0,   0,   36,  47,
      95,  96,  97,  98,  99,  100, 49,  72, 74,  75, 76, 77, 78,  48,  73,  51,  83,  85,  86,
      87,  88,  89,  50,  84,  0,   35,  0,  101, 0,  0,  0,  53,  0,   0,   0,   0,   0,   0,
      0,   121, 21,  124, 0,   125, 0,   0,  19,  0,  20, 37, 39,  41,  43,  40,  22,  23,  24,
      25,  26,  27,  56,  17,  18,  10,  11, 13,  14, 15, 16, 116, 117, 119, 42,  118, 120, 123,
      0,   128, 104, 0,   103, 0,   102, 55, 0,   52, 0,  54, 34,  0,   29,  0,   33,  95,  118,
      0,   0,   0,   0,   0,   0,   0,   0,  0,   0,  0,  0,  0,   0,   62,  127, 129, 130, 0,
      107, 105, 106, 0,   109, 0,   0,   61, 0,   57, 31, 30, 32,  0,   0,   0,   0,   0,   0,
      0,   131, 108, 111, 0,   110, 0,   60, 0,   59, 0,  0,  0,   0,   113, 112, 0,   0,   0,
      0,   114, 0,   58,  0,   115};

  const signed char Parser::yypgoto_[] = {-34, -34, -34, -18, 95, 81, -4};

  const short Parser::yydefgoto_[] = {-1, 1, 6, 27, 28, 66, 161};

  const unsigned short Parser::yytable_[] = {
      30,  65,  2,   3,   84,  85,  86,  87,  88,  31,  32,  33,  34,  35,  36,  67,  68,  7,   69,
      71,  72,  4,   85,  86,  87,  88,  47,  114, 115, 116, 117, 118, 119, 37,  38,  121, 122, 123,
      124, 125, 126, 87,  88,  128, 130, 131, 132, 133, 134, 135, 56,  5,   137, 138, 140, 144, 148,
      151, 89,  59,  153, 60,  57,  111, 160, 162, 192, 58,  61,  57,  90,  91,  92,  93,  94,  95,
      81,  159, 82,  83,  83,  140, 140, 165, 167, 29,  96,  96,  85,  86,  87,  88,  212, 241, 175,
      176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 187, 188, 189, 0,   154, 0,   0,   113, 193,
      73,  194, 74,  75,  70,  76,  120, 81,  0,   82,  83,  0,   112, 0,   127, 129, 0,   90,  91,
      92,  93,  94,  95,  136, 0,   0,   143, 147, 150, 0,   0,   236, 0,   96,  207, 0,   0,   0,
      208, 139, 142, 146, 149, 0,   0,   0,   0,   0,   158, 90,  91,  92,  93,  94,  95,  0,   197,
      168, 169, 170, 171, 172, 173, 174, 0,   96,  163, 164, 77,  166, 78,  79,  211, 80,  0,   214,
      0,   216, 85,  86,  87,  88,  219, 0,   0,   221, 0,   222, 0,   0,   186, 0,   224, 0,   226,
      0,   8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  243, 20,  0,   0,   244,
      0,   245, 108, 156, 110, 113, 0,   111, 0,   248, 0,   198, 0,   199, 0,   0,   0,   22,  23,
      257, 203, 258, 204, 24,  0,   25,  26,  0,   265, 90,  91,  92,  93,  94,  95,  202, 270, 0,
      90,  91,  92,  93,  94,  95,  0,   96,  213, 227, 215, 228, 0,   217, 0,   218, 96,  0,   220,
      85,  86,  87,  88,  0,   223, 0,   0,   90,  91,  92,  93,  94,  95,  0,   90,  91,  92,  93,
      94,  95,  0,   0,   242, 96,  0,   0,   0,   0,   0,   246, 96,  247, 8,   9,   10,  11,  12,
      13,  14,  15,  16,  17,  18,  0,   19,  0,   20,  0,   0,   21,  259, 0,   0,   0,   0,   264,
      0,   0,   266, 106, 107, 108, 156, 110, 0,   0,   111, 22,  23,  0,   0,   0,   0,   24,  0,
      25,  26,  8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  141, 20,  8,   9,
      10,  11,  12,  13,  14,  15,  16,  17,  18,  0,   19,  145, 20,  0,   0,   0,   0,   0,   22,
      23,  0,   0,   0,   0,   24,  0,   25,  26,  0,   0,   0,   0,   0,   22,  23,  0,   0,   0,
      0,   24,  263, 25,  26,  8,   9,   62,  11,  12,  13,  14,  63,  16,  17,  0,   0,   19,  0,
      20,  90,  91,  92,  93,  94,  95,  0,   0,   0,   0,   0,   0,   234, 0,   235, 0,   96,  0,
      0,   0,   22,  64,  255, 0,   256, 0,   24,  0,   25,  26,  90,  91,  92,  93,  94,  95,  0,
      0,   0,   0,   90,  91,  92,  93,  94,  95,  96,  230, 0,   0,   0,   0,   0,   238, 0,   0,
      96,  0,   0,   0,   0,   0,   0,   90,  91,  92,  93,  94,  95,  90,  91,  92,  93,  94,  95,
      249, 0,   0,   0,   96,  0,   254, 0,   0,   0,   96,  0,   0,   0,   0,   0,   90,  91,  92,
      93,  94,  95,  90,  91,  92,  93,  94,  95,  267, 0,   0,   0,   96,  0,   269, 0,   0,   0,
      96,  0,   0,   0,   0,   0,   90,  91,  92,  93,  94,  95,  90,  91,  92,  93,  94,  95,  231,
      0,   232, 0,   96,  0,   0,   233, 0,   0,   96,  0,   0,   0,   0,   0,   98,  99,  100, 101,
      102, 103, 104, 105, 106, 107, 108, 156, 110, 190, 0,   111, 0,   191, 0,   0,   0,   0,   39,
      40,  41,  42,  43,  44,  0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110,
      0,   195, 111, 196, 45,  46,  0,   0,   0,   0,   48,  49,  50,  51,  52,  53,  0,   98,  99,
      100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 0,   200, 111, 201, 54,  55,  0,   0,
      0,   0,   152, 32,  33,  34,  35,  36,  0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107,
      108, 109, 110, 0,   205, 111, 206, 37,  38,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 0,   251, 111, 252, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104,
      105, 106, 107, 108, 156, 110, 0,   260, 111, 261, 0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 253, 0,
      111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102,
      103, 104, 105, 106, 107, 108, 156, 110, 268, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 155,
      0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103,
      104, 105, 106, 107, 108, 156, 110, 229, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 237, 0,   111, 0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106,
      107, 108, 156, 110, 239, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,
      99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 240, 0,   111, 0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156,
      110, 250, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101,
      102, 103, 104, 105, 106, 107, 108, 156, 110, 262, 0,   111, 0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 271, 0,
      111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104,
      105, 106, 107, 108, 156, 110, 157, 0,   111, 0,   0,   0,   0,   0,   0,   0,   0,   0,   98,
      99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 225, 0,   111, 0,   0,   0,   0,
      0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 97,
      0,   111, 0,   0,   0,   0,   0,   0,   0,   98,  99,  100, 101, 102, 103, 104, 105, 106, 107,
      108, 109, 110, 209, 0,   111, 90,  91,  92,  93,  94,  95,  0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   96,  210, 98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 0,
      0,   111, 98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 0,   0,   111, 98,
      99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110, 0,   0,   111, 98,  99,  100, 101,
      102, 103, 104, 105, 0,   0,   0,   0,   110, 0,   0,   111, 99,  100, 101, 102, 103, 104, 105,
      106, 107, 108, 156, 110, 0,   0,   111, 100, 101, 102, 103, 104, 105, 106, 107, 108, 156, 110,
      0,   0,   111};

  const short Parser::yycheck_[] = {
      4,   19,  0,   1,   20,  38,  39,  40,  41,  22,  23,  24,  25,  26,  27,  19,  20,  20,  22,
      23,  24,  19,  38,  39,  40,  41,  22,  31,  32,  33,  34,  35,  36,  46,  47,  39,  40,  41,
      42,  43,  44,  40,  41,  47,  48,  49,  50,  51,  52,  53,  22,  49,  56,  57,  58,  59,  60,
      61,  20,  15,  64,  15,  17,  45,  82,  83,  16,  22,  15,  17,  32,  33,  34,  35,  36,  37,
      28,  81,  30,  31,  31,  85,  86,  87,  88,  4,   48,  48,  38,  39,  40,  41,  22,  22,  98,
      99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, -1,  16,  -1,  -1,  31,  14,
      3,   16,  5,   6,   23,  8,   39,  28,  -1,  30,  31,  -1,  31,  -1,  47,  48,  -1,  32,  33,
      34,  35,  36,  37,  56,  -1,  -1,  59,  60,  61,  -1,  -1,  14,  -1,  48,  152, -1,  -1,  -1,
      156, 58,  59,  60,  61,  -1,  -1,  -1,  -1,  -1,  81,  32,  33,  34,  35,  36,  37,  -1,  16,
      90,  91,  92,  93,  94,  95,  96,  -1,  48,  85,  86,  3,   88,  5,   6,   190, 8,   -1,  193,
      -1,  195, 38,  39,  40,  41,  200, -1,  -1,  203, -1,  205, -1,  -1,  109, -1,  210, -1,  212,
      -1,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,  227, 17,  -1,  -1,  231,
      -1,  233, 40,  41,  42,  152, -1,  45,  -1,  241, -1,  14,  -1,  16,  -1,  -1,  -1,  38,  39,
      251, 14,  253, 16,  44,  -1,  46,  47,  -1,  260, 32,  33,  34,  35,  36,  37,  16,  268, -1,
      32,  33,  34,  35,  36,  37,  -1,  48,  193, 14,  195, 16,  -1,  198, -1,  200, 48,  -1,  203,
      38,  39,  40,  41,  -1,  209, -1,  -1,  32,  33,  34,  35,  36,  37,  -1,  32,  33,  34,  35,
      36,  37,  -1,  -1,  227, 48,  -1,  -1,  -1,  -1,  -1,  234, 48,  236, 3,   4,   5,   6,   7,
      8,   9,   10,  11,  12,  13,  -1,  15,  -1,  17,  -1,  -1,  20,  255, -1,  -1,  -1,  -1,  260,
      -1,  -1,  263, 38,  39,  40,  41,  42,  -1,  -1,  45,  38,  39,  -1,  -1,  -1,  -1,  44,  -1,
      46,  47,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,  16,  17,  3,   4,
      5,   6,   7,   8,   9,   10,  11,  12,  13,  -1,  15,  16,  17,  -1,  -1,  -1,  -1,  -1,  38,
      39,  -1,  -1,  -1,  -1,  44,  -1,  46,  47,  -1,  -1,  -1,  -1,  -1,  38,  39,  -1,  -1,  -1,
      -1,  44,  14,  46,  47,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  -1,  -1,  15,  -1,
      17,  32,  33,  34,  35,  36,  37,  -1,  -1,  -1,  -1,  -1,  -1,  14,  -1,  16,  -1,  48,  -1,
      -1,  -1,  38,  39,  14,  -1,  16,  -1,  44,  -1,  46,  47,  32,  33,  34,  35,  36,  37,  -1,
      -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  48,  16,  -1,  -1,  -1,  -1,  -1,  16,  -1,  -1,
      48,  -1,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  32,  33,  34,  35,  36,  37,
      16,  -1,  -1,  -1,  48,  -1,  16,  -1,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,
      35,  36,  37,  32,  33,  34,  35,  36,  37,  16,  -1,  -1,  -1,  48,  -1,  16,  -1,  -1,  -1,
      48,  -1,  -1,  -1,  -1,  -1,  32,  33,  34,  35,  36,  37,  32,  33,  34,  35,  36,  37,  14,
      -1,  16,  -1,  48,  -1,  -1,  21,  -1,  -1,  48,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,
      34,  35,  36,  37,  38,  39,  40,  41,  42,  14,  -1,  45,  -1,  18,  -1,  -1,  -1,  -1,  22,
      23,  24,  25,  26,  27,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
      -1,  14,  45,  16,  46,  47,  -1,  -1,  -1,  -1,  22,  23,  24,  25,  26,  27,  -1,  30,  31,
      32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  46,  47,  -1,  -1,
      -1,  -1,  22,  23,  24,  25,  26,  27,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
      40,  41,  42,  -1,  14,  45,  16,  46,  47,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,
      37,  38,  39,  40,  41,  42,  -1,  14,  45,  16,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  14,  -1,
      45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,
      35,  36,  37,  38,  39,  40,  41,  42,  14,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,
      -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,
      36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,
      39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
      42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,
      34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  16,  -1,
      45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,
      37,  38,  39,  40,  41,  42,  18,  -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  18,  -1,  45,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  20,
      -1,  45,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
      40,  41,  42,  29,  -1,  45,  32,  33,  34,  35,  36,  37,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
      -1,  -1,  -1,  48,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,
      -1,  45,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  30,
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  -1,  -1,  45,  30,  31,  32,  33,
      34,  35,  36,  37,  -1,  -1,  -1,  -1,  42,  -1,  -1,  45,  31,  32,  33,  34,  35,  36,  37,
      38,  39,  40,  41,  42,  -1,  -1,  45,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
      -1,  -1,  45};

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
      16, 14, 16, 21, 14, 16, 14, 16, 16, 16, 16, 22, 55, 56, 56, 56, 55, 55, 56, 16, 16, 14, 16,
      14, 16, 14, 16, 56, 56, 55, 14, 16, 16, 14, 55, 56, 55, 16, 14, 16, 56, 16};

  const unsigned char Parser::yyr1_[] = {
      0,  50, 51, 51, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
      53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56};

  const unsigned char Parser::yyr2_[] = {
      0, 2, 0, 2, 1,  3,  3,  3,  2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
      3, 1, 4, 6, 6,  6,  4,  4,  3, 3, 3, 2, 3, 3, 3, 3, 3, 1, 1, 1, 3, 3, 3, 3, 3, 4, 3,
      4, 4, 3, 6, 12, 8,  8,  6,  5, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2,
      2, 2, 3, 3, 3,  3,  3,  3,  3, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 6, 6, 6,
      8, 6, 8, 8, 10, 10, 12, 14, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 1, 5, 4, 6, 6, 8};

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

#if YYDEBUG
  const unsigned short Parser::yyrline_[] = {
      0,   146, 146, 147, 150, 151, 158, 162, 163, 164, 167, 168, 169, 170, 171, 172, 173, 174, 175,
      176, 177, 178, 181, 182, 183, 184, 185, 186, 188, 189, 195, 201, 207, 213, 219, 225, 228, 230,
      238, 240, 248, 249, 250, 251, 260, 261, 262, 263, 265, 268, 272, 273, 274, 280, 286, 292, 298,
      299, 305, 311, 317, 323, 329, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 342, 345, 346,
      347, 348, 349, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 366, 368, 371, 374, 377,
      380, 382, 385, 388, 391, 394, 401, 408, 415, 422, 429, 436, 443, 450, 457, 463, 469, 475, 481,
      487, 493, 499, 500, 501, 502, 509, 516, 517, 518, 521, 522, 525, 526, 527, 528, 529, 551};

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
#endif // YYDEBUG

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
#line 2471 "apr_parser.cc" // lalr1.cc:1242
#line 574 "aprepro.yy"     // lalr1.cc:1243

void SEAMS::Parser::error(const std::string &m) { aprepro.error(m); }
