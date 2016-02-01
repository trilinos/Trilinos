/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 1 "aprepro.y" /* yacc.c:339  */

/* 
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include "my_aprepro.h"
#include "ap_array.h"
#include <stdlib.h>

void immutable_modify(symrec* var);
void undefined_warning(char* var);
void redefined_warning(symrec* var);
void set_type(symrec *var, int type);
void yyerror(char* var);
int  yylex(void);
extern int echo;
extern int state_immutable;
extern aprepro_options ap_options;

symrec *format;

#line 118 "y.tab.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NUM = 258,
    QSTRING = 259,
    UNDVAR = 260,
    VAR = 261,
    SVAR = 262,
    IMMVAR = 263,
    IMMSVAR = 264,
    AVAR = 265,
    FNCT = 266,
    SFNCT = 267,
    AFNCT = 268,
    EQ_PLUS = 269,
    EQ_MINUS = 270,
    EQ_TIME = 271,
    EQ_DIV = 272,
    EQ_POW = 273,
    LOR = 274,
    LAND = 275,
    LE = 276,
    GE = 277,
    EQ = 278,
    NE = 279,
    UNARY = 280,
    NOT = 281,
    POW = 282,
    INC = 283,
    DEC = 284,
    CONCAT = 285
  };
#endif
/* Tokens.  */
#define NUM 258
#define QSTRING 259
#define UNDVAR 260
#define VAR 261
#define SVAR 262
#define IMMVAR 263
#define IMMSVAR 264
#define AVAR 265
#define FNCT 266
#define SFNCT 267
#define AFNCT 268
#define EQ_PLUS 269
#define EQ_MINUS 270
#define EQ_TIME 271
#define EQ_DIV 272
#define EQ_POW 273
#define LOR 274
#define LAND 275
#define LE 276
#define GE 277
#define EQ 278
#define NE 279
#define UNARY 280
#define NOT 281
#define POW 282
#define INC 283
#define DEC 284
#define CONCAT 285

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 53 "aprepro.y" /* yacc.c:355  */

  double  val;          /* For returning numbers.               */
  symrec *tptr;         /* For returning symbol-table pointers  */
  char   *string;       /* For returning quoted strings         */
  array  *array;

#line 225 "y.tab.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 242 "y.tab.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1052

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  50
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  7
/* YYNRULES -- Number of rules.  */
#define YYNRULES  121
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  255

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   285

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      41,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,    34,     2,     2,
      44,    45,    33,    31,    46,    30,     2,    32,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    21,    47,
      24,    14,    25,    20,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    48,     2,    49,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    42,     2,    43,     2,     2,     2,     2,
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
       5,     6,     7,     8,     9,    10,    11,    12,    13,    15,
      16,    17,    18,    19,    22,    23,    26,    27,    28,    29,
      35,    36,    37,    38,    39,    40
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    94,    94,    95,    98,    99,   103,   105,   106,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   119,
     120,   123,   124,   125,   126,   127,   128,   130,   131,   132,
     134,   136,   138,   141,   143,   151,   153,   161,   162,   163,
     164,   173,   174,   175,   176,   178,   181,   185,   186,   187,
     188,   189,   190,   191,   196,   198,   200,   202,   204,   206,
     207,   208,   209,   210,   211,   212,   213,   214,   215,   217,
     220,   221,   222,   223,   224,   229,   230,   231,   232,   233,
     234,   235,   236,   237,   238,   239,   241,   243,   246,   249,
     252,   255,   257,   260,   263,   266,   269,   275,   276,   277,
     278,   279,   281,   283,   285,   287,   289,   291,   293,   294,
     295,   296,   303,   310,   311,   312,   315,   316,   321,   322,
     324,   325
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUM", "QSTRING", "UNDVAR", "VAR",
  "SVAR", "IMMVAR", "IMMSVAR", "AVAR", "FNCT", "SFNCT", "AFNCT", "'='",
  "EQ_PLUS", "EQ_MINUS", "EQ_TIME", "EQ_DIV", "EQ_POW", "'?'", "':'",
  "LOR", "LAND", "'<'", "'>'", "LE", "GE", "EQ", "NE", "'-'", "'+'", "'/'",
  "'*'", "'%'", "UNARY", "NOT", "POW", "INC", "DEC", "CONCAT", "'\\n'",
  "'{'", "'}'", "'('", "')'", "','", "';'", "'['", "']'", "$accept",
  "input", "line", "bool", "aexp", "sexp", "exp", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,    61,   269,   270,   271,   272,   273,
      63,    58,   274,   275,    60,    62,   276,   277,   278,   279,
      45,    43,    47,    42,    37,   280,   281,   282,   283,   284,
     285,    10,   123,   125,    40,    41,    44,    59,    91,    93
};
# endif

#define YYPACT_NINF -40

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-40)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     -40,     7,   -40,   -39,   -40,   297,   -40,   -40,   -40,   -40,
     255,   349,    -8,   376,    -1,    -9,   -27,   -26,   -22,   297,
     165,   165,   180,   202,   165,   165,    40,   197,   925,   849,
     297,   165,   165,   165,   165,   165,   -40,   -40,   165,   165,
     165,   165,   165,   165,   -40,   -40,   165,   165,   165,   165,
     165,   165,   165,   -40,   -40,   165,   297,   165,   119,   251,
     297,   -40,   210,   -16,   386,   -25,   165,   -16,   -16,   -40,
     -40,   -40,   -40,   -40,   -40,   -40,   -40,    44,   629,   289,
     165,   165,   165,   297,   297,   165,   297,   -40,   165,   165,
     165,   165,   165,   165,   165,   -40,   165,   165,   165,   165,
     165,   165,   165,   165,   165,   165,   165,   297,   165,   165,
     -40,    38,   210,   970,   986,   986,   986,   986,   986,   210,
     986,   986,   986,   986,   986,   986,   210,   986,   210,   986,
     986,   986,   986,   986,   986,   210,   986,    38,   970,   529,
     -40,   -21,   116,   404,   -40,   121,   763,   429,   149,   869,
     454,   165,   -16,   165,   -40,   -40,   -40,   934,   954,     8,
     986,   -40,   -30,   -30,   384,   -40,   384,     0,     0,     0,
       0,     0,     0,   -40,  1001,  1015,   246,   246,   246,   246,
     246,   246,   101,   101,   -16,   -40,   -16,   -16,   -16,   165,
     -40,   -40,   165,   -40,   165,   -40,   -40,   165,   -40,   165,
     -40,   -40,   -40,   165,   986,   -16,   165,   165,   324,   875,
     653,   352,   794,   554,   801,   677,   210,   986,    42,   -40,
     -40,   -40,   165,   165,   165,   165,   165,   -40,   165,   479,
     579,   897,   701,   771,   986,   -40,   165,   165,   -40,   -40,
     -40,   165,   504,   725,   824,   -40,   165,   -40,   165,   604,
     903,   165,   -40,   749,   -40
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,     4,     0,     3,     8,    59,    41,
      86,    62,    42,    63,    43,    27,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   118,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    89,    90,     0,     0,
       0,     0,     0,     0,    66,    67,     0,     0,     0,     0,
       0,     0,     0,    77,    78,     0,     0,     0,     0,     0,
       0,    35,     0,   113,    86,     0,     0,   114,    11,    60,
      87,    64,    75,    61,    88,    65,    76,   118,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     7,     0,     0,
       0,     0,     0,     0,     0,     6,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       5,    33,    44,    91,    92,    93,    94,    95,    96,    46,
      68,    70,    71,    72,    73,    74,    45,    69,    48,    79,
      81,    82,    83,    84,    85,    47,    80,    32,     0,     0,
      97,     0,     0,     0,    50,     0,     0,     0,     0,     0,
       0,     0,   113,     0,    20,   116,   117,     0,     0,    18,
       0,    19,    36,    34,    38,    40,    37,    21,    22,    23,
      24,    25,    26,    53,    16,    17,     9,    10,    12,    13,
      14,    15,   109,   108,   111,    39,   110,   112,   115,     0,
     100,    99,     0,    98,     0,    51,    49,     0,    52,     0,
      31,    28,    30,     0,    91,   110,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    58,   119,   120,   102,
     101,   103,     0,     0,     0,     0,     0,    29,     0,     0,
       0,     0,     0,     0,   121,   104,     0,     0,    56,    57,
      55,     0,     0,     0,     0,   106,     0,   105,     0,     0,
       0,     0,    54,     0,   107
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -40,   -40,   -40,   -23,    53,    27,    -5
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,     6,    26,    27,    62,   160
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      29,    77,    85,    86,     7,    56,    46,     2,     3,    83,
      84,    85,    86,    55,    63,    67,    68,    58,    59,    78,
      79,   109,    60,    57,   190,   113,   114,   115,   116,   117,
     118,    82,    28,   120,   121,   122,   123,   124,   125,    57,
      94,   127,   129,   130,   131,   132,   133,   134,     4,     5,
     136,   138,   139,   143,   147,   150,   228,   112,   159,   161,
      80,   152,    81,    82,    80,   119,    81,    82,    83,    84,
      85,    86,    61,   126,   128,   158,     0,     0,   138,   138,
     164,   166,   135,   111,     0,   142,   146,   149,     0,   154,
       0,   174,   175,   176,   177,   178,   179,   180,   181,   182,
     183,   184,   186,   187,   188,     0,     0,   157,     0,   137,
       0,   141,   145,   148,     0,   167,   168,   169,   170,   171,
     172,   173,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,   106,   153,   108,   162,   163,   109,   165,
      88,    89,    90,    91,    92,    93,   204,     0,   205,    19,
      20,    83,    84,    85,    86,    21,    94,    22,    23,     0,
     185,   191,   192,    24,   140,     0,   195,    25,     8,     9,
      64,    11,    12,    13,    14,    65,    16,    17,   112,    83,
      84,    85,    86,    69,   208,    70,    71,   210,    72,   211,
       0,     0,   213,     0,   200,    66,    20,     0,   215,     0,
       0,    21,   217,    22,    23,    73,     0,    74,    75,    24,
      76,     0,     0,    25,     0,     0,     0,   229,   230,   209,
     232,     0,     0,   234,   212,     0,   214,    83,    84,    85,
      86,   242,   243,   216,    88,    89,    90,    91,    92,    93,
      87,   249,     0,     0,     0,     0,   253,     0,     0,     0,
      94,   231,     0,   233,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,     0,     0,     0,   244,    30,
      31,    32,    33,    34,    35,   250,   104,   105,   106,   153,
     108,    19,    20,   109,     0,     0,     0,    21,     0,    22,
      23,     0,     0,    36,    37,    24,   144,     0,     0,    25,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   153,   108,     0,     0,   109,    19,    20,     0,
       0,     0,     0,    21,     0,    22,    23,     0,   156,     0,
       0,    24,     0,     0,     0,    25,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   153,   108,     0,
       0,   109,     0,    38,    39,    40,    41,    42,    43,     0,
       0,     0,     0,   218,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   153,   108,    44,    45,   109,
      47,    48,    49,    50,    51,    52,     0,   221,   222,   223,
     151,    31,    32,    33,    34,    35,    96,    97,    98,    99,
     100,   101,   102,   103,    53,    54,     0,     0,   108,     0,
       0,   109,     0,     0,    36,    37,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,     0,
       0,   109,     0,     0,     0,     0,     0,     0,     0,   193,
     194,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,     0,     0,   109,     0,     0,     0,
       0,     0,     0,     0,   198,   199,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,     0,
       0,   109,     0,     0,     0,     0,     0,     0,     0,   202,
     203,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   153,   108,     0,     0,   109,     0,     0,     0,
       0,     0,     0,     0,   235,   236,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   153,   108,     0,
       0,   109,     0,     0,     0,     0,     0,     0,     0,   245,
     246,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   153,   108,     0,     0,   109,     0,     0,     0,
       0,     0,     0,     0,     0,   189,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   153,   108,     0,
       0,   109,     0,     0,     0,     0,     0,     0,     0,     0,
     225,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   153,   108,     0,     0,   109,     0,     0,     0,
       0,     0,     0,     0,     0,   237,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   153,   108,     0,
       0,   109,     0,     0,     0,     0,     0,     0,     0,     0,
     251,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   153,   108,     0,     0,   109,     0,     0,     0,
       0,     0,     0,     0,   155,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,   106,   153,   108,     0,     0,
     109,     0,     0,     0,     0,     0,     0,     0,   220,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     153,   108,     0,     0,   109,     0,     0,     0,     0,     0,
       0,     0,   227,    96,    97,    98,    99,   100,   101,   102,
     103,   104,   105,   106,   153,   108,     0,     0,   109,     0,
       0,     0,     0,     0,     0,     0,   239,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   153,   108,
       0,     0,   109,     0,     0,     0,     0,     0,     0,     0,
     247,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   153,   108,     0,     0,   109,    88,    89,    90,
      91,    92,    93,     0,   254,    88,    89,    90,    91,    92,
      93,     0,     0,    94,     0,     0,     0,     0,   196,   197,
       0,    94,     0,     0,     0,     0,   240,   241,    88,    89,
      90,    91,    92,    93,     0,    88,    89,    90,    91,    92,
      93,     0,     0,     0,    94,     0,     0,     0,     0,     0,
     224,    94,     0,     0,     0,     0,     0,   226,    88,    89,
      90,    91,    92,    93,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    94,     0,     0,     0,     0,     0,
     248,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,     0,     0,   109,     0,     0,     0,
       0,     0,   110,    88,    89,    90,    91,    92,    93,    88,
      89,    90,    91,    92,    93,     0,     0,     0,     0,    94,
       0,     0,     0,     0,   201,    94,     0,     0,     0,     0,
     219,    88,    89,    90,    91,    92,    93,    88,    89,    90,
      91,    92,    93,     0,     0,     0,     0,    94,     0,     0,
       0,     0,   238,    94,     0,     0,     0,     0,   252,    88,
      89,    90,    91,    92,    93,   206,     0,     0,    88,    89,
      90,    91,    92,    93,     0,    94,     0,     0,    95,     0,
       0,     0,     0,     0,    94,   207,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   153,   108,     0,
       0,   109,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,     0,     0,   109,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,   106,   153,
     108,     0,     0,   109,    97,    98,    99,   100,   101,   102,
     103,   104,   105,   106,   153,   108,     0,     0,   109,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   153,   108,
       0,     0,   109
};

static const yytype_int16 yycheck[] =
{
       5,    24,    32,    33,    43,    14,    14,     0,     1,    30,
      31,    32,    33,    14,    19,    20,    21,    44,    44,    24,
      25,    37,    44,    48,    45,    30,    31,    32,    33,    34,
      35,    23,     5,    38,    39,    40,    41,    42,    43,    48,
      40,    46,    47,    48,    49,    50,    51,    52,    41,    42,
      55,    56,    57,    58,    59,    60,    14,    30,    81,    82,
      20,    66,    22,    23,    20,    38,    22,    23,    30,    31,
      32,    33,    19,    46,    47,    80,    -1,    -1,    83,    84,
      85,    86,    55,    30,    -1,    58,    59,    60,    -1,    45,
      -1,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,    -1,    -1,    80,    -1,    56,
      -1,    58,    59,    60,    -1,    88,    89,    90,    91,    92,
      93,    94,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,    32,    33,    34,    83,    84,    37,    86,
      24,    25,    26,    27,    28,    29,   151,    -1,   153,    30,
      31,    30,    31,    32,    33,    36,    40,    38,    39,    -1,
     107,    45,    46,    44,    45,    -1,    45,    48,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,   151,    30,
      31,    32,    33,     3,   189,     5,     6,   192,     8,   194,
      -1,    -1,   197,    -1,    45,    30,    31,    -1,   203,    -1,
      -1,    36,   207,    38,    39,     3,    -1,     5,     6,    44,
       8,    -1,    -1,    48,    -1,    -1,    -1,   222,   223,   192,
     225,    -1,    -1,   228,   197,    -1,   199,    30,    31,    32,
      33,   236,   237,   206,    24,    25,    26,    27,    28,    29,
      43,   246,    -1,    -1,    -1,    -1,   251,    -1,    -1,    -1,
      40,   224,    -1,   226,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    -1,    -1,    -1,   241,    14,
      15,    16,    17,    18,    19,   248,    30,    31,    32,    33,
      34,    30,    31,    37,    -1,    -1,    -1,    36,    -1,    38,
      39,    -1,    -1,    38,    39,    44,    45,    -1,    -1,    48,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    -1,    -1,    37,    30,    31,    -1,
      -1,    -1,    -1,    36,    -1,    38,    39,    -1,    49,    -1,
      -1,    44,    -1,    -1,    -1,    48,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    -1,
      -1,    37,    -1,    14,    15,    16,    17,    18,    19,    -1,
      -1,    -1,    -1,    49,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    38,    39,    37,
      14,    15,    16,    17,    18,    19,    -1,    45,    46,    47,
      14,    15,    16,    17,    18,    19,    22,    23,    24,    25,
      26,    27,    28,    29,    38,    39,    -1,    -1,    34,    -1,
      -1,    37,    -1,    -1,    38,    39,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    -1,
      -1,    37,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    45,
      46,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    -1,    -1,    37,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    45,    46,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    -1,
      -1,    37,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    45,
      46,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    -1,    -1,    37,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    45,    46,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    -1,
      -1,    37,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    45,
      46,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    -1,    -1,    37,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    46,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    -1,
      -1,    37,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      46,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    -1,    -1,    37,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    46,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    -1,
      -1,    37,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      46,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    -1,    -1,    37,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    45,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    -1,    -1,
      37,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    45,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    -1,    -1,    37,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    45,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    -1,    -1,    37,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    45,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      -1,    -1,    37,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      45,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    -1,    -1,    37,    24,    25,    26,
      27,    28,    29,    -1,    45,    24,    25,    26,    27,    28,
      29,    -1,    -1,    40,    -1,    -1,    -1,    -1,    45,    46,
      -1,    40,    -1,    -1,    -1,    -1,    45,    46,    24,    25,
      26,    27,    28,    29,    -1,    24,    25,    26,    27,    28,
      29,    -1,    -1,    -1,    40,    -1,    -1,    -1,    -1,    -1,
      46,    40,    -1,    -1,    -1,    -1,    -1,    46,    24,    25,
      26,    27,    28,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    40,    -1,    -1,    -1,    -1,    -1,
      46,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    -1,    -1,    37,    -1,    -1,    -1,
      -1,    -1,    43,    24,    25,    26,    27,    28,    29,    24,
      25,    26,    27,    28,    29,    -1,    -1,    -1,    -1,    40,
      -1,    -1,    -1,    -1,    45,    40,    -1,    -1,    -1,    -1,
      45,    24,    25,    26,    27,    28,    29,    24,    25,    26,
      27,    28,    29,    -1,    -1,    -1,    -1,    40,    -1,    -1,
      -1,    -1,    45,    40,    -1,    -1,    -1,    -1,    45,    24,
      25,    26,    27,    28,    29,    21,    -1,    -1,    24,    25,
      26,    27,    28,    29,    -1,    40,    -1,    -1,    43,    -1,
      -1,    -1,    -1,    -1,    40,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    -1,
      -1,    37,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    -1,    -1,    37,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    -1,    -1,    37,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    -1,    -1,    37,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      -1,    -1,    37
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    51,     0,     1,    41,    42,    52,    43,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    30,
      31,    36,    38,    39,    44,    48,    53,    54,    55,    56,
      14,    15,    16,    17,    18,    19,    38,    39,    14,    15,
      16,    17,    18,    19,    38,    39,    14,    14,    15,    16,
      17,    18,    19,    38,    39,    14,    14,    48,    44,    44,
      44,    54,    55,    56,     5,    10,    30,    56,    56,     3,
       5,     6,     8,     3,     5,     6,     8,    53,    56,    56,
      20,    22,    23,    30,    31,    32,    33,    43,    24,    25,
      26,    27,    28,    29,    40,    43,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    37,
      43,    54,    55,    56,    56,    56,    56,    56,    56,    55,
      56,    56,    56,    56,    56,    56,    55,    56,    55,    56,
      56,    56,    56,    56,    56,    55,    56,    54,    56,    56,
      45,    54,    55,    56,    45,    54,    55,    56,    54,    55,
      56,    14,    56,    33,    45,    45,    49,    55,    56,    53,
      56,    53,    54,    54,    56,    54,    56,    55,    55,    55,
      55,    55,    55,    55,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    54,    56,    56,    56,    46,
      45,    45,    46,    45,    46,    45,    45,    46,    45,    46,
      45,    45,    45,    46,    56,    56,    21,    21,    56,    55,
      56,    56,    55,    56,    55,    56,    55,    56,    49,    45,
      45,    45,    46,    47,    46,    46,    46,    45,    14,    56,
      56,    55,    56,    55,    56,    45,    46,    46,    45,    45,
      45,    46,    56,    56,    55,    45,    46,    45,    46,    56,
      55,    46,    45,    56,    45
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    50,    51,    51,    52,    52,    52,    52,    52,    53,
      53,    53,    53,    53,    53,    53,    53,    53,    53,    53,
      53,    53,    53,    53,    53,    53,    53,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    55,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    55,    55,    55,    55,    55,    55,    55,    55,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     3,     3,     3,     2,     3,
       3,     2,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     4,     6,
       4,     4,     3,     3,     3,     2,     3,     3,     3,     3,
       3,     1,     1,     1,     3,     3,     3,     3,     3,     4,
       3,     4,     4,     3,    12,     8,     8,     8,     5,     1,
       2,     2,     1,     1,     2,     2,     2,     2,     3,     3,
       3,     3,     3,     3,     3,     2,     2,     2,     2,     3,
       3,     3,     3,     3,     3,     3,     1,     2,     2,     2,
       2,     3,     3,     3,     3,     3,     3,     3,     4,     4,
       4,     6,     6,     6,     8,    10,    10,    14,     3,     3,
       3,     3,     3,     2,     2,     3,     3,     3,     1,     5,
       6,     8
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
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
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 4:
#line 98 "aprepro.y" /* yacc.c:1646  */
    { if (echo) fprintf(yyout,"\n");        }
#line 1640 "y.tab.c" /* yacc.c:1646  */
    break;

  case 5:
#line 99 "aprepro.y" /* yacc.c:1646  */
    { if (echo) {
                                     format = getsym("_FORMAT");
                                     fprintf(yyout, format->value.svar, (yyvsp[-1].val));
                                   }                                    }
#line 1649 "y.tab.c" /* yacc.c:1646  */
    break;

  case 6:
#line 103 "aprepro.y" /* yacc.c:1646  */
    { if (echo && (yyvsp[-1].string) != NULL)
                                    fprintf(yyout, "%s", (yyvsp[-1].string));   }
#line 1656 "y.tab.c" /* yacc.c:1646  */
    break;

  case 7:
#line 105 "aprepro.y" /* yacc.c:1646  */
    { ; }
#line 1662 "y.tab.c" /* yacc.c:1646  */
    break;

  case 8:
#line 106 "aprepro.y" /* yacc.c:1646  */
    { yyerrok;                              }
#line 1668 "y.tab.c" /* yacc.c:1646  */
    break;

  case 9:
#line 109 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) < (yyvsp[0].val);                         }
#line 1674 "y.tab.c" /* yacc.c:1646  */
    break;

  case 10:
#line 110 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) > (yyvsp[0].val);                         }
#line 1680 "y.tab.c" /* yacc.c:1646  */
    break;

  case 11:
#line 111 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = !((yyvsp[0].val));                           }
#line 1686 "y.tab.c" /* yacc.c:1646  */
    break;

  case 12:
#line 112 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) <= (yyvsp[0].val);                        }
#line 1692 "y.tab.c" /* yacc.c:1646  */
    break;

  case 13:
#line 113 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) >= (yyvsp[0].val);                        }
#line 1698 "y.tab.c" /* yacc.c:1646  */
    break;

  case 14:
#line 114 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) == (yyvsp[0].val);                        }
#line 1704 "y.tab.c" /* yacc.c:1646  */
    break;

  case 15:
#line 115 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) != (yyvsp[0].val);                        }
#line 1710 "y.tab.c" /* yacc.c:1646  */
    break;

  case 16:
#line 116 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) || (yyvsp[0].val);                        }
#line 1716 "y.tab.c" /* yacc.c:1646  */
    break;

  case 17:
#line 117 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) && (yyvsp[0].val);                        }
#line 1722 "y.tab.c" /* yacc.c:1646  */
    break;

  case 18:
#line 118 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) || (yyvsp[0].val);                        }
#line 1728 "y.tab.c" /* yacc.c:1646  */
    break;

  case 19:
#line 119 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) && (yyvsp[0].val);                        }
#line 1734 "y.tab.c" /* yacc.c:1646  */
    break;

  case 20:
#line 120 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-1].val);                              }
#line 1740 "y.tab.c" /* yacc.c:1646  */
    break;

  case 21:
#line 123 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (strcmp((yyvsp[-2].string),(yyvsp[0].string)) <  0 ? 1 : 0);    }
#line 1746 "y.tab.c" /* yacc.c:1646  */
    break;

  case 22:
#line 124 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (strcmp((yyvsp[-2].string),(yyvsp[0].string)) >  0 ? 1 : 0);    }
#line 1752 "y.tab.c" /* yacc.c:1646  */
    break;

  case 23:
#line 125 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (strcmp((yyvsp[-2].string),(yyvsp[0].string)) <= 0 ? 1 : 0);    }
#line 1758 "y.tab.c" /* yacc.c:1646  */
    break;

  case 24:
#line 126 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (strcmp((yyvsp[-2].string),(yyvsp[0].string)) >= 0 ? 1 : 0);    }
#line 1764 "y.tab.c" /* yacc.c:1646  */
    break;

  case 25:
#line 127 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (strcmp((yyvsp[-2].string),(yyvsp[0].string)) == 0 ? 1 : 0);    }
#line 1770 "y.tab.c" /* yacc.c:1646  */
    break;

  case 26:
#line 128 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (strcmp((yyvsp[-2].string),(yyvsp[0].string)) != 0 ? 1 : 0);    }
#line 1776 "y.tab.c" /* yacc.c:1646  */
    break;

  case 27:
#line 130 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = (yyvsp[0].tptr)->value.avar;}
#line 1782 "y.tab.c" /* yacc.c:1646  */
    break;

  case 28:
#line 131 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = (*((yyvsp[-3].tptr)->value.arrfnct))((yyvsp[-1].string));      }
#line 1788 "y.tab.c" /* yacc.c:1646  */
    break;

  case 29:
#line 133 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = (*((yyvsp[-5].tptr)->value.arrfnct))((yyvsp[-3].val),(yyvsp[-1].val));   }
#line 1794 "y.tab.c" /* yacc.c:1646  */
    break;

  case 30:
#line 135 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = (*((yyvsp[-3].tptr)->value.arrfnct))((yyvsp[-1].val));      }
#line 1800 "y.tab.c" /* yacc.c:1646  */
    break;

  case 31:
#line 137 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = (*((yyvsp[-3].tptr)->value.arrfnct))((yyvsp[-1].array));      }
#line 1806 "y.tab.c" /* yacc.c:1646  */
    break;

  case 32:
#line 138 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = (yyvsp[0].array); (yyvsp[-2].tptr)->value.avar = (yyvsp[0].array); 
                                  redefined_warning((yyvsp[-2].tptr));
                                  set_type((yyvsp[-2].tptr), AVAR); }
#line 1814 "y.tab.c" /* yacc.c:1646  */
    break;

  case 33:
#line 141 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = (yyvsp[0].array); (yyvsp[-2].tptr)->value.avar = (yyvsp[0].array); 
                                  set_type((yyvsp[-2].tptr), AVAR); }
#line 1821 "y.tab.c" /* yacc.c:1646  */
    break;

  case 34:
#line 143 "aprepro.y" /* yacc.c:1646  */
    { if ((yyvsp[-2].array)->cols == (yyvsp[0].array)->cols && (yyvsp[-2].array)->rows == (yyvsp[0].array)->rows ) {
                                     (yyval.array) = array_add((yyvsp[-2].array), (yyvsp[0].array)); 
                                  }
                                  else {
                                    yyerror("Arrays do not have same row and column count"); 
                                    yyerrok;
                                  }
                                }
#line 1834 "y.tab.c" /* yacc.c:1646  */
    break;

  case 35:
#line 151 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = array_scale((yyvsp[0].array), -1.0);           }
#line 1840 "y.tab.c" /* yacc.c:1646  */
    break;

  case 36:
#line 153 "aprepro.y" /* yacc.c:1646  */
    { if ((yyvsp[-2].array)->cols == (yyvsp[0].array)->cols && (yyvsp[-2].array)->rows == (yyvsp[0].array)->rows ) {
                                     (yyval.array) = array_sub((yyvsp[-2].array), (yyvsp[0].array)); 
                                  }
                                  else {
                                    yyerror("Arrays do not have same row and column count"); 
                                    yyerrok;
                                  }
                                }
#line 1853 "y.tab.c" /* yacc.c:1646  */
    break;

  case 37:
#line 161 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = array_scale((yyvsp[-2].array), (yyvsp[0].val));             }
#line 1859 "y.tab.c" /* yacc.c:1646  */
    break;

  case 38:
#line 162 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = array_scale((yyvsp[-2].array), 1.0/(yyvsp[0].val));         }
#line 1865 "y.tab.c" /* yacc.c:1646  */
    break;

  case 39:
#line 163 "aprepro.y" /* yacc.c:1646  */
    { (yyval.array) = array_scale((yyvsp[0].array), (yyvsp[-2].val));             }
#line 1871 "y.tab.c" /* yacc.c:1646  */
    break;

  case 40:
#line 164 "aprepro.y" /* yacc.c:1646  */
    { if ((yyvsp[-2].array)->cols == (yyvsp[0].array)->rows) {
                                    (yyval.array) = array_mult((yyvsp[-2].array), (yyvsp[0].array));
                                  }
                                  else {
                                    yyerror("Column count of first array does not match row count of second array"); 
                                    yyerrok;
                                  }
				}
#line 1884 "y.tab.c" /* yacc.c:1646  */
    break;

  case 41:
#line 173 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (yyvsp[0].string);                              }
#line 1890 "y.tab.c" /* yacc.c:1646  */
    break;

  case 42:
#line 174 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (yyvsp[0].tptr)->value.svar;                  }
#line 1896 "y.tab.c" /* yacc.c:1646  */
    break;

  case 43:
#line 175 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (yyvsp[0].tptr)->value.svar;                  }
#line 1902 "y.tab.c" /* yacc.c:1646  */
    break;

  case 44:
#line 176 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (yyvsp[0].string); (yyvsp[-2].tptr)->value.svar = (yyvsp[0].string);
                                  set_type((yyvsp[-2].tptr), SVAR);                   }
#line 1909 "y.tab.c" /* yacc.c:1646  */
    break;

  case 45:
#line 178 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (yyvsp[0].string); 
                                  (yyvsp[-2].tptr)->value.svar = (yyvsp[0].string);
                                  redefined_warning((yyvsp[-2].tptr));          }
#line 1917 "y.tab.c" /* yacc.c:1646  */
    break;

  case 46:
#line 181 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (yyvsp[0].string); 
                                  (yyvsp[-2].tptr)->value.svar= (yyvsp[0].string);
                                  redefined_warning((yyvsp[-2].tptr));          
                                  (yyvsp[-2].tptr)->type = SVAR;              }
#line 1926 "y.tab.c" /* yacc.c:1646  */
    break;

  case 47:
#line 185 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-2].tptr)); YYERROR; }
#line 1932 "y.tab.c" /* yacc.c:1646  */
    break;

  case 48:
#line 186 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-2].tptr)); YYERROR; }
#line 1938 "y.tab.c" /* yacc.c:1646  */
    break;

  case 49:
#line 187 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (*((yyvsp[-3].tptr)->value.strfnct))((yyvsp[-1].string));      }
#line 1944 "y.tab.c" /* yacc.c:1646  */
    break;

  case 50:
#line 188 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (*((yyvsp[-2].tptr)->value.strfnct))();        }
#line 1950 "y.tab.c" /* yacc.c:1646  */
    break;

  case 51:
#line 189 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (*((yyvsp[-3].tptr)->value.strfnct))((yyvsp[-1].array));      }
#line 1956 "y.tab.c" /* yacc.c:1646  */
    break;

  case 52:
#line 190 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (*((yyvsp[-3].tptr)->value.strfnct))((yyvsp[-1].val));      }
#line 1962 "y.tab.c" /* yacc.c:1646  */
    break;

  case 53:
#line 191 "aprepro.y" /* yacc.c:1646  */
    { int len1 = strlen((yyvsp[-2].string));
                                  int len3 = strlen((yyvsp[0].string));
                                  ALLOC((yyval.string), len1+len3+1, char *);
                                  memcpy((yyval.string), (yyvsp[-2].string), len1+1);
                                  strcat((yyval.string), (yyvsp[0].string)); }
#line 1972 "y.tab.c" /* yacc.c:1646  */
    break;

  case 54:
#line 197 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (*((yyvsp[-11].tptr)->value.strfnct))((yyvsp[-9].val), (yyvsp[-7].string), (yyvsp[-5].string), (yyvsp[-3].string), (yyvsp[-1].string)); }
#line 1978 "y.tab.c" /* yacc.c:1646  */
    break;

  case 55:
#line 199 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (*((yyvsp[-7].tptr)->value.strfnct))((yyvsp[-5].val), (yyvsp[-3].string), (yyvsp[-1].string)); }
#line 1984 "y.tab.c" /* yacc.c:1646  */
    break;

  case 56:
#line 201 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (*((yyvsp[-7].tptr)->value.strfnct))((yyvsp[-5].string), (yyvsp[-3].string), (yyvsp[-1].string)); }
#line 1990 "y.tab.c" /* yacc.c:1646  */
    break;

  case 57:
#line 203 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = (*((yyvsp[-7].tptr)->value.strfnct))((yyvsp[-5].string), (yyvsp[-3].val), (yyvsp[-1].val)); }
#line 1996 "y.tab.c" /* yacc.c:1646  */
    break;

  case 58:
#line 204 "aprepro.y" /* yacc.c:1646  */
    { (yyval.string) = ((yyvsp[-4].val)) ? ((yyvsp[-2].string)) : ((yyvsp[0].string));              }
#line 2002 "y.tab.c" /* yacc.c:1646  */
    break;

  case 59:
#line 206 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].val);                              }
#line 2008 "y.tab.c" /* yacc.c:1646  */
    break;

  case 60:
#line 207 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].val) + 1;                          }
#line 2014 "y.tab.c" /* yacc.c:1646  */
    break;

  case 61:
#line 208 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].val) - 1;                          }
#line 2020 "y.tab.c" /* yacc.c:1646  */
    break;

  case 62:
#line 209 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].tptr)->value.var;                   }
#line 2026 "y.tab.c" /* yacc.c:1646  */
    break;

  case 63:
#line 210 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].tptr)->value.var;                   }
#line 2032 "y.tab.c" /* yacc.c:1646  */
    break;

  case 64:
#line 211 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = ++((yyvsp[0].tptr)->value.var);               }
#line 2038 "y.tab.c" /* yacc.c:1646  */
    break;

  case 65:
#line 212 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = --((yyvsp[0].tptr)->value.var);               }
#line 2044 "y.tab.c" /* yacc.c:1646  */
    break;

  case 66:
#line 213 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = ((yyvsp[-1].tptr)->value.var)++;               }
#line 2050 "y.tab.c" /* yacc.c:1646  */
    break;

  case 67:
#line 214 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = ((yyvsp[-1].tptr)->value.var)--;               }
#line 2056 "y.tab.c" /* yacc.c:1646  */
    break;

  case 68:
#line 215 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].val); (yyvsp[-2].tptr)->value.var = (yyvsp[0].val);
                                  redefined_warning((yyvsp[-2].tptr));          }
#line 2063 "y.tab.c" /* yacc.c:1646  */
    break;

  case 69:
#line 217 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].val); (yyvsp[-2].tptr)->value.var = (yyvsp[0].val);
                                  redefined_warning((yyvsp[-2].tptr));          
                                  (yyvsp[-2].tptr)->type = VAR;                       }
#line 2071 "y.tab.c" /* yacc.c:1646  */
    break;

  case 70:
#line 220 "aprepro.y" /* yacc.c:1646  */
    { (yyvsp[-2].tptr)->value.var += (yyvsp[0].val); (yyval.val) = (yyvsp[-2].tptr)->value.var; }
#line 2077 "y.tab.c" /* yacc.c:1646  */
    break;

  case 71:
#line 221 "aprepro.y" /* yacc.c:1646  */
    { (yyvsp[-2].tptr)->value.var -= (yyvsp[0].val); (yyval.val) = (yyvsp[-2].tptr)->value.var; }
#line 2083 "y.tab.c" /* yacc.c:1646  */
    break;

  case 72:
#line 222 "aprepro.y" /* yacc.c:1646  */
    { (yyvsp[-2].tptr)->value.var *= (yyvsp[0].val); (yyval.val) = (yyvsp[-2].tptr)->value.var; }
#line 2089 "y.tab.c" /* yacc.c:1646  */
    break;

  case 73:
#line 223 "aprepro.y" /* yacc.c:1646  */
    { (yyvsp[-2].tptr)->value.var /= (yyvsp[0].val); (yyval.val) = (yyvsp[-2].tptr)->value.var; }
#line 2095 "y.tab.c" /* yacc.c:1646  */
    break;

  case 74:
#line 224 "aprepro.y" /* yacc.c:1646  */
    { errno = 0;
                                  (yyvsp[-2].tptr)->value.var = pow((yyvsp[-2].tptr)->value.var,(yyvsp[0].val)); 
                                  (yyval.val) = (yyvsp[-2].tptr)->value.var; 
                                  MATH_ERROR("Power");
                                }
#line 2105 "y.tab.c" /* yacc.c:1646  */
    break;

  case 75:
#line 229 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[0].tptr)); YYERROR; }
#line 2111 "y.tab.c" /* yacc.c:1646  */
    break;

  case 76:
#line 230 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[0].tptr)); YYERROR; }
#line 2117 "y.tab.c" /* yacc.c:1646  */
    break;

  case 77:
#line 231 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-1].tptr)); YYERROR; }
#line 2123 "y.tab.c" /* yacc.c:1646  */
    break;

  case 78:
#line 232 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-1].tptr)); YYERROR; }
#line 2129 "y.tab.c" /* yacc.c:1646  */
    break;

  case 79:
#line 233 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-2].tptr)); YYERROR; }
#line 2135 "y.tab.c" /* yacc.c:1646  */
    break;

  case 80:
#line 234 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-2].tptr)); YYERROR; }
#line 2141 "y.tab.c" /* yacc.c:1646  */
    break;

  case 81:
#line 235 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-2].tptr)); YYERROR; }
#line 2147 "y.tab.c" /* yacc.c:1646  */
    break;

  case 82:
#line 236 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-2].tptr)); YYERROR; }
#line 2153 "y.tab.c" /* yacc.c:1646  */
    break;

  case 83:
#line 237 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-2].tptr)); YYERROR; }
#line 2159 "y.tab.c" /* yacc.c:1646  */
    break;

  case 84:
#line 238 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-2].tptr)); YYERROR; }
#line 2165 "y.tab.c" /* yacc.c:1646  */
    break;

  case 85:
#line 239 "aprepro.y" /* yacc.c:1646  */
    { immutable_modify((yyvsp[-2].tptr)); YYERROR; }
#line 2171 "y.tab.c" /* yacc.c:1646  */
    break;

  case 86:
#line 241 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].tptr)->value.var;
                                  undefined_warning((yyvsp[0].tptr)->name);          }
#line 2178 "y.tab.c" /* yacc.c:1646  */
    break;

  case 87:
#line 243 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = ++((yyvsp[0].tptr)->value.var);               
                                  set_type((yyvsp[0].tptr), VAR);                       
                                  undefined_warning((yyvsp[0].tptr)->name);          }
#line 2186 "y.tab.c" /* yacc.c:1646  */
    break;

  case 88:
#line 246 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = --((yyvsp[0].tptr)->value.var);               
                                  set_type((yyvsp[0].tptr), VAR);                       
                                  undefined_warning((yyvsp[0].tptr)->name);          }
#line 2194 "y.tab.c" /* yacc.c:1646  */
    break;

  case 89:
#line 249 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = ((yyvsp[-1].tptr)->value.var)++;               
                                  set_type((yyvsp[-1].tptr), VAR);                       
                                  undefined_warning((yyvsp[-1].tptr)->name);          }
#line 2202 "y.tab.c" /* yacc.c:1646  */
    break;

  case 90:
#line 252 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = ((yyvsp[-1].tptr)->value.var)--;               
                                  set_type((yyvsp[-1].tptr), VAR);                       
                                  undefined_warning((yyvsp[-1].tptr)->name);          }
#line 2210 "y.tab.c" /* yacc.c:1646  */
    break;

  case 91:
#line 255 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].val); (yyvsp[-2].tptr)->value.var = (yyvsp[0].val);
                                  set_type((yyvsp[-2].tptr), VAR);                    }
#line 2217 "y.tab.c" /* yacc.c:1646  */
    break;

  case 92:
#line 257 "aprepro.y" /* yacc.c:1646  */
    { (yyvsp[-2].tptr)->value.var += (yyvsp[0].val); (yyval.val) = (yyvsp[-2].tptr)->value.var; 
                                  set_type((yyvsp[-2].tptr), VAR);                       
                                  undefined_warning((yyvsp[-2].tptr)->name);          }
#line 2225 "y.tab.c" /* yacc.c:1646  */
    break;

  case 93:
#line 260 "aprepro.y" /* yacc.c:1646  */
    { (yyvsp[-2].tptr)->value.var -= (yyvsp[0].val); (yyval.val) = (yyvsp[-2].tptr)->value.var; 
                                  set_type((yyvsp[-2].tptr), VAR);                       
                                  undefined_warning((yyvsp[-2].tptr)->name);          }
#line 2233 "y.tab.c" /* yacc.c:1646  */
    break;

  case 94:
#line 263 "aprepro.y" /* yacc.c:1646  */
    { (yyvsp[-2].tptr)->value.var *= (yyvsp[0].val); (yyval.val) = (yyvsp[-2].tptr)->value.var; 
                                  set_type((yyvsp[-2].tptr), VAR);                       
                                  undefined_warning((yyvsp[-2].tptr)->name);          }
#line 2241 "y.tab.c" /* yacc.c:1646  */
    break;

  case 95:
#line 266 "aprepro.y" /* yacc.c:1646  */
    { (yyvsp[-2].tptr)->value.var /= (yyvsp[0].val); (yyval.val) = (yyvsp[-2].tptr)->value.var; 
                                  set_type((yyvsp[-2].tptr), VAR);                       
                                  undefined_warning((yyvsp[-2].tptr)->name);          }
#line 2249 "y.tab.c" /* yacc.c:1646  */
    break;

  case 96:
#line 269 "aprepro.y" /* yacc.c:1646  */
    { errno = 0;
                                  (yyvsp[-2].tptr)->value.var = pow((yyvsp[-2].tptr)->value.var,(yyvsp[0].val)); 
                                  (yyval.val) = (yyvsp[-2].tptr)->value.var; 
                                  set_type((yyvsp[-2].tptr), VAR);                       
                                  MATH_ERROR("Power");
                                  undefined_warning((yyvsp[-2].tptr)->name);          }
#line 2260 "y.tab.c" /* yacc.c:1646  */
    break;

  case 97:
#line 275 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-2].tptr)->value.fnctptr))();        }
#line 2266 "y.tab.c" /* yacc.c:1646  */
    break;

  case 98:
#line 276 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-3].tptr)->value.fnctptr))((yyvsp[-1].val));      }
#line 2272 "y.tab.c" /* yacc.c:1646  */
    break;

  case 99:
#line 277 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-3].tptr)->value.fnctptr))((yyvsp[-1].string));      }
#line 2278 "y.tab.c" /* yacc.c:1646  */
    break;

  case 100:
#line 278 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-3].tptr)->value.fnctptr))((yyvsp[-1].array));      }
#line 2284 "y.tab.c" /* yacc.c:1646  */
    break;

  case 101:
#line 280 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-5].tptr)->value.fnctptr))((yyvsp[-3].string), (yyvsp[-1].val));  }
#line 2290 "y.tab.c" /* yacc.c:1646  */
    break;

  case 102:
#line 282 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-5].tptr)->value.fnctptr))((yyvsp[-3].string), (yyvsp[-1].string));  }
#line 2296 "y.tab.c" /* yacc.c:1646  */
    break;

  case 103:
#line 284 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-5].tptr)->value.fnctptr))((yyvsp[-3].val), (yyvsp[-1].val));  }
#line 2302 "y.tab.c" /* yacc.c:1646  */
    break;

  case 104:
#line 286 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-7].tptr)->value.fnctptr))((yyvsp[-5].val), (yyvsp[-3].val), (yyvsp[-1].val)); }
#line 2308 "y.tab.c" /* yacc.c:1646  */
    break;

  case 105:
#line 288 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-9].tptr)->value.fnctptr))((yyvsp[-7].val), (yyvsp[-5].val), (yyvsp[-3].val), (yyvsp[-1].val)); }
#line 2314 "y.tab.c" /* yacc.c:1646  */
    break;

  case 106:
#line 290 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-9].tptr)->value.fnctptr))((yyvsp[-7].val), (yyvsp[-5].val), (yyvsp[-3].val), (yyvsp[-1].val)); }
#line 2320 "y.tab.c" /* yacc.c:1646  */
    break;

  case 107:
#line 292 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (*((yyvsp[-13].tptr)->value.fnctptr))((yyvsp[-11].val), (yyvsp[-9].val), (yyvsp[-7].val), (yyvsp[-5].val), (yyvsp[-3].val), (yyvsp[-1].val)); }
#line 2326 "y.tab.c" /* yacc.c:1646  */
    break;

  case 108:
#line 293 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) + (yyvsp[0].val);                         }
#line 2332 "y.tab.c" /* yacc.c:1646  */
    break;

  case 109:
#line 294 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) - (yyvsp[0].val);                         }
#line 2338 "y.tab.c" /* yacc.c:1646  */
    break;

  case 110:
#line 295 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) * (yyvsp[0].val);                         }
#line 2344 "y.tab.c" /* yacc.c:1646  */
    break;

  case 111:
#line 296 "aprepro.y" /* yacc.c:1646  */
    { if ((yyvsp[0].val) == 0.)
                                    {
                                      yyerror("Zero divisor"); 
                                      yyerrok;
                                    }
                                  else
                                    (yyval.val) = (yyvsp[-2].val) / (yyvsp[0].val);                       }
#line 2356 "y.tab.c" /* yacc.c:1646  */
    break;

  case 112:
#line 303 "aprepro.y" /* yacc.c:1646  */
    { if ((yyvsp[0].val) == 0.)
                                    {
                                      yyerror("Zero divisor");
                                      yyerrok;
                                    }
                                  else
                                    (yyval.val) = (int)(yyvsp[-2].val) % (int)(yyvsp[0].val);             }
#line 2368 "y.tab.c" /* yacc.c:1646  */
    break;

  case 113:
#line 310 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = -(yyvsp[0].val);                             }
#line 2374 "y.tab.c" /* yacc.c:1646  */
    break;

  case 114:
#line 311 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) =  (yyvsp[0].val);                             }
#line 2380 "y.tab.c" /* yacc.c:1646  */
    break;

  case 115:
#line 312 "aprepro.y" /* yacc.c:1646  */
    { errno = 0;
                                  (yyval.val) = pow((yyvsp[-2].val), (yyvsp[0].val)); 
                                  MATH_ERROR("Power");                  }
#line 2388 "y.tab.c" /* yacc.c:1646  */
    break;

  case 116:
#line 315 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-1].val);                              }
#line 2394 "y.tab.c" /* yacc.c:1646  */
    break;

  case 117:
#line 316 "aprepro.y" /* yacc.c:1646  */
    { errno = 0;
                                  (yyval.val) = (double)((yyvsp[-1].val) < 0 ? 
                                        -floor(-((yyvsp[-1].val))): floor((yyvsp[-1].val)) );
                                  MATH_ERROR("floor (int)");            }
#line 2403 "y.tab.c" /* yacc.c:1646  */
    break;

  case 118:
#line 321 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = ((yyvsp[0].val)) ? 1 : 0; }
#line 2409 "y.tab.c" /* yacc.c:1646  */
    break;

  case 119:
#line 322 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = ((yyvsp[-4].val)) ? ((yyvsp[-2].val)) : ((yyvsp[0].val));              }
#line 2415 "y.tab.c" /* yacc.c:1646  */
    break;

  case 120:
#line 324 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = array_value((yyvsp[-5].tptr)->value.avar, (yyvsp[-3].val), (yyvsp[-1].val)); }
#line 2421 "y.tab.c" /* yacc.c:1646  */
    break;

  case 121:
#line 326 "aprepro.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].val);
				    array *arr = (yyvsp[-7].tptr)->value.avar;
                                    int cols = arr->cols;
                                    int rows = arr->rows;
				    int row = (yyvsp[-5].val);
				    int col = (yyvsp[-3].val);
				    if (ap_options.one_based_index == True) {
				      row--;
				      col--;
				    }
				    if (row < rows && col < cols) {
                                      int offset = row*cols+col;
                                      (yyvsp[-7].tptr)->value.avar->data[offset] = (yyvsp[0].val);
                                    }
                                    else {
                                      yyerror("Row or Column index out of range"); 
                                      yyerrok;
                                    }
                                  }
#line 2445 "y.tab.c" /* yacc.c:1646  */
    break;


#line 2449 "y.tab.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
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
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
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

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 348 "aprepro.y" /* yacc.c:1906  */

# include "lex.yy.c"


