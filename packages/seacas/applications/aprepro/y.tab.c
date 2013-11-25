/* A Bison parser, made by GNU Bison 2.7.12-4996.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.
   
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
#define YYBISON_VERSION "2.7.12-4996"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
/* Line 371 of yacc.c  */
#line 1 "aprepro.y"

/* 
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
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

/* Line 371 of yacc.c  */
#line 119 "y.tab.c"

# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
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
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUM = 258,
     QSTRING = 259,
     UNDVAR = 260,
     VAR = 261,
     SVAR = 262,
     IMMVAR = 263,
     IMMSVAR = 264,
     FNCT = 265,
     SFNCT = 266,
     EQ_MINUS = 267,
     EQ_PLUS = 268,
     EQ_DIV = 269,
     EQ_TIME = 270,
     EQ_POW = 271,
     LOR = 272,
     LAND = 273,
     NE = 274,
     EQ = 275,
     GE = 276,
     LE = 277,
     NOT = 278,
     UNARY = 279,
     POW = 280,
     DEC = 281,
     INC = 282,
     CONCAT = 283
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
#define FNCT 265
#define SFNCT 266
#define EQ_MINUS 267
#define EQ_PLUS 268
#define EQ_DIV 269
#define EQ_TIME 270
#define EQ_POW 271
#define LOR 272
#define LAND 273
#define NE 274
#define EQ 275
#define GE 276
#define LE 277
#define NOT 278
#define UNARY 279
#define POW 280
#define DEC 281
#define INC 282
#define CONCAT 283



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 387 of yacc.c  */
#line 52 "aprepro.y"

  double  val;		/* For returning numbers.		*/
  symrec *tptr;		/* For returning symbol-table pointers	*/
  char   *string;	/* For returning quoted strings		*/


/* Line 387 of yacc.c  */
#line 225 "y.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

/* Line 390 of yacc.c  */
#line 253 "y.tab.c"

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
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
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
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
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

#ifndef __attribute__
/* This feature is available in gcc versions 2.5 and later.  */
# if (! defined __GNUC__ || __GNUC__ < 2 \
      || (__GNUC__ == 2 && __GNUC_MINOR__ < 5))
#  define __attribute__(Spec) /* empty */
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif


/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(N) (N)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
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
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
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
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
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
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
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
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

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
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   805

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  48
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  6
/* YYNRULES -- Number of rules.  */
#define YYNRULES  102
/* YYNRULES -- Number of states.  */
#define YYNSTATES  207

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   283

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      39,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,    32,     2,     2,
      42,    43,    31,    29,    44,    28,     2,    30,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    19,    45,
      22,    12,    23,    18,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    46,     2,    47,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    40,     2,    41,     2,     2,     2,     2,
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
       5,     6,     7,     8,     9,    10,    11,    13,    14,    15,
      16,    17,    20,    21,    24,    25,    26,    27,    33,    34,
      35,    36,    37,    38
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    13,    17,    20,    24,
      28,    31,    35,    39,    43,    47,    51,    55,    59,    63,
      67,    71,    75,    79,    83,    87,    91,    93,    95,    97,
     101,   105,   109,   113,   117,   122,   126,   131,   135,   148,
     157,   166,   175,   181,   183,   186,   189,   191,   193,   196,
     199,   202,   205,   209,   213,   217,   221,   225,   229,   233,
     236,   239,   242,   245,   249,   253,   257,   261,   265,   269,
     273,   275,   278,   281,   284,   287,   291,   295,   299,   303,
     307,   311,   315,   320,   325,   332,   339,   346,   355,   366,
     377,   392,   396,   400,   404,   408,   412,   415,   418,   422,
     426,   430,   432
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      49,     0,    -1,    -1,    49,    50,    -1,    39,    -1,    40,
      53,    41,    -1,    40,    52,    41,    -1,     1,    41,    -1,
      53,    22,    53,    -1,    53,    23,    53,    -1,    33,    53,
      -1,    53,    27,    53,    -1,    53,    26,    53,    -1,    53,
      25,    53,    -1,    53,    24,    53,    -1,    53,    20,    53,
      -1,    53,    21,    53,    -1,    51,    20,    51,    -1,    51,
      21,    51,    -1,    42,    51,    43,    -1,    52,    22,    52,
      -1,    52,    23,    52,    -1,    52,    27,    52,    -1,    52,
      26,    52,    -1,    52,    25,    52,    -1,    52,    24,    52,
      -1,     4,    -1,     7,    -1,     9,    -1,     5,    12,    52,
      -1,     7,    12,    52,    -1,     6,    12,    52,    -1,     9,
      12,    52,    -1,     8,    12,    52,    -1,    11,    42,    52,
      43,    -1,    11,    42,    43,    -1,    11,    42,    53,    43,
      -1,    52,    38,    52,    -1,    11,    42,    53,    44,    52,
      44,    52,    44,    52,    44,    52,    43,    -1,    11,    42,
      53,    44,    52,    44,    52,    43,    -1,    11,    42,    52,
      44,    52,    44,    52,    43,    -1,    11,    42,    52,    44,
      53,    44,    53,    43,    -1,    51,    18,    52,    19,    52,
      -1,     3,    -1,    37,     3,    -1,    36,     3,    -1,     6,
      -1,     8,    -1,    37,     6,    -1,    36,     6,    -1,     6,
      37,    -1,     6,    36,    -1,     6,    12,    53,    -1,     7,
      12,    53,    -1,     6,    14,    53,    -1,     6,    13,    53,
      -1,     6,    16,    53,    -1,     6,    15,    53,    -1,     6,
      17,    53,    -1,    37,     8,    -1,    36,     8,    -1,     8,
      37,    -1,     8,    36,    -1,     8,    12,    53,    -1,     9,
      12,    53,    -1,     8,    14,    53,    -1,     8,    13,    53,
      -1,     8,    16,    53,    -1,     8,    15,    53,    -1,     8,
      17,    53,    -1,     5,    -1,    37,     5,    -1,    36,     5,
      -1,     5,    37,    -1,     5,    36,    -1,     5,    12,    53,
      -1,     5,    14,    53,    -1,     5,    13,    53,    -1,     5,
      16,    53,    -1,     5,    15,    53,    -1,     5,    17,    53,
      -1,    10,    42,    43,    -1,    10,    42,    53,    43,    -1,
      10,    42,    52,    43,    -1,    10,    42,    52,    44,    53,
      43,    -1,    10,    42,    52,    44,    52,    43,    -1,    10,
      42,    53,    44,    53,    43,    -1,    10,    42,    53,    44,
      53,    44,    53,    43,    -1,    10,    42,    53,    44,    53,
      45,    53,    44,    53,    43,    -1,    10,    42,    53,    44,
      53,    44,    53,    44,    53,    43,    -1,    10,    42,    53,
      44,    53,    44,    53,    44,    53,    44,    53,    44,    53,
      43,    -1,    53,    29,    53,    -1,    53,    28,    53,    -1,
      53,    31,    53,    -1,    53,    30,    53,    -1,    53,    32,
      53,    -1,    28,    53,    -1,    29,    53,    -1,    53,    35,
      53,    -1,    42,    53,    43,    -1,    46,    53,    47,    -1,
      51,    -1,    51,    18,    53,    19,    53,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    89,    89,    90,    93,    94,    98,   100,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     117,   118,   119,   120,   121,   122,   124,   125,   126,   127,
     129,   132,   136,   137,   138,   139,   140,   141,   146,   148,
     150,   152,   154,   156,   157,   158,   159,   160,   161,   162,
     163,   164,   165,   167,   170,   171,   172,   173,   174,   179,
     180,   181,   182,   183,   184,   185,   186,   187,   188,   189,
     191,   193,   196,   199,   202,   205,   207,   210,   213,   216,
     219,   225,   226,   227,   228,   230,   232,   234,   236,   238,
     240,   242,   243,   244,   245,   252,   259,   260,   261,   264,
     265,   270,   271
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUM", "QSTRING", "UNDVAR", "VAR",
  "SVAR", "IMMVAR", "IMMSVAR", "FNCT", "SFNCT", "'='", "EQ_MINUS",
  "EQ_PLUS", "EQ_DIV", "EQ_TIME", "EQ_POW", "'?'", "':'", "LOR", "LAND",
  "'<'", "'>'", "NE", "EQ", "GE", "LE", "'-'", "'+'", "'/'", "'*'", "'%'",
  "NOT", "UNARY", "POW", "DEC", "INC", "CONCAT", "'\\n'", "'{'", "'}'",
  "'('", "')'", "','", "';'", "'['", "']'", "$accept", "input", "line",
  "bool", "sexp", "exp", YY_NULL
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,    61,   267,   268,   269,   270,   271,    63,    58,
     272,   273,    60,    62,   274,   275,   276,   277,    45,    43,
      47,    42,    37,   278,   279,   280,   281,   282,   283,    10,
     123,   125,    40,    41,    44,    59,    91,    93
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    48,    49,    49,    50,    50,    50,    50,    51,    51,
      51,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    51,    51,    51,    51,    52,    52,    52,    52,
      52,    52,    52,    52,    52,    52,    52,    52,    52,    52,
      52,    52,    52,    53,    53,    53,    53,    53,    53,    53,
      53,    53,    53,    53,    53,    53,    53,    53,    53,    53,
      53,    53,    53,    53,    53,    53,    53,    53,    53,    53,
      53,    53,    53,    53,    53,    53,    53,    53,    53,    53,
      53,    53,    53,    53,    53,    53,    53,    53,    53,    53,
      53,    53,    53,    53,    53,    53,    53,    53,    53,    53,
      53,    53,    53
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     3,     3,     2,     3,     3,
       2,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     1,     1,     1,     3,
       3,     3,     3,     3,     4,     3,     4,     3,    12,     8,
       8,     8,     5,     1,     2,     2,     1,     1,     2,     2,
       2,     2,     3,     3,     3,     3,     3,     3,     3,     2,
       2,     2,     2,     3,     3,     3,     3,     3,     3,     3,
       1,     2,     2,     2,     2,     3,     3,     3,     3,     3,
       3,     3,     4,     4,     6,     6,     6,     8,    10,    10,
      14,     3,     3,     3,     3,     3,     2,     2,     3,     3,
       3,     1,     5
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,     4,     0,     3,     7,    43,    26,
      70,    46,    27,    47,    28,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   101,     0,     0,     0,     0,     0,
       0,     0,     0,    74,    73,     0,     0,     0,     0,     0,
       0,    51,    50,     0,     0,     0,     0,     0,     0,     0,
      62,    61,     0,     0,     0,     0,    96,    97,    10,    45,
      72,    49,    60,    44,    71,    48,    59,   101,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       6,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     5,    29,    75,    77,    76,
      79,    78,    80,    31,    52,    55,    54,    57,    56,    58,
      30,    53,    33,    63,    66,    65,    68,    67,    69,    32,
      64,    81,     0,     0,    35,     0,     0,    19,    99,   100,
       0,     0,    17,     0,    18,    20,    21,    25,    24,    23,
      22,    37,    15,    16,     8,     9,    14,    13,    12,    11,
      92,    91,    94,    93,    95,    98,    83,     0,    82,     0,
      34,     0,    36,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    42,   102,    85,    84,    86,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    87,     0,     0,
      40,    41,    39,     0,     0,     0,     0,    89,     0,    88,
       0,     0,     0,     0,    38,     0,    90
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,     6,    24,    55,   133
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -40
static const yytype_int16 yypact[] =
{
     -40,     6,   -40,   -39,   -40,   261,   -40,   -40,   -40,   -40,
      57,   238,    -3,   296,    -2,   -31,   -26,   261,   261,   261,
     142,   213,   261,   261,    79,   696,   648,   261,   261,   261,
     261,   261,   261,   -40,   -40,   261,   261,   261,   261,   261,
     261,   -40,   -40,   261,   261,   261,   261,   261,   261,   261,
     -40,   -40,   261,   158,   203,   235,   -20,   -20,   -20,   -40,
     -40,   -40,   -40,   -40,   -40,   -40,   -40,    92,    32,    99,
     261,   261,   261,   261,   261,   261,   261,   261,   261,   261,
     -40,   261,   261,   261,   261,   261,   261,   261,   261,   261,
     261,   261,   261,   261,   261,   -40,   235,   741,   741,   741,
     741,   741,   741,   235,   741,   741,   741,   741,   741,   741,
     235,   741,   235,   741,   741,   741,   741,   741,   741,   235,
     741,   -40,   555,   294,   -40,   563,   319,   -40,   -40,   -40,
     705,   725,    -1,   741,   -40,   -19,   -19,   -19,   -19,   -19,
     -19,   -40,   756,   770,    73,    73,    73,    73,    73,    73,
     -27,   -27,   -20,   -20,   -20,   -20,   -40,   261,   -40,   261,
     -40,   261,   -40,   261,   261,   261,   154,   469,   256,   593,
     394,   616,   235,   741,   -40,   -40,   -40,   261,   261,   261,
     261,   261,   344,   419,   668,   493,   586,   -40,   261,   261,
     -40,   -40,   -40,   261,   369,   517,   623,   -40,   261,   -40,
     261,   444,   674,   261,   -40,   541,   -40
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -40,   -40,   -40,   -21,    63,    -5
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      26,    67,     7,    91,    92,    93,     2,     3,    94,    43,
      52,    53,    56,    57,    58,    94,    54,    68,    69,    79,
      72,     0,    97,    98,    99,   100,   101,   102,     0,     0,
     104,   105,   106,   107,   108,   109,     0,     0,   111,   113,
     114,   115,   116,   117,   118,     4,     5,   120,   123,   126,
     132,   134,    81,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    91,    92,    93,   131,     0,    94,    25,    27,
      28,    29,    30,    31,    32,   128,   142,   143,   144,   145,
     146,   147,   148,   149,   150,   151,   152,   153,   154,   155,
      96,     0,     0,    33,    34,     0,     0,    70,   103,    71,
      72,    89,    90,    91,    92,    93,   110,   112,    94,     0,
      70,     0,    71,    72,     0,   119,   122,   125,     0,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      92,    93,     0,   130,    94,   127,   135,   136,   137,   138,
     139,   140,   141,     0,     0,    59,   129,    60,    61,     0,
      62,     0,   167,     0,   168,     0,   170,     0,     0,     0,
     173,     8,     9,    10,    11,    12,    13,    14,    15,    16,
       0,     0,   182,   183,     0,   185,    73,    74,    75,    76,
      77,    78,     0,   194,   195,     0,    17,    18,     0,     0,
       0,    19,    79,   201,    20,    21,     0,   174,   205,     0,
      22,   121,     0,     0,    23,     0,     8,     9,    10,    11,
      12,    13,    14,    15,    16,     0,    63,     0,    64,    65,
     166,    66,     0,     0,   169,     0,   171,   172,     0,     0,
       0,    17,    18,     0,     0,     0,    19,     0,     0,    20,
      21,     0,   184,     0,   186,    22,   124,     0,     0,    23,
      35,    36,    37,    38,    39,    40,   196,    73,    74,    75,
      76,    77,    78,   202,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    79,    41,    42,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    17,
      18,    94,     0,     0,    19,     0,     0,    20,    21,   176,
     177,   178,     0,    22,     0,     0,     0,    23,    44,    45,
      46,    47,    48,    49,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,     0,     0,    94,
       0,     0,    50,    51,     0,     0,     0,   158,   159,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      92,    93,     0,     0,    94,     0,     0,     0,     0,     0,
       0,     0,   162,   163,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,     0,     0,    94,
       0,     0,     0,     0,     0,     0,     0,   187,   188,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      92,    93,     0,     0,    94,     0,     0,     0,     0,     0,
       0,     0,   197,   198,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,     0,     0,    94,
       0,     0,     0,     0,     0,     0,     0,     0,   180,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      92,    93,     0,     0,    94,     0,     0,     0,     0,     0,
       0,     0,     0,   189,    81,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,     0,     0,    94,
       0,     0,     0,     0,     0,     0,     0,     0,   203,    81,
      82,    83,    84,    85,    86,    87,    88,    89,    90,    91,
      92,    93,     0,     0,    94,     0,     0,     0,     0,     0,
       0,     0,   175,    81,    82,    83,    84,    85,    86,    87,
      88,    89,    90,    91,    92,    93,     0,     0,    94,     0,
       0,     0,     0,     0,     0,     0,   191,    81,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
       0,     0,    94,     0,     0,     0,     0,     0,     0,     0,
     199,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,     0,     0,    94,    73,    74,    75,
      76,    77,    78,     0,   206,    73,    74,    75,    76,    77,
      78,     0,     0,    79,     0,     0,     0,     0,   156,   157,
       0,    79,     0,     0,     0,     0,   160,   161,    73,    74,
      75,    76,    77,    78,     0,    73,    74,    75,    76,    77,
      78,     0,     0,     0,    79,     0,     0,     0,     0,   192,
     193,    79,     0,     0,     0,     0,     0,   179,    73,    74,
      75,    76,    77,    78,     0,    73,    74,    75,    76,    77,
      78,     0,     0,     0,    79,     0,     0,     0,     0,     0,
     181,    79,     0,     0,     0,     0,     0,   200,    81,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    92,
      93,     0,     0,    94,     0,     0,     0,     0,     0,    95,
      73,    74,    75,    76,    77,    78,    73,    74,    75,    76,
      77,    78,     0,     0,     0,     0,    79,     0,     0,     0,
       0,   190,    79,     0,     0,     0,     0,   204,    73,    74,
      75,    76,    77,    78,   164,     0,     0,    73,    74,    75,
      76,    77,    78,     0,    79,     0,     0,    80,     0,     0,
       0,     0,     0,    79,   165,    81,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    92,    93,     0,     0,
      94,    81,    82,    83,    84,    85,    86,    87,    88,    89,
      90,    91,    92,    93,     0,     0,    94,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,     0,
       0,    94,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    92,    93,     0,     0,    94
};

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-40)))

#define yytable_value_is_error(Yytable_value) \
  YYID (0)

static const yytype_int16 yycheck[] =
{
       5,    22,    41,    30,    31,    32,     0,     1,    35,    12,
      12,    42,    17,    18,    19,    35,    42,    22,    23,    38,
      21,    -1,    27,    28,    29,    30,    31,    32,    -1,    -1,
      35,    36,    37,    38,    39,    40,    -1,    -1,    43,    44,
      45,    46,    47,    48,    49,    39,    40,    52,    53,    54,
      71,    72,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    70,    -1,    35,     5,    12,
      13,    14,    15,    16,    17,    43,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      27,    -1,    -1,    36,    37,    -1,    -1,    18,    35,    20,
      21,    28,    29,    30,    31,    32,    43,    44,    35,    -1,
      18,    -1,    20,    21,    -1,    52,    53,    54,    -1,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    -1,    70,    35,    43,    73,    74,    75,    76,
      77,    78,    79,    -1,    -1,     3,    47,     5,     6,    -1,
       8,    -1,   157,    -1,   159,    -1,   161,    -1,    -1,    -1,
     165,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      -1,    -1,   177,   178,    -1,   180,    22,    23,    24,    25,
      26,    27,    -1,   188,   189,    -1,    28,    29,    -1,    -1,
      -1,    33,    38,   198,    36,    37,    -1,    43,   203,    -1,
      42,    43,    -1,    -1,    46,    -1,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    -1,     3,    -1,     5,     6,
     157,     8,    -1,    -1,   161,    -1,   163,   164,    -1,    -1,
      -1,    28,    29,    -1,    -1,    -1,    33,    -1,    -1,    36,
      37,    -1,   179,    -1,   181,    42,    43,    -1,    -1,    46,
      12,    13,    14,    15,    16,    17,   193,    22,    23,    24,
      25,    26,    27,   200,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    38,    36,    37,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    28,
      29,    35,    -1,    -1,    33,    -1,    -1,    36,    37,    43,
      44,    45,    -1,    42,    -1,    -1,    -1,    46,    12,    13,
      14,    15,    16,    17,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    -1,    35,
      -1,    -1,    36,    37,    -1,    -1,    -1,    43,    44,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    -1,    -1,    35,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    43,    44,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    -1,    35,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    43,    44,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    -1,    -1,    35,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    43,    44,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    -1,    35,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    44,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    -1,    -1,    35,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    44,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    -1,    35,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    44,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    -1,    -1,    35,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    43,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    -1,    -1,    35,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    43,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      -1,    -1,    35,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      43,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    -1,    -1,    35,    22,    23,    24,
      25,    26,    27,    -1,    43,    22,    23,    24,    25,    26,
      27,    -1,    -1,    38,    -1,    -1,    -1,    -1,    43,    44,
      -1,    38,    -1,    -1,    -1,    -1,    43,    44,    22,    23,
      24,    25,    26,    27,    -1,    22,    23,    24,    25,    26,
      27,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,    43,
      44,    38,    -1,    -1,    -1,    -1,    -1,    44,    22,    23,
      24,    25,    26,    27,    -1,    22,    23,    24,    25,    26,
      27,    -1,    -1,    -1,    38,    -1,    -1,    -1,    -1,    -1,
      44,    38,    -1,    -1,    -1,    -1,    -1,    44,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    -1,    -1,    35,    -1,    -1,    -1,    -1,    -1,    41,
      22,    23,    24,    25,    26,    27,    22,    23,    24,    25,
      26,    27,    -1,    -1,    -1,    -1,    38,    -1,    -1,    -1,
      -1,    43,    38,    -1,    -1,    -1,    -1,    43,    22,    23,
      24,    25,    26,    27,    19,    -1,    -1,    22,    23,    24,
      25,    26,    27,    -1,    38,    -1,    -1,    41,    -1,    -1,
      -1,    -1,    -1,    38,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    -1,    -1,
      35,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    -1,    -1,    35,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    -1,
      -1,    35,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    -1,    -1,    35
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    49,     0,     1,    39,    40,    50,    41,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    28,    29,    33,
      36,    37,    42,    46,    51,    52,    53,    12,    13,    14,
      15,    16,    17,    36,    37,    12,    13,    14,    15,    16,
      17,    36,    37,    12,    12,    13,    14,    15,    16,    17,
      36,    37,    12,    42,    42,    52,    53,    53,    53,     3,
       5,     6,     8,     3,     5,     6,     8,    51,    53,    53,
      18,    20,    21,    22,    23,    24,    25,    26,    27,    38,
      41,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    35,    41,    52,    53,    53,    53,
      53,    53,    53,    52,    53,    53,    53,    53,    53,    53,
      52,    53,    52,    53,    53,    53,    53,    53,    53,    52,
      53,    43,    52,    53,    43,    52,    53,    43,    43,    47,
      52,    53,    51,    53,    51,    52,    52,    52,    52,    52,
      52,    52,    53,    53,    53,    53,    53,    53,    53,    53,
      53,    53,    53,    53,    53,    53,    43,    44,    43,    44,
      43,    44,    43,    44,    19,    19,    52,    53,    53,    52,
      53,    52,    52,    53,    43,    43,    43,    44,    45,    44,
      44,    44,    53,    53,    52,    53,    52,    43,    44,    44,
      43,    43,    43,    44,    53,    53,    52,    43,    44,    43,
      44,    53,    52,    44,    43,    53,    43
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

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
      YYERROR;							\
    }								\
while (YYID (0))

/* Error token number */
#define YYTERROR	1
#define YYERRCODE	256


/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */
#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

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
#ifndef	YYINITDEPTH
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
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
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
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
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
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
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
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
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

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YYUSE (yytype);
}




/* The lookahead symbol.  */
int yychar;


#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval YY_INITIAL_VALUE(yyval_default);

/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

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
      yychar = YYLEX;
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
     `$$ = $1'.

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
/* Line 1787 of yacc.c  */
#line 93 "aprepro.y"
    { if (echo) fprintf(yyout,"\n");	}
    break;

  case 5:
/* Line 1787 of yacc.c  */
#line 94 "aprepro.y"
    { if (echo) {
	                             format = getsym("_FORMAT");
	                             fprintf(yyout, format->value.svar, (yyvsp[(2) - (3)].val));
				   }                                    }
    break;

  case 6:
/* Line 1787 of yacc.c  */
#line 98 "aprepro.y"
    { if (echo && (yyvsp[(2) - (3)].string) != NULL)
				    fprintf(yyout, "%s", (yyvsp[(2) - (3)].string));	}
    break;

  case 7:
/* Line 1787 of yacc.c  */
#line 100 "aprepro.y"
    { yyerrok;				}
    break;

  case 8:
/* Line 1787 of yacc.c  */
#line 103 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) < (yyvsp[(3) - (3)].val);                         }
    break;

  case 9:
/* Line 1787 of yacc.c  */
#line 104 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) > (yyvsp[(3) - (3)].val);                         }
    break;

  case 10:
/* Line 1787 of yacc.c  */
#line 105 "aprepro.y"
    { (yyval.val) = !((yyvsp[(2) - (2)].val));                           }
    break;

  case 11:
/* Line 1787 of yacc.c  */
#line 106 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) <= (yyvsp[(3) - (3)].val);                        }
    break;

  case 12:
/* Line 1787 of yacc.c  */
#line 107 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) >= (yyvsp[(3) - (3)].val);                        }
    break;

  case 13:
/* Line 1787 of yacc.c  */
#line 108 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) == (yyvsp[(3) - (3)].val);                        }
    break;

  case 14:
/* Line 1787 of yacc.c  */
#line 109 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) != (yyvsp[(3) - (3)].val);                        }
    break;

  case 15:
/* Line 1787 of yacc.c  */
#line 110 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) || (yyvsp[(3) - (3)].val);                        }
    break;

  case 16:
/* Line 1787 of yacc.c  */
#line 111 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) && (yyvsp[(3) - (3)].val);                        }
    break;

  case 17:
/* Line 1787 of yacc.c  */
#line 112 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) || (yyvsp[(3) - (3)].val);                        }
    break;

  case 18:
/* Line 1787 of yacc.c  */
#line 113 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) && (yyvsp[(3) - (3)].val);                        }
    break;

  case 19:
/* Line 1787 of yacc.c  */
#line 114 "aprepro.y"
    { (yyval.val) = (yyvsp[(2) - (3)].val);                              }
    break;

  case 20:
/* Line 1787 of yacc.c  */
#line 117 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) <  0 ? 1 : 0);	}
    break;

  case 21:
/* Line 1787 of yacc.c  */
#line 118 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) >  0 ? 1 : 0);	}
    break;

  case 22:
/* Line 1787 of yacc.c  */
#line 119 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) <= 0 ? 1 : 0);	}
    break;

  case 23:
/* Line 1787 of yacc.c  */
#line 120 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) >= 0 ? 1 : 0);	}
    break;

  case 24:
/* Line 1787 of yacc.c  */
#line 121 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) == 0 ? 1 : 0);	}
    break;

  case 25:
/* Line 1787 of yacc.c  */
#line 122 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) != 0 ? 1 : 0);	}
    break;

  case 26:
/* Line 1787 of yacc.c  */
#line 124 "aprepro.y"
    { (yyval.string) = (yyvsp[(1) - (1)].string);				}
    break;

  case 27:
/* Line 1787 of yacc.c  */
#line 125 "aprepro.y"
    { (yyval.string) = (yyvsp[(1) - (1)].tptr)->value.svar;			}
    break;

  case 28:
/* Line 1787 of yacc.c  */
#line 126 "aprepro.y"
    { (yyval.string) = (yyvsp[(1) - (1)].tptr)->value.svar;			}
    break;

  case 29:
/* Line 1787 of yacc.c  */
#line 127 "aprepro.y"
    { (yyval.string) = (yyvsp[(3) - (3)].string); (yyvsp[(1) - (3)].tptr)->value.svar = (yyvsp[(3) - (3)].string);
	                          set_type((yyvsp[(1) - (3)].tptr), SVAR);			}
    break;

  case 30:
/* Line 1787 of yacc.c  */
#line 129 "aprepro.y"
    { (yyval.string) = (yyvsp[(3) - (3)].string); 
                                  (yyvsp[(1) - (3)].tptr)->value.svar = (yyvsp[(3) - (3)].string);
                                  redefined_warning((yyvsp[(1) - (3)].tptr));          }
    break;

  case 31:
/* Line 1787 of yacc.c  */
#line 132 "aprepro.y"
    { (yyval.string) = (yyvsp[(3) - (3)].string); 
                                  (yyvsp[(1) - (3)].tptr)->value.svar= (yyvsp[(3) - (3)].string);
                                  redefined_warning((yyvsp[(1) - (3)].tptr));          
                                  (yyvsp[(1) - (3)].tptr)->type = SVAR;              }
    break;

  case 32:
/* Line 1787 of yacc.c  */
#line 136 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (3)].tptr)); YYERROR; }
    break;

  case 33:
/* Line 1787 of yacc.c  */
#line 137 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (3)].tptr)); YYERROR; }
    break;

  case 34:
/* Line 1787 of yacc.c  */
#line 138 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (4)].tptr)->value.strfnct))((yyvsp[(3) - (4)].string));	}
    break;

  case 35:
/* Line 1787 of yacc.c  */
#line 139 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (3)].tptr)->value.strfnct))();	}
    break;

  case 36:
/* Line 1787 of yacc.c  */
#line 140 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (4)].tptr)->value.strfnct))((yyvsp[(3) - (4)].val));	}
    break;

  case 37:
/* Line 1787 of yacc.c  */
#line 141 "aprepro.y"
    { int len1 = strlen((yyvsp[(1) - (3)].string));
				  int len3 = strlen((yyvsp[(3) - (3)].string));
				  ALLOC((yyval.string), len1+len3+1, char *);
				  memcpy((yyval.string), (yyvsp[(1) - (3)].string), len1+1);
				  strcat((yyval.string), (yyvsp[(3) - (3)].string)); }
    break;

  case 38:
/* Line 1787 of yacc.c  */
#line 147 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (12)].tptr)->value.strfnct))((yyvsp[(3) - (12)].val), (yyvsp[(5) - (12)].string), (yyvsp[(7) - (12)].string), (yyvsp[(9) - (12)].string), (yyvsp[(11) - (12)].string)); }
    break;

  case 39:
/* Line 1787 of yacc.c  */
#line 149 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (8)].tptr)->value.strfnct))((yyvsp[(3) - (8)].val), (yyvsp[(5) - (8)].string), (yyvsp[(7) - (8)].string)); }
    break;

  case 40:
/* Line 1787 of yacc.c  */
#line 151 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (8)].tptr)->value.strfnct))((yyvsp[(3) - (8)].string), (yyvsp[(5) - (8)].string), (yyvsp[(7) - (8)].string)); }
    break;

  case 41:
/* Line 1787 of yacc.c  */
#line 153 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (8)].tptr)->value.strfnct))((yyvsp[(3) - (8)].string), (yyvsp[(5) - (8)].val), (yyvsp[(7) - (8)].val)); }
    break;

  case 42:
/* Line 1787 of yacc.c  */
#line 154 "aprepro.y"
    { (yyval.string) = ((yyvsp[(1) - (5)].val)) ? ((yyvsp[(3) - (5)].string)) : ((yyvsp[(5) - (5)].string));              }
    break;

  case 43:
/* Line 1787 of yacc.c  */
#line 156 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); 				}
    break;

  case 44:
/* Line 1787 of yacc.c  */
#line 157 "aprepro.y"
    { (yyval.val) = (yyvsp[(2) - (2)].val) + 1;				}
    break;

  case 45:
/* Line 1787 of yacc.c  */
#line 158 "aprepro.y"
    { (yyval.val) = (yyvsp[(2) - (2)].val) - 1;				}
    break;

  case 46:
/* Line 1787 of yacc.c  */
#line 159 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (1)].tptr)->value.var;			}
    break;

  case 47:
/* Line 1787 of yacc.c  */
#line 160 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (1)].tptr)->value.var;			}
    break;

  case 48:
/* Line 1787 of yacc.c  */
#line 161 "aprepro.y"
    { (yyval.val) = ++((yyvsp[(2) - (2)].tptr)->value.var);		}
    break;

  case 49:
/* Line 1787 of yacc.c  */
#line 162 "aprepro.y"
    { (yyval.val) = --((yyvsp[(2) - (2)].tptr)->value.var);		}
    break;

  case 50:
/* Line 1787 of yacc.c  */
#line 163 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (2)].tptr)->value.var)++;		}
    break;

  case 51:
/* Line 1787 of yacc.c  */
#line 164 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (2)].tptr)->value.var)--;		}
    break;

  case 52:
/* Line 1787 of yacc.c  */
#line 165 "aprepro.y"
    { (yyval.val) = (yyvsp[(3) - (3)].val); (yyvsp[(1) - (3)].tptr)->value.var = (yyvsp[(3) - (3)].val);
				  redefined_warning((yyvsp[(1) - (3)].tptr));          }
    break;

  case 53:
/* Line 1787 of yacc.c  */
#line 167 "aprepro.y"
    { (yyval.val) = (yyvsp[(3) - (3)].val); (yyvsp[(1) - (3)].tptr)->value.var = (yyvsp[(3) - (3)].val);
				  redefined_warning((yyvsp[(1) - (3)].tptr));          
				  (yyvsp[(1) - (3)].tptr)->type = VAR;			}
    break;

  case 54:
/* Line 1787 of yacc.c  */
#line 170 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var += (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; }
    break;

  case 55:
/* Line 1787 of yacc.c  */
#line 171 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var -= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; }
    break;

  case 56:
/* Line 1787 of yacc.c  */
#line 172 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var *= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; }
    break;

  case 57:
/* Line 1787 of yacc.c  */
#line 173 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var /= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; }
    break;

  case 58:
/* Line 1787 of yacc.c  */
#line 174 "aprepro.y"
    { errno = 0;
				  (yyvsp[(1) - (3)].tptr)->value.var = pow((yyvsp[(1) - (3)].tptr)->value.var,(yyvsp[(3) - (3)].val)); 
				  (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
				  MATH_ERROR("Power");
				}
    break;

  case 59:
/* Line 1787 of yacc.c  */
#line 179 "aprepro.y"
    { immutable_modify((yyvsp[(2) - (2)].tptr)); YYERROR; }
    break;

  case 60:
/* Line 1787 of yacc.c  */
#line 180 "aprepro.y"
    { immutable_modify((yyvsp[(2) - (2)].tptr)); YYERROR; }
    break;

  case 61:
/* Line 1787 of yacc.c  */
#line 181 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (2)].tptr)); YYERROR; }
    break;

  case 62:
/* Line 1787 of yacc.c  */
#line 182 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (2)].tptr)); YYERROR; }
    break;

  case 63:
/* Line 1787 of yacc.c  */
#line 183 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (3)].tptr)); YYERROR; }
    break;

  case 64:
/* Line 1787 of yacc.c  */
#line 184 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (3)].tptr)); YYERROR; }
    break;

  case 65:
/* Line 1787 of yacc.c  */
#line 185 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (3)].tptr)); YYERROR; }
    break;

  case 66:
/* Line 1787 of yacc.c  */
#line 186 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (3)].tptr)); YYERROR; }
    break;

  case 67:
/* Line 1787 of yacc.c  */
#line 187 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (3)].tptr)); YYERROR; }
    break;

  case 68:
/* Line 1787 of yacc.c  */
#line 188 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (3)].tptr)); YYERROR; }
    break;

  case 69:
/* Line 1787 of yacc.c  */
#line 189 "aprepro.y"
    { immutable_modify((yyvsp[(1) - (3)].tptr)); YYERROR; }
    break;

  case 70:
/* Line 1787 of yacc.c  */
#line 191 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (1)].tptr)->value.var;
				  undefined_warning((yyvsp[(1) - (1)].tptr)->name);          }
    break;

  case 71:
/* Line 1787 of yacc.c  */
#line 193 "aprepro.y"
    { (yyval.val) = ++((yyvsp[(2) - (2)].tptr)->value.var);		
	                          set_type((yyvsp[(2) - (2)].tptr), VAR);                       
				  undefined_warning((yyvsp[(2) - (2)].tptr)->name);          }
    break;

  case 72:
/* Line 1787 of yacc.c  */
#line 196 "aprepro.y"
    { (yyval.val) = --((yyvsp[(2) - (2)].tptr)->value.var);		
	                          set_type((yyvsp[(2) - (2)].tptr), VAR);                       
				  undefined_warning((yyvsp[(2) - (2)].tptr)->name);          }
    break;

  case 73:
/* Line 1787 of yacc.c  */
#line 199 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (2)].tptr)->value.var)++;		
	                          set_type((yyvsp[(1) - (2)].tptr), VAR);                       
				  undefined_warning((yyvsp[(1) - (2)].tptr)->name);          }
    break;

  case 74:
/* Line 1787 of yacc.c  */
#line 202 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (2)].tptr)->value.var)--;		
	                          set_type((yyvsp[(1) - (2)].tptr), VAR);                       
				  undefined_warning((yyvsp[(1) - (2)].tptr)->name);          }
    break;

  case 75:
/* Line 1787 of yacc.c  */
#line 205 "aprepro.y"
    { (yyval.val) = (yyvsp[(3) - (3)].val); (yyvsp[(1) - (3)].tptr)->value.var = (yyvsp[(3) - (3)].val);
	                          set_type((yyvsp[(1) - (3)].tptr), VAR);                    }
    break;

  case 76:
/* Line 1787 of yacc.c  */
#line 207 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var += (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
	                          set_type((yyvsp[(1) - (3)].tptr), VAR);                       
				  undefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 77:
/* Line 1787 of yacc.c  */
#line 210 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var -= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
	                          set_type((yyvsp[(1) - (3)].tptr), VAR);                       
				  undefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 78:
/* Line 1787 of yacc.c  */
#line 213 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var *= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
	                          set_type((yyvsp[(1) - (3)].tptr), VAR);                       
				  undefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 79:
/* Line 1787 of yacc.c  */
#line 216 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var /= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
	                          set_type((yyvsp[(1) - (3)].tptr), VAR);                       
				  undefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 80:
/* Line 1787 of yacc.c  */
#line 219 "aprepro.y"
    { errno = 0;
				  (yyvsp[(1) - (3)].tptr)->value.var = pow((yyvsp[(1) - (3)].tptr)->value.var,(yyvsp[(3) - (3)].val)); 
				  (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
	                          set_type((yyvsp[(1) - (3)].tptr), VAR);                       
				  MATH_ERROR("Power");
				  undefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 81:
/* Line 1787 of yacc.c  */
#line 225 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (3)].tptr)->value.fnctptr))();	}
    break;

  case 82:
/* Line 1787 of yacc.c  */
#line 226 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (4)].tptr)->value.fnctptr))((yyvsp[(3) - (4)].val)); 	}
    break;

  case 83:
/* Line 1787 of yacc.c  */
#line 227 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (4)].tptr)->value.fnctptr))((yyvsp[(3) - (4)].string)); 	}
    break;

  case 84:
/* Line 1787 of yacc.c  */
#line 229 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (6)].tptr)->value.fnctptr))((yyvsp[(3) - (6)].string), (yyvsp[(5) - (6)].val)); 	}
    break;

  case 85:
/* Line 1787 of yacc.c  */
#line 231 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (6)].tptr)->value.fnctptr))((yyvsp[(3) - (6)].string), (yyvsp[(5) - (6)].string)); 	}
    break;

  case 86:
/* Line 1787 of yacc.c  */
#line 233 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (6)].tptr)->value.fnctptr))((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));	}
    break;

  case 87:
/* Line 1787 of yacc.c  */
#line 235 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (8)].tptr)->value.fnctptr))((yyvsp[(3) - (8)].val), (yyvsp[(5) - (8)].val), (yyvsp[(7) - (8)].val)); }
    break;

  case 88:
/* Line 1787 of yacc.c  */
#line 237 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (10)].tptr)->value.fnctptr))((yyvsp[(3) - (10)].val), (yyvsp[(5) - (10)].val), (yyvsp[(7) - (10)].val), (yyvsp[(9) - (10)].val)); }
    break;

  case 89:
/* Line 1787 of yacc.c  */
#line 239 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (10)].tptr)->value.fnctptr))((yyvsp[(3) - (10)].val), (yyvsp[(5) - (10)].val), (yyvsp[(7) - (10)].val), (yyvsp[(9) - (10)].val)); }
    break;

  case 90:
/* Line 1787 of yacc.c  */
#line 241 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (14)].tptr)->value.fnctptr))((yyvsp[(3) - (14)].val), (yyvsp[(5) - (14)].val), (yyvsp[(7) - (14)].val), (yyvsp[(9) - (14)].val), (yyvsp[(11) - (14)].val), (yyvsp[(13) - (14)].val)); }
    break;

  case 91:
/* Line 1787 of yacc.c  */
#line 242 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) + (yyvsp[(3) - (3)].val); 			}
    break;

  case 92:
/* Line 1787 of yacc.c  */
#line 243 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) - (yyvsp[(3) - (3)].val); 			}
    break;

  case 93:
/* Line 1787 of yacc.c  */
#line 244 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) * (yyvsp[(3) - (3)].val); 			}
    break;

  case 94:
/* Line 1787 of yacc.c  */
#line 245 "aprepro.y"
    { if ((yyvsp[(3) - (3)].val) == 0.)
				    {
				      yyerror("Zero divisor"); 
				      yyerrok;
				    }
				  else
				    (yyval.val) = (yyvsp[(1) - (3)].val) / (yyvsp[(3) - (3)].val); 			}
    break;

  case 95:
/* Line 1787 of yacc.c  */
#line 252 "aprepro.y"
    { if ((yyvsp[(3) - (3)].val) == 0.)
				    {
				      yyerror("Zero divisor");
				      yyerrok;
				    }
				  else
				    (yyval.val) = (int)(yyvsp[(1) - (3)].val) % (int)(yyvsp[(3) - (3)].val);		}
    break;

  case 96:
/* Line 1787 of yacc.c  */
#line 259 "aprepro.y"
    { (yyval.val) = -(yyvsp[(2) - (2)].val);				}
    break;

  case 97:
/* Line 1787 of yacc.c  */
#line 260 "aprepro.y"
    { (yyval.val) =  (yyvsp[(2) - (2)].val);				}
    break;

  case 98:
/* Line 1787 of yacc.c  */
#line 261 "aprepro.y"
    { errno = 0;
				  (yyval.val) = pow((yyvsp[(1) - (3)].val), (yyvsp[(3) - (3)].val)); 
				  MATH_ERROR("Power");			}
    break;

  case 99:
/* Line 1787 of yacc.c  */
#line 264 "aprepro.y"
    { (yyval.val) = (yyvsp[(2) - (3)].val);				}
    break;

  case 100:
/* Line 1787 of yacc.c  */
#line 265 "aprepro.y"
    { errno = 0;
				  (yyval.val) = (double)((yyvsp[(2) - (3)].val) < 0 ? 
					-floor(-((yyvsp[(2) - (3)].val))): floor((yyvsp[(2) - (3)].val)) );
				  MATH_ERROR("floor (int)");		}
    break;

  case 101:
/* Line 1787 of yacc.c  */
#line 270 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (1)].val)) ? 1 : 0; }
    break;

  case 102:
/* Line 1787 of yacc.c  */
#line 271 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (5)].val)) ? ((yyvsp[(3) - (5)].val)) : ((yyvsp[(5) - (5)].val));              }
    break;


/* Line 1787 of yacc.c  */
#line 2398 "y.tab.c"
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

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
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

  /* Do not reclaim the symbols of the rule which action triggered
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
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

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
  /* Do not reclaim the symbols of the rule which action triggered
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
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


/* Line 2050 of yacc.c  */
#line 275 "aprepro.y"

# include "lex.yy.c"


