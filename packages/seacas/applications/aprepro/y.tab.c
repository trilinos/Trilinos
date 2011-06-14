/* A Bison parser, made by GNU Bison 2.4.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006,
   2009, 2010 Free Software Foundation, Inc.
   
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
#define YYBISON_VERSION "2.4.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
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

void undefined_warning(char* var);
void redefined_warning(char* var);
void yyerror(char* var);
void warning(char *string);
int  yylex(void);

int   echo = True;

symrec *format;


/* Line 189 of yacc.c  */
#line 122 "y.tab.c"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
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
     FNCT = 263,
     SFNCT = 264,
     EQ_MINUS = 265,
     EQ_PLUS = 266,
     EQ_DIV = 267,
     EQ_TIME = 268,
     EQ_POW = 269,
     LOR = 270,
     LAND = 271,
     NE = 272,
     EQ = 273,
     GE = 274,
     LE = 275,
     NOT = 276,
     UNARY = 277,
     POW = 278,
     DEC = 279,
     INC = 280,
     CONCAT = 281
   };
#endif
/* Tokens.  */
#define NUM 258
#define QSTRING 259
#define UNDVAR 260
#define VAR 261
#define SVAR 262
#define FNCT 263
#define SFNCT 264
#define EQ_MINUS 265
#define EQ_PLUS 266
#define EQ_DIV 267
#define EQ_TIME 268
#define EQ_POW 269
#define LOR 270
#define LAND 271
#define NE 272
#define EQ 273
#define GE 274
#define LE 275
#define NOT 276
#define UNARY 277
#define POW 278
#define DEC 279
#define INC 280
#define CONCAT 281




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 50 "aprepro.y"

  double  val;		/* For returning numbers.		*/
  symrec *tptr;		/* For returning symbol-table pointers	*/
  char   *string;	/* For returning quoted strings		*/



/* Line 214 of yacc.c  */
#line 218 "y.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 230 "y.tab.c"

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
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
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
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
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
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
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

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

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

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   752

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  46
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  6
/* YYNRULES -- Number of rules.  */
#define YYNRULES  86
/* YYNRULES -- Number of states.  */
#define YYNSTATES  181

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   281

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      37,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,    30,     2,     2,
      40,    41,    29,    27,    42,    26,     2,    28,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    17,    43,
      20,    10,    21,    16,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    44,     2,    45,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    38,     2,    39,     2,     2,     2,     2,
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
       5,     6,     7,     8,     9,    11,    12,    13,    14,    15,
      18,    19,    22,    23,    24,    25,    31,    32,    33,    34,
      35,    36
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    13,    17,    20,    24,
      28,    31,    35,    39,    43,    47,    51,    55,    59,    63,
      67,    71,    75,    79,    83,    87,    91,    93,    95,    99,
     103,   107,   112,   116,   121,   125,   138,   147,   156,   162,
     164,   167,   170,   172,   175,   178,   181,   184,   188,   192,
     196,   200,   204,   208,   212,   214,   217,   220,   223,   226,
     230,   234,   238,   242,   246,   250,   254,   259,   264,   271,
     278,   285,   294,   305,   316,   331,   335,   339,   343,   347,
     351,   354,   357,   361,   365,   369,   371
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      47,     0,    -1,    -1,    47,    48,    -1,    37,    -1,    38,
      51,    39,    -1,    38,    50,    39,    -1,     1,    39,    -1,
      51,    20,    51,    -1,    51,    21,    51,    -1,    31,    51,
      -1,    51,    25,    51,    -1,    51,    24,    51,    -1,    51,
      23,    51,    -1,    51,    22,    51,    -1,    51,    18,    51,
      -1,    51,    19,    51,    -1,    49,    18,    49,    -1,    49,
      19,    49,    -1,    40,    49,    41,    -1,    50,    20,    50,
      -1,    50,    21,    50,    -1,    50,    25,    50,    -1,    50,
      24,    50,    -1,    50,    23,    50,    -1,    50,    22,    50,
      -1,     4,    -1,     7,    -1,     5,    10,    50,    -1,     7,
      10,    50,    -1,     6,    10,    50,    -1,     9,    40,    50,
      41,    -1,     9,    40,    41,    -1,     9,    40,    51,    41,
      -1,    50,    36,    50,    -1,     9,    40,    51,    42,    50,
      42,    50,    42,    50,    42,    50,    41,    -1,     9,    40,
      51,    42,    50,    42,    50,    41,    -1,     9,    40,    50,
      42,    50,    42,    50,    41,    -1,    49,    16,    50,    17,
      50,    -1,     3,    -1,    35,     3,    -1,    34,     3,    -1,
       6,    -1,    35,     6,    -1,    34,     6,    -1,     6,    35,
      -1,     6,    34,    -1,     6,    10,    51,    -1,     7,    10,
      51,    -1,     6,    12,    51,    -1,     6,    11,    51,    -1,
       6,    14,    51,    -1,     6,    13,    51,    -1,     6,    15,
      51,    -1,     5,    -1,    35,     5,    -1,    34,     5,    -1,
       5,    35,    -1,     5,    34,    -1,     5,    10,    51,    -1,
       5,    12,    51,    -1,     5,    11,    51,    -1,     5,    14,
      51,    -1,     5,    13,    51,    -1,     5,    15,    51,    -1,
       8,    40,    41,    -1,     8,    40,    51,    41,    -1,     8,
      40,    50,    41,    -1,     8,    40,    50,    42,    51,    41,
      -1,     8,    40,    50,    42,    50,    41,    -1,     8,    40,
      51,    42,    51,    41,    -1,     8,    40,    51,    42,    51,
      42,    51,    41,    -1,     8,    40,    51,    42,    51,    43,
      51,    42,    51,    41,    -1,     8,    40,    51,    42,    51,
      42,    51,    42,    51,    41,    -1,     8,    40,    51,    42,
      51,    42,    51,    42,    51,    42,    51,    42,    51,    41,
      -1,    51,    27,    51,    -1,    51,    26,    51,    -1,    51,
      29,    51,    -1,    51,    28,    51,    -1,    51,    30,    51,
      -1,    26,    51,    -1,    27,    51,    -1,    51,    33,    51,
      -1,    40,    51,    41,    -1,    44,    51,    45,    -1,    49,
      -1,    49,    16,    51,    17,    51,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    85,    85,    86,    89,    90,    94,    96,    99,   100,
     101,   102,   103,   104,   105,   106,   107,   108,   109,   110,
     113,   114,   115,   116,   117,   118,   120,   121,   122,   124,
     127,   131,   132,   133,   134,   139,   141,   143,   145,   147,
     148,   149,   150,   151,   152,   153,   154,   155,   157,   160,
     161,   162,   163,   164,   169,   171,   174,   177,   180,   183,
     185,   188,   191,   194,   197,   203,   204,   205,   206,   208,
     210,   212,   214,   216,   218,   220,   221,   222,   223,   230,
     237,   238,   239,   242,   243,   248,   249
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUM", "QSTRING", "UNDVAR", "VAR",
  "SVAR", "FNCT", "SFNCT", "'='", "EQ_MINUS", "EQ_PLUS", "EQ_DIV",
  "EQ_TIME", "EQ_POW", "'?'", "':'", "LOR", "LAND", "'<'", "'>'", "NE",
  "EQ", "GE", "LE", "'-'", "'+'", "'/'", "'*'", "'%'", "NOT", "UNARY",
  "POW", "DEC", "INC", "CONCAT", "'\\n'", "'{'", "'}'", "'('", "')'",
  "','", "';'", "'['", "']'", "$accept", "input", "line", "bool", "sexp",
  "exp", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
      61,   265,   266,   267,   268,   269,    63,    58,   270,   271,
      60,    62,   272,   273,   274,   275,    45,    43,    47,    42,
      37,   276,   277,   278,   279,   280,   281,    10,   123,   125,
      40,    41,    44,    59,    91,    93
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    46,    47,    47,    48,    48,    48,    48,    49,    49,
      49,    49,    49,    49,    49,    49,    49,    49,    49,    49,
      49,    49,    49,    49,    49,    49,    50,    50,    50,    50,
      50,    50,    50,    50,    50,    50,    50,    50,    50,    51,
      51,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    51,    51,    51,    51,    51
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     3,     3,     2,     3,     3,
       2,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     1,     1,     3,     3,
       3,     4,     3,     4,     3,    12,     8,     8,     5,     1,
       2,     2,     1,     2,     2,     2,     2,     3,     3,     3,
       3,     3,     3,     3,     1,     2,     2,     2,     2,     3,
       3,     3,     3,     3,     3,     3,     4,     4,     6,     6,
       6,     8,    10,    10,    14,     3,     3,     3,     3,     3,
       2,     2,     3,     3,     3,     1,     5
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,     4,     0,     3,     7,    39,    26,
      54,    42,    27,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    85,     0,     0,     0,     0,     0,     0,     0,
       0,    58,    57,     0,     0,     0,     0,     0,     0,    46,
      45,     0,     0,     0,     0,    80,    81,    10,    41,    56,
      44,    40,    55,    43,    85,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     6,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     5,    28,    59,    61,    60,    63,    62,    64,
      30,    47,    50,    49,    52,    51,    53,    29,    48,    65,
       0,     0,    32,     0,     0,    19,    83,    84,     0,     0,
      17,     0,    18,    20,    21,    25,    24,    23,    22,    34,
      15,    16,     8,     9,    14,    13,    12,    11,    76,    75,
      78,    77,    79,    82,    67,     0,    66,     0,    31,     0,
      33,     0,     0,     0,     0,     0,     0,     0,     0,    38,
      86,    69,    68,    70,     0,     0,     0,     0,     0,     0,
       0,     0,    71,     0,     0,    37,    36,     0,     0,     0,
       0,    73,     0,    72,     0,     0,     0,     0,    35,     0,
      74
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,     6,    22,    44,   111
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -38
static const yytype_int16 yypact[] =
{
     -38,     4,   -38,   -37,   -38,   194,   -38,   -38,   -38,   -38,
      69,    76,    -7,   -33,   -27,   194,   194,   194,     3,    55,
     194,   194,    89,    23,   597,   194,   194,   194,   194,   194,
     194,   -38,   -38,   194,   194,   194,   194,   194,   194,   -38,
     -38,   194,    93,   139,   132,   -16,   -16,   -16,   -38,   -38,
     -38,   -38,   -38,   -38,    37,   425,   221,   194,   194,   194,
     194,   194,   194,   194,   194,   194,   194,   -38,   194,   194,
     194,   194,   194,   194,   194,   194,   194,   194,   194,   194,
     194,   194,   -38,   132,   690,   690,   690,   690,   690,   690,
     132,   690,   690,   690,   690,   690,   690,   132,   690,   -38,
     191,   275,   -38,   511,   300,   -38,   -38,   -38,   654,   674,
      -1,   690,   -38,   -22,   -22,   -22,   -22,   -22,   -22,   -38,
     705,   719,   229,   229,   229,   229,   229,   229,    21,    21,
     -16,   -16,   -16,   -16,   -38,   194,   -38,   194,   -38,   194,
     -38,   194,   194,   194,   617,   449,   249,   542,   549,   132,
     690,   -38,   -38,   -38,   194,   194,   194,   194,   325,   375,
     623,   519,   -38,   194,   194,   -38,   -38,   194,   350,   473,
     572,   -38,   194,   -38,   194,   400,   645,   194,   -38,   497,
     -38
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -38,   -38,   -38,   -19,    52,    -5
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      24,    54,     7,    41,     2,     3,    48,    42,    49,    50,
      45,    46,    47,    43,    66,    55,    56,    81,    59,     0,
      84,    85,    86,    87,    88,    89,     0,     0,    91,    92,
      93,    94,    95,    96,     0,     0,    98,   101,   104,   110,
     112,     4,     5,    60,    61,    62,    63,    64,    65,    78,
      79,    80,   109,    57,    81,    58,    59,    23,    51,    66,
      52,    53,    67,   120,   121,   122,   123,   124,   125,   126,
     127,   128,   129,   130,   131,   132,   133,    83,   105,    25,
      26,    27,    28,    29,    30,    90,    33,    34,    35,    36,
      37,    38,     0,    97,   100,   103,     8,     9,    10,    11,
      12,    13,    14,    31,    32,    57,     0,    58,    59,   108,
      39,    40,   113,   114,   115,   116,   117,   118,   119,    15,
      16,     0,     0,     0,    17,     0,     0,    18,    19,     0,
     145,     0,   146,    20,    99,     0,     0,    21,   150,     0,
       0,     0,     8,     9,    10,    11,    12,    13,    14,   158,
     159,     0,    60,    61,    62,    63,    64,    65,   168,   169,
       0,     0,     0,     0,     0,    15,    16,   175,    66,     0,
      17,     0,   179,    18,    19,     0,     0,     0,     0,    20,
     102,     0,     0,    21,     0,     0,     0,   144,     0,     0,
       0,   147,     0,   148,   149,     0,     0,     8,     9,    10,
      11,    12,    13,    14,     0,     0,     0,     0,   160,   161,
       0,    60,    61,    62,    63,    64,    65,     0,     0,   170,
      15,    16,     0,     0,     0,    17,   176,    66,    18,    19,
       0,     0,   134,   135,    20,     0,     0,     0,    21,    68,
      69,    70,    71,    72,    73,    74,    75,    76,    77,    78,
      79,    80,     0,     0,    81,    76,    77,    78,    79,    80,
       0,     0,    81,     0,     0,     0,   107,    68,    69,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
       0,     0,    81,     0,     0,     0,     0,     0,     0,     0,
     153,   154,   155,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,     0,     0,    81,     0,
       0,     0,     0,     0,     0,     0,   136,   137,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      80,     0,     0,    81,     0,     0,     0,     0,     0,     0,
       0,   140,   141,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,     0,     0,    81,     0,
       0,     0,     0,     0,     0,     0,   162,   163,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      80,     0,     0,    81,     0,     0,     0,     0,     0,     0,
       0,   171,   172,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,     0,     0,    81,     0,
       0,     0,     0,     0,     0,     0,     0,   164,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      80,     0,     0,    81,     0,     0,     0,     0,     0,     0,
       0,     0,   177,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,     0,     0,    81,     0,
       0,     0,     0,     0,     0,     0,   106,    68,    69,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
       0,     0,    81,     0,     0,     0,     0,     0,     0,     0,
     152,    68,    69,    70,    71,    72,    73,    74,    75,    76,
      77,    78,    79,    80,     0,     0,    81,     0,     0,     0,
       0,     0,     0,     0,   173,    68,    69,    70,    71,    72,
      73,    74,    75,    76,    77,    78,    79,    80,     0,     0,
      81,    60,    61,    62,    63,    64,    65,     0,   180,    60,
      61,    62,    63,    64,    65,     0,     0,    66,     0,     0,
       0,     0,   138,   139,     0,    66,     0,     0,     0,     0,
     166,   167,    60,    61,    62,    63,    64,    65,     0,    60,
      61,    62,    63,    64,    65,     0,     0,     0,    66,     0,
       0,     0,     0,     0,   156,    66,     0,     0,     0,     0,
       0,   157,    60,    61,    62,    63,    64,    65,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    66,     0,
       0,     0,     0,     0,   174,    68,    69,    70,    71,    72,
      73,    74,    75,    76,    77,    78,    79,    80,     0,     0,
      81,     0,     0,     0,     0,     0,    82,    60,    61,    62,
      63,    64,    65,    60,    61,    62,    63,    64,    65,     0,
       0,     0,     0,    66,     0,     0,     0,     0,   151,    66,
       0,     0,     0,     0,   165,    60,    61,    62,    63,    64,
      65,   142,     0,     0,    60,    61,    62,    63,    64,    65,
       0,    66,     0,     0,     0,     0,   178,     0,     0,     0,
      66,   143,    68,    69,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    80,     0,     0,    81,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      80,     0,     0,    81,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,     0,     0,    81,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
       0,     0,    81
};

static const yytype_int16 yycheck[] =
{
       5,    20,    39,    10,     0,     1,     3,    40,     5,     6,
      15,    16,    17,    40,    36,    20,    21,    33,    19,    -1,
      25,    26,    27,    28,    29,    30,    -1,    -1,    33,    34,
      35,    36,    37,    38,    -1,    -1,    41,    42,    43,    58,
      59,    37,    38,    20,    21,    22,    23,    24,    25,    28,
      29,    30,    57,    16,    33,    18,    19,     5,     3,    36,
       5,     6,    39,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    25,    41,    10,
      11,    12,    13,    14,    15,    33,    10,    11,    12,    13,
      14,    15,    -1,    41,    42,    43,     3,     4,     5,     6,
       7,     8,     9,    34,    35,    16,    -1,    18,    19,    57,
      34,    35,    60,    61,    62,    63,    64,    65,    66,    26,
      27,    -1,    -1,    -1,    31,    -1,    -1,    34,    35,    -1,
     135,    -1,   137,    40,    41,    -1,    -1,    44,   143,    -1,
      -1,    -1,     3,     4,     5,     6,     7,     8,     9,   154,
     155,    -1,    20,    21,    22,    23,    24,    25,   163,   164,
      -1,    -1,    -1,    -1,    -1,    26,    27,   172,    36,    -1,
      31,    -1,   177,    34,    35,    -1,    -1,    -1,    -1,    40,
      41,    -1,    -1,    44,    -1,    -1,    -1,   135,    -1,    -1,
      -1,   139,    -1,   141,   142,    -1,    -1,     3,     4,     5,
       6,     7,     8,     9,    -1,    -1,    -1,    -1,   156,   157,
      -1,    20,    21,    22,    23,    24,    25,    -1,    -1,   167,
      26,    27,    -1,    -1,    -1,    31,   174,    36,    34,    35,
      -1,    -1,    41,    42,    40,    -1,    -1,    -1,    44,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    -1,    -1,    33,    26,    27,    28,    29,    30,
      -1,    -1,    33,    -1,    -1,    -1,    45,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      -1,    -1,    33,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      41,    42,    43,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    -1,    -1,    33,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    41,    42,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    -1,    -1,    33,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    41,    42,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    -1,    -1,    33,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    41,    42,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    -1,    -1,    33,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    41,    42,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    -1,    -1,    33,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    42,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    -1,    -1,    33,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    42,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    -1,    -1,    33,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    41,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      -1,    -1,    33,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      41,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    -1,    -1,    33,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    41,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    -1,    -1,
      33,    20,    21,    22,    23,    24,    25,    -1,    41,    20,
      21,    22,    23,    24,    25,    -1,    -1,    36,    -1,    -1,
      -1,    -1,    41,    42,    -1,    36,    -1,    -1,    -1,    -1,
      41,    42,    20,    21,    22,    23,    24,    25,    -1,    20,
      21,    22,    23,    24,    25,    -1,    -1,    -1,    36,    -1,
      -1,    -1,    -1,    -1,    42,    36,    -1,    -1,    -1,    -1,
      -1,    42,    20,    21,    22,    23,    24,    25,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    36,    -1,
      -1,    -1,    -1,    -1,    42,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    -1,    -1,
      33,    -1,    -1,    -1,    -1,    -1,    39,    20,    21,    22,
      23,    24,    25,    20,    21,    22,    23,    24,    25,    -1,
      -1,    -1,    -1,    36,    -1,    -1,    -1,    -1,    41,    36,
      -1,    -1,    -1,    -1,    41,    20,    21,    22,    23,    24,
      25,    17,    -1,    -1,    20,    21,    22,    23,    24,    25,
      -1,    36,    -1,    -1,    -1,    -1,    41,    -1,    -1,    -1,
      36,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    -1,    -1,    33,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    -1,    -1,    33,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    -1,    -1,    33,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      -1,    -1,    33
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    47,     0,     1,    37,    38,    48,    39,     3,     4,
       5,     6,     7,     8,     9,    26,    27,    31,    34,    35,
      40,    44,    49,    50,    51,    10,    11,    12,    13,    14,
      15,    34,    35,    10,    11,    12,    13,    14,    15,    34,
      35,    10,    40,    40,    50,    51,    51,    51,     3,     5,
       6,     3,     5,     6,    49,    51,    51,    16,    18,    19,
      20,    21,    22,    23,    24,    25,    36,    39,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    33,    39,    50,    51,    51,    51,    51,    51,    51,
      50,    51,    51,    51,    51,    51,    51,    50,    51,    41,
      50,    51,    41,    50,    51,    41,    41,    45,    50,    51,
      49,    51,    49,    50,    50,    50,    50,    50,    50,    50,
      51,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    51,    51,    41,    42,    41,    42,    41,    42,
      41,    42,    17,    17,    50,    51,    51,    50,    50,    50,
      51,    41,    41,    41,    42,    43,    42,    42,    51,    51,
      50,    50,    41,    42,    42,    41,    41,    42,    51,    51,
      50,    41,    42,    41,    42,    51,    50,    42,    41,    51,
      41
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

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
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
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
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

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
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

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
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


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

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

       Refer to the stacks thru separate pointers, to allow yyoverflow
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
  int yytoken;
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

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

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
  if (yyn == YYPACT_NINF)
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
      if (yyn == 0 || yyn == YYTABLE_NINF)
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
  *++yyvsp = yylval;

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

/* Line 1464 of yacc.c  */
#line 89 "aprepro.y"
    { if (echo) fprintf(yyout,"\n");	}
    break;

  case 5:

/* Line 1464 of yacc.c  */
#line 90 "aprepro.y"
    { if (echo) {
	                             format = getsym("_FORMAT");
	                             fprintf(yyout, format->value.svar, (yyvsp[(2) - (3)].val));
				   }                                    }
    break;

  case 6:

/* Line 1464 of yacc.c  */
#line 94 "aprepro.y"
    { if (echo && (yyvsp[(2) - (3)].string) != NULL)
				    fprintf(yyout, "%s", (yyvsp[(2) - (3)].string));	}
    break;

  case 7:

/* Line 1464 of yacc.c  */
#line 96 "aprepro.y"
    { yyerrok;				}
    break;

  case 8:

/* Line 1464 of yacc.c  */
#line 99 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) < (yyvsp[(3) - (3)].val);                         }
    break;

  case 9:

/* Line 1464 of yacc.c  */
#line 100 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) > (yyvsp[(3) - (3)].val);                         }
    break;

  case 10:

/* Line 1464 of yacc.c  */
#line 101 "aprepro.y"
    { (yyval.val) = !((yyvsp[(2) - (2)].val));                           }
    break;

  case 11:

/* Line 1464 of yacc.c  */
#line 102 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) <= (yyvsp[(3) - (3)].val);                        }
    break;

  case 12:

/* Line 1464 of yacc.c  */
#line 103 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) >= (yyvsp[(3) - (3)].val);                        }
    break;

  case 13:

/* Line 1464 of yacc.c  */
#line 104 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) == (yyvsp[(3) - (3)].val);                        }
    break;

  case 14:

/* Line 1464 of yacc.c  */
#line 105 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) != (yyvsp[(3) - (3)].val);                        }
    break;

  case 15:

/* Line 1464 of yacc.c  */
#line 106 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) || (yyvsp[(3) - (3)].val);                        }
    break;

  case 16:

/* Line 1464 of yacc.c  */
#line 107 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) && (yyvsp[(3) - (3)].val);                        }
    break;

  case 17:

/* Line 1464 of yacc.c  */
#line 108 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) || (yyvsp[(3) - (3)].val);                        }
    break;

  case 18:

/* Line 1464 of yacc.c  */
#line 109 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) && (yyvsp[(3) - (3)].val);                        }
    break;

  case 19:

/* Line 1464 of yacc.c  */
#line 110 "aprepro.y"
    { (yyval.val) = (yyvsp[(2) - (3)].val);                              }
    break;

  case 20:

/* Line 1464 of yacc.c  */
#line 113 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) <  0 ? 1 : 0);	}
    break;

  case 21:

/* Line 1464 of yacc.c  */
#line 114 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) >  0 ? 1 : 0);	}
    break;

  case 22:

/* Line 1464 of yacc.c  */
#line 115 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) <= 0 ? 1 : 0);	}
    break;

  case 23:

/* Line 1464 of yacc.c  */
#line 116 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) >= 0 ? 1 : 0);	}
    break;

  case 24:

/* Line 1464 of yacc.c  */
#line 117 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) == 0 ? 1 : 0);	}
    break;

  case 25:

/* Line 1464 of yacc.c  */
#line 118 "aprepro.y"
    { (yyval.val) = (strcmp((yyvsp[(1) - (3)].string),(yyvsp[(3) - (3)].string)) != 0 ? 1 : 0);	}
    break;

  case 26:

/* Line 1464 of yacc.c  */
#line 120 "aprepro.y"
    { (yyval.string) = (yyvsp[(1) - (1)].string);				}
    break;

  case 27:

/* Line 1464 of yacc.c  */
#line 121 "aprepro.y"
    { (yyval.string) = (yyvsp[(1) - (1)].tptr)->value.svar;			}
    break;

  case 28:

/* Line 1464 of yacc.c  */
#line 122 "aprepro.y"
    { (yyval.string) = (yyvsp[(3) - (3)].string); (yyvsp[(1) - (3)].tptr)->value.svar = (yyvsp[(3) - (3)].string);
				  (yyvsp[(1) - (3)].tptr)->type = SVAR;			}
    break;

  case 29:

/* Line 1464 of yacc.c  */
#line 124 "aprepro.y"
    { (yyval.string) = (yyvsp[(3) - (3)].string); 
				  (yyvsp[(1) - (3)].tptr)->value.svar = (yyvsp[(3) - (3)].string);
				  redefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 30:

/* Line 1464 of yacc.c  */
#line 127 "aprepro.y"
    { (yyval.string) = (yyvsp[(3) - (3)].string); 
				  (yyvsp[(1) - (3)].tptr)->value.svar= (yyvsp[(3) - (3)].string);
				  redefined_warning((yyvsp[(1) - (3)].tptr)->name);          
				  (yyvsp[(1) - (3)].tptr)->type = SVAR; 		}
    break;

  case 31:

/* Line 1464 of yacc.c  */
#line 131 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (4)].tptr)->value.strfnct))((yyvsp[(3) - (4)].string));	}
    break;

  case 32:

/* Line 1464 of yacc.c  */
#line 132 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (3)].tptr)->value.strfnct))();	}
    break;

  case 33:

/* Line 1464 of yacc.c  */
#line 133 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (4)].tptr)->value.strfnct))((yyvsp[(3) - (4)].val));	}
    break;

  case 34:

/* Line 1464 of yacc.c  */
#line 134 "aprepro.y"
    { int len1 = strlen((yyvsp[(1) - (3)].string));
				  int len3 = strlen((yyvsp[(3) - (3)].string));
				  ALLOC((yyval.string), len1+len3+1, char *);
				  memcpy((yyval.string), (yyvsp[(1) - (3)].string), len1+1);
				  strcat((yyval.string), (yyvsp[(3) - (3)].string)); }
    break;

  case 35:

/* Line 1464 of yacc.c  */
#line 140 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (12)].tptr)->value.strfnct))((yyvsp[(3) - (12)].val), (yyvsp[(5) - (12)].string), (yyvsp[(7) - (12)].string), (yyvsp[(9) - (12)].string), (yyvsp[(11) - (12)].string)); }
    break;

  case 36:

/* Line 1464 of yacc.c  */
#line 142 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (8)].tptr)->value.strfnct))((yyvsp[(3) - (8)].val), (yyvsp[(5) - (8)].string), (yyvsp[(7) - (8)].string)); }
    break;

  case 37:

/* Line 1464 of yacc.c  */
#line 144 "aprepro.y"
    { (yyval.string) = (*((yyvsp[(1) - (8)].tptr)->value.strfnct))((yyvsp[(3) - (8)].string), (yyvsp[(5) - (8)].string), (yyvsp[(7) - (8)].string)); }
    break;

  case 38:

/* Line 1464 of yacc.c  */
#line 145 "aprepro.y"
    { (yyval.string) = ((yyvsp[(1) - (5)].val)) ? ((yyvsp[(3) - (5)].string)) : ((yyvsp[(5) - (5)].string));              }
    break;

  case 39:

/* Line 1464 of yacc.c  */
#line 147 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (1)].val); 				}
    break;

  case 40:

/* Line 1464 of yacc.c  */
#line 148 "aprepro.y"
    { (yyval.val) = (yyvsp[(2) - (2)].val) + 1;				}
    break;

  case 41:

/* Line 1464 of yacc.c  */
#line 149 "aprepro.y"
    { (yyval.val) = (yyvsp[(2) - (2)].val) - 1;				}
    break;

  case 42:

/* Line 1464 of yacc.c  */
#line 150 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (1)].tptr)->value.var;			}
    break;

  case 43:

/* Line 1464 of yacc.c  */
#line 151 "aprepro.y"
    { (yyval.val) = ++((yyvsp[(2) - (2)].tptr)->value.var);		}
    break;

  case 44:

/* Line 1464 of yacc.c  */
#line 152 "aprepro.y"
    { (yyval.val) = --((yyvsp[(2) - (2)].tptr)->value.var);		}
    break;

  case 45:

/* Line 1464 of yacc.c  */
#line 153 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (2)].tptr)->value.var)++;		}
    break;

  case 46:

/* Line 1464 of yacc.c  */
#line 154 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (2)].tptr)->value.var)--;		}
    break;

  case 47:

/* Line 1464 of yacc.c  */
#line 155 "aprepro.y"
    { (yyval.val) = (yyvsp[(3) - (3)].val); (yyvsp[(1) - (3)].tptr)->value.var = (yyvsp[(3) - (3)].val);
				  redefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 48:

/* Line 1464 of yacc.c  */
#line 157 "aprepro.y"
    { (yyval.val) = (yyvsp[(3) - (3)].val); (yyvsp[(1) - (3)].tptr)->value.var = (yyvsp[(3) - (3)].val);
				  redefined_warning((yyvsp[(1) - (3)].tptr)->name);          
				  (yyvsp[(1) - (3)].tptr)->type = VAR;			}
    break;

  case 49:

/* Line 1464 of yacc.c  */
#line 160 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var += (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; }
    break;

  case 50:

/* Line 1464 of yacc.c  */
#line 161 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var -= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; }
    break;

  case 51:

/* Line 1464 of yacc.c  */
#line 162 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var *= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; }
    break;

  case 52:

/* Line 1464 of yacc.c  */
#line 163 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var /= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; }
    break;

  case 53:

/* Line 1464 of yacc.c  */
#line 164 "aprepro.y"
    { errno = 0;
				  (yyvsp[(1) - (3)].tptr)->value.var = pow((yyvsp[(1) - (3)].tptr)->value.var,(yyvsp[(3) - (3)].val)); 
				  (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
				  MATH_ERROR("Power");
				}
    break;

  case 54:

/* Line 1464 of yacc.c  */
#line 169 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (1)].tptr)->value.var;
				  undefined_warning((yyvsp[(1) - (1)].tptr)->name);          }
    break;

  case 55:

/* Line 1464 of yacc.c  */
#line 171 "aprepro.y"
    { (yyval.val) = ++((yyvsp[(2) - (2)].tptr)->value.var);		
				  (yyvsp[(2) - (2)].tptr)->type = VAR;                       
				  undefined_warning((yyvsp[(2) - (2)].tptr)->name);          }
    break;

  case 56:

/* Line 1464 of yacc.c  */
#line 174 "aprepro.y"
    { (yyval.val) = --((yyvsp[(2) - (2)].tptr)->value.var);		
				  (yyvsp[(2) - (2)].tptr)->type = VAR;                       
				  undefined_warning((yyvsp[(2) - (2)].tptr)->name);          }
    break;

  case 57:

/* Line 1464 of yacc.c  */
#line 177 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (2)].tptr)->value.var)++;		
				  (yyvsp[(1) - (2)].tptr)->type = VAR;                       
				  undefined_warning((yyvsp[(1) - (2)].tptr)->name);          }
    break;

  case 58:

/* Line 1464 of yacc.c  */
#line 180 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (2)].tptr)->value.var)--;		
				  (yyvsp[(1) - (2)].tptr)->type = VAR;                       
				  undefined_warning((yyvsp[(1) - (2)].tptr)->name);          }
    break;

  case 59:

/* Line 1464 of yacc.c  */
#line 183 "aprepro.y"
    { (yyval.val) = (yyvsp[(3) - (3)].val); (yyvsp[(1) - (3)].tptr)->value.var = (yyvsp[(3) - (3)].val);
				  (yyvsp[(1) - (3)].tptr)->type = VAR;                       }
    break;

  case 60:

/* Line 1464 of yacc.c  */
#line 185 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var += (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
				  (yyvsp[(1) - (3)].tptr)->type = VAR;                       
				  undefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 61:

/* Line 1464 of yacc.c  */
#line 188 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var -= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
				  (yyvsp[(1) - (3)].tptr)->type = VAR;                       
				  undefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 62:

/* Line 1464 of yacc.c  */
#line 191 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var *= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
				  (yyvsp[(1) - (3)].tptr)->type = VAR;                       
				  undefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 63:

/* Line 1464 of yacc.c  */
#line 194 "aprepro.y"
    { (yyvsp[(1) - (3)].tptr)->value.var /= (yyvsp[(3) - (3)].val); (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
				  (yyvsp[(1) - (3)].tptr)->type = VAR;                       
				  undefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 64:

/* Line 1464 of yacc.c  */
#line 197 "aprepro.y"
    { errno = 0;
				  (yyvsp[(1) - (3)].tptr)->value.var = pow((yyvsp[(1) - (3)].tptr)->value.var,(yyvsp[(3) - (3)].val)); 
				  (yyval.val) = (yyvsp[(1) - (3)].tptr)->value.var; 
				  (yyvsp[(1) - (3)].tptr)->type = VAR;                       
				  MATH_ERROR("Power");
				  undefined_warning((yyvsp[(1) - (3)].tptr)->name);          }
    break;

  case 65:

/* Line 1464 of yacc.c  */
#line 203 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (3)].tptr)->value.fnctptr))();	}
    break;

  case 66:

/* Line 1464 of yacc.c  */
#line 204 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (4)].tptr)->value.fnctptr))((yyvsp[(3) - (4)].val)); 	}
    break;

  case 67:

/* Line 1464 of yacc.c  */
#line 205 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (4)].tptr)->value.fnctptr))((yyvsp[(3) - (4)].string)); 	}
    break;

  case 68:

/* Line 1464 of yacc.c  */
#line 207 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (6)].tptr)->value.fnctptr))((yyvsp[(3) - (6)].string), (yyvsp[(5) - (6)].val)); 	}
    break;

  case 69:

/* Line 1464 of yacc.c  */
#line 209 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (6)].tptr)->value.fnctptr))((yyvsp[(3) - (6)].string), (yyvsp[(5) - (6)].string)); 	}
    break;

  case 70:

/* Line 1464 of yacc.c  */
#line 211 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (6)].tptr)->value.fnctptr))((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));	}
    break;

  case 71:

/* Line 1464 of yacc.c  */
#line 213 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (8)].tptr)->value.fnctptr))((yyvsp[(3) - (8)].val), (yyvsp[(5) - (8)].val), (yyvsp[(7) - (8)].val)); }
    break;

  case 72:

/* Line 1464 of yacc.c  */
#line 215 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (10)].tptr)->value.fnctptr))((yyvsp[(3) - (10)].val), (yyvsp[(5) - (10)].val), (yyvsp[(7) - (10)].val), (yyvsp[(9) - (10)].val)); }
    break;

  case 73:

/* Line 1464 of yacc.c  */
#line 217 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (10)].tptr)->value.fnctptr))((yyvsp[(3) - (10)].val), (yyvsp[(5) - (10)].val), (yyvsp[(7) - (10)].val), (yyvsp[(9) - (10)].val)); }
    break;

  case 74:

/* Line 1464 of yacc.c  */
#line 219 "aprepro.y"
    { (yyval.val) = (*((yyvsp[(1) - (14)].tptr)->value.fnctptr))((yyvsp[(3) - (14)].val), (yyvsp[(5) - (14)].val), (yyvsp[(7) - (14)].val), (yyvsp[(9) - (14)].val), (yyvsp[(11) - (14)].val), (yyvsp[(13) - (14)].val)); }
    break;

  case 75:

/* Line 1464 of yacc.c  */
#line 220 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) + (yyvsp[(3) - (3)].val); 			}
    break;

  case 76:

/* Line 1464 of yacc.c  */
#line 221 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) - (yyvsp[(3) - (3)].val); 			}
    break;

  case 77:

/* Line 1464 of yacc.c  */
#line 222 "aprepro.y"
    { (yyval.val) = (yyvsp[(1) - (3)].val) * (yyvsp[(3) - (3)].val); 			}
    break;

  case 78:

/* Line 1464 of yacc.c  */
#line 223 "aprepro.y"
    { if ((yyvsp[(3) - (3)].val) == 0.)
				    {
				      yyerror("Zero divisor"); 
				      yyerrok;
				    }
				  else
				    (yyval.val) = (yyvsp[(1) - (3)].val) / (yyvsp[(3) - (3)].val); 			}
    break;

  case 79:

/* Line 1464 of yacc.c  */
#line 230 "aprepro.y"
    { if ((yyvsp[(3) - (3)].val) == 0.)
				    {
				      yyerror("Zero divisor");
				      yyerrok;
				    }
				  else
				    (yyval.val) = (int)(yyvsp[(1) - (3)].val) % (int)(yyvsp[(3) - (3)].val);		}
    break;

  case 80:

/* Line 1464 of yacc.c  */
#line 237 "aprepro.y"
    { (yyval.val) = -(yyvsp[(2) - (2)].val);				}
    break;

  case 81:

/* Line 1464 of yacc.c  */
#line 238 "aprepro.y"
    { (yyval.val) =  (yyvsp[(2) - (2)].val);				}
    break;

  case 82:

/* Line 1464 of yacc.c  */
#line 239 "aprepro.y"
    { errno = 0;
				  (yyval.val) = pow((yyvsp[(1) - (3)].val), (yyvsp[(3) - (3)].val)); 
				  MATH_ERROR("Power");			}
    break;

  case 83:

/* Line 1464 of yacc.c  */
#line 242 "aprepro.y"
    { (yyval.val) = (yyvsp[(2) - (3)].val);				}
    break;

  case 84:

/* Line 1464 of yacc.c  */
#line 243 "aprepro.y"
    { errno = 0;
				  (yyval.val) = (double)((yyvsp[(2) - (3)].val) < 0 ? 
					-floor(-((yyvsp[(2) - (3)].val))): floor((yyvsp[(2) - (3)].val)) );
				  MATH_ERROR("floor (int)");		}
    break;

  case 85:

/* Line 1464 of yacc.c  */
#line 248 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (1)].val)) ? 1 : 0; }
    break;

  case 86:

/* Line 1464 of yacc.c  */
#line 249 "aprepro.y"
    { (yyval.val) = ((yyvsp[(1) - (5)].val)) ? ((yyvsp[(3) - (5)].val)) : ((yyvsp[(5) - (5)].val));              }
    break;



/* Line 1464 of yacc.c  */
#line 2343 "y.tab.c"
      default: break;
    }
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
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
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
      if (yyn != YYPACT_NINF)
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

  *++yyvsp = yylval;


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

#if !defined(yyoverflow) || YYERROR_VERBOSE
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
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
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



/* Line 1684 of yacc.c  */
#line 253 "aprepro.y"

# include "lex.yy.c"



