/* A Bison parser, made by GNU Bison 2.4.2.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2006, 2009-2010 Free Software
   Foundation, Inc.
   
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

/* Line 1749 of yacc.c  */
#line 52 "aprepro.y"

  double  val;		/* For returning numbers.		*/
  symrec *tptr;		/* For returning symbol-table pointers	*/
  char   *string;	/* For returning quoted strings		*/



/* Line 1749 of yacc.c  */
#line 115 "y.tab.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;


