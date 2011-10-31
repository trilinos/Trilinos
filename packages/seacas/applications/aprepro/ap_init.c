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
/***
   NAME
     init
   PURPOSE
     Initialize variables and functions Aprepro
***/
#include <stdlib.h> /* for malloc */
#include "my_aprepro.h"
#include "y.tab.h"
#include <sys/types.h>
#include "init_structs.h"

void init_table(char comment);
char comm_string[2];

extern double
  do_fabs(double x), do_acos(double x), do_acosd(double x), do_asin(double x), do_asind(double x),
  do_atan(double x), do_atan2(double x, double y), do_atan2d(double x, double y),
  do_atand(double x), do_ceil(double x), do_cos(double x), do_cosd(double x), do_cosh(double x),
  do_d2r(double x), do_dim(double x, double y), do_dist(double x1, double y1, double x2, double y2),
  do_exp(double x), do_floor(double x), do_fmod(double x, double y), do_int(double x),
  do_log(double x), do_log10(double x), do_max(double x, double y), do_min(double x, double y),
  do_r2d(double x), do_rand(double xl, double xh), do_sign(double x, double y), do_sin(double x),
  do_rand_normal(double mean, double stddev), do_rand_weibull(double alpha, double beta),
  do_rand_lognormal(double mean, double stddev), do_srand(double seed),
  do_sind(double x), do_sinh(double x), do_sqrt(double x), do_tan(double x), do_tand(double x),
  do_tanh(double x), do_hypot(double x, double y), do_polarX(double rad, double ang),
  do_polarY(double rad, double ang), do_angle(double x1, double y1, double x2, double y2),
  do_angled(double x1, double y1, double x2, double y2), do_lgamma(double val),
  do_julday(double mon, double day, double year),
  do_juldayhms(double mon, double day, double year, double h, double mi, double se),
  do_log1p(double x), do_acosh(double x), do_asinh(double x), do_atanh(double x),
  do_word_count(char *string, char *delm), do_strtod(char *string),
  do_nint(double x), do_option(char *option, double value);

struct init arith_fncts[] =
{
  {"Vangle",         do_angle,  "Vangle(x1,y1,x2,y2)","Angle (radians) between vector x1_i+y1_j and x2_i+y2_j."},
  {"Vangled",        do_angled, "Vangled(x1,y1,x2,y2)","Angle (degrees) between vector x1_i+y1_j and x2_i+y2_j."},
  {"abs",            do_fabs,   "abs(x)",   "Absolute value of x. |x|."},
  {"acos",           do_acos,   "acos(x)",  "Inverse cosine of x, returns radians."},
  {"acosd",          do_acosd,  "acosd(x)", "Inverse cosine of x, returns degrees."},
  {"acosh",          do_acosh,  "acosh(x)", "Inverse hyperbolic cosine of x."},
  {"asin",           do_asin,   "asin(x)",  "Inverse sine of x, returns radians."},
  {"asind",          do_asind,  "asin(x)",  "Inverse sine of x, returns degrees."},
  {"asinh",          do_asinh,  "asinh(x)", "Inverse hyperbolic sine of x."},
  {"atan",           do_atan,   "atan(x)",  "Inverse tangent of x, returns radians."},
  {"atan2",          do_atan2,  "atan2(x,y)","Inverse tangent of y/x, returns radians."},
  {"atan2d",         do_atan2d, "atan2d(x,y)","Inverse tangent of y/x, returns degrees."},
  {"atand",          do_atand,  "atand(x)",  "Inverse tangent of x, returns degrees."},
  {"atanh",          do_atanh,  "atanh(x)",  "Inverse hyperbolic tangent of x."},
  {"ceil",           do_ceil,   "ceil(x)",   "Smallest integer not less than x."},
  {"cos",            do_cos,    "cos(x)",    "Cosine of x, with x in radians"},
  {"cosd",           do_cosd,   "cosd(x)",   "Cosine of x, with x in degrees"},
  {"cosh",           do_cosh,   "cosh(x)",   "Hyperbolic cosine of x."},
  {"d2r",            do_d2r,    "d2r(x)",    "Degrees to radians."},
  {"dim",            do_dim,    "dim(x,y)",  "x - min(x,y)"},
  {"dist",           do_dist,   "dist(x1,y1, x2,y2)", "sqrt((x1-x2)^2 + (y1-y2)^2)"},
  {"exp",            do_exp,    "exp(x)",    "Exponential: e^x"},
  {"floor",          do_floor,  "floor(x)",  "Largest integer not greater than x."},
  {"fmod",           do_fmod,   "fmod(x,y)", "Floating-point remainder of x/y."},
  {"hypot",          do_hypot,  "hypot(x,y)", "sqrt(x^2+y^2)."},
  {"int",            do_int,    "int(x), [x]","Integer part of x truncated toward 0."},
  {"julday",         do_julday, "julday(mm, dd, yy)","Julian day corresponding to mm/dd/yy. "},
  {"juldayhms",      do_juldayhms,"juldayhms(mm, dd, yy, hh, mm, ss)","Julian day corresponding to mm/dd/yy at hh:mm:ss "},
  {"lgamma",         do_lgamma, "lgamma(x)", "log(Gamma(x))."},
  {"ln",             do_log,    "ln(x)",     "Natural (base e) logarithm of x."},
  {"log",            do_log,    "log(x)",    "Natural (base e) logarithm of x."},
  {"log10",          do_log10,  "log10(x)",  "Base 10 logarithm of x. "},
  {"log1p",          do_log1p,  "log1p(x)",  "log(1+x) "},
  {"max",            do_max,    "max(x,y)",  "Maximum of x and y. "},
  {"min",            do_min,    "min(x,y)",  "Minimum of x and y. "},
  {"nint",           do_nint,   "nint(x)",   "Rounds x to nearest integer. <0.5 down; >=0.5 up."},
  {"polarX",         do_polarX, "polarX(r,a)","r * cos(a), a is in degrees "},
  {"polarY",         do_polarY, "polarY(r,a)","r * sin(a), a is in degrees "},
  {"r2d",            do_r2d,    "r2d(x)",     "Radians to degrees. "},
  {"srand",          do_srand,  "srand(seed)","Seed the random number generator with the given integer value. "},
  {"rand",           do_rand,   "rand(xl,xh)","Random value between xl and xh; uniformly distributed. "},
  {"rand_normal",    do_rand_normal, "rand_normal(m,s)","Random value normally distributed with mean m and stddev s."},
  {"rand_lognormal", do_rand_lognormal,"rand_lognormal(m,s)","Random value with lognormal distribution with mean m and stddev s."},
  {"rand_weibull",   do_rand_weibull,"rand_weibull(a, b)","Random value with weibull distribution with alpha=a and beta=b. "},
  {"sign",           do_sign,   "sign(x,y)","x * sgn(y)"},
  {"sin",            do_sin,    "sin(x)","Sine of x, with x in radians. "},
  {"sind",           do_sind,   "sind(x)","Sine of x, with x in degrees. "},
  {"sinh",           do_sinh,   "sinh(x)","Hyperbolic sine of x "},
  {"sqrt",           do_sqrt,   "sqrt(x)","Square root of x. "},
  {"tan",            do_tan,    "tan(x)","Tangent of x, with x in radians. "},
  {"tand",           do_tand,   "tand(x)","Tangent of x, with x in radians. "},
  {"tanh",           do_tanh,   "tanh(x)","Hyperbolic tangent of x. "},
  {"word_count",     do_word_count,"word_count(svar,del)","Number of words in svar. Words are separated by one or more of the characters in the string variable del."},
  {"strtod",         do_strtod, "strtod(svar)","Returns a double-precision floating-point number equal to the value represented by the character string pointed to by svar."},
  {"option",         do_option, "option(?,?)","Internal"},
  {0, 0, 0, 0}				/* Last line must be 0, 0 */
};

extern char *do_tolower(char *string), *do_toupper(char *string), *do_tostring(double x),
  *do_output(char *filename), *do_get_word(double n, char *string, char *delm),
  *do_execute(char *string), *do_getenv(char *env), *do_error(char *error_string),
  *do_rescan(char *string),  *do_Units(char *type), *do_dumpsym(void), *do_dumpfunc(void), *do_help(void),
  *do_file_to_string(char *filename), *do_append(char *filename), *do_dumpvar(void),
#if !defined(NO_EXODUSII)
  *do_exodus_info(char *filename), *do_exodus_meta(char *filename),
#endif
  *do_include_path(char *new_path), *do_intout(double intval), *do_get_date(void), *do_get_time(void),
  *do_extract(char *string, char *begin, char *end);
  

struct str_init string_fncts[] =
{
  {"tolower",        do_tolower,     "tolower(svar)","Translates all uppercase characters in svar to lowercase. It modifies svar and returns the resulting string.  "},
  {"toupper",        do_toupper,     "toupper(svar)","Translates all lowercase character in svar to uppercase. It modifies svar and returns the resulting string. "},
  {"tostring",       do_tostring,    "tostring(x)","Returns a string representation of the numerical varaible x. The variable x is unchanged. "},
  {"to_lower",       do_tolower,     "to_lower(svar)","Translates all uppercase characters in svar to lowercase. It modifies svar and returns the resulting string.  "},
  {"to_upper",       do_toupper,     "toupper(svar)","Translates all lowercase character in svar to uppercase. It modifies svar and returns the resulting string. "},
  {"to_string",      do_tostring,    "tostring(x)","Returns a string representation of the numerical varaible x. The variable x is unchanged. "},
  {"getenv",         do_getenv,      "getenv(svar)","Returns a string containing the value of the environment variable svar. If the environment variable is not defined, an empty string is returned. "},
  {"error",          do_error,       "error(svar)","Outputs the string svar to stderr and then terminates the code with an error exit status."},
  {"output",         do_output,      "output(filename)","Creates the file specified by filename and sends all subsequent output from aprepro to that file."},
  {"output_append",  do_append,      "output_append(fn)","If file with name fn exists, append output to it; otherwise create the file and send all subsequent output from aprepro to that file."},
  {"get_word",       do_get_word,    "get_word(n,svar,del)","Returns a string containing the nth word of svar. The words are separated by one or more of the characters in the string variable del "},
  {"file_to_string", do_file_to_string,    "file_to_string(fn)","Opens the file specified by fn and returns the contents as a multi-line string."},
  {"execute",        do_execute,     "execute(svar)","svar is parsed and executed as if it were a line read from the input file."},
  {"rescan",         do_rescan,      "rescan(svar)","The difference between execute(sv1) and rescan(sv2) is that sv1 must be a valid expression, but sv2 can contain zero or more expressions. "},
  {"Units",          do_Units,       "Units(svar)","See manual. svar is one of the defined units systems:\n\t\t\t'si', 'cgs', 'cgs-ev', 'shock', 'swap', 'ft-lbf-s', 'ft-lbm-s', 'in-lbf-s'"},
  {"help",           do_help,        "help()","Tell how to get help on variables, functions, ..."},
  {"DUMP",           do_dumpsym,     "DUMP()","Output a list of all defined variables and their value."},
  {"DUMP_FUNC",      do_dumpfunc,    "DUMP_FUNC()","Output a list of all double and string functions recognized by aprepro."},
  {"DUMP_PREVAR",    do_dumpvar,    "DUMP_PREVAR()","Output a list of all predefined variables and their value."},
  {"include_path",   do_include_path,"include_path(path)","Specify an optional path to be prepended to a filename when opening a file. Can also be specified via the -I command line option when executing aprepro."},
  {"IO",             do_intout,      "IO(x)","Convert x to an integer and then to a string. "},
  {"get_date",       do_get_date,    "get_date()","Returns a string representing the current date in the form YYYY/MM/DD."},
  {"get_time",       do_get_time,    "get_time()","Returns a string representing the current time in the form HH:MM:SS."},
  {"extract",        do_extract,     "extract(s, b, e)","Return substring [b,e). 'b' is included; 'e' is not. If 'b' not found, return empty; If 'e' not found, return rest of string. If 'b' empty, start at beginning; if 'e' empty, return rest of string."},
#if !defined(NO_EXODUSII)
  {"exodus_info",    do_exodus_info, "exodus_info(ex_fn)","Parses the info records extracted from the exodus file 'ex_fn'"},
  {"exodus_meta",    do_exodus_meta, "exodus_meta(ex_fn)","Creates several variables related to the exodusII metadata in the specified file. Experimental."},
#endif
  {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
};

struct var_init variables[] =
{
  {"DEG",  57.29577951308232087680},	/* 180/pi, degrees per radian */
  {"RAD",   0.01745329251994329576},	/* pi/180, radians per degree */
  {"E",     2.71828182845904523536},	/* e, base of natural log     */
  {"GAMMA", 0.57721566490153286060},	/* euler-mascheroni constant  */
  {"PHI",   1.61803398874989484820},	/* golden ratio               */
  {"PI",    3.14159265358979323846},	/* pi                         */
  {"PI_2",  1.57079632679489661923},	/* pi / 2			 */
  {"SQRT2", 1.41421356237309504880},	/* square root of 2		 */
  {"TRUE",  1},
  {"FALSE", 0},
  {0, 0}				/* Last line must be 0, 0 */
};

struct svar_init svariables[] =
{
  {"_FORMAT", "%.10g"},	/* Default output format */
  {0, 0}		/* Last line must be 0, 0 */
};
/* NOTE: The current comment is stored in "_C_"
 *	 Since it can be changed by user on command line, we
 *	 initialize is differently than the other string variables.
 */

void init_table(char comment)
{
  int i;
  symrec *ptr;
  for (i = 0; arith_fncts[i].fname != 0; i++)
    {
      ptr = putsym(arith_fncts[i].fname, FNCT, 1);
      ptr->value.fnctptr = arith_fncts[i].fnct;
      ptr->info = arith_fncts[i].description;
      ptr->syntax = arith_fncts[i].syntax;
    }
  for (i = 0; string_fncts[i].fname != 0; i++)
    {
      ptr = putsym(string_fncts[i].fname, SFNCT, 1);
      ptr->value.strfnct = string_fncts[i].fnct;
      ptr->info = string_fncts[i].description;
      ptr->syntax = string_fncts[i].syntax;
    }
  for (i = 0; variables[i].vname != 0; i++)
    {
      ptr = putsym(variables[i].vname, VAR, 1);
      ptr->value.var = variables[i].value;
    }
  for (i = 0; svariables[i].vname != 0; i++)
    {
      ptr = putsym(svariables[i].vname, SVAR, 1);
      ptr->value.svar = svariables[i].value;
    }
  sprintf(comm_string, "%c", comment);
  ptr = putsym("_C_", SVAR, 1);
  ptr->value.svar = comm_string;
  {
    char *version_string = (char *)malloc(8);
    version(version_string);
    ptr = putsym("VERSION", SVAR, 1);
    ptr->value.svar = version_string;
  }
}

