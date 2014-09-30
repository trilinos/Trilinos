// Copyright (c) 2014, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

/***
   NAME
     init
   PURPOSE
     Initialize variables and functions Aprepro
***/
#include "apr_builtin.h"
#include <cstring>                      // for strncpy
#include <string>                       // for string
#include "aprepro.h"                    // for symrec, Aprepro, etc
#include "init_structs.h"               // for array_a_init, array_c_init, etc

namespace SEAMS {
extern SEAMS::Aprepro *aprepro;

void init_table(const char *comment);
char comm_string[32];
char vers_string[32];
  
struct init_d arith_fncts[] =
{
  {"abs",            do_fabs,   "abs(x)",   "Absolute value of x. |x|."},
  {"acos",           do_acos,   "acos(x)",  "Inverse cosine of x, returns radians."},
  {"acosd",          do_acosd,  "acosd(x)", "Inverse cosine of x, returns degrees."},
  {"acosh",          do_acosh,  "acosh(x)", "Inverse hyperbolic cosine of x."},
  {"asin",           do_asin,   "asin(x)",  "Inverse sine of x, returns radians."},
  {"asind",          do_asind,  "asin(x)",  "Inverse sine of x, returns degrees."},
  {"asinh",          do_asinh,  "asinh(x)", "Inverse hyperbolic sine of x."},
  {"atan",           do_atan,   "atan(x)",  "Inverse tangent of x, returns radians."},
  {"atand",          do_atand,  "atand(x)",  "Inverse tangent of x, returns degrees."},
  {"atanh",          do_atanh,  "atanh(x)",  "Inverse hyperbolic tangent of x."},
  {"ceil",           do_ceil,   "ceil(x)",   "Smallest integer not less than x."},
  {"cos",            do_cos,    "cos(x)",    "Cosine of x, with x in radians"},
  {"cosd",           do_cosd,   "cosd(x)",   "Cosine of x, with x in degrees"},
  {"cosh",           do_cosh,   "cosh(x)",   "Hyperbolic cosine of x."},
  {"d2r",            do_d2r,    "d2r(x)",    "Degrees to radians."},
  {"exp",            do_exp,    "exp(x)",    "Exponential: e^x"},
  {"floor",          do_floor,  "floor(x)",  "Largest integer not greater than x."},
  {"int",            do_int,    "int(x), [x]","Integer part of x truncated toward 0."},
  {"lgamma",         do_lgamma, "lgamma(x)", "log(Gamma(x))."},
  {"ln",             do_log,    "ln(x)",     "Natural (base e) logarithm of x."},
  {"log",            do_log,    "log(x)",    "Natural (base e) logarithm of x."},
  {"log10",          do_log10,  "log10(x)",  "Base 10 logarithm of x. "},
  {"log1p",          do_log1p,  "log1p(x)",  "log(1+x) "},
  {"nint",           do_nint,   "nint(x)",   "Rounds x to nearest integer. <0.5 down; >=0.5 up."},
  {"r2d",            do_r2d,    "r2d(x)",     "Radians to degrees. "},
  {"sin",            do_sin,    "sin(x)","Sine of x, with x in radians. "},
  {"sind",           do_sind,   "sind(x)","Sine of x, with x in degrees. "},
  {"sinh",           do_sinh,   "sinh(x)","Hyperbolic sine of x "},
  {"srand",          do_srand,  "srand(seed)","Seed the random number generator with the given integer value. "},
  {"sqrt",           do_sqrt,   "sqrt(x)","Square root of x. "},
  {"tan",            do_tan,    "tan(x)","Tangent of x, with x in radians. "},
  {"tand",           do_tand,   "tand(x)","Tangent of x, with x in radians. "},
  {"tanh",           do_tanh,   "tanh(x)","Hyperbolic tangent of x. "},
  {0, 0, 0, 0}				/* Last line must be 0, 0 */
};

struct init_a  arith_a_fncts[] =
{
  {"rows",           do_rows,   "rows(array)","Returns the number of rows in the array. "},
  {"cols",           do_cols,   "cols(array)","Returns the number of columns in the array. "},
  {0, 0, 0, 0}				/* Last line must be 0, 0 */
};

struct init_dd arith_dd_fncts[] =
{
  {"atan2",          do_atan2,  "atan2(x,y)","Inverse tangent of y/x, returns radians."},
  {"atan2d",         do_atan2d, "atan2d(x,y)","Inverse tangent of y/x, returns degrees."},
  {"dim",            do_dim,    "dim(x,y)",  "x - min(x,y)"},
  {"fmod",           do_fmod,   "fmod(x,y)", "Floating-point remainder of x/y."},
  {"hypot",          do_hypot,  "hypot(x,y)", "sqrt(x^2+y^2)."},
  {"max",            do_max,    "max(x,y)",  "Maximum of x and y. "},
  {"min",            do_min,    "min(x,y)",  "Minimum of x and y. "},
  {"polarX",         do_polarX, "polarX(r,a)","r * cos(a), a is in degrees "},
  {"polarY",         do_polarY, "polarY(r,a)","r * sin(a), a is in degrees "},
  {"rand",           do_rand,   "rand(xl,xh)","Random value between xl and xh; uniformly distributed. "},
  {"rand_normal",    do_rand_normal, "rand_normal(m,s)","Random value normally distributed with mean m and stddev s."},
  {"rand_lognormal", do_rand_lognormal,"rand_lognormal(m,s)","Random value with lognormal distribution with mean m and stddev s."},
  {"rand_weibull",   do_rand_weibull,"rand_weibull(a, b)","Random value with weibull distribution with alpha=a and beta=b. "},
  {"sign",           do_sign,   "sign(x,y)","x * sgn(y)"},
  {0, 0, 0, 0}				/* Last line must be 0, 0 */
};

struct init_dddd arith_dddd_fncts[] =
{
  {"Vangle",         do_angle,  "Vangle(x1,y1,x2,y2)", "Angle (radians) between vector x1_i+y1_j and x2_i+y2_j."},
  {"Vangled",        do_angled, "Vangled(x1,y1,x2,y2)","Angle (degrees) between vector x1_i+y1_j and x2_i+y2_j."},
  {"dist",           do_dist,   "dist(x1,y1, x2,y2)", "sqrt((x1-x2)^2 + (y1-y2)^2)"},
  {0, 0, 0, 0}				/* Last line must be 0, 0 */
};

struct init_cc arith_cc_fncts[] =
{
  {"word_count",     do_word_count,"word_count(svar,del)","Number of words in svar. Words are separated by one or more of the characters in the string variable 'del'."},
  {0, 0, 0, 0}				/* Last line must be 0, 0 */
};

struct init_c arith_c_fncts[] =
{
  {"strtod",         do_strtod, "strtod(svar)","Returns a double-precision floating-point number equal to the value represented by the character string pointed to by svar."},
  {0, 0, 0, 0}				/* Last line must be 0, 0 */
};

struct init_cd arith_cd_fncts[] =
{
  {"option", do_option, "option(?,?)", "Internal"},
  {0, 0, 0, 0} /* Last line must be 0, 0, 0, 0*/
};

struct init_ddd arith_ddd_fncts[] =
{
  {"julday", do_julday, "julday(mm, dd, yy)", "Julian day corresponding to mm/dd/yy. "},
  {0, 0, 0, 0} /* Last line must be 0, 0, 0, 0*/
};

struct init_dddddd arith_dddddd_fncts[] =
{
  {"juldayhms", do_juldayhms, "juldayhms(mm, dd, yy, hh, mm, ss)", "Julian day corresponding to mm/dd/yy at hh:mm:ss "},
  {0, 0, 0, 0} /* Last line must be 0, 0, 0, 0*/
};

struct str_init string_fncts[] =
{
  {"DUMP",           do_dumpsym,    "DUMP()","Output a list of all defined variables and their value."},
  {"DUMP_FUNC",      do_dumpfunc,   "DUMP_FUNC()","Output a list of all double and string functions recognized by aprepro."},
  {"DUMP_PREVAR",    do_dumpvar,    "DUMP_PREVAR()","Output a list of all predefined variables and their value."},
  {"get_date",       do_get_date,   "get_date()","Returns a string representing the current date in the form YYYY/MM/DD."},
  {"get_iso_date",   do_get_iso_date,"get_iso_date()","Returns a string representing the current date in the form YYYYMMDD."},
  {"get_time",       do_get_time,   "get_time()","Returns a string representing the current time in the form HH:MM:SS."},
  {"get_temp_filename", do_get_temp_filename, "get_temp_filename()", "Returns a string which can be used for a temporary filename without conflicting with any other filenames."},
  {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
};

struct str_c_init string_c_fncts[] =
{
  {"tolower",        do_tolower,     "tolower(svar)","Translates all uppercase characters in svar to lowercase. It modifies svar and returns the resulting string.  "},
  {"toupper",        do_toupper,     "toupper(svar)","Translates all lowercase character in svar to uppercase. It modifies svar and returns the resulting string. "},
  {"to_lower",       do_tolower,     "to_lower(svar)","Translates all uppercase characters in svar to lowercase. It modifies svar and returns the resulting string.  "},
  {"to_upper",       do_toupper,     "toupper(svar)","Translates all lowercase character in svar to uppercase. It modifies svar and returns the resulting string. "},
  {"getenv",         do_getenv,      "getenv(svar)","Returns a string containing the value of the environment variable svar. If the environment variable is not defined, an empty string is returned. "},
  {"file_to_string", do_file_to_string,    "file_to_string(fn)","Opens the file specified by fn and returns the contents as a multi-line string."},
  {"error",          do_error,       "error(svar)","Outputs the string svar to stderr and then terminates the code with an error exit status."},
  {"execute",        do_execute,     "execute(svar)","svar is parsed and executed as if it were a line read from the input file."},
  {"output",         do_output,      "output(filename)","Creates the file specified by filename and sends all subsequent output from aprepro to that file."},
  {"output_append",  do_append,      "output_append(fn)","If file with name fn exists, append output to it; otherwise create the file and send all subsequent output from aprepro to that file."},
  {"rescan",         do_rescan,      "rescan(svar)","The difference between execute(sv1) and rescan(sv2) is that sv1 must be a valid expression, but sv2 can contain zero or more expressions. "},
  {"include_path",   do_include_path,"include_path(path)","Specify an optional path to be prepended to a filename when opening a file. Can also be specified via the -I command line option when executing aprepro."},
  {"Units",          do_Units,       "Units(svar)","See manual. svar is one of the defined units systems:\n\t\t\t'si', 'cgs', 'cgs-ev', 'shock', 'swap', 'ft-lbf-s', 'ft-lbm-s', 'in-lbf-s'"},
#if defined(EXODUSII)
  {"exodus_info",    do_exodus_info, "exodus_info(ex_fn)","Parses the info records extracted from the exodus file 'ex_fn'"},
  {"exodus_meta",    do_exodus_meta, "exodus_meta(ex_fn)","Creates several variables related to the exodusII metadata in the specified file. Experimental."},
#endif
  {"delete",         do_delete,      "delete(var_name)", "Delete the variable with name 'var_name'."},
  {"if",             do_str_if,      "if(x)", "Handles the if statements. x can be any valid expression; nonzero is true"},
  {"If",             do_str_if,      "If(x)", "Handles the if statements. x can be any valid expression; nonzero is true"},
  {"elseif",         do_str_elseif,  "elseif(x)", "Handles the if statements. x can be any valid expression; nonzero is true"},
  {"Elseif",         do_str_elseif,  "Elseif(x)", "Handles the if statements. x can be any valid expression; nonzero is true"},
  {"ifdef",          do_str_if,      "ifdef(x)", "Handles the if statements. x can be any valid expression; nonzero is true (deprecated, use if)"},
  {"Ifdef",          do_str_if,      "Ifdef(x)", "Handles the if statements. x can be any valid expression; nonzero is true (deprecated, use if)"},
  {"ifndef",         do_str_notif,   "ifndef(x)", "Handles the if statements. x can be any valid expression; nonzero is true (deprecated, use if)"},
  {"Ifndef",         do_str_notif,   "Ifndef(x)", "Handles the if statements. x can be any valid expression; nonzero is true (deprecated, use if)"},
  {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
};
  
struct str_d_init string_d_fncts[] =
  {
  {"IO",             do_intout,      "IO(x)","Convert x to an integer and then to a string. "},
  {"to_string",      do_tostring,    "to_string(x)","Returns a string representation of the numerical variable x. The variable x is unchanged."},
  {"tostring",       do_tostring,    "tostring(x)","Returns a string representation of the numerical variable x. The variable x is unchanged."},
  {"if",             do_if,          "if(x)", "Handles the if statements. x can be any valid expression; nonzero is true"},
  {"If",             do_if,          "If(x)", "Handles the if statements. x can be any valid expression; nonzero is true"},
  {"elseif",         do_elseif,      "elseif(x)", "Handles the if statements. x can be any valid expression; nonzero is true"},
  {"Elseif",         do_elseif,      "Elseif(x)", "Handles the if statements. x can be any valid expression; nonzero is true"},
  {"ifdef",          do_if,          "ifdef(x)", "Handles the if statements. x can be any valid expression; nonzero is true (deprecated, use if)"},
  {"Ifdef",          do_if,          "Ifdef(x)", "Handles the if statements. x can be any valid expression; nonzero is true (deprecated, use if)"},
  {"ifndef",         do_notif,       "ifndef(x)", "Handles the if statements. x can be any valid expression; nonzero is true (deprecated, use if)"},
  {"Ifndef",         do_notif,       "Ifndef(x)", "Handles the if statements. x can be any valid expression; nonzero is true (deprecated, use if)"},
  {"switch",         do_switch,      "switch(x)", "Switch statement. Select from the following case statements which matches 'x' and execute that one. End with endswitch"},
  {"Switch",         do_switch,      "Switch(x)", "Switch statement. Select from the following case statements which matches 'x' and execute that one. End with endswitch"},
  {"case",           do_case,        "case(x)", "Switch statement. A case used in a containing switch statement."},
  {"Case",           do_case,        "Case(x)", "Switch statement. A case used in a containing switch statement."},
  {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
  };

struct str_dcc_init string_dcc_fncts[] =
  {
  {"get_word",       do_get_word,    "get_word(n,svar,del)","Returns a string containing the nth word of svar. The words are separated by one or more of the characters in the string variable del "},
  {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
  };

struct str_ccc_init string_ccc_fncts[] =
  {
  {"extract",        do_extract,     "extract(s, b, e)","Return substring [b,e). 'b' is included; 'e' is not. If 'b' not found, return empty; If 'e' not found, return rest of string. If 'b' empty, start at beginning; if 'e' empty, return rest of string."},
  {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
  };

struct str_a_init string_a_fncts[] =
  {
    {"print_array",  do_print_array, "print_array(array)","Prints the data in the array."},
    {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
  };

struct array_c_init array_c_fncts[] =
  {
    {"csv_array",         do_csv_array1,      "csv_array(filename)",
     "Create a 2D array from the data in a csv file."},
    {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
  };

struct array_cd_init array_cd_fncts[] =
  {
    {"csv_array",         do_csv_array,      "csv_array(filename, [skip])",
     "Create a 2D array from the data in a csv file optionally skipping rows."
     " If skip is integer skip that many rows; if skip is a character, skip lines beginning with that character"},
    {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
  };

struct array_cc_init array_cc_fncts[] =
  {
    {"csv_array",         do_csv_array2,      "csv_array(filename, [skip])",
     "Create a 2D array from the data in a csv file optionally skipping rows."
     " If skip is integer skip that many rows; if skip is a character, skip lines beginning with that character"},
    {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
  };

struct array_dd_init array_dd_fncts[] =
  {
    {"make_array",        do_make_array,     "make_array(rows, cols)",
     "Create a 2D array of size 'rows' by 'cols' initialized to zero."},
    {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
  };

struct array_d_init array_d_fncts[] =
  {
    {"identity",          do_identity,     "identity(size)",
     "Create a 2D identity array with 'size' rows and columns. Diagonal = 1.0"},
    {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
  };

struct array_a_init array_a_fncts[] =
  {
    {"transpose",         do_transpose,      "transpose(array)",
     "Return the transpose of input array"},
    {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
  };

struct var_init variables[] =
{
  {"DEG",  57.29577951308232087680},	/* 180/pi, degrees per radian */
  {"RAD",   0.01745329251994329576},	/* pi/180, radians per degree */
  {"E",     2.71828182845904523536},	/* e, base of natural log     */
  {"GAMMA", 0.57721566490153286060},	/* euler-mascheroni constant  */
  {"PHI",   1.61803398874989484820},	/* golden ratio               */
  {"TAU",   6.28318530717958623200},    /* 2*PI see Tau Manifesto, http://tauday.com */
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

  void Aprepro::init_table(const char *comment)
  {
    for (int i = 0; arith_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(arith_fncts[i].fname, FUNCTION, 1);
      ptr->value.fnctptr_d = arith_fncts[i].fnct;
      ptr->info = arith_fncts[i].description;
      ptr->syntax = arith_fncts[i].syntax;
    }

    for (int i = 0; arith_dd_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(arith_dd_fncts[i].fname, FUNCTION, 1);
      ptr->value.fnctptr_dd = arith_dd_fncts[i].fnct;
      ptr->info = arith_dd_fncts[i].description;
      ptr->syntax = arith_dd_fncts[i].syntax;
    }

    for (int i = 0; arith_a_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(arith_a_fncts[i].fname, FUNCTION, 1);
      ptr->value.fnctptr_a = arith_a_fncts[i].fnct;
      ptr->info = arith_a_fncts[i].description;
      ptr->syntax = arith_a_fncts[i].syntax;
    }

    for (int i = 0; arith_dddd_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(arith_dddd_fncts[i].fname, FUNCTION, 1);
      ptr->value.fnctptr_dddd = arith_dddd_fncts[i].fnct;
      ptr->info = arith_dddd_fncts[i].description;
      ptr->syntax = arith_dddd_fncts[i].syntax;
    }

    for (int i = 0; arith_cc_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(arith_cc_fncts[i].fname, FUNCTION, 1);
      ptr->value.fnctptr_cc = arith_cc_fncts[i].fnct;
      ptr->info = arith_cc_fncts[i].description;
      ptr->syntax = arith_cc_fncts[i].syntax;
    }

    for (int i = 0; arith_c_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(arith_c_fncts[i].fname, FUNCTION, 1);
      ptr->value.fnctptr_c = arith_c_fncts[i].fnct;
      ptr->info = arith_c_fncts[i].description;
      ptr->syntax = arith_c_fncts[i].syntax;
    }

    for (int i = 0; arith_cd_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(arith_cd_fncts[i].fname, FUNCTION, 1);
      ptr->value.fnctptr_cd = arith_cd_fncts[i].fnct;
      ptr->info = arith_cd_fncts[i].description;
      ptr->syntax = arith_cd_fncts[i].syntax;
    }

    for (int i = 0; arith_ddd_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(arith_ddd_fncts[i].fname, FUNCTION, 1);
      ptr->value.fnctptr_ddd = arith_ddd_fncts[i].fnct;
      ptr->info = arith_ddd_fncts[i].description;
      ptr->syntax = arith_ddd_fncts[i].syntax;
    }

    for (int i = 0; arith_dddddd_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(arith_dddddd_fncts[i].fname, FUNCTION, 1);
      ptr->value.fnctptr_dddddd = arith_dddddd_fncts[i].fnct;
      ptr->info = arith_dddddd_fncts[i].description;
      ptr->syntax = arith_dddddd_fncts[i].syntax;
    }

    for (int i = 0; string_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(string_fncts[i].fname, STRING_FUNCTION, 1);
      ptr->value.strfnct = string_fncts[i].fnct;
      ptr->info = string_fncts[i].description;
      ptr->syntax = string_fncts[i].syntax;
    }

    for (int i = 0; string_c_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(string_c_fncts[i].fname, STRING_FUNCTION, 1);
      ptr->value.strfnct_c = string_c_fncts[i].fnct;
      ptr->info = string_c_fncts[i].description;
      ptr->syntax = string_c_fncts[i].syntax;
    }

    for (int i = 0; string_d_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(string_d_fncts[i].fname, STRING_FUNCTION, 1);
      ptr->value.strfnct_d = string_d_fncts[i].fnct;
      ptr->info = string_d_fncts[i].description;
      ptr->syntax = string_d_fncts[i].syntax;
    }

    for (int i = 0; string_dcc_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(string_dcc_fncts[i].fname, STRING_FUNCTION, 1);
      ptr->value.strfnct_dcc = string_dcc_fncts[i].fnct;
      ptr->info = string_dcc_fncts[i].description;
      ptr->syntax = string_dcc_fncts[i].syntax;
    }

    for (int i = 0; string_ccc_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(string_ccc_fncts[i].fname, STRING_FUNCTION, 1);
      ptr->value.strfnct_ccc = string_ccc_fncts[i].fnct;
      ptr->info = string_ccc_fncts[i].description;
      ptr->syntax = string_ccc_fncts[i].syntax;
    }

    for (int i = 0; string_a_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(string_a_fncts[i].fname, STRING_FUNCTION, 1);
      ptr->value.strfnct_a = string_a_fncts[i].fnct;
      ptr->info = string_a_fncts[i].description;
      ptr->syntax = string_a_fncts[i].syntax;
    }

    for (int i = 0; array_c_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(array_c_fncts[i].fname, ARRAY_FUNCTION, 1);
      ptr->value.arrfnct_c = array_c_fncts[i].fnct;
      ptr->info = array_c_fncts[i].description;
      ptr->syntax = array_c_fncts[i].syntax;
    }

    for (int i = 0; array_cc_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(array_cc_fncts[i].fname, ARRAY_FUNCTION, 1);
      ptr->value.arrfnct_cc = array_cc_fncts[i].fnct;
      ptr->info = array_cc_fncts[i].description;
      ptr->syntax = array_cc_fncts[i].syntax;
    }

    for (int i = 0; array_cd_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(array_cd_fncts[i].fname, ARRAY_FUNCTION, 1);
      ptr->value.arrfnct_cd = array_cd_fncts[i].fnct;
      ptr->info = array_cd_fncts[i].description;
      ptr->syntax = array_cd_fncts[i].syntax;
    }

    for (int i = 0; array_dd_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(array_dd_fncts[i].fname, ARRAY_FUNCTION, 1);
      ptr->value.arrfnct_dd = array_dd_fncts[i].fnct;
      ptr->info = array_dd_fncts[i].description;
      ptr->syntax = array_dd_fncts[i].syntax;
    }

    for (int i = 0; array_d_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(array_d_fncts[i].fname, ARRAY_FUNCTION, 1);
      ptr->value.arrfnct_d = array_d_fncts[i].fnct;
      ptr->info = array_d_fncts[i].description;
      ptr->syntax = array_d_fncts[i].syntax;
    }

    for (int i = 0; array_a_fncts[i].fname != 0; i++) {
      symrec *ptr = putsym(array_a_fncts[i].fname, ARRAY_FUNCTION, 1);
      ptr->value.arrfnct_a = array_a_fncts[i].fnct;
      ptr->info = array_a_fncts[i].description;
      ptr->syntax = array_a_fncts[i].syntax;
    }

    for (int i = 0; variables[i].vname != 0; i++) {
      symrec *ptr = putsym(variables[i].vname, VARIABLE, 1);
      ptr->value.var = variables[i].value;
    }

    for (int i = 0; svariables[i].vname != 0; i++) {
      symrec *ptr = putsym(svariables[i].vname, STRING_VARIABLE, 1);
      ptr->value.svar = svariables[i].value;
    }

    comm_string[0] = '#';
    
    symrec *ptr = putsym("_C_", STRING_VARIABLE, 1);
    ptr->value.svar = comm_string;

    {
      std::strncpy(vers_string, aprepro->version().c_str(), 32);
      vers_string[31] = '\0';
      ptr = putsym("VERSION", STRING_VARIABLE, 1);
      ptr->value.svar = vers_string;
    }
  }
} // namespace SEAMS
