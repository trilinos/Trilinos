/***
   NAME
     init
   PURPOSE
     Initialize variables and functions Aprepro
***/
#include "aprepro.h"
#include <cstdio>
#include <sys/types.h>
#include "init_structs.h"

namespace SEAMS {
extern SEAMS::Aprepro *aprepro;

void init_table(char comment);
char comm_string[2];

extern double
  do_fabs(double x),
  do_acos(double x),
  do_acosd(double x),
  do_asin(double x),
  do_asind(double x),
  do_atan(double x),
  do_atand(double x),
  do_ceil(double x),
  do_cos(double x),
  do_cosd(double x),
  do_cosh(double x),
  do_d2r(double x),
  do_exp(double x),
  do_floor(double x),
  do_int(double x),
  do_log(double x),
  do_log10(double x),
  do_r2d(double x),
  do_sin(double x),
  do_sind(double x),
  do_sinh(double x),
  do_sqrt(double x),
  do_tan(double x),
  do_tand(double x),
  do_tanh(double x),
  do_lgamma(double val),
  do_log1p(double x),
  do_acosh(double x),
  do_asinh(double x),
  do_atanh(double x),
  do_nint(double x),
  do_atan2(double x, double y),
  do_atan2d(double x, double y),
  do_dim(double x, double y),
  do_fmod(double x, double y),
  do_max(double x, double y),
  do_min(double x, double y),
  do_srand(double seed),
  do_rand(double xl, double xh),
  do_sign(double x, double y),
  do_rand_normal(double mean, double stddev),
  do_rand_weibull(double alpha, double beta),
  do_rand_lognormal(double mean, double stddev), 
  do_hypot(double x, double y),
  do_polarX(double rad, double ang),
  do_polarY(double rad, double ang),
  do_dist(double x1, double y1, double x2, double y2),
  do_angle(double x1, double y1, double x2, double y2),
  do_angled(double x1, double y1, double x2, double y2),
  do_word_count(char *string, char *delm),
  do_strtod(char *string);

#if 0
  do_julday(double mon, double day, double year),
  do_juldayhms(double mon, double day, double year, double h, double mi, double se),
  do_option(char *option, double value);
#endif

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

//TODO
#if 0  
  {"option",         do_option, "option(?,?)","Internal"},
  {"julday",         do_julday, "julday(mm, dd, yy)","Julian day corresponding to mm/dd/yy. "},
  {"juldayhms",      do_juldayhms,"juldayhms(mm, dd, yy, hh, mm, ss)","Julian day corresponding to mm/dd/yy at hh:mm:ss "},
  {0, 0, 0, 0}				/* Last line must be 0, 0 */
};
#endif

extern char
  *do_dumpsym(),
  *do_dumpfunc(),
  *do_dumpvar(),
  *do_get_date(),
  *do_get_iso_date(),
  *do_get_time(),
  *do_get_temp_filename(),
  
  *do_tolower(char *string),
  *do_toupper(char *string),
  *do_Units(char *type),
  *do_file_to_string(char *filename),
  *do_error(char *error_string),
  *do_include_path(char *new_path),
  *do_getenv(char *env),
  *do_output(char *filename),
  *do_append(char *filename),
  *do_execute(char *string),
  *do_rescan(char *string),

  *do_intout(double intval),
  *do_tostring(double x),

  *do_get_word(double n, char *string, char *delm),
  *do_extract(char *string, char *begin, char *end);
  
#if defined(EXODUSII)
  *do_exodus_info(char *filename),
  *do_exodus_meta(char *filename),
#endif


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
  {0, 0, 0, 0}				/* Last line must be 0, 0, 0, 0 */
};
  
struct str_d_init string_d_fncts[] =
  {
  {"IO",             do_intout,      "IO(x)","Convert x to an integer and then to a string. "},
  {"to_string",      do_tostring,    "to_string(x)","Returns a string representation of the numerical variable x. The variable x is unchanged."},
  {"tostring",       do_tostring,    "tostring(x)","Returns a string representation of the numerical variable x. The variable x is unchanged."},
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

  void Aprepro::init_table(char comment)
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

    for (int i = 0; variables[i].vname != 0; i++) {
      symrec *ptr = putsym(variables[i].vname, VARIABLE, 1);
      ptr->value.var = variables[i].value;
    }

    for (int i = 0; svariables[i].vname != 0; i++) {
      symrec *ptr = putsym(svariables[i].vname, STRING_VARIABLE, 1);
      ptr->value.svar = svariables[i].value;
    }

    std::sprintf(comm_string, "%c", comment);

    symrec *ptr = putsym("_C_", STRING_VARIABLE, 1);
    ptr->value.svar = comm_string;

    {
      ptr = putsym("VERSION", STRING_VARIABLE, 1);
      ptr->value.svar = aprepro->version().c_str();
    }
  }
} // namespace SEAMS
