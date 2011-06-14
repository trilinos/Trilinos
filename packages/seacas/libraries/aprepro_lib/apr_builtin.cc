#include <cmath>
#include <cctype>
#include <errno.h>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <assert.h>
#include "aprepro.h"
#include "apr_scanner.h"
#include "aprepro_parser.h"
#include "apr_util.h"
#include "apr_tokenize.h"

#ifndef RAND_MAX
#include <limits.h>
#define RAND_MAX INT_MAX
#endif

#ifndef PI
#define PI  3.141592653589793238462643
#endif

namespace SEAMS {

extern SEAMS::Aprepro *aprepro;

#if defined(VMS) || defined(_hpux_) || defined(sun) || defined(__linux__)
#define HYPOT(x,y) hypot(x,y)
#else
#define HYPOT(x,y) do_hypot(x,y)
#endif

#define d2r(x)  ((x)*PI/180.)
#define r2d(x)  ((x)*180./PI)

#ifndef max
#define max(x,y) (x) > (y) ? (x) : (y)
#define min(x,y) (x) < (y) ? (x) : (y)
#endif

#if defined(sun) || defined(__linux__)
#define LOG1P(x)	log1p(x)
#else
#define LOG1P(x)	std::log(1.0 + (x))
#endif

double do_acos(double x);
double do_acosd(double x);
double do_acosh(double x);
double do_angle(double x1, double y1, double x2, double y2);
double do_angled(double x1, double y1, double x2, double y2);
double do_asin(double x);
double do_asind(double x);
double do_asinh(double x);
double do_atan(double x);
double do_atan2(double x, double y);
double do_atan2d(double x, double y);
double do_atand(double x);
double do_atanh(double x);
double do_ceil(double x);
double do_cos(double x);
double do_cosd(double x);
double do_cosh(double x);
double do_d2r(double x);
double do_dim(double x, double y);
double do_dist(double x1, double y1, double x2, double y2);
double do_exp(double x);
double do_fabs(double x);
double do_floor(double x);
double do_nint(double x);
double do_fmod(double x, double y);
double do_hypot(double x, double y);
double do_int(double x);
double do_log(double x);
double do_log10(double x);
double do_max(double x, double y);
double do_min(double x, double y);
double do_r2d(double x);
double do_rand(double xl, double xh);
double do_srand(double seed);
double do_rand_normal(double mean, double stddev);
double do_rand_lognormal(double mean, double stddev);
double do_rand_weibull(double alpha, double beta);
double do_sign(double x, double y);
double do_sin(double x);
double do_sind(double x);
double do_sinh(double x);
double do_sqrt(double x);
double do_tan(double x);
double do_tand(double x);
double do_tanh(double x);
double do_polarX(double rad, double ang);
double do_polarY(double rad, double ang);
double do_strtod(char *string);

const char  *do_execute(char *string);
const char  *do_getenv(char *string);
const char  *do_tolower(char *string);
const char  *do_toupper(char *string);
const char  *do_tostring(double x);
const char  *do_output(char *newfile);
const char  *do_append(char *newfile);
const char  *do_error(char *error_string);
const char  *do_get_date(void);
const char  *do_get_iso_date(void);
const char  *do_get_time(void);
double do_option(char *option, double value);
double do_word_count(char *string, char *delm );
const char  *do_get_word(double n, char *string, char *delm);
const char  *do_file_to_string(char *filename);
const char  *do_extract(char *string, char *begin, char *end);
double do_Material(double id, char *type, char *name, char *model, char *code, FILE * yyout);
double do_lgamma(double val);
double do_juldayhms(double mon, double day, double year,
			    double h, double mi, double se);
double do_julday(double mon, double day, double year);
double do_log1p(double mag);
const char  *do_include_path(char *newpath);
const char  *do_intout(double intval);

/* DO_INT:  Calculate integer nearest to zero from value */
double do_int(double x)
{
  double temp;
  errno = 0;
  temp = (double) (x < 0 ? -std::floor (-(x)) : std::floor (x));
  SEAMS::math_error ("int");
  return (temp);
}

/* DO_NINT:  Calculate integer nearest value */
double do_nint(double x)
{
  double temp;
  errno = 0;
  temp = (double) (x < 0 ? -std::floor(0.5-x) : std::floor(x+0.5));
  SEAMS::math_error ("nint");
  return (temp);
}

/* DO_DIST: Calculate distance between point 1 at (x1,y1) and
 *          point 2 at (x2,y2).
 */
double do_dist(double x1, double y1, double x2, double y2)
{
  double temp;
  errno = 0;
  temp = HYPOT((x1 - x2), (y1 - y2));
  SEAMS::math_error("hypot");
  return (temp);
}

/* DO_ANGLE: Calculate angle (radians) between vector 1 at (0,0; x1,y1) and
 *          vector 2 at (0,0; x2,y2).
 */
double do_angle(double x1, double y1, double x2, double y2)
{
  double temp;
  temp = ((x1 * x2) + (y1 * y2)) / (HYPOT(x1, y1) * HYPOT(x2, y2));
  errno = 0;
  temp = acos(temp);
  SEAMS::math_error("angle");
  return (temp);
}

/* DO_ANGLE: Calculate angle (degrees) between vector 1 at (0,0; x1,y1) and
 *          vector 2 at (0,0; x2,y2).
 */
double do_angled(double x1, double y1, double x2, double y2)
{
  double temp;
  temp = ((x1 * x2) + (y1 * y2)) / (HYPOT(x1, y1) * HYPOT(x2, y2));
  errno = 0;
  temp = r2d(acos(temp));
  SEAMS::math_error("angled");
  return (temp);
}

/* DO_HYPOT: calcluate sqrt(p^2 + q^2)     */
/* Algorithm from "More Programming Pearls," Jon Bentley */
/* Accuracy: 6.5 digits after 2 iterations,
 *           20  digits after 3 iterations,
 *           62  digits after 4 iterations.
 */

double do_hypot(double x, double y)
{
  double r;
  int i;

  x = fabs(x);
  y = fabs(y);
  if (x < y)
    {
      r = y;
      y = x;
      x = r;
    }
  if (x == 0.0)
    return (y);

  for (i = 0; i < 3; i++)
    {
      r = y / x;
      r *= r;
      r /= (4.0 + r);
      x += (2.0 * r * x);
      y *= r;
    }
  return (x);
}

double 
do_max(double x, double y)
{
  double temp;
  errno = 0;
  temp = max(x, y);
  SEAMS::math_error("max");
  return (temp);
}

double 
do_min(double x, double y)
{
  double temp;
  errno = 0;
  temp = min(x, y);
  SEAMS::math_error("min");
  return (temp);
}

double 
do_d2r(double x)
{
  return (d2r(x));
}

double 
do_r2d(double x)
{
  return (r2d(x));
}

double 
do_sind(double x)
{
  double temp;
  errno = 0;
  temp = sin(d2r(x));
  SEAMS::math_error("sind");
  return (temp);
}

double 
do_sin(double x)
{
  double temp;
  errno = 0;
  temp = sin(x);
  SEAMS::math_error("sin");
  return (temp);
}

double 
do_cosd(double x)
{
  double temp;
  errno = 0;
  temp = cos(d2r(x));
  SEAMS::math_error("cosd");
  return (temp);
}

double 
do_cos(double x)
{
  double temp;
  errno = 0;
  temp = cos(x);
  SEAMS::math_error("cos");
  return (temp);
}

double 
do_tand(double x)
{
  double temp;
  errno = 0;
  temp = tan(d2r(x));
  SEAMS::math_error("tand");
  return (temp);
}

double 
do_tan(double x)
{
  double temp;
  errno = 0;
  temp = tan(x);
  SEAMS::math_error("tan");
  return (temp);
}

double 
do_atan2d(double x, double y)
{
  double temp;
  errno = 0;
  temp = r2d(atan2(x, y));
  SEAMS::math_error("atan2d");
  return (temp);
}

double 
do_atan2(double x, double y)
{
  double temp;
  errno = 0;
  temp = atan2(x, y);
  SEAMS::math_error("atan2");
  return (temp);
}

double 
do_atand(double x)
{
  double temp;
  errno = 0;
  temp = r2d(atan(x));
  SEAMS::math_error("atand");
  return (temp);
}

double 
do_atan(double x)
{
  double temp;
  errno = 0;
  temp = atan(x);
  SEAMS::math_error("atan");
  return (temp);
}

double 
do_asind(double x)
{
  double temp;
  errno = 0;
  temp = r2d(asin(x));
  SEAMS::math_error("asind");
  return (temp);
}

double 
do_asin(double x)
{
  double temp;
  errno = 0;
  temp = asin(x);
  SEAMS::math_error("asin");
  return (temp);
}

double 
do_acosd(double x)
{
  double temp;
  errno = 0;
  temp = r2d(acos(x));
  SEAMS::math_error("acosd");
  return (temp);
}

double 
do_acos(double x)
{
  double temp;
  errno = 0;
  temp = acos(x);
  SEAMS::math_error("acos");
  return (temp);
}

/* do_srand(x) Seed the random generator with the specified integer value */
double
do_srand(double seed)
{
  srand((unsigned)seed);
  return (0);
}

/* do_rand(x) returns a random double in the range 0<= do_rand <= x */
double 
do_rand(double xl, double xh)
{
  double temp;
  errno = 0;
  temp = xl + (xh - xl) * ((double) rand() / (double) RAND_MAX);
  SEAMS::math_error("rand");
  return (temp);
}

double 
do_rand_normal(double mean, double stddev)
{
  /* boxmuller.c
     Implements the Polar form of the Box-Muller Transformation
   
     (c) Copyright 1994, Everett F. Carter Jr.  Permission is granted by
     the author to use this software for any application provided this
     copyright notice is preserved.
  */
  double x1, x2, w, y1;
  static double y2;
  static int use_last = 0;
  
  if (use_last) {
    y1 = y2;
    use_last = 0;
  }
  else {
    do {
      x1 = 2.0 * ((double)rand()/(double)RAND_MAX) - 1.0;
      x2 = 2.0 * ((double)rand()/(double)RAND_MAX) - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = std::sqrt( (-2.0 * std::log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }
  return ( mean + y1 * stddev);
}

double 
do_rand_lognormal(double mean, double stddev)
{
  double x;
  double logstd = std::log(1.0 + (stddev/mean)*(stddev/mean));
  double logmean = std::log(mean) - 0.5 * logstd;
  logstd = std::sqrt(logstd);

  x = do_rand_normal(logmean, logstd);

  return exp(x);
}

double 
do_rand_weibull(double alpha, double beta)
{
  double temp = (double) rand() / (double) RAND_MAX;
  errno = 0;
  temp = pow( (-1.0 / alpha * std::log(1.0 - temp)), (1.0/beta) );
  SEAMS::math_error("weibull");
  return (temp);
}

double 
do_sign(double x, double y)
{
  double temp;
  errno = 0;
  temp = (y) >= 0 ? fabs(x) : -fabs(x);
  SEAMS::math_error("sign");
  return (temp);
}

double 
do_dim(double x, double y)
{
  double temp;
  errno = 0;
  temp = x - (min(x, y));
  SEAMS::math_error("dim");
  return (temp);
}

double 
do_fabs(double x)
{
  double temp;
  errno = 0;
  temp = fabs(x);
  SEAMS::math_error("fabs");
  return (temp);
}

double 
do_ceil(double x)
{
  double temp;
  errno = 0;
  temp = ceil(x);
  SEAMS::math_error("ceil");
  return (temp);
}

double 
do_cosh(double x)
{
  double temp;
  errno = 0;
  temp = cosh(x);
  SEAMS::math_error("cosh");
  return (temp);
}

double 
do_exp(double x)
{
  double temp;
  errno = 0;
  temp = exp(x);
  SEAMS::math_error("exp");
  return (temp);
}

double 
do_floor(double x)
{
  double temp;
  errno = 0;
  temp = floor(x);
  SEAMS::math_error("floor");
  return (temp);
}

double 
do_fmod(double x, double y)
{
  double temp;
  errno = 0;
  temp = fmod(x, y);
  SEAMS::math_error("fmod");
  return (temp);
}

double 
do_log(double x)
{
  double temp;
  errno = 0;
  temp = std::log(x);
  SEAMS::math_error("log");
  return (temp);
}

double 
do_log10(double x)
{
  double temp;
  errno = 0;
  temp = std::log10(x);
  SEAMS::math_error("log10");
  return (temp);
}

double 
do_sinh(double x)
{
  double temp;
  errno = 0;
  temp = sinh(x);
  SEAMS::math_error("sinh");
  return (temp);
}

double 
do_sqrt(double x)
{
  double temp;
  errno = 0;
  temp = std::sqrt(x);
  SEAMS::math_error("sqrt");
  return (temp);
}

double 
do_tanh(double x)
{
  double temp;
  errno = 0;
  temp = tanh(x);
  SEAMS::math_error("tanh");
  return (temp);
}

double do_polarX(double rad, double ang)
{
  return (rad * cos(d2r(ang)));
}

double do_polarY(double rad, double ang)
{
  return (rad * sin(d2r(ang)));
}

double cof[] =
{76.18009173, -86.50532033, 24.01409822,
 -1.231739516, 0.120858003e-2, -0.536382e-5};
double do_lgamma(double val)
{
#define STP	2.50662827465
  double x, tmp, ser;
  int j;

  x = val - 1.0;
  tmp = x + 5.5;
  tmp = (x + 0.5) * std::log(tmp) - tmp;
  ser = 1.0;
  for (j = 0; j < 6; j++)
    {
      x += 1.0;
      ser += (cof[j] / x);
    }
  return (tmp + std::log(STP * ser));
}

double do_juldayhms(double mon, double day, double year,
		    double h, double mi, double se)
{
  long m = (long)mon, d = (long)day, y = (long)year;
  long c, ya, j;
  double seconds = h * 3600.0 + mi * 60 + se;

  if (m > 2)
    m -= 3;
  else
    {
      m += 9;
      --y;
    }
  c = y / 100L;
  ya = y - (100L * c);
  j = (146097L * c) / 4L + (1461L * ya) / 4L + (153L * m + 2L) / 5L + d + 1721119L;
  if (seconds < 12 * 3600.0)
    {
      j--;
      seconds += 12.0 * 3600.0;
    }
  else
    {
      seconds = seconds - 12.0 * 3600.0;
    }
  return (j + (seconds / 3600.0) / 24.0);
}

double do_julday(double mon, double day, double year)
{
  return do_juldayhms(mon, day, year, 0.0, 0.0, 0.0);
}

double do_log1p(double x)
{
  return LOG1P(x);
}

double do_acosh(double x)
{
  double t;
  if (x > 1.0e20)
    return (LOG1P(x) + std::log(2.0));
  else
    {
      t = std::sqrt(x - 1.0);
      return (LOG1P(t * (t + std::sqrt(x + 1))));
    }
}

double do_asinh(double x)
{
  double s, t;
  if (1.0 + x * x == 1.0)
    return (x);
  if (std::sqrt(1.0 + x * x) == 1.0)
    return (do_sign(1.0, x) * (LOG1P(x) + std::log(2.0)));
  else
    {
      t = fabs(x);
      s = 1.0 / t;
      return (do_sign(1.0, x) * LOG1P(t + t / (s + std::sqrt(1 + s * s))));
    }
}

double do_atanh(double x)
{
  double z;
  z = do_sign(0.5, x);
  x = do_sign(x, 1.0);
  x = x / (1.0 - x);
  return (z * LOG1P(x + x));
}
/*
  --------------------------STRING FUNCTIONS------------------------
 */

const char *do_get_date(void)
{
  char *tmp;
  const size_t bufsize = 32;
  static char tmpstr[32];

  time_t timer = time(NULL);
  struct tm *timeptr = localtime(&timer);
  
  /* First  the date in the form CCYY/MM/DD */
  strftime(tmpstr, bufsize, "%Y/%m/%d", timeptr);
  new_string(tmpstr, &tmp);
  return(tmp);
}

const char *do_get_iso_date(void)
{
  char *tmp;
  const size_t bufsize = 32;
  static char tmpstr[32];

  time_t timer = time(NULL);
  struct tm *timeptr = localtime(&timer);
  
  /* First  the date in the form CCYY/MM/DD */
  strftime(tmpstr, bufsize, "%Y%m%d", timeptr);
  new_string(tmpstr, &tmp);
  return(tmp);
}

const char *do_get_time(void)
{
  char *tmp;
  const size_t bufsize = 32;
  static char tmpstr[32];

  time_t timer = time(NULL);
  struct tm *timeptr = localtime(&timer);
  
  /* Now the time in the form HH:MM:SS where 0 <= HH < 24 */
  strftime(tmpstr, bufsize, "%H:%M:%S", timeptr);
  new_string(tmpstr, &tmp);
  return(tmp);
}

const char *do_tolower(char *string)
{
  char *p = string;
  while (*p != '\0')
    {
      if (isupper((int)*p))
	*p = tolower((int)*p);
      p++;
    }
  return (string);
}

const char *do_toupper(char *string)
{
  char *p = string;
  while (*p != '\0')
    {
      if (islower((int)*p))
	*p = toupper((int)*p);
      p++;
    }
  return (string);
}

const char *do_tostring(double x)
{
  char *tmp;
  static char tmpstr[128];
  if (x == 0.0)
    {
      new_string("0", &tmp);
      return (tmp);
    }
  else
    {
      SEAMS::symrec *format;
      format = aprepro->getsym("_FORMAT");
      (void) sprintf(tmpstr, format->value.svar, x);
      new_string(tmpstr, &tmp);
      return (tmp);
    }
}

const char *do_output(char *filename)
{
  aprepro->outputStream.top()->flush();

  if (std::strcmp(filename, "pop") == 0 ||
      std::strcmp(filename, "stdout") == 0) {
    while (aprepro->outputStream.size() > 1) {
      std::ostream* output = aprepro->outputStream.top();
      aprepro->outputStream.pop();
      delete output;
    }

    if (aprepro->ap_options.info_msg == true) {
      std::cerr << "Aprepro: INFO: Output now redirected to original output stream.\n";
    }
  }
  else {
    std::ostream* output = new std::ofstream(filename);
    if (output != NULL) {
      aprepro->outputStream.push(output);

      if (aprepro->ap_options.info_msg == true) {
	std::cerr << "Aprepro: INFO: Output now redirected to file '"
		  << filename << "'.\n";
      }
    } else {
	std::cerr << "Aprepro: ERROR: Could not open output file '"
		  << filename << "'.\n";
    }
  }
  return (NULL);
}

const char *do_append(char *filename)
{
  aprepro->outputStream.top()->flush();

  if (std::strcmp(filename, "pop") == 0 ||
      std::strcmp(filename, "stdout") == 0) {
    while (aprepro->outputStream.size() > 1) {
      std::ostream* output = aprepro->outputStream.top();
      aprepro->outputStream.pop();
      delete output;
    }

    if (aprepro->ap_options.info_msg == true) {
      std::cerr << "Aprepro: INFO: Output now redirected to original output stream.\n";
    }
  }
  else {
    std::ofstream* output = new std::ofstream(filename, std::ios_base::app); // Append
    if (output != NULL) {
      aprepro->outputStream.push(output);

      if (aprepro->ap_options.info_msg == true) {
	std::cerr << "Aprepro: INFO: Output now redirected to file '"
		  << filename << "'\n";
      }
    } else {
	std::cerr << "Aprepro: ERROR: Could not open output file '"
		  << filename << "' for appending.\n";
    }
  }
  return (NULL);
}

double do_word_count(char *string, char *delm)
{
  std::string temp = string;
  std::vector<std::string> tokens;
  tokenize(temp, delm, tokens);
  return (double)tokens.size();
}

const char *do_get_word(double n, char *string, char *delm)
{
  size_t in = (size_t)n;
  std::string temp = string;
  std::vector<std::string> tokens;
  tokenize(temp, delm, tokens);

  if (tokens.size() >= in) {
    char *word = NULL;
    new_string(tokens[in-1].c_str(), &word);
    return word;
  } else {
    return NULL;
  }
}

const char *do_file_to_string(char *filename)
{
  std::fstream *file = aprepro->open_file(filename, "r");

  std::ostringstream lines;
  std::string line;
  while (std::getline(*file, line)) {
    lines << line << '\n';
  }

  char *ret_string;
  new_string(lines.str().c_str(), &ret_string);
  return ret_string;
}

const char *do_getenv(char *env)
{
  char *tmp;
  char *ret_string;
  if (env == NULL)
    return (NULL);
  
  tmp = (char *)getenv(env);
  if (tmp != NULL) {
    new_string(tmp, &ret_string);
    return (ret_string);
  } else {
    return (NULL);
  }
}

double do_strtod(char *string)
{
  double x;
  errno = 0;
  x = atof(string);
  SEAMS::math_error("strtod");
  return x;
}

const char *
do_dumpsym(void)
{
  aprepro->dumpsym(SEAMS::Parser::token::VAR, 0);
  return(NULL);
}

const char *
do_dumpfunc(void)
{
  aprepro->dumpsym(SEAMS::Parser::token::FNCT, 1);
  return(NULL);
}

const char *
do_dumpvar(void)
{
  aprepro->dumpsym(SEAMS::Parser::token::VAR, 1);
  return(NULL);
}

double do_option(char *option, double value)
{
  double current;
  current = -1;
  
  if (std::strcmp(option, "warning") == 0) {
    current = aprepro->ap_options.warning_msg;
    aprepro->ap_options.warning_msg = (value == 0.0) ? false : true;
  }

  else if (std::strcmp(option, "info") == 0) {
    current = aprepro->ap_options.info_msg;
    aprepro->ap_options.info_msg = (value == 0.0) ? false : true;
  }
  
  else if (std::strcmp(option, "debugging") == 0) {
    current = aprepro->ap_options.debugging;
    aprepro->ap_options.debugging = (value == 0.0) ? false : true;
  }
  
  else if (std::strcmp(option, "statistics") == 0) {
    aprepro->statistics();
  }

  else {
    fprintf(stderr, "Valid arguments to option are: 'warning', 'info', 'debugging', and 'statistics'\n");
  }
  return current;
}

const char *do_include_path(char *new_path)
{
  aprepro->ap_options.include_path = new_path;
  return (NULL);
}

const char *do_intout(double intval)
{
  /* convert 'intval' to a string using an integer format
   * This can be used if you need to set the default output
   * format to force the decimal point.  In that case, integers
   * also have a decimal point which is usually not wanted.
   * Using 'intout(val)', val will be converted to a string
   * using an integer format
   */
  
  char *tmp;
  static char tmpstr[128];
  if (intval == 0.0) {
    new_string("0", &tmp);
    return (tmp);
  }
  else {
    (void) sprintf(tmpstr, "%d", (int)intval);
    new_string(tmpstr, &tmp);
    return (tmp);
  }
}

const char *do_execute(char *string)
{
  aprepro->lexer->execute(string);
  return NULL;
}

const char *do_rescan(char *string)
{
  aprepro->lexer->rescan(string);
  return NULL;
}

const char *do_extract(char *string, char *begin, char *end)
{
  /* From 'string' return a substring delimited by 'begin' and 'end'.
   *  'begin' is included in the string, but 'end' is not. If
   *  'begin' does not appear in the string, return NULL; If 'end'
   *  does not appear, then return the remainder of the string. If
   *  'begin' == "", then start at beginning; if 'end' == "", then
   *  return remainder of the string.
   */
  
  char *start = string;
  char *tmp;
  int len = 0;
  
  if (std::strlen(begin) > 0) {
    start = std::strstr(string, begin);
    if (start == NULL)
      return NULL;
  }
  
  len = std::strlen(start);
  if (std::strlen(end) > 0) {
    char *finish = std::strstr(start, end);
    if (finish != NULL) {
      len = finish-start;
    }
  }

  char *tmpstr = new char[len+1];
  std::strncpy(tmpstr, start, len);
  tmpstr[len] = '\0';
  new_string(tmpstr, &tmp);
  delete[] tmpstr;
  return tmp;
}

const char *do_get_temp_filename()
{
  char *filename;
  new_string(get_temp_filename(), &filename);
  return filename;
}

const char *do_error (char *error_string)
{
  /* Print error message (to stderr) and exit */
  yyerror(*aprepro, error_string);
  exit(EXIT_FAILURE);
  /* NOTREACHED */
  return(NULL);
}
} // namespace SEAMS
