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

#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include "my_aprepro.h"
#include "getline.h"
#include "y.tab.h"

#ifndef RAND_MAX
#include <limits.h>
#define RAND_MAX INT_MAX
#endif

#ifndef PI
#define PI  3.141592653589793238462643
#endif

#if defined(VMS) || defined(_hpux_) || defined(sun)
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

#if defined(sun)
#define LOG1P(x)	log1p(x)
#else
#define LOG1P(x)	log(1.0 + (x))
#endif

extern aprepro_options ap_options;
extern FILE *open_file(char *file, char *mode);

void check_math_error(double);
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
double do_srand(double seed);
double do_rand(double xl, double xh);
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

char  *do_getenv(char *string);
char  *do_tolower(char *string);
char  *do_toupper(char *string);
char  *do_tostring(double x);
char  *do_output(char *newfile);
char  *do_append(char *newfile);
char  *do_error(char *error_string);
char  *do_help(void);
char  *do_dumpsym(void);
char  *do_get_date(void);
char  *do_get_iso_date(void);
char  *do_get_time(void);
double do_option(char *option, double value);
extern void dumpsym(int type, int doInternal);  /* in hash.c */
double do_word_count(char *string, char *delm );
char  *do_get_word(double n, char *string, char *delm);
char  *do_file_to_string(char *filename);
char  *do_extract(char *string, char *begin, char *end);
double do_Material(double id, char *type, char *name, char *model, char *code, FILE * yyout);
double do_lgamma(double val);
double do_juldayhms(double mon, double day, double year,
			    double h, double mi, double se);
double do_julday(double mon, double day, double year);
double do_log1p(double mag);
char  *do_include_path(char *newpath);
char  *do_intout(double intval);

/* DO_INT:  Calculate integer nearest to zero from value */
double do_int(double x)
{
  double temp;
  errno = 0;
  temp = (double) (x < 0 ? -floor (-(x)) : floor (x));
  MATH_ERROR ("int");
  return (temp);
}

/* DO_NINT:  Calculate integer nearest value */
double do_nint(double x)
{
  double temp;
  errno = 0;
  temp = (double) (x < 0 ? -floor(0.5-x) : floor(x+0.5));
  MATH_ERROR ("nint");
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
  MATH_ERROR("hypot");
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
  MATH_ERROR("angle");
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
  MATH_ERROR("angled");
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
  MATH_ERROR("max");
  return (temp);
}

double 
do_min(double x, double y)
{
  double temp;
  errno = 0;
  temp = min(x, y);
  MATH_ERROR("min");
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
  MATH_ERROR("sind");
  return (temp);
}

double 
do_sin(double x)
{
  double temp;
  errno = 0;
  temp = sin(x);
  MATH_ERROR("sin");
  return (temp);
}

double 
do_cosd(double x)
{
  double temp;
  errno = 0;
  temp = cos(d2r(x));
  MATH_ERROR("cosd");
  return (temp);
}

double 
do_cos(double x)
{
  double temp;
  errno = 0;
  temp = cos(x);
  MATH_ERROR("cos");
  return (temp);
}

double 
do_tand(double x)
{
  double temp;
  errno = 0;
  temp = tan(d2r(x));
  MATH_ERROR("tand");
  return (temp);
}

double 
do_tan(double x)
{
  double temp;
  errno = 0;
  temp = tan(x);
  MATH_ERROR("tan");
  return (temp);
}

double 
do_atan2d(double x, double y)
{
  double temp;
  errno = 0;
  temp = r2d(atan2(x, y));
  MATH_ERROR("atan2d");
  return (temp);
}

double 
do_atan2(double x, double y)
{
  double temp;
  errno = 0;
  temp = atan2(x, y);
  MATH_ERROR("atan2");
  return (temp);
}

double 
do_atand(double x)
{
  double temp;
  errno = 0;
  temp = r2d(atan(x));
  MATH_ERROR("atand");
  return (temp);
}

double 
do_atan(double x)
{
  double temp;
  errno = 0;
  temp = atan(x);
  MATH_ERROR("atan");
  return (temp);
}

double 
do_asind(double x)
{
  double temp;
  errno = 0;
  temp = r2d(asin(x));
  MATH_ERROR("asind");
  return (temp);
}

double 
do_asin(double x)
{
  double temp;
  errno = 0;
  temp = asin(x);
  MATH_ERROR("asin");
  return (temp);
}

double 
do_acosd(double x)
{
  double temp;
  errno = 0;
  temp = r2d(acos(x));
  MATH_ERROR("acosd");
  return (temp);
}

double 
do_acos(double x)
{
  double temp;
  errno = 0;
  temp = acos(x);
  MATH_ERROR("acos");
  return (temp);
}

/* do_rand(x) returns a random double in the range 0<= do_rand <= x */
double 
do_rand(double xl, double xh)
{
  double temp;
  errno = 0;
  temp = xl + (xh - xl) * ((double) rand() / (double) RAND_MAX);
  MATH_ERROR("rand");
  return (temp);
}

/* do_srand(x) Seed the random generator with the specified integer value */
double 
do_srand(double seed)
{
  srand((unsigned)seed);
  return (0);
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

    w = sqrt( (-2.0 * log(w)) / w);
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
  double logstd = log(1.0 + (stddev/mean)*(stddev/mean));
  double logmean = log(mean) - 0.5 * logstd;
  logstd = sqrt(logstd);

  x = do_rand_normal(logmean, logstd);

  return exp(x);
}

double 
do_rand_weibull(double alpha, double beta)
{
  double temp = (double) rand() / (double) RAND_MAX;
  errno = 0;
  temp = pow( (-1.0 / alpha * log(1.0 - temp)), (1.0/beta) );
  MATH_ERROR("weibull");
  return (temp);
}

double 
do_sign(double x, double y)
{
  double temp;
  errno = 0;
  temp = (y) >= 0 ? fabs(x) : -fabs(x);
  MATH_ERROR("sign");
  return (temp);
}

double 
do_dim(double x, double y)
{
  double temp;
  errno = 0;
  temp = x - (min(x, y));
  MATH_ERROR("dim");
  return (temp);
}

double 
do_fabs(double x)
{
  double temp;
  errno = 0;
  temp = fabs(x);
  MATH_ERROR("fabs");
  return (temp);
}

double 
do_ceil(double x)
{
  double temp;
  errno = 0;
  temp = ceil(x);
  MATH_ERROR("ceil");
  return (temp);
}

double 
do_cosh(double x)
{
  double temp;
  errno = 0;
  temp = cosh(x);
  MATH_ERROR("cosh");
  return (temp);
}

double 
do_exp(double x)
{
  double temp;
  errno = 0;
  temp = exp(x);
  MATH_ERROR("exp");
  return (temp);
}

double 
do_floor(double x)
{
  double temp;
  errno = 0;
  temp = floor(x);
  MATH_ERROR("floor");
  return (temp);
}

double 
do_fmod(double x, double y)
{
  double temp;
  errno = 0;
  temp = fmod(x, y);
  MATH_ERROR("fmod");
  return (temp);
}

double 
do_log(double x)
{
  double temp;
  errno = 0;
  temp = log(x);
  MATH_ERROR("log");
  return (temp);
}

double 
do_log10(double x)
{
  double temp;
  errno = 0;
  temp = log10(x);
  MATH_ERROR("log10");
  return (temp);
}

double 
do_sinh(double x)
{
  double temp;
  errno = 0;
  temp = sinh(x);
  MATH_ERROR("sinh");
  return (temp);
}

double 
do_sqrt(double x)
{
  double temp;
  errno = 0;
  temp = sqrt(x);
  MATH_ERROR("sqrt");
  return (temp);
}

double 
do_tanh(double x)
{
  double temp;
  errno = 0;
  temp = tanh(x);
  MATH_ERROR("tanh");
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
  tmp = (x + 0.5) * log(tmp) - tmp;
  ser = 1.0;
  for (j = 0; j < 6; j++)
    {
      x += 1.0;
      ser += (cof[j] / x);
    }
  return (tmp + log(STP * ser));
}

double do_juldayhms(double mon, double day, double year,
	      double h, double mi, double se)
{
  long m = mon, d = day, y = year;
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
    return (LOG1P(x) + log(2.0));
  else
    {
      t = sqrt(x - 1.0);
      return (LOG1P(t * (t + sqrt(x + 1))));
    }
}

double do_asinh(double x)
{
  double s, t;
  if (1.0 + x * x == 1.0)
    return (x);
  if (sqrt(1.0 + x * x) == 1.0)
    return (do_sign(1.0, x) * (LOG1P(x) + log(2.0)));
  else
    {
      t = fabs(x);
      s = 1.0 / t;
      return (do_sign(1.0, x) * LOG1P(t + t / (s + sqrt(1 + s * s))));
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

char *do_get_date(void)
{
  char *tmp;
  const size_t bufsize = 32;
  static char tmpstr[32];

  time_t timer = time(NULL);
  struct tm *timeptr = localtime(&timer);
  
  /* First  the date in the form CCYY/MM/DD */
  strftime(tmpstr, bufsize, "%Y/%m/%d", timeptr);
  NEWSTR(tmpstr, tmp);
  return(tmp);
}

char *do_get_iso_date(void)
{
  char *tmp;
  const size_t bufsize = 32;
  static char tmpstr[32];

  time_t timer = time(NULL);
  struct tm *timeptr = localtime(&timer);
  
  /* First  the date in the form CCYY/MM/DD */
  strftime(tmpstr, bufsize, "%Y%m%d", timeptr);
  NEWSTR(tmpstr, tmp);
  return(tmp);
}

char *do_get_time(void)
{
  char *tmp;
  const size_t bufsize = 32;
  static char tmpstr[32];

  time_t timer = time(NULL);
  struct tm *timeptr = localtime(&timer);
  
  /* Now the time in the form HH:MM:SS where 0 <= HH < 24 */
  strftime(tmpstr, bufsize, "%H:%M:%S", timeptr);
  NEWSTR(tmpstr, tmp);
  return(tmp);
}

char *do_tolower(char *string)
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

char *do_toupper(char *string)
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

char *do_tostring(double x)
{
  char *tmp;
  static char tmpstr[128];
  if (x == 0.0)
    {
      NEWSTR("0", tmp);
      return (tmp);
    }
  else
    {
      symrec *format;
      format = getsym("_FORMAT");
      (void) sprintf(tmpstr, format->value.svar, x);
      NEWSTR(tmpstr, tmp);
      return (tmp);
    }
}

char *do_output(char *filename)
{
  errno = 0;
  fflush(yyout);
  MATH_ERROR("output (fflush)");
  if (yyout != stdout) 
    fclose(yyout);
  MATH_ERROR("output (fclose)");
  if (strcmp(filename, "stdout") == 0)
    yyout = stdout;
  else
    yyout = open_file(filename, "w");
  return (NULL);
}

char *do_append(char *filename)
{
  errno = 0;
  fflush(yyout);
  MATH_ERROR("append (fflush)");
  if (yyout != stdout) 
    fclose(yyout);
  MATH_ERROR("append (fclose)");
  if (strcmp(filename, "stdout") == 0)
    yyout = stdout;
  else
    yyout = open_file(filename, "a");
  return (NULL);
}

double do_word_count(char *string, char *delm)
{
   char *temp ;
   double i = 0;

  NEWSTR(string, temp);
  
  if( strtok(temp,delm)) {
      i++;
      while( strtok(NULL,delm) ) {
      i++;
        }
      }
  free(temp);
  return (i) ;
}

char *do_get_word(double n, char *string, char *delm)
{
   char *temp, *token, *word;
   int i;

    NEWSTR(string, temp);
    token = strtok(temp,delm);
    if( n == 1 )
     {
      NEWSTR(token,word);
      free(temp);
      return(word);
     }
    for(i=1; i<n; i++)
	{
        if( (token = strtok(NULL,delm)) == NULL )
	    {
	        free(temp);
            return(NULL);
            }
	 }
     NEWSTR(token,word);
     free(temp);
     return(word);
}

char *do_file_to_string(char *filename)
{
  FILE * fp;
  int size = 0;
  char *line = NULL;
  char *ret_string = NULL; 
  size_t len = 0;
  int error = 0;
  struct stat st;
  
  char *lines = NULL;
  errno = 0;

  error = stat(filename, &st);
  if (error < 0) {
    char tmpstr[128];
    sprintf(tmpstr, "Aprepro: ERR:  Can't open '%s'",filename); 
    perror(tmpstr);
    exit(EXIT_FAILURE);
  }

  /* Add extra characters to size in case file does not end in newline
   * getline add the newline at the end which results in overrunning
   * the memory.
   */
  size = st.st_size+2;

  lines = malloc(size * sizeof(char)+1);
  lines[0] = '\0';
  
  fp = open_file(filename, "r");

  while (getline(&line, &len, fp) != -1) {
    strcat ( lines, line );
    assert(strlen(lines) <= size);
  }

  assert(strlen(lines) <= size);
  NEWSTR(lines, ret_string);
  if (line) free(line);
  if (lines) free(lines);
  return ret_string;
}

char *do_getenv(char *env)
{
  char *tmp;
  char *ret_string;
  if (env == NULL)
    return (NULL);
  
  tmp = (char *)getenv(env);
  if (tmp != NULL) {
    NEWSTR(tmp, ret_string);
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
  MATH_ERROR("strtod");
  return x;
}

char *
do_help(void)
{
  char comment = getsym("_C_")->value.svar[0];
  printf ("\n%c   Enter {DUMP()}        to list defined variables\n", comment);
  printf ("%c         {DUMP_FUNC()}   to list of all double and string functions\n", comment);
  printf ("%c         {DUMP_PREVAR()} to list all predefined variables\n", comment);
  return("");
}

char *
do_dumpsym(void)
{
  dumpsym(VAR, 0);
  return("");
}

char *
do_dumpfunc(void)
{
  dumpsym(FNCT, 1);
  return("");
}

char *
do_dumpvar(void)
{
  dumpsym(VAR, 1);
  return("");
}

double do_option(char *option, double value)
{
  double current;
  current = -1;
  
  if (strcmp(option, "warning") == 0) {
    current = ap_options.warning_msg;
    ap_options.warning_msg = (value == 0.0) ? False : True;
  }

  else if (strcmp(option, "info") == 0) {
    current = ap_options.info_msg;
    ap_options.info_msg = (value == 0.0) ? False : True;
  }
  
  else if (strcmp(option, "debugging") == 0) {
    current = ap_options.debugging;
    ap_options.debugging = (value == 0.0) ? False : True;
  }
  
  else if (strcmp(option, "statistics") == 0) {
    current = ap_options.statistics;
    ap_options.statistics = (value == 0.0) ? False : True;
  }

  else {
    fprintf(stderr, "Valid arguments to option are: 'warning', 'info', 'debugging', and 'statistics'\n");
  }
  return current;
}

char *do_include_path(char *new_path)
{
  if (ap_options.include_path != NULL) {
    free(ap_options.include_path);
  }
  NEWSTR(new_path, ap_options.include_path);
  return (NULL);
}

char *do_intout(double intval)
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
  if (intval == 0.0)
    {
      NEWSTR("0", tmp);
      return (tmp);
    }
  else
    {
      (void) sprintf(tmpstr, "%d", (int)intval);
      NEWSTR(tmpstr, tmp);
      return (tmp);
    }
}

char *do_extract(char *string, char *begin, char *end)
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
  
  if (strlen(begin) > 0) {
    start = strstr(string, begin);
    if (start == NULL)
      return "";
  }
  
  len = strlen(start);
  if (strlen(end) > 0) {
    char *finish = strstr(start, end);
    if (finish != NULL) {
      len = finish-start;
    }
  }

  {
    char *tmpstr = malloc(len+1);
    strncpy(tmpstr, start, len);
    tmpstr[len] = '\0';
    NEWSTR(tmpstr, tmp);
    free(tmpstr);
  }      
  return tmp;
}
