// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "apr_builtin.h"

#include <cctype>
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include <sys/stat.h>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif
#include "apr_scanner.h"
#include "apr_tokenize.h"
#include "apr_util.h"
#include "aprepro.h"
#include "aprepro_parser.h"
#include <cassert>
#include <cstring>

#include <random>

#ifndef PI
#define PI 3.141592653589793238462643
#endif

namespace {
  std::unordered_map<size_t, std::vector<std::string>> tokenized_strings;

  std::vector<std::string> &get_tokenized_strings(const char *string, const char *delm)
  {
    // key is address of string + hash of delimiter
    size_t key = std::hash<std::string>{}(string) + std::hash<std::string>{}(delm);
    if (tokenized_strings.find(key) == tokenized_strings.end()) {
      std::string temp       = string;
      auto        tokens     = SEAMS::tokenize(temp, delm);
      tokenized_strings[key] = tokens;
    }
    return tokenized_strings[key];
  }

  std::random_device rd;
  std::mt19937_64    rng(rd());

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

  extern SEAMS::Aprepro *aprepro;

#define d2r(x) ((x)*PI / 180.)
#define r2d(x) ((x)*180. / PI)

#ifndef max
#define max(x, y) (x) > (y) ? (x) : (y)
#define min(x, y) (x) < (y) ? (x) : (y)
#endif

  double do_time()
  {
    time_t timer = time(nullptr);
    return timer;
  }

  // DO_INT:  Calculate integer nearest to zero from value
  double do_int(double x)
  {
    reset_error();
    double temp = (x < 0 ? -std::floor(-(x)) : std::floor(x));
    SEAMS::math_error("int");
    return (temp);
  }

  // DO_NINT:  Calculate integer nearest value
  double do_nint(double x)
  {
    reset_error();
    double temp = (x < 0 ? -std::floor(0.5 - x) : std::floor(x + 0.5));
    SEAMS::math_error("nint");
    return (temp);
  }

  // DO_DIST: Calculate distance between point 1 at (x1,y1) and
  //          point 2 at (x2,y2).

  double do_dist(double x1, double y1, double x2, double y2)
  {
    reset_error();
    double temp = std::hypot((x1 - x2), (y1 - y2));
    SEAMS::math_error("hypot");
    return (temp);
  }

  // DO_ANGLE: Calculate angle (radians) between vector 1 at (0,0; x1,y1) and
  //           vector 2 at (0,0; x2,y2).

  double do_angle(double x1, double y1, double x2, double y2)
  {
    double temp = ((x1 * x2) + (y1 * y2)) / (std::hypot(x1, y1) * std::hypot(x2, y2));
    reset_error();
    temp = acos(temp);
    SEAMS::math_error("angle");
    return (temp);
  }

  // DO_ANGLE: Calculate angle (degrees) between vector 1 at (0,0; x1,y1) and
  //           vector 2 at (0,0; x2,y2).

  double do_angled(double x1, double y1, double x2, double y2)
  {
    double temp = ((x1 * x2) + (y1 * y2)) / (std::hypot(x1, y1) * std::hypot(x2, y2));
    reset_error();
    temp = r2d(acos(temp));
    SEAMS::math_error("angled");
    return (temp);
  }

  // DO_HYPOT: calculate sqrt(p^2 + q^2)
  double do_hypot(double x, double y) { return std::hypot(x, y); }

  double do_pow(double x, double y)
  {
    // X^Y
    return std::pow(x, y);
  }

  double do_max(double x, double y)
  {
    reset_error();
    double temp = max(x, y);
    SEAMS::math_error("max");
    return (temp);
  }

  double do_min(double x, double y)
  {
    reset_error();
    double temp = min(x, y);
    SEAMS::math_error("min");
    return (temp);
  }

  double do_d2r(double x) { return (d2r(x)); }

  double do_r2d(double x) { return (r2d(x)); }

  double do_sind(double x)
  {
    reset_error();
    double temp = sin(d2r(x));
    SEAMS::math_error("sind");
    return (temp);
  }

  double do_sin(double x)
  {
    reset_error();
    double temp = sin(x);
    SEAMS::math_error("sin");
    return (temp);
  }

  double do_cosd(double x)
  {
    reset_error();
    double temp = cos(d2r(x));
    SEAMS::math_error("cosd");
    return (temp);
  }

  double do_cos(double x)
  {
    reset_error();
    double temp = cos(x);
    SEAMS::math_error("cos");
    return (temp);
  }

  double do_tand(double x)
  {
    reset_error();
    double temp = tan(d2r(x));
    SEAMS::math_error("tand");
    return (temp);
  }

  double do_tan(double x)
  {
    reset_error();
    double temp = tan(x);
    SEAMS::math_error("tan");
    return (temp);
  }

  double do_atan2d(double x, double y)
  {
    reset_error();
    double temp = r2d(atan2(x, y));
    SEAMS::math_error("atan2d");
    return (temp);
  }

  double do_atan2(double x, double y)
  {
    reset_error();
    double temp = atan2(x, y);
    SEAMS::math_error("atan2");
    return (temp);
  }

  double do_atand(double x)
  {
    reset_error();
    double temp = r2d(atan(x));
    SEAMS::math_error("atand");
    return (temp);
  }

  double do_atan(double x)
  {
    reset_error();
    double temp = atan(x);
    SEAMS::math_error("atan");
    return (temp);
  }

  double do_asind(double x)
  {
    reset_error();
    double temp = r2d(asin(x));
    SEAMS::math_error("asind");
    return (temp);
  }

  double do_asin(double x)
  {
    reset_error();
    double temp = asin(x);
    SEAMS::math_error("asin");
    return (temp);
  }

  double do_acosd(double x)
  {
    reset_error();
    double temp = r2d(acos(x));
    SEAMS::math_error("acosd");
    return (temp);
  }

  double do_acos(double x)
  {
    reset_error();
    double temp = acos(x);
    SEAMS::math_error("acos");
    return (temp);
  }

  // do_srand(x) Seed the random generator with the specified integer value
  double do_srand(double seed)
  {
    rng.seed(static_cast<size_t>(seed));
    return (0);
  }

  // do_rand(x) returns a random double in the range 0<= do_rand <= x
  double do_rand(double xl, double xh)
  {
    std::uniform_real_distribution<double> dist(xl, xh);
    return dist(rng);
  }

  double do_rand_normal(double mean, double stddev)
  {
    std::normal_distribution<> dist(mean, stddev);
    return dist(rng);
  }

  double do_rand_lognormal(double mean, double stddev)
  {
    std::lognormal_distribution<> dist(mean, stddev);
    return dist(rng);
  }

  double do_rand_weibull(double alpha, double beta)
  {
    std::weibull_distribution<> dist(alpha, beta);
    return dist(rng);
  }

  double do_sign(double x, double y)
  {
    reset_error();
    double temp = (y) >= 0 ? fabs(x) : -fabs(x);
    SEAMS::math_error("sign");
    return (temp);
  }

  double do_dim(double x, double y)
  {
    reset_error();
    double temp = x - (min(x, y));
    SEAMS::math_error("dim");
    return (temp);
  }

  double do_fabs(double x)
  {
    reset_error();
    double temp = fabs(x);
    SEAMS::math_error("fabs");
    return (temp);
  }

  double do_ceil(double x)
  {
    reset_error();
    double temp = ceil(x);
    SEAMS::math_error("ceil");
    return (temp);
  }

  double do_cosh(double x)
  {
    reset_error();
    double temp = cosh(x);
    SEAMS::math_error("cosh");
    return (temp);
  }

  double do_exp(double x)
  {
    reset_error();
    double temp = exp(x);
    SEAMS::math_error("exp");
    return (temp);
  }

  double do_expm1(double x)
  {
    reset_error();
    double temp = std::expm1(x);
    SEAMS::math_error("exp");
    return (temp);
  }

  double do_erf(double x) { return std::erf(x); }

  double do_erfc(double x) { return std::erfc(x); }

  double do_floor(double x)
  {
    reset_error();
    double temp = floor(x);
    SEAMS::math_error("floor");
    return (temp);
  }

  double do_fmod(double x, double y)
  {
    reset_error();
    double temp = fmod(x, y);
    SEAMS::math_error("fmod");
    return (temp);
  }

  double do_log(double x)
  {
    reset_error();
    double temp = std::log(x);
    SEAMS::math_error("log");
    return (temp);
  }

  double do_log10(double x)
  {
    reset_error();
    double temp = std::log10(x);
    SEAMS::math_error("log10");
    return (temp);
  }

  double do_sinh(double x)
  {
    reset_error();
    double temp = sinh(x);
    SEAMS::math_error("sinh");
    return (temp);
  }

  double do_sqrt(double x)
  {
    feclearexcept(FE_ALL_EXCEPT);
    reset_error();
    double temp = std::sqrt(x);
    if (fetestexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO) != 0) {
      SEAMS::math_error("sqrt");
    }
    return (temp);
  }

  double do_cbrt(double x)
  {
    feclearexcept(FE_ALL_EXCEPT);
    reset_error();
    double temp = std::cbrt(x);
    if (fetestexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO) != 0) {
      SEAMS::math_error("sqrt");
    }
    return (temp);
  }

  double do_tanh(double x)
  {
    reset_error();
    double temp = tanh(x);
    SEAMS::math_error("tanh");
    return (temp);
  }

  double do_polarX(double rad, double ang) { return (rad * cos(d2r(ang))); }

  double do_polarY(double rad, double ang) { return (rad * sin(d2r(ang))); }

  double do_lgamma(double val) { return std::lgamma(val); }

  double do_tgamma(double val) { return std::tgamma(val); }

  double do_juldayhms(double mon, double day, double year, double h, double mi, double se)
  {
    long   m = static_cast<long>(mon), d = static_cast<long>(day), y = static_cast<long>(year);
    long   c, ya, j;
    double seconds = h * 3600.0 + mi * 60 + se;

    if (m > 2) {
      m -= 3;
    }
    else {
      m += 9;
      --y;
    }
    c  = y / 100L;
    ya = y - (100L * c);
    j  = (146097L * c) / 4L + (1461L * ya) / 4L + (153L * m + 2L) / 5L + d + 1721119L;
    if (seconds < 12 * 3600.0) {
      j--;
      seconds += 12.0 * 3600.0;
    }
    else {
      seconds = seconds - 12.0 * 3600.0;
    }
    return (j + (seconds / 3600.0) / 24.0);
  }

  double do_julday(double mon, double day, double year)
  {
    return do_juldayhms(mon, day, year, 0.0, 0.0, 0.0);
  }

  double do_log1p(double x) { return std::log1p(x); }

  double do_acosh(double x) { return std::acosh(x); }

  double do_asinh(double x) { return std::asinh(x); }

  double do_atanh(double x) { return std::atanh(x); }

  double do_rows(const array *arr) { return arr->rows; }

  double do_cols(const array *arr) { return arr->cols; }

  // --------------------------STRING FUNCTIONS------------------------
  const char *do_get_date()
  {
    char *       tmp;
    const size_t bufsize = 32;
    static char  tmpstr[32];

    time_t     timer   = time(nullptr);
    struct tm *timeptr = localtime(&timer);

    // First  the date in the form CCYY/MM/DD
    strftime(tmpstr, bufsize, "%Y/%m/%d", timeptr);
    new_string(tmpstr, &tmp);
    return (tmp);
  }

  const char *do_get_iso_date()
  {
    char *       tmp;
    const size_t bufsize = 32;
    static char  tmpstr[32];

    time_t     timer   = time(nullptr);
    struct tm *timeptr = localtime(&timer);

    // First  the date in the form CCYY/MM/DD
    strftime(tmpstr, bufsize, "%Y%m%d", timeptr);
    new_string(tmpstr, &tmp);
    return (tmp);
  }

  const char *do_get_time()
  {
    char *       tmp;
    const size_t bufsize = 32;
    static char  tmpstr[32];

    time_t     timer   = time(nullptr);
    struct tm *timeptr = localtime(&timer);

    // Now the time in the form HH:MM:SS where 0 <= HH < 24
    strftime(tmpstr, bufsize, "%H:%M:%S", timeptr);
    new_string(tmpstr, &tmp);
    return (tmp);
  }

  const char *do_tolower(char *string)
  {
    char *p = string;
    while (*p != '\0') {
      if (isupper(static_cast<int>(*p)) != 0) {
        *p = tolower(static_cast<int>(*p));
      }
      p++;
    }
    return (string);
  }

  const char *do_toupper(char *string)
  {
    char *p = string;
    while (*p != '\0') {
      if (islower(static_cast<int>(*p)) != 0) {
        *p = toupper(static_cast<int>(*p));
      }
      p++;
    }
    return (string);
  }

  const char *do_tostring(double x)
  {
    char *      tmp;
    static char tmpstr[128];
    if (x == 0.0) {
      new_string("0", &tmp);
      return (tmp);
    }

    SEAMS::symrec *format;
    format = aprepro->getsym("_FORMAT");
    (void)sprintf(tmpstr, format->value.svar.c_str(), x);
    new_string(tmpstr, &tmp);
    return (tmp);
  }

  const char *do_output(char *filename)
  {
    aprepro->outputStream.top()->flush();

    if (std::strcmp(filename, "pop") == 0 || std::strcmp(filename, "stdout") == 0) {
      while (aprepro->outputStream.size() > 1) {
        std::ostream *output = aprepro->outputStream.top();
        aprepro->outputStream.pop();
        delete output;
      }

      aprepro->info("Output now redirected to original output stream.\n");
    }
    else {
      std::ostream *output = new std::ofstream(filename);
      aprepro->outputStream.push(output);
      aprepro->info("Output now redirected to file'" + std::string(filename) + "'.\n");
    }
    return (nullptr);
  }

  const char *do_append(char *filename)
  {
    aprepro->outputStream.top()->flush();

    if (std::strcmp(filename, "pop") == 0 || std::strcmp(filename, "stdout") == 0) {
      while (aprepro->outputStream.size() > 1) {
        std::ostream *output = aprepro->outputStream.top();
        aprepro->outputStream.pop();
        delete output;
      }

      aprepro->info("Output now redirected to original output stream.\n");
    }
    else {
      auto output = new std::ofstream(filename, std::ios_base::app); // Append
      aprepro->outputStream.push(output);

      aprepro->info("Output now redirected to file '" + std::string(filename) + "'\n");
    }
    return (nullptr);
  }

  double do_word_count(char *string, char *delm)
  {
    return static_cast<double>(get_tokenized_strings(string, delm).size());
  }

  double do_find_word(char *word, char *string, char *delm)
  {
    auto &      tokens = get_tokenized_strings(string, delm);
    std::string sword{word};
    for (size_t i = 0; i < tokens.size(); i++) {
      if (tokens[i] == sword) {
        return i + 1;
      }
    }
    return 0;
  }

  const char *do_get_word(double n, char *string, char *delm)
  {
    auto &tokens = get_tokenized_strings(string, delm);

    size_t in = static_cast<size_t>(n);
    if (tokens.size() >= in) {
      char *word = nullptr;
      new_string(tokens[in - 1], &word);
      return word;
    }
    return "";
  }

  const char *do_file_to_string(char *filename)
  {
    char *        ret_string = nullptr;
    std::fstream *file       = aprepro->open_file(filename, "r");

    if (file != nullptr) {
      std::ostringstream lines;
      std::string        line;
      while (std::getline(*file, line)) {
        lines << line << '\n';
      }

      new_string(lines.str(), &ret_string);
    }
    return ret_string;
  }

  const char *do_getenv(char *env)
  {
    if (env == nullptr) {
      return "";
    }
    char *tmp = getenv(env);
    if (tmp != nullptr) {
      char *ret_string;
      new_string(tmp, &ret_string);
      return (ret_string);
    }

    return "";
  }

  double do_strtod(char *string)
  {
    reset_error();
    double x = atof(string);
    SEAMS::math_error("strtod");
    return x;
  }

  const char *do_dumpsym()
  {
    aprepro->dumpsym(SEAMS::Parser::token::VAR, false);
    return (nullptr);
  }

  const char *do_dumpfunc()
  {
    aprepro->dumpsym(SEAMS::Parser::token::FNCT, true);
    return (nullptr);
  }

  const char *do_dumpvar()
  {
    aprepro->dumpsym(SEAMS::Parser::token::VAR, true);
    return (nullptr);
  }

  const char *do_dumpsym1(char *pre)
  {
    aprepro->dumpsym(SEAMS::Parser::token::VAR, pre, false);
    return (nullptr);
  }

  const char *do_dumpfunc1(char *pre)
  {
    aprepro->dumpsym(SEAMS::Parser::token::FNCT, pre, true);
    return (nullptr);
  }

  const char *do_dumpvar1(char *pre)
  {
    aprepro->dumpsym(SEAMS::Parser::token::VAR, pre, true);
    return (nullptr);
  }

  double do_option(char *option, double value)
  {
    double current = -1;

    if (std::strcmp(option, "warning") == 0) {
      current                         = static_cast<double>(aprepro->ap_options.warning_msg);
      aprepro->ap_options.warning_msg = value != 0.0;
    }

    else if (std::strcmp(option, "info") == 0) {
      current                      = static_cast<double>(aprepro->ap_options.info_msg);
      aprepro->ap_options.info_msg = value != 0.0;
    }

    else if (std::strcmp(option, "debugging") == 0) {
      current                       = static_cast<double>(aprepro->ap_options.debugging);
      aprepro->ap_options.debugging = value != 0.0;
    }

    else if (std::strcmp(option, "statistics") == 0) {
      aprepro->statistics();
    }

    else {
      aprepro->error(
          "Valid arguments to option are: 'warning', 'info', 'debugging', and 'statistics'\n",
          false);
    }
    return current;
  }

  const char *do_include_path(char *new_path)
  {
    aprepro->ap_options.include_path = new_path;
    return (nullptr);
  }

  const char *do_intout(double intval)
  {
    // convert 'intval' to a string using an integer format
    // This can be used if you need to set the default output
    // format to force the decimal point.  In that case, integers
    // also have a decimal point which is usually not wanted.
    // Using 'intout(val)', val will be converted to a string
    // using an integer format

    char *      tmp;
    static char tmpstr[128];
    if (intval == 0.0) {
      new_string("0", &tmp);
      return (tmp);
    }

    (void)sprintf(tmpstr, "%d", static_cast<int>(intval));
    new_string(tmpstr, &tmp);
    return (tmp);
  }

  const char *do_execute(char *string)
  {
    aprepro->lexer->execute(string);
    return nullptr;
  }

  const char *do_rescan(char *string)
  {
    aprepro->lexer->rescan(string);
    return nullptr;
  }

  const char *do_if(double x)
  {
    aprepro->inIfdefGetvar = false;
    aprepro->lexer->if_handler(x);
    return nullptr;
  }

  const char *do_notif(double x)
  {
    aprepro->lexer->if_handler(!x);
    return nullptr;
  }

  const char *do_elseif(double x)
  {
    aprepro->lexer->elseif_handler(x);
    return nullptr;
  }

  const char *do_str_if(char *string)
  {
    std::string test(string);
    aprepro->lexer->if_handler(static_cast<double>(!test.empty()));

    return nullptr;
  }

  const char *do_str_notif(char *string)
  {
    std::string test(string);
    aprepro->lexer->if_handler(static_cast<double>(test.empty()));

    return nullptr;
  }

  const char *do_str_elseif(char *string)
  {
    std::string test(string);
    aprepro->lexer->elseif_handler(static_cast<double>(!test.empty()));

    return nullptr;
  }

  const char *do_switch(double x)
  {
    aprepro->lexer->switch_handler(x);
    return nullptr;
  }

  const char *do_case(double x)
  {
    aprepro->lexer->case_handler(x);
    return nullptr;
  }

  const char *do_extract(char *string, char *begin, char *end)
  {
    // From 'string' return a substring delimited by 'begin' and 'end'.
    // 'begin' is included in the string, but 'end' is not. If
    // 'begin' does not appear in the string, return nullptr; If 'end'
    // does not appear, then return the remainder of the string. If
    // 'begin' == "", then start at beginning; if 'end' == "", then
    // return remainder of the string.
    char *start = string;

    if (std::strlen(begin) > 0) {
      start = std::strstr(string, begin);
      if (start == nullptr) {
        return "";
      }
    }

    int len = std::strlen(start);
    if (std::strlen(end) > 0) {
      char *finish = std::strstr(start, end);
      if (finish != nullptr) {
        len = finish - start;
      }
    }

    std::string tmpstr(start, 0, len);
    char *      tmp;
    new_string(tmpstr, &tmp);
    return tmp;
  }

  const char *do_get_temp_filename()
  {
    char *filename;
    new_string(get_temp_filename(), &filename);
    return filename;
  }

  const char *do_error(char *error_string)
  {
    // Print error message (to stderr) and exit
    yyerror(*aprepro, error_string);
    throw std::runtime_error(std::string(error_string));
    // NOTREACHED
    return (nullptr);
  }

  const char *do_print_array(const array *my_array_data)
  {
    if (my_array_data != nullptr) {
      std::ostringstream lines;

      int rows = my_array_data->rows;
      int cols = my_array_data->cols;

      int idx = 0;

      for (int ir = 0; ir < rows; ir++) {
        if (ir > 0) {
          lines << "\n";
        }
        lines << "\t";
        for (int ic = 0; ic < cols; ic++) {
          lines << my_array_data->data[idx++];
          if (ic < cols - 1) {
            lines << "\t";
          }
        }
      }
      char *ret_string;
      new_string(lines.str(), &ret_string);
      return ret_string;
    }

    return "";
  }

  const char *do_delete(char *string)
  {
    aprepro->remove_variable(string);

    return nullptr;
  }

  array *do_make_array(double rows, double cols)
  {
    auto array_data = new array(rows, cols);
    return array_data;
  }

  array *do_make_array_init(double rows, double cols, double init)
  {
    auto array_data = new array(rows, cols);
    int  isize      = (int)rows * int(cols);
    for (int i = 0; i < isize; i++) {
      array_data->data[i] = init;
    }
    return array_data;
  }

  array *do_identity(double size)
  {
    int  i;
    int  isize      = size;
    auto array_data = new array(size, size);

    for (i = 0; i < isize; i++) {
      array_data->data[i * isize + i] = 1.0;
    }
    return array_data;
  }

  array *do_linear_array(double init, double final, double count)
  {
    // Create 1D array with `count` rows and 1 column.
    // Values are linearly spaced from `init` to `final`
    int  isize      = count;
    auto array_data = new array(count, 1);

    double inc = (final - init) / (count - 1);
    for (int i = 0; i < isize; i++) {
      array_data->data[i] = init + (double)i * inc;
    }
    return array_data;
  }

  array *do_transpose(const array *a)
  {
    int  i, j;
    auto array_data = new array(a->cols, a->rows);

    for (i = 0; i < a->rows; i++) {
      for (j = 0; j < a->cols; j++) {
        array_data->data[j * a->rows + i] = a->data[i * a->cols + j];
      }
    }
    return array_data;
  }

  array *do_csv_array1(const char *filename) { return do_csv_array(filename, 0.0); }

  array *do_csv_array(const char *filename, double skip)
  {
    size_t rows_to_skip = static_cast<size_t>(skip);

    const char *  delim = ",\t ";
    std::fstream *file  = aprepro->open_file(filename, "r");
    if (file != nullptr) {

      size_t rows = 0;
      size_t cols = 0;

      std::string line;
      while (std::getline(*file, line)) {
        rows++;
        if (rows > rows_to_skip) {
          auto tokens = tokenize(line, delim);
          cols        = tokens.size() > cols ? tokens.size() : cols;
        }
      }

      auto array_data = new array(rows - rows_to_skip, cols);

      // Read file again storing entries in array_data->data
      file->clear();
      file->seekg(0);

      int idx = 0;
      rows    = 0;
      while (std::getline(*file, line)) {
        if (++rows > rows_to_skip) {
          auto tokens = tokenize(line, delim);
          for (size_t i = 0; i < static_cast<size_t>(array_data->cols); i++) {
            if (i < tokens.size()) {
              array_data->data[idx++] = std::stod(tokens[i]);
            }
            else {
              array_data->data[idx++] = 0.0;
            }
          }
        }
      }
      assert(rows - rows_to_skip == (size_t)array_data->rows);
      delete file;
      return array_data;
    }

    return nullptr;
  }

  array *do_csv_array2(const char *filename, const char *comment)
  {
    std::fstream *file = aprepro->open_file(filename, "r");
    if (file != nullptr) {

      size_t      rows  = 0;
      size_t      cols  = 0;
      const char *delim = ",\t ";

      std::string line;
      while (std::getline(*file, line)) {
        if (line[0] != comment[0]) {
          rows++;
          auto tokens = tokenize(line, delim);
          cols        = tokens.size() > cols ? tokens.size() : cols;
        }
      }

      auto array_data = new array(rows, cols);

      // Read file again storing entries in array_data->data
      file->clear();
      file->seekg(0);

      int idx = 0;
      rows    = 0;
      while (std::getline(*file, line)) {
        if (line[0] != comment[0]) {
          rows++;
          auto tokens = tokenize(line, delim);
          for (size_t i = 0; i < static_cast<size_t>(array_data->cols); i++) {
            if (i < tokens.size()) {
              array_data->data[idx++] = std::stod(tokens[i]);
            }
            else {
              array_data->data[idx++] = 0.0;
            }
          }
        }
      }
      assert((int)rows == array_data->rows);
      delete file;
      return array_data;
    }
    return nullptr;
  }

  array *do_array_from_string(const char *string, const char *delm)
  {
    auto tokens     = SEAMS::tokenize(string, delm);
    auto array_data = new array(tokens.size(), 1);

    int idx = 0;
    for (const auto &token : tokens) {
      array_data->data[idx++] = std::stod(token);
    }
    assert(idx == array_data->rows);
    return array_data;
  }
} // namespace SEAMS
