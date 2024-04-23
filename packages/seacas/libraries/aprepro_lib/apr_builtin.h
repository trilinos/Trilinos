/*
 * Copyright(C) 1999-2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include <cstdio>

namespace SEAMS {
  struct array;

  double do_FtoC(double F);
  double do_CtoF(double C);
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
  double do_cbrt(double x);
  double do_cos(double x);
  double do_cosd(double x);
  double do_cosh(double x);
  double do_d2r(double x);
  double do_dim(double x, double y);
  double do_dist(double x1, double y1, double x2, double y2);
  double do_erf(double x);
  double do_erfc(double x);
  double do_exp(double x);
  double do_expm1(double x);
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
  double do_pow(double x, double y);
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
  double do_option(char *option, double value);
  double do_word_count(char *string, char *delm);
  double do_find_word(char *word, char *string, char *delm);
  double do_lgamma(double val);
  double do_tgamma(double val);
  double do_juldayhms(double mon, double day, double year, double h, double mi, double se);
  double do_julday(double mon, double day, double year);
  double do_log1p(double x);
  double do_rows(const array *arr);
  double do_cols(const array *arr);
  double do_time();

  const char *do_dumpsym();
  const char *do_dumpsym_json();
  const char *do_dumpfunc();
  const char *do_dumpvar();
  const char *do_get_date();
  const char *do_get_iso_date();
  const char *do_get_time();
  const char *do_get_temp_filename();
  const char *do_version();

  const char *do_dumpsym1(char *pre);
  const char *do_dumpfunc1(char *pre);
  const char *do_dumpvar1(char *pre);
  const char *do_Units(char *type);
  const char *do_file_to_string(char *filename);
  const char *do_error(char *error_string);
  const char *do_include_path(char *new_path);
  const char *do_getenv(char *env);
  const char *do_output(char *filename);
  const char *do_append(char *filename);
  const char *do_rescan(char *string);

  const char *do_if(double x);
  const char *do_notif(double x);
  const char *do_elseif(double x);
  const char *do_switch(double x);
  const char *do_case(double x);

  const char *do_get_word(double n, char *string, char *delm);
  const char *do_print_array(const array *my_array_data);

  const char *do_import(char *string);
  const char *do_execute(char *string);
  const char *do_tolower(char *string);
  const char *do_toupper(char *string);
  const char *do_tostring(double x);
  const char *do_extract(char *string, char *begin, char *end);
  const char *do_intout(double intval);
  const char *do_print_array(array *my_array_data);
  const char *do_str_if(char *string);
  const char *do_str_notif(char *string);
  const char *do_str_elseif(char *string);
  const char *do_delete(char *string);

#if defined(EXODUS_SUPPORT)
  const char *do_exodus_info_range(char *filename, char *beg, char *end);
  const char *do_exodus_info(char *filename, char *prefix);
  const char *do_exodus_meta(char *filename);
#endif

  array *do_csv_array(const char *filename, double skip);
  array *do_csv_array1(const char *filename);
  array *do_csv_array2(const char *filename, const char *comment);
  array *do_make_array(double rows, double cols);
  array *do_make_array_init(double rows, double cols, double init);
  array *do_linear_array(double init, double final, double count);
  array *do_identity(double size);
  array *do_transpose(const array *a);
  array *do_principal(const array *a);
  array *do_array_from_string(const char *string, const char *delm);
  array *do_sym_tensor_from_string(const char *string, const char *delm);
} // namespace SEAMS
