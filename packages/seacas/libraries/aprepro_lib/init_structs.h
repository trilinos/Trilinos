/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef INIT_STRUCTS_H
#define INIT_STRUCTS_H

namespace SEAMS {
  struct array;
} // namespace SEAMS

struct init
{
  const char *fname;
  double (*fnct)();
  const char *syntax;
  const char *description;
};

struct init_d
{
  const char *fname;
  double (*fnct)(double);
  const char *syntax;
  const char *description;
};

struct init_dd
{
  const char *fname;
  double (*fnct)(double, double);
  const char *syntax;
  const char *description;
};

struct init_cd
{
  const char *fname;
  double (*fnct)(char *, double);
  const char *syntax;
  const char *description;
};

struct init_ddd
{
  const char *fname;
  double (*fnct)(double, double, double);
  const char *syntax;
  const char *description;
};

struct init_dddd
{
  const char *fname;
  double (*fnct)(double, double, double, double);
  const char *syntax;
  const char *description;
};

struct init_dddddd
{
  const char *fname;
  double (*fnct)(double, double, double, double, double, double);
  const char *syntax;
  const char *description;
};

struct init_ccc
{
  const char *fname;
  double (*fnct)(char *, char *, char *);
  const char *syntax;
  const char *description;
};

struct init_cc
{
  const char *fname;
  double (*fnct)(char *, char *);
  const char *syntax;
  const char *description;
};

struct init_c
{
  const char *fname;
  double (*fnct)(char *);
  const char *syntax;
  const char *description;
};

struct init_a
{
  const char *fname;
  double (*fnct)(const SEAMS::array *);
  const char *syntax;
  const char *description;
};

struct str_init
{
  const char *fname;
  const char *(*fnct)();
  const char *syntax;
  const char *description;
};

struct str_c_init
{
  const char *fname;
  const char *(*fnct)(char *);
  const char *syntax;
  const char *description;
};

struct str_d_init
{
  const char *fname;
  const char *(*fnct)(double);
  const char *syntax;
  const char *description;
};

struct str_a_init
{
  const char *fname;
  const char *(*fnct)(const SEAMS::array *);
  const char *syntax;
  const char *description;
};

struct str_dcc_init
{
  const char *fname;
  const char *(*fnct)(double, char *, char *);
  const char *syntax;
  const char *description;
};

struct str_cc_init
{
  const char *fname;
  const char *(*fnct)(char *, char *);
  const char *syntax;
  const char *description;
};

struct str_ccc_init
{
  const char *fname;
  const char *(*fnct)(char *, char *, char *);
  const char *syntax;
  const char *description;
};

struct array_c_init
{
  const char *fname;
  SEAMS::array *(*fnct)(const char *);
  const char *syntax;
  const char *description;
};

struct array_cc_init
{
  const char *fname;
  SEAMS::array *(*fnct)(const char *, const char *);
  const char *syntax;
  const char *description;
};

struct array_cd_init
{
  const char *fname;
  SEAMS::array *(*fnct)(const char *, double);
  const char *syntax;
  const char *description;
};

struct array_ddd_init
{
  const char *fname;
  SEAMS::array *(*fnct)(double, double, double);
  const char *syntax;
  const char *description;
};

struct array_dd_init
{
  const char *fname;
  SEAMS::array *(*fnct)(double, double);
  const char *syntax;
  const char *description;
};

struct array_d_init
{
  const char *fname;
  SEAMS::array *(*fnct)(double);
  const char *syntax;
  const char *description;
};

struct array_a_init
{
  const char *fname;
  SEAMS::array *(*fnct)(const SEAMS::array *);
  const char *syntax;
  const char *description;
};

struct var_init
{
  const char *vname;
  double      value;
};

struct svar_init
{
  const char *vname;
  const char *value;
};

#endif
