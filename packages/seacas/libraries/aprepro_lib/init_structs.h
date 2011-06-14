#ifndef INIT_STRUCTS_H
#define INIT_STRUCTS_H
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

struct init_dddd
  {
    const char *fname;
    double (*fnct)(double, double, double, double);
    const char *syntax;
    const char *description;
  };

struct init_cc
  {
    const char *fname;
    double (*fnct)(char*, char*);
    const char *syntax;
    const char *description;
  };

struct init_c
  {
    const char *fname;
    double (*fnct)(char*);
    const char *syntax;
    const char *description;
  };

struct str_init
  {
    const char *fname;
    char *(*fnct)();
    const char *syntax;
    const char *description;
  };

struct str_c_init
  {
    const char *fname;
    char *(*fnct)(char*);
    const char *syntax;
    const char *description;
  };

struct str_d_init
  {
    const char *fname;
    char *(*fnct)(double);
    const char *syntax;
    const char *description;
  };

struct str_dcc_init
  {
    const char *fname;
    char *(*fnct)(double, char*, char*);
    const char *syntax;
    const char *description;
  };

struct str_ccc_init
  {
    const char *fname;
    char *(*fnct)(char*, char*, char*);
    const char *syntax;
    const char *description;
  };

struct var_init
  {
    const char *vname;
    double value;
  };

struct svar_init
  {
    const char *vname;
    const char *value;
  };

#endif
