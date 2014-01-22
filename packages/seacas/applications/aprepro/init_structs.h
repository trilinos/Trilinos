#ifndef INIT_STRUCTS_H
#define INIT_STRUCTS_H
#include "my_aprepro.h"

struct init
  {
    char *fname;
    double (*fnct)();
    char *syntax;
    char *description;
  };

struct str_init
  {
    char *fname;
    char *(*fnct)();
    char *syntax;
    char *description;
  };

struct array_init
  {
    char *fname;
    array *(*fnct)();
    char *syntax;
    char *description;
  };

struct var_init
  {
    char *vname;
    double value;
  };

struct svar_init
  {
    char *vname;
    char *value;
  };

#endif
