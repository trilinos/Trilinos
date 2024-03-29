// Copyright(C) 1999-, 20212021, , ,  National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

// Might be good to add a callback function which would be called
// when there was output -- In LexerOutput for example.  Default
// could be to just write to std::cout or to resultsOutput stringstream...
#pragma once

#include <string>
#include <vector>

namespace SEAMS {
  struct array
  {
    std::vector<double> data{};
    int                 rows{0};
    int                 cols{0};

    array(int r, int c) : rows(r), cols(c) { data.resize(r * c); }
  };

  struct symrec
  {
    std::string name{};
    std::string syntax{};
    std::string info{};
    int         type;
    bool        isInternal;
    struct value
    {
      double var{0};
      double (*fnctptr)(){nullptr};
      double (*fnctptr_d)(double){nullptr};
      double (*fnctptr_c)(char *){nullptr};
      double (*fnctptr_dc)(double, char *){nullptr};
      double (*fnctptr_cd)(char *, double){nullptr};
      double (*fnctptr_cc)(char *, char *){nullptr};
      double (*fnctptr_dd)(double, double){nullptr};
      double (*fnctptr_ddd)(double, double, double){nullptr};
      double (*fnctptr_ccc)(char *, char *, char *){nullptr};
      double (*fnctptr_ccd)(char *, char *, double){nullptr};
      double (*fnctptr_dddd)(double, double, double, double){nullptr};
      double (*fnctptr_ddddc)(double, double, double, double, char *){nullptr};
      double (*fnctptr_dddddd)(double, double, double, double, double, double){nullptr};
      double (*fnctptr_a)(const array *){nullptr};
      std::string svar{};
      const char *(*strfnct)(){nullptr};
      const char *(*strfnct_c)(char *){nullptr};
      const char *(*strfnct_d)(double){nullptr};
      const char *(*strfnct_a)(const array *){nullptr};
      const char *(*strfnct_dd)(double, double){nullptr};
      const char *(*strfnct_cc)(char *, char *){nullptr};
      const char *(*strfnct_ccc)(char *, char *, char *){nullptr};
      const char *(*strfnct_dcc)(double, char *, char *){nullptr};
      const char *(*strfnct_dcccc)(double, char *, char *, char *, char *){nullptr};
      array *avar{nullptr}; /* Array Variable */
      array *(*arrfnct_c)(const char *){nullptr};
      array *(*arrfnct_cc)(const char *, const char *){nullptr};
      array *(*arrfnct_cd)(const char *, double){nullptr};
      array *(*arrfnct_ddd)(double, double, double){nullptr};
      array *(*arrfnct_dd)(double, double){nullptr};
      array *(*arrfnct_d)(double){nullptr};
      array *(*arrfnct_a)(const array *){nullptr};
    } value;
    symrec *next;

    symrec(const char *my_name, int my_type, bool is_internal = false)
        : name(my_name), syntax(my_name), info("UNDEFINED"), type(my_type), isInternal(is_internal),
          next(nullptr)
    {
      value.var = 0;
    }

    symrec(const std::string &my_name, int my_type, bool is_internal = false)
        : name(my_name), syntax(my_name), info("UNDEFINED"), type(my_type), isInternal(is_internal),
          next(nullptr)
    {
      value.var = 0;
    }
  };
} // namespace SEAMS
