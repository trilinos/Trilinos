/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef SEAMS_UTIL_H
#define SEAMS_UTIL_H

#include "aprepro.h"
#include <string>

namespace SEAMS {
  bool        arg_check(SEAMS::symrec *symbol, bool is_null);
  void        conv_string(char *string);
  void        new_string(const std::string &from, char **to);
  void        new_string(const char *from, char **to);
  void        concat_string(const char *from1, const char *from2, char **to);
  const char *get_temp_filename();
  void        math_error(const SEAMS::Aprepro &aprepro, const char *function);
  void        math_error(const char *function);
  void        yyerror(const SEAMS::Aprepro &aprepro, const std::string &s);
  void        undefined_error(const SEAMS::Aprepro &aprepro, const std::string &var);
  void        redefined_warning(const SEAMS::Aprepro &aprepro, const SEAMS::symrec *var);
  void        warning(const SEAMS::Aprepro &aprepro, const std::string &var);
  void        immutable_modify(const SEAMS::Aprepro &aprepro, const SEAMS::symrec *var);
  void        set_type(const SEAMS::Aprepro &apr, SEAMS::symrec *var, int type);
  void        cleanup_memory();
  bool        is_directory(const std::string &filepath);
  bool        check_valid_var(const char *s);
} // namespace SEAMS
#endif
