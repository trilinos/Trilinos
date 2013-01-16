#ifndef SEAMS_UTIL_H
#define SEAMS_UTIL_H

#include <string>
#include "aprepro.h"

namespace SEAMS {
  void conv_string(char *string);
  void new_string(const char *from, char **to);
  void concat_string(const char *from1, const char *from2, char **to);
  char *get_temp_filename();
  void math_error(const SEAMS::Aprepro &aprepro, const char *function);
  void math_error(const char *function);
  void yyerror (const SEAMS::Aprepro &aprepro, const std::string &s);
  void undefined_warning (const SEAMS::Aprepro &aprepro, const std::string &var);
  void redefined_warning (const SEAMS::Aprepro &aprepro, const std::string &var);
  void warning (const SEAMS::Aprepro &aprepro, const std::string &var);
  void immutable_modify(const SEAMS::Aprepro &aprepro, const SEAMS::symrec* var);
  void set_type(const SEAMS::Aprepro &apr, SEAMS::symrec* var, int type);
  void cleanup_memory();
  bool is_directory(const std::string &filepath);
}
#endif
