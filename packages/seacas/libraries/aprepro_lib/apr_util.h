/*
 * Copyright (c) 2014, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Governement retains certain rights in this software.
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
#ifndef SEAMS_UTIL_H
#define SEAMS_UTIL_H

#include <string>
#include "aprepro.h"

namespace SEAMS {
  bool arg_check(SEAMS::symrec* symbol, bool is_null);
  void conv_string(char *string);
  void new_string(const char *from, char **to);
  void concat_string(const char *from1, const char *from2, char **to);
  char *get_temp_filename();
  void math_error(const SEAMS::Aprepro &aprepro, const char *function);
  void math_error(const char *function);
  void yyerror (const SEAMS::Aprepro &aprepro, const std::string &s);
  void undefined_warning (const SEAMS::Aprepro &aprepro, const std::string &var);
  void redefined_warning (const SEAMS::Aprepro &aprepro, const SEAMS::symrec* var);
  void warning (const SEAMS::Aprepro &aprepro, const std::string &var);
  void immutable_modify(const SEAMS::Aprepro &aprepro, const SEAMS::symrec* var);
  void set_type(const SEAMS::Aprepro &apr, SEAMS::symrec* var, int type);
  void cleanup_memory();
  bool is_directory(const std::string &filepath);
  bool check_valid_var(const char *s);
}
#endif
