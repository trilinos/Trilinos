/*
 * Copyright (c) 2014, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
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
#ifndef INIT_STRUCTS_H
#define INIT_STRUCTS_H

namespace SEAMS {
  struct array;
}

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
  double (*fnct)(char*, double);
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

struct init_a
  {
    const char *fname;
    double (*fnct)(const SEAMS::array*);
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
    const char *(*fnct)(char*);
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
    const char *(*fnct)(const SEAMS::array*);
    const char *syntax;
    const char *description;
  };

struct str_dcc_init
  {
    const char *fname;
    const char *(*fnct)(double, char*, char*);
    const char *syntax;
    const char *description;
  };

struct str_ccc_init
  {
    const char *fname;
    const char *(*fnct)(char*, char*, char*);
    const char *syntax;
    const char *description;
  };

struct array_c_init
  {
    const char *fname;
    SEAMS::array *(*fnct)(const char*);
    const char *syntax;
    const char *description;
  };

struct array_cc_init
  {
    const char *fname;
    SEAMS::array *(*fnct)(const char*, const char*);
    const char *syntax;
    const char *description;
  };

struct array_cd_init
  {
    const char *fname;
    SEAMS::array *(*fnct)(const char*, double);
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
    SEAMS::array *(*fnct)(const SEAMS::array*);
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
