/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
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
 *     * Neither the name of NTESS nor the names of its
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

#include "elb_err.h"
#include <cstddef> // for size_t
#include <cstdio>  // for fprintf, stderr
#include <vector>  // for vector

const int MAX_ERR_MSG = 1024;
int       error_lev   = 1;

static std::vector<error_message> error_info;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function adds the specified error message to the array of error
 * structures.
 *
 * A level 0 error indicates a fatal error, otherwise it's a warning.
 *****************************************************************************/
void error_add(int level, const std::string &message, const std::string &filename, int line_no)
{
  if (error_info.size() < static_cast<size_t>(MAX_ERR_MSG)) {
    error_info.emplace_back(level, message, line_no, filename);
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function outputs the accumulated error messages to stderr
 *****************************************************************************/
void error_report()
{
  int iflag = 0;

  if (error_lev >= 1) {
    size_t error_cnt = error_info.size();
    for (size_t i = 0; i < error_cnt; i++) {
      if (error_info[i].level == 0 || error_info[i].level >= error_lev) {
        if (iflag == 0) {
          fprintf(stderr, "================================");
          fprintf(stderr, "messages");
          fprintf(stderr, "================================\n");
          iflag = 1;
        }

        fprintf(stderr, "\t%s\n", error_info[i].err_mesg.c_str());
        if (error_lev >= 2) {
          fprintf(stderr, "\t\tin file %s\n", error_info[i].filename.c_str());
        }

        if (error_lev >= 3) {
          fprintf(stderr, "\t\t\tat line %d\n", error_info[i].line_no);
        }
      }
    }
  }
  return;
}
