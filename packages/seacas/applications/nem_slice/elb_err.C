/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "elb_err.h"
#include <cstddef> // for size_t
#include <cstdio>
#include <fmt/ostream.h>
#include <vector> // for vector

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
          fmt::print(stderr, "================================");
          fmt::print(stderr, "messages");
          fmt::print(stderr, "================================\n");
          iflag = 1;
        }

        fmt::print(stderr, "\t{}\n", error_info[i].err_mesg);
        if (error_lev >= 2) {
          fmt::print(stderr, "\t\tin file {}\n", error_info[i].filename);
        }

        if (error_lev >= 3) {
          fmt::print(stderr, "\t\t\tat line {}\n", error_info[i].line_no);
        }
      }
    }
  }
}
