/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef EXODIFF_UTIL_H
#define EXODIFF_UTIL_H

#include <fmt/color.h>
#include <sstream>
#include <string>

char **get_name_array(int size, int length);
void   free_name_array(char **names, int size);

void Error(const std::string &x);
void ERR_OUT(std::ostringstream &buf);
void DIFF_OUT(std::ostringstream &buf, fmt::internal::color_type color = fmt::color::red);
void DIFF_OUT(const std::string &buf, fmt::internal::color_type color = fmt::color::red);

#endif
