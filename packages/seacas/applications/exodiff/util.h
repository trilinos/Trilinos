/*
 * Copyright(C) 1999-2025 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once
#include "vector_data.h"
#include <fmt/color.h>
#include <sstream>
#include <string>

using NameList = std::vector<std::string>;

int    name_length();
char **get_name_array(size_t size, size_t length);
void   free_name_array(char **names, size_t size);

[[noreturn]] void Error(const std::string &x);
[[noreturn]] void Error(std::ostringstream &buf);
void              Warning(const std::string &x);
void              ERR_OUT(std::ostringstream &buf);
void              WARN_OUT(std::ostringstream &buf);
void DIFF_OUT(std::ostringstream &buf, fmt::detail::color_type color = fmt::color::red);
void DIFF_OUT(const std::string &buf, fmt::detail::color_type color = fmt::color::red);
