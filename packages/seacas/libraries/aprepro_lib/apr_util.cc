// Copyright (c) 2014-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "aprepro.h"        // for symrec, Aprepro, etc
#include "aprepro_parser.h" // for Parser, Parser::token, etc
#include <cctype>           // for isalnum, isalpha, isupper, etc
#include <cerrno>           // for errno, EDOM, ERANGE
#include <cfenv>            // for fetestexcept, FE_DIVBYZERO, etc
#include <cmath>            // for math_errhandling, etc
#include <cstdio>           // for perror
#include <cstdlib>          // for mkstemp
#include <cstring>          // for strlen, strcpy, memcpy, etc
#include <iostream>         // for operator<<, cerr, ostream
#include <string>           // for allocator, operator+, etc
#include <sys/stat.h>       // for stat, S_ISDIR
#include <unistd.h>         // for close
#include <vector>           // for vector

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#include <windows.h>
#else
#include <unistd.h> // for close
#endif

#if !defined(S_ISDIR) && defined(_WIN32)
#define S_ISDIR(mode) (((mode)&S_IFMT) == S_IFDIR)
#endif

namespace {
  std::vector<char *> allocations;

  void new_string_int(const char *from, char **to)
  {
    int len = strlen(from);
    *to     = new char[len + 1];
    std::memcpy(*to, from, len + 1);
    allocations.push_back(*to);
  }

} // namespace

namespace SEAMS {
  extern Aprepro *aprepro;

  bool arg_check(SEAMS::symrec *symbol, bool is_null)
  {
    if (is_null) {
      std::string err_msg = "Incorrect argument count/type for function '" + symbol->name;
      err_msg += "'.\n                 ";
      err_msg += "The correct syntax is " + symbol->syntax;
      aprepro->error(err_msg);
      return false;
    }
    return true;
  }

  void set_type(const SEAMS::Aprepro &apr, SEAMS::symrec *var, int type)
  {
    if (var->name[0] == '_' || !apr.state_is_immutable()) {
      var->type = type;
    }
    else {
      if (type == Parser::token::VAR) {
        var->type = Parser::token::IMMVAR;
      }
      else if (type == Parser::token::SVAR) {
        var->type = Parser::token::IMMSVAR;
      }
      else {
        var->type = type;
      }
    }
  }

  void new_string(const std::string &from, char **to) { new_string_int(from.c_str(), to); }
  void new_string(const char *from, char **to) { new_string_int(from, to); }

  void concat_string(const char *from1, const char *from2, char **to)
  {
    int len = strlen(from1) + strlen(from2);
    *to     = new char[len + 1];
    std::strcpy(*to, from1);
    std::strcat(*to, from2);
    allocations.push_back(*to);
  }

  /* This function returns a pointer to a static character array.
   * If called multiple times, the name will change on each call,
   * so don't save the pointer; copy the data it points to if you
   * need to keep it for awhile.
   */
  char *get_temp_filename()
  {
    static char tmp_name[] = "./aprepro_temp_XXXXXX";
    int         fd;

    std::strcpy(tmp_name, "./aprepro_temp_XXXXXX");
#if defined(__CYGWIN__) && defined(__NO_CYGWIN_OPTION__)
    fd = mkstemps(tmp_name, 0);
    if (fd >= 0)
      close(fd);
#elif defined(_WIN32)
    std::strcpy(tmp_name, _mktemp(tmp_name));
#else
    fd = mkstemp(tmp_name);
    if (fd >= 0) {
      close(fd);
    }
#endif
    return tmp_name;
  }

  void yyerror(const SEAMS::Aprepro &apr, const std::string &s) { apr.error(s); }

  void immutable_modify(const SEAMS::Aprepro &apr, const SEAMS::symrec *var)
  {
    apr.error("(IMMUTABLE) Variable " + var->name + " is immutable and cannot be modified", true,
              false);
  }

  void undefined_error(const SEAMS::Aprepro &apr, const std::string &var)
  {
    if (!apr.inIfdefGetvar) {
      apr.warning("Undefined variable '" + var + "'");
    }
    else {
      apr.inIfdefGetvar = false;
    }
  }

  void redefined_warning(const SEAMS::Aprepro &apr, const SEAMS::symrec *var)
  {
    if (var->name[0] != '_' && apr.ap_options.warning_msg) {
      // See if internal or user-defined variable...
      std::string type;
      if (var->isInternal) {
        type = "Pre";
      }
      else {
        type = "User";
      }
      apr.warning(type + "-defined Variable '" + var->name + "' redefined");
    }
  }

  void warning(const SEAMS::Aprepro &apr, const std::string &s) { apr.warning(s); }

  void math_error(const SEAMS::Aprepro &apr, const char *function)
  {
    if (math_errhandling & MATH_ERRNO) {
      if (errno != 0) {
        yyerror(apr, function);
      }
      if (errno == EDOM) {
        perror("	DOMAIN error");
      }
      else if (errno == ERANGE) {
        perror("	RANGE error");
      }
      else if (errno != 0) {
        perror("	Unknown error");
      }
      errno = 0;
    }
    else if (math_errhandling & MATH_ERREXCEPT) {
      if (std::fetestexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO) != 0) {
        yyerror(apr, function);
      }

      if (std::fetestexcept(FE_INVALID) != 0) {
        std::cerr << "	DOMAIN error\n";
      }
      else if (std::fetestexcept(FE_OVERFLOW) != 0) {
        std::cerr << "  RANGE error -- overflow\n";
      }
      else if (std::fetestexcept(FE_DIVBYZERO) != 0) {
        std::cerr << "	RANGE error -- divide by zero\n";
      }
    }
  }

  void math_error(const char *function) { math_error(*SEAMS::aprepro, function); }

  /* Convert string to all lower-case and replace all spaces with '_' */
  void conv_string(char *string)
  {
    char *p = string;
    while (*p != '\0') {
      if (*p == ' ') {
        *p = '_';
      }
      else if (isupper(static_cast<int>(*p)) != 0) {
        *p = tolower(static_cast<int>(*p));
      }
      p++;
    }
  }

  void cleanup_memory()
  {
    for (auto &allocation : allocations) {
      delete[] allocation;
    }

    // Clear the vector to avoid stale pointers.
    allocations.clear();
  }

  bool is_directory(const std::string &filepath)
  {
    struct stat s
    {
    };
    int ok = stat(filepath.c_str(), &s);
    if (ok == 0) {
      return S_ISDIR(s.st_mode);
    }
    else {
      return false;
    }
  }

  bool check_valid_var(const char *var)
  {
    /* Check that 'var' meets the restriction for a variable name
     * L(:|L|D)*
     * D [0-9]
     * L [A-Za-z_]
     */

    int length = strlen(var);
    if (length == 0) {
      return false;
    }

    if ((isalpha(var[0]) == 0) && var[0] != '_') {
      return false;
    }

    for (int i = 1; i < length; i++) {
      char c = var[i];
      if ((isalnum(c) == 0) && c != ':' && c != '_') {
        return false;
      }
    }
    return true;
  }
} // namespace SEAMS
