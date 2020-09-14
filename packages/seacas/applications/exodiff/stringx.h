// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef STRINGX_H
#define STRINGX_H

#include "util.h" // for free_name_array, etc

//! Compare a string against another "master" string, where the string, str,
//! can be abbreiviated to as little as min_length characters.  Returns true
//! only if str has at least min_length characters and those that it does
//! have match the master string exactly.
bool abbreviation(const std::string &s, const std::string &master, unsigned min_length);

//! Compares two string ignoring letter case.  Returns true if they are equal.
bool no_case_equals(const std::string &s1, const std::string &s2);

//! Removes whitespace from the end of the string.  Returns the string given
//! as the argument.
std::string &chop_whitespace(std::string &s);

//! Separates the next token from the given string.  The next token is
//! returned and the given string has the next token removed (so it is
//! modified in place).
std::string extract_token(std::string &s, const char *delimiters = " \t\n\r");

//! Counts how many tokens are contained in the given string.
int count_tokens(const std::string &s, const char *delimiters = " \t\n\r");

//! Runs each string in the vector and returns the maximum size.
int max_string_length(const std::vector<std::string> &names);

//! Replaces each character of the string with its lower case equivalent.
void to_lower(std::string &s);

//! Returns the first non-white space character of the string.  If the string
//! is empty or all white space, a nullptr char is returned.
char first_character(const std::string &s);

//! Searches the list of strings for a particular string value.  Letter case
//! will be ignored if the last argument is true.  Returns the index of the
//! string in the vector if found, otherwise returns -1.
int find_string(const std::vector<std::string> &lst, const std::string &s, bool nocase);

#endif
