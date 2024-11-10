// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _TOKENIZER_H
#define _TOKENIZER_H

#include "RTC_commonRTC.hh"

#include <string>
#include <vector>
#include <stack>
#include <cassert>

namespace PG_RuntimeCompiler {

/**
 *
 */
class Tokenizer {
 public:
  Tokenizer(const std::string& body, std::string& errors);

  Tokenizer(std::vector<Token> tokens);

  static void test();

  void print();

  void printCurrLine();

  int lineNum();

  void nextLine();

  void nextToken();

  void previousToken();

  Token token() {assert(!eol()); return *_currToken;}

  bool eof() { return _currLine == _tokenizedLines.end();}

  bool eol() { return _currToken == _currLine->end();}

  bool isArg() { return _tokenizedLines.size() == 1 && _currLine->size() == 1;}

private:
  Tokenizer(const Tokenizer&) {} //no copying allowed

  bool checkStack(const std::string& op);

  bool check(const std::string& testName, const std::string& errs,
             unsigned int expectedNumLines,
             unsigned int* expectedNumTokens,
             Token** expectedTokens);

  std::vector<std::vector<Token> >           _tokenizedLines;
  std::vector<std::string>                   _stack;
  std::vector<std::vector<Token> >::iterator _currLine;
  std::vector<Token>::iterator               _currToken;
};

}
#endif
