// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "RTC_ForBlockRTC.hh"
#include "RTC_WhileBlockRTC.hh"
#include "RTC_LineRTC.hh"
#include "RTC_BlockRTC.hh"
#include "RTC_IfElseifElseBlockRTC.hh"
#include "RTC_commonRTC.hh"

#include <string>
#include <map>
#include <stack>

using namespace std;
using namespace PG_RuntimeCompiler;

/*****************************************************************************/
ForBlock::ForBlock(
    map<string, Variable*> vars, Tokenizer& lines,
    string& errs
    ):
  Block(vars),
  _init(NULL),
  _condition(NULL),
  _postloop(NULL)
{
  if (errs != "") return;

  lines.nextToken(); //move past "for" token

  _init      = new Line(lines, this, errs, false);
  if (errs != "") return;

  _condition = new Line(lines, this, errs, true);
  if (errs != "") return;

  _postloop  = new Line(lines, this, errs, false);
  if (errs != "") return;

  createSubStatements(lines, errs);
}

/*****************************************************************************/
ForBlock::~ForBlock()
{
  if (_init != NULL)
    delete _init;
  if (_condition != NULL)
    delete _condition;
  if (_postloop != NULL)
    delete _postloop;
}

/*****************************************************************************/
Value* ForBlock::execute()
{
  //execute the initialization statement (ex: int i = 0)
  _init->execute();

  //check the conditional statement (ex: i < 10)
  while (_condition->execute()->getValue()) {
    list<Executable*>::iterator itr = _statements.begin();

    //execute all sub statements in the block
    while(itr != _statements.end()) {
      (*itr)->execute();
      ++itr;
    }

    //execute the postloop statement (ex: ++i)
    _postloop->execute();
  }
  return NULL;
}
