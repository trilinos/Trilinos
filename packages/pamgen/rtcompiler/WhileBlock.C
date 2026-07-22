// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "RTC_WhileBlockRTC.hh"
#include "RTC_ForBlockRTC.hh"
#include "RTC_LineRTC.hh"
#include "RTC_BlockRTC.hh"
#include "RTC_IfElseifElseBlockRTC.hh"
#include "RTC_commonRTC.hh"

#include <string>
#include <map>

using namespace std;
using namespace PG_RuntimeCompiler;

/*****************************************************************************/
WhileBlock::WhileBlock(
    map<string, Variable*> vars,
    Tokenizer& lines,
    string& errs
    ):
  Block(vars),
  _condition(NULL)
{
  if (errs != "") return;

  lines.nextToken(); //move past "while"

  //create the condition statement
  _condition = new Line(lines, this, errs, true);
  if (errs != "") return;

  createSubStatements(lines, errs);
}

/*****************************************************************************/
WhileBlock::~WhileBlock()
{
  if (_condition != NULL)
    delete _condition;
}

/*****************************************************************************/
Value* WhileBlock::execute()
{
  while (_condition->execute()->getValue()) {
    list<Executable*>::iterator itr = _statements.begin();

    while(itr != _statements.end()) {
      (*itr)->execute();
      ++itr;
    }
  }
  return NULL;
}
