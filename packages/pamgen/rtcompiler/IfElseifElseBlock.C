// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "RTC_IfElseifElseBlockRTC.hh"
#include "RTC_BlockRTC.hh"
#include "RTC_ConditionalBlockRTC.hh"
#include "RTC_NormalBlockRTC.hh"
#include "RTC_commonRTC.hh"

#include <string>
#include <list>
#include <map>

using namespace std;
using namespace PG_RuntimeCompiler;

/*****************************************************************************/
IfElseifElseBlock::
IfElseifElseBlock(
    map<string, Variable*> vars,
    Tokenizer& lines, string& errs
    ):
  Block(vars),
  _if(NULL),
  _else(NULL)
{
  if (errs != "") return;

  //create the conditionalBlock for the if part
  _if = new ConditionalBlock(_vars, lines, errs);
  if (errs != "" || lines.eof()) return;

  //create the conditional blocks for the else-if parts
  while (lines.token() == Token(BLOCKOPENER, "else if", 0)) {
    _elseifs.push_back(new ConditionalBlock(_vars, lines, errs));
    if (errs != "" || lines.eof()) return;
  }

  //create the block for the else part
  if (lines.token() == Token(BLOCKOPENER, "else", 0)) {
    lines.nextLine(); //skip else statement and treat it as a normal block
    _else = new NormalBlock(_vars, lines, errs);
  }
}

/*****************************************************************************/
IfElseifElseBlock::~IfElseifElseBlock()
{
  delete _if;
  list<ConditionalBlock*>::iterator itr = _elseifs.begin();

  while (itr != _elseifs.end()) {
    delete (*itr);
    ++itr;
  }

  if (_else != NULL)
    delete _else;
}

/*****************************************************************************/
Value* IfElseifElseBlock::execute()
{
  _if->execute();

  if (_if->wasExecuted()) {
    _if->reset();
    return NULL;
  }

  list<ConditionalBlock*>::iterator itr = _elseifs.begin();

  while(itr != _elseifs.end()) {
    (*itr)->execute();

    if ( (*itr)->wasExecuted()) {
      (*itr)->reset();
      return NULL;
    }

    ++itr;
  }

  if (_else != NULL)
    _else->execute();

  return NULL;
}

/*****************************************************************************/
ostream& IfElseifElseBlock::operator<<(ostream& os) const
{
  os << *_if;

  for (list<ConditionalBlock*>::const_iterator itr = _elseifs.begin();
      itr != _elseifs.end(); itr++) {
    os << *(*itr);
  }

  if (_else != NULL) {
    os << *_else;
  }
  return os;
}
