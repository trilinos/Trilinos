// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "RTC_BlockRTC.hh"
#include "RTC_VariableRTC.hh"
#include "RTC_LineRTC.hh"
#include "RTC_WhileBlockRTC.hh"
#include "RTC_ForBlockRTC.hh"
#include "RTC_IfElseifElseBlockRTC.hh"
#include "RTC_TokenizerRTC.hh"

#include <string>
#include <list>
#include <map>
#include <iostream>

using namespace std;
using namespace PG_RuntimeCompiler;

int Block::indent = 0;

/*****************************************************************************/
Block::Block(const map<string, Variable*>& vars)
  : _vars(vars)
/*****************************************************************************/
{

}

/*****************************************************************************/
Block::~Block()
/*****************************************************************************/
{
  list<Executable*>::iterator itr = _statements.begin();

  //delete all sub-blocks/sub-lines
  while(itr != _statements.end()) {
    delete (*itr);
    ++itr;
  }

  //delele all owned variables
  list<Variable*>::iterator myVars = _varsIOwn.begin();
  for ( ; myVars != _varsIOwn.end(); ++myVars) {
    delete (*myVars);
  }
}

/*****************************************************************************/
void Block::addStatement(Executable* statement)
/*****************************************************************************/
{
  _statements.push_back(statement);
}

/*****************************************************************************/
void Block::addVariable(Variable* var)
/*****************************************************************************/
{
  _vars[var->getName()] = var;
  _varsIOwn.push_back(var);
}

/*****************************************************************************/
Variable* Block::getVar(const string& name)
/*****************************************************************************/
{
  map<string, Variable*>::iterator itr = _vars.find(name);
  if (itr != _vars.end())
    return itr->second;
  else
    return NULL;
}

/*****************************************************************************/
void Block::createSubStatements(Tokenizer& lines, string& errs)
/*****************************************************************************/
{
  if (errs != "") return;

  while (!lines.eof() && lines.token().type() != CLOSEBRACE) {
    if (lines.token().type() == BLOCKOPENER) {
      if (lines.token().value() == "if")
	addStatement(new IfElseifElseBlock(_vars, lines, errs));
      else if (lines.token().value() == "while")
	addStatement(new WhileBlock(_vars, lines, errs));
      else if (lines.token().value() == "for")
	addStatement(new ForBlock(_vars, lines, errs));
      else
	errs += "Unexpected block opener: " + lines.token().value() +
	  " at line: " + intToString(lines.lineNum());
    }
    //otherwise we are dealing with an ordinary line
    else
      addStatement(new Line(lines, this, errs, false));

    if (errs != "") return;
  }

  //move past the closing brace
  if (!lines.eof())
    lines.nextLine();
}

/*****************************************************************************/
ostream& Block::operator<<(ostream& os) const
/*****************************************************************************/
{
  string indent_str(indent, ' ');
  os << "PRINTING BLOCK" << endl;
  os << indent_str << "My variables: ";
  for (list<Variable*>::const_iterator itr = _varsIOwn.begin();
       itr != _varsIOwn.end(); itr++) {
    os << *(*itr) << ", ";
  }
  os << endl;

  os << indent_str << "Block begin:" << endl;
  for (list<Executable*>::const_iterator itr = _statements.begin();
       itr != _statements.end(); itr++) {
    indent += 2;
    string new_indent(indent, ' ');
    os << new_indent << *(*itr) << endl;
    indent -= 2;
  }
  return os;
}
