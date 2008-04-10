#include <string>
#include <list>
#include <map>
#include <iostream>
#include <math.h>
#include "BlockRTC.hh"
#include "VariableRTC.hh"
#include "LineRTC.hh"
#include "WhileBlockRTC.hh"
#include "ForBlockRTC.hh"
#include "IfElseifElseBlockRTC.hh"
#include "TokenizerRTC.hh"

using namespace std;
using namespace PG_RuntimeCompiler;

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
  for ( ; myVars != _varsIOwn.end(); ++myVars)
    delete (*myVars);
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
