#include "WhileBlockRTC.hh"
#include "ForBlockRTC.hh"
#include "LineRTC.hh"
#include "BlockRTC.hh"
#include "IfElseifElseBlockRTC.hh"
#include "commonRTC.hh"

#include <string>
#include <map>

using namespace std;
using namespace PG_RuntimeCompiler;

/*****************************************************************************/
WhileBlock::WhileBlock(map<string, Variable*> vars, Tokenizer& lines, 
                       string& errs) : Block(vars) 
/*****************************************************************************/
{
  if (errs != "") return;

  _condition = NULL;

  lines.nextToken(); //move past "while"

  //create the condition statement
  _condition = new Line(lines, this, errs, true);
  if (errs != "") return;

  createSubStatements(lines, errs);
}

/*****************************************************************************/
WhileBlock::~WhileBlock() 
/*****************************************************************************/
{
  if (_condition != NULL)
    delete _condition;
}

/*****************************************************************************/
Value* WhileBlock::execute() 
/*****************************************************************************/
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
