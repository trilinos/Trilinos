#include "BlockRTC.hh"
#include "LineRTC.hh"
#include "commonRTC.hh"
#include "VariableRTC.hh"
#include "ConditionalBlockRTC.hh"

#include <string>
#include <map>

using namespace std;
using namespace PG_RuntimeCompiler;

/*****************************************************************************/
ConditionalBlock::
ConditionalBlock(map<string, Variable*> vars, Tokenizer& lines, string& errs) 
  : Block(vars) 
/*****************************************************************************/
{
  if (errs != "") return;

  _executed  = false;
  _condition = NULL;
  
  lines.nextToken(); //move past "if"

  //create the conditional statement
  _condition = new Line(lines, this, errs, true);
  if (errs != "") return;

  //create this block's substatements
  createSubStatements(lines, errs);
}

/*****************************************************************************/
ConditionalBlock::~ConditionalBlock() 
/*****************************************************************************/
{
  if (_condition != NULL)
    delete _condition;
}

/*****************************************************************************/
Value* ConditionalBlock::execute() 
/*****************************************************************************/
{
  //evaluate condition
  if (_condition->execute()->getValue()) {
    list<Executable*>::iterator itr = _statements.begin();
      
    //execute every substatement in the block
    while(itr != _statements.end()) {
      (*itr)->execute();
      ++itr;
    }
    _executed = true;
  }

  return NULL;
}
