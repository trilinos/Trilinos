#include <iostream>
#include <string>
#include <map>
#include "NormalBlockRTC.hh"
#include "BlockRTC.hh"
#include "commonRTC.hh"

using namespace std;
using namespace PG_RuntimeCompiler;

/*****************************************************************************/
NormalBlock::NormalBlock(map<string, Variable*> vars, Tokenizer& lines, 
                         string& errs) : Block(vars) 
/*****************************************************************************/
{
  createSubStatements(lines, errs);
}

/*****************************************************************************/
Value* NormalBlock::execute() 
/*****************************************************************************/
{
  list<Executable*>::iterator itr = _statements.begin();
  
  //execute all substatements
  while(itr != _statements.end()) {
    (*itr)->execute();
    ++itr;
  }
  return NULL;
}  
