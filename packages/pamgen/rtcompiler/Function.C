#include <string>
#include <list>
#include <stack>
#include <math.h>
#include "FunctionRTC.hh"
#include "BlockRTC.hh"
#include "NormalBlockRTC.hh"
#include "VariableRTC.hh"
#include "commonRTC.hh"
#include "ArrayVarRTC.hh"
#include "ScalarVarRTC.hh"
#include "TokenizerRTC.hh"

using namespace std;
using namespace PG_RuntimeCompiler;

/*****************************************************************************/
Function::Function(const unsigned varCount,const string & name):
/*****************************************************************************/
_name(name),
_errors(""),
_mainBlock(0),
_vars(),
_firstRun(true),
_compiled(false),
_clean(false)
{
  _vars.reserve(varCount);
}


/*****************************************************************************/
bool Function::addVar(const string& type, const string& name) 
/*****************************************************************************/
{
  //check to make sure the name is valid
  if (name == "") {
    _errors += "Illegal name given for argument\n";
    for (unsigned int i = 0; i < name.length(); ++i) {
      if (!isValidVariableChar(name[i])) {
	_errors += "Illegal name given for argument. Make sure you did not include whitespace in the strings passed to addVar.\n";
	return false;
      }
    }
  }


  //create a new variable
  Variable* ptr;
  if (type == "int") 
    ptr = new ScalarVar<int>(name, IntT);
  else if (type == "char")
    ptr = new ScalarVar<char>(name, CharT);
  else if (type == "float")
    ptr = new ScalarVar<float>(name, FloatT);
  else if (type == "double")
    ptr = new ScalarVar<double>(name, DoubleT);
  else if (type == "int[]")
    ptr = new ArrayVar<int>(name, IntT);
  else if (type == "char[]")
    ptr = new ArrayVar<char>(name, CharT);
  else if (type == "float[]")
    ptr = new ArrayVar<float>(name, FloatT);
  else if (type == "double[]")
    ptr = new ArrayVar<double>(name, DoubleT);
  else {
    _errors += "Illegal type provided for argument. Make sure you did not have whitespace in the strings passed to addVar.\n";
    return false;
  }

  //add the variable to the array of argument variables
  ptr->init();
  _vars.push_back(ptr);

  return true;
}

/*****************************************************************************/
bool Function::addBody(const string& body) 
/*****************************************************************************/
{
  Tokenizer tokens(body, _errors);
  
  if (_errors == "") {
    map<string, Variable*> temp; //a map of variables that are in scope
    
    //add arguments to map of available variables
    for (unsigned int i = 0; i < _vars.size(); ++i) 
      temp[_vars[i]->getName()] = _vars[i];
    
    //Note: additional errors can occur here
    _mainBlock = new NormalBlock(temp, tokens, _errors);
  }
  
  _compiled = (_errors == "");

  return _compiled;
}

/*****************************************************************************/
bool Function::varValueFill(unsigned int index, double value) 
/*****************************************************************************/
{
  if (index >= _vars.size()) {
    _errors += "Index is too large.\n";
    return false;
  }
  if (_vars[index]->getObjectType() != ScalarVarOT) {
    _errors += "Must use the varValueFill on scalar variables.\n";
    return false;
  }

  _vars[index]->setValue(value);
  _vars[index]->init();
  return true;
}

/*****************************************************************************/
bool Function::execute() 
/*****************************************************************************/
{
  if (!_compiled) {
    _errors += "Cannot run. Program has not compiled successfully.\n";
    return false;
  }

  if (_firstRun) {
    for (unsigned int i = 0; i < _vars.size(); ++i) {
      if (!_vars[i]->isSet()) {
	_errors += "Not all arguments were set.\n";
	return false;
      }
    }
    _firstRun = 0;
  }

  _mainBlock->execute();
  
  return true;
}

/*****************************************************************************/
double Function::getValueOfVar(const std::string& var) 
/*****************************************************************************/
{
  if (!_clean)
    return (_mainBlock->getVar(var))->getValue();
  else
    return 0;
}

/*****************************************************************************/
bool Function::cleanup() 
/*****************************************************************************/
{
  if (_clean) return true;

  if (_mainBlock != 0)
  {
    delete _mainBlock;
    _mainBlock=0;
  }

  for (unsigned int i = 0; i < _vars.size(); ++i) {
    delete _vars[i];
  }

  _vars.clear();

  _clean = true;

  return true;
}

/*****************************************************************************/
Function::~Function() 
/*****************************************************************************/
{
  cleanup();
}

void Function::checkType(unsigned int index, int size, int* addr, string& errs)
{
  CHECKARGERR((_vars[index]->getType() != IntT  || 
               (_vars[index]->getObjectType() == ScalarVarOT && size != 0) || 
               (_vars[index]->getObjectType() == ArrayVarOT  && size == 0)),
              type_err(index))
}
void Function::checkType(unsigned int index, int size, float* addr, string& errs) 
{
  CHECKARGERR((_vars[index]->getType() != FloatT  || 
               (_vars[index]->getObjectType() == ScalarVarOT && size != 0) || 
               (_vars[index]->getObjectType() == ArrayVarOT  && size == 0)),
              type_err(index))
}
void Function::checkType(unsigned int index, int size, double* addr, string& errs) 
{
  CHECKARGERR((_vars[index]->getType() != DoubleT  || 
               (_vars[index]->getObjectType() == ScalarVarOT && size != 0) || 
               (_vars[index]->getObjectType() == ArrayVarOT  && size == 0)),
              type_err(index))
}
void Function::checkType(unsigned int index, int size, char* addr,string& errs)
{
  CHECKARGERR((_vars[index]->getType() != CharT  || 
               (_vars[index]->getObjectType() == ScalarVarOT && size != 0) || 
               (_vars[index]->getObjectType() == ArrayVarOT  && size == 0)),
              type_err(index))
}
