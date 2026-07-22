// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "RTC_BlockRTC.hh"
#include "RTC_ObjectRTC.hh"
#include "RTC_ScalarNumberRTC.hh"
#include "RTC_ArrayNumberRTC.hh"
#include "RTC_commonRTC.hh"
#include "RTC_ArrayIndexRTC.hh"
#include "RTC_OperatorRTC.hh"
#include "RTC_RegistrarRTC.hh"
#include "RTC_TokenizerRTC.hh"
#include "RTC_VariableRTC.hh"
#include "RTC_ScalarVarRTC.hh"
#include "RTC_ArrayVarRTC.hh"
#include "RTC_LineRTC.hh"
#include "RTC_NormalBlockRTC.hh"

#include <string>
#include <list>
#include <vector>
#include <stack>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace PG_RuntimeCompiler;

/*****************************************************************************/
Line::Line(Tokenizer& tokens, Block* parent, string& errs, bool needReturn)
/*****************************************************************************/
{
  _parent        = parent;
  _needReturnVal = needReturn;
  _tempSize      = 0;
  _curr          = tokens.lineNum();
  _tempHolder    = NULL;

  infToPost(tokens, errs);

  tokens.nextLine();
}

/*****************************************************************************/
Line::~Line()
/*****************************************************************************/
{
  set<Object*>::iterator itr = _objsToDelete.begin();
  if (_objsToDelete.size() != 0) {
    while(itr != _objsToDelete.end()) {
        delete (*itr);
        ++itr;
    }
  }
  if (_tempHolder != NULL)
    delete[] _tempHolder;
}

/*****************************************************************************/
void Line::stackOp(stack<Operator*>& ops, Operator* op)
/*****************************************************************************/
{
  //open parens always go on the stack
  if (op->getName() == "(") {
    ops.push(op);
  }
  //if we see a close paren, pop every item on the stack until we see an open
  //paren. Popped items go to the end of the post fix list
  else if (op->getName() == ")") {
    while (true) {
      if (ops.empty()) {
	break;
      }
      else if (ops.top()->getName() == "(") {
	ops.pop();
	break;
      }
      else {
	addNewObject(ops.top());
	ops.pop();
      }
    }
  }
  //if the stack is empty, operators will always be pushed onto stack
  else if (ops.empty()) {
    ops.push(op);
  }
  //if we are dealing with a non-paren operator and the stack is not empty,
  //we must pop from the stack based on precedence and associativity
  else {
    while (!ops.empty() && ops.top()->getPrecedence() >= op->getPrecedence()) {
      //note: unary operators are right associative
      if (!ops.top()->isLeftAssociative() && !op->isLeftAssociative() &&
          ops.top()->getPrecedence() == op->getPrecedence())
        break;

      addNewObject(ops.top());
      ops.pop();
    }

    ops.push(op);
  }
}

/*****************************************************************************/
void Line::process_function(Tokenizer& tokens, string& errs)
/*****************************************************************************/
{
  //create a function call object
  string temp = tokens.token().value();
  FunctionCall* funcCall = Registrar::generateCall(temp);
  CHECKERR((funcCall == NULL), func_err(temp))

  tokens.nextToken(); //move past funcname
  tokens.nextToken(); //move past open paren

  //put tokens into the argument lists.. We loop until we have seen the paren
  //that terminates this function call or until we have run out of tokens on
  //this line
  list< vector<Token> > args;
  vector<Token> currArg;
  int depth = 0;
  while(!((tokens.token().value() == ")" && depth == 0) || tokens.eol())) {

    //if we see a comma at paren depth zero, we have just reached the end of an
    //argument
    if (tokens.token().type() == COMMA && depth == 0) {
      assert(!currArg.empty());
      args.push_back(currArg);
      currArg.clear();
    }
    else {
      currArg.push_back(tokens.token());

      if (tokens.token() == Token(OPERATOR, "(", 0))
        ++depth;
      if (tokens.token() == Token(OPERATOR, ")", 0))
        --depth;
    }
    tokens.nextToken();
  }
  if (!currArg.empty())
    args.push_back(currArg);

  CHECKERR(
      tokens.eol() || tokens.token().value() != ")",
      arg_err(temp)
      )

  if (funcCall->hasVariableArgs()) {
    CHECKERR (
        args.size() < funcCall->getNumArgs(),
        arg_err(temp)
        )
  } else {
    CHECKERR(
        args.size() != funcCall->getNumArgs(),
        arg_err(temp)
        )
  }

  //Construct a Line for each argument
  list< vector<Token> >::iterator arg_itr = args.begin();
  for ( ; arg_itr != args.end(); ++arg_itr) {
    CHECKERR (((*arg_itr).size() == 0), arg_err(temp))
    Tokenizer tempToken(*arg_itr);
    Line* arg = new Line(tempToken, _parent, errs, true);
    funcCall->fillArg(arg);
  }

  addNewObject(funcCall);
}

/*****************************************************************************/
void Line::add_newvar(string& type, string& newVar, Line* sizePtr,bool isArray)
/*****************************************************************************/
{
  Object* ptr;
  if (type == "int") {
    if (isArray)
      ptr = new ArrayVar<int>(newVar, sizePtr);
    else
      ptr = new ScalarVar<int>(newVar);
  }
  else if (type == "long") {
    if (isArray)
      ptr = new ArrayVar<long>(newVar, sizePtr);
    else
      ptr = new ScalarVar<long>(newVar);
  }
  else if (type == "float") {
    if (isArray)
      ptr = new ArrayVar<float>(newVar, sizePtr);
    else
      ptr = new ScalarVar<float>(newVar);
  }
  else if (type == "double") {
    if (isArray)
      ptr = new ArrayVar<double>(newVar, sizePtr);
    else
      ptr = new ScalarVar<double>(newVar);
  }
  else  {
    if (isArray)
      ptr = new ArrayVar<char>(newVar, sizePtr);
    else
      ptr = new ScalarVar<char>(newVar);
  }

  addNewObject(ptr);
  _parent->addVariable( (Variable*) ptr);

  if (isArray)
    addNewObject(Operator::getOp("init"));
}

/*****************************************************************************/
void Line::process_newvar(Tokenizer& tokens, string& errs)
/*****************************************************************************/
{
  string type = tokens.token().value();
  tokens.nextToken();

  CHECKERR(tokens.eol() || tokens.token().type() != VAR, var_err())

  string name = tokens.token().value();
  tokens.nextToken();

  bool isArray = false;
  Line* sizePtr = NULL;

  //if the following token is an open index, we know that our new variable
  //is an array
  if (!tokens.eol() && tokens.token().type() == OPENINDEX) {
    vector<Token> size_expr;
    tokens.nextToken(); //move past openindex

    //get all the tokens that are part of the array's size expression
    while (!tokens.eol() && tokens.token().type() != CLOSEINDEX) {
      size_expr.push_back(tokens.token());
      tokens.nextToken();
    }
    CHECKERR ((size_expr.size() == 0), ara_err(name))
    CHECKERR ((tokens.token().type() != CLOSEINDEX), ara_err(name))

    isArray = true;
    Tokenizer tempToken(size_expr);
    sizePtr = new Line(tempToken, _parent, errs, true);
  }
  else
    tokens.previousToken();

  if (_parent->getVar(name) == NULL)
    add_newvar(type, name, sizePtr, isArray);
  else
    { CHECKERR(true, dec_err(name)) }
}

/*****************************************************************************/
void Line::process_existing_var(Tokenizer& tokens, string& errs)
/*****************************************************************************/
{
  string temp = tokens.token().value();
  Variable* v = _parent->getVar(temp);
  CHECKERR ((v == NULL), und_err(temp))

  //Note: we must allow for arrays to be passed to RTBoundFuncs without
  //having to use braces [].
  if (tokens.isArg()) {
    addNewObject(v);
    return;
  }

  //When we see an array variable, it must be followed by an index
  if (v->getObjectType() == ArrayVarOT) {
    tokens.nextToken();
    CHECKERR((tokens.eol() || tokens.token().type()!=OPENINDEX),ara_err(temp))
    tokens.nextToken(); //move past OPENINDEX
    vector<Token> index_list;

    //get all the tokens that are part of the array's index expression
    while (!tokens.eol() && tokens.token().type() != CLOSEINDEX) {
      index_list.push_back(tokens.token());
      tokens.nextToken();
    }
    CHECKERR ((index_list.size() == 0), ara_err(temp))
    CHECKERR ((tokens.eol()||tokens.token().type()!=CLOSEINDEX), ara_err(temp))

    Tokenizer tempToken(index_list);
    Line* indexPtr = new Line(tempToken, _parent, errs, true);
    ArrayIndex* ai = new ArrayIndex(v, indexPtr);
    addNewObject(ai);
  }
  else {
    addNewObject(v);
  }
}

/*****************************************************************************/
void Line::process_number(Tokenizer& tokens)
/*****************************************************************************/
{
  Object* ptr;
  string temp = tokens.token().value();
  if (isInt(temp)) {
    long value = stol(temp);
    ptr = new ScalarNumber<long>(value);
  }
  else if (isDouble(temp)) {
    double value = stod(temp);
    ptr = new ScalarNumber<double>(value);
  }
  else if (isString(temp)) {
    ptr = new ArrayNumber<char>(temp.c_str(), temp.size());
  }
  else {
    assert(isChar(temp));
    char value = (char) temp[1];
    ptr = new ScalarNumber<char>(value);
  }
  addNewObject(ptr);
  _objsToDelete.insert(ptr);
}

/*****************************************************************************/
void Line::process_operator(Tokenizer& tokens, stack<Operator*>& ops)
/*****************************************************************************/
{
  Operator* op = Operator::getOp(tokens.token().value());

  stackOp(ops, op);

  if (op->getName() != ")" && op->getName() != "(")
    ++_tempSize;
}

/*****************************************************************************/
void Line::infToPost(Tokenizer& tokens, string& errs)
/*****************************************************************************/
{
  stack<Operator*> ops;

  //loop through all tokens and handle them
  for ( ; !tokens.eol(); tokens.nextToken()) {
    if (tokens.token().type() == FUNC) {
      process_function(tokens, errs);
      ++_tempSize; //for the return value
    }
    else if (tokens.token().type() == DECL) {
      process_newvar(tokens, errs);
    }
    else if (tokens.token().type() == CONSTANT) {
      process_number(tokens);
    }
    else if (tokens.token().type() == VAR) {
      process_existing_var(tokens, errs);
    }
    else if (tokens.token().type() == OPERATOR) {
      process_operator(tokens, ops);
    }
    else if (tokens.token().type() == SEMICOLON ||
	     tokens.token().type() == OPENBRACE) {
      tokens.nextToken();
      assert(tokens.eol());
      break;
    }
    else {
      CHECKERR(true, syntax_err(tokens.token().value()))
    }
    if (errs != "") return;
  }

  //put remaining opps at end of postfixLine
  while (!ops.empty()) {
    addNewObject(ops.top());
    ops.pop();
  }

  compile(errs, tokens);
  performNumericOps();
}

/*****************************************************************************/
void Line::compile(string& errs, Tokenizer& tokens)
/*****************************************************************************/
{
  int assignments = 0;
  int numOps = 0;
  int numArgs = 0;
  list<Object*>::iterator itr = _postfixLine.begin();

  while (itr != _postfixLine.end()) {
    ObjectType currType = (*itr)->getObjectType();
    if (currType == OperatorOT) {
      Operator* op = (Operator*)(*itr);
      if (op->getName() != "init" && !op->isUnary())
	++numOps;
      if (op->getName() == "=")
	++assignments;
    }
    else {
      ++ numArgs;
      if (isVariable(currType) && itr != _postfixLine.begin()) {
	Variable* v = (Variable*)(*itr);
	CHECKERR((!v->isInit()),uninit_err(v->getName()))
      }
    }
    ++itr;
  }

  CHECKERR((numArgs != (numOps+1)), opp_err())
  CHECKERR((assignments > 1), ass_err())

  if (assignments == 1) {
    Object* ptr = *(_postfixLine.begin());
    CHECKERR ((!isAssignable(ptr->getObjectType())), nonv_err())

    if(ptr->getObjectType() == ScalarVarOT) {
      Variable* v = (Variable*)ptr;
      v->init();
    }
  }
}

/*****************************************************************************/
void Line::performNumericOps()
/*****************************************************************************/
{
  list<Object*>::iterator itr = _postfixLine.begin();
  list<Object*>::iterator arg1;
  list<Object*>::iterator arg2;
  list<Object*>::iterator op;

  while (itr != _postfixLine.end()) {
    if ( (*itr)->getObjectType() == OperatorOT &&
         ((Operator*)*itr)->getName() != "=" &&
	 ((Operator*)*itr)->getName() != "init") {
      op = itr;
      if (itr != _postfixLine.begin()) {
	arg2 = --itr;
	if (((Operator*)(*op))->isUnary()) {
	  if ( (*arg2)->getObjectType() == ScalarNumberOT ) {
	    //if we have a unary operator operating on a constant, we can
            //precompute the result

            --_tempSize; //we no longer need a temporary for this operation

            //create a new constant containing the result of the operation
	    ScalarNumber<double>* d = new ScalarNumber<double>(0);
	    ((Operator*)(*op))->performOp(((Value*)*arg2), *d);

            //put this new value in postfix
	    _postfixLine.insert(arg2, d);

            //put this new value in set of deletable items
            //(Lines delete their own constants)
            _objsToDelete.insert(d);

            //remove the old constant from set of items to delete
            assert(_objsToDelete.find(*arg2) != _objsToDelete.end());
            _objsToDelete.erase(*arg2);

            delete *arg2; //delete old constant

            //remove old constant and operator from postfix
	    _postfixLine.erase(arg2);
	    itr = _postfixLine.erase(op);

	    continue;
	  }
	  else
	    ++itr;
	}
	else if (itr != _postfixLine.begin()) {
	  arg1 = --itr;
	  if ( (*arg1)->getObjectType() == ScalarNumberOT &&
               (*arg2)->getObjectType() == ScalarNumberOT) {
	    //if we have a binary operator operating on two constants, we can
            //precompute the result

            --_tempSize; //we no longer need a temporary for this operation

            //create a new constant containing the result of the operation
	    ScalarNumber<double>* d = new ScalarNumber<double>(0);
	    ((Operator*)(*op))->performOp(((Value*)*arg1),((Value*)*arg2),*d);

            //put this new value in postfix
	    _postfixLine.insert(arg1, d);

            //put this new value in set of deletable items
            _objsToDelete.insert(d);

            //remove the old constants from set of items to delete
            assert(_objsToDelete.find(*arg1) != _objsToDelete.end());
            assert(_objsToDelete.find(*arg2) != _objsToDelete.end());
            _objsToDelete.erase(*arg1);
            _objsToDelete.erase(*arg2);

            //delete old constant objects
            delete *arg1;
            delete *arg2;

            //remove old constants and operator from postfix
	    _postfixLine.erase(arg1);
	    _postfixLine.erase(arg2);
	    itr = _postfixLine.erase(op);

	    continue;
	  }
	  else {
	    ++itr;
	    ++itr;
	  }
	}
	else
	  ++itr;
      }
    }
    else if ((*itr)->getObjectType() == FunctionOT) {
      if ( ((FunctionCall*)(*itr))->canGoEarly() ) {
        //if we have an optimizable function with constant arguments, we can
        //precompute the result

        --_tempSize; //we no longer need a temporary for this operation

        //create a new constant containing the result of the operation
        ScalarNumber<double>* d =
          new ScalarNumber<double>(((FunctionCall*)(*itr))->execute());

        //put this new value in postfix
        _postfixLine.insert(itr, d);

        //put this new value in set of deletable items
        _objsToDelete.insert(d);

        //remove the old constant from set of items to delete
        assert(_objsToDelete.find(*itr) != _objsToDelete.end());
        _objsToDelete.erase(*itr);

        //delete the function call
        delete *itr;

        //remove old constant and operator from postfix
        itr = _postfixLine.erase(itr);

        continue;
      }
    }
    ++itr;
  }

  if (_tempSize > 0)
    _tempHolder = new ScalarNumber<double>[_tempSize];
}

/*****************************************************************************/
Value* Line::execute()
/*****************************************************************************/
{
  stack<Value*> vals;
  list<Object*>::iterator itr = _postfixLine.begin ();
  int curr = 0;

  while(itr != _postfixLine.end()) {
    if ((*itr)->getObjectType() == OperatorOT) {
      Operator* op = ((Operator*)(*itr));
      if (op->getName() == "init") {
	((Variable*)vals.top())->evaluateSizeExpr();
      }
      else if (op->isUnary()) {
	op->performOp(vals.top(), _tempHolder[curr]);
	vals.pop();
	vals.push(&_tempHolder[curr++]);
      }
      else {
	Value* arg = vals.top();
	vals.pop();
	op->performOp(vals.top(), arg, _tempHolder[curr]);
	vals.pop();
	vals.push(&_tempHolder[curr++]);
      }
    }
    else if ((*itr)->getObjectType() == FunctionOT) {
      FunctionCall* fc = ((FunctionCall*)(*itr));
      _tempHolder[curr].setValue(fc->execute());
      vals.push(&_tempHolder[curr++]);
    }
    else {
      vals.push((Value*) (*itr));
    }
    ++itr;
  }

  if (!_needReturnVal)
    return NULL;
  else
    return vals.top();
}

/*****************************************************************************/
bool Line::can_go() const
/*****************************************************************************/
{
  //jgfouca: why can't multi-argument functions be pre-computed if all the
  //         arguments are constants?
  return (_postfixLine.size() == 1 &&
          _postfixLine.front()->getObjectType() == ScalarNumberOT);
}

/*****************************************************************************/
ostream& Line::operator<<(ostream& os) const
/*****************************************************************************/
{
  list<Object*>::const_iterator itr = _postfixLine.begin();

  for ( ; itr != _postfixLine.end(); ++itr) {
    os << *(*itr) << " ";
  }
  return os;
}

/*****************************************************************************/
bool Line::indvTest(const string& line, Block* parent, double expectedResult,
		    bool examineResult)
/*****************************************************************************/
{
  string errs = "";
  Tokenizer tokens(line, errs);
  if (errs != "") {
    cerr << "Line: " << line << " did not make it past tokenizing phase."
	 << " Tokenizer returned the error: " << errs << endl;
    return false;
  }

  Line lineObj(tokens, parent, errs, examineResult);

  if (errs != "") {
    cerr << "Line: " << line << " had errors: " << errs << endl;
    cerr << lineObj;
    return false;
  }

  if (examineResult) {
    double result = lineObj.execute()->getValue();

    if (expectedResult != result) {
      cerr << "Line: " << line << " was expected to generate result: "
	   << expectedResult << " instead, it generated: " << result << endl;
      cerr << lineObj;
      return false;
    }
  }
  else
    lineObj.execute();

  cout << "Line test passed: " << line << endl;

  return true;
}

/*****************************************************************************/
void Line::test()
/*****************************************************************************/
{
  map<string, Variable*> dummy;
  string errs = "";
  Tokenizer dummyTokens("int i = 0;", errs);

  NormalBlock emptyBlock(dummy, dummyTokens, errs);

  if (!indvTest("int temp = 5;", &emptyBlock, 5.0, true)) return;

  if (!indvTest("double array[5];", &emptyBlock, 0, false)) return;

  if (!indvTest("array[4-1-2] = 4.34;", &emptyBlock, 4.34, true)) return;

  if (!indvTest("double res = array[1] * sin(3);", &emptyBlock, (4.34*sin(3.0)), true)) return;

  if (!indvTest("(temp == 0)", &emptyBlock, 0.0, true)) return;

  if (!indvTest("temp = -7;", &emptyBlock, -7.0, true)) return;

  if (!indvTest("temp = -6*-7;", &emptyBlock, 42, true)) return;

  if (!indvTest("temp = --7;", &emptyBlock, 7.0, true)) return;

  if (!indvTest("!0 && !0;", &emptyBlock, 1.0, true)) return;

  if (!indvTest("arrayFill(array, -3);", &emptyBlock, -3.0, true)) return;

  if (!indvTest("array[0];", &emptyBlock, -3.0, true)) return;
  if (!indvTest("array[1];", &emptyBlock, -3.0, true)) return;
  if (!indvTest("array[2];", &emptyBlock, -3.0, true)) return;
  if (!indvTest("array[3];", &emptyBlock, -3.0, true)) return;
  if (!indvTest("array[4];", &emptyBlock, -3.0, true)) return;

  //Test sci notation
  if (!indvTest("double d = 1e10;", &emptyBlock, 1e10, true)) return;
  if (!indvTest(" d = 1e+10;", &emptyBlock, 1e10, true)) return;
  if (!indvTest(" d = 1.15e-10;", &emptyBlock, 1.15e-10, true)) return;
  if (!indvTest(" d = -1E-10;", &emptyBlock, -1e-10, true)) return;

  //Test comments
  if (!indvTest("/*comment*/ double c =  /*comment*/ 23;",
                &emptyBlock, 23, true)) return;

  //Need to test a variable argument function
  if (!indvTest("printf(\"One:% Two:% Three:% \", 5-4, 2.0e0, 'c');",
                &emptyBlock, 22, true)) return;

  //welcoming ideas for more
}

/*****************************************************************************/
void Line::addNewObject(Object* newObj)
/*****************************************************************************/
{
  _postfixLine.push_back(newObj);
  if (!isVariable(newObj->getObjectType()) &&
      newObj->getObjectType() != OperatorOT) {
    _objsToDelete.insert(newObj);
  }
}
