// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _LINERTC_H
#define _LINERTC_H

#include "RTC_BlockRTC.hh"
#include "RTC_ObjectRTC.hh"
#include "RTC_commonRTC.hh"
#include "RTC_ExecutableRTC.hh"
#include "RTC_TokenizerRTC.hh"
#include "RTC_ScalarNumberRTC.hh"
#include "RTC_OperatorRTC.hh"

#include <string>
#include <list>
#include <stack>
#include <set>

namespace PG_RuntimeCompiler {

/**
 * A Line object represents a single line of code. This class is reposible
 * for parsing single lines of code and turning them into post-fix lists
 * of operands and operators. Lines of code can be executed.
 */

class Line: public Executable
{
 public:

  /**
   * Constructor -> The constructor will parse line and turn it into a postfix
   *                list.
   *
   * @param tokens     - A tokenized line
   * @param parent     - The block containing this line of code
   * @param errs       - A string containing errors
   * @param needReturn - Whether or not this line needs to return a value to
   *                     the entity that executed it
   */
  Line(Tokenizer& tokens, Block* parent, std::string& errs,
       bool needReturn);

  /**
   * Destructor -> The destructor will delete all objects in _postfixLine
   *               except for Variables. The Blocks are reponsible for deleting
   *               Variables. It will also delete _tempHolder.
   */
  ~Line();

  /**
   * stackop -> This method is used to sort operators into postfix.
   *
   * @param ops - A stack used to help produce postfix ordering
   * @param op  - The op we are adding to the stack
   */
  void stackOp(std::stack<Operator*>& ops, Operator* op);

  /**
   * infToPost -> This method is used by the constructor to convert a line
   *              into a postfix list
   *
   * @param tokens - A list of all the tokens in the list
   * @param err    - A string to hold errors
   */
  void infToPost(Tokenizer& tokens, std::string& errs);

  /**
   * compile -> This method does some additional checks to make sure that there
   *            are no errors in the Line. If any errors are found, they will
   *            be appended onto errs.
   *
   * @param err - A string to hold errors
   */
  void compile(std::string& errs, Tokenizer& tokens);

  /**
   * execute -> This method will execute this Line object.
   */
  Value* execute();

  /**
   * performNumericOps -> This method checks the operands at compile time. If
   *                      all the operands involved in an operation are
   *                      constants, it will perform this operation at compile
   *                      time.
   */
  void performNumericOps();

  /**
   * process_function -> Processes a function call that occurs in the user's
   *                     code. It parses the argument list and creates
   *                     executable argument statements for each argument.
   *
   * @param tokens - A list of tokens in this line
   * @param errs   - A string for holding errors
   */
  void process_function(Tokenizer& tokens, std::string& errs);

  /**
   * process_newvar -> Handles the creation of a new variable. Determines if
   *                   the variable is a scalar or an array.
   *
   * @param tokens - A list of tokens in this line
   * @param errs   - A string for holding errors
   */
  void process_newvar(Tokenizer& tokens, std::string& errs);

  /**
   * process_existing_var -> Handles the case where we come across a variable
   *                         that has been declared ealier. Array variables
   *                         must be indexed.
   *
   * @param tokens - A list of tokens in this line
   * @param itr    - An iterator that points to the token currently being
   *                 examined
   * @param errs   - A string for holding errors
   */
  void process_existing_var(Tokenizer& tokens, std::string& errs);

  /**
   * process_number -> Determines what type of number itr is pointing to and
   *                   creates the corresponding object.
   *
   * @param itr - An iterator that points to the token currently being
   *              examined
   */
  void process_number(Tokenizer& tokens);

  /**
   * process_operator -> Handles operators by delegating them to the stackOp
   *                     method.
   *
   * @param itr - An iterator that points to the token currently being
   *              examined
   * @param ops - A stack used to help produce postfix ordering of operations
   */
  void process_operator(Tokenizer& tokens, std::stack<Operator*>& ops);

  /**
   * add_newvar -> Adds a new variable into our system
   *
   * @param type    - The type of the new variable
   * @param newvar  - The name of the new variable
   * @param sizePtr - A pointer to a size expression (only applies to arrays)
   * @param isArray - Tells us if we are dealing with an array or not
   */
  void add_newvar(std::string& type, std::string& newvar, Line* sizePtr, bool isArray);

  /**
   * can_go -> This method is used by RTCBoundFunc objects to see check to see
   *           if the expressions (line objects) that make up their arguments
   *           were evaluated at compile time.
   */
  bool can_go() const;

  std::string func_err(const std::string& currVar) const {
    return "Function: " + currVar + " does not exist.\n";
  }

  std::string arg_err(const std::string& currVar) const {
    return "Incorrect number of args given to Function: "+currVar+"\n";
  }

  std::string var_err() const {
    return "Tried to declare a variable with no name. \n";
  }

  std::string ara_err(const std::string& currVar) const {
    return "Incorrect array syntax for array: " + currVar + '\n';
  }

  std::string dec_err(const std::string& newVar) const {
    return "Tried to declare a variable: " +newVar+ " that already exists.\n";
  }

  std::string und_err(const std::string& currVar) const {
    return "Use of undeclared variable: " + currVar + '\n';
  }

  std::string uninit_err(const std::string& var) const {
    return "Tried to use an uninitialized variable: " + var + '\n';
  }

  std::string opp_err() const {
    return "Number of operators does not match number of values. \n";
  }

  std::string ass_err() const {
    return "Only one assignment(=) per line of code is allowed. \n";
  }

  std::string nonv_err() const {
    return "Tried to assign to a non-assignable entity.\n";
  }

  std::string syntax_err(const std::string& value) const {
    return "Unexpected token " + value + '\n';
  }

  void static test(); //unit test

  std::ostream& operator<<(std::ostream& os) const;

 private:
  void addNewObject(Object* newObj);
  void removeObject(Object* obj);

  bool static indvTest(const std::string& line, Block* parent,
		       double expectedResult, bool examineResult);

  std::list<Object*> _postfixLine; //!< A post-fix list of operands/operators

  std::set<Object*> _objsToDelete; /**!< A list of objects this obj should
                                    *    delete. Operators are reused, so
                                    *    we don't have to delete those.
                                    *    Variables are deleted by the
                                    *    block that owns them, so we don't
                                    *    have to delete those either.
                                    */

  Block* _parent; //!< The Block of code that this Line belongs to

  ScalarNumber<double>* _tempHolder; /**!< Used to allocate memory needed for
                                      *    temporaries at compile time. This
                                      *    boosts runtime performance.
                                      */

  int _tempSize; //!< The size of the _tempHolder array

  bool _needReturnVal; /**!< Tells us if this Line needs to return result of
                        *    its execution. This is necessary because sometimes
                        *    Lines represent condition statements for loops or
                        *    ConditionalBlocks.
                        */

  int _curr; //The line number of this line
};

}

#endif
