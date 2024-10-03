// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _BLOCKRTC_H
#define _BLOCKRTC_H

#include "RTC_ExecutableRTC.hh"
#include "RTC_TokenizerRTC.hh"
#include "RTC_VariableRTC.hh"

#include <string>
#include <map>
#include <list>
#include <iostream>

namespace PG_RuntimeCompiler {

/**
 * The Block class represents a block of code. A block of code begins with
 * a { and ends with a }. Blocks of code can contain Lines of code and other
 * blocks of code which are its sub-blocks. Sub-blocks of code will have
 * access to all the variables that its parent has access to. Block extends
 * Executable because a Block of code can be executed.
 */

class Block: public Executable
{
 public:

  /**
   * Constructor -> Initializes instance variables
   *
   * @param vars - A list of vars that are already active (in scope) at the
   *               time this Block was created
   */
  Block(const std::map<std::string, Variable*>& vars);

  /**
   * Destructor -> The destructor will call the delete every statement in the
   *               Block. It will also delete all of the variables that were
   *               declared inside this Block.
   */
  virtual ~Block();

  /**
   * addStatement -> This method adds a statement object to statement list
   *
   * @param statement - The statement we are adding
   */
  void addStatement(Executable* statement);

  /**
   * addVariable -> This method adds the a variable to the variable list
   *
   * @param var - The variable we are adding
   */
  void addVariable(Variable* var);

  /**
   * getVar -> The method returns a variable with matching name.
   *           If no such variable exists, NULL is returned.
   *
   * @param name - The name of the variable we are looking for
   */
  Variable* getVar(const std::string& name);

  /**
   * createSubStatements -> Looks at each line until it sees the closing }. It
   *                        will create and add the sublines and substatements
   *                        of this block
   *
   * @param lines - The lines of code
   * @param errs  - A string containing errors
   */
  void createSubStatements(Tokenizer& lines, std::string& errs);

  std::ostream& operator<<(std::ostream& os) const;

 private:

  std::list<Variable*> _varsIOwn; /**!< A list of variables that this block
                                   *    created and therefore is responsible
                                   *    for deleting.
                                   */

  static int indent;

 protected:

  std::list<Executable*> _statements; /**!< A list of executable objects
                                       *    contained by this block. These
                                       *    objects may be Lines of code or
                                       *    other Blocks. When a Block of code
                                       *    is executed, it will execute all
                                       *    the objects in its _statements list
                                       */

  std::map<std::string, Variable*> _vars; /**!< The map of available variables
                                           *    for this Block. When a variable
                                           *    is declared, it's added to this
                                           *    list. It maps names to Variable
                                           *    objects.
                                           */
};

}
#endif
