// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _EXECUTABLERTC_H
#define _EXECUTABLERTC_H

#include <iostream>

namespace PG_RuntimeCompiler {

class Value;

/**
 * The Executable class is the parent of all classes of object that can
 * be executed (run). They are executed by calling the execute() method
 * which they must implement.
 */
class Executable
{
 public:

  /**
   * Destructor -> Is a no-op.
   */
  virtual ~Executable() {}

  /**
   * execute -> A pure virtual method. It will run the code inside of an
   *            executable object.
   */
  virtual Value* execute() = 0;

  virtual std::ostream& operator<<(std::ostream& os) const = 0;
};

std::ostream& operator<<(std::ostream& os, const Executable& obj);

}

#endif
