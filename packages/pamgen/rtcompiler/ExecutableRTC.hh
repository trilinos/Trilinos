#ifndef _EXECUTABLERTC_H
#define _EXECUTABLERTC_H

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

  virtual void print() = 0;
};

}

#endif
