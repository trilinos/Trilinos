#ifndef _FUNCTIONRTC_H
#define _FUNCTIONRTC_H

#include <string>
#include <typeinfo>
#include <vector>
#include "commonRTC.hh"
#include "BlockRTC.hh"
#include "VariableRTC.hh"
#include "ExecutableRTC.hh"

/**************** IMPORTANT ***********************
 * FunctionRTC.h is the only file you need to #include
 * in order to use the RTC tool. The methods below are 
 * the interface you will use RTC with.
 * Please read README.txt for more information.
 *************************************************/

namespace PG_RuntimeCompiler {

/**
 * A Function object represents the function that you have written and
 * compiled. It contains all the lines of code, and all the arguments 
 * and variables. A Function will not execute unless it is compiled 
 * successfully first.
 */

class Function 
{
 public:

  /**
   * Constructor
   *
   * @param varCount - The number of args this function will take (optional, used to reserve memory
   *                   and also used for backwards compatibility with a version which required
   *                   an integer as the first argument of the constructor).
   * @param name     - The name of the function (optional)
   */
  Function(const unsigned varCount=0,const std::string & name = "NONAME");

  /**
   * Destructor -> If the Function is not clean, the destructor deletes 
   *               _mainBlock, all of the arguments, and safely destroys _vars.
   */
  ~Function();


  /**
   * getName -> Returns the name of this Function.
   */
  std::string getName() const { return _name; }

  /** 
   * addVar -> Adds a variable to the argument list for this function. Note, 
   *           you should remember the order in which you add your variables 
   *           because otherwise you will have no idea what index to give to 
   *           the filler functions.   
   *
   * @param type - The type of the variable being added
   * @param name - The name of the variable being added
   */
  bool addVar(const std::string& type, const std::string& name);

  /**
   * addBody -> This method compiles the string body. Body should not include 
   *            any code related to the signature of the function.
   */
  bool addBody(const std::string& body);

  /** 
   * varValueFill -> This method fills the indexth variable of the signature 
   *                 with a value. Note!: the order of the variables in the 
   *                 signature is equal to the order in which the variables 
   *                 were added using addVar.
   *
   * @param index - The index of the argument we are filling
   * @param value - The value we are setting the argument to
   */
  bool varValueFill(unsigned int index, double value);

  /** 
   * varAddrFill -> This method is overloaded to take a pointer from any of the
   *                four types as the second argument. This method has similar 
   *                functionality to varValueFill except use varAddrFill when 
   *                you want to pass a value by reference instead of by copy.\
   *
   * @param index   - The index of the argument we are filling
   * @param address - The address of the variable
   * @param size    - Only used when calling on an array, its the array's size
   */  
  template <class T>
  bool varAddrFill(unsigned int index, T* address, int size = 0) {
    if (index >= _vars.size()) {
      _errors += "Index is too large.\n";
      return false;
    }

    checkType(index, size, address, _errors);
    if (_errors != "") return false;
    
    _vars[index]->setAddress(address);
    _vars[index]->setSize(size);
    _vars[index]->init();
    return true; 
  }

  /**
   * arrayAddrFill -> Left in for backwards compatibility, simply forwards the
   *                  call the varAddrFill.
   */
  template<class T>
  bool arrayAddrFill(unsigned int index, T* address, int size) {
    return varAddrFill(index, address, size);
  }

  /**
   * getErrors -> This method returns any compiler errors that have occured
   */
  std::string getErrors() const { return _errors;}

  /** 
   * getValueOfVar -> This methods searches the function for a variable named 
   *                  var, if it finds this variable, it returns its value. 
   *                  Note, this method can be used for arguments AND for 
   *                  variables declared inside the function.
   *
   * @param var - The variable we want the value of
   */
  double getValueOfVar(const std::string& var);

  /** 
   * execute -> This method executes your function. Make sure you have added 
   *            and filled all the arguments before calling this.
   */
  bool execute();

  /** 
   * cleanup -> This method cleans up all the heap memory created by 
   *            the program.
   */
  bool cleanup();

 private:

  /**
   * type_err -> This method creates a type error.
   */ 
  std::string type_err(unsigned int index) const {
    return "Wrong type passed to varAddrFill for variable: " +
      _vars[index]->getName() + " at index: " + intToString(index) + " \n " 
      "Note: both the basic type and scalar/nonscalar modifier must match up \n";
  }

  void checkType(unsigned int index, int size, int* addr, std::string& errs);
  void checkType(unsigned int index, int size, float* addr, std::string& errs);
  void checkType(unsigned int index, int size,double* addr, std::string& errs);
  void checkType(unsigned int index, int size, char* addr, std::string& errs);

  std::string _name; //!< The name of this Function

  std::string _errors; //!< String containing errors generated by the compiler.

  Block* _mainBlock; /**!< The main block of the Function, this is where 
                      *    execution begins.
                      */

  std::vector<Variable *> _vars; //!< Vector of pointers to Variables.

  bool _firstRun; /**!< A bool that tells us if the Function has been run yet 
                   *    or not.
                   */

  bool _compiled; /**!< A bool that tells us if the Function has been 
                   *    successfully compiled. 
                   */

  bool _clean; /**!< Tells us if the memory newed by the Function has been 
                *    cleaned up yet.
                */
};

}
  
#endif
