#ifndef _VARIABLERTC_H
#define _VARIABLERTC_H

#include <assert.h>
#include <string>
#include "commonRTC.hh"
#include "ValueRTC.hh"

namespace PG_RuntimeCompiler {

/**
 * A Varaible object represents the variables in the code the user gives us. 
 */

class Variable: public Value
{
 public:

  /**
   * Constructor -> Trivial
   *
   * @param name    - The name of the variable
   * @param type    - The data type of the variable
   * @param objType - The object type of the variable
   */
  Variable(const std::string& name, Type type, ObjectType objType) 
    : Value(type, objType) {
    _name    = name;
    _address = NULL;

    _isSet                 = false;
    _willBeInitAtTimeOfUse = false;
  }
  
  /**
   * getName -> This method returns the name of the variable
   */
  std::string getName() const {return _name;}

  /** 
   * setAddress -> This method sets the memory address of a variable. This 
   *               only applies to Variables that have been passed in to 
   *               Function by-reference.
   *
   * @param addr - The new address
   */
  virtual void setAddress(void* addr) = 0;

  /**
   * setSize - Applies to arrays only, no-op here
   */
  virtual void setSize(int size) {}
  
  /**
   * evaluateSizeExpr - Applies to arrays only, no-op here
   */
  virtual void evaluateSizeExpr() {}

  /** 
   * isSet -> This method returns _isSet. It is only used for argument 
   *          variables.
   */
  bool isSet() const {return _isSet;};
  
  /** 
   * init -> This method sets _willBeInitAtTimeOfUse to true
   */ 
  void init() { _willBeInitAtTimeOfUse = true;}

  /** 
   * isInit -> This method returns _willBeInitAtTimeOfUse
   */ 
  bool isInit() const { return _willBeInitAtTimeOfUse;}

 protected:

  std::string _name; //!< The name of the variable

  bool _isSet; /**!< Tells us if an argument has been set yet. If someone tries
                *    to run a Function that has not had all its arguments set, 
                *    it will generate an error.
                */

  bool _willBeInitAtTimeOfUse; /**!< helps us find errors where the user is 
                                *    trying to use an uninitialized variable.
                                */

  void* _address; //!< The address location of the variable
};

}

#endif
