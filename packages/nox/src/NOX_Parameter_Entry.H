// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_PARAMETER_ENTRY_H
#define _NOX_PARAMETER_ENTRY_H

#include<string>	    // one parameter type

namespace NOX {

//! NOX Parameter support.
namespace Parameter {

class List; // anonther parameter type (forward declaration)

//! Manipulating single parameters, including sublists.
class Entry {

public:
  
  //! Default Constructor
  Entry();

  //! Bool constructor
  Entry(bool value);

  //! Integer constructor
  Entry(int value);

  //! Double constructor
  Entry(double value);

  //! String constructor
  Entry(const string& value);

  //! Destructor
  ~Entry();

  /** @name ParameterList parameters
   *
   * Functions for handling parameters that are themselves lists.  */
  //@{
  List& setList();
  List& getListValue();
  const List& getListValue() const;
  //@}

  /** @name Set functions. 
   *
   * The input value type determines the type of parameter
   * stored. Invalidates any previous values stored by this object,
   * although it doesn't necessarily erase them. Resets 'isused'
   * functionality. */
  //@{ 
  void setValue(bool value);
  void setValue(int value);
  void setValue(double value);
  void setValue(const char* value);
  void setValue(const string& value);
  //@}

  /** @name Is functions. 
   *
   * Return true if the parameter is of the specified type; otherwise,
   * return false. */
  //@{ 
  bool isBool() const;
  bool isInt() const;
  bool isDouble() const;
  bool isString() const;
  bool isList() const;
  //@}

  
  /** @name Get functions. 
   *
   * Returns value of parameter. The value is nonsense if we do not
   * request the correct type of value. We cannot name all of these
   * functions the same since the language does not allow us to
   * overload functions based solely on return value. */
  //@{ 
  bool getBoolValue() const;
  int getIntValue() const;
  double getDoubleValue() const;
  const string& getStringValue() const;
  //@}

  //! Return whether or not the value is used; i.e., whether or not
  //! the value has been retrieved via a get function.
  bool isUsed() const;

  //! Output the parameter to the given stream. 
  /*! Formats the output as "<type,value>", except in the case of a
    list which just outputs "<sublist>". If the parameter has not yet
    been set, it outputs "<NONE>". This is the function called by the
    ostream operator<<. */
  ostream& leftshift(ostream& stream) const;

private:
  
  //! All possible parameter types that this class can store
  enum EntryType { NONE, BOOL, INT, DOUBLE, STRING, LIST };

  //! Type of parameter stored in this object.
  EntryType type;

  //! Boolean value, if this is of type BOOL
  bool bval;

  //! Integer value, if this is of type INT
  int ival;

  //! Double value, if this is of type DOUBLE
  double dval;

  //! String value, if this is of type STRING
  string sval;

  //! Pointer to list, if this is of type LIST
  List* lval;		

  //! Has this parameter been accessed by a "get" function?
  mutable bool isused;

};

}
}

//! Output the parameter. Relies of leftshift operator defined in the class.
ostream& operator<<(ostream& stream, const NOX::Parameter::Entry& e);

#endif
