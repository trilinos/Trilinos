// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_PARAMETER_ENTRY_H
#define NOX_PARAMETER_ENTRY_H

#include "NOX_Common.H"		// class data element (string)
#include "NOX_Parameter_Arbitrary.H"

namespace NOX {

//! %NOX %Parameter support.
namespace Parameter {

class List; // anonther parameter type (forward declaration)

//! Manipulating single parameters, including sublists.
class Entry {

public:
  
  //! Default Constructor
  Entry();
  
  //! Copy constructor
  Entry(const Entry& source);

  //! Copy
  Entry& operator=(const Entry& source);

  //! Bool constructor
  Entry(bool value, bool isCreatedByGet = false);

  //! Integer constructor
  Entry(int value, bool isCreatedByGet = false);

  //! Double constructor
  Entry(double value, bool isCreatedByGet = false);

  //! String constructor (creates its own copy of the string)
  Entry(const string& value, bool isCreatedByGet = false);

  //! Arbitrary constructor (creates its own copy of the Arbitrary object)
  Entry(const Arbitrary& value, bool isCreatedByGet = false);

  //! Destructor
  ~Entry();

  /** @name ParameterList parameters
   *
   * Functions for handling parameters that are themselves lists.  */
  //@{
  List& setList(bool isCreatedByGet = false);
  List& getListValue();
  const List& getListValue() const;
  //@}

  /** @name Set functions. 
   
    The input value type determines the type of parameter
    stored. Invalidates any previous values stored by this object,
    although it doesn't necessarily erase them. Resets 'isused'
    functionality. 
  */
  //@{ 
  void setValue(bool value, bool isCreatedByGet = false);
  void setValue(int value, bool isCreatedByGet = false);
  void setValue(double value, bool isCreatedByGet = false);
  void setValue(const char* value, bool isCreatedByGet = false);
  void setValue(const string& value, bool isCreatedByGet = false);
  void setValue(const Arbitrary& value, bool isCreatedByGet = false);
  //@}

  /** @name Is functions. 
   
    Return true if the parameter is of the specified type; otherwise,
    return false.
  */
  //@{ 
  bool isBool() const;
  bool isInt() const;
  bool isDouble() const;
  bool isString() const;
  bool isList() const;
  bool isArbitrary() const;
  //@}

  
  //! Return whether or not the value is used; i.e., whether or not
  //! the value has been retrieved via a get function.
  bool isUsed() const;

  /** @name Get functions. 
   
    Returns value of parameter. The value is nonsense if we do not
    request the correct type of value. We cannot name all of these
    functions the same since the language does not allow us to
    overload functions based solely on return value. 
  */
  //@{ 
  bool getBoolValue() const;
  int getIntValue() const;
  double getDoubleValue() const;
  const string& getStringValue() const;
  const Arbitrary& getArbitraryValue() const;
  //@}

  //! Output the parameter to the given stream. 
  /*! 
    Formats the output as "<type,value>", except in the case of a
    list which just outputs "\<sublist\>". If the parameter has not yet
    been set, it outputs "\<NONE\>". This is the function called by the
    ostream operator<<. 
  */
  ostream& leftshift(ostream& stream) const;

private:

  //! Reset the entry
  void reset();
  
  //! All possible parameter types that this class can store
  enum EntryType { 
    //! No entry type set yet (will be set later by setValue()
    NOX_NONE, 
    //! Boolean
    NOX_BOOL, 
    //! Integer
    NOX_INT, 
    //! Double
    NOX_DOUBLE, 
    //! String
    NOX_STRING,
    //! NOX::Parameter::Arbitrary
    NOX_ARBITRARY,
    //! Sublist (NOX::Parameter::List)
    NOX_LIST 
  };

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

  //! Pointer to Arbitrary object, is this is of type ARBITRARY
  Arbitrary* aval;

  //! Pointer to list, if this is of type LIST
  List* lval;		

  //! Has this parameter been accessed by a "get" function?
  mutable bool isGotten;

  //! Was this parameter a nominal value assigned by a "get" function?
  mutable bool isSetByGet;

};

} // namespace Parameter
} // namespace NOX

//! Output the parameter. Relies of leftshift operator defined in the class.
ostream& operator<<(ostream& stream, const NOX::Parameter::Entry& e);

#endif
