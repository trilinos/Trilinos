// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef AMESOS_PARAMETER_ENTRY_H
#define AMESOS_PARAMETER_ENTRY_H

//  #include "NOX_Common.h"		// class data element (string)
#include <string>
#include <vector>
#include "Teuchos_ParameterList.hpp"

//#include "Amesos_Parameter_Arbitrary.h"

namespace AMESOS {

#ifdef OLDLIST

//! %AMESOS %Parameter support.
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
    list which just outputs "<sublist>". If the parameter has not yet
    been set, it outputs "<NONE>". This is the function called by the
    ostream operator<<. 
  */
  ostream& leftshift(ostream& stream) const;

private:

  //! Reset the entry
  void reset();
  
  //! All possible parameter types that this class can store
  enum EntryType { 
    //! No entry type set yet (will be set later by setValue()
    AMESOS_NONE, 
    //! Boolean
    AMESOS_BOOL, 
    //! Integer
    AMESOS_INT, 
    //! Double
    AMESOS_DOUBLE, 
    //! String
    AMESOS_STRING,
    //! AMESOS::Parameter::Arbitrary
    AMESOS_ARBITRARY,
    //! Sublist (AMESOS::Parameter::List)
    AMESOS_LIST 
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

#endif

} // namespace AMESOS

#ifdef OLD_LIST
//! Output the parameter. Relies of leftshift operator defined in the class.
ostream& operator<<(ostream& stream, const AMESOS::Parameter::Entry& e);
#endif

#endif
