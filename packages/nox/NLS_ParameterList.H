// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_PARAMETERLIST_H
#define _NLS_PARAMETERLIST_H

#include<string>
#include<vector>
#include<map>

#include "NLS_Parameter.H"

//! Manipulating lists of parameters
class NLS_ParameterList {

public:

  //! Constructor
  NLS_ParameterList();

  //! Deconstructor
  ~NLS_ParameterList();

  //! List unused parameters
  void unused() const;

  /** @name Setting Parameters 
   *
   * Sets different types of parameters. The type depends on the
   * second entry. Be sure to use static_cast<type>() when the type is
   * ambiguous. Both char* and string map to are stored as strings
   * internally. Sets the parameter as "unused". */
  //@{
  void setParameter(const string& name, bool value);
  void setParameter(const string& name, int value);
  void setParameter(const string& name, double value);
  void setParameter(const string& name, const char* value);
  void setParameter(const string& name, const string& value);

  //! Note that this only retains a reference to the list, not a copy.
  //! Returns false is the "value" cannot be added. 
  //! This happens when "value" contains "this".
  bool setParameter(const string& name, const List& value);
  //@}

  /** @name Getting Parameters 
   *
   * Get different types of parameters. The type depends on the second
   * entry. Returns the nominal value if that parameter has not been *
   * specified. Be sure to use static_cast<type>() when the type is *
   * ambiguous. Both char* and string map return string values.  Sets
   * the parameters as "used". */
  //@{
  bool getParameter(const string& name, bool nominal) const;
  int getParameter(const string& name, int nominal) const;
  double getParameter(const string& name, double nominal) const;
  const string& getParameter(const string& name, const char* nominal) const;
  const string& getParameter(const string& name, const string& nominal) const;
  const List& getParameter(const string& name, const List& nominal) const;
  //@}

  //! Printing 
  ostream& print(ostream& stream, int indent = 0) const;

private:

  //! Check to see if "l" or any of its sublists is "this"
  bool isRecursive(const List& l) const;


private:

  //! Parameter list
  map<string, NLS_Parameter> params;
 
  //! Used to create a string when the getParameter is called with a
  //! char* nominal value. A new string is created for each such
  //! argument. The whole group of strings is destroyed when this object
  //! is destroyed. This is really annoying, but I don't know a better 
  //! way.
  mutable vector<string> tmpstrings;
};

ostream& operator<<(ostream& stream, const NLS_ParameterList& p);

#endif


