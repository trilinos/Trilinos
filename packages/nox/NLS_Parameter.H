// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_PARAMETER_H
#define _NLS_PARAMETER_H

#include<string>
#include<vector>
#include<map>

//forward declaration
class NLS_ParameterList;
typedef NLS_ParameterList List;

//! Manipulating single parameters
class NLS_Parameter {

public:

  NLS_Parameter();
  NLS_Parameter(bool value);
  NLS_Parameter(int value);
  NLS_Parameter(double value);
  NLS_Parameter(const string& value);
  NLS_Parameter(const List& p);

  ~NLS_Parameter();

  void setValue(bool value);
  void setValue(int value);
  void setValue(double value);
  void setValue(const string& value);
  void setValue(const List& p);

  bool isBool() const;
  bool isInt() const;
  bool isDouble() const;
  bool isString() const;
  bool isList() const;

  bool getBoolValue() const;
  int getIntValue() const;
  double getDoubleValue() const;
  const string& getStringValue() const;
  const List& getListValue() const;

  bool isUsed() const;

  ostream& leftshift(ostream& stream) const;

private:

  enum NLS_ParameterType { NONE, BOOL, INT, DOUBLE, STRING, LIST };
  NLS_ParameterType type;
  bool bval;
  int ival;
  double dval;
  string sval;
  const List* lval;		
  mutable bool isused;

};

ostream& operator<<(ostream& stream, const NLS_Parameter& e);

#endif
