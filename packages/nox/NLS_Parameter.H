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

//! Manipulating single parameters
class NLS_Parameter {

public:

  NLS_Parameter();
  NLS_Parameter(bool value);
  NLS_Parameter(int value);
  NLS_Parameter(double value);
  NLS_Parameter(const string& value);

  ~NLS_Parameter();

  void setValue(bool value);
  void setValue(int value);
  void setValue(double value);
  void setValue(const string& value);

  bool isBool() const;
  bool isInt() const;
  bool isDouble() const;
  bool isString() const;

  bool getBoolValue();
  int getIntValue();
  double getDoubleValue();
  string& getStringValue();

  bool isUsed() const;

  ostream& leftshift(ostream& stream) const;

private:

  enum NLS_ParameterType { NONE, BOOL, INT, DOUBLE, STRING };
  NLS_ParameterType type;
  bool bval;
  int ival;
  double dval;
  string sval;
  bool isused;
};

ostream& operator<<(ostream& stream, const NLS_Parameter& e);

#endif
