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
