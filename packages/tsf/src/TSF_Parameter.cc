
#include "TSF_Parameter.h"

namespace TSF {

//=============================================================================
Parameter::Parameter(void) 
{
  initializeDefaults();
}
//=============================================================================
Parameter::Parameter(const TSF::Parameter& parameter) 
{
  assign(parameter);
}

//=============================================================================
void Parameter::assign(const TSF::Parameter& parameter) 
{
  label_ = 0; // Must initialize this specially here  

  initializeDefaults(parameter.getLabel());
  
  // Call each get function.  One of them will eventually work!
  if (parameter.isChar()) parameter.getChar(char_);
  if (parameter.isCharString()) parameter.getCharString(charString_);
  if (parameter.isInt()) parameter.getInt(int_);
  if (parameter.isIntArray()) parameter.getIntArray(length_, intArray_);
  if (parameter.isDouble()) parameter.getDouble(double_);
  if (parameter.isDoubleArray()) parameter.getDoubleArray(length_, doubleArray_);

  return;
}

//=============================================================================
Parameter::Parameter(char const * const & label, char parameter)  
{
  set(label, parameter);
}
//=============================================================================
Parameter::Parameter(char const * const & label, char * parameter)  
{
  set(label, parameter);
}
//=============================================================================
Parameter::Parameter(char const * const & label, int parameter)  
{
  set(label, parameter);
}
//=============================================================================
Parameter::Parameter(char const * const & label, int length, int * parameter)  
{
  set(label, length, parameter);
}
//=============================================================================
Parameter::Parameter(char const * const & label, double parameter)  
{
  set(label, parameter);
}
//=============================================================================
Parameter::Parameter(char const * const & label, int length, double * parameter)  
{
  set(label, length, parameter);
}

//=============================================================================
Parameter::~Parameter(void)  
{

  if (label_!=0) delete [] label_;
  if (charString_!=0) delete [] charString_;
  if (intArray_!=0) delete [] intArray_;
  if (doubleArray_!=0) delete [] doubleArray_;
  
}
//=============================================================================
void Parameter::initializeDefaults(char const * const & label)
{
  if (label_!=0) delete [] label_;
  initializeDefaults();
  int length = strlen(label);
  if (length>0) label_ = new char[length];
  strcpy(label_, label);
  return;
}
//=============================================================================
void Parameter::initializeDefaults(void)
{
  label_ = 0;
  char_ = 0;
  charString_ = 0;
  int_ = 0;
  intArray_ = 0;
  double_ = 0.0;
  doubleArray_ = 0;
  isChar_ = false;
  isCharString_ = false;
  isInt_ = false;
  isIntArray_ = false;
  isDouble_ = false;
  isDoubleArray_ = false;
  length_ = 0;
  return;
}
//=============================================================================
Parameter::set(char const * const & label, char parameter)  
{
  initializeDefaults(label);
  char_ = parameter;
  isChar_ = true;
}
//=============================================================================
Parameter::set(char const * const & label, char * parameter)  
{
  initializeDefaults(label);
  length_ = strlen(parameter);
  if (length_>0) charString_ = new char[length_];
  strcpy(charString_, parameter);
  isCharString_ = true;
}
//=============================================================================
Parameter::set(char const * const & label, int parameter)  
{
  initializeDefaults(label);
  int_ = parameter;
  isInt_ = true;
}
//=============================================================================
Parameter::set(char const * const & label, int length, int * parameter)  
{
  initializeDefaults(label);
  if (length_>0) intArray_ = new int[length_];
  for (int i=0; i<length_; i++) intArray_[i] = parameter[i];
  isIntArray_ = true;
}
//=============================================================================
Parameter::set(char const * const & label, double parameter)  
{
  initializeDefaults(label);
  double_ = parameter;
  isDouble_ = true;
}
//=============================================================================
Parameter::set(char const * const & label, int length, double * parameter)  
{
  initializeDefaults(label);
  if (length_>0) doubleArray_ = new double[length_];
  for (int i=0; i<length_; i++) doubleArray_[i] = parameter[i];
  isDoubleArray_ = true;
}

} // TSF namespace
