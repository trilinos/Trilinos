
#include "TSF_ParameterList.h"

namespace TSF {
//=============================================================================
ParameterList::ParameterList(void) 
{
  initializeDefaults();
}
//=============================================================================
ParameterList::ParameterList(int numParameters, const TSF::Parameter * & parameters) 
{
  initializeDefaults();
  for (int i=0; i<numParameters; i++) addParameter(parameters[i]);
  
}
//=============================================================================
ParameterList::ParameterList(const TSF::ParameterList& ParameterList) 
{
  initializeDefaults();
  int numParameters = ParameterList.getNumParameters();
  for (int i=0; i<numParameters; i++) addParameter(ParameterList.parameters_[i]);
}

//=============================================================================
ParameterList::~ParameterList(void)  
{
  deleteArrays();
  numParameters_ = 0;
}
//=============================================================================
void ParameterList::initializeDefaults(void)
{
  maxNumParameters_ = 50;
  numParameters_ = 0;
  parameters_ = new TSF::Parameter[maxNumParameters_];
  return;
}
//=============================================================================
void ParameterList::deleteArrays(void)  
{
  //for (int i=0; i<getNumParameters(); i++) delete parameters_[i];
  delete [] parameters_;
}
//=============================================================================
int ParameterList::setKeyLoc(char const * const & key)  
{
  // First see if this key is already in the list
  int keyLoc = getKeyLoc(key);
  if (KeyLoc>-1) return(keyLoc); // Return if it is
  if (numParameters_ >= maxNumParameters_) {

    int tmp_maxNumParameters = maxNumParameters_ + 10;
    TSF::Parameter * tmp_parameters = new TSF::Parameter[tmp_maxNumParameters];
    for (int i=0; i<numParameters_; i++) tmp_parameters[i] = parameters_[i];
    deleteArrays();
    parameters_ = tmp_parameters;
    maxNumParameters_ = tmp_maxNumParameters;
  }
  keyLoc = numParameters_++;
  return(keyLoc);
}
//=============================================================================
void ParameterList::setParameter(const TSF::Parameter parameter)  
{
  parameters_[setKeyLoc()] = parameter;
}
//=============================================================================
void ParameterList::setParameter(char const * const & key, char parameter)  
{
  parameters_[setKeyLoc()].set(key,parameter);
}
//=============================================================================
void ParameterList::setParameter(char const * const & key, char * parameter)  
{
  parameters_[setKeyLoc()].set(key,parameter);
}
//=============================================================================
void ParameterList::setParameter(char const * const & key, int parameter)  
{
  parameters_[setKeyLoc()].set(key,parameter);
}
//=============================================================================
void ParameterList::setParameter(char const * const & key, int length, int * parameter)  
{
  parameters_[setKeyLoc()].set(key,length,parameter);
}
//=============================================================================
void ParameterList::setParameter(char const * const & key, double parameter)  
{
  parameters_[setKeyLoc()].set(key,parameter);
}
//=============================================================================
void ParameterList::setParameter(char const * const & key, int length, double * parameter)  
{
  parameters_[setKeyLoc()].set(key,length,parameter);
}
//=============================================================================
bool ParameterList::hasParameter(char const * const & key)  const {
  int i = getKeyLoc(key);
  return (i>-1);
}
//=============================================================================
bool ParameterList::hasChar(char const * const & key) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].isChar());
  else return(false);
}
//=============================================================================
bool ParameterList::hasCharString(char const * const & key, char * & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].CharString(parameter));
  else return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, int & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].getInt(parameter));
  else return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, int & length, int * & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].getIntArray(length, parameter));
  else return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, double & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].getDouble(parameter));
  else return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, int & length, double * & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].getDoubleArray(length, parameter));
  else return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, TSF::Parameter & parameter)  const {
  int i = getKeyLoc(key);
  if (i>-1) {
    parameter = parameters_[i];
    return(true);
  }

  return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, char & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].getChar(parameter));
  else return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, char * & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].getCharString(parameter));
  else return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, int & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].getInt(parameter));
  else return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, int & length, int * & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].getIntArray(length, parameter));
  else return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, double & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].getDouble(parameter));
  else return(false);
}
//=============================================================================
bool ParameterList::getParameter(char const * const & key, int & length, double * & parameter) const {
  int i = getKeyLoc(key);
  if (i>-1) return(parameters_[i].getDoubleArray(length, parameter));
  else return(false);
}
//=============================================================================
int ParameterList::getKeyLoc(char const * const & key) const { 

  // Note:  This is a very quick and dirty search.  Shouldn't be a problem for a small list of parameters

  for (int i=0; i<numParameters_; i++) if (strcmp(key,parameters_[i].getLabel())==0) return(i);

  return(-1);
}
}// TSF namespace
