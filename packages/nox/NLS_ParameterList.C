#include "NLS_ParameterList.H"

//----------------------------------------------------------------------
// NLS_ParameterList
//----------------------------------------------------------------------

// A useful typedef for manipulating the params map
typedef map<string, NLS_Parameter>::const_iterator paramsiter;

NLS_ParameterList::NLS_ParameterList() {}

NLS_ParameterList::~NLS_ParameterList() 
{
  unused();
}

void NLS_ParameterList::unused() const
{
  // In deconstructor, warn about any unused parameters
  for (paramsiter it = params.begin(); it != params.end(); it ++) {
    if (!(it->second.isUsed())) {
      cout << "Parameter \"" << it->first << "\" " << it->second
	   << " is unused" << endl;
    }
  }
}

void NLS_ParameterList::setParameter(const string& name, bool value)
{
  params[name].setValue(value);
}

void NLS_ParameterList::setParameter(const string& name, int value)
{
  params[name].setValue(value);
}

void NLS_ParameterList::setParameter(const string& name, double value)
{
  params[name].setValue(value);
}

void NLS_ParameterList::setParameter(const string& name, const char* value)
{
  params[name].setValue(value);
}

void NLS_ParameterList::setParameter(const string& name, const string& value)
{
  params[name].setValue(value);
}

bool NLS_ParameterList::getParameter(const string& name, bool nominal)
{
  if ((params.find(name) == params.end()) || 
      (!(params[name].isBool())))
    return nominal;
  
  return params[name].getBoolValue();
}

int NLS_ParameterList::getParameter(const string& name, int nominal)
{
  if ((params.find(name) == params.end()) ||   
      (!(params[name].isInt())))
    return nominal;
  
  return params[name].getIntValue();
}

double NLS_ParameterList::getParameter(const string& name, double nominal)
{
  if ((params.find(name) == params.end()) ||   
      (!(params[name].isDouble())))
    return nominal;
  
  return params[name].getDoubleValue();
}

const string& NLS_ParameterList::getParameter(const string& name, const char* nominal)
{
  tmpstrings.push_back(nominal);
  return getParameter(name, tmpstrings[tmpstrings.size() - 1]);
}

const string& NLS_ParameterList::getParameter(const string& name, const string& nominal)
{
  if ((params.find(name) == params.end()) ||   
      (!(params[name].isString())))
    return nominal;

  return params[name].getStringValue();
}
  
