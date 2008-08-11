/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_ParameterSet.hpp>

fei::ParameterSet::ParameterSet()
  : params_(NULL)
{
  params_ = new std::vector<const Param*>;
}

fei::ParameterSet::~ParameterSet()
{
  const_iterator iter = begin(), iter_end = end();
  for(; iter != iter_end; ++iter) {
    delete &(*iter);
  }

  delete params_; params_ = NULL;
}

int fei::ParameterSet::getIntParamValue(const char* name,
					int& paramValue) const
{
  const fei::Param* param = get(name);
  if (param == NULL) return(-1);

  if (param->getType() != fei::Param::INT) return(-1);

  paramValue = param->getIntValue();
  return(0);
}

int fei::ParameterSet::getDoubleParamValue(const char* name,
					   double& paramValue) const
{
  const fei::Param* param = get(name);
  if (param == NULL) return(-1);

  if (param->getType() == fei::Param::INT) {
    paramValue = param->getIntValue();
  }
  else if (param->getType() == fei::Param::DOUBLE) {
    paramValue = param->getDoubleValue();
  }
  else return(-1);

  return(0);
}

int fei::ParameterSet::getStringParamValue(const char* name,
					   std::string& paramValue) const
{
  const fei::Param* param = get(name);
  if (param == NULL) return(-1);

  if (param->getType() != fei::Param::STRING) return(-1);

  paramValue = param->getStringValue();
  return(0);
}

int fei::ParameterSet::getBoolParamValue(const char* name,
					 bool& paramValue) const
{
  const fei::Param* param = get(name);
  if (param == NULL) return(-1);

  if (param->getType() != fei::Param::BOOL) return(-1);

  paramValue = param->getBoolValue();
  return(0);
}

int fei::ParameterSet::getVoidParamValue(const char* name,
					 const void*& paramValue) const
{
  const fei::Param* param = get(name);
  if (param == NULL) return(-1);

  if (param->getType() != fei::Param::VOID) return(-1);

  paramValue = param->getVoidValue();
  return(0);
}
