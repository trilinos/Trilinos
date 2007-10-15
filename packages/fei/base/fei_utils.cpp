/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include "fei_fstream.hpp"

#include "fei_utils.hpp"
#include "snl_fei_Utils.hpp"
#include "fei_MatrixReducer.hpp"
#include "fei_Matrix_Impl.hpp"
#include "fei_LinearSystemCore.hpp"
#include "fei_ParameterSet.hpp"

#include "fei_version.h"

#ifdef FEI_HAVE_TIME_H
#include <time.h>
#endif

//----------------------------------------------------------------------------
const char* fei_VERSION::version()
{
  static int times_called = 0;

  static std::string static_fei_version_string;

  if (times_called == 0) {
    FEI_OSTRINGSTREAM osstr;

    osstr << FEI_MAJOR_VERSION << "."
          << FEI_MINOR_VERSION << "."
          << FEI_PATCH_VERSION;

    static_fei_version_string = osstr.str();

    times_called = 1;
  }

  return( static_fei_version_string.c_str() );
}

//----------------------------------------------------------------------------
double fei::utils::cpu_time()
{
  double cpu_seconds = 0.0;

#ifdef FEI_HAVE_TIME_H
  cpu_seconds = (1.0*clock())/CLOCKS_PER_SEC;
#endif

  return(cpu_seconds);
}

//----------------------------------------------------------------------------
LinearSystemCore*
fei::utils::get_LinearSystemCore(fei::Matrix* matrix)
{
  fei::Matrix* matptr = matrix;
  fei::MatrixReducer* matred = dynamic_cast<fei::MatrixReducer*>(matptr);
  if (matred != NULL) matptr = matred->getTargetMatrix().get();

  fei::Matrix_Impl<LinearSystemCore>* mat_lsc =
    dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*>(matptr);

  if (mat_lsc != NULL) {
    return mat_lsc->getMatrix().get();
  }

  return(NULL);
}

//----------------------------------------------------------------------------
void fei::utils::char_ptrs_to_strings(int numStrings,
                                     const char*const* charstrings,
                                     std::vector<std::string>& stdstrings)
{
  stdstrings.resize(0);
  for(int i=0; i<numStrings; ++i) {
    if (charstrings[i] != NULL) {
      std::string tempstr(charstrings[i]);
      stdstrings.push_back(tempstr);
    }
  }
}

//----------------------------------------------------------------------------
void fei::utils::strings_to_char_ptrs(std::vector<std::string>& stdstrings,
                                      int& numStrings,
                                      const char**& charPtrs)
{
  numStrings = stdstrings.size();
  charPtrs = numStrings > 0 ? new const char*[numStrings] : NULL;

  for(int i=0; i<numStrings; ++i) {
    charPtrs[i] = stdstrings[i].c_str();
  }
}

//----------------------------------------------------------------------------
void fei::utils::parse_strings(std::vector<std::string>& stdstrings,
                              const char* separator_string,
                              fei::ParameterSet& paramset)
{
  std::vector<std::string>::iterator
    s_iter = stdstrings.begin(),
    s_end = stdstrings.end();

  int intval = 0;
  double doubleval = 0.0;
  const char* charstring_key = NULL;
  const char* charstring_val = NULL;
  int key_len = 0;
  int val_len = 0;

  for(; s_iter != s_end; ++s_iter) {
    snl_fei::separate_string((*s_iter).c_str(), separator_string,
                    charstring_key, key_len, charstring_val, val_len);
    if (key_len == 0) {
      continue;
    }

    std::string keystr(charstring_key, key_len);

    if (val_len == 0) {
      fei::Param vparam(keystr.c_str(), true);
      paramset.add(vparam, false);
      continue;
    }

    std::string valstr(charstring_val);

    if (valstr == "true" || valstr == "True" || valstr == "TRUE") {
      fei::Param bparam(keystr.c_str(), true);
      paramset.add(bparam, false);
      continue;
    }

    if (valstr == "false" || valstr == "False" || valstr == "FALSE") {
      fei::Param bparam(keystr.c_str(), false);
      paramset.add(bparam, false);
      continue;
    }

    FEI_ISTRINGSTREAM isstr(valstr);

    //Does charstring_val contain a floating-point value?
    //If so, we'll store it as a double.
    std::string::size_type i = valstr.find(".");
    std::string::size_type valstrsize = valstr.size();

    if (i < valstrsize) {
      isstr >> doubleval;
      if (!isstr.fail()) {
        fei::Param dparam(keystr.c_str(), doubleval);
        paramset.add(dparam, false);
        continue;
      }
      isstr.clear();
    }

    //Does charstring_val contain an int?
    isstr >> intval;
    if (!isstr.fail()) {
      fei::Param iparam(keystr.c_str(), intval);
      paramset.add(iparam, false);
      continue;
    }
    isstr.clear();

    //If charstring_val doesn't contain an int or a real, we'll just store
    //it as a string.
    fei::Param sparam(keystr.c_str(), charstring_val);
    paramset.add(sparam, false);
  }
}

//----------------------------------------------------------------------------
void
fei::utils::convert_ParameterSet_to_strings(const fei::ParameterSet* paramset,
                                        std::vector<std::string>& paramStrings)
{
  paramStrings.resize(0);

  fei::ParameterSet::const_iterator
    iter = paramset->begin(),
    iter_end = paramset->end();

  for(; iter != iter_end; ++iter) {
    const fei::Param& param = *iter;
    fei::Param::ParamType ptype = param.getType();

    FEI_OSTRINGSTREAM osstr;
    osstr << param.getName();

    switch(ptype) {
    case fei::Param::STRING:
      osstr << " " << param.getStringValue();
      break;
    case fei::Param::DOUBLE:
      osstr << " " << param.getDoubleValue();
      break;
    case fei::Param::INT:
      osstr << " " << param.getIntValue();
      break;
    case fei::Param::VOID:
      break;
    case fei::Param::BOOL:
      if (param.getBoolValue()) osstr << " true";
      else osstr << " false";
      break;
    default:
      break;
    }

    paramStrings.push_back(osstr.str());
  }
}

