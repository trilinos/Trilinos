#ifndef _Aztec_ESI_Translator_cpp_
#define _Aztec_ESI_Translator_cpp_

#include "Aztec_ESI_Translator.h"

#include "Epetra_ESI_utils.h"

aztecoo_esi::Translator::Translator()
{
}

int aztecoo_esi::Translator::stringsToAztecSettings(
                            int numParams, char** paramStrings,
                            int* options, double* params,
                            bool quiet)
{
  //First, let's get those string-arrays in Epetra_Array objects...
  //(This is a very light-weight operation.)

  Epetra_Array<const char*>& azDefStrings = get_az_def_map();

  Epetra_Array<const char*>& azOptionStrings = get_az_option_strs();

  Epetra_Array<const char*>& azParamStrings = get_az_param_strs();

  //Now let's set up some work-strings...

  Epetra_Array<char> keyString(128), valString(128);
  char* keyStr = keyString.dataPtr();
  char* valStr = valString.dataPtr();

  //Now loop over the input parameters, and do the deciphering.
  for(int i=0; i<numParams; i++) {
    //
    //first see if this input parameter is a space-separated pair of strings.
    //if it isn't, we aren't interested in it.
    //
    int num = sscanf(paramStrings[i], "%s %s", keyStr, valStr);
    if (num < 2) continue;

    //Now we need to determine whether the key-string is an aztec-option
    //or an aztec-param. (Note the implicit assumption: it can't be both
    //an aztec-option AND an aztec-param. So, we jump out of this loop-
    //iteration if it's an aztec-option.)

    int optIndx = epetra_esi::findString(azOptionStrings, keyStr);
    if (optIndx >= 0) {
      options[optIndx] = azOptionValue(valStr, azDefStrings, quiet);
      continue;
    }

    int prmIndx = epetra_esi::findString(azParamStrings, keyStr);
    if (prmIndx >= 0) {
      params[prmIndx] = azParamValue(valStr, azDefStrings, quiet);
    }
  }

  return(0);
}

Epetra_Array<const char*>& aztecoo_esi::Translator::get_az_def_map()
{
  //
  //Note from ABW: The following is ugly, but for now I can't think of
  //a clean way to translate string parameters into Aztec options/params....
  //The string-lists need to be updated whenever changes are
  //made in az_aztec_defs.h.
  //
  //az_def_map.h contains a list of symbol-value pairs gleaned directly
  //from az_aztec_defs.h by a script that uses awk. This is what we'll use to
  //figure out that, for example, "AZ_gmres" has the value 1.
  //

#include "az_def_map.h"
  static Epetra_Array<const char*> azDefMap(num_def_strs, num_def_strs,
                                            (const char**)az_def_map);
  return(azDefMap);
}

Epetra_Array<const char*>& aztecoo_esi::Translator::get_az_option_strs()
{
  //
  //And here's the REALLY bad part: the order of the strings in az_option_strs
  //and az_param_strs matters, because their positions are used for the
  //index into the Aztec options_ and params_ arrays... In other words, don't
  //mess with these string lists.
  //Don't try this at home, folks.
  //
  static const char* az_option_strs[] = {
  "AZ_solver",
  "AZ_scaling",
  "AZ_precond",
  "AZ_conv",
  "AZ_output",
  "AZ_pre_calc",
  "AZ_max_iter",
  "AZ_poly_ord",
  "AZ_overlap",
  "AZ_type_overlap",
  "AZ_kspace",
  "AZ_orthog",
  "AZ_aux_vec",
  "AZ_reorder",
  "AZ_keep_info",
  "AZ_recursion_level",
  "AZ_print_freq",
  "AZ_graph_fill",
  "AZ_subdomain_solve",
  "AZ_init_guess"
  };
  int num_option_strs = 20;
  static Epetra_Array<const char*> azOptionStrs(num_option_strs,
                                                num_option_strs,
                                              (const char**)az_option_strs);
  return(azOptionStrs);
}

Epetra_Array<const char*>& aztecoo_esi::Translator::get_az_param_strs()
{
  static const char* az_param_strs[] = {
    "AZ_tol",
    "AZ_drop",
    "AZ_ilut_fill",
    "AZ_omega",
    "AZ_rthresh",
    "AZ_athresh",
    "AZ_weights"
  };
  int num_param_strs = 7;
  static Epetra_Array<const char*> azParamStrs(num_param_strs,
                                               num_param_strs,
                                        (const char**)az_param_strs);
  return(azParamStrs);
}

int aztecoo_esi::Translator::azOptionValue(const char* valStr,
                                       Epetra_Array<const char*>& azDefStrings,
                                       bool quiet)
{
  static char tmpStr[128];

  int optionValue;
  int num = sscanf(valStr, "%d", &optionValue);

  if (num == 0) {
    int indx = epetra_esi::findHasSubString(azDefStrings, valStr);
    if (indx >= 0) {
      num = sscanf(azDefStrings[indx], "%s %d", tmpStr, &optionValue);
    }
  }

  if (num > 0) {
    return(optionValue);
  }

  if (!quiet) {
    cerr << "aztecoo_esi::Translator WARNING: azOptionValue failed to"
       << " convert '"<< valStr << "' to an integer." << endl;
  }

  return(0);
}

float aztecoo_esi::Translator::azParamValue(const char* valStr,
                                       Epetra_Array<const char*>& azDefStrings,
                                       bool quiet)
{
  static char tmpStr[128];

  float paramValue;
  int num = sscanf(valStr, "%e", &paramValue);

  if (num == 0) {
    int indx = epetra_esi::findHasSubString(azDefStrings, valStr);
    if (indx >= 0) {
      num = sscanf(azDefStrings[indx], "%s %e", tmpStr, &paramValue);
    }
  }

  if (num > 0) {
    return(paramValue);
  }

  if (!quiet) {
    cerr << "aztecoo_esi::Solver warning: azParamValue failed to"
         << " convert '"<< valStr << "' to a double." << endl;
  }

  return(0.0);
}

#endif

