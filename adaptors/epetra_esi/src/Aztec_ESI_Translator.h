#ifndef _Aztec_ESI_Translator_h_
#define _Aztec_ESI_Translator_h_

namespace aztecoo_esi {

class Translator {
 public:
  Translator();
  ~Translator() {}

  /** Function for translating string parameters to Aztec options
      and params settings.
     @param numParams Input, number of parameter strings.
     @param paramStrings Input, list of strings.
     @param options Input/Output, user-allocated Aztec options array (of length
               AZ_OPTIONS_SIZE).
     @param params Input/Output, user-allocated Aztec params array (of length
               AZ_PARAMS_SIZE).
     @param quiet Determines whether a warning should be printed when a
              string is encountered that is not successfully translated to
              an Aztec option or param value. By default, a warning is printed.
     @return error-code 0 if successful.
  */
  static int stringsToAztecSettings(int numParams, char** paramStrings,
                                    int* options, double* params,
                                    bool quiet=false);

  static Epetra_Array<const char*>& get_az_def_map();

  static Epetra_Array<const char*>& get_az_option_strs();

  static Epetra_Array<const char*>& get_az_param_strs();

  static int azOptionValue(const char* valStr,
                           Epetra_Array<const char*>& azDefStrings,
                           bool quiet=false);

  static float azParamValue(const char* valStr,
                            Epetra_Array<const char*>& azDefStrings,
                            bool quiet=false);
};

}; //close aztecoo_esi namespace

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Aztec_ESI_Translator.cpp"
#endif

#endif

