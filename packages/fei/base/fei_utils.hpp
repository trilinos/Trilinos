#ifndef _fei_utils_hpp_
#define _fei_utils_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_fwd.hpp>
#include <fei_version.h>

#include <Teuchos_ParameterList.hpp>

#include <string>
#include <vector>



namespace fei_VERSION {
//  Mangled 'version' function. The name of this namespace, 'fei_VERSION'
//  is a macro that is defined (in fei_version.h) to include the current
//  version number.
//  This fei_VERSION::version() function is not intended for public use.
//  There is another 'version' function below, which is for public use,
//  but internally it calls this function in the fei_VERSION namespace. This
//  prevents header-mismatch errors, where a user application accidentally
//  includes headers from a different fei version than the fei libraries that
//  are being linked.
//  (In that scenario, unresolved symbol errors will occur, since the value
//  of fei_VERSION in the headers won't match what's been compiled into the
//  library.)
const char* version();

}//namespace fei_VERSION



/** The fei namespace contains public functions, classes and interfaces.
*/
namespace fei {

/** The utils namespace contains general utility functions.
*/
namespace utils {

/** Return a const char-ptr containing the fei version.
*/
inline
const char* version()
{
  return( fei_VERSION::version() );
}

/** Return CPU time. To measure an elapsed time, take the difference
    between two returned values.
*/
double cpu_time();

/** Convert a string to an fei::OutputLevel enum value.
   Valid string values are strings that match one of the enum names 
   in fei_fwd.hpp.
   If an invalid string is given, then fei::NONE will be returned.
*/
fei::OutputLevel string_to_output_level(const std::string& str);

/** Attempt to extract a LinearSystemCore from a fei::Matrix.
   Returns NULL if unsuccessful.
*/
LinearSystemCore* get_LinearSystemCore(fei::Matrix* matrix);

/** Return element-node connectivity in a pair of arrays, as follows:
  The 'nodes' array holds all of the node-identifiers (for nodes connected
  to local elements). The 'elem_offsets' array holds offsets into the
  'nodes' array at which the nodes for a given element can be found.

  Thus:

   num-elems = elem_offsets.size()-1;
   nodes for i-th element lie in these positions:
     nodes[elem_offsets[i] .. elem_offsets[i+1]-1]
*/
void getConnectivityArrays(fei::MatrixGraph& matrixGraph,
                           std::vector<int>& nodes,
                           std::vector<int>& elem_offsets);

/** Given an integer length 'numStrings' and a list of pointers-to-char-pointer,
    wrap them in a std::vector of std::string objects.
*/
void char_ptrs_to_strings(int numStrings,
                         const char*const* charstrings,
                         std::vector<std::string>& stdstrings);

/** populate an array of raw char-ptrs with the 'c_str()' pointers
    from the specified std::string objects.
    The caller is responsible for deleting the char** array, but *NOT*
    the individual char* pointers in the array.
*/
void strings_to_char_ptrs(std::vector<std::string>& stdstrings,
                          int& numStrings,
                          const char**& charPtrs);

/** Populate a fei::ParameterSet object, taking input from the
    specified list of strings. Each string is assumed to contain
    a key-value pair, separated by the specified 'separator_string'.
*/
void parse_strings(std::vector<std::string>& stdstrings,
                  const char* separator_string,
                  fei::ParameterSet& paramset);

/** Convert the contents of a fei::ParameterSet object to
    a collection of strings. Each string will contain a space-separated
    key-value pair.
*/
void convert_ParameterSet_to_strings(const fei::ParameterSet* paramset,
                                     std::vector<std::string>& paramStrings);

}//namespace utils
}//namespace fei

#endif

