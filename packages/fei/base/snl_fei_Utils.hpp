/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_Utils_hpp_
#define _snl_fei_Utils_hpp_

#include <fei_fwd.hpp>
#include <fei_mpi.h>
#include <fei_SharedPtr.hpp>

namespace fei {
  class CSVec;
  class Matrix;
}

#include <vector>
#include <map>
#include <string>

namespace snl_fei {

  /** Given a search key, and a list of strings, search the list of strings
      looking for a string that starts with the search key. i.e., a string
      such that key is a leading substring. If found, return a pointer to
      that string.<br>
      Important Note: The returned pointer should be treated as read-only. It
      is not a separate copy of the string being pointed to.

      @param key String to be searched for.
      @param numParams Number of strings to be searched.
      @param paramStrings List of strings to be searched.

      @return search-result Pointer to the entry in paramStrings that has key
      as a leading substring, or NULL if key is not found.
  */
  const char* getParam(const char* key,
		       int numParams,
		       const char* const* paramStrings);

  /** Given a search key, and a list of strings containing key-value pairs,
      search the list of strings looking for one that starts with the search
      key. i.e., a string such that key is a leading substring. If found,
      return a pointer to the value portion of that key-value pair.<br>
      Important Note: The returned pointer should be treated as read-only. It
      is not a separate copy of the string being pointed to.

      @param key String to be searched for.
      @param numParams Number of strings to be searched.
      @param paramStrings List of strings to be searched.
      @param separator Optional argument, defaults to ' ' (space). This is
      the character that is the separator between keys and values in the
      parameter strings.
      @return search-result Pointer to the entry in paramStrings that has key
      as a leading substring, or NULL if key is not found.
  */
  const char* getParamValue(const char* key,
			    int numParams,
			    const char* const* paramStrings,
			    char separator=' ');

  /** Given a search key, and a vector of strings containing key-value pairs,
      search the strings looking for one that starts with the search
      key. i.e., a string such that key is a leading substring. If found,
      return a pointer to the value portion of that key-value pair.<br>
      Important Note: The returned pointer should be treated as read-only. It
      is not a separate copy of the string being pointed to.

      @param key String to be searched for.
      @param params vector of strings to be searched.
      @param separator Optional separator char, defaults to space
      @return search-result Pointer to the entry in paramStrings that has key
      as a leading substring, or NULL if key is not found.
  */
  const char* getParamValue(const char* key,
			    std::vector<std::string>& params,
			    char separator=' ');

  /** Given a search key, and a list of strings containing key-value pairs,
      search the list of strings looking for one that starts with the search
      key. i.e., a string such that key is a leading substring. If found,
      return a pointer to the value portion of that key-value pair.<br>
      Important Note: The returned pointer should be treated as read-only. It
      is not a separate copy of the string being pointed to.

      @param key String to be searched for.
      @param numParams Number of strings to be searched.
      @param paramStrings List of strings to be searched.
      @param foundOffset offset at which key was found in paramStrings. If not
      found, then this parameter is set to -1.
      @param separator Optional argument, defaults to ' ' (space). This is
      the character that is the separator between keys and values in the
      parameter strings.
      @return search-result Pointer to the entry in paramStrings that has key
      as a leading substring, or NULL if key is not found.
  */
  const char* getParamValue(const char* key,
			    int numParams,
			    const char* const* paramStrings,
			    int& foundOffset,
			    char separator=' ');

  /** Get the integer value of a named key-value pair from an array of strings.
   */
  int getIntParamValue(const char* key,
		       int numParams,
		       const char*const* params,
		       int& paramValue);

  /** Get the double-precision value of a named key-value pair from an array of
      strings.
   */
  int getDoubleParamValue(const char* key,
			  int numParams,
			  const char*const* params,
			  double& paramValue);

  /** Get the double-precision value of a named key-value pair from a vector of
      strings.
   */
  int getDoubleParamValue(const char* key,
			  std::vector<std::string>& params,
			  double& paramValue);

  /** Given a search key, and a list of strings, search the list of strings
      looking for a string that starts with the search key. i.e., a string
      such that key is a leading substring. If found, return a pointer to
      that string, and also provide the offset of the location in
      'paramStrings' at which it was found.<br>
      Important Note: The returned pointer should be treated as read-only. It
      is not a separate copy of the string being pointed to.

      This method is case-sensitive.

      @param key String to be searched for.
      @param numParams Number of strings to be searched.
      @param paramStrings List of strings to be searched.
      @param foundOffset offset at which key was found in paramStrings. If not
      found, then this parameter is set to -1.
      @return search-result Pointer to the entry in paramStrings that has key
      as a leading substring, or NULL if key is not found.
  */
  const char* getParam(const char* key,
		       int numParams,
		       const char* const* paramStrings,
		       int& foundOffset);

  /** Given a search key, and a vector of strings, search the strings
      looking for a string that starts with the search key. i.e., a string
      such that key is a leading substring. If found, return a pointer to
      that string, and also provide the offset of the location in
      'paramStrings' at which it was found.<br>
      Important Note: The returned pointer should be treated as read-only. It
      is not a separate copy of the string being pointed to.

      This method is case-sensitive.

      @param key String to be searched for.
      @param paramStrings vector of strings to be searched.
      @param foundOffset offset at which key was found in paramStrings. If not
      found, then this parameter is set to -1.
      @return search-result Pointer to the entry in paramStrings that has key
      as a leading substring, or NULL if key is not found.
  */
  const char* getParam(const char* key,
		       std::vector<std::string>& paramStrings,
		       int& foundOffset);

  /** stored a named void pointer in a vector. (poor-man's map) */
  int storeNamedAttribute(const char* name,
			  void* attribute,
			  std::vector<char*>& attributeNames,
			  std::vector<void*>& attributes);

  /** retrived a named void pointer from a vector. (poor-man's map) */
  void* retrieveNamedAttribute(const char* name,
			       std::vector<char*>& attributeNames,
			       std::vector<void*>& attributes);

  /** separate a string into pieces. unsupported, for power-users only */
  void separate_string(const char* input_string,
		       const char* substring,
		       const char*& before_substring,
		       int& len_before_substring,
		       const char*& after_substring,
		       int& len_after_substring);

  /** Given a string, return the length of the leading substring, which is the
      characters that precede any space, tab, equals sign, or null character.
  */
  unsigned leading_substring_length(const char* string);

  /** Given a parameter-string, which is assumed to contain a pair of
      substrings (key-value pair) separated by a specified separator character,
      return a pointer to the second substring. In other words, return a pointer
      to the first character after the separator. Allows for the possibility that
      the substrings are separated by more than one copy of the separator.

      @param paramString String to be searched.
      @param separator Optional parameter, defaults to ' ' (space).
      @return result First character after separator, or NULL if separator
      is not found in paramString.
  */
  const char* skipSeparator(const char* paramString,
			    char separator=' ');

  /** Merge a list of strings into a list of other strings. Maintain
      uniqueness, don't merge strings that are already present.
  */
  int mergeStringLists(char**& strings, int& numStrings,
		       const char*const* stringsToMerge, int numStringsToMerge);

  /** Resolve conflicts between constraint-relations and essential (dirichlet)
      boundary conditions.
  */
  int resolveConflictingCRs(fei::MatrixGraph& matrixGraph,
			    fei::Matrix& bcEqns,
                            const std::vector<int>& bcEqnNumbers);

  /** Do appropriate communications to gather column-portions of remotely-held
      essential BCs onto local processor.
  */
  int gatherRemoteEssBCs(fei::CSVec& essBCs,
			 fei::SparseRowGraph* remoteGraph,
			 fei::Matrix& matrix);

  fei::SharedPtr<fei::SparseRowGraph>
    mergeSparseRowGraphs(const fei::SparseRowGraph* srg1,
                         const fei::SparseRowGraph* srg2);

  /** Copy values from a 2D "C-style" block-diagonal table (list of pointers) to
      a "fortran-style" flat column-contiguous array.
  */
  void copy2DBlockDiagToColumnContig(int numBlocks,
				     const int* blockSizes,
				     const double*const* values2d,
				     int format,
				     double* colcontigvalues);

  /** Copy values from a 2D "C-style" table (list of pointers) to a "fortran-style"
      flat column-contiguous array.
  */
  void copy2DToColumnContig(int numrows,
			    int numcols,
			    const double*const* values2d,
			    int format,
			    double* colcontigvalues);
} //namespace snl_fei

#endif // _snl_fei_Utils_hpp_

