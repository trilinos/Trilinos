#ifndef _Epetra_ESI_utils_h_
#define _Epetra_ESI_utils_h_

namespace epetra_esi {

  /** check whether two strings are the same. */
  bool stringsMatch(const char* str1, const char* str2);


  /** if 'string' is contained in 'strings', return its index. return -1 if it
      is not present.
  */
  int findString(Epetra_Array<const char*>& strings, const char* string);

  /** check whether one string occurs as a substring of another string.
    @param string Input, the string to be searched.
    @param sub Input, the substring to be searched for.
    @return true if sub occurs within string, false if it doesn't.
  */
  bool hasSubString(const char* string, const char* sub);

  /** if 'sub' is a substring of any string in 'strings', return the index of
      that string.  return -1 if sub is not a substring of any of them.
  */
  int findHasSubString(Epetra_Array<const char*>& strings, const char* strng);

}; //namespace epetra_esi

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_ESI_utils.cpp"
#endif

#endif

