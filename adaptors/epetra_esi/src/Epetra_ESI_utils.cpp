#ifndef _Epetra_ESI_utils_cpp_
#define _Epetra_ESI_utils_cpp_

bool epetra_esi::stringsMatch(const char* str1, const char* str2)
{
  int len1 = strlen(str1);   int len2 = strlen(str2);

  if (len1 != len2) return(false);

  if (strncmp(str1, str2, len1) == 0) return(true);
  else return(false);
}

int epetra_esi::findString(Epetra_Array<const char*>& strings,
                                 const char* string)
{
  int index = -1;
  for(int i=0; i<strings.length(); i++) {
    if (stringsMatch(strings[i], string)) { index = i; break; }
  }

  return(index);
}

bool epetra_esi::hasSubString(const char* string, const char* sub)
{
  int len1 = strlen(string);  int len2 = strlen(sub);
  if (len1 < len2) return(false);

  char* strptr = (char*)string;
  for(int i=0; i<(len1-len2); i++) {
    if (strncmp(strptr, sub, len2) == 0) return(true);
    strptr++;
  }

  return(false);
}

int epetra_esi::findHasSubString(Epetra_Array<const char*>& strings,
                                       const char* string)
{
  int index = -1;
  for(int i=0; i<strings.length(); i++) {
    if (epetra_esi::hasSubString(strings[i], string)) { index = i; break; }
  }

  return(index);
}

#endif

