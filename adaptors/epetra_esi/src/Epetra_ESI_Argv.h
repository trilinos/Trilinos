#ifndef _Epetra_ESI_Argv_h_
#define _Epetra_ESI_Argv_h_

#include "Epetra_ESI_CHK_ERR.h"
#include "Epetra_ESI_platforms.h"

#include "Epetra_Array.h"

#include "esi/ESI.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

namespace epetra_esi {

/** Simple container-of-string-arguments class. This is basically an
    alternative to storing a collection of strings as a char** with an
    accompanying 'int numParams' descriptor. This class allows strings
    to be added one at a time, retrieved one at a time, and there's a
    function for querying how many strings are currently stored. That's
    about it. (STL vector could be used to provide the same functionality...)
*/

class Argv : public virtual esi::Argv {
 public:
  /** Default constructor. */
  Argv();

  /** Destructor. */
  virtual ~Argv(){ for(int i=0; i<argv_.length(); i++) delete[] argv_[i]; };

  /** @return the ith element of argv (not to be modified) if exists.
   *       returns 0 otherwise. */
  const char * get(int index)
    { if (index < 0 || index >= argv_.length()) return(NULL);
      return(argv_[index]);
    }

  /** @return the number of strings. */ 
  int getArgCount() { return( argv_.length() ); }

  /** Set a string on this Argv object. This object will take
      a copy of the string argument. The caller can destroy the input argument
      after this function returns.
      @param str Input. 
      @return error-code 0 if successful
  */
  int appendArg(const char* str)
    {
      char* newstr = new char[strlen(str)+1];
      if (newstr == NULL) return(-1);
      sprintf(newstr, str);
      return(argv_.append(newstr));
    }

 private:
  Epetra_Array<char*> argv_;
};

}; // namespace epetra_esi

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_ESI_Argv.cpp"
#endif

#endif
