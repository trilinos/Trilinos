#ifndef _fei_Data_hpp_
#define _fei_Data_hpp_

#include <string.h>
#include <stdlib.h>

/**
  This is a very simple class for passing stuff around
  in a void pointer. It has the ability to store and query
  a type name, so at least there can be user-enforced
  type safety.

  When setTypeName is called, a char* is created and a copy
  of the input argument is taken. This char* is later destroyed
  by the Data destructor. The void* dataPtr_ member is not
  destroyed, it is just a copy of a pointer.
*/

class Data {
 public:
  /** Default constructor. */
   Data() {typeName_ = NULL; dataPtr_ = NULL;};

   /** Default destructor. */
   virtual ~Data() {if (typeName_) delete [] typeName_;};

   /** Set a string representing the type of the object stored in
       'getDataPtr()'. */
   void setTypeName(const char* name) {if (typeName_) delete [] typeName_;
                                 size_t len = strlen(name);
                                 typeName_ = new char[len+1];
                                 strcpy(typeName_, name);
                                 typeName_[len] = '\0';};

   /** Query the string representing the type of the object stored in
       'getDataPtr()'. */
   char* getTypeName() const {return(typeName_);};

   /** Set the contents of the data pointer. */
   void setDataPtr(void* ptr) {dataPtr_ = ptr;};

  /** Retrieve the contents of the data pointer. */
   void* getDataPtr() const {return(dataPtr_);};

 private:
   char* typeName_;
   void* dataPtr_;
};

#endif

