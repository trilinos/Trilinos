#ifndef _Epetra_ESI_Object_h_
#define _Epetra_ESI_Object_h_

#include <stdlib.h>
#include <string.h>

#include "Epetra_ESI_CHK_ERR.h"
#include "Epetra_ESI_platforms.h"

#include "Epetra_Array.h"

#include "esi/ESI.h"

#ifdef EPETRA_MPI
#include <mpi.h>
#endif

namespace epetra_esi {

/** Petra's ESI Object implementation.
This class implements this ESI interface:
<ul>
<li>   esi::Object
</ul>
*/

class Object : public virtual esi::Object 
{
 public:

  /** Default constructor. */
  Object();

  /** Destructor. */
  virtual ~Object();


  /** Given an interface name, output a pointer to an instance of that 
      interface, if available.
    @param name Name of the interface being requested.
    @param object Output, double-pointer contains a pointer to the requested
         interface.
    @return error-code 0 if successful, -1 if interface 'name' is not available.
  */
  virtual esi::ErrorCode getInterface(const char* name, void *& object);


   /** Obtain a list of the names of the interfaces currently stored by this
     object. The esi::Argv argument is allocated/owned by the caller. It is
     basically just a list of strings. See Epetra_ESI_Argv documentation.
     @param names Input Petra's implementation of esi::Argv is called
     Epetra_ESI_Argv. Caller is responsible for allocating/deleting this object.
     @return error-code 0 if successful, -1 if any error was encountered.
   */
   virtual esi::ErrorCode getInterfacesSupported(esi::Argv* list);


   /** Store a 'Run-Time-Model' pointer and its name. (This is typically
     intended to be a pointer to an MPI_Comm, but can be anything else that
     the user provides.) This function takes a copy of the name (allocates
     a char* and copies name into it), and stores the comm pointer. Note that
     a RTModel of a given name may only be stored once. E.g., if this object
     already has an MPI_Comm stored with the name "MPI", then trying to store
     another one using "MPI" will return an error. The Petra ESI objects 
     Epetra_ESI_IndexSpace, Epetra_ESI_Vector and Epetra_ESI_RowMatrix store an MPI_Comm
     automatically upon construction, using the name "MPI", if built with
     -DEPETRA_MPI.
     @param name Name by which the RTModel should be referred to.
     @param comm Pointer to the RTModel being stored.
     @return error-code 0 if successful, -1 if the allocation failed when 
        making a copy of name.
   */
   virtual esi::ErrorCode setRunTimeModel(const char* name, void *comm);


   /** Obtain a named 'Run-Time-Model'. 
     @param name Name of the RTModel being requested.
     @param comm Output, double pointer contains a pointer to the RTModel that
       was requested, if available. This argument is not referenced if no
       RTModel by that name has been stored in this object.
     @return error-code 0 if successful, -1 if named RTModel is not available.
   */
   virtual esi::ErrorCode getRunTimeModel(const char* name, void*& comm);


   /** Obtain a list of the names of the RTModels currently stored by this
     object. The ugly triple-pointer argument represents the address
     of a list of strings, each string being a char*. This function will 
     allocate the list of strings, which the caller is then responsible for
     deleting. The user should NOT delete the contents of those strings. This
     function is probably the ugliest function in any ESI interface, and needs
     to be changed/improved!
     @param numbersupported pointer-to-int. This function sets the int to be
       the number of RTModels supported. (This is the C 'pass-by-reference'
       mechanism, should be replaced by int-reference (int&).)
     @param names Triple-pointer. Address of double-pointer (allocated by this
       function if 'numbersupported' is greater than 0), which is a list of
       strings. Caller is responsible for deleting the double-pointer, but not
       the contents of the strings.
     @return error-code 0 if successful, -1 if any error was encountered.
   */
   virtual esi::ErrorCode getRunTimeModelsSupported(esi::Argv* list);


   /** Reference-counting mechanism.  Calling this function increments the
       object's internal reference-count. The reference-count starts at 1
       when this object is constructed. Users are not required to use this
       function, or its companion function deleteReference. However, if this
       mechanism is used, addReference and deleteReference MUST be called in
       matching pairs. */
   virtual esi::ErrorCode addReference() { ++refCount_; return(0);}

   /** Reference-counting mechanism.  Calling this function decrements the
       object's internal reference-count. If the reference-count becomes
       0, the 'this' pointer is deleted, destroying the object.
       Users are not required to use this
       function, or its companion function addReference. However, if this
       mechanism is used, addReference and deleteReference MUST be called in
       matching pairs. */
   virtual esi::ErrorCode deleteReference()
     { --refCount_; if (refCount_ <= 0) delete this; return(0);}

   /** Function for serializing an object's state. Should be over-ridden by
     classes that derive from this one. Not implemented here.
   */
   virtual esi::ErrorCode saveState(FILE* file, int* format)
   {(void)file; (void)format; return(-1);};

   /** Function for reading an object's state from a file. Should be
      over-ridden by classes that derive from this one. Not implemented here.
   */
   virtual esi::ErrorCode loadState(FILE* file)
   {(void)file; return(-1);};


 protected:

  /** check whether one string occurs as a substring of another string.
    @param string Input, the string to be searched.
    @param sub Input, the substring to be searched for.
    @return true if sub occurs within string, false if it doesn't.
  */
  static bool hasSubString(const char* string, const char* sub);


  //check whether two strings are the same.
  static bool stringsMatch(const char* str1, const char* str2);


  //if 'string' is contained in 'strings', return its index. return -1 if it
  //is not present.
  static int findString(Epetra_Array<const char*>& strings, const char* string);


  //if 'sub' is a substring of any string in 'strings', return the index of
  //that string.  return -1 if sub is not a substring of any of them.
  static int findHasSubString(Epetra_Array<const char*>& strings, const char* string);


  //add an interface that this object is capable of returning a pointer to.
  int addInterface(const char* ifName, void* ifPtr);


  /**
    A function for printing a message to cerr and then aborting. This is
    used by objects that derive from epetra_esi::Object, ONLY when they
    encounter an error in their constructor. In that case, there's no
    way to return an error code to the user. Also, the only kinds of things
    that get done in the constructors are things that must succeed in order
    to bring the object up in a usable state. I.e., we only abort when we
    REALLY have to, because it's a very unfriendly thing to do to a user.
  */
  void msgAbort(const char* msg);

  int refCount_;

 private:
  //RTModel-names, and interface-names.
  Epetra_Array<const char*> rtmNames_, ifNames_;

  //RTModel-pointers, and interface-pointers.
  Epetra_Array<void*> rtmPtrs_,  ifPtrs_;
};

}; //namespace epetra_esi

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_ESI_Object.cpp"
#endif

#endif

