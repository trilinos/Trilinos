#ifndef _Epetra_ESI_Object_cpp_
#define _Epetra_ESI_Object_cpp_

#include "Epetra_Object.h"

#ifndef EPETRA_ESI_INCLUDE_IMPLEMENTATION
//we're going to turn on the flag to include implementations, but if it isn't
//already on, we'll need to turn it back off afterwards... is that clear?
#define Epet_ESI_Obj_UNDEF
#endif

#define EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_ESI_Object.h"

#ifdef Epet_ESI_Obj_UNDEF
#undef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#undef Epet_ESI_Obj_UNDEF
#endif

//------------------------------------------------------------------------------
epetra_esi::Object::Object()
    : refCount_(1), rtmNames_(0,1), rtmPtrs_(0,1), ifNames_(0,1), ifPtrs_(0,1)
{
   int err = addInterface("esi::Object", (void*)(((esi::Object*)this)) );
   if (err) {
      cerr << "Epetra_ESI_Object: ERROR in constructor. aborting."<<endl;
      abort();
   }
}

//------------------------------------------------------------------------------
epetra_esi::Object::~Object()
{
  for(int i=0; i<rtmNames_.length(); i++) delete [] rtmNames_[i];
  for(int j=0; j<ifNames_.length(); j++) delete [] ifNames_[j];
}

//------------------------------------------------------------------------------
esi::ErrorCode epetra_esi::Object::getInterface(const char* name,
				       void *& object)
{
  int index = findString(ifNames_, name);

  if (index >= 0) { object = ifPtrs_[index]; return(0); }

  EPETRA_ESI_ERR_BEHAVIOR(-1);
}

//------------------------------------------------------------------------------
esi::ErrorCode epetra_esi::Object::getInterfacesSupported(esi::Argv* list)
{
  int len = ifNames_.length();
  if (len == 0) return(0);

  int err = 0;
  for(int i=0; i<len; i++) err += list->appendArg( ifNames_[i] );

  EPETRA_ESI_ERR_BEHAVIOR(err);
}

//------------------------------------------------------------------------------
esi::ErrorCode epetra_esi::Object::setRunTimeModel(const char* name, void *comm)
{
  //check whether an RTModel of this name has already been stored.
  int index = findString(rtmNames_, name);

  if (index >= 0) { rtmPtrs_[index] = comm; return(0); }

  char* newName = new char[strlen(name)+1]; strcpy(newName, name);
  if (newName == NULL) return(-1);
  rtmNames_.append(newName);
  rtmPtrs_.append(comm);

  return(0);
}

//------------------------------------------------------------------------------
esi::ErrorCode epetra_esi::Object::getRunTimeModel(const char* name,
					  void*& comm)
{
  int index = findString(rtmNames_, name);

  if (index >= 0) { comm = rtmPtrs_[index]; return(0); }

  EPETRA_ESI_ERR_BEHAVIOR(-1);
}

//------------------------------------------------------------------------------
esi::ErrorCode epetra_esi::Object::getRunTimeModelsSupported(esi::Argv* list)
{
  int len = rtmNames_.length();
  if (len == 0) return(0);

  int err = 0;
  for(int i=0; i<len; i++) err += list->appendArg( rtmNames_[i] );

  EPETRA_ESI_ERR_BEHAVIOR(err);
}

//------------------------------------------------------------------------------
bool epetra_esi::Object::hasSubString(const char* string, const char* sub)
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

//------------------------------------------------------------------------------
bool epetra_esi::Object::stringsMatch(const char* str1, const char* str2)
{
  int len1 = strlen(str1);   int len2 = strlen(str2);

  if (len1 != len2) return(false);

  if (strncmp(str1, str2, len1) == 0) return(true);
  else return(false);
}

//------------------------------------------------------------------------------
int epetra_esi::Object::findString(Epetra_Array<const char*>& strings,
				 const char* string)
{
  int index = -1;
  for(int i=0; i<strings.length(); i++) {
    if (stringsMatch(strings[i], string)) { index = i; break; }
  }

  return(index);
}

//------------------------------------------------------------------------------
int epetra_esi::Object::findHasSubString(Epetra_Array<const char*>& strings,
				       const char* string)
{
  int index = -1;
  for(int i=0; i<strings.length(); i++) {
    if (hasSubString(strings[i], string)) { index = i; break; }
  }

  return(index);
}

//------------------------------------------------------------------------------
int epetra_esi::Object::addInterface(const char* ifName, void* ifPtr)
{
  int index = findString(ifNames_, ifName);
  if (index >= 0) { ifPtrs_[index] = ifPtr; return(0); }

  char* newName = new char[strlen(ifName)+1];
  if (newName == NULL) EPETRA_ESI_ERR_BEHAVIOR(-1);
  strcpy(newName, ifName);
  ifNames_.append(newName);
  ifPtrs_.append(ifPtr);
  return(0);
}

//------------------------------------------------------------------------------
void epetra_esi::Object::msgAbort(const char* msg)
{
  cerr << msg << ". aborting." << endl;
#ifdef EPETRA_MPI
  //Since we have MPI, let's try to get a communicator and bring the
  //whole group down.
  MPI_Comm* mpicomm;
  int err = getRunTimeModel("MPI", (void*&)mpicomm);

  //if we couldn't get the communicator, then just use regular abort. 
  //Otherwise, use MPI_Abort...
  if (err) abort(); else MPI_Abort(*mpicomm, -1);
#else
  //we don't have MPI, so just abort the regular way...
  abort();
#endif
}

#endif

