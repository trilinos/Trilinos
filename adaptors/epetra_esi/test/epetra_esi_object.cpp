#include <iostream.h>

#include "Epetra_ESI.h"

//========= utility function prototype ========

int create_ESI_Object(esi::Object*& esi_obj);


//========= main program =======================

int main(int argc, char** argv) {

#ifdef EPETRA_MPI
   if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      cerr << "error in MPI_Init." << endl; return(-1);
   }
#endif

   //This program is a simple test driver for Epetra's implementation
   //of the esi::Object interface.
   //First, instantiate the object. Note that the function create_ESI_Object
   //returns an instance of the abstract ESI interface esi::Object, but the
   //run-time type of the object is epetra_esi::Object.
   //
   esi::Object* esiObj = NULL;
   CHK_ERR( create_ESI_Object(esiObj) );

   //An esi::Object can store a collection of "RunTimeModels", which in
   //practice are usually MPI communicators. For testing purposes, let's just
   //store some named NULL pointers, and make sure we can retrieve them by
   //name later.
   //
   CHK_ERR( esiObj->setRunTimeModel("test1", NULL) );
   CHK_ERR( esiObj->setRunTimeModel("test2", (void*)0x002) );
   CHK_ERR( esiObj->setRunTimeModel("test3", NULL) );

   //The esi::Argv interface is a very simple encapsulation of a
   //collection of strings.
   //
   esi::Argv* argvObj = new epetra_esi::Argv();
   CHK_ERR( esiObj->getRunTimeModelsSupported(argvObj) );

   //At this point 'argvObj' holds the names of all of the RunTimeModels
   //we stored in 'esiObj' above.
   //
   if (argvObj->getArgCount() != 3) {
      cerr << "esiObj->getRunTimeModelsSupported returned wrong number."
	   << endl;
      return(-1);
   }

   //Print out those names...
   //
   for(int i=0; i<argvObj->getArgCount(); i++) {
      cout << "rtNames["<<i<<"]: " << argvObj->get(i) << endl;
   }

   //Query for a particular RunTimeModel.
   void* rtmodel = NULL;
   CHK_ERR( esiObj->getRunTimeModel("test2", rtmodel) );

   if (rtmodel != (void*)0x002) {
      cerr << "esiObj->getRunTimeModel output the wrong pointer." << endl; 
      return(-1);
   }

   delete argvObj;
   delete esiObj;

   esiObj = NULL;

   //Now let's instantiate a concrete epetra_esi::Object and make sure that its
   //getInterface function works, by requesting an abstract esi::Object
   //interface from it.
   //
   epetra_esi::Object pobj;

   CHK_ERR( pobj.getInterface("esi::Object", (void*&)esiObj) );

   CHK_ERR( esiObj->setRunTimeModel("BOGUS", NULL) );

#ifdef EPETRA_MPI
   cout << "calling MPI_Finalize."<<endl;
   MPI_Finalize();
#endif
   return(0);
}

//----------------------------------------------
int create_ESI_Object(esi::Object*& esi_obj)
{
   esi_obj = new epetra_esi::Object;

   //return error-code -1 if the allocation failed
   if (esi_obj == NULL) return(-1);

   return(0);
}

