#include <iostream.h>

//
//Epetra_ESI.h includes all of the Epetra_ESI_* headers, which in turn
//include the Epetra headers that they depend on.
//
#include "Epetra_ESI.h"

//========= utility function prototypes ========

int create_ESI_Stuff(int numGlobal, int numLocal, Epetra_Comm& comm,
		     esi::IndexSpace<int>*& ispace,
                     esi::Vector<double,int>*& vec);

//========= main program =======================

int main(int argc, char** argv) {
   Epetra_Comm* comm = NULL;

#ifdef EPETRA_MPI
   if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      cerr << "error in MPI_Init." << endl; return(-1);
   }
   comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
   comm = new Epetra_SerialComm();
#endif

   if (comm == NULL) {
      cerr << "error allocating Epetra_Comm." << endl; return(-1);
   }

   int numProcs = comm->NumProc();
   int localProc = comm->MyPID();
   int numLocal = 50;
   int numGlobal = numProcs*numLocal;

   //This program is a simple test driver for Epetra's implementation of the
   //abstract esi::Vector interface. Epetra's implementation is called 
   //epetra_esi::Vector.

   //We'll first declare a index-space, because a vector needs to be
   //constructed with an index-space, and we'll need to destroy it at the
   //end of the program.
   //
   esi::IndexSpace<int>* esiIndexSpace = NULL;
   esi::Vector<double,int>* esiVec = NULL;

   //Now create the esi::IndexSpace, along with an esi::Vector.

   CHK_ERR( create_ESI_Stuff(numGlobal, numLocal, *comm, esiIndexSpace, esiVec) );

   //Let's use the 'put' function to put 4.0's in esiVec.
   CHK_ERR( esiVec->put(4.0) );

   //let's declare a couple more esi::Vectors.
   esi::Vector<double,int> *x = NULL, *y = NULL;

   CHK_ERR( esiVec->clone(x) );
   CHK_ERR( esiVec->clone(y) );

   //fill x with 2's and y with 3's.
   CHK_ERR( x->put(2.0) );
   CHK_ERR( y->put(3.0) );

   //form the dot product of x and y.
   double dotResult = 0.0;
   CHK_ERR( x->dot(*y, dotResult) );

   cout << "x dot y: " << dotResult << endl;

   //we know that the dot product should be 6 * numGlobal...
   if (dotResult != (double)(6*numGlobal)) return(-1);

   //let's create a vector z = x + y
   esi::Vector<double,int>* z = NULL;
   CHK_ERR( y->clone(z) );
   CHK_ERR( z->axpby(1.0, *x, 1.0, *y) );

   //z should be full of 5's now.

   double nrm1 = 0.0;
   CHK_ERR( z->norm1(nrm1) );

   cout << "1-norm of z: " << nrm1 << endl;

   //the 1-norm of z should be 5*numGlobal...
   if (nrm1 != (double)(5*numGlobal)) return(-1);

   //let's try out the epetra_esi::Vector constructor that accepts a
   //Epetra_Vector. First we'll cast esiVec.
   Epetra_Vector* petra_vec = dynamic_cast<Epetra_Vector*>(esiVec);
   if (petra_vec == NULL) return(-1);

   epetra_esi::Vector<double,int>* pesivec =
         new epetra_esi::Vector<double,int>(*petra_vec);

   //pesivec now holds a "view" of esiVec. No data copy should have occurred.
   //4.0's were put in esiVec above, so the 1-norm of pesivec should be
   //4*numGlobal...
   CHK_ERR( pesivec->norm1(nrm1) );
   if (nrm1 != (double)(4*numGlobal)) return(-1);

   //now, use the axpby function to make pesivec hold x + y
   CHK_ERR( pesivec->axpby(1.0, *x, 1.0, *y) );

   //finally, the 1-norm of pesivec should now be 5*numGlobal.
   CHK_ERR( pesivec->norm1(nrm1) );
   if (nrm1 != (double)(5*numGlobal)) return(-1);

   delete (esi::Object*)pesivec;
   delete esiVec;
   delete x;
   delete y;
   delete z;
   delete esiIndexSpace;
   delete comm;

#ifdef EPETRA_MPI
   MPI_Finalize();
#endif

   return(0);
}

//----------------------------------------------
int create_ESI_Stuff(int numGlobal, int numLocal, Epetra_Comm& comm,
		   esi::IndexSpace<int>*& esiIndexSpace, esi::Vector<double,int>*& esiVec)
{
  //the 0 means we're using an indexBase of 0.
  esiIndexSpace = new epetra_esi::IndexSpace<int>(numGlobal, numLocal, 0, comm);

  //return error-code -1 if the allocation failed
  if (esiIndexSpace == NULL) return(-1);

  epetra_esi::IndexSpace<int>* pispace = NULL;
  int err = esiIndexSpace->getInterface("epetra_esi::IndexSpace", (void*&)pispace);
  if (err != 0) return(-1);
 
  esiVec = new epetra_esi::Vector<double,int>(*pispace);
  if (esiVec == NULL) return(-1);

  return(0);
}

