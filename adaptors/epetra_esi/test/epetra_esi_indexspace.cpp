#include <iostream.h>

#include "Epetra_ESI.h"

//========= utility function prototypes ========

int create_ESI_IndexSpace(int numGlobal, int numLocal,
                   Epetra_Comm& comm, esi::IndexSpace<int>*& esi_ispc);


//========= main program =======================

int main(int argc, char** argv) {

   //This is a simple test driver to exercise Epetra's epetra_esi::IndexSpace, which is
   //an implementation of the abstract esi::IndexSpace interface.
   //The epetra_esi::IndexSpace object is actually built on top of the
   //Epetra_Map class, which requires a Epetra_Comm object at construction time.
   //So first, create a Epetra_Comm.
   //
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

   //Set up a trivially simple decomposition of an equation space, where each
   //processor has 50 equations. This is the decomposition that will be
   //described by the index-space object.
   //
   int numProcs = comm->NumProc();
   int localProc = comm->MyPID();
   int numLocal = 50;
   int numGlobal = numProcs*numLocal;
   int localOffset = numLocal*localProc; //localOffset is first-local-eqn.

   //Now we're ready to instantiate an esi::IndexSpace<int> instance, the run-time type
   //of which is actually epetra_esi::IndexSpace.
   //
   esi::IndexSpace<int>* esiIndexSpace = NULL;
   CHK_ERR( create_ESI_IndexSpace(numGlobal, numLocal, *comm, esiIndexSpace) );

#ifdef EPETRA_MPI
   //If we're using MPI, let's test getting the MPI RunTimeModel from the esiIndexSpace
   //and using it to find out how many procs there are.

   MPI_Comm* mpicomm;
   int mpiNumProcs;
   CHK_ERR( esiIndexSpace->getRunTimeModel("MPI", (void*&)mpicomm) );
   if (MPI_Comm_size(*mpicomm, &mpiNumProcs) != MPI_SUCCESS) {
      cerr << "err in MPI_Comm_size." << endl; return(-1);
   }

   if (mpiNumProcs != numProcs) {
      cerr << "error, numProcs from esiIndexSpace->getRunTimeModel comm is wrong."
	   <<endl;
      return(-1);
   }
#endif

   //Can we get a regular (native) Epetra_Map from this esi::IndexSpace?
   Epetra_Map* petraMap = NULL;
   CHK_ERR( esiIndexSpace->getInterface("Epetra_Map", (void*&)petraMap) );

   int petraMapNumGlobal = petraMap->NumGlobalElements();
   cout << "petraMap->NumGlobalElements(): " << petraMapNumGlobal << endl;

   int esiIndSpcNumGlobal = 0;
   CHK_ERR( esiIndexSpace->getGlobalSize(esiIndSpcNumGlobal) );

   //We've now obtained the global problem size from a couple of different
   //interfaces, which should all be "windows" on the same object. Make sure
   //they gave the right answer.
   //
   if (esiIndSpcNumGlobal != numGlobal ||
       petraMapNumGlobal != numGlobal) {
     cerr << "globalSize inconsistency." << endl; return(-1);
   }

   //Make sure that the esi::MapPartition interface gives the right answer for
   //number-of-processors.
   //
   int esiIndSpcPartNumProcs = 0;
   CHK_ERR( esiIndexSpace->getGlobalPartitionSetSize(esiIndSpcPartNumProcs) );

   if (esiIndSpcPartNumProcs != numProcs) {
      cerr << "esi::IndexSpace returned wrong numProcs." << endl; return(-1);
   }

   //esi::IndexSpace can supply a list, of length number-of-processors, which
   //contains each processor's local offset (each processor's first local eqn).
   //
   int* globalOffsets = new int[esiIndSpcPartNumProcs];
   if (globalOffsets == NULL) {
      cerr << "allocation failure, globalOffsets." << endl; return(-1);
   }

   CHK_ERR( esiIndexSpace->getGlobalPartitionOffsets(globalOffsets) );

   //esi::IndexSpace's localPartitionRank is our processor's 'MPI rank'.
   int localRank = -1;
   CHK_ERR( esiIndexSpace->getLocalPartitionRank(localRank) );

   //Make sure that "our" entry in the globalOffsets list matches our previously
   //calculated value for localOffset.
   if (globalOffsets[localRank] != localOffset) {
     cerr << "esi::IndexSpace's globalOffsets are not correct." << endl;
     return(-1);
   }

   for(int i=0; i<esiIndSpcPartNumProcs; i++) {
      cout << "proc " << localRank << ", globalOffsets["<<i<<"]: "
          << globalOffsets[i] << endl;
   }

   delete [] globalOffsets;
   delete esiIndexSpace;
   delete comm;

#ifdef EPETRA_MPI
   MPI_Finalize();
#endif

   return(0);
}

//----------------------------------------------
int create_ESI_IndexSpace(int numGlobal, int numLocal,
                   Epetra_Comm& comm, esi::IndexSpace<int>*& esi_ispc)
{
   esi_ispc = new epetra_esi::IndexSpace<int>(numGlobal, numLocal, 0, comm);

   //return error-code -1 if the allocation failed
   if (esi_ispc == NULL) return(-1);

   return(0);
}

