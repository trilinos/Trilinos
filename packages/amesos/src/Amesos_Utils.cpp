#include "Amesos_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "Amesos_Utils.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"

void Amesos_BreakForDebugger(Epetra_Comm& Comm)
{
  char hostname[80];
  char buf[80];
  if (Comm.MyPID()  == 0) cout << "Host and Process Ids for tasks" << endl;
  for (int i = 0; i <Comm.NumProc() ; i++) {
    if (i == Comm.MyPID() ) {
#ifdef COUGAR
      sprintf(buf, "Host: %s   PID: %d", "janus", getpid());
#else
      gethostname(hostname, sizeof(hostname));
      sprintf(buf, "Host: %s\tComm.MyPID(): %d\tPID: %d",
	      hostname, Comm.MyPID(), getpid());
#endif
      printf("%s\n",buf);
      fflush(stdout);
      sleep(1);
    }
  }
  if(Comm.MyPID() == 0) {
    printf("\n");
    printf("** Pausing to attach debugger...\n");
    printf("** You may now attach debugger to the processes listed above.\n");
    printf( "**\n");
    printf( "** Enter a character to continue > "); fflush(stdout);
    char go;
    scanf("%c",&go);
  }

  Comm.Barrier();

}

//============================================================================
Epetra_CrsMatrix* CreateOverlappingCrsMatrix(Epetra_CrsMatrix* Matrix,
					     const int OverlappingLevel)
{

  if (OverlappingLevel == 0) 
    return(0); // All done
  if (Matrix->Comm().NumProc() == 1) 
    return(0); // All done

  Epetra_CrsMatrix* OverlappingMatrix;
  OverlappingMatrix = const_cast<Epetra_CrsMatrix*>(Matrix);
  Epetra_Map* OverlappingMap;
  OverlappingMap = const_cast<Epetra_Map*>(&(Matrix->RowMap()));

  Epetra_CrsMatrix* OldMatrix;
  const Epetra_Map* DomainMap = &(Matrix->DomainMap());
  const Epetra_Map* RangeMap = &(Matrix->RangeMap());

  for (int level = 1; level <= OverlappingLevel ; ++level) {

    OldMatrix = OverlappingMatrix;

    Epetra_Import* OverlappingImporter;
    OverlappingImporter = const_cast<Epetra_Import*>(OldMatrix->Importer());
    int NumMyElements = OverlappingImporter->TargetMap().NumMyElements();
    int* MyGlobalElements = OverlappingImporter->TargetMap().MyGlobalElements();

    // need to build an Epetra_Map in this way because Epetra_CrsMatrix
    // requires Epetra_Map and not Epetra_BlockMap

    OverlappingMap = new Epetra_Map(-1,NumMyElements,MyGlobalElements,
				    0, Matrix->Comm());

    if (level < OverlappingLevel)
      OverlappingMatrix = new Epetra_CrsMatrix(Copy, *OverlappingMap, 0);
    else
      // On last iteration, we want to filter out all columns except 
      // those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlappingMatrix = new Epetra_CrsMatrix(Copy, *OverlappingMap, 
					       *OverlappingMap, 0);

    OverlappingMatrix->Import(*OldMatrix, *OverlappingImporter, Insert);
    if (level < OverlappingLevel) {
      OverlappingMatrix->TransformToLocal(DomainMap, RangeMap);
    }
    else {
      OverlappingMatrix->TransformToLocal(DomainMap, RangeMap);
    }

    delete OverlappingMap;

    if (level > 1) {
      delete OldMatrix;
    }
  }

  return(OverlappingMatrix);
}

