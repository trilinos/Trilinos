// Epetra_Object Test routine

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Object.h"

int main(int argc, char *argv[]) {

  int ierr = 0;
#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Epetra_MpiComm comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

  Epetra_SerialComm comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


  // I'm alive !!!
  if (verbose) cout << comm <<endl;

  Epetra_Object obj;

  // Test Epetra_Object label and the method to get the label attribute
  char* ObjLabel = obj.Label();
  char* ObjLabel1 = "Epetra::Object";
  if (verbose) cout << endl << endl << "This should say " << ObjLabel1 << ": " << ObjLabel << endl << endl << endl;
  ierr = strcmp(ObjLabel1,ObjLabel);
  assert (ierr==0);

  // Test Epetra_Object SetLabel attribute set method
  char* NewObjLabel = "New name for Epetra_Object";
  obj.SetLabel(NewObjLabel);
  char* NewObjLabel1 = obj.Label(); 
  if (verbose) cout << endl << "This should say " << NewObjLabel << ": " << NewObjLabel1 << endl << endl << endl;
  ierr = strcmp(NewObjLabel1,NewObjLabel);
  assert (ierr==0);

  // Test GetRacebackMode and SetTracebackMode methods
  assert(obj.GetTracebackMode()==DefaultTracebackMode);
  if (verbose) cout << endl <<"Default Traceback Mode value = " << obj.GetTracebackMode() << endl;

  obj.SetTracebackMode(DefaultTracebackMode-1);
  if (verbose) cout << "Set Traceback Mode value to one less than default = " << obj.GetTracebackMode() << endl << endl;
  Epetra_Object obj0;
  assert(obj0.GetTracebackMode()==DefaultTracebackMode-1);

  // Test constructors other than the default
  Epetra_Object obj1(1); // pass only TracebackMode
  int TbM = obj1.GetTracebackMode();
  if (verbose) cout << endl << endl << "This should say 1: " << TbM << endl << endl;
  assert (1==TbM);

  Epetra_Object obj2(NewObjLabel); // pass only a label
  char* NewObjLabel2 = obj2.Label();
  if (verbose) cout << endl << endl << "This should say " << NewObjLabel << ": " << NewObjLabel2 << endl << endl << endl;
  ierr = strcmp(NewObjLabel2,NewObjLabel);
  assert (ierr==0);

  Epetra_Object obj3(NewObjLabel,1); // pass a label and a TracebackMode
  char* NewObjLabel3 = obj3.Label();
  int TbM1 = obj3.GetTracebackMode();
  if (verbose) cout << endl << "This should say " << NewObjLabel << "," << "1: " << NewObjLabel3 << "," << TbM1 << endl << endl << endl;
  ierr = strcmp(NewObjLabel3,NewObjLabel);
  assert (ierr==0);
  assert (1==TbM1);
  
  Epetra_Object obj4(obj3); // copy constructor
  char* NewObjLabel4 = obj4.Label();
  int TbM2 = obj4.GetTracebackMode();
  if (verbose) cout << endl << "This should say " << NewObjLabel << "," << "1: " << NewObjLabel4 << "," << TbM2 << endl << endl << endl;
  ierr = strcmp(NewObjLabel4,NewObjLabel);
  assert (ierr==0);
  assert (1==TbM2);
  

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  return 0;
}

/*
  end of file main.cc
*/
