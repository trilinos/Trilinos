// $Id$
//
//  C++ example of Zoltan library, including test of
//   copy constructor and copy operator.
//
//  MPICPP - Define this if your C++ implementation of MPI works.
//  NAMESPACES_OK - Define this if your system uses namespaces.
//  TEST_COPY - Define this if you want to test the Zoltan C++
//              copy operator and copy constructor.
//

#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <zoltan_cpp.h>
#include <iostream>
#include "exzoltan.h"

//#define MPICPP
//#define NAMESPACES_OK
#define TEST_COPY

#ifdef NAMESPACES_OK
 using namespace std;

 #ifdef MPICPP
 using namespace MPI;
 #endif
#else
 #define cout std::cout
 #define endl std::endl
#endif

static int MyNumPts=0;
static int *Gids;
static float *Pts;
static int rank, size;

#define RCB 1
#define RIB 2
#define HSFC 3

static void FreePoints();
static void MPIExit();
static Zoltan_Object *makeZoltanObject(int method);
static void SetRCB_Parameters(Zoltan_Object &);
static void SetRIB_Parameters(Zoltan_Object &);
static void SetHSFC_Parameters(Zoltan_Object &);

int main(int argc, char *argv[])
{
#ifdef MPICPP
  MPI::Init(argc, argv);
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();  
#else
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  if (getenv("WAIT_FOR_DEBUGGER")) sleep(60);

  // Create a rectilinear test mesh and divide it across processors.
  // exSetDivisions(div) sets mesh dimensions to div X div X div

  exSetDivisions(32);

  MyNumPts = exInitializePoints(&Pts, &Gids, rank, size);

  // Initialize Zoltan. 

  float version;

  Zoltan_Initialize(argc, argv, &version); 

  // Create Zoltan object.  

  if (rank == 0) cout << "\nRecursive Coordinate Bisection" << endl;
  Zoltan_Object *zz = makeZoltanObject(RCB);

#ifdef TEST_COPY
  // Create a copy - this will use the copy constructor

  Zoltan_Object zzCopy = *zz;

  // Check that the two look similar:
  if (0 == rank){
    cout << "RCB structure for original:" << endl;
    zz->PrintRCB(10);

    cout << "RCB structure for copy:" << endl;
    zzCopy.PrintRCB(10);
  }
  delete zz;

  // Verify copy is still valid after original was deleted

  if (0 == rank){
    cout << "RCB structure for copy, after deleting original:" << endl;
    zzCopy.PrintRCB(10);
  }
#else
  if (0 == rank){
    zz->PrintRCB(10);
  }
  delete zz;
#endif

  // Try RIB

  if (rank == 0) cout << "\nRecursive Inertial Bisection" << endl;
  zz = makeZoltanObject(RIB);

#ifdef TEST_COPY
  // Create a copy - this will use the copy constructor

  Zoltan_Object RIBCopy = *zz;

  // Create a copy using the copy operator

  Zoltan_Object zzCopy2;
  zzCopy2 = *zz;

  if (0 == rank){
    cout << "(" << rank << ") RIB structure for original:" << endl;
    zz->PrintRIB(10);

    cout << "(" << rank << ") RIB structure for copy from copy constructor:" << endl;
    RIBCopy.PrintRIB(10);

    cout << "(" << rank << ") RIB structure for copy from copy operator:" << endl;
    zzCopy2.PrintRIB(10);
  }
#else
  if (0 == rank){
    zz->PrintRIB(10);
  }
#endif

  delete zz;

  // Try HSFC

  if (rank == 0) cout << "\nHilbert Space Filling Curve" << endl;
  zz = makeZoltanObject(HSFC);

#ifdef TEST_COPY
  // Create a copy - this will use the copy constructor

  Zoltan_Object HSFCCopy = *zz;

  // Check that the two look similar:
  if (0 == rank){
    cout << "(" << rank << ") HSFC structure for original:" << endl;
    zz->PrintHSFC(10);

    cout << "(" << rank << ") HSFC structure for copy:" << endl;
    HSFCCopy.PrintHSFC(10);
  }
#else
  if (rank == 0){
    zz->PrintHSFC(10);
  }
#endif

  delete zz;

  FreePoints();

  MPIExit();
}
static void FreePoints()
{
  // These were allocated by a C-function, so free them
  // with "free()".

  if (MyNumPts > 0)
    {
    free(Gids);
    free(Pts);
    }
}
static void MPIExit()
{
#ifdef MPICPP
  MPI::Finalize();
#else
  MPI_Finalize();
#endif
}

static Zoltan_Object *makeZoltanObject(int method)
{
#ifdef MPICPP
  Zoltan_Object *zz = new Zoltan_Object(MPI::COMM_WORLD);
#else
  Zoltan_Object *zz = new Zoltan_Object(MPI_COMM_WORLD);
#endif

  if (method == RCB){
    SetRCB_Parameters(*zz);
  }
  else if (method == RIB){
    SetRIB_Parameters(*zz);
  }
  else if (method == HSFC){
    SetHSFC_Parameters(*zz);
  }
  else{
    cout << "ERROR" << endl;
    exit(0);
  }

  // Call backs:

  zz->Set_Num_Obj_Fn(exGetNumberOfAssignedObjects, NULL);
  zz->Set_Obj_List_Fn(exGetObjectList, NULL);
  zz->Set_Num_Geom_Fn(exGetObjectSize, NULL);
  zz->Set_Geom_Multi_Fn(exGetObject, NULL);

  // Perform the load balancing partitioning

  int changes;
  int numGidEntries;
  int numLidEntries;
  int numImport;
  ZOLTAN_ID_PTR importGlobalIds;
  ZOLTAN_ID_PTR importLocalIds;
  int *importProcs;
  int *importToPart;
  int numExport;
  ZOLTAN_ID_PTR exportGlobalIds;
  ZOLTAN_ID_PTR exportLocalIds;
  int *exportProcs;
  int *exportToPart;

  int rc = zz->LB_Partition(changes, numGidEntries, numLidEntries, 
    numImport, importGlobalIds, importLocalIds, importProcs, importToPart,
    numExport, exportGlobalIds, exportLocalIds, exportProcs, exportToPart);

  rc = exGlobalSuccess(rc, size, rank, (rank == 0));
  
  if (!rc){
    exPrintGlobalResult("==============Result==============", size, rank,
                           MyNumPts, numImport, numExport, changes);
  }
  else{
    MPIExit();
    FreePoints();
    exit(0);
  }
  
  // Free the memory allocated for lists returned by LB_Parition()

  zz->LB_Free_Part(importGlobalIds, importLocalIds, importProcs, importToPart);
  zz->LB_Free_Part(exportGlobalIds, exportLocalIds, exportProcs, exportToPart);

  return zz;
}
static void SetRCB_Parameters(Zoltan_Object &zz)
{
  /* General parameters */

  zz.Set_Param("DEBUG_LEVEL", "0");
  zz.Set_Param("LB_METHOD", "RCB");
  zz.Set_Param("NUM_GID_ENTRIES", "1");
  zz.Set_Param("NUM_LID_ENTRIES", "1");
  zz.Set_Param("RETURN_LISTS", "ALL");

  /* RCB parameters */

  zz.Set_Param("KEEP_CUTS", "1");
  zz.Set_Param("RCB_OUTPUT_LEVEL", "0");
  zz.Set_Param("RCB_RECTILINEAR_BLOCKS", "1");
}
static void SetHSFC_Parameters(Zoltan_Object &zz)
{
  /* General parameters */

  zz.Set_Param("DEBUG_LEVEL", "0");
  zz.Set_Param("LB_METHOD", "HSFC");
  zz.Set_Param("NUM_GID_ENTRIES", "1");
  zz.Set_Param("NUM_LID_ENTRIES", "1");
  zz.Set_Param("RETURN_LISTS", "ALL");

  /* HSFC parameters */

  zz.Set_Param("KEEP_CUTS", "1");
}
static void SetRIB_Parameters(Zoltan_Object &zz)
{
  /* General parameters */

  zz.Set_Param("DEBUG_LEVEL", "0");
  zz.Set_Param("LB_METHOD", "RIB");
  zz.Set_Param("NUM_GID_ENTRIES", "1");
  zz.Set_Param("NUM_LID_ENTRIES", "1");
  zz.Set_Param("RETURN_LISTS", "ALL");

  /* RIB parameters */

  zz.Set_Param("KEEP_CUTS", "1");
  zz.Set_Param("RIB_OUTPUT_LEVEL", "0");
  zz.Set_Param("AVERAGE_CUTS", "1");
}
