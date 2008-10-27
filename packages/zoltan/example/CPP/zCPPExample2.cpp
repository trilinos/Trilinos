// $Id$
//
//  C++ example of Zoltan library
//
//  MPICPP - Define this if your C++ implementation of MPI works.
//  NAMESPACES_OK - Define this if your system uses namespaces.
//

#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <zoltan_cpp.h>
#include <iostream>
#include "exzoltan.h"

//#define MPICPP
//#define NAMESPACES_OK

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
static Zoltan *makeZoltanObject(int method);
static void SetRCB_Parameters(Zoltan &);
static void SetRIB_Parameters(Zoltan &);
static void SetHSFC_Parameters(Zoltan &);

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
  Zoltan *zz = makeZoltanObject(RCB);
  delete zz;

  // Try RIB

  if (rank == 0) cout << "\nRecursive Inertial Bisection" << endl;
  zz = makeZoltanObject(RIB);
  delete zz;

  // Try HSFC

  if (rank == 0) cout << "\nHilbert Space Filling Curve" << endl;
  zz = makeZoltanObject(HSFC);
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

static Zoltan *makeZoltanObject(int method)
{
#ifdef MPICPP
  Zoltan *zz = new Zoltan(MPI::COMM_WORLD);
#else
  Zoltan *zz = new Zoltan(MPI_COMM_WORLD);
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
  zz->Set_Geom_Multi_Fn(exGetObjectCoords, NULL);

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

  Zoltan::LB_Free_Part(&importGlobalIds, &importLocalIds, &importProcs, 
                   &importToPart);
  Zoltan::LB_Free_Part(&exportGlobalIds, &exportLocalIds, &exportProcs, 
                   &exportToPart);

  return zz;
}
static void SetRCB_Parameters(Zoltan &zz)
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
static void SetHSFC_Parameters(Zoltan &zz)
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
static void SetRIB_Parameters(Zoltan &zz)
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
