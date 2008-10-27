// $Id$
//
//  Simple C++ example of Zoltan library
//
//  MPICPP - Define this if your C++ bindings for MPI work.
//  NAMESPACES_OK - Define this if your system uses namespaces.
//

#include <mpi.h>
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
#endif

static int MyNumPts=0;
static int *Gids;
static float *Pts;
static int rank, size;

static void FreePoints();
static void MPIExit();

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

  // Create a rectilinear test mesh and divide it across processors.
  // exSetDivisions(div) sets mesh dimensions to div X div X div

  exSetDivisions(32);    /* rectilinear mesh is div X div X div */

  MyNumPts = exInitializePoints(&Pts, &Gids, rank, size);

  // Initialize Zoltan.  This is a C call.  The simple C++ code 
  // that creates Zoltan objects does not keep track of whether 
  // Zoltan_Initialize has been called.

  float version;

  Zoltan_Initialize(argc, argv, &version); 

  // Create Zoltan object.  This calls Zoltan_Create.  

#ifdef MPICPP
  Zoltan *zz = new Zoltan(MPI::COMM_WORLD);
#else
  Zoltan *zz = new Zoltan(MPI_COMM_WORLD);
#endif

  // General parameters:

  zz->Set_Param("DEBUG_LEVEL", "0");     // no debug messages
  zz->Set_Param("LB_METHOD", "RCB");     // perform recursive coordinate bisection
  zz->Set_Param("NUM_GID_ENTRIES", "1"); // global IDs are 1 integer
  zz->Set_Param("NUM_LID_ENTRIES", "1"); // local IDs are 1 integer
  zz->Set_Param("RETURN_LISTS", "ALL");  // return all lists in LB_Partition

  // RCB parameters:

  zz->Set_Param("KEEP_CUTS", "1");              // save decomposition
  zz->Set_Param("RCB_OUTPUT_LEVEL", "0");       // provide extra info
  zz->Set_Param("RCB_RECTILINEAR_BLOCKS", "1"); // don't split point on boundary

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
    exPrintGlobalResult("Recursive Coordinate Bisection", size, rank,
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

  FreePoints();

  // Implementation note:  A Zoltan object contains an MPI communicator.
  //   When the Zoltan object is destroyed, it uses it's MPI communicator.
  //   So it is important that the Zoltan object is destroyed before
  //   the MPI communicator is destroyed.  To ensure this, dynamically
  //   allocate the Zoltan object, so you can explicitly destroy it.
  //   If you create a Zoltan object on the stack, it's destructor will
  //   be invoked atexit, possibly after the communicator's
  //   destructor.

  delete zz;

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
