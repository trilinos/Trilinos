////////////////////////////////////////////////////////////////////
// Basic C++ example of using Zoltan to compute an RCB partitioning
// of a very simple mesh.
////////////////////////////////////////////////////////////////////

#include <rectangularMesh.h>

static void MPIExit()
{
#ifdef MPICPP
  MPI::Finalize();
#else
  MPI_Finalize();
#endif
}

int main(int argc, char *argv[])
{
  /////////////////////////////////
  // Initialize MPI and Zoltan
  /////////////////////////////////

  int rank, size;
  float version;

#ifdef MPICPP
  MPI::Init(argc, argv);
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
#else
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  Zoltan_Initialize(argc, argv, &version);

  /////////////////////////////////
  // Create a Zoltan object 
  /////////////////////////////////
  
#ifdef MPICPP
  Zoltan *zz = new Zoltan(MPI::COMM_WORLD);
#else
  Zoltan *zz = new Zoltan(MPI_COMM_WORLD);
#endif

  if (zz == NULL){
    MPIExit();
    exit(0);
  }

  ///////////////////////////////////////////////////////////
  // Create a simple rectangular mesh.  It's vertices are the
  // objects to be partitioned.
  ///////////////////////////////////////////////////////////

  rectangularMesh *mesh = new rectangularMesh();

  mesh->set_x_dim(10);
  mesh->set_y_dim(8);
  mesh->set_x_stride(1);
  mesh->set_y_stride(4);

  ///////////////////////////////////////////////////////////
  // Supply parameters to Zoltan to govern the partitioning
  // of the vertices across the processes.  Also supply
  // functions that Zoltan can call to obtain mesh information.
  // 
  // Many more parameters can be found in the User's Guide.
  ///////////////////////////////////////////////////////////

  // General parameters:

  zz->Set_Param("DEBUG_LEVEL", "0");     // amount of debug messages desired
  zz->Set_Param("LB_METHOD", "RCB");     // recursive coordinate bisection
  zz->Set_Param("NUM_GID_ENTRIES", "1"); // number of integers in a global ID
  zz->Set_Param("NUM_LID_ENTRIES", "1"); // number of integers in a local ID
  zz->Set_Param("OBJ_WEIGHT_DIM", "1");  // dimension of a vertex weight
  zz->Set_Param("RETURN_LISTS", "ALL");  // return all lists in LB_Partition

  // RCB parameters:
  //   Note that relaxing RCB_RECTILINEAR_BLOCKS (the requirement that
  //   all partitions be rectangles) may provide a better balance, but
  //   more communication.

  zz->Set_Param("KEEP_CUTS", "1");              // save decomposition
  zz->Set_Param("RCB_OUTPUT_LEVEL", "0");       // amount of output desired
  zz->Set_Param("RCB_RECTILINEAR_BLOCKS", "1"); // create rectilinear regions

  // Query functions:

  zz->Set_Num_Obj_Fn(rectangularMesh::get_number_of_objects, (void *)mesh);
  zz->Set_Obj_List_Fn(rectangularMesh::get_object_list, (void *)mesh);
  zz->Set_Num_Geom_Fn(rectangularMesh::get_num_geometry, NULL);
  zz->Set_Geom_Multi_Fn(rectangularMesh::get_geometry_list, (void *)mesh);
  
  ////////////////////////////////////////////////////////////////
  // Zoltan can now partition the vertices in the simple mesh.
  // In this simple example, we assume the number of partitions is
  // equal to the number of processes.  Process rank 0 will own
  // partition 0, process rank 1 will own partition 1, and so on.
  ////////////////////////////////////////////////////////////////

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
  
  if (rc != ZOLTAN_OK){
    printf("Partitioning failed on process %d\n",rank);
    MPIExit();
    delete zz;
    delete mesh;
    exit(0);
  }

  ////////////////////////////////////////////////////////////////
  // In a real application, you would rebalance the problem now by
  // sending the objects to their new partitions.
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // Visualize the new partitioning.
  // Create a list of GIDs now assigned to my partition, let
  // process zero display the partitioning.
  ////////////////////////////////////////////////////////////////

  int nMyGids = rectangularMesh::get_number_of_objects((void *)mesh, &rc);
  int nGids   = 0;
  int includeVertexWeights = 1; /* if OBJ_WEIGHT_DIM is 0, this should be 0 */

#ifdef MPICPP
  MPI::COMM_WORLD.Allreduce(&nMyGids, &nGids, 1, MPI::INT, MPI::SUM);
#else
  MPI_Allreduce(&nMyGids, &nGids, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  int *gid_list = new int[nMyGids];
  int *lid_list = new int[nMyGids];

  rectangularMesh::get_object_list((void *)mesh, nMyGids, nMyGids,
                  (ZOLTAN_ID_PTR)gid_list, (ZOLTAN_ID_PTR)lid_list,
                  0, NULL, &rc);

  mesh->draw_partitions("initial", nMyGids, gid_list, includeVertexWeights);

  int *gid_flags = new int[nGids];
  memset(gid_flags, 0, sizeof(int) * nGids);
  for (int i=0; i<nMyGids; i++){
    gid_flags[gid_list[i]-1] = 1;    // my original vertices 
  }
  for (int i=0; i<numImport; i++){
    gid_flags[importGlobalIds[i] - 1] = 1;  // my imports 
  }
  for (int i=0; i<numExport; i++){
    gid_flags[exportGlobalIds[i] - 1] = 0;  // my exports 
  }
  int nextIdx = 0;
  for (int i=0; i<nGids; i++){
    if (gid_flags[i]){
      gid_flags[nextIdx] = i+1; // my new GID list
      nextIdx++;
    }
  }
  mesh->draw_partitions("new partitioning", 
                        nextIdx, gid_flags, includeVertexWeights);

  delete [] gid_flags;
  delete [] gid_list;
  delete [] lid_list;

  ////////////////////////////////////////////////////////////////
  // Free the arrays allocated by LB_Partition, and free
  // the storage allocated for the Zoltan structure and the mesh.
  ////////////////////////////////////////////////////////////////

  zz->LB_Free_Part(&importGlobalIds, &importLocalIds, &importProcs,
                   &importToPart);
  zz->LB_Free_Part(&exportGlobalIds, &exportLocalIds, &exportProcs,
                   &exportToPart);

  delete zz;
  delete mesh;

  ////////////////////////////////////////////////////////////////
  // all done ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  MPIExit();

  return 0;
}
