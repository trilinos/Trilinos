// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// Basic C++ example of using Zoltan to compute a quick partitioning
// of a set of objects.
//

#ifdef MPICPP
#undef MPICPP
#endif /* MPICPP */

//#define MPICPP // Uncomment to use C++ interface for MPI.

#include <mpi.h>
#include <stdio.h>
#include <zoltan_cpp.h>

// Class representing collection of objects to be partitioned.

class objectCollection {

private:

  int numGlobalObjects;
  int numMyObjects;
  int *myGlobalIDs;

public:

  // constructor

  objectCollection():numGlobalObjects(0), numMyObjects(0), myGlobalIDs(NULL){}

  // destructor

  ~objectCollection(){ if (myGlobalIDs) delete [] myGlobalIDs;}

  void set_num_global_objects(int n) {numGlobalObjects = n;}
  int get_num_global_objects() {return numGlobalObjects;}

  void set_num_my_objects(int n) {numMyObjects = n;}
  int get_num_my_objects() {return numMyObjects;}

  void set_my_global_ids(int *p) {myGlobalIDs = p;}
  int *get_my_global_ids() {return myGlobalIDs;}

  // query functions that respond to requests from Zoltan 

  static int get_number_of_objects(void *data, int *ierr){

    objectCollection *objs = (objectCollection *)data;
    *ierr = ZOLTAN_OK;

     return objs->numMyObjects;
  }

  static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr){

    objectCollection *objs = (objectCollection *)data;
    *ierr = ZOLTAN_OK;

    // In this example, return the IDs of our objects, but no weights.
    // Zoltan will assume equally weighted objects.

    for (int i=0; i<objs->get_num_my_objects(); i++){
      globalID[i] = objs->get_my_global_ids()[i];
      localID[i] = i;
    }
    return;
  }

};

static const char *global_fname="objects.txt"; // File containing objects to be partitioned.

static int get_next_line(FILE *fp, char *buf, int bufsize);
static void input_file_error(int numProcs, int tag, int startProc);
static void showSimpleMeshPartitions(int myProc, int numIDs, int *GIDs, int *parts);
static void read_input_objects(int myRank, int numProcs, const char *fname, objectCollection &myData);

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

  //////////////////////////////////////////////////////////////////
  // Read objects from input file and distribute them unevenly
  //////////////////////////////////////////////////////////////////

  FILE *fp = fopen(global_fname, "r");
  if (!fp){
    if (rank == 0) fprintf(stderr,"ERROR: Can not open %s\n",global_fname);
    MPIExit();
    exit(1);
  }
  fclose(fp);

  objectCollection objects;

  read_input_objects(rank, size, global_fname, objects);

  ///////////////////////////////////////////////////////////////////
  // Set the Zoltan parameters, and the names of the query functions
  ///////////////////////////////////////////////////////////////////

  // General parameters 

  zz->Set_Param( "LB_METHOD", "BLOCK");    /* Zoltan method: "BLOCK" */
  zz->Set_Param( "NUM_GID_ENTRIES", "1");  /* global ID is 1 integer */
  zz->Set_Param( "NUM_LID_ENTRIES", "1");  /* local ID is 1 integer */
  zz->Set_Param( "OBJ_WEIGHT_DIM", "0");   /* we omit object weights */

  // Query functions 

  zz->Set_Num_Obj_Fn(objectCollection::get_number_of_objects, &objects);
  zz->Set_Obj_List_Fn(objectCollection::get_object_list, &objects);

  ////////////////////////////////////////////////////////////////
  // Zoltan can now partition the objects in this collection.
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
    exit(0);
  }

  /////////////////////////////////////////////////////////////////
  // Visualize the partitioning before and after calling Zoltan.
  //
  // In this example, partition number equals process rank.
  /////////////////////////////////////////////////////////////////

  int *parts = new int [objects.get_num_my_objects()];

  for (int i=0; i < objects.get_num_my_objects(); i++){
    parts[i] = rank;
  }

  if (rank == 0){
    printf("\nObject partition assignments before calling Zoltan\n");
  }

  showSimpleMeshPartitions(rank, objects.get_num_my_objects(), 
                           objects.get_my_global_ids(), parts);

  for (int i=0; i < numExport; i++){
    parts[exportLocalIds[i]] = exportToPart[i];
  }

  if (rank == 0){
    printf("Object partition assignments after calling Zoltan\n");
  }

  showSimpleMeshPartitions(rank, objects.get_num_my_objects(), 
                           objects.get_my_global_ids(), parts);

  delete [] parts;

  ////////////////////////////////////////////////////////////////
  // Free the arrays allocated by LB_Partition, and free
  // the storage allocated for the Zoltan structure and the mesh.
  ////////////////////////////////////////////////////////////////

  Zoltan::LB_Free_Part(&importGlobalIds, &importLocalIds, &importProcs,
                   &importToPart);
  Zoltan::LB_Free_Part(&exportGlobalIds, &exportLocalIds, &exportProcs,
                   &exportToPart);

  delete zz;

  ////////////////////////////////////////////////////////////////
  // all done ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  MPIExit();

  return 0;
}

/* Function to find next line of information in input file */
 
static int get_next_line(FILE *fp, char *buf, int bufsize)
{
int i, cval, len;
char *c;

  while (1){

    c = fgets(buf, bufsize, fp);

    if (c == NULL)
      return 0;  /* end of file */

    len = strlen(c);

    for (i=0, c=buf; i < len; i++, c++){
      cval = (int)*c; 
      if (isspace(cval) == 0) break;
    }
    if (i == len) continue;   /* blank line */
    if (*c == '#') continue;  /* comment */

    if (c != buf){
      strcpy(buf, c);
    }
    break;
  }

  return strlen(buf);  /* number of characters */
}

// Proc 0 notifies others of error and exits 

static void input_file_error(int numProcs, int tag, int startProc)
{
int i, val;

  val = -1;

  fprintf(stderr,"ERROR in input file.\n");

  for (i=startProc; i < numProcs; i++){
    // these procs have posted receive for "tag" 
#ifdef MPICPP
    MPI::COMM_WORLD.Send(&val, 1, MPI::INT, i, tag);
#else
    MPI_Send(&val, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
#endif
    
  }
  for (i=1; i < startProc; i++){
    // these procs are done 
#ifdef MPICPP
    MPI::COMM_WORLD.Send(&val, 1, MPI::INT, i, 0);
#else
    MPI_Send(&val, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
#endif
  }

  MPIExit();
  exit(0);
}

/* Draw the partition assignments of the objects */

void showSimpleMeshPartitions(int myProc, int numIDs, int *GIDs, int *parts)
{
int partAssign[25], allPartAssign[25];
int i, j, part;

  memset(partAssign, 0, sizeof(int) * 25);

  for (i=0; i < numIDs; i++){
    partAssign[GIDs[i]-1] = parts[i];
  }

#ifdef MPICPP
  MPI::COMM_WORLD.Reduce(partAssign, allPartAssign, 25, MPI::INT, MPI::MAX, 0);
#else
  MPI_Reduce(partAssign, allPartAssign, 25, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
#endif

  if (myProc == 0){

    for (i=20; i >= 0; i-=5){
      for (j=0; j < 5; j++){
        part = allPartAssign[i + j];
        if (j < 4)
          printf("%d-----",part);
        else
          printf("%d\n",part);
      }
      if (i > 0)
        printf("|     |     |     |     |\n");
    }
    printf("\n");
  }
}

//
// Proc 0 reads the objects in the input file and divides them across processes 
//

void read_input_objects(int myRank, int numProcs, const char *fname, objectCollection &myData)
{
int val, nobj, remainingObj;
int obj_ack_tag = 5, obj_count_tag = 10, obj_id_tag = 15;
#ifdef MPICPP
MPI::Status status;
#else
MPI_Status status;
#endif

  if (myRank == 0){

    char *buf = new char [512];
    FILE *fp = fopen(fname, "r");

    int num = get_next_line(fp, buf, 512);
    if (num == 0) input_file_error(numProcs, obj_count_tag, 1);
    num = sscanf(buf, "%d", &val);
    myData.set_num_global_objects(val);
    if (num != 1) input_file_error(numProcs, obj_count_tag, 1);

    if (numProcs > 1){
      nobj = myData.get_num_global_objects() / 2;
      remainingObj = myData.get_num_global_objects() - nobj;
    }
    else{
      nobj = myData.get_num_global_objects();
      remainingObj = 0;
    }

    int *mygids = new int [nobj];
    myData.set_num_my_objects(nobj);
    myData.set_my_global_ids(mygids);

    for (int i=0; i < nobj; i++){

      num = get_next_line(fp, buf, 512);
      if (num == 0) input_file_error(numProcs, obj_count_tag, 1);
      num = sscanf(buf, "%d", &val);
      if (num != 1) input_file_error(numProcs, obj_count_tag, 1);
      mygids[i] = val;
  
    }

    int *gids = new int [nobj + 1];
    int ack = 0;

    for (int i=1; i < numProcs; i++){
    
      if (remainingObj > 1){
        nobj = remainingObj / 2;
        remainingObj -= nobj;
      }
      else if (remainingObj == 1){
        nobj = 1;
        remainingObj = 0;
      }
      else{
        nobj = 0;
      }

      if ((i == numProcs - 1) && (remainingObj > 0))
        nobj += remainingObj;

      if (nobj > 0){
        for (int j=0; j < nobj; j++){
          num = get_next_line(fp, buf, 512);
          if (num == 0) input_file_error(numProcs, obj_count_tag, i);
          num = sscanf(buf, "%d", &val);
          if (num != 1) input_file_error(numProcs, obj_count_tag, i);
          gids[j] = val;;
        }
      }

#ifdef MPICPP
      MPI::COMM_WORLD.Send(&nobj, 1, MPI::INT, i, obj_count_tag);
      MPI::COMM_WORLD.Recv(&ack, 1, MPI::INT, i, obj_ack_tag, status);
#else
      MPI_Send(&nobj, 1, MPI_INT, i, obj_count_tag, MPI_COMM_WORLD);
      MPI_Recv(&ack, 1, MPI_INT, i, obj_ack_tag, MPI_COMM_WORLD, &status);
#endif

      if (nobj > 0)
        MPI_Send(gids, nobj, MPI_INT, i, obj_id_tag, MPI_COMM_WORLD);
      
    }

    delete [] gids;
    delete [] buf;
    
    fclose(fp);

    /* signal all procs it is OK to go on */
    ack = 0;
    for (int i=1; i < numProcs; i++){
#ifdef MPICPP
      MPI::COMM_WORLD.Send(&ack, 1, MPI::INT, i, 0);
#else
      MPI_Send(&ack, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
#endif
    }
  }
  else{

#ifdef MPICPP
    MPI::COMM_WORLD.Recv(&val, 1, MPI::INT, 0, obj_count_tag, status);
#else
    MPI_Recv(&val, 1, MPI_INT, 0, obj_count_tag, MPI_COMM_WORLD, &status);
#endif

    myData.set_num_my_objects(val);

    int ack = 0;

    if (myData.get_num_my_objects() > 0){
      int *mygids = new int [myData.get_num_my_objects()];
#ifdef MPICPP
      MPI::COMM_WORLD.Send(&ack, 1, MPI::INT, 0, obj_ack_tag);
      MPI::COMM_WORLD.Recv(mygids, myData.get_num_my_objects(), MPI::INT, 0, 
               obj_id_tag, status);
#else
      MPI_Send(&ack, 1, MPI_INT, 0, obj_ack_tag, MPI_COMM_WORLD);
      MPI_Recv(mygids, myData.get_num_my_objects(), MPI_INT, 0, 
               obj_id_tag, MPI_COMM_WORLD, &status);
#endif

      myData.set_my_global_ids(mygids);
    }
    else if (myData.get_num_my_objects() == 0){
#ifdef MPICPP
      MPI::COMM_WORLD.Send(&ack, 1, MPI::INT, 0, obj_ack_tag);
#else
      MPI_Send(&ack, 1, MPI_INT, 0, obj_ack_tag, MPI_COMM_WORLD);
#endif
    }
    else{
      MPIExit();
      exit(1);
    }

#ifdef MPICPP
    MPI::COMM_WORLD.Recv(&ack, 1, MPI::INT, 0, 0, status);
#else
    MPI_Recv(&ack, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
#endif
    if (ack < 0){
      MPIExit();
      exit(1);
    }
  }
}
