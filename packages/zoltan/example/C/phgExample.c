/*
** $Id$
**
** Basic example of using Zoltan to partition a hypergraph.
**
** Also a test of Zoltan library's ability to handle very general
** hypergraph query functions: 
**   each process returns rows (use COMPLETE_ROWS in exSetHGDivisions)
**   each process returns columns                  (use COMPLETE_COLS)
**   each returns an arbitrary set of pins                (use BLOCKS)
**   some processes return no pins    (run on odd number of processes)
**
** (Note: "rows" are hyperedges, "columns" are vertices.)
** Zoltan_* functions are functions in the Zoltan library.
**
** ex* functions are functions in the example library, a library
**  created to make it possible to write simple examples.  Your
**  application would replace these with functions that access
**  your own data.
*/

#include "zoltan.h"
#include "exzoltan.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
  int rc, me, nprocs;
  float ver;
  struct Zoltan_Struct *zz;
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
  int *phgHandle;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /*
  ** The example library defines a simple hypergraph, and an initial
  ** partioning of vertices across the processes.
  **
  ** Also, based on certain options, the example library determines
  ** which pins will be supplied by each process in the application
  ** supplied query functions.  The query function interface is very
  ** general.  Each process can supply any subset of hypergraph
  ** pins (in a compressed row storage format or a compressed column
  ** storage format), but no two processes may supply the same
  ** pin.
  */

  phgHandle = 
    exSetHGDivisions(
      BLOCKS,          /* each process has some rows of hg */
      ALL_HAVE_EDGE_WEIGHTS,  /* each process has some of the edge weights */ 
      1,                  /* do edge weights supplied by processes overlap */ 
      ZOLTAN_COMPRESSED_EDGE); /* compressed pin format for query function */

  if (!phgHandle){
    printf("error in initializing example library\n");
    exit(0);
  }

  /*
  ** All applications which use the Zoltan library must call Zoltan_Initialize
  */
  rc = Zoltan_Initialize(argc, argv, &ver);

  if (rc != ZOLTAN_OK){
    printf("sorry...\n");
    exit(0);
  }

  /******************************************************************
  ** Prepare to partition the example hypergraph using the Zoltan 
  ** parallel hypergraph package.
  **
  ** Set parameter values, and supply to Zoltan the query functions 
  ** defined in the example library which will provide to the Zoltan 
  ** library the hypergraph data.
  ******************************************************************/

  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
  Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); /* test lib gives 1 vtx weight */
  Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");/* test lib gives 1 edge weight*/ 

  /* Graph parameters */

#if 0
  Zoltan_Set_Param(zz, "ADD_OBJ_WEIGHT", "unit"); /* Add a unit weight to each vtx */
  Zoltan_Set_Param(zz, "ADD_OBJ_WEIGHT", "pins"); /* Add a weight equal to # pins  */
#else
  Zoltan_Set_Param(zz, "ADD_OBJ_WEIGHT", "none"); /* Don't calculate extra weights */
#endif

  /* Parallel hypergraph parameters */

  Zoltan_Set_Param(zz, "PHG_OUTPUT_LEVEL", "0"); /* in general, be quiet   */
  Zoltan_Set_Param(zz, "PHG_VERTEX_VISIT_ORDER", "0");  /* see User's Guide */

#if 0
  /* 
   * If multiple processes supply a weight for the same edge, do
   * we "add" them or do we take the "max"imum weight?  Or do we
   * return an "error" if the weights do not agree.
   */
  Zoltan_Set_Param(zz, "PHG_EDGE_WEIGHT_OPERATION", "error");
  Zoltan_Set_Param(zz, "PHG_EDGE_WEIGHT_OPERATION", "add");
#else
  Zoltan_Set_Param(zz, "PHG_EDGE_WEIGHT_OPERATION", "max");
#endif

  /* Application defined query functions (defined in exphg.c) */

  Zoltan_Set_Num_Obj_Fn(zz, exGetHgNumVertices, phgHandle);
  Zoltan_Set_Obj_List_Fn(zz, exGetHgVerticesAndWeights, phgHandle);
  Zoltan_Set_HG_Size_CS_Fn(zz, exGetHgSizeAndFormat, phgHandle);
  Zoltan_Set_HG_CS_Fn(zz, exGetHg, phgHandle);
  Zoltan_Set_HG_Size_Edge_Weights_Fn(zz, exGetHgEdgeWeightSize, phgHandle);
  Zoltan_Set_HG_Edge_Weights_Fn(zz, exGetHgEdgeWeights, phgHandle);

  /* Parallel hypergraph partitioning occurs now */

  rc = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries,
    &numImport, &importGlobalIds, &importLocalIds, 
    &importProcs, &importToPart,
    &numExport, &exportGlobalIds, &exportLocalIds, 
    &exportProcs, &exportToPart);

  /* Processes compare rc for global success or failure */

  rc = exGlobalSuccess(rc, nprocs, me, (me == 0));

  if (!rc){
    if (me == 0){
      printf("Partitioning succeeded on all processes\n");
    }
  }
  else{
    goto End;
  }

  /* update list of objects assigned to this process under new partitioning */

  exUpdateDivisions(phgHandle, numExport, numImport, 
                   exportGlobalIds, importGlobalIds);

  /* evaluate the quality of the partitioning */

  Zoltan_LB_Eval(zz, 1, NULL, NULL, NULL, NULL, NULL, NULL);


  /* Print out the partitioning to a text file */

  Zoltan_Generate_Files(zz, "example", 0, 0, 0, 1);

End:
  /*
  ** Clean up
  */

  Zoltan_LB_Free_Part(&importGlobalIds, &importLocalIds, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalIds, &exportLocalIds, 
                      &exportProcs, &exportToPart);

  Zoltan_Destroy(&zz);

  exFreeDivisions(phgHandle);

  /**********************
  ** all done ***********
  **********************/

  MPI_Finalize();

  return 0;
}
