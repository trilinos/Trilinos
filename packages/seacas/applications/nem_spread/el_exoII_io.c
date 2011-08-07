/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sort_utils.h"

#include "rf_salsa.h"

#include "pe_common.h"
#include "pe_str_util_const.h"

#include "el_elm.h"
#include "el_geom_const.h"
#include "el_exoII_io.h"

#include "rf_fem_const.h"
#include "rf_allo.h"
#include "rf_message.h"
#include "rf_io_const.h"
#include "rf_mp_const.h"
#include "rf_util.h"

#include "netcdf.h"
#include "exodusII.h"
#include "ne_nemesisI.h"

#include "ps_pario_const.h"

/************* R O U T I N E S   I N   T H I S   F I L E **********************

  Name_of_Routine               type                 Called by
  ---------------------------------------------------------------

  load_mesh ()                                  main:rf_salsa.c
  read_coord ()                                 load_mesh
  read_elem_blk_ids ()                          load_mesh
  read_elem_blk ()                              load_mesh
  extract_elem_connect ()                         read_elem_blk
  extract_elem_attr ()                            read_elem_blk
  find_elem_block ()                              read_elem_blk
  read_node_set_ids ()                          load_mesh
  read_side_set_ids ()                          load_mesh
  read_node_sets ()                             load_mesh
  read_side_sets ()                             load_mesh
  read_nodal_vars ()                            load_mesh
  read_aux_nodal_vars()                         load_mesh
  ****************************************************************************/

/******** P R O T O T Y P E S   O F    F U N C T I O N S *********************/

extern int  check_monot       (int vector[], int length);
extern void check_exodus_error(
                               int  error,           /* error code         */
                               char function_name[]
                                    /* EXODUS function returning the error   */
                               );

static void read_coord(int mesh_exoid, int io_ws, int max_name_length);
static void read_elem_blk_ids(int mesh_exoid, int max_name_length);
static void read_elem_blk(int mesh_exoid, int io_ws);
static void extract_elem_blk(void);
static void extract_global_element_ids(int global_ids[], int Num_Elem, int iproc);
static void extract_global_node_ids(int global_ids[], int Num_Node, int iproc);
static int extract_elem_connect(int elem_blk[], int icurrent_elem_blk,
                                int istart_elem, int iend_elem,
                                int *local_ielem_blk, int indx);
static void extract_elem_attr(void *elem_attr, int icurrent_elem_blk,
                              int istart_elem, int iend_elem,
                              int natt_p_elem, int indx, int io_ws);
static void find_elem_block(int *proc_elem_blk, int indx, int proc_for);
static void read_node_set_ids(int mesh_exoid, int [], int [], int max_name_length);
static void read_side_set_ids(int mesh_exoid, int [], int [], int max_name_length);
static void read_node_sets(int mesh_exoid, int *, int *, int);
static void read_side_sets(int mesh_exoid, int *, int *, int);

void read_nodal_vars(int mesh_exoid);
void read_aux_nodal_vars(int mesh_exoid);
void create_elem_types(void);
void read_exo_nv(int exoid, int itotal_nodes, int n_glob_nodes, int proc_id,
                 int num_proc, int n_gnv_to_read, int *gnode, int time_index,
                 int *gnv_index, float **glob_var_arr);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void load_mesh(int io_ws)

/* Function which reads the EXODUS II database which contains the mesh.  It
 * does this under the control of another EXODUS II database file which
 * contains the load-balance information.  The load balance information must
 * have been already read.
 *
 *----------------------------------------------------------------------------
 *
 * Functions called:
 *
 * check_exodus_error -- function which handles the error code returned by
 *                        calls to EXODUS II API routines.
 * construct_lb_filename --  function which appends the string '-#' where
 *                            '#' is the number of processors to the name of
 *                            the EXODUS II filename.
 * read_lb          -- function which reads the load-balance information
 *                     from an EXODUS II database for a given processor.
 * read_coord       -- function which reads the nodal coordinates information
 *                     from an EXODUS II database for a given processor.
 * read_elem_blk_ids-- Function which read the element block ids and other
 *                     global information specific to the element block level.
 * read_elem_blk    -- Function which reads the element level information
 *                     for each element block in the mesh
 * read_node_sets   -- function which reads the node sets information from an
 *                     EXODUS II database for a given processor.
 * read_side_sets   -- function which reads the side sets information from an
 *                     EXODUS II database for a given processor.
 *
 *----------------------------------------------------------------------------
 */

{

  /* Local variables */

  int   *num_nodes_in_node_set = NULL, mesh_exoid;
  int   *num_elem_in_ssets=NULL, *num_df_in_ssets=NULL, *num_df_in_nsets=NULL;
  int    cpu_ws;
  int    i1, i2, iproc, index, length_qa, glob_pindx;
  char  *yo = "load_mesh: ";
  float  version;
  double start_time;
  int    max_name_length;

  char   cTemp[512];

  /*************************** execution begins ******************************/

  /* computing precision should be the same as the database precision
   *
   * EXCEPTION: if the io_ws is smaller than the machine precision,
   * ie - database with io_ws == 4 on a Cray (sizeof(float) == 8),
   * then the cpu_ws must be the machine precision.
   */
  if (io_ws < sizeof(float)) cpu_ws = sizeof(float);
  else                       cpu_ws = io_ws;

  /* Allocate some memory for each processor read by this processor */
  Proc_Num_Elem_Blk   = (int *)array_alloc(__FILE__, __LINE__, 1,
                                           5 * Proc_Info[2], sizeof(int));
  Proc_Num_Node_Sets       = Proc_Num_Elem_Blk   + Proc_Info[2];
  Proc_NS_List_Length      = Proc_Num_Node_Sets  + Proc_Info[2];
  Proc_Num_Side_Sets       = Proc_NS_List_Length + Proc_Info[2];
  Proc_SS_Elem_List_Length = Proc_Num_Side_Sets  + Proc_Info[2];

  /* Initialize */
  for(iproc=0; iproc < 5*Proc_Info[2]; iproc++)
    Proc_Num_Elem_Blk[iproc] = 0;

  GElem_Blks           = (int **)array_alloc(__FILE__, __LINE__, 1,
                                             26 * Proc_Info[2], sizeof(int *));
  Proc_Nodes_Per_Elem   = GElem_Blks            + Proc_Info[2];
  Proc_Elem_Blk_Ids     = Proc_Nodes_Per_Elem   + Proc_Info[2];
  Proc_Elem_Blk_Types   = Proc_Elem_Blk_Ids     + Proc_Info[2];
  Proc_Num_Attr         = Proc_Elem_Blk_Types   + Proc_Info[2];
  Proc_Num_Elem_In_Blk  = Proc_Num_Attr         + Proc_Info[2];
  Proc_Connect_Ptr      = Proc_Num_Elem_In_Blk  + Proc_Info[2];
  Proc_Elem_Connect     = Proc_Connect_Ptr      + Proc_Info[2];
  Proc_NS_Ids           = Proc_Elem_Connect     + Proc_Info[2];
  Proc_NS_Count         = Proc_NS_Ids           + Proc_Info[2];
  Proc_NS_DF_Count      = Proc_NS_Count         + Proc_Info[2];
  Proc_NS_Pointers      = Proc_NS_DF_Count      + Proc_Info[2];
  Proc_NS_List          = Proc_NS_Pointers      + Proc_Info[2];
  GNode_Sets            = Proc_NS_List          + Proc_Info[2];
  Proc_NS_GNMap_List    = GNode_Sets            + Proc_Info[2];
  Proc_SS_Ids           = Proc_NS_GNMap_List    + Proc_Info[2];
  Proc_SS_Elem_Count    = Proc_SS_Ids           + Proc_Info[2];
  Proc_SS_DF_Count      = Proc_SS_Elem_Count    + Proc_Info[2];
  Proc_SS_Elem_Pointers = Proc_SS_DF_Count      + Proc_Info[2];
  Proc_SS_Elem_List     = Proc_SS_Elem_Pointers + Proc_Info[2];
  Proc_SS_Side_List     = Proc_SS_Elem_List     + Proc_Info[2];
  Proc_SS_DF_Pointers   = Proc_SS_Side_List     + Proc_Info[2];
  GSide_Sets            = Proc_SS_DF_Pointers   + Proc_Info[2];
  Proc_SS_GEMap_List    = GSide_Sets            + Proc_Info[2];
  Proc_Global_Elem_Id_Map= Proc_SS_GEMap_List   + Proc_Info[2];
  Proc_Global_Node_Id_Map=Proc_Global_Elem_Id_Map+ Proc_Info[2];

  /* Initialize */
  for(iproc=0; iproc < 25*Proc_Info[2]; iproc++)
    GElem_Blks[iproc] = NULL;

  if (io_ws <= sizeof(float)) {
    Coor_sp = (float ***)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2],
                                     sizeof(float **));

    /* Initialize */
    for(iproc=0; iproc < Proc_Info[2]; iproc++)
      Coor_sp[iproc] = NULL;

    Proc_Elem_Attr_sp  = (float **) array_alloc(__FILE__, __LINE__, 1,
                                                3 * Proc_Info[2],
                                                sizeof(float *));
    Proc_NS_Dist_Fact_sp = Proc_Elem_Attr_sp    + Proc_Info[2];
    Proc_SS_Dist_Fact_sp = Proc_NS_Dist_Fact_sp + Proc_Info[2];

    /* Initialize */
    for(iproc=0; iproc < 3*Proc_Info[2]; iproc++)
      Proc_Elem_Attr_sp[iproc] = NULL;
  }
  else {
    Coor_dp = (double ***)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2],
                                      sizeof(double **));

    /* Initialize */
    for(iproc=0; iproc < Proc_Info[2]; iproc++)
      Coor_dp[iproc] = NULL;

    Proc_Elem_Attr_dp  = (double **) array_alloc(__FILE__, __LINE__, 1,
                                                 3 * Proc_Info[2],
                                                 sizeof(double *));
    Proc_NS_Dist_Fact_dp = Proc_Elem_Attr_dp    + Proc_Info[2];
    Proc_SS_Dist_Fact_dp = Proc_NS_Dist_Fact_dp + Proc_Info[2];

    /* Initialize */
    for(iproc=0; iproc < 3*Proc_Info[2]; iproc++)
      Proc_Elem_Attr_dp[iproc] = NULL;
  }

  if (Proc == 0) {

    /* Check for a problem which has too many processors for a given mesh */

    if (Num_Node/Proc_Info[0] < 1) {
      fprintf(stderr, "%sERROR: Problem divided among too many "
              "processors.\n", yo);
      exit(1);
    }
    else if (Num_Elem/Proc_Info[0] < 1) {
      fprintf(stderr, "%sERROR: Problem divided among too many "
              "processors.\n", yo);
      exit(1);
    }

    /* Open the EXODUS II mesh file */

    mesh_exoid = ex_open(ExoFile, EX_READ, &cpu_ws, &io_ws, &version);
    if (mesh_exoid < 0) {
      fprintf(stderr, "%sExodus returned error opening mesh file, %s\n",
              yo, ExoFile);
      exit(1);
    }

    max_name_length = ex_inquire_int(mesh_exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
    ex_set_max_name_length(mesh_exoid, max_name_length);

    if (ex_inquire(mesh_exoid, EX_INQ_QA, &Num_QA_Recs, (float *) NULL,
                   (char *) NULL) < 0) {
      fprintf(stderr,
              "%sERROR, could not determine number of QA records\n", yo);
      exit(1);
    }
  }

  /* Broadcast a count of QA Records to every processor */
  brdcst(Proc, Num_Proc, (char *)&Num_QA_Recs, sizeof(int), 0);

  /* Add a QA record for the spreader */
  length_qa = 4*(Num_QA_Recs+1);
  QA_Record = (char **) array_alloc(__FILE__, __LINE__, 1, length_qa,
                                    sizeof(char *));
  for (i1 = 0, index = 0; i1 < Num_QA_Recs+1; i1++) {
    for (i2 = 0; i2 < 4; i2++) {
      QA_Record[index++] = (char *)array_alloc(__FILE__, __LINE__, 1,
                                               MAX_STR_LENGTH+1, sizeof(char));
    }
  }

  if(Proc == 0) {
    if (ex_get_qa(mesh_exoid, (char *(*)[]) &QA_Record[0]) < 0) {
      fprintf(stderr, "%sERROR, could not get QA record(s)\n", yo);
      exit(1);
    }

  } /* End "if (Proc == 0)" */

  /* Now broadcast each of the QA Records to all processors */
  for (i1=0, index=0; i1 < Num_QA_Recs; i1++) {
    for (i2=0; i2 < 4; i2++) {
      brdcst(Proc, Num_Proc, QA_Record[index++],
             (MAX_STR_LENGTH+1)*sizeof(char), 0);
    }
  }

  /* Read in the information records */
  if (Proc == 0) {

    if (ex_inquire(mesh_exoid, EX_INQ_INFO, &Num_Info_Recs, (float *) NULL,
                   (char *) NULL) < 0) {
      fprintf(stderr,
              "%sERROR, could not determine number of Info records\n", yo);
      exit(1);
    }
  }

  /* Broadcast a count of Info Records to every processor */
  brdcst(Proc, Num_Proc, (char *)&Num_Info_Recs, sizeof(int), 0);

  if(Num_Info_Recs > 0) {
    /* Allocate the Information records */
    Info_Record = (char **)array_alloc(__FILE__, __LINE__, 2, Num_Info_Recs,
                                       MAX_LINE_LENGTH + 1, sizeof(char));
    if(!Info_Record) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n", Proc, yo);
      exit(1);
    }

    if(Proc == 0) {
      if(ex_get_info(mesh_exoid, Info_Record) < 0) {
        fprintf(stderr, "%sERROR, could not get Info record(s)\n", yo);
        exit(1);
      }
    }

    /* Now broadcast the info records */
    for(i1=0; i1 < Num_Info_Recs; i1++) {
      brdcst(Proc, Num_Proc, Info_Record[i1],
             (MAX_LINE_LENGTH+1)*sizeof(char), 0);
    }
  }

  /* Read in the coordinate frame information */
  if (Proc == 0) {
    Num_Coordinate_Frames = ex_inquire_int(mesh_exoid, EX_INQ_COORD_FRAMES);
  }

  /* Broadcast a count of Coordinate Frames to every processor */
  brdcst(Proc, Num_Proc, (char *)&Num_Coordinate_Frames, sizeof(int), 0);

  if(Num_Coordinate_Frames > 0) {
    void *Coordinate_Frame_Coordinates = NULL;
    
    /* Allocate the Coordinate Frame records */
    Coordinate_Frame_Ids = (int *)array_alloc(__FILE__, __LINE__, 1, Num_Coordinate_Frames, sizeof(int));
    if(!Coordinate_Frame_Ids) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n", Proc, yo);
      exit(1);
    }

    if (io_ws <= sizeof(float)) {
      Coordinate_Frame_Coordinates_sp = (float  *)array_alloc(__FILE__, __LINE__, 1, 9*Num_Coordinate_Frames,
							      sizeof(float));
      if(!Coordinate_Frame_Coordinates_sp) {
	fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n", Proc, yo);
	exit(1);
      }
      Coordinate_Frame_Coordinates = Coordinate_Frame_Coordinates_sp;
    } else {
      Coordinate_Frame_Coordinates_dp = (double *)array_alloc(__FILE__, __LINE__, 1, 9*Num_Coordinate_Frames,
							      sizeof(double));
      if(!Coordinate_Frame_Coordinates_dp) {
	fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n", Proc, yo);
	exit(1);
      }
      Coordinate_Frame_Coordinates = Coordinate_Frame_Coordinates_dp;
    }
    Coordinate_Frame_Tags = (char *)array_alloc(__FILE__, __LINE__, 1, Num_Coordinate_Frames, sizeof(char));
    if(!Coordinate_Frame_Tags) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n", Proc, yo);
      exit(1);
    }

    if(Proc == 0) {
      int num_frames = 0;
      if(ex_get_coordinate_frames(mesh_exoid, &num_frames, Coordinate_Frame_Ids,
				  Coordinate_Frame_Coordinates, Coordinate_Frame_Tags) < 0) {
        fprintf(stderr, "%sERROR, could not get Coordinate Frame record(s)\n", yo);
        exit(1);
      }
      if (num_frames != Num_Coordinate_Frames) {
        fprintf(stderr, "%sERROR, frame count inconsistency\n", yo);
        exit(1);
      }
    }

    /* Now broadcast the coordinate frame records */
    brdcst(Proc, Num_Proc, (char *)Coordinate_Frame_Ids,           Num_Coordinate_Frames*sizeof(int),  0);
    brdcst(Proc, Num_Proc, (char *)Coordinate_Frame_Coordinates, 9*Num_Coordinate_Frames*io_ws,        0);
    brdcst(Proc, Num_Proc, (char *)Coordinate_Frame_Tags,          Num_Coordinate_Frames*sizeof(char), 0);
  }

  /*
   * Allocate Temporary Arrays that are only used in el_exoII_io.c (Note: Calls
   * to array_alloc are minimized by using pointer arithmetic. Resulting memory
   * structure may be assumed during later broadcast routines)
   */

  if (Num_Elem_Blk > 0) {
    Num_Elem_In_Blk    = (int *) array_alloc (__FILE__, __LINE__, 1,
                                              (4 * Num_Elem_Blk),
                                              sizeof(int));
    Num_Nodes_Per_Elem = Num_Elem_In_Blk    + Num_Elem_Blk;
    Num_Attr_Per_Elem  = Num_Nodes_Per_Elem + Num_Elem_Blk;
    Elem_Blk_Ids       = Num_Attr_Per_Elem  + Num_Elem_Blk;
    Elem_Blk_Types     =
      (char **) array_alloc (__FILE__, __LINE__, 2, Num_Elem_Blk,
                             MAX_STR_LENGTH + 1, sizeof(char));
    Elem_Blk_Names     =
      (char **) array_alloc (__FILE__, __LINE__, 2, Num_Elem_Blk,
                             max_name_length + 1, sizeof(char));
    Elem_Blk_Attr_Names = (char***) array_alloc(__FILE__, __LINE__, 1, Num_Elem_Blk,
						sizeof(char**));
  } else {
    fprintf(stderr,"ERROR, proc %d, Num_Elem_Blk = %d\n", Proc,
            Num_Elem_Blk);
    exit(1);
  }

  if (Num_Node_Set > 0)  {
    Node_Set_Ids          = (int *)array_alloc(__FILE__, __LINE__, 1,
                                               (3 * Num_Node_Set),
                                               sizeof(int));
    num_nodes_in_node_set = Node_Set_Ids          + Num_Node_Set;
    num_df_in_nsets       = num_nodes_in_node_set + Num_Node_Set;
    Node_Set_Names     =
      (char **) array_alloc (__FILE__, __LINE__, 2, Num_Node_Set,
                             max_name_length + 1, sizeof(char));
  } else {
    Node_Set_Ids = NULL;
    Node_Set_Names = NULL;
  }

  if (Num_Side_Set > 0) {
    Side_Set_Ids      = (int *)array_alloc(__FILE__, __LINE__, 1,
                                           (3 * Num_Side_Set), sizeof(int));
    num_elem_in_ssets = Side_Set_Ids      + Num_Side_Set;
    num_df_in_ssets   = num_elem_in_ssets + Num_Side_Set;
    Side_Set_Names     =
      (char **) array_alloc (__FILE__, __LINE__, 2, Num_Side_Set,
                             max_name_length + 1, sizeof(char));
  } else {
    Side_Set_Ids = NULL;
    Side_Set_Names = NULL;
  }

  /*
   * Process the Global Element Block IDs and associated information related to
   * element blocks (i.e., read them, broadcast them, and check them against
   * the input file).
   */
  if(Proc == 0)
    start_time = second();

  read_elem_blk_ids(mesh_exoid, max_name_length);

  if(Proc == 0) {
    printf("\tTime to read element block IDs: %.2f\n",
           second() - start_time);
  }

  /*
   * Process the Node Set IDs and associated information related to node sets
   * (i.e., read them, broadcast them, and check them against the input file).
   */
  if(Proc == 0)
    start_time = second();

  read_node_set_ids(mesh_exoid, num_nodes_in_node_set, num_df_in_nsets, max_name_length);

  if(Proc == 0) {
    printf("\tTime to read node set IDs: %.2f\n",
           second() - start_time);
  }

  /*
   * Process the Side Set IDs and associated information related to side sets
   * (i.e., read them, broadcast them, and check them against the input file).
   */
  if(Proc == 0)
    start_time = second();

  read_side_set_ids(mesh_exoid, num_elem_in_ssets, num_df_in_ssets, max_name_length);

  if(Proc == 0) {
    printf("\tTime to read side set IDs: %.2f\n",
           second() - start_time);
  }

  /*
   * Process the element block information.  Find out which element blocks have
   * elements needed by which processors.  Set-up and fill in vectors which
   * have length equal to the number of element blocks on the processor,
   * Proc_Num_Elem_Blk.
   */
  if(Proc == 0)
    start_time = second();

  extract_elem_blk();

  if(Proc == 0) {
    printf("\tTime to extract element block information: %.2f\n",
           second() - start_time);
  }

  /*
   * Read the mesh information from the exodus II file and broadcast it to the
   * processors.  NOTE: this will only read as much information as will fit in
   * the communication buffer.  It will then send the information to all the
   * processors so they can extract the information they need. Then, the next
   * block of mesh data will be read and sent, etc., until the entire mesh has
   * been read and distributed to the processors.  Note, that the load balance
   * information has already been distributed by this point.  Therefore, each
   * node already nodes which global nodes it owns, which elements it has
   * global nodes in that it owns, and which ghost nodes it needs information
   * about.
   *
   *        The following information is defined as the mesh information:
   *
   *            1 - Coodinate values of global nodes
   *            2 - Element connectivity graph
   *            3 - Attributes of the elements
   *            4 - Composition of Node Sets
   *            5 - Composition of Side Sets
   */

  /*      Read the coordinates, broadcast, and sort */

  if(Proc == 0)
    start_time = second();

  read_coord(mesh_exoid, io_ws, max_name_length);

  if(Proc == 0) {
    printf("\tTime to read nodal coordinates: %.2f\n",
           second() - start_time);
  }

  /*
   * Read the element connectivity and attributes information.  Broadcast it
   * and extract only the information that the current processor needs.
   */

  if(Proc == 0)
    start_time = second();

  read_elem_blk(mesh_exoid, io_ws);

  if(Proc == 0) {
    printf("\tTime to read element blocks: %.2f\n",
           second() - start_time);
  }

  /* Read the node sets, broadcast, and sort */

  if(Num_Node_Set > 0 && Proc == 0)
    start_time = second();

  if (Num_Node_Set > 0)
    read_node_sets(mesh_exoid, num_nodes_in_node_set, num_df_in_nsets, io_ws);

  if(Num_Node_Set > 0 && Proc == 0) {
    printf("\tTime to read node sets: %.2f\n",
           second() - start_time);
  }

  /* Assign the element types. */

  if(Proc == 0)
    start_time = second();

  create_elem_types();

  if(Proc == 0) {
    printf("\tTime to categorize element types: %.2f\n",
           second() - start_time);
  }

  /* Read the side sets, broadcast, and sort */

  if(Num_Side_Set > 0 && Proc == 0)
    start_time = second();

  if (Num_Side_Set > 0)
    read_side_sets(mesh_exoid, num_elem_in_ssets, num_df_in_ssets, io_ws);

  if(Num_Side_Set > 0 && Proc == 0) {
    printf("\tTime to read side sets: %.2f\n",
           second() - start_time);
  }

  if (Proc == 0) {

    /* Close the EXODUS II  mesh file */
    check_exodus_error(ex_close(mesh_exoid), "ex_close");
  }

/************* Output the information to the parallel disks *****************/
/*==========================================================================*/
/*==========================================================================*/

  /* Generate the processor to disk map */

  gen_disk_map(&PIO_Info, Proc_Info, Proc, Num_Proc);

  /* Generate the parallel exodus II file name */

  /* See if any '/' in the name.  IF present, isolate the basename of the file */
  if (strrchr(PIO_Info.Scalar_LB_File_Name, '/') != NULL) {
    /* There is a path separator.  Get the portion after the
     * separator
     */
    strcpy(cTemp, strrchr(PIO_Info.Scalar_LB_File_Name, '/')+1);
  } else {
    /* No separator; this is already just the basename... */
    strcpy(cTemp, PIO_Info.Scalar_LB_File_Name);
  }

  if (strlen(PIO_Info.Exo_Extension) == 0)
    add_fname_ext(cTemp, ".par");
  else
    add_fname_ext(cTemp, PIO_Info.Exo_Extension);

  for(iproc=0; iproc < Proc_Info[2]; iproc++) {
    gen_par_filename(cTemp, Par_Nem_File_Name, Proc_Ids[iproc], Proc_Info[0]);

    /* Stage the writes if specified in the input file */
    pdisk_stage_begin(&PIO_Info, Proc, Proc_Info, Proc_Ids, iproc,
                      &glob_pindx);

    /* Create the parallel Exodus II file for writing */
    if (Debug_Flag >= 7)
      printf("%sParallel mesh file name is %s\n", yo, Par_Nem_File_Name);
    else {
      if (iproc%10 == 0 || iproc == Proc_Info[2]-1)
        fprintf(stderr, "%d", iproc);
      else
        fprintf(stderr, ".");
    }

#if !defined(sun)
    if ((mesh_exoid=ex_create(Par_Nem_File_Name, EX_CLOBBER, &cpu_ws,
                              &io_ws)) == -1) {
#else
    if ((mesh_exoid=ex_create(Par_Nem_File_Name, EX_CLOBBER|EX_SHARE, &cpu_ws,
                              &io_ws)) == -1) {
#endif
      fprintf(stderr,"[%d] %sCould not create parallel Exodus II file\n",
              Proc, yo);
      exit(1);
    }

    ex_set_max_name_length(mesh_exoid, max_name_length);
    
    /* Set fill mode off... */
    {
    int old_fill;
    nc_set_fill(mesh_exoid, NC_NOFILL, &old_fill);
    }
    
    if (Debug_Flag >= 7)
      printf("%sParallel mesh file id is %d\n", yo, mesh_exoid);

    /* Write out a parallel mesh file local to each processor */
    write_parExo_data(mesh_exoid, max_name_length, iproc, io_ws, Num_QA_Recs+1, QA_Record,
                      Num_Info_Recs, Info_Record, Elem_Blk_Types,
                      Node_Set_Ids, Side_Set_Ids, Elem_Blk_Ids,
		      num_nodes_in_node_set, num_elem_in_ssets, Num_Elem_In_Blk,
		      E_Comm_Map[iproc], N_Comm_Map[iproc],
		      Node_Set_Names, Side_Set_Names, Elem_Blk_Names,
		      Elem_Blk_Attr_Names);

    /* Close the parallel exodus II file */
    if (ex_close(mesh_exoid) == -1) {
      fprintf(stderr, "%sCould not close the parallel Exodus II file\n",
              yo);
      exit(1);
    }

    /*
     * If staged writes are enabled then the current processor tells the next
     * processor in the chain to perform it's writes.
     */

    pdisk_stage_end(&PIO_Info, Proc, Proc_Info, Proc_Ids, iproc, glob_pindx);

  } /* End "for(iproc=0; iproc < Proc_Info[2]; iproc++)" */

  /*
   * The corresponding print_sync_start() is contained toward the bottom of
   * the file pe_write_parExo_info.c.
   */
  if (Debug_Flag >= 4) {
    print_sync_start(Proc, Num_Proc, TRUE);
    if (Proc == 0) {
      printf("\n\n\t\tTIMING TABLE FOR PROCESSORS\n");
      printf("===========================================================\n");
    }
    printf("Processor %d:\t[0]: %6.2f\n", Proc, PIO_Time_Array[0]);
    for (i1 = 1; i1 < 25; i1++) {
      printf("\t\t[%d]: %6.2f\n", i1, PIO_Time_Array[i1]);
    }
    printf("\nOutput rate: %6.2f kB/s\n", PIO_Time_Array[25]);
    print_sync_end(Proc, Num_Proc, TRUE);
  }
  else if (Debug_Flag >= 1)
    printf("\n\nProcessor %d output rate: %6.2f kB/s\n", Proc,
           PIO_Time_Array[25]);

 /*--------------------------------------------------------------------------*/
 /*--------------------------------------------------------------------------*/

  /* Free QA records */
  if (Proc == 0) {
    for (i1 = 0; i1 < 4*(Num_QA_Recs+1); i1++)
      safe_free((void **) &(QA_Record[i1]));

    safe_free((void **) &QA_Record);
  }

  if (Proc == 0) {
    for(i1=0; i1 < Num_Elem_Blk; i1++) {
      if (Num_Attr_Per_Elem[i1] > 0) {
	safe_free((void **) &(Elem_Blk_Attr_Names[i1]));
      }
    }
    safe_free((void **) &Elem_Blk_Attr_Names);
  }
  
  /* Free Coordinate Frames */
  if (Proc == 0 && Num_Coordinate_Frames > 0) {
    safe_free((void **) &Coordinate_Frame_Ids);
    if (io_ws <= sizeof(float)) {
      safe_free((void **) &Coordinate_Frame_Coordinates_sp);
    } else {
      safe_free((void **) &Coordinate_Frame_Coordinates_dp);
    }
    safe_free((void **) &Coordinate_Frame_Tags);
  }

  /* done with the Coordinate names */
  for(i1=0; i1 < Num_Dim; i1++)
    safe_free((void **) &(Coord_Name[i1]));

  /* Free some local arrays */
  safe_free((void **) &Num_Elem_In_Blk);
  safe_free((void **) &Elem_Blk_Types);
  safe_free((void **) &Node_Set_Ids);
  safe_free((void **) &Side_Set_Ids);

  safe_free((void **) &Node_Set_Names);
  safe_free((void **) &Side_Set_Names);
  safe_free((void **) &Elem_Blk_Names);

  /*
   * free up some other memory so that there is more
   * memory available for reading restart variables
   */
  for(iproc=0; iproc < Proc_Info[2]; iproc++) {
    safe_free((void **) &(Proc_Connect_Ptr[iproc]));
    safe_free((void **) &(Proc_Elem_Connect[iproc]));

    if (io_ws <= sizeof(float)) {
      safe_free((void **) &(Coor_sp[iproc]));
      if (Proc_Elem_Attr_sp[iproc])
        safe_free((void **) &(Proc_Elem_Attr_sp[iproc]));
      if (Proc_NS_Dist_Fact_sp[iproc])
        safe_free((void **) &(Proc_NS_Dist_Fact_sp[iproc]));
      if (Proc_SS_Dist_Fact_sp[iproc])
        safe_free((void **) &(Proc_SS_Dist_Fact_sp[iproc]));
    }
    else {
      safe_free((void **) &(Coor_dp[iproc]));
      if (Proc_Elem_Attr_dp[iproc])
        safe_free((void **) &(Proc_Elem_Attr_dp[iproc]));
      if (Proc_NS_Dist_Fact_dp[iproc])
        safe_free((void **) &(Proc_NS_Dist_Fact_dp[iproc]));
      if (Proc_SS_Dist_Fact_dp[iproc])
        safe_free((void **) &(Proc_SS_Dist_Fact_dp[iproc]));
    }

    if (N_Comm_Map[iproc])
      safe_free((void **) &(N_Comm_Map[iproc]));
    if (E_Comm_Map[iproc])
      safe_free((void **) &(E_Comm_Map[iproc]));
  }
  free(N_Comm_Map);
  free(E_Comm_Map);
  printf("\n");

  return;

} /* END of routine load_mesh () *********************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_elem_blk_ids(int mesh_exoid, int max_name_length)

/* This function reads part of the element block info from the EXODUS II file.
 * It reads all information having a length equal to the number of elements
 * blocks, specifically.
 *
 * The function then broadcasts this information to all processors.
 *
 * The function also calls check_matrl_blks, which checks the element block ids
 * in the input file and in the mesh file for consistency (and prints out a
 * table for a large enough value of Debug_Flag).
 */

{
  int   i;

  if (Proc == 0) {

    /* Get the Element Block IDs from the input file */
    check_exodus_error (ex_get_elem_blk_ids (mesh_exoid, Elem_Blk_Ids),
                        "ex_get_elem_blk_ids");

    check_exodus_error (ex_get_names(mesh_exoid, EX_ELEM_BLOCK, Elem_Blk_Names),
			"ex_get_elem_blk_ids");

    /*
     *     Get from the input file:
     *         Number of Elements               in each element block
     *         Number of nodes per element      in each element block
     *         Number of attributes per element in each element block
     *         The element type for elements    in each element block
     */
    for (i = 0; i < Num_Elem_Blk; i++) {

      check_exodus_error(ex_get_elem_block(mesh_exoid, Elem_Blk_Ids[i],
                                           Elem_Blk_Types[i],
                                           &Num_Elem_In_Blk[i],
                                           &Num_Nodes_Per_Elem[i],
                                           &Num_Attr_Per_Elem[i]),
                         "ex_get_elem_block");

      check_exodus_error(ne_get_elem_type(mesh_exoid, Elem_Blk_Ids[i],
                                          Elem_Blk_Types[i]),
                         "ne_get_elem_type");

      /* Convert element block types to lower case here */
      string_to_lower(Elem_Blk_Types[i], '\0');

      /* Allocate space for attribute names (if any), read and store. */
      if (Num_Attr_Per_Elem[i] > 0) {
	Elem_Blk_Attr_Names[i] = 
	  (char **) array_alloc (__FILE__, __LINE__, 2, Num_Attr_Per_Elem[i],
				 max_name_length + 1, sizeof(char));
	check_exodus_error(ex_get_attr_names(mesh_exoid, EX_ELEM_BLOCK,
					     Elem_Blk_Ids[i], Elem_Blk_Attr_Names[i]),
			   "ex_get_attr_names");
      } else {
	Elem_Blk_Attr_Names[i] = NULL;
      }
    }

  } /* End "if (Proc == 0)" */

  /*
   *             Broadcast the element block information using one
   *             broadcast, because these are in contiguous memory locations.
   *                       Num_Elem_In_Blk,
   *                       Num_Nodes_Per_Elem,
   *                       Num_Attr_Per_Elem
   *                       Elem_Blk_Ids
   */
  brdcst(Proc, Num_Proc, (char *)Num_Elem_In_Blk,
         (4*Num_Elem_Blk*sizeof(int)), 0);

  /*
   * (Note: arrays with dimension greater than 1, that were allocated with the
   * array_alloc function must be communicated with either a series of
   * broadcasts due to potential memory positioning conflicts or a broadcast
   * starting at an address after the pointer information.
   */
  brdcst(Proc, Num_Proc, Elem_Blk_Types[0],
         (Num_Elem_Blk *(MAX_STR_LENGTH+1)*sizeof(char)), 0);

  brdcst(Proc, Num_Proc, Elem_Blk_Names[0],
         (Num_Elem_Blk *(max_name_length+1)*sizeof(char)), 0);

  /*
   * It would be more efficient to pack this into a single array of
   * the correct size and do a single broadcast, but since we very
   * seldom run this on multiple processors anyway, lets just do it
   * the simple way for now and fix if/when we start running
   * nem_spread on multiple processors routinely and this becomes a
   * bottleneck. 
   */
  for (i = 0; i < Num_Elem_Blk; i++) {
    if (Num_Attr_Per_Elem[i] > 0) {
      brdcst(Proc, Num_Proc, (char *)Elem_Blk_Attr_Names[i],
	     (Num_Attr_Per_Elem[i] *(max_name_length+1)*sizeof(char)), 0);
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_node_set_ids(int mesh_exoid, int num_nodes_in_node_set[],
                              int num_df_in_nsets[], int max_name_length)

/*
 * This function reads part of the node set info from the EXODUS II file.  It
 * reads all information having a length equal to the number of node sets,
 * specifically.
 *
 * The function then broadcasts this information to all processors.
 *
 */

{
  int  i, error;

  if (Proc == 0) {

    if (Num_Node_Set > 0) {
      error = ex_get_node_set_ids(mesh_exoid, Node_Set_Ids);
      check_exodus_error(error, "ex_get_node_set_ids");

      error = ex_get_names(mesh_exoid, EX_NODE_SET, Node_Set_Names);
      check_exodus_error(error, "ex_get_node_set_ids");

      for (i = 0; i < Num_Node_Set; i++) {
        check_exodus_error(ex_get_node_set_param(mesh_exoid, Node_Set_Ids[i],
                                                 &num_nodes_in_node_set[i],
                                                 &num_df_in_nsets[i]),
                           "ex_get_node_set_param");
      }
    }

    /* Output debug info */
    if ((Debug_Flag > 1) && (Proc == 0)) {
      (void) printf("\n\n");
      print_line("=", 79);
      (void) printf("\tTABLE OF NODE SET ID\'s\n\n");
      (void) printf("Node_Set_Num   ID  Num_Nodes\n" );
      print_line("-", 79);

      if (Num_Node_Set > 0) {
        for (i = 0; i < Num_Node_Set; i++)
          (void) printf("%6d%11d%9d\n", i, Node_Set_Ids[i],
                        num_nodes_in_node_set[i]);
      }

      else {
        (void) printf("\tNO NODE SETS ARE DEFINED IN THE MESH FILE\n");
      }
      print_line("=", 79); printf("\n");
    }
  }

  /* Broadcast node set information to all processors (Note: Node_Set_Ids,
   * num_nodes_in_node_set and num_df_in_nsets are continguous in memory so
   * only one broadcast is necessary for communication of both arrays)
   */
  if (Num_Node_Set > 0) {
    brdcst(Proc, Num_Proc, (char *)Node_Set_Ids,
           3*Num_Node_Set*sizeof(int), 0);

    brdcst(Proc, Num_Proc, Node_Set_Names[0],
	   (Num_Node_Set *(max_name_length+1)*sizeof(char)), 0);
  }

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_side_set_ids(int mesh_exoid, int num_elem_in_ssets[],
                              int num_df_in_ssets[], int max_name_length)

/* This function reads part of the side set info from the EXODUS II file.  It
 * reads all information having a length equal to the number of side sets,
 * specifically.
 *
 * The function then broadcasts this information to all processors.
 *
 */

{
  int i, error;

  if (Proc == 0) {

    if (Num_Side_Set > 0) {
      error = ex_get_side_set_ids(mesh_exoid, Side_Set_Ids);
      check_exodus_error(error, "ex_get_side_set_ids");

      error = ex_get_names(mesh_exoid, EX_SIDE_SET, Side_Set_Names);
      check_exodus_error(error, "ex_get_side_set_ids");

      for (i = 0; i < Num_Side_Set; i++) {
        check_exodus_error(ex_get_side_set_param(mesh_exoid, Side_Set_Ids[i],
                                                 &num_elem_in_ssets[i],
                                                 &num_df_in_ssets[i]),
                           "ex_get_side_set_param");
      }
    }

    /* Output debug information */
    if ((Debug_Flag > 1) && (Proc == 0)) {
      printf("\n\n");
      print_line("=", 79);
      printf("\tTABLE OF SIDE SET ID\'s\n\n");
      printf("Side_Set_Num   ID   Number Elements\n");
      print_line ("-", 79);

      if (Num_Side_Set > 0) {
        for (i = 0; i < Num_Side_Set; i++)
          printf("%6d%11d  %11d\n", i, Side_Set_Ids[i],
                 num_elem_in_ssets[i]);
      }

      else {
        printf("\tNO SIDE SETS ARE DEFINED IN THE MESH FILE\n");
      }
      print_line("=", 79);
      printf("\n");
    }

  }

  /*
   * Broadcast the side-set information (Note: Side_Set_Ids,
   * num_elem_in_ssets and num_df_in_ssets are continguous in memory so
   * only one broadcast is necessary for communication of the three arrays)
   */
  if (Num_Side_Set > 0) {
    brdcst(Proc, Num_Proc, (char *)Side_Set_Ids, (3*Num_Side_Set*sizeof(int)),
           0);
    brdcst(Proc, Num_Proc, Side_Set_Names[0],
	   (Num_Side_Set *(max_name_length+1)*sizeof(char)), 0);

  }

  return;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_coord(int exoid, int io_ws, int max_name_length)

/* Function which reads the nodal coordinates information from an * EXODUS II
 * database for a given processor.
 *
 * Author: Scott Hutchinson (1421)
 *
 * ------------------------------------------------------------------------
 *
 *     Functions called:
 *     check_exodus_error -- function which handles the error code returned by
 *     calls to EXODUS II API routines.
 *
 * ------------------------------------------------------------------------
 */

{

  /* Local Variables */
  void   *x_coor = NULL, *y_coor = NULL, *z_coor = NULL, *coord_vector;
  int     i, j, icoord_size, num_coord_in_mesg, num_messages, ioffset;
  int     itotal_nodes, istart_node, iend_node, inode, *mesg_start;
  int     max_coord_per_mesg;
  int    *global_node_ids = NULL;
  float  *coord_sp;
  double *coord_dp;

  int    iproc;

  /* Function declarations */
  extern void  print_line(char *, int);

  /*************************** execution begins ******************************/

  /* check to see if the io_ws is smaller than the machine precision */
  if (io_ws < sizeof(float)) io_ws = sizeof(float);

  /*
   * Calculate the size of the coordinate space for each processors
   * coordinate matrix.
   */
  for(iproc=0; iproc < Proc_Info[2]; iproc++) {
    itotal_nodes = Num_Internal_Nodes[iproc] + Num_Border_Nodes[iproc] +
                   Num_External_Nodes[iproc];

    /* Allocate permament storage for the coordinates */
    if (io_ws == sizeof(float))
      Coor_sp[iproc] = (float **) array_alloc(__FILE__, __LINE__, 2, Num_Dim,
                                              itotal_nodes, sizeof(float));
    else
      Coor_dp[iproc] = (double **) array_alloc(__FILE__, __LINE__, 2, Num_Dim,
                                               itotal_nodes, sizeof(double));
  }

  /* Calculate the size of a single node's coordinates */
  icoord_size = Num_Dim * io_ws;

  /*
   * Break coordinate messages up into chunks that can fit in a single
   * processor's memory
   */
  num_messages = break_message_up(icoord_size, Num_Node, MAX_CHUNK_SIZE,
                                  &mesg_start);
  max_coord_per_mesg = mesg_start[1] - mesg_start[0];

  if (Debug_Flag > 1 && Proc == 0) {
    printf("\nNumber of coordinate messages: %d\nMaximum coordinates "
           "per message: %d\n\n", num_messages, max_coord_per_mesg);
  }

  /* Allocate temporary space to hold the longest message */

  coord_vector = array_alloc(__FILE__, __LINE__, 1, Num_Dim*max_coord_per_mesg,
                             io_ws);

  if (io_ws == sizeof(float)) coord_sp = (float *) coord_vector;
  else                        coord_dp = (double *) coord_vector;

  /* Read in the coordinates and broadcast to the processors */

  for (i = 0; i < num_messages; i++) {
    istart_node = mesg_start[i];
    iend_node   = mesg_start[i+1];
    num_coord_in_mesg = iend_node - istart_node;

    if (io_ws == sizeof(float)) {

      switch (Num_Dim) {

      case (3): z_coor = coord_sp + (2 * num_coord_in_mesg);
        /*FALLTHROUGH*/
      case (2): y_coor = coord_sp + num_coord_in_mesg;
        /*FALLTHROUGH*/
      case (1): x_coor = coord_sp;
      }
    }
    else {

      switch (Num_Dim) {

      case (3): z_coor = coord_dp + (2 * num_coord_in_mesg);
        /*FALLTHROUGH*/
      case (2): y_coor = coord_dp + num_coord_in_mesg;
        /*FALLTHROUGH*/
      case (1): x_coor = coord_dp;
      }
    }

    /* Read a slab of coordinate values from the Exodus II mesh file */
    if (Proc == 0) {
      check_exodus_error(ne_get_n_coord(exoid, (istart_node + 1),
                                        num_coord_in_mesg, x_coor,
                                        y_coor, z_coor),
                         "ne_get_n_coord");
    }

    /* Broadcast the slab of values to all of the processors */
    if (Num_Proc > 1) {
      brdcst_maxlen(Proc, Num_Proc, (char *) coord_vector,
		    (Num_Dim*num_coord_in_mesg*io_ws), 0);
    }

    for(iproc=0; iproc < Proc_Info[2]; iproc++) {

      itotal_nodes = Num_Internal_Nodes[iproc] + Num_Border_Nodes[iproc] +
                     Num_External_Nodes[iproc];

      /*
       * On each processor, extract the coordinates that the processor
       * needs.
       */

      switch(Num_Dim) {

      case 3:
        for (j = 0; j < itotal_nodes; j++) {
          inode = GNodes[iproc][j];
          if (inode >= istart_node && inode < iend_node) {
            ioffset = inode - istart_node;
            if (io_ws == sizeof(float)) {
              Coor_sp[iproc][0][j] = coord_sp[ioffset];
              Coor_sp[iproc][1][j] = coord_sp[ioffset + num_coord_in_mesg];
              Coor_sp[iproc][2][j] = coord_sp[ioffset + 2*num_coord_in_mesg];
            }
            else {
              Coor_dp[iproc][0][j] = coord_dp[ioffset];
              Coor_dp[iproc][1][j] = coord_dp[ioffset + num_coord_in_mesg];
              Coor_dp[iproc][2][j] = coord_dp[ioffset + 2*num_coord_in_mesg];
            }
          }
        }
        break;

      case 2:
        for (j = 0; j < itotal_nodes; j++) {
          inode = GNodes[iproc][j];
          if (inode >= istart_node && inode < iend_node) {
            ioffset = inode - istart_node;
            if (io_ws == sizeof(float)) {
              Coor_sp[iproc][0][j] = coord_sp[ioffset];
              Coor_sp[iproc][1][j] = coord_sp[ioffset + num_coord_in_mesg];
            }
            else {
              Coor_dp[iproc][0][j] = coord_dp[ioffset];
              Coor_dp[iproc][1][j] = coord_dp[ioffset + num_coord_in_mesg];
            }
          }
        }
        break;

      case 1:
        for (j = 0; j < itotal_nodes; j++) {
          inode = GNodes[iproc][j];
          if (inode >= istart_node && inode < iend_node) {
            if (io_ws == sizeof(float))
              Coor_sp[iproc][0][j] = coord_sp[ioffset];
            else
              Coor_dp[iproc][0][j] = coord_dp[ioffset];
          }
        }
        break;
      }
    }

  } /* END of for (i = 0; i < num_messages; i++) */

  for(i=0; i < Num_Dim; i++)
    Coord_Name[i] = (char *)array_alloc(__FILE__, __LINE__, 1,
                                        max_name_length + 1, sizeof(char));

  /* Get the coordinate names */
  if(Proc == 0) {
    if(ex_get_coord_names(exoid, Coord_Name) < 0) {
      fprintf(stderr, "ERROR:Unable to obtain coordinate names\n");
      exit(1);
    }
  }

  /* Broadcast the coordinate names to each processor */
  for(i=0; i < Num_Dim; i++)
    brdcst(Proc, Num_Proc, Coord_Name[i], (max_name_length+1)*sizeof(char), 0);

  /*  Output the Coordinate Information if Debug_Flag is large enough */

  if (Debug_Flag > 6) {
    print_sync_start (Proc, Num_Proc, TRUE);
    for(iproc=0; iproc < Proc_Info[2]; iproc++) {

      itotal_nodes = Num_Internal_Nodes[iproc] + Num_Border_Nodes[iproc] +
                     Num_External_Nodes[iproc];

      print_line ("=", 79);
      printf("Coordinates of Nodes on Proc %d for Proc %d\n", Proc,
             Proc_Ids[iproc]);

      print_line ("-", 79);
      printf("\tNode_Num | Node_Type |    Coordinates\n");
      print_line ("-", 79);

      for (i = 0; i < itotal_nodes; i++) {
        printf("\t %6d  | ", i);
        if (i < Num_Internal_Nodes[iproc])  printf("Internal  |");
        else if (i < (Num_Internal_Nodes[iproc]+Num_Border_Nodes[iproc]))
          printf("Border    |");
        else                         printf("External  |");
        for (j = 0; j < Num_Dim; j++) {
          if (io_ws == sizeof(float))
            printf(" %11.3e", Coor_sp[iproc][j][i]);
          else
            printf(" %11.3e", Coor_dp[iproc][j][i]);
        }
        printf("\n");
      }
      print_line ("=", 79);
    }
    print_sync_end (Proc, Num_Proc, TRUE);
  }

  safe_free ((void **) &coord_vector);
  safe_free ((void **) &mesg_start);

  /* Handle global node ids... */
  if (Proc == 0) {
    global_node_ids = (int *) array_alloc(__FILE__, __LINE__, 1, Num_Node,
					  sizeof(int));

    check_exodus_error(ex_get_node_num_map(exoid, global_node_ids),
		       "ex_get_node_num_map");
  }
  if (Num_Proc > 1) {
    brdcst_maxlen(Proc, Num_Proc, (char *)global_node_ids,
		  Num_Node*sizeof(int), 0);
  }

  /*
   * Check whether map is sequential (1..Num_Node). If it is, then it
   * provides no information and we don't need to store it in the
   * output databases.
   */
  {
    int sequential = 1;
    for (i=0; i < Num_Node; i++) {
      if (global_node_ids[i] != i+1) {
	sequential = 0;
	break;
      }
    }
    if (sequential == 0) {
      for(iproc=0; iproc < Proc_Info[2]; iproc++) {

	itotal_nodes = Num_Internal_Nodes[iproc] + Num_Border_Nodes[iproc] +  Num_External_Nodes[iproc];
	Proc_Global_Node_Id_Map[iproc] = (int *) array_alloc(__FILE__, __LINE__, 1, itotal_nodes,
							     sizeof(int));

	extract_global_node_ids(global_node_ids, Num_Node, iproc);
      }
    } else {
      /* Should be NULL already, but make it more clear */
      for(iproc=0; iproc < Proc_Info[2]; iproc++) {
	Proc_Global_Node_Id_Map[iproc] = NULL;
      }
    }
    safe_free((void **) &global_node_ids);
  }
} /* END of routine read_coord */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void extract_elem_blk (void)

/* Function which calculates the element block information for the current
 * processor, given the global element block information.  If an error is
 * tripped, the program will set a flag, and abort at the end of the routine.
 *
 * The following global variables are allocated and calculated in this routine:
 *
 *           Proc_Num_Elem_Blk = Number of element blocks on the current
 *                               processor
 *           GElem_Blks [Proc_Num_Elem_Blk]
 *                             = Map from the local element block number
 *                               to the global element block number.
 *           Proc_Nodes_Per_Elem [Proc_Num_Elem_Blk]
 *                             = Number of nodes per element for each
 *                               block on the current processor
 *           Proc_Elem_Blk_Ids [Proc_Num_Elem_Blk]
 *                             = Element block id's for the processor's
 *                               element blocks
 *           Proc_Elem_Blk_Types [Proc_Num_Elem_Blk]
 *                             = Element block types for the processor's
 *                               element blocks
 *                               (this is a unique integer number)
 *           Proc_Num_Attr [Proc_Num_Elem_Blk]
 *                             = Number of attributes for each block on the
 *                               current processor
 *           Proc_Num_Elem_In_Blk [Proc_Num_Elem_Blk]
 *                             = Number of elements in the processor's
 *                               element blocks
 *
 *    ------------------------------------------------------------------------
 *
 *       Functions called:
 *
 *       find_elem_block -- function which finds the element block which
 *                          owns each element in this processor.  In addition,
 *                          it determines a map from the current element to an
 *                          element block.
 *
 *    ------------------------------------------------------------------------
 */

{
  int     i, j, iglobal_blk, iproc;
  int    *proc_elem_blk = NULL;

  /* Element blocks local to the current processor */

#ifdef DEBUG
  extern void print_line(char *, int);
#endif

/**************************** execution begins ******************************/

  /*
   * Allocate temporary array for a mapping between the local element number
   * and the corresponding element block id for the element block that the
   * element belongs to.
   */

  for(iproc=0; iproc < Proc_Info[2]; iproc++) {

    proc_elem_blk = (int *) array_alloc(__FILE__, __LINE__, 1,
                                        Num_Internal_Elems[iproc] +
                                        Num_Border_Elems[iproc],
                                        sizeof(int));

    /* Find out which element block each element in this processor belongs to.
     * Fill this information into the temporary vector, proc_elem_blk.  Also,
     * calculate:
     *
     *              Proc_Num_Elem_Blk = Number of element blocks defined on
     *                                  the current processor.
     *              GElem_Blks        = Map from the local element block number
     *                                  to the global element block number.
     */
    find_elem_block(proc_elem_blk, iproc, Proc_Ids[iproc]);

    /* Allocate Permament integer arrays, which have lengths equal to the
     * number of element blocks defined on the processor.  This is done with a
     * single malloc call.  Initialize this space to zero.
     */
    if (Num_Elem_Blk > 0) {
      Proc_Nodes_Per_Elem[iproc]  = (int *)
            array_alloc(__FILE__, __LINE__, 1,
                        (4 * Num_Elem_Blk + Proc_Num_Elem_Blk[iproc]),
                        sizeof(int));
      Proc_Elem_Blk_Ids[iproc]    =  Proc_Nodes_Per_Elem[iproc] + Num_Elem_Blk;
      Proc_Elem_Blk_Types[iproc]  =  Proc_Elem_Blk_Ids[iproc]   + Num_Elem_Blk;
      Proc_Num_Attr[iproc]        =  Proc_Elem_Blk_Types[iproc] +
                                     Proc_Num_Elem_Blk[iproc];
      Proc_Num_Elem_In_Blk[iproc] =  Proc_Num_Attr[iproc]       + Num_Elem_Blk;

      /* Initialize */
      for (i = 0; i < (4*Num_Elem_Blk+Proc_Num_Elem_Blk[iproc]); i++)
        Proc_Nodes_Per_Elem[iproc][i] = 0;

    } else {
      fprintf(stderr, "ERROR Proc %d, Num_Elem_Blk = %d\n", Proc,
              Num_Elem_Blk);
      exit(1);
    }

    /*
     * Fill in the local processor element block arrays from a lookup from the
     * global element block arrays.  Note that the global element block arrays
     * will be freed at the end of the file.
     */
    for (i = 0; i < Proc_Num_Elem_Blk[iproc]; i++) {

      iglobal_blk                   = GElem_Blks[iproc][i];
      Proc_Nodes_Per_Elem[iproc][i] = Num_Nodes_Per_Elem[iglobal_blk];
      Proc_Elem_Blk_Ids[iproc][i]   = Elem_Blk_Ids[iglobal_blk];
      Proc_Num_Attr[iproc][i]       = Num_Attr_Per_Elem[iglobal_blk];

      /* Determine the element type integer id for this current element block.
       * The element type integer ID is a unique identifier (unique for all
       * element types that the program knows about).  It is determined from
       * the character string name, Elem_Blk_Types, and the number of nodes in
       * the element.
       */
      Proc_Elem_Blk_Types[iproc][i] = get_type(Elem_Blk_Types[iglobal_blk],
                                               Proc_Nodes_Per_Elem[iproc][i]);
    }

    /*
     * Determine the number of elements in each of the element blocks defined
     * on the local processor.
     */

    for (i = 0; i < Proc_Num_Elem_Blk[iproc]; i++) {
      for (j = 0; j < Num_Internal_Elems[iproc]+Num_Border_Elems[iproc]; j++) {
        if (proc_elem_blk[j] == Proc_Elem_Blk_Ids[iproc][i]) {
          (Proc_Num_Elem_In_Blk[iproc][i])++;
        }
      }
    }

    /* Sort GElems so that each element block is monotonic */
    j = 0;
    for(i=0; i < Proc_Num_Elem_Blk[iproc]; i++) {
      gds_qsort((GElems[iproc])+j, Proc_Num_Elem_In_Blk[iproc][i]);
      j += Proc_Num_Elem_In_Blk[iproc][i];
    }

    /* Free temporary vectors */
    safe_free((void **) &proc_elem_blk);

  } /* End "for(iproc=0; iproc < Proc_Info[2]; iproc++)" */

  if(Debug_Flag >= 5) {
    print_sync_start(Proc, Num_Proc, FALSE);
    for(iproc=0; iproc < Proc_Info[2]; iproc++) {

      /* Printout the Element Block Information defined on each Processor */
      print_line ("=", 79);
      printf("\t\tLocal Element Block information for Proc = %d on Proc %d\n",
             Proc_Ids[iproc], Proc);
      printf("\t\tNumber of Elem blocks on processor = %d\n",
             Proc_Num_Elem_Blk[iproc]);
      printf("%s%s\n",
             "Local_Block_Num  Global_Block_Num  Block_ID Nodes_Per_Elem ",
             "Num_Attributes  Elem_Blk_Type  Num_Elem_In_Blk "
             "Glb_Elm_In_Blk");
      print_line("-", 79);
      for (i = 0; i < Proc_Num_Elem_Blk[iproc]; i++ )
        printf("%4d\t\t%5d\t%8d\t%8d\t%8d\t%8d\t%8d\t%8d\n",
               i, GElem_Blks[iproc][i], Proc_Elem_Blk_Ids[iproc][i],
               Proc_Nodes_Per_Elem[iproc][i], Proc_Num_Attr[iproc][i],
               Proc_Elem_Blk_Types[iproc][i], Proc_Num_Elem_In_Blk[iproc][i],
               Num_Elem_In_Blk[GElem_Blks[iproc][i]]);
      print_line ("=", 79);
    }
    print_sync_end(Proc, Num_Proc, FALSE);
  }

} /* END extract_elem_blk() */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_elem_blk(int exoid, int io_ws)

/* Function which reads the element block information from an EXODUS II
 * database.  If an error is tripped, the program will set a flag, and abort at
 * the end of the routine.
 *
 *    ------------------------------------------------------------------------
 *
 *       extract_elem_connect
 *                        - Function which sorts through the element block
 *                          connectivity information and extracts the
 *                          information required by each processor.
 *
 *    ------------------------------------------------------------------------
 */

{

  /* Local variables */

  void   *elem_attr = NULL;
  float  *attr_sp;
  double *attr_dp;
  int     i, j, k, ielem, ielem_blk, iconnect_length;
#ifdef DEBUG
  int     ielem_count;
#endif
  int    *elem_blk = NULL, iattr_length = 0, ipos, iptr_count;
  int    *global_ids = NULL;
  int     icount, error_flag = 0;
  int     iconnect_size, iattr_size, num_elem_per_message;
  int     num_attr_per_message, iproc;
  int     num_elem_messages, num_attr_messages, num_elem_left_over;
  int     num_attr_left_over, istart_elem, iend_elem, istart_attr, iend_attr;
  int     local_ielem_blk;

  /* Function declarations */

  extern void  find_message_info(int , int , int *, int *, int *);
  extern void  print_line(char *, int);

/**************************** execution begins ******************************/

  /* check to see if the io_ws is smaller than the machine precision */
  if (io_ws < sizeof(float)) io_ws = sizeof(float);

  /*
   * Allocate memory for the vector of global element types in the entire
   * mesh. Used in sideset manipulations
   */
  GM_Elem_Types = (int *) array_alloc(__FILE__, __LINE__, 1, Num_Elem,
                                      sizeof(int));

  global_ids = (int *) array_alloc(__FILE__, __LINE__, 1, Num_Elem,
				   sizeof(int));

  for(iproc=0; iproc < Proc_Info[2]; iproc++) {

    /*
     * Build Proc_Connect_Ptr, a vector of pointers to the start of each
     * element's connectivity list in Proc_Elem_Connect.  The last entry in
     * Proc_Connect_Ptr is the length of Proc_Elem_Connect.
     */
    Proc_Connect_Ptr[iproc] = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                  Num_Internal_Elems[iproc] +
                                                  Num_Border_Elems[iproc] + 1,
                                                  sizeof(int));

    for (i = 0, iconnect_length = 0, iptr_count = 0;
         i < Proc_Num_Elem_Blk[iproc]; i++) {
      for (j = 0; j < Proc_Num_Elem_In_Blk[iproc][i]; j++) {
        Proc_Connect_Ptr[iproc][iptr_count++] = iconnect_length;
        iconnect_length += Proc_Nodes_Per_Elem[iproc][i];
      }
    }
    Proc_Connect_Ptr[iproc][iptr_count] = iconnect_length;

    /*
     * Allocate the processor's connectivity list vector
     * - This is a global vector
     */
    Proc_Elem_Connect[iproc] = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                   iconnect_length,
                                                   sizeof(int));

#ifdef DEBUG
    for (i = 0; i < iconnect_length; i++)
      Proc_Elem_Connect[iproc][i] = -1111111;
#endif

    /*
     * Allocate the processor's attribute list vector, if its length
     * is nonzero
     */
    for (i = 0, iattr_length = 0; i < Proc_Num_Elem_Blk[iproc]; i++)
      iattr_length += Proc_Num_Attr[iproc][i]*Proc_Num_Elem_In_Blk[iproc][i];
    if (iattr_length > 0) {
      if (io_ws == sizeof(float))
        Proc_Elem_Attr_sp[iproc] = (float *) array_alloc(__FILE__, __LINE__, 1,
                                                         iattr_length,
                                                         sizeof(float));
      else
        Proc_Elem_Attr_dp[iproc] = (double *) array_alloc(__FILE__, __LINE__,
                                                          1, iattr_length,
                                                          sizeof(double));
    }

  } /* End "for(iproc=0; iproc < Proc_Info[2]; iproc)" */

  icount = 0;

  /*
   * READ THE ELEMENT CONNECTIVITY AND ATTRIBUTE INFORMATION
   * FROM THE EXODUS II FILE FOR EACH ELEMENT BLOCK
   */
  for(ielem_blk=0; ielem_blk < Num_Elem_Blk; ielem_blk++) {

    if (Num_Elem_In_Blk[ielem_blk] > 0) {

      /*
       * Calculate the size of a single element's connectivity list and a
       * single element's attribute list
       */
      iconnect_size = Num_Nodes_Per_Elem[ielem_blk]*sizeof(int);
      iattr_size    = Num_Attr_Per_Elem[ielem_blk]*io_ws;

      /*
       * Determine how many element connectivity list and how many element
       * attributes can be send in single message
       */
      find_message_info(iconnect_size, Num_Elem_In_Blk[ielem_blk],
                        &num_elem_per_message, &num_elem_messages,
                        &num_elem_left_over);

      if (iattr_size > 0)
        find_message_info(iattr_size, Num_Elem_In_Blk[ielem_blk],
                          &num_attr_per_message, &num_attr_messages,
                          &num_attr_left_over);
      else
        num_attr_messages = 0;

      if (Debug_Flag > 1 && Proc == 0) {
        printf("\n\nMessage summary for Element Block number %d, ",
               ielem_blk);
        printf("having a block id of %d:\n", Elem_Blk_Ids[ielem_blk]);
        printf("\tNumber of messages needed for the element connectivity "
               "vector = %d\n", num_elem_messages);
        printf("\tNumber of elements per message = %d\n",
               num_elem_per_message);
        printf("\tNumber of nodes per element = %d\n",
               Num_Nodes_Per_Elem[ielem_blk]);
        printf("\tLength of each message = %d bytes\n",
               (int)(Num_Nodes_Per_Elem[ielem_blk] * num_elem_per_message *
               sizeof(int)));
        if (num_attr_messages > 0)
          printf("\tNumber of attribute messages: %d\n\tNumber "
                 "of attributes per message: %d\n\n",
                 num_attr_messages, num_attr_per_message);
      }

      /*
       * Allocate the arrays for reading from the ExodusII file and
       * broadcasting
       */
      elem_blk = (int *) array_alloc(__FILE__, __LINE__, 1,
                                     Num_Nodes_Per_Elem[ielem_blk] *
                                     num_elem_per_message, sizeof(int));
      if (num_attr_messages > 0)
        elem_attr = (void *) array_alloc(__FILE__, __LINE__, 1,
                                         Num_Attr_Per_Elem[ielem_blk]*
                                         num_attr_per_message, io_ws);
      else
        elem_attr = NULL;

      if (io_ws == sizeof(float)) attr_sp = (float *) elem_attr;
      else                        attr_dp = (double *) elem_attr;

      /*
       * Read in the element connectivity list for the current element block,
       * ielem_blk, from Proc 0, and then, broadcast it to all of the processors
       */
      for(i=0; i < num_elem_messages; i++) {

        if(Debug_Flag >= 2 && Proc == 0)
          printf("\telem block message: %d of %d\n", i+1, num_elem_messages);

        /* Initialize the element connectivity list to a value of -1.0 */
        for (j = 0; j < Num_Nodes_Per_Elem[ielem_blk]*num_elem_per_message;
             elem_blk[j++] = -1);

        istart_elem = i*num_elem_per_message;
        if (num_elem_left_over == 0 || i < num_elem_messages - 1)
          iend_elem = istart_elem + num_elem_per_message;
        else
          iend_elem = istart_elem + num_elem_left_over;

        if (Proc == 0) {

          if (num_elem_left_over == 0 || i < num_elem_messages - 1) {
            int el_type;
            check_exodus_error(ne_get_n_elem_conn(exoid,
                                                  Elem_Blk_Ids[ielem_blk],
                                                  (istart_elem + 1),
                                                  num_elem_per_message,
                                                  elem_blk),
                               "ne_get_n_elem_conn");
            if(Debug_Flag >= 2 && Proc == 0)
              printf("\t\tread connectivity\n");

            el_type = get_type(Elem_Blk_Types[ielem_blk],
                               Num_Nodes_Per_Elem[ielem_blk]);
            for(ielem=0; ielem < num_elem_per_message; ielem++) {
              GM_Elem_Types[icount++] = el_type;
            }

            if(Debug_Flag >= 2 && Proc == 0)
              printf("\t\tgot element types\n");

          }
          else {
            check_exodus_error(ne_get_n_elem_conn(exoid,
                                                  Elem_Blk_Ids[ielem_blk],
                                                  (istart_elem + 1),
                                                  num_elem_left_over,
                                                  elem_blk),
                               "ne_get_n_elem_conn");
            if(Debug_Flag >= 2 && Proc == 0)
              printf("\t\tread connectivity\n");

            for(ielem=0; ielem < num_elem_left_over; ielem++) {

              GM_Elem_Types[icount++] = get_type(Elem_Blk_Types[ielem_blk],
                                                 Num_Nodes_Per_Elem[ielem_blk]);
            }
            if(Debug_Flag >= 2 && Proc == 0)
              printf("\t\tgot element types\n");
          }

        }

        /* PRINT OUT THE ELEMENT CONNECTIVITY TABLE IF IN DEBUGGING MODE */
        if (Debug_Flag >= 6 && Proc == 0) {
          printf("\n\n\n");
          print_line("=", 79);
          printf("Printout of Element connectivity list obtained from "
                 "Exodus II file:\n");
          printf("\tGlobal element block number = %d\n", ielem_blk);
          printf("\tElement ID number     = %d\n",
                 Elem_Blk_Ids[ielem_blk]);
          printf("\tMessage number        = %d\n", i);
          print_line("-", 79);
          ipos = 0;
          if(num_elem_left_over == 0 || i < num_elem_messages - 1) {
            for (j = 0; j <num_elem_per_message; j++) {
              printf("\t elem: %d, nodes:", j);
              for (k = 0; k < Num_Nodes_Per_Elem[ielem_blk]; k++)
                printf(" %d", elem_blk[ipos++]);
              printf("\n");
            }
          }
          else {
            for (j = 0; j < num_elem_left_over; j++) {
              printf("\t elem: %d, nodes:", j);
              for (k = 0; k < Num_Nodes_Per_Elem[ielem_blk]; k++)
                printf(" %d", elem_blk[ipos++]);
              printf("\n");
            }
          }
          print_line("=", 79);
        }

        /* Broadcast the Element connectivity matrix to all of the processors */
	if (Num_Proc > 1) {
	  brdcst_maxlen(Proc, Num_Proc, (char *) elem_blk,
			Num_Nodes_Per_Elem[ielem_blk]*
			num_elem_per_message*sizeof(int), 0);
	  if(Debug_Flag >= 2 && Proc == 0)
	    printf("\t\tbroadcast connectivity\n");
	}

        psync(Proc, Num_Proc);

        /*
         * On each processor, extract the element connectivity lists that the
         * processor needs
         */
        for(iproc=0; iproc < Proc_Info[2]; iproc++) {
          extract_elem_connect(elem_blk, ielem_blk, istart_elem,
                               iend_elem, &local_ielem_blk, iproc);
        }
        if(Debug_Flag >= 2 && Proc == 0)
          printf("\t\textract connectivity\n");

      } /* End "for(i=0; i < num_elem_messages; i++)" */

      /* Read in the element attribute lists and broadcast to the processors */
      for (i = 0; i < num_attr_messages; i++) {

        if(Debug_Flag >= 2 && Proc == 0)
          printf("\tattribute message: %d of %d\n", i+1, num_attr_messages);

        /* Initialize */
        if (io_ws == sizeof(float))
          for (j = 0; j < Num_Attr_Per_Elem[ielem_blk]*num_attr_per_message;
               attr_sp[j++] = 0.0);
        else
          for (j = 0; j < Num_Attr_Per_Elem[ielem_blk]*num_attr_per_message;
               attr_dp[j++] = 0.0);

        istart_attr = i*num_attr_per_message;

        if(num_attr_left_over == 0 || i < (num_attr_messages-1))
          iend_attr = istart_attr + num_attr_per_message;
        else
          iend_attr = istart_attr + num_attr_left_over;

        if (Proc == 0) {
          if (num_attr_left_over == 0 || i < (num_attr_messages - 1)) {
            check_exodus_error(ne_get_n_elem_attr(exoid,
                                                  Elem_Blk_Ids[ielem_blk],
                                                  (istart_attr + 1),
                                                  num_attr_per_message,
                                                  elem_attr),
                               "ne_get_n_elem_attr");

          }
          else {
            check_exodus_error(ne_get_n_elem_attr(exoid,
                                                  Elem_Blk_Ids[ielem_blk],
                                                  (istart_attr + 1),
                                                  num_attr_left_over,
                                                  elem_attr),
                               "ne_get_n_elem_attr");

          }

        }

        /* Broadcast to processors */
	if (Num_Proc > 1) {
	  brdcst_maxlen(Proc, Num_Proc, (char *)elem_attr,
			Num_Attr_Per_Elem[ielem_blk]*
			num_attr_per_message*io_ws, 0);
	  
	  psync(Proc, Num_Proc);
	}

        for(iproc=0; iproc < Proc_Info[2]; iproc++) {

          if(Debug_Flag > 6 && Proc == 0)
            printf("\t\tExtract attributes for processor %d\n",
                   Proc_Ids[iproc]);
          /*
           * On each processor, extract the element attributes that the
           * processor needs
           */
          extract_elem_attr(elem_attr, ielem_blk, istart_attr, iend_attr,
                            Num_Attr_Per_Elem[ielem_blk], iproc, io_ws);
        }
      }

      /*       Free vectors */
      safe_free((void **) &elem_blk);
      safe_free((void **) &elem_attr);

    } /* END "if (Num_Elem_In_Blk[ielem_blk] > 0)" */

  } /* End "for(ielem_blk=0; ielem_blk < Num_Elem_Blk; ielem_blk++)" */

  /* Broadcast the global element type array */
  if (Num_Proc > 1) {
    brdcst_maxlen(Proc, Num_Proc, (char *)GM_Elem_Types,
		  Num_Elem*sizeof(int), 0);
  }

  /* Handle global element ids... */
  if (Proc == 0) {
    check_exodus_error(ex_get_elem_num_map(exoid, global_ids),
		       "ex_get_elem_num_map");
  }
  if (Num_Proc > 1) {
    brdcst_maxlen(Proc, Num_Proc, (char *)global_ids,
		  Num_Elem*sizeof(int), 0);
  }

  /*
   * Check whether map is sequential (1..Num_Elem). If it is, then it
   * provides no information and we don't need to store it in the
   * output databases.
   */
  {
    int sequential = 1;
    for (i=0; i < Num_Elem; i++) {
      if (global_ids[i] != i+1) {
	sequential = 0;
	break;
      }
    }
    if (sequential == 0) {
      for(iproc=0; iproc < Proc_Info[2]; iproc++) {

	Proc_Global_Elem_Id_Map[iproc] = (int *) array_alloc(__FILE__, __LINE__, 1,
							Num_Internal_Elems[iproc] +
							Num_Border_Elems[iproc],
							sizeof(int));

	extract_global_element_ids(global_ids, Num_Elem, iproc);
      }
    } else {
      /* Should be NULL already, but make it more clear */
      for(iproc=0; iproc < Proc_Info[2]; iproc++) {
	Proc_Global_Elem_Id_Map[iproc] = NULL;
      }
    }
    safe_free((void **) &global_ids);
  }
#ifdef DEBUG
  if (Debug_Flag > 6) {
    print_sync_start(Proc, Num_Proc, TRUE);

    for(iproc=0; iproc < Proc_Info[2]; iproc++) {
      ipos = ielem_count = 0;
      printf("\n\n\n");
      print_line("=", 79);
      printf("Printout of Element connectivity lists on proc %d "
             "for proc %d\n", Proc, Proc_Ids[iproc]);
      printf("\t Number of element blocks on the current processor = %d\n",
             Proc_Num_Elem_Blk[iproc]);
      printf("\t\tLocal_block_ID Num_Elem_In_Blk Nodes_Per_Elem_In_Blk\n");
      printf("\t\t----------------------------------------------------\n");
      for (i = 0; i < Proc_Num_Elem_Blk[iproc]; i++)
        printf("\t\t\t %d   \t   %d \t\t   %d\n",
               i, Proc_Num_Elem_In_Blk[iproc][i],
               Proc_Nodes_Per_Elem[iproc][i]);
      printf("\t\t----------------------------------------------------\n");
      print_line("-", 79);
      for (i = 0; i < Proc_Num_Elem_Blk[iproc]; i++) {
        printf("\n\n\tOutput of local Element block %d\n", i);
        print_line("-", 79);
        for (j = 0; j < Proc_Num_Elem_In_Blk[iproc][i]; j++) {
          printf("\t elem: %d (%d), nodes:", ielem_count,
                 GElems[iproc][ielem_count]);
	  ielem_count++;
          for (k = 0; k < Proc_Nodes_Per_Elem[iproc][i]; k++)
            printf(" %d", Proc_Elem_Connect[iproc][ipos++]);
          printf("\n");
        }
        print_line ("-", 79);
      }
      print_line ("=", 79);
    }
    print_sync_end(Proc, Num_Proc, TRUE);
  }
#endif

  if (error_flag) {
    fprintf(stderr, "read_elem_blk FAILURE (Proc = %d): error flag was "
            "set in routine\n", Proc);
    exit(1);
  }

} /* read_elem_blk */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void extract_global_element_ids(int global_ids[], int Num_Elem, int iproc)
{
  /*
   * Num_Elem -- number of elements in serial mesh.
   * global_ids[] (size Num_Elem) = global id corresponding to index
   *                                 of element in serial file 
   * Num_Internal_Elems[iproc] + Num_Border_Elems[iproc] = number of element on this processor 
   * GElems[iproc][i] = index of element in serial file corresponding
                        to i'th element on this processor. Both are zero-based?
  */
  int i;
  int num_local_element = Num_Internal_Elems[iproc] + Num_Border_Elems[iproc];
  for (i=0; i < num_local_element; i++) {
    Proc_Global_Elem_Id_Map[iproc][i] = global_ids[GElems[iproc][i]];
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void extract_global_node_ids(int global_ids[], int Num_Node, int iproc)
{
  /*
   * Num_Elem -- number of elements in serial mesh.
   * global_ids[] (size Num_Elem) = global id corresponding to index
   *                                 of element in serial file 
   * Num_Internal_Elems[iproc] + Num_Border_Elems[iproc] = number of element on this processor 
   * GElems[iproc][i] = index of element in serial file corresponding
                        to i'th element on this processor. Both are zero-based?
  */
  int i;
  int num_local_node = Num_Internal_Nodes[iproc] + Num_Border_Nodes[iproc] +
                       Num_External_Nodes[iproc];
  for (i=0; i < num_local_node; i++) {
    Proc_Global_Node_Id_Map[iproc][i] = global_ids[GNodes[iproc][i]];
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int extract_elem_connect(int elem_blk[], int icurrent_elem_blk,
                                int istart_elem, int iend_elem,
                                int *local_ielem_blk, int indx)

     /*
      * Function which extracts the element connectivity list from the current
      * element connectivity message.  The results are put into the global
      * vector, Proc_Elem_Connect.  This function is run on all processors,
      * include Proc = 0
      *
      *       Author:          Scott Hutchinson (1421)
      *
      *      ------------------------------------------------------------------
      *
      *      input
      *    -----------
      *       elem_blk[] = Vector containing the element connectivity message.
      *                    This message consists of the either the whole or
      *                    part of the element connectivity array for a single
      *                    element block.
      *
      *      istart_elem = This is the starting value of the element number for
      *                    the element connectivity message.  i.e., the first
      *                    message in an element block will have
      *                    istart_elem = 0.  The second message for that
      *                    element block will have the value:
      *                         istart_elem = num_elem_per_message.
      *      iend_elem   = This is the value + 1 of the last element number for
      *                    the element connectivity message.  i.e., the first
      *                    message in an element block will have
      *                          iend_elem = num_elem_per_message
      *      icurrent_elem_blk = Value of the global element block number.
      *
      *     Return Value
      *   -------------------
      *     The function returns the number of elements in the current message
      *    that are defined on the current processor.
      *      ------------------------------------------------------------------
      */

{

  /* Local variables */

  int i, j, iglobal_elem;
  int iglobal_begin, iglobal_end;
  int iglobal_pos, ielem_blk, ipos, ibegin, iend;
  int iglobal_offset, iproc_offset, iproc_start, num_elem;
  int count_hits = 0, found = FALSE;

  extern int in_list_mono(int, int *, int, int vector[]);

  /***************************** Execution begins *****************************/

  /* Match the Element Block Id of the global block with the local block number
   * The end result is the currect value of ielem_blk, Store this value in the
   * output variable, local_ielem_blk.
   */
  for (ielem_blk = 0; ielem_blk < Proc_Num_Elem_Blk[indx]; ielem_blk++)
    if (Proc_Elem_Blk_Ids[indx][ielem_blk] ==
        Elem_Blk_Ids[icurrent_elem_blk]) {
      *local_ielem_blk = ielem_blk;
      found = TRUE;

      /* Some elements in the current element block passed in from
         read_elem_blk may be on this processor */

      /* Calculate the number of elements in the current message */
      num_elem = iend_elem - istart_elem;

      /* Calculate iglobal_offset - The sum of the elements in all the global
       *                            element blocks preceding this one
       */
      iglobal_offset = 0;
      for (i = 0; i < GElem_Blks[indx][ielem_blk]; i++)
        iglobal_offset += Num_Elem_In_Blk[i];

      /* Calculate iproc_offset - The sum of the elements in all the
       *                          processor's
       *                          element blocks preceding this one.
       * Calculate iproc_start  - The starting position in the connectivity
       *                          vector for this element block,
       *                          Proc_Elem_Connect[]
       */
      iproc_offset = iproc_start = 0;
      for (i = 0; i < ielem_blk; i++) {
        iproc_offset += Proc_Num_Elem_In_Blk[indx][i];
        iproc_start  += Proc_Num_Elem_In_Blk[indx][i] *
          Proc_Nodes_Per_Elem[indx][i];
      }

      iglobal_begin = iglobal_offset + istart_elem;
      iglobal_end   = iglobal_begin  + num_elem;

      ibegin       = iproc_offset;
      iend         = iproc_offset + Proc_Num_Elem_In_Blk[indx][ielem_blk];

      /* OPTIMIZE: Check that max GElems is >min < max...*/
      /* Note that GElems is sorted, so we can check that the
       * ranges overlap and skip the loop if they don't */
      if (GElems[indx][ibegin] < iglobal_end  &&
          GElems[indx][iend-1] >= iglobal_begin) {

        for (i=ibegin; i < iend; i++) {
          iglobal_elem = GElems[indx][i];

          /* GElems is sorted, so if we are outside the iglobal range,
           * break out of this loop...
           */
          if (iglobal_elem > iglobal_end)
            break;

          if (iglobal_begin <= iglobal_elem && iglobal_elem < iglobal_end) {
            count_hits ++;

            ipos = iproc_start + (i - iproc_offset) *
              Proc_Nodes_Per_Elem[indx][ielem_blk];

            iglobal_pos = (iglobal_elem - iglobal_begin) *
              Proc_Nodes_Per_Elem[indx][ielem_blk];

            assert(iglobal_pos < num_elem*Proc_Nodes_Per_Elem[indx][ielem_blk]);
            /*
             * Store the connectivity information for the current element into
             * Proc_Elem_Connect [].  At the same time, Decrement by one the
             * value of all nodes in the connectivity matrix.
             */
            for (j = 0; j < Proc_Nodes_Per_Elem[indx][ielem_blk]; j++) {
              Proc_Elem_Connect[indx][ipos++] = elem_blk[iglobal_pos++] - 1;
            }
          }
        }
      }
    } /* END if (Proc_Elem_Blk_Ids[ielem_blk]== ... */

  if (!found)
    *local_ielem_blk = -1;

  return count_hits;

} /* extract_elem_connect */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void extract_elem_attr(void *elem_attr, int icurrent_elem_blk,
                              int istart_elem, int iend_elem,
                              int natt_p_elem, int indx, int io_ws)
     /*
      * Function which extracts the element connectivity list from the current
      * element connectivity message.  The results are put into the global
      * vector, Proc_Elem_Connect.  This function is run on all processors,
      * include Proc = 0
      *
      *       Author:          Scott Hutchinson (1421)
      *
      *      ------------------------------------------------------------------
      *
      *      input
      *    -----------
      *       elem_attr  = Vector containing the element attributes for this message.
      *                    This message consists of the either the whole or
      *                    part of the element attribute array for a single
      *                    element block.
      *
      *      istart_elem = This is the starting value of the element number for
      *                    the element attribute message.  i.e., the first
      *                    message in an element block will have
      *                    istart_elem = 0.  The second message for that
      *                    element block will have the value:
      *                         istart_elem = num_elem_per_message.
      *      iend_elem   = This is the value + 1 of the last element number for
      *                    the element attribute message.  i.e., the first
      *                    message in an element block will have
      *                          iend_elem = num_elem_per_message
      */

{

  /* Local variables */

  int i, j, iglobal_elem;
  int iglobal_begin, iglobal_end;
  int iglobal_pos, ielem_blk, ipos, ibegin, iend;
  int iglobal_offset, iproc_offset, iproc_start, num_elem;
  float  *attr_sp;
  double *attr_dp;

  extern int in_list_mono(int, int *, int, int vector[]);

  /***************************** Execution begins *****************************/

  /* check to see if the io_ws is smaller than the machine precision */
  if (io_ws < sizeof(float)) io_ws = sizeof(float);

  if (io_ws == sizeof(float)) attr_sp = (float *) elem_attr;
  else                        attr_dp = (double *) elem_attr;

  /* Match the Element Block Id of the global block with the local block number
   */
  for (ielem_blk = 0; ielem_blk < Proc_Num_Elem_Blk[indx]; ielem_blk++)
    if (Proc_Elem_Blk_Ids[indx][ielem_blk] ==
        Elem_Blk_Ids[icurrent_elem_blk]) {

      /* Some elements in the current element block passed in from
         read_elem_blk may be on this processor */

      /* Calculate the number of elements in the current message */
      num_elem = iend_elem - istart_elem;

      /* Calculate iglobal_offset - The sum of the elements in all the global
       *                            element blocks preceding this one
       */
      iglobal_offset = 0;
      for (i = 0; i < GElem_Blks[indx][ielem_blk]; i++)
        iglobal_offset += Num_Elem_In_Blk[i];

      /* Calculate iproc_offset - The sum of the elements in all the
       *                          processor's
       *                          element blocks preceding this one.
       * Calculate iproc_start  - The starting position in the connectivity
       *                          vector for this element block,
       *                          Proc_Elem_Connect[]
       */
      iproc_offset = iproc_start = 0;
      for (i = 0; i < ielem_blk; i++) {
        iproc_offset += Proc_Num_Elem_In_Blk[indx][i];
        iproc_start  += Proc_Num_Elem_In_Blk[indx][i] * Proc_Num_Attr[indx][i];
      }

      iglobal_begin = iglobal_offset + istart_elem;
      iglobal_end   = iglobal_begin  + num_elem;

      ibegin       = iproc_offset;
      iend         = iproc_offset + Proc_Num_Elem_In_Blk[indx][ielem_blk];

      /* OPTIMIZE: Check that max GElems is >min < max...*/
      /* Note that GElems is sorted, so we can check that the
       * ranges overlap and skip the loop if they don't */
      if (GElems[indx][ibegin] < iglobal_end  &&
          GElems[indx][iend-1] >= iglobal_begin) {

        for (i=ibegin; i < iend; i++) {
          iglobal_elem = GElems[indx][i];

          /* GElems is sorted, so if we are outside the iglobal range,
           * break out of this loop...
           */
          if (iglobal_elem > iglobal_end)
            break;

          if (iglobal_begin <= iglobal_elem && iglobal_elem < iglobal_end) {
            ipos = iproc_start + (i - iproc_offset) * natt_p_elem;

            iglobal_pos = (iglobal_elem - iglobal_begin) * natt_p_elem;

            assert(iglobal_pos < num_elem*natt_p_elem);
            /*
             * Store the connectivity information for the current element into
             * Proc_Elem_Connect [].  At the same time, Decrement by one the
             * value of all nodes in the connectivity matrix.
             */
	    if (io_ws == sizeof(float)) {
	      for (j = 0; j < natt_p_elem; j++)
		Proc_Elem_Attr_sp[indx][ipos++] = attr_sp[iglobal_pos++];
	    }
	    else {
	      for (j = 0; j < natt_p_elem; j++)
		Proc_Elem_Attr_dp[indx][ipos++] = attr_dp[iglobal_pos++];
	    }
          }
        }
      }
    } /* END if (Proc_Elem_Blk_Ids[ielem_blk]== ... */
} /* extract_elem_connect */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void find_elem_block(int *proc_elem_blk, int indx, int proc_for)

/* Function which finds the element block which owns each element on the
 * current processor.  In addition, a map from the local element block number
 * to the global element block number, *GElem_Blks, is created.
 * Proc_Num_Elem_Blk, the number of element blocks on the current processor,
 * is calculated This function is called by every processor.
 *
 * Author(s):          Scott Hutchinson (9221)
 *                     Gary Hennigan (9221)
 *
 *       Output:
 *
 *         Proc_Num_Elem_Blk  = Number of element blocks on the processor
 *                              (Global int variable)
 *
 *        *proc_elem_blk      = Vector of element block ids for the local
 *                              elements defined on the processor.
 *                              (Local int vector of length
 *                              Num_Internal_Elems + Num_Border_Elems)
 *
 *        *GElem_Blks         = Map from the local element block number to
 *                              the global element block number.
 *                              (Global int vector of length Proc_Num_Elem_Blk)
 */

{

  /* Local variables */

  int  i, j, found, icount = 0;
  int *elem_in_blk;           /* Boolean vector of length Num_Elem_Blk
                                 If the ith element block exists on the
                                 current processor, the ith entry is set to
                                 TRUE      */
  int *elem_blk_point;        /* Vector of integer offsets into the vector
                                 GElems.
                                 It has a length of Num_Elem_Blk+1.
                                 The ith entry points to the beginning of
                                 of the element map information for the
                                 first element in the ith element block.
                                 The (i+1)th entry points to the last element
                                 for the ith block in the vector GElem. */
  /***************************** execution begins ****************************/

  /* Allocate array and initialize  */
  elem_in_blk    = (int *) array_alloc(__FILE__, __LINE__, 1,
                                       (2 * Num_Elem_Blk + 1), sizeof(int));
  elem_blk_point = elem_in_blk + Num_Elem_Blk;

  for (i = 0; i < Num_Elem_Blk; elem_in_blk[i++] = FALSE)
    ;

  /*
   * Construct a vector of pointers to the beginning of element blocks in an
   * element connectivity list (vector)
   */
  elem_blk_point[0] = 0;
  for (i = 0; i < Num_Elem_Blk; i++)
    elem_blk_point[i+1] = elem_blk_point[i] + Num_Elem_In_Blk[i];

  /*
   * Find out which elements belonging to the current processor are in which
   * local element block.  Store this in *proc_elem_blk.
   */

  /* Internal Elements */
  if (check_monot(&GElems[indx][0], Num_Internal_Elems[indx])) {
    int tmp_cnt = Num_Internal_Elems[indx];
    j = 0;
    i = 0;
    while (i < tmp_cnt && j < Num_Elem_Blk) {
      while (i < tmp_cnt && GElems[indx][i] < elem_blk_point[j+1]) {
        assert(GElems[indx][i] >= elem_blk_point[j]);
        proc_elem_blk[i++] = j;
        elem_in_blk[j] = TRUE;
      }
      j++;
    }
  } else {
    for (i = 0; i < Num_Internal_Elems[indx]; i++) {
      found = FALSE;
      for (j = 0; j < Num_Elem_Blk && !found; j++) {
        if (GElems[indx][i] <  elem_blk_point[j+1] &&
            GElems[indx][i] >= elem_blk_point[j]         )  {
          proc_elem_blk[i] = j;
          elem_in_blk[j]   = found = TRUE;
        }
      }
      if (!found) {
        fprintf(stderr, "find_elem_block: Error!  (Proc = %d):\n", Proc);
        fprintf(stderr, "\tElement %d not found in any element "
                "block.\n", i);
        exit(1);
      }
    }
  }

  /* Border Elements */
  if (check_monot(&GElems[indx][Num_Internal_Elems[indx]],
                               Num_Border_Elems[indx])) {
    int tmp_cnt = Num_Internal_Elems[indx] + Num_Border_Elems[indx];
    j = 0;
    i = Num_Internal_Elems[indx];
    while (i < tmp_cnt && j < Num_Elem_Blk) {
      while (i < tmp_cnt && GElems[indx][i] < elem_blk_point[j+1]) {
        assert(GElems[indx][i] >= elem_blk_point[j]);
        proc_elem_blk[i++] = j;
        elem_in_blk[j] = TRUE;
      }
      j++;
    }
  } else {
    int tmp_cnt = Num_Internal_Elems[indx] + Num_Border_Elems[indx];
    for (i = Num_Internal_Elems[indx]; i < tmp_cnt; i++) {
      found = FALSE;
      for (j = 0; j < Num_Elem_Blk && !found; j++) {
        if (GElems[indx][i] <  elem_blk_point[j+1] &&
            GElems[indx][i] >= elem_blk_point[j]         )  {
          proc_elem_blk[i] = j;
          elem_in_blk[j]   = found = TRUE;
        }
      }
      if (!found) {
        fprintf(stderr, "find_elem_block: Error!  (Proc = %d):\n", Proc);
        fprintf(stderr, "\tElement %d not found in any element "
                "block.\n", i);
        exit(1);
      }
    }
  }


  /*
   * Reorder GElems based on the element block information. This is necessary
   * to insure that internal and border elements that are in the same element
   * block are sequentially contained in GElems, which is assumed in some
   * later operations.
   */

  sort_int_int(Num_Internal_Elems[indx]+Num_Border_Elems[indx],
                proc_elem_blk, GElems[indx]);

  /* Now change proc_elem_blk to be a list of global element block IDs */

  for (i = 0; i < Num_Internal_Elems[indx]+Num_Border_Elems[indx]; i++) {
    j = proc_elem_blk[i];
    proc_elem_blk[i] = Elem_Blk_Ids[j];
  }

  /* Count the number of element blocks defined on this processor */
  Proc_Num_Elem_Blk[indx] = 0;
  for (i = 0; i < Num_Elem_Blk; i++)
    if (elem_in_blk[i]) Proc_Num_Elem_Blk[indx]++;

  /*
   * Create a map from the current processor's element block number to the
   * 'global' element block number, GElem_Blks.  This is a permament Global
   * Map.
   */
  GElem_Blks[indx] = (int *) array_alloc(__FILE__, __LINE__, 1, Num_Elem_Blk,
                                         sizeof(int));

  for (i = 0, icount = 0; i < Num_Elem_Blk; i++)
    if (elem_in_blk[i]) GElem_Blks[indx][icount++] = i;

  /* Free temporary arrays */
  safe_free ((void **) &elem_in_blk);

  return;

} /* END of routine find_elem_block() ****************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_node_sets(int exoid, int *num_nodes_in_node_set,
                           int *num_df_in_nsets, int io_ws)

/* Function which reads the node sets information from an EXODUS II database.
 * It read information in chunks.  Then, it broadcasts the chunk to all
 * processors.  Each processor, then searches for the chunk for what it is
 * responsible for.  This function is called by all processors.
 *
 *     Global Variables which are set by this routine
 *     -----------------------------------------------
 *
 *       Proc_Num_Node_Sets = Number of node sets on the current processor
 *
 *       GNode_Sets [Proc_Num_Node_Sets]
 *                          = Mapping between the local node set number and the
 *                            global node set number.
 *
 *       Proc_NS_Ids [Proc_Num_Node_Sets]
 *                          = Set IDs for the node sets defined on the current
 *                            processor
 *
 *       Proc_NS_Count [Proc_Num_Node_Sets]
 *                          = Number of nodes in the node sets defined on
 *                            the current processor
 *
 *       Proc_NS_Pointers [Proc_Num_Node_Sets]
 *                          = Pointers into the node set list for the node sets
 *                            defined on the current processor.
 *
 *       Proc_NS_List_Length = Length of the local node set list.
 *
 *       Proc_NS_List [Proc_NS_List_Length]
 *                          = Concatenated node set list for the current
 *                            processor.
 */

{
  char *yo="read_node_sets";

  int      i, j, k, kbeg, ipos, iproc, iss_size, brdcst_size;
  int      num_messages, num_left_over, num_node_per_message, imess, istart_ns;
  int     *ns_cntr, *first_message, *ns_on_proc, *proc_num_ns;

  int     *node_set;          /* the current node set's node number list     */
  void    *node_set_df;       /* the current node set dist. factors.         */
  int     *proc_ns_pointer;   /* pointers into 'node_set' for the nodes that
                                 are common to both 'node_set' and the internal
                                 and border nodes for the current processor  */
  int   ***proc_list;         /* the node number lists for the current node
                                 set for the current processor               */

  int   ***proc_gmap_list;    /*  Vector of pointers to nodeset global nodeeset node map lists.
                                  There is one pointer for each nodeset     */

  float  ***proc_list_df_sp;  /* the dist. factors for the current processor */
  double ***proc_list_df_dp;  /* must handle both single/double precision    */

  int    **proc_list_pointer; /* pointers into the concatenated node set list
                                 which locate the start of the node sets     */
  int    **ns_proc2glob;      /* mapping from the processor node set numbers
                                 to the global node set numbers              */
  int    **proc_list_count;   /* the number of nodes in each processor's node
                                 set                                         */
  int **proc_list_df_count;   /* the number of df's in each processor's node
                                 set (usually either zero or proc_list_count)*/
  int     *list_length = NULL; /* length of the processor's node-set list    */

  /* generic pointers for single/double precision */
  void    *ptr;
  float   *ptr_sp;
  double  *ptr_dp;

  int *gmap;
  int **gmap_dups = NULL;
  int  dups = 0;
  
  extern void find_message_info (int , int , int *, int *, int *);
  extern int  int_cmp           (int *i1, int *i2);

  /*************************** Execution begins ******************************/

  /* check to see if the io_ws is smaller than the machine precision */
  if (io_ws < sizeof(float)) io_ws = sizeof(float);

  /* Allocate arrays */
  list_length       = (int *)array_alloc(__FILE__, __LINE__, 1, 5*Proc_Info[2],
                                         sizeof(int));
  proc_num_ns       = list_length   + Proc_Info[2];
  ns_cntr           = proc_num_ns   + Proc_Info[2];
  first_message     = ns_cntr       + Proc_Info[2];
  ns_on_proc        = first_message + Proc_Info[2];

  ns_proc2glob       = (int **) array_alloc(__FILE__, __LINE__, 1,
                                            4*Proc_Info[2], sizeof(int *));
  proc_list_pointer  = ns_proc2glob      + Proc_Info[2];
  proc_list_count    = proc_list_pointer + Proc_Info[2];
  proc_list_df_count = proc_list_count   + Proc_Info[2];

  proc_list         = (int ***)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2], sizeof(int **));
  proc_gmap_list    = (int ***)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2], sizeof(int **));

  if (io_ws == sizeof(float))
    proc_list_df_sp    = (float ***)array_alloc(__FILE__, __LINE__, 1,
                                                Proc_Info[2],
                                                sizeof(float **));
  else
    proc_list_df_dp    = (double ***)array_alloc(__FILE__, __LINE__, 1,
                                                 Proc_Info[2],
                                                 sizeof(double **));

  /* Initialize */
  for(iproc=0; iproc < Proc_Info[2]; iproc++) {
    ns_proc2glob[iproc] = (int *)array_alloc(__FILE__, __LINE__, 1,
                                             4*Num_Node_Set, sizeof(int));
    proc_list_pointer[iproc]  = ns_proc2glob[iproc]      + Num_Node_Set;
    proc_list_count[iproc]    = proc_list_pointer[iproc] + Num_Node_Set;
    proc_list_df_count[iproc] = proc_list_count[iproc]   + Num_Node_Set;

    proc_list[iproc]      = (int **) array_alloc(__FILE__, __LINE__, 1, Num_Node_Set, sizeof(int *));
    proc_gmap_list[iproc] = (int **) array_alloc(__FILE__, __LINE__, 1, Num_Node_Set, sizeof(int *));

    if (io_ws == sizeof(float))
      proc_list_df_sp[iproc] = (float **)array_alloc(__FILE__, __LINE__, 1,
                                                     Num_Node_Set,
                                                     sizeof(float *));
    else
      proc_list_df_dp[iproc] = (double **)array_alloc(__FILE__, __LINE__, 1,
                                                      Num_Node_Set,
                                                      sizeof(double *));

    Proc_Num_Node_Sets[iproc] = 0;

    list_length[iproc] = 0;
  }

  /* Create a map from 'global node id' to processor */
  {
    gmap = (int *) array_alloc(__FILE__, __LINE__, 1, Num_Node, sizeof(int));
    
    /* Initialize to -1 */
    for (i=0; i < Num_Node; i++)
      gmap[i] = -1;
    
    /* Fill it in.... If run in parallel, only fills with the
     * nodes being processed on this processor, not entire model
     * Note that since nodes can be on multiple processors, we can't
     * normally do a simple one-to-one map of node to it
     * processor. Instead, we try to do a map from node->processor.
     * If there are multiple nodes, then we flag it in gmap by setting
     * the value to num_procs+first_value where first_value is the
     * lowest processor that the node appears on.  We also count the
     * number of "duplicates" or nodes on multiple processors.
     *
     * We then allocate the "gmap_dups" array to hold all nodes that
     * are on multiple processors.  This is a 2D array which holds the
     * node id (0-based) and the processor(s) it is on.
     */
    for (iproc=0; iproc < Proc_Info[2]; iproc++) {
      int size = Num_Internal_Nodes[iproc]+Num_Border_Nodes[iproc];
      for (i=0; i < size; i++) {
	int ind = GNodes[iproc][i];
	if (gmap[ind] == -1)
	  gmap[ind] = iproc;
	else {
	  if (gmap[ind] < Proc_Info[2]) {
	    gmap[ind] += Proc_Info[2];
	  }
	  dups ++;
	}
      }
    }
    if (dups > 0) {
      gmap_dups = (int **) array_alloc(__FILE__, __LINE__, 2, 2, dups, sizeof(int));
      dups = 0;
      for (iproc=1; iproc < Proc_Info[2]; iproc++) { /* Don't need to worry about proc 0 */
	int size = Num_Internal_Nodes[iproc]+Num_Border_Nodes[iproc];
	for (i=0; i < size; i++) {
	  int ind = GNodes[iproc][i];
	  if (gmap[ind] >= Proc_Info[2] && gmap[ind] != Proc_Info[2]+iproc) {
	    gmap_dups[0][dups] = ind;
	    gmap_dups[1][dups] = iproc;
	    dups++;
	  }
	}
      }
    }
  }

  /*-------------------------------------------------------------------------*/
  /*                    LOOP OVER THE NODE SETS                              */
  /*-------------------------------------------------------------------------*/

  for(i=0; i < Num_Node_Set; i++) {

    if (num_nodes_in_node_set[i] > 0) {

      for(iproc=0; iproc < Proc_Info[2]; iproc++) {
        first_message[iproc] = TRUE;
        proc_num_ns[iproc] = 0;
        ns_cntr[iproc] = 0;
        ns_on_proc[iproc] = FALSE;
      }

      if(num_df_in_nsets[i] > 0)
        iss_size = sizeof(int) + io_ws;
      else
        iss_size = sizeof(int);

      find_message_info(iss_size, num_nodes_in_node_set[i],
                        &num_node_per_message, &num_messages, &num_left_over);

      if (Debug_Flag > 1 && Proc == 0) {
        printf("\nMessage summary for Node Set number %d, with an ID of %d:\n",
               i, Node_Set_Ids[i]);
        printf("\tNumber of messages need for node set = %d\n",
               num_messages);
        printf("\tNumber of node IDs and dist. factors per message = %d\n",
               num_node_per_message);
        printf("\tLength of each message = %d\n",
               num_node_per_message*iss_size);
      }

      /*
       * Allocate temporary storage for the current message of the current node
       * set on each processor.
       */
      if(num_df_in_nsets[i] > 0)
        proc_ns_pointer = malloc(num_node_per_message*(2*sizeof(int) + io_ws));
      else
        proc_ns_pointer = malloc(num_node_per_message*(2*sizeof(int)));

      if(!proc_ns_pointer) {
        fprintf(stderr, "[%d-%s]: ERROR, insufficient memory!\n", Proc,
                yo);
        exit(1);
      }

      node_set = proc_ns_pointer + num_node_per_message;

      if(num_df_in_nsets[i] > 0)
        node_set_df     =  (void *) (node_set + num_node_per_message);

      /*---------------------------------------------------------------------*/
      /*                  LOOP OVER THE MESSAGES                             */
      /*---------------------------------------------------------------------*/

      for(imess=0; imess < num_messages; imess++) {

        istart_ns = imess*num_node_per_message;
        if (num_left_over != 0 && imess == num_messages - 1) {
          num_node_per_message = num_left_over;

          /* Pointers must also be adjusted for the last message */
          if(num_df_in_nsets[i] > 0)
            node_set_df = (void *) (node_set + num_node_per_message);
        }

        /* Read in the part of the node set that will fit in the message */
        if (Proc == 0) {
          check_exodus_error(ne_get_n_node_set(exoid, Node_Set_Ids[i],
                                               (istart_ns + 1),
                                               num_node_per_message,
                                               node_set),
                             "ne_get_n_node_set");

          if(num_df_in_nsets[i] > 0) {
            check_exodus_error(ne_get_n_node_set_df(exoid, Node_Set_Ids[i],
                                                    (istart_ns + 1),
                                                    num_node_per_message,
                                                    node_set_df),
                               "ne_get_n_node_set_df");
          }
        }

        if(num_df_in_nsets[i] > 0)
          brdcst_size = num_node_per_message*(sizeof(int) + io_ws);
        else
          brdcst_size = num_node_per_message*sizeof(int);

        /* Broadcast the node_set vector to all procs */
	if (Num_Proc > 1) {
	  brdcst_maxlen(Proc, Num_Proc, (char *) node_set, brdcst_size, 0);
	}

        /* Renumber nodes to start at node '0' instead of node '1' */
        for (j = 0; j < num_node_per_message; node_set[j++]--) ;

        /* Loop over the number or processors being handled */
	kbeg = 0;
        for(iproc=0; iproc < Proc_Info[2]; iproc++) {

	  if (dups > 0) {
	    for (; kbeg < dups; kbeg++) {
	      if (gmap_dups[1][kbeg] >= iproc) 
		break;
	    }
	  }

          /*
           * Find the intersection between the node set and Internal nodes and
           * Border nodes.  See above description of 'gmap' array
           * which explains the duplicates (nodes on multiple
           * proceessors). 
           */

	  ipos = 0;
	  for (j=0; j < num_node_per_message; j++) {
            if (gmap[node_set[j]] == iproc)
              proc_ns_pointer[ipos++] = j;
	    else if (gmap[node_set[j]] == iproc + Proc_Info[2])
              proc_ns_pointer[ipos++] = j;
	    else if (gmap[node_set[j]] >= Proc_Info[2]) {
	      /* Node is on multiple processors. See if on this
		 one... */
	      for (k = kbeg; k < dups && gmap_dups[1][k] == iproc; k++) {
		if (gmap_dups[0][k] == node_set[j] && gmap_dups[1][k] == iproc) {
		  proc_ns_pointer[ipos++] = j;
		  break;
		}
	      }
	    }
	  }

          /*
           * If the message node set and either the Internal node or the Border
           * node sets intersect, then add this intersection to the list
           */
          if (ipos > 0) {
            ns_on_proc[iproc] = TRUE;
            proc_num_ns[iproc] += ipos;

            ptr = NULL;

            /* Allocate and store node information in a temporary vector */
            if (first_message[iproc]) {

              proc_list[iproc][Proc_Num_Node_Sets[iproc]] =
                (int *)array_alloc(__FILE__, __LINE__, 1, proc_num_ns[iproc],
                                   sizeof(int));

              proc_gmap_list[iproc][Proc_Num_Node_Sets[iproc]] =
                (int *) array_alloc(__FILE__, __LINE__, 1,   proc_num_ns[iproc], sizeof(int));

              if(num_df_in_nsets[i] > 0)
                ptr = array_alloc(__FILE__, __LINE__, 1, proc_num_ns[iproc],
                                  io_ws);

              first_message[iproc] = FALSE;
            }
            else {

              proc_list[iproc][Proc_Num_Node_Sets[iproc]] =
                (int *) realloc(
                  proc_list[iproc][Proc_Num_Node_Sets[iproc]],
                  (proc_num_ns[iproc] * sizeof(int)));

              if(!(proc_list[iproc][Proc_Num_Node_Sets[iproc]])) {
                fprintf(stderr, "[%d-%s]: ERROR, insufficient memory\n",
                        Proc, yo);
                exit(1);
              }

              if((proc_gmap_list[iproc][Proc_Num_Node_Sets[iproc]] = (int *)
                  realloc(proc_gmap_list[iproc][Proc_Num_Node_Sets[iproc]],
                          (proc_num_ns[iproc]*sizeof(int)))) == NULL ) {
                fprintf(stderr, "ERROR: unable to realloc proc_gmap_list in "
			"%s\n", yo);
                exit(1);
              }

              if(num_df_in_nsets[i] > 0) {

                ptr = realloc(ptr, (proc_num_ns[iproc] * io_ws));

                if(!(ptr)) {
                  fprintf(stderr, "[%d-%s]: ERROR, insufficient memory\n",
                          Proc, yo);
                  exit(1);
                }

              }
            }

            if (io_ws == sizeof(float)) {
              proc_list_df_sp[iproc][Proc_Num_Node_Sets[iproc]] = (float *) ptr;
              ptr_sp = (float *) node_set_df;
            }
            else {
              proc_list_df_dp[iproc][Proc_Num_Node_Sets[iproc]] =
                (double *) ptr;
              ptr_dp = (double *) node_set_df;
            }

            for (j = 0; j < ipos; j++) {

              proc_list[iproc][Proc_Num_Node_Sets[iproc]][ns_cntr[iproc]] =
                node_set[proc_ns_pointer[j]];

	      proc_gmap_list[iproc][Proc_Num_Node_Sets[iproc]][ns_cntr[iproc]] =
		proc_ns_pointer[j]+istart_ns;

              if(num_df_in_nsets[i] > 0) {
                if (io_ws == sizeof(float)) {
                  proc_list_df_sp[iproc][Proc_Num_Node_Sets[iproc]]
                    [ns_cntr[iproc]] = ptr_sp[proc_ns_pointer[j]];
                }
                else {
                  proc_list_df_dp[iproc][Proc_Num_Node_Sets[iproc]]
                    [ns_cntr[iproc]] = ptr_dp[proc_ns_pointer[j]];
                }
              }

              (ns_cntr[iproc])++;
            }

          } /* End "if(intersection)" */

        } /* End "for(iproc=0; iproc < Proc_Info[2]; iproc++)" */

      } /* End "for(imess=0; imess < num_messages; imess++)" */

      /*
       * If any part of this node-set is on the processor, update the various
       * pointers, lengths, etc.
       */
      for(iproc=0; iproc < Proc_Info[2]; iproc++) {
        if (ns_on_proc[iproc] ) {
          ns_proc2glob[iproc][Proc_Num_Node_Sets[iproc]]      = i;

          proc_list_pointer[iproc][Proc_Num_Node_Sets[iproc]] =
            list_length[iproc];

          proc_list_count[iproc][Proc_Num_Node_Sets[iproc]]   =
            proc_num_ns[iproc];

          if(num_df_in_nsets[i] > 0) {
            proc_list_df_count[iproc][Proc_Num_Node_Sets[iproc]] =
              proc_num_ns[iproc];
          }
          else
            proc_list_df_count[iproc][Proc_Num_Node_Sets[iproc]] = 0;

          (Proc_Num_Node_Sets[iproc])++;

          list_length[iproc] += proc_num_ns[iproc];
        }
      }

      /* Free arrays */
      safe_free((void **) &proc_ns_pointer);

    } /* END "if (num_nodes_in_node_set[i] > 0)" */

  } /* END "for (i = 0; i < Num_Node_Set; i++)" */
  safe_free((void **) &gmap);
  safe_free((void **) &gmap_dups);

 /*--------------------------------------------------------------------------*/
 /*             WRITE PERMAMENT ARRAYS FOR NODE SET INFO                     */
 /*--------------------------------------------------------------------------*/

  for(iproc=0; iproc < Proc_Info[2]; iproc++) {

    /* Allocate Permament Arrays for node sets in one long malloc */
    if (Num_Node_Set > 0) {

      /*
       * Note that the alloc is done based on Num_Node_Set, rather than
       * Proc_Num_Node_Sets[] due to the fact that NULL entities are
       * stored on processors not having a particular node set.
       */
      Proc_NS_Ids[iproc]      =  (int *) array_alloc(__FILE__, __LINE__, 1,
                                                 (3*Num_Node_Set +
                                                  2*Proc_Num_Node_Sets[iproc] +
                                                  list_length[iproc]),
                                                     sizeof(int));

      Proc_NS_Count[iproc]     =  Proc_NS_Ids[iproc]      + Num_Node_Set;
      Proc_NS_DF_Count[iproc]  =  Proc_NS_Count[iproc]    + Num_Node_Set;
      Proc_NS_Pointers[iproc]  =  Proc_NS_DF_Count[iproc] + Num_Node_Set;
      GNode_Sets[iproc]        =  Proc_NS_Pointers[iproc] +
                                  Proc_Num_Node_Sets[iproc];
      Proc_NS_List[iproc]      =  GNode_Sets[iproc]       +
                                  Proc_Num_Node_Sets[iproc];

      Proc_NS_GNMap_List[iproc] = (int *) array_alloc(__FILE__, __LINE__, 1, list_length[iproc],
						     sizeof(int));

      if(list_length[iproc] > 0) {
        if (io_ws == sizeof(float)) {
          Proc_NS_Dist_Fact_sp[iproc] = (float *) array_alloc(__FILE__,
                                                              __LINE__, 1,
                                                          list_length[iproc],
                                                          sizeof(float));
        }
        else {
          Proc_NS_Dist_Fact_dp[iproc] = (double *) array_alloc(__FILE__,
                                                               __LINE__, 1,
                                                           list_length[iproc],
                                                           sizeof(double));
        }
      }
      else {
        if(io_ws == sizeof(float))
          Proc_NS_Dist_Fact_sp[iproc] = NULL;
        else
          Proc_NS_Dist_Fact_dp[iproc] = NULL;
      }
    }
    else {
      Proc_NS_Ids[iproc]      = NULL;
      Proc_NS_Count[iproc]    = NULL;
      Proc_NS_Pointers[iproc] = NULL;
      Proc_NS_List[iproc]     = NULL;
    }

    /*
     * Fill in the permament node set arrays which have length,
     * Proc_Num_Node_Sets, the total number of node sets defined on the
     * current processor.
     */
    Proc_NS_List_Length[iproc] = 0;
    for (i = 0; i < Proc_Num_Node_Sets[iproc]; i++) {
      GNode_Sets[iproc][i]        = ns_proc2glob[iproc][i];
      Proc_NS_Ids[iproc][i]       = Node_Set_Ids[ns_proc2glob[iproc][i]];
      Proc_NS_Count[iproc][i]     = proc_list_count[iproc][i];
      Proc_NS_Pointers[iproc][i]  = proc_list_pointer[iproc][i];
      Proc_NS_List_Length[iproc] += Proc_NS_Count[iproc][i];
      Proc_NS_DF_Count[iproc][i]  = proc_list_df_count[iproc][i];
    }

    /* Construct the concatenated node list */
    if (io_ws == sizeof(float)) {

      for (i = 0; i < Proc_Num_Node_Sets[iproc]; i++) {
        for (j = 0; j < Proc_NS_Count[iproc][i]; j++) {
          Proc_NS_List[iproc][Proc_NS_Pointers[iproc][i]+j] =
            proc_list[iproc][i][j];
	  Proc_NS_GNMap_List[iproc][Proc_NS_Pointers[iproc][i]+j] =
	    proc_gmap_list[iproc][i][j];

          if(Proc_NS_DF_Count[iproc][i] > 0) {
            Proc_NS_Dist_Fact_sp[iproc][Proc_NS_Pointers[iproc][i]+j] =
              proc_list_df_sp[iproc][i][j];
          }
        }
      }
    }
    else {

      for (i = 0; i < Proc_Num_Node_Sets[iproc]; i++) {
        for (j = 0; j < Proc_NS_Count[iproc][i]; j++) {
          Proc_NS_List[iproc][Proc_NS_Pointers[iproc][i]+j] =
            proc_list[iproc][i][j];
	  Proc_NS_GNMap_List[iproc][Proc_NS_Pointers[iproc][i]+j] =
	    proc_gmap_list[iproc][i][j];

          if(Proc_NS_DF_Count[iproc][i] > 0) {
            Proc_NS_Dist_Fact_dp[iproc][Proc_NS_Pointers[iproc][i]+j] =
              proc_list_df_dp[iproc][i][j];
          }
        }
      }
    }

  } /* End "for(iproc=0; iproc < Proc_Info[2]; iproc++)" */

  /*
   * Print Out a Table Showing the Distribution of Node Sets Across the
   * Processors
   */
#ifdef DEBUG
  if ((Debug_Flag >= 4) && (Num_Node_Set > 0)) {
    print_sync_start (Proc, Num_Proc, TRUE);
    if (Proc == 0) {
      printf("\n\tPRINT OUT OF THE DISTRIBUTION OF NODE SETS ACROSS THE "
             "PROCESSORS\n");
      printf("\n\n");
    }
    for(iproc=0; iproc < Proc_Info[2]; iproc++) {
      if (Proc_Num_Node_Sets[iproc] > 0) {
        printf("\nNode Sets Defined on Proc %d for Proc %d:\n\n", Proc,
               Proc_Ids[iproc]);
        printf("%s%s\n  ",
               " Loc_Node_# Glob_Node_# NS_ID  Ptr_Val Loc_Num_in_NS",
               " Glob_Num_in_NS");
        print_line ("-", 76);
        for (i = 0; i < Proc_Num_Node_Sets[iproc]; i++) {
          printf(" %4d %11d %9d %9d %9d %14d\n", i, GNode_Sets[iproc][i],
                 Proc_NS_Ids[iproc][i], Proc_NS_Pointers[iproc][i],
                 Proc_NS_Count[iproc][i],
                 num_nodes_in_node_set[GNode_Sets[iproc][i]]);
        }
        printf("  "); print_line ("-", 76); printf("\n");
      }
      if ((Debug_Flag >= 5) && (Proc_Num_Node_Sets[iproc] > 0)) {
        printf("\tDump of Node_Set Nodes on Proc %d for Proc %d:\n\n",
               Proc, Proc_Ids[iproc]);
        for (i = 0; i < Proc_Num_Node_Sets[iproc]; i++) {
          printf("\t\tLoc_NS_# = %5d, Glob_NS_# = %5d, NS_ID = %5d:\n\n",
                 i, GNode_Sets[iproc][i], Proc_NS_Ids[iproc][i]);
          printf("\t Nodes_In_Node_Set\t| Distribution factor\n\t");
          print_line ("-", 70);
          for (j = 0; j < Proc_NS_Count[iproc][i]; j++) {
            printf("\t %7d\t\t|     ",
                   Proc_NS_List[iproc][Proc_NS_Pointers[iproc][i]+j]);

            if(Proc_NS_DF_Count[iproc][i] > 0) {
              if (io_ws == sizeof(float))
                printf("%7.5f\n",
                     Proc_NS_Dist_Fact_sp[iproc][Proc_NS_Pointers[iproc][i]+j]);
              else
                printf("%7.5lf\n",
                     Proc_NS_Dist_Fact_dp[iproc][Proc_NS_Pointers[iproc][i]+j]);
            }
            else
              printf("none\n");

          }
          printf("\t");
          print_line ("-", 70);  printf("\n");
        }
      }
    }
    print_sync_end(Proc, Num_Proc, TRUE);
  }
#endif

  /* Free temporary arrays */
  safe_free((void **) &list_length);
  for(iproc=0; iproc < Proc_Info[2]; iproc++) {
    safe_free((void **) &(ns_proc2glob[iproc]));
    for(i=0; i < Proc_Num_Node_Sets[iproc]; i++) {
      safe_free((void **) &(proc_list[iproc][i]));
      safe_free((void **) &(proc_gmap_list[iproc][i]));

      if(Proc_NS_DF_Count[iproc][i] > 0) {
        if (io_ws == sizeof(float))
          safe_free((void **) &(proc_list_df_sp[iproc][i]));
        else
          safe_free((void **) &(proc_list_df_dp[iproc][i]));
      }
    }

    safe_free((void **) &(proc_list[iproc]));
    safe_free((void **) &(proc_gmap_list[iproc]));
    if (io_ws == sizeof(float)) safe_free((void **) &(proc_list_df_sp[iproc]));
    else                        safe_free((void **) &(proc_list_df_dp[iproc]));
  }

  safe_free((void **) &ns_proc2glob);
  safe_free((void **) &proc_list);
  safe_free((void **) &proc_gmap_list);
  if (io_ws == sizeof(float)) safe_free((void **) &proc_list_df_sp);
  else                        safe_free((void **) &proc_list_df_dp);

} /* END of routine read_node_sets () ****************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_side_sets(int exoid, int *num_elem_in_ssets,
                           int *num_df_in_ssets, int io_ws)

/*
 * Function which reads the side sets information from an EXODUS II database
 * for a given processor. It then broadcasts all information to every
 * processor. Each processor extracts and keeps only the information that it
 * needs to run its local problem.
 *
 *     Global Variables which are set by this routine
 *     -----------------------------------------------
 *
 *       Proc_Num_Side_Sets = Number of side sets on the current processor
 *
 *       GSide_Sets [Proc_Num_Side_Sets]
 *                          = Mapping between the local side set number and the
 *                            global side set number.
 *
 *       Proc_SS_Ids [Proc_Num_Side_Sets]
 *                          = Set IDs for the side sets defined on the current
 *                            processor
 *
 *       Proc_SS_Elem_Count [Proc_Num_Side_Sets]
 *                          = Number of elements in the side set for the
 *                            current processor
 *
 *       Proc_SS_Side_Count [Proc_Num_Side_Sets]
 *                          = Number of nodes in the node set list for the side
 *                            set for the current processor.
 *
 *       Proc_SS_Elem_Pointers [Proc_Num_Side_Sets]
 *                          = Side sets pointer record for elements for
 *                            the side sets in a given processor
 *
 *       Proc_SS_Elem_List_Length
 *                          = Length of the element list length for side sets
 *                            on the current processor
 *
 *       Proc_SS_Elem_List [Proc_SS_Elem_List_Length]
 *                          = Concatenated vector of elements that comprise
 *                            the side set definitions on the current
 *                            processor.
 *
 *       Proc_SS_Side_List [Proc_SS_Elem_List_Length]
 *                          = Concatenated vector of sides that comprise
 *                            the side set definitions on the current
 *                            processor.
 */

{

  /* Local variables */

  char *yo = "read_side_sets";

  int i, j, k, iproc, ipos_elem;
  int iss_size, num_messages, num_left_over, num_elem_per_message;
  int imess, istart_ss, ntotal_s, df_loc, idx_begin, idx_end;
  int *ss_elem_cntr, *first_message;

  int *ss_on_proc;            /* Flag to indicate that the ss is on the proc */
  int *proc_num_sides_ss;     /* Number of sides in the current side set that
                                 exist for the current processor             */
  int *ss_elem_list;          /* side-set element list                       */
  int *ss_side_list;          /* side-set node list                          */
  int *ss_df_ptr;             /* pointer into the dist. factor list          */
  int *proc_es_pointer;       /* pointers into the global element side-set list
                                 for elements which are common to the processor
                                 and the side-set                            */
  void *ss_dist_fact;         /* side-set distribution factors                     */
  int *proc_ss_pointer;       /* pointers into the global node side-set list
                                 for nodes which are common to the processor
                                 and the side-set                            */
  int ***proc_elem_list;      /*  Vector of pointers to side-set element lists.
                                  There is one pointer for each side-set     */
  int ***proc_side_list;      /*  Vector of pointers to side-set side lists.
                                  There is one pointer for each side-set     */
  int ***proc_gmap_list;      /*  Vector of pointers to side-set global sideset element map lists.
                                  There is one pointer for each side-set     */
  int ***proc_df_ptr;        /*  Vector of pointers to side-set df pointers
				 There is one pointer for each side-set     */
  int ***proc_df_indx;
  int ***proc_df_indx_cnt;
  int  **proc_df_map;        /* maps the sorted list of df pointer back to
                                their original location for each proc        */
  float  ***proc_ss_df_sp;   /* Array of the values of the dist. factors.    */
  double ***proc_ss_df_dp;   /* Array of the values of the dist. factors.    */

  int **proc_ss_df_cnt;      /* Array of counts of the dist. factors.        */

  int **ss_proc2glob;         /* mapping from the processor side-set numbers
                                 to the global side-set numbers              */
  int **proc_elem_list_ptr;   /* pointers into the concatenated element list
                                 which locate the start of the element lists */
  int **proc_elem_list_cnt;   /* the number of elements in each processor's
                                 element set                                 */
  int *elem_list_length;      /* length of the element side-set list for the
                                 current processor                           */

  int *gmap;

  /* generic pointers for single/double precision */
  void   *ptr;
  float  *ptr_sp;
  double *ptr_dp;

  int ilast, ilast_side, ss_indx;
  int loc_elem, loc_side, *ss_num, *ntotal, indx;

  extern void find_message_info(int , int , int *, int *, int *);
  extern int  bin_search2(int value, int num, int List[]);

  /*************************** execution begins ******************************/

  /* check to see if the io_ws is smaller than the machine precision */
  if (io_ws < sizeof(float)) io_ws = sizeof(float);

  proc_elem_list   = (int ***)array_alloc(__FILE__, __LINE__, 1,
                                          5*Proc_Info[2], sizeof(int **));
  proc_side_list   = proc_elem_list + Proc_Info[2];
  proc_df_ptr      = proc_side_list + Proc_Info[2];
  proc_df_indx     = proc_df_ptr    + Proc_Info[2];
  proc_df_indx_cnt = proc_df_indx   + Proc_Info[2];

  proc_gmap_list   = (int ***)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2], sizeof(int **));

  if (io_ws == sizeof(float)) {
    proc_ss_df_sp = (float ***)array_alloc(__FILE__, __LINE__, 1,
                                           3*Proc_Info[2], sizeof(float **));
  }
  else {
    proc_ss_df_dp = (double ***)array_alloc(__FILE__, __LINE__, 1,
                                            3*Proc_Info[2], sizeof(double **));
  }

  ss_elem_cntr      = (int *)array_alloc(__FILE__, __LINE__, 1, 7*Proc_Info[2],
                                         sizeof(int));
  ss_on_proc        = ss_elem_cntr      + Proc_Info[2];
  elem_list_length  = ss_on_proc        + Proc_Info[2];
  proc_num_sides_ss = elem_list_length  + Proc_Info[2];
  first_message     = proc_num_sides_ss + Proc_Info[2];
  ntotal            = first_message     + Proc_Info[2];
  ss_num            = ntotal            + Proc_Info[2];

  /* Allocate temporary arrays */
  ss_proc2glob       = (int **)array_alloc(__FILE__, __LINE__, 1,
                                           4*Proc_Info[2], sizeof(int *));
  proc_elem_list_ptr = ss_proc2glob       + Proc_Info[2];
  proc_elem_list_cnt = proc_elem_list_ptr + Proc_Info[2];
  proc_ss_df_cnt     = proc_elem_list_cnt + Proc_Info[2];

  for(iproc=0; iproc < Proc_Info[2]; iproc++) {
    proc_elem_list[iproc]   = (int **) array_alloc(__FILE__, __LINE__, 1,
                                                   5*Num_Side_Set,
                                                   sizeof(int *));
    proc_side_list[iproc]   = proc_elem_list[iproc] + Num_Side_Set;
    proc_df_ptr[iproc]      = proc_side_list[iproc] + Num_Side_Set;
    proc_df_indx[iproc]     = proc_df_ptr[iproc]    + Num_Side_Set;
    proc_df_indx_cnt[iproc] = proc_df_indx[iproc]   + Num_Side_Set;
    proc_gmap_list[iproc]   = (int **) array_alloc(__FILE__, __LINE__, 1, Num_Side_Set,
                                                   sizeof(int *));

    if (io_ws == sizeof(float))
      proc_ss_df_sp[iproc]       = (float **)array_alloc(__FILE__, __LINE__, 1,
                                                         Num_Side_Set,
                                                         sizeof(float *));
    else
      proc_ss_df_dp[iproc]       = (double **)array_alloc(__FILE__, __LINE__,
                                                          1, Num_Side_Set,
                                                          sizeof(double *));

    ss_proc2glob[iproc]       = (int *)  array_alloc(__FILE__, __LINE__, 1,
                                                     4*Num_Side_Set,
                                                     sizeof(int));
    proc_elem_list_ptr[iproc] = ss_proc2glob[iproc]       + Num_Side_Set;
    proc_elem_list_cnt[iproc] = proc_elem_list_ptr[iproc] + Num_Side_Set;
    proc_ss_df_cnt[iproc]     = proc_elem_list_cnt[iproc] + Num_Side_Set;

    for(i=0; i < 2*Num_Side_Set; i++)
      proc_elem_list[iproc][i] = NULL;

    for(i=0; i < Num_Side_Set; i++) {
      proc_ss_df_cnt[iproc][i] = 0;
    }

  }

  /* Initialize */
  for(iproc=0; iproc < Proc_Info[2]; iproc++) {
    Proc_Num_Side_Sets[iproc] = 0;
    elem_list_length[iproc] = 0;
  }

  /* Create a map from 'global element id' to processor */
  {
    gmap = (int *) array_alloc(__FILE__, __LINE__, 1, Num_Elem, sizeof(int));

    /* Initialize to -1 */
    for (i=0; i < Num_Elem; i++)
      gmap[i] = -1;

    /* Fill it in.... If run in parallel, only fills with the
     * elements being processed on this processor, not entire model
     */
    for (iproc=0; iproc < Proc_Info[2]; iproc++) {
      int size = Num_Internal_Elems[iproc]+Num_Border_Elems[iproc];
      for (i=0; i < size; i++) {
	int ind = GElems[iproc][i];
	gmap[ind] = iproc;
      }
    }
  }

  /*-------------------------------------------------------------------------*/
  /*                    LOOP OVER THE SIDE SETS                              */
  /*-------------------------------------------------------------------------*/
  for(i=0; i < Num_Side_Set; i++) {

    if (num_elem_in_ssets[i] > 0) {

      ilast_side = 0;

      for(iproc=0; iproc < Proc_Info[2]; iproc++) {
        ss_on_proc[iproc] = FALSE;
        proc_num_sides_ss[iproc] = 0;
        first_message[iproc]     = TRUE;
        ss_elem_cntr[iproc] = 0;
      }

      /* One element ID + One side ID */
      iss_size = 3*sizeof(int);

      find_message_info(iss_size, num_elem_in_ssets[i], &num_elem_per_message,
                        &num_messages, &num_left_over);

      if(Debug_Flag >= 2) {
        if(Proc == 0) {
          printf("Message summary for Side Set number %d, with an ID of %d:\n",
                 i, Side_Set_Ids[i]);
          printf("\tNumber of messages needed for element and side list = %d\n",
                 num_messages);
          printf("\tNumber of element and side IDs per message = %d\n",
                 num_elem_per_message);
          printf("\tLength of each message = %d\n",
                 iss_size*num_elem_per_message);
        }
      }

      /* Allocate temporary storage for the current message for the
       * current side set on each proc
       */
      proc_es_pointer = (int *) array_alloc(__FILE__, __LINE__, 1,
                                            5*num_elem_per_message,
                                            sizeof(int));
      ss_elem_list    = proc_es_pointer + num_elem_per_message;
      ss_side_list    = ss_elem_list    + num_elem_per_message;
      ss_df_ptr       = ss_side_list    + num_elem_per_message;
      proc_ss_pointer = ss_df_ptr       + num_elem_per_message;

      /*---------------------------------------------------------------------*/
      /*                  LOOP OVER THE MESSAGES                             */
      /*---------------------------------------------------------------------*/
      for(imess=0; imess < num_messages; imess++) {

        if(Debug_Flag >= 2 && Proc == 0)
          printf("\tside set message: %d of %d\n", imess+1, num_messages);

        istart_ss = imess*num_elem_per_message;

        if (num_left_over != 0 && imess == num_messages - 1) {
          num_elem_per_message = num_left_over;

          /* Pointers need to be adjusted */
          ss_side_list = ss_elem_list + num_elem_per_message;
          ss_df_ptr    = ss_side_list + num_elem_per_message;
        }

        /* Read in the part of the side set that will fit in the message. */

        if (Proc == 0) {
          check_exodus_error(ne_get_n_side_set(exoid, Side_Set_Ids[i],
                                               (istart_ss + 1),
                                               num_elem_per_message,
                                               ss_elem_list, ss_side_list),
                             "ne_get_n_side_set");

          /* Fill in the distribution factor pointer vector */
          if(imess == 0) {
            ilast = 0;
            ilast_side = 0;
          }

          ss_df_ptr[0] = ilast + ilast_side;

          for(j=1; j < num_elem_per_message; j++) {

	    /*
	     * Kluge for "special" sidesest which has a single distribution
	     * factor for each element instead of 1 per face node..
	     */
	    if (num_elem_in_ssets[i] == num_df_in_ssets[i]) {
	      ilast_side = 1;
	    } else {

	      ilast_side = elem_info(NN_SIDE,
				     GM_Elem_Types[(ss_elem_list[j-1])-1],
				     ss_side_list[j-1]);
	      
	      /*
	       * kludge for HEXSHELL's
	       * where distribution factors are concerned, only count
	       * 4 nodes on the 6 node faces (they are just like HEX's)
	       */
	      if (GM_Elem_Types[(ss_elem_list[j-1])-1] == HEXSHELL)
		ilast_side = 4;
	    }
	    
            ss_df_ptr[j] = ss_df_ptr[j-1] + ilast_side;
	    
            ilast = ss_df_ptr[j];
          }
        }

        /* Broadcast the side set information to all of the processors */
	if (Num_Proc > 1) {
	  brdcst_maxlen(Proc, Num_Proc, (char *) ss_elem_list,
			3*num_elem_per_message*sizeof(int), 0);

	  psync(Proc, Num_Proc);
	}

        /* Renumber elements to start at '0' instead of '1' */
        for (j = 0; j < num_elem_per_message; ss_elem_list[j++]--);

        for(iproc=0; iproc < Proc_Info[2]; iproc++) {

          /*
           * Find the intersection between the elements in the side set and the
           * elements owned by the current processor.  
           */

	  /* Note that an element may be in ss_elem_list multiple times */
          ipos_elem = 0;
          for (j=0; j < num_elem_per_message; j++) {
            if (gmap[ss_elem_list[j]] == iproc)
              proc_es_pointer[ipos_elem++] = j;
          }

          if(ipos_elem > 0) {
            proc_num_sides_ss[iproc] += ipos_elem;
            ss_on_proc[iproc] = TRUE;

	    /*
	     * Store the information for the current side set in temporary arrays.
	     *
	     * This part of the routine only gets executed if the side set is
	     * defined on the current processor.
	     */
	    
	    
            /* Allocate and store element information in a temporary vector */
            if(first_message[iproc]) {

              proc_elem_list[iproc][Proc_Num_Side_Sets[iproc]] =
                (int *) array_alloc(__FILE__, __LINE__, 1,   proc_num_sides_ss[iproc], sizeof(int));

              proc_side_list[iproc][Proc_Num_Side_Sets[iproc]] =
                (int *) array_alloc(__FILE__, __LINE__, 1,   proc_num_sides_ss[iproc], sizeof(int));

              proc_gmap_list[iproc][Proc_Num_Side_Sets[iproc]] =
                (int *) array_alloc(__FILE__, __LINE__, 1,   proc_num_sides_ss[iproc], sizeof(int));

              proc_df_ptr[iproc][Proc_Num_Side_Sets[iproc]] =
                (int *) array_alloc(__FILE__, __LINE__, 1,   proc_num_sides_ss[iproc], sizeof(int));

              /* This could probably be allocated in one-shot later */
              proc_df_indx[iproc][Proc_Num_Side_Sets[iproc]] = 
                (int *) array_alloc(__FILE__, __LINE__, 1, 2*proc_num_sides_ss[iproc], sizeof(int));

              proc_df_indx_cnt[iproc][Proc_Num_Side_Sets[iproc]] =
                proc_df_indx[iproc][Proc_Num_Side_Sets[iproc]] +
                proc_num_sides_ss[iproc];

              first_message[iproc] = FALSE;
            }
            else {

              if((proc_elem_list[iproc][Proc_Num_Side_Sets[iproc]] = (int *)
                  realloc(proc_elem_list[iproc][Proc_Num_Side_Sets[iproc]],
                          (proc_num_sides_ss[iproc]*sizeof(int)))) == NULL ) {
                fprintf(stderr, "ERROR: unable to realloc proc_elem_list in "
			"%s\n", yo);
                exit(1);
              }

              if((proc_side_list[iproc][Proc_Num_Side_Sets[iproc]] = (int *)
                  realloc(proc_side_list[iproc][Proc_Num_Side_Sets[iproc]],
                          (proc_num_sides_ss[iproc]*sizeof(int)))) == NULL ) {
                fprintf(stderr, "ERROR: unable to realloc proc_side_list in "
			"%s\n", yo);
                exit(1);
              }

              if((proc_gmap_list[iproc][Proc_Num_Side_Sets[iproc]] = (int *)
                  realloc(proc_gmap_list[iproc][Proc_Num_Side_Sets[iproc]],
                          (proc_num_sides_ss[iproc]*sizeof(int)))) == NULL ) {
                fprintf(stderr, "ERROR: unable to realloc proc_gmap_list in "
			"%s\n", yo);
                exit(1);
              }

              if((proc_df_ptr[iproc][Proc_Num_Side_Sets[iproc]] = (int *)
                  realloc(proc_df_ptr[iproc][Proc_Num_Side_Sets[iproc]],
                          (proc_num_sides_ss[iproc]*sizeof(int)))) == NULL ) {
                fprintf(stderr, "ERROR: unable to realloc proc_df_ptr in "
			"%s\n", yo);
                exit(1);
              }

              /* This could probably be allocated in one-shot later */
              if((proc_df_indx[iproc][Proc_Num_Side_Sets[iproc]] = (int *)
                  realloc(proc_df_indx[iproc][Proc_Num_Side_Sets[iproc]],
                          2*(proc_num_sides_ss[iproc]*sizeof(int)))) == NULL ) {
                fprintf(stderr, "ERROR: unable to realloc proc_df_indx in "
			"%s\n", yo);
                exit(1);
              }

              proc_df_indx_cnt[iproc][Proc_Num_Side_Sets[iproc]] =
                proc_df_indx[iproc][Proc_Num_Side_Sets[iproc]] +
                proc_num_sides_ss[iproc];

            }

	    /* Transfer the element, side number, and df pointers  into "permanent" storage... */
            for(j=0; j < ipos_elem; j++) {
	      proc_elem_list[iproc][Proc_Num_Side_Sets[iproc]][ss_elem_cntr[iproc]] =
		ss_elem_list[proc_es_pointer[j]];
	      
	      proc_side_list[iproc][Proc_Num_Side_Sets[iproc]][ss_elem_cntr[iproc]] =
		ss_side_list[proc_es_pointer[j]];
	      
	      proc_gmap_list[iproc][Proc_Num_Side_Sets[iproc]][ss_elem_cntr[iproc]] =
		proc_es_pointer[j]+istart_ss;
	      
	      proc_df_ptr[iproc][Proc_Num_Side_Sets[iproc]][(ss_elem_cntr[iproc])++] =
		ss_df_ptr[proc_es_pointer[j]];
	      
            }
          } /* End "if(ipos_elem > 0)" */
        } /* End "for(iproc=0; iproc < Proc_Info[2]; iproc++)" */
      } /* End "for(imess=0; imess < num_messages; imess++)" */
      
      /* Entire sideset has been read at this point */
      /*
       * If any part of this side-set is on the processor, update the various
       * pointers, lengths, etc.
       */

      for(iproc=0; iproc < Proc_Info[2]; iproc++) {
        if (ss_on_proc[iproc]) {
          ss_proc2glob[iproc][Proc_Num_Side_Sets[iproc]] = i;

          proc_elem_list_ptr[iproc][Proc_Num_Side_Sets[iproc]] =
            elem_list_length[iproc];

          proc_elem_list_cnt[iproc][Proc_Num_Side_Sets[iproc]] =
            proc_num_sides_ss[iproc];

          (Proc_Num_Side_Sets[iproc])++;
          elem_list_length[iproc] += proc_num_sides_ss[iproc];
        }
      }

      /* Process any distribution factors in the side set */
      if(num_df_in_ssets[i] > 0) {

        iss_size = sizeof(float);
        find_message_info(iss_size, num_df_in_ssets[i], &num_elem_per_message,
			  &num_messages, &num_left_over);

        if(Debug_Flag >= 4) {
          if(Proc == 0) {
            printf("Message summary for Side Set number %d, with ID of %d:\n",
                   i, Side_Set_Ids[i]);
            printf("\tNumber of messages needed for distribution "
                   "factors = %d\n", num_messages);
            printf("\tNumber of dist. factors in each message = %d\n",
                   num_elem_per_message);
            printf("\tLength of each message = %d\n",
                   (int)(num_elem_per_message * sizeof(float)));
          }
        }

        ss_dist_fact = array_alloc(__FILE__, __LINE__, 1, num_elem_per_message,
                                   io_ws);

        /* set up an array to help speed up the searches below */
        proc_df_map = (int **) array_alloc(__FILE__, __LINE__, 1,
                                           Proc_Info[2], sizeof(int *));

        /* Loop over the messages needed for the distribution factors */
        for(imess=0; imess < num_messages; imess++) {

          istart_ss = imess*num_elem_per_message;
          if(num_left_over != 0 && imess == (num_messages-1))
            num_elem_per_message = num_left_over;

          /* Read in the part of the side set df's that will fit in the msg. */
          if(Proc == 0) {
            check_exodus_error(ne_get_n_side_set_df(exoid, Side_Set_Ids[i],
                                                    (istart_ss + 1),
                                                    num_elem_per_message,
                                                    ss_dist_fact),
                               "ne_get_n_side_set_df");
          }

          /* Broadcast the side set df information to all processors */
	  if (Num_Proc > 1) {
	    brdcst_maxlen(Proc, Num_Proc, (char *)ss_dist_fact,
			  num_elem_per_message*io_ws, 0);
	  }

          /*
           * At this point a processor has the list of global element IDs
           * that it owns that are contained in the side set. It also has
           * a vector of pointers into the distribution factor vector telling
           * where the distribution factors for a particular element in
           * the side set begin. Now need to find out if the information
           * just broadcast is needed by this processor.
           */

          /*
           * First figure out how many df's there are for this side set
           * on this processor.
           */
          for(iproc=0; iproc < Proc_Info[2]; iproc++) {

            if(imess == 0)
              ntotal[iproc] = 0;

            if(imess == 0 && ss_on_proc[iproc]) {

              ss_num[iproc] = Proc_Num_Side_Sets[iproc];
              ntotal[iproc] = 0;

              indx = 0;

              proc_df_map[iproc] = (int *) array_alloc(__FILE__, __LINE__, 1,
						       proc_elem_list_cnt[iproc][ss_num[iproc]-1],
                                                       sizeof(int));

              for(j=0; j < proc_elem_list_cnt[iproc][ss_num[iproc]-1]; j++) {

                loc_elem = proc_elem_list[iproc][ss_num[iproc]-1][j];
                loc_side = proc_side_list[iproc][ss_num[iproc]-1][j];

		/*
		 * Kluge for "special" sidesest which has a single distribution
		 * factor for each element instead of 1 per face node..
		 */
		if (num_elem_in_ssets[i] == num_df_in_ssets[i]) {
		  ntotal_s = 1;
		} else {

		  ntotal_s = elem_info(NN_SIDE, GM_Elem_Types[loc_elem],
				       loc_side);
		  
		  /*
		   * kludge for HEXSHELL's
		   * where distribution factors are concerned, only count
		   * 4 nodes on the 6 node faces (they are just like HEX's)
		   */
		  if (GM_Elem_Types[loc_elem] == HEXSHELL)
		    ntotal_s = 4;
		}
		
                ntotal[iproc] += ntotal_s;

                proc_df_indx[iproc][ss_num[iproc]-1][j] = indx;
                proc_df_indx_cnt[iproc][ss_num[iproc]-1][j] = 0;

                indx += ntotal_s;

                /* and set up map for sort in a minute */
                proc_df_map[iproc][j] = j;

              }

              /*
               * now sort the proc_dr_ptr array so that it is monotonic,
               * and can be searched more easily
               */
              sortN_int_int(proc_elem_list_cnt[iproc][ss_num[iproc]-1],
                            proc_df_ptr[iproc][ss_num[iproc]-1], 1,
                            proc_df_map[iproc]);

            }

            /* Allocate memory for the dfs */
            if(ntotal[iproc] > 0) {

              if(imess == 0) {
                ptr = array_alloc(__FILE__, __LINE__, 1, ntotal[iproc], io_ws);
                if (io_ws == sizeof(float))
                  proc_ss_df_sp[iproc][ss_num[iproc]-1] = (float *) ptr;
                else
                  proc_ss_df_dp[iproc][ss_num[iproc]-1] = (double *) ptr;
              }

              if (io_ws == sizeof(float))
                ptr_sp = (float *) ss_dist_fact;
              else
                ptr_dp = (double *) ss_dist_fact;

              k = 0;
              ss_indx = ss_num[iproc]-1;

              for(j=0; j < num_elem_per_message; j++) {

                indx = istart_ss + j;

                /*
                 * Find out if this index fits in the list of this
                 * processors indices.
                 */

                /* 'istart_ss' <= indx <= 'istart_ss+num_elem_per_message' */
                /* 'indx' increments by 1 each time through loop */

                /* 'proc_df_ptr[iproc][ss_indx]' is the list
                   and it is sorted....
                   ... So, don't need a binary search, can just increment
                   as we proceed through the loops...
                */
                if (indx >= proc_df_ptr[iproc][ss_indx][0]) {

                  while (k < proc_elem_list_cnt[iproc][ss_indx]-1 &&
                         indx >= proc_df_ptr[iproc][ss_indx][k+1])
                    k++;

                  /* now get the proper location from the map */
                  ipos_elem = proc_df_map[iproc][k];

                  idx_begin = proc_df_indx[iproc][ss_indx][ipos_elem];
                  if ((ipos_elem+1) < proc_elem_list_cnt[iproc][ss_indx])
                    idx_end = proc_df_indx[iproc][ss_indx][ipos_elem+1];
                  else
                    idx_end = ntotal[iproc];

                  ntotal_s = idx_end - idx_begin;

                  /*
                   * check to see if indx really fits in this spot
                   *
                   * need to make sure that indx > the k-th value since
                   * bin_search_min will return 0 for anything less than
                   * the first entry in the list
                   */
                  if (indx >= proc_df_ptr[iproc][ss_indx][k] &&
                      indx < proc_df_ptr[iproc][ss_indx][k] + ntotal_s) {

                    df_loc = proc_df_indx[iproc][ss_indx][ipos_elem] +
                      proc_df_indx_cnt[iproc][ss_indx][ipos_elem];

                    if (io_ws == sizeof(float))
                      proc_ss_df_sp[iproc][ss_indx][df_loc] =
                        ptr_sp[indx-istart_ss];
                    else
                      proc_ss_df_dp[iproc][ss_indx][df_loc] =
                        ptr_dp[indx-istart_ss];

                    proc_ss_df_cnt[iproc][ss_indx]++;

                    proc_df_indx_cnt[iproc][ss_indx][ipos_elem]++;
                  }
                }

              } /* End "for(j=0; j < num_elem_per_message; j++)" */

            } /* End "if(ntotal[iproc] > 0)" */

          } /* End "for(iproc=0; iproc < Proc_Info[2]; iproc++)" */

        } /* End "for(imess=0; imess < num_messages; imess++)" */

        /* Free temporary memory */
        safe_free((void **) &ss_dist_fact);
        for (iproc=0; iproc < Proc_Info[2]; iproc++)
          if (ss_on_proc[iproc])
            safe_free((void **) &(proc_df_map[iproc]));
        safe_free((void **) &proc_df_map);

      } /* End "if(num_df_in_ssets[i] > 0)" */

      /* Free some of the temporary arrays */
      safe_free((void **) &proc_es_pointer);

    } /* END "if (num_elem_in_ssets[i] > 0)" */

  } /* END "for (i = 0; i < Num_Side_Set; i++)" */

  safe_free((void **) &gmap);
  /*-------------------------------------------------------------------------*/
  /*---------------- Store Structures back into a packed form ---------------*/
  /*-------------------------------------------------------------------------*/

  /*
   * Allocate storage for permament side set integer info vectors using one
   * long malloc statement.
   */
  for(iproc=0; iproc < Proc_Info[2]; iproc++) {
    if (Num_Side_Set > 0) {

      /*
       * Note that the alloc is done based on Num_Side_Set, rather than
       * Proc_Num_Side_Sets[] due to the fact that NULL entities are
       * stored on processors not having a particular side set.
       */
      Proc_SS_Ids[iproc] =  (int *) array_alloc(__FILE__, __LINE__, 1,
                                                (3*Num_Side_Set  +
						 3*Proc_Num_Side_Sets[iproc] +
						 2*elem_list_length[iproc] + 1),
                                                sizeof(int));
      Proc_SS_Elem_Count[iproc]    = Proc_SS_Ids[iproc]           + Num_Side_Set;
      Proc_SS_Elem_Pointers[iproc] = Proc_SS_Elem_Count[iproc]    + Num_Side_Set;
      Proc_SS_DF_Count[iproc]      = Proc_SS_Elem_Pointers[iproc] + Proc_Num_Side_Sets[iproc];
      GSide_Sets[iproc]            = Proc_SS_DF_Count[iproc]      + Num_Side_Set;
      Proc_SS_Elem_List[iproc]     = GSide_Sets[iproc]            + Proc_Num_Side_Sets[iproc];
      Proc_SS_Side_List[iproc]     = Proc_SS_Elem_List[iproc]     + elem_list_length[iproc];
      Proc_SS_DF_Pointers[iproc]   = Proc_SS_Side_List[iproc]     + elem_list_length[iproc];

      Proc_SS_GEMap_List[iproc] = (int *) array_alloc(__FILE__, __LINE__, 1, elem_list_length[iproc],
						     sizeof(int));
    }

    /* Construct the side sets global array information */
    Proc_SS_Elem_List_Length[iproc] = 0;
    for (i = 0; i < Proc_Num_Side_Sets[iproc]; i++) {
      GSide_Sets[iproc][i]             = ss_proc2glob[iproc][i];
      Proc_SS_Ids[iproc][i]            = Side_Set_Ids[ss_proc2glob[iproc][i]];
      Proc_SS_Elem_Count[iproc][i]     = proc_elem_list_cnt[iproc][i];
      Proc_SS_Elem_Pointers[iproc][i]  = proc_elem_list_ptr[iproc][i];
      Proc_SS_Elem_List_Length[iproc] += Proc_SS_Elem_Count[iproc][i];
      Proc_SS_DF_Count[iproc][i]       = proc_ss_df_cnt[iproc][i];

      if(i == 0)
        Proc_SS_DF_Pointers[iproc][i]  = 0;
      else
        Proc_SS_DF_Pointers[iproc][i]  = Proc_SS_DF_Pointers[iproc][i-1] +
          proc_ss_df_cnt[iproc][i-1];
    }

    if(Proc_Num_Side_Sets[iproc] > 0) {
      Proc_SS_DF_Pointers[iproc][i] = Proc_SS_DF_Pointers[iproc][i-1] +
        proc_ss_df_cnt[iproc][i-1];

      /* Allocate memory for the DF vector */
      if(Proc_SS_DF_Pointers[iproc][i] > 0) {
        ptr = array_alloc(__FILE__, __LINE__, 1, Proc_SS_DF_Pointers[iproc][i],
                          io_ws);
        if (io_ws == sizeof(float))
          Proc_SS_Dist_Fact_sp[iproc] = (float *) ptr;
        else
          Proc_SS_Dist_Fact_dp[iproc] = (double *) ptr;
      }
      else {
        if(io_ws == sizeof(float))
          Proc_SS_Dist_Fact_sp[iproc] = NULL;
        else
          Proc_SS_Dist_Fact_dp[iproc] = NULL;
      }
    }

    /* Construct the concatenated element and side list */
    for (i = 0; i < Proc_Num_Side_Sets[iproc]; i++) {
      for (j = 0; j < Proc_SS_Elem_Count[iproc][i]; j++) {
        Proc_SS_Elem_List[iproc][Proc_SS_Elem_Pointers[iproc][i]+j] =
          proc_elem_list[iproc][i][j];
        Proc_SS_Side_List[iproc][Proc_SS_Elem_Pointers[iproc][i]+j] =
          proc_side_list[iproc][i][j];
	Proc_SS_GEMap_List[iproc][Proc_SS_Elem_Pointers[iproc][i]+j] =
	  proc_gmap_list[iproc][i][j];
      }
      if (io_ws == sizeof(float)) {
        for (j = 0; j < Proc_SS_DF_Count[iproc][i]; j++) {
          Proc_SS_Dist_Fact_sp[iproc][Proc_SS_DF_Pointers[iproc][i]+j] =
            proc_ss_df_sp[iproc][i][j];
        }
      }
      else {
        for (j = 0; j < Proc_SS_DF_Count[iproc][i]; j++) {
          Proc_SS_Dist_Fact_dp[iproc][Proc_SS_DF_Pointers[iproc][i]+j] =
            proc_ss_df_dp[iproc][i][j];
        }
      }
    }

  } /* End "for(iproc)" */

  /*
   * Print Out a Table Showing the Distribution of Side Sets Across the
   * Processors
   */
#ifdef DEBUG
  if (Debug_Flag > 3 && (Num_Side_Set > 0)) {

    /*
     * Sync the processors, printing out a long line, because this will be a
     * long print-out.
     */
    print_sync_start (Proc, Num_Proc, TRUE);
    if (Proc == 0) {
      printf("\n  PRINT OUT OF THE DISTRIBUTION OF SIDE SETS ACROSS THE "
             "PROCESSORS\n\n\n");
    }
    for(iproc=0; iproc < Proc_Info[2]; iproc++) {
      if (Proc_Num_Side_Sets[iproc] > 0) {
        printf("\nSide Sets Defined on Proc %d for Proc %d:\n\n", Proc,
               Proc_Ids[iproc]);
        printf("%s%s\n ",
               " Loc_SS# Glb_SS# SS_ID Elem_Ptr ",
               "Loc_#Sides Loc_#DF Glob_#Sides");
        print_line ("-", 76);
        for (i = 0; i < Proc_Num_Side_Sets[iproc]; i++) {
          printf("%6d %6d %7d %6d %9d %9d %9d\n", i,
                 GSide_Sets[iproc][i], Proc_SS_Ids[iproc][i],
                 Proc_SS_Elem_Pointers[iproc][i],
                 Proc_SS_Elem_Count[iproc][i],
                 Proc_SS_DF_Count[iproc][i],
                 num_elem_in_ssets[GSide_Sets[iproc][i]]);
        }
        printf(" "); print_line ("-", 76); printf("\n");

        /*
         * If Debug_Flag is large enough, dump out a complete listing of the
         * definition of the side set on the current processor
         */
        if ((Debug_Flag > 4) && (Proc_Num_Side_Sets[iproc] > 0)) {
          printf("  Dump of Side_Set Sides on Proc %d for Proc %d:\n\n",
                 Proc, Proc_Ids[iproc]);
          indx = 0;
          for (i = 0; i < Proc_Num_Side_Sets[iproc]; i++) {
            printf("\tLoc_SS_# = %5d, Glob_SS_# = %5d, SS_ID = %5d:\n\n",
                   i, GSide_Sets[iproc][i], Proc_SS_Ids[iproc][i]);
            printf("\t\tLoc_Elem_Num |    Glob_Elem_Num  |  Side Ids  | "
                   "Dist. Factors\n\t\t");
            print_line ("-", 60);
            for (j = 0; j < Proc_SS_Elem_Count[iproc][i]; j++) {

              elem_num = Proc_SS_Elem_List[iproc][
						  Proc_SS_Elem_Pointers[iproc][i]+j];

              /* Find the local element ID */
              for(loc_elem=0; loc_elem < Num_Internal_Elems[iproc]+
                    Num_Border_Elems[iproc];
                  loc_elem++) {
                if(GElems[iproc][loc_elem] == elem_num)
                  break;
              }

              loc_side = Proc_SS_Side_List[iproc][
						  Proc_SS_Elem_Pointers[iproc][i]+j];

              printf("\t\t %6d      | ", loc_elem);

              printf("    %6d        | ", elem_num+1);

              printf("%6d     | ",
                     Proc_SS_Side_List[iproc][
					      Proc_SS_Elem_Pointers[iproc][i]+j]);

              if(Proc_SS_DF_Count[iproc][i] > 0) {
		if (num_elem_in_ssets[i] == num_df_in_ssets[i]) {
		  imess = 1;
		} else {
		  imess = elem_info(NN_SIDE, Elem_Type[iproc][loc_elem],
				    loc_side);
		  
		  /*
		   * kludge for HEXSHELL's
		   * where distribution factors are concerned, only count
		   * 4 nodes on the 6 node faces (they are just like HEX's)
		   */
		  if (Elem_Type[iproc][loc_elem] == HEXSHELL)
		    imess = 4;
		}
                printf(" ");
                for(k=0; k < imess-1; k++) {
                  if (io_ws == sizeof(float))
                    printf("%5.2f ", Proc_SS_Dist_Fact_sp[iproc][indx++]);
                  else
                    printf("%5.2lf ", Proc_SS_Dist_Fact_dp[iproc][indx++]);
                }

                if (io_ws == sizeof(float))
                  printf("%5.2f\n", Proc_SS_Dist_Fact_sp[iproc][indx++]);
                else
                  printf("%5.2lf\n", Proc_SS_Dist_Fact_dp[iproc][indx++]);
              }
              else
                printf("none\n");

            }

            printf("\t\t"); print_line ("-", 60); printf("\n");

          }
          printf("\n");
        }
      }
    }
    print_sync_end(Proc, Num_Proc, TRUE);
  }
#endif

  /* Free temporary arrays */
  for(iproc=0; iproc < Proc_Info[2]; iproc++) {
    for (i = 0 ; i < Proc_Num_Side_Sets[iproc] ; i++) {
      safe_free((void **) &(proc_elem_list[iproc][i]));
      safe_free((void **) &(proc_side_list[iproc][i]));
      safe_free((void **) &(proc_gmap_list[iproc][i]));

      if(Proc_SS_DF_Count[iproc][i] > 0) {

        if (io_ws == sizeof(float))
          safe_free((void **) &(proc_ss_df_sp[iproc][i]));
        else
          safe_free((void **) &(proc_ss_df_dp[iproc][i]));
      }

      safe_free((void **) &(proc_df_indx[iproc][i]));
      safe_free((void **) &(proc_df_ptr[iproc][i]));
    }
    safe_free((void **) &(proc_elem_list[iproc]));
    safe_free((void **) &(proc_gmap_list[iproc]));
    if (io_ws == sizeof(float)) safe_free((void **) &(proc_ss_df_sp[iproc]));
    else                        safe_free((void **) &(proc_ss_df_dp[iproc]));
  }

  for(iproc=0; iproc < Proc_Info[2]; iproc++)
    safe_free((void **) &(ss_proc2glob[iproc]));

  safe_free((void **) &ss_proc2glob);
  safe_free((void **) &proc_elem_list);
  safe_free((void **) &proc_gmap_list);
  safe_free((void **) &ss_elem_cntr);
  if (io_ws == sizeof(float)) safe_free((void **) &proc_ss_df_sp);
  else                        safe_free((void **) &proc_ss_df_dp);

  /* Free list of global element types */
  safe_free((void **) &GM_Elem_Types);

} /*      END of routine read_side_sets                                      */
