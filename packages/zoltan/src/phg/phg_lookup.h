// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef __ZOLTAN_PHG_LOOKUP_H
#define __ZOLTAN_PHG_LOOKUP_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "third_library_tools.h"

/*****************************************************************************/

/* 
 * Structures to hold hypergraph data returned by query functions,
 * and hypergraph data gathered by processes to which edges/vertices
 * map to via a hash function.
 */
 
typedef struct _myObj{  /* Vertices returned in Get_Obj_List queries */
  int    *vtxHash;      /* Process to which GID hashes, temporary owner */
}zoltan_objects;

typedef struct _myPin{      /* Pins returned by hypergraph query functions */
  int           nHedges;    /* number of (partial) hyperedges */
  ZOLTAN_ID_PTR edgeGID;    /* edge global IDs */
  int           *esizes;    /* local size in pins of each hyperedge */
  ZOLTAN_ID_PTR pinGID;     /* global ID of pin vertex */
  int           numPins;    /* sum of esizes array */
  int           *edgeHash;  /* process assigned edgeGID by hash function */
}zoltan_pins;

typedef struct _myEW{     /* Values returned by edge weight query functions */
  int           size;       /* number of edges */
  ZOLTAN_ID_PTR edgeGID;   /* edge global IDs */
  int           *edgeHash;  /* process assigned this edge by hash function */
  float         *wgt;       /* weights supplied by query function for edge */
}zoltan_ews;

typedef struct _hshEdge{ /* Edges assigned to this process with hash func */
  ZOLTAN_ID_PTR edgeGID;    /* edge global IDs */
  ZOLTAN_ID_PTR pinGID;     /* vertex ID of each pin*/
  int           *pinHash;   /* process to which pin vertex is hashed */
}zoltan_temp_edges;

typedef struct _hshVtx{ /* Vertices assigned to this process with hash func */
  int           size;      /* number of vertices assigned to this process */
  ZOLTAN_ID_PTR vtxGID;    /* vertex global IDs  */
  int           *vtxOwner; /* process that returned vtx in Get_Obj_List  */
  ZOLTAN_GNO_TYPE *vtxGNO;   /* vertex global number */
}zoltan_temp_vertices;

/* 
 * A search structure, to find the index of a global ID in any of the
 * above structures.
 */

typedef struct _GID_lookup{
  struct Hash_Node *htTop;
  struct Hash_Node **ht;
  int table_size;
  int numGIDs;
  int lenGID;
}phg_GID_lookup;

/*****************************************************************************/

void phg_free_objects(zoltan_objects *zo);
void phg_free_pins(zoltan_pins *zp);
void phg_free_ews(zoltan_ews *zew);
void phg_free_temp_edges(zoltan_temp_edges *zte);
void phg_free_temp_vertices(zoltan_temp_vertices *ztv);

int phg_map_GIDs_to_processes(ZZ *zz, ZOLTAN_ID_PTR eid, int size, int lenGID, 
                              int **hashedProc, int nprocs);

phg_GID_lookup *phg_create_GID_lookup_table(ZOLTAN_ID_PTR gids, int size, int lenGID);
phg_GID_lookup *phg_create_GID_lookup_table2(ZOLTAN_ID_PTR gids, int ngids, int lenGID);
int phg_lookup_GID(phg_GID_lookup *lu, ZOLTAN_ID_PTR gid);
void phg_free_GID_lookup_table(phg_GID_lookup **lu);

#endif
