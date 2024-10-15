// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef _G2L_HASH_H_
#define _G2L_HASH_H_

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Structure used for hashing */
struct G2L_Hash_Node {
    ZOLTAN_GNO_TYPE gno;           /* Global number */
    int lno;           /* Mapped id of gno*/
    struct G2L_Hash_Node * next;
};

typedef struct G2L_Hash_Node G2LHashNode;

struct G2L_Hash {
    int   maxsize;
    int   size;          /* number of ids stored in the hash */
    ZOLTAN_GNO_TYPE base, baseend; /* base and baseend are inclusive gno's of local vertices */
    int   nlvtx;         /* it is #localy owened vertices: simply equal to "baseend-base+1" */
    int   num_gid_entries;   /* multiple of ZOLTAN_ID_TYPEs in a key */
    
    G2LHashNode **table;
    G2LHashNode *nodes;
};

typedef struct G2L_Hash G2LHash;

int Zoltan_G2LHash_Create(G2LHash *hash, int maxsize, ZOLTAN_GNO_TYPE base, int nlvtx);
int Zoltan_G2LHash_Destroy(G2LHash *hash);
int Zoltan_G2LHash_G2L(G2LHash *hash, ZOLTAN_GNO_TYPE gno);
/*
  if gno exist it returns lno, if it does not exist,
  it inserts andr returns newly assigned lno */
int Zoltan_G2LHash_Insert(G2LHash *hash, ZOLTAN_GNO_TYPE gno);
    
#define Zoltan_G2LHash_L2G(hash, lno) ((lno<(hash)->nlvtx) ? (hash)->base+lno : (hash)->nodes[lno-(hash)->nlvtx].gno)


/* Key&Value hash functions using same data structure above
   the only difference will be the insert function */
typedef struct G2L_Hash KVHash;

int Zoltan_KVHash_Create(KVHash *hash, int maxsize);
int Zoltan_KVHash_Destroy(KVHash *hash);

int Zoltan_KVHash_Insert(KVHash *hash, ZOLTAN_GNO_TYPE key, int value);
int Zoltan_KVHash_GetValue(KVHash *hash, ZOLTAN_GNO_TYPE key);

    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
    

#endif
