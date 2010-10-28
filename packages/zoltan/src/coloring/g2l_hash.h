/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
#ifndef _G2L_HASH_H_
#define _G2L_HASH_H_

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Structure used for hashing */
struct G2L_Hash_Node {
    int gno;           /* Global number */
    int lno;           /* Mapped id of gno*/
    struct G2L_Hash_Node * next;
};

typedef struct G2L_Hash_Node G2LHashNode;

struct G2L_Hash {
    int   maxsize;
    int   size;          /* number of ids stored in the hash */
    int   base, baseend; /* base and baseend are inclusive gno's of local vertices */
    int   nlvtx;         /* it is #localy owened vertices: simply equal to "baseend-base+1" */
    
    G2LHashNode **table;
    G2LHashNode *nodes;
};

typedef struct G2L_Hash G2LHash;

int Zoltan_G2LHash_Create(G2LHash *hash, int maxsize, int base, int nlvtx);
int Zoltan_G2LHash_Destroy(G2LHash *hash);
int Zoltan_G2LHash_G2L(G2LHash *hash, int gno);
/*
  if gno exist it returns lno, if it does not exist,
  it inserts andr returns newly assigned lno */
int Zoltan_G2LHash_Insert(G2LHash *hash, int gno);
    
#define Zoltan_G2LHash_L2G(hash, lno) ((lno<(hash)->nlvtx) ? (hash)->base+lno : (hash)->nodes[lno-(hash)->nlvtx].gno)


/* Key&Value hash functions using same data structure above
   the only difference will be the insert function */
typedef struct G2L_Hash KVHash;

int Zoltan_KVHash_Create(KVHash *hash, int maxsize);
int Zoltan_KVHash_Destroy(KVHash *hash);

int Zoltan_KVHash_Insert(KVHash *hash, int key, int value);
int Zoltan_KVHash_GetValue(KVHash *hash, int key);

    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
    

#endif

