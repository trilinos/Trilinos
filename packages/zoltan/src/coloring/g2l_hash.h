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

#ifndef _COLOR_HASH_H_
#define _COLOR_HASH_H_

/* Structure used for hashing */
struct G2L_Hash_Node {
    int gno;           /* Global number */
    int lno;           /* Mapped id of gno*/
    struct G2L_Hash_Node * next;
};

typedef struct G2L_Hash_Node G2LHashNode;

struct G2L_Hash {
    int   size;
    int   lastlno;
    
    G2LHashNode **table;
    G2LHashNode *nodes;
};

typedef struct G2L_Hash G2LHash;


int Zoltan_G2LHash_Create(G2LHash *hash, int size);
int Zoltan_G2LHash_Destroy(G2LHash *hash);
int Zoltan_G2LHash_G2L(G2LHash *hash, int gno);
#define Zoltan_G2LHash_L2G(hash, lno) ((hash)->nodes[lno].gno)


#endif

