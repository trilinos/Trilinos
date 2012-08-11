/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */
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

