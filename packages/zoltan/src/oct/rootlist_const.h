
#ifndef __OCT_ROOTLIST_CONST_H
#define __OCT_ROOTLIST_CONST_H

#include "octant_const.h"

/* Structure used to store root octants for on and off processor roots */
/* The first and last nodes are sentinals */
/* Sentinals are signified by NULL oct fields */

typedef struct RL_Node *pRList;         /* typedef for a pointer to rootlist */
typedef struct RL_Node {                /* an entry in the local root list   */
  struct RL_Node *next;                 /* pointer to next node in the list  */
  pOctant oct;                          /* pointer to the root octant        */
} RList;

extern pRList  RL_initRootList();
extern pOctant RL_nextRootOctant(pRList *rlist);
extern int     RL_addRootOctant(pRList rlist, pOctant oct);
extern int     RL_delRootOctant(pRList rlist, pOctant oct);
extern int     RL_freeList(pRList rlist);
extern int     RL_numRootOctants(pRList);
extern int     RL_clearRootOctants(pRList rlist);
extern int     RL_printRootOctants(pRList rlist);
#endif /*__OCT_ROOTLIST_CONST_H*/
