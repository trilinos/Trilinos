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
#include "g2l_hash.h"
#include "coloring.h"

int Zoltan_G2LHash_Create(G2LHash *hash, int size)
{
    hash->table = NULL;
    hash->nodes = NULL;
    hash->size = size;
    hash->lastlno = 0;
    hash->table = (G2LHashNode **) ZOLTAN_CALLOC(size, sizeof(G2LHashNode *));
    hash->nodes = (G2LHashNode *) ZOLTAN_MALLOC(size * sizeof(G2LHashNode));
    if (!hash->table || !hash->nodes) {
        Zoltan_G2LHash_Destroy(hash);
        return ZOLTAN_MEMERR;
    }
    return ZOLTAN_OK;
}

int Zoltan_G2LHash_Destroy(G2LHash *hash)
{
    ZOLTAN_FREE(&hash->table);
    ZOLTAN_FREE(&hash->nodes);

    return ZOLTAN_OK;
}

int Zoltan_G2LHash_Insert(G2LHash *hash, int gno)
{
    int i, lno;
    G2LHashNode *ptr;

    i = Zoltan_Hash((ZOLTAN_ID_PTR) &gno, 1, (unsigned int) hash->size);
    for (ptr=hash->table[i]; ptr && ptr->gno!=gno; ptr = ptr->next);
    if (!ptr) {
        lno = hash->lastlno++;

        if (lno >= hash->size) {
            ZOLTAN_PRINT_ERROR(-1, "Zoltan_G2LHash_G2L", "Hash is full!");
            return -1;
        }
        
        ptr = &(hash->nodes[lno]);
        ptr->gno = gno;
        ptr->lno = lno;
        ptr->next = hash->table[i];
        hash->table[i] = ptr;
    } else
        lno = ptr->lno;

    return lno;
}

int Zoltan_G2LHash_G2L(G2LHash *hash, int gno)
{
    int i;
    G2LHashNode *ptr;

    i = Zoltan_Hash(&gno, 1, (unsigned int) hash->size);
    for (ptr=hash->table[i]; ptr && ptr->gno!=gno; ptr = ptr->next);
    if (!ptr)
        return -1;
    else
        return ptr->lno;
}
