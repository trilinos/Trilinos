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
#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zoltan_util.h"
#include "coloring.h"    
#include "g2l_hash.h"



/*****************************************************************************/
/* Returns the "next" closest prime number after (or equal to) stopafter */
int Zoltan_GenPrime(int stopafter, int *prime_num)
{
    int nap, cont=1;
    int num, c, i;
    int *prime;

    prime = (int *) ZOLTAN_MALLOC(((stopafter/2) + 7) * sizeof(int));
    if (!prime)
        return ZOLTAN_MEMERR;
    
    prime[0] = 2; prime[1] = 3;
    c = 2; /* initial primes */
    
    /* only have to check odd numbers */
    for (num=5; cont; num = num + 2) {
        nap = 0;  /* set not-a-prime false */        
        /* cycle through list of known primes */
        for (i=0; i < c; i++) { 
            /* check if a previous prime divides evenly */
            /* if so the number is not a prime */
            if ((num % prime[i]) == 0) {
                nap = 1;
                break;
            }            
            /* stop if prime squared is bigger than the number */
            if ((prime[i] * prime[i]) > num)
                break;
        }        
        /* if not-a-prime, then we found a prime */
        if (nap != 1) {
           /* add prime to list of known primes */
            prime[c] = num;
            c++;
            if (num >= stopafter)
                cont = 0;
        }
    }
    *prime_num = prime[c-1];

    ZOLTAN_FREE(&prime);
    return ZOLTAN_OK;
}



int Zoltan_G2LHash_Create(G2LHash *hash, int maxsize, ZOLTAN_GNO_TYPE base, int nlvtx)
{
    if (maxsize == 0) /* to avoid memory allocation errors */
        maxsize = 1;

    if (Zoltan_GenPrime(maxsize , &(hash->maxsize))==ZOLTAN_MEMERR)
      return ZOLTAN_MEMERR;
    
    hash->table = NULL;
    hash->nodes = NULL;
    hash->base = base;
    hash->baseend = base+nlvtx-1;
    hash->nlvtx = nlvtx;
    hash->size = 0;
    hash->num_gid_entries = sizeof(ZOLTAN_GNO_TYPE) / sizeof(ZOLTAN_ID_TYPE);
    hash->table = (G2LHashNode **) ZOLTAN_CALLOC(hash->maxsize, sizeof(G2LHashNode *));
    hash->nodes = (G2LHashNode *) ZOLTAN_MALLOC(hash->maxsize * sizeof(G2LHashNode));
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

int Zoltan_G2LHash_Insert(G2LHash *hash, ZOLTAN_GNO_TYPE gno)
{
    int i, lno;
    G2LHashNode *ptr;

    if (gno<hash->base || gno>hash->baseend) {
        i = Zoltan_Hash((ZOLTAN_ID_PTR) &gno, hash->num_gid_entries, (unsigned int) hash->maxsize);
        for (ptr=hash->table[i]; ptr && ptr->gno!=gno; ptr = ptr->next);
        if (!ptr) {
            if (hash->size >= hash->maxsize) {
                char st[2048];
                sprintf(st, "Hash is full! #entries=%d  maxsize=%d", hash->size, hash->maxsize);
                ZOLTAN_PRINT_ERROR(-1, "Zoltan_G2LHash_G2L", st);
                return -1;
            }
            ptr = &(hash->nodes[hash->size]);
            ptr->gno = gno;
            lno = ptr->lno = hash->nlvtx + hash->size;
            ptr->next = hash->table[i];
            hash->table[i] = ptr;
            ++hash->size;
        } else
            lno = ptr->lno;
    } else
        return gno-hash->base;

    return lno;
}

int Zoltan_G2LHash_G2L(G2LHash *hash, ZOLTAN_GNO_TYPE gno)
{
    int i;
    G2LHashNode *ptr;

    if (gno<hash->base || gno>hash->baseend) {
        i = Zoltan_Hash((ZOLTAN_ID_PTR) &gno, hash->num_gid_entries, (unsigned int) hash->maxsize);
        for (ptr=hash->table[i]; ptr && ptr->gno!=gno; ptr = ptr->next);
        if (!ptr)
            return -1;
        else
            return ptr->lno;
    } else
        return gno-hash->base;
}


/* --------------------------- KVHash -------------------------------- */

int Zoltan_KVHash_Create(KVHash *hash, int maxsize)
{
    if (maxsize == 0) /* to avoid memory allocation errors */
        maxsize = 1;

    if (Zoltan_GenPrime(maxsize , &(hash->maxsize))==ZOLTAN_MEMERR)
      return ZOLTAN_MEMERR;
    
    hash->table = NULL;
    hash->nodes = NULL;
    hash->size = 0;
    hash->num_gid_entries = sizeof(ZOLTAN_GNO_TYPE) / sizeof(ZOLTAN_ID_TYPE);
    hash->table = (G2LHashNode **) ZOLTAN_CALLOC(hash->maxsize, sizeof(G2LHashNode *));
    hash->nodes = (G2LHashNode *) ZOLTAN_MALLOC(hash->maxsize * sizeof(G2LHashNode));
    if (!hash->table || !hash->nodes) {
        Zoltan_G2LHash_Destroy(hash);
        return ZOLTAN_MEMERR;
    }
    return ZOLTAN_OK;
}

int Zoltan_KVHash_Destroy(KVHash *hash)
{
    ZOLTAN_FREE(&hash->table);
    ZOLTAN_FREE(&hash->nodes);

    return ZOLTAN_OK;
}    

int Zoltan_KVHash_Insert(KVHash *hash, ZOLTAN_GNO_TYPE key, int value)
{
    int i;
    
    G2LHashNode *ptr;

    i = Zoltan_Hash((ZOLTAN_ID_PTR) &key, hash->num_gid_entries, (unsigned int) hash->maxsize);
    for (ptr=hash->table[i]; ptr && ptr->gno!=key; ptr = ptr->next);
    if (!ptr) {
        if (hash->size >= hash->maxsize) {
            ZOLTAN_PRINT_ERROR(-1, "Zoltan_KVHash_Insert", "Hash is full!");
            return -1;
        }
        
        ptr = &(hash->nodes[hash->size]);
        ptr->gno = key;
        ptr->lno = value;
        ptr->next = hash->table[i];
        hash->table[i] = ptr;
        ++hash->size;
    } else
        value = ptr->lno;

    return value;   
}


int Zoltan_KVHash_GetValue(KVHash *hash, ZOLTAN_GNO_TYPE key)    
{
    int i;
    G2LHashNode *ptr;

    i = Zoltan_Hash((ZOLTAN_ID_PTR) &key, hash->num_gid_entries, (unsigned int) hash->maxsize);
    for (ptr=hash->table[i]; ptr && ptr->gno!=key; ptr = ptr->next);
    if (!ptr)
        return -1;
    else
        return ptr->lno;
}
    


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
