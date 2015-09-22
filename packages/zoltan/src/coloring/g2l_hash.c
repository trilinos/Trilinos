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
#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zoltan_util.h"
#include "coloring.h"    
#include "g2l_hash.h"



/*****************************************************************************/
/* Returns a prime number larger than (or equal to) stopafter,
   hopefully close to stopafter. If stopafter is larger than
   2147483647. 2147483647 will be returned.
 */

#define MAX_PRIME 193
    
int Zoltan_GenPrime(int stopafter, int *prime_num)
{
    static const int primes[MAX_PRIME]=
        {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 47, 53, 59, 67, 79, 89,
         101, 113, 127, 149, 167, 191, 211, 233, 257, 283, 313, 347, 383, 431,
         479, 541, 599, 659, 727, 809, 907, 1009, 1117, 1229, 1361, 1499, 1657,
         1823, 2011, 2213, 2437, 2683, 2953, 3251, 3581, 3943, 4339, 4783,
         5273, 5801, 6389, 7039, 7753, 8537, 9391, 10331, 11369, 12511, 13763,
         15149, 16673, 18341, 20177, 22229, 24469, 26921, 29629, 32603, 35869,
         39461, 43411, 47777, 52561, 57829, 63617, 69991, 76991, 84691, 93169,
         102497, 112757, 124067, 136481, 150131, 165161, 181693, 199873,
         219871, 241861, 266051, 292661, 321947, 354143, 389561, 428531,
         471389, 518533, 570389, 627433, 690187, 759223, 835207, 918733,
         1010617, 1111687, 1222889, 1345207, 1479733, 1627723, 1790501,
         1969567, 2166529, 2383219, 2621551, 2883733, 3172123, 3489347,
         3838283, 4222117, 4644329, 5108767, 5619667, 6181639, 6799811,
         7479803, 8227787, 9050599, 9955697, 10951273, 12046403, 13251047,
         14576161, 16033799, 17637203, 19400929, 21341053, 23475161, 25822679,
         28404989, 31245491, 34370053, 37807061, 41587807, 45746593, 50321261,
         55353391, 60888739, 66977621, 73675391, 81042947, 89147249, 98061979,
         107868203, 118655027, 130520531, 143572609, 157929907, 173722907,
         191095213, 210204763, 231225257, 254347801, 279782593, 307760897,
         338536987, 372390691, 409629809, 450592801, 495652109, 545217341,
         599739083, 659713007, 725684317, 798252779, 878078057, 965885863,
         1062474559, 1168722059, 1285594279, 1414153729, 1555569107,
         1711126033, 1882238639, 2070462533, 2147483647};

    int uplimit=MAX_PRIME-1;
    int botlimit=0;
    int j;
    int result=primes[uplimit];
    while(1) {
        if (uplimit-botlimit < 5) {
            for (j = botlimit; primes[j]<= stopafter && j <= uplimit; j++);
            result = (j==MAX_PRIME) ? j-1 : j;
	    break;
	}
	if (primes[botlimit + (uplimit - botlimit) / 2] < stopafter)
	    botlimit = botlimit + (uplimit-botlimit) / 2;
	else
            uplimit = uplimit - (uplimit-botlimit) / 2;
    }
    *prime_num = primes[result];
    return ZOLTAN_OK;
}

#undef MAX_PRIME
    

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
    hash->table = (G2LHashNode **) ZOLTAN_CALLOC((size_t) hash->maxsize, sizeof(G2LHashNode *));
    hash->nodes = (G2LHashNode *) ZOLTAN_MALLOC((size_t) hash->maxsize * (size_t) sizeof(G2LHashNode));
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
        i = Zoltan_Hash((ZOLTAN_ID_PTR) (void *)&gno, hash->num_gid_entries, (unsigned int) hash->maxsize);
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
        i = Zoltan_Hash((ZOLTAN_ID_PTR) (void *)&gno, hash->num_gid_entries, (unsigned int) hash->maxsize);
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

    i = Zoltan_Hash((ZOLTAN_ID_PTR) (void *)&key, hash->num_gid_entries, (unsigned int) hash->maxsize);
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

    i = Zoltan_Hash((ZOLTAN_ID_PTR) (void *) &key, hash->num_gid_entries, (unsigned int) hash->maxsize);
    for (ptr=hash->table[i]; ptr && ptr->gno!=key; ptr = ptr->next);
    if (!ptr)
        return -1;
    else
        return ptr->lno;
}
    


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
