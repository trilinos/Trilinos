/*BHEADER**********************************************************************
 * (c) 1999   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision$
 *********************************************************************EHEADER*/
/******************************************************************************
 *
 * Hash - Open addressing hash table with linear probing.  The keys are stored
 * in the table (also known as a closed table).  Conflicts are resolved with
 * linear probing, which may be faster in an environment with cache memory.
 *
 * We allow rehashing the data into a larger or smaller table, and thus
 * allow a data item (an integer, but a pointer would be more general)
 * to be stored with each key in the table.  (If we only return the 
 * storage location of the key in the table (the implied index), then 
 * rehashing would change the implied indices.)
 *
 * The modulus function is used as the hash function.
 * The keys must not equal HASH_EMPTY, which is -1.
 * The integer data associated with a key must not equal HASH_NOTFOUND,
 * which is -1.
 *
 *****************************************************************************/

#include <stdlib.h>
#include <assert.h>
#include "Common.h"
#include "Hash.h"

/*--------------------------------------------------------------------------
 * HashCreate - Return (a pointer to) a hash table of size "size".
 * "size" should be prime, if possible.
 *--------------------------------------------------------------------------*/

Hash *HashCreate(int size)
{
    int i, *p;

    Hash *h = (Hash *) malloc(sizeof(Hash));

    h->size  = size;
    h->num   = 0;
    h->keys  = (int *) malloc(size * sizeof(int));
    h->table = (int *) malloc(size * sizeof(int));
    h->data  = (int *) malloc(size * sizeof(int));

    /* Initialize the table to empty */
    p = h->table;
    for (i=0; i<size; i++)
        *p++ = HASH_EMPTY;

    return h;
}

/*--------------------------------------------------------------------------
 * HashDestroy - Destroy a hash table object "h".
 *--------------------------------------------------------------------------*/

void HashDestroy(Hash *h)
{
    free(h->keys);
    free(h->table);
    free(h->data);
    free(h);
}

/*--------------------------------------------------------------------------
 * HashLookup - Look up the "key" in hash table "h" and return the data 
 * associated with the key, or return HASH_NOTFOUND.
 *--------------------------------------------------------------------------*/

int HashLookup(Hash *h, int key)
{
    int loc;

    /* loc = key % h->size; */
    double keyd = key * 0.6180339887;
    loc = (int) (h->size * (keyd - (int) keyd));

    while (h->table[loc] != key)
    {
        if (h->table[loc] == HASH_EMPTY)
            return HASH_NOTFOUND;

        loc = (loc + 1) % h->size;
    }

    return h->data[loc];
}

/*--------------------------------------------------------------------------
 * HashInsert - Insert "key" with data "data" into hash table "h".
 * If the key is already in the hash table, the data item is replaced.
 *--------------------------------------------------------------------------*/

void HashInsert(Hash *h, int key, int data)
{
    int loc;

    /* loc = key % h->size; */
    double keyd = (double) key * 0.6180339887;
    loc = (int) ((double) h->size * (keyd - (int) keyd));

    while (h->table[loc] != key)
    {
        if (h->table[loc] == HASH_EMPTY)
        {
            assert(h->num < h->size);

	    h->keys[h->num++] = key;
            h->table[loc] = key;
            break;
        }

        loc = (loc + 1) % h->size;
    }

    h->data[loc] = data;
}

/*--------------------------------------------------------------------------
 * HashRehash - Given two hash tables, put the entries in one table into
 * the other.
 *--------------------------------------------------------------------------*/

void HashRehash(Hash *old, Hash *new)
{
    int i, data;

    for (i=0; i<old->num; i++)
    {
	data = HashLookup(old, old->keys[i]);
	HashInsert(new, old->keys[i], data);
    }
}

/*--------------------------------------------------------------------------
 * HashReset - Reset the hash table to all empty.
 *--------------------------------------------------------------------------*/

void HashReset(Hash *h)
{
    int i, *p;

    h->num = 0;
    p = h->table;
    for (i=0; i<h->size; i++)
	*p++ = HASH_EMPTY;
}

/*--------------------------------------------------------------------------
 * HashPrint - Print hash table to stdout.
 *--------------------------------------------------------------------------*/

void HashPrint(Hash *h)
{
    int i, j, *p;
    int lines = h->size/38;

    printf("Hash size: %d\n", h->size);

    p = h->table;
    for (i=0; i<lines; i++)
    {
	for (j=0; j<38; j++)
	    printf("%d ", ((*p++ == HASH_EMPTY) ? 0 : 1));
	    /*printf("%d ", *p++);*/
	printf("\n");
    }
}

