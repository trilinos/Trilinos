// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include "zoltan_mem.h"
#include "bucket.h"    


void Zoltan_Bucket_Insert(Bucket* bs, int id, int value)
{
#if 0
    assert (bs != NULL);
    assert (value >= 0);
    assert (value <= bs->max_value);    
    assert (id >= 0);
    assert (id < bs->nb_elements);
#endif
    
    bs->values[id] = value;
    
    bs->elements[id].prev = NULL;
    bs->elements[id].next = bs->buckets[value];
    
    if (bs->buckets[value] != NULL) 
        bs->buckets[value]->prev = &(bs->elements[id]);
    else if (bs->current_min_value > value)
        bs->current_min_value = value;

    bs->buckets[value] = &(bs->elements[id]);

}

void Zoltan_Bucket_Update(Bucket* bs, int id, int new_value)
{
    int old_value = bs->values[id];

    if (old_value == INT_MAX)
        return;
    
#if 0  
    assert (bs != NULL);
    assert (new_value >= 0);
    assert (new_value <= bs->max_value);        
    assert (id >= 0);
    assert (id < bs->nb_elements);
#endif
  
    bs->values[id] = new_value;


    if (bs->elements[id].prev == NULL)
        bs->buckets[old_value] = bs->elements[id].next;
    else
        bs->elements[id].prev->next = bs->elements[id].next;
  
    if (bs->elements[id].next != NULL)
        bs->elements[id].next->prev = bs->elements[id].prev;

    Zoltan_Bucket_Insert(bs, id, new_value);
}

int Zoltan_Bucket_PopMin(Bucket* bs)
{
    int id;

#if 0  
    assert (bs != NULL);
    assert (bs->current_min_value >= 0);
#endif

    for (; bs->current_min_value<=bs->max_value; bs->current_min_value++) {
        if (bs->buckets[bs->current_min_value] != NULL) {
            id = bs->buckets[bs->current_min_value] - bs->elements;
            bs->buckets[bs->current_min_value] = bs->buckets[bs->current_min_value]->next;
            if (bs->buckets[bs->current_min_value] != NULL) {
                bs->buckets[bs->current_min_value]->prev = NULL;
            }
            bs->values[id] = INT_MAX;
            return id;
        }
    }
    return -1;
}

Bucket Zoltan_Bucket_Initialize(int max_value, int nb_element)
{
    Bucket bs;
    int i;

#if 0  
    assert (max_value>=0);
    assert (nb_element>=0);
#endif

    bs.buckets = (Bucket_element **)ZOLTAN_MALLOC(sizeof(Bucket_element *)*(max_value+1));
    bs.elements = (Bucket_element *)ZOLTAN_MALLOC(sizeof(Bucket_element)*nb_element);
    bs.values = (int *) ZOLTAN_MALLOC(sizeof(int)*nb_element);
    bs.max_value = max_value;
    bs.nb_elements = nb_element;

    if (bs.buckets == NULL || bs.elements == NULL || bs.values == NULL) {
        ZOLTAN_FREE(&(bs.values));
        ZOLTAN_FREE(&(bs.buckets));
        ZOLTAN_FREE(&(bs.elements));
    } else {
        for (i=0; i<=max_value; i++)
            bs.buckets[i] = NULL;

        for (i=0; i<nb_element; i++) {
            bs.elements[i].prev = NULL;
            bs.elements[i].next = NULL;
        }
    }
	
  
    bs.current_min_value = max_value+1;
  
    return bs;
}
    
void Zoltan_Bucket_Free(Bucket* bs)
{
    ZOLTAN_FREE(&(bs->values));
    ZOLTAN_FREE(&(bs->buckets));
    ZOLTAN_FREE(&(bs->elements));
} 


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
