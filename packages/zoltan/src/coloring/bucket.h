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
#ifndef __BUCKET__H
#define __BUCKET__H


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


    /* This is ID-less bucket datastructure to save memory
       and hence speed up bucket updates. It assumes IDs are
       0-based indices without any gap */
typedef struct S
{
  struct S* prev;
  struct S* next;
} Bucket_element;


typedef struct arg
{
    Bucket_element **buckets; /* actual pointers to bucket heads */
    Bucket_element *elements; /* for direct access to bucket elements
                                 elements[id] is the id-th element */
    int          nb_elements; 
    int            max_value;
    int              *values; /* needed for update, incase bucket head
                                 changed. */

  int current_min_value;
} Bucket;

/* value == INT_MAX means not present in bucket */
void Zoltan_Bucket_Insert(Bucket* bs, int id, int value);

void Zoltan_Bucket_Update(Bucket* bs, int id, int new_value);

#define Zoltan_Bucket_DecVal(bs, id) Zoltan_Bucket_Update(bs, id, (bs)->values[id]-1)

/*returns -1 if empty*/
int Zoltan_Bucket_PopMin(Bucket* bs);

Bucket Zoltan_Bucket_Initialize(int max_value, int nb_element);

void Zoltan_Bucket_Free(Bucket* bs);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
    


#endif
