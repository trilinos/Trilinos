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
