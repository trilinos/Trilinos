/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "hypergraph.h"

#define INT_SWAP(A,B)         {int    _C_=(A);(A)=(B);(B)=_C_;}

/* This module implements a binary (max-) heap.
 * Three arrays are associated with a heap:
 *   ele   - the elements (ints between 0 and n; for example, 
 *           vertex numbers) that are stored in the heap. 
 *   pos   - gives the current position in the heap for each element
 *   value - key values (floats) by which the heap are arranged.
 *           Not in arranged order.
 */

static void heapify (HEAP*, int);

int heap_init (ZZ *zz, HEAP *h, int space)
{ char *yo = "heap_init";
  int i;

  h->space = space;
  h->n = 0;
  if ((space>0) &&
      (!(h->ele   = (int *)  ZOLTAN_CALLOC(space,sizeof(int))) ||
       !(h->pos   = (int *)  ZOLTAN_CALLOC(space,sizeof(int))) ||
       !(h->value = (float *)ZOLTAN_CALLOC(space,sizeof(float))) ))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<space; i++)
    h->pos[i] = -1; 
  return ZOLTAN_OK;
}

void heap_free (HEAP *h)
{ ZOLTAN_FREE ((void **) &(h->ele));
  ZOLTAN_FREE ((void **) &(h->pos));
  ZOLTAN_FREE ((void **) &(h->value));
  h->space = 0;
  h->n = 0;
}

int heap_check (HEAP *h)
{ int i, left, right;
  static char * yo = "heap_check";

  for (i=0; i<h->n; i++)
  { left = 2*i+1;
    right = 2*i+2;
    if ((left <h->n && h->value[h->ele[left ]]>h->value[h->ele[i]]) ||
        (right<h->n && h->value[h->ele[right]]>h->value[h->ele[i]]))
    { ZOLTAN_PRINT_ERROR(0, yo, "No heap property!\n");
      return ZOLTAN_FATAL;
  } }
  return ZOLTAN_OK;
}

/* heap_input adds one item to the heap but does NOT
 * rearrange the heap! We need a function heap_insert
 * to add an item and preserve the heap property. 
 */
int heap_input (HEAP *h, int element, float value)
{
  static char *yo = "heap_input";

  if (element >= h->space){
    ZOLTAN_PRINT_ERROR(0, yo, "Inserted heap element out of range!\n");
    return ZOLTAN_FATAL;
  }
  if (h->n >= h->space){
    ZOLTAN_PRINT_ERROR(0, yo, "Heap is full!\n");
    return ZOLTAN_FATAL;
  }
  h->value[element] = value;
  h->pos[element] = h->n;
  h->ele[(h->n)++] = element;
  return ZOLTAN_OK;
}

int heap_make (HEAP *h)
{ int i;
  
  for (i=h->n/2; i>=0; i--)
    heapify(h, i);
  return ZOLTAN_OK;
}

static void heapify (HEAP *h, int root)
{ int	left=root*2+1, right=root*2+2, largest=root; 

  if ((left <h->n) && (h->value[h->ele[left ]]>h->value[h->ele[largest]]))
    largest = left;
  if ((right<h->n) && (h->value[h->ele[right]]>h->value[h->ele[largest]]))
    largest = right;
  if (largest != root)
  { h->pos[h->ele[root]] = largest;
    h->pos[h->ele[largest]] = root;
    INT_SWAP(h->ele[root],h->ele[largest]);
    heapify(h,largest);
  }
}

int heap_change_value (HEAP *h, int element, float value)
{ int position, father;

  if ((element<0) || (element>=h->space)){
    return ZOLTAN_FATAL; /* Error */
  }

  position = h->pos[element];
  if (position >= 0)
  { if (value < h->value[element])
    { h->value[element] = value;
      heapify(h,position);
    }
    else if (value > h->value[element])
    { h->value[element] = value;
      father = (position-1)/2;
      while (position>0 && h->value[element]>h->value[h->ele[father]])
      { h->pos[h->ele[position]] = father;
        h->pos[h->ele[father]] = position;
        INT_SWAP(h->ele[father],h->ele[position]);
        position = father;
        father = (father-1)/2;
  } } }
  return ZOLTAN_OK;
}

int heap_extract_max (HEAP *h)
{ int max = h->ele[0];

  h->value[max] = 0.0;
  h->pos[max] = -1;
  h->pos[h->ele[0]=h->ele[--(h->n)]] = 0;
  heapify(h,0);
  return max;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

