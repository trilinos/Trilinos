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

#define INT_CHANGE(A,B)         {int    _C_=(A);(A)=(B);(B)=_C_;}

int heap_init (HEAP *h, int space)
{ char *yo = "heap_init";

  h->space = space;
  h->n = 0;
  if (!(h->ele   = (int *)  ZOLTAN_CALLOC(space,sizeof(int))) ||
      !(h->pos   = (int *)  ZOLTAN_CALLOC(space,sizeof(int))) ||
      !(h->value = (float *)ZOLTAN_CALLOC(space,sizeof(float))) )
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  return ZOLTAN_OK;
}

void heap_free (HEAP *h)
{ ZOLTAN_FREE ((void **) &(h->ele));
  ZOLTAN_FREE ((void **) &(h->pos));
  ZOLTAN_FREE ((void **) &(h->value));
}

int heap_check (HEAP *h)
{ int i, left, right;

  for (i=0; i<h->n; i++)
  { left = 2*i+1;
    right = 2*i+2;
    if ((left <h->n && h->value[h->ele[left ]]>h->value[h->ele[i]]) ||
        (right<h->n && h->value[h->ele[right]]>h->value[h->ele[i]]))
    { printf("ERROR...no heap property!\n");
      return ZOLTAN_FATAL;
  } }
  return ZOLTAN_OK;
}

void heap_input (HEAP *h, int element, float value)
{
  h->ele[h->n] = element;
  h->pos[element] = (h->n)++;
  h->value[element] = value;
}

void heap_make (HEAP *h)
{ int i;
  
  for (i=h->n/2; i>=0; i--)
    heapify(h, i);
}

void heapify (HEAP *h, int root)
{ int	left=root*2+1, right=root*2+2, largest=root; 

  if ((left <h->n) && (h->value[h->ele[left ]]>h->value[h->ele[largest]]))
    largest = left;
  if ((right<h->n) && (h->value[h->ele[right]]>h->value[h->ele[largest]]))
    largest = right;
  if (largest != root)
  { h->pos[h->ele[root]] = largest;
    h->pos[h->ele[largest]] = root;
    INT_CHANGE(h->ele[root],h->ele[largest]);
    heapify(h,largest);
  }
}

void heap_change_value (HEAP *h, int element, float value)
{ int position=h->pos[element], father;

  if (h->value[element] > value)
  { h->value[element] = value;
    heapify(h,position);
  }
  else if (h->value[element] < value)
  { h->value[element] = value;
    father = (position-1)/2;
    while (position>0 && h->value[element]>h->value[h->ele[father]])
    { h->pos[h->ele[position]] = father;
      h->pos[h->ele[father]] = position;
      INT_CHANGE(h->ele[father],h->ele[position]);
      position = father;
      father = (father-1)/2;
  } }
}

int heap_extract_max (HEAP *h)
{ int max = h->ele[0];

  h->value[max] = 0.0;
  h->pos[h->ele[0]=h->ele[--(h->n)]] = 0;
  heapify(h,0);
  return max;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

