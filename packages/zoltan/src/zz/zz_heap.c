#include "hypergraph.h"

#define INT_CHANGE(A,B)         {int    _C_=(A);(A)=(B);(B)=_C_;}

int heap_init (HEAP *heap, int n)
{ char *yo = "heap_init";

  heap->space = n;
  heap->n = 0;

  if (!(heap->ele = (int *)  ZOLTAN_CALLOC(n,sizeof(int))) ||
      !(heap->pos = (int *)  ZOLTAN_CALLOC(n,sizeof(int))) ||
      !(heap->d   = (float *)ZOLTAN_CALLOC(n,sizeof(float))) )
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  return ZOLTAN_OK;
}

void heap_free (HEAP *heap)
{ ZOLTAN_FREE ((void **) &(heap->ele));
  ZOLTAN_FREE ((void **) &(heap->pos));
  ZOLTAN_FREE ((void **) &(heap->d));
}

int heap_check (HEAP *heap)
{ int i, left, right;

  for (i=0; i<heap->n; i++)
  { left = 2*i+1;
    right = 2*i+2;
    if ((left<heap->n && heap->d[heap->ele[left]]>heap->d[heap->ele[i]]) ||
        (right<heap->n && heap->d[heap->ele[right]]>heap->d[heap->ele[i]]))
    { printf("ERROR...no heap property!\n");
      return ZOLTAN_FATAL;
    }
  }
  return ZOLTAN_OK;
}

void heap_input (HEAP *heap, int element, float d)
{
  heap->ele[heap->n] = element;
  heap->pos[element] = (heap->n)++;
  heap->d[element] = d;
}

void heap_make (HEAP *heap)
{ int i;
  
  for (i=heap->n/2; i>=0; i--)
    heapify(heap, i);
}

void heapify (HEAP *heap, int root)
{ int	left=root*2+1, right=root*2+2, largest=root; 

  if ((left < heap->n) && (heap->d[heap->ele[left]] > heap->d[heap->ele[largest]]))
    largest = left;
  if((right < heap->n) && (heap->d[heap->ele[right]] > heap->d[heap->ele[largest]]))
    largest = right;
  if (largest != root)
  { heap->pos[heap->ele[root]] = largest;
    heap->pos[heap->ele[largest]] = root;
    INT_CHANGE(heap->ele[root],heap->ele[largest]);
    heapify(heap,largest);
  }
}

void heap_change_key (HEAP *heap, int element, float d)
{ int position=heap->pos[element], father;

  if (heap->d[element] > d)
  { heap->d[element] = d;
    heapify(heap,position);
  }
  else if (heap->d[element] < d)
  { heap->d[element] = d;
    father = (position-1)/2;
    while (position>0 && heap->d[element]>heap->d[heap->ele[father]])
    { heap->pos[heap->ele[position]] = father;
      heap->pos[heap->ele[father]] = position;
      INT_CHANGE(heap->ele[father],heap->ele[position]);
      position = father;
      father = (father-1)/2;
    } 
  }
}

int heap_extract_max (HEAP *heap)
{ int max = heap->ele[0];
  heap->d[max] = 0.0;
  heap->pos[heap->ele[0]=heap->ele[--(heap->n)]] = 0;
  heapify(heap,0);
  return max;
}
