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

#include "zz_heap.h"



#define INT_SWAP(A,B)         {int    _C_=(A);(A)=(B);(B)=_C_;}

/* This module implements a binary (max-) heap.
 * Three arrays are associated with a heap:
 *   ele   - ele[0] is the max element in the heap, ele[1] and
 *           ele[2] are its children, and so on.
 *   pos   - gives the current position in the heap for each element,
 *           where the elements are usually integers from 0 to n,
 *           for example, vertex numbers.
 *   value - key values (floats) by which the heap are arranged.
 *           Not in (heap) arranged order.
 *
 * EBEB: Idea for saving memory: Typically, the value array is stored
 * both in the application and a copy within the heap. We could
 * let the heap use an external array, or alternatively, make the
 * (internal) heap array visible to the application. Both approaches are
 * dangerous because the application developer must take care to ensure
 * the heap is properly updated when values change.
 */


/* prototypes */
static void heapify (HEAP*, int);


/* Inititializes the heap values and allocates space. */
/* Currently, we allocate space for as many heap entries as */
/* the range of the input elements. These two parameters may be */
/* decoupled in the future. */

int Zoltan_heap_init (ZZ *zz, HEAP *h, int space)
{
char *yo = "Zoltan_heap_init";
int i;

  h->space = space;
  h->n = 0;
  if ((space > 0)
   && (!(h->ele   = (int*)   ZOLTAN_CALLOC(space, sizeof(int)))
   ||  !(h->pos   = (int*)   ZOLTAN_CALLOC(space, sizeof(int)))
   ||  !(h->value = (float*) ZOLTAN_CALLOC(space, sizeof(float))) )) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
  for (i = 0; i < space; i++)
     h->pos[i] = -1;
  return ZOLTAN_OK;
}



void Zoltan_heap_clear (HEAP *h)
{
int i;

  h->n = 0;
  for (i = 0; i < h->space; i++) {
     h->ele[i]   = 0;
     h->value[i] = 0.0;
     h->pos[i]   = -1;
     }
}



/* Frees all memory and sets the heap value back to default */
void Zoltan_heap_free (HEAP *h)
{
  if (h->space != 0){
     Zoltan_Multifree(__FILE__, __LINE__, 3, &h->ele, &h->pos, &h->value);
     h->space = 0;
     }
  h->n = 0;
}



/* Checks wheather the heap has the Max-Heap property */
int Zoltan_heap_check (HEAP *h)
{
int i, left, right;
static char *yo = "Zoltan_heap_check";

  for (i = 0; i < h->n; i++) {
     left  = 2*i + 1;
     right = 2*i + 2;
     if ((left < h->n && h->value[h->ele[left ]] > h->value[h->ele[i]])
      || (right< h->n && h->value[h->ele[right]] > h->value[h->ele[i]])) {
         ZOLTAN_PRINT_ERROR(0, yo, "No heap property!\n");
         return ZOLTAN_FATAL;
         }
     }
  return ZOLTAN_OK;
}



/* Zoltan_heap_input adds one item to the heap but does NOT rearrange the
   heap! Constant time. We might want to write a function Zoltan_heap_insert
   that adds an item and preserves the heap property. */

int Zoltan_heap_input (HEAP *h, int element, float value)
{
static char *yo = "Zoltan_heap_input";

  if (element >= h->space) {
     ZOLTAN_PRINT_ERROR(0, yo, "Inserted heap element out of range!\n");
     return ZOLTAN_FATAL;
     }
  if (h->n >= h->space) {
     ZOLTAN_PRINT_ERROR(0, yo, "Heap is full!\n");
     return ZOLTAN_FATAL;
     }
  h->value[element] = value;
  h->pos[element]   = h->n;
  h->ele[(h->n)++]  = element;
  return ZOLTAN_OK;
}



/* Moves the values in the heap to gain the Max-Heap property. Linear time */
int Zoltan_heap_make (HEAP *h)
{ int i;

  for (i = h->n / 2; i >= 0;  i--)
     heapify(h, i);
  return ZOLTAN_OK;
}



/* Subroutine which gets the heap property if both subtrees already
   have the heap property. */
static void heapify (HEAP *h, int root)
{
int left = root*2 + 1, right = root*2 + 2, largest = root;

  if ((left < h->n)  && (h->value[h->ele[left ]] > h->value[h->ele[largest]]))
     largest = left;
  if ((right < h->n) && (h->value[h->ele[right]] > h->value[h->ele[largest]]))
     largest = right;
  if (largest != root) {
     h->pos[h->ele[root]]    = largest;
     h->pos[h->ele[largest]] = root;
     INT_SWAP(h->ele[root], h->ele[largest]);
     heapify(h, largest);
     }
}



/* Changes the value of an element in the heap and restores the
   heap property. This can take O(log(n)) time */
int Zoltan_heap_change_value (HEAP *h, int element, float value)
{
int position, father;

  if (element < 0 || element >= h->space)
     return ZOLTAN_FATAL;                           /* Error */

  if ((position = h->pos[element]) >= 0) {
     if (value < h->value[element]) {
        h->value[element] = value;
        heapify(h,position);
        }
     else if (value > h->value[element]) {
        h->value[element] = value;
        father = (position-1)/2;
        while (position > 0 && h->value[element] > h->value[h->ele[father]]) {
           h->pos[h->ele[position]] = father;
           h->pos[h->ele[father]]   = position;
           INT_SWAP(h->ele[father], h->ele[position]);
           position = father;
           father   = (father-1)/2;
           }
        }
     }
  return ZOLTAN_OK;
}




/* Extracts the maximum element & restores the heap property. Time O(log(n))*/
int Zoltan_heap_extract_max (HEAP *h)
{
int max;

  if (h->n == 0)
     return -1;           /* No elements in heap. */

  max     = h->ele[0];

  h->value[max] = 0.0;
  h->pos[max] = -1;
  h->pos[h->ele[0] = h->ele[--(h->n)]] = 0;

  heapify(h,0);
  return max;
}


/* Extracts an element from the heap */
int Zoltan_heap_extract (HEAP *h, int element)
{
  int position;

  if (element < 0 || element >= h->space)
     return ZOLTAN_FATAL;
  if (h->n == 0)
     return ZOLTAN_FATAL;

  position = h->pos[element];

  h->value[element] = 0.0;
  h->pos[element] = -1;
  h->pos[h->ele[position] = h->ele[--(h->n)]] = position;

  heapify(h,position);

  return ZOLTAN_OK;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

