/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_gr_list_c = "$Id$";
#endif

#include "all_const.h"
#include "gr_list_const.h"
#include "id_util_const.h"
#include "vx_const.h"

/*****************************************************************************/
/*
 *  Routines for maintaining linked lists.  For now, we do NOT sort the lists.
 */
/*****************************************************************************/

static void delete_list_entry(LIST_ENTRY **ptr);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_add_to_bucket(LIST_ENTRY **head_list, VERTEX *vertex)
{
/*
 *  This function creates a new list entry structure and inserts it as the 
 *  first entry in the doubly-linked list headed by head_list.
 *
 *  Input:
 *    head_list:    Pointer to the pointer to the first entry in the list.
 *    vertex:       Pointer to the data to be entered into the list.
 */

char *yo = "LB_add_to_bucket";
LIST_ENTRY *new;

  new = (LIST_ENTRY *) LB_MALLOC(sizeof(LIST_ENTRY));
  if (!new) {
    fprintf(stderr, "Error from %s: Insufficient memory\n", yo);
    exit(-1);
  }

  new->Next = *head_list;
  new->Prev = NULL;
  new->Vertex = vertex;

  if (*head_list != NULL) {
    (*head_list)->Prev = new;
  }
  *head_list = new;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_remove_from_bucket(LIST_ENTRY **head_list, VERTEX *vertex)
{
/*
 *  Routine to remove a vertex from a hash bucket.  It does NOT free the
 *  memory for the vertex.
 */

char *yo = "LB_remove_from_bucket";
LIST_ENTRY *ptr = *head_list;

  while (ptr != NULL && ptr->Vertex != vertex) 
    ptr = ptr->Next;

  if (ptr == NULL) {
    fprintf(stderr, "Error in %s:  Item %x not found in bucket.\n");
    exit(-1);
  }

  if (ptr == *head_list) {
    *head_list = ptr->Next;
  }

  delete_list_entry(&ptr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

VERTEX *LB_search_bucket(LIST_ENTRY *head_list, ID *id)
{
/*
 *  Given a bucket, return a pointer to a vertex based on its ID.
 */
LIST_ENTRY *ptr = head_list;
VERTEX *vertex = NULL;

  while (ptr != NULL && !BL_ID_Util.Compare(&(ptr->Vertex->Id), id))
    ptr = ptr->Next;
  
  if (ptr != NULL) 
    vertex = ptr->Vertex;

  return (vertex);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_free_bucket(LIST_ENTRY **head_list) 
{
/*
 *  Frees the list representing a hash bucket.  Also frees the 
 *  vertex storage for vertices in the bucket.
 */

LIST_ENTRY *ptr, *next_ptr;

  next_ptr = *head_list;

  while ((ptr = next_ptr) != NULL) {
    next_ptr = ptr->Next;
    LB_free_vertex(&(ptr->Vertex));
    delete_list_entry(&ptr);
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void delete_list_entry(LIST_ENTRY **ptr) 
{
LIST_ENTRY *prev, *next;

  prev = (*ptr)->Prev;
  next = (*ptr)->Next;

  if (prev != NULL) {
    prev->Next = next;
  }

  if (next != NULL) {
    next->Prev = prev;
  }

  LB_FREE(ptr);
}
