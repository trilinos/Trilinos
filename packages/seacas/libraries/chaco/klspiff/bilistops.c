/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h" // for bilist
#include <stdio.h>   // for NULL

/* Note: bi-directional lists aren't assumed to be sorted. */

void add2bilist(                      /* add val to unsorted list */
                struct bilist * lptr, /* element to add */
                struct bilist **list  /* list added to */
)
{
  lptr->next = *list;
  if (*list != NULL) {
    (*list)->prev = lptr;
  }
  lptr->prev = NULL;
  *list      = lptr;
}

void removebilist(struct bilist * lptr, /* ptr to element to remove */
                  struct bilist **list  /* head of list to remove it from */
)

/* Remove an element from a bidirectional list. */
{
  if (lptr->next != NULL) {
    lptr->next->prev = lptr->prev;
  }
  if (lptr->prev != NULL) {
    lptr->prev->next = lptr->next;
  }
  else {
    *list = lptr->next;
  }
}

void movebilist(struct bilist * lptr,    /* ptr to element to move */
                struct bilist **oldlist, /* head of list to remove it from */
                struct bilist **newlist  /* head of list to add it to */
)

/* Move an element from a old bidirectional list to new one. */
{
  removebilist(lptr, oldlist);

  add2bilist(lptr, newlist);
}
