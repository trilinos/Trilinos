/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "defs.h"
#include "structs.h"
#include "smalloc.h"

struct scanlink *mkscanlist(depth)
int       depth;
{
    struct scanlink *prevlnk;
    struct scanlink *newlnk;
    int       i;


    prevlnk = smalloc(sizeof(struct scanlink));
    prevlnk->pntr = NULL;
    newlnk = prevlnk;		/* in case the list is one long */
    for (i = 1; i <= (depth - 1); i++) {
	newlnk = smalloc(sizeof(struct scanlink));
	newlnk->pntr = prevlnk;
	prevlnk = newlnk;
    }
    return (newlnk);
}
