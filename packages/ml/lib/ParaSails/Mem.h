/*BHEADER**********************************************************************
 * (c) 1999   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision$
 *********************************************************************EHEADER*/
/******************************************************************************
 *
 * Mem.h header file.
 *
 *****************************************************************************/

#include <stdio.h>

#ifndef _MEM_H
#define _MEM_H

#define MEM_BLOCKSIZE (2*1024*1024)
#define MEM_MAXBLOCKS 1024

typedef struct
{
    int   num_blocks;
    int   bytes_left;

    long  total_bytes;
    long  bytes_alloc;
    int   num_over;

    char *avail;
    char *blocks[MEM_MAXBLOCKS];
}
Mem;

Mem  *MemCreate();
void  MemDestroy(Mem *m);
char *MemAlloc(Mem *m, int size);
void  MemStat(Mem *m, FILE *stream, char *msg);

#endif /* _MEM_H */
