/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef ZOLTAN_HEAP_H
#define ZOLTAN_HEAP_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"


/* Heap datastructure */
typedef struct {
   int    space;
   int    n;
   int   *ele;
   int   *pos;
   float *value;
   } HEAP;

#define Zoltan_heap_empty(H)         (((H)->n)==0)
#define Zoltan_heap_not_empty(H)     (((H)->n)!=0)
#define Zoltan_heap_max_value(H)     ((H)->value[(H)->ele[0]])
#define Zoltan_heap_peek_max(H)      ((H)->ele[0])
#define Zoltan_heap_count(H)         ((H)->n)

int  Zoltan_heap_init         (ZZ*, HEAP*, int);
void Zoltan_heap_clear        (HEAP*);
void Zoltan_heap_free         (HEAP*);
int  Zoltan_heap_check        (HEAP*);
int  Zoltan_heap_input        (HEAP*, int, float);
int  Zoltan_heap_make         (HEAP*);
int  Zoltan_heap_change_value (HEAP*, int, float);
int  Zoltan_heap_extract_max  (HEAP*);
int  Zoltan_heap_extract      (HEAP*, int);



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* ZOLTAN_HEAP_H_ */
