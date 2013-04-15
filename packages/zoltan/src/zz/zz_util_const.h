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


#ifndef __ZOLTAN_UTIL_CONST_H
#define __ZOLTAN_UTIL_CONST_H

#include <mpi.h>
#include "zz_const.h"
#include "zoltan_types.h"
#include "zz_hash.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef ZOLTAN_PURIFY
extern int Zoltan_strcasecmp(const char *, const char *);
extern int Zoltan_strncasecmp(const char *, const char *, size_t n);
#endif
extern int Zoltan_Clean_String(const char *, char **);
extern char *Zoltan_Strdup(const char *);
void Zoltan_Transform_Point( double *p, double (*m)[3], int *a, int d,
  int ndims, double *v);
void Zoltan_Transform_Box(double *lo, double *hi, double (*m)[3], int *a, 
  int d, int ndims);
void Zoltan_Transform_Box_Points(double *lo, double *hi, double (*m)[3], 
  int *a, int d, int ndims, double (*v)[3]);
int Zoltan_AllReduceInPlace(void *, int , MPI_Datatype , MPI_Op , MPI_Comm );
void Zoltan_write_linux_meminfo(int append, char *msg, int committedOnly);
long Zoltan_get_process_kilobytes();
int Zoltan_get_global_id_type(char **name);
int Zoltan_overflow_test(size_t val);

MPI_Datatype Zoltan_mpi_gno_type();
char *Zoltan_mpi_gno_name();

/* A Zoltan_Map is like a C++ STL map.  It uses Zoltan_Hash.
 */

struct Zoltan_Map_Entry{
  char *key;             /* pointer to arbitrary length key */
  intptr_t data;         /* value */
  struct Zoltan_Map_Entry *next;
};

typedef struct Zoltan_Map_Entry ZOLTAN_ENTRY;

struct Zoltan_Map_List{
  ZOLTAN_ENTRY **entries; /* hash array, length max_index + 1 */

  ZOLTAN_ENTRY *top;      /* if dynamicEntries==0, entries are here */
  char *keys;             /* If copyKeys and !dynamicEntries, keys are here */

  int key_size;             /* size in bytes of key */
  int num_zoltan_id_types;  /* number of ZOLTAN_ID_TYPES required to hold key */
  int max_index;        /* hash number range */
  int max_entries;      /* size of top array, or 0 if dynamicEntries == 1 */

  int prev_index;       /* index of top element returned in iterator */
  int prev_hash_index;  /* hash index of previous returned element */
  ZOLTAN_ENTRY *prev;   /* pointer to previous returned element */

  int dynamicEntries;   /* 1 - entries allocated as they are added */
			/* 0 - entries allocated at the start in top array */

  int copyKeys;         /* 1 - We create a copy of the added keys */
			/* 0 - We keep a pointer to the caller's copy of the key */

  int used;             /* 1 - this map is being used, 0 - it's free */
  int entry_count;      /* how many entries have been added to the map */
  ZOLTAN_ID_PTR zid;    /* buffer to hold key if it is not a multiple of ZOLTAN_ID_TYPEs */
};

typedef struct Zoltan_Map_List ZOLTAN_MAP;

ZOLTAN_MAP* Zoltan_Map_Create(ZZ *zz, int hash_range, int key_size_in_bytes, int store_keys, int num_entries);
int Zoltan_Map_Destroy(ZZ *zz, ZOLTAN_MAP **map);
int Zoltan_Map_Add(ZZ *zz, ZOLTAN_MAP *map, char *key, intptr_t data);
int Zoltan_Map_Find(ZZ *zz, ZOLTAN_MAP *map, char *key, intptr_t *data);
int Zoltan_Map_Find_Add(ZZ *zz, ZOLTAN_MAP* map, char *key, intptr_t datain, intptr_t *dataout);
int Zoltan_Map_Size(ZZ *zz, ZOLTAN_MAP *map);
int Zoltan_Map_First(ZZ *zz, ZOLTAN_MAP *map, char **key, intptr_t *data);
int Zoltan_Map_Next(ZZ *zz, ZOLTAN_MAP *map, char **key, intptr_t *data);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
