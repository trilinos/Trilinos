// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdlib.h>
#include <zz_util_const.h>
#include <zoltan_mem.h>

#define key_match(m, k1, k2) (memcmp(k1, k2, m->key_size) == 0);
#define current_key(m, key) (m->zid ? memcpy(m->zid, key, m->key_size) : key);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 * A Zoltan_Map is like a C++ STL map.  It uses Zoltan_Hash to maintain a
 * a table mapping keys values.  Keys are unique and of arbitrary length.
 * The value must fit in an intptr_t.
 *
 * It supports:
 *
 *  Zoltan_Map_Create       Initialize the map
 *  Zoltan_Map_Add          Add an element to the map
 *  Zoltan_Map_Find         Find an element in the map and return it.
 *  Zoltan_Map_Find_Add     Try to find an element in the map, if not present, add it.
 *  Zoltan_Map_First        Return the first element in the map (like an iterator)
 *  Zoltan_Map_Next         Return the next element in the map
 *  Zoltan_Map_Size         Return the number of elements in the map
 *  Zoltan_Map_Destroy      Destroy the memory allocated for the map
 *
 * There's no delete.
 *
 * If you know the maximum number of unique items that will be added to the
 * map, include that as an argument in Zoltan_Map_Create.  Else if this value
 * is zero, the entries will be created dynamically (more time, less memory).
 *
 * The Zoltan_Map can store a copy of the key, or just store a pointer to
 * the caller's copy of the key.
 *
 * NOTE:
 * Zoltan_Map will not be efficient when num_entries = 0, in Zoltan_Map_Create
 * and later large number of entries are added. The hash table size is hard
 * coded to 1000. num_entries in the same order should be ok.
 */

/*****************************************************************
 * Create a new map.
 * Return the map number to be used in subsequent calls.
 * (Return NULL on error)
 *
 * Originally, Zoltan_Map_* assumed the key was a Zoltan global ID.  There
 * is Zoltan code that used Zoltan_Map with a key that is an integer.  This
 * worked when ZOLTAN_ID_TYPE was always an integer.
 *
 * But now ZOLTAN_ID_TYPE can be set at compile time.  So the third parameter,
 * instead of being zz->Num_GID, should be the number of bytes in the key.
 */

ZOLTAN_MAP* Zoltan_Map_Create(ZZ *zz,     /* just need this for error messages */
		      int hash_range_max, /* > 0: maximum hash value */
					  /*   0: Zoltan_Map_Create will choose value */
		      int num_bytes,       /* length of a key */
		      int store_keys,     /* 1 - keep a copy of each key */
					  /* 0 - keep a pointer to caller's key */
		      int num_entries)    /* > 0: number of keys that will be added */
					  /*   0: dynamically allocate entries */
{
  char *yo = "Zoltan_Map_Create";
  ZOLTAN_ENTRY *top = NULL;
  ZOLTAN_ENTRY **entries = NULL;
  ZOLTAN_MAP* map = NULL;
  char *keys = NULL;
  int num, rem;
  int my_range ;

  if ((num_bytes < 1) || (num_entries < 0) || (hash_range_max < 0)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Bad parameters\n");
    return NULL;
  }

  map = (ZOLTAN_MAP*) ZOLTAN_MALLOC(sizeof(ZOLTAN_MAP));

  if (num_entries > 0){

    /* we create storage for the entries in advance */

    top = (ZOLTAN_ENTRY *)ZOLTAN_CALLOC(num_entries, sizeof(ZOLTAN_ENTRY));
    if (!top){
      ZOLTAN_FREE(&map);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Out of memory\n");
      return NULL;
    }
    if (store_keys) {
      keys = (char *)ZOLTAN_CALLOC(num_entries, num_bytes);
      if (!keys) {
        ZOLTAN_FREE(&top);
        ZOLTAN_FREE(&map);
	ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Out of memory\n");
	return NULL;
      }
    }
  }

  /* hash table */

  if (hash_range_max == 0){
    if (num_entries > 0){
      /* SRSR : This hash range maximum will be a prime number closer to
      num_entries for smaller problems and closer to num_entries/2 for larger
      problems. For very small problems hash_range_max can be larger than
      num_entries, but choosing hash_range_max = num_entries may lead to
      too many collisions. */
      my_range = (int) pow((double)num_entries, 0.90) ;
      hash_range_max = Zoltan_Recommended_Hash_Size(2*my_range) ;
    }
    else{
      /* This could be a performance bottleneck depending on num_entries.  */
      hash_range_max = 1000;
    }
  }

  entries = (ZOLTAN_ENTRY **)ZOLTAN_CALLOC(hash_range_max+1, sizeof(ZOLTAN_ENTRY*));

  if (!entries){
    ZOLTAN_FREE(&top);
    ZOLTAN_FREE(&map);
    ZOLTAN_FREE(&keys);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Out of memory\n");
    return NULL;
  }

  /* If the key is a multiple of ZOLTAN_ID_TYPEs, then we can use the key passed in
   * by the caller.  Otherwise we need to write it to a buffer that is a multiple
   * of ZOLTAN_ID_TYPEs, because Zoltan_Hash takes a tuple of ZOLTAN_ID_TYPEs.
   */

  if (num_bytes < (int)(sizeof(ZOLTAN_ID_TYPE))){
    map->num_zoltan_id_types = 1;
    map->zid = (ZOLTAN_ID_PTR)calloc(sizeof(ZOLTAN_ID_TYPE), 1);
  }
  else if (num_bytes == (int)(sizeof(ZOLTAN_ID_TYPE))){
    map->num_zoltan_id_types = 1;
    map->zid = NULL;
  }
  else {
    num = num_bytes / sizeof(ZOLTAN_ID_TYPE);
    rem = num_bytes % sizeof(ZOLTAN_ID_TYPE);

    if (rem == 0){
      map->num_zoltan_id_types = num;
      map->zid = NULL;
    }
    else{
      map->num_zoltan_id_types = num+1;
      map->zid = (ZOLTAN_ID_PTR)calloc(sizeof(ZOLTAN_ID_TYPE), num+1);
    }
  }


  map->entries     = entries;
  map->top         = top;
  map->key_size     = num_bytes;
  map->max_index   = hash_range_max;
  map->max_entries = num_entries;
  map->prev_index      = -1;
  map->prev_hash_index = -1;
  map->prev            = NULL;
  map->used            = 1;
  map->dynamicEntries  = (num_entries == 0);
  map->copyKeys        = store_keys;
  map->entry_count     = 0;
  map->keys            = keys;

  return map;
}

/*****************************************************************
 * Free the memory used by a map.  It is not an error to call this for
 * a map that has already been destroyed, or for a map that was never created.
 * (Return ZOLTAN_OK, etc.)
 */
int Zoltan_Map_Destroy(ZZ *zz, ZOLTAN_MAP** map)
{
  int i;
  ZOLTAN_ENTRY *nextEntry, *tmpEntry;

  if (!map || !*map){
    return ZOLTAN_OK;  /* OK to call Destroy more than once */
  }

  if ((*map)->zid != NULL){
    free((*map)->zid);
  }

  if ((*map)->copyKeys){

    /* free our copies of the keys */

    if (!(*map)->dynamicEntries){
      ZOLTAN_FREE(&((*map)->keys));
    }
    else{
      for (i=0; i <= (*map)->max_index; i++){
	nextEntry = (*map)->entries[i];
	while (nextEntry){
	  tmpEntry = nextEntry->next;
	  ZOLTAN_FREE(&(nextEntry->key));
	  nextEntry = tmpEntry;
	}
      }
    }
  }

  /* free the map entries */

  if (!(*map)->dynamicEntries){
    ZOLTAN_FREE(&(*map)->entries);
    ZOLTAN_FREE(&(*map)->top);
  }
  else{
    for (i=0; i <= (*map)->max_index; i++){
      nextEntry = (*map)->entries[i];
      while (nextEntry){
	tmpEntry = nextEntry->next;
	ZOLTAN_FREE(&nextEntry);
	nextEntry = tmpEntry;
      }
    }
    ZOLTAN_FREE(&(*map)->entries);
  }

  ZOLTAN_FREE(map);
  return ZOLTAN_OK;
}

/*****************************************************************
 * Add a key/value pair to the map.  Duplicate keys are ignored.
 * (Return ZOLTAN_OK, etc.)
 */

int Zoltan_Map_Add(ZZ *zz, ZOLTAN_MAP* map, char *key, intptr_t data) {
  return (Zoltan_Map_Find_Add(zz, map, key, data, NULL));
}

int Zoltan_Map_Find_Add(ZZ *zz, ZOLTAN_MAP* map, char *key, intptr_t datain, intptr_t *dataout)
{
  char *yo = "Zoltan_Map_Add";
  int index, match, i;
  ZOLTAN_ENTRY *element;
  char *hash_key;

  if (!map){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Map specified does not exist\n");
    return ZOLTAN_FATAL;
  }

  hash_key = current_key(map, key);

  index = Zoltan_Hash((ZOLTAN_ID_PTR)hash_key, map->num_zoltan_id_types, map->max_index);

  /* If this key is not found in the map, then add it */

  element = map->entries[index];
  match = 0;

  while (element != NULL){
    match = key_match(map, element->key, key);
    if (match){
      break;
    }
    element = element->next;
  }

  if (!match){
    i = map->entry_count;

    if (!map->dynamicEntries && (i == map->max_entries)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Fixed size map overflow\n");
      return ZOLTAN_FATAL;
    }

    if (map->dynamicEntries){
      element = (ZOLTAN_ENTRY *)ZOLTAN_MALLOC(sizeof(ZOLTAN_ENTRY));
      if (!element)
	return ZOLTAN_MEMERR;
    }
    else{
      element = map->top + i;
    }

    if (map->copyKeys){
      if(map->dynamicEntries) {
	element->key = (char *)ZOLTAN_MALLOC(map->key_size);
	if (!element->key) {
          ZOLTAN_FREE(&element);
	  return ZOLTAN_MEMERR;
        }
      }
      else
	element->key = (char *)map->keys + (map->entry_count * map->key_size);

      memcpy(element->key, key, map->key_size);
    }
    else{
      element->key = key;
    }

    element->data = datain;
    element->next = map->entries[index];

    map->entries[index] = element;

    map->entry_count++;
  }
  if (dataout != NULL) {
    *dataout = element->data;
  }

  return ZOLTAN_OK;
}

/*****************************************************************
 * Find the key in the map.  If found, set data to the key value.  
 * If not found, the data is set to an invalid value.
 * (Return ZOLTAN_OK, etc.)
 */

int Zoltan_Map_Find(ZZ *zz, ZOLTAN_MAP* map, char *key, intptr_t *data)
{
  char *yo = "Zoltan_Map_Find";
  int index, match;
  ZOLTAN_ENTRY *element;
  char *hash_key;

  *data = ZOLTAN_NOT_FOUND;

  if (!map){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Map specified does not exist\n");
    return ZOLTAN_FATAL;
  }

  hash_key = current_key(map, key);

  index = Zoltan_Hash((ZOLTAN_ID_PTR)hash_key, map->num_zoltan_id_types, map->max_index);

  element = map->entries[index];
  match = 0;

  while (element != NULL){
    match = key_match(map, element->key, key);
    if (match){
      *data = element->data;
      break;
    }
    element = element->next;
  }
  return ZOLTAN_OK;
}

/*****************************************************************
 * Return the number of keys in the map.
 */

int Zoltan_Map_Size(ZZ *zz, ZOLTAN_MAP* map)
{
  if (!map){
    return 0;
  }

  return map->entry_count;
}

/*****************************************************************
 * Begin iterating through the key/value pairs.
 *
 * *key is set to point to the first key
 * *data is set to the value that is associated with the first key
 * (Return ZOLTAN_OK, etc)
 */

int Zoltan_Map_First(ZZ *zz, ZOLTAN_MAP* map, char **key, intptr_t *data)
{
  char *yo = "Zoltan_Map_First";
  ZOLTAN_ENTRY *entry = NULL;
  int i;

  *key = NULL;
  *data = ZOLTAN_NOT_FOUND;

  if (map){

    if (map->entry_count == 0){
      map->prev_index = -1;
      map->prev_hash_index = -1;
      map->prev = NULL;
    }
    else{
      if (!map->dynamicEntries){
	map->prev_index = 0;
	entry = map->top;
      }
      else{

	/* find the first entry in the map */

	for (i=0; i <= map->max_index; i++){
	  if (map->entries[i]){

	    map->prev_hash_index = i;
	    entry = map->prev = map->entries[i];
	    break;
	  }
	}
	if (!entry){
	  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Entry not found\n");
	  return ZOLTAN_FATAL;
	}
      }

      *key = entry->key;
      *data = entry->data;
    }
  }
  else {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid map\n");
    return ZOLTAN_FATAL;
  }

  return ZOLTAN_OK;
}

/***********************************************************************
 * Return the next key/data pair.
 *
 * *key is set to point to the next key
 * *data is set to the pointer that is associated with that key
 *
 * *key and *data are NULL if there are no more key/data pairs
 *
 * (Return ZOLTAN_OK, etc)
 */

int Zoltan_Map_Next(ZZ *zz, ZOLTAN_MAP* map, char **key, intptr_t *data)
{
  ZOLTAN_ENTRY *next = NULL;
  int i;

  *key = (char *)NULL;
  *data = (intptr_t)ZOLTAN_NOT_FOUND;

  if (!map){
    return ZOLTAN_OK;
  }

  if (!map->dynamicEntries){
    if ((map->prev_index + 1 >= map->entry_count)){
      return ZOLTAN_OK;
    }

    map->prev_index++;

    next = map->top + map->prev_index;
  }
  else{
    if (map->prev == NULL){
      return Zoltan_Map_First(zz, map, key, data);
    }
    if (map->prev->next != NULL){
      next = map->prev->next;
      map->prev = next;
    }
    else{
      for (i=map->prev_hash_index+1; i <= map->max_index; i++){
	if (map->entries[i]){
	  map->prev_hash_index = i;
	  next = map->prev = map->entries[i];
	  break;
	}
      }
      if (!next){
	return ZOLTAN_OK;   /* at the end of the entries */
      }
    }
  }

  *key = next->key;
  *data = next->data;

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
