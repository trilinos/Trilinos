/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile $
 *    $Author $
 *    $Date $
 *    $Revision $
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdlib.h>
#include <zz_util_const.h>
#include <zoltan_mem.h>


static int key_match(int key_size, int *key1, int *key2);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 * A Zoltan_Map is like a C++ STL map.  It uses Zoltan_Hash to maintain a
 * a table mapping integer tuples (keys) to pointers (values).  Keys are unique.
 *
 * It supports:
 *
 *  Zoltan_Map_Create       Initialize the map
 *  Zoltan_Map_Add          Add an element to the map
 *  Zoltan_Map_Find         Find an element in the map, return the associated pointer
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
 */
ZOLTAN_MAP* Zoltan_Map_Create(ZZ *zz,     /* just need this for error messages */
		      int hash_range_max, /* > 0: maximum hash value */
					  /*   0: Zoltan_Map_Create will choose value */
		      int num_id_entries, /* length of a key */
		      int store_keys,     /* 1 - keep a copy of each key */
					  /* 0 - keep a pointer to caller's key */
		      int num_entries)    /* > 0: number of keys that will be added */
					  /*   0: dynamically allocate entries */
{
  char *yo = "Zoltan_Map_Create";
  ZOLTAN_ENTRY *top = NULL;
  ZOLTAN_ENTRY **entries = NULL;
  ZOLTAN_MAP* map = NULL;
  int *keys = NULL;
  int my_range ;

  if ((num_id_entries < 1) || (num_entries < 0) || (hash_range_max < 0)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Bad parameters\n");
    return NULL;
  }

  map = (ZOLTAN_MAP*) ZOLTAN_MALLOC(sizeof(ZOLTAN_MAP));

  if (num_entries > 0){

    /* we create storage for the entries in advance */

    top = (ZOLTAN_ENTRY *)ZOLTAN_CALLOC(num_entries, sizeof(ZOLTAN_ENTRY));
    if (!top){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Out of memory\n");
      return NULL;
    }
    if (store_keys) {
      keys = (int*)ZOLTAN_CALLOC(num_entries, sizeof(int)*num_id_entries);
      if (!keys) {
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
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Out of memory\n");
    return NULL;
  }


  map->entries     = entries;
  map->top         = top;
  map->id_size     = num_id_entries;
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

int Zoltan_Map_Add(ZZ *zz, ZOLTAN_MAP* map, int *key, int data) {
  return (Zoltan_Map_Find_Add(zz, map, key, data, NULL));
}

int Zoltan_Map_Find_Add(ZZ *zz, ZOLTAN_MAP* map, int *key, int datain, int *dataout)
{
  char *yo = "Zoltan_Map_Add";
  int index, match, i;
  ZOLTAN_ENTRY *element;
  ZOLTAN_ID_PTR zkey = (ZOLTAN_ID_PTR)key;

  if (!map){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Map specified does not exist\n");
    return ZOLTAN_FATAL;
  }

  index = Zoltan_Hash(zkey, map->id_size, map->max_index);

  /* If this key is not found in the map, then add it */

  element = map->entries[index];
  match = 0;

  while (element != NULL){
    match = key_match(map->id_size, element->key, key);
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
	element->key = (int *)ZOLTAN_MALLOC(sizeof(int) * map->id_size);
	if (!element->key)
	  return ZOLTAN_MEMERR;
      }
      else
	element->key = map->keys + map->entry_count;
      memcpy(element->key, key, sizeof(int)*map->id_size);
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
 * Find the key in the map.  If found, the data pointer is set to
 * the key value.  If not found, the data pointer is set to ZOLTAN_NOT_FOUND.
 * (Return ZOLTAN_OK, etc.)
 */

int Zoltan_Map_Find(ZZ *zz, ZOLTAN_MAP* map, int *key, int *data)
{
  char *yo = "Zoltan_Map_Find";
  int index, match;
  ZOLTAN_ENTRY *element;
  ZOLTAN_ID_PTR zkey = (ZOLTAN_ID_PTR)key;

  *data = ZOLTAN_NOT_FOUND;

  if (!map){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Map specified does not exist\n");
    return ZOLTAN_FATAL;
  }

  index = Zoltan_Hash(zkey, map->id_size, map->max_index);

  element = map->entries[index];
  match = 0;

  while (element != NULL){
    match = key_match(map->id_size, element->key, key);
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
 * *data is set to the integer that is associated with the first key
 * (Return ZOLTAN_OK, etc)
 */

int Zoltan_Map_First(ZZ *zz, ZOLTAN_MAP* map, int **key, int *data)
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
 * *data is set to the integer that is associated with that key
 *
 * *key is NULL if there are no more key/data pairs; *data is ZOLTAN_NOT_FOUND
 *
 * (Return ZOLTAN_OK, etc)
 */

int Zoltan_Map_Next(ZZ *zz, ZOLTAN_MAP* map, int **key, int *data)
{
  ZOLTAN_ENTRY *next = NULL;
  int i;

  *key = NULL;
  *data = ZOLTAN_NOT_FOUND;

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

/***********************************************************************
 * static helper functions
 */

static int key_match(int key_size, int *key1, int *key2)
{
  int i;

  for (i=0; i < key_size; i++){
    if (key1[i] != key2[i])
      return 0;           /* no match */
  }
  return 1;               /* it's a match */
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
