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

struct Zoltan_Map_Entry{
  int *key;          /* a copy of or a pointer to callers key */
  void *data;        /* a pointer provided by caller */
  struct Zoltan_Map_Entry *next;
};

typedef struct Zoltan_Map_Entry ZOLTAN_ENTRY;

struct Zoltan_Map_List{
  ZOLTAN_ENTRY **entries; /* hash array, length max_index + 1 */

  ZOLTAN_ENTRY *top;      /* if dynamicEntries==0, entries are here */

  int id_size;          /* size of integer tuple */
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
};

typedef struct Zoltan_Map_List ZOLTAN_MAP;

static ZOLTAN_MAP map_array[ZOLTAN_MAX_MAP];

static int key_match(int key_size, int *key1, int *key2);

static ZOLTAN_MAP *get_used_map(int map_num);

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
 */

/*****************************************************************
 * Create a new map.  
 * Return the map number to be used in subsequent calls.
 * (Return -1 on error)
 */
int Zoltan_Map_Create(ZZ *zz,     /* just need this for error messages */
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
  int map_num = 0;

  if ((num_id_entries < 1) || (num_entries < 0) || (hash_range_max < 0)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Bad parameters\n");
    return -1;
  }

  while ((map_num < ZOLTAN_MAX_MAP) && (map_array[map_num].used != 0)) map_num++;

  if (map_num == ZOLTAN_MAX_MAP){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Maximum number of maps already allocated\n");
    return -1;
  }

  if (num_entries > 0){

    /* we create storage for the entries in advance */

    top = (ZOLTAN_ENTRY *)ZOLTAN_CALLOC(num_entries, sizeof(ZOLTAN_ENTRY));
    if (!top){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Out of memory\n");
      return -1;
    }
  }

  /* hash table */

  if (hash_range_max == 0){  /* we should give this more thought */
    if (num_entries > 0){
      if (num_entries > 10000)
        hash_range_max = 10000;
      else
        hash_range_max = num_entries;
    }
    else{
      hash_range_max = 1000;
    }
  }

  entries = (ZOLTAN_ENTRY **)ZOLTAN_CALLOC(hash_range_max+1, sizeof(ZOLTAN_ENTRY*));

  if (!entries){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Out of memory\n");
    return -1;
  }

  map_array[map_num].entries     = entries;
  map_array[map_num].top         = top;
  map_array[map_num].id_size     = num_id_entries;
  map_array[map_num].max_index   = hash_range_max;
  map_array[map_num].max_entries = num_entries;
  map_array[map_num].prev_index      = -1;
  map_array[map_num].prev_hash_index = -1;
  map_array[map_num].prev            = NULL;
  map_array[map_num].used            = 1;
  map_array[map_num].dynamicEntries  = (num_entries == 0); 
  map_array[map_num].copyKeys        = store_keys;
  map_array[map_num].entry_count     = 0;

  return map_num;
}

/*****************************************************************
 * Free the memory used by a map.  It is not an error to call this for
 * a map that has already been destroyed, or for a map that was never created.
 * (Return ZOLTAN_OK, etc.)
 */
int Zoltan_Map_Destroy(ZZ *zz, int map_num)
{
  int i;
  ZOLTAN_ENTRY *nextEntry, *tmpEntry;
  ZOLTAN_MAP *map = get_used_map(map_num);

  if (!map){
    return ZOLTAN_OK;  /* OK to call Destroy more than once */
  }

  if (map->copyKeys){

    /* free our copies of the keys */

    if (!map->dynamicEntries){
      for (i=0; i < map->entry_count; i++){
        ZOLTAN_FREE(&(map->top[i].key));
      }
    }
    else{
      for (i=0; i <= map->max_index; i++){
        nextEntry = map->entries[i];
        while (nextEntry){
          tmpEntry = nextEntry->next;
          ZOLTAN_FREE(&(nextEntry->key));
          nextEntry = tmpEntry;
        }
      }
    }
  }

  /* free the map entries */

  if (!map->dynamicEntries){
    ZOLTAN_FREE(&map->entries);
    ZOLTAN_FREE(&map->top);
  }
  else{
    for (i=0; i <= map->max_index; i++){
      nextEntry = map->entries[i];
      while (nextEntry){
        tmpEntry = nextEntry->next;
        ZOLTAN_FREE(&nextEntry);
        nextEntry = tmpEntry;
      }
    }
    ZOLTAN_FREE(&map->entries);
  }

  memset(map, 0, sizeof(ZOLTAN_MAP));

  return ZOLTAN_OK;
}

/*****************************************************************
 * Add a key/value pair to the map.  Duplicate keys are ignored.
 * (Return ZOLTAN_OK, etc.)
 */

int Zoltan_Map_Add(ZZ *zz, int map_num, int *key, void *data)
{
  char *yo = "Zoltan_Map_Add";
  int index, match, i;
  ZOLTAN_ENTRY *element;
  ZOLTAN_ID_PTR zkey = (ZOLTAN_ID_PTR)key;

  ZOLTAN_MAP *map = get_used_map(map_num);

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
      element->key = (int *)ZOLTAN_MALLOC(sizeof(int) * map->id_size);
      if (!element->key) 
        return ZOLTAN_MEMERR;

      for (i=0; i < map->id_size; i++){
        element->key[i] = key[i];
      }
    }
    else{
      element->key = key;
    }

    element->data = data;
    element->next = map->entries[index];

    map->entries[index] = element;

    map->entry_count++;
  }

  return ZOLTAN_OK;
}

/*****************************************************************
 * Find the key in the map.  If found, the data pointer is set to
 * the key value.  If not found, the data pointer is set to NULL.
 * (Return ZOLTAN_OK, etc.)
 */

int Zoltan_Map_Find(ZZ *zz, int map_num, int *key, void **data)
{
  char *yo = "Zoltan_Map_Find";
  int index, match;
  ZOLTAN_ENTRY *element;
  ZOLTAN_ID_PTR zkey = (ZOLTAN_ID_PTR)key;

  ZOLTAN_MAP *map = get_used_map(map_num);

  *data = NULL;

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

int Zoltan_Map_Size(ZZ *zz, int map_num)
{
  ZOLTAN_MAP *map = get_used_map(map_num);

  if (!map){
    return 0;
  }

  return map->entry_count;
}

/*****************************************************************
 * Begin iterating through the key/value pairs.
 *
 * *key is set to point to the first key
 * *data is set to the pointer that is associated with the first key
 * (Return ZOLTAN_OK, etc)
 */

int Zoltan_Map_First(ZZ *zz, int map_num, int **key, void **data)
{
  char *yo = "Zoltan_Map_First";
  ZOLTAN_MAP *map = NULL;
  ZOLTAN_ENTRY *entry = NULL;
  int i;

  *key = NULL;
  *data = NULL;

  map = get_used_map(map_num);

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

int Zoltan_Map_Next(ZZ *zz, int map_num, int **key, void **data)
{
  ZOLTAN_MAP *map = NULL;
  ZOLTAN_ENTRY *next = NULL; 
  int i;

  *key = NULL;
  *data = NULL;

  map = get_used_map(map_num);

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
      return Zoltan_Map_First(zz, map_num, key, data);
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

static ZOLTAN_MAP *get_used_map(int map_num)
{
  ZOLTAN_MAP *map = NULL;

  if ((map_num >= 0) &&
      (map_num < ZOLTAN_MAX_MAP) &&
      (map_array[map_num].used != 0)){
  
    map = map_array + map_num;
  }
  return map;
}

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

