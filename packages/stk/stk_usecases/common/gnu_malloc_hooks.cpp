/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifdef USE_GNU_MALLOC_HOOKS

#include <malloc.h>

#include <cstdlib>
#include <cmath>
#include <map>

//if USE_GNU_MALLOC_HOOKS is not defined, then
//the compiler won't see the rest of this file.

enum { ALLOC_BYTES=0,
       FREED_BYTES=1,
       ALLOC_BLKS=2,
       FREED_BLKS=3,
       NUM_STATS=4};

static bool gnu_malloc_hooks_disabled = false;
static double my_stats[NUM_STATS];

typedef std::map<void*,unsigned> size_map;

size_map* my_sizes = NULL;

const unsigned ONEMB = 1024*1024;

/* Prototypes for our hooks.  */
static void my_init_hook (void);
static void *my_malloc_hook (size_t, const void *);
static void *(*old_malloc_hook)(size_t, const void*);
static void *my_realloc_hook (void*, size_t, const void *);
static void *(*old_realloc_hook)(void*, size_t, const void*);
static void my_free_hook (void*, const void *);
static void (*old_free_hook)(void*, const void*);

/* Override initializing hook from the C library. */
void (*__malloc_initialize_hook) (void) = my_init_hook;

void
my_init_hook (void)
{
  my_sizes = new size_map;
  old_malloc_hook = __malloc_hook;
  old_realloc_hook = __realloc_hook;
  old_free_hook = __free_hook;
  __malloc_hook = my_malloc_hook;
  __realloc_hook = my_realloc_hook;
  __free_hook = my_free_hook;
  for(unsigned i=0; i<NUM_STATS; ++i) my_stats[i] = 0.0;
}

void *
my_malloc_hook (size_t size, const void *caller)
{
  void *result;
  /* Restore all old hooks */
  __malloc_hook = old_malloc_hook;
  __free_hook = old_free_hook;

  /* Call recursively */
  result = malloc (size);
  my_stats[ALLOC_BYTES] += size;
  my_stats[ALLOC_BLKS] += 1;

  (*my_sizes)[result] = size;

  /* Restore our own hooks */
  __malloc_hook = my_malloc_hook;
  __free_hook = my_free_hook;
  return result;
}

void *
my_realloc_hook (void* ptr, size_t size, const void *caller)
{
  void *result;
  /* Restore all old hooks */
  __malloc_hook = old_malloc_hook;
  __realloc_hook = old_realloc_hook;
  __free_hook = old_free_hook;

  size_map::iterator iter = my_sizes->find(ptr);
  if (iter != my_sizes->end()) {
    unsigned oldsize = iter->second;
    my_stats[FREED_BYTES] += oldsize;
    my_stats[FREED_BLKS] += 1;
    my_sizes->erase(iter);
  }

  /* Call recursively */
  result = realloc(ptr, size);
  my_stats[ALLOC_BYTES] += size;
  my_stats[ALLOC_BLKS] += 1;
  (*my_sizes)[result] = size;

  /* Restore our own hooks */
  __malloc_hook = my_malloc_hook;
  __realloc_hook = my_realloc_hook;
  __free_hook = my_free_hook;
  return result;
}

void
my_free_hook (void *ptr, const void *caller)
{
  /* Restore all old hooks */
  __malloc_hook = old_malloc_hook;
  __free_hook = old_free_hook;

  size_map::iterator iter = my_sizes->find(ptr);
  if (iter != my_sizes->end()) {
    unsigned size = iter->second;
    my_stats[FREED_BYTES] += size;
    my_stats[FREED_BLKS] += 1;
    my_sizes->erase(iter);
  }
 
  /* Call recursively */
  free(ptr);

  /* Restore our own hooks */
  __malloc_hook = my_malloc_hook;
  __free_hook = my_free_hook;
}

void disable_gnu_malloc_hooks()
{
  /* Restore all old hooks */
  __malloc_hook = old_malloc_hook;
  __free_hook = old_free_hook;
  gnu_malloc_hooks_disabled = true;
}


void reset_malloc_stats()
{
  for(unsigned i=0; i<NUM_STATS; ++i) my_stats[i] = 0.0;
}

double alloc_MB()
{
  return my_stats[ALLOC_BYTES]/ONEMB;
}

unsigned alloc_blks()
{
  return (int)std::floor(my_stats[ALLOC_BLKS]);
}

double freed_MB()
{
  return my_stats[FREED_BYTES]/ONEMB;
}

unsigned freed_blks()
{
  return(int)std::floor(my_stats[FREED_BLKS]);
}

#endif

//USE_GNU_MALLOC_HOOKS

