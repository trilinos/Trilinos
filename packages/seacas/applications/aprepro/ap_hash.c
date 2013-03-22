/* 
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
/***
   NAME
     hash
   PURPOSE
     Symbol table manipulation for Aprepro
***/
#include <stdlib.h>
#include "my_aprepro.h"
#include "y.tab.h"
#include <sys/types.h>
#include <limits.h>

void dumpsym(int type, int doInternal);
void pstats(void);
unsigned hash(char *symbol);
extern double mean(void);
extern void newsample(int);
extern double deviation(void);

#define HASHSIZE 5939
unsigned hash (char *symbol)
{
  unsigned hashval;
  for (hashval = 0; *symbol != '\0'; symbol++)
    hashval = *symbol + 65599 * hashval;
  return (hashval % HASHSIZE);
}

static symrec *sym_table[HASHSIZE];

symrec *putsym (char *sym_name, int sym_type, int is_internal)
{
  symrec *ptr;
  unsigned hashval;

  ptr = (symrec *) malloc (sizeof (symrec));
  if (ptr == NULL)
    return NULL;
  NEWSTR (sym_name, ptr->name);
  ptr->type = sym_type;
  ptr->value.var = 0;
  ptr->isInternal = is_internal;
  
  hashval = hash (ptr->name);
  ptr->next = sym_table[hashval];
  sym_table[hashval] = ptr;
  return ptr;
}

symrec *getsym (char *sym_name)
{
  symrec *ptr;
  for (ptr = sym_table[hash (sym_name)]; ptr != (symrec *) 0; ptr = ptr->next)
    if (strcmp (ptr->name, sym_name) == 0)
      return ptr;
  return (symrec *) 0;
}

void dumpsym (int type, int doInternal)
{
  symrec *ptr;
  unsigned hashval;

  char comment = getsym("_C_")->value.svar[0];
  
  if (type == VAR || type == SVAR) {
    printf ("\n%c   Variable = Value\n", comment);
    for (hashval = 0; hashval < HASHSIZE; hashval++) {
      for (ptr = sym_table[hashval]; ptr != (symrec *) 0; ptr = ptr->next) {
	if ((doInternal && ptr->isInternal) || (!doInternal && !ptr->isInternal)) {
	  if (ptr->type == VAR)
	    printf ("%c  {%-10s\t= %.10g}\n", comment, ptr->name, ptr->value.var);
	  else if (ptr->type == IMMVAR)
	    printf ("%c  {%-10s\t= %.10g} (immutable)\n", comment, ptr->name, ptr->value.var);
	  else if (ptr->type == SVAR)
	    printf ("%c  {%-10s\t= \"%s\"}\n", comment, ptr->name, ptr->value.svar);
	  else if (ptr->type == IMMSVAR)
	    printf ("%c  {%-10s\t= \"%s\"} (immutable)\n", comment, ptr->name, ptr->value.svar);
	}
      }
    }
  }
  else if (type == FNCT || type == SFNCT) {
    printf ("\nFunctions returning double:\n");
    for (hashval = 0; hashval < HASHSIZE; hashval++) {
      for (ptr = sym_table[hashval]; ptr != (symrec *) 0; ptr = ptr->next) {
	if (ptr->type == FNCT) {
	  printf ("%-20s:  %s\n", ptr->syntax, ptr->info);
	}
      }
    }
    printf ("\nFunctions returning string:\n");
    for (hashval = 0; hashval < HASHSIZE; hashval++) {
      for (ptr = sym_table[hashval]; ptr != (symrec *) 0; ptr = ptr->next) {
	if (ptr->type == SFNCT) {
	  printf ("%-20s:  %s\n", ptr->syntax, ptr->info);
	}
      }
    }
  }
}

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

#define MAXLEN 16
void pstats (void)
{
  symrec *ptr;
  unsigned hashval;
  unsigned entries = 0;
  int i;
  int chain_len;
  int maxlen = 0;
  int minlen = INT_MAX;
  int lengths[MAXLEN];
  int longer = 0;
  double hash_ratio = 0.0;

  memset ((void *) lengths, 0, sizeof (lengths));

  for (hashval = 0; hashval < HASHSIZE; hashval++)
    {
      chain_len = 0;
      for (ptr = sym_table[hashval]; ptr != (symrec *) 0; ptr = ptr->next)
	chain_len++;

      hash_ratio += chain_len * (chain_len + 1.0);
      entries += chain_len;
      if (chain_len >= MAXLEN)
	++longer;
      else
	++lengths[chain_len];

      minlen = min (minlen, chain_len);
      maxlen = max (maxlen, chain_len);

      if (chain_len > 0)
	newsample (chain_len);
    }

  hash_ratio = hash_ratio / ((float) entries / HASHSIZE *
			     (float) (entries + 2.0 * HASHSIZE - 1.0));
  printf ("%d entries in %d element hash table, ", (int)entries, HASHSIZE);
  printf ("%d (%1.0f%%) empty.\n",
	  lengths[0], ((double) lengths[0] / HASHSIZE) * 100.0);
  printf ("Hash ratio = %f\n", hash_ratio);
  printf ("Mean (nonempty) chain length = %4.1f, max = %d, min = %d, deviation = %4.2f\n",
	  mean (), maxlen, minlen, deviation ());

  for (i = 0; i < MAXLEN; i++)
    if (lengths[i])
      printf ("%3d chains of length %d\n", lengths[i], i);

  if (longer)
    printf ("%3d chains of length %d or longer\n", longer, MAXLEN);
}
