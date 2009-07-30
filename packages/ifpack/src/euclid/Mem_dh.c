/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Parser_dh.h"
#include "Mem_dh.h"

/* TODO: error checking is not complete; memRecord_dh  need to
         be done in Mem_dhMalloc() and Mem_dhFree()l
*/


  /* a memRecord_dh is pre and post-pended to every 
   * piece of memory obtained by calling MALLOC_DH 
   */
typedef struct {
    double size;
    double cookie;
} memRecord_dh;

struct _mem_dh {
  double maxMem;        /* max allocated at any point in time */
  double curMem;        /* total currently allocated */
  double totalMem;      /* total cumulative malloced */
  double mallocCount;  /* number of times mem_dh->malloc has been called. */
  double freeCount;    /* number of times mem_dh->free has been called. */
};


#undef __FUNC__
#define __FUNC__ "Mem_dhCreate"
void Mem_dhCreate(Mem_dh *m)
{
  START_FUNC_DH
  struct _mem_dh *tmp = (struct _mem_dh*)PRIVATE_MALLOC(sizeof(struct _mem_dh)); CHECK_V_ERROR;
  *m = tmp;
  tmp->maxMem = 0.0;
  tmp->curMem = 0.0;
  tmp->totalMem = 0.0;
  tmp->mallocCount = 0.0;
  tmp->freeCount = 0.0;
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "Mem_dhDestroy"
void Mem_dhDestroy(Mem_dh m)
{
  START_FUNC_DH
  if (Parser_dhHasSwitch(parser_dh, "-eu_mem")) {
    Mem_dhPrint(m, stdout, false); CHECK_V_ERROR;
  }

  PRIVATE_FREE(m);
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "Mem_dhMalloc"
void* Mem_dhMalloc(Mem_dh m, size_t size)
{
  START_FUNC_DH_2
  void *retval;
  memRecord_dh *tmp;
  size_t s = size + 2*sizeof(memRecord_dh);
  void *address;

  address = PRIVATE_MALLOC(s);

  if (address == NULL) {
    sprintf(msgBuf_dh, "PRIVATE_MALLOC failed; totalMem = %g; requested additional = %i", m->totalMem, (int)s);
    SET_ERROR(NULL, msgBuf_dh);
  }

  retval = (char*)address + sizeof(memRecord_dh);

  /* we prepend and postpend a private record to the 
   * requested chunk of memory; this permits tracking the
   * sizes of freed memory, along with other rudimentary
   * error checking.  This is modeled after the PETSc code.
   */
  tmp = (memRecord_dh*)address;
  tmp->size = (double) s;

  m->mallocCount += 1;
  m->totalMem += (double)s;
  m->curMem += (double)s;
  m->maxMem = MAX(m->maxMem, m->curMem);

  END_FUNC_VAL_2( retval )
}


#undef __FUNC__
#define __FUNC__ "Mem_dhFree"
void Mem_dhFree(Mem_dh m, void *ptr)
{
  START_FUNC_DH_2
  double size;
  char *tmp = (char*)ptr;
  memRecord_dh *rec;
  tmp -= sizeof(memRecord_dh);
  rec = (memRecord_dh*)tmp;
  size = rec->size;   

  mem_dh->curMem -= size;
  mem_dh->freeCount += 1;

  PRIVATE_FREE(tmp);
  END_FUNC_DH_2
}


#undef __FUNC__
#define __FUNC__ "Mem_dhPrint"
void  Mem_dhPrint(Mem_dh m, FILE* fp, bool allPrint)
{
  START_FUNC_DH_2
  if (fp == NULL) SET_V_ERROR("fp == NULL");
  if (myid_dh == 0 || allPrint) {
    double tmp;
    fprintf(fp, "---------------------- Euclid memory report (start)\n");
    fprintf(fp, "malloc calls = %g\n", m->mallocCount);
    fprintf(fp, "free   calls = %g\n", m->freeCount);
    fprintf(fp, "curMem          = %g Mbytes (should be zero)\n", 
                                                   m->curMem/1000000);
    tmp = m->totalMem / 1000000;
    fprintf(fp, "total allocated = %g Mbytes\n", tmp);
    fprintf(fp, "max malloc      = %g Mbytes (max allocated at any point in time)\n", m->maxMem/1000000);
    fprintf(fp, "\n");
    fprintf(fp, "---------------------- Euclid memory report (end)\n");
  } 
  END_FUNC_DH_2
}
