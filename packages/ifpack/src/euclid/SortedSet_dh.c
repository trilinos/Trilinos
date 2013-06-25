/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "SortedSet_dh.h"
#include "shellSort_dh.h"
#include "Mem_dh.h"


#undef __FUNC__
#define __FUNC__ "SortedSet_dhCreate"
void
SortedSet_dhCreate (SortedSet_dh * ss, int size)
{
  START_FUNC_DH
    struct _sortedset_dh *tmp =
    (struct _sortedset_dh *) MALLOC_DH (sizeof (struct _sortedset_dh));
  CHECK_V_ERROR;
  *ss = tmp;

  tmp->n = size;
  tmp->list = (int *) MALLOC_DH (size * sizeof (int));
  CHECK_V_ERROR;
  tmp->count = 0;
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "SortedSet_dhDestroy"
void
SortedSet_dhDestroy (SortedSet_dh ss)
{
  START_FUNC_DH if (ss->list != NULL)
    {
      FREE_DH (ss->list);
      CHECK_V_ERROR;
    }
  FREE_DH (ss);
  CHECK_V_ERROR;
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "SortedSet_dhInsert"
void
SortedSet_dhInsert (SortedSet_dh ss, int idx)
{
  START_FUNC_DH bool isInserted = false;
  int ct = ss->count;
  int *list = ss->list;
  int i, n = ss->n;

  /* determine if item was already inserted */
  for (i = 0; i < ct; ++i)
    {
      if (list[i] == idx)
	{
	  isInserted = true;
	  break;
	}
    }

  /* is we need to insert the item, first check for overflow
     and reallocate if necessary, then append the index to the
     end of the list.
   */
  if (!isInserted)
    {
      if (ct == n)
	{
	  int *tmp = (int *) MALLOC_DH (n * 2 * sizeof (int));
	  CHECK_V_ERROR;
	  memcpy (tmp, list, n * sizeof (int));
	  FREE_DH (list);
	  CHECK_V_ERROR;
	  list = ss->list = tmp;
	  ss->n *= 2;
	}

      list[ct] = idx;
      ss->count += 1;
    }
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "SortedSet_dhGetList"
void
SortedSet_dhGetList (SortedSet_dh ss, int **list, int *count)
{
  START_FUNC_DH shellSort_int (ss->count, ss->list);
  *list = ss->list;
  *count = ss->count;
END_FUNC_DH}
