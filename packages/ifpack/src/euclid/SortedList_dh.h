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

#ifndef SORTEDLIST_DH_H
#define SORTEDLIST_DH_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

/* for private use by mpi factorization algorithms */

#include "euclid_common.h"
#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct _srecord
  {
    int col;
    int level;
    double val;
    int next;
  } SRecord;


  extern void SortedList_dhCreate (SortedList_dh * sList);
  extern void SortedList_dhDestroy (SortedList_dh sList);
  extern void SortedList_dhInit (SortedList_dh sList, SubdomainGraph_dh sg);
  extern void SortedList_dhEnforceConstraint (SortedList_dh sList,
					      SubdomainGraph_dh sg);

  extern void SortedList_dhReset (SortedList_dh sList, int row);

  extern int SortedList_dhReadCount (SortedList_dh sList);
  /* returns number of records inserted since last reset */

  extern void SortedList_dhResetGetSmallest (SortedList_dh sList);
  /* resets index used for SortedList_dhGetSmallestLowerTri().
   */

  extern SRecord *SortedList_dhGetSmallest (SortedList_dh sList);
  /* returns record with smallest column value that hasn't been
     retrieved via this method since last call to SortedList_dhReset()
     or SortedList_dhResetGetSmallest().
     If all records have been retrieved, returns NULL.
   */

  extern SRecord *SortedList_dhGetSmallestLowerTri (SortedList_dh sList);
  /* returns record with smallest column value that hasn't been
     retrieved via this method since last call to reset.  
     Only returns records where SRecord sr.col < row (per Init).
     If all records have been retrieved, returns NULL.
   */

  extern void SortedList_dhInsert (SortedList_dh sList, SRecord * sr);
  /* unilateral insert (does not check to see if item is already
     in list); does not permute sr->col; used in numeric
     factorization routines.
   */

  extern void SortedList_dhInsertOrUpdateVal (SortedList_dh sList,
					      SRecord * sr);
  /* unilateral insert: does not check to see if already
     inserted; does not permute sr->col; used in numeric 
     factorization routines.
   */

  extern bool SortedList_dhPermuteAndInsert (SortedList_dh sList,
					     SRecord * sr, double thresh);
  /* permutes sr->col, and inserts record in sorted list.
     Note: the contents of the passed variable "sr" may be changed.
     Note: this performs sparsification 
   */


  extern void SortedList_dhInsertOrUpdate (SortedList_dh sList, SRecord * sr);
  /* if a record with identical sr->col was inserted, updates sr->level
     to smallest of the two values; otherwise, inserts the record.
     Unlike SortedList_dhPermuteAndInsert, does not permute sr->col.
     Note: the contents of the passed variable "sr" may be changed.
     Warning: do not call SortedList_dhGetSmallestLowerTri() again
     until reset is called.
   */

  extern SRecord *SortedList_dhFind (SortedList_dh sList, SRecord * sr);
  /* returns NULL if no record is found containing sr->col 
   */

  extern void SortedList_dhUpdateVal (SortedList_dh sList, SRecord * sr);
#ifdef __cplusplus
}
#endif
#endif
