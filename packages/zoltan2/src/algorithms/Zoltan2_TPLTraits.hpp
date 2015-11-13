// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef _ZOLTAN2_TPLTRAITS_HPP_
#define _ZOLTAN2_TPLTRAITS_HPP_

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_as.hpp>
#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Environment.hpp>

#include <zoltan_types.h>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_TPLTraits.hpp
//! \brief Traits class to handle conversions between gno_t/lno_t 
//! and TPL data types (e.g., ParMETIS's idx_t, Scotch's SCOTCH_NUM,
//! Zoltan's ZOLTAN_ID_PTR)

namespace Zoltan2 {

//////////////////
// General case //
//////////////////

template <typename tpl_t, typename zno_t>
struct TPL_Traits {

  static inline bool OK_TO_CAST_TPL_T() 
  {
    // Return true if pointer to tpl_t can be safely used as pointer to zno_t
    return ((sizeof(tpl_t) == sizeof(zno_t)) && 
            (std::numeric_limits<tpl_t>::is_signed == 
             std::numeric_limits<zno_t>::is_signed));
  }

  static inline void ASSIGN_TPL_T(tpl_t &a, zno_t b)
  {
    // Assign a = b; make sure tpl_t is large enough to accept zno_t.
    try {
      a = Teuchos::asSafe<tpl_t, zno_t>(b);
    }
    catch (std::exception &e) {
      throw std::runtime_error(
       "TPL_Traits:  Value too large for TPL index type. "
       "Rebuild TPL with larger index type or rebuild without the TPL.");
    }
  }

  static inline void ASSIGN_TPL_T_ARRAY(tpl_t **a, ArrayView<zno_t> &b)
  {
    // Allocate array a; copy b values into a.
    size_t size = b.size();
    if (size > 0) {
      *a = new tpl_t[size];
      for (size_t i = 0; i < size; i++) ASSIGN_TPL_T((*a)[i], b[i]);
    }
    else {
      *a = NULL;
      // Note:  the Scotch manual says that if any rank has a non-NULL array,
      //        every process must have a non-NULL array.  In practice, 
      //        however, this condition is not needed for the arrays we use.
      //        For now, we'll set these arrays to NULL.  We could allocate
      //        a dummy value here if needed.  KDD 1/23/14
      // Note:  ParMETIS would likely prefer a dummy value as well.  It does
      //        not like NULL adjcny array.  KDD 10/7/14
    }
  }

  static inline void DELETE_TPL_T_ARRAY(tpl_t **a)
  {
    // Delete the copy made in ASSIGN_TPL_T_ARRAY.
    delete [] *a;
  }
};

///////////////////////////////////
// Special case:  zno_t == tpl_t //
// No error checking or copies   //
///////////////////////////////////

template <typename tpl_t>
struct TPL_Traits<tpl_t, tpl_t> {

  static inline bool OK_TO_CAST_TPL_T() {return true;}

  static inline void ASSIGN_TPL_T(tpl_t &a, tpl_t b)
  { a = b; }

  static inline void ASSIGN_TPL_T_ARRAY(tpl_t **a, ArrayView<tpl_t> &b)
  {
    if (b.size() > 0)
      *a = const_cast<tpl_t *> (b.getRawPtr());
    else
      *a = NULL;
      // Note:  the Scotch manual says that if any rank has a non-NULL array,
      //        every process must have a non-NULL array.  In practice, 
      //        however, this condition is not needed for the arrays we use.
      //        For now, we'll set these arrays to NULL.  We could allocate
      //        a dummy value here if needed.  KDD 1/23/14
  }
  static inline void DELETE_TPL_T_ARRAY(tpl_t **a) { }
};

////////////////////////////////////////////////////////////
// Special case:  Zoltan ZOLTAN_ID_PTR 
// Copy the data bitwise into the array of ZOLTAN_ID_TYPE   
// We assume that any memory pointed to by ZOLTAN_ID_PTR is
// big enough to hold the bits of zno_t -- that is, that  
// the ZOLTAN_ID_PTR's memory correctly accomodates Zoltan's
// num_gid_entries or num_lid_entries                     
////////////////////////////////////////////////////////////

template <typename zno_t>
struct TPL_Traits<ZOLTAN_ID_PTR, zno_t> {

  static const int NUM_ID = ((sizeof(zno_t) / sizeof(ZOLTAN_ID_TYPE) > 0)
                           ? (sizeof(zno_t) / sizeof(ZOLTAN_ID_TYPE))
                           : 1);

  static inline bool OK_TO_CAST_TPL_T() 
  {
    // There may be cases where if it OK to cast a pointer to a
    // zno_t to a ZOLTAN_ID_PTR, but the semantics of this 
    // function ask whether a pointer to a zno_t can be cast
    // to a pointer to a ZOLTAN_ID_PTR.  Thus, the answer is "no."
    return false;
  }

  static inline void ASSIGN_TPL_T(ZOLTAN_ID_PTR &a, zno_t b)
  {
    switch (NUM_ID) {
    case 1:
      a[0] = static_cast<ZOLTAN_ID_TYPE>(b);
      break;
    case 2: {
      ZOLTAN_ID_TYPE *ptr = (ZOLTAN_ID_TYPE *)(&b);
      a[0] = ptr[0]; 
      a[1] = ptr[1];
      break;
    }
    default: {
      ZOLTAN_ID_TYPE *ptr = (ZOLTAN_ID_TYPE *)(&b);
      for (int i = 0; i < NUM_ID; i++) a[i] = ptr[i];
    }
    }
  }

  static inline void ASSIGN_TPL_T_ARRAY(ZOLTAN_ID_PTR *a, ArrayView<zno_t> &b)
  {
    // Allocate array a; copy b values into a.
    size_t size = b.size();
    if (size > 0) {
      if (NUM_ID == 1) {
        *a = b.getRawPtr();  // Don't have to make a new copy
      }
      else {
        *a = new ZOLTAN_ID_TYPE[size*NUM_ID];
        for (size_t i = 0; i < size; i++) ASSIGN_TPL_T((*a)[i], b[i]);
      }
    }
    else {
      *a = NULL;
    }
  }

  static inline void DELETE_TPL_T_ARRAY(ZOLTAN_ID_PTR *a)
  {
    // Delete the copy made in ASSIGN_TPL_T_ARRAY.
    if (NUM_ID != 1)
      delete [] *a;
  }
};
} // namespace Zoltan2

#endif
