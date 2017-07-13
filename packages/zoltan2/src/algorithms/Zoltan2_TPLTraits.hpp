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

/////////////////////////////////////////////////
// General case:  first_t and second_t differ  //
/////////////////////////////////////////////////

template <typename first_t, typename second_t>
struct TPL_Traits {

  static inline bool OK_TO_CAST() 
  {
    // Return true if pointer to first_t can be used as pointer to second_t
    return ((sizeof(first_t) == sizeof(second_t)) && 
            (std::numeric_limits<first_t>::is_signed == 
             std::numeric_limits<second_t>::is_signed));
  }

  static inline void ASSIGN(first_t &a, second_t b)
  {
    // Assign a = b; make sure first_t is large enough to accept second_t.
    try {
      a = Teuchos::asSafe<first_t, second_t>(b);
    }
    catch (std::exception &e) {
      throw std::runtime_error(
       "TPL_Traits:  Value too large for TPL index type. "
       "Rebuild TPL with larger index type or rebuild without the TPL.");
    }
  }

  static inline void ASSIGN_ARRAY(first_t **a, ArrayView<second_t> &b)
  {
    // Allocate array a; copy b values into a.
    size_t size = b.size();
    if (size > 0) {
      *a = new first_t[size];
      for (size_t i = 0; i < size; i++) ASSIGN((*a)[i], b[i]);
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

  static inline void SAVE_ARRAYRCP(ArrayRCP<first_t> *a, second_t *b, 
                                   size_t size)
  {
    // Allocate array tmp; copy size values of b into tmp; 
    // save tmp as ArrayRCP a
    if (size > 0) {
      first_t *tmp = new first_t[size];
      for (size_t i = 0; i < size; i++) ASSIGN(tmp[i], b[i]);
      *a = arcp(tmp, 0, size, true);
    }
    else {
      *a = Teuchos::null;
    }
  }

  static inline void DELETE_ARRAY(first_t **a)
  {
    // Delete the copy made in ASSIGN_ARRAY.
    delete [] *a;
  }
};

////////////////////////////////////////
// Special case:  second_t == first_t //
// No error checking or copies        //
////////////////////////////////////////

template <typename first_t>
struct TPL_Traits<first_t, first_t> {

  static inline bool OK_TO_CAST() {return true;}

  static inline void ASSIGN(first_t &a, first_t b) { a = b; }

  static inline void ASSIGN_ARRAY(first_t **a, ArrayView<first_t> &b)
  {
    if (b.size() > 0)
      *a = b.getRawPtr();
    else
      *a = NULL;
      // Note:  the Scotch manual says that if any rank has a non-NULL array,
      //        every process must have a non-NULL array.  In practice, 
      //        however, this condition is not needed for the arrays we use.
      //        For now, we'll set these arrays to NULL.  We could allocate
      //        a dummy value here if needed.  KDD 1/23/14
  }

  static inline void SAVE_ARRAYRCP(ArrayRCP<first_t> *a, first_t *b, 
                                         size_t size)
  {
    // Return in a an ArrayRCP of b.
    if (size > 0)
      *a = arcp(b, 0, size, true);
    else
      *a = Teuchos::null;
  }

  static inline void DELETE_ARRAY(first_t **a) { }
};

/////////////////////////////////////////////////////
// Special case:  first_t == Zoltan ZOLTAN_ID_PTR  //
/////////////////////////////////////////////////////

template <typename second_t>
struct TPL_Traits<ZOLTAN_ID_PTR, second_t> {

  // Copy the data bitwise INTO the array of ZOLTAN_ID_TYPE   
  // We assume that any memory pointed to by ZOLTAN_ID_PTR is
  // big enough to hold the bits of second_t -- that is, that  
  // the ZOLTAN_ID_PTR's memory correctly accomodates Zoltan's
  // num_gid_entries or num_lid_entries                     

  static const int NUM_ID = ((sizeof(second_t) / sizeof(ZOLTAN_ID_TYPE) > 0)
                           ? (sizeof(second_t) / sizeof(ZOLTAN_ID_TYPE))
                           : 1);

  static inline bool OK_TO_CAST() 
  {
    // There may be cases where if it OK to cast a pointer to a
    // second_t to a ZOLTAN_ID_PTR, but the semantics of this 
    // function ask whether a pointer to a second_t can be cast
    // to a pointer to a ZOLTAN_ID_PTR.  Thus, the answer is "no."
    return false;
  }

  static inline void ASSIGN(ZOLTAN_ID_PTR &a, second_t b)
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

  static inline void ASSIGN_ARRAY(ZOLTAN_ID_PTR *a, ArrayView<second_t> &b)
  {
    // Allocate array a; copy b values into a.
    size_t size = b.size();
    if (size > 0) {
      if (NUM_ID == 1) {
        // Don't have to make a new copy
        *a = reinterpret_cast<ZOLTAN_ID_PTR> (b.getRawPtr());  
      }
      else {
        *a = new ZOLTAN_ID_TYPE[size*NUM_ID];
        for (size_t i = 0; i < size; i++) {
          ZOLTAN_ID_PTR tmp = &((*a)[i*NUM_ID]);
          ASSIGN(tmp, b[i]);
        }
      }
    }
    else {
      *a = NULL;
    }
  }

  static inline void SAVE_ARRAYRCP(ArrayRCP<ZOLTAN_ID_PTR> *a, second_t *b, 
                                         size_t size)
  {
    throw std::runtime_error(
               "TPL_Traits::SAVE_ARRAYRCP<ZOLTAN_ID_PTR,second_t> "
               "is not implemented.");
  }

  static inline void DELETE_ARRAY(ZOLTAN_ID_PTR *a)
  {
    // Delete the copy made in ASSIGN_ARRAY.
    if (NUM_ID != 1)
      delete [] *a;
  }
};

/////////////////////////////////////////////////////
// Special case:  second_t == Zoltan ZOLTAN_ID_PTR  //
/////////////////////////////////////////////////////

template <typename first_t>
struct TPL_Traits<first_t, ZOLTAN_ID_PTR> {

  // Copy the data bitwise FROM the array of ZOLTAN_ID_TYPE   

  static const int NUM_ID = TPL_Traits<ZOLTAN_ID_PTR, first_t>::NUM_ID;

  static inline bool OK_TO_CAST()
  {
    // There may be cases where if it OK to cast a pointer to a
    // first_t to a ZOLTAN_ID_PTR, but the semantics of this
    // function ask whether a pointer to a first_t can be cast
    // to a pointer to a ZOLTAN_ID_PTR.  Thus, the answer is "no."
    return false;
  }

  static inline void ASSIGN(first_t &a, ZOLTAN_ID_PTR b)
  {
    switch (NUM_ID) {
    case 1:
      a = static_cast<first_t>(b[0]);
      break;
    default:
      first_t *tmp = (first_t *) b;
      a = *tmp;
    }
  }

  static inline void ASSIGN_ARRAY(first_t *a, ArrayView<ZOLTAN_ID_PTR> &b)
  {
    throw std::runtime_error(
               "TPL_Traits::ASSIGN_ARRAY<first_t,ZOLTAN_ID_PTR> "
               "is not implemented.");
  }

  static inline void SAVE_ARRAYRCP(ArrayRCP<first_t> *a, ZOLTAN_ID_PTR b, 
                                         size_t size)
  {
    // Here, size * NUM_ID == length of b; that is, size is the number of
    // objects in b.
    // Return in a an ArrayRCP of b.
    if (size > 0) {
      if (NUM_ID == 1)
        *a = arcp(b, 0, size, true);  // Don't have to make a new copy
      else {
        first_t *tmp = new first_t[size];
        for (size_t i = 0; i < size; i++) ASSIGN(tmp[i], &(b[i*NUM_ID]));
        *a = arcp(tmp, 0, size, true);
      }
    }
    else {
      *a = Teuchos::null;
    }
  }

  static inline void DELETE_ARRAY(first_t **a)
  {
    // Delete the copy made in ASSIGN_ARRAY.
    if (NUM_ID != 1)
      delete [] *a;
  }
};

} // namespace Zoltan2

#endif
