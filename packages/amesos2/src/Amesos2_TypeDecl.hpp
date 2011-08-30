// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER


/**
 * \file   Amesos2_TypeDecl.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Mon Jul 18 11:57:14 2011
 * 
 * \brief  Enum and other types declarations for Amesos2
 * 
 * 
 */

#ifndef AMESOS2_TYPEDECL_HPP
#define AMESOS2_TYPEDECL_HPP

namespace Amesos2 {

  /**
   * \brief Used to indicate a phase in the direct solution.
   *
   * \ingroup amesos2_enums
   */
  typedef enum {
    CLEAN,
    PREORDERING,
    SYMBFACT,
    NUMFACT,
    SOLVE
  } EPhase;

  /** \internal
   *
   * \brief Indicates that the concrete class has a special
   * implementation that should be called.
   *
   * Matrix Adapters \c typedef this as \c {get_ccs|get_crs}_spec
   * indicating that the concrete adapter has a special
   * implementation for either the \c getCrs or \c getCcs functions.
   */
  struct has_special_impl {};

  /** \internal
   *
   * \brief Indicates that the concrete class can use the generic
   * getC{c|r}s methods implemented in MatrixAdapter.
   */
  struct no_special_impl {};

  /** \internal
   * 
   * \brief Indicates that the object of an adapter provides row
   * access to its data.
   */
  struct row_access {};
    
  /** \internal
   *
   * \brief Indicates that the object of an adapter provides column
   * access to its data.
   */
  struct col_access {};

  /** \internal
   *
   * \enum EDistribution
   *
   * An enum of this type is expected by the Matrix adapters' getCrs
   * and getCcs functions to describe the layout of the
   * representation on the calling processors.
   *
   * \ingroup amesos2_enums
   */
  typedef enum {
    DISTRIBUTED,                /**< no processor has a view of the entire matrix, only local pieces */
    DISTRIBUTED_NO_OVERLAP,     /**< no row or column may be present on more than one processor */
    GLOBALLY_REPLICATED,        /**< each processor has a view of the entire matrix */
    ROOTED                      /**< only \c rank=0 has a full view, all others have nothing. */
  } EDistribution;

  /** \internal
   *
   * \enum EStorage_Ordering
   *
   * This enum also used by the matrix adapters to indicate whether
   * the indices of the representation must be in sorted order or
   * can have an arbitrary order.
   *
   * \ingroup amesos2_enums
   */
  typedef enum {
    SORTED_INDICES,             /**< row/col indices need to appear in sorted order */
    ARBITRARY                   /**< index order can be arbitrary */
  } EStorage_Ordering;

}

#endif	// AMESOS2_TYPEDECL_HPP
