// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  struct has_special_impl {
    bool set;
  };

  /** \internal
   *
   * \brief Indicates that the concrete class can use the generic
   * getC{c|r}s methods implemented in MatrixAdapter.
   */
  struct no_special_impl {
    bool set;
  };

  /** \internal
   *
   * \brief Indicates that the object of an adapter provides row
   * access to its data.
   */
  struct row_access {
    bool set;
  };

  /** \internal
   *
   * \brief Indicates that the object of an adapter provides column
   * access to its data.
   */
  struct col_access {
    bool set;
  };

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
    ROOTED,                     /**< only \c rank=0 has a full view, all others have nothing. */
    CONTIGUOUS_AND_ROOTED       /**< the global GIDs of the map are not contiguous */
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
