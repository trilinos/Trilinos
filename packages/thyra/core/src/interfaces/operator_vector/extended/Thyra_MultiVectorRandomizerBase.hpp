// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_RANDOMIZER_BASE_HPP
#define THYRA_MULTI_VECTOR_RANDOMIZER_BASE_HPP

#include "Thyra_OperatorVectorTypes.hpp"

namespace Thyra {


/** \brief Base interface for a strategy object for randomizing a
 * multi-vector.
 *
 * This object is *not* stateless in its use!  Every time it generates a new
 * random multi-vector its behavior changes.
 *
 * A single <tt>MultiVectorRandomizerBase</tt> object may be compatible with
 * many different types of concrete vector space implementations or may
 * compatible with only a specific instantiation of a concrete vector space
 * subclass.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class MultiVectorRandomizerBase {
public:

  /** \brief . */
  virtual ~MultiVectorRandomizerBase() {}

  /** \brief Determines if <tt>*this</tt> is compatible with multi-vectors
   * from the <tt>VectorSpace</tt> <tt>space</tt>.
   */
  virtual bool isCompatible( const VectorSpaceBase<Scalar> &space ) const = 0;

  /** \brief Randomize a "compatible" multi-vector.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>mv!=NULL</tt>
   * <li><tt>this->isCompatible(*mv->range()) == true</tt>
   * </ul>
   */
  void randomize(const Ptr<MultiVectorBase<Scalar> > &mv)
    { randomizeImpl(mv); }

private:

  /** \brief . */
  virtual void randomizeImpl(const Ptr<MultiVectorBase<Scalar> > &mv) = 0;
  
};


} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_RANDOMIZER_BASE_HPP
