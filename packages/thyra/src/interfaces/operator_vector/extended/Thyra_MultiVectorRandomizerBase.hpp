// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#ifndef THYRA_MULTI_VECTOR_RANDOMIZER_BASE_HPP
#define THYRA_MULTI_VECTOR_RANDOMIZER_BASE_HPP

#include "Thyra_OperatorVectorTypes.hpp"

namespace Thyra {

/** \brief Base interface for a strategy object for randomizing a multi-vector.
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
  virtual void randomize( MultiVectorBase<Scalar> *mv ) = 0;
  
};

} // namespace Thyra

#endif // THYRA_MULTI_VECTOR_RANDOMIZER_BASE_HPP
