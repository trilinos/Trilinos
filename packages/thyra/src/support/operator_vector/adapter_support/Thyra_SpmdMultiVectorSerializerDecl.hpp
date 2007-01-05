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

#ifndef THYRA_SPMD_MULTI_VECTOR_SERIALIZER_DECL_HPP
#define THYRA_SPMD_MULTI_VECTOR_SERIALIZER_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace Thyra {

/** \brief Concrete utility class for reading and writing SPMD-based
 * MultiVectorBase objects to and from standard streams.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_concrete_std_grp
 */
template<class Scalar>
class SpmdMultiVectorSerializer {
public:

  /// Set to true if to use binary IO and to false if using ASCII.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, binaryMode );

  /** \brief . */
  SpmdMultiVectorSerializer(
    const bool  binaryMode = false
    );

  /** \brief Determine if the multi-vector is compatible or not. */
  bool isCompatible( const MultiVectorBase<Scalar> &mv ) const;

  /** \brief Write to a stream.
   *
   * ToDo: Finish documentation!
   */
  void serialize( const MultiVectorBase<Scalar>& mv, std::ostream& out ) const;

  /** \brief Read from a stream
   *
   * ToDo: Finish documentation!
   */
  void deserialize( std::istream& in, MultiVectorBase<Scalar>* mv ) const;

}; // end class SpmdMultiVectorSerializer

} // end namespace Thyra

#endif // THYRA_SPMD_MULTI_VECTOR_SERIALIZER_DECL_HPP
