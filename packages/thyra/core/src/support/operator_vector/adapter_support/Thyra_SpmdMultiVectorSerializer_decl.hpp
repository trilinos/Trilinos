// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
