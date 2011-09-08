// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
