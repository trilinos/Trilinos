//@HEADER
// ************************************************************************
//
//             KrylovSupport: Krylov Methods Support Package
//                  Copyright 2010 Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef __Krylov_OperatorTraits_hpp
#define __Krylov_OperatorTraits_hpp

/// \file Krylov_OperatorTraits.hpp
/// \brief Traits class for operators in Trilinos iterative methods
///

#include "Krylov_ConfigDefs.hpp"

namespace Krylov {

  template< class ScalarType, class MV, class OP >
  struct UndefinedOperatorTraits
  {
    /// This function should not compile if there is an attempt to
    /// instantiate!
    ///
    /// \note Any attempt to compile this function results in a
    ///   compile time error.  This means that the template
    ///   specialization of Krylov::OperatorTraits class does not
    ///   exist for type <tt>OP</tt>, or is not complete.
    static inline void notDefined() { OP::this_type_is_missing_a_specialization(); };
  };
 
  /// \class OperatorTraits
  /// \brief Traits class for operators in Trilinos iterative methods
  ///
  /// \note An adapter for this traits class must exist for the
  ///   <tt>MV</tt> and <tt>OP</tt> types.  If not, this class will
  ///   produce a compile-time error.
  ///
  template< class ScalarType, class MV, class OP >
  class OperatorTraits 
  {
  public:
    
    /// \fn Apply
    /// \brief y := A * x
    ///
    /// This class method takes the MV x and applies the operator to
    /// it, resulting in the MV y.
    ///
    /// \param x [in] Input multivector
    /// \param y [out] Output multivector
    ///
    /// \param trans [in] Whether to apply the operator (NOTRANS, the
    ///   default) or the (Hermitian) transpose of the operator.  If
    ///   the operator does not support the (Hermitian) transpose
    ///   operation, an exception is thrown.
    ///
    static void 
    Apply (const OP& Op,
	   const MV& x,
	   MV& y,
	   const Teuchos::ETrans trans = NOTRANS)
    { 
      // This won't compile unless OperatorTraits is specialized for
      // ScalarType, MV, and OP.
      UndefinedOperatorTraits<ScalarType, MV, OP>::notDefined(); 
    }
  };
  
} // namespace Krylov

#endif // __Krylov_OperatorTraits_hpp
