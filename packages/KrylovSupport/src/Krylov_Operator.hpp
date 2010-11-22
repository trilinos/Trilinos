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

#ifndef __Krylov_Operator_hpp
#define __Krylov_Operator_hpp

/// \file Krylov_Operator.hpp
///
/// Virtual base class which defines the operator interface required
/// by iterative linear solvers.

#include "Krylov_ConfigDefs.hpp"
#include "Krylov_OperatorTraits.hpp"
#include "Krylov_MultiVec.hpp"

/// \class Krylov::Operator
///
/// Templated pure virtual class for constructing the operator that is
/// used by the linear solver.
///
/// This operator is used as the interface to the matrix (<tt>A</tt>),
/// solution (<tt>X</tt>), and right-hand side (<tt>B</tt>) of the
/// linear system <tt>AX = B</tt>.  Furthermore, it is also the
/// interface to left/right preconditioning and left/right scaling of
/// the linear system.
///
/// A concrete implementation of this class is necessary.  The user
/// can create their own implementation if those supplied are not
/// suitable for their needs.
///
/// \author Original authors: Michael Heroux and Heidi Thornquist.
/// Moved to the KrylovSupport package and modified slightly by Mark
/// Hoemmen.

namespace Krylov {
  
  template< class ScalarType >
  class Operator {
  public:
    //! @name Constructor/Destructor
    //@{ 
    
    //! Default constructor
    Operator() {};
    
    //! Destructor.
    virtual ~Operator() {};
    //@}
    
    /// \fn Apply
    /// \brief y := A * x
    ///
    /// This method takes the Krylov::MultiVec x and applies the
    /// operator to it, resulting in the Belos::MultiVec y.
    ///
    /// \param x [in] Input multivector
    /// \param y [out] Output multivector
    ///
    /// \param trans [in] Whether to apply the operator (NOTRANS, the
    ///   default) or the (Hermitian) transpose of the operator.  If
    ///   the operator does not support the (Hermitian) transpose
    ///   operation, an exception is thrown.
    ///
    /// \note Implementations should report errors by throwing a
    ///   subclass of std::exception.
    virtual void 
    Apply (const MultiVec< ScalarType >& x, 
	   MultiVec< ScalarType >& y, 
	   const Teuchos::ETrans trans = NOTRANS) const = 0;
  };
  

  /// Template specialization of the Krylov::OperatorTraits class,
  /// using the Krylov::Operator and Krylov::MultiVec abstract base
  /// classes.
  /// 
  /// \note Any class that inherits from Krylov::Operator will be
  ///   accepted by the Anasazi and Belos templated solvers, due to
  ///   this interface to the Krylov::OperatorTraits class.
  template< class ScalarType > 
  class OperatorTraits< ScalarType, MultiVec< ScalarType >, Operator< ScalarType > > 
  {
  public:

    /// \fn Apply
    /// \brief y := A * x
    ///
    /// This class method takes the Krylov::MultiVec x and applies the
    /// operator to it, resulting in the Belos::MultiVec y.
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
    Apply (const Operator< ScalarType >& Op, 
	   const MultiVec< ScalarType >& x, 
	   MultiVec< ScalarType >& y,
	   const Teuchos::ETrans trans = NOTRANS)
    { 
      Op.Apply (x, y, trans); 
    }
    
  };
  
} // namespace Krylov

#endif // __Krylov_Operator_hpp

