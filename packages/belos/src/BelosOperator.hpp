//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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

#ifndef BELOS_OPERATOR_HPP
#define BELOS_OPERATOR_HPP

/*!     \file BelosOperator.hpp
        \brief Virtual base class which defines the operator interface 
	required by the iterative linear solver.
*/

#include "BelosOperatorTraits.hpp"
#include "BelosMultiVec.hpp"
#include "BelosConfigDefs.hpp"

/*!	\class Belos::Operator

	\brief Belos's templated pure virtual class for constructing the operator that is
	used by the linear solver.  

	This operator is used as the interface to the matrix (<tt>A</tt>), 
	solution (<tt>X</tt>), and right-hand side (<tt>B</tt>) of the linear system <tt>AX = B</tt>.
	Furthermore, it is also the interface to left/right preconditioning and left/right scaling of the
	linear system.

	A concrete implementation of this class is necessary.  The user can create their own implementation
	if those supplied are not suitable for their needs.

	\author Michael Heroux and Heidi Thornquist
*/

namespace Belos {
  
  template <class ScalarType>
  class Operator {
  public:
    
    //! @name Constructor/Destructor
    //@{ 
    
    //! Default constructor
    Operator() {};
    
    //! Destructor.
    virtual ~Operator() {};
    //@}
    
    //! @name Operator application method
    //@{ 
    
    /*! \brief This routine takes the Belos::MultiVec \c x and applies the operator
      to it resulting in the Belos::MultiVec \c y, which is returned.
        \note It is expected that any problem with applying this operator to \c x will be
	indicated by an std::exception being thrown.
    */
    virtual void Apply ( const MultiVec<ScalarType>& x, 
			 MultiVec<ScalarType>& y, ETrans trans=NOTRANS ) const = 0;
  };
  
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for Belos::Operator 
  //                                               and Belos::MultiVec.
  //
  ////////////////////////////////////////////////////////////////////  
  
  /*!  \brief Template specialization of Belos::OperatorTraits class using Belos::Operator and Belos::MultiVec virtual
    base classes.
    
    Any class that inherits from Belos::Operator will be accepted by the Belos templated solvers due to this
    interface to the Belos::OperatorTraits class.
  */

  template <class ScalarType> 
  class OperatorTraits < ScalarType, MultiVec<ScalarType>, Operator<ScalarType> > 
  {
  public:
    
    ///
    static void Apply ( const Operator<ScalarType>& Op, 
			const MultiVec<ScalarType>& x, 
			MultiVec<ScalarType>& y,
			ETrans trans=NOTRANS )
    { Op.Apply( x, y, trans ); }
    
  };
  
} // end Belos namespace

#endif

// end of file BelosOperator.hpp
