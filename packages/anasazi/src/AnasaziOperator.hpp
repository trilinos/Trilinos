// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// @HEADER

/*! \file AnasaziOperator.hpp
  \brief  Templated virtual class for creating operators that can interface with the Anasazi::OperatorTraits class
*/

#ifndef ANASAZI_OPERATOR_HPP
#define ANASAZI_OPERATOR_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziMultiVec.hpp"
#include "Teuchos_ScalarTraits.hpp"


namespace Anasazi {
  
/*!	
	\brief Anasazi's templated virtual class for constructing an operator that can interface with the 
	OperatorTraits class used by the eigensolvers.
	
	A concrete implementation of this class is necessary.  The user can create their own implementation
	if those supplied are not suitable for their needs.

	\author Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/
  template <class ScalarType>
  class Operator {
  public:
    //! @name Constructor/Destructor
    //@{ 
    //! Default constructor.
    Operator() {};
    
    //! Destructor.
    virtual ~Operator() {};
    //@}
    
    //! @name Operator application method
    //@{ 
    
    /*! \brief This method takes the Anasazi::MultiVec \c x and
      applies the operator to it resulting in the Anasazi::MultiVec \c y.
    */
    virtual void Apply ( const MultiVec<ScalarType>& x, MultiVec<ScalarType>& y ) const = 0;

    //@}
  };
  
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::OperatorTraits for Anasazi::Operator 
  //                                               and Anasazi::MultiVec.
  //
  ////////////////////////////////////////////////////////////////////  
  
  /*! 
    \brief Template specialization of Anasazi::OperatorTraits class using Anasazi::Operator and Anasazi::MultiVec virtual
    base classes.

    Any class that inherits from Anasazi::Operator will be accepted by the Anasazi templated solvers due to this
    interface to the Anasazi::OperatorTraits class.
  */

  template <class ScalarType> 
  class OperatorTraits < ScalarType, MultiVec<ScalarType>, Operator<ScalarType> > 
  {
  public:
  
    //! @name Operator application method
    //@{ 

    /*! \brief This method takes the Anasazi::MultiVec \c x and
      applies the Anasazi::Operator \c Op to it resulting in the Anasazi::MultiVec \c y.
    */
    static void Apply ( const Operator<ScalarType>& Op, 
			      const MultiVec<ScalarType>& x, 
			      MultiVec<ScalarType>& y )
    { Op.Apply( x, y ); }

    //@}
    
  };
  
} // end of Anasazi namespace

#endif

// end of file AnasaziOperator.hpp
