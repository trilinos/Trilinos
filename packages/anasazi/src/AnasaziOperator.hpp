// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
