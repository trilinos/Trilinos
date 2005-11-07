// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef ANASAZI_LOCA_MATRIX_H
#define ANASAZI_LOCA_MATRIX_H

#include "Teuchos_RefCountPtr.hpp"

// Anasazi includes
#include "AnasaziOperator.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziReturnType.hpp"

// LOCA includes
#include "LOCA_AnasaziOperator_AbstractStrategy.H"

// forward declarations
namespace LOCA {
  class GlobalData;
}

namespace Anasazi {

  /*!
   * \brief Implementation of the Anasazi::Matrix for LOCA operators.
   */
  /*!
   * This class provides an implementation of the Anasazi::Matrix class
   * using LOCA operators.  The operator is represented by the 
   * LOCA::AnasaziOperator::AbstractStrategy class which serves as a base class
   * for many different types of operators, including shift-invert and
   * Generalized Cayley transformations.
   */
  class LOCAMatrix : public Anasazi::Operator<double> {
  public:

    //! Constructor
    LOCAMatrix(
    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
    const Teuchos::RefCountPtr< ::LOCA::AnasaziOperator::AbstractStrategy>& op);
    
    //! Destructor
    ~LOCAMatrix();
    
    //! Apply operator
    Anasazi::ReturnType Apply(const Anasazi::MultiVec<double>& x, 
			      Anasazi::MultiVec<double>& y ) const;
  protected:

    //! Global data
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData;
    
    //! LOCA operator
    Teuchos::RefCountPtr< ::LOCA::AnasaziOperator::AbstractStrategy> locaOp;
    
  }; // class Matrix
} // namespace Anasazi

#endif 
